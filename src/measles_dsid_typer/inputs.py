"""Input handling: FASTA, FASTQ (consensus), and BAM (consensus).

Design goals:
  * FASTA works with zero third-party dependencies (stdlib + our parser).
  * FASTQ consensus works dependency-free via a simple N-450-anchored pileup:
    every read is anchored to the N-450 frame using the same reference k-mer
    seeds as assembly input, and a per-position majority consensus (with a
    configurable minimum depth) is called. This avoids requiring an aligner for
    the common case of amplicon N-450 reads.
  * BAM consensus uses ``pysam`` when available (it is the correct tool for
    aligned data). If ``pysam`` is not installed we raise a clear, actionable
    error rather than guessing.

All record yields are (sample_id, sequence) tuples.
"""

from __future__ import annotations

import gzip
import io
import os
from collections import Counter
from typing import Dict, Iterator, List, Optional, Tuple

from . import sequtils
from .n450 import SEED_K, _SeedIndex
from .reference import Reference, N450_LENGTH


def _open_text(path: str):
    """Open plain or gzip-compressed text transparently."""
    if path.endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "rb"))
    return open(path, "r")


def detect_format(path: str) -> str:
    """Guess input format from extension. Returns 'fasta'|'fastq'|'bam'."""
    low = path.lower()
    base = low[:-3] if low.endswith(".gz") else low
    if base.endswith((".bam", ".sam", ".cram")):
        return "bam"
    if base.endswith((".fastq", ".fq")):
        return "fastq"
    if base.endswith((".fasta", ".fa", ".fna", ".fas")):
        return "fasta"
    # fall back to sniffing the first non-empty byte
    try:
        with _open_text(path) as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    return "fasta"
                if line.startswith("@"):
                    return "fastq"
                break
    except Exception:
        pass
    return "fasta"


# --------------------------------------------------------------------------- #
# FASTA
# --------------------------------------------------------------------------- #
def iter_fasta_file(path: str) -> Iterator[Tuple[str, str]]:
    with _open_text(path) as fh:
        for header, seq in sequtils.iter_fasta(fh):
            sample = header.split()[0] if header else header
            yield sample, seq


# --------------------------------------------------------------------------- #
# FASTQ (dependency-free N-450 amplicon consensus)
# --------------------------------------------------------------------------- #
def iter_fastq_reads(path: str) -> Iterator[Tuple[str, str, str]]:
    """Yield (read_id, sequence, qualities) from a FASTQ file."""
    with _open_text(path) as fh:
        while True:
            h = fh.readline()
            if not h:
                break
            seq = fh.readline().rstrip("\n").rstrip("\r")
            plus = fh.readline()
            qual = fh.readline().rstrip("\n").rstrip("\r")
            if not qual and not seq:
                break
            read_id = h[1:].strip().split()[0] if h.startswith("@") else h.strip()
            yield read_id, seq, qual


def fastq_consensus(
    path: str,
    reference: Reference,
    seed_index: Optional[_SeedIndex] = None,
    min_depth: int = 3,
    min_base_quality: int = 20,
    min_reads: int = 1,
) -> Tuple[str, Dict[str, int]]:
    """Build an N-450 consensus from amplicon FASTQ reads.

    Each read is anchored to the 450-nt N-450 frame via reference k-mer seeds
    (forward and reverse-complement are both tried). A per-position base count
    is accumulated, filtering bases below ``min_base_quality`` (Phred+33). The
    consensus calls the majority base where depth >= ``min_depth`` and 'N'
    elsewhere.

    Returns (consensus_sequence, stats).
    """
    if seed_index is None:
        seed_index = _SeedIndex(reference)

    # columns[pos] -> Counter of bases
    columns: List[Counter] = [Counter() for _ in range(N450_LENGTH)]
    n_reads = 0
    n_anchored = 0

    for _rid, seq, qual in iter_fastq_reads(path):
        n_reads += 1
        norm = sequtils.normalize(seq)
        norm = "".join(c for c in norm if c not in sequtils.GAP_CHARS)
        if len(norm) < SEED_K:
            continue

        best = None  # (oriented_seq, oriented_qual, start, hits)
        q_fwd = qual
        q_rev = qual[::-1] if qual else ""
        for oriented, oq, _name in (
            (norm, q_fwd, "fwd"),
            (sequtils.reverse_complement(norm), q_rev, "rev"),
        ):
            loc = seed_index.locate(oriented)
            if loc is None:
                continue
            start, hits = loc
            if best is None or hits > best[3]:
                best = (oriented, oq, start, hits)
        if best is None:
            continue
        n_anchored += 1
        oriented, oq, start, _hits = best

        # ``start`` from seed_index.locate() is the inferred start of the full
        # 450-nt N-450 frame *within the read*. For a short read that is itself
        # a sub-window of the frame this is negative; the read's own index ``i``
        # therefore maps to N-450 coordinate (i - start).
        for i, base in enumerate(oriented):
            pos = i - start  # N-450 coordinate
            if pos < 0 or pos >= N450_LENGTH:
                continue
            if base not in "ACGT":
                continue
            if oq and i < len(oq):
                if (ord(oq[i]) - 33) < min_base_quality:
                    continue
            columns[pos][base] += 1

    consensus_chars: List[str] = []
    covered = 0
    for col in columns:
        depth = sum(col.values())
        if depth >= min_depth:
            # majority; deterministic tie-break on base order A<C<G<T
            top = max(col.items(), key=lambda kv: (kv[1], "TGCA".index(kv[0])))
            consensus_chars.append(top[0])
            covered += 1
        else:
            consensus_chars.append("N")
    consensus = "".join(consensus_chars)

    stats = {
        "reads_total": n_reads,
        "reads_anchored": n_anchored,
        "positions_covered": covered,
        "positions_total": N450_LENGTH,
    }
    if n_anchored < min_reads:
        stats["warning"] = 1  # too few anchored reads for a trustworthy call
    return consensus, stats


# --------------------------------------------------------------------------- #
# BAM (consensus via pysam)
# --------------------------------------------------------------------------- #
def bam_consensus(
    path: str,
    min_depth: int = 3,
    min_base_quality: int = 20,
    region: Optional[str] = None,
) -> Tuple[str, Dict[str, int]]:
    """Build a consensus from an aligned BAM/CRAM using pysam pileup.

    Requires ``pysam``. The consensus spans the covered reference positions;
    the N-450 window is then located downstream by the typer's extractor, so
    this returns the full covered consensus rather than a hard 450-nt slice.
    """
    try:
        import pysam  # type: ignore
    except ImportError as exc:  # pragma: no cover - environment dependent
        raise RuntimeError(
            "BAM input requires the optional dependency 'pysam'. "
            "Install with `pip install measles-dsid-typer[bam]` or "
            "`conda install -c bioconda pysam`."
        ) from exc

    af = pysam.AlignmentFile(path, "rb")
    # map reference position -> Counter of bases
    columns: Dict[int, Counter] = {}
    for col in af.pileup(
        region=region, truncate=bool(region), min_base_quality=min_base_quality
    ):
        pos = col.reference_pos
        counter = columns.setdefault(pos, Counter())
        for pread in col.pileups:
            if pread.is_del or pread.is_refskip:
                continue
            qpos = pread.query_position
            if qpos is None:
                continue
            base = pread.alignment.query_sequence[qpos].upper()
            if base in "ACGT":
                counter[base] += 1
    af.close()

    if not columns:
        return "", {"positions_covered": 0}

    lo = min(columns)
    hi = max(columns)
    chars: List[str] = []
    covered = 0
    for pos in range(lo, hi + 1):
        col = columns.get(pos)
        depth = sum(col.values()) if col else 0
        if col and depth >= min_depth:
            top = max(col.items(), key=lambda kv: (kv[1], "TGCA".index(kv[0])))
            chars.append(top[0])
            covered += 1
        else:
            chars.append("N")
    consensus = "".join(chars)
    stats = {
        "positions_covered": covered,
        "span": hi - lo + 1,
        "ref_start": lo,
        "ref_end": hi,
    }
    return consensus, stats
