"""Extraction of the 450-nt N-450 genotyping window from an input sequence.

Inputs to the typer can be:
  * already-trimmed N-450 fragments (exactly or ~450 nt), or
  * longer sequences (full N gene, or whole-genome assemblies ~15.9 kb)

To stay dependency-free and deterministic we anchor the N-450 window using
k-mer seeds drawn from the curated reference alleles, then take a fixed-length
450-nt window around the best-supported offset. Orientation (forward vs reverse
complement) is resolved automatically.

The approach is coordinate-anchored to the reference allele set itself rather
than to an external WHO genotype-A genome, so no network access or extra files
are needed. For genotypes present in the reference (B3, D8, and any added
later) this reliably recovers the same 450-nt frame the references use.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

from . import sequtils
from .reference import Reference, N450_LENGTH

# k-mer size used to seed the N-450 window. 21 is long enough to be specific
# within the ~15.9 kb measles genome yet short enough to tolerate SNPs.
SEED_K = 21
# minimum number of distinct anchoring seed hits to trust an extraction
MIN_SEED_HITS = 3


class N450ExtractionError(Exception):
    """Raised when a 450-nt N-450 window cannot be confidently located."""


@dataclass(frozen=True)
class N450Extraction:
    sequence: str          # the 450-nt window (normalized, forward orientation)
    orientation: str       # "forward" or "reverse_complement"
    source_length: int     # length of the input sequence
    method: str            # "verbatim", "anchored"
    offset: Optional[int]  # start offset in the (oriented) source, if anchored


class _SeedIndex:
    """Maps reference k-mers -> the offset they occur at within a 450-nt frame.

    Because every reference allele is itself a 450-nt N-450 frame, a reference
    k-mer observed at position ``p`` within an allele implies that, in a query,
    the same k-mer marks N-450 offset ``p``. We take the modal offset across
    all seed hits to place the window.
    """

    def __init__(self, reference: Reference, k: int = SEED_K):
        self.k = k
        self._offsets: Dict[str, int] = {}
        for allele in reference.alleles:
            s = allele.sequence
            for i in range(len(s) - k + 1):
                sub = s[i : i + k]
                if all(c in "ACGT" for c in sub):
                    # keep first-seen offset; identical k-mers across alleles
                    # share the same N-450 coordinate in practice
                    self._offsets.setdefault(sub, i)

    def locate(self, seq: str) -> Optional[Tuple[int, int]]:
        """Return (window_start, n_hits) for the modal N-450 placement, or None.

        window_start is the inferred start of the 450-nt frame within ``seq``.
        """
        k = self.k
        votes: Dict[int, int] = {}
        for i in range(len(seq) - k + 1):
            sub = seq[i : i + k]
            ref_off = self._offsets.get(sub)
            if ref_off is None:
                continue
            start = i - ref_off
            votes[start] = votes.get(start, 0) + 1
        if not votes:
            return None
        # modal start; tie-break on smallest start for determinism
        best_start = max(sorted(votes), key=lambda s: votes[s])
        return best_start, votes[best_start]


def extract_n450(
    seq: str,
    reference: Reference,
    seed_index: Optional[_SeedIndex] = None,
) -> N450Extraction:
    """Locate and return the 450-nt N-450 window from ``seq``.

    Fast path: if the input is already exactly 450 nt we take it verbatim
    (still resolving orientation against the reference if needed).

    Otherwise we anchor via reference k-mer seeds in both orientations and
    carve out the 450-nt frame.
    """
    norm = sequtils.normalize(seq)
    norm = "".join(c for c in norm if c not in sequtils.GAP_CHARS)
    if not norm:
        raise N450ExtractionError("empty sequence after normalization")

    if seed_index is None:
        seed_index = _SeedIndex(reference)

    # ---- verbatim fast path -------------------------------------------- #
    if len(norm) == N450_LENGTH:
        # confirm orientation: does forward anchor at offset 0 better than RC?
        fwd = seed_index.locate(norm)
        rc = sequtils.reverse_complement(norm)
        rc_hit = seed_index.locate(rc)
        fwd_hits = fwd[1] if fwd else 0
        rc_hits = rc_hit[1] if rc_hit else 0
        if rc_hits > fwd_hits:
            return N450Extraction(rc, "reverse_complement", len(norm), "verbatim", 0)
        return N450Extraction(norm, "forward", len(norm), "verbatim", 0)

    # ---- anchored path ------------------------------------------------- #
    best: Optional[Tuple[str, str, int, int]] = None  # (window, orient, start, hits)
    for oriented, orient_name in (
        (norm, "forward"),
        (sequtils.reverse_complement(norm), "reverse_complement"),
    ):
        loc = seed_index.locate(oriented)
        if loc is None:
            continue
        start, hits = loc
        if best is None or hits > best[3]:
            best = (oriented, orient_name, start, hits)

    if best is None or best[3] < MIN_SEED_HITS:
        hits = best[3] if best else 0
        raise N450ExtractionError(
            f"could not anchor N-450 window (max seed hits={hits}, "
            f"need >= {MIN_SEED_HITS}); input length {len(norm)}"
        )

    oriented, orient_name, start, hits = best
    # clamp the window to the available sequence
    lo = max(0, start)
    hi = lo + N450_LENGTH
    window = oriented[lo:hi]
    if len(window) != N450_LENGTH:
        # try to recover by shifting left if we ran off the 3' end
        if hi > len(oriented) and len(oriented) >= N450_LENGTH:
            window = oriented[len(oriented) - N450_LENGTH : len(oriented)]
            lo = len(oriented) - N450_LENGTH
        else:
            raise N450ExtractionError(
                f"anchored window is {len(window)} nt, expected {N450_LENGTH} "
                f"(input too short to contain a full N-450)"
            )
    return N450Extraction(window, orient_name, len(norm), "anchored", lo)
