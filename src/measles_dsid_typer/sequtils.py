"""Low-level nucleotide sequence utilities with IUPAC ambiguity handling.

All functions here are pure and deterministic. The IUPAC logic is the single
source of truth for how ambiguous bases are compared throughout the typer.
"""

from __future__ import annotations

from typing import Dict, FrozenSet, Iterable, Tuple

# IUPAC nucleotide ambiguity codes -> the set of concrete bases they represent.
IUPAC_CODES: Dict[str, FrozenSet[str]] = {
    "A": frozenset("A"),
    "C": frozenset("C"),
    "G": frozenset("G"),
    "T": frozenset("T"),
    "U": frozenset("T"),  # RNA uracil treated as thymine
    "R": frozenset("AG"),
    "Y": frozenset("CT"),
    "S": frozenset("GC"),
    "W": frozenset("AT"),
    "K": frozenset("GT"),
    "M": frozenset("AC"),
    "B": frozenset("CGT"),
    "D": frozenset("AGT"),
    "H": frozenset("ACT"),
    "V": frozenset("ACG"),
    "N": frozenset("ACGT"),
}

# Bases considered "ambiguous" (anything expanding to >1 concrete base, plus N).
AMBIGUOUS_BASES: FrozenSet[str] = frozenset(
    code for code, bases in IUPAC_CODES.items() if len(bases) > 1
)

# Gap / missing-data characters.
GAP_CHARS: FrozenSet[str] = frozenset("-.")

_COMPLEMENT = str.maketrans(
    "ACGTUNRYSWKMBDHVacgtunryswkmbdhv",
    "TGCAANYRSWMKVHDBtgcaanyrswmkvhdb",
)


def normalize(seq: str) -> str:
    """Upper-case a sequence and map U->T. Non-letter whitespace is stripped."""
    return seq.strip().upper().replace("U", "T")


def reverse_complement(seq: str) -> str:
    """Reverse complement, preserving IUPAC codes."""
    return seq.translate(_COMPLEMENT)[::-1]


def is_ambiguous(base: str) -> bool:
    return base in AMBIGUOUS_BASES


def bases_compatible(a: str, b: str) -> bool:
    """True if two (possibly ambiguous) bases could represent the same base.

    Two IUPAC codes are compatible when the base sets they expand to share at
    least one concrete nucleotide. Gaps are never compatible with anything.
    This is the IUPAC-aware definition of a non-mismatch.
    """
    if a in GAP_CHARS or b in GAP_CHARS:
        return False
    sa = IUPAC_CODES.get(a)
    sb = IUPAC_CODES.get(b)
    if sa is None or sb is None:
        return False
    return not sa.isdisjoint(sb)


def count_mismatches(query: str, reference: str) -> Tuple[int, int, int]:
    """Compare two equal-length sequences position by position.

    Returns ``(mismatches, ambiguous_positions, called_positions)`` where:
      * ``mismatches``  - positions whose base sets are *incompatible*
                          (a true substitution). An ambiguous query base that
                          is still compatible with the reference is NOT a
                          mismatch.
      * ``ambiguous_positions`` - query positions that carried an IUPAC
                          ambiguity code (incl. N) but were handled.
      * ``called_positions`` - positions where the query had a real base
                          (not a gap / missing); the denominator for identity.

    Raises ValueError if lengths differ (callers must align/trim first).
    """
    if len(query) != len(reference):
        raise ValueError(
            f"length mismatch: query={len(query)} reference={len(reference)}"
        )
    mismatches = 0
    ambiguous = 0
    called = 0
    for q, r in zip(query, reference):
        if q in GAP_CHARS:
            # missing query base: not called, not a mismatch
            continue
        called += 1
        if is_ambiguous(q):
            ambiguous += 1
        if not bases_compatible(q, r):
            mismatches += 1
    return mismatches, ambiguous, called


def kmer_set(seq: str, k: int) -> FrozenSet[str]:
    """All length-k substrings of seq that contain only unambiguous ACGT."""
    out = set()
    n = len(seq)
    for i in range(n - k + 1):
        sub = seq[i : i + k]
        if all(c in "ACGT" for c in sub):
            out.add(sub)
    return frozenset(out)


def iter_fasta(lines: Iterable[str]):
    """Minimal streaming FASTA parser. Yields (header, sequence) tuples.

    Header is the text after '>' up to the first whitespace is kept as the id,
    but the full header line (minus '>') is returned so callers can parse it.
    """
    header = None
    chunks = []
    for raw in lines:
        line = raw.rstrip("\n").rstrip("\r")
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(chunks)
            header = line[1:].strip()
            chunks = []
        else:
            chunks.append(line.strip())
    if header is not None:
        yield header, "".join(chunks)
