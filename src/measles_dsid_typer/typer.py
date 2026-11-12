"""Core DSID typing logic.

Given an input sequence, the typer:
  1. extracts the 450-nt N-450 window (see n450.py),
  2. looks for a byte-exact reference match (coordinate-exact N-450 identity),
  3. otherwise finds the nearest reference allele by IUPAC-aware mismatch count,
  4. reports genotype, DSID sub-lineage, mismatch count, % identity and the
     ambiguity/IUPAC handling that was applied.

Determinism: ties between equally-good reference alleles are broken by a stable
sort on (mismatches, dsid, accession), so the same input always yields the same
assignment regardless of reference file ordering.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, Tuple

from . import sequtils
from .n450 import N450Extraction, N450ExtractionError, _SeedIndex, extract_n450
from .reference import Reference, ReferenceAllele, N450_LENGTH


@dataclass
class TyperResult:
    """The result of typing a single input sequence."""

    sample: str
    status: str                       # "typed", "no_n450", "no_reference"
    genotype: Optional[str] = None
    dsid: Optional[str] = None
    collision_dsids: Tuple[str, ...] = field(default_factory=tuple)
    match_type: Optional[str] = None  # "exact" | "nearest"
    n450_coordinate_match: bool = False   # True iff byte-exact 450-nt identity
    mismatches: Optional[int] = None
    pident: Optional[float] = None    # percent identity over called positions
    aligned_length: Optional[int] = None
    ambiguous_positions: int = 0      # IUPAC/N positions handled in the query
    tied_dsids: Tuple[str, ...] = field(default_factory=tuple)  # all DSIDs at min mm
    dsid_ambiguous: bool = False      # True iff >1 DSID tied (ambiguous call)
    n450_length: Optional[int] = None
    orientation: Optional[str] = None
    extraction_method: Optional[str] = None
    best_reference: Optional[str] = None  # reference allele label
    n450_sequence: Optional[str] = None
    note: Optional[str] = None

    def as_row(self, include_sequence: bool = False) -> dict:
        """Flat, deterministic dict for CSV/TSV/JSON emission."""
        row = {
            "sample": self.sample,
            "status": self.status,
            "genotype": self.genotype or "",
            "dsid": self.dsid or "",
            "collision_dsids": ";".join(self.collision_dsids),
            "match_type": self.match_type or "",
            "n450_coordinate_match": self.n450_coordinate_match,
            "mismatches": self.mismatches if self.mismatches is not None else "",
            "pident": f"{self.pident:.3f}" if self.pident is not None else "",
            "aligned_length": self.aligned_length
            if self.aligned_length is not None
            else "",
            "ambiguous_positions": self.ambiguous_positions,
            "tied_dsids": ";".join(self.tied_dsids),
            "dsid_ambiguous": self.dsid_ambiguous,
            "n450_length": self.n450_length if self.n450_length is not None else "",
            "orientation": self.orientation or "",
            "extraction_method": self.extraction_method or "",
            "best_reference": self.best_reference or "",
            "note": self.note or "",
        }
        if include_sequence:
            row["n450_sequence"] = self.n450_sequence or ""
        return row


# Ordered field list so CSV/TSV/JSON columns are stable across runs/platforms.
RESULT_FIELDS: Tuple[str, ...] = (
    "sample",
    "status",
    "genotype",
    "dsid",
    "collision_dsids",
    "match_type",
    "n450_coordinate_match",
    "mismatches",
    "pident",
    "aligned_length",
    "ambiguous_positions",
    "tied_dsids",
    "dsid_ambiguous",
    "n450_length",
    "orientation",
    "extraction_method",
    "best_reference",
    "note",
)


def _best_allele(
    n450: str, reference: Reference
) -> Tuple[ReferenceAllele, int, int, int, Tuple[str, ...]]:
    """Nearest reference allele by IUPAC-aware mismatch count.

    Returns (allele, mismatches, ambiguous_positions, called_positions,
    tied_dsids) where ``tied_dsids`` is every distinct DSID achieving the
    minimal mismatch count (including the chosen one). Ambiguous (IUPAC/N) query
    positions can make several DSIDs equally compatible; surfacing the tie set
    keeps that ambiguity visible instead of hiding it behind one pick.

    Deterministic tie-break for the single chosen allele: fewer mismatches,
    then lexicographic dsid, then accession.
    """
    best: Optional[Tuple[int, str, str, ReferenceAllele, int, int]] = None
    min_mm: Optional[int] = None
    tied: set = set()
    for allele in reference.alleles:
        mm, amb, called = sequtils.count_mismatches(n450, allele.sequence)
        if min_mm is None or mm < min_mm:
            min_mm = mm
            tied = {allele.dsid}
        elif mm == min_mm:
            tied.add(allele.dsid)
        key = (mm, allele.dsid, allele.accession)
        if best is None or key < (best[0], best[1], best[2]):
            best = (mm, allele.dsid, allele.accession, allele, amb, called)
    assert best is not None
    return best[3], best[0], best[4], best[5], tuple(sorted(tied))


def type_sequence(
    sample: str,
    sequence: str,
    reference: Reference,
    seed_index: Optional[_SeedIndex] = None,
    keep_sequence: bool = False,
) -> TyperResult:
    """Type a single already-assembled nucleotide sequence."""
    try:
        extraction: N450Extraction = extract_n450(sequence, reference, seed_index)
    except N450ExtractionError as exc:
        return TyperResult(
            sample=sample,
            status="no_n450",
            note=str(exc),
        )

    n450 = extraction.sequence

    # Always compute the nearest allele + full tie set (IUPAC-aware). This
    # surfaces DSID ambiguity caused by N / IUPAC positions in the query.
    nearest_allele, nearest_mm, amb, called, tied = _best_allele(n450, reference)
    dsid_ambiguous = len(tied) > 1

    # 1) exact coordinate match (byte-identical 450-nt N-450)
    exact = reference.exact_lookup(n450)
    if exact:
        # deterministic pick among identical references
        exact_sorted = sorted(exact, key=lambda a: (a.dsid, a.accession))
        allele = exact_sorted[0]
        # union of all DSIDs covered by identical references + recorded collisions
        coll = set()
        for a in exact:
            coll.update(a.collision_dsids)
            coll.add(a.dsid)
        coll.discard(allele.dsid)
        return TyperResult(
            sample=sample,
            status="typed",
            genotype=allele.genotype,
            dsid=allele.dsid,
            collision_dsids=tuple(sorted(coll)),
            match_type="exact",
            n450_coordinate_match=True,
            mismatches=0,
            pident=100.0,
            aligned_length=N450_LENGTH,
            ambiguous_positions=sum(1 for c in n450 if sequtils.is_ambiguous(c)),
            tied_dsids=tied,
            dsid_ambiguous=dsid_ambiguous,
            n450_length=len(n450),
            orientation=extraction.orientation,
            extraction_method=extraction.method,
            best_reference=allele.label,
            n450_sequence=n450 if keep_sequence else None,
        )

    # 2) nearest neighbour by IUPAC-aware mismatch count
    allele, mm = nearest_allele, nearest_mm
    pident = 100.0 * (called - mm) / called if called else 0.0
    coll = set(allele.collision_dsids)
    coll.discard(allele.dsid)
    return TyperResult(
        sample=sample,
        status="typed",
        genotype=allele.genotype,
        dsid=allele.dsid,
        collision_dsids=tuple(sorted(coll)),
        match_type="nearest",
        n450_coordinate_match=False,
        mismatches=mm,
        pident=pident,
        aligned_length=called,
        ambiguous_positions=amb,
        tied_dsids=tied,
        dsid_ambiguous=dsid_ambiguous,
        n450_length=len(n450),
        orientation=extraction.orientation,
        extraction_method=extraction.method,
        best_reference=allele.label,
        n450_sequence=n450 if keep_sequence else None,
    )


def type_record(
    sample: str,
    sequence: str,
    reference: Reference,
    seed_index: Optional[_SeedIndex] = None,
    keep_sequence: bool = False,
) -> TyperResult:
    """Alias kept for API clarity when iterating over parsed records."""
    return type_sequence(
        sample, sequence, reference, seed_index, keep_sequence
    )
