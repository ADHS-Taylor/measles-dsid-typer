"""Loading and indexing of the curated WHO/MeaNS N-450 DSID reference alleles.

Reference FASTA headers follow the MeaNS-derived convention:

    >DSId_<GENOTYPE>-<DSID>_<ACCESSION>
    e.g. >DSId_D8-9171_PZ023615.1

Some reference alleles are byte-identical across several DSIDs; the
``dsid_collisions.csv`` file records, for each kept reference, every DSID that
the identical sequence is known to represent. We surface those as
``collision_dsids`` so downstream users can see the ambiguity instead of
silently collapsing it.
"""

from __future__ import annotations

import csv
import os
import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from . import sequtils

_HEADER_RE = re.compile(
    r"^DSId_(?P<genotype>[A-Za-z0-9]+)-(?P<dsid>[0-9]+)_(?P<accession>\S+)$"
)

# Reference files shipped inside the package data dir.
_PKG_DIR = os.path.dirname(os.path.abspath(__file__))
_DATA_DIR = os.path.join(_PKG_DIR, "data")
DEFAULT_REFERENCE_FASTA = os.path.join(_DATA_DIR, "dsid_n450_reference.fasta")
DEFAULT_COLLISIONS_CSV = os.path.join(_DATA_DIR, "dsid_collisions.csv")

N450_LENGTH = 450


@dataclass(frozen=True)
class ReferenceAllele:
    """A single curated N-450 reference allele."""

    label: str            # full header, e.g. DSId_D8-9171_PZ023615.1
    genotype: str         # e.g. D8
    dsid: str             # full DSID label incl. genotype, e.g. D8-9171
    accession: str        # e.g. PZ023615.1
    sequence: str         # normalized 450-nt N-450 sequence
    collision_dsids: Tuple[str, ...] = field(default_factory=tuple)


@dataclass
class Reference:
    """An indexed collection of N-450 reference alleles."""

    alleles: List[ReferenceAllele]
    # exact-sequence index for O(1) exact matching
    _by_sequence: Dict[str, List[ReferenceAllele]] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not self._by_sequence:
            idx: Dict[str, List[ReferenceAllele]] = {}
            for a in self.alleles:
                idx.setdefault(a.sequence, []).append(a)
            self._by_sequence = idx

    # ------------------------------------------------------------------ #
    def exact_lookup(self, seq: str) -> List[ReferenceAllele]:
        """Return alleles whose sequence is byte-identical to ``seq``."""
        return list(self._by_sequence.get(seq, ()))

    @property
    def n_alleles(self) -> int:
        return len(self.alleles)

    @property
    def genotypes(self) -> Tuple[str, ...]:
        return tuple(sorted({a.genotype for a in self.alleles}))

    @property
    def dsids(self) -> Tuple[str, ...]:
        return tuple(sorted({a.dsid for a in self.alleles}))

    # ------------------------------------------------------------------ #
    @classmethod
    def load(
        cls,
        fasta_path: Optional[str] = None,
        collisions_path: Optional[str] = None,
    ) -> "Reference":
        fasta_path = fasta_path or DEFAULT_REFERENCE_FASTA
        collisions_path = collisions_path or DEFAULT_COLLISIONS_CSV

        collisions = _load_collisions(collisions_path)

        alleles: List[ReferenceAllele] = []
        with open(fasta_path, "r") as fh:
            for header, seq in sequtils.iter_fasta(fh):
                m = _HEADER_RE.match(header)
                if not m:
                    # Skip headers that do not match the DSID convention rather
                    # than guessing; a malformed reference should be loud.
                    raise ValueError(
                        f"reference header does not match DSId convention: {header!r}"
                    )
                genotype = m.group("genotype")
                dsid = f"{genotype}-{m.group('dsid')}"
                accession = m.group("accession")
                norm = sequtils.normalize(seq)
                coll = tuple(collisions.get(header, ()))
                alleles.append(
                    ReferenceAllele(
                        label=header,
                        genotype=genotype,
                        dsid=dsid,
                        accession=accession,
                        sequence=norm,
                        collision_dsids=coll,
                    )
                )
        if not alleles:
            raise ValueError(f"no reference alleles loaded from {fasta_path}")
        return cls(alleles=alleles)


def _load_collisions(path: str) -> Dict[str, Tuple[str, ...]]:
    """Parse dsid_collisions.csv -> {reference_label: (dsid, dsid, ...)}."""
    out: Dict[str, Tuple[str, ...]] = {}
    if not os.path.exists(path):
        return out
    with open(path, "r", newline="") as fh:
        reader = csv.DictReader(fh)
        # tolerate either Kept_Reference/Associated_DSIds naming
        for row in reader:
            ref = row.get("Kept_Reference") or row.get("kept_reference")
            raw = row.get("Associated_DSIds") or row.get("associated_dsids")
            if not ref or not raw:
                continue
            dsids = tuple(d.strip() for d in raw.split(",") if d.strip())
            out[ref.strip()] = dsids
    return out
