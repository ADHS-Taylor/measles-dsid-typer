"""measles-dsid-typer: deterministic N-450 DSID sub-lineage typing for measles virus.

Assigns WHO/MeaNS-style DSID sub-lineages to measles sequences from the N-450
region (the 450 nucleotides encoding the carboxyl-terminal 150 amino acids of
the nucleoprotein, the WHO-standard genotyping window).

Public API:
    from measles_dsid_typer import type_sequence, TyperResult, Reference
"""

from .reference import Reference, ReferenceAllele
from .typer import TyperResult, type_sequence, type_record
from .n450 import extract_n450, N450ExtractionError

__all__ = [
    "Reference",
    "ReferenceAllele",
    "TyperResult",
    "type_sequence",
    "type_record",
    "extract_n450",
    "N450ExtractionError",
]

__version__ = "0.1.0"
