"""Command-line interface for measles-dsid-typer.

Usage examples
--------------
    # FASTA (one or many sequences), CSV to stdout
    measles-dsid-typer sample.fasta

    # explicit format, TSV to a file
    measles-dsid-typer reads.fastq.gz --format fastq --output-format tsv -o out.tsv

    # BAM input, JSON out, keep the extracted N-450 sequence
    measles-dsid-typer aligned.bam --output-format json --emit-sequence -o out.json

Output is deterministic: rows are sorted by sample id and columns are fixed.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from typing import List, Optional

from . import __version__
from .inputs import (
    bam_consensus,
    detect_format,
    fastq_consensus,
    iter_fasta_file,
)
from .n450 import _SeedIndex
from .reference import Reference
from .typer import RESULT_FIELDS, TyperResult, type_sequence


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="measles-dsid-typer",
        description="Deterministic N-450 DSID sub-lineage typing for measles virus.",
    )
    p.add_argument("input", help="Input file (FASTA / FASTQ / BAM; .gz ok).")
    p.add_argument(
        "--format",
        choices=["auto", "fasta", "fastq", "bam"],
        default="auto",
        help="Input format (default: auto-detect from extension/content).",
    )
    p.add_argument(
        "--output-format",
        choices=["csv", "tsv", "json"],
        default="csv",
        help="Output format (default: csv).",
    )
    p.add_argument(
        "-o", "--output", default="-", help="Output path ('-' for stdout, default)."
    )
    p.add_argument(
        "--reference",
        default=None,
        help="Override reference FASTA (defaults to the packaged N-450 set).",
    )
    p.add_argument(
        "--collisions",
        default=None,
        help="Override DSID collisions CSV (defaults to the packaged file).",
    )
    p.add_argument(
        "--emit-sequence",
        action="store_true",
        help="Include the extracted 450-nt N-450 sequence in the output.",
    )
    # consensus tuning (fastq / bam)
    p.add_argument("--min-depth", type=int, default=3,
                   help="Min read depth to call a consensus base (fastq/bam).")
    p.add_argument("--min-base-quality", type=int, default=20,
                   help="Min Phred base quality for consensus (fastq/bam).")
    p.add_argument("--sample-name", default=None,
                   help="Sample id to use for fastq/bam consensus (default: file stem).")
    p.add_argument("--version", action="version",
                   version=f"measles-dsid-typer {__version__}")
    return p


def _sample_stem(path: str) -> str:
    import os

    base = os.path.basename(path)
    for ext in (".gz",):
        if base.endswith(ext):
            base = base[: -len(ext)]
    base = os.path.splitext(base)[0]
    return base


def run(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)

    reference = Reference.load(args.reference, args.collisions)
    seed_index = _SeedIndex(reference)

    fmt = args.format
    if fmt == "auto":
        fmt = detect_format(args.input)

    results: List[TyperResult] = []

    if fmt == "fasta":
        for sample, seq in iter_fasta_file(args.input):
            results.append(
                type_sequence(
                    sample, seq, reference, seed_index,
                    keep_sequence=args.emit_sequence,
                )
            )
    elif fmt == "fastq":
        sample = args.sample_name or _sample_stem(args.input)
        consensus, stats = fastq_consensus(
            args.input, reference, seed_index,
            min_depth=args.min_depth,
            min_base_quality=args.min_base_quality,
        )
        res = type_sequence(
            sample, consensus, reference, seed_index,
            keep_sequence=args.emit_sequence,
        )
        note = (
            f"fastq consensus: {stats.get('reads_anchored', 0)}/"
            f"{stats.get('reads_total', 0)} reads anchored, "
            f"{stats.get('positions_covered', 0)}/{stats.get('positions_total', 0)} "
            f"positions covered"
        )
        res.note = (res.note + "; " + note) if res.note else note
        results.append(res)
    elif fmt == "bam":
        sample = args.sample_name or _sample_stem(args.input)
        consensus, stats = bam_consensus(
            args.input,
            min_depth=args.min_depth,
            min_base_quality=args.min_base_quality,
        )
        res = type_sequence(
            sample, consensus, reference, seed_index,
            keep_sequence=args.emit_sequence,
        )
        note = (
            f"bam consensus: {stats.get('positions_covered', 0)} positions covered"
        )
        res.note = (res.note + "; " + note) if res.note else note
        results.append(res)
    else:  # pragma: no cover
        print(f"unsupported format: {fmt}", file=sys.stderr)
        return 2

    # deterministic ordering
    results.sort(key=lambda r: r.sample)

    _emit(results, args.output_format, args.output, args.emit_sequence)
    return 0


def _emit(
    results: List[TyperResult],
    output_format: str,
    output: str,
    include_sequence: bool,
) -> None:
    fields = list(RESULT_FIELDS)
    if include_sequence:
        fields.append("n450_sequence")

    close = False
    if output == "-":
        fh = sys.stdout
    else:
        fh = open(output, "w", newline="")
        close = True
    try:
        if output_format == "json":
            payload = [r.as_row(include_sequence) for r in results]
            json.dump(payload, fh, indent=2, sort_keys=False)
            fh.write("\n")
        else:
            delimiter = "\t" if output_format == "tsv" else ","
            writer = csv.DictWriter(
                fh, fieldnames=fields, delimiter=delimiter, lineterminator="\n"
            )
            writer.writeheader()
            for r in results:
                writer.writerow(r.as_row(include_sequence))
    finally:
        if close:
            fh.close()


def main() -> None:  # console_scripts entry point
    sys.exit(run())


if __name__ == "__main__":
    main()
