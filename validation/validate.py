#!/usr/bin/env python3
"""Empirical validation harness for measles-dsid-typer.

Compares the typer's DSID assignments against published ground truth and
reports a concordance matrix.

Two validation modes are supported:

1. ``--mode public-wgs`` (default)
   Types a set of public measles D8 whole-genome sequences (Pathoplexus /
   GenBank accessions) and compares the assigned DSID to a published
   ground-truth TSV (accession -> DSID). This exercises the full pipeline:
   N-450 extraction from ~15.9 kb genomes + DSID assignment.

2. ``--mode reference-selfcheck``
   Types every curated N-450 reference allele back through the typer and
   confirms it recovers its own published DSID. This validates exact-match
   integrity and the collision handling across all 77 DSIDs.

Outputs:
  * a per-sample results TSV
  * a concordance summary (printed + JSON)
  * a mismatch-analysis table (where typer != ground truth)

Usage:
  python validation/validate.py \
      --sequences path/to/sequences.fasta \
      --truth path/to/ground_truth.tsv \
      --out-dir validation/out
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys
from collections import Counter, defaultdict
from typing import Dict, List, Optional, Tuple

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
SRC = os.path.join(REPO, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from measles_dsid_typer.reference import Reference  # noqa: E402
from measles_dsid_typer.n450 import _SeedIndex  # noqa: E402
from measles_dsid_typer.inputs import iter_fasta_file  # noqa: E402
from measles_dsid_typer.typer import type_sequence  # noqa: E402


def load_truth(path: str) -> Dict[str, str]:
    """Load accession -> published DSID. Accepts TSV or CSV with a header."""
    truth: Dict[str, str] = {}
    delim = "\t" if path.endswith((".tsv", ".txt")) else ","
    with open(path, newline="") as fh:
        reader = csv.reader(fh, delimiter=delim)
        header = next(reader, None)
        # find the accession and dsid columns heuristically
        acc_i, dsid_i = 0, 1
        if header:
            lower = [h.lower() for h in header]
            for i, h in enumerate(lower):
                if h in ("accession", "query_id", "id", "sample", "accession_id"):
                    acc_i = i
                if "dsid" in h or "lineage" in h or "assigned" in h:
                    dsid_i = i
        for row in reader:
            if len(row) <= max(acc_i, dsid_i):
                continue
            truth[row[acc_i].strip()] = row[dsid_i].strip()
    return truth


def run_validation(
    sequences: str,
    truth_path: str,
    out_dir: str,
    reference: Reference,
) -> dict:
    os.makedirs(out_dir, exist_ok=True)
    seed_index = _SeedIndex(reference)
    truth = load_truth(truth_path)

    rows: List[dict] = []
    n_total = 0
    n_concordant = 0
    n_no_truth = 0
    n_untyped = 0
    mismatches: List[Tuple[str, str, str, int]] = []
    ambiguity_total = 0
    exact_matches = 0
    n_ambiguous_calls = 0
    n_concordant_unambiguous = 0
    n_concordant_via_tie = 0

    for sample, seq in iter_fasta_file(sequences):
        n_total += 1
        r = type_sequence(sample, seq, reference, seed_index)
        published = truth.get(sample)
        ambiguity_total += r.ambiguous_positions
        if r.n450_coordinate_match:
            exact_matches += 1
        if r.dsid_ambiguous:
            n_ambiguous_calls += 1

        if r.status != "typed":
            n_untyped += 1
            concordant = False
        elif published is None:
            n_no_truth += 1
            concordant = None
        else:
            # concordant if the published DSID is the top pick, appears in the
            # collision set (identical reference sequences), OR is among the
            # IUPAC-tied best candidates (genuine ambiguity from N/IUPAC bases).
            top_set = {r.dsid, *r.collision_dsids}
            tie_set = set(r.tied_dsids)
            if published in top_set:
                concordant = True
                n_concordant += 1
                n_concordant_unambiguous += 1
            elif published in tie_set:
                concordant = True
                n_concordant += 1
                n_concordant_via_tie += 1
            else:
                concordant = False
                mismatches.append(
                    (sample, published, r.dsid or "NA", r.mismatches or -1)
                )

        rows.append(
            {
                "sample": sample,
                "published_dsid": published or "",
                "assigned_dsid": r.dsid or "",
                "collision_dsids": ";".join(r.collision_dsids),
                "tied_dsids": ";".join(r.tied_dsids),
                "dsid_ambiguous": r.dsid_ambiguous,
                "concordant": "" if concordant is None else str(concordant),
                "match_type": r.match_type or "",
                "n450_coordinate_match": r.n450_coordinate_match,
                "mismatches": r.mismatches if r.mismatches is not None else "",
                "pident": f"{r.pident:.3f}" if r.pident is not None else "",
                "ambiguous_positions": r.ambiguous_positions,
                "status": r.status,
            }
        )

    # write per-sample results
    results_path = os.path.join(out_dir, "validation_results.tsv")
    fields = list(rows[0].keys()) if rows else []
    with open(results_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for row in sorted(rows, key=lambda x: x["sample"]):
            w.writerow(row)

    n_with_truth = n_total - n_no_truth
    concordance = (n_concordant / n_with_truth) if n_with_truth else 0.0

    summary = {
        "sequences_file": os.path.abspath(sequences),
        "truth_file": os.path.abspath(truth_path),
        "n_sequences": n_total,
        "n_with_published_dsid": n_with_truth,
        "n_without_published_dsid": n_no_truth,
        "n_untyped": n_untyped,
        "n_concordant": n_concordant,
        "n_concordant_unambiguous": n_concordant_unambiguous,
        "n_concordant_via_iupac_tie": n_concordant_via_tie,
        "n_ambiguous_calls": n_ambiguous_calls,
        "n_discordant": len(mismatches),
        "concordance_rate": round(concordance, 5),
        "exact_n450_matches": exact_matches,
        "total_ambiguous_positions_handled": ambiguity_total,
        "discordant_examples": [
            {
                "sample": m[0],
                "published": m[1],
                "assigned": m[2],
                "mismatches": m[3],
            }
            for m in mismatches[:25]
        ],
    }
    with open(os.path.join(out_dir, "validation_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)

    return summary


def run_reference_selfcheck(reference: Reference, out_dir: str) -> dict:
    os.makedirs(out_dir, exist_ok=True)
    seed_index = _SeedIndex(reference)
    per_dsid = defaultdict(lambda: {"n": 0, "ok": 0})
    n_total = 0
    n_ok = 0
    for allele in reference.alleles:
        n_total += 1
        r = type_sequence(allele.label, allele.sequence, reference, seed_index)
        covered = {r.dsid, *r.collision_dsids}
        ok = r.n450_coordinate_match and (allele.dsid in covered)
        per_dsid[allele.dsid]["n"] += 1
        if ok:
            per_dsid[allele.dsid]["ok"] += 1
            n_ok += 1
    summary = {
        "n_alleles": n_total,
        "n_recovered": n_ok,
        "recovery_rate": round(n_ok / n_total, 5) if n_total else 0.0,
        "per_dsid": {k: v for k, v in sorted(per_dsid.items())},
    }
    with open(os.path.join(out_dir, "reference_selfcheck.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    return summary


def main(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--mode",
        choices=["public-wgs", "reference-selfcheck"],
        default="public-wgs",
    )
    p.add_argument("--sequences", help="FASTA of sequences to type.")
    p.add_argument("--truth", help="Ground-truth TSV/CSV (accession -> DSID).")
    p.add_argument("--out-dir", default=os.path.join(HERE, "out"))
    p.add_argument("--reference", default=None)
    p.add_argument("--collisions", default=None)
    args = p.parse_args(argv)

    reference = Reference.load(args.reference, args.collisions)

    if args.mode == "reference-selfcheck":
        summary = run_reference_selfcheck(reference, args.out_dir)
        print("== Reference self-check ==")
        print(f"alleles        : {summary['n_alleles']}")
        print(f"recovered      : {summary['n_recovered']}")
        print(f"recovery rate  : {summary['recovery_rate'] * 100:.2f}%")
        return 0 if summary["n_recovered"] == summary["n_alleles"] else 1

    if not args.sequences or not args.truth:
        p.error("--sequences and --truth are required for public-wgs mode")

    summary = run_validation(args.sequences, args.truth, args.out_dir, reference)
    print("== Empirical validation (public WGS) ==")
    print(f"sequences typed       : {summary['n_sequences']}")
    print(f"with published DSID   : {summary['n_with_published_dsid']}")
    print(f"concordant            : {summary['n_concordant']}")
    print(f"  exact/collision     : {summary['n_concordant_unambiguous']}")
    print(f"  via IUPAC tie       : {summary['n_concordant_via_iupac_tie']}")
    print(f"ambiguous calls       : {summary['n_ambiguous_calls']}")
    print(f"discordant            : {summary['n_discordant']}")
    print(f"concordance rate      : {summary['concordance_rate'] * 100:.2f}%")
    print(f"exact N-450 matches   : {summary['exact_n450_matches']}")
    print(f"ambiguous pos handled : {summary['total_ambiguous_positions_handled']}")
    if summary["discordant_examples"]:
        print("\n-- discordant examples --")
        for ex in summary["discordant_examples"]:
            print(
                f"  {ex['sample']}: published={ex['published']} "
                f"assigned={ex['assigned']} (mm={ex['mismatches']})"
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
