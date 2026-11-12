#!/usr/bin/env bash
# benchmark_run.sh — end-to-end demonstration + benchmark for measles-dsid-typer.
#
# Runs the typer on the bundled example (FASTA), reports timing on the full
# public validation set, and executes both validation modes. No network access
# required: everything uses data shipped in the repo.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"
export PYTHONPATH="${REPO_ROOT}/src:${PYTHONPATH:-}"

PY="${PYTHON:-python3}"
TYPER="$PY -m measles_dsid_typer.cli"

echo "=============================================================="
echo " measles-dsid-typer benchmark"
echo "=============================================================="

echo
echo ">> 1. Example FASTA (CSV output)"
$TYPER examples/example.fasta

echo
echo ">> 2. Example FASTA (JSON output, N-450 sequence included)"
$TYPER examples/example.fasta --output-format json --emit-sequence | head -40

echo
echo ">> 3. Throughput on the public validation set (660 WGS genomes)"
SECONDS=0
$TYPER validation/data/public_d8_wgs.fasta --output-format tsv \
    > /tmp/bench_out.tsv
N=$(( $(wc -l < /tmp/bench_out.tsv) - 1 ))
echo "    typed ${N} whole-genome sequences in ${SECONDS}s"

echo
echo ">> 4. Determinism check (byte-identical across two runs)"
$TYPER validation/data/public_d8_wgs.fasta > /tmp/run_a.csv
$TYPER validation/data/public_d8_wgs.fasta > /tmp/run_b.csv
if diff -q /tmp/run_a.csv /tmp/run_b.csv >/dev/null; then
    echo "    OK: output is deterministic"
else
    echo "    FAIL: output differs between runs" >&2
    exit 1
fi

echo
echo ">> 5. Empirical validation vs published DSIDs"
$PY validation/validate.py --mode public-wgs \
    --sequences validation/data/public_d8_wgs.fasta \
    --truth validation/data/ground_truth_blast.tsv \
    --out-dir validation/out

echo
echo ">> 6. Reference self-check (all 813 curated alleles)"
$PY validation/validate.py --mode reference-selfcheck --out-dir validation/out

echo
echo "Done. Detailed outputs in validation/out/."
