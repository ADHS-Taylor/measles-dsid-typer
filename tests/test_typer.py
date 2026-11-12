"""Unit + regression tests for measles-dsid-typer.

Run with:  pytest -q  (from the repo root, after `pip install -e .[dev]`)
or:        PYTHONPATH=src pytest -q
"""

import json
import os
import subprocess
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
SRC = os.path.join(REPO, "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from measles_dsid_typer import sequtils  # noqa: E402
from measles_dsid_typer.reference import Reference  # noqa: E402
from measles_dsid_typer.n450 import _SeedIndex, extract_n450  # noqa: E402
from measles_dsid_typer.typer import type_sequence  # noqa: E402


@pytest.fixture(scope="module")
def reference():
    return Reference.load()


@pytest.fixture(scope="module")
def seed_index(reference):
    return _SeedIndex(reference)


# --------------------------------------------------------------------------- #
# sequtils: IUPAC handling
# --------------------------------------------------------------------------- #
def test_normalize_uracil():
    assert sequtils.normalize("acguACGU") == "ACGTACGT"


def test_reverse_complement_iupac():
    assert sequtils.reverse_complement("ACGT") == "ACGT"
    assert sequtils.reverse_complement("AARYN") == "NRYTT"


def test_bases_compatible_concrete():
    assert sequtils.bases_compatible("A", "A")
    assert not sequtils.bases_compatible("A", "C")


def test_bases_compatible_ambiguous():
    # R = A or G, compatible with A and with G, not with C/T
    assert sequtils.bases_compatible("R", "A")
    assert sequtils.bases_compatible("R", "G")
    assert not sequtils.bases_compatible("R", "C")
    # N is compatible with everything
    for b in "ACGT":
        assert sequtils.bases_compatible("N", b)


def test_gap_never_compatible():
    assert not sequtils.bases_compatible("-", "A")
    assert not sequtils.bases_compatible("A", "-")


def test_count_mismatches_iupac_not_counted():
    ref = "ACGTACGT"
    # query with N and R at compatible positions -> 0 mismatches, 2 ambiguous
    qry = "ACGTNCGR"  # pos4 N (ref A) ok; pos7 R (ref T) -> R!=T => mismatch
    mm, amb, called = sequtils.count_mismatches(qry, ref)
    assert amb == 2
    assert called == 8
    assert mm == 1  # only the R-vs-T position


def test_count_mismatches_length_guard():
    with pytest.raises(ValueError):
        sequtils.count_mismatches("ACGT", "ACG")


# --------------------------------------------------------------------------- #
# reference loading
# --------------------------------------------------------------------------- #
def test_reference_loads(reference):
    assert reference.n_alleles == 813
    assert set(reference.genotypes) == {"B3", "D8"}
    assert len(reference.dsids) == 77


def test_reference_all_450(reference):
    assert all(len(a.sequence) == 450 for a in reference.alleles)


def test_reference_collisions_present(reference):
    # D8-9171 reference is known to collide with D8-9421/D8-9274/D8-9461
    colliding = [a for a in reference.alleles if a.collision_dsids]
    assert colliding, "expected at least one allele with recorded collisions"


# --------------------------------------------------------------------------- #
# exact / nearest typing
# --------------------------------------------------------------------------- #
def test_exact_match_every_dsid(reference, seed_index):
    # Every reference allele must type back to its own DSID with 0 mismatches.
    seen = set()
    for allele in reference.alleles:
        if allele.dsid in seen:
            continue
        seen.add(allele.dsid)
        r = type_sequence("t", allele.sequence, reference, seed_index)
        assert r.status == "typed"
        assert r.n450_coordinate_match is True
        assert r.mismatches == 0
        assert r.pident == 100.0
        # the exact-matched allele should itself represent this DSID
        # (either as primary dsid or via identical-sequence collision set)
        covered = {r.dsid, *r.collision_dsids}
        assert allele.dsid in covered


def test_single_snp_is_nearest(reference, seed_index):
    a = reference.alleles[1]
    mut = list(a.sequence)
    mut[100] = "A" if mut[100] != "A" else "C"
    r = type_sequence("snp", "".join(mut), reference, seed_index)
    assert r.match_type == "nearest"
    assert r.n450_coordinate_match is False
    assert r.mismatches >= 1
    assert r.pident < 100.0


def test_iupac_ambiguity_not_mismatch(reference, seed_index):
    a = reference.alleles[1]
    amb = list(a.sequence)
    amb[5] = "N"   # compatible with anything
    r = type_sequence("amb", "".join(amb), reference, seed_index)
    # still exact-compatible: N does not create a mismatch
    assert r.mismatches == 0
    assert r.ambiguous_positions >= 1


def test_reverse_complement_detected(reference, seed_index):
    a = reference.alleles[1]
    rc = sequtils.reverse_complement(a.sequence)
    r = type_sequence("rc", rc, reference, seed_index)
    assert r.orientation == "reverse_complement"
    assert r.mismatches == 0


def test_wgs_extraction(reference, seed_index):
    a = reference.alleles[1]
    flank = "ATGCATGCGT" * 50
    wgs = flank + a.sequence + flank
    r = type_sequence("wgs", wgs, reference, seed_index)
    assert r.extraction_method == "anchored"
    assert r.n450_length == 450
    assert r.mismatches == 0


def test_too_short_input(reference, seed_index):
    r = type_sequence("short", "ACGT" * 5, reference, seed_index)
    assert r.status == "no_n450"


# --------------------------------------------------------------------------- #
# determinism
# --------------------------------------------------------------------------- #
def test_determinism(reference, seed_index):
    a = reference.alleles[10]
    mut = list(a.sequence)
    mut[50] = "N"
    mut[300] = "A" if mut[300] != "A" else "T"
    s = "".join(mut)
    r1 = type_sequence("d", s, reference, seed_index)
    r2 = type_sequence("d", s, reference, seed_index)
    assert r1.as_row() == r2.as_row()


# --------------------------------------------------------------------------- #
# FASTQ consensus
# --------------------------------------------------------------------------- #
def test_fastq_consensus(reference, seed_index, tmp_path):
    from measles_dsid_typer.inputs import fastq_consensus

    a = reference.alleles[1]
    seq = a.sequence
    fq = tmp_path / "reads.fastq"
    rl = 100
    with open(fq, "w") as fh:
        for cov in range(8):
            for start in range(0, len(seq) - rl + 1, 25):
                # ensure the 3' end is covered too
                sub = seq[start : start + rl]
                fh.write(f"@r{cov}_{start}\n{sub}\n+\n{'I' * len(sub)}\n")
            tail = seq[-rl:]
            fh.write(f"@r{cov}_tail\n{tail}\n+\n{'I' * len(tail)}\n")
    cons, stats = fastq_consensus(str(fq), reference, seed_index, min_depth=2)
    assert stats["positions_covered"] == 450
    r = type_sequence("fq", cons, reference, seed_index)
    assert r.mismatches == 0


# --------------------------------------------------------------------------- #
# CLI integration
# --------------------------------------------------------------------------- #
def _run_cli(args):
    env = dict(os.environ)
    env["PYTHONPATH"] = SRC + os.pathsep + env.get("PYTHONPATH", "")
    return subprocess.run(
        [sys.executable, "-m", "measles_dsid_typer.cli", *args],
        capture_output=True,
        text=True,
        env=env,
    )


def test_cli_csv(reference, tmp_path):
    a = reference.alleles[1]
    fa = tmp_path / "s.fasta"
    fa.write_text(f">s1\n{a.sequence}\n")
    res = _run_cli([str(fa)])
    assert res.returncode == 0
    lines = res.stdout.strip().splitlines()
    assert lines[0].startswith("sample,status,genotype,dsid")
    assert "exact" in lines[1]
    assert "True" in lines[1]


def test_cli_json(reference, tmp_path):
    a = reference.alleles[1]
    fa = tmp_path / "s.fasta"
    fa.write_text(f">s1\n{a.sequence}\n")
    res = _run_cli([str(fa), "--output-format", "json"])
    assert res.returncode == 0
    payload = json.loads(res.stdout)
    assert payload[0]["dsid"].startswith("D8") or payload[0]["dsid"].startswith("B3")
    assert payload[0]["n450_coordinate_match"] is True


def test_cli_deterministic(reference, tmp_path):
    a = reference.alleles[2]
    b = reference.alleles[3]
    fa = tmp_path / "multi.fasta"
    fa.write_text(f">zzz\n{a.sequence}\n>aaa\n{b.sequence}\n")
    r1 = _run_cli([str(fa)])
    r2 = _run_cli([str(fa)])
    assert r1.stdout == r2.stdout
    # rows sorted by sample id -> 'aaa' before 'zzz'
    lines = r1.stdout.strip().splitlines()
    assert lines[1].startswith("aaa")
    assert lines[2].startswith("zzz")
