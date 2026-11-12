# measles-dsid-typer

Deterministic **N-450 DSID sub-lineage typing** for measles virus, from
FASTA, FASTQ, or BAM input.

Standard measles genotyping resolves sequences to a genotype (e.g. **D8**,
**B3**) — often too coarse to separate co-circulating transmission chains. This
tool adds the finer **DSID sub-lineage** layer (the WHO/MeaNS-style
distinct-sequence identifier) on top of the genotype, using the **N-450**
region: the 450 nucleotides encoding the carboxyl-terminal 150 amino acids of
the nucleoprotein, which is the WHO-standard window for measles molecular
surveillance.

- **Deterministic.** The same input always produces byte-identical output.
- **IUPAC-aware.** Ambiguous bases (N, R, Y, …) are handled explicitly and
  never silently counted as mismatches; genuinely ambiguous DSID calls are
  flagged, not hidden.
- **Zero-dependency core.** FASTA and FASTQ typing need only the Python
  standard library. BAM/CRAM support is an optional `pysam` extra.
- **Validated.** 100 % recovery across all 813 curated reference alleles and
  **98.9 % concordance** against published DSIDs on 660 public whole-genome
  sequences (see [Validation](#validation)).

---

## Quick start

### Install

```bash
# from a clone
pip install .

# with BAM/CRAM support
pip install ".[bam]"
```

Or run in a container (no local Python needed):

```bash
docker build -t measles-dsid-typer .
docker run --rm -v "$PWD":/data measles-dsid-typer /data/sample.fasta
```

### Type sequences

```bash
# FASTA (one or many sequences) -> CSV on stdout
measles-dsid-typer sample.fasta

# TSV to a file
measles-dsid-typer sample.fasta --output-format tsv -o results.tsv

# JSON for downstream pipelines, keeping the extracted N-450 sequence
measles-dsid-typer sample.fasta --output-format json --emit-sequence -o results.json

# FASTQ amplicon reads -> consensus -> typing
measles-dsid-typer reads.fastq.gz --format fastq --min-depth 5

# aligned BAM -> consensus -> typing (needs the [bam] extra)
measles-dsid-typer aligned.bam --format bam
```

Try it on the bundled example:

```bash
measles-dsid-typer examples/example.fasta
```

```
sample,status,genotype,dsid,...,match_type,n450_coordinate_match,mismatches,pident,...
sample_2snp_D8-9171,typed,D8,D8-9171,...,nearest,False,2,99.556,...
sample_exact_D8-10887,typed,D8,D8-10887,...,exact,True,0,100.000,...
sample_iupac_B3-6418,typed,B3,B3-6418,...,nearest,False,0,100.000,...   # IUPAC handled
sample_revcomp_D8-9461,typed,D8,D8-9461,...,exact,True,0,100.000,...     # reverse-complement auto-detected
```

### Python API

```python
from measles_dsid_typer import Reference, type_sequence
from measles_dsid_typer.n450 import _SeedIndex

ref = Reference.load()              # bundled curated N-450 reference
idx = _SeedIndex(ref)               # optional: build once, reuse across calls

result = type_sequence("my_sample", my_sequence, ref, idx)
print(result.genotype, result.dsid, result.mismatches, result.n450_coordinate_match)
```

---

## Output columns

| column | meaning |
|--------|---------|
| `sample` | sequence id (FASTA header token, or file stem for FASTQ/BAM) |
| `status` | `typed`, `no_n450` (window not locatable) |
| `genotype` | WHO genotype, e.g. `D8`, `B3` |
| `dsid` | assigned DSID sub-lineage, e.g. `D8-9171` |
| `collision_dsids` | other DSIDs represented by a **byte-identical** reference allele |
| `match_type` | `exact` (byte-identical N-450) or `nearest` (closest allele) |
| `n450_coordinate_match` | `True` iff the 450-nt N-450 is byte-identical to a reference |
| `mismatches` | IUPAC-aware substitution count vs the best reference allele |
| `pident` | percent identity over called (non-gap) positions |
| `aligned_length` | number of called positions compared |
| `ambiguous_positions` | IUPAC/N positions in the query that were handled |
| `tied_dsids` | every DSID tied at the minimum mismatch count |
| `dsid_ambiguous` | `True` iff more than one DSID is tied (ambiguous call) |
| `n450_length` | length of the extracted N-450 window (450 when complete) |
| `orientation` | `forward` or `reverse_complement` (auto-detected) |
| `extraction_method` | `verbatim` (input was ~450 nt) or `anchored` (carved from a longer sequence) |
| `best_reference` | the reference allele header that produced the call |
| `note` | free-text diagnostics (e.g. FASTQ/BAM consensus coverage) |

---

## Architecture

```
input (FASTA / FASTQ / BAM)
        │
        ▼
 ┌─────────────────┐   FASTQ/BAM → per-position majority consensus
 │  inputs.py      │   (depth + base-quality filtered)
 └─────────────────┘
        │  nucleotide sequence(s)
        ▼
 ┌─────────────────┐   locate the 450-nt N-450 window by reference k-mer
 │  n450.py        │   seeding; resolve forward vs reverse-complement
 └─────────────────┘
        │  450-nt N-450
        ▼
 ┌─────────────────┐   1. byte-exact reference lookup  → coordinate match
 │  typer.py       │   2. else nearest allele by IUPAC-aware mismatch count
 │  (sequtils.py)  │   3. record genotype, DSID, ties, collisions, identity
 └─────────────────┘
        │  TyperResult
        ▼
 ┌─────────────────┐   deterministic CSV / TSV / JSON
 │  cli.py         │   (rows sorted by sample, fixed column order)
 └─────────────────┘
```

**N-450 extraction.** Inputs may be pre-trimmed N-450 fragments, full N-gene
amplicons, or whole-genome assemblies (~15.9 kb). The window is located by
anchoring reference-derived k-mer seeds and taking the modal 450-nt frame;
orientation is resolved automatically. No external genome or network access is
required — extraction is anchored to the reference allele set itself.

**IUPAC handling.** Every base comparison expands IUPAC codes to their concrete
base sets. Two bases are a *match* when their sets intersect (so `N` matches
anything, `R` matches `A`/`G`) and a *mismatch* only when the sets are disjoint.
Gaps never match. When ambiguity leaves several DSIDs equally compatible, all of
them are reported in `tied_dsids` and `dsid_ambiguous=True`.

**Determinism.** Reference scanning uses a stable tie-break
(`mismatches`, then DSID, then accession), rows are sorted by sample id, and
columns are fixed — so output is reproducible across runs, machines, and
reference-file orderings.

---

## Reference data

`data/reference/dsid_n450_reference.fasta` — **813 curated N-450 alleles**
(all exactly 450 nt), spanning **77 distinct DSIDs** across genotypes **B3**
(122 alleles) and **D8** (691 alleles). Headers follow the MeaNS-derived
convention `>DSId_<GENOTYPE>-<DSID>_<ACCESSION>` (e.g.
`>DSId_D8-9171_PZ023615.1`). The accessions are public GenBank / Pathoplexus
records (`data/reference/reference_accessions.txt`).

`data/reference/dsid_collisions.csv` — records, for each kept reference, every
DSID that a byte-identical N-450 is known to represent (surfaced as
`collision_dsids`).

To use your own reference set, pass `--reference my_refs.fasta`
(and optionally `--collisions my_collisions.csv`).

---

## Benchmark

Measured on the bundled public validation set (660 measles D8 whole-genome
sequences, ~15.9 kb each), Python 3.10, single core.

| metric | value |
|--------|-------|
| sequences typed | 660 whole genomes |
| wall-clock time | ~57 s (~0.09 s / genome) |
| N-450 extracted from | full ~15.9 kb assemblies (anchored) |
| exact N-450 coordinate matches | 629 / 660 (95.3 %) |
| nearest-allele calls | 31 / 660 |
| core runtime dependencies | **none** (stdlib only) |

Reproduce with:

```bash
bash examples/benchmark_run.sh
```

---

## Validation

Two independent checks, both reproducible offline from bundled data.

### 1. Reference self-check (integrity)

Every curated reference allele is typed back through the full pipeline and must
recover its own published DSID by exact coordinate match.

| alleles | recovered | recovery rate |
|---------|-----------|---------------|
| 813 | 813 | **100.00 %** |

```bash
python validation/validate.py --mode reference-selfcheck
```

### 2. Empirical concordance vs published DSIDs (660 public WGS)

The typer's assignments are compared to published DSIDs for 660 public measles
D8 whole-genome sequences (Pathoplexus / GenBank accessions; ground truth in
`validation/data/ground_truth_blast.tsv`).

| metric | count | rate |
|--------|-------|------|
| sequences with a published DSID | 660 | — |
| **concordant** | **653** | **98.94 %** |
| &nbsp;&nbsp;• via exact / collision match | 644 | 97.58 % |
| &nbsp;&nbsp;• via IUPAC-tie (published DSID among tied candidates) | 9 | 1.36 % |
| flagged as ambiguous (`dsid_ambiguous=True`) | 14 | 2.12 % |
| discordant | 7 | 1.06 % |
| total IUPAC/N positions handled | 1 287 | — |

```bash
python validation/validate.py --mode public-wgs \
    --sequences validation/data/public_d8_wgs.fasta \
    --truth     validation/data/ground_truth_blast.tsv
```

#### Mismatch analysis

All 7 discordant sequences carry heavy N-masking in the N-450 window (27–189
ambiguous positions out of 450). In each case the sequence is a **0-mismatch,
IUPAC-compatible** match to one DSID but **1–2 mismatches** from the published
DSID. The published labels come from a BLAST-bitscore pipeline
(`max_target_seqs 1`), which can prefer a slightly more-mismatched reference
when ambiguity codes inflate the alignment score. By strict N-450 identity the
typer's call is the closer match; these cases are the expected divergence
between heuristic bitscore assignment and exact N-450 coordinate identity on
low-quality (N-rich) sequences. They are surfaced — not hidden — via
`dsid_ambiguous` and `tied_dsids`.

#### IUPAC handling in the validation matrix

Of the 660 sequences, 14 produced genuinely ambiguous N-450 windows (multiple
DSIDs compatible at 0 mismatches). The typer flags every one
(`dsid_ambiguous=True`) and lists all compatible DSIDs in `tied_dsids`, so a
reviewer can see exactly where N-450 quality limits sub-lineage resolution
rather than receiving a false-precision single call.

### Refreshing / expanding the validation set

The bundled set reproduces the numbers above without network access. To pull
additional public sequences from NCBI GenBank:

```bash
python validation/fetch_ncbi.py \
    --accessions validation/data/extra_accessions.txt \
    --out validation/data/extra.fasta \
    --email you@example.org
```

---

## Testing

```bash
pip install ".[dev]"
pytest -q           # 21 unit + regression + CLI integration tests
```

The suite covers IUPAC compatibility logic, reference loading, exact/nearest
assignment, reverse-complement detection, WGS extraction, FASTQ consensus,
determinism, and the CSV/TSV/JSON CLI.

---

## Install via Bioconda

A Bioconda recipe is provided in `conda-recipe/meta.yaml`. After the package is
published to a tagged GitHub release (and the recipe submitted to
bioconda-recipes), installation will be:

```bash
conda install -c bioconda measles-dsid-typer
```

---

## How DSID typing relates to the broader workflow

This typer is the sub-lineage layer referenced in the measles genomic
investigation toolkit: genotype-level typing (D8/B3) is too coarse to separate
transmission chains, whereas N-450 DSID sub-lineages (e.g. `D8-9171`,
`D8-10574`, `D8-10707`) resolve distinct introductions vs. one diversifying
chain. The DSID assignments produced here feed directly into downstream
sub-lineage SNP-distance, monophyly, and import-vs-local analyses.

---

## License

MIT — see [`LICENSE`](LICENSE). Copyright © 2026 Taylor Martins.

## Citing / method background

- WHO. *Measles virus nomenclature update* and the N-450 genotyping window
  (450 nt encoding the C-terminal 150 aa of the nucleoprotein).
- Bankamp B, et al. *Improving molecular tools for global surveillance of
  measles virus.* J Clin Virol. 2013.

> Reference alleles are derived from public GenBank / Pathoplexus records. No
> case-level or personally identifying data are included in this repository.
