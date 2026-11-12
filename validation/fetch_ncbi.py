#!/usr/bin/env python3
"""Pull public measles N-450 / WGS sequences from NCBI GenBank for validation.

This is a convenience fetcher. It uses NCBI E-utilities over HTTPS (stdlib only,
no Biopython needed) to download sequences by accession. Supply a file of
accessions (one per line) and it writes a single multi-FASTA.

The repo already ships a bundled public validation set
(``validation/data/public_d8_wgs.fasta`` + ground truth), so network access is
NOT required to reproduce the published concordance numbers. Use this script
only if you want to expand or refresh the validation set.

Example:
    python validation/fetch_ncbi.py \
        --accessions validation/data/extra_accessions.txt \
        --out validation/data/extra.fasta

Respect NCBI usage policy: set --email, keep batch sizes modest, and avoid
hammering the service. An API key (--api-key) raises the rate limit.
"""

from __future__ import annotations

import argparse
import sys
import time
import urllib.parse
import urllib.request
from typing import List

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


def fetch_batch(accessions: List[str], email: str, api_key: str = "") -> str:
    params = {
        "db": "nuccore",
        "id": ",".join(accessions),
        "rettype": "fasta",
        "retmode": "text",
        "tool": "measles-dsid-typer",
        "email": email,
    }
    if api_key:
        params["api_key"] = api_key
    url = EUTILS + "?" + urllib.parse.urlencode(params)
    with urllib.request.urlopen(url, timeout=60) as resp:
        return resp.read().decode("utf-8", errors="replace")


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--accessions", required=True,
                   help="File with one GenBank accession per line.")
    p.add_argument("--out", required=True, help="Output multi-FASTA path.")
    p.add_argument("--email", required=True, help="Contact email (NCBI policy).")
    p.add_argument("--api-key", default="", help="NCBI API key (optional).")
    p.add_argument("--batch-size", type=int, default=100)
    p.add_argument("--sleep", type=float, default=0.4,
                   help="Seconds between batches (be polite to NCBI).")
    args = p.parse_args(argv)

    with open(args.accessions) as fh:
        accs = [ln.strip() for ln in fh if ln.strip() and not ln.startswith("#")]
    if not accs:
        print("no accessions found", file=sys.stderr)
        return 1

    n_written = 0
    with open(args.out, "w") as out:
        for i in range(0, len(accs), args.batch_size):
            batch = accs[i : i + args.batch_size]
            fasta = fetch_batch(batch, args.email, args.api_key)
            out.write(fasta)
            if not fasta.endswith("\n"):
                out.write("\n")
            n_written += fasta.count(">")
            print(f"fetched {i + len(batch)}/{len(accs)} accessions...",
                  file=sys.stderr)
            time.sleep(args.sleep)
    print(f"wrote {n_written} sequences to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
