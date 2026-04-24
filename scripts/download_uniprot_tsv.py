"""Download UniProt TSV slices for each Pfam accession in the whitelist.

Run via:
    python scripts/download_uniprot_tsv.py \
        --whitelist genome_atlas/data/pfam_whitelist.yaml \
        --output data/raw/uniprot/
"""
from __future__ import annotations

import argparse
import gzip
import sys
import time
from pathlib import Path

import requests
import yaml

UNIPROT_STREAM = "https://rest.uniprot.org/uniprotkb/stream"
FIELDS = (
    "accession,organism_name,organism_id,length,xref_pfam,xref_pdb,"
    "xref_alphafolddb,sequence,reviewed,protein_name,lineage_ids"
)


def download_one(acc: str, output_dir: Path, max_retries: int = 5) -> Path:
    """Download TSV for one Pfam accession, gzipped on disk."""
    out = output_dir / f"{acc}.tsv.gz"
    if out.exists() and out.stat().st_size > 0:
        print(f"  [skip] {acc} already downloaded ({out.stat().st_size / 1e6:.1f} MB)")
        return out

    params = {
        "format": "tsv",
        "query": f"xref:pfam-{acc}",
        "fields": FIELDS,
        "compressed": "true",
    }

    for attempt in range(max_retries):
        try:
            with requests.get(UNIPROT_STREAM, params=params, stream=True, timeout=300) as r:
                r.raise_for_status()
                with open(out, "wb") as f:
                    for chunk in r.iter_content(chunk_size=65536):
                        f.write(chunk)
            size_mb = out.stat().st_size / 1e6
            print(f"  [ok]   {acc} -> {out.name} ({size_mb:.1f} MB)")
            return out
        except (requests.exceptions.RequestException, IOError) as e:
            wait = 2 ** attempt
            print(f"  [retry {attempt+1}/{max_retries}] {acc}: {e} -- sleeping {wait}s")
            time.sleep(wait)
    raise RuntimeError(f"Failed to download {acc} after {max_retries} retries")


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--whitelist", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()

    args.output.mkdir(parents=True, exist_ok=True)
    cfg = yaml.safe_load(args.whitelist.read_text())
    accs = [e["accession"] for e in cfg["domains"]]

    print(f"Downloading UniProt TSV for {len(accs)} Pfam accessions -> {args.output}")
    for acc in accs:
        download_one(acc, args.output)
    print("\nAll downloads complete.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
