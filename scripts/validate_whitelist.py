"""Verify every Pfam accession in the whitelist resolves to a real family.

Run via:
    docker run --rm -v ~/pen-stack/code/repos/genome-atlas:/pkg \
               -w /pkg pen-stack/data:0.1.0 \
               python scripts/validate_whitelist.py
"""
import sys
from pathlib import Path
import yaml
import requests

WHITELIST = Path("genome_atlas/data/pfam_whitelist.yaml")


def check_pfam(acc: str) -> dict:
    """Query InterPro for a Pfam accession. Returns {ok, name, count, url}."""
    url = f"https://www.ebi.ac.uk/interpro/api/entry/pfam/{acc}/"
    r = requests.get(url, timeout=15)
    if r.status_code != 200:
        return {"ok": False, "error": f"HTTP {r.status_code}", "url": url}
    data = r.json()
    meta = data.get("metadata", {})
    return {
        "ok": True,
        "accession": meta.get("accession"),
        "name": meta.get("name", {}).get("name"),
        "count": meta.get("counters", {}).get("proteins", 0),
        "url": url,
    }


def main() -> int:
    cfg = yaml.safe_load(WHITELIST.read_text())
    failures = []
    total_proteins = 0
    for entry in cfg["domains"] + cfg.get("auxiliary", []):
        result = check_pfam(entry["accession"])
        if not result["ok"]:
            print(f"FAIL  {entry['accession']}  {entry['name']}: {result.get('error')}")
            failures.append(entry["accession"])
        else:
            print(
                f"OK    {entry['accession']}  {entry['name']:40s}  "
                f"→ {result['count']:>8,} proteins in UniProt"
            )
            total_proteins += result["count"]
    if failures:
        print(f"\n{len(failures)} accession(s) failed to resolve: {failures}")
        return 1
    print(f"\nAll whitelist accessions validated.")
    print(f"Combined protein count (upper bound, before dedup): {total_proteins:,}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
