import sys
import yaml
import requests
from pathlib import Path

WHITELIST = Path("genome_atlas/data/pfam_whitelist.yaml")

def check_pfam(acc: str) -> dict:
    r = requests.get(f"https://www.ebi.ac.uk/interpro/api/entry/pfam/{acc}/", timeout=15)
    if r.status_code != 200:
        return {"ok": False, "error": f"HTTP {r.status_code}"}
    meta = r.json().get("metadata", {})
    return {"ok": True, "count": meta.get("counters", {}).get("proteins", 0)}

def main() -> int:
    cfg = yaml.safe_load(WHITELIST.read_text())
    failures = []
    for entry in cfg["domains"] + cfg.get("auxiliary", []):
        res = check_pfam(entry["accession"])
        if not res["ok"]:
            failures.append(entry["accession"])
            print(f"FAIL {entry['accession']}: {res['error']}")
        else:
            print(f"OK   {entry['accession']}  {entry['name']:40s} -> {res['count']:>8,} proteins")
    return 1 if failures else 0

if __name__ == "__main__":
    sys.exit(main())
