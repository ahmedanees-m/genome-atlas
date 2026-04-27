"""Upgraded whitelist validator - catches name/description mismatches as well as existence.

Run:
    python scripts/validate_whitelist.py --strict
"""
from __future__ import annotations
import argparse
import sys
import time
from pathlib import Path

import yaml
import requests

WHITELIST = Path("genome_atlas/data/pfam_whitelist.yaml")


def fetch_interpro_pfam(acc: str) -> dict | None:
    url = f"https://www.ebi.ac.uk/interpro/api/entry/pfam/{acc}/"
    r = requests.get(url, timeout=30)
    if r.status_code != 200:
        return None
    return r.json().get("metadata", {})


def fetch_uniprot_has_domain(uniprot_acc: str, pfam_acc: str) -> bool:
    """Verify the example UniProt protein actually carries this Pfam."""
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_acc}.json"
    r = requests.get(url, timeout=30)
    if r.status_code != 200:
        return False
    data = r.json()
    xrefs = data.get("uniProtKBCrossReferences", [])
    return any(x.get("database") == "Pfam" and x.get("id") == pfam_acc for x in xrefs)


def validate_entry(entry: dict, strict: bool) -> list[str]:
    """Return list of issue strings; empty list means entry passes."""
    issues: list[str] = []
    acc = entry["accession"]
    meta = fetch_interpro_pfam(acc)
    if meta is None:
        return [f"{acc}: NOT FOUND in InterPro"]

    # Name/description coherence check
    ipr_name = (meta.get("name", {}) or {}).get("name", "")
    ipr_short = (meta.get("name", {}) or {}).get("short", "")
    declared_name = entry.get("name", "")
    declared_ipr = entry.get("interpro_name", "")

    # Soft check: declared name should share keywords with InterPro name
    if declared_ipr and declared_ipr.lower() != ipr_name.lower():
        issues.append(
            f"{acc}: declared interpro_name '{declared_ipr}' != InterPro's '{ipr_name}'"
        )

    # Hard check: InterPro name or short name must contain at least one keyword from declared name
    declared_keywords = [w.lower() for w in declared_name.replace("_", " ").split()
                         if len(w) > 2]
    ipr_text = (ipr_name + " " + ipr_short).lower()
    if declared_keywords and not any(k in ipr_text for k in declared_keywords):
        issues.append(
            f"{acc}: declared name '{declared_name}' has no keyword overlap with "
            f"InterPro '{ipr_name}' (short: '{ipr_short}')"
        )

    # Example UniProt protein sanity check (only in strict mode - slow)
    if strict:
        for up_acc in entry.get("example_uniprot", []):
            if not fetch_uniprot_has_domain(up_acc, acc):
                issues.append(
                    f"{acc}: example protein {up_acc} does not carry this Pfam"
                )
            time.sleep(0.3)  # be kind to UniProt

    return issues


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--whitelist", type=Path, default=WHITELIST)
    p.add_argument("--strict", action="store_true",
                   help="Also verify example UniProt proteins carry the domain (slower)")
    args = p.parse_args()

    cfg = yaml.safe_load(args.whitelist.read_text())
    entries = cfg["domains"] + cfg.get("auxiliary", [])
    print(f"Validating {len(entries)} entries...")
    print(f"Strict mode: {args.strict}")
    print()

    total_issues = 0
    for entry in entries:
        issues = validate_entry(entry, strict=args.strict)
        if not issues:
            print(f"  OK   {entry['accession']:10s} {entry.get('name','?'):30s}")
        else:
            for iss in issues:
                print(f"  FAIL {iss}")
                total_issues += 1

    print(f"\nTotal issues: {total_issues}")
    return 1 if total_issues else 0


if __name__ == "__main__":
    sys.exit(main())
