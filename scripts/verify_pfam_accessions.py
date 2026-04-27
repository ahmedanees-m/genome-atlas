"""Independent verification of all Pfam accessions in whitelist against live APIs."""
import requests
import yaml
from pathlib import Path

WHITELIST = Path("genome_atlas/data/pfam_whitelist.yaml")

def check(acc: str, claimed_name: str):
    results = {}
    # InterPro API
    try:
        r = requests.get(f"https://www.ebi.ac.uk/interpro/api/entry/pfam/{acc}/", timeout=15)
        if r.status_code == 200:
            data = r.json().get("metadata", {})
            results["interpro"] = {
                "name": data.get("name", {}).get("name"),
                "proteins": data.get("counters", {}).get("proteins", 0),
                "type": data.get("type"),
            }
        else:
            results["interpro"] = {"error": f"HTTP {r.status_code}"}
    except Exception as e:
        results["interpro"] = {"error": str(e)}
    
    # Pfam website
    try:
        r = requests.get(f"https://pfam.xfam.org/family/{acc}", timeout=15, allow_redirects=False)
        results["pfam_status"] = r.status_code
    except Exception as e:
        results["pfam_status"] = str(e)
    
    return acc, claimed_name, results

def main():
    cfg = yaml.safe_load(WHITELIST.read_text())
    all_entries = cfg["domains"] + cfg.get("auxiliary", [])
    
    print("=" * 110)
    print(f"Verifying {len(all_entries)} Pfam accessions against InterPro & Pfam APIs")
    print("=" * 110)
    print(f"{'Acc':<10} {'Claimed Name':<35} {'API Name':<40} {'Proteins':>10} {'Status':<10}")
    print("-" * 110)
    
    issues = []
    for entry in all_entries:
        acc = entry["accession"]
        claimed = entry["name"]
        _, _, results = check(acc, claimed)
        ip = results.get("interpro", {})
        
        if "error" in ip:
            api_name = f"ERROR: {ip['error']}"
            proteins = "N/A"
            status = "FAIL"
            issues.append((acc, claimed, f"InterPro error: {ip['error']}"))
        else:
            api_name = ip.get("name", "N/A") or "N/A"
            proteins = ip.get("proteins", 0)
            
            # Check name match
            claimed_lc = claimed.lower().replace("_", " ").replace("-", " ")
            api_lc = api_name.lower().replace("_", " ").replace("-", " ")
            
            if claimed_lc in api_lc or api_lc in claimed_lc or \
               any(word in api_lc for word in claimed_lc.split() if len(word) > 3):
                status = "OK"
            else:
                status = "NAME_MISMATCH"
                issues.append((acc, claimed, f"Name mismatch: API says '{api_name}'"))
        
        print(f"{acc:<10} {claimed:<35} {str(api_name):<40} {str(proteins):>10} {status:<10}")
    
    print("=" * 110)
    if issues:
        print(f"\nISSUES FOUND ({len(issues)}):")
        for acc, claimed, issue in issues:
            print(f"  {acc} ({claimed}): {issue}")
    else:
        print("\nAll accessions verified successfully — names match API records.")
    print("=" * 110)
    
    return len(issues)

if __name__ == "__main__":
    import sys
    sys.exit(main())
