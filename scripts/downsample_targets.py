"""Downsample targets_v1.parquet to a curated ~9,500-protein atlas."""
import argparse
from pathlib import Path
import pandas as pd


def score_protein(r) -> int:
    s = 0
    # Structural evidence (highest priority)
    pdb = r.get("xref_pdb")
    if pd.notna(pdb) and str(pdb).strip() not in ("", "None", "nan"):
        s += 100

    af = r.get("xref_alphafolddb")
    if pd.notna(af) and str(af).strip() not in ("", "None", "nan"):
        s += 50

    # Domain diversity = composite systems
    pfams = r.get("xref_pfam")
    if pd.notna(pfams) and str(pfams).strip() not in ("", "None", "nan"):
        pfam_list = [p.strip() for p in str(pfams).split(";") if p.strip()]
        s += 30 * len(set(pfam_list))

    # Length sweet spot for therapeutic deliverability
    length = r.get("length", 0)
    if 200 <= length <= 1200:
        s += 20

    # Priority model organisms (CRISPR/CAST/Fanzor sources)
    priority_taxa = {83333, 1314, 1280, 9606, 287, 208964, 666, 1133849, 93061, 99287}
    if r.get("organism_id") in priority_taxa:
        s += 10

    return s


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--target-n", type=int, default=9500,
                   help="Total proteins in curated atlas")
    args = p.parse_args()

    df = pd.read_parquet(args.input)
    print(f"Input: {len(df):,} proteins")

    # Tier 1: All SwissProt-reviewed entries
    swissprot = df[df["reviewed"] == "reviewed"].copy()
    print(f"Tier 1 (SwissProt): {len(swissprot):,}")

    # Tier 2: Top-scoring TrEMBL to fill remaining slots
    trembl = df[df["reviewed"] != "reviewed"].copy()
    trembl["curiosity_score"] = trembl.apply(score_protein, axis=1)

    slots = max(0, args.target_n - len(swissprot))
    print(f"Slots for TrEMBL: {slots:,}")

    if slots > 0 and len(trembl) > 0:
        top_trembl = trembl.nlargest(slots, "curiosity_score")
        final = pd.concat([swissprot, top_trembl])
    else:
        final = swissprot

    final = final.drop_duplicates(subset=["accession"]).reset_index(drop=True)
    print(f"Final targets_v2: {len(final):,} proteins")
    print("By mechanism bucket:")
    print(final["primary_mechanism_bucket"].value_counts())
    print("By review status:")
    print(final["reviewed"].value_counts())
    print("With PDB:", final["xref_pdb"].notna().sum())
    print("With AlphaFold:", final["xref_alphafolddb"].notna().sum())

    args.output.parent.mkdir(parents=True, exist_ok=True)
    final.to_parquet(args.output, compression="zstd", compression_level=9)
    print(f"Saved -> {args.output}")


if __name__ == "__main__":
    main()
