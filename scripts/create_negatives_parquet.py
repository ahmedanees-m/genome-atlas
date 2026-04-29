"""Create targets_v2_with_negatives.parquet for ESM-2 re-embedding.

Reads the 500 negative-control accessions saved by add_negative_controls.py,
fetches their sequences from targets_v1.parquet, and concatenates with
targets_v2.parquet to produce a single 10,000-row input file for embed_esm2.py.

Usage (on VM):
    python3 scripts/create_negatives_parquet.py
    # or with explicit paths:
    python3 scripts/create_negatives_parquet.py \
        --accessions /data/graphs/negative_control_accessions.txt \
        --targets-v1 /data/processed/targets_v1.parquet \
        --targets-v2 /data/processed/targets_v2.parquet \
        --output     /data/processed/targets_v2_with_negatives.parquet
"""
import argparse
import os
from pathlib import Path

import duckdb
import pandas as pd


def main(
    accessions_path: Path,
    t1_path: Path,
    t2_path: Path,
    output_path: Path,
) -> None:
    print("=" * 60)
    print("Create targets_v2_with_negatives.parquet")
    print("=" * 60)

    # Read negative accessions
    neg_accs = [l.strip() for l in accessions_path.read_text().splitlines() if l.strip()]
    print(f"  Negative accessions: {len(neg_accs)}")

    con = duckdb.connect()

    # Load targets_v2 (confirmed genome editors)
    v2 = con.execute(f"SELECT * FROM read_parquet('{t1_path.as_posix()}')"
                     ).fetchdf()  # will override below; just get schema
    v2 = con.execute(f"SELECT * FROM read_parquet('{t2_path.as_posix()}')").fetchdf()
    print(f"  targets_v2 rows: {len(v2)}")

    # Fetch negative sequences from targets_v1
    # Build a values list for IN clause (safe: accessions are UniProt IDs, alphanumeric)
    accs_sql = ", ".join(f"'{a}'" for a in neg_accs)
    negatives = con.execute(f"""
        SELECT
            accession,
            'NEGATIVE_CONTROL'       AS primary_mechanism_bucket,
            'unreviewed'             AS reviewed,
            protein_name,
            organism_name,
            organism_id,
            length,
            xref_pfam,
            xref_pdb,
            xref_alphafolddb,
            sequence,
            lineage_ids,
            0.0                      AS curiosity_score
        FROM read_parquet('{t1_path.as_posix()}')
        WHERE accession IN ({accs_sql})
    """).fetchdf()
    print(f"  Negatives with sequences fetched: {len(negatives)}")

    if len(negatives) < len(neg_accs):
        missing = len(neg_accs) - len(negatives)
        print(f"  WARNING: {missing} negatives not found in targets_v1 (sequences missing)")

    # Concatenate — v2 first, negatives appended
    combined = pd.concat([v2, negatives], ignore_index=True)
    print(f"  Combined rows: {len(combined)}")
    assert len(combined) == len(v2) + len(negatives), "Row count mismatch!"

    output_path.parent.mkdir(parents=True, exist_ok=True)
    combined.to_parquet(str(output_path), compression="zstd", index=False)
    sz = output_path.stat().st_size / (1024 * 1024)
    print(f"  Saved: {output_path}  ({sz:.1f} MB)")
    print("\nDone. Next: python3 archive/embed_esm2.py --input <this file>")


if __name__ == "__main__":
    _base = Path(os.environ.get("PENSTACK_DATA", "/data"))
    ap = argparse.ArgumentParser()
    ap.add_argument("--accessions", type=Path,
                    default=_base / "graphs/negative_control_accessions.txt")
    ap.add_argument("--targets-v1", type=Path,
                    default=_base / "processed/targets_v1.parquet")
    ap.add_argument("--targets-v2", type=Path,
                    default=_base / "processed/targets_v2.parquet")
    ap.add_argument("--output",     type=Path,
                    default=_base / "processed/targets_v2_with_negatives.parquet")
    args = ap.parse_args()
    main(args.accessions, args.targets_v1, args.targets_v2, args.output)
