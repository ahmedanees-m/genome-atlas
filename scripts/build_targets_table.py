"""Merge Pfam-sliced UniProt TSVs into a deduplicated, filtered targets table."""
from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path

import pandas as pd
import yaml
from rich.progress import Progress, SpinnerColumn, BarColumn, TextColumn, TimeElapsedColumn


def load_whitelist(path: Path) -> dict:
    """Return {accession: mechanism_bucket} mapping from whitelist."""
    cfg = yaml.safe_load(path.read_text())
    return {
        e["accession"]: e["mechanism_bucket"]
        for e in cfg["domains"]
    }


def parse_xref_pfam(val: str | float) -> list[str]:
    """Parse semicolon-separated Pfam accessions from xref_pfam field."""
    if pd.isna(val) or not str(val).strip():
        return []
    return [x.strip().split(";")[0].split(".")[0] for x in str(val).split(";") if x.strip()]


def pct_x(seq: str) -> float:
    """Return percentage of X residues in sequence."""
    if not seq:
        return 100.0
    return (seq.upper().count("X") / len(seq)) * 100


def chunk_filter(path: Path, whitelist_accs: set[str], chunksize: int = 50_000):
    """Yield filtered DataFrame chunks from a gzipped TSV."""
    for chunk in pd.read_csv(
        path,
        sep="\t",
        compression="gzip",
        chunksize=chunksize,
        low_memory=False,
        on_bad_lines="skip",
    ):
        # Standardise column names (UniProt TSV field names -> our names)
        col_map = {
            "entry": "accession",
            "organism": "organism_name",
            "organism (id)": "organism_id",
            "pfam": "xref_pfam",
            "pdb": "xref_pdb",
            "alphafolddb": "xref_alphafolddb",
            "protein names": "protein_name",
            "taxonomic lineage (ids)": "lineage_ids",
        }
        chunk.columns = [c.lower().strip() for c in chunk.columns]
        chunk = chunk.rename(columns=col_map)

        # Ensure required columns exist
        required = {"accession", "length", "sequence", "reviewed", "xref_pfam"}
        missing = required - set(chunk.columns)
        if missing:
            # Some files may have slightly different columns; skip if critical ones missing
            if {"accession", "sequence"} & missing:
                continue
            for col in missing:
                chunk[col] = None

        # Filter by length
        chunk["length"] = pd.to_numeric(chunk["length"], errors="coerce")
        chunk = chunk[(chunk["length"] >= 100) & (chunk["length"] <= 2000)]
        if chunk.empty:
            continue

        # Filter by X residue content
        chunk["pct_x"] = chunk["sequence"].astype(str).apply(pct_x)
        chunk = chunk[chunk["pct_x"] < 50]
        if chunk.empty:
            continue

        # Filter: must have at least one whitelisted Pfam
        chunk["pfam_list"] = chunk["xref_pfam"].apply(parse_xref_pfam)
        chunk = chunk[chunk["pfam_list"].apply(lambda lst: any(a in whitelist_accs for a in lst))]
        if chunk.empty:
            continue

        yield chunk


def deduplicate(best: dict, chunk: pd.DataFrame) -> dict:
    """Update best-row dict with chunk, preferring reviewed entries."""
    for _, row in chunk.iterrows():
        acc = str(row["accession"]).strip()
        if not acc or acc == "nan":
            continue
        is_reviewed = str(row.get("reviewed", "")).strip().upper() == "REVIEWED"
        if acc not in best:
            best[acc] = row
        else:
            existing_reviewed = str(best[acc].get("reviewed", "")).strip().upper() == "REVIEWED"
            if is_reviewed and not existing_reviewed:
                best[acc] = row
    return best


def assign_primary_bucket(row, whitelist_map: dict[str, str]) -> str | None:
    """Assign the first matching mechanism bucket from whitelist."""
    pfams = row.get("pfam_list", [])
    for p in pfams:
        if p in whitelist_map:
            return whitelist_map[p]
    return None


def build_targets(
    raw_dir: Path,
    whitelist_path: Path,
    output_path: Path,
) -> int:
    whitelist_map = load_whitelist(whitelist_path)
    whitelist_accs = set(whitelist_map.keys())

    files = sorted(raw_dir.glob("*.tsv.gz"))
    if not files:
        print("ERROR: No .tsv.gz files found in", raw_dir)
        return 1

    best_rows: dict[str, pd.Series] = {}

    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
        TimeElapsedColumn(),
    ) as progress:
        task = progress.add_task("Processing TSVs", total=len(files))
        for fpath in files:
            progress.update(task, description=f"Processing {fpath.name}")
            try:
                for chunk in chunk_filter(fpath, whitelist_accs):
                    best_rows = deduplicate(best_rows, chunk)
            except Exception as e:
                print(f"WARNING: error processing {fpath.name}: {e}")
            progress.advance(task)

    print(f"Total unique accessions after deduplication: {len(best_rows)}")

    if not best_rows:
        print("ERROR: No rows passed filters.")
        return 1

    # Build final dataframe
    df = pd.DataFrame(best_rows.values())

    # Assign mechanism bucket
    df["primary_mechanism_bucket"] = df.apply(lambda r: assign_primary_bucket(r, whitelist_map), axis=1)

    # Validate all rows have a bucket
    missing_bucket = df["primary_mechanism_bucket"].isna().sum()
    if missing_bucket:
        print(f"WARNING: {missing_bucket} rows missing mechanism bucket; dropping them.")
        df = df[df["primary_mechanism_bucket"].notna()]

    # Drop helper columns not needed in output
    drop_cols = [c for c in ["pct_x", "pfam_list"] if c in df.columns]
    if drop_cols:
        df = df.drop(columns=drop_cols)

    # Reorder: put key cols first
    first_cols = ["accession", "primary_mechanism_bucket", "reviewed", "protein_name",
                  "organism_name", "organism_id", "length", "xref_pfam", "xref_pdb",
                  "xref_alphafolddb", "sequence", "lineage_ids"]
    first_cols = [c for c in first_cols if c in df.columns]
    other_cols = [c for c in df.columns if c not in first_cols]
    df = df[first_cols + other_cols]

    # Write parquet with zstd compression
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(output_path, compression="zstd", index=False)

    print(f"Wrote {len(df)} rows to {output_path}")
    print(f"  Reviewed: {(df['reviewed'].astype(str).str.upper() == 'REVIEWED').sum()}")
    print(f"  Mechanism bucket counts:")
    for bucket, count in df["primary_mechanism_bucket"].value_counts().items():
        print(f"    {bucket}: {count}")

    return 0


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--raw-dir", type=Path, required=True)
    p.add_argument("--whitelist", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    return build_targets(args.raw_dir, args.whitelist, args.output)


if __name__ == "__main__":
    sys.exit(main())
