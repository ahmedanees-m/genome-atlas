#!/usr/bin/env python3
"""Compute ESM-2 embeddings for every protein in a targets parquet.

Loads esm2_t30_150M_UR50D, runs each sequence through the model, and mean-pools
the layer-30 residue representations into a single 640-dim vector per protein.
The output parquet (accession, embedding, seq_length) is what the GNN uses to
initialise Protein node features.

Usage:
    python scripts/embed_esm2.py \
        --targets ~/pen-stack/data/processed/targets_v2_with_negatives.parquet \
        --output  ~/pen-stack/data/embeddings/esm2_150M_v6.parquet

Needs a GPU for anything beyond a few hundred proteins. Falls back to CPU.
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import torch

# ESM-2 caps out around 1022 residues plus the BOS/EOS tokens. Longer proteins
# are truncated; the catalytic domains we care about sit well under this.
MAX_LEN = 1022
REPR_LAYER = 30


def embed(targets_path: Path, output_path: Path, batch_size: int):
    import esm

    df = pd.read_parquet(targets_path, columns=["accession", "sequence"])
    df = df.dropna(subset=["sequence"]).drop_duplicates("accession").reset_index(drop=True)
    print(f"Proteins to embed: {len(df)}")

    model, alphabet = esm.pretrained.esm2_t30_150M_UR50D()
    model.eval()
    device = "cuda" if torch.cuda.is_available() else "cpu"
    model = model.to(device)
    batch_converter = alphabet.get_batch_converter()
    print(f"Device: {device}")

    records = []
    for start in range(0, len(df), batch_size):
        chunk = df.iloc[start:start + batch_size]
        data = [(r.accession, r.sequence[:MAX_LEN]) for r in chunk.itertuples()]
        _, _, tokens = batch_converter(data)
        tokens = tokens.to(device)

        with torch.no_grad():
            out = model(tokens, repr_layers=[REPR_LAYER], return_contacts=False)
        reps = out["representations"][REPR_LAYER]

        for i, (acc, seq) in enumerate(data):
            # Skip BOS (index 0) and any padding/EOS past the sequence length.
            vec = reps[i, 1:len(seq) + 1].mean(0).cpu().numpy()
            records.append({"accession": acc, "embedding": vec.tolist(), "seq_length": len(seq)})

        done = min(start + batch_size, len(df))
        print(f"  {done}/{len(df)}")

    out_df = pd.DataFrame(records)
    dim = len(out_df["embedding"].iloc[0]) if len(out_df) else 0
    output_path.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_parquet(output_path, index=False, compression="zstd")
    print(f"Wrote {len(out_df)} embeddings ({dim}-dim) to {output_path}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets", type=Path, required=True,
                    help="Parquet with 'accession' and 'sequence' columns")
    ap.add_argument("--output", type=Path, required=True,
                    help="Destination parquet (accession, embedding, seq_length)")
    ap.add_argument("--batch-size", type=int, default=8)
    args = ap.parse_args()
    embed(args.targets, args.output, args.batch_size)
