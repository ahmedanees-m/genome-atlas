#!/usr/bin/env python3
import pandas as pd, numpy as np, torch, os
from pathlib import Path

DATA = Path('/data')
ESM_PATH = DATA / 'embeddings/esm2_150M_v6.parquet'
TARGETS_PATH = DATA / 'processed/targets_v2_with_negatives.parquet'

esm_df = pd.read_parquet(ESM_PATH)
print(f'ESM-2 rows before: {len(esm_df)}')

if 'D2TGM5' in esm_df['accession'].values:
    print('D2TGM5 already in ESM-2, nothing to do')
    exit(0)

tgt = pd.read_parquet(TARGETS_PATH, columns=['accession', 'sequence'])
row = tgt[tgt['accession'] == 'D2TGM5']
sequence = row['sequence'].iloc[0]
print(f'Sequence length: {len(sequence)}')

import esm as esm_lib
print('Loading ESM-2 150M...')
model, alphabet = esm_lib.pretrained.esm2_t30_150M_UR50D()
model.eval()
device = 'cuda' if torch.cuda.is_available() else 'cpu'
model = model.to(device)
print(f'Device: {device}')

batch_converter = alphabet.get_batch_converter()
data = [('D2TGM5', sequence)]
_, _, batch_tokens = batch_converter(data)
batch_tokens = batch_tokens.to(device)

with torch.no_grad():
    results = model(batch_tokens, repr_layers=[30], return_contacts=False)

token_rep = results['representations'][30]
embedding = token_rep[0, 1:len(sequence)+1].mean(0).cpu().numpy()
print(f'Embedding dim: {embedding.shape[0]}')

new_row = pd.DataFrame([{'accession': 'D2TGM5', 'embedding': embedding.tolist(), 'seq_length': len(sequence)}])
updated = pd.concat([esm_df, new_row], ignore_index=True)
updated.to_parquet(ESM_PATH, index=False, compression='zstd')
print(f'ESM-2 rows after: {len(updated)}')
print('D2TGM5 embedding added successfully')
