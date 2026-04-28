# genome-atlas

[![CI](https://github.com/ahmedanees-m/genome-atlas/workflows/CI/badge.svg)](https://github.com/ahmedanees-m/genome-atlas/actions)
[![codecov](https://codecov.io/gh/ahmedanees-m/genome-atlas/branch/main/graph/badge.svg)](https://codecov.io/gh/ahmedanees-m/genome-atlas)
[![PyPI](https://img.shields.io/pypi/v/genome-atlas.svg)](https://pypi.org/project/genome-atlas/)
[![Docs](https://img.shields.io/badge/docs-readthedocs-blue)](https://genome-atlas.readthedocs.io)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

A unified knowledge graph and embedding space for programmable genome-writing enzymes — CRISPR, CAST, bridge recombinases, Fanzor, and more — with cross-system selection decision support.

Part of the [PEN-STACK](https://github.com/ahmedanees-m/pen-stack) infrastructure for non-destructive genome engineering.

## Install

```bash
pip install genome-atlas
```

## Quickstart

```python
from genome_atlas.api import Atlas

atlas = Atlas(
    graph_path="atlas.gpickle",
    embeddings_path="embeddings.parquet",
    targets_path="targets.parquet",
)
recs = atlas.select_editor(
    cell_type="HEK293T",
    edit_type="insertion",
    cargo_size_bp=1500,
    delivery="AAV",
    top_k=5,
)
for r in recs:
    print(r.system, r.pen_score, r.dsb_free)
```

## CLI

```bash
genome-atlas query-system System_SpCas9
genome-atlas select --cell HEK293T --edit deletion --top-k 5
```

## Benchmarks

Primary benchmark: Protein→Domain (`HAS_DOMAIN`) link prediction, 20% hold-out test set.
GNN confidence intervals via Mann-Whitney SE; Node2Vec via 1000× bootstrap on test set.

| Model | AUROC | AUPRC |
|-------|-------|-------|
| GAT (1 head, residual) | 0.9705 [0.9446–0.9964] | 0.9421 |
| GraphSAGE | 0.9664 [0.9405–0.9923] | 0.9184 |
| Classical RBF (node classification) | 0.9331 | — |
| Quantum Kernel (node classification) | 0.8731 | — |
| Node2Vec | 0.8411 [0.8052–0.8342] | 0.8905 |

## Citation

See [CITATION.cff](CITATION.cff) or cite the preprint:
> Ahmed, A. (2026). *GENOME-ATLAS: A unified knowledge graph for programmable genome-writing enzymes.* bioRxiv.

## License

MIT — see [LICENSE](LICENSE).
