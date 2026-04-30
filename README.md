# GENOME-ATLAS

[![CI](https://github.com/ahmedanees-m/genome-atlas/workflows/CI/badge.svg)](https://github.com/ahmedanees-m/genome-atlas/actions)
[![codecov](https://codecov.io/gh/ahmedanees-m/genome-atlas/branch/main/graph/badge.svg)](https://codecov.io/gh/ahmedanees-m/genome-atlas)
[![PyPI](https://img.shields.io/pypi/v/genome-atlas.svg)](https://pypi.org/project/genome-atlas/)
[![Docs](https://img.shields.io/badge/docs-readthedocs-blue)](https://genome-atlas.readthedocs.io)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Version](https://img.shields.io/badge/version-0.6.0-green)](CHANGELOG.md)

**Part of [PEN-STACK](https://github.com/ahmedanees-m/pen-stack)** (Programmable Enzyme Networks — Systematic Tool for Atlas and Knowledge)

Heterogeneous GNN knowledge graph for DNA-modifying enzyme annotation. GENOME-ATLAS builds a knowledge graph of programmable genome-writing enzymes — CRISPR, CAST, bridge recombinases, Fanzor, and more — trains GraphSAGE and GAT graph neural networks for link prediction, and provides a decision-support API for editor selection.

---

## Architecture

```
Raw sources                Build pipeline               Downstream
────────────               ───────────────              ──────────
UniProt TSVs  ──┐
PDB structures ─┤  materialize_graph.py ──► atlas.gpickle ──► build_pyg_hetero()
Pfam whitelist ─┤                                                     │
AlphaFold DB  ──┘                                                     ▼
                            ESM-2 (esm2_t30_150M_UR50D)   HeteroData (PyG)
                            embed_esm2.py                       │         │
                            640-dim per protein                 │         │
                                   │                     GraphSAGE   GAT
                                   └───────────────────► train_gnn.py (--model sage/gat)
                                                              │
                           Node2Vec (inductive)               │
                           node2vec_inductive_v6.py           │
                                   │                          │
                                   └──────────────────────────┤
                                                              ▼
                                              bootstrap_cis_v6.py
                                              95% CIs (1000× bootstrap)
                                                              │
                                                              ▼
                                                    Atlas API / CLI
                                                    atlas.select_editor(...)
```

### Knowledge Graph (v0.6.0)

| Statistic | Value |
|-----------|-------|
| Nodes | 12,267 |
| Edges | 13,645 |
| Node types | 6 (Protein, Domain, Structure, Mechanism, RNA, System) |
| Edge types | 7 (HAS_DOMAIN, HAS_PROTEIN, STRUCTURE_OF, SIMILAR_TO, HAS_MECHANISM, HAS_RNA, PART_OF) |
| Protein embeddings | ESM-2 esm2_t30_150M_UR50D (640-dim) |

---

## Quick Start

### Install

```bash
pip install genome-atlas
# or from source:
pip install -e ".[dev]"
```

### Python API

```python
from genome_atlas.api import Atlas

atlas = Atlas(
    graph_path="atlas.gpickle",
    embeddings_path="embeddings.parquet",
    targets_path="targets.parquet",
)

# Query a specific system
result = atlas.query_system("SpCas9")
print(result["name"], result["mechanism_bucket"])

# Select best editor for a use case
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

### CLI

```bash
genome-atlas query-system System_SpCas9
genome-atlas select --cell HEK293T --edit deletion --top-k 5
genome-atlas --help
```

---

## Benchmark Results (v0.6.0)

**Task**: Protein→Domain (`HAS_DOMAIN`) link prediction, 10% hold-out test set (20% combined val+test).
**Split**: 80/10/10 train/val/test per edge type, seed=42, applied per edge type independently.
**Test size**: n=1908 per model (positive + matched negatives; type-consistent negative sampling).
**CIs**: 1000× bootstrap resampling on held-out test predictions (seed=42).

### Table 1 — Primary Benchmark (Inductive, link prediction)

| Model | AUROC | 95% CI | AUPRC |
|-------|-------|--------|-------|
| GraphSAGE | **0.9707** | [0.9618, 0.9790] | **0.9456** |
| GAT (4-head, residual) | 0.9685 | [0.9588, 0.9771] | 0.9325 |

GraphSAGE and GAT are statistically tied (AUROC delta = 0.0022, confidence intervals fully overlap). Neither model is definitively superior on this dataset.

### Node2Vec (excluded from Table 1)

| Variant | AUROC | 95% CI | AUPRC | Note |
|---------|-------|--------|-------|------|
| Node2Vec (inductive) | 0.9867 | [0.9806, 0.9921] | — | Topology only; no ESM-2 features |
| Node2Vec (transductive) | 0.9965 | — | — | Supplementary only; inflated |

Node2Vec is excluded from the primary table for two reasons:

1. **Inductive variant**: random walks use only graph topology (no ESM-2 protein features). The fair comparison with GNNs requires equivalent input modalities; topology-only models have an inherent structural advantage on this graph.
2. **Transductive variant**: random walks are trained on the full graph including test edges, making the AUROC incomparable to inductive GNN results. Reported in supplementary as an upper bound.

---

## Directory Structure

```
genome-atlas/
├── genome_atlas/              # Python package (importable)
│   ├── __init__.py
│   ├── _version.py            # Version string (0.6.0)
│   ├── api.py                 # Atlas API class
│   ├── cli.py                 # Click CLI entry point
│   ├── selection.py           # Editor selection engine
│   ├── graph/
│   │   ├── __init__.py
│   │   └── build.py           # build_pyg_hetero(), add_train_val_test_split()
│   ├── models/
│   │   ├── __init__.py
│   │   └── graphsage.py       # HeteroGNN (SAGEConv/GATConv) + LinkPredictor
│   └── utils/
│       ├── __init__.py
│       └── size.py
│
├── scripts/                   # Training and analysis scripts
│   ├── train_gnn.py           # GNN training (--model sage|gat)
│   ├── bootstrap_cis_v6.py    # 95% CIs for GNN + transductive Node2Vec
│   ├── node2vec_inductive_v6.py  # Inductive Node2Vec evaluation
│   ├── materialize_graph.py   # Build NetworkX gpickle from DuckDB
│   ├── create_train_gpickle.py   # Remove isolated nodes
│   ├── create_negatives_parquet.py
│   ├── build_edges_v8.py
│   ├── add_rna_nodes.py
│   ├── add_structural_edges.py
│   ├── add_negative_controls.py
│   └── figures/               # Figure generation scripts (fig1–fig5)
│
├── tests/
│   ├── conftest.py            # Shared fixtures (skip if atlas not present)
│   ├── test_placeholder.py    # Version smoke test
│   ├── unit/                  # Fast, no data required
│   │   ├── test_api.py
│   │   ├── test_cli.py
│   │   ├── test_selection.py
│   │   └── test_selection_scoring.py
│   ├── integration/           # Skipped unless atlas.gpickle on disk
│   │   └── test_api_graph.py
│   └── regression/            # Checks AUROC does not drift
│       └── test_benchmark_ci.py
│
├── notebooks/
│   ├── 01_validation_selection.py
│   ├── benchmark_results.json # Recorded baseline (used by regression tests)
│   └── validation_scenarios.yaml
│
├── docs/                      # Sphinx docs + manuscript
│   ├── figures/               # fig1–fig5 (PNG + PDF)
│   └── manuscript/
│       ├── main.tex
│       └── references.bib
│
├── data/
│   ├── raw/
│   │   ├── uniprot/           # PF*.tsv.gz downloads
│   │   └── pdb/               # *.pdb.gz structures
│   └── processed/
│       └── targets_v1.parquet
│
├── containers/                # Docker build contexts (empty — images pre-built on VM)
├── dist/                      # Built wheels (genome_atlas-0.5.1-*)
├── pyproject.toml
├── setup.py
├── CHANGELOG.md
├── CITATION.cff
└── LICENSE
```

---

## Docker Images

All training steps run inside Docker on the training VM (GPU: NVIDIA V100 16.7 GB).

| Image | Purpose | Key packages |
|-------|---------|--------------|
| `pen-stack/data:0.1.0` | Data ingestion, chown | pandas, pyarrow, duckdb |
| `pen-stack/graph:0.1.0` | GNN training, graph ops | torch, torch-geometric, networkx, node2vec, scikit-learn |
| `pen-stack/plm:0.1.0` | ESM-2 inference (GPU) | transformers, esm, torch+CUDA |

Images are pre-built on the VM. To rebuild from scratch, Dockerfiles are stored in the VM at `/home/anees_22phd0670/pen-stack/dockerfiles/`.

---

## Reproducibility — v0.6.0 Pipeline

All steps use fixed seeds throughout. The full pipeline is `scratch/pipeline_v6.sh` on the local G: drive, copied to the VM before execution.

### Seeds

| Component | Seed | Where set |
|-----------|------|-----------|
| Train/val/test split | 42 | `add_train_val_test_split(seed=42)` |
| Test-set negative sampling (GNN eval) | 42 | `evaluate(..., seed=42)` in `train_gnn.py` |
| Bootstrap resampling | 42 | `bootstrap_ci(..., seed=42)` |
| Node2Vec training | 42 | `_Node2Vec(..., seed=42)` |
| Node2Vec probe negatives (train) | 42 | `build_Xy(..., rng_seed=42)` |
| Node2Vec probe negatives (test) | 43 | `build_Xy(..., rng_seed=43)` |
| LogisticRegression probes | 42 | `LogisticRegression(random_state=42)` |

### Pipeline Steps (estimated wall-clock times on V100)

| Step | Script | Time | Docker image |
|------|--------|------|--------------|
| 0 | `git pull` | <1 min | — |
| 1 | `create_train_gpickle.py` | ~5 min | graph |
| 2 | `create_negatives_parquet.py` | ~5 min | data |
| 3A | `embed_esm2.py` (GPU, 10k proteins) | ~60 min | plm |
| 3B | `embed_node2vec.py` (CPU, 16 workers) | ~30 min | graph (parallel with 3A) |
| 4 | `train_gnn.py --model sage --epochs 100` | <1 min¹ | graph (GPU) |
| 5 | `train_gnn.py --model gat --epochs 300 --lr 0.0001` | <1 min¹ | graph (GPU) |
| 6 | `benchmark_link_prediction.py` | ~10 min | graph |
| 7 | `bootstrap_cis_v6.py` (1000× bootstrap) | ~30 min | graph |
| 8 | `chown -R 1005:1005` (fix Docker ownership) | <1 min | data |

**Total**: ~75 min end-to-end (dominated by ESM-2 inference in Step 3A; GNN training <1 min on this graph).

¹ This graph has 12,267 nodes and 13,645 edges — a small graph. On a V100 GPU, 100 epochs of GraphSAGE complete in ~7 seconds and 300 epochs of GAT in ~25 seconds. Larger graphs would scale accordingly.

### Exact Reproduction Commands

```bash
# On the VM — assumes repo is at $REPO and data at $DATA
export REPO=/home/anees_22phd0670/pen-stack/code/repos/genome-atlas
export DATA=/home/anees_22phd0670/pen-stack/data

# Step 4: GraphSAGE
docker run --rm --gpus all \
    -v "$DATA":/data -v "$REPO":/repo -e PYTHONPATH=/repo \
    pen-stack/graph:0.1.0 \
    python3 /repo/scripts/train_gnn.py \
        --model sage \
        --graph             /data/graphs/atlas_train.gpickle \
        --esm-embeddings    /data/embeddings/esm2_150M_v6.parquet \
        --output-embeddings /data/embeddings/graphsage_v6.parquet \
        --output-metrics    /data/embeddings/graphsage_v6_metrics.parquet \
        --output-test-preds /data/embeddings/graphsage_v6_test_preds.parquet \
        --epochs 100

# Step 5: GAT
docker run --rm --gpus all \
    -v "$DATA":/data -v "$REPO":/repo -e PYTHONPATH=/repo \
    pen-stack/graph:0.1.0 \
    python3 /repo/scripts/train_gnn.py \
        --model gat \
        --graph             /data/graphs/atlas_train.gpickle \
        --esm-embeddings    /data/embeddings/esm2_150M_v6.parquet \
        --output-embeddings /data/embeddings/gat_v6.parquet \
        --output-metrics    /data/embeddings/gat_v6_metrics.parquet \
        --output-test-preds /data/embeddings/gat_v6_test_preds.parquet \
        --epochs 300 --lr 0.0001

# Step 7: Bootstrap CIs
docker run --rm \
    -v "$DATA":/data -v "$REPO":/repo -e PYTHONPATH=/repo \
    pen-stack/graph:0.1.0 \
    python3 /repo/scripts/bootstrap_cis_v6.py \
        --graph           /data/graphs/atlas_train.gpickle \
        --node2vec        /data/embeddings/node2vec_v6.parquet \
        --graphsage-preds /data/embeddings/graphsage_v6_test_preds.parquet \
        --gat-preds       /data/embeddings/gat_v6_test_preds.parquet \
        --output          /data/embeddings/bootstrap_cis_v6.parquet
```

Full pipeline: `nohup bash scratch/pipeline_v6.sh > ~/pen-stack/logs/pipeline_v6.log 2>&1 &`

---

## Running Tests

```bash
# Unit tests only (no data files required)
pytest tests/unit/ tests/test_placeholder.py -v

# All tests (integration and regression tests skip gracefully if atlas not on disk)
pytest -v

# With coverage
pytest --cov=genome_atlas --cov-report=term-missing
```

The unit test suite passes locally without any VM data. Integration and regression tests are marked `skipif` and skip cleanly in CI.

**Structural note**: `tests/unit/`, `tests/integration/`, and `tests/regression/` do not have `__init__.py` files. pytest discovers them correctly via `testpaths = ["tests"]` in `pyproject.toml`, but adding `__init__.py` to each subdirectory would make them importable as packages if needed.

---

## Key Findings

1. **GNNs are statistically tied**: GraphSAGE (AUROC 0.9707) and GAT (AUROC 0.9685) have overlapping 95% bootstrap CIs (delta = 0.0022). Neither is definitively superior; both achieve high link-prediction performance on the Protein→Domain task.

2. **Graph topology is the dominant signal**: Node2Vec (inductive, topology only, AUROC 0.9867) outperforms GNNs that use ESM-2 sequence embeddings (GraphSAGE 0.9707, GAT 0.9685). For Protein→Domain prediction on this graph, the connectivity structure is more informative than sequence features alone. ESM-2 embeddings provide complementary signal within the GNN message-passing framework but do not surpass the topology-only baseline.

3. **Transductive Node2Vec is inflated**: Node2Vec with full-graph random walks (AUROC 0.9965) is not comparable to inductive GNN results because test edges participate in training. It serves only as a topology upper-bound in supplementary material.

4. **GAT residual connections prevent embedding collapse**: Without residual connections, GAT produces only 4 unique embeddings across 10,000 Protein nodes because most nodes receive zero messages from their System neighbours. The residual (`x_dict[nt] = new_x[nt] + x`) restores ESM-2 diversity. GraphSAGE does not require this fix (self-features are concatenated internally by SAGEConv).

5. **Negative sampling is type-consistent**: All link prediction negatives are sampled from (src_type, dst_type) pools matching the positive edge type. This prevents the trivial classification pattern where a model learns "is the destination a Domain node?" rather than learning embedding similarity.

---

## Citation

See [CITATION.cff](CITATION.cff) for the machine-readable citation. Preprint forthcoming on bioRxiv.

---

## License

MIT — see [LICENSE](LICENSE).
