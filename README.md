# GENOME-ATLAS

[![CI](https://github.com/ahmedanees-m/genome-atlas/workflows/CI/badge.svg)](https://github.com/ahmedanees-m/genome-atlas/actions)
[![codecov](https://codecov.io/gh/ahmedanees-m/genome-atlas/branch/main/graph/badge.svg)](https://codecov.io/gh/ahmedanees-m/genome-atlas)
[![PyPI](https://img.shields.io/pypi/v/genome-atlas.svg)](https://pypi.org/project/genome-atlas/)
[![Docs](https://img.shields.io/badge/docs-readthedocs-blue)](https://genome-atlas.readthedocs.io)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Version](https://img.shields.io/badge/version-0.7.2-green)](CHANGELOG.md)

Part of [PEN-STACK](https://github.com/ahmedanees-m/pen-stack), a set of tools for programmable genome writing.

---

## What is GENOME-ATLAS?

GENOME-ATLAS is a **knowledge graph and machine-learning benchmark** for programmable genome-writing enzymes - the molecular tools scientists use to make precise, permanent changes to DNA.

In plain terms: it is a structured database that knows *which proteins belong to which editing system*, *which protein domains give them their function*, *what 3D structures they adopt*, and *what kind of edit each tool is best suited for*. On top of this database, it trains graph neural networks (GNNs) to predict new relationships - and provides a Python API that researchers can query to pick the right editor for their experiment.

Think of it as a smart, interconnected map of the genome-editing toolkit, built to be queried by both humans and machine-learning models.

---

## Why was it built?

Gene editing is advancing rapidly. Over the past decade, researchers have developed dozens of distinct molecular tools - CRISPR nucleases, base editors, prime editors, bridge recombinases, transposases, and more - each with different strengths, limitations, sizes, and ideal use cases.

**The problem**: choosing the right tool for a given experiment is non-trivial. The decision depends on many factors simultaneously - the cell type you are editing, whether a DNA break is acceptable, how much cargo you need to deliver, whether AAV packaging size is a constraint, and more. There is no single authoritative, machine-readable resource that connects all of this information in a way that supports computational reasoning.

**The goal of GENOME-ATLAS**:

1. **Build a unified knowledge graph** that links proteins, protein domains, 3D structures, editing mechanisms, RNA guides, and organisms into a single queryable graph for all major genome-writing systems.
2. **Benchmark graph neural networks** on this graph, demonstrating that GNNs can learn biologically meaningful relationships (e.g. predicting which protein domains a new enzyme is likely to carry).
3. **Provide a decision-support API** (`select_editor`) that researchers can use to receive ranked editor recommendations for a specific editing scenario.
4. **Ship a reproducible, versioned, citable data artifact** that downstream tools and analyses can build on directly.

---

## How does it work?

GENOME-ATLAS is built in three layers that work together:

### Layer 1 - The Knowledge Graph

The graph is a **heterogeneous network**: it has multiple node types (System, Protein, Domain, Structure, RNA, Organism) and multiple edge types connecting them. Each node stores metadata specific to its type, and each edge represents a biologically meaningful relationship.

```
           HAS_PROTEIN          HAS_DOMAIN
  System ────────────► Protein ──────────► Domain
    │                     │
    │ HAS_RNA              │ STRUCTURE_OF
    ▼                     ▼
   RNA               Structure ──────────► Structure
                               SIMILAR_TO
```

For example: `SpCas9` (System) **HAS_PROTEIN** `Q99ZW2` (Protein) **HAS_DOMAIN** `PF09650` (HNH nuclease domain). A structure entry for `7S7W` is linked back to the protein via **STRUCTURE_OF**, and similar structures are linked by **SIMILAR_TO**.

The graph is built from public databases (UniProt, RCSB PDB, AlphaFold DB, CRISPRCasDB) using a reproducible ingestion pipeline, stored in DuckDB, then exported as a NetworkX graph for training.

**Graph statistics (v0.7.1):**

| Property | Value |
|----------|-------|
| Total nodes | 13,401 |
| Training nodes (non-isolated) | 11,763 |
| Edges | 11,817 |
| Node types | 7 (System, Protein, Domain, Structure, RNA, Mechanism, Organism) |
| Primary edge types | 4 (HAS_DOMAIN, HAS_PROTEIN, STRUCTURE_OF, USES_MECHANISM) |
| Full edge types | 7 (+ HAS_RNA, PART_OF, SIMILAR_TO via `graph_view='full'`) |
| Foundational systems | 28 |
| Protein embeddings | 10,001 proteins at 640 dimensions (ESM-2) |

### Layer 2 - Graph Neural Networks

Each protein node is initialised with a **640-dimensional ESM-2 embedding** - a numerical fingerprint of its amino acid sequence, computed by a protein language model (`esm2_t30_150M_UR50D`, 150 million parameters). This captures sequence-level biochemistry before any graph training begins.

Two GNN architectures then learn to propagate information through the graph:

- **GraphSAGE** - samples and aggregates neighbourhood features; fast and robust on sparse graphs
- **GAT (Graph Attention Network)** - uses attention heads to weight the contribution of each neighbour differently; requires residual connections on this graph (see [Key Findings](#key-findings))

Both models are trained on a **link prediction task**: given a protein node and a domain node, predict whether a `HAS_DOMAIN` edge should exist between them. This mirrors a real biological question - *"does this protein carry this functional domain?"*

Training uses an 80/10/10 split (train/validation/test), all seeds fixed at 42 for full reproducibility.

### Layer 3 - The Editor Selection API

The `select_editor()` function combines the graph's structured metadata with a multi-factor scoring model (PEN-SCORE) to rank editors for a specific experimental scenario. It considers:

- Edit type (deletion, SNV, small insertion, large insertion)
- Cargo size and AAV packaging constraints
- Cell type and delivery method
- Whether a DNA double-strand break is acceptable
- Whether RNA guidance is required

The result is a ranked list of editors with per-factor scores, so researchers can understand *why* a particular tool was recommended.

---

## Quick Start

### Install

```bash
pip install genome-atlas
```

Or to install from source with development dependencies:

```bash
git clone https://github.com/ahmedanees-m/genome-atlas.git
cd genome-atlas
pip install -e ".[dev]"
```

To also install GNN training dependencies (requires PyTorch):

```bash
pip install -e ".[gnn]"
```

### Python API

```python
from genome_atlas import Atlas

# Load the knowledge graph and associated data
atlas = Atlas(
    graph_path="path/to/atlas.gpickle",
    embeddings_path="path/to/esm2_150M_v6.parquet",
    targets_path="path/to/targets_v2.parquet",
)

# --- Query a system ---
system = atlas.query_system("SpCas9")
print(system["name"], system["mechanism_bucket"])
# -> SpCas9  DSB_NUCLEASE

# --- List all systems filtered by mechanism ---
df = atlas.systems(mechanism_bucket="DSB_FREE_PRIME_EDITOR")
print(df[["name", "mechanism_bucket"]])

# --- Find all proteins in the graph for a system ---
protein = atlas.query_protein("Q99ZW2")

# --- Get editor recommendations for an experiment ---
recs = atlas.select_editor(
    cell_type="HEK293T",
    edit_type="insertion",
    cargo_size_bp=1500,
    delivery="AAV",
    prefer_dsb_free=True,
    top_k=5,
)

for r in recs:
    print(f"{r.system:20s}  score={r.pen_score:.3f}  dsb_free={r.dsb_free}")
```

### CLI

```bash
# Query a system by node ID
genome-atlas query-system System_SpCas9

# Get editor recommendations from the command line
genome-atlas select --cell HEK293T --edit deletion --top-k 5

# Get help
genome-atlas --help
```

### Load system metadata (no graph file required)

```python
from genome_atlas import load_systems, resolve_system_name

# Get all 28 foundational systems as structured dataclasses
systems = load_systems()
iscro4 = systems["ISCro4"]
print(iscro4.uniprot)        # D2TGM5
print(iscro4.mechanism_bucket)  # DSB_FREE_BRIDGE_RECOMBINASE

# Resolve a deprecated name (e.g. from an old paper) to the canonical one
canonical = resolve_system_name("IS622")  # -> "ISCro4" + DeprecationWarning
```

---

## What's Covered - The 28 Foundational Systems

GENOME-ATLAS v0.7.1 includes 28 curated systems representing **57% of a reference set of 49 well-characterised therapeutic genome editors**. Only systems with published mechanistic characterisation, structural evidence, and demonstrated therapeutic relevance are included - this is a deliberate quality filter, not an omission.

Coverage by mechanism class:

| Mechanism | Systems included | Coverage |
|-----------|-----------------|----------|
| DSB nuclease (Cas9, Cas12) | SpCas9, SaCas9, CjCas9, NmCas9, AsCas12a, ... | 100% of well-characterised families |
| Base editor | ABE8e, BE4max, ... | 100% of mechanism classes |
| Prime editor | PE2, PE3, twinPE | 100% |
| Bridge recombinase | IS621, ISCro4 | 66.7% (2/3); ISCro4 added in v0.7.1 |
| Serine integrase | Bxb1, phiC31, TP901-1 | 100% |
| Transposase / CAST | CAST-I-F, ShCAST | 100% of families |
| Fanzor / OMEGA | SpuFz1, ... | 100% |

**ISCro4** (formerly called IS622 in preprints; UniProt D2TGM5; *Citrobacter rodentium* ICC168) was added in v0.7.1 as the 28th system. It is an IS110-family bridge recombinase achieving ~20% insertion efficiency in HEK293T cells (Perry 2025 *bioRxiv*; Pelea 2026 *Science* DOI:10.1126/science.adz1884). See [DATA_PROVENANCE.md](DATA_PROVENANCE.md) for the full system list with references.

---

## Benchmark Results

The primary benchmark tests how well GNNs can predict **Protein->Domain** relationships (`HAS_DOMAIN` edges) on the held-out test set - a biologically meaningful task that checks whether the model can infer what functional domain a protein carries.

**Task setup:**
- Split: 80% train / 10% validation / 10% test, per edge type, seed = 42
- Test set: n = 1,908 samples (positive edges + type-consistent negative samples)
- Confidence intervals: 1,000x bootstrap resampling on held-out predictions, seed = 42

### Table 1 - Primary Benchmark (HAS_DOMAIN link prediction)

| Model | AUROC | 95% CI | AUPRC |
|-------|-------|--------|-------|
| GraphSAGE | **0.9714** | [0.9625, 0.9797] | **0.9451** |
| GAT (4-head, residual) | 0.9690 | [0.9590, 0.9778] | 0.9331 |

Both models remain within v0.7.0 confidence intervals. GraphSAGE and GAT are **statistically tied** (AUROC delta = 0.0024; CIs substantially overlap).

### Secondary task - Structure->Protein (STRUCTURE_OF)

| Model | AUROC | 95% CI |
|-------|-------|--------|
| GraphSAGE | 0.9971 | [0.9913, 1.0000] |
| GAT | 0.9971 | [0.9913, 1.0000] |

### Node2Vec (topology-only baseline, not in Table 1)

| Variant | AUROC | 95% CI | Note |
|---------|-------|--------|------|
| Node2Vec (inductive) | 0.9890 | [0.9825, 0.9940] | No ESM-2 features; topology only |
| Node2Vec (transductive) | 0.9965 | [0.9924, 0.9997] | Supplementary only - inflated; test edges seen during training |

Node2Vec is excluded from the primary table because (a) it uses only graph topology without protein sequence features, making it not comparable to GNNs on equal footing, and (b) the transductive variant is an upper bound, not an inductive result.

Authoritative CI file: `reproduction/bootstrap_ci_v7.json` (Parquet format; read with `pd.read_parquet()`).

---

## Editor Selection Validation

`select_editor()` was validated against **10 published therapeutic editing scenarios** from peer-reviewed literature (2016-2025). For each scenario, the ground truth is the editor that was actually used in the published experiment. GENOME-ATLAS was asked to rank candidates without seeing the answer in advance.

**Result: 7 out of 10 correct (70%) in top-3 recommendations.**

| Scenario | Target / Disease | Published Editor | ATLAS Rank | Correct? |
|----------|-----------------|-----------------|------------|----------|
| 1 | Sickle cell, BCL11A enhancer | SpCas9 | #1 | Yes |
| 2 | DMD exon skipping (AAV) | SpuFz1 Fanzor | n/a | No |
| 3 | CAR-T 5 kb knockin | CAST-I-F evoCAST | n/a | No |
| 4 | Point mutation SNV (AAV) | PE2 prime editor | n/a | No |
| 5 | Megabase rearrangement | IS621 bridge recombinase | #3 | Yes |
| 6 | Liver PCSK9 knockdown | SpCas9 | #1 | Yes |
| 7 | Retinal CEP290 restoration | PE2 prime editor | #1 | Yes |
| 8 | TRAC knockout (T cells) | SpCas9 | #1 | Yes |
| 9 | Landing pad integration | Bxb1 integrase | #3 | Yes |
| 10 | Compact AAV deletion | Cas12f | #3 | Yes |

Each of the three misses points to a specific gap in the current selection heuristic. Full details: [VALIDATION.md](VALIDATION.md).

---

## Key Findings

1. **GraphSAGE and GAT are statistically tied.** AUROC delta of 0.0024 with substantially overlapping confidence intervals. Both achieve strong link-prediction performance, so architecture choice is secondary to graph quality and initialisation.

2. **Graph topology is the dominant signal.** Node2Vec (topology-only, AUROC 0.9890) outperforms both GNNs that use 640-dimensional ESM-2 embeddings (0.9714 and 0.9690). For Protein->Domain prediction on this graph, *who your neighbours are* tells you more than *what your sequence looks like* - the connectivity structure is highly informative. ESM-2 embeddings provide complementary signal but do not surpass topology.

3. **GAT requires residual connections on sparse bipartite graphs.** Without them, GAT produces only 4 unique embeddings across ~9,500 Protein nodes. The reason: most Protein nodes have no System neighbours, so they receive zero attention messages and collapse to the same output. The fix (`x_dict[nt] = new_x[nt] + x`) preserves ESM-2 diversity by adding back the original features. GraphSAGE does not need this fix because SAGEConv concatenates self-features internally.

4. **Selection accuracy is 70% top-3 on published scenarios.** The three misses each expose a specific boundary: DSB penalisation is too aggressive for moderate-cargo AAV scenarios, wild-type and evolved CAST variants score identically without efficiency data, and component size conflates total system size with per-AAV payload. These inform the PEN-SCORE v2 design.

---

## How the Pipeline Works End to End

```
 Data sources            Build                  Train                  Use
 ────────────            ─────                  ─────                  ───
 UniProt TSVs ─┐
 PDB structs  ─┤──► DuckDB ──► atlas.gpickle ──► HeteroData (PyG) ──► GraphSAGE/GAT
 Pfam list    ─┤    (atlas.duckdb)              (build_pyg_hetero)        │
 AlphaFold DB ─┘                                                          │ 128-dim
 CRISPRCasDB ──┘                                                          ▼
                                                              bootstrap_cis_v6.py
 ESM-2 model ──────────────────────────────────────────────► 640-dim embeddings
 (GPU, 150M params)      esm2_150M_v6.parquet                    (init for GNN)
                                                                          │
 foundational_systems.yaml ─────────────────────────────────────► Atlas API
 (28 curated systems)                                             select_editor()
```

**Step-by-step:**

1. **Ingest** - `ingest_to_duckdb_v2.py` reads UniProt TSVs, AlphaFold metadata, CRISPRCasDB entries, and the curated `foundational_systems.yaml` into a structured DuckDB database with tables for proteins, domains, structures, mechanisms, systems, and edges.

2. **Materialise** - `materialize_graph.py` queries DuckDB and builds a NetworkX `MultiDiGraph` with typed nodes and edges, saved as `atlas.gpickle`.

3. **Prune** - `create_train_gpickle.py` removes isolated nodes (nodes with no edges) to produce `atlas_train.gpickle`, which is used for GNN message-passing.

4. **Embed** - `embed_esm2.py` runs the ESM-2 protein language model (GPU) over all protein sequences, producing 640-dimensional embeddings saved to `esm2_150M_v6.parquet`.

5. **Train** - `train_gnn.py --model sage` and `--model gat` load the training graph and ESM-2 embeddings, build a PyTorch Geometric `HeteroData` object, and train for 100 epochs. Best model by validation AUROC is kept.

6. **Evaluate** - `bootstrap_cis_v6.py` takes held-out test predictions and computes 1,000x bootstrap confidence intervals for AUROC and AUPRC.

7. **Use** - The `Atlas` API wraps the graph and embeddings with query and selection methods for downstream research use.

All training runs inside Docker containers on a NVIDIA V100 GPU. See [Reproducibility](#reproducibility) for exact commands.

---

## Installation & Requirements

**Minimum (API only - no GPU required):**

```bash
pip install genome-atlas
```

Core dependencies installed automatically: `networkx`, `pandas`, `numpy`, `pyarrow`, `pyyaml`, `scikit-learn`, `click`.

**For GNN training (requires PyTorch + GPU recommended):**

```bash
pip install "genome-atlas[gnn]"
# installs: torch>=2.0, torch-geometric>=2.3
```

**For development (tests, linting, docs):**

```bash
pip install "genome-atlas[dev]"
# installs: pytest, pytest-cov, black, flake8, sphinx
```

**Python version:** 3.10 or 3.11 (tested in CI on both).

**Data files** (not included in the pip package; provided as raw data files):
- `graphs/atlas.gpickle` - the full knowledge graph
- `embeddings/esm2_150M_v6.parquet` - protein embeddings (10,001 rows, 640-dim)
- `graphs/atlas.duckdb` - source relational database

---

## Running Tests

```bash
# Unit tests only - fast, no data files needed
pytest tests/unit/ tests/test_placeholder.py -v

# Full suite - integration and regression tests skip automatically if atlas.gpickle is absent
pytest -v

# With coverage report
pytest --cov=genome_atlas --cov-report=term-missing
```

**Test suite:** 133 tests pass in CI (Python 3.10 and 3.11), 5 skipped (torch-dependent GNN tests). Coverage: **97%** across all core modules. Torch-dependent files (`models/`, `graph/build.py`, `graph/view.py`) are omitted from CI coverage since they require the optional `[gnn]` extra - this is standard practice for optional-dependency code and does not affect the core API coverage.

---

## Reproducibility

All training steps use fixed random seeds throughout.

| Component | Seed |
|-----------|------|
| Train/val/test split | 42 |
| GNN test-set negative sampling | 42 |
| Bootstrap resampling | 42 |
| Node2Vec training | 42 |
| Node2Vec probe negatives (train) | 42 |
| Node2Vec probe negatives (test) | 43 |
| LogisticRegression probes | 42 |

The build runs on Python 3.10 or 3.11 with the `[gnn]` extra installed. ESM-2
inference needs a GPU; the graph and API steps run on CPU. The ordered rebuild
sequence is documented in [UPDATE_STRATEGY.md](UPDATE_STRATEGY.md).

Approximate runtimes on a single GPU: GraphSAGE 100 epochs ~7 s; GAT 100 epochs
~10 s; Node2Vec ~25 min. ESM-2 inference dominates a full rebuild.

---

## Project Structure

```
genome-atlas/
├── genome_atlas/                      # Installable Python package
│   ├── __init__.py                    # Public API: Atlas, get_graph, load_systems
│   ├── _version.py                    # Version string ("0.7.2")
│   ├── api.py                         # Atlas class - query, select, embed
│   ├── cli.py                         # Command-line interface (Click)
│   ├── selection.py                   # PEN-SCORE editor selection engine
│   ├── systems.py                     # load_systems(), SystemEntry, resolve_system_name
│   ├── data/
│   │   ├── foundational_systems.yaml  # 28 curated systems (ISCro4 canonical)
│   │   └── pfam_whitelist.yaml        # 18 Pfam family whitelist
│   ├── graph/
│   │   ├── build.py                   # PyG HeteroData builder + train/val/test split
│   │   └── view.py                    # get_graph(graph_view='primary'|'full')
│   ├── models/
│   │   └── graphsage.py               # HeteroGNN (SAGEConv / GATConv) + LinkPredictor
│   └── utils/
│       └── size.py                    # system_total_size_aa helper
│
├── scripts/                           # Data pipeline and training scripts
│   ├── train_gnn.py                   # Train GraphSAGE or GAT (--model sage|gat)
│   ├── embed_esm2.py                  # ESM-2 protein embedding (GPU)
│   ├── materialize_graph.py           # DuckDB -> NetworkX gpickle
│   ├── ingest_to_duckdb_v2.py         # Raw data -> DuckDB
│   ├── create_train_gpickle.py        # Remove isolated nodes
│   ├── bootstrap_cis_v6.py            # 95% bootstrap CIs (1,000x)
│   ├── node2vec_inductive_v6.py       # Inductive Node2Vec baseline
│   ├── audit_coverage.py              # Coverage vs 49-system reference
│   └── add_structural_edges.py        # Foldseek TM-score similarity edges
│
├── tests/
│   ├── unit/                          # Fast - no data files required
│   │   ├── test_aliases.py            # ISCro4/IS622 alias resolution (22 tests)
│   │   ├── test_api.py                # Atlas API basics
│   │   ├── test_api_extended.py       # Full graph traversal coverage (34 tests)
│   │   ├── test_cli.py                # CLI basics
│   │   ├── test_cli_extended.py       # CLI command integration (7 tests)
│   │   ├── test_graph_views.py        # get_graph views (18 tests; skipped without torch)
│   │   ├── test_init.py               # Package-level imports (6 tests)
│   │   ├── test_selection.py          # Editor selection basics
│   │   ├── test_selection_extended.py # All scoring branches (32 tests)
│   │   ├── test_selection_scoring.py  # PEN-SCORE scoring
│   │   └── test_utils.py              # Utility functions (5 tests)
│   ├── integration/                   # Skipped if atlas.gpickle not on disk
│   │   └── test_api_graph.py
│   └── regression/                    # AUROC drift check vs benchmark_results.json
│       └── test_benchmark_ci.py
│
├── notebooks/
│   ├── 01_validation_selection.py     # 10-scenario validation script
│   ├── benchmark_results.json         # Authoritative v0.7.1 CI file (regression baseline)
│   └── validation_scenarios.yaml      # The 10 validation scenarios with ground truth
│
├── reproduction/
│   └── bootstrap_ci_v7.json           # 95% CIs, all models (Parquet; use pd.read_parquet)
│
├── docs/                              # Sphinx documentation
│   ├── GRAPH_SCHEMA.md                # Node/edge type reference
│   ├── conf.py                        # release = "0.7.2"
│   └── ...
│
├── pyproject.toml                     # Package config; version = "0.7.2"
├── CHANGELOG.md                       # Version history
├── DATA_PROVENANCE.md                 # Source databases, system list, coverage stats
├── UPDATE_STRATEGY.md                 # Versioning cadence; v0.8/v1.0 roadmap
├── VALIDATION.md                      # Full validation report (10 scenarios)
└── LICENSE                            # MIT
```

---

## Versioning

| Version | What changed |
|---------|-------------|
| **v0.7.2** (2026-05-25) | ISCro4 canonical naming; `load_systems()` API; 141 unit tests; 97% coverage |
| **v0.7.1** (2026-05-23) | ISCro4 (28th system) added; ESM-2 embedding for D2TGM5; graph rebuilt |
| **v0.7.0** | 27 systems; 11 new systems added (SaCas9, CjCas9, ABE8e, PE3, twinPE, ...) |
| **v0.6.0** | 16 systems; initial GraphSAGE + GAT benchmark; ESM-2 initialisation |

The data files (graph, embeddings, DuckDB) correspond to the **v0.7.1 graph rebuild**. The code and API correspond to **v0.7.2**. No graph rebuild was needed for v0.7.2 - only the system naming and API changed.

Full history: [CHANGELOG.md](CHANGELOG.md) | Update policy: [UPDATE_STRATEGY.md](UPDATE_STRATEGY.md)

---

## Part of PEN-STACK

GENOME-ATLAS is the knowledge-graph component of **PEN-STACK**, a set of interconnected tools for programmable genome writing:

| Component | Role |
|-----------|------|
| **GENOME-ATLAS** (this repo) | Knowledge graph + GNN benchmark; editor selection API |
| **PEN-SCORE** | Multi-objective scoring model for editor recommendation |
| **PEN-COMPARE** | Comparative analysis across editing systems |
| **PEN-ASSEMBLE** | Assembly and annotation of new editing systems |
| **mech-class** | Mechanism classification module |

Each component is versioned independently and uses GENOME-ATLAS as its data source via the `load_systems()` and `Atlas` APIs.

---

## Citation

If you use GENOME-ATLAS in your research, please cite the software:

```
Ahmed, A. (2026). GENOME-ATLAS: A unified knowledge graph and GNN benchmark
for programmable genome-writing enzymes. https://github.com/ahmedanees-m/genome-atlas
```

---

## License

MIT - see [LICENSE](LICENSE).

---

## Contact

**Anees Ahmed**
GitHub: [@ahmedanees-m](https://github.com/ahmedanees-m) | Email: ahmedaneesm@gmail.com
