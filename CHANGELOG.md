# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.6.0] - 2026-04-29

### Added

- **RNA nodes** (`nodes_rna` table): 6 guide/scaffold RNA types — `sgRNA`, `crRNA`, `tracrRNA`, `bridge_RNA`, `omegaRNA`, `pegRNA` — added to `foundational_systems.yaml` with length and provenance metadata.
- **`HAS_RNA` edges**: `System → RNA` edges automatically wired from YAML `rna_components` field (14 edges across 10 RNA-guided systems).
- **Negative control proteins**: `scripts/add_negative_controls.py` samples 500 random unreviewed TrEMBL proteins (no whitelisted Pfam domain) and adds them as a fourth mechanism bucket (`NEGATIVE_CONTROL`) for specificity validation in UMAP and benchmarks.
- **Structural similarity edges** (`SIMILAR_TO`): `scripts/run_foldseek.sh` runs Foldseek all-vs-all on 2,284 PDB structures; `scripts/add_structural_edges.py` ingests TM-score ≥ 0.5 hits (top-10 per node) as weighted `Structure → Structure` edges.
- **API methods**: `rna_guides_of_system(name)` returns RNA nodes for a system; `structurally_similar(accession, top_k)` returns structurally homologous proteins via SIMILAR_TO edges and TM-score ranking.
- **`materialize_graph.py`**: gracefully handles missing `nodes_rna` table for backwards compatibility with pre-v0.6.0 databases.

### Changed

- `genome_atlas/data/foundational_systems.yaml`: added `rnas:` section documenting all 6 RNA types.
- `.gitignore`: added `logs/`, `*.log`, and one-off `scripts/vm_*`, `scripts/check_*` patterns.

## [0.5.2] - 2026-04-28

### Fixed

- **API graph queries**: `query_system()`, `proteins_with_domain()`, `domains_of_protein()`, and `structures_of_protein()` now resolve nodes via attribute lookup caches at init time, supporting both named-ID graphs (tests) and the production numeric-ID graph (`System_1`, `Protein_1`, etc.). Previously all graph-based queries silently raised `KeyError` against the real atlas.
- **CI workflow**: Replaced unset shell variable `${pkg}` with literal `genome_atlas` in the import smoke test step.
- **README benchmarks**: Node2Vec AUROC corrected to 0.8411 (single-run test-set value, consistent with GNN reporting); CI range [0.8052–0.8342] retained from bootstrap. Quantum Kernel and Classical RBF rows annotated as node-classification tasks, not link prediction.

## [0.5.1] - 2026-04-27

### Fixed

- **C1 — Cargo-tiered DSB penalty**: Reduced DSB penalty for nucleases on insertions ≤ 1 kb (HDR-efficient) and SNVs (0.55/0.60 vs flat 0.20).
- **C2 — SNV-specific prime editor boost**: Prime editors now receive `s_cargo = 1.0` for SNV corrections, reflecting their purpose-built design.
- **C3 — Unknown-size AAV guard**: Systems with unlinked proteins (`total_aa == 0`) now score `s_aav = 0.1` instead of falsely receiving ultra-compact AAV scores.

### Added

- `SpuFz1_V4` and `enNlovFz2` engineered variants to `foundational_systems.yaml`.
- Validation scenarios updated to expect engineered variants where appropriate.

## [0.5.0] - 2026-04-27

### Added

- Heterogeneous GNN models: GraphSAGE and GAT with residual connections.
- Protein language model embeddings via ESM-2 (640-dim).
- Quantum kernel link prediction baseline (qiskit-machine-learning).
- Selection decision support engine with 70% top-3 validation accuracy.
- CLI entry point (`genome-atlas`) for querying and editor selection.
- Full test suite (unit, integration, regression) with pytest; 64% coverage.
- Sphinx documentation scaffold.
- GitHub Actions CI workflow.

## [0.1.0] - 2026-04-22

### Added

- Initial scaffold and repository setup.
