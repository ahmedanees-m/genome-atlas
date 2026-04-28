# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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
