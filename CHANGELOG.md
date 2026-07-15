# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.7.2] - 2026-05-25

### Added

- **`genome_atlas.load_systems()`** - new public API function; parses
  `genome_atlas/data/foundational_systems.yaml` and returns a
  `dict[str, SystemEntry]` keyed by canonical system name.
- **`genome_atlas.resolve_system_name(name)`** - alias-aware lookup; resolves
  deprecated names (e.g. "IS622") to canonical (e.g. "ISCro4") with a
  `DeprecationWarning`.
- **`genome_atlas.SystemEntry`** - frozen dataclass exposing per-system metadata
  fields (`name`, `aliases`, `uniprot`, `pfam`, `organism`, `mechanism_bucket`,
  `tier_a_gate`, and others).
- **`tests/unit/test_aliases.py`** - 22 new unit tests covering ISCro4 canonical
  naming, IS622 alias resolution, DeprecationWarning emission, and top-level
  package import surface.
- `foundational_systems.yaml`: new fields on the ISCro4 entry - `aliases`,
  `organism`, `pfam`, `renamed_version` - enabling programmatic metadata access.

### Changed

- **Renamed foundational system `IS622` -> `ISCro4`** in
  `genome_atlas/data/foundational_systems.yaml`.  `IS622` is retained as an
  `aliases` entry for backward compatibility.
- `DATA_PROVENANCE.md`: updated all references to use ISCro4 as the primary
  name; added *Naming convention* section explaining ISCro4 vs IS622 vs IS621.
- `docs/benchmark.rst`, `docs/GRAPH_SCHEMA.md`: updated IS622 -> ISCro4.
- `notebooks/benchmark_results.json`: `_graph_version` string updated.
- `scripts/add_d2tgm5_targets_v2.py`: `protein_name` updated to canonical.

### Rationale

The bridge recombinase from *Citrobacter rodentium* ICC168 (UniProt D2TGM5) has
two names in the literature: **ISCro4** (canonical - UniProt gene name +
Pelea 2026 *Science* DOI:10.1126/science.adz1884) and **IS622** (deprecated -
Perry 2025 *bioRxiv* preprint label, retired upon publication).
This patch adopts the canonical name everywhere and exposes it via the API.

### Compatibility

- Fully backward-compatible: existing call sites that load the YAML and match
  on `"IS622"` will continue to find it in `SystemEntry.aliases`; use
  `resolve_system_name("IS622")` for the canonical name.
- All downstream packages updated in lockstep: mech-class v0.5.4, pen-score
  v0.1.3, pen-assemble v0.5.2.
- No graph rebuild required - data files are unchanged; only the system label
  in the YAML and the new `load_systems()` API are new.

## [0.7.1] - 2026-05-23

### Added

- **IS622 / ISCro4** (`D2TGM5`, 326 aa, *Citrobacter rodentium* ICC168) added to
  `foundational_systems.yaml` as the 28th foundational system.
  - PF01548 + PF02371 dual-domain IS110-family bridge recombinase
  - 20% insertion efficiency in human cells (Perry 2025 bioRxiv; 10.1101/2025.05.14.653916)
  - Highest-profile IS110 human-cell result after IS621; Pelea 2026 *Science* (10.1126/science.adz1884)
  - `tier_a_gate: true` - mech-class Tier-A gate fires for D2TGM5
- **`genome_atlas.graph.view` module** - `get_graph()` function with `graph_view` parameter:
  - `graph_view='primary'` (default): 4-edge ML-training schema; backward-compatible with v0.7.0
  - `graph_view='full'`: primary + three secondary edge types (HAS_RNA, PART_OF, SIMILAR_TO) that
    were present in v0.6.0 but removed from the v0.7.0 training graph
  - `HAS_RNA` derived from `rna_components` in `foundational_systems.yaml`
  - `PART_OF` derived as reverse of `HAS_PROTEIN` (which proteins belong to which system)
  - `SIMILAR_TO` derived from ESM-2 cosine similarity >= configurable threshold (default 0.90)
- **`docs/GRAPH_SCHEMA.md`** - full schema reference: node/edge types, primary vs full view,
  version history, downstream package compatibility table
- **18 new unit tests** in `tests/unit/test_graph_views.py` - all graph_view scenarios,
  edge derivation, self-loop exclusion, and error handling

### Changed

- `pyproject.toml` + `setup.py`: version bumped `0.7.0` -> `0.7.1`
- `docs/conf.py`: release bumped `"0.7.0"` -> `"0.7.1"`
- `genome_atlas/graph/__init__.py`: exports `get_graph` from new `view` module

## [0.7.0] - 2026-05-21

### Added

- **11 new foundational systems** - expands the reference universe coverage from 16 to 27 systems:
  - `SaCas9` (J7RUA5, *S. aureus*; PDB 5CZZ/5CZY) - compact Cas9 for AAV delivery; PAM NNGRRT
  - `CjCas9` (Q0P897, *C. jejuni*; PDB 6GFX) - smallest natural Cas9 (~2.9 kb); PAM NNNNRYAC
  - `NmCas9` (A1IQ68, *N. meningitidis*; PDB 6JDQ) - orthogonal Cas9 for multiplexing; PAM NNNNGATT
  - `evoCas9` - PACE-evolved SpCas9 with expanded PAM compatibility (Thuronyi 2019 *Nat Chem Biol*)
  - `HiFi_Cas9` - high-fidelity SpCas9 R691A variant; >99% off-target reduction (Vakulskas 2018 *Nat Med*)
  - `ABE8e` - adenine base editor; TadA-8e + SpCas9-D10A nickase; A->G editing (Richter 2020 *Nat Biotechnol*)
  - `BE4max` - cytosine base editor; APOBEC1-nSpCas9-UGI; C->T editing (Koblan 2018 *Nat Biotechnol*)
  - `PE3` - third-generation prime editor; adds nicking sgRNA for 3x efficiency boost (Anzalone 2019 *Nature*)
  - `twinPE` - twin prime editing with two pegRNAs; enables large-segment replacement (Anzalone 2021 *Nat Biotechnol*)
  - `phiC31_integrase` - Streptomyces phage serine integrase; clinical landing-pad integration (Groth 2000 *PNAS*)
  - `TP901-1_integrase` - Lactococcus phage serine integrase; orthogonal to phiC31 (Stoll 2002 *Mol Cell*)
- **1 new RNA node** - `nicking_sgRNA` (PE3 non-edited-strand nick guide)
- SaCas9 (J7RUA5), CjCas9 (Q0P897), NmCas9 (A1IQ68) proteins already present in the 10,000-protein graph; `CONTAINS` edges auto-wired on rebuild
- `UPDATE_STRATEGY.md` v0.7.0 rows marked complete

### Changed

- `genome_atlas/data/foundational_systems.yaml`: 16 -> 27 systems, 6 -> 7 RNA types
- `docs/conf.py`: release bumped to `"0.7.0"`
- `docs/benchmark.rst`: benchmark table updated with v0.7.0 AUROC/CI numbers (GraphSAGE 0.9717, GAT 0.9690, Node2Vec 0.9889)
- `README.md`: version badge, graph stats, benchmark table, coverage section, GPU reference updated to v0.7.0
- `DATA_PROVENANCE.md`: 27-system table, coverage 55.1%, new family breakdown
- `notebooks/benchmark_results.json`: all v0.7.0 numbers (_version="0.7.0", _rebuild_date="2026-05-21")

### Fixed

- `scripts/bootstrap_cis_v6.py`: auto-detect node2vec parquet format (v0.7.0 saves 2-row metrics summary; v0.6.0 saved per-node embeddings); prevents `KeyError: 'node_id'` on rebuild

## [0.6.0] - 2026-04-29

### Added

- **RNA nodes** (`nodes_rna` table): 6 guide/scaffold RNA types - `sgRNA`, `crRNA`, `tracrRNA`, `bridge_RNA`, `omegaRNA`, `pegRNA` - added to `foundational_systems.yaml` with length and provenance metadata.
- **`HAS_RNA` edges**: `System -> RNA` edges automatically wired from YAML `rna_components` field (14 edges across 10 RNA-guided systems).
- **Negative control proteins**: `scripts/add_negative_controls.py` samples 500 random unreviewed TrEMBL proteins (no whitelisted Pfam domain) and adds them as a fourth mechanism bucket (`NEGATIVE_CONTROL`) for specificity validation in UMAP and benchmarks.
- **Structural similarity edges** (`SIMILAR_TO`): `scripts/run_foldseek.sh` runs Foldseek all-vs-all on 2,284 PDB structures; `scripts/add_structural_edges.py` ingests TM-score >= 0.5 hits (top-10 per node) as weighted `Structure -> Structure` edges.
- **API methods**: `rna_guides_of_system(name)` returns RNA nodes for a system; `structurally_similar(accession, top_k)` returns structurally homologous proteins via SIMILAR_TO edges and TM-score ranking.
- **`materialize_graph.py`**: gracefully handles missing `nodes_rna` table for backwards compatibility with pre-v0.6.0 databases.

### Changed

- `genome_atlas/data/foundational_systems.yaml`: added `rnas:` section documenting all 6 RNA types.
- `.gitignore`: added `logs/`, `*.log`, and one-off `scripts/vm_*`, `scripts/check_*` patterns.

## [0.5.2] - 2026-04-28

### Fixed

- **API graph queries**: `query_system()`, `proteins_with_domain()`, `domains_of_protein()`, and `structures_of_protein()` now resolve nodes via attribute lookup caches at init time, supporting both named-ID graphs (tests) and the production numeric-ID graph (`System_1`, `Protein_1`, etc.). Previously all graph-based queries silently raised `KeyError` against the real atlas.
- **CI workflow**: Replaced unset shell variable `${pkg}` with literal `genome_atlas` in the import smoke test step.
- **README benchmarks**: Node2Vec AUROC corrected to 0.8411 (single-run test-set value, consistent with GNN reporting); CI range [0.8052-0.8342] retained from bootstrap. Quantum Kernel and Classical RBF rows annotated as node-classification tasks, not link prediction.

## [0.5.1] - 2026-04-27

### Fixed

- **C1 - Cargo-tiered DSB penalty**: Reduced DSB penalty for nucleases on insertions <= 1 kb (HDR-efficient) and SNVs (0.55/0.60 vs flat 0.20).
- **C2 - SNV-specific prime editor boost**: Prime editors now receive `s_cargo = 1.0` for SNV corrections, reflecting their purpose-built design.
- **C3 - Unknown-size AAV guard**: Systems with unlinked proteins (`total_aa == 0`) now score `s_aav = 0.1` instead of falsely receiving ultra-compact AAV scores.

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
