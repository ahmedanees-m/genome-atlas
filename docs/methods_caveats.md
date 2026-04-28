# Data Source Caveats (ATLAS v0.5.2)

## CRISPRCasdb
No public bulk JSON dump was available via the documented API at the time of construction (April 2026).
The 13 systems and 19 proteins encoded in `foundational_systems.yaml` were curated manually from the primary literature and cross-referenced against UniProt/CRISPRCasdb web records.
This is sufficient for the v0.5 atlas scope; full CRISPRCasdb bulk ingestion is planned for v1.0.

## ISfinder
The canonical ISfinder FTP (Toulouse) and NCBI mirror were unreachable during data acquisition.
We obtained 8,433 IS elements from a community GitHub mirror and parsed them with `scripts/parse_isfinder.py`.
These are used for family-level annotation; individual protein accessions were mapped to UniProt via sequence identity.

## AlphaFold DB
Species-tar downloads covered SwissProt + 9 priority proteomes.
Structural coverage for the curated 9,500-protein atlas is ~15% (PDB + AFDB); the remaining 85% are sequence-only nodes.
The full 802,559-protein extended catalog is retained for future releases.

## Downsampled Atlas Scope
The v0.5 atlas uses a curated 9,500-protein subset (targets_v2.parquet) selected from the full 802,559-protein extended catalog (targets_v1.parquet).
Selection criteria: all 863 SwissProt-reviewed entries preserved, plus top-scoring TrEMBL entries by structural evidence, domain diversity, length, and priority organism membership.
This scope was chosen for the NAR Database Issue reviewer envelope (5k-10k proteins).

## Benchmark Reporting Notes

### Confidence Intervals
GNN (GraphSAGE, GAT) confidence intervals are computed via the Mann-Whitney SE formula applied to the test-set AUC, not bootstrap resampling. The `bootstrap_cis.parquet` file in `data/embeddings/` reflects bootstrap resampling of the **training** graph embeddings for a separate internal diagnostic and reports inflated AUROCs (~0.999); it should not be cited for primary benchmark confidence intervals.

### `USES_MECHANISM` Edge Type (GNN Benchmark)
The `System_USES_MECHANISM_Mechanism` edge type has only 14 edges in the graph (one per foundational system). With a 20% test split, the expected number of positive test edges is ~2.8, making statistical AUC estimation unreliable. GraphSAGE reports AUROC = 0.000 (inverted random) and GAT reports 0.500 (chance) for this edge type — both are artifacts of the trivial test set size. These results are excluded from the primary Protein→Domain benchmark and are reported in Supplementary Table S3 with this note.

### Graph Node ID Convention
The production atlas graph (`atlas.gpickle`) uses numeric node IDs (`System_1`, `Protein_1`, etc.) with named attributes. The `Atlas` API resolves queries by attribute lookup (system `name`, domain `accession`, protein `accession`) and is compatible with both the production graph and test fixtures using named-ID convention.
