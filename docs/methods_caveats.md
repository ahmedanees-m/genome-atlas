# Data Source Caveats (ATLAS v0.5)

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
