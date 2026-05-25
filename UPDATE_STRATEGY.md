# Update Strategy for GENOME-ATLAS

Gene editing tools evolve quickly, so ATLAS is built as a living resource rather
than a static snapshot. This document describes the versioning policy, the
criteria for adding new systems, and the pipeline for rebuilding the graph.

---

## Version Policy

| Version series | Scope | Trigger |
|----------------|-------|---------|
| v0.6.x | Bug fixes, documentation, minor script changes | Continuous |
| v0.7.0 | 11 new foundational systems (compact Cas9s, base editors, PE3/twinPE, large serine integrases) | Released 2026-05-21 |
| v0.7.x | Bug fixes and benchmark refresh post-rebuild | Continuous |
| v0.8.0 | Cas13 RNA-targeting family; IscB/Fanzor3; updated scoring for RNA editors | Planned |
| v1.0.0 | Metagenomic mining of the extended catalog; mechanism-class integration | Planned |

Semantic versioning: `MAJOR.MINOR.PATCH`.

- PATCH: no graph changes, no API changes.
- MINOR: new systems added, graph rebuilt, benchmark re-run.
- MAJOR: architecture change (new node types, new scoring model, new API).

---

## Criteria for Adding a New System

A system is included when it satisfies all three of the following:

1. **Mechanistic characterisation:** peer-reviewed paper with biochemical or
   structural evidence for the editing mechanism, not just an activity assay.
2. **Structural evidence:** resolved cryo-EM or X-ray structure deposited in
   PDB, or an AlphaFold model with mean pLDDT >= 70 for the catalytic subunit.
3. **Therapeutic or biotechnological relevance:** demonstrated activity in
   human cells, a mammalian animal model, or active clinical development.

Systems meeting criteria 1 and 3 but not yet criterion 2 are tracked in the
"pending structure" queue below and added once structural data is available.

---

## Systems Tracked for Inclusion

| System | Publication | Status | Version |
|--------|-------------|--------|---------|
| SaCas9 | Ran 2015 *Nat Biotechnol* | Added (J7RUA5, PDB 5CZZ) | v0.7.0 |
| CjCas9 | Kim 2017 *Nat Commun* | Added (Q0P897, PDB 6GFX) | v0.7.0 |
| NmCas9 | Hou 2013 *PNAS* | Added (A1IQ68, PDB 6JDQ) | v0.7.0 |
| evoCas9 (PACE-evolved SpCas9) | Thuronyi 2019 *Nat Chem Biol* | Added (engineered proteins) | v0.7.0 |
| HiFi Cas9 | Vakulskas 2018 *Nat Med* | Added (engineered proteins) | v0.7.0 |
| ABE8e (adenine base editor) | Richter 2020 *Nat Biotechnol* | Added (replaces generic ABE) | v0.7.0 |
| BE4max (cytosine base editor) | Koblan 2018 *Nat Biotechnol* | Added (replaces generic CBE) | v0.7.0 |
| PE3 | Anzalone 2019 *Nature* | Added (nicking_sgRNA added) | v0.7.0 |
| twinPE | Anzalone 2021 *Nat Biotechnol* | Added | v0.7.0 |
| phiC31 integrase | Groth 2000 *PNAS* | Added (PDB 3G1R) | v0.7.0 |
| TP901-1 integrase | Stoll 2002 *Mol Cell* | Added (PDB 3BVP) | v0.7.0 |
| CAST-V-E | Metagenomi 2024 | Structure pending deposition | v0.8.0 |
| IscB | Altae-Tran 2021 *Science* | Literature confirmed; protein ingestion pending | v0.8.0 |
| Cas13d (RNA editor) | Konermann 2018 *Cell* | RNA-targeting scoring logic | v0.8.0 |
| dCas12a-ABE fusion | Xu 2024 *Nat Biotechnol* | Fusion architecture | v0.8.0 |

---

## Rebuild Pipeline

When a qualifying new system is published:

```
1. Add an entry to genome_atlas/data/foundational_systems.yaml
2. Download UniProt TSVs for the Pfam families of the new catalytic domain
   (scripts/download_uniprot_tsv.py)
3. Fetch PDB or AlphaFold structures
   (scripts/download_pdb.py / scripts/download_alphafold.py)
4. Re-run the build:
   python scripts/ingest_to_duckdb_v2.py
   python scripts/materialize_graph.py
   python scripts/create_train_gpickle.py
   python scripts/embed_esm2.py           # GPU recommended
   python scripts/train_gnn.py --model sage
   python scripts/train_gnn.py --model gat
   python scripts/bootstrap_cis_v6.py
5. Update notebooks/benchmark_results.json
6. Bump the MINOR version in pyproject.toml and genome_atlas/_version.py
7. Tag the release and push
8. Refresh the raw data files
```

ESM-2 inference dominates the wall-clock time of a full rebuild; GNN training
takes under a minute on the current graph size.

---

## Community Contributions

Pull requests for new systems are welcome. Each submission should include:

1. A peer-reviewed paper or preprint with a deposited structure (PDB ID or
   AlphaFold accession).
2. The UniProt accession for the catalytic subunit.
3. A `mechanism_bucket` annotation:
   `DSB_NUCLEASE | DSB_FREE_TRANSEST_RECOMBINASE | TRANSPOSASE`.
4. At least one `validation_scenarios.yaml` entry (cell type, edit type,
   delivery vector, published efficiency, and reference DOI).
