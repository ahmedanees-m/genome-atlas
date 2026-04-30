# Update Strategy for GENOME-ATLAS

Gene editing tools evolve rapidly. ATLAS is designed as a **living resource**,
not a static snapshot. This document describes the versioning policy, the
criteria for adding new systems, and the automated pipeline that keeps the
companion PEN-COMPARE tool current.

---

## Version Policy

| Version series | Scope | Trigger |
|----------------|-------|---------|
| v0.6.x | Bug fixes, documentation, minor script changes | Continuous (GitHub push) |
| v0.7.0 | New foundational systems with published structures; base-editor split-node model | Q3 2026 |
| v0.8.0 | Cas13 RNA-targeting family; updated scoring for RNA editors | Q4 2026 |
| v1.0.0 | Full metagenomic mining of the 802k extended catalog; MECH-CLASS integration | After Paper 2 acceptance |

Semantic versioning: `MAJOR.MINOR.PATCH`.  
- PATCH: no graph changes, no API changes.  
- MINOR: new systems added, graph rebuilt, benchmark re-run.  
- MAJOR: architecture change (new node types, new scoring model, new API).

---

## Criteria for Adding a New System

A system is included when it satisfies **all three** of the following:

1. **Mechanistic characterisation** — peer-reviewed paper with biochemical or
   structural evidence for the editing mechanism (not just activity assay).
2. **Structural evidence** — resolved cryo-EM/X-ray structure deposited in
   PDB, or AlphaFold model with mean pLDDT ≥ 70 for the catalytic subunit.
3. **Therapeutic or biotechnological relevance** — demonstrated activity in
   human cells, a mammalian animal model, or active clinical development (IND
   filed or equivalent).

Systems meeting criteria 1 and 3 but not yet criterion 2 are tracked in the
"pending structure" queue below and added once structural data is available.

---

## Systems Tracked for Inclusion

| System | Publication | Blocking criterion | Planned version |
|--------|-------------|-------------------|-----------------|
| SaCas9 | Ran 2015 *Nat Biotechnol* | Ingestion pending (structure available) | v0.7.0 |
| CjCas9 | Kim 2017 *Nat Commun* | Ingestion pending | v0.7.0 |
| evoCas9 (PACE-evolved SpCas9) | Thuronyi 2019 *Nat Chem Biol* | Ingestion pending | v0.7.0 |
| HiFi Cas9 | Vakulskas 2018 *Nat Med* | Ingestion pending | v0.7.0 |
| CBE / ABE (base editors) | Komor 2016 / Gaudelli 2017 | Split-node architecture for fusion proteins | v0.7.0 |
| PE3 / PE4 / twinPE | Anzalone 2021 *Nat Biotechnol* | Extension of PE2 (already in graph) | v0.7.0 |
| phiC31 / TP901-1 | Groth 2000 / Stoll 2002 | Ingestion pending | v0.7.0 |
| CAST-V-E | Metagenomi 2024 | Structure pending deposition | v0.7.0 or v1.0 |
| Cas13d (RNA editor) | Konermann 2018 *Cell* | RNA-targeting scoring logic | v0.8.0 |
| dCas12a-ABE fusion | Xu 2024 *Nat Biotechnol* | Fusion architecture | v0.8.0 |
| Engineered Fanzor3 | Zhao 2026 (expected) | Awaiting publication | v0.8.0 |

---

## Update Pipeline

### ATLAS Core (graph + GNN, manual)

When a qualifying new system is published:

```
1. Add entry to genome_atlas/data/foundational_systems.yaml
2. Download UniProt TSVs for all Pfam families of the new catalytic domain
   (scripts/download_uniprot_tsv.py)
3. Fetch PDB or AlphaFold structures
   (scripts/download_pdb.py / scripts/download_alphafold.py)
4. Re-run build pipeline:
   python scripts/materialize_graph.py
   python scripts/create_train_gpickle.py
   python scripts/embed_esm2.py          # GPU required
   python scripts/train_gnn.py --model sage
   python scripts/train_gnn.py --model gat
   python scripts/bootstrap_cis_v6.py
5. Update notebooks/benchmark_results.json
6. Bump MINOR version in pyproject.toml / genome_atlas/_version.py
7. Tag release and push to GitHub
8. Update Zenodo deposit (DOI forthcoming)
```

Estimated wall-clock time for a full rebuild: ~75 minutes on V100 GPU
(dominated by ESM-2 inference; GNN training <1 minute on current graph size).

### PEN-COMPARE Companion Tool (automated weekly)

The companion web tool PEN-COMPARE (Deliverable 5 of the PEN-STACK roadmap)
runs a **weekly automated refresh**:

1. **Discovery** — Query bioRxiv API and PubMed RSS for new manuscripts
   mentioning: `"CRISPR" OR "CAST" OR "bridge recombinase" OR "Fanzor"
   OR "prime editor" OR "base editor"`.
2. **Triage** — Flag new candidate systems for manual review (~5–10 min/week).
3. **Score** — Run PEN-SCORE on confirmed entries during the next HPC
   short-queue slot (no GPU required for scoring alone).
4. **Deploy** — Auto-commit updated scores to the `pen-compare` Parquet
   database; Vercel rebuilds the frontend automatically.

---

## Community Contributions

Pull requests for new systems are welcome. Requirements:

1. Peer-reviewed paper or preprint with deposited structure (PDB ID or
   AlphaFold accession).
2. UniProt accession for the catalytic subunit.
3. `mechanism_bucket` annotation:
   `DSB_NUCLEASE | DSB_FREE_TRANSEST_RECOMBINASE | TRANSPOSASE`.
4. At least one `validation_scenarios.yaml` entry (cell type, edit type,
   delivery vector, published efficiency and reference DOI).

See `CONTRIBUTING.md` for the submission template and CI requirements.

---

## Long-Term Vision: Metagenomic Discovery Engine

Year-2 stretch goal: systematic ATLAS-driven mining of the 802k extended
protein catalog (currently archived at `data/raw/` on the training VM) for
overlooked RNA-guided systems in bacterial and archaeal genomes. Target:
identify 1–2 new families of programmable recombinases and generate
AlphaFold structure predictions, with wet-lab validation in collaboration
with an experimental partner lab.

This work is tracked as Paper 2 (MECH-CLASS) in the PEN-STACK roadmap.
