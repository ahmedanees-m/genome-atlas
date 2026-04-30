# Data Provenance and Annotation Coverage

This document describes what is in the GENOME-ATLAS knowledge graph, what is
deliberately excluded, and how the coverage compares to the broader literature.

---

## Knowledge Graph Statistics (v0.6.0)

| Metric | Value |
|--------|-------|
| Total nodes | 12,267 |
| Total edges | 13,645 |
| Node types | 6 (Protein, Domain, Structure, Mechanism, RNA, System) |
| Edge types | 7 (HAS_DOMAIN, HAS_PROTEIN, STRUCTURE_OF, SIMILAR_TO, HAS_MECHANISM, HAS_RNA, PART_OF) |
| Foundational systems | 16 |
| Protein sequences with ESM-2 embeddings | 10,000 |
| AlphaFold structures (high-confidence, pLDDT ≥ 70) | 823 |
| UniProt Pfam families included | 18 |

---

## Foundational Systems (16 total)

The 16 systems in `genome_atlas/data/foundational_systems.yaml` were selected
by three criteria applied simultaneously:

1. **Mechanistic characterisation** — a peer-reviewed paper with biochemical
   or structural evidence for the editing mechanism.
2. **Structural evidence** — a resolved cryo-EM/X-ray structure or a
   high-confidence AlphaFold model (pLDDT ≥ 70 for the catalytic subunit).
3. **Therapeutic or biotechnological relevance** — demonstrated activity in
   human cells, animal models, or active clinical development.

| System | Mechanism Class | RNA-guided | Included Since |
|--------|----------------|------------|----------------|
| SpCas9 | DSB nuclease | Yes (sgRNA) | v0.1 |
| Cas12a | DSB nuclease | Yes (crRNA) | v0.1 |
| Cas12f | DSB nuclease | Yes (crRNA) | v0.3 |
| CAST-V-K | DSB-free integrase | Yes (sgRNA) | v0.3 |
| CAST-I-F_evoCAST | DSB-free integrase | Yes (crRNA) | v0.5 |
| IS621_bridge_recombinase | DSB-free recombinase | Yes (bridge RNA) | v0.5 |
| SpuFz1_Fanzor | DSB nuclease (OMEGA) | Yes (omegaRNA) | v0.5 |
| SpuFz1_V4 | DSB nuclease (OMEGA) | Yes (omegaRNA) | v0.5 |
| MmeFz2_Fanzor | DSB nuclease (OMEGA) | Yes (omegaRNA) | v0.5 |
| enNlovFz2 | DSB nuclease (OMEGA) | Yes (omegaRNA) | v0.5 |
| Cre_recombinase | DSB-free recombinase | No | v0.4 |
| Bxb1_integrase | DSB-free integrase | No | v0.4 |
| Tn5_transposase | Transposase | No | v0.4 |
| PE2_prime_editor | DSB nuclease (fusion) | Yes (pegRNA) | v0.5 |
| eePASSIGE | DSB-free recombinase (evolved) | No | v0.6 |
| TnsABC_CAST | DSB-free integrase | Yes (crRNA) | v0.4 |

---

## Coverage Relative to the Literature

Running `scripts/audit_coverage.py` against a curated reference set of 49
well-characterised, therapeutically relevant genome-writing systems (April 2026)
gives the following coverage:

```
Atlas systems:         16
Reference universe:    49 characterised systems
Covered by ATLAS:      19  (38.8%)

Family                                   Coverage
-------------------------------------------------
CRISPR-Cas nucleases                      4/17  (23.5%)
Base and prime editors                     1/8  (12.5%)
CRISPR-associated transposases (CAST)      4/6  (66.7%)
Bridge recombinases                        1/3  (33.3%)
Fanzor / OMEGA                             4/4 (100.0%)
Site-specific recombinases                 2/5  (40.0%)
Transposases                               1/4  (25.0%)
Engineered / evolved systems               2/2 (100.0%)
```

**The 38.8% coverage figure should be interpreted carefully:**

- **All three mechanism classes** (DSB nuclease, DSB-free integrase/recombinase,
  transposase) are represented — 100% class coverage.
- **All five RNA-guided programmable families** (CRISPR-Cas, CAST, bridge RNA,
  Fanzor/OMEGA, prime editing) are represented — 100% family coverage.
- Low coverage in "CRISPR-Cas nucleases" is because the 17 reference entries
  include 8 SpCas9 PAM variants (SaCas9, CjCas9, evoCas9, HiFi Cas9, etc.)
  that differ primarily in PAM sequence, not mechanism. These are planned for
  v0.7 as they become structurally characterised.
- Low coverage in "Base editors" is deliberate: base editors (CBE, ABE, CGBE)
  are fusions of Cas9 + deaminase. The Cas9 component is in ATLAS; the
  deaminase fusion variants are planned for v0.7.
- The full 802k-protein metagenomic catalog (archived) covers thousands of
  candidate systems; ATLAS v0.6.0 exposes only the characterised, validated
  subset.

**Notably absent systems and planned versions:**

| System | Reason absent | Planned |
|--------|--------------|---------|
| SaCas9 / CjCas9 / NmCas9 | Structural characterisation complete; ingestion pending | v0.7 |
| CBE / ABE (base editors) | Fusion architecture requires split-node modelling | v0.7 |
| PE3 / PE4 / twinPE | Extension of PE2 with additional nicking; PE2 present | v0.7 |
| Cas13 variants | RNA-targeting; different scoring logic needed | v0.8 |
| evoCas9 / HiFi Cas9 | PAM-variants of SpCas9; minor additions | v0.7 |
| CAST-V-E / CAST-I-B | Structural data emerging; pending structure confirmation | v0.7 |
| phiC31 / TP901-1 | Site-specific integrases; well-characterised; pending ingestion | v0.7 |

---

## Protein Sequence Coverage

The 10,000 Protein nodes in the graph represent all UniProt accessions from
the 18 included Pfam families (see `data/raw/uniprot/PF*.tsv.gz`). Each
protein is embedded with ESM-2 (esm2_t30_150M_UR50D, 640-dim). Proteins
without a UniProt entry (e.g., synthetic or engineered systems not yet
deposited) are represented as zero-feature nodes.

The 18 Pfam families were chosen to cover the catalytic domains of all 16
foundational systems. Coverage is intentionally narrow — including all Pfam
families would add millions of non-relevant proteins and dilute the graph.

---

## Structural Coverage

- **823 AlphaFold structures** for target enzymes (pLDDT ≥ 70 for at least
  one catalytic domain chain).
- **PDB structures** included where available (via `data/raw/pdb/`).
- Structural similarity edges (SIMILAR_TO) computed by Foldseek with
  TM-score ≥ 0.5.

---

## Raw Data Sources

| Source | Files | Version / Date |
|--------|-------|----------------|
| UniProt | `data/raw/uniprot/PF*.tsv.gz` (18 files) | April 2026 release |
| AlphaFold DB | `data/raw/alphafold_targets/*.pdb.gz` (823 files) | v4 models |
| RCSB PDB | `data/raw/pdb/*.pdb.gz` | April 2026 snapshot |
| CRISPRCasDB | `data/processed/crisprcasdb_*.parquet` | April 2026 |
| ISfinder | `data/processed/isfinder.parquet` | April 2026 |

All source data is deposited in the Zenodo archive (DOI forthcoming) alongside
the processed graphs and embeddings.
