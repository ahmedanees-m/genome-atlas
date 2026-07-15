# Data Provenance and Annotation Coverage

This document describes what is in the GENOME-ATLAS knowledge graph, what is
deliberately excluded, and how the coverage compares to the broader literature.

---

## Knowledge Graph Statistics (v0.7.2)

| Metric | Value |
|--------|-------|
| Total nodes (full graph) | ~13,400+ (rebuild pending) |
| Total nodes (training graph) | ~11,762+ (rebuild pending) |
| Total edges | ~11,820+ (rebuild pending) |
| Node types | 7 (Protein, Domain, Structure, RNA, Mechanism, Organism, System) |
| Edge types | 4 primary (HAS_DOMAIN, HAS_PROTEIN, USES_MECHANISM, STRUCTURE_OF) + 3 secondary via `graph_view='full'` |
| Foundational systems | 28 (ISCro4 added in v0.7.1; canonical name adopted in v0.7.2) |
| Proteins | 9,501+ (D2TGM5 added) |
| ESM-2 embedding file rows | 10,001 (D2TGM5 embedding added 2026-05-23) |
| Protein nodes in graph with embeddings | 9,501 (subset of ESM-2 file) |
| AlphaFold structures (high-confidence, pLDDT >= 70) | 823 (unchanged) |
| UniProt Pfam families included | 18 (unchanged - ISCro4 uses PF01548+PF02371, both already in whitelist) |
| HAS_DOMAIN edges | ~9,535+ (2 new edges for D2TGM5->PF01548, D2TGM5->PF02371) |
| HAS_PROTEIN edges | 18 (1 new: ISCro4->D2TGM5) |
| USES_MECHANISM edges | 28 (1 new: ISCro4->DSB_FREE_TRANSEST_RECOMBINASE) |
| STRUCTURE_OF edges | 2,239 (unchanged) |
| PDB structures | 2,284 (unchanged) |
| Organisms | 745+ (637910 Citrobacter rodentium ICC168 added if not present) |

*v0.7.1 graph rebuilt 2026-05-23. GraphSAGE AUROC=0.9714 [0.9625, 0.9797]; GAT AUROC=0.9690 [0.9590, 0.9778]. Bootstrap CIs in `reproduction/bootstrap_ci_v7.json`.*
*v0.7.0 graph rebuilt from scratch 2026-05-21.*

---

## Foundational Systems (28 total)

The 28 systems in `genome_atlas/data/foundational_systems.yaml` were selected
by three criteria applied simultaneously:

1. **Mechanistic characterisation** - a peer-reviewed paper with biochemical
   or structural evidence for the editing mechanism.
2. **Structural evidence** - a resolved cryo-EM/X-ray structure or a
   high-confidence AlphaFold model (pLDDT >= 70 for the catalytic subunit).
3. **Therapeutic or biotechnological relevance** - demonstrated activity in
   human cells, animal models, or active clinical development.

### Original 16 systems (v0.1 - v0.6.0)

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

### New 11 systems added in v0.7.0

| System | Description | PAM / Key Feature | Mechanism Class |
|--------|-------------|-------------------|----------------|
| SaCas9 | *Staphylococcus aureus* Cas9 | NNGRRT; compact for AAV | DSB nuclease |
| CjCas9 | *Campylobacter jejuni* Cas9 | NNNNRYAC; smallest natural Cas9 | DSB nuclease |
| NmCas9 | *Neisseria meningitidis* Cas9 | NNNNGATT; orthogonal for multiplexing | DSB nuclease |
| evoCas9 | PACE-evolved SpCas9 | Expanded PAM compatibility | DSB nuclease |
| HiFi_Cas9 | SpCas9 R691A variant | >99% off-target reduction | DSB nuclease |
| ABE8e | Adenine base editor | TadA-8e + SpCas9 nickase; A->G | DSB-free (base edit) |
| BE4max | Cytosine base editor | APOBEC1-nSpCas9-UGI; C->T | DSB-free (base edit) |
| PE3 | Prime editor 3 | PE2 + nicking sgRNA; ~3x efficiency vs PE2 | DSB-free (prime edit) |
| twinPE | Twin prime editing | Two pegRNAs; large-segment replacement | DSB-free (prime edit) |
| phiC31 | Streptomyces phage integrase | Serine integrase; clinical landing-pad | DSB-free integrase |
| TP901-1 | Lactococcus phage integrase | Serine integrase; orthogonal to phiC31 | DSB-free integrase |

New RNA node added in v0.7.0: **nicking_sgRNA** (PE3 non-edited-strand nicking guide).

### New 1 system added in v0.7.1 (canonical name adopted in v0.7.2)

| System | Description | PAM / Key Feature | Mechanism Class |
|--------|-------------|-------------------|----------------|
| **ISCro4** | *Citrobacter rodentium* ICC168 IS110 bridge recombinase | PF01548+PF02371 dual-domain; >6% insertion efficiency in human cells | DSB-free recombinase |

D2TGM5 (326 aa) - primary reference: Pelea 2026 *Science* 10.1126/science.adz1884; secondary: Perry 2025 bioRxiv 10.1101/2025.05.14.653916. `tier_a_gate: true`. The deprecated preprint label "IS622" is retained as an alias in `foundational_systems.yaml` for backward compatibility.

---

## Coverage Relative to the Literature

Running `scripts/audit_coverage.py` against a curated reference set of 49
well-characterised, therapeutically relevant genome-writing systems (April 2026)
gives the following coverage for v0.7.1:

```
Atlas systems:         28  (was 27 in v0.7.0)
Reference universe:    49 characterised systems
Covered by ATLAS:      28  (~57.1%)

Family                                   Coverage
-------------------------------------------------
CRISPR-Cas nucleases                      9/17  (52.9%)
Base and prime editors                     4/8  (50.0%)
CRISPR-associated transposases (CAST)      4/6  (66.7%)
Bridge recombinases                        2/3  (66.7%)  <- +1 from ISCro4 (D2TGM5)
Fanzor / OMEGA                             4/4 (100.0%)
Site-specific recombinases                 4/5  (80.0%)
Transposases                               1/4  (25.0%)
Engineered / evolved systems               2/2 (100.0%)
```

**The ~57.1% coverage figure should be interpreted carefully:**

- **All three mechanism classes** (DSB nuclease, DSB-free integrase/recombinase,
  transposase) are represented - 100% class coverage.
- **All five RNA-guided programmable families** (CRISPR-Cas, CAST, bridge RNA,
  Fanzor/OMEGA, prime editing) are represented - 100% family coverage.
- Bridge recombinase coverage increased from 33.3% (v0.7.0) to **66.7%** (v0.7.1)
  with the addition of ISCro4 (D2TGM5; *Citrobacter rodentium* ICC168).
- CRISPR-Cas nuclease coverage increased from 23.5% (v0.6.0) to 52.9% (v0.7.0)
  with the addition of SaCas9, CjCas9, NmCas9, evoCas9, and HiFi_Cas9.
- Base/prime editor coverage increased from 12.5% to 50.0% with ABE8e, BE4max,
  PE3, and twinPE.
- Site-specific recombinase coverage increased from 40.0% to 80.0% with phiC31
  and TP901-1.
- The full 802k-protein metagenomic catalog (archived) covers thousands of
  candidate systems; ATLAS v0.7.1 exposes only the characterised, validated
  subset.

**Remaining absent systems and planned versions:**

| System | Reason absent | Planned |
|--------|--------------|---------|
| CAST-V-E / CAST-I-B | Structural data emerging; pending structure confirmation | v0.8 |
| Cas13 variants | RNA-targeting; different scoring logic needed | v0.8 |
| Additional Tn7-based CAST | Partial structural data | v0.8 |
| SpRY / xCas9 | Expanded-PAM SpCas9 variants; minor additions | v0.8 |
| Additional site-specific recombinases | Well-characterised but pending ingestion | v1.0 |

---

## Protein Sequence Coverage

The 9,500 Protein nodes in the v0.7.0 graph represent UniProt accessions from
the 18 included Pfam families after deduplication. Each protein is embedded
with ESM-2 (esm2_t30_150M_UR50D, 640-dim). Proteins without a UniProt entry
(e.g., synthetic or engineered systems not yet deposited) are represented as
zero-feature nodes.

The ESM-2 embedding file (`esm2_150M_v6.parquet`) contains 10,000 protein
embeddings (the full v0.6.0 universe). Of these, 9,500 appear as Protein nodes
in the v0.7.0 graph; the remaining 500 are in the file but not the graph.

The 18 Pfam families were chosen to cover the catalytic domains of all 28
foundational systems. The new v0.7.0 systems (Cas9 variants, base editors,
serine integrases) share Pfam families already in the whitelist - no new
families were required and no new ESM-2 inference was needed.

---

## Structural Coverage

- **823 AlphaFold structures** for target enzymes (pLDDT >= 70 for at least
  one catalytic domain chain). Unchanged from v0.6.0.
- **2,284 PDB structures** included where available.
- Structural data unchanged from v0.6.0; new v0.7.0 systems reuse existing
  structural entries where PDB/AlphaFold coverage overlaps with SpCas9 family.

---

## Raw Data Sources

| Source | Files | Version / Date |
|--------|-------|----------------|
| UniProt | `data/raw/uniprot/PF*.tsv.gz` (18 files) | April 2026 release (unchanged from v0.6.0) |
| AlphaFold DB | `data/raw/alphafold_targets/*.pdb.gz` (823 files) | v4 models (unchanged from v0.6.0) |
| RCSB PDB | `data/raw/pdb/*.pdb.gz` | April 2026 snapshot (unchanged) |
| CRISPRCasDB | `data/processed/crisprcasdb_*.parquet` | April 2026 (unchanged) |
| New systems | 11 systems from peer-reviewed literature | 2000-2021 publications |

The raw source data, processed graphs, and embeddings are available as raw data files.

---

## Naming convention: ISCro4 (canonical) vs IS622 (deprecated)

The bridge recombinase from *Citrobacter rodentium* ICC168 (UniProt D2TGM5, 326 aa,
PF01548 + PF02371 domain architecture) is referred to by two names in the literature:

| Name | Source | Status |
|------|--------|--------|
| **ISCro4** | UniProt D2TGM5 gene name; Pelea et al. 2026 *Science* DOI:10.1126/science.adz1884 | **Canonical - use this** |
| IS622 | Perry et al. 2025 *bioRxiv* preprint (Hsu lab) DOI:10.1101/2025.05.14.653916 | Deprecated preprint label |

As of v0.7.2, GENOME-ATLAS uses **ISCro4** as the canonical system identifier
in `foundational_systems.yaml`, all API responses, and all downstream packages.
The string "IS622" is preserved in the `aliases` field of the ISCro4 entry for
backward compatibility and to ensure searches against either name resolve correctly
via `genome_atlas.resolve_system_name("IS622")`.

Both names refer to the same gene and protein - there is no biological difference.

**Key distinctions** (frequently confused):

| Name | Organism | UniProt | Relationship |
|------|----------|---------|--------------|
| **ISCro4** | *Citrobacter rodentium* ICC168 | D2TGM5 | Same protein; canonical published name |
| IS622 | *Citrobacter rodentium* ICC168 | D2TGM5 | Same protein; deprecated preprint label |
| IS621 | *Escherichia coli* | A0A1H8KQQ4 | **Different protein** - do not conflate |

IS621 and ISCro4 are different proteins from different organisms.
IS622 and ISCro4 are identical proteins - only the label changed.
