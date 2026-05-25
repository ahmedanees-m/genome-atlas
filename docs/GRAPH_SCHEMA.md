# GENOME-ATLAS Graph Schema

This document describes the heterogeneous knowledge graph schema, the two
available views (`primary` and `full`), and how to select between them.

---

## Node Types

| Node type   | Description                                            | Feature vector        |
|-------------|--------------------------------------------------------|-----------------------|
| `Protein`   | UniProt protein accession                              | ESM-2 640-dim         |
| `Domain`    | Pfam domain family (e.g. PF01548)                      | zeros (topological)   |
| `Structure` | PDB or AlphaFold structure ID                          | zeros (topological)   |
| `Mechanism` | Mechanism class (DSB_NUCLEASE / DSB_FREE / TRANSPOSASE)| zeros (topological)   |
| `Organism`  | NCBI taxonomy ID                                       | zeros (topological)   |
| `System`    | Named foundational system (e.g. IS621_bridge_recombinase) | zeros (topological)|
| `RNA`       | RNA guide/scaffold component (e.g. bridge_RNA, sgRNA) | zeros (topological)   |

Protein nodes with no ESM-2 embedding (e.g. synthetic or not-yet-deposited
proteins) have zero feature vectors. This is the **IS110 OOD scenario** that
prompted the Tier-A gate in `mech-class` v0.5.2.

---

## Edge Types

### Primary View (ML training - v0.7.0 schema)

The **primary view** contains only the four edge types used for GNN benchmark
training. This is the default for `get_graph()` and is backward-compatible with
v0.7.0:

| Edge type (src, rel, dst)                      | Count (v0.7.1) | Description                                   |
|------------------------------------------------|----------------|-----------------------------------------------|
| `Protein -> HAS_DOMAIN -> Domain`                | ~9,533         | Protein has a Pfam domain annotation          |
| `Structure -> STRUCTURE_OF -> Protein`           | ~2,239         | PDB / AlphaFold structure for this protein    |
| `System -> USES_MECHANISM -> Mechanism`          | 28             | System belongs to mechanism class             |
| `System -> HAS_PROTEIN -> Protein`               | 18             | System contains this protein subunit          |

The primary benchmark task is **link prediction on HAS_DOMAIN** (80/10/10 split,
seed=42, n_test=~1,908 pairs including matched negatives).

### Full View (annotation graph - v0.7.1+)

The **full view** adds three secondary edge types:

| Edge type (src, rel, dst)                | Derivation                                                       |
|------------------------------------------|------------------------------------------------------------------|
| `System -> HAS_RNA -> RNA`                 | From `rna_components` in `foundational_systems.yaml`             |
| `Protein -> PART_OF -> System`             | Reverse of `HAS_PROTEIN` (which proteins belong to which system) |
| `Protein -> SIMILAR_TO -> Protein`         | Cosine similarity >= threshold on ESM-2 embeddings                |

These edge types were present in v0.6.0 but were removed from the training
graph in v0.7.0 to simplify the ML benchmark. They are restored in v0.7.1
via the `graph_view='full'` parameter without changing the training graph.

---

## Selecting a View

Use the `get_graph()` function from `genome_atlas.graph`:

```python
from genome_atlas.graph import get_graph

# ML training - compact primary view (default, backward-compatible)
data = get_graph(
    gpickle_path="data/graphs/atlas_train.gpickle",
    esm_emb_path="data/embeddings/esm2_150M_v6.parquet",
    graph_view="primary",   # default
)

# Full annotation graph (all edge types including HAS_RNA, PART_OF, SIMILAR_TO)
data = get_graph(
    gpickle_path="data/graphs/atlas_train.gpickle",
    esm_emb_path="data/embeddings/esm2_150M_v6.parquet",
    graph_view="full",
    similarity_threshold=0.90,   # cosine sim >= 0.90 for SIMILAR_TO edges
)
```

The legacy `build_pyg_hetero()` function is still available and still includes
all edge types present in the gpickle. If you call it directly, use
`_filter_primary(data)` from `genome_atlas.graph.view` to strip secondary edges.

---

## Version History

| Version | Edge types in training graph | Notes |
|---------|------------------------------|-------|
| v0.6.0  | HAS_DOMAIN, STRUCTURE_OF, USES_MECHANISM, HAS_PROTEIN, HAS_RNA, PART_OF, SIMILAR_TO | All 7 types in gpickle |
| v0.7.0  | HAS_DOMAIN, STRUCTURE_OF, USES_MECHANISM, HAS_PROTEIN | HAS_RNA / PART_OF / SIMILAR_TO dropped for ML-benchmark clarity |
| v0.7.1  | HAS_DOMAIN, STRUCTURE_OF, USES_MECHANISM, HAS_PROTEIN | + `graph_view='full'` restores all 7 types at load time; ISCro4 (D2TGM5) added |

---

## Downstream Package Compatibility

| Package        | Pin in v0.7.0              | Pin after v0.7.1           | Requires       |
|----------------|----------------------------|----------------------------|----------------|
| `mech-class`   | `genome-atlas<0.7.0`       | `genome-atlas<0.8.0`       | primary view   |
| `pen-score`    | `genome-atlas<0.7.0`       | `genome-atlas<0.8.0`       | primary view   |
| `pen-assemble` | `genome-atlas<0.7.0`       | `genome-atlas<0.8.0`       | full view      |

After v0.7.1 is released, downstream packages can safely bump their pin from
`<0.7.0` to `<0.8.0`.

---

## Statistics (v0.7.1)

| Metric                  | Value     |
|-------------------------|-----------|
| Foundational systems    | 28 (ISCro4 added; v0.7.2 canonical name) |
| Proteins                | 9,500+    |
| ESM-2 embedding rows    | 10,000    |
| HAS_DOMAIN edges        | ~9,533    |
| HAS_PROTEIN edges       | 18        |
| USES_MECHANISM edges    | 28        |
| STRUCTURE_OF edges      | ~2,239    |
| HAS_RNA edges (full)    | derived from YAML |
| PART_OF edges (full)    | = HAS_PROTEIN count (reversed) |
| SIMILAR_TO edges (full) | varies by threshold |
