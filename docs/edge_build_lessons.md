# Edge Build Engineering Notes

## Lesson 1: Drop indexes before bulk insert
DuckDB `INSERT INTO ... SELECT` on a table with 3 indexes + PRIMARY KEY caused >5 min hangs at 100% CPU for 2M rows.
Fix: `DROP INDEX` -> bulk insert -> `CREATE INDEX`. Runtime dropped from ~300s to ~1s.

## Lesson 2: Column-name drift between Parquet and DuckDB
The Parquet file used `xref_pdb` / `xref_pfam` / `primary_mechanism_bucket`, but DuckDB `nodes_protein` did not have these columns.
Edge builders that joined against `nodes_protein` silently failed or produced empty sets.
Fix: Join against the original Parquet file for edge construction, not the normalized DuckDB node table, or add the columns to the node schema.

## Lesson 3: Validate edge cardinality by source_type
A bug in v8 assigned `edge_type='HAS_PROTEIN'` to `Organism->Protein` edges, producing 802k false System-protein links.
Fix: Always `GROUP BY source_type, edge_type` after ingestion and assert expected cardinalities (e.g., `HAS_PROTEIN` may only originate from `System`).

## Lesson 4: Correct edge semantics matter more than speed
The original edge builder created `Protein->Mechanism` edges (410k) instead of `System->Mechanism` edges (~14).
This broke the selection API which traverses `System -> USES_MECHANISM -> Mechanism`.
Fix: Map `mechanism_bucket` from `foundational_systems.yaml` to `nodes_mechanism`, creating System->Mechanism edges.
