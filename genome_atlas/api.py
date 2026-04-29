"""Public API for the GENOME-ATLAS."""
from __future__ import annotations

import pickle
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd


DEFAULT_GPICKLE = Path.home() / "pen-stack/data/graphs/atlas.gpickle"
DEFAULT_EMBEDDINGS = Path.home() / "pen-stack/data/embeddings/graphsage.parquet"
DEFAULT_TARGETS = Path.home() / "pen-stack/data/processed/targets_v2.parquet"


@dataclass
class EditorRecommendation:
    system: str
    pen_score: float
    aav_fit: bool
    mechanism: str
    reasoning: list[str]


class Atlas:
    """Main entry point for querying the ATLAS.

    Usage::
        >>> atlas = Atlas.load()
        >>> cas9 = atlas.query_system("SpCas9")
        >>> recs = atlas.select_editor(cell_type="HEK293T",
        ...                            edit_type="insertion",
        ...                            cargo_size_bp=1500,
        ...                            delivery="AAV")
    """

    def __init__(self, gpickle_path: Optional[Path] = None,
                 embeddings_path: Optional[Path] = None,
                 targets_path: Optional[Path] = None):
        self._G = None
        self._embeddings = None
        self._length_map: dict[str, int] = {}
        self._protein_map: dict[str, dict] = {}

        # Lookup caches for attribute-based queries (handles both named-ID and
        # numeric-ID graph conventions, e.g. System_SpCas9 vs System_1)
        self._system_by_name: dict[str, str] = {}
        self._domain_by_accession: dict[str, str] = {}
        self._protein_by_accession: dict[str, str] = {}
        self._rna_by_name: dict[str, str] = {}        # v0.6.0

        if targets_path and targets_path.exists():
            targets = pd.read_parquet(targets_path, columns=["accession", "length", "protein_name", "organism_name"])
            self._length_map = {
                f"Protein_{r['accession']}": int(r["length"])
                for _, r in targets.iterrows()
            }
            self._protein_map = {
                f"Protein_{r['accession']}": {
                    "accession": r["accession"],
                    "length": int(r["length"]),
                    "protein_name": r.get("protein_name", ""),
                    "organism_name": r.get("organism_name", ""),
                }
                for _, r in targets.iterrows()
            }

        if gpickle_path and gpickle_path.exists():
            with open(gpickle_path, "rb") as f:
                self._G = pickle.load(f)

        if self._G is not None:
            for n, d in self._G.nodes(data=True):
                nt = d.get("node_type", "")
                if nt == "System" and "name" in d:
                    self._system_by_name[d["name"]] = n
                elif nt == "Domain" and "accession" in d:
                    self._domain_by_accession[d["accession"]] = n
                elif nt == "Protein" and "accession" in d:
                    self._protein_by_accession[d["accession"]] = n
                elif nt == "RNA" and "name" in d:
                    self._rna_by_name[d["name"]] = n

        if embeddings_path and embeddings_path.exists():
            self._embeddings = pd.read_parquet(embeddings_path).set_index("node_id")

    @classmethod
    def load(cls, gpickle_path: Path = DEFAULT_GPICKLE,
             embeddings_path: Path = DEFAULT_EMBEDDINGS,
             targets_path: Path = DEFAULT_TARGETS) -> "Atlas":
        return cls(gpickle_path, embeddings_path, targets_path)

    # ------- Basic queries -------

    def _resolve_node(self, prefix: str, key: str,
                      cache: dict[str, str]) -> Optional[str]:
        """Resolve a node ID by trying the named-ID convention then the cache."""
        direct = f"{prefix}_{key}"
        if self._G is not None and direct in self._G:
            return direct
        return cache.get(key)

    def query_system(self, name: str) -> dict:
        """Return system metadata by name."""
        node_id = self._resolve_node("System", name, self._system_by_name)
        if self._G is not None and node_id is not None and node_id in self._G:
            return {"node_id": node_id, **self._G.nodes[node_id]}
        raise KeyError(f"System '{name}' not found")

    def query_protein(self, accession: str) -> dict:
        """Return protein metadata by UniProt accession."""
        node_id = self._resolve_node("Protein", accession, self._protein_by_accession)
        map_key = f"Protein_{accession}"
        if node_id is not None and self._G is not None and node_id in self._G:
            meta = self._protein_map.get(map_key, {})
            return {"node_id": node_id, **self._G.nodes[node_id], **meta}
        if map_key in self._protein_map:
            return {"node_id": map_key, **self._protein_map[map_key]}
        raise KeyError(f"Protein '{accession}' not found")

    def systems(self, mechanism_bucket: Optional[str] = None) -> pd.DataFrame:
        """Return all systems, optionally filtered by mechanism bucket."""
        if self._G is None:
            raise RuntimeError("Graph not loaded")
        rows = []
        for n, d in self._G.nodes(data=True):
            if d.get("node_type") == "System":
                if mechanism_bucket is None or d.get("mechanism_bucket") == mechanism_bucket:
                    rows.append({"node_id": n, **d})
        return pd.DataFrame(rows)

    def proteins_with_domain(self, domain_accession: str) -> pd.DataFrame:
        """Return proteins that have a given Pfam domain."""
        if self._G is None:
            raise RuntimeError("Graph not loaded")
        domain_node = self._resolve_node("Domain", domain_accession, self._domain_by_accession)
        if domain_node is None or domain_node not in self._G:
            raise KeyError(f"Domain '{domain_accession}' not found")
        proteins = []
        for u, v, edge_data in self._G.edges(data=True):
            if edge_data.get("edge_type") == "HAS_DOMAIN" and v == domain_node:
                acc = self._G.nodes[u].get("accession", "")
                info = self._protein_map.get(f"Protein_{acc}", {"node_id": u, **self._G.nodes[u]})
                proteins.append(info)
        return pd.DataFrame(proteins)

    def domains_of_protein(self, accession: str) -> pd.DataFrame:
        """Return domains of a given protein."""
        if self._G is None:
            raise RuntimeError("Graph not loaded")
        protein_node = self._resolve_node("Protein", accession, self._protein_by_accession)
        if protein_node is None or protein_node not in self._G:
            raise KeyError(f"Protein '{accession}' not found")
        domains = []
        for u, v, edge_data in self._G.edges(data=True):
            if edge_data.get("edge_type") == "HAS_DOMAIN" and u == protein_node:
                domains.append({"node_id": v, **self._G.nodes[v]})
        return pd.DataFrame(domains)

    def rna_guides_of_system(self, name: str) -> pd.DataFrame:
        """Return RNA guide/scaffold nodes for a given system.

        Example::
            >>> atlas.rna_guides_of_system("SpCas9")
               node_id      name    rna_type  length_nt
            0  RNA_1    sgRNA   guide_RNA       100
        """
        if self._G is None:
            raise RuntimeError("Graph not loaded")
        sys_node = self._resolve_node("System", name, self._system_by_name)
        if sys_node is None or sys_node not in self._G:
            raise KeyError(f"System '{name}' not found")
        rnas = []
        for u, v, edge_data in self._G.edges(data=True):
            if edge_data.get("edge_type") == "HAS_RNA" and u == sys_node:
                rnas.append({"node_id": v, **self._G.nodes[v]})
        return pd.DataFrame(rnas)

    def structurally_similar(self, accession: str,
                             top_k: int = 5) -> pd.DataFrame:
        """Return structures similar to those of a given protein (via SIMILAR_TO edges).

        Returns a DataFrame of Structure nodes sorted by TM-score descending.

        Example::
            >>> atlas.structurally_similar("Q99ZW2", top_k=3)
        """
        if self._G is None:
            raise RuntimeError("Graph not loaded")
        protein_node = self._resolve_node("Protein", accession, self._protein_by_accession)
        if protein_node is None or protein_node not in self._G:
            raise KeyError(f"Protein '{accession}' not found")

        # Collect Structure nodes for this protein
        struct_nodes = [
            v for u, v, d in self._G.edges(data=True)
            if d.get("edge_type") == "STRUCTURE_OF" and u == protein_node
        ]

        hits = []
        for s_node in struct_nodes:
            for u, v, d in self._G.edges(data=True):
                if d.get("edge_type") == "SIMILAR_TO" and u == s_node:
                    hits.append({
                        "node_id": v,
                        "tmscore": d.get("weight", 0.0),
                        **self._G.nodes[v],
                    })

        if not hits:
            return pd.DataFrame()
        df = pd.DataFrame(hits).drop_duplicates("node_id")
        return df.sort_values("tmscore", ascending=False).head(top_k).reset_index(drop=True)

    def structures_of_protein(self, accession: str) -> pd.DataFrame:
        """Return structures (PDB / AlphaFold) of a given protein."""
        if self._G is None:
            raise RuntimeError("Graph not loaded")
        protein_node = self._resolve_node("Protein", accession, self._protein_by_accession)
        if protein_node is None or protein_node not in self._G:
            raise KeyError(f"Protein '{accession}' not found")
        structures = []
        for u, v, edge_data in self._G.edges(data=True):
            if edge_data.get("edge_type") == "STRUCTURE_OF" and v == protein_node:
                structures.append({"node_id": u, **self._G.nodes[u]})
        return pd.DataFrame(structures)

    # ------- Embedding-based similarity -------

    def get_embedding(self, node_id: str) -> np.ndarray:
        """Return 128-dim embedding for a node."""
        if self._embeddings is None:
            raise RuntimeError("No embeddings loaded.")
        if node_id not in self._embeddings.index:
            raise KeyError(f"No embedding for {node_id}")
        return np.array(self._embeddings.loc[node_id, "embedding"])

    def similar_nodes(self, node_id: str, node_type: str = "System",
                      top_k: int = 5) -> pd.DataFrame:
        """Top-k most similar nodes by cosine similarity in embedding space."""
        if self._embeddings is None:
            raise RuntimeError("No embeddings loaded.")
        query_emb = self.get_embedding(node_id)
        subset = self._embeddings[self._embeddings["node_type"] == node_type]
        if subset.empty:
            return pd.DataFrame()
        embs = np.stack([np.array(e) for e in subset["embedding"]])
        query_norm = query_emb / (np.linalg.norm(query_emb) + 1e-10)
        embs_norm = embs / (np.linalg.norm(embs, axis=1, keepdims=True) + 1e-10)
        sims = embs_norm @ query_norm
        top_idx = np.argsort(-sims)[:top_k + 1]
        results = subset.iloc[top_idx].copy()
        results["similarity"] = sims[top_idx]
        results = results[results.index != node_id].head(top_k)
        return results

    # ------- Selection Decision Support -------

    def select_editor(self, cell_type: str = "HEK293T",
                      edit_type: str = "insertion",
                      cargo_size_bp: int = 0,
                      delivery: str = "AAV",
                      prefer_dsb_free: bool = True,
                      top_k: int = 5) -> list[EditorRecommendation]:
        """Rank candidate editors for a use case."""
        from genome_atlas.selection import SelectionEngine
        engine = SelectionEngine(self)
        return engine.rank(
            cell_type=cell_type, edit_type=edit_type,
            cargo_size_bp=cargo_size_bp, delivery=delivery,
            prefer_dsb_free=prefer_dsb_free, top_k=top_k,
        )
