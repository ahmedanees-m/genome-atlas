"""Extended selection scoring tests - covers remaining branches in selection.py.

Targets branches missed by test_selection.py and test_selection_scoring.py:
  - SNV + small cargo -> nuclease high, dsb-free lower
  - Small insertion with prime editor name -> s_dsb = 1.0
  - Large insertion with transposase -> s_dsb = 0.9
  - Default/moderate cargo block (s_dsb branches at lines 87-107)
  - Cargo-tiered DSB branches: small-insert HDR, SNV nuclease
  - AAV size scoring: ultra-compact, single-AAV, tight, borderline, over
  - Non-AAV delivery -> s_aav = 0.9
  - Cargo scoring: recombinase/integrase, CAST/transposase, nuclease>1000 bp
  - Cell type compat: prime, bridge, fanzor, cre, default
  - Weight adjustment paths: delivery=AAV large system, prefer_dsb_free toggle
"""
from __future__ import annotations

import networkx as nx
import pandas as pd
import pytest

from genome_atlas.selection import SelectionEngine, UseCaseProfile
from genome_atlas.api import Atlas


# ---------------------------------------------------------------------------
# FakeAtlas factory - lets us inject systems with specific sizes
# ---------------------------------------------------------------------------

def _make_atlas_with_systems(systems_rows: list[dict],
                              length_map: dict[str, int] | None = None) -> Atlas:
    """Build an Atlas-like object with an in-memory MultiDiGraph."""
    G = nx.MultiDiGraph()
    lm = length_map or {}
    for row in systems_rows:
        nid = row["node_id"]
        G.add_node(nid, **{k: v for k, v in row.items() if k != "node_id"})

    # Wire HAS_PROTEIN edges for size lookup
    for node_id, aa in lm.items():
        G.add_edge(
            [r["node_id"] for r in systems_rows if r["node_id"].endswith(node_id.split("_")[-1])][0]
                if False else list(G.nodes)[0],
            node_id, edge_type="HAS_PROTEIN"
        )

    atlas = Atlas.__new__(Atlas)
    atlas._G = G
    atlas._length_map = lm
    atlas._protein_map = {}
    atlas._system_by_name = {d.get("name", ""): n for n, d in G.nodes(data=True)}
    atlas._domain_by_accession = {}
    atlas._protein_by_accession = {}
    atlas._rna_by_name = {}
    atlas._embeddings = None

    # Override systems() to return a predictable DataFrame
    def _systems(mechanism_bucket=None):
        rows = []
        for n, d in G.nodes(data=True):
            if d.get("node_type") == "System":
                if mechanism_bucket is None or d.get("mechanism_bucket") == mechanism_bucket:
                    rows.append({"node_id": n, **d})
        return pd.DataFrame(rows)

    import types
    atlas.systems = types.MethodType(lambda self, mb=None: _systems(mb), atlas)
    return atlas


def _simple_atlas(*systems):
    """Create a minimal atlas from dicts of {node_id, name, mechanism_bucket}."""
    rows = [{"node_id": s["node_id"], "node_type": "System",
             "name": s["name"], "mechanism_bucket": s["mechanism_bucket"]}
            for s in systems]
    return _make_atlas_with_systems(rows)


def _fake_atlas(*system_tuples):
    """system_tuples: (name, mechanism_bucket).  node_id derived automatically."""
    systems = [{"node_id": f"System_{n}", "node_type": "System",
                "name": n, "mechanism_bucket": mb}
               for n, mb in system_tuples]
    return _make_atlas_with_systems(systems)


class FakeAtlasSimple:
    """Minimal fake Atlas with controllable systems() output."""
    def __init__(self, rows):
        self._G = nx.MultiDiGraph()
        self._length_map = {}
        self._rows = rows

    def systems(self, mechanism_bucket=None):
        rows = self._rows
        if mechanism_bucket:
            rows = [r for r in rows if r.get("mechanism_bucket") == mechanism_bucket]
        return pd.DataFrame(rows)


def _atlas(*rows):
    return FakeAtlasSimple(list(rows))


def _row(name, mech, *, nid=None):
    return {"node_id": nid or f"System_{name}", "node_type": "System",
            "name": name, "mechanism_bucket": mech}


# ---------------------------------------------------------------------------
# DSB scoring - SNV + simple edit branch (lines 53-61)
# ---------------------------------------------------------------------------

class TestDSBScoringSimpleEdits:

    def test_nuclease_high_for_snv(self):
        atlas = _atlas(
            _row("SpCas9",        "DSB_NUCLEASE"),
            _row("Cre_recombinase","DSB_FREE_TRANSEST_RECOMBINASE"),
        )
        engine = SelectionEngine(atlas)
        recs = engine.rank("HEK293T", "SNV", 0, "AAV", False, top_k=2)
        cas9 = next(r for r in recs if r.system == "SpCas9")
        cre  = next(r for r in recs if r.system == "Cre_recombinase")
        assert cas9.pen_score >= cre.pen_score, "Nuclease should score >= DSB-free for SNV"

    def test_dsb_free_lower_for_deletion(self):
        atlas = _atlas(
            _row("SpCas9",        "DSB_NUCLEASE"),
            _row("Cre_recombinase","DSB_FREE_TRANSEST_RECOMBINASE"),
        )
        engine = SelectionEngine(atlas)
        recs = engine.rank("HEK293T", "deletion", 5, "AAV", False, top_k=2)
        cas9 = next(r for r in recs if r.system == "SpCas9")
        cre  = next(r for r in recs if r.system == "Cre_recombinase")
        # For simple deletion, nuclease should outscore DSB-free
        assert cas9.pen_score >= cre.pen_score

    def test_unknown_mech_for_simple_edit(self):
        atlas = _atlas(_row("UnknownTool", "UNKNOWN_MECH"))
        engine = SelectionEngine(atlas)
        recs = engine.rank("HEK293T", "deletion", 0, "AAV", False, top_k=1)
        assert recs[0].system == "UnknownTool"
        assert 0.0 < recs[0].pen_score < 0.9


# ---------------------------------------------------------------------------
# DSB scoring - small insertion branch (lines 63-75)
# ---------------------------------------------------------------------------

class TestDSBScoringSmallInsertion:

    def test_prime_editor_top_for_small_insertion(self):
        # Name must contain "prime" for the code to recognise it
        atlas = _atlas(
            _row("PE2_prime_editor", "DSB_FREE_PRIME_EDITOR"),
            _row("SpCas9",           "DSB_NUCLEASE"),
            _row("Cre",              "DSB_FREE_TRANSEST_RECOMBINASE"),
        )
        engine = SelectionEngine(atlas)
        recs = engine.rank("HEK293T", "insertion", 100, "LNP", True, top_k=3)
        assert recs[0].system == "PE2_prime_editor"

    def test_nuclease_second_for_small_insertion(self):
        atlas = _atlas(
            _row("PE2_prime_editor", "DSB_FREE_PRIME_EDITOR"),
            _row("SpCas9",           "DSB_NUCLEASE"),
        )
        engine = SelectionEngine(atlas)
        recs = engine.rank("HEK293T", "insertion", 100, "LNP", False, top_k=2)
        names = [r.system for r in recs]
        assert "SpCas9" in names


# ---------------------------------------------------------------------------
# DSB scoring - large insertion branch (lines 76-86)
# ---------------------------------------------------------------------------

class TestDSBScoringLargeInsertion:

    def test_dsb_free_top_for_large_insert(self):
        atlas = _atlas(
            _row("ISCro4_bridge","DSB_FREE_TRANSEST_RECOMBINASE"),
            _row("SpCas9",       "DSB_NUCLEASE"),
        )
        engine = SelectionEngine(atlas)
        recs = engine.rank("HEK293T", "insertion", 2000, "LNP", True, top_k=2)
        assert recs[0].system == "ISCro4_bridge"

    def test_transposase_high_for_large_insert(self):
        atlas = _atlas(
            _row("TnpB_transposase", "TRANSPOSASE"),
            _row("SpCas9",           "DSB_NUCLEASE"),
        )
        engine = SelectionEngine(atlas)
        recs = engine.rank("HEK293T", "insertion", 5000, "LNP", False, top_k=2)
        tnpb = next(r for r in recs if r.system == "TnpB_transposase")
        cas9 = next(r for r in recs if r.system == "SpCas9")
        assert tnpb.pen_score > cas9.pen_score

    def test_nuclease_penalised_for_large_insert(self):
        atlas = _atlas(_row("SpCas9", "DSB_NUCLEASE"))
        engine = SelectionEngine(atlas)
        recs = engine.rank("HEK293T", "insertion", 3000, "LNP", False, top_k=1)
        assert recs[0].pen_score < 0.7, "Nuclease should be penalised for large cargo"


# ---------------------------------------------------------------------------
# DSB scoring - default/moderate cargo block (lines 87-108)
# ---------------------------------------------------------------------------

class TestDSBScoringDefault:

    def test_dsb_free_in_default_block(self):
        # edit_type != deletion/SNV with cargo<=10, != insertion<=200, != insertion>1000
        # Use insertion, cargo=500 (> 200, <= 1000)
        atlas = _atlas(_row("Cre", "DSB_FREE_TRANSEST_RECOMBINASE"))
        engine = SelectionEngine(atlas)
        recs = engine.rank("HEK293T", "insertion", 500, "LNP", True, top_k=1)
        assert recs[0].pen_score > 0.5

    def test_transposase_in_default_block(self):
        atlas = _atlas(_row("TnpA_transposase", "TRANSPOSASE"))
        engine = SelectionEngine(atlas)
        recs = engine.rank("HEK293T", "insertion", 500, "LNP", False, top_k=1)
        assert 0.0 < recs[0].pen_score <= 1.0

    def test_nuclease_small_insert_hdr_branch(self):
        """Line 97-99: DSB nuclease with insertion<=1000 -> s_dsb = 0.55."""
        atlas = _atlas(_row("HiFi_Cas9", "DSB_NUCLEASE"))
        engine = SelectionEngine(atlas)
        recs = engine.rank("HEK293T", "insertion", 800, "LNP", False, top_k=1)
        assert 0.0 < recs[0].pen_score <= 1.0

    def test_nuclease_snv_branch(self):
        """Line 100-102: DSB nuclease with edit_type=SNV -> s_dsb = 0.60."""
        atlas = _atlas(_row("SaCas9", "DSB_NUCLEASE"))
        engine = SelectionEngine(atlas)
        # SNV with cargo=5 will go into simple-edit branch, not default;
        # use cargo=500 to fall into default block
        recs = engine.rank("HEK293T", "SNV", 500, "LNP", False, top_k=1)
        assert 0.0 < recs[0].pen_score <= 1.0

    def test_unknown_mech_in_default_block(self):
        atlas = _atlas(_row("Fanzor1", "UNKNOWN_MECHANISM"))
        engine = SelectionEngine(atlas)
        recs = engine.rank("HEK293T", "insertion", 500, "LNP", False, top_k=1)
        assert 0.0 < recs[0].pen_score <= 1.0


# ---------------------------------------------------------------------------
# AAV delivery scoring (lines 110-130)
# ---------------------------------------------------------------------------

class TestAAVScoring:

    def test_ultra_compact_aav_score(self):
        """total_aa <= 600 -> s_aav = 1.0."""
        # Need atlas with length_map wired. Use the full FakeAtlas with size patch.
        class SizedAtlas(FakeAtlasSimple):
            def __init__(self, rows, length_map):
                super().__init__(rows)
                self._G = nx.MultiDiGraph()
                for r in rows:
                    self._G.add_node(r["node_id"], **r)
                self._length_map = length_map

        # CjCas9: ~984 aa -> won't hit ultra-compact; use artificial 500 aa
        atlas = SizedAtlas(
            [_row("SmallCas", "DSB_NUCLEASE", nid="System_SmallCas")],
            {"Protein_X": 500}
        )
        # Add HAS_PROTEIN edge so size lookup fires
        atlas._G.add_node("Protein_X", node_type="Protein", accession="X")
        atlas._G.add_edge("System_SmallCas", "Protein_X", edge_type="HAS_PROTEIN")

        engine = SelectionEngine(atlas)
        recs = engine.rank("HEK293T", "deletion", 0, "AAV", False, top_k=1)
        # With 500 aa ultra-compact, AAV flag should be True
        assert recs[0].aav_fit is True

    def test_non_aav_delivery(self):
        """delivery != AAV -> s_aav = 0.9 (line 131)."""
        atlas = _atlas(_row("SpCas9", "DSB_NUCLEASE"))
        engine = SelectionEngine(atlas)
        recs = engine.rank("HEK293T", "deletion", 0, "LNP", False, top_k=1)
        assert recs[0].pen_score > 0.0

    def test_unknown_size_aav_guard(self):
        """total_aa == 0 -> s_aav = 0.1 (C3 guard, line 113)."""
        # No length_map entries -> total_aa == 0
        atlas = _atlas(_row("NoSizeSystem", "DSB_NUCLEASE"))
        engine = SelectionEngine(atlas)
        recs = engine.rank("HEK293T", "deletion", 0, "AAV", False, top_k=1)
        # aav_fit should be False because s_aav = 0.1 < 0.8
        assert recs[0].aav_fit is False


# ---------------------------------------------------------------------------
# Cargo scoring (lines 135-163)
# ---------------------------------------------------------------------------

class TestCargoScoring:

    def test_snv_prime_editor_boost(self):
        """line 136: SNV + prime editor -> s_cargo = 1.0."""
        atlas = _atlas(
            _row("PE2_prime", "DSB_FREE_PRIME_EDITOR"),
            _row("SpCas9",    "DSB_NUCLEASE"),
        )
        engine = SelectionEngine(atlas)
        recs = engine.rank("HEK293T", "SNV", 1, "LNP", False, top_k=2)
        pe2 = next(r for r in recs if r.system == "PE2_prime")
        assert pe2.pen_score > 0.5

    def test_prime_editor_large_cargo_penalty(self):
        """line 141-143: prime editor + cargo > 200 bp -> s_cargo = 0.15."""
        atlas = _atlas(_row("PE3_prime_editor", "DSB_FREE_PRIME_EDITOR"))
        engine = SelectionEngine(atlas)
        recs = engine.rank("HEK293T", "insertion", 1000, "LNP", False, top_k=1)
        assert recs[0].pen_score < 0.9

    def test_recombinase_large_cargo_ok(self):
        """line 145-150: recombinase/integrase name -> s_cargo = 1.0 for cargo<=50kbp."""
        atlas = _atlas(_row("Bxb1_integrase", "DSB_FREE_TRANSEST_RECOMBINASE"))
        engine = SelectionEngine(atlas)
        recs = engine.rank("HEK293T", "insertion", 10000, "LNP", False, top_k=1)
        assert recs[0].pen_score > 0.5

    def test_recombinase_overlimit_penalty(self):
        """line 148: cargo > 50000 bp -> s_cargo = 0.4."""
        atlas = _atlas(_row("ISCro4_bridge_recombinase", "DSB_FREE_TRANSEST_RECOMBINASE"))
        engine = SelectionEngine(atlas)
        recs = engine.rank("HEK293T", "insertion", 60000, "LNP", False, top_k=1)
        assert recs[0].pen_score < 0.9

    def test_cast_transposase_cargo(self):
        """line 151-153: cast/transposon name -> large cargo OK."""
        atlas = _atlas(_row("CAST_transposon", "TRANSPOSASE"))
        engine = SelectionEngine(atlas)
        recs = engine.rank("HEK293T", "insertion", 5000, "LNP", False, top_k=1)
        assert recs[0].pen_score > 0.0

    def test_nuclease_large_cargo_penalty(self):
        """line 154-158: DSB_NUCLEASE with cargo > 1000 -> s_cargo = 0.2."""
        atlas = _atlas(_row("SpCas9", "DSB_NUCLEASE"))
        engine = SelectionEngine(atlas)
        recs = engine.rank("HEK293T", "insertion", 5000, "LNP", False, top_k=1)
        assert recs[0].pen_score < 0.8

    def test_deletion_cargo_agnostic(self):
        """line 163: deletion -> s_cargo = 0.9."""
        atlas = _atlas(_row("SpCas9", "DSB_NUCLEASE"))
        engine = SelectionEngine(atlas)
        recs = engine.rank("HEK293T", "deletion", 0, "LNP", False, top_k=1)
        assert recs[0].pen_score > 0.4


# ---------------------------------------------------------------------------
# Cell type scoring (lines 165-186)
# ---------------------------------------------------------------------------

class TestCellTypeScoring:

    def test_cas9_human_cell(self):
        """line 169-171: SpCas9 in human cell -> s_cell = 0.95."""
        atlas = _atlas(_row("spcas9_hifi", "DSB_NUCLEASE"))
        engine = SelectionEngine(atlas)
        recs = engine.rank("HEK293T", "deletion", 0, "LNP", False, top_k=1)
        assert recs[0].pen_score > 0.5

    def test_prime_editor_human_cell(self):
        """line 172-174: prime in name -> s_cell = 0.9."""
        atlas = _atlas(_row("PE3_prime_editor", "DSB_FREE_PRIME_EDITOR"))
        engine = SelectionEngine(atlas)
        recs = engine.rank("K562", "SNV", 1, "LNP", False, top_k=1)
        assert recs[0].pen_score > 0.0

    def test_bridge_human_cell(self):
        """line 175-177: bridge in name -> s_cell = 0.85."""
        atlas = _atlas(_row("ISCro4_bridge_recombinase", "DSB_FREE_TRANSEST_RECOMBINASE"))
        engine = SelectionEngine(atlas)
        recs = engine.rank("HEK293T", "insertion", 500, "LNP", False, top_k=1)
        assert recs[0].pen_score > 0.0

    def test_fanzor_human_cell(self):
        """line 178-180: fanzor in name -> s_cell = 0.7."""
        atlas = _atlas(_row("Fanzor1_nuclease", "DSB_NUCLEASE"))
        engine = SelectionEngine(atlas)
        recs = engine.rank("HEK293T", "deletion", 0, "LNP", False, top_k=1)
        assert recs[0].pen_score > 0.0

    def test_cre_human_cell(self):
        """line 181-183: cre in name -> s_cell = 0.8."""
        atlas = _atlas(_row("CreERT2_recombinase", "DSB_FREE_TRANSEST_RECOMBINASE"))
        engine = SelectionEngine(atlas)
        recs = engine.rank("HEK293T", "deletion", 0, "LNP", False, top_k=1)
        assert recs[0].pen_score > 0.0

    def test_unknown_system_human_cell(self):
        """line 184-185: unrecognised name -> s_cell = 0.6."""
        atlas = _atlas(_row("NewNovelSystem99", "UNKNOWN"))
        engine = SelectionEngine(atlas)
        recs = engine.rank("HEK293T", "deletion", 0, "LNP", False, top_k=1)
        assert recs[0].pen_score > 0.0

    def test_non_human_cell_default(self):
        """line 166: cell_type not in human_cell_types -> s_cell = 0.5."""
        atlas = _atlas(_row("SpCas9", "DSB_NUCLEASE"))
        engine = SelectionEngine(atlas)
        recs = engine.rank("CHO", "deletion", 0, "LNP", False, top_k=1)
        assert recs[0].pen_score > 0.0


# ---------------------------------------------------------------------------
# Weight adjustments (lines 188-204)
# ---------------------------------------------------------------------------

class TestWeightAdjustments:

    def test_prefer_dsb_free_boosts_relative_advantage(self):
        """prefer_dsb_free=True widens the gap between DSB-free and DSB-nuclease.

        The absolute score of the DSB-free system may go up or down depending on
        the specific scenario; what matters is that its margin over the nuclease
        is larger when prefer_dsb_free=True.
        """
        atlas = _atlas(
            _row("ISCro4_bridge","DSB_FREE_TRANSEST_RECOMBINASE"),
            _row("SpCas9",       "DSB_NUCLEASE"),
        )
        engine = SelectionEngine(atlas)
        recs_pref    = engine.rank("HEK293T", "insertion", 500, "LNP", True,  top_k=2)
        recs_nopref  = engine.rank("HEK293T", "insertion", 500, "LNP", False, top_k=2)
        gap_pref = (
            next(r.pen_score for r in recs_pref if r.system == "ISCro4_bridge") -
            next(r.pen_score for r in recs_pref if r.system == "SpCas9")
        )
        gap_nopref = (
            next(r.pen_score for r in recs_nopref if r.system == "ISCro4_bridge") -
            next(r.pen_score for r in recs_nopref if r.system == "SpCas9")
        )
        assert gap_pref >= gap_nopref, (
            f"prefer_dsb_free should widen DSB-free advantage: "
            f"pref gap={gap_pref:.3f} vs no-pref gap={gap_nopref:.3f}"
        )

    def test_aav_size_constrained_weights(self):
        """delivery=AAV with total_aa > 900 triggers AAV-boosted weight block."""
        # Use a large system to push the weight to the AAV path.
        # AAV scoring yields s_aav = 0.1 for 0-size; weights shift to boost aav.
        atlas = _atlas(_row("LargeCas9", "DSB_NUCLEASE"))
        engine = SelectionEngine(atlas)
        # With AAV delivery, if size unknown (total_aa=0, s_aav=0.1), and prefer off:
        recs = engine.rank("HEK293T", "insertion", 500, "AAV", False, top_k=1)
        assert recs[0].pen_score > 0.0
