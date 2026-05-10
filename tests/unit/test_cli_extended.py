"""Extended CLI tests - exercises query_system and select commands.

Covers cli.py lines missed by test_cli.py:
  - lines 17-19: no-subcommand help path
  - lines 25, 29: default graph/targets path assignment
  - lines 31: embeddings Path conversion
  - lines 40-41: query_system result output
  - lines 54-57: select command result formatting
"""
from __future__ import annotations

import json
import pickle
from pathlib import Path

import networkx as nx
import pandas as pd
import pytest
from click.testing import CliRunner

from genome_atlas.cli import cli


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_mock_graph() -> nx.MultiDiGraph:
    G = nx.MultiDiGraph()
    G.add_node("System_SpCas9", node_type="System", name="SpCas9",
               mechanism_bucket="DSB_NUCLEASE")
    G.add_node("Protein_Q99ZW2", node_type="Protein", accession="Q99ZW2",
               name="Cas9")
    G.add_edge("System_SpCas9", "Protein_Q99ZW2", edge_type="HAS_PROTEIN")
    return G


def _setup_atlas_files(tmp_path: Path):
    gp = tmp_path / "atlas.gpickle"
    tp = tmp_path / "targets.parquet"
    with open(gp, "wb") as f:
        pickle.dump(_make_mock_graph(), f)
    df = pd.DataFrame([{
        "accession": "Q99ZW2", "length": 1368,
        "protein_name": "Cas9", "organism_name": "S. pyogenes"
    }])
    df.to_parquet(tp, index=False)
    return gp, tp


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestCliNoSubcommand:
    def test_no_subcommand_prints_help(self):
        """Hitting cli with no subcommand prints help (lines 17-18)."""
        runner = CliRunner()
        result = runner.invoke(cli, [])
        assert result.exit_code == 0
        assert "GENOME-ATLAS" in result.output


class TestQuerySystemCommand:

    def test_query_system_outputs_json(self, tmp_path):
        gp, tp = _setup_atlas_files(tmp_path)
        runner = CliRunner()
        result = runner.invoke(cli, [
            "--graph", str(gp),
            "--targets", str(tp),
            "query-system", "SpCas9",
        ])
        assert result.exit_code == 0, result.output
        data = json.loads(result.output)
        assert data["name"] == "SpCas9"

    def test_query_system_with_default_paths(self, tmp_path, monkeypatch):
        """cli uses Path.home() defaults when --graph/--targets not supplied (lines 25, 29).

        We patch genome_atlas.cli.Atlas (where it was already imported) to avoid
        needing real home-dir files.
        """
        from unittest.mock import patch, MagicMock
        import genome_atlas.cli as cli_module

        mock_atlas = MagicMock()
        mock_atlas.query_system.return_value = {
            "node_id": "System_SpCas9", "name": "SpCas9"
        }

        with patch.object(cli_module, "Atlas", return_value=mock_atlas) as MockAtlas:
            runner = CliRunner()
            result = runner.invoke(cli, ["query-system", "SpCas9"])
            # Atlas was called without explicit graph/targets
            MockAtlas.assert_called_once()
            call_args = MockAtlas.call_args[0]
            # first positional arg is the graph path - should be the default home path
            assert "pen-stack" in str(call_args[0])

    def test_query_system_with_embeddings_path(self, tmp_path):
        """--embeddings flag triggers Path(embeddings) conversion (line 31)."""
        gp, tp = _setup_atlas_files(tmp_path)
        # Create a minimal embeddings parquet so Click --path exists check passes
        ep = tmp_path / "embeddings.parquet"
        pd.DataFrame(columns=["node_id", "embedding"]).to_parquet(ep, index=False)

        runner = CliRunner()
        result = runner.invoke(cli, [
            "--graph", str(gp),
            "--targets", str(tp),
            "--embeddings", str(ep),
            "query-system", "SpCas9",
        ])
        # Just check it ran (embeddings may be empty so no embedding lookup needed)
        assert result.exit_code == 0


class TestSelectCommand:

    def test_select_outputs_ranked_lines(self, tmp_path):
        gp, tp = _setup_atlas_files(tmp_path)
        runner = CliRunner()
        result = runner.invoke(cli, [
            "--graph", str(gp),
            "--targets", str(tp),
            "select",
            "--cell", "HEK293T",
            "--edit", "deletion",
            "--cargo", "0",
            "--delivery", "AAV",
            "--top-k", "1",
        ])
        assert result.exit_code == 0, result.output
        assert "pen_score=" in result.output

    def test_select_shows_mechanism(self, tmp_path):
        gp, tp = _setup_atlas_files(tmp_path)
        runner = CliRunner()
        result = runner.invoke(cli, [
            "--graph", str(gp),
            "--targets", str(tp),
            "select",
            "--top-k", "2",
        ])
        assert result.exit_code == 0, result.output
        assert "mechanism=" in result.output

    def test_select_no_prefer_dsb_free(self, tmp_path):
        gp, tp = _setup_atlas_files(tmp_path)
        runner = CliRunner()
        result = runner.invoke(cli, [
            "--graph", str(gp),
            "--targets", str(tp),
            "select",
            "--no-prefer-dsb-free",
            "--top-k", "1",
        ])
        assert result.exit_code == 0, result.output
