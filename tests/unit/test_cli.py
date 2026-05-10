"""Unit tests for CLI."""
from click.testing import CliRunner

from genome_atlas.cli import cli


def test_cli_help():
    runner = CliRunner()
    result = runner.invoke(cli, ["--help"])
    assert result.exit_code == 0
    assert "GENOME-ATLAS" in result.output


def test_query_system_help():
    runner = CliRunner()
    result = runner.invoke(cli, ["query-system", "--help"])
    assert result.exit_code == 0


def test_select_help():
    runner = CliRunner()
    result = runner.invoke(cli, ["select", "--help"])
    assert result.exit_code == 0
