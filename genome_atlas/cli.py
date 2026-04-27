"""CLI entry point for genome-atlas."""
import json
from pathlib import Path

import click

from genome_atlas.api import Atlas


@click.group()
@click.option("--graph", "-g", type=click.Path(exists=True), help="Path to atlas.gpickle")
@click.option("--embeddings", "-e", type=click.Path(exists=True), help="Path to embeddings parquet")
@click.option("--targets", "-t", type=click.Path(exists=True), help="Path to targets parquet")
@click.pass_context
def cli(ctx, graph, embeddings, targets):
    """GENOME-ATLAS: Programmable Genome-Writing Knowledge Graph."""
    if ctx.invoked_subcommand is None:
        click.echo(ctx.get_help())
        return
    ctx.ensure_object(dict)
    # Defaults
    if graph is None:
        graph = Path.home() / "pen-stack/data/graphs/atlas.gpickle"
    else:
        graph = Path(graph)
    if targets is None:
        targets = Path.home() / "pen-stack/data/processed/targets_v2.parquet"
    else:
        targets = Path(targets)
    if embeddings is not None:
        embeddings = Path(embeddings)
    ctx.obj["atlas"] = Atlas(graph, embeddings, targets)


@cli.command()
@click.argument("system_id")
@click.pass_context
def query_system(ctx, system_id):
    """Query a system by node_id (e.g. System_SpCas9)."""
    result = ctx.obj["atlas"].query_system(system_id)
    click.echo(json.dumps(result, indent=2))


@cli.command()
@click.option("--cell", default="HEK293T", help="Cell type")
@click.option("--edit", default="deletion", help="Edit type")
@click.option("--cargo", default=0, type=int, help="Cargo size in bp")
@click.option("--delivery", default="AAV", help="Delivery vector")
@click.option("--prefer-dsb-free/--no-prefer-dsb-free", default=True)
@click.option("--top-k", default=5, type=int)
@click.pass_context
def select(ctx, cell, edit, cargo, delivery, prefer_dsb_free, top_k):
    """Rank genome-writing systems for a therapeutic scenario."""
    recs = ctx.obj["atlas"].select_editor(cell, edit, cargo, delivery, prefer_dsb_free, top_k)
    for i, r in enumerate(recs, 1):
        flag = "✓" if r.published_hit else " "
        click.echo(f"{flag} {i}. {r.system:20s}  pen_score={r.pen_score:.3f}  dsb_free={r.dsb_free}")


if __name__ == "__main__":
    cli()
