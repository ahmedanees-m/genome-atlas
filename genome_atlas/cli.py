"""CLI for genome-atlas."""
from __future__ import annotations

import click


@click.group()
@click.version_option()
def main() -> int:
    """GENOME-ATLAS CLI."""
    return 0


if __name__ == "__main__":
    main()
