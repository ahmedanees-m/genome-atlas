"""GENOME-ATLAS: knowledge graph for programmable genome-writing enzymes."""

try:
    from genome_atlas._version import __version__
except ImportError:
    __version__ = "unknown"

__all__ = ["__version__"]
