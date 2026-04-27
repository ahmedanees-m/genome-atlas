"""GENOME-ATLAS: knowledge graph for programmable genome-writing enzymes.

Public API surface:
    __version__          - package version string
    Atlas                - main query/selection class (genome_atlas.api)
    get_graph            - PyG HeteroData builder with graph_view= parameter
                           (genome_atlas.graph.view)
    load_systems()       - load foundational systems YAML -> dict[str, SystemEntry]
    resolve_system_name  - resolve deprecated alias -> canonical name (e.g. IS622->ISCro4)
    SystemEntry          - frozen dataclass for per-system metadata
    select_editor        - convenience wrapper: Atlas().select_editor(...)
                           (genome_atlas.api.Atlas.select_editor)
"""

try:
    from genome_atlas._version import __version__
except ImportError:
    __version__ = "unknown"

# Lightweight imports (no torch dependency) - always available
from genome_atlas.systems import load_systems, resolve_system_name, SystemEntry

# Lazy public-API re-exports - imported here so users can write:
#   from genome_atlas import Atlas, get_graph
# Heavy torch / torch-geometric deps are only pulled in when these are used.
def __getattr__(name: str):  # PEP 562 lazy module attribute
    if name == "Atlas":
        from genome_atlas.api import Atlas
        return Atlas
    if name == "get_graph":
        from genome_atlas.graph.view import get_graph
        return get_graph
    raise AttributeError(f"module 'genome_atlas' has no attribute {name!r}")


__all__ = [
    "__version__",
    "Atlas",
    "get_graph",
    "load_systems",
    "resolve_system_name",
    "SystemEntry",
]
