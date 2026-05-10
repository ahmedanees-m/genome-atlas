"""Tests for genome_atlas.__init__ - covers lazy import paths.

Covers __init__.py lines missed by other test modules:
  - lines 28-29: PEP-562 __getattr__ returning Atlas class
  - lines 31-32: PEP-562 __getattr__ returning get_graph function
  - line 33: __getattr__ raising AttributeError for unknown attribute
"""
from __future__ import annotations

import pytest


class TestLazyImports:

    def test_atlas_accessible_from_package(self):
        """Accessing genome_atlas.Atlas triggers the __getattr__ lazy import."""
        import genome_atlas
        Atlas = genome_atlas.Atlas
        assert Atlas.__name__ == "Atlas"

    def test_atlas_importable_directly(self):
        from genome_atlas import Atlas
        assert callable(Atlas)

    def test_get_graph_accessible_from_package(self):
        """Accessing genome_atlas.get_graph triggers the __getattr__ lazy import.

        genome_atlas/graph/view.py uses lazy function-level torch imports,
        so the module itself is importable without torch installed.
        """
        import genome_atlas
        get_graph = genome_atlas.get_graph
        assert callable(get_graph)

    def test_unknown_attribute_raises(self):
        """__getattr__ must raise AttributeError for unknown names."""
        import genome_atlas
        with pytest.raises(AttributeError):
            _ = genome_atlas.this_does_not_exist_anywhere

    def test_version_string_present(self):
        import genome_atlas
        assert isinstance(genome_atlas.__version__, str)
        assert genome_atlas.__version__ != ""

    def test_all_exports_present(self):
        import genome_atlas
        for name in ["__version__", "load_systems", "resolve_system_name", "SystemEntry"]:
            assert hasattr(genome_atlas, name), f"Missing export: {name}"
