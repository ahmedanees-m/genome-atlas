"""v0.7.2: ISCro4 canonical naming + IS622 alias resolution.

Verifies:
  - ISCro4 is the canonical top-level key (not IS622)
  - ISCro4 carries the correct UniProt accession and Pfam architecture
  - IS622 is retained as a deprecated alias in SystemEntry.aliases
  - resolve_system_name("IS622") returns "ISCro4" with a DeprecationWarning
  - load_systems() and resolve_system_name are importable from genome_atlas directly
"""
from __future__ import annotations

import warnings

import pytest

import genome_atlas
from genome_atlas.systems import load_systems, resolve_system_name


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def systems():
    """Return the bundled foundational-systems dict (loaded once per module)."""
    return load_systems()


# ---------------------------------------------------------------------------
# Tests - ISCro4 canonical entry
# ---------------------------------------------------------------------------

class TestISCro4Canonical:

    def test_iscro4_present_in_foundational_systems(self, systems):
        assert "ISCro4" in systems, "ISCro4 must be a top-level canonical key"

    def test_iscro4_uniprot(self, systems):
        assert systems["ISCro4"].uniprot == "D2TGM5"

    def test_iscro4_pfam_architecture(self, systems):
        pfam = systems["ISCro4"].pfam
        assert "PF01548" in pfam, "IS110 catalytic domain PF01548 missing"
        assert "PF02371" in pfam, "IS110 catalytic domain PF02371 missing"

    def test_iscro4_organism(self, systems):
        assert "Citrobacter rodentium" in systems["ISCro4"].organism

    def test_iscro4_tier_a_gate(self, systems):
        assert systems["ISCro4"].tier_a_gate is True

    def test_iscro4_mechanism_bucket(self, systems):
        assert systems["ISCro4"].mechanism_bucket == "DSB_FREE_TRANSEST_RECOMBINASE"


# ---------------------------------------------------------------------------
# Tests - IS622 alias
# ---------------------------------------------------------------------------

class TestIS622Alias:

    def test_is622_not_in_canonical_keys(self, systems):
        """IS622 must not be a top-level canonical name - only as alias inside ISCro4."""
        assert "IS622" not in systems, (
            "IS622 should not be a canonical key. "
            "It is a deprecated alias stored in ISCro4.aliases."
        )

    def test_iscro4_has_is622_alias(self, systems):
        assert "IS622" in systems["ISCro4"].aliases

    def test_resolve_is622_returns_iscro4(self, systems):
        with warnings.catch_warnings(record=True):
            warnings.simplefilter("always")
            canonical = resolve_system_name("IS622", systems, warn=False)
        assert canonical == "ISCro4"

    def test_resolve_is622_emits_deprecation_warning(self, systems):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            resolve_system_name("IS622", systems, warn=True)
        deprecation_warnings = [
            w for w in caught if issubclass(w.category, DeprecationWarning)
        ]
        assert deprecation_warnings, "Expected a DeprecationWarning for IS622 alias"
        assert "deprecated" in str(deprecation_warnings[0].message).lower()

    def test_resolve_canonical_name_no_warning(self, systems):
        """Resolving the canonical name ISCro4 must not emit any warnings."""
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            canonical = resolve_system_name("ISCro4", systems, warn=True)
        assert canonical == "ISCro4"
        deprecation_warnings = [
            w for w in caught if issubclass(w.category, DeprecationWarning)
        ]
        assert not deprecation_warnings, "No DeprecationWarning for canonical name"

    def test_unknown_system_raises_key_error(self, systems):
        with pytest.raises(KeyError):
            resolve_system_name("NotARealSystem_XYZ", systems)


# ---------------------------------------------------------------------------
# Tests - top-level genome_atlas import
# ---------------------------------------------------------------------------

class TestTopLevelImport:

    def test_load_systems_importable_from_package(self):
        """genome_atlas.load_systems is a callable re-exported from __init__."""
        assert callable(genome_atlas.load_systems)

    def test_resolve_system_name_importable_from_package(self):
        assert callable(genome_atlas.resolve_system_name)

    def test_system_entry_importable_from_package(self):
        from genome_atlas import SystemEntry
        assert SystemEntry is not None

    def test_package_load_systems_returns_iscro4(self):
        systems = genome_atlas.load_systems()
        assert "ISCro4" in systems
        assert systems["ISCro4"].uniprot == "D2TGM5"


# ---------------------------------------------------------------------------
# Tests - SystemEntry.uniprot edge cases
# ---------------------------------------------------------------------------

class TestSystemEntryUniprot:

    def test_uniprot_returns_none_for_multi_protein(self):
        """SystemEntry.uniprot returns None when proteins has >1 entry (line 114)."""
        from genome_atlas.systems import SystemEntry
        entry = SystemEntry(
            name="MultiProteinSystem",
            proteins=("P00001", "P00002"),  # two proteins
        )
        assert entry.uniprot is None

    def test_uniprot_returns_accession_for_single_protein(self):
        """SystemEntry.uniprot returns proteins[0] for single-protein systems."""
        from genome_atlas.systems import SystemEntry
        entry = SystemEntry(name="SpCas9", proteins=("Q99ZW2",))
        assert entry.uniprot == "Q99ZW2"

    def test_uniprot_returns_none_for_no_proteins(self):
        from genome_atlas.systems import SystemEntry
        entry = SystemEntry(name="NoProtein", proteins=())
        assert entry.uniprot is None


# ---------------------------------------------------------------------------
# Tests - resolve_system_name called without systems arg (line 220-221)
# ---------------------------------------------------------------------------

class TestResolveWithoutSystemsArg:

    def test_resolve_canonical_without_systems_arg(self):
        """resolve_system_name(name) with no systems= calls load_systems() internally."""
        canonical = resolve_system_name("ISCro4")
        assert canonical == "ISCro4"

    def test_resolve_alias_without_systems_arg(self):
        """Alias resolution via internal load_systems() call (line 220-221)."""
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            canonical = resolve_system_name("IS622")
        assert canonical == "ISCro4"
        dep_warnings = [w for w in caught if issubclass(w.category, DeprecationWarning)]
        assert dep_warnings

    def test_unknown_raises_without_systems_arg(self):
        with pytest.raises(KeyError):
            resolve_system_name("TOTALLY_UNKNOWN_SYSTEM_XYZ")
