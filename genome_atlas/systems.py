"""Loader for the GENOME-ATLAS foundational systems YAML.

Public API
----------
load_systems() -> dict[str, SystemEntry]
    Parse ``genome_atlas/data/foundational_systems.yaml`` and return a dict
    keyed by canonical system name (e.g. "ISCro4", "SpCas9").

SystemEntry
    Frozen dataclass exposing the per-system metadata fields used by
    PEN-COMPARE v3.2 and downstream packages.

Alias resolution
----------------
``load_systems()`` returns the *canonical* name as the dict key.  Aliases are
stored in ``SystemEntry.aliases`` but are **not** top-level keys.  Callers that
need alias-aware lookup should use :func:`resolve_system_name`.

Example::

    >>> from genome_atlas.systems import load_systems, resolve_system_name
    >>> systems = load_systems()
    >>> "ISCro4" in systems
    True
    >>> systems["ISCro4"].uniprot
    'D2TGM5'
    >>> systems["ISCro4"].pfam
    ['PF01548', 'PF02371']
    >>> "IS622" in systems["ISCro4"].aliases
    True
    >>> resolve_system_name("IS622", systems)
    'ISCro4'
"""
from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import yaml

_YAML_PATH = Path(__file__).parent / "data" / "foundational_systems.yaml"


# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class SystemEntry:
    """Immutable record for one foundational genome-writing system.

    Attributes
    ----------
    name:
        Canonical system identifier (e.g. ``"ISCro4"``).
    aliases:
        Deprecated / alternative names for the same system (e.g. ``["IS622"]``).
    type:
        High-level system family (``"Bridge_Recombinase"``, ``"CRISPR-Cas"``, ...).
    subtype:
        Subfamily within the type (e.g. ``"IS110"``, ``"Type II-A"``).
    mechanism_bucket:
        Canonical mechanism class string (e.g. ``"DSB_FREE_TRANSEST_RECOMBINASE"``).
    proteins:
        List of UniProt accession strings for the catalytic subunit(s).
    organism:
        Source organism string (e.g. ``"Citrobacter rodentium ICC168"``).
    pfam:
        List of Pfam family IDs covering the catalytic domain(s) of the primary
        protein.  Populated for systems where the Pfam architecture is
        diagnostically significant (e.g. IS110 bridge recombinases).
    rna_components:
        RNA guide/scaffold names referenced in ``rnas:`` block of the YAML.
    canonical_structures:
        PDB / AlphaFold accessions for the primary structural model.
    reference_doi:
        DOI of the primary canonical publication.
    reference_doi_2:
        DOI of a secondary reference (e.g. preprint superseded by reference_doi).
    notes:
        Free-text provenance / naming notes.
    tier_a_gate:
        If ``True``, mech-class Tier-A rule fires unconditionally for this system's
        primary protein based on Pfam architecture.
    """

    name: str
    aliases: tuple[str, ...] = field(default_factory=tuple)
    type: Optional[str] = None
    subtype: Optional[str] = None
    mechanism_bucket: Optional[str] = None
    proteins: tuple[str, ...] = field(default_factory=tuple)
    organism: Optional[str] = None
    pfam: tuple[str, ...] = field(default_factory=tuple)
    rna_components: tuple[str, ...] = field(default_factory=tuple)
    canonical_structures: tuple[str, ...] = field(default_factory=tuple)
    reference_doi: Optional[str] = None
    reference_doi_2: Optional[str] = None
    notes: Optional[str] = None
    tier_a_gate: bool = False

    @property
    def uniprot(self) -> Optional[str]:
        """Primary UniProt accession - returns ``proteins[0]`` if exactly one entry.

        Returns ``None`` for multi-protein systems (e.g. CAST complexes) where a
        single accession is ambiguous.
        """
        if len(self.proteins) == 1:
            return self.proteins[0]
        return None


# ---------------------------------------------------------------------------
# Loader
# ---------------------------------------------------------------------------


def load_systems(yaml_path: Optional[Path] = None) -> dict[str, SystemEntry]:
    """Load foundational systems from YAML, keyed by canonical system name.

    Parameters
    ----------
    yaml_path:
        Path to a ``foundational_systems.yaml``-formatted file.  Defaults to
        the bundled data file inside the installed ``genome_atlas`` package.

    Returns
    -------
    dict[str, SystemEntry]
        Maps canonical name -> :class:`SystemEntry`.  Aliases are **not**
        top-level keys; use :func:`resolve_system_name` for alias-aware lookup.

    Raises
    ------
    FileNotFoundError
        If ``yaml_path`` does not exist.
    """
    path = yaml_path or _YAML_PATH
    with open(path, "r", encoding="utf-8") as fh:
        doc = yaml.safe_load(fh)

    result: dict[str, SystemEntry] = {}
    for raw in doc.get("systems", []):
        name = raw.get("name")
        if not name:
            continue

        # Aliases: stored as YAML list, exposed as a frozen tuple
        raw_aliases = raw.get("aliases", [])
        if isinstance(raw_aliases, str):
            raw_aliases = [raw_aliases]

        entry = SystemEntry(
            name=name,
            aliases=tuple(raw_aliases),
            type=raw.get("type"),
            subtype=raw.get("subtype"),
            mechanism_bucket=raw.get("mechanism_bucket"),
            proteins=tuple(raw.get("proteins") or []),
            organism=raw.get("organism"),
            pfam=tuple(raw.get("pfam") or []),
            rna_components=tuple(raw.get("rna_components") or []),
            canonical_structures=tuple(raw.get("canonical_structures") or []),
            reference_doi=raw.get("reference_doi"),
            reference_doi_2=raw.get("reference_doi_2"),
            notes=raw.get("notes"),
            tier_a_gate=bool(raw.get("tier_a_gate", False)),
        )
        result[name] = entry

    return result


# ---------------------------------------------------------------------------
# Alias resolution
# ---------------------------------------------------------------------------


def resolve_system_name(
    name: str,
    systems: Optional[dict[str, SystemEntry]] = None,
    *,
    warn: bool = True,
) -> str:
    """Resolve a (possibly deprecated) system name to its canonical form.

    Parameters
    ----------
    name:
        System name to resolve (canonical or deprecated alias).
    systems:
        Pre-loaded systems dict from :func:`load_systems`.  If ``None``,
        :func:`load_systems` is called internally (slightly slower on repeated calls).
    warn:
        If ``True`` (default), emit a :class:`DeprecationWarning` when ``name``
        is a deprecated alias rather than the canonical name.

    Returns
    -------
    str
        Canonical system name.

    Raises
    ------
    KeyError
        If ``name`` is neither a canonical name nor a known alias.

    Examples
    --------
    >>> resolve_system_name("IS622")
    # DeprecationWarning: 'IS622' is a deprecated alias for 'ISCro4'. Use 'ISCro4'.
    'ISCro4'
    >>> resolve_system_name("ISCro4")
    'ISCro4'
    """
    if systems is None:
        systems = load_systems()

    # Direct canonical match
    if name in systems:
        return name

    # Alias scan
    for canonical, entry in systems.items():
        if name in entry.aliases:
            if warn:
                warnings.warn(
                    f"'{name}' is a deprecated alias for '{canonical}'. "
                    f"Use '{canonical}' in new code.",
                    DeprecationWarning,
                    stacklevel=2,
                )
            return canonical

    raise KeyError(
        f"System '{name}' not found in foundational_systems.yaml "
        f"(neither as canonical name nor alias)."
    )
