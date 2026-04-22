"""Placeholder test — ensures CI passes until real tests are written."""
import genome_atlas


def test_version():
    """Verify package has a version attribute."""
    assert hasattr(genome_atlas, "__version__")
    assert genome_atlas.__version__ == "0.0.1"
