import os
import sys

sys.path.insert(0, os.path.abspath(".."))

project = "GENOME-ATLAS"
copyright = "2026, Anees Ahmed"
author = "Anees Ahmed"
release = "0.6.0"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "sphinx_autodoc_typehints",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
html_theme = "furo"
html_static_path = []          # _static dir not used; avoids sphinx warning

# Mock heavy optional deps so autodoc doesn't fail if they're unavailable
autodoc_mock_imports = ["torch_geometric"]

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "torch": ("https://pytorch.org/docs/stable", None),
}
