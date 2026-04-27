import os
import sys

sys.path.insert(0, os.path.abspath(".."))

project = "GENOME-ATLAS"
copyright = "2026, PEN-STACK Consortium"
author = "PEN-STACK Consortium"
release = "0.5.0"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
html_theme = "alabaster"
html_static_path = ["_static"]
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "torch": ("https://pytorch.org/docs/stable", None),
}
