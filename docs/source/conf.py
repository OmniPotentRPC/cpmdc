#!/usr/bin/env python3
import os
import sys

project = "cpmdc"
copyright = "2026-present, cpmdc developers"
author = "Rohit Goswami"
release = "0.1.0"
version = "0.1.0"

doxyrest_prefix = os.environ.get("CONDA_PREFIX")
extensions = [
    "myst_parser",
    "sphinx.ext.intersphinx",
    "sphinx_sitemap",
]
if doxyrest_prefix:
    sys.path.insert(0, os.path.join(doxyrest_prefix, "share", "doxyrest", "sphinx"))
    try:
        import doxyrest  # noqa: F401
        extensions = ["doxyrest", "cpplexer"] + extensions
    except ImportError:
        pass

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "xml"]
source_suffix = {".rst": "restructuredtext", ".md": "markdown"}

intersphinx_mapping = {"python": ("https://docs.python.org/3", None)}

html_theme = "shibuya"
html_static_path = ["_static"]
html_baseurl = "https://cpmdc.rgoswami.me/"

html_theme_options = {
    "github_url": "https://github.com/OmniPotentRPC/cpmdc",
    "accent_color": "violet",
    "dark_code": True,
    "globaltoc_expand_depth": 1,
}

html_context = {
    "source_type": "github",
    "source_user": "OmniPotentRPC",
    "source_repo": "cpmdc",
    "source_version": "main",
    "source_docs_path": "/docs/source/",
}

html_css_files = ["custom.css"]
html_favicon = "_static/mark.svg"
