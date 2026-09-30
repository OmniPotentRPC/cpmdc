#!/usr/bin/env python3
import os
import re
import sys


def _meson_project_version():
    """The project version is the version field of the root meson.build."""
    meson_build = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir, "meson.build"
    )
    with open(meson_build, encoding="ascii") as handle:
        for line in handle:
            match = re.match(r"\s*version:\s*'([^']+)'\s*,?\s*$", line)
            if match:
                return match.group(1)
    raise RuntimeError("project version is missing from meson.build")


project = "cpmdc"
copyright = "2026-present, cpmdc developers"
author = "Rohit Goswami"
release = _meson_project_version()
version = release

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
