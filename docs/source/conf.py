#!/usr/bin/env python3
import os
import sys


def _project_version():
    """The project version is the single line of the root version.txt."""
    version_txt = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir, "version.txt"
    )
    with open(version_txt, encoding="ascii") as handle:
        return handle.read().strip()


project = "cpmdc"
copyright = "2026-present, cpmdc developers"
author = "Rohit Goswami"
release = _project_version()
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
# A warning is a failed book.
warningiserror = True


def _reference_well():
    """Draw the reference-well figure before Sphinx reads the how-to."""
    import importlib.util

    script = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), os.pardir, "render_reference_well.py"
    )
    spec = importlib.util.spec_from_file_location("render_reference_well", script)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    dest = os.path.join(os.path.dirname(os.path.abspath(__file__)), "_static", "reference-well.svg")
    if module.main(["--out", dest, "--check"]) != 0:
        raise SystemExit("reference well was not drawn")


def _api_stub():
    """Keep the C API toctree resolvable before Doxygen has run."""
    api = os.path.join(os.path.dirname(os.path.abspath(__file__)), "api", "index.rst")
    if os.path.isfile(api):
        return
    os.makedirs(os.path.dirname(api), exist_ok=True)
    with open(api, "w", encoding="utf-8") as handle:
        handle.write(
            "C API\n"
            "=====\n"
            "\n"
            "These pages are filled by Doxygen and Doxyrest.\n"
            "Until that step has run, the headers under ``include/`` are the API.\n"
        )


_reference_well()
_api_stub()
source_suffix = {".rst": "restructuredtext", ".md": "markdown"}

intersphinx_mapping = {"python": ("https://docs.python.org/3", None)}

html_theme = "shibuya"
html_static_path = ["_static"]
html_baseurl = "https://cpmdc.rgoswami.me/"
# Sphinx sets language to "en". The default sitemap scheme is
# {lang}{version}{link}, so every page would be published under /en/.
sitemap_url_scheme = "{link}"
# Navbar falls back to html_logo when light_logo/dark_logo are unset.
# The badge reads on both headers; the wordmarks are the page hero.
html_logo = "_static/icon.svg"

html_theme_options = {
    "github_url": "https://github.com/OmniPotentRPC/cpmdc",
    "accent_color": "violet",
    "dark_code": True,
    "globaltoc_expand_depth": 1,
    # Shibuya toggles .light-logo / .dark-logo with color mode.
    "light_logo": "_static/icon.svg",
    "dark_logo": "_static/icon.svg",
    "nav_links": [
        {
            "title": "Ecosystem",
            "children": [
                {
                    "title": "eOn",
                    "url": "https://eondocs.org",
                    "summary": "Saddle point finding on potential energy surfaces",
                },
                {
                    "title": "rgpot",
                    "url": "https://rgpot.rgoswami.me",
                    "summary": "Potential-energy library and Cap'n Proto RPC server",
                },
                {
                    "title": "rgmin",
                    "url": "https://omnipotentrpc.github.io/rgmin",
                    "summary": "Local minimization over eindir objectives",
                },
                {
                    "title": "rgsaddle",
                    "url": "https://omnipotentrpc.github.io/rgsaddle",
                    "summary": "Band and minimum-mode saddle mechanics over rgmin",
                },
                {
                    "title": "rgpycrumbs",
                    "url": "https://rgpycrumbs.rgoswami.me",
                    "summary": "Pure-Python computational library and CLI for chemical physics",
                },
                {
                    "title": "minimage",
                    "url": "https://lode-org.github.io/minimage",
                    "summary": "Minimum-image convention on a general parallelepiped",
                },
                {
                    "title": "linkcell",
                    "url": "https://docs.rs/linkcell",
                    "summary": "Periodic linked-cell k-nearest neighbour search",
                },
                {
                    "title": "readcon",
                    "url": "https://lode-org.github.io/readcon-core",
                    "summary": "Reference implementation of versioned CON checkpoints",
                },
                {
                    "title": "readcon-db",
                    "url": "https://lode-org.github.io/readcon-db",
                    "summary": "Mmap-backed CON campaign corpus store",
                },
                {
                    "title": "readcon-ops",
                    "url": "https://github.com/lode-org/readcon-ops",
                    "summary": "Comparisons and rigid alignment of coordinate frames",
                },
            ],
        },
    ],
}

html_context = {
    "source_type": "github",
    "source_user": "OmniPotentRPC",
    "source_repo": "cpmdc",
    "source_version": "main",
    "source_docs_path": "/docs/source/",
}

html_sidebars = {
    "**": [
        "sidebars/localtoc.html",
        "sidebars/repo-stats.html",
        "sidebars/edit-this-page.html",
    ],
}

html_css_files = ["custom.css"]
html_favicon = "_static/favicons/favicon.svg"
