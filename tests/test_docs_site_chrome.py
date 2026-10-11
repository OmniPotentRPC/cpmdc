#!/usr/bin/env python3
"""Check the published-site chrome: logos, ecosystem menu, and docs CI."""

from __future__ import annotations

import re
import sys
from pathlib import Path


ECOSYSTEM = [
    (
        "eOn",
        "https://eondocs.org",
        "Saddle point finding on potential energy surfaces",
    ),
    (
        "rgpot",
        "https://rgpot.rgoswami.me",
        "Potential-energy library and Cap'n Proto RPC server",
    ),
    (
        "rgmin",
        "https://omnipotentrpc.github.io/rgmin",
        "Local minimization over eindir objectives",
    ),
    (
        "rgsaddle",
        "https://omnipotentrpc.github.io/rgsaddle",
        "Band and minimum-mode saddle mechanics over rgmin",
    ),
    (
        "rgpycrumbs",
        "https://rgpycrumbs.rgoswami.me",
        "Pure-Python computational library and CLI for chemical physics",
    ),
    (
        "minimage",
        "https://lode-org.github.io/minimage",
        "Minimum-image convention on a general parallelepiped",
    ),
    (
        "linkcell",
        "https://docs.rs/linkcell",
        "Periodic linked-cell k-nearest neighbour search",
    ),
    (
        "readcon",
        "https://lode-org.github.io/readcon-core",
        "Reference implementation of versioned CON checkpoints",
    ),
    (
        "readcon-db",
        "https://lode-org.github.io/readcon-db",
        "Mmap-backed CON campaign corpus store",
    ),
    (
        "readcon-ops",
        "https://github.com/lode-org/readcon-ops",
        "Comparisons and rigid alignment of coordinate frames",
    ),
]


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: test_docs_site_chrome.py SOURCE_ROOT", file=sys.stderr)
        return 2

    root = Path(sys.argv[1])
    failures: list[str] = []

    conf = (root / "docs/source/conf.py").read_text(encoding="utf-8")
    for needle in (
        '"light_logo": "_static/icon.svg"',
        '"dark_logo": "_static/icon.svg"',
        'html_favicon = "_static/favicons/favicon.svg"',
        "sidebars/localtoc.html",
        "sidebars/repo-stats.html",
        "sidebars/edit-this-page.html",
        'sitemap_url_scheme = "{link}"',
    ):
        if needle not in conf:
            failures.append(f"conf.py missing {needle}")
    menu = re.findall(
        r'\{\s*"title": "([^"]+)",\s*"url": "([^"]+)",\s*"summary": "([^"]+)",',
        conf,
    )
    if menu != ECOSYSTEM:
        failures.append(f"ecosystem menu is {menu}, expected {ECOSYSTEM}")
    if '"title": "cpmdc"' in conf:
        failures.append("conf.py lists cpmdc in its own ecosystem menu")
    if "HaoZeke/readcon" in conf:
        failures.append("conf.py still points readcon at HaoZeke")

    for relative in (
        "docs/source/_static/icon.svg",
        "docs/source/_static/logo-light.svg",
        "docs/source/_static/logo-dark.svg",
        "docs/source/_static/logo-light.webp",
        "docs/source/_static/logo-dark.webp",
        "branding/logo/icon.svg",
        "branding/logo/icon-16.svg",
        "branding/logo/icon-32.svg",
        "branding/logo/logo-light.svg",
        "branding/logo/logo-dark.svg",
        "branding/logo/logo-light.webp",
        "branding/logo/logo-dark.webp",
        "branding/logo/favicon/favicon.svg",
        "docs/source/_static/favicons/favicon.svg",
        "docs/source/_static/favicons/favicon-16x16.png",
        "docs/source/_static/favicons/favicon-32x32.png",
        "docs/source/_static/favicons/apple-touch-icon.png",
        "docs/source/_static/favicons/android-chrome-192x192.png",
        "docs/source/_static/favicons/android-chrome-512x512.png",
        "docs/source/_static/favicons/site.webmanifest",
    ):
        if not (root / relative).is_file():
            failures.append(f"missing logo {relative}")

    head = (root / "docs/source/_templates/partials/extra-head.html").read_text(
        encoding="utf-8"
    )
    for needle in (
        "favicon-16x16.png",
        "favicon-32x32.png",
        "android-chrome-192x192.png",
        "android-chrome-512x512.png",
        "apple-touch-icon.png",
        "site.webmanifest",
    ):
        if needle not in head:
            failures.append(f"extra-head.html missing {needle}")

    footer = (root / "docs/source/_templates/partials/site-foot.html").read_text(
        encoding="utf-8"
    )
    for needle in ("antics.turtletech.us", "TurtleTech ehf"):
        if needle not in footer:
            failures.append(f"site-foot.html missing {needle}")

    gitignore = (root / ".gitignore").read_text(encoding="utf-8")
    if "/docs/source/**/*.rst" not in gitignore:
        failures.append("gitignore does not ignore generated rst")

    export_el = (root / "docs/export.el").read_text(encoding="utf-8")
    for needle in ("ox-rst", "CPMDC_ORG_EXPORT", "cpmdc-publish-with-pandoc"):
        if needle not in export_el:
            failures.append(f"export.el missing {needle}")
    if 'getenv "CPMDC_ORG_EXPORT"' not in export_el.replace("(", " ").replace(")", " "):
        # The gate is (getenv "CPMDC_ORG_EXPORT").
        if '(getenv "CPMDC_ORG_EXPORT")' not in export_el:
            failures.append("export.el does not gate pandoc on CPMDC_ORG_EXPORT")

    workflow = (root / ".github/workflows/ci_docs.yml").read_text(encoding="utf-8")
    build, marker, deploy = workflow.partition("deploy_docs:")
    if not marker:
        failures.append("ci_docs.yml has no deploy_docs job")
    else:
        if "contents: write" in build:
            failures.append("ci_docs.yml grants contents: write on the build job")
        if "contents: write" not in deploy:
            failures.append("ci_docs.yml deploy job lacks contents: write")
        if "github.event_name == 'push'" not in deploy:
            failures.append("ci_docs.yml deploy job is not push-only")
    for use in re.findall(r"uses:\s*(\S+)", workflow):
        if not re.search(r"@[0-9a-f]{40}\b", use):
            failures.append(f"ci_docs.yml action is not SHA-pinned: {use}")
    if "cname: cpmdc.rgoswami.me" not in workflow:
        failures.append("ci_docs.yml lost the cpmdc CNAME")

    commenter = (root / ".github/workflows/ci_doc_commenter.yml").read_text(
        encoding="utf-8"
    )
    for needle in (
        "Build and Deploy Documentation",
        "HaoZeke/doc-previewer@",
        "artifact_name: documentation",
    ):
        if needle not in commenter:
            failures.append(f"ci_doc_commenter.yml missing {needle}")
    for use in re.findall(r"uses:\s*(\S+)", commenter):
        if not re.search(r"@[0-9a-f]{40}\b", use):
            failures.append(f"ci_doc_commenter.yml action is not SHA-pinned: {use}")

    if failures:
        print("docs site chrome check failed:", file=sys.stderr)
        print("\n".join(failures), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
