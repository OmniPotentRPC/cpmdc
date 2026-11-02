The version is written by hand only in `version.txt`.
Meson and Sphinx read it, and `pixi run version-sync` regenerates the other copies with cog.
CI fails when a copy differs.
