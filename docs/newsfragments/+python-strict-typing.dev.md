`python/cpmdc` passes `mypy --strict` with the `[tool.mypy]` table of `pyproject.toml`, and the `python-mypy-strict` Meson test runs it (it is skipped when mypy is not installed; the sanitizer and wheel jobs of the workflow install it).
`cpmdc.gromacs.write_qmi` now lists the keyword options of `qmi_input` in its signature instead of accepting `**options`; the accepted keywords are unchanged.
