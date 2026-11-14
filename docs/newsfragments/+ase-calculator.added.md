`tests/python` is a pytest suite for the ASE calculator on the reference evaluator: forces against finite differences, stress against numerical stress, units and sign, translation, caching, boundaries, `set()`, copies, and BFGS with a trajectory.
Meson runs it as `python-pytest`, the wheel builds run it against the installed package, and the `test` extra installs pytest.
The Python how-to gains a property table, a section on units, boundaries, and caching, and examples that run as written.
