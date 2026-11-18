The sdist workflow installs the `test` extra before the installed-wheel check, so pytest is in that environment. The `ase` extra does not install pytest, and the check runs `python -m pytest`.
