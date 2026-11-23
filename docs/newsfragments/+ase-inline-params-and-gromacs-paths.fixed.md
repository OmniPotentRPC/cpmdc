ASE and GROMACS adapter ergonomics were tightened for drop-in use:

- `CPMDCCalculator.todict()` now keeps in-memory `params` bytes by emitting a
  JSON-safe `{"cpmdparams_base64": ...}` mapping, and the calculator accepts
  that mapping in `params` for a lossless rebuild path.
- `cpmdc.gromacs.qmi_input()` and `write_qmi()` now validate `params` and output
  paths up front, with typed errors for binary or empty path inputs.
