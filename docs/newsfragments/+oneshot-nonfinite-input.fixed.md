`cpmdc_energy()`, `cpmdc_energy_gradient()`, and `cpmdc_energy_forces()` reject
a position that is not finite before the evaluation.
`ok` is 0, the energy is 0, and
`cpmdc_last_error_code()` is `CPMDC_ERROR_INVALID_ARGUMENT`.
The message is `positions contain a non-finite component`.
