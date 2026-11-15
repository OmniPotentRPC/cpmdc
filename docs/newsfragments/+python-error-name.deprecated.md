`cpmdc.CpmdcError` is deprecated in favour of `cpmdc.CPMDCError`, which spells the acronym as `CPMDCCalculator` and the C types (`CPMDCSession`, `CPMDCResult`) do.
Both names are the same class, so `except cpmdc.CpmdcError` still catches every error the package raises; reading the old name warns with `DeprecationWarning`, and it is no longer in `cpmdc.__all__`.
