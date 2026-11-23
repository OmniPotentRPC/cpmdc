An older `libcpmdc` without `cpmdc_params_method` still loads in Python;
calling `params_method()` against it raises `CPMDCError` reporting the
missing symbol.
