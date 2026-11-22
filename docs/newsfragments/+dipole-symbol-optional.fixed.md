An older `libcpmdc` without `cpmdc_session_last_dipole` still loads at the
same ABI. Energy and forces are unchanged, and the dipole stays unavailable.
Asking ASE for the dipole names the missing symbol.
