A repeated session step keeps the gradient snapshot in the caller's atom order.
`cpmdc_session_last_gradient()` and a later cache hit match the forces when the
atoms are reordered.
