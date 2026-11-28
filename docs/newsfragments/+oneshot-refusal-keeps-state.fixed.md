`cpmdc_energy()`, `cpmdc_energy_gradient()`, and `cpmdc_energy_forces()` check
the geometry before installing parameters.
A call rejected before that install leaves the parameters, the active session,
the selected orbitals, and the previous result in place.
