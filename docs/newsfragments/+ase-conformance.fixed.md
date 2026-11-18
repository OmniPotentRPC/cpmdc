`CPMDCCalculator.calculate()` called directly replaces every stored result, so the stress of an earlier periodic step is not left beside the energy of a cluster, and a step that fails stores nothing instead of leaving the previous step under the new structure.
A step whose energy, forces, or stress is NaN or infinite raises `cpmdc.ase.CalculationFailed`, and so does a failure to open the session.
The Voigt stress averages the two off-diagonal elements of each shear component, as ASE's `full_3x3_to_voigt_6_stress` does, instead of taking the upper triangle.
