`cpmdc.gromacs.qmi_input` and `write_qmi` reject a non-finite or non-positive `ewald_rcut` or `ewald_precision`, a non-finite `qm_charges` entry, and a `charge` or `multiplicity` that is not an integer, with a message that names the option.
A `multiplicity` below 1 raises `ValueError`.
They used to write `nan` or `inf` into the file and truncate `charge=1.7` to 1.
