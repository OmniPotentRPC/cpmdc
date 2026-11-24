`libcpmdc_gromacs` refuses a `-qmi` input whose `&DFT` `CHARGE` or
`MULTIPLICITY` differs from the charge or multiplicity of the
`CPMDParams` message in `&CPMDC PARAMS`, as it already did for the input
that GROMACS generates. The explicit route ran the message's system and
ignored the values that GROMACS used for the QM charge. A value that is
not an integer is refused too. `qmi_input()` and `write_qmi()` take the
deck's `charge` and `multiplicity` for a deck that is not neutral and
singlet. A link hydrogen is not counted into either value.
