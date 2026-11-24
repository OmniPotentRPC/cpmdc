`libcpmdc_gromacs` refuses a `-qmi` input whose first `CHARGE` or
`MULTIPLICITY` anywhere in the file differs from the charge or
multiplicity of the `CPMDParams` message in `&CPMDC PARAMS`, matching
GROMACS' own token-stream reader just as it does for the input that
GROMACS generates. The explicit route ran the message's system and
ignored the values that GROMACS used for the QM charge. Later
`CHARGE` or `MULTIPLICITY` lines do not override the first occurrence,
and a line without a single value or a non-integer first value is
refused too. `qmi_input()` and `write_qmi()` take the deck's `charge`
and `multiplicity` for a deck that is not neutral and singlet. A link
hydrogen is not counted into either value.
