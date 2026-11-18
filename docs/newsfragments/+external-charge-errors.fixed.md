External-charge force reads now set `cpmdc_last_error()` and
`cpmdc_last_error_code()`. A deck that already defines `EXTERNAL POTENTIAL`
is rejected with `CPMDC_ERROR_DECK_EXTERNAL_POTENTIAL` instead of a generic
evaluation failure. Python raises `CPMDCInvalidArgument`,
`CPMDCUnavailable`, `CPMDCUnsupportedPeriodicEmbedding`,
`CPMDCEvaluationFailed`, or `CPMDCDeckExternalPotential`.
