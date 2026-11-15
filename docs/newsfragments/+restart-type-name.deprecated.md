`include/cpmdc_restart.h` names its opaque image `CPMDCRestart`, as `cpmdc.h` names `CPMDCSession` and `CPMDCResult`, and its prototypes use that name.
`cpmdc_restart`, the 0.4.0 name, is deprecated; it stays a typedef of the same `struct cpmdc_restart`, so existing code compiles and links unchanged.
