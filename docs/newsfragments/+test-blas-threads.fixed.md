Meson tests now use a default setup that pins `OPENBLAS_NUM_THREADS`,
`OMP_NUM_THREADS`, and `MKL_NUM_THREADS` to 1 so parallel `meson test`
workers do not oversubscribe OpenBLAS. OpenCPMD live targets
`extended-surface`, `session-cached-geometry`, and `configure-overlay`
get 300 s timeouts when `with_cpmd` is enabled. Meson replaces those
three variables on test children; parent-shell exports do not apply (see
`docs/orgmode/reference/build-options.org`).
