The GROMACS CI job sets `GMX_SIMD` to `AVX2_256`, includes the runner operating
system, architecture, and that setting in the cache key, and the build script
passes `GMX_SIMD` through as `-DGMX_SIMD`.
A restored binary that fails its version command is removed and built again.
