`CPMDAtomsPseudopotential.raggio` and `nonlinearCore` reach CPMD: cpmdc wrote
`RAGGIO=` and `NLCC` on the `LMAX` line, where CPMD reads neither, so the ionic
charge kept the default Gaussian width and the core charge was dropped.
They now follow the path on the `*file` line, as `KLEINMAN-BYLANDER` does, in
the deck cpmdc renders and in the `&ATOMS` block that the OpenCPMD path
rebuilds from each geometry, which also keeps the other options of a deck's own
`*file` line.
A path that ends in `.upf`, in any case, gains `UPF`, so that CPMD reads the
file as UPF; the rebuilt block used to drop that option, and OpenCPMD stopped
with `RECPNEW| &ATOM SECTION NOT FOUND` and a segmentation fault.
A deck line or `*file` option that does not fit its buffer fails the call
instead of dropping the tail of the line.
