The instruction set architecture (ISA) allowlist names the minimage kernels built for Advanced Vector Extensions (AVX), AVX2 or AVX-512.
minimage calls each kernel only after the Rust standard library reports its features at run time.
The allowlist also names three functions that minimage links but does not define.
One rule matches `find_impl` and `find_prefilter_impl` of the AVX2 packed-pair searcher in the `memchr` crate, which the standard library depends on through the `object` crate.
The standard library build links both functions and never calls them.
The other rule names `_xgetbv` from the Rust core library.
It runs only after the processor identification instruction (CPUID) reports that the processor and the operating system support extended state saving.
The portable ISA test rejects an allow rule that is not anchored, uses a wildcard, matches a cpmdc symbol, or states no run-time check.
