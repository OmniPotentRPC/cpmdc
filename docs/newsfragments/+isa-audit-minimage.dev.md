The instruction set architecture (ISA) allowlist names the minimage kernels built for Advanced Vector Extensions (AVX), AVX2 or AVX-512.
minimage calls each kernel only after the Rust standard library reports its features at run time.
The allowlist also names two functions of the Rust standard library.
One is the `memchr` AVX2 searcher, which that build links and never calls.
The other is `_xgetbv`, which runs only after the processor identification instruction (CPUID) reports that the processor and the operating system support extended state saving.
The portable ISA test rejects an allow rule that is not anchored, uses a wildcard, matches a cpmdc symbol, or states no run-time check.
