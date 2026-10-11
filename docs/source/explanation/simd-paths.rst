Paths
=====

``cpmdc`` has no hot numeric kernel. The C sources compile a deck and
read Cap'n Proto messages. The Fortran sources are the embedding bridge
into OpenCPMD. There is no AVX-512 path, no AVX2 path, and no NEON path
in this repository. Every translation unit compiles for x86-64 (v1)
on x86_64 and armv8-a on aarch64. The Meson build passes
``-march=x86-64`` or ``-march=armv8-a``. conda-forge compilers export
``CFLAGS=-march=nocona``, and Meson appends that environment after its
own arguments, so the later ``-march`` would win.
``scripts/with_baseline_flags.sh`` drops ``-march`` and the other ISA
selectors from ``CFLAGS``, ``CXXFLAGS``, ``FFLAGS``, and ``FCFLAGS``
before ``meson setup``. The Pixi test tasks source that script. The
build does not pass
``-march=native``, ``-mavx*``, or ``-march=x86-64-v3``.

OpenCPMD is Fortran. This wave does not add runtime dispatch inside it.
The portable configure files in ``tools/`` are the build target. See
:doc:`the archive how-to <../howto/opencpmd-archive>`.

How selection works
===================

There is no path table and no ``CPMDC_SIMD`` variable. ``cpmdc`` does
not read ``CPMDC_SIMD``. An environment value cannot select a kernel
this library does not ship. The scalar compiler default is the only
path, and it is the path the tests run.

OpenBLAS and FFTW, which an OpenCPMD link pulls in, select kernels at
run time when built with ``DYNAMIC_ARCH`` and with FFTW's SIMD
codelets. Those instructions sit in those libraries, not in a
``cpmdc`` symbol. The ISA audit allowlist names no ``cpmdc`` symbol.
It names OpenBLAS or libgfortran only when a build links that library
statically into an audited artefact. The default link is dynamic, so
the allowlist shipped here has no symbol rule and no vendored path.

Dispatched kernels elsewhere in this set of libraries follow the
cluster-pair layout of Páll and Hess, Comput. Phys. Commun. 184, 2641
(2013), doi:10.1016/j.cpc.2013.06.003. One binary per microarchitecture
belongs in EESSI (Dröge et al., Softw. Pract. Exp. 53, 176 (2023),
doi:10.1002/spe.3075).

Bit-identical results
=====================

The same input gives the same output bits on every baseline build,
because there is one code path. C and Fortran compile with
``-ffp-contract=off``. On MSVC the matching flag is ``/fp:precise``,
which does not contract. The compiler does not fuse a multiply and an
add, so a result does not depend on whether the build host would have
contracted those operations. Reductions in this repository do not
split a sum into partials. Sources do not call a fused multiply-add
intrinsic.

Two exceptions sit outside that guarantee. Cell angles in
``src/cpmdc_params.c`` call libm ``cos`` when a six-element cell is
written from degrees. That ``cos`` is libm, not an in-house polynomial,
and glibc can pick an implementation at run time, so those characters
of a rendered deck can differ across libm builds. OpenBLAS, FFTW, and
OpenCPMD's own math are not compiled by ``cpmdc``. Their kernels are
not bit-identical across CPUs. glibc ``exp``, ``log``, and ``pow``
have the same kind of runtime choice. ``cpmdc`` does not call them for
a value it defines.
