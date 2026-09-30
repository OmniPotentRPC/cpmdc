Goal
====

Produce an OpenCPMD build tree whose ``lib/libcpmd.a`` can be linked
into the shared ``libcpmdc.so`` and evaluated repeatedly from one host
process. A stock OpenCPMD archive does not work inside a host:
``libcpmdc`` imports routines that only the patches add, the shared link
stops on code built without ``-fPIC``, and the unpatched SCF driver
drops the forces and the orbitals that ``cpmdc`` reads after each call.
Seven patches in ``tools/`` and one compiler flag fix those.
``opencpmd_c_mem_addrs.patch`` makes ``cGetMemAddrs`` return the address
as ``size_t``. GCC 14 rejects the stock return of a pointer from that
function.

Patches and what each one is for
================================

+--------------------------------------+--------------------------------+-----------------------------+
| Patch                                | File patched                   | Without it                  |
+======================================+================================+=============================+
| ``opencpmd_embed_geometry.patch``    | ``src/embed_ctrl.mod.F90``,    | ``initrun`` still writes    |
|                                      | ``src/initrun_driver.mod.F90``,| ``GEOMETRY``, and           |
|                                      | ``src/SOURCES``                | ``rwfopt`` has no           |
|                                      |                                | ``embed_ctrl`` to read      |
+--------------------------------------+--------------------------------+-----------------------------+
| ``opencpmd_embed_rwfopt.patch``      | ``src/rwfopt_utils.mod.F90``   | ``rwfopt`` drops ``fion``   |
|                                      |                                | and starts every SCF from   |
|                                      |                                | a fresh guess, and the      |
|                                      |                                | embed bridge does not       |
|                                      |                                | compile                     |
+--------------------------------------+--------------------------------+-----------------------------+
| ``opencpmd_converged_state.patch``   | ``src/updwf_utils.mod.F90``    | the optimiser applies one   |
|                                      |                                | more update after           |
|                                      |                                | convergence, so the saved   |
|                                      |                                | orbitals no longer match    |
|                                      |                                | the reported energy and     |
|                                      |                                | forces; steepest descent    |
|                                      |                                | with ``iproj <= 1`` needs   |
|                                      |                                | its ``gemax`` reset before  |
|                                      |                                | the convergence check       |
+--------------------------------------+--------------------------------+-----------------------------+
| ``opencpmd_kpoints_inputfile.patch`` | ``src/rkpnt_utils.mod.F90``    | the k-point report reads    |
|                                      |                                | the host's ``argv[1]`` as   |
|                                      |                                | the input file name and     |
|                                      |                                | stops with ``STOP 12345``   |
|                                      |                                | when it is missing or       |
|                                      |                                | longer than 80 characters   |
+--------------------------------------+--------------------------------+-----------------------------+
| ``opencpmd_stopgm_return.patch``     | ``src/error_handling.mod.F90`` | a CPMD stop inside an       |
|                                      |                                | evaluation takes the        |
|                                      |                                | ``my_stopall`` path, which  |
|                                      |                                | OpenCPMD writes to end the  |
|                                      |                                | run                         |
+--------------------------------------+--------------------------------+-----------------------------+
| ``opencpmd_tistopgm.patch``          | ``src/timer.mod.F90``          | ``tistopgm`` indexes        |
|                                      |                                | ``trace_depth`` before      |
|                                      |                                | ``tistart``. The depth      |
|                                      |                                | is still ``HUGE(0)``,       |
|                                      |                                | so the process faults       |
|                                      |                                | before ``LocalError``       |
+--------------------------------------+--------------------------------+-----------------------------+
| ``opencpmd_c_mem_addrs.patch``       | ``src/c_mem_utils.c``          | ``cGetMemAddrs`` returns a  |
|                                      |                                | ``size_t`` pointer. GCC 14  |
|                                      |                                | stops on that conversion    |
+--------------------------------------+--------------------------------+-----------------------------+

The rwfopt patch publishes ``embed_set_warm_orbitals``,
``embed_set_need_forces``, and ``embed_reset_warm_orbitals``.
``embed_warm_orbitals`` and ``embed_need_forces`` default to false, which
is the ``cpmd.x`` behaviour. ``rwfopt`` leaves ``fion`` allocated only when
``embed_need_forces`` is true, and frees a surviving ``fion`` before its
own ``ALLOCATE``. It saves ``c0`` only when ``embed_warm_orbitals`` is
true and the SCF converged. An unconverged pass leaves the previous copy.
A failed store allocation calls ``stopgm``. On the next call it restores
that copy after ``initrun`` when the shape still matches, and a different
shape calls ``stopgm``. It sets ``tfor`` when ``embed_need_forces`` is
true. The stop patch publishes ``cpmd_stopgm_hook``, a C function
pointer. ``stopgm`` calls it with the stop code passed by value, after
writing ``LocalError-*.log``. A nonzero return makes ``stopgm`` return to
its caller, and a null pointer makes ``stopgm`` call ``my_stopall``. The
caller then continues past the failed check. The wavefunction, ``fion``,
and module state of that call are undefined, so ``cpmdc`` discards the
result and sets CPMD up again before the next call. ``stopgm`` runs only
on the ranks that hit the error. A host with more than one rank aborts
the others, because a rank that returns while the rest wait in an MPI
call leaves them waiting. A standalone ``cpmd.x`` leaves the
pointer null and links without ``libcpmdc``. While an evaluation is armed,
``libcpmdc`` stores a catch in the pointer. The catch records the code and
returns 1, and ``cpmdc`` reports the stop as a failed call (see
:doc:`troubleshooting <troubleshooting>`). ``cpmdc_embed_catch`` and
``cpmdc_note_stop`` are exported from ``libcpmdc``. The calculator split
does not need a patch: ``mp_start`` assigns ``mp_comm_world`` only when
CPMD itself calls ``MPI_Init``, and ``cpmdc_embed_bind_calculator``
initialises MPI and stores its communicator before that. A module flag,
``embed_calculator_bound``, makes a second bind return the same index.
The timer patch changes ``tistopgm``. ``trace_depth`` stays ``HUGE(0)``
until ``tistart``, and without the guard the call-stack walk indexes
``trace_names`` from that value, so a ``stopgm`` before ``tistart`` faults
instead of writing ``LocalError``. The walk runs only when the depth is
inside ``SIZE(tname%trace_names)``.

``tests/test_opencpmd_patch_integrity.py`` checks that the patches are
portable unified diffs against the files named above. Apply
``opencpmd_embed_geometry.patch`` first: ``rwfopt`` reads
``embed_ctrl``. The other five are the commits on
`OpenCPMD pull request 9 <https://github.com/OpenCPMD/CPMD/pull/9>`__
(branch ``embedding-hooks``, ``43cf4e7``) and apply in sequence to
OpenCPMD commit ``062582b``. Apply ``opencpmd_c_mem_addrs.patch`` with
them. ``cuda_get_address`` stores the ``size_t`` that function returns.

Build a patched tree
====================

Start from an OpenCPMD source checkout and apply the patches before the
first build:

.. code:: bash

   git clone https://github.com/OpenCPMD/CPMD.git opencpmd
   cd opencpmd
   for p in embed_geometry embed_rwfopt converged_state \
            kpoints_inputfile stopgm_return tistopgm c_mem_addrs; do
     patch -p1 < /path/to/cpmdc/tools/opencpmd_$p.patch
   done

The OpenCPMD sample ``LINUX-X86_64-GFORTRAN-MPI`` sets ``LIBS`` to one
person's library path. Copying that file and adding ``-fPIC`` leaves
the path in the build. ``cpmdc`` ships
``tools/LINUX-X86_64-GFORTRAN-MPI-PIC`` for the command below. It calls
``mpif90`` and ``gcc`` from ``PATH``, adds ``-fPIC`` to ``FFLAGS`` and
``CFLAGS``, and sets ``LIBS`` to ``-lopenblas``. When ``mpif90`` accepts
it, ``FFLAGS`` also gains ``-fallow-argument-mismatch``.
``mp_bcast_byte`` broadcasts the bytes of its first argument, and one
procedure calls it with several derived types. gfortran 10 and later
stop on those calls without the flag. That OpenBLAS build has to
export the LAPACK symbols. ``libcpmdc`` is a shared library; a
member of ``libcpmd.a`` compiled without ``-fPIC`` fails the link with
a ``R_X86_64_PC32`` relocation error. Copy the shipped file into the
clone, then configure and build:

.. code:: bash

   cp /path/to/cpmdc/tools/LINUX-X86_64-GFORTRAN-MPI-PIC configure/
   ./configure.sh -DEST=/path/to/opencpmd-build LINUX-X86_64-GFORTRAN-MPI-PIC
   make -C /path/to/opencpmd-build -j 8
   ls /path/to/opencpmd-build/lib/libcpmd.a /path/to/opencpmd-build/obj/timetag.o

``/path/to/opencpmd-build`` is the ``cpmd_root`` passed to Meson. The
``cpmdc`` build reads ``lib/libcpmd.a``, the module files in ``obj/``,
and the headers in ``src/``.

Patch a tree that is already built
==================================

When the tree was built before the patches were applied, patch the
source and rebuild only the affected members:

.. code:: bash

   cd /path/to/opencpmd-source
   for p in embed_geometry embed_rwfopt converged_state \
            kpoints_inputfile stopgm_return tistopgm c_mem_addrs; do
     patch -p1 < /path/to/cpmdc/tools/opencpmd_$p.patch
   done
   /path/to/cpmdc/tools/rebuild_opencpmd_embed.sh /path/to/opencpmd-build

``rebuild_opencpmd_embed.sh`` recompiles ``error_handling.mod.o``,
``scex_utils.mod.o``, ``rwfopt_utils.mod.o``,
``updwf_utils.mod.o``, ``rkpnt_utils.mod.o``, and ``timer.mod.o`` with one Make job, then
replaces those members with ``ar r`` and runs ``ranlib``. It runs one
job because gfortran stores a copy of ``Scex_t`` inside
``rwfopt_utils.mod``: a parallel rebuild can leave the two module files
naming different components. It calls ``ar`` directly because the
generated Makefile sets ``AR`` without an operation letter, so its
archive rule does not replace members. ``MAKE`` and ``AR`` in the
environment override the programs it runs.

Link and check
==============

OpenCPMD leaves ``timetag`` out of ``libcpmd.a`` and links
``obj/timetag.o`` into ``cpmd.x``. ``header`` calls ``timetag``, so the
embed link needs that object too. Build it from the configured tree
before ``meson setup``:

.. code:: bash

   make -C /path/to/opencpmd-build/obj -f /path/to/opencpmd-build/Makefile timetag.o

.. code:: bash

   meson setup build-cpmd -Dwith_cpmd=true -Dcpmd_root=/path/to/opencpmd-build
   meson compile -C build-cpmd
   nm -D build-cpmd/libcpmdc.so | grep -c ' T cpmdc_'

The last command counts the exported ``cpmdc_*`` symbols. The first
evaluation needs the pseudopotential directory in ``CPMDC_PSEUDO_DIR``;
the :doc:`first energy tutorial <../tutorials/first-energy>` runs one.
Stress tensors, which ``PotentialResult.stress`` carries for periodic
cells, need no extra patch: ``cpmdc`` sets ``cntl%tpres`` before the SCF
and copies ``paiu/omega`` afterwards. An isolated cell does not set the
snapshot. The box has a volume, and CPMD does not compute the tensor.
