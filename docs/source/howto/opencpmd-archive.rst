Goal
====

Produce an OpenCPMD build tree whose ``lib/libcpmd.a`` can be linked
into the shared ``libcpmdc.so`` and evaluated repeatedly from one host
process. A stock OpenCPMD archive does not work inside a host:
``libcpmdc`` imports routines that only the patches add, the shared link
stops on code built without ``-fPIC``, and the unpatched SCF driver
drops the forces and the orbitals that ``cpmdc`` reads after each call.
Four patches in ``tools/`` and one compiler flag fix those.

Patches and what each one is for
================================

+--------------------------------------+--------------------------------+-----------------------------+
| Patch                                | File patched                   | Without it                  |
+======================================+================================+=============================+
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

The rwfopt patch publishes ``embed_set_warm_orbitals``,
``embed_set_need_forces``, and ``embed_reset_warm_orbitals``.
``embed_warm_orbitals`` and ``embed_need_forces`` default to false, which
is the ``cpmd.x`` behaviour. ``rwfopt`` leaves ``fion`` allocated only when
``embed_need_forces`` is true, and frees a surviving ``fion`` before its
own ``ALLOCATE``. It saves ``c0`` only when ``embed_warm_orbitals`` is
true. A failed store allocation calls ``stopgm``. On the next call it
restores that copy after ``initrun`` when the dimensions still match, and
sets ``tfor`` when ``embed_need_forces`` is true. The stop patch
publishes ``cpmd_stopgm_hook``, a C function pointer. ``stopgm`` calls it
with the stop code passed by value, after writing ``LocalError-*.log``.
A nonzero return makes ``stopgm`` return to its caller, and a null pointer
makes ``stopgm`` call ``my_stopall``. The caller then continues past the
failed check, so ``cpmdc`` treats every result of that call as invalid and
sets CPMD up again before the next one. A standalone ``cpmd.x`` leaves the
pointer null and links without ``libcpmdc``. While an evaluation is armed,
``libcpmdc`` stores a catch in the pointer. The catch records the code and
returns 1, and ``cpmdc`` reports the stop as a failed call (see
:doc:`troubleshooting <troubleshooting>`). ``cpmdc_embed_catch`` and
``cpmdc_note_stop`` are exported from ``libcpmdc``. The calculator split
does not need a patch: ``mp_start`` assigns ``mp_comm_world`` only when
CPMD itself calls ``MPI_Init``, and ``cpmdc_embed_bind_calculator``
initialises MPI and stores its communicator before that. A module flag,
``embed_calculator_bound``, makes a second bind return the same index.

``tests/test_opencpmd_patch_integrity.py`` checks that the patches are
portable unified diffs against the files named above. All four apply in
sequence to OpenCPMD commit ``062582b``.

Build a patched tree
====================

Start from an OpenCPMD source checkout and apply the patches before the
first build:

.. code:: bash

   git clone https://github.com/OpenCPMD/CPMD.git opencpmd
   cd opencpmd
   for p in embed_rwfopt converged_state \
            kpoints_inputfile stopgm_return; do
     patch -p1 < /path/to/cpmdc/tools/opencpmd_$p.patch
   done

Pick a configuration from ``./configure.sh -help`` that matches your
compilers, for example ``LINUX-X86_64-GFORTRAN-MPI``. Copy it under a
new name in ``configure/`` and add ``-fPIC`` to both ``FFLAGS`` and
``CFLAGS``. ``libcpmdc`` is a shared library; a member of ``libcpmd.a``
compiled without ``-fPIC`` fails the link with a ``R_X86_64_PC32``
relocation error. Then configure into a separate build directory and
build:

.. code:: bash

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
   for p in embed_rwfopt converged_state \
            kpoints_inputfile stopgm_return; do
     patch -p1 < /path/to/cpmdc/tools/opencpmd_$p.patch
   done
   /path/to/cpmdc/tools/rebuild_opencpmd_embed.sh /path/to/opencpmd-build

``rebuild_opencpmd_embed.sh`` recompiles ``error_handling.mod.o``,
``scex_utils.mod.o``, ``rwfopt_utils.mod.o``,
``updwf_utils.mod.o``, and ``rkpnt_utils.mod.o`` with one Make job, then
replaces those members with ``ar r`` and runs ``ranlib``. It runs one
job because gfortran stores a copy of ``Scex_t`` inside
``rwfopt_utils.mod``: a parallel rebuild can leave the two module files
naming different components. It calls ``ar`` directly because the
generated Makefile sets ``AR`` without an operation letter, so its
archive rule does not replace members. ``MAKE`` and ``AR`` in the
environment override the programs it runs.

Link and check
==============

.. code:: bash

   meson setup build-cpmd -Dwith_cpmd=true -Dcpmd_root=/path/to/opencpmd-build
   meson compile -C build-cpmd
   nm -D build-cpmd/libcpmdc.so | grep -c ' T cpmdc_'

The last command counts the exported ``cpmdc_*`` symbols. The first
evaluation needs the pseudopotential directory in ``CPMDC_PSEUDO_DIR``;
the :doc:`first energy tutorial <../tutorials/first-energy>` runs one.
Stress tensors, which ``PotentialResult.stress`` carries for periodic
cells, need no extra patch: ``cpmdc`` sets ``cntl%tpres`` before the SCF
and copies ``paiu/omega`` afterwards.
