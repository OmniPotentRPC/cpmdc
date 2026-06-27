Stub And Parser Build
=====================

The default build is the fastest way to check the public ABI. It builds
the parser, generated Cap'n Proto readers, shared ``libcpmdc``, and a
deterministic reference evaluator. No OpenCPMD checkout is needed for
this path.

Install Meson, Ninja, Cap'n Proto, cmocka, C and Fortran compilers, and
pkg-config with your system package manager or use the checked-in Pixi
environment.

+---------------------+----------------------------------------------------------+-----------------------+
| Path                | Commands                                                 | Expected coverage     |
+=====================+==========================================================+=======================+
| Full default suite  | ``meson setup build -Dwith_tests=true``;                 | parser, ABI, feature  |
|                     | ``meson compile -C build``;                              | inventory, sessions,  |
|                     | ``meson test -C build --print-errorlogs``                | E2E reference         |
|                     |                                                          | evaluator             |
+---------------------+----------------------------------------------------------+-----------------------+
| cmocka/parser focus | ``meson test -C build --suite cmocka --print-errorlogs`` | Cap'n Proto decode,   |
|                     |                                                          | deck rendering,       |
|                     |                                                          | feature table checks  |
+---------------------+----------------------------------------------------------+-----------------------+
| session E2E focus   | ``meson test -C build --suite e2e --print-errorlogs``    | one-shot calls,       |
|                     |                                                          | session calls,        |
|                     |                                                          | topology rejection,   |
|                     |                                                          | unit conversion       |
+---------------------+----------------------------------------------------------+-----------------------+
| docs                | ``pixi run -e docs docbld``                              | Org export, generated |
|                     |                                                          | C API pages, Sphinx   |
|                     |                                                          | HTML                  |
+---------------------+----------------------------------------------------------+-----------------------+

.. code:: bash

   meson setup build -Dwith_tests=true
   meson compile -C build
   meson test -C build --print-errorlogs

This configuration compiles the vendored ``capnp-c`` runtime, generates
C readers for ``schema/Potentials.capnp``, and runs cmocka suites:

- ``stub`` — ``tests/test_stub_abi.c``
- ``abi`` — exported header and feature-discovery surface
- ``cmocka`` / ``capnp`` — structured ``CPMDParams`` deck rendering
- ``params`` / ``inventory`` — schema field coverage and CPMD catalog
  feature IDs
- ``socket`` — ``ForceInput`` sizing and session ``PotentialResult``
  buffer contract
- ``e2e`` — one-shot and session calls through the deterministic
  reference evaluator

cmocka is resolved through pkg-config. The parser tests use Cap'n Proto
text fixtures encoded by the ``capnp`` CLI (``tests/encode_capnp.py``).

The standalone stub target is different from the default shared engine:
``libcpmdc_stub.a`` only provides linkable symbols and reports
``cpmdc_available() == 0``. The default shared ``libcpmdc`` reports
availability through the reference evaluator and can run the session
tests.

Run The Host Example
====================

``examples/host_step.c`` is a minimal C host for the public ABI. It
reads one serialized ``CPMDParams`` file and one serialized
``ForceInput`` file, creates a ``CPMDCSession``, sizes the
``PotentialResult`` output, evaluates one step, and prints the result
summary.

.. code:: bash

   meson test -C build example-host-step --print-errorlogs

The example output has this shape:

.. code:: text

   energy_h=...
   potential_result_size_bytes=...
   message=...

The Meson test feeds the example from the same generated fixture
binaries used by the E2E suites. Those fixture binaries are produced
from ``tests/cpmd_params_parser.capnp.txt`` or
``tests/water_blyp_params.capnp.txt`` and
``tests/force_input_step_a.capnp.txt``, depending on whether the build
links OpenCPMD archives.

Test Selection
==============

Use the smallest suite that proves the layer you changed:

+----------------------+-------------------------------------------------------------------------------------+-----------------------+
| Change               | Command                                                                             | Pass condition        |
+======================+=====================================================================================+=======================+
| Parser, feature      | ``meson test -C build --print-errorlogs``                                           | all configured tests  |
| inventory, or C ABI  |                                                                                     | pass                  |
+----------------------+-------------------------------------------------------------------------------------+-----------------------+
| Deck rendering or    | ``meson test -C build --suite cmocka --print-errorlogs``                            | cmocka render suites  |
| Cap'n Proto fixtures |                                                                                     | pass                  |
+----------------------+-------------------------------------------------------------------------------------+-----------------------+
| Session lifecycle,   | ``meson test -C build --suite e2e --print-errorlogs``                               | single-point and      |
| topology, or units   |                                                                                     | optimizer-session     |
|                      |                                                                                     | tests pass            |
+----------------------+-------------------------------------------------------------------------------------+-----------------------+
| Public documentation | ``pixi run -e docs sphinxbld``                                                      | Sphinx build succeeds |
+----------------------+-------------------------------------------------------------------------------------+-----------------------+
| Org-to-RST refresh   | ``pixi run -e docs mkrst``                                                          | checked-in            |
|                      |                                                                                     | ``docs/source/*.rst`` |
|                      |                                                                                     | matches the org       |
|                      |                                                                                     | sources               |
+----------------------+-------------------------------------------------------------------------------------+-----------------------+
| OpenCPMD archive     | ``CPMDC_PSEUDO_DIR=/path/to/PP_LIBRARY meson test -C build-cpmd --print-errorlogs`` | default tests plus    |
| link                 |                                                                                     | live water parity     |
|                      |                                                                                     | pass                  |
+----------------------+-------------------------------------------------------------------------------------+-----------------------+

``meson test --print-errorlogs`` is preferred because cmocka assertion
output and OpenCPMD diagnostics stay attached to the failed test.

End-to-End Suites
=================

The embed shell ships a deterministic harmonic evaluator so single-point
and multi-step session paths work without linking OpenCPMD archives.
Real PW-DFT evaluation replaces that evaluator when archives are wired.

.. code:: bash

   meson test -C build --suite e2e --print-errorlogs
   # or
   pixi run test-e2e

- ``e2e-single-point`` — ``cpmdc_calculate_result`` +
  ``cpmdc_energy_forces`` on one ``ForceInput`` (HO-like geometry).
- ``e2e-optimizer-session`` — one ``CPMDCSession``, three geometry steps
  (A → B → A), forces on B, topology rejection, then a second session
  for species change.

OpenCPMD Archive Build
======================

The default build emits ``libcpmdc.so`` for dlopen consumers such as
rgpot. Passing a tree path links that same C ABI surface against
``libcpmd.a``.

.. code:: bash

   export CPMD_ROOT=/path/to/OpenCPMD/CPMD
   export CPMDC_PSEUDO_DIR=/path/to/CPMD-Regtests/tests/PP_LIBRARY
   meson setup build-cpmd \
     -Dwith_cpmd=true \
     -Dcpmd_root="$CPMD_ROOT" \
     -Dwith_tests=true
   meson compile -C build-cpmd
   meson test -C build-cpmd --print-errorlogs

The OpenCPMD tree must contain ``lib/libcpmd.a`` and ``obj/timetag.o``
from a completed executable build. If the params use pseudopotential
library tokens such as ``O_MT_BLYP.psp``, set ``CPMDC_PSEUDO_DIR`` to
the directory containing those files before running real backend tests
or embedding the library.

The live build runs the same public C ABI as the default build. When
``with_cpmd=true`` is enabled, water Cap'n Proto parity tests are added
and the long E2E tests exercise the linked OpenCPMD archive.

Common Failures
===============

These checks usually identify environment problems before code changes
are needed:

+------------------------------------------+---------------------------------------------------+
| Symptom                                  | Check                                             |
+==========================================+===================================================+
| ``with_cpmd requires ... lib/libcpmd.a`` | ``CPMD_ROOT`` points at a completed OpenCPMD      |
|                                          | build tree with ``lib/libcpmd.a``                 |
+------------------------------------------+---------------------------------------------------+
| OpenCPMD cannot find ``O_MT_BLYP.psp``   | ``CPMDC_PSEUDO_DIR`` points at the regtest        |
| or ``H_CVB_BLYP.psp``                    | ``PP_LIBRARY`` or the params use absolute         |
|                                          | pseudopotential paths                             |
+------------------------------------------+---------------------------------------------------+
| A session reports a topology change      | the same ``CPMDCSession`` received a different    |
|                                          | atom count or ordered atomic-number list; create  |
|                                          | a new session                                     |
+------------------------------------------+---------------------------------------------------+
| ``PotentialResult buffer too small``     | call                                              |
|                                          | ``cpmdc_potential_result_size_for_force_input()`` |
|                                          | and allocate at least that many bytes             |
+------------------------------------------+---------------------------------------------------+
| Unit parsing fails                       | use supported length/energy strings such as       |
|                                          | ``angstrom``, ``bohr``, ``eV``, and ``hartree``   |
+------------------------------------------+---------------------------------------------------+
| ``capnp encode`` fails on a fixture      | re-run the command against                        |
|                                          | ``schema/Potentials.capnp`` and the intended root |
|                                          | type, either ``CPMDParams`` or ``ForceInput``     |
+------------------------------------------+---------------------------------------------------+

Cap'n Proto Fixtures
====================

Encode a ``CPMDParams`` or ``ForceInput`` text message the same way
tests do:

.. code:: bash

   capnp encode schema/Potentials.capnp CPMDParams \
     < tests/cpmd_params_parser.capnp.txt > /tmp/cpmd_params.bin
   capnp encode schema/Potentials.capnp ForceInput \
     < tests/force_input_step_a.capnp.txt > /tmp/force_input.bin

Pass those buffers to ``cpmdc_session_create()`` and
``cpmdc_session_calculate_result()`` from C, Python (``pycapnp``), or
``rgpot``.
