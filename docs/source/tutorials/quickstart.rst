Stub And Parser Build
=====================

The default build is the fastest way to check the public ABI. It builds
the parser, generated Cap'n Proto readers, shared ``libcpmdc``, and a
deterministic reference evaluator. No OpenCPMD checkout is needed for
this path.

.. code:: bash

   meson setup build -Dwith_tests=true
   meson compile -C build
   meson test -C build --print-errorlogs

This configuration compiles the vendored ``capnp-c`` runtime, generates
C readers for ``schema/Potentials.capnp``, and runs cmocka suites:

-  ``stub`` — ``tests/test_stub_abi.c``
-  ``cmocka`` / ``capnp`` — ``tests/cmocka/test_params_render_cmocka.c``
   (deck render)
-  ``socket`` — ForceInput sizing and session ``PotentialResult`` buffer
   contract

cmocka is resolved through pkg-config. The parser tests use Cap'n Proto
text fixtures encoded by the ``capnp`` CLI (``tests/encode_capnp.py``).

The standalone stub target is different from the default shared engine:
``libcpmdc_stub.a`` only provides linkable symbols and reports
``cpmdc_available() == 0``. The default shared ``libcpmdc`` reports
availability through the reference evaluator and can run the session
tests.

End-to-End Suites
=================

The embed shell ships a deterministic harmonic evaluator so single-point
and multi-step session paths work without linking OpenCPMD archives.
Real PW-DFT evaluation replaces that evaluator when archives are wired.

.. code:: bash

   meson test -C build --suite e2e --print-errorlogs
   # or
   pixi run test-e2e

-  ``e2e-single-point`` — ``cpmdc_calculate_result`` +
   ``cpmdc_energy_forces`` on one ``ForceInput`` (HO-like geometry).
-  ``e2e-optimizer-session`` — one ``CPMDCSession``, three geometry
   steps (A → B → A), forces on B, topology rejection, then a second
   session for species change.

OpenCPMD Archive Build
======================

The default build emits ``libcpmdc.so`` for dlopen consumers such as
rgpot. Passing a tree path links that same C ABI surface against
``libcpmd.a``.

.. code:: bash

   export CPMD_ROOT=/path/to/OpenCPMD/CPMD
   meson setup build-cpmd \
     -Dwith_cpmd=true \
     -Dcpmd_root="$CPMD_ROOT" \
     -Dwith_tests=true
   meson compile -C build-cpmd

The OpenCPMD tree must contain ``lib/libcpmd.a`` and ``obj/timetag.o``
from a completed executable build. If the params use pseudopotential
library tokens such as ``O_MT_BLYP.psp``, set ``CPMDC_PSEUDO_DIR`` to
the directory containing those files before running real backend tests
or embedding the library.

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
