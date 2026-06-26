Stub And Parser Build
=====================

The default build exercises the stable C ABI and the Cap'n Proto parser
without linking a local OpenCPMD tree.

.. code:: bash

   meson setup build -Dwith_tests=true
   meson compile -C build
   meson test -C build --print-errorlogs

This configuration builds the stub library, the ISO\ :sub:`C` embed
shell (evaluation unavailable until OpenCPMD is wired), compiles the
vendored ``capnp-c`` runtime, generates C readers for
``schema/Potentials.capnp``, and runs cmocka suites:

-  ``stub`` — ``tests/test_stub_abi.c``
-  ``cmocka`` / ``capnp`` — ``tests/cmocka/test_params_render_cmocka.c``
   (deck render)
-  ``socket`` — ForceInput sizing and session ``PotentialResult`` buffer
   contract

cmocka is resolved through pkg-config. The parser test uses a Cap'n
Proto text fixture encoded by the ``capnp`` CLI
(``tests/encode_capnp.py``).

End-to-End Suites
=================

The embed shell ships a deterministic **reference PEF** (harmonic,
Z-scaled) so single-point and multi-step session paths work without
linking OpenCPMD archives. Real PW-DFT evaluation replaces that PEF when
archives are wired.

.. code:: bash

   meson test -C build --suite e2e --print-errorlogs
   # or
   pixi run test-e2e

-  ``e2e-single-point`` — ``cpmdc_calculate_result`` +
   ``cpmdc_energy_forces`` on one ``ForceInput`` (HO-like geometry).
-  ``e2e-optimizer-session`` — one ``CPMDCSession``, three geometry
   steps (A → B → A), forces on B, topology rejection, then a second
   session for species change.

OpenCPMD Embed Build (0.1 shell)
================================

The v0.1 shared library path installs the C ABI and Fortran
ISO\ :sub:`C` symbols when a tree path is provided. Full archive linking
and energy/force evaluation are the next embed milestone;
``cpmdc_available()`` remains ``0`` until that lands.

.. code:: bash

   export CPMD_ROOT=/path/to/OpenCPMD/CPMD
   meson setup build-cpmd \
     -Dwith_cpmd=true \
     -Dcpmd_root="$CPMD_ROOT" \
     -Dwith_tests=true
   meson compile -C build-cpmd

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
