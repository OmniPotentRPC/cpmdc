Development Checks
==================

Run the default parser, socket, and stub checks before changing the ABI,
schema, renderer, or session runtime:

.. code:: bash

   meson setup build -Dwith_tests=true
   meson compile -C build
   meson test -C build --print-errorlogs

Use focused suites while iterating, then run the default suite before
pushing a public change:

+----------------------------------+-------------------------------------------------------------------------------------+
| Contract                         | Command                                                                             |
+==================================+=====================================================================================+
| Full default ABI and parser      | ``meson test -C build --print-errorlogs``                                           |
| suite                            |                                                                                     |
+----------------------------------+-------------------------------------------------------------------------------------+
| Cap'n Proto decode and deck      | ``meson test -C build --suite cmocka --print-errorlogs``                            |
| rendering                        |                                                                                     |
+----------------------------------+-------------------------------------------------------------------------------------+
| Session, topology, result        | ``meson test -C build --suite e2e --print-errorlogs``                               |
| sizing, and unit conversion      |                                                                                     |
+----------------------------------+-------------------------------------------------------------------------------------+
| Generated feature inventory      | ``meson test -C build feature-inventory --print-errorlogs``                         |
+----------------------------------+-------------------------------------------------------------------------------------+
| OpenCPMD archive build           | ``CPMDC_PSEUDO_DIR=/path/to/PP_LIBRARY meson test -C build-cpmd --print-errorlogs`` |
+----------------------------------+-------------------------------------------------------------------------------------+

Params render coverage lives under ``tests/cmocka/`` and uses cmocka
assertions on generated Cap'n Proto types (``read_CPMDParams``, deck
substrings) only. Session tests use encoded ``CPMDParams`` and
``ForceInput`` fixtures from ``tests/*.capnp.txt``.

Source Layout
=============

+------------------------------+---------------------------------------+
| Path                         | Ownership                             |
+==============================+=======================================+
| ``schema/Potentials.capnp``  | shared RPC and direct-call wire       |
|                              | contract                              |
+------------------------------+---------------------------------------+
| ``include/cpmdc.h``          | public C ABI and session result entry |
|                              | points                                |
+------------------------------+---------------------------------------+
| ``include/cpmdc_features.h`` | runtime feature discovery ABI         |
+------------------------------+---------------------------------------+
| ``src/cpmdc_params.c``       | Cap'n Proto decode, effective config, |
|                              | CPMD deck rendering                   |
+------------------------------+---------------------------------------+
| ``src/cpmdc.c``              | session lifecycle, topology guard,    |
|                              | unit conversion, result writing       |
+------------------------------+---------------------------------------+
| ``src/cpmd_embed_c_api.F90`` | Fortran ``iso_c_binding`` bridge and  |
|                              | OpenCPMD archive path                 |
+------------------------------+---------------------------------------+
| ``tests/``                   | encoded fixtures and C/cmocka         |
|                              | contract tests                        |
+------------------------------+---------------------------------------+

Documentation
=============

Documentation source lives under ``docs/orgmode`` (authoritative prose),
with light/dark logos in ``docs/source/_static/`` and
``branding/logo/``. Build docs (exports Org → RST, Doxygen XML, Sphinx
HTML with Antics analytics):

.. code:: bash

   pixi run -e docs docbld

Regenerate RST after editing Org sources:

.. code:: bash

   pixi run -e docs mkrst

The Sphinx theme injects ``antics-api.turtletech.us/antics.js`` via
``docs/source/_templates/partials/extra-head.html`` (same pattern as
``nwchemc`` / ``rgpot``).

Release Process
===============

.. toctree::
   :maxdepth: 2

   developer/release
