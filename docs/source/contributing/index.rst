Development Checks
==================

Run the default parser, socket, and stub checks before changing the ABI
or schema:

.. code:: bash

   meson setup build -Dwith_tests=true
   meson compile -C build
   meson test -C build --print-errorlogs

Params render coverage lives under ``tests/cmocka/`` and uses cmocka
assertions on generated Cap'n Proto types (``read_CPMDParams``, deck
substrings) only.

Documentation
=============

Documentation source lives under ``docs/orgmode`` (authoritative prose),
with light/dark logos in ``docs/source/_static/`` and
``branding/logo/``. Build docs (exports Org → RST, Doxygen XML, Sphinx
HTML with Antics analytics):

.. code:: bash

   pixi run -e docs docbld

The Sphinx theme injects ``antics-api.turtletech.us/antics.js`` via
``docs/source/_templates/partials/extra-head.html`` (same pattern as
``nwchemc`` / ``rgpot``).

Release Process
===============

.. toctree::
   :maxdepth: 2

   developer/release
