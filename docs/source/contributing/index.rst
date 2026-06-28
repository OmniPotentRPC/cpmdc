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
| Base ``&CPMD`` keyword inventory | ``meson test -C build cpmd-base-keyword-inventory --print-errorlogs``               |
+----------------------------------+-------------------------------------------------------------------------------------+
| Top-level ``CPMDParams`` feature | ``meson test -C build cpmd-params-field-inventory --print-errorlogs``               |
| rows and duplicate inventory     |                                                                                     |
| lists                            |                                                                                     |
+----------------------------------+-------------------------------------------------------------------------------------+
| ``CPMDCpmdSection`` render       | ``meson test -C build cpmd-schema-render-coverage --print-errorlogs``               |
| mappings                         |                                                                                     |
+----------------------------------+-------------------------------------------------------------------------------------+
| Inline option token coverage     | ``meson test -C build cpmd-option-token-coverage --print-errorlogs``                |
+----------------------------------+-------------------------------------------------------------------------------------+
| Typed render field coverage      | ``meson test -C build cpmd-typed-render-field-coverage --print-errorlogs``          |
+----------------------------------+-------------------------------------------------------------------------------------+
| OpenCPMD archive build           | ``CPMDC_PSEUDO_DIR=/path/to/PP_LIBRARY meson test -C build-cpmd --print-errorlogs`` |
+----------------------------------+-------------------------------------------------------------------------------------+

Params render coverage lives under ``tests/cmocka/`` and uses cmocka
assertions on generated Cap'n Proto types (``read_CPMDParams``, deck
substrings) only. Session tests use encoded ``CPMDParams`` and
``ForceInput`` fixtures from ``tests/*.capnp.txt``.
``feature-inventory`` ties the generated JSON inventory to the schema,
public feature table, docs, public headers, ABI symbol list, and
optional ``CPMD_ROOT`` parser section probe.
``cpmd-base-keyword-inventory`` requires every base ``&CPMD`` keyword in
``schema/inventory/cpmd_cp_keywords.txt`` to resolve to a
``catalog.cpmd.*`` row. ``cpmd-params-field-inventory`` requires every
top-level ``CPMDParams`` schema field to resolve to a ``params.*``
feature row and rejects duplicate inventory IDs/lists before set
comparisons can hide them. ``cpmd-schema-render-coverage`` and
``cpmd-option-token-coverage`` keep typed ``CPMDCpmdSection`` fields and
fixture inline tokens tied to render coverage.
``cpmd-typed-render-field-coverage`` requires typed ``cpmd``,
``system``, ``dft``, and ``atoms`` fields to appear in render fixtures
or render assertions.

Add a typed CPMD keyword
========================

Use a typed field when a CPMD control has stable structure, affects
ordinary host inputs, or should be discoverable through
``cpmdc_feature_find()``. Use ``directives``, ``set``, ``generic``,
``raw``, or ``inputBlocks`` for rare literal deck fragments that should
stay text.

Start from the test surface:

+---------------------------------------+----------------------------------+
| File                                  | Required change                  |
+=======================================+==================================+
| ``tests/test_params_cp_dft_render.c`` | Add catalog coverage and         |
|                                       | rendered-token checks            |
+---------------------------------------+----------------------------------+
| ``tests/test_features_cmocka.c``      | Add ``catalog.cpmd.*`` and       |
|                                       | ``params.inputSections.cpmd.*``  |
|                                       | expectations                     |
+---------------------------------------+----------------------------------+
| ``tests/params_*.capnp.txt``          | Add a fixture value that         |
|                                       | exercises the new field          |
+---------------------------------------+----------------------------------+

Run the focused failure before adding production code:

.. code:: bash

   CPMD_ROOT=/tmp/OpenCPMD-cpmdc meson test -C build \
     params-cp-dft-render features-cmocka feature-inventory \
     cpmd-base-keyword-inventory cpmd-params-field-inventory \
     cpmd-schema-render-coverage cpmd-option-token-coverage \
     cpmd-typed-render-field-coverage \
     --print-errorlogs

Then update the production and docs surface in one change:

+---------------------------------------------+----------------------------------+
| File                                        | Required change                  |
+=============================================+==================================+
| ``schema/Potentials.capnp``                 | Add the ``CPMDCpmdSection``      |
|                                             | field with the next field number |
+---------------------------------------------+----------------------------------+
| ``src/cpmdc_params.c``                      | Render the exact ``&CPMD`` deck  |
|                                             | spelling accepted by OpenCPMD    |
+---------------------------------------------+----------------------------------+
| ``src/cpmdc_features.c``                    | Add catalog and parameter        |
|                                             | feature table rows               |
+---------------------------------------------+----------------------------------+
| ``schema/inventory/cpmd_features.json``     | Add matching inventory metadata  |
+---------------------------------------------+----------------------------------+
| ``schema/inventory/cpmd_cp_keywords.txt``   | Add a base ``&CPMD`` keyword     |
|                                             | when the parser accepts it as a  |
|                                             | top-level control                |
+---------------------------------------------+----------------------------------+
| ``docs/orgmode/reference/cpmd-options.org`` | Document every writable field ID |
+---------------------------------------------+----------------------------------+

Re-run the focused command, then run:

.. code:: bash

   meson test -C build --print-errorlogs
   pixi run -e docs sphinxbld

For parser coverage audits, compare OpenCPMD ``control_utils.mod.F90``
``keyword_contains(line, ...)`` branches against ``catalog.cpmd.*``
rows. Treat context-specific branch tokens as covered by their
structured parent when the schema field renders the accepted deck form,
such as ``LBFGS_NTRUST`` for ``LBFGS NTRUST`` or ``NOSE_IONS_LOCAL`` for
``NOSE IONS LOCAL``.

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
