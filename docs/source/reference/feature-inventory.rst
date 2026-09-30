Two copies of one list
======================

``cpmdc`` keeps its capability list twice. ``src/cpmdc_features.c``
compiles it into the library, where ``cpmdc_feature_count()``,
``cpmdc_feature_table()``, and ``cpmdc_feature_find()`` serve it at run
time. ``schema/inventory/cpmd_features.json`` holds the same rows with
metadata for the checks. The Meson suite fails when the two copies, the
schema, the headers, or the docs drift apart.

Feature IDs
===========

+----------------------------------------+---------------------------+----------------------------------------+----------------------+
| Namespace                              | Kind                      | Example                                | Means                |
+========================================+===========================+========================================+======================+
| ``abi.*``                              | ``CPMDC_FEATURE_ABI``     | ``abi.cpmdc_session_calculate_result`` | an exported C        |
|                                        |                           |                                        | function             |
+----------------------------------------+---------------------------+----------------------------------------+----------------------+
| ``params.*``                           | ``CPMDC_FEATURE_PARAMS``  | ``params.cutOffRy``                    | a top-level          |
|                                        |                           |                                        | ``CPMDParams`` field |
+----------------------------------------+---------------------------+----------------------------------------+----------------------+
| ``params.inputSections.<arm>.<field>`` | ``CPMDC_FEATURE_PARAMS``  | ``params.inputSections.cpmd.maxIter``  | a typed field inside |
|                                        |                           |                                        | a section arm        |
+----------------------------------------+---------------------------+----------------------------------------+----------------------+
| ``catalog.section.*``                  | ``CPMDC_FEATURE_SECTION`` | ``catalog.section.PIMD``               | a deck section       |
|                                        |                           |                                        | OpenCPMD's parser    |
|                                        |                           |                                        | reads with           |
|                                        |                           |                                        | ``inscan('&NAME')``  |
+----------------------------------------+---------------------------+----------------------------------------+----------------------+
| ``catalog.cpmd.*``                     | ``CPMDC_FEATURE_KEYWORD`` | ``catalog.cpmd.ODIIS``                 | a ``&CPMD`` keyword  |
|                                        |                           |                                        | or inline option     |
|                                        |                           |                                        | spelling             |
+----------------------------------------+---------------------------+----------------------------------------+----------------------+
| ``catalog.dft.*``                      | ``CPMDC_FEATURE_KEYWORD`` | ``catalog.dft.PCGC``                   | a ``&DFT`` keyword   |
+----------------------------------------+---------------------------+----------------------------------------+----------------------+
| ``section.*``                          | inventory only            | ``section.generic``                    | a                    |
|                                        |                           |                                        | ``CPMDInputSection`` |
|                                        |                           |                                        | union arm            |
+----------------------------------------+---------------------------+----------------------------------------+----------------------+

A ``params.*`` row names a field you can write. A ``catalog.*`` row says
that the renderer can emit a keyword or section, often through a text
field of a typed section rather than a field of its own.
``stub_applicable`` and ``embed_applicable`` say whether the stub build
and the OpenCPMD build expose the row.

Query it from C
===============

.. code:: c

   const CPMDCFeatureEntry *e = cpmdc_feature_find("params.inputSections.cpmd.odiisVectors");
   if (e != NULL && e->embed_applicable) {
     /* this library renders ODIIS with a vector count */
   }
   for (size_t i = 0; i < cpmdc_feature_count(); ++i) {
     const CPMDCFeatureEntry *row = &cpmdc_feature_table()[i];
     /* row->feature_id, row->kind, row->stub_applicable, row->embed_applicable */
   }

The inventory file
==================

Top-level keys of ``schema/inventory/cpmd_features.json``:

+----------------------------------+----------------------------------+
| Key                              | Holds                            |
+==================================+==================================+
| ``version``                      | inventory format version         |
+----------------------------------+----------------------------------+
| ``schema``                       | the schema file it describes     |
+----------------------------------+----------------------------------+
| ``opencpmd_reference``           | the OpenCPMD source it was       |
|                                  | probed against, as a             |
|                                  | description, never a local path  |
+----------------------------------+----------------------------------+
| ``opencpmd_inscan_sections``,    | the 25 sections OpenCPMD parses  |
| ``cpmd_sections``                | with ``inscan``; the two lists   |
|                                  | must be equal                    |
+----------------------------------+----------------------------------+
| ``section_kinds``                | the ``CPMDSectionKind`` enum     |
|                                  | values with their ``section.*``  |
|                                  | IDs                              |
+----------------------------------+----------------------------------+
| ``params_fields``                | the top-level ``CPMDParams``     |
|                                  | fields                           |
+----------------------------------+----------------------------------+
| ``abi_symbols``                  | every public ``cpmdc_*``         |
|                                  | function                         |
+----------------------------------+----------------------------------+
| ``cpmd_keywords``,               | keyword and functional catalogs  |
| ``dft_keywords``,                |                                  |
| ``dft_functionals``              |                                  |
+----------------------------------+----------------------------------+
| ``features``                     | one object per row:              |
|                                  | ``feature_id``, ``kind``,        |
|                                  | ``name``, ``stub_applicable``,   |
|                                  | ``embed_applicable``, ``role``   |
+----------------------------------+----------------------------------+

Plain-text companions in ``schema/inventory/``:
``opencpmd_sections.txt`` (the section allowlist),
``cpmd_cp_keywords.txt`` (base ``&CPMD`` keywords),
``dft_functionals.txt``, and ``surface_map_abi.txt``.

The checks
==========

All run in the default Meson suite; the ``inventory`` suite selects most
of them.

+--------------------------------------+----------------------------------+
| Test                                 | Fails when                       |
+======================================+==================================+
| ``feature-inventory``                | the JSON inventory disagrees     |
|                                      | with                             |
|                                      | ``schema/Potentials.capnp``, the |
|                                      | C table, the public headers, the |
|                                      | README, or the option docs; with |
|                                      | ``CPMD_ROOT`` set, also when the |
|                                      | live ``inscan`` probe finds a    |
|                                      | section the inventory lacks      |
+--------------------------------------+----------------------------------+
| ``cpmd-base-keyword-inventory``      | a base ``&CPMD`` keyword in      |
|                                      | ``cpmd_cp_keywords.txt`` has no  |
|                                      | ``catalog.cpmd.*`` row           |
+--------------------------------------+----------------------------------+
| ``cpmd-params-field-inventory``      | a top-level ``CPMDParams`` field |
|                                      | has no ``params.*`` row, or the  |
|                                      | inventory or the C table has a   |
|                                      | duplicate                        |
+--------------------------------------+----------------------------------+
| ``cpmd-public-abi-inventory``        | a public ``cpmdc_*`` function is |
|                                      | missing from ``abi_symbols``,    |
|                                      | the table, or its implementation |
|                                      | in ``src/cpmdc.c``,              |
|                                      | ``src/cpmdc_stub.c``, or         |
|                                      | ``src/cpmdc_features.c``         |
+--------------------------------------+----------------------------------+
| ``stub-abi-symbol-coverage``         | the stub test does not assert    |
|                                      | every ``abi_symbols`` row        |
+--------------------------------------+----------------------------------+
| ``shared-dlopen-symbol-coverage``    | the ``dlopen`` test does not     |
|                                      | load every ``abi_symbols`` entry |
|                                      | from ``libcpmdc.so``             |
+--------------------------------------+----------------------------------+
| ``cpmd-schema-render-coverage``      | a typed ``CPMDCpmdSection``      |
|                                      | field has no render mapping      |
+--------------------------------------+----------------------------------+
| ``cpmd-option-token-coverage``       | an inline option token in a      |
|                                      | fixture has no render coverage   |
+--------------------------------------+----------------------------------+
| ``cpmd-typed-render-field-coverage`` | a typed ``cpmd``, ``system``,    |
|                                      | ``dft``, or ``atoms`` field is   |
|                                      | absent from the render fixtures  |
|                                      | and assertions                   |
+--------------------------------------+----------------------------------+
| ``cpmd-long-tail-section-coverage``  | a ``CPMDDirectiveSection`` arm   |
|                                      | is missing from the long-tail    |
|                                      | fixture, the render assertions,  |
|                                      | or the option reference          |
+--------------------------------------+----------------------------------+
| ``cpmd-escape-hatch-coverage``       | ``inputBlocks``, ``generic``,    |
|                                      | ``set``, or ``raw`` lose schema, |
|                                      | fixture, render, or docs         |
|                                      | coverage                         |
+--------------------------------------+----------------------------------+

Run them together:

.. code:: bash

   meson test -C build --suite inventory --print-errorlogs

To include the live probe, point ``CPMD_ROOT`` at an OpenCPMD source
tree:

.. code:: bash

   CPMD_ROOT=/path/to/opencpmd meson test -C build feature-inventory -v

A passing live probe prints
``probed /path/to/opencpmd: 25 inscan sections``. Without ``CPMD_ROOT``
the test prints a warning and skips that part.

Add a row
=========

Adding a typed keyword touches the schema, the renderer, the C table,
the inventory file, and the option docs in one change; the
:doc:`contributing guide <../contributing/index>` lists the files and
the focused test command.
