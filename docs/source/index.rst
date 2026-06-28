.. image:: /_static/logo-light.svg
   :class: light-only
   :width: 460
   :align: center

.. image:: /_static/logo-dark.svg
   :class: dark-only
   :width: 460
   :align: center

Overview
========

``cpmdc`` is the C ABI layer for embedding
`OpenCPMD <https://github.com/OpenCPMD/CPMD>`__ in OmniPotentRPC tools.
It keeps the public boundary small and language-neutral. Hosts pass
method setup once, pass geometry for each calculation step, and get
native C results or a serialized result message back.

+---------------------+-----------------------+-----------------------+
| Message             | Owns                  | Typical lifetime      |
+=====================+=======================+=======================+
| ``CPMDParams``      | method setup, CPMD    | one session           |
|                     | sections,             |                       |
|                     | pseudopotentials,     |                       |
|                     | backend hints         |                       |
+---------------------+-----------------------+-----------------------+
| ``ForceInput``      | positions, atomic     | one calculation step  |
|                     | numbers, optional     |                       |
|                     | cell, requested units |                       |
+---------------------+-----------------------+-----------------------+
| ``PotentialResult`` | energy, forces,       | one calculation       |
|                     | status, message       | result                |
+---------------------+-----------------------+-----------------------+

The ABI does not expose C++ types, Rust types, or a second JSON/TOML
options language. All portable configuration lives in
``schema/Potentials.capnp``; the C library decodes those bytes with
generated ``capnp-c`` readers and renders the CPMD ``INPUT`` deck
internally.

Long-running drivers create one ``CPMDCSession`` from ``CPMDParams``,
then pass ``ForceInput`` per step. The result-carrier call,
``cpmdc_session_calculate_result()``, writes the same
``PotentialResult`` message used by ``rgpot`` RPC clients.

The first accepted session evaluation fixes atom count and ordered
atomic numbers. Later steps may change coordinates, units, and the 3x3
cell; topology changes require a new session.

Start Here
==========

.. list-table::
   :header-rows: 1

   * - You are
     - Start with
     - Page
   * - Building the project for the first time
     - Run the default Meson suite and learn what it covers
     - :doc:`Quickstart <tutorials/quickstart>`
   * - Calling cpmdc from a host process
     - Use the session result path and size output from ``ForceInput``
     - :doc:`Embedding cpmdc <howto/embedding>`
   * - Constructing ``CPMDParams``
     - Pick typed section fields before raw deck text
     - :doc:`CPMD option mapping <reference/cpmd-options>`
   * - Linking a real OpenCPMD archive
     - Check the C/Fortran layer boundary and runtime data ownership
     - :doc:`Architecture <reference/architecture>`
   * - Checking exported C symbols
     - Inspect generated declarations and feature structs
     - :doc:`API reference <api/index>`

First Command
=============

.. code:: bash

   meson setup build -Dwith_tests=true
   meson compile -C build
   meson test -C build --print-errorlogs

The default build does not need an OpenCPMD checkout. It uses a
deterministic reference evaluator to test the public ABI, parser,
feature table, sessions, result sizing, and unit conversion.

Runtime Contracts
=================

+-------------------+------------------------+-------------------------------------+
| Concern           | Carrier                | Public entry points                 |
+===================+========================+=====================================+
| Method setup      | serialized             | ``cpmdc_session_create``,           |
|                   | ``CPMDParams``         | ``cpmdc_session_set_params``,       |
|                   |                        | ``cpmdc_set_params``                |
+-------------------+------------------------+-------------------------------------+
| Geometry step     | serialized             | ``cpmdc_session_calculate_result``, |
|                   | ``ForceInput`` or C    | ``cpmdc_session_calculate_forces``, |
|                   | arrays                 | ``cpmdc_session_energy_forces``     |
+-------------------+------------------------+-------------------------------------+
| Result carrier    | serialized             | ``cpmdc_session_calculate_result``, |
|                   | ``PotentialResult`` or | ``cpmdc_calculate_result``,         |
|                   | native C values        | ``CPMDCResult``                     |
+-------------------+------------------------+-------------------------------------+
| Runtime discovery | stable feature IDs     | ``cpmdc_feature_count``,            |
|                   |                        | ``cpmdc_feature_table``,            |
|                   |                        | ``cpmdc_feature_find``              |
+-------------------+------------------------+-------------------------------------+

An OpenCPMD archive build links the same C ABI surface against
``libcpmd.a``.

.. toctree::
   :maxdepth: 2
   :caption: Tutorials

   tutorials/quickstart

.. toctree::
   :maxdepth: 2
   :caption: How-To Guides

   howto/embedding

.. toctree::
   :maxdepth: 2
   :caption: Reference

   reference/architecture
   reference/cpmd-options
   api/index

.. toctree::
   :maxdepth: 2
   :caption: Contributing

   contributing/index

.. toctree::
   :maxdepth: 1
   :caption: Project

   changelog
