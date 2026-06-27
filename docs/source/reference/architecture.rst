Layers
======

``cpmdc`` mirrors ``nwchemc``: Cap'n Proto on the wire, stable C ABI,
Fortran ``iso_c_binding`` embed that talks to the engine **as a
library**.

- ``include/cpmdc.h`` — language-neutral C ABI (``CPMDCResult``,
  ``CPMDCSession``, socket entry points, feature discovery via
  ``cpmdc_feature_count`` / ``cpmdc_feature_table`` /
  ``cpmdc_feature_find``).
- ``schema/Potentials.capnp`` — ``CPMDParams``, ``ForceInput``,
  ``PotentialResult`` (method knobs vs per-step geometry). Host tools
  may also use
  `readcon-core <https://github.com/HaoZeke/readcon-core>`__ / ``.con``
  landscapes; those become ``ForceInput`` before ``cpmdc``.
- ``src/cpmdc_params.c`` — decode Cap'n Proto only (no CPMD types).
- ``src/cpmdc.c`` — sessions, topology, unit conversion, calls embed.
- ``src/cpmd_embed_c_api.F90`` — ``bind(C)`` surface with a
  deterministic reference PEF; with ``CPMDC_HAS_CPMD`` it links
  ``libcpmd.a`` and sets OpenCPMD **module state** in memory.
- ``src/cpmdc_stub.c`` — same ABI symbols with
  ``cpmdc_available() == 0``.

No CPMD INPUT / ``control`` I/O for configuration
=================================================

Embed does **not** write ``INPUT`` decks or ``CALL control`` to parse
method options.

+----------------------+----------------------+----------------------------+
| Concern              | Carrier              | Where it is applied        |
+======================+======================+============================+
| Method (XC, cutoff,  | Cap'n Proto          | C decode →                 |
| charge, …)           | ``CPMDParams``       | ``cpmdc_embed_set_config`` |
|                      |                      | → ``control_def`` + typed  |
|                      |                      | module writes              |
|                      |                      | (``cntr%ecut``, ``func1``, |
|                      |                      | spin, …)                   |
+----------------------+----------------------+----------------------------+
| Geometry / cell /    | Cap'n Proto          | ``tau0`` / ``parm%a1..a3`` |
| species Z            | ``ForceInput`` (or C | / in-memory species table  |
|                      | arrays)              |                            |
+----------------------+----------------------+----------------------------+
| PP data files        | Paths inside         | PP readers only (data      |
|                      | structured ``atoms`` | files), not method config  |
|                      | fields               |                            |
+----------------------+----------------------+----------------------------+

Upstream CPMD CLI still uses ``&SECTION`` files;
``cpmdc_params_render_input_deck`` exists for **debugging and CLI
parity**, not for the embed hot path. Structured ``CPMDInputSection``
arms include ``cpmd``, ``system``, ``dft``, ``atoms``, ``generic``,
``set``, and ``raw``; ``set`` accepts dotted ``SECTION.KEYWORD`` keys
and merges the keyword into the named section.

Parameter flow
==============

#. Host builds ``CPMDParams`` (pycapnp, rgpot, readcon-core → geometry
   step as ``ForceInput``).
#. ``cpmdc_session_create`` copies bytes and calls
   ``cpmdc_embed_set_config`` with extracted scalars (functional, cutoff
   Ry, charge, multiplicity).
#. Each ``cpmdc_session_calculate_result`` parses ``ForceInput``,
   updates ionic positions in ``coor%tau0`` (Angstrom → a.u.), runs
   library ``wfopts`` / forces, reads ``ener_com%etot`` and
   ``coor%fion``, writes ``PotentialResult``.

Long-running sessions
=====================

First successful eval fixes topology (atom count + ordered Z). Later
steps only change coordinates, units, and cell vectors. Species or count
changes need a new session. Method knobs are fixed at session create
unless ``cpmdc_session_set_params`` runs **before** topology is
accepted.

Build
=====

Without ``libcpmd.a``: stub, Cap'n Proto, sizing, shared dlopen, and
reference PEF tests; the embed shell reports ``available()==1`` and the
stub reports ``available()==0``. With
``-Dwith_cpmd=true -Dcpmd_root=...`` and ``$cpmd_root/lib/libcpmd.a``:
full embed symbols link against the OpenCPMD archive.
