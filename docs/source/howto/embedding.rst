Pick The Entry Point
====================

+----------------+------------------------------------+----------------+-----------------------+
| Need           | Entry point                        | Input          | Output                |
+================+====================================+================+=======================+
| Repeated RPC   | ``cpmdc_session_calculate_result`` | one            | serialized            |
| or optimizer   |                                    | ``CPMDParams`` | ``PotentialResult``   |
| steps          |                                    | session, one   |                       |
|                |                                    | ``ForceInput`` |                       |
|                |                                    | per step       |                       |
+----------------+------------------------------------+----------------+-----------------------+
| Repeated       | ``cpmdc_session_energy_forces``    | one            | Hartree plus          |
| native C       |                                    | ``CPMDParams`` | Hartree/Bohr forces   |
| forces         |                                    | session, C     |                       |
|                |                                    | coordinate     |                       |
|                |                                    | arrays         |                       |
+----------------+------------------------------------+----------------+-----------------------+
| Single         | ``cpmdc_calculate_result``         | ``CPMDParams`` | serialized            |
| serialized     |                                    | and            | ``PotentialResult``   |
| calculation    |                                    | ``ForceInput`` |                       |
|                |                                    | bytes          |                       |
+----------------+------------------------------------+----------------+-----------------------+
| Scalar         | ``cpmdc_energy_forces``            | C coordinate   | Hartree plus          |
| compatibility  |                                    | arrays plus    | Hartree/Bohr forces   |
| call           |                                    | ``CPMDParams`` |                       |
|                |                                    | bytes          |                       |
+----------------+------------------------------------+----------------+-----------------------+
| Capability     | ``cpmdc_feature_find``             | stable feature | ``CPMDCFeatureEntry`` |
| discovery      |                                    | ID             | or ``NULL``           |
+----------------+------------------------------------+----------------+-----------------------+

Use the session result path for new drivers. It keeps method setup and
topology state in one object while preserving the same
``PotentialResult`` carrier used by RPC frontends.

The compatibility array calls are useful when a host already owns native
C arrays and only needs Hartree/Hartree-per-Bohr values. The serialized
result calls are the better boundary for RPC frontends, cross-language
bindings, and drivers that already use the shared OmniPotentRPC schema.

C ABI
=====

Include ``cpmdc.h`` and pass serialized Cap'n Proto message bytes to the
ABI functions. Wire structs are generated from
``schema/Potentials.capnp``; the header does not redefine ``CPMDParams``
fields as a second C struct language.

.. code:: c

   #include <cpmdc.h>
   #include <stdlib.h>

   int rc = cpmdc_set_params(params_bytes, params_size);
   CPMDCResult result = cpmdc_energy_gradient(
       n_atoms, positions_ang, atomic_numbers, params_bytes, params_size,
       gradient_h_bohr);
   CPMDCResult forces = cpmdc_energy_forces(
       n_atoms, positions_ang, atomic_numbers, params_bytes, params_size,
       forces_h_bohr);

``params_bytes`` is an unpacked flat Cap'n Proto message whose root is
``CPMDParams``. It can come from pycapnp, a memory-mapped file, rgpot,
or another Cap'n Proto binding that writes the standard flat stream
format. Use the scalar calls when the caller already owns arrays in
Angstrom and wants Hartree or Hartree/Bohr output directly.

For a complete host-side skeleton, see ``examples/host_step.c``. The
checked ``example-host-step`` Meson test runs that program with
generated fixture bytes and expects output shaped like:

.. code:: text

   energy_h=...
   potential_result_size_bytes=...
   message=...

Message Flow
============

There are two wire messages in a normal embedded driver:

+----------------+----------------------+------------------------------+
| Message        | Lifetime             | Contents                     |
+================+======================+==============================+
| ``CPMDParams`` | Session setup        | method, structured CPMD      |
|                |                      | sections, pseudopotentials,  |
|                |                      | engine hints                 |
+----------------+----------------------+------------------------------+
| ``ForceInput`` | One calculation step | positions, atomic numbers,   |
|                |                      | optional cell, requested     |
|                |                      | output units                 |
+----------------+----------------------+------------------------------+

A host should build ``CPMDParams`` once, serialize it as an unpacked
flat Cap'n Proto message, and create a ``CPMDCSession``. Each geometry
step is a separate ``ForceInput`` message passed to
``cpmdc_session_calculate_result()`` or the lower-level force calls. The
``PotentialResult`` path converts output to ``ForceInput.energyUnit``
and ``ForceInput.energyUnit / ForceInput.lengthUnit``.

Feature Discovery
=================

``cpmdc.h`` also exposes a small feature table. Embedders can inspect it
before building inputs, hiding unsupported controls, or deciding whether
a stub build is sufficient for a workflow.

.. code:: c

   size_t feature_count = cpmdc_feature_count();
   const CPMDCFeatureEntry *features = cpmdc_feature_table();
   const CPMDCFeatureEntry *pimd =
       cpmdc_feature_find("catalog.section.PIMD");
   const CPMDCFeatureEntry *result_call =
       cpmdc_feature_find("abi.cpmdc_session_calculate_result");

   if (pimd != NULL && pimd->embed_applicable) {
     /* The embedded OpenCPMD build can render a &PIMD section. */
   }

Feature IDs use namespaces such as ``catalog.section.*``,
``catalog.cpmd.*``, ``catalog.dft.*``, ``params.*``, and ``abi.*``. Each
entry reports whether it applies to stub builds, embedded OpenCPMD
builds, or both. Structured field controls are discoverable with
``params.inputSections.<section>.<field>`` IDs, for example
``params.inputSections.cpmd.maxIter``,
``params.inputSections.dft.hfxScreening``, and
``params.inputSections.atoms.pseudopotentials``.

Session Step Calls (direct-call socket)
=======================================

Callers that drive multiple geometry steps keep method state in a
session and pass a serialized ``ForceInput`` for each step. The
result-carrier entry point writes an unpacked flat ``PotentialResult``
for RPC-style or in-process loops.

.. code:: c

   CPMDCSession *session = cpmdc_session_create(params_bytes, params_size);

   if (session != NULL) {
     size_t forces_len = n_atoms * 3;
     double *forces_h_bohr = calloc(forces_len, sizeof(*forces_h_bohr));
     size_t potential_result_capacity =
         cpmdc_potential_result_size_for_force_input(
             force_input_bytes, force_input_size);
     unsigned char *potential_result_bytes =
         malloc(potential_result_capacity);
     size_t potential_result_size = 0;

     CPMDCResult step = cpmdc_session_calculate_forces(
         session, force_input_bytes, force_input_size, forces_h_bohr,
         forces_len);

     CPMDCResult rpc_step = cpmdc_session_calculate_result(
         session, force_input_bytes, force_input_size,
         potential_result_bytes, potential_result_capacity,
         &potential_result_size);

     free(forces_h_bohr);
     free(potential_result_bytes);
     cpmdc_session_destroy(session);
   }

When ``potential_result_capacity`` is too small, ``ok == 0``, the
required byte count is written to ``potential_result_size``, and
OpenCPMD is not evaluated (same contract as
``nwchemc_session_calculate_result``).

Use the session path for optimizers, molecular dynamics drivers, and RPC
frontends. It avoids reparsing method setup on every step and gives the
runtime one place to enforce topology consistency.

Result Buffer Contract
======================

Serialized result calls write an unpacked flat ``PotentialResult``. Size
the output buffer from the step message before evaluating:

.. code:: c

   size_t needed = cpmdc_potential_result_size_for_force_input(
       force_input_bytes, force_input_size);
   unsigned char *out = malloc(needed);
   size_t wrote = 0;

   CPMDCResult r = cpmdc_session_calculate_result(
       session, force_input_bytes, force_input_size, out, needed, &wrote);

``needed == 0`` means the ``ForceInput`` message is invalid or too large
for the C ABI. If ``needed`` is positive but the supplied capacity is
smaller, the call returns ``ok == 0``, writes the required byte count to
``wrote``, and skips evaluation.

Session Lifetime
================

``cpmdc_session_create()`` copies the serialized ``CPMDParams`` buffer.
Callers may free or reuse their original input bytes after the session
is created.

The first successful session evaluation fixes the topology: atom count
and ordered atomic numbers. Later steps may change coordinates, the
optional 3x3 cell, and requested units. A species change or atom-count
change requires a new ``CPMDCSession``.

``cpmdc_session_set_params()`` can replace method setup only before
topology is accepted. Once a step has succeeded, method changes also
need a new session.

CPMD Input Ownership
====================

``CPMDParams`` owns method and backend setup. ``ForceInput`` owns
coordinates, atomic numbers, unit strings, and the optional 3x3 cell for
a single evaluation. The runtime merges those carriers into a CPMD
``INPUT`` deck for each accepted step.

Structured ``inputSections`` should be preferred over raw blocks when a
typed arm exists:

- ``system`` for ``CELL``, ``CUTOFF``, ``SCALE``, ``CHARGE``, and spin
  multiplicity
- ``cpmd`` for wavefunction optimization, MD, convergence, restart, and
  trajectory controls
- ``dft`` for functional and spin-polarized DFT controls
- ``atoms`` for pseudopotential grouping and fixed non-coordinate
  ``&ATOMS`` directives

Use ``generic``, ``set``, or ``raw`` for CPMD sections that do not have
a typed arm. Typed ``atoms`` sections must cover every element present
in the step geometry. Without an explicit ``atoms`` section, built-in
BLYP defaults cover H and O only.

Choose the least lossy structured carrier:

+----------------------------------+----------------------------------+
| Need                             | Carrier                          |
+==================================+==================================+
| A field exists in                | typed section arm                |
| ``CPMDCpmdSection``,             |                                  |
| ``CPMDDftSection``,              |                                  |
| ``CPMDSystemSection``, or        |                                  |
| ``CPMDAtomsSection``             |                                  |
+----------------------------------+----------------------------------+
| One unsupported keyword belongs  | ``set`` with ``SECTION.KEYWORD`` |
| inside a typed section           |                                  |
+----------------------------------+----------------------------------+
| A whole unsupported section can  | ``generic``                      |
| be expressed as keywords and     |                                  |
| argument lists                   |                                  |
+----------------------------------+----------------------------------+
| A complete deck fragment must be | ``raw`` or top-level             |
| preserved as text                | ``inputBlocks``                  |
+----------------------------------+----------------------------------+

Units
=====

Native evaluation units on ``CPMDCResult.energy_h`` and force buffers
are Hartree and Hartree/Bohr. ``PotentialResult`` energy and forces are
converted to ``ForceInput.energyUnit`` and ``energyUnit / lengthUnit``
(for example ``eV`` and ``eV/angstrom``).

``ForceInput.pos`` and ``ForceInput.box`` use ``ForceInput.lengthUnit``.
The native array calls take positions in Angstrom regardless of the
schema defaults.

OpenCPMD Runtime Inputs
=======================

Archive builds require ``-Dwith_cpmd=true`` and
``-Dcpmd_root=/path/to/OpenCPMD``. That tree must contain
``lib/libcpmd.a`` and the OpenCPMD object/include files from a completed
executable build.

Pseudopotential names in ``atoms.pseudopotentials`` can be absolute
paths or library-style names. For library-style names, set
``CPMDC_PSEUDO_DIR`` to the directory containing the files. The embed
layer also checks common pseudopotential directories under ``cpmdRoot``.
