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
format.

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
``catalog.cpmd.*``, ``catalog.dft.*``, ``params.*``, and ``abi.*``.
Each entry reports whether it applies to stub builds, embedded OpenCPMD
builds, or both.

Session Step Calls (direct-call socket)
=======================================

Callers that drive multiple geometry steps keep method state in a
session and pass a serialized ``ForceInput`` for each step. The
result-carrier entry point writes an unpacked flat ``PotentialResult``
for RPC-style or in-process loops.

.. code:: c

   CPMDCSession *session = cpmdc_session_create(params_bytes, params_size);

   if (session != NULL) {
     double forces_h_bohr[/* n_atoms * 3 */];
     size_t potential_result_capacity =
         cpmdc_potential_result_size_for_force_input(
             force_input_bytes, force_input_size);
     unsigned char *potential_result_bytes =
         malloc(potential_result_capacity);
     size_t potential_result_size = 0;

     CPMDCResult step = cpmdc_session_calculate_forces(
         session, force_input_bytes, force_input_size, forces_h_bohr,
         /* forces_len */ 0 /* set to n_atoms * 3 */);

     CPMDCResult rpc_step = cpmdc_session_calculate_result(
         session, force_input_bytes, force_input_size,
         potential_result_bytes, potential_result_capacity,
         &potential_result_size);

     free(potential_result_bytes);
     cpmdc_session_destroy(session);
   }

When ``potential_result_capacity`` is too small, ``ok == 0``, the
required byte count is written to ``potential_result_size``, and
OpenCPMD is not evaluated (same contract as
``nwchemc_session_calculate_result``).

Units
=====

Native evaluation units on ``CPMDCResult.energy_h`` and force buffers
are Hartree and Hartree/Bohr. ``PotentialResult`` energy and forces are
converted to ``ForceInput.energyUnit`` and ``energyUnit / lengthUnit``
(for example ``eV`` and ``eV/angstrom``).

Topology
========

The first accepted session evaluation fixes atom count and ordered
atomic numbers. Later steps may change coordinates, units, and the 3×3
cell in ``ForceInput.box``. Topology changes require a new session.
