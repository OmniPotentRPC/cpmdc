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

The public C entry points are feature rows too. Use those rows when a
host needs to negotiate the serialized result path before allocating
Cap'n Proto buffers:

+----------------------+----------------------------------------+----------------------+
| Host action          | ABI feature ID                         | Carrier              |
+======================+========================================+======================+
| Accept a             | ``abi.cpmdc_set_params``               | ``CPMDParams`` bytes |
| process-wide method  |                                        |                      |
| buffer               |                                        |                      |
+----------------------+----------------------------------------+----------------------+
| Create a reusable    | ``abi.cpmdc_session_create``           | ``CPMDParams`` bytes |
| method session       |                                        |                      |
+----------------------+----------------------------------------+----------------------+
| Replace method setup | ``abi.cpmdc_session_set_params``       | ``CPMDParams`` bytes |
| before topology is   |                                        |                      |
| accepted             |                                        |                      |
+----------------------+----------------------------------------+----------------------+
| Evaluate one session | ``abi.cpmdc_session_calculate_result`` | ``ForceInput`` to    |
| step into a result   |                                        | ``PotentialResult``  |
| message              |                                        |                      |
+----------------------+----------------------------------------+----------------------+
| Evaluate one         | ``abi.cpmdc_calculate_result``         | ``CPMDParams`` plus  |
| one-shot serialized  |                                        | ``ForceInput`` to    |
| calculation          |                                        | ``PotentialResult``  |
+----------------------+----------------------------------------+----------------------+

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
optional 3 by 3 cell, and requested units. A species change or
atom-count change requires a new ``CPMDCSession``.

``cpmdc_session_set_params()`` can replace method setup only before
topology is accepted. Once a step has succeeded, method changes also
need a new session.

CPMD Input Ownership
====================

``CPMDParams`` owns method and backend setup. ``ForceInput`` owns
coordinates, atomic numbers, unit strings, and the optional 3 by 3 cell
for a single evaluation. The runtime renders a CPMD ``INPUT`` deck from
``CPMDParams`` and merges the step geometry into it on the first
evaluation of a session.

Structured ``inputSections`` should be preferred over raw blocks when a
typed arm exists:

- ``system`` for cells, cutoff/grid/mesh controls, state occupation,
  external fields, CDFT (constrained DFT) Gaussian controls, pressure
  and stress controls, Poisson settings, isolated-shape controls,
  charge, and spin multiplicity
- ``cpmd`` for wavefunction optimization, MD, convergence, restart, and
  trajectory controls
- ``dft`` for functional and spin-polarized DFT controls
- ``atoms`` for pseudopotential grouping and fixed non-coordinate
  ``&ATOMS`` directives
- named OpenCPMD section arms such as ``pimd``, ``vdw``, ``linres``, and
  ``tddft`` for keyword/value lines and nested subsection blocks

Use ``generic`` for non-catalog aliases, ``set`` for merge-only keywords
inside a named section, and ``raw`` or ``inputBlocks`` for
text-preserving deck fragments. On a geometry render that is not the
OpenCPMD merge, typed ``atoms`` sections must cover every element
present in the step, and a missing ``atoms`` section covers H and O
only. On the OpenCPMD path the ``&ATOMS`` block CPMD reads is rebuilt
from each step. Elements listed in ``atoms.pseudopotentials`` supply
the file, ``LMAX``, ``LOC``, and ``KLEINMAN-BYLANDER`` on the ``*file``
line. Elements the message omits use the built-in table for H, C, N, O,
Si, and Ge. An element in neither is an error. See
:doc:`writing a CPMDParams message <write-cpmdparams>`.

Choose the least lossy structured carrier:

+-----------------------+----------------------+-------------------------------------------------+
| Need                  | Carrier              | Feature ID to discover                          |
+=======================+======================+=================================================+
| A field exists in     | typed ``cpmd`` arm   | ``params.inputSections.cpmd.maxIter``           |
| ``CPMDCpmdSection``   |                      |                                                 |
+-----------------------+----------------------+-------------------------------------------------+
| A field exists in     | typed ``system`` arm | ``params.inputSections.system.cell``            |
| ``CPMDSystemSection`` |                      |                                                 |
+-----------------------+----------------------+-------------------------------------------------+
| A field exists in     | typed ``dft`` arm    | ``params.inputSections.dft.hfxScreening``       |
| ``CPMDDftSection``    |                      |                                                 |
+-----------------------+----------------------+-------------------------------------------------+
| Pseudopotentials and  | typed ``atoms`` arm  | ``params.inputSections.atoms.pseudopotentials`` |
| fixed non-coordinate  |                      |                                                 |
| ``&ATOMS`` controls   |                      |                                                 |
+-----------------------+----------------------+-------------------------------------------------+
| A catalog section     | matching typed       | ``catalog.section.PIMD`` plus                   |
| contains              | long-tail arm        | ``params.inputSections.pimd.directives``        |
| keyword/value lines   |                      |                                                 |
+-----------------------+----------------------+-------------------------------------------------+
| A catalog section     | matching typed       | ``catalog.section.VDW`` plus                    |
| contains nested       | long-tail arm        | ``params.inputSections.vdw.subsections``        |
| blocks                |                      |                                                 |
+-----------------------+----------------------+-------------------------------------------------+
| One merge-only        | ``set`` with         | ``params.inputSections.set.key``                |
| keyword belongs       | ``SECTION.KEYWORD``  |                                                 |
| inside a typed        |                      |                                                 |
| section               |                      |                                                 |
+-----------------------+----------------------+-------------------------------------------------+
| A non-catalog alias   | ``generic``          | ``params.inputSections.generic.directives``     |
| can be expressed as   |                      |                                                 |
| keyword/argument      |                      |                                                 |
| lines                 |                      |                                                 |
+-----------------------+----------------------+-------------------------------------------------+
| A complete section    | ``raw``              | ``params.inputSections.raw``                    |
| fragment must be      |                      |                                                 |
| preserved as text     |                      |                                                 |
+-----------------------+----------------------+-------------------------------------------------+
| A complete deck       | top-level            | ``params.inputBlocks``                          |
| fragment must be      | ``inputBlocks``      |                                                 |
| prepended before      |                      |                                                 |
| structured sections   |                      |                                                 |
+-----------------------+----------------------+-------------------------------------------------+

RESTART files
=============

The binary ``RESTART`` file is a separate object from the deck keyword.
Read and rewrite it with ``cpmdc-restart`` or ``cpmdc_restart.h``
(``libcpmdc_restart``). Coordinates in that file are Bohr in species
order. A coordinate patch copies the wavefunction records unchanged,
which is how a host keeps orbitals while it owns the next geometry. Put
``RESTART WAVEFUNCTION`` in the deck to start the first SCF of a session
from the patched file; the OpenCPMD path ignores restarted coordinates,
velocities, and geometry, because each step supplies them.

Units
=====

After a successful in-process evaluation
(``cpmdc_session_energy_forces``, ``cpmdc_session_calculate_result``, or
one-shot ``cpmdc_energy*``), call ``cpmdc_last_energy_components`` to
read the OpenCPMD ``ener_com`` snapshot (``etot``, ``ekin``, ``epseu``,
``enl``, ``eht``, ``exc``, and the rest) in Hartree without parsing CLI
``ENERGY`` files. ``etot`` matches ``CPMDCResult.energy_h`` for that
step.

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
executable build with the patches in ``tools/``; see
:doc:`building the archive <opencpmd-archive>`.

Set ``CPMDC_PSEUDO_DIR`` to the directory holding the pseudopotential
files. The first evaluation of a session changes into that directory
while CPMD reads them (``ratom`` and ``recpnew``), then returns to the
host's directory. ``CPMD_PP_LIBRARY_PATH`` is the fallback when
``CPMDC_PSEUDO_DIR`` is unset; with neither set, the first evaluation
fails.

``CPMDC_DECK_OUT`` names a file that receives the method deck cpmdc
renders from ``CPMDParams``; see :doc:`debugging a deck <debug-deck>`
for what the OpenCPMD path adds before CPMD parses it.

Nuclear forces from repeated SCF calls
======================================

A live OpenCPMD link evaluates each geometry with fixed nuclei
(``OPTIMIZE WAVEFUNCTION``) and then exports nuclear forces from
``coor%fion``, the same Born-Oppenheimer force definition as a
single-point force evaluation after the SCF (self-consistent field).

On the embed path the Fortran bridge always calls
``cpmdc_set_need_forces(.TRUE.)`` before ``wfopts``. OpenCPMD must be
patched so ``rwfopt`` sets
``tfor = (iprint_force == 1) .OR. cpmdc_need_forces``. Without that,
OpenCPMD zeros ``fion`` after ``forcedr`` when ``tfor`` is false and the
C force buffer stays all zeros even though the energy is finite.

Warm calls (same process, same session, same cell, new ``ForceInput``
positions):

#. First force call: cold setup once, from the in-memory deck.
#. Later calls: update ``tau0`` from C arrays, ``phfac``, full SCF from
   the retained orbitals (``cpmdc_set_warm_orbitals``), no second setup.
#. Every call requests forces through ``cpmdc_set_need_forces`` and
   keeps the deck's ``MAXITER``;
   ``tests/test_embed_warm_no_nomore_clamp.py`` and
   ``tests/test_embed_bomd_force_export.py`` guard both.

:doc:`Wavefunction state <../explanation/wavefunction-state>` explains
what makes a call cold or warm.

Prefer an off-equilibrium water or similar system when checking forces.
Some tightly packed cluster geometries can report a zero nuclear
gradient even under native ``cpmd.x``; an energy that changes between
calls does not show that the forces are right.

Example environment for a live driver:

.. code:: bash

   export CPMDC_LIBRARY=/path/to/libcpmdc.so
   export CPMDC_PSEUDO_DIR=/path/to/pseudopotentials
   # one session, many ForceInput steps via cpmdc_session_calculate_forces
