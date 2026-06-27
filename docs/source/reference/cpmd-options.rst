Contract
========

All method knobs for the OpenCPMD backend are fields of ``CPMDParams``
in ``schema/Potentials.capnp``. Geometry is never part of
``CPMDParams``; it is always a ``ForceInput`` step message (or the raw C
array entry points that take positions and atomic numbers).

Top-level ``CPMDParams`` fields
===============================

+----------------------------------+----------------------------------+
| Field                            | Role                             |
+----------------------------------+----------------------------------+
| ``functional``                   | Default DFT functional (also     |
|                                  | used when no ``dft`` section is  |
|                                  | present)                         |
+----------------------------------+----------------------------------+
| ``cutOffRy``                     | Plane-wave cutoff in Rydberg     |
+----------------------------------+----------------------------------+
| ``charge`` / ``multiplicity``    | System charge and 2S+1           |
|                                  | (multiplicity ``> 1`` implies    |
|                                  | LSD defaults)                    |
+----------------------------------+----------------------------------+
| ``task``                         | Frontend hint: ``energy``,       |
|                                  | ``gradient``, ``md``,            |
|                                  | ``optimize``                     |
+----------------------------------+----------------------------------+
| ``title``                        | Optional comment header in       |
|                                  | rendered decks                   |
+----------------------------------+----------------------------------+
| ``memoryMb``                     | Frontend memory hint             |
+----------------------------------+----------------------------------+
| ``scratchDir``, ``permanentDir`` | ``&CPMD FILEPATH`` runtime file  |
|                                  | placement; ``permanentDir`` wins |
+----------------------------------+----------------------------------+
| ``cpmdRoot``                     | OpenCPMD source/build tree       |
+----------------------------------+----------------------------------+
| ``enginePath``                   | Frontend engine selection hint   |
+----------------------------------+----------------------------------+
| ``inputBlocks``                  | Raw ``&SECTION`` text blocks     |
|                                  | prepended to the deck            |
+----------------------------------+----------------------------------+
| ``inputSections``                | Structured ``CPMDInputSection``  |
|                                  | list                             |
+----------------------------------+----------------------------------+

``CPMDInputSection`` kinds
==========================

+-------------+------------------------+--------------------------+
| Kind        | Cap'n Proto type       | Deck effect              |
+-------------+------------------------+--------------------------+
| ``cpmd``    | ``CPMDCpmdSection``    | ``&CPMD`` (optimize      |
|             |                        | wavefunction, MD,        |
|             |                        | convergence, timestep,   |
|             |                        | …)                       |
+-------------+------------------------+--------------------------+
| ``system``  | ``CPMDSystemSection``  | ``&SYSTEM`` (symmetry,   |
|             |                        | angstrom, cell, cutoff,  |
|             |                        | charge, …)               |
+-------------+------------------------+--------------------------+
| ``dft``     | ``CPMDDftSection``     | ``&DFT`` (functional,    |
|             |                        | LSD, extra directives)   |
+-------------+------------------------+--------------------------+
| ``atoms``   | ``CPMDAtomsSection``   | ``&ATOMS``               |
|             |                        | (pseudopotential lines;  |
|             |                        | coordinates come from    |
|             |                        | ``ForceInput``)          |
+-------------+------------------------+--------------------------+
| ``generic`` | ``CPMDGenericSection`` | Arbitrary ``&NAME`` …    |
|             |                        | ``&END``                 |
+-------------+------------------------+--------------------------+
| ``set``     | ``CPMDSetDirective``   | Embed-path logical       |
|             |                        | ``SECTION.KEYWORD``      |
|             |                        | merged into the section  |
+-------------+------------------------+--------------------------+
| ``raw``     | ``Text``               | Full section text        |
|             |                        | inserted as-is           |
+-------------+------------------------+--------------------------+

Missing core sections receive minimal defaults from the renderer so a
sparse ``CPMDParams`` still produces a valid-looking deck for debugging.
When structured ``dft`` or ``system`` sections provide scalar method
values (functional, cutoff, charge, multiplicity), those section values
are also the effective embed configuration passed through the C ABI.

Escape hatches
==============

Long-tail CPMD keywords go in ``CPMDDirective`` lists inside a section,
``CPMDSetDirective`` dotted keys, ``CPMDInputSection.raw``, or
``inputBlocks``. ``set.key`` uses ``SECTION.KEYWORD`` form, for example
``CPMD.PRINT FORCES ON`` or ``SYSTEM.POISSON SOLVER HOCKNEY``; non-empty
``set.value`` is emitted on the following indented line. Prefer new structured
fields in the schema over inventing a second config file format.

``PotentialConfig``
===================

For rgpot / multi-backend configure RPCs, ``PotentialConfig`` is a
tagged union with ``cpmd @2 :CPMDParams`` (ordinal aligned with the
shared schema; ``nwchem @1`` is reserved void in this package).
