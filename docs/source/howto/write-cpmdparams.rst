Goal
====

Turn a CPMD input deck into a ``CPMDParams`` message file that
``cpmdc_session_create()``, ``rgpot``, or eOn's ``params_path`` can
load. The message carries method setup only. Coordinates, species, and
the cell of each step travel in ``ForceInput``.

Write the message as Cap'n Proto text
=====================================

Cap'n Proto text is the shortest way to write a message by hand, and the
``capnp`` tool encodes it without any code. This message describes a
closed-shell Si and N cluster in a 14 Angstrom cube, BLYP at 70 Ry, with
Kleinman-Bylander pseudopotentials:

.. code:: capnp

   (
     functional = "BLYP",
     cutOffRy = 70.0,
     inputSections = [
       ( cpmd = ( optimizeWavefunction = true,
                  convergenceOrbitals = 1.0e-5,
                  maxIter = 400,
                  odiisVectors = 10,
                  centerMoleculeOff = true ) ),
       ( system = ( symmetry = 1, angstrom = true, cutOffRy = 70.0,
                    cell = [14.0, 1.0, 1.0, 0.0, 0.0, 0.0] ) ),
       ( dft = ( functional = "BLYP", gcCutoff = 1.0e-7, newCode = true ) ),
       ( atoms = ( pseudopotentials = [
           ( element = "Si", path = "Si_MT_BLYP.psp", lmax = 2, loc = 2,
             kleinmanBylander = true ),
           ( element = "N", path = "N_MT_BLYP.psp", lmax = 1, loc = 1,
             kleinmanBylander = true )
         ] ) )
     ]
   )

Save it as ``cluster.params.txt`` and encode it against the schema in
the ``cpmdc`` checkout:

.. code:: bash

   capnp encode schema/Potentials.capnp CPMDParams \
     < cluster.params.txt > cluster.params.bin

``capnp encode`` writes the unpacked flat stream that every ``cpmdc``
entry point reads. It rejects a field name that is not in the schema, so
a typo fails here and not inside CPMD.

Map deck lines to fields
========================

Each deck line has one preferred carrier. Look for a typed field first;
the :doc:`option mapping <../reference/cpmd-options>` lists every one of
them with its feature ID.

+---------------------------------------+------------+---------------------------------------------------------------------------------------------+
| Deck line                             | Section    | Field                                                                                       |
+=======================================+============+=============================================================================================+
| ``OPTIMIZE WAVEFUNCTION``             | ``cpmd``   | ``optimizeWavefunction = true``                                                             |
+---------------------------------------+------------+---------------------------------------------------------------------------------------------+
| ``CONVERGENCE ORBITALS`` / ``1.0d-5`` | ``cpmd``   | ``convergenceOrbitals = 1.0e-5``                                                            |
+---------------------------------------+------------+---------------------------------------------------------------------------------------------+
| ``MAXITER`` / ``400``                 | ``cpmd``   | ``maxIter = 400``                                                                           |
+---------------------------------------+------------+---------------------------------------------------------------------------------------------+
| ``ODIIS`` / ``10``                    | ``cpmd``   | ``odiisVectors = 10``                                                                       |
+---------------------------------------+------------+---------------------------------------------------------------------------------------------+
| ``PCG MINIMIZE``                      | ``cpmd``   | ``directives = [ ( keyword = "PCG MINIMIZE", args = [] ) ]``                                |
+---------------------------------------+------------+---------------------------------------------------------------------------------------------+
| ``CENTER MOLECULE OFF``               | ``cpmd``   | ``centerMoleculeOff = true``                                                                |
+---------------------------------------+------------+---------------------------------------------------------------------------------------------+
| ``SYMMETRY`` / ``1``                  | ``system`` | ``symmetry = 1``                                                                            |
+---------------------------------------+------------+---------------------------------------------------------------------------------------------+
| ``CELL`` /                            | ``system`` | ``cell = [14.0, 1.0, 1.0, 0.0, 0.0, 0.0]``                                                  |
| ``14.0 1.0 1.0 0.0 0.0 0.0``          |            |                                                                                             |
+---------------------------------------+------------+---------------------------------------------------------------------------------------------+
| ``POISSON SOLVER HOCKNEY``            | ``system`` | ``poissonSolver = "HOCKNEY"``                                                               |
+---------------------------------------+------------+---------------------------------------------------------------------------------------------+
| ``CUTOFF`` / ``70``                   | ``system`` | ``cutOffRy = 70.0``                                                                         |
+---------------------------------------+------------+---------------------------------------------------------------------------------------------+
| ``FUNCTIONAL BLYP``                   | ``dft``    | ``functional = "BLYP"``                                                                     |
+---------------------------------------+------------+---------------------------------------------------------------------------------------------+
| ``NEWCODE``                           | ``dft``    | ``newCode = true``                                                                          |
+---------------------------------------+------------+---------------------------------------------------------------------------------------------+
| ``GC-CUTOFF`` / ``1.0d-7``            | ``dft``    | ``gcCutoff = 1.0e-7``                                                                       |
+---------------------------------------+------------+---------------------------------------------------------------------------------------------+
| ``*Si_MT_BLYP.psp KLEINMAN-BYLANDER`` | ``atoms``  | ``( element = "Si", path = "Si_MT_BLYP.psp", lmax = 2, loc = 2, kleinmanBylander = true )`` |
| / ``LMAX=D LOC=D``                    |            |                                                                                             |
+---------------------------------------+------------+---------------------------------------------------------------------------------------------+

``lmax``, ``loc``, and ``skip`` take the channel as an integer: 0
renders ``S``, 1 ``P``, 2 ``D``, 3 ``F``. A keyword without a typed
field goes into that section's ``directives`` list: ``keyword`` is the
line, and each entry of ``args`` becomes one indented line below it.

Defaults that change the physics
================================

Four defaults decide what CPMD computes when the message leaves them
out.

+----------------------------+--------------------------------------------------------------+--------------------------------------------------------------+
| Left out                   | What ``cpmdc`` renders                                       | Consequence                                                  |
+============================+==============================================================+==============================================================+
| the ``system`` section     | ``ANGSTROM`` and ``CUTOFF`` from ``cutOffRy``, no            | the force call chooses the symmetry from the box: 1 for a    |
|                            | ``SYMMETRY`` line and no ``CELL``                            | cube, 8 for unequal orthogonal edges, ``CELL VECTORS`` with  |
|                            |                                                              | no symmetry line when the box is tilted; a zero box is       |
|                            |                                                              | ``SYMMETRY 0`` with the Hockney solver                       |
+----------------------------+--------------------------------------------------------------+--------------------------------------------------------------+
| ``symmetry`` and           | ``SYMMETRY 0`` plus ``POISSON SOLVER HOCKNEY``               | the Hockney isolated-system solver, even for a periodic box  |
| ``poissonSolver``          |                                                              |                                                              |
+----------------------------+--------------------------------------------------------------+--------------------------------------------------------------+
| ``maxIter``                | no ``MAXITER`` line; the OpenCPMD path then inserts          | an SCF that needs more than 40 steps fails the call          |
|                            | ``MAXITER 40`` into the first deck                           |                                                              |
+----------------------------+--------------------------------------------------------------+--------------------------------------------------------------+
| ``kleinmanBylander``       | the pseudopotential line without ``KLEINMAN-BYLANDER``       | Gauss-Hermite integration of the nonlocal projectors, a      |
|                            |                                                              | different energy                                             |
+----------------------------+--------------------------------------------------------------+--------------------------------------------------------------+

Set each of them on purpose.
A raw ``inputBlocks`` entry or a ``raw`` section that already opens
``&CPMD``, ``&SYSTEM``, ``&DFT``, or ``&ATOMS`` is that section.
The default writer does not open the same name again, so a periodic
``SYMMETRY 1`` block is not followed by ``SYMMETRY 0`` and
``POISSON SOLVER HOCKNEY``.
``ForceInput.box`` supplies ``CELL`` only when the deck has no ``CELL``
line.

Check the rendered deck
=======================

``CPMDC_DECK_OUT`` writes the deck ``cpmdc`` renders from the message,
so the check needs no OpenCPMD:

.. code:: bash

   CPMDC_DECK_OUT=cluster.inp build/example_host_step \
     cluster.params.bin step.bin
   cat cluster.inp

``step.bin`` is any encoded ``ForceInput``; the default build evaluates
it with the reference evaluator after the deck is written. Compare
``cluster.inp`` line by line against the deck you started from. On the
OpenCPMD path a second deck, the one CPMD parses with the geometry
merged in, overwrites the first; the
:doc:`deck debugging how-to <debug-deck>` shows both.

Where the geometry-dependent lines go
=====================================

The OpenCPMD path writes ``&ATOMS`` itself from each step's atomic
numbers and positions. An element listed in ``atoms.pseudopotentials``
takes its pseudopotential file, ``LMAX``, ``LOC``, and
``KLEINMAN-BYLANDER`` from that entry. ``KLEINMAN-BYLANDER`` is written
on the ``*file`` line, which is where CPMD reads it, and only when that
entry sets ``kleinmanBylander = true``. The rendered method deck marks
each listed entry with a ``!SPECIES SYM`` comment so the geometry merge
can match the file to an atomic number. An element the message does not
list uses this table:

======= ================== ================
Element File               Channels
======= ================== ================
H       ``H_CVB_BLYP.psp`` ``LMAX=S``
C       ``C_MT_BLYP.psp``  ``LMAX=P``
N       ``N_MT_BLYP.psp``  ``LMAX=P``
O       ``O_MT_BLYP.psp``  ``LMAX=P``
Si      ``Si_MT_BLYP.psp`` ``LMAX=D LOC=D``
Ge      ``Ge_MT_BLYP.psp`` ``LMAX=P``
======= ================== ================

An element in neither the message nor the table fails the evaluation.
The message names the atomic number. A method deck whose ``&ATOMS``
block already has coordinate lines, for example from ``inputBlocks``,
is passed to CPMD unchanged; each step then overwrites only the
positions.
