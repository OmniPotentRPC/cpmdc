Goal
====

See the CPMD input deck that ``cpmdc`` renders from a ``CPMDParams``
message, and work out what the OpenCPMD path changes before CPMD
(Car-Parrinello Molecular Dynamics) parses it.

Write the deck to a file
========================

Set ``CPMDC_DECK_OUT`` to a file name before the process configures
``cpmdc``:

.. code:: bash

   export CPMDC_DECK_OUT=$PWD/method.inp
   build/example_host_step params.bin step.bin
   cat method.inp

``cpmdc`` writes the file each time it renders a message: in
``cpmdc_session_create()``, ``cpmdc_session_set_params()``,
``cpmdc_set_params()``, the ``PotentialConfig`` entry points, and every
one-shot call. Each write replaces the previous content, so the file
holds the deck of the most recent configuration. The write happens in
the default build too, so a message can be checked on a machine without
OpenCPMD. An empty ``CPMDC_DECK_OUT``, or a path that cannot be opened
for writing, skips the write without an error.

Under ``mpirun`` every rank renders the same deck to the same path.
Point each rank at its own file when the ranks share a file system, for
example through a wrapper script that appends ``$OMPI_COMM_WORLD_RANK``.

What the OpenCPMD path changes
==============================

The file holds the method deck. On the first evaluation of a session,
the OpenCPMD path builds the deck that CPMD parses from it in three
steps:

#. If the method deck has an ``&ATOMS`` block with coordinate lines, it
   is used as it is. Otherwise every ``&ATOMS`` block is dropped and a
   new one is written from the step's atomic numbers and positions. A
   ``!SPECIES SYM`` comment in the method deck selects that element's
   file, ``LMAX``, ``LOC``, and, when the ``*file`` line has it,
   ``KLEINMAN-BYLANDER``. An element with no such comment uses the
   built-in table (see :doc:`writing a CPMDParams message
   <write-cpmdparams>`). An element in neither fails before CPMD starts,
   and the error names the atomic number. When the method deck has no
   ``&CPMD``, ``&SYSTEM``, or ``&DFT`` section at all, a minimal deck
   replaces it: ``OPTIMIZE WAVEFUNCTION``, ``CONVERGENCE ORBITALS
   1.0d-5``, ``MAXITER 40``, ``ANGSTROM``, and the functional, cutoff,
   charge, and multiplicity of the message. The symmetry comes from the
   step's box. A cube is ``SYMMETRY 1``. Unequal orthogonal edges are
   ``SYMMETRY 8``. A tilted box is ``CELL VECTORS`` with no ``SYMMETRY``
   line. A missing or zero box is ``SYMMETRY 0`` with ``POISSON SOLVER
   HOCKNEY``.
#. If the deck has no ``CELL`` keyword, a ``CELL`` line from
   ``ForceInput.box`` is inserted after ``ANGSTROM`` in ``&SYSTEM``, or
   a 12 Angstrom cube when the step has no box. A tilted box is written
   as ``CELL VECTORS``.
#. If the deck names no ``SYMMETRY``, no ``POISSON SOLVER``, no
   ``CLUSTER``, and no isolated-molecule keyword, the same symmetry
   choice is written into ``&SYSTEM``. An explicit symmetry line, Hockney
   solver, ``CLUSTER``, or isolated-molecule keyword is left as written.
   ``CELL VECTORS`` is not given a ``SYMMETRY`` line.
#. If the deck has no ``MAXITER`` keyword, ``MAXITER 40`` is inserted
   after ``&CPMD``.

The method text is kept in full, including a block longer than a few
thousand characters. A preview buffer shorter than the composed deck
returns an error and does not keep a shortened copy.

So read ``method.inp`` with those edits in mind. The ``&ATOMS``
block in the file lists the pseudopotential lines with a placeholder
atom count of 0; the coordinates come from each step.

Checks worth making on the file
===============================

+----------------------------------+--------------------------------------------------------+
| Look for                         | Why                                                    |
+==================================+========================================================+
| ``SYMMETRY`` and                 | ``SYMMETRY 0`` with ``POISSON SOLVER HOCKNEY`` is an   |
| ``POISSON SOLVER`` in ``&SYSTEM``| isolated system: a zero box, or a deck that asked for  |
|                                  | isolation. A positive-volume box with no such keyword  |
|                                  | is symmetry 1, symmetry 8, or ``CELL VECTORS`` with no |
|                                  | symmetry line                                          |
+----------------------------------+--------------------------------------------------------+
| ``CELL``                         | present means the message fixes the cell; absent means |
|                                  | each step's box supplies it                            |
+----------------------------------+--------------------------------------------------------+
| ``KLEINMAN-BYLANDER`` on the     | on the ``LMAX`` line CPMD ignores it and integrates    |
| ``*file`` line                   | the projectors by Gauss-Hermite                        |
+----------------------------------+--------------------------------------------------------+
| ``MAXITER``                      | absent means the OpenCPMD path caps the SCF at 40      |
|                                  | steps                                                  |
+----------------------------------+--------------------------------------------------------+
| the optimiser keyword            | decides the cost per SCF step; see                     |
| (``ODIIS``, ``PCG``,             | :doc:`choosing the optimiser <wavefunction-optimiser>` |
| ``PCG MINIMIZE``)                |                                                        |
+----------------------------------+--------------------------------------------------------+
| ``FUNCTIONAL`` in ``&DFT``       | the functional CPMD uses; the ``dft`` section wins     |
|                                  | over the top-level ``functional``                      |
+----------------------------------+--------------------------------------------------------+

Compare with a deck that works
==============================

When a hand-written CPMD deck gives the expected energy under ``cpmd.x``
and the library does not, compare the two after the edits above:

.. code:: bash

   diff -u reference.inp method.inp

Differences in ``&ATOMS`` are expected: coordinates, atom counts, and
the species order come from the step. Every other difference is a field
that the message sets differently or leaves to a default.
