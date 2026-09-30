A force call on a geometry close to the previous one should not start
its SCF from scratch. The converged orbitals of the last call, ``c0`` in
OpenCPMD, are a far better starting point than atomic orbitals, and
keeping them in memory is the main reason to run CPMD inside the host.
This page follows ``c0`` from call to call.

Cold and warm calls
===================

The process keeps the converged orbitals, the cell, and the method that
produced them. Each session keeps a counter of its own successful calls.

.. list-table::
   :header-rows: 1
   :widths: 8 46 46

   * - Call
     - Condition
     - What runs
   * - cold
     - first call of the process, or a changed functional, deck, cell, cutoff, charge, multiplicity, atom count, or elemental composition
     - the full CPMD setup on the composed deck, then the SCF from CPMD's initial guess, or from ``RESTART`` when the deck has ``RESTART WAVEFUNCTION``
   * - warm
     - a previous success in this process kept the same functional, deck, cell, cutoff, charge, multiplicity, and elemental composition, including a new session or a reordering of those atomic numbers
     - new positions into ``tau0``, ``phfac``, then the SCF from the stored ``c0``

A cell counts as unchanged when every component of ``ForceInput.box``
matches the stored cell within 1e-8 Angstrom, and when both calls agree
on having a box at all. A new session has no cell of its own yet, so
that comparison uses the cell stored with the orbitals. The cutoff
matches within 1e-8 Rydberg. The elemental composition is the count of
each atomic number. Reordering the same atoms keeps that composition.
The charge and the multiplicity belong to the same match. A failed warm
call leaves the counter
where it was, so the next call is warm again; a failed cold call clears
the stored results. An SCF that does not converge does not replace the
stored orbitals. When ODIIS exhausts ``MAXITER``, the same call
continues once with ``PCG MINIMIZE``. That pass restores the previous
converged ``c0`` when a converged copy is already stored, and it still
counts as warm. A ``stopgm`` during the call is not that failure: the
wavefunction, forces, and module state of the call are undefined,
``cpmdc`` drops the stored orbitals and the warm cell, and the next call
runs setup again. With more than one MPI rank, ``cpmdc`` aborts the
other ranks.

Where the orbitals are kept
===========================

``opencpmd_embed_rwfopt.patch`` adds a module-level copy of ``c0`` to
``rwfopt_utils``. At the end of ``rwfopt``, every rank saves its own
slice of ``c0`` when ``embed_warm_orbitals`` is true and
``ropt_mod%convwf`` is true. A call that does not converge leaves the
previous copy in place. A failed store allocation calls ``stopgm``.
``cpmdc`` sets ``embed_warm_orbitals`` on every SCF, including the
first. Restore does nothing until a converged copy exists, so the first
call still starts from CPMD's guess and stores ``c0`` only if that SCF
converges. On a later call, ``rwfopt`` first runs ``initrun``, which
sets up the iteration state and scratch arrays as for any SCF, and then
overwrites the generated starting orbitals with the saved copy when its
shape still matches. A different shape calls ``stopgm`` and names
``embed_reset_warm_orbitals``. ``cpmdc`` makes that call when the
functional, the deck, the cutoff, the cell, the charge, the
multiplicity, the atom count, or the elemental composition changes,
before ``rwfopt`` runs.
The SCF then runs to the orbital convergence threshold with the full
``MAXITER`` budget of the deck; a warm call is a complete SCF from a
better start, not a shortened one.

``opencpmd_converged_state.patch`` makes the saved copy trustworthy. In
stock OpenCPMD, ``updwf`` applies an optimiser step to ``c0`` before it
checks convergence, so the orbitals left behind after the last iteration
have moved one step past the ones that produced the reported energy and
forces. The patch checks convergence first and updates only unconverged
orbitals, so the saved ``c0`` is the state that ``forcedr`` used.
Steepest descent with ``iproj <= 1`` sets ``gemax`` above the threshold
before that check: those steps do not converge on ``gemax``.

What clears the orbitals
========================

.. list-table::
   :header-rows: 1
   :widths: 50 50

   * - Event
     - Why
   * - a cell change
     - the plane-wave basis depends on the cell, so the stored coefficients no longer fit
   * - a cutoff change
     - the plane-wave basis changes, so the stored coefficients no longer fit
   * - a change of atom count or elemental composition, or of charge or multiplicity
     - those inputs set the number of states, and a restored ``c0`` of another shape calls ``stopgm``
   * - ``cpmdc_session_set_params()`` or ``cpmdc_session_configure()`` when the functional, cutoff, charge, multiplicity, or deck differs
     - the next evaluation drops the orbitals; the call itself keeps them
   * - a new session with the same functional, deck, cutoff, cell, charge, multiplicity, and elemental composition
     - reuses the stored orbitals and skips setup
   * - a new session whose functional, deck, cutoff, cell, charge, multiplicity, or elemental composition differs
     - the next evaluation drops the stored orbitals and runs setup
   * - a one-shot or global call (``cpmdc_set_params()``, ``cpmdc_energy*()``, ``cpmdc_calculate_result()``) whose basis differs from the stored one
     - that call runs setup; a matching basis reuses the orbitals
   * - ``stopgm`` returns through ``cpmd_stopgm_hook``
     - the call continued past a failed check, so the stored orbitals are not a result

A topology change is refused outright rather than treated as cold: the
atom count and the elemental composition belong to the session. A
reordering of those atomic numbers stays on the session. A new session
with that composition keeps the stored orbitals. An evaluation whose
positions match the previous successful call on that session, including
that reordering, returns the stored energy and forces.

The RESTART file
================

A library force call writes no ``RESTART.1``, ``LATEST``, ``GEOMETRY``,
or ``GEOMETRY.xyz``. The orbitals for the next call stay in memory.
``cpmd.x`` writes those files. A cold library call reads a ``RESTART.1``
only when the deck asks for ``RESTART WAVEFUNCTION``. ``cpmdc`` turns off
CPMD's restart of coordinates, velocities, and the stored geometry
before every SCF. The host owns the geometry, and a
``RESTART COORDINATES`` in the deck would otherwise overwrite the step's
positions.

To seed a cold start from an earlier run with a new geometry, patch the
positions of a ``RESTART.1`` with ``cpmdc-restart patch-positions``. It
rewrites only the coordinate records and copies the wavefunction records
unchanged; the deck then needs ``RESTART WAVEFUNCTION``.

The optimiser on a warm start
=============================

A warm start changes how far the SCF has to go, not how it gets there.
CPMD's optimiser choice then decides the cost: ODIIS reuses its history
of residuals, while ``PCG MINIMIZE`` runs a line minimisation per step.
On the Si3N4 cluster of the :doc:`route comparison <routes>`, the same
minimisation took 515 SCF steps at 1.26 s with ``PCG MINIMIZE`` and 296
steps at 0.64 s with ODIIS, to the same minimum.
