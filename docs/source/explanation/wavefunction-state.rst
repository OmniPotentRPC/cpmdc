A force call on a geometry close to the previous one should not start
its SCF from scratch. The converged orbitals of the last call, ``c0`` in
OpenCPMD, are a far better starting point than atomic orbitals, and
keeping them in memory is the main reason to run CPMD inside the host.
This page follows ``c0`` from call to call.

Cold and warm calls
===================

Each session keeps a counter of successful warm-capable calls and the
cell of the call that started them.

+------+-------------------------------+-------------------------------+
| Call | Condition                     | What runs                     |
+======+===============================+===============================+
| cold | first call of a session, a    | the full CPMD setup on the    |
|      | changed cell, or the session  | composed deck, then the SCF   |
|      | was re-configured             | from CPMD's initial guess, or |
|      |                               | from ``RESTART`` when the     |
|      |                               | deck has                      |
|      |                               | ``RESTART WAVEFUNCTION``      |
+------+-------------------------------+-------------------------------+
| warm | the previous call of this     | new positions into ``tau0``,  |
|      | session succeeded and the     | ``phfac``, then the SCF from  |
|      | cell is unchanged             | the stored ``c0``             |
+------+-------------------------------+-------------------------------+

A cell counts as unchanged when every component of ``ForceInput.box``
matches the stored cell within 1e-8 Angstrom, and when both calls agree
on having a box at all. A failed warm call leaves the counter where it
was, so the next call is warm again; a failed cold call clears the
stored results.

Where the orbitals are kept
===========================

``opencpmd_warm_orbitals.patch`` adds a module-level copy of ``c0`` to
``rwfopt_utils``. At the end of each ``rwfopt``, every rank saves its
own slice of ``c0``. On a warm call, ``rwfopt`` first runs ``initrun``,
which sets up the iteration state and scratch arrays as for any SCF, and
then overwrites the generated starting orbitals with the saved copy when
its dimensions still match. The SCF then runs to the orbital convergence
threshold with the full ``MAXITER`` budget of the deck; a warm call is a
complete SCF from a better start, not a shortened one.

``opencpmd_converged_state.patch`` makes the saved copy trustworthy. In
stock OpenCPMD, ``updwf`` applies an optimiser step to ``c0`` before it
checks convergence, so the orbitals left behind after the last iteration
have moved one step past the ones that produced the reported energy and
forces. The patch checks convergence first and updates only unconverged
orbitals, so the saved ``c0`` is the state that ``forcedr`` used.

What clears the orbitals
========================

+----------------------------------+----------------------------------+
| Event                            | Why                              |
+==================================+==================================+
| a cell change                    | the plane-wave basis depends on  |
|                                  | the cell, so the stored          |
|                                  | coefficients no longer fit       |
+----------------------------------+----------------------------------+
| ``cpmdc_session_set_params()``   | a new method invalidates the     |
| or ``cpmdc_session_configure()`` | orbitals                         |
+----------------------------------+----------------------------------+
| evaluating a different session   | the stored copy is process-wide; |
| in the same process              | the other session re-applies its |
|                                  | configuration and starts cold    |
+----------------------------------+----------------------------------+
| any global call                  | each applies its message afresh  |
| (``cpmdc_set_params()``,         |                                  |
| ``cpmdc_energy*()``) or a        |                                  |
| one-shot call                    |                                  |
| (``cpmdc_calculate_result()``)   |                                  |
+----------------------------------+----------------------------------+

A topology change is refused outright rather than treated as cold: the
atom count and the ordered atomic numbers belong to the session.

The RESTART file
================

The warm path never reads or writes ``RESTART`` on the host's behalf.
CPMD still writes ``RESTART.1`` at the end of each SCF, in the working
directory of that moment, and reads one on a cold call when the deck
asks for ``RESTART WAVEFUNCTION``. ``cpmdc`` turns off CPMD's restart of
coordinates, velocities, and the stored geometry before every SCF: the
host owns the geometry, and a ``RESTART COORDINATES`` in the deck would
otherwise overwrite the step's positions.

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
