A force call on a geometry close to the previous one should not start
its SCF from scratch. The converged orbitals of the last call, ``c0`` in
OpenCPMD, are a far better starting point than atomic orbitals, and
keeping them in memory is the main reason to run CPMD inside the host.
This page follows ``c0`` from call to call.

Cold and warm calls
===================

The process keeps the converged orbitals, the cell, and the method that
produced them. Each session keeps a counter of its own successful calls.

+------+-------------------------------+-------------------------------+
| Call | Condition                     | What runs                     |
+======+===============================+===============================+
| cold | first call of the process, or | the full CPMD setup on the    |
|      | a changed functional, deck,   | composed deck, then the SCF   |
|      | cell, cutoff, charge,         | from CPMD's initial guess, or |
|      | multiplicity, atom count, or  | from ``RESTART`` when the     |
|      | elemental composition         | deck has                      |
|      |                               | ``RESTART WAVEFUNCTION``      |
+------+-------------------------------+-------------------------------+
| warm | a previous success in this    | new positions into ``tau0``,  |
|      | process kept the same         | ``phfac``, then the SCF from  |
|      | functional, deck, cell,       | the stored ``c0``             |
|      | cutoff, charge, multiplicity, |                               |
|      | and elemental composition,    |                               |
|      | including a new session or a  |                               |
|      | reordering of those atomic    |                               |
|      | numbers                       |                               |
+------+-------------------------------+-------------------------------+

A cell counts as unchanged when every component of ``ForceInput.box``
matches the stored cell within 1e-8 Angstrom, and when both calls agree
on having a box at all. A new session has no cell of its own yet, so
that comparison uses the cell stored with the orbitals. The cutoff
matches within 1e-8 Rydberg. The elemental composition is the count of
each atomic number. Reordering the same atoms keeps that composition.
The charge and the multiplicity belong to the same match. A failed warm
call leaves the counter where it was, so the next call is warm again; a
failed cold call clears the stored results. An SCF that does not
converge does not replace the stored orbitals. When ODIIS exhausts
``MAXITER``, the same call continues once with ``PCG MINIMIZE``. That
pass restores the previous converged ``c0`` when a converged copy is
already stored, and it still counts as warm. A ``stopgm`` during the
call is not that failure: the wavefunction, forces, and module state of
the call are undefined, ``cpmdc`` drops the stored orbitals and the warm
cell, and the next call runs setup again. With more than one MPI rank,
``cpmdc`` aborts the other ranks.

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
converges. On a later call, ``rwfopt`` copies the saved ``c0`` in before
``initrun``. ``rinitwf`` then returns after the phase factors, so it
does not build a guess, and the density ``initrun`` builds comes from
those orbitals. A restart read inside ``initrun`` can replace ``c0``, so
``rwfopt`` copies the saved orbitals in again after ``initrun`` when the
shape still matches. A different shape calls ``stopgm`` and names
``embed_reset_warm_orbitals``. ``cpmdc`` makes that call when the
functional, the deck, the cutoff, the cell, the charge, the
multiplicity, the atom count, or the elemental composition changes,
before ``rwfopt`` runs. The SCF then runs to the orbital convergence
threshold with the full ``MAXITER`` budget of the deck; a warm call is a
complete SCF from a better start, not a shortened one.

``opencpmd_converged_state.patch`` makes the saved copy trustworthy. In
stock OpenCPMD, ``updwf`` applies an optimiser step to ``c0`` before it
checks convergence, so the orbitals left behind after the last iteration
have moved one step past the ones that produced the reported energy and
forces. The patch checks convergence first and updates only unconverged
orbitals, so the saved ``c0`` is the state that ``forcedr`` used.
Steepest descent with ``iproj <= 1`` sets ``gemax`` above the threshold
before that check: those steps do not converge on ``gemax``.

Orbitals per key
================

One calculator often evaluates several calculations in turn: the images
of a band on fewer calculators than images, or the beads of a ring
polymer. With one stored copy, each SCF starts from the orbitals of
whichever calculation ran last on that calculator, which for a ring is a
neighbouring bead rather than the same bead one step earlier.
``cpmdc_session_select_orbitals()`` names the calculation the next
evaluations belong to. ``opencpmd_embed_rwfopt.patch`` keeps a table of
stored ``c0`` copies, one per key, beside the store ``rwfopt`` restores
from. Selecting another key parks the converged ``c0`` under the key
being left and moves the copy parked under the new key into the store. A
partial copy from an unconverged pass is not parked. Every rank of a
calculator selects the same key, and each rank parks its own slice of
``c0``. A basis change frees the table with the store, since every copy
in it has the shape of the old basis. A session that never selects a key
keeps one stored copy.

Distance limit
==============

Every stored and parked ``c0`` carries the geometry it converged at, the
Angstrom positions the caller passed. Before each SCF,
``embed_choose_warm_start`` measures the root-mean-square distance per
atom from the new geometry to each of them, without alignment. The
nearest one starts the SCF. It is either the store, which holds the
selected key's own copy or, when that key has none, the copy of the key
last left, or another key's parked copy, which is copied so the donor
keeps its own. When that distance exceeds ``CPMDC_WARM_RMSD`` (default
0.15 Angstrom), the SCF starts from the deck's guess instead, and the
stored copy of that call is dropped. A negative limit starts from the
store at any distance. Periodic systems use minimage's Euclidean closest
lattice image in the native cell, including skew cells. Isolated cluster
and Hockney/Tuckerman systems use absolute Cartesian distances even when
a cell is supplied. Coordinates are compared in species order without
rotation or alignment. A cell or boundary-model change invalidates every
stored orbital state, including when the distance limit is negative.

The default comes from the Si3N4 BLYP 70 Ry band (7 atoms, 14 Angstrom
cell, ODIIS, 4 ranks per calculation). One process makes two calls. The
second sits a set distance along the band from the first and starts from
the first call's orbitals or from the deck's atomic guess. Five starting
images, means over them per call:

+----------------+----------------+----------------+----------------+
| RMSD           | SCF steps,     | SCF steps,     | wall saved by  |
| (Angstrom)     | warm           | guess          | the warm start |
|                |                |                | (s)            |
+================+================+================+================+
| 0.047          | 24.0           | 32.4           |                |
+----------------+----------------+----------------+----------------+
| 0.116          | 29.2           | 33.2           | 4.5            |
+----------------+----------------+----------------+----------------+
| 0.149          | 32.4           | 33.0           | 2.1            |
+----------------+----------------+----------------+----------------+
| 0.167          | 33.8           | 33.4           | 1.2            |
+----------------+----------------+----------------+----------------+
| 0.186          | 35.4           | 33.0           | -0.4           |
+----------------+----------------+----------------+----------------+
| 0.233          | 37.6           | 32.8           | -2.4           |
+----------------+----------------+----------------+----------------+
| 0.46 (3)       | 47.7           | 37.7           |                |
+----------------+----------------+----------------+----------------+
| 0.9 to 3.1     | 36 to 283      | 19 to 52       |                |
+----------------+----------------+----------------+----------------+

Next to an endpoint the guess is already close (19 to 26 steps): a warm
start there gains at most 3 steps and loses from 0.12 Angstrom on.
Mid-path images gain up to 27 steps at 0.05 Angstrom. The 0.46 row has
three starting images; the last row covers the warm starts from images
0, 4 and 8 to every image that far away. The deck's guess, not the
simple guess (``inwfun 3``) a process that holds orbitals otherwise
takes, is the start past the limit: the simple guess took 65 to 205
steps on the same images where the deck's took 37 to 52. A band step
capped at 0.1 Angstrom per atom stays inside the limit, so an image that
returns to its own key is warm. ``CPMDC_SCF_STEPS=1`` prints the start
of every call.

What clears the orbitals
========================

+----------------------------------+----------------------------------+
| Event                            | Why                              |
+==================================+==================================+
| a cell change                    | the plane-wave basis depends on  |
|                                  | the cell, so the stored          |
|                                  | coefficients no longer fit       |
+----------------------------------+----------------------------------+
| a cutoff change                  | the plane-wave basis changes, so |
|                                  | the stored coefficients no       |
|                                  | longer fit                       |
+----------------------------------+----------------------------------+
| a change of atom count or        | those inputs set the number of   |
| elemental composition, or of     | states, and a restored ``c0`` of |
| charge or multiplicity           | another shape calls ``stopgm``   |
+----------------------------------+----------------------------------+
| ``cpmdc_session_set_params()``   | the next evaluation drops the    |
| or ``cpmdc_session_configure()`` | orbitals; the call itself keeps  |
| when the functional, cutoff,     | them                             |
| charge, multiplicity, or deck    |                                  |
| differs                          |                                  |
+----------------------------------+----------------------------------+
| a new session with the same      | reuses the stored orbitals and   |
| functional, deck, cutoff, cell,  | skips setup                      |
| charge, multiplicity, and        |                                  |
| elemental composition            |                                  |
+----------------------------------+----------------------------------+
| a new session whose functional,  | the next evaluation drops the    |
| deck, cutoff, cell, charge,      | stored orbitals and runs setup   |
| multiplicity, or elemental       |                                  |
| composition differs              |                                  |
+----------------------------------+----------------------------------+
| a one-shot or global call        | that call runs setup; a matching |
| (``cpmdc_set_params()``,         | basis reuses the orbitals        |
| ``cpmdc_energy*()``,             |                                  |
| ``cpmdc_calculate_result()``)    |                                  |
| whose basis differs from the     |                                  |
| stored one                       |                                  |
+----------------------------------+----------------------------------+
| ``stopgm`` returns through       | the call continued past a failed |
| ``cpmd_stopgm_hook``             | check, so the stored orbitals    |
|                                  | are not a result                 |
+----------------------------------+----------------------------------+

A topology change is refused outright rather than treated as cold: the
atom count and the elemental composition belong to the session. A
reordering of those atomic numbers stays on the session. A new session
with that composition keeps the stored orbitals. An evaluation whose
positions match the previous successful call on that session, including
that reordering, returns the stored energy and forces.

An evaluation that needs a second setup, for a basis that differs from
the first or after a ``stopgm``, tears the first one down in the same
process. ``cpmdc_embed_teardown`` frees the module arrays,
communicators, and first-call state of that setup
(``tools/opencpmd_embed_teardown.patch``), and the new setup gives the
energies of a fresh process. ``cpmdc_embed_reset_state`` tears down too.
``tistart`` runs once per process: OpenCPMD executes ``STOP`` on a
second call.

The RESTART file
================

A library force call writes no ``RESTART.1``, ``LATEST``, ``GEOMETRY``,
or ``GEOMETRY.xyz``. The orbitals for the next call stay in memory.
``cpmd.x`` writes those files. A cold library call reads a ``RESTART.1``
only when the deck asks for ``RESTART WAVEFUNCTION``. ``cpmdc`` turns
off CPMD's restart of coordinates, velocities, and the stored geometry
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
