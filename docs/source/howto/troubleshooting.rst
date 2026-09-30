Each entry names the symptom, the cause in the code, and the fix. Set
``CPMDC_DECK_OUT`` first when the symptom involves energies (see
:doc:`debugging a deck <debug-deck>`).

Energies
========

The energy is off by about 2 eV from a hand-written deck
--------------------------------------------------------

**Cause:** ``KLEINMAN-BYLANDER`` belongs on the pseudopotential file
line, ``*Si_MT_BLYP.psp KLEINMAN-BYLANDER``. On the ``LMAX=`` line below
it CPMD ignores the keyword without a warning and integrates the
nonlocal projectors by Gauss-Hermite quadrature. On a 7-atom Si3N4
cluster that alone moved the energy from -1396.2695 eV to -1398.3638 eV,
2.09 eV.

**Fix:** set ``kleinmanBylander = true`` on each
``CPMDAtomsPseudopotential`` entry, which renders the keyword on the
file line. In ``inputBlocks`` or ``raw`` text, write it on the ``*file``
line yourself.

The system is treated as isolated
---------------------------------

**Cause:** a message without a ``system`` section renders ``SYMMETRY 0``
with no ``CELL``, and ``SYMMETRY 0`` without a Poisson solver adds
``POISSON SOLVER HOCKNEY``. That is the isolated-cluster setup, with the
cell taken from ``ForceInput.box`` (or a 12 Angstrom cube without a
box). A host that fills only the top-level scalars, ``functional``,
``cutOffRy``, ``charge``, and ``multiplicity``, gets this deck.

**Fix:** write a ``system`` section with ``symmetry`` and ``cell``, and
``poissonSolver`` when the system is isolated on purpose. In eOn, pass
the whole message through ``[RgpotPot] params_path`` instead of the
scalar keys (see :doc:`the eOn tutorial <../tutorials/eon-rgpot>`).

The energy changes when only the cell changes
---------------------------------------------

**Cause:** a cell that differs from the previous call's cell makes the
next call cold: ``cpmdc`` clears the stored orbitals and runs CPMD's
full setup again. The energy is then that of a fresh SCF in the new
cell.

**Fix:** none needed when the cell change is intended. Keep
``ForceInput.box`` identical between steps when the cell is fixed; a
difference above 1e-8 Angstrom in any component counts as a change.

Failed calls
============

``orbitals not converged within MAXITER, or no energy``
-------------------------------------------------------

**Cause:** the SCF stopped at ``MAXITER`` before reaching the orbital
convergence threshold. OpenCPMD computes ionic forces only for converged
orbitals, so the forces of such a call are zero, which a caller would
read as a stationary point; ``cpmdc`` reports the call as failed
instead, with ``ok == 0``. Without ``maxIter`` in the message the
OpenCPMD path inserts ``MAXITER 40``.

The same message has three other causes, all visible before the SCF
starts:

+----------------------------------+----------------------------------+
| Check                            | Cause                            |
+==================================+==================================+
| ``CPMDC_PSEUDO_DIR`` and         | the first call returns before    |
| ``CPMD_PP_LIBRARY_PATH`` are     | CPMD's setup                     |
| both unset or name a missing     |                                  |
| directory                        |                                  |
+----------------------------------+----------------------------------+
| an element absent from the       | no pseudopotential in the        |
| message and outside H, C,        | message or the built-in          |
| N, O, Si, and Ge                 | table; the error names the       |
|                                  | atomic number                    |
+----------------------------------+----------------------------------+
| atoms of one element are not     | positions are copied species by  |
| contiguous in ``ForceInput``, or | species into CPMD's ``tau0``,    |
| the species order differs from   | and a mismatch stops the copy    |
| the deck's ``&ATOMS`` block      |                                  |
+----------------------------------+----------------------------------+

**Fix:** raise ``maxIter``, or switch the optimiser (see
:doc:`choosing the optimiser <wavefunction-optimiser>`); otherwise fix
the item from the table.

``CPMD stopgm during embed SCF``
--------------------------------

**Cause:** CPMD called ``stopgm``, its fatal-error routine, inside the
call. With ``opencpmd_stopgm_return.patch`` in the archive, ``stopgm``
calls ``cpmd_stopgm_hook``. While an evaluation is armed the catch records
the stop code and returns 1, so ``stopgm`` returns to its caller, and
``cpmdc`` returns ``ok == 0`` instead of ending the host. CPMD writes its reason
to a ``LocalError-*.log`` file in the working directory of that moment.
During the first call that is the pseudopotential directory, because
``cpmdc`` changes into it while CPMD reads the pseudopotentials.

**Fix:** read the ``LocalError`` file in ``CPMDC_PSEUDO_DIR`` or in the
host's working directory. The stopped call has no result. ``cpmdc`` drops
the stored orbitals and the warm cell, and the next call on the same
session sets CPMD up again.

``topology change requires a new session``
------------------------------------------

**Cause:** the first successful step of a session fixes the atom count
and the ordered atomic numbers. A later step with a different count or
order is refused before CPMD runs.

**Fix:** create a new ``CPMDCSession`` for the new composition.

``PotentialResult buffer too small``
------------------------------------

**Cause:** the output buffer is smaller than the serialized result.
``cpmdc`` writes the required size to the size argument and does not
evaluate.

**Fix:** allocate ``cpmdc_potential_result_size_for_force_input()``
bytes.

``CPMD embed not available``
----------------------------

**Cause:** the library was finalized with ``cpmdc_finalize()``, or the
session could not apply its configuration.

**Fix:** call ``cpmdc_finalize()`` only at shutdown; check the message
with ``CPMDC_DECK_OUT`` for a rendering failure.

Forces
======

Forces are all zero with a finite energy
----------------------------------------

**Cause:** the OpenCPMD archive lacks ``opencpmd_embed_rwfopt.patch``, so
``fion`` is deallocated before ``cpmdc`` reads it.

**Fix:** apply the patches in ``tools/`` and rebuild the archive (see
:doc:`building the archive <opencpmd-archive>`). A symmetric cluster at
a stationary geometry also gives zero forces; displace an atom to tell
the two apart.

Ranks disagree about forces under mpirun
----------------------------------------

**Cause:** only CPMD's parent rank holds the calculator's forces.

**Fix:** broadcast the parent rank's result (see
:doc:`running under mpirun <mpi>`).

mpirun reports an abnormal termination at exit
----------------------------------------------

**Cause:** CPMD initialised MPI and nothing finalized it; Open MPI 5
kills the ranks still running when one exits without ``MPI_Finalize``.

**Fix:** call ``MPI_Finalize`` on every rank before exit.

Build and link
==============

``with_cpmd requires DIR/lib/libcpmd.a``
----------------------------------------

**Cause:** ``cpmd_root`` does not point at a completed OpenCPMD build
tree.

**Fix:** pass the ``-DEST`` directory of the OpenCPMD build, which holds
``lib/``, ``obj/``, and ``src/``.

Link fails with ``R_X86_64_PC32`` against ``libcpmd.a``
-------------------------------------------------------

**Cause:** OpenCPMD was compiled without ``-fPIC``.

**Fix:** add ``-fPIC`` to ``FFLAGS`` and ``CFLAGS`` of the OpenCPMD
configuration and rebuild the archive.

``cpmd_embed_c_api.F90`` fails to compile on ``embed_set_warm_orbitals``
---------------------------------------------------------------------

**Cause:** the archive's module files come from an unpatched tree.

**Fix:** apply ``opencpmd_embed_rwfopt.patch``, then rebuild with
``tools/rebuild_opencpmd_embed.sh``.

CPMD cannot find ``O_MT_BLYP.psp``
----------------------------------

**Cause:** the pseudopotential directory is not set. CPMD takes
``argv[2]`` as its pseudopotential library whenever the process has more
than one argument, so a host with command-line arguments cannot rely on
``CPMD_PP_LIBRARY_PATH``. ``cpmdc`` works around that by changing into
the directory for the first call.

**Fix:** ``export CPMDC_PSEUDO_DIR=/path/to/pseudopotentials``.

CPMD files appear in the pseudopotential directory
--------------------------------------------------

**Cause:** the first call of a session runs with its working directory
in ``CPMDC_PSEUDO_DIR``, so the ``RESTART.1``, ``LATEST``, ``GEOMETRY``,
and ``GEOMETRY.xyz`` that CPMD writes at the end of that SCF land there.
Later calls run in the host's working directory.

**Fix:** give each run a writable copy of the pseudopotential directory,
or link the needed files into a per-run directory and point
``CPMDC_PSEUDO_DIR`` at it.

A long method deck loses its end
--------------------------------

**Cause:** the OpenCPMD path stores the rendered method deck in a
4096-character buffer; text past that is dropped.

**Fix:** move long literal blocks into typed fields, or check the deck
size in ``CPMDC_DECK_OUT`` (``wc -c``).
