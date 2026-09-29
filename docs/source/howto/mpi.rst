Goal
====

Run a host that embeds ``cpmdc`` on several MPI (Message Passing
Interface) ranks, so that CPMD parallelises each self-consistent field
(SCF) calculation over the ranks, and keep every rank's view of the
result consistent. This page needs an OpenCPMD build of ``cpmdc``; the
default build runs the reference evaluator on each rank independently.

The rules
=========

#. Start every rank of the host under ``mpirun``. CPMD initialises MPI
   on the first evaluation when the host has not.
#. Every rank makes the same ``cpmdc`` calls with the same messages in
   the same order. An evaluation is collective over the ranks of its
   calculator.
#. Take the energy and forces from the first rank of each calculator,
   CPMD's parent rank, and broadcast them to the other ranks. The other
   ranks return from the call with values of their own. In a 4-rank
   water single point all four printed the parent's energy and forces,
   but in a 4-rank dimer search, where every rank also ran the host's
   optimizer, rank 1 took 3 steps on different curvatures while rank 0
   took 5, and only rank 0's results matched an independent ``cpmd.x``
   run.
#. Keep one ``CPMDCSession`` per calculator.
#. Call ``MPI_Finalize`` on every rank before the host exits.
   ``cpmdc_finalize()`` does not finalize MPI.

One calculator over all ranks
=============================

Without further setup, CPMD uses ``MPI_COMM_WORLD``, so all ranks form
one calculator:

.. code:: bash

   export CPMDC_PSEUDO_DIR=/path/to/pseudopotentials
   mpirun -np 16 ./host params.bin step.bin

Rank 0 of ``MPI_COMM_WORLD`` is then CPMD's parent rank.

Several calculators in one job
==============================

``cpmdc_bind_calculator(ranks_per_calc)`` splits ``MPI_COMM_WORLD`` into
groups of ``ranks_per_calc`` consecutive ranks. Each group becomes one
CPMD calculator with its own communicator and its own wavefunction,
running the same deck. Call it on every rank, once, before the first
evaluation:

.. code:: c

   #include <cpmdc.h>
   #include <mpi.h>

   int main(int argc, char **argv) {
     MPI_Init(&argc, &argv);
     int group = cpmdc_bind_calculator(4); /* 16 ranks -> 4 calculators */
     if (group < 0) {
       /* world size not a multiple of 4, or no CPMD backend */
       MPI_Abort(MPI_COMM_WORLD, 1);
     }
     int world_rank = 0;
     MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
     int is_parent = (world_rank % 4) == 0;

     /* group 0 evaluates image 0, group 1 image 1, ... */
     /* ... create one session, evaluate, share the parent's result ... */

     cpmdc_finalize();
     MPI_Finalize();
     return 0;
   }

+-------------------------+--------------------------------------------+
| Argument or result      | Meaning                                    |
+=========================+============================================+
| ``ranks_per_calc > 0``  | ranks per calculator; the world size must  |
|                         | be a multiple of it                        |
+-------------------------+--------------------------------------------+
| ``ranks_per_calc <= 0`` | one calculator over the whole world        |
+-------------------------+--------------------------------------------+
| return value ``>= 0``   | index of this rank's calculator,           |
|                         | ``world_rank / ranks_per_calc``            |
+-------------------------+--------------------------------------------+
| return value ``-1``     | the split was refused, or the library has  |
|                         | no CPMD backend (the default build and the |
|                         | stub)                                      |
+-------------------------+--------------------------------------------+

A second call returns the same index without splitting again. The
group's first rank, world rank ``group * ranks_per_calc``, is its parent
rank. The split needs ``opencpmd_mp_comm_set.patch`` in the OpenCPMD
archive: it keeps CPMD's ``mp_start`` from resetting the communicator
back to ``MPI_COMM_WORLD``.

Share the parent's result
=========================

The host owns the broadcast. With one calculator over the whole world,
broadcast from world rank 0:

.. code:: c

   double buf[1 + 3 * N_ATOMS];
   if (world_rank == 0) {
     buf[0] = result.energy_h;
     memcpy(buf + 1, forces, sizeof(double) * 3 * N_ATOMS);
   }
   MPI_Bcast(buf, 1 + 3 * N_ATOMS, MPI_DOUBLE, 0, MPI_COMM_WORLD);

With several calculators, each group's parent broadcasts inside its
group, and a host that needs every group's result on every rank
broadcasts once per group over ``MPI_COMM_WORLD``. An optimizer that
runs on every rank must step from the broadcast values; ranks that step
from their own copies drift apart after the first step.

One session per calculator
==========================

The OpenCPMD module state and the stored orbitals belong to the process,
not to a ``CPMDCSession``. When a process evaluates a second session,
``cpmdc`` applies that session's configuration and clears the stored
orbitals; the next call of either session starts cold. Run independent
calculations, such as the images of a nudged elastic band, on separate
calculators, one session each.

Finalize MPI
============

CPMD calls ``MPI_Init`` during its setup, and nothing in ``cpmdc`` calls
``MPI_Finalize``. Open MPI 5 treats a rank that exits without
``MPI_Finalize`` as an abnormal termination and kills the ranks still
running; ``mpirun`` then reports a rank "exiting improperly" and lists a
missing finalize among the reasons. Call ``MPI_Finalize`` on every rank,
after ``cpmdc_finalize()``, as the example above does. When the host
cannot restructure its exit path, register a handler with ``atexit``
that calls ``MPI_Finalize`` when ``MPI_Initialized`` reports true and
``MPI_Finalized`` reports false.

Threads
=======

CPMD also runs OpenMP threads inside each rank when OpenCPMD was built
with ``-omp``. Set ``OMP_NUM_THREADS`` per rank so that ranks times
threads matches the cores the job holds.
