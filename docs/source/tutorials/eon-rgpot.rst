This tutorial relaxes the water molecule of the
:doc:`first tutorial <first-energy>` with eOn's minimiser, with every
force call going through ``rgpot`` into ``libcpmdc`` inside the
``eonclient`` process. It reuses ``water.params.bin`` from that tutorial
and the OpenCPMD build of ``libcpmdc`` in ``$CPMDC/build-cpmd``.

What connects to what
=====================

eOn's ``RGPOT`` potential loads rgpot's ``CPMDPot``, which opens
``libcpmdc.so`` with ``dlopen`` and keeps one ``CPMDCSession`` for the
life of the client. eOn hands each geometry to ``CPMDPot``, ``CPMDPot``
turns it into a ``ForceInput``, and ``cpmdc`` runs the SCF in the same
process, starting every call after the first from the previous orbitals.
The method comes from a ``CPMDParams`` message file that eOn reads
through ``[RgpotPot] params_path``.

You need an ``eonclient`` built with RGPOT support whose ``[RgpotPot]``
section reads ``params_path``; eOn lists the key under ``RgpotPot`` in
``eon/config.yaml``. Without ``params_path``, eOn builds ``CPMDParams``
from four scalar keys (``functional``, ``cutoff_ry``, ``charge``,
``multiplicity``). A positive cell is then periodic: symmetry 1 when the
edges are equal, symmetry 8 when they are not, and ``CELL VECTORS`` with
no symmetry line when the box is tilted. A zero cell stays ``SYMMETRY 0``
with the Hockney solver. To keep a molecule isolated, the message has to
ask for it. This run uses ``water.params.bin``, which sets
``poissonSolver = "HOCKNEY"``, and stays isolated.

Prepare the run directory
=========================

.. code:: bash

   mkdir -p eon-water
   cd eon-water
   cp ../cpmdc-tutorial/water.params.bin .

Save the starting geometry as ``pos.con``, eOn's structure format. It is
the water molecule of the first tutorial, in the same 10 Angstrom box,
oxygen first:

.. code:: text

   Generated for the cpmdc eOn tutorial

   10.000000 10.000000 10.000000
   90.000000 90.000000 90.000000


   2
   1 2
   15.999000 1.008000
   O
   Coordinates of Component 1
   5.000000 5.000000 5.117300 0 1
   H
   Coordinates of Component 2
   5.000000 5.757200 4.530800 0 2
   5.000000 4.242800 4.530800 0 3

eOn lists each element in one component. ``cpmdc`` also accepts a step
whose atoms of one element are split: it places positions into CPMD's
species blocks and returns forces in the order of the input.

Save the eOn settings as ``config.ini``:

.. code:: ini

   [Main]
   job = minimization

   [Potential]
   potential = rgpot

   [RgpotPot]
   backend = cpmdc
   params_path = water.params.bin

``backend = cpmdc`` selects rgpot's ``CPMDPot``. ``params_path`` is read
relative to the directory ``eonclient`` runs in. ``RGPOT_PARAMS_PATH``
in the environment overrides it.

Run the minimisation
====================

Point rgpot at the library and CPMD at the pseudopotentials, then start
the client:

.. code:: bash

   export CPMDC_LIBRARY=$CPMDC/build-cpmd/libcpmdc.so
   export CPMDC_PSEUDO_DIR=$PWD/../cpmdc-tutorial/Regtests/tests/PP_LIBRARY
   export OMP_NUM_THREADS=1
   eonclient

rgpot tries ``[RgpotPot] engine_path`` first, then ``CPMDC_LIBRARY``,
``RGPOT_CPMDC_ENGINE``, and ``RGPOT_CPMD_ENGINE``, then ``libcpmdc.so``
on the loader path. CPMD writes its output for each force call to the
client's standard output. The first call starts from atomic orbitals;
later calls start from the previous call's orbitals and need fewer SCF
iterations.

When the run ends, the directory holds:

+-----------------+----------------------------------------------------+
| File            | Contents                                           |
+=================+====================================================+
| ``min.con``     | the relaxed structure                              |
+-----------------+----------------------------------------------------+
| ``results.dat`` | the termination status, the potential name, the    |
|                 | number of force calls, and the final energy in eV  |
+-----------------+----------------------------------------------------+

eOn works in eV and Angstrom; rgpot converts from the Hartree and
Hartree/Bohr values that ``cpmdc`` returns.

Change the method
=================

Edit ``water.params.txt``, encode it again, and rerun; eOn needs no
other change:

.. code:: bash

   capnp encode "$CPMDC/schema/Potentials.capnp" CPMDParams \
     < ../cpmdc-tutorial/water.params.txt > water.params.bin

To see what CPMD receives, set ``CPMDC_DECK_OUT=$PWD/deck.inp`` before
``eonclient`` (see :doc:`debugging a deck <../howto/debug-deck>`). For a
periodic system or a larger cluster, write the message as the
:doc:`CPMDParams how-to <../howto/write-cpmdparams>` shows, and choose
the optimiser for warm starts with the
:doc:`optimiser how-to <../howto/wavefunction-optimiser>`.

Run on several ranks
====================

Start ``eonclient`` under ``mpirun`` to let CPMD spread each SCF over
the ranks:

.. code:: bash

   mpirun -np 4 eonclient

Every rank runs the whole client. Only CPMD's parent rank holds the
calculator's result, so the client must take each result from that rank
and must finalize MPI at exit; check that your eOn and rgpot builds do
both, as :doc:`running under mpirun <../howto/mpi>` explains, before
trusting a multi-rank run.
