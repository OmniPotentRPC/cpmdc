The schema is ``schema/Potentials.capnp``, a vendored copy of the shared
``potentials-schema`` release pinned in
``subprojects/potentials-schema.wrap``. This page lists the fields that
decide what CPMD computes. The :doc:`option mapping <cpmd-options>`
lists every typed field with its feature ID and deck spelling.

``CPMDParams``
==============

+----------------------+----------------------------+-----------------------+
| Field                | Type, default              | Role                  |
+======================+============================+=======================+
| ``functional``       | ``Text``, ``"BLYP"``       | functional when no    |
|                      |                            | ``dft`` section sets  |
|                      |                            | one                   |
+----------------------+----------------------------+-----------------------+
| ``cutOffRy``         | ``Float64``, 70.0          | plane-wave cutoff in  |
|                      |                            | Rydberg when no       |
|                      |                            | ``system`` section    |
|                      |                            | sets one              |
+----------------------+----------------------------+-----------------------+
| ``charge``           | ``Int32``, 0               | system charge when    |
|                      |                            | the ``system``        |
|                      |                            | section does not set  |
|                      |                            | one                   |
+----------------------+----------------------------+-----------------------+
| ``multiplicity``     | ``Int32``, 1               | 2S+1; above 1 selects |
|                      |                            | spin-polarized        |
|                      |                            | defaults              |
+----------------------+----------------------------+-----------------------+
| ``task``             | ``Text``, ``"gradient"``   | host hint:            |
|                      |                            | ``energy``,           |
|                      |                            | ``gradient``, ``md``, |
|                      |                            | ``optimize``          |
+----------------------+----------------------------+-----------------------+
| ``title``            | ``Text``                   | comment header of the |
|                      |                            | rendered deck         |
+----------------------+----------------------------+-----------------------+
| ``memoryMb``         | ``UInt32``, 0              | host memory hint      |
+----------------------+----------------------------+-----------------------+
| ``scratchDir``,      | ``Text``                   | ``FILEPATH`` in       |
| ``permanentDir``     |                            | ``&CPMD``;            |
|                      |                            | ``permanentDir`` wins |
+----------------------+----------------------------+-----------------------+
| ``cpmdRoot``         | ``Text``                   | OpenCPMD tree hint    |
|                      |                            | for hosts             |
+----------------------+----------------------------+-----------------------+
| ``enginePath``       | ``Text``                   | library path hint for |
|                      |                            | hosts such as         |
|                      |                            | ``rgpot``             |
+----------------------+----------------------------+-----------------------+
| ``inputBlocks``      | ``List(Text)``             | literal               |
|                      |                            | ``&SECTION ... &END`` |
|                      |                            | blocks placed before  |
|                      |                            | the rendered sections |
+----------------------+----------------------------+-----------------------+
| ``inputSections``    | ``List(CPMDInputSection)`` | structured sections,  |
|                      |                            | rendered in list      |
|                      |                            | order                 |
+----------------------+----------------------------+-----------------------+

A section missing from ``inputSections`` gets a default: ``&CPMD`` with
``OPTIMIZE WAVEFUNCTION`` and ``CONVERGENCE ORBITALS`` 1e-6, ``&SYSTEM``
with ``SYMMETRY 0``, ``ANGSTROM``, and the top-level cutoff, charge, and
multiplicity, and ``&DFT`` with the top-level functional.

``CPMDInputSection``
====================

A union; each list entry sets exactly one arm.

+----------------------+------------------------------+-----------------------+
| Arm                  | Payload                      | Renders               |
+======================+==============================+=======================+
| ``cpmd``             | ``CPMDCpmdSection``          | ``&CPMD``             |
+----------------------+------------------------------+-----------------------+
| ``system``           | ``CPMDSystemSection``        | ``&SYSTEM``           |
+----------------------+------------------------------+-----------------------+
| ``dft``              | ``CPMDDftSection``           | ``&DFT``              |
+----------------------+------------------------------+-----------------------+
| ``atoms``            | ``CPMDAtomsSection``         | ``&ATOMS`` without    |
|                      |                              | coordinates           |
+----------------------+------------------------------+-----------------------+
| ``vdwParams``,       | typed long-tail sections     | the matching section  |
| ``propParams``,      |                              |                       |
| ``linresParams``,    |                              |                       |
| ``pimdParams``,      |                              |                       |
| ``pathParams``,      |                              |                       |
| ``tddftParams``,     |                              |                       |
| ``respParams``,      |                              |                       |
| ``exteParams``,      |                              |                       |
| ``vectorsParams``    |                              |                       |
+----------------------+------------------------------+-----------------------+
| ``atom``, ``basis``, | ``CPMDDirectiveSection``:    | the matching section, |
| ``clas``, ``eam``,   | ``directives`` and nested    | one line per          |
| ``exte``,            | ``subsections``              | directive             |
| ``hardness``,        |                              |                       |
| ``info``,            |                              |                       |
| ``linres``,          |                              |                       |
| ``molstates``,       |                              |                       |
| ``mts``, ``nlcc``,   |                              |                       |
| ``path``, ``pimd``,  |                              |                       |
| ``potential``,       |                              |                       |
| ``prop``,            |                              |                       |
| ``ptddft``,          |                              |                       |
| ``resp``, ``tddft``, |                              |                       |
| ``vdw``,             |                              |                       |
| ``vectors``,         |                              |                       |
| ``wavefunction``     |                              |                       |
+----------------------+------------------------------+-----------------------+
| ``set``              | ``CPMDSetDirective``:        | one keyword merged    |
|                      | ``key = "SECTION.KEYWORD"``, | into a named section  |
|                      | ``value``                    |                       |
+----------------------+------------------------------+-----------------------+
| ``generic``          | ``CPMDGenericSection``:      | a section by name     |
|                      | ``name``, ``directives``     |                       |
+----------------------+------------------------------+-----------------------+
| ``raw``              | ``Text``                     | a full                |
|                      |                              | ``&SECTION ... &END`` |
|                      |                              | block, verbatim       |
+----------------------+------------------------------+-----------------------+

``CPMDDirective``
=================

+-------------+--------------------------------------------------------+
| Field       | Meaning                                                |
+=============+========================================================+
| ``keyword`` | the line, for example ``"PCG MINIMIZE"``               |
+-------------+--------------------------------------------------------+
| ``args``    | ``List(Text)``; each entry becomes one indented line   |
|             | below the keyword                                      |
+-------------+--------------------------------------------------------+

Every typed section has a ``directives`` list for keywords without a
typed field.

``CPMDCpmdSection``, the fields a force call uses
=================================================

+----------------------------+------------+----------------------------+
| Field                      | Default    | Deck                       |
+============================+============+============================+
| ``optimizeWavefunction``   | ``true``   | ``OPTIMIZE WAVEFUNCTION``  |
+----------------------------+------------+----------------------------+
| ``convergenceOrbitals``    | 1e-6       | ``CONVERGENCE ORBITALS``   |
+----------------------------+------------+----------------------------+
| ``maxIter``                | 0, omitted | ``MAXITER``; the OpenCPMD  |
|                            |            | path inserts 40 when       |
|                            |            | absent                     |
+----------------------------+------------+----------------------------+
| ``maxStep``                | 0, omitted | ``MAXSTEP``, the geometry  |
|                            |            | step limit of CPMD's own   |
|                            |            | drivers                    |
+----------------------------+------------+----------------------------+
| ``odiis``,                 | off        | ``ODIIS`` and its vector   |
| ``odiisVectors``,          |            | count                      |
| ``odiisOptions``           |            |                            |
+----------------------------+------------+----------------------------+
| ``pcg``                    | off        | ``PCG``; ``PCG MINIMIZE``  |
|                            |            | goes in ``directives``     |
+----------------------------+------------+----------------------------+
| ``diis``                   | off        | ``DIIS``                   |
+----------------------------+------------+----------------------------+
| ``centerMoleculeOff``,     | off        | ``CENTER MOLECULE OFF`` /  |
| ``centerMoleculeOn``       |            | ``ON``                     |
+----------------------------+------------+----------------------------+
| ``restartWavefunction``    | off        | ``RESTART WAVEFUNCTION``;  |
|                            |            | read on the cold first     |
|                            |            | call only                  |
+----------------------------+------------+----------------------------+
| ``isolatedMolecule``       | off        | ``ISOLATED MOLECULE``      |
+----------------------------+------------+----------------------------+

The OpenCPMD path ignores restarted coordinates, velocities, and
geometry: the host owns them.

``CPMDSystemSection``, the fields that set the cell
===================================================

+-----------------------------+----------+-----------------------------------------------------+
| Field                       | Default  | Deck                                                |
+=============================+==========+=====================================================+
| ``symmetry``                | 0        | ``SYMMETRY``; 0 is isolated                         |
+-----------------------------+----------+-----------------------------------------------------+
| ``angstrom``                | ``true`` | ``ANGSTROM``, the unit of ``CELL`` and ``&ATOMS``   |
+-----------------------------+----------+-----------------------------------------------------+
| ``cell``                    | empty    | ``CELL``: six values                                |
|                             |          | ``a, b/a, c/a, cos(alpha), cos(beta), cos(gamma)``, |
|                             |          | or nine with ``cellVectors``                        |
+-----------------------------+----------+-----------------------------------------------------+
| ``cellAbsolute``,           | off      | ``CELL ABSOLUTE``, ``DEGREE``, ``VECTORS``          |
| ``cellDegree``,             |          |                                                     |
| ``cellVectors``             |          |                                                     |
+-----------------------------+----------+-----------------------------------------------------+
| ``cutOffRy``                | 70.0     | ``CUTOFF``                                          |
+-----------------------------+----------+-----------------------------------------------------+
| ``poissonSolver``,          | empty    | ``POISSON SOLVER``; ``SYMMETRY 0`` without one adds |
| ``poissonParameter``        |          | ``HOCKNEY``                                         |
+-----------------------------+----------+-----------------------------------------------------+
| ``charge``,                 | 0, 1     | ``CHARGE``, ``MULTIPLICITY``                        |
| ``multiplicity``            |          |                                                     |
+-----------------------------+----------+-----------------------------------------------------+
| ``kpointsMonkhorstPack``    | empty    | ``KPOINTS MONKHORST-PACK``                          |
+-----------------------------+----------+-----------------------------------------------------+

``CPMDDftSection``
==================

========================= ========== =========================
Field                     Default    Deck
========================= ========== =========================
``functional``            ``"BLYP"`` ``FUNCTIONAL``
``lsd``                   ``false``  ``LSD``
``gcCutoff``              0, omitted ``GC-CUTOFF``
``oldCode``, ``newCode``  off        ``OLDCODE`` / ``NEWCODE``
``hfx``, ``hfxScreening`` off        hybrid exchange controls
========================= ========== =========================

``CPMDAtomsSection`` and ``CPMDAtomsPseudopotential``
=====================================================

+---------------------------+-------------+---------------------------+
| Field                     | Default     | Deck                      |
+===========================+=============+===========================+
| ``element``               |             | element symbol, such as   |
|                           |             | ``"Si"``                  |
+---------------------------+-------------+---------------------------+
| ``path``                  |             | the ``*file`` line, a     |
|                           |             | library name or a path    |
+---------------------------+-------------+---------------------------+
| ``lmax``, ``loc``,        | -1, omitted | ``LMAX=``, ``LOC=``,      |
| ``skip``                  |             | ``SKIP=``; 0 is ``S``, 1  |
|                           |             | ``P``, 2 ``D``, 3 ``F``   |
+---------------------------+-------------+---------------------------+
| ``kleinmanBylander``      | ``false``   | ``KLEINMAN-BYLANDER`` on  |
|                           |             | the ``*file`` line        |
+---------------------------+-------------+---------------------------+
| ``raggio``                | 0, omitted  | ``RAGGIO=``               |
+---------------------------+-------------+---------------------------+
| ``nonlinearCore``         | ``false``   | ``NLCC``                  |
+---------------------------+-------------+---------------------------+

The ``atoms`` section also carries ``constraints``, ``isotopes``,
``velocities``, ``dummyAtoms``, ``changeBonds``, ``generate``, and
``directives``. Coordinates never go here. On the OpenCPMD path, the
``&ATOMS`` block CPMD reads is rebuilt from each step's geometry.
Entries in ``pseudopotentials`` supply the file, ``LMAX``, ``LOC``, and
``KLEINMAN-BYLANDER`` for the elements they name. Elements left out of
that list use a built-in table for H, C, N, O, Si, and Ge; see
:doc:`writing a CPMDParams message <../howto/write-cpmdparams>`.

``ForceInput``
==============

+----------------+----------------+------------------------------------+
| Field          | Default        | Meaning                            |
+================+================+====================================+
| ``pos``        |                | ``3 * natoms`` coordinates in      |
|                |                | ``lengthUnit``                     |
+----------------+----------------+------------------------------------+
| ``atmnrs``     |                | atomic numbers                     |
+----------------+----------------+------------------------------------+
| ``box``        | empty          | nine values, row-major cell        |
|                |                | vectors in ``lengthUnit``          |
+----------------+----------------+------------------------------------+
| ``lengthUnit`` | ``"angstrom"`` | unit of ``pos`` and ``box``        |
+----------------+----------------+------------------------------------+
| ``energyUnit`` | ``"eV"``       | unit of the ``PotentialResult``    |
|                |                | energy                             |
+----------------+----------------+------------------------------------+

``PotentialConfig``
===================

A union of backend arms plus a ``common`` overlay
(``CommonMethodSpec``). ``cpmdc`` accepts the ``cpmd`` arm, which holds
a ``CPMDParams``, or the overlay alone; see the
:doc:`C ABI reference <c-abi>` for the lowering rules.
