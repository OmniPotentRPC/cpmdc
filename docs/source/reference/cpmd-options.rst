Contract
========

All method knobs for the OpenCPMD backend are fields of ``CPMDParams``
in ``schema/Potentials.capnp``. Geometry is never part of
``CPMDParams``; it is always a ``ForceInput`` step message (or the raw C
array entry points that take positions and atomic numbers).

Top-level ``CPMDParams`` fields
===============================

+--------------------------+-------------------+------------------------+
| Feature ID               | Field             | Role                   |
+==========================+===================+========================+
| ``params.functional``    | ``functional``    | Default DFT functional |
|                          |                   | (also used when no     |
|                          |                   | ``dft`` section is     |
|                          |                   | present)               |
+--------------------------+-------------------+------------------------+
| ``params.cutOffRy``      | ``cutOffRy``      | Plane-wave cutoff in   |
|                          |                   | Rydberg                |
+--------------------------+-------------------+------------------------+
| ``params.charge``        | ``charge``        | System charge          |
+--------------------------+-------------------+------------------------+
| ``params.multiplicity``  | ``multiplicity``  | 2S+1; values ``> 1``   |
|                          |                   | imply LSD defaults     |
+--------------------------+-------------------+------------------------+
| ``params.task``          | ``task``          | Frontend hint:         |
|                          |                   | ``energy``,            |
|                          |                   | ``gradient``, ``md``,  |
|                          |                   | ``optimize``           |
+--------------------------+-------------------+------------------------+
| ``params.title``         | ``title``         | Optional comment       |
|                          |                   | header in rendered     |
|                          |                   | decks                  |
+--------------------------+-------------------+------------------------+
| ``params.memoryMb``      | ``memoryMb``      | Frontend memory hint   |
+--------------------------+-------------------+------------------------+
| ``params.scratchDir``    | ``scratchDir``    | ``&CPMD FILEPATH``     |
|                          |                   | scratch placement      |
+--------------------------+-------------------+------------------------+
| ``params.permanentDir``  | ``permanentDir``  | ``&CPMD FILEPATH``     |
|                          |                   | permanent placement;   |
|                          |                   | wins over              |
|                          |                   | ``scratchDir``         |
+--------------------------+-------------------+------------------------+
| ``params.cpmdRoot``      | ``cpmdRoot``      | OpenCPMD source/build  |
|                          |                   | tree                   |
+--------------------------+-------------------+------------------------+
| ``params.enginePath``    | ``enginePath``    | Frontend engine        |
|                          |                   | selection hint         |
+--------------------------+-------------------+------------------------+
| ``params.inputBlocks``   | ``inputBlocks``   | Raw ``&SECTION`` text  |
|                          |                   | blocks prepended to    |
|                          |                   | the deck               |
+--------------------------+-------------------+------------------------+
| ``params.inputSections`` | ``inputSections`` | Structured             |
|                          |                   | ``CPMDInputSection``   |
|                          |                   | list                   |
+--------------------------+-------------------+------------------------+

The feature-discovery table exposes these as ``params.<field>`` IDs.
Structured fields use ``params.inputSections.<section>.<field>``, for
example ``params.inputSections.cpmd.maxIter``,
``params.inputSections.dft.hfxScreening``,
``params.inputSections.system.cell``, and
``params.inputSections.atoms.pseudopotentials``.

``CPMDInputSection`` kinds
==========================

+-------------+------------------------+-------------------------------+
| Kind        | Cap'n Proto type       | Deck effect                   |
+=============+========================+===============================+
| ``cpmd``    | ``CPMDCpmdSection``    | ``&CPMD``                     |
|             |                        | (wavefunction/geometry        |
|             |                        | optimization, MD,             |
|             |                        | convergence, iteration        |
|             |                        | limits, timestep, electron    |
|             |                        | mass, …)                      |
+-------------+------------------------+-------------------------------+
| ``system``  | ``CPMDSystemSection``  | ``&SYSTEM`` (symmetry,        |
|             |                        | angstrom, cell, cutoff,       |
|             |                        | charge, …)                    |
+-------------+------------------------+-------------------------------+
| ``dft``     | ``CPMDDftSection``     | ``&DFT`` (functional, LSD, GC |
|             |                        | cutoff, XC driver,            |
|             |                        | hybrid/Hubbard scalars, extra |
|             |                        | directives)                   |
+-------------+------------------------+-------------------------------+
| ``atoms``   | ``CPMDAtomsSection``   | ``&ATOMS`` (pseudopotential   |
|             |                        | lines; coordinates come from  |
|             |                        | ``ForceInput``)               |
+-------------+------------------------+-------------------------------+
| ``generic`` | ``CPMDGenericSection`` | Arbitrary ``&NAME`` …         |
|             |                        | ``&END``                      |
+-------------+------------------------+-------------------------------+
| ``set``     | ``CPMDSetDirective``   | Embed-path logical            |
|             |                        | ``SECTION.KEYWORD`` merged    |
|             |                        | into the section              |
+-------------+------------------------+-------------------------------+
| ``raw``     | ``Text``               | Full section text inserted    |
|             |                        | as-is                         |
+-------------+------------------------+-------------------------------+

Missing core sections receive minimal defaults from the renderer so a
sparse ``CPMDParams`` still produces a valid-looking deck for debugging.
When structured ``dft`` or ``system`` sections provide scalar method
values (functional, cutoff, charge, multiplicity), those section values
are also the effective embed configuration passed through the C ABI. For
geometry deck rendering, ``atoms.pseudopotentials`` groups
``ForceInput`` coordinates by element symbol. Explicit pseudopotential
entries accept standard element symbols and must cover every atomic
number present in the step geometry; otherwise rendering fails. When no
``atoms`` section is supplied, the renderer only provides built-in BLYP
pseudopotential defaults for H and O. An explicit ``atoms`` section may
rely on those H/O defaults while still carrying ``CPMDDirective``
entries for additional ``&ATOMS`` keywords.

Escape hatches
==============

Long-tail CPMD keywords go in ``CPMDDirective`` lists inside a section,
``CPMDSetDirective`` dotted keys, ``CPMDInputSection.raw``, or
``inputBlocks``. ``set.key`` uses ``SECTION.KEYWORD`` form, for example
``CPMD.PRINT FORCES ON`` or ``SYSTEM.POISSON SOLVER HOCKNEY``; non-empty
``set.value`` is emitted on the following indented line. Prefer new
structured fields in the schema over inventing a second config file
format.

Catalog section exposure
========================

The public feature table advertises every OpenCPMD
``inscan('&SECTION')`` name as a ``catalog.section.<NAME>`` feature ID.
Typed sections have dedicated schema arms; every section can still be
rendered through ``set``, ``generic``, ``raw``, or ``inputBlocks`` when
a dedicated field is unnecessary.

+----------------------------------+-------------------+------------------------+
| Feature ID                       | Deck section      | ABI carrier            |
+==================================+===================+========================+
| ``catalog.section.ATOM``         | ``&ATOM``         | ``set``, ``generic``,  |
|                                  |                   | ``raw``, or            |
|                                  |                   | ``inputBlocks``        |
+----------------------------------+-------------------+------------------------+
| ``catalog.section.ATOMS``        | ``&ATOMS``        | typed ``atoms``,       |
|                                  |                   | ``set``, ``generic``,  |
|                                  |                   | ``raw``, or            |
|                                  |                   | ``inputBlocks``        |
+----------------------------------+-------------------+------------------------+
| ``catalog.section.BASIS``        | ``&BASIS``        | ``set``, ``generic``,  |
|                                  |                   | ``raw``, or            |
|                                  |                   | ``inputBlocks``        |
+----------------------------------+-------------------+------------------------+
| ``catalog.section.CLAS``         | ``&CLAS``         | ``set``, ``generic``,  |
|                                  |                   | ``raw``, or            |
|                                  |                   | ``inputBlocks``        |
+----------------------------------+-------------------+------------------------+
| ``catalog.section.CPMD``         | ``&CPMD``         | typed ``cpmd``,        |
|                                  |                   | ``set``, ``generic``,  |
|                                  |                   | ``raw``, or            |
|                                  |                   | ``inputBlocks``        |
+----------------------------------+-------------------+------------------------+
| ``catalog.section.DFT``          | ``&DFT``          | typed ``dft``,         |
|                                  |                   | ``set``, ``generic``,  |
|                                  |                   | ``raw``, or            |
|                                  |                   | ``inputBlocks``        |
+----------------------------------+-------------------+------------------------+
| ``catalog.section.EAM``          | ``&EAM``          | ``set``, ``generic``,  |
|                                  |                   | ``raw``, or            |
|                                  |                   | ``inputBlocks``        |
+----------------------------------+-------------------+------------------------+
| ``catalog.section.EXTE``         | ``&EXTE``         | ``set``, ``generic``,  |
|                                  |                   | ``raw``, or            |
|                                  |                   | ``inputBlocks``        |
+----------------------------------+-------------------+------------------------+
| ``catalog.section.HARDNESS``     | ``&HARDNESS``     | ``set``, ``generic``,  |
|                                  |                   | ``raw``, or            |
|                                  |                   | ``inputBlocks``        |
+----------------------------------+-------------------+------------------------+
| ``catalog.section.INFO``         | ``&INFO``         | ``set``, ``generic``,  |
|                                  |                   | ``raw``, or            |
|                                  |                   | ``inputBlocks``        |
+----------------------------------+-------------------+------------------------+
| ``catalog.section.LINRES``       | ``&LINRES``       | ``set``, ``generic``,  |
|                                  |                   | ``raw``, or            |
|                                  |                   | ``inputBlocks``        |
+----------------------------------+-------------------+------------------------+
| ``catalog.section.MOLSTATES``    | ``&MOLSTATES``    | ``set``, ``generic``,  |
|                                  |                   | ``raw``, or            |
|                                  |                   | ``inputBlocks``        |
+----------------------------------+-------------------+------------------------+
| ``catalog.section.MTS``          | ``&MTS``          | ``set``, ``generic``,  |
|                                  |                   | ``raw``, or            |
|                                  |                   | ``inputBlocks``        |
+----------------------------------+-------------------+------------------------+
| ``catalog.section.NLCC``         | ``&NLCC``         | ``set``, ``generic``,  |
|                                  |                   | ``raw``, or            |
|                                  |                   | ``inputBlocks``        |
+----------------------------------+-------------------+------------------------+
| ``catalog.section.PATH``         | ``&PATH``         | ``set``, ``generic``,  |
|                                  |                   | ``raw``, or            |
|                                  |                   | ``inputBlocks``        |
+----------------------------------+-------------------+------------------------+
| ``catalog.section.PIMD``         | ``&PIMD``         | ``set``, ``generic``,  |
|                                  |                   | ``raw``, or            |
|                                  |                   | ``inputBlocks``        |
+----------------------------------+-------------------+------------------------+
| ``catalog.section.POTENTIAL``    | ``&POTENTIAL``    | ``set``, ``generic``,  |
|                                  |                   | ``raw``, or            |
|                                  |                   | ``inputBlocks``        |
+----------------------------------+-------------------+------------------------+
| ``catalog.section.PROP``         | ``&PROP``         | ``set``, ``generic``,  |
|                                  |                   | ``raw``, or            |
|                                  |                   | ``inputBlocks``        |
+----------------------------------+-------------------+------------------------+
| ``catalog.section.PTDDFT``       | ``&PTDDFT``       | ``set``, ``generic``,  |
|                                  |                   | ``raw``, or            |
|                                  |                   | ``inputBlocks``        |
+----------------------------------+-------------------+------------------------+
| ``catalog.section.RESP``         | ``&RESP``         | ``set``, ``generic``,  |
|                                  |                   | ``raw``, or            |
|                                  |                   | ``inputBlocks``        |
+----------------------------------+-------------------+------------------------+
| ``catalog.section.SYSTEM``       | ``&SYSTEM``       | typed ``system``,      |
|                                  |                   | ``set``, ``generic``,  |
|                                  |                   | ``raw``, or            |
|                                  |                   | ``inputBlocks``        |
+----------------------------------+-------------------+------------------------+
| ``catalog.section.TDDFT``        | ``&TDDFT``        | ``set``, ``generic``,  |
|                                  |                   | ``raw``, or            |
|                                  |                   | ``inputBlocks``        |
+----------------------------------+-------------------+------------------------+
| ``catalog.section.VDW``          | ``&VDW``          | ``set``, ``generic``,  |
|                                  |                   | ``raw``, or            |
|                                  |                   | ``inputBlocks``        |
+----------------------------------+-------------------+------------------------+
| ``catalog.section.VECTORS``      | ``&VECTORS``      | ``set``, ``generic``,  |
|                                  |                   | ``raw``, or            |
|                                  |                   | ``inputBlocks``        |
+----------------------------------+-------------------+------------------------+
| ``catalog.section.WAVEFUNCTION`` | ``&WAVEFUNCTION`` | ``set``, ``generic``,  |
|                                  |                   | ``raw``, or            |
|                                  |                   | ``inputBlocks``        |
+----------------------------------+-------------------+------------------------+

Choosing an input carrier
=========================

Use the most structured carrier that represents the deck intent:

+----------------------+----------------------+----------------------+
| Requirement          | Carrier              | Reason               |
+======================+======================+======================+
| Supported ``&CPMD``, | typed section field  | Schema readers,      |
| ``&DFT``,            |                      | tests, and feature   |
| ``&SYSTEM``, or      |                      | discovery can name   |
| ``&ATOMS`` keyword   |                      | the field directly   |
+----------------------+----------------------+----------------------+
| One extra keyword    | ``set``              | The renderer merges  |
| inside a generated   |                      | it into the matching |
| section              |                      | section without      |
|                      |                      | duplicating defaults |
+----------------------+----------------------+----------------------+
| Unsupported section  | ``generic``          | The section remains  |
| with ordinary        |                      | structured and       |
| keyword/argument     |                      | discoverable         |
| lines                |                      |                      |
+----------------------+----------------------+----------------------+
| Existing deck text   | ``raw`` or           | The renderer does    |
| that must survive    | ``inputBlocks``      | not reinterpret the  |
| byte-for-byte apart  |                      | fragment             |
| from surrounding     |                      |                      |
| output               |                      |                      |
+----------------------+----------------------+----------------------+

For example, a compact typed setup with one extra ``&SYSTEM`` keyword
looks like:

.. code:: capnp

   (
     functional = "PBE0",
     cutOffRy = 85.0,
     multiplicity = 2,
     task = "gradient",
     inputSections = [
       ( system = (
           symmetry = 0,
           angstrom = true,
           cutOffRy = 85.0,
           cell = [12.0, 1.0, 1.0, 0.0, 0.0, 0.0]
         ) ),
       ( cpmd = (
           optimizeWavefunction = true,
           convergenceOrbitals = 1.0e-6,
           maxStep = 50,
           electronMass = 450.0
         ) ),
       ( dft = (
           functional = "PBE0",
           lsd = true,
           hfx = true,
           hfxScreening = "0.2"
         ) ),
       ( atoms = ( pseudopotentials = [
           ( element = "O", path = "O_MT_BLYP.psp", lmax = 1 ),
           ( element = "H", path = "H_CVB_BLYP.psp", lmax = 0 )
         ] ) ),
       ( set = ( key = "SYSTEM.POISSON SOLVER", value = "HOCKNEY" ) )
     ]
   )

Other structured parameter features
===================================

These ``params.inputSections.*`` IDs cover carriers and fields that are
not single CPMD/DFT catalog keyword rows in the typed tables below.

+-------------------------------------------------+----------------------------+--------------------------+
| Parameter feature ID                            | Field                      | Deck effect              |
+=================================================+============================+==========================+
| ``params.inputSections.generic.name``           | ``generic.name``           | Section name for an      |
|                                                 |                            | arbitrary ``&NAME``      |
|                                                 |                            | block                    |
+-------------------------------------------------+----------------------------+--------------------------+
| ``params.inputSections.generic.directives``     | ``generic.directives``     | Keyword/value lines      |
|                                                 |                            | inside the generic       |
|                                                 |                            | section                  |
+-------------------------------------------------+----------------------------+--------------------------+
| ``params.inputSections.system.symmetry``        | ``system.symmetry``        | ``&SYSTEM SYMMETRY``     |
+-------------------------------------------------+----------------------------+--------------------------+
| ``params.inputSections.system.angstrom``        | ``system.angstrom``        | ``&SYSTEM ANGSTROM``     |
+-------------------------------------------------+----------------------------+--------------------------+
| ``params.inputSections.system.cell``            | ``system.cell``            | ``&SYSTEM CELL``         |
+-------------------------------------------------+----------------------------+--------------------------+
| ``params.inputSections.system.cutOffRy``        | ``system.cutOffRy``        | ``&SYSTEM CUTOFF``       |
+-------------------------------------------------+----------------------------+--------------------------+
| ``params.inputSections.system.scale``           | ``system.scale``           | ``&SYSTEM SCALE``        |
+-------------------------------------------------+----------------------------+--------------------------+
| ``params.inputSections.system.charge``          | ``system.charge``          | ``&SYSTEM CHARGE``       |
+-------------------------------------------------+----------------------------+--------------------------+
| ``params.inputSections.system.multiplicity``    | ``system.multiplicity``    | ``&SYSTEM MULTIPLICITY`` |
+-------------------------------------------------+----------------------------+--------------------------+
| ``params.inputSections.system.directives``      | ``system.directives``      | Additional ``&SYSTEM``   |
|                                                 |                            | keyword/value lines      |
+-------------------------------------------------+----------------------------+--------------------------+
| ``params.inputSections.cpmd.directives``        | ``cpmd.directives``        | Additional ``&CPMD``     |
|                                                 |                            | keyword/value lines      |
+-------------------------------------------------+----------------------------+--------------------------+
| ``params.inputSections.dft.directives``         | ``dft.directives``         | Additional ``&DFT``      |
|                                                 |                            | keyword/value lines      |
+-------------------------------------------------+----------------------------+--------------------------+
| ``params.inputSections.atoms.pseudopotentials`` | ``atoms.pseudopotentials`` | Pseudopotential entries  |
|                                                 |                            | grouped with             |
|                                                 |                            | ``ForceInput``           |
|                                                 |                            | coordinates              |
+-------------------------------------------------+----------------------------+--------------------------+
| ``params.inputSections.atoms.directives``       | ``atoms.directives``       | Additional               |
|                                                 |                            | non-coordinate           |
|                                                 |                            | ``&ATOMS`` keyword/value |
|                                                 |                            | lines                    |
+-------------------------------------------------+----------------------------+--------------------------+
| ``params.inputSections.set.key``                | ``set.key``                | Dotted                   |
|                                                 |                            | ``SECTION.KEYWORD``      |
|                                                 |                            | merge target             |
+-------------------------------------------------+----------------------------+--------------------------+
| ``params.inputSections.set.value``              | ``set.value``              | Optional value emitted   |
|                                                 |                            | under the dotted key     |
+-------------------------------------------------+----------------------------+--------------------------+
| ``params.inputSections.raw``                    | ``raw``                    | Full section text        |
|                                                 |                            | inserted as-is           |
+-------------------------------------------------+----------------------------+--------------------------+

Typed ``&CPMD`` controls
========================

``CPMDCpmdSection`` exposes these driver controls directly. The feature
IDs are the stable discovery keys returned by ``cpmdc_feature_table()``.

+-----------------------------------------------+----------------------------------------------------------+--------------------------------+----------------------------------+
| Catalog feature ID                            | Parameter feature ID                                     | Field                          | Deck keyword                     |
+===============================================+==========================================================+================================+==================================+
| ``catalog.cpmd.OPTIMIZE_WAVEFUNCTION``        | ``params.inputSections.cpmd.optimizeWavefunction``       | ``optimizeWavefunction``       | ``OPTIMIZE WAVEFUNCTION``        |
+-----------------------------------------------+----------------------------------------------------------+--------------------------------+----------------------------------+
| ``catalog.cpmd.OPTIMIZE_GEOMETRY``            | ``params.inputSections.cpmd.optimizeGeometry``           | ``optimizeGeometry``           | ``OPTIMIZE GEOMETRY``            |
+-----------------------------------------------+----------------------------------------------------------+--------------------------------+----------------------------------+
| ``catalog.cpmd.MOLECULAR_DYNAMICS``           | ``params.inputSections.cpmd.molecularDynamics``          | ``molecularDynamics``          | ``MOLECULAR DYNAMICS``           |
+-----------------------------------------------+----------------------------------------------------------+--------------------------------+----------------------------------+
| ``catalog.cpmd.MOLECULAR_DYNAMICS_CP``        | ``params.inputSections.cpmd.molecularDynamicsCp``        | ``molecularDynamicsCp``        | ``MOLECULAR DYNAMICS CP``        |
+-----------------------------------------------+----------------------------------------------------------+--------------------------------+----------------------------------+
| ``catalog.cpmd.MOLECULAR_DYNAMICS_BO``        | ``params.inputSections.cpmd.molecularDynamicsBo``        | ``molecularDynamicsBo``        | ``MOLECULAR DYNAMICS BO``        |
+-----------------------------------------------+----------------------------------------------------------+--------------------------------+----------------------------------+
| ``catalog.cpmd.MOLECULAR_DYNAMICS_EH``        | ``params.inputSections.cpmd.molecularDynamicsEh``        | ``molecularDynamicsEh``        | ``MOLECULAR DYNAMICS EH``        |
+-----------------------------------------------+----------------------------------------------------------+--------------------------------+----------------------------------+
| ``catalog.cpmd.MOLECULAR_DYNAMICS_PT``        | ``params.inputSections.cpmd.molecularDynamicsPt``        | ``molecularDynamicsPt``        | ``MOLECULAR DYNAMICS PT``        |
+-----------------------------------------------+----------------------------------------------------------+--------------------------------+----------------------------------+
| ``catalog.cpmd.MOLECULAR_DYNAMICS_CLASSICAL`` | ``params.inputSections.cpmd.molecularDynamicsClassical`` | ``molecularDynamicsClassical`` | ``MOLECULAR DYNAMICS CLASSICAL`` |
+-----------------------------------------------+----------------------------------------------------------+--------------------------------+----------------------------------+
| ``catalog.cpmd.MOLECULAR_DYNAMICS_FILE``      | ``params.inputSections.cpmd.molecularDynamicsFile``      | ``molecularDynamicsFile``      | ``MOLECULAR DYNAMICS FILE``      |
+-----------------------------------------------+----------------------------------------------------------+--------------------------------+----------------------------------+
| ``catalog.cpmd.EMASS``                        | ``params.inputSections.cpmd.electronMass``               | ``electronMass``               | ``EMASS``                        |
+-----------------------------------------------+----------------------------------------------------------+--------------------------------+----------------------------------+
| ``catalog.cpmd.TIMESTEP``                     | ``params.inputSections.cpmd.timestep``                   | ``timestep``                   | ``TIMESTEP``                     |
+-----------------------------------------------+----------------------------------------------------------+--------------------------------+----------------------------------+
| ``catalog.cpmd.MAXSTEP``                      | ``params.inputSections.cpmd.maxStep``                    | ``maxStep``                    | ``MAXSTEP``                      |
+-----------------------------------------------+----------------------------------------------------------+--------------------------------+----------------------------------+
| ``catalog.cpmd.MAXITER``                      | ``params.inputSections.cpmd.maxIter``                    | ``maxIter``                    | ``MAXITER``                      |
+-----------------------------------------------+----------------------------------------------------------+--------------------------------+----------------------------------+
| ``catalog.cpmd.CONVERGENCE_ORBITALS``         | ``params.inputSections.cpmd.convergenceOrbitals``        | ``convergenceOrbitals``        | ``CONVERGENCE ORBITALS``         |
+-----------------------------------------------+----------------------------------------------------------+--------------------------------+----------------------------------+
| ``catalog.cpmd.CONVERGENCE_GEOMETRY``         | ``params.inputSections.cpmd.convergenceGeometry``        | ``convergenceGeometry``        | ``CONVERGENCE GEOMETRY``         |
+-----------------------------------------------+----------------------------------------------------------+--------------------------------+----------------------------------+
| ``catalog.cpmd.NOSE``                         | ``params.inputSections.cpmd.nose``                       | ``nose``                       | ``NOSE``                         |
+-----------------------------------------------+----------------------------------------------------------+--------------------------------+----------------------------------+
| ``catalog.cpmd.NOSE_IONS``                    | ``params.inputSections.cpmd.noseIons``                   | ``noseIons``                   | ``NOSE IONS``                    |
+-----------------------------------------------+----------------------------------------------------------+--------------------------------+----------------------------------+
| ``catalog.cpmd.NOSE_ELECTRONS``               | ``params.inputSections.cpmd.noseElectrons``              | ``noseElectrons``              | ``NOSE ELECTRONS``               |
+-----------------------------------------------+----------------------------------------------------------+--------------------------------+----------------------------------+
| ``catalog.cpmd.BERENDSEN``                    | ``params.inputSections.cpmd.berendsen``                  | ``berendsen``                  | ``BERENDSEN``                    |
+-----------------------------------------------+----------------------------------------------------------+--------------------------------+----------------------------------+
| ``catalog.cpmd.LANGEVIN``                     | ``params.inputSections.cpmd.langevin``                   | ``langevin``                   | ``LANGEVIN``                     |
+-----------------------------------------------+----------------------------------------------------------+--------------------------------+----------------------------------+
| ``catalog.cpmd.ANNEALING``                    | ``params.inputSections.cpmd.annealing``                  | ``annealing``                  | ``ANNEALING``                    |
+-----------------------------------------------+----------------------------------------------------------+--------------------------------+----------------------------------+
| ``catalog.cpmd.QUENCH``                       | ``params.inputSections.cpmd.quench``                     | ``quench``                     | ``QUENCH``                       |
+-----------------------------------------------+----------------------------------------------------------+--------------------------------+----------------------------------+
| ``catalog.cpmd.RATTLE``                       | ``params.inputSections.cpmd.rattle``                     | ``rattle``                     | ``RATTLE``                       |
+-----------------------------------------------+----------------------------------------------------------+--------------------------------+----------------------------------+
| ``catalog.cpmd.SHAKE``                        | ``params.inputSections.cpmd.shake``                      | ``shake``                      | ``SHAKE``                        |
+-----------------------------------------------+----------------------------------------------------------+--------------------------------+----------------------------------+
| ``catalog.cpmd.CONSTRAINT``                   | ``params.inputSections.cpmd.constraint``                 | ``constraint``                 | ``CONSTRAINT``                   |
+-----------------------------------------------+----------------------------------------------------------+--------------------------------+----------------------------------+
| ``catalog.cpmd.TROTTER``                      | ``params.inputSections.cpmd.trotter``                    | ``trotter``                    | ``TROTTER``                      |
+-----------------------------------------------+----------------------------------------------------------+--------------------------------+----------------------------------+
| ``catalog.cpmd.RESTART``                      | ``params.inputSections.cpmd.restart``                    | ``restart``                    | ``RESTART``                      |
+-----------------------------------------------+----------------------------------------------------------+--------------------------------+----------------------------------+
| ``catalog.cpmd.RESTART_WAVEFUNCTION``         | ``params.inputSections.cpmd.restartWavefunction``        | ``restartWavefunction``        | ``RESTART WAVEFUNCTION``         |
+-----------------------------------------------+----------------------------------------------------------+--------------------------------+----------------------------------+
| ``catalog.cpmd.TRAJECTORY``                   | ``params.inputSections.cpmd.trajectory``                 | ``trajectory``                 | ``TRAJECTORY``                   |
+-----------------------------------------------+----------------------------------------------------------+--------------------------------+----------------------------------+
| ``catalog.cpmd.PRINT``                        | ``params.inputSections.cpmd.printOptions``               | ``printOptions``               | ``PRINT``                        |
+-----------------------------------------------+----------------------------------------------------------+--------------------------------+----------------------------------+
| ``catalog.cpmd.STORE``                        | ``params.inputSections.cpmd.storeOptions``               | ``storeOptions``               | ``STORE``                        |
+-----------------------------------------------+----------------------------------------------------------+--------------------------------+----------------------------------+
| ``catalog.cpmd.CENTER_MOLECULE_OFF``          | ``params.inputSections.cpmd.centerMoleculeOff``          | ``centerMoleculeOff``          | ``CENTER MOLECULE OFF``          |
+-----------------------------------------------+----------------------------------------------------------+--------------------------------+----------------------------------+
| ``catalog.cpmd.CENTER_MOLECULE_ON``           | ``params.inputSections.cpmd.centerMoleculeOn``           | ``centerMoleculeOn``           | ``CENTER MOLECULE ON``           |
+-----------------------------------------------+----------------------------------------------------------+--------------------------------+----------------------------------+
| ``catalog.cpmd.DIIS``                         | ``params.inputSections.cpmd.diis``                       | ``diis``                       | ``DIIS``                         |
+-----------------------------------------------+----------------------------------------------------------+--------------------------------+----------------------------------+
| ``catalog.cpmd.ODIIS``                        | ``params.inputSections.cpmd.odiis``                      | ``odiis``                      | ``ODIIS``                        |
+-----------------------------------------------+----------------------------------------------------------+--------------------------------+----------------------------------+
| ``catalog.cpmd.PCG``                          | ``params.inputSections.cpmd.pcg``                        | ``pcg``                        | ``PCG``                          |
+-----------------------------------------------+----------------------------------------------------------+--------------------------------+----------------------------------+
| ``catalog.cpmd.DIAGONALIZATION``              | ``params.inputSections.cpmd.diagonalization``            | ``diagonalization``            | ``DIAGONALIZATION``              |
+-----------------------------------------------+----------------------------------------------------------+--------------------------------+----------------------------------+
| ``catalog.cpmd.FREE_ENERGY``                  | ``params.inputSections.cpmd.freeEnergy``                 | ``freeEnergy``                 | ``FREE-ENERGY``                  |
+-----------------------------------------------+----------------------------------------------------------+--------------------------------+----------------------------------+
| ``catalog.cpmd.INTERFACE``                    | ``params.inputSections.cpmd.interface``                  | ``interface``                  | ``INTERFACE``                    |
+-----------------------------------------------+----------------------------------------------------------+--------------------------------+----------------------------------+
| ``catalog.cpmd.QMMM``                         | ``params.inputSections.cpmd.qmmm``                       | ``qmmm``                       | ``QMMM``                         |
+-----------------------------------------------+----------------------------------------------------------+--------------------------------+----------------------------------+
| ``catalog.cpmd.BICANONICAL_ENSEMBLE``         | ``params.inputSections.cpmd.bicanonicalEnsemble``        | ``bicanonicalEnsemble``        | ``BICANONICAL ENSEMBLE``         |
+-----------------------------------------------+----------------------------------------------------------+--------------------------------+----------------------------------+
| ``catalog.cpmd.CDFT``                         | ``params.inputSections.cpmd.cdft``                       | ``cdft``                       | ``CDFT``                         |
+-----------------------------------------------+----------------------------------------------------------+--------------------------------+----------------------------------+
| ``catalog.cpmd.PROPERTIES``                   | ``params.inputSections.cpmd.properties``                 | ``properties``                 | ``PROPERTIES``                   |
+-----------------------------------------------+----------------------------------------------------------+--------------------------------+----------------------------------+

Typed ``&DFT`` controls
=======================

``CPMDDftSection`` exposes these DFT scalars directly:

+-------------------------------+-------------------------------------------+------------------+-------------------+
| Catalog feature ID            | Parameter feature ID                      | Field            | Deck keyword      |
+===============================+===========================================+==================+===================+
| ``catalog.dft.FUNCTIONAL``    | ``params.inputSections.dft.functional``   | ``functional``   | ``FUNCTIONAL``    |
+-------------------------------+-------------------------------------------+------------------+-------------------+
| ``catalog.dft.LSD``           | ``params.inputSections.dft.lsd``          | ``lsd``          | ``LSD``           |
+-------------------------------+-------------------------------------------+------------------+-------------------+
| ``catalog.dft.GC_CUTOFF``     | ``params.inputSections.dft.gcCutoff``     | ``gcCutoff``     | ``GC-CUTOFF``     |
+-------------------------------+-------------------------------------------+------------------+-------------------+
| ``catalog.dft.OLDCODE``       | ``params.inputSections.dft.oldCode``      | ``oldCode``      | ``OLDCODE``       |
+-------------------------------+-------------------------------------------+------------------+-------------------+
| ``catalog.dft.NEWCODE``       | ``params.inputSections.dft.newCode``      | ``newCode``      | ``NEWCODE``       |
+-------------------------------+-------------------------------------------+------------------+-------------------+
| ``catalog.dft.XC_DRIVER``     | ``params.inputSections.dft.xcDriver``     | ``xcDriver``     | ``XC_DRIVER``     |
+-------------------------------+-------------------------------------------+------------------+-------------------+
| ``catalog.dft.LIBXC``         | ``params.inputSections.dft.libxc``        | ``libxc``        | ``LIBXC``         |
+-------------------------------+-------------------------------------------+------------------+-------------------+
| ``catalog.dft.LR_KERNEL``     | ``params.inputSections.dft.lrKernel``     | ``lrKernel``     | ``LR KERNEL``     |
+-------------------------------+-------------------------------------------+------------------+-------------------+
| ``catalog.dft.REFUNCT``       | ``params.inputSections.dft.refunct``      | ``refunct``      | ``REFUNCT``       |
+-------------------------------+-------------------------------------------+------------------+-------------------+
| ``catalog.dft.MTS_HIGH_FUNC`` | ``params.inputSections.dft.mtsHighFunc``  | ``mtsHighFunc``  | ``MTS_HIGH_FUNC`` |
+-------------------------------+-------------------------------------------+------------------+-------------------+
| ``catalog.dft.MTS_LOW_FUNC``  | ``params.inputSections.dft.mtsLowFunc``   | ``mtsLowFunc``   | ``MTS_LOW_FUNC``  |
+-------------------------------+-------------------------------------------+------------------+-------------------+
| ``catalog.dft.HFX``           | ``params.inputSections.dft.hfx``          | ``hfx``          | ``HFX``           |
+-------------------------------+-------------------------------------------+------------------+-------------------+
| ``catalog.dft.HFX_SCREENING`` | ``params.inputSections.dft.hfxScreening`` | ``hfxScreening`` | ``HFX-SCREENING`` |
+-------------------------------+-------------------------------------------+------------------+-------------------+
| ``catalog.dft.HUBBARD``       | ``params.inputSections.dft.hubbard``      | ``hubbard``      | ``HUBBARD``       |
+-------------------------------+-------------------------------------------+------------------+-------------------+
| ``catalog.dft.CORRELATION``   | ``params.inputSections.dft.correlation``  | ``correlation``  | ``CORRELATION``   |
+-------------------------------+-------------------------------------------+------------------+-------------------+
| ``catalog.dft.EXCHANGE``      | ``params.inputSections.dft.exchange``     | ``exchange``     | ``EXCHANGE``      |
+-------------------------------+-------------------------------------------+------------------+-------------------+
| ``catalog.dft.BECKE88``       | ``params.inputSections.dft.becke88``      | ``becke88``      | ``BECKE88``       |
+-------------------------------+-------------------------------------------+------------------+-------------------+
| ``catalog.dft.ALPHA``         | ``params.inputSections.dft.alpha``        | ``alpha``        | ``ALPHA``         |
+-------------------------------+-------------------------------------------+------------------+-------------------+
| ``catalog.dft.BETA``          | ``params.inputSections.dft.beta``         | ``beta``         | ``BETA``          |
+-------------------------------+-------------------------------------------+------------------+-------------------+

``catalog.dft.FUNCTIONAL_*`` IDs name accepted ``functional`` values:

======================================= ====================
Feature ID                              ``functional`` value
======================================= ====================
``catalog.dft.FUNCTIONAL_LDA``          ``LDA``
``catalog.dft.FUNCTIONAL_BONLY``        ``BONLY``
``catalog.dft.FUNCTIONAL_BP``           ``BP``
``catalog.dft.FUNCTIONAL_BP86``         ``BP86``
``catalog.dft.FUNCTIONAL_BLYP``         ``BLYP``
``catalog.dft.FUNCTIONAL_GGA``          ``GGA``
``catalog.dft.FUNCTIONAL_PW91``         ``PW91``
``catalog.dft.FUNCTIONAL_PBE``          ``PBE``
``catalog.dft.FUNCTIONAL_REVPBE``       ``REVPBE``
``catalog.dft.FUNCTIONAL_HCTH``         ``HCTH``
``catalog.dft.FUNCTIONAL_OPTX``         ``OPTX``
``catalog.dft.FUNCTIONAL_OLYP``         ``OLYP``
``catalog.dft.FUNCTIONAL_XLYP``         ``XLYP``
``catalog.dft.FUNCTIONAL_B3LYP``        ``B3LYP``
``catalog.dft.FUNCTIONAL_B1LYP``        ``B1LYP``
``catalog.dft.FUNCTIONAL_PBE0``         ``PBE0``
``catalog.dft.FUNCTIONAL_REVPBE0``      ``REVPBE0``
``catalog.dft.FUNCTIONAL_PBES``         ``PBES``
``catalog.dft.FUNCTIONAL_PBE1W``        ``PBE1W``
``catalog.dft.FUNCTIONAL_TPSS``         ``TPSS``
``catalog.dft.FUNCTIONAL_HARTREE``      ``HARTREE``
``catalog.dft.FUNCTIONAL_HARTREE_FOCK`` ``HARTREE-FOCK``
``catalog.dft.FUNCTIONAL_HSE06``        ``HSE06``
``catalog.dft.FUNCTIONAL_ACM0``         ``ACM0``
``catalog.dft.FUNCTIONAL_ACM1``         ``ACM1``
``catalog.dft.FUNCTIONAL_ACM3``         ``ACM3``
======================================= ====================

``PotentialConfig``
===================

For rgpot / multi-backend configure RPCs, ``PotentialConfig`` is a
tagged union with ``cpmd @2 :CPMDParams`` (ordinal aligned with the
shared schema; ``nwchem @1`` is reserved void in this package).
