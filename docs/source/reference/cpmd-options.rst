Contract
========

All method knobs for the OpenCPMD backend are fields of ``CPMDParams``
in ``schema/Potentials.capnp``. Geometry is never part of
``CPMDParams``; it is always a ``ForceInput`` step message (or the raw C
array entry points that take positions and atomic numbers).

Top-level ``CPMDParams`` fields
===============================

+----------------------------------+-----------------------------------+
| Field                            | Role                              |
+==================================+===================================+
| ``functional``                   | Default DFT functional (also used |
|                                  | when no ``dft`` section is        |
|                                  | present)                          |
+----------------------------------+-----------------------------------+
| ``cutOffRy``                     | Plane-wave cutoff in Rydberg      |
+----------------------------------+-----------------------------------+
| ``charge`` / ``multiplicity``    | System charge and 2S+1            |
|                                  | (multiplicity ``> 1`` implies LSD |
|                                  | defaults)                         |
+----------------------------------+-----------------------------------+
| ``task``                         | Frontend hint: ``energy``,        |
|                                  | ``gradient``, ``md``,             |
|                                  | ``optimize``                      |
+----------------------------------+-----------------------------------+
| ``title``                        | Optional comment header in        |
|                                  | rendered decks                    |
+----------------------------------+-----------------------------------+
| ``memoryMb``                     | Frontend memory hint              |
+----------------------------------+-----------------------------------+
| ``scratchDir``, ``permanentDir`` | ``&CPMD FILEPATH`` runtime file   |
|                                  | placement; ``permanentDir`` wins  |
+----------------------------------+-----------------------------------+
| ``cpmdRoot``                     | OpenCPMD source/build tree        |
+----------------------------------+-----------------------------------+
| ``enginePath``                   | Frontend engine selection hint    |
+----------------------------------+-----------------------------------+
| ``inputBlocks``                  | Raw ``&SECTION`` text blocks      |
|                                  | prepended to the deck             |
+----------------------------------+-----------------------------------+
| ``inputSections``                | Structured ``CPMDInputSection``   |
|                                  | list                              |
+----------------------------------+-----------------------------------+

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

Typed ``&CPMD`` controls
========================

``CPMDCpmdSection`` exposes the common driver controls directly:

============================== ================================
Field                          Deck keyword
============================== ================================
``optimizeWavefunction``       ``OPTIMIZE WAVEFUNCTION``
``optimizeGeometry``           ``OPTIMIZE GEOMETRY``
``molecularDynamics``          ``MOLECULAR DYNAMICS``
``molecularDynamicsCp``        ``MOLECULAR DYNAMICS CP``
``molecularDynamicsBo``        ``MOLECULAR DYNAMICS BO``
``molecularDynamicsEh``        ``MOLECULAR DYNAMICS EH``
``molecularDynamicsPt``        ``MOLECULAR DYNAMICS PT``
``molecularDynamicsClassical`` ``MOLECULAR DYNAMICS CLASSICAL``
``molecularDynamicsFile``      ``MOLECULAR DYNAMICS FILE``
``convergenceOrbitals``        ``CONVERGENCE ORBITALS``
``convergenceGeometry``        ``CONVERGENCE GEOMETRY``
``maxStep``                    ``MAXSTEP``
``maxIter``                    ``MAXITER``
``timestep``                   ``TIMESTEP``
``electronMass``               ``EMASS``
``nose``                       ``NOSE``
``noseIons``                   ``NOSE IONS``
``noseElectrons``              ``NOSE ELECTRONS``
``berendsen``                  ``BERENDSEN``
``langevin``                   ``LANGEVIN``
``annealing``                  ``ANNEALING``
``quench``                     ``QUENCH``
``rattle``                     ``RATTLE``
``shake``                      ``SHAKE``
``constraint``                 ``CONSTRAINT``
``trotter``                    ``TROTTER``
``restart``                    ``RESTART``
``printOptions``               ``PRINT``
``storeOptions``               ``STORE``
``centerMoleculeOff``          ``CENTER MOLECULE OFF``
``centerMoleculeOn``           ``CENTER MOLECULE ON``
``diis``                       ``DIIS``
``odiis``                      ``ODIIS``
``pcg``                        ``PCG``
``diagonalization``            ``DIAGONALIZATION``
``freeEnergy``                 ``FREE-ENERGY``
``interface``                  ``INTERFACE``
``qmmm``                       ``QMMM``
``bicanonicalEnsemble``        ``BICANONICAL ENSEMBLE``
``cdft``                       ``CDFT``
``properties``                 ``PROPERTIES``
``restartWavefunction``        ``RESTART WAVEFUNCTION``
``trajectory``                 ``TRAJECTORY``
============================== ================================

Typed ``&DFT`` controls
=======================

``CPMDDftSection`` exposes common DFT scalars directly:

================ =================
Field            Deck keyword
================ =================
``functional``   ``FUNCTIONAL``
``lsd``          ``LSD``
``gcCutoff``     ``GC-CUTOFF``
``xcDriver``     ``XC_DRIVER``
``libxc``        ``LIBXC``
``lrKernel``     ``LR KERNEL``
``refunct``      ``REFUNCT``
``mtsHighFunc``  ``MTS_HIGH_FUNC``
``mtsLowFunc``   ``MTS_LOW_FUNC``
``hfx``          ``HFX``
``hfxScreening`` ``HFX-SCREENING``
``hubbard``      ``HUBBARD``
``alpha``        ``ALPHA``
``beta``         ``BETA``
``oldCode``      ``OLDCODE``
``newCode``      ``NEWCODE``
``correlation``  ``CORRELATION``
``exchange``     ``EXCHANGE``
``becke88``      ``BECKE88``
================ =================

``PotentialConfig``
===================

For rgpot / multi-backend configure RPCs, ``PotentialConfig`` is a
tagged union with ``cpmd @2 :CPMDParams`` (ordinal aligned with the
shared schema; ``nwchem @1`` is reserved void in this package).
