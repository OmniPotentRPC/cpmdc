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
``params.inputSections.pimd.directives``.

Feature ID namespaces
=====================

Feature IDs separate writable schema fields from parser/catalog
capability rows.

+-----------------------------------+----------------------------------+
| Namespace                         | Meaning                          |
+===================================+==================================+
| ``abi.*``                         | Exported C ABI symbols           |
+-----------------------------------+----------------------------------+
| ``params.*``                      | Top-level ``CPMDParams`` fields  |
+-----------------------------------+----------------------------------+
| ``params.inputSections.*``        | Typed fields and directive       |
|                                   | carriers inside                  |
|                                   | ``CPMDInputSection``             |
+-----------------------------------+----------------------------------+
| ``params.inputSections.system.*`` | Structured ``&SYSTEM`` controls  |
|                                   | such as cell, cutoff, Poisson,   |
|                                   | KPOINTS, CDFT, pressure, and     |
|                                   | stress fields                    |
+-----------------------------------+----------------------------------+
| ``catalog.cpmd.*``                | Named ``&CPMD`` keyword          |
|                                   | capability rows                  |
+-----------------------------------+----------------------------------+
| ``catalog.dft.*``                 | Named ``&DFT`` keyword           |
|                                   | capability rows                  |
+-----------------------------------+----------------------------------+
| ``catalog.section.*``             | OpenCPMD section names           |
|                                   | discovered from parser           |
|                                   | ``inscan('&SECTION')`` calls     |
+-----------------------------------+----------------------------------+

``catalog.section.SYSTEM`` means the ``&SYSTEM`` deck section is
recognized as a section kind. Individual ``&SYSTEM`` controls stay under
``params.inputSections.system.*`` rather than a separate
``catalog.system.*`` namespace.

Manual keyword names do not always map to a ``catalog.cpmd.*`` feature
ID. Use the field namespace when the renderer owns a typed section
field:

+-----------------------------------+-------------------------------------------------------------+-----------------------------------+
| Manual name                       | Feature or field ID                                         | Rendered deck form                |
+===================================+=============================================================+===================================+
| ``BOX WALLS``                     | ``params.inputSections.system.boxWalls``                    | ``&SYSTEM`` / ``BOX WALLS``       |
+-----------------------------------+-------------------------------------------------------------+-----------------------------------+
| ``MODIFIED GOEDECKER``            | ``params.inputSections.system.modifiedGoedecker``           | ``&SYSTEM`` /                     |
|                                   |                                                             | ``MODIFIED GOEDECKER``            |
+-----------------------------------+-------------------------------------------------------------+-----------------------------------+
| ``MODIFIED GOEDECKER PARAMETERS`` | ``params.inputSections.system.modifiedGoedeckerParameters`` | ``&SYSTEM`` /                     |
|                                   |                                                             | ``MODIFIED GOEDECKER PARAMETERS`` |
+-----------------------------------+-------------------------------------------------------------+-----------------------------------+
| ``SOC``                           | ``catalog.cpmd.SOC`` and                                    | ``&CPMD`` /                       |
|                                   | ``params.inputSections.cpmd.spinOrbitCouplingStates``       | ``SPIN-ORBIT COUPLING``           |
+-----------------------------------+-------------------------------------------------------------+-----------------------------------+

``FFTW WISDOM`` is not a typed feature row because the OpenCPMD manual
source comments out that keyword and the parser tree used for inventory
probing has no matching branch. Use ``inputSections.raw`` or
``inputBlocks`` only for a downstream tree that still accepts it.

``CPMDInputSection`` kinds
==========================

+------------------+--------------------------+--------------------------+
| Kind             | Cap'n Proto type         | Deck effect              |
+==================+==========================+==========================+
| ``atom``         | ``CPMDDirectiveSection`` | ``&ATOM`` keyword/value  |
|                  |                          | lines and nested         |
|                  |                          | subsections              |
+------------------+--------------------------+--------------------------+
| ``cpmd``         | ``CPMDCpmdSection``      | ``&CPMD``                |
|                  |                          | (wavefunction/geometry   |
|                  |                          | optimization, MD,        |
|                  |                          | isolated/centered        |
|                  |                          | molecule controls,       |
|                  |                          | convergence, iteration   |
|                  |                          | limits, timestep,        |
|                  |                          | electron mass, …)        |
+------------------+--------------------------+--------------------------+
| ``system``       | ``CPMDSystemSection``    | ``&SYSTEM`` (symmetry,   |
|                  |                          | cell/reference/classical |
|                  |                          | cell, cutoff/grid/mesh   |
|                  |                          | controls, state          |
|                  |                          | occupation, external     |
|                  |                          | fields, CDFT Gaussian    |
|                  |                          | controls,                |
|                  |                          | pressure/stress          |
|                  |                          | controls, Poisson        |
|                  |                          | solver,                  |
|                  |                          | surface/polymer/cluster, |
|                  |                          | charge, …)               |
+------------------+--------------------------+--------------------------+
| ``dft``          | ``CPMDDftSection``       | ``&DFT`` (functional,    |
|                  |                          | LSD, GC cutoff, XC       |
|                  |                          | driver, hybrid/Hubbard   |
|                  |                          | scalars, extra           |
|                  |                          | directives)              |
+------------------+--------------------------+--------------------------+
| ``atoms``        | ``CPMDAtomsSection``     | ``&ATOMS``               |
|                  |                          | (pseudopotential lines;  |
|                  |                          | coordinates come from    |
|                  |                          | ``ForceInput``)          |
+------------------+--------------------------+--------------------------+
| ``basis``        | ``CPMDDirectiveSection`` | ``&BASIS`` keyword/value |
|                  |                          | lines and nested         |
|                  |                          | subsections              |
+------------------+--------------------------+--------------------------+
| ``clas``         | ``CPMDDirectiveSection`` | ``&CLAS`` keyword/value  |
|                  |                          | lines and nested         |
|                  |                          | subsections              |
+------------------+--------------------------+--------------------------+
| ``eam``          | ``CPMDDirectiveSection`` | ``&EAM`` keyword/value   |
|                  |                          | lines and nested         |
|                  |                          | subsections              |
+------------------+--------------------------+--------------------------+
| ``exte``         | ``CPMDDirectiveSection`` | ``&EXTE`` keyword/value  |
|                  |                          | lines and nested         |
|                  |                          | subsections              |
+------------------+--------------------------+--------------------------+
| ``hardness``     | ``CPMDDirectiveSection`` | ``&HARDNESS``            |
|                  |                          | keyword/value lines and  |
|                  |                          | nested subsections       |
+------------------+--------------------------+--------------------------+
| ``info``         | ``CPMDDirectiveSection`` | ``&INFO`` keyword/value  |
|                  |                          | lines and nested         |
|                  |                          | subsections              |
+------------------+--------------------------+--------------------------+
| ``linres``       | ``CPMDDirectiveSection`` | ``&LINRES``              |
|                  |                          | keyword/value lines and  |
|                  |                          | nested subsections       |
+------------------+--------------------------+--------------------------+
| ``molstates``    | ``CPMDDirectiveSection`` | ``&MOLSTATES``           |
|                  |                          | keyword/value lines and  |
|                  |                          | nested subsections       |
+------------------+--------------------------+--------------------------+
| ``mts``          | ``CPMDDirectiveSection`` | ``&MTS`` keyword/value   |
|                  |                          | lines and nested         |
|                  |                          | subsections              |
+------------------+--------------------------+--------------------------+
| ``nlcc``         | ``CPMDDirectiveSection`` | ``&NLCC`` keyword/value  |
|                  |                          | lines and nested         |
|                  |                          | subsections              |
+------------------+--------------------------+--------------------------+
| ``path``         | ``CPMDDirectiveSection`` | ``&PATH`` keyword/value  |
|                  |                          | lines and nested         |
|                  |                          | subsections              |
+------------------+--------------------------+--------------------------+
| ``pimd``         | ``CPMDDirectiveSection`` | ``&PIMD`` keyword/value  |
|                  |                          | lines and nested         |
|                  |                          | subsections              |
+------------------+--------------------------+--------------------------+
| ``potential``    | ``CPMDDirectiveSection`` | ``&POTENTIAL``           |
|                  |                          | keyword/value lines and  |
|                  |                          | nested subsections       |
+------------------+--------------------------+--------------------------+
| ``prop``         | ``CPMDDirectiveSection`` | ``&PROP`` keyword/value  |
|                  |                          | lines and nested         |
|                  |                          | subsections              |
+------------------+--------------------------+--------------------------+
| ``ptddft``       | ``CPMDDirectiveSection`` | ``&PTDDFT``              |
|                  |                          | keyword/value lines and  |
|                  |                          | nested subsections       |
+------------------+--------------------------+--------------------------+
| ``resp``         | ``CPMDDirectiveSection`` | ``&RESP`` keyword/value  |
|                  |                          | lines and nested         |
|                  |                          | subsections              |
+------------------+--------------------------+--------------------------+
| ``tddft``        | ``CPMDDirectiveSection`` | ``&TDDFT`` keyword/value |
|                  |                          | lines and nested         |
|                  |                          | subsections              |
+------------------+--------------------------+--------------------------+
| ``vdw``          | ``CPMDDirectiveSection`` | ``&VDW`` keyword/value   |
|                  |                          | lines and nested         |
|                  |                          | subsections such as      |
|                  |                          | ``EMPIRICAL CORRECTION`` |
+------------------+--------------------------+--------------------------+
| ``vectors``      | ``CPMDDirectiveSection`` | ``&VECTORS``             |
|                  |                          | keyword/value lines and  |
|                  |                          | nested subsections       |
+------------------+--------------------------+--------------------------+
| ``wavefunction`` | ``CPMDDirectiveSection`` | ``&WAVEFUNCTION``        |
|                  |                          | keyword/value lines and  |
|                  |                          | nested subsections       |
+------------------+--------------------------+--------------------------+
| ``generic``      | ``CPMDGenericSection``   | Arbitrary ``&NAME`` …    |
|                  |                          | ``&END``                 |
+------------------+--------------------------+--------------------------+
| ``set``          | ``CPMDSetDirective``     | Embed-path logical       |
|                  |                          | ``SECTION.KEYWORD``      |
|                  |                          | merged into the section  |
+------------------+--------------------------+--------------------------+
| ``raw``          | ``Text``                 | Full section text        |
|                  |                          | inserted as-is           |
+------------------+--------------------------+--------------------------+

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

Long-tail CPMD keywords go in ``CPMDDirective`` lists or ``subsections``
inside the matching typed section, ``CPMDSetDirective`` dotted keys,
``CPMDInputSection.raw``, or ``inputBlocks``. ``&VDW`` uses
``subsections`` for nested blocks such as ``EMPIRICAL CORRECTION``.
``set.key`` uses ``SECTION.KEYWORD`` form, for example
``CPMD.PRINT FORCES ON`` or ``SYSTEM.POISSON SOLVER HOCKNEY``; non-empty
``set.value`` is emitted on the following indented line. Prefer new
structured fields in the schema over inventing a second config file
format.

Choosing an input carrier
=========================

Use the most structured carrier that represents the deck intent. The
usual path is:

#. Put functional, cutoff, charge, multiplicity, task, and embed hints
   in the top-level ``CPMDParams`` fields.
#. Put named OpenCPMD section controls in typed ``inputSections`` arms
   such as ``system``, ``cpmd``, ``dft``, and ``atoms``.
#. Put extra section keywords in that section's ``directives`` list when
   they naturally belong to a typed section.
#. Use ``set`` for one merge-only keyword, ``generic`` for a non-catalog
   section alias, and ``raw`` or ``inputBlocks`` for text that must
   remain literal.

+------------------------+----------------------+----------------------+
| Requirement            | Carrier              | Reason               |
+========================+======================+======================+
| Supported OpenCPMD     | typed section field  | Schema readers,      |
| ``inscan('&SECTION')`` |                      | tests, and feature   |
| block                  |                      | discovery can name   |
|                        |                      | the field directly   |
+------------------------+----------------------+----------------------+
| One extra keyword      | ``set``              | The renderer merges  |
| inside a generated     |                      | it into the matching |
| section                |                      | section without      |
|                        |                      | duplicating defaults |
+------------------------+----------------------+----------------------+
| Non-catalog section    | ``generic``          | The section remains  |
| alias with ordinary    |                      | structured even when |
| keyword/argument lines |                      | it is not an         |
|                        |                      | inventory ID         |
+------------------------+----------------------+----------------------+
| Existing deck text     | ``raw`` or           | The renderer does    |
| that must survive      | ``inputBlocks``      | not reinterpret the  |
| byte-for-byte apart    |                      | fragment             |
| from surrounding       |                      |                      |
| output                 |                      |                      |
+------------------------+----------------------+----------------------+

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
| ``catalog.section.ATOM``         | ``&ATOM``         | typed ``atom``,        |
|                                  |                   | ``set``, ``generic``,  |
|                                  |                   | ``raw``, or            |
|                                  |                   | ``inputBlocks``        |
+----------------------------------+-------------------+------------------------+
| ``catalog.section.ATOMS``        | ``&ATOMS``        | typed ``atoms``,       |
|                                  |                   | ``set``, ``generic``,  |
|                                  |                   | ``raw``, or            |
|                                  |                   | ``inputBlocks``        |
+----------------------------------+-------------------+------------------------+
| ``catalog.section.BASIS``        | ``&BASIS``        | typed ``basis``,       |
|                                  |                   | ``set``, ``generic``,  |
|                                  |                   | ``raw``, or            |
|                                  |                   | ``inputBlocks``        |
+----------------------------------+-------------------+------------------------+
| ``catalog.section.CLAS``         | ``&CLAS``         | typed ``clas``,        |
|                                  |                   | ``set``, ``generic``,  |
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
| ``catalog.section.EAM``          | ``&EAM``          | typed ``eam``,         |
|                                  |                   | ``set``, ``generic``,  |
|                                  |                   | ``raw``, or            |
|                                  |                   | ``inputBlocks``        |
+----------------------------------+-------------------+------------------------+
| ``catalog.section.EXTE``         | ``&EXTE``         | typed ``exte``,        |
|                                  |                   | ``set``, ``generic``,  |
|                                  |                   | ``raw``, or            |
|                                  |                   | ``inputBlocks``        |
+----------------------------------+-------------------+------------------------+
| ``catalog.section.HARDNESS``     | ``&HARDNESS``     | typed ``hardness``,    |
|                                  |                   | ``set``, ``generic``,  |
|                                  |                   | ``raw``, or            |
|                                  |                   | ``inputBlocks``        |
+----------------------------------+-------------------+------------------------+
| ``catalog.section.INFO``         | ``&INFO``         | typed ``info``,        |
|                                  |                   | ``set``, ``generic``,  |
|                                  |                   | ``raw``, or            |
|                                  |                   | ``inputBlocks``        |
+----------------------------------+-------------------+------------------------+
| ``catalog.section.LINRES``       | ``&LINRES``       | typed ``linres``,      |
|                                  |                   | ``set``, ``generic``,  |
|                                  |                   | ``raw``, or            |
|                                  |                   | ``inputBlocks``        |
+----------------------------------+-------------------+------------------------+
| ``catalog.section.MOLSTATES``    | ``&MOLSTATES``    | typed ``molstates``,   |
|                                  |                   | ``set``, ``generic``,  |
|                                  |                   | ``raw``, or            |
|                                  |                   | ``inputBlocks``        |
+----------------------------------+-------------------+------------------------+
| ``catalog.section.MTS``          | ``&MTS``          | typed ``mts``,         |
|                                  |                   | ``set``, ``generic``,  |
|                                  |                   | ``raw``, or            |
|                                  |                   | ``inputBlocks``        |
+----------------------------------+-------------------+------------------------+
| ``catalog.section.NLCC``         | ``&NLCC``         | typed ``nlcc``,        |
|                                  |                   | ``set``, ``generic``,  |
|                                  |                   | ``raw``, or            |
|                                  |                   | ``inputBlocks``        |
+----------------------------------+-------------------+------------------------+
| ``catalog.section.PATH``         | ``&PATH``         | typed ``path``,        |
|                                  |                   | ``set``, ``generic``,  |
|                                  |                   | ``raw``, or            |
|                                  |                   | ``inputBlocks``        |
+----------------------------------+-------------------+------------------------+
| ``catalog.section.PIMD``         | ``&PIMD``         | typed ``pimd``,        |
|                                  |                   | ``set``, ``generic``,  |
|                                  |                   | ``raw``, or            |
|                                  |                   | ``inputBlocks``        |
+----------------------------------+-------------------+------------------------+
| ``catalog.section.POTENTIAL``    | ``&POTENTIAL``    | typed ``potential``,   |
|                                  |                   | ``set``, ``generic``,  |
|                                  |                   | ``raw``, or            |
|                                  |                   | ``inputBlocks``        |
+----------------------------------+-------------------+------------------------+
| ``catalog.section.PROP``         | ``&PROP``         | typed ``prop``,        |
|                                  |                   | ``set``, ``generic``,  |
|                                  |                   | ``raw``, or            |
|                                  |                   | ``inputBlocks``        |
+----------------------------------+-------------------+------------------------+
| ``catalog.section.PTDDFT``       | ``&PTDDFT``       | typed ``ptddft``,      |
|                                  |                   | ``set``, ``generic``,  |
|                                  |                   | ``raw``, or            |
|                                  |                   | ``inputBlocks``        |
+----------------------------------+-------------------+------------------------+
| ``catalog.section.RESP``         | ``&RESP``         | typed ``resp``,        |
|                                  |                   | ``set``, ``generic``,  |
|                                  |                   | ``raw``, or            |
|                                  |                   | ``inputBlocks``        |
+----------------------------------+-------------------+------------------------+
| ``catalog.section.SYSTEM``       | ``&SYSTEM``       | typed ``system``,      |
|                                  |                   | ``set``, ``generic``,  |
|                                  |                   | ``raw``, or            |
|                                  |                   | ``inputBlocks``        |
+----------------------------------+-------------------+------------------------+
| ``catalog.section.TDDFT``        | ``&TDDFT``        | typed ``tddft``,       |
|                                  |                   | ``set``, ``generic``,  |
|                                  |                   | ``raw``, or            |
|                                  |                   | ``inputBlocks``        |
+----------------------------------+-------------------+------------------------+
| ``catalog.section.VDW``          | ``&VDW``          | typed ``vdw``,         |
|                                  |                   | ``set``, ``generic``,  |
|                                  |                   | ``raw``, or            |
|                                  |                   | ``inputBlocks``        |
+----------------------------------+-------------------+------------------------+
| ``catalog.section.VECTORS``      | ``&VECTORS``      | typed ``vectors``,     |
|                                  |                   | ``set``, ``generic``,  |
|                                  |                   | ``raw``, or            |
|                                  |                   | ``inputBlocks``        |
+----------------------------------+-------------------+------------------------+
| ``catalog.section.WAVEFUNCTION`` | ``&WAVEFUNCTION`` | typed                  |
|                                  |                   | ``wavefunction``,      |
|                                  |                   | ``set``, ``generic``,  |
|                                  |                   | ``raw``, or            |
|                                  |                   | ``inputBlocks``        |
+----------------------------------+-------------------+------------------------+

Other structured parameter features
===================================

These ``params.inputSections.*`` IDs cover carriers and fields that are
not single CPMD/DFT catalog keyword rows in the typed tables below.

+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| Parameter feature ID                                                  | Field                                            | Deck effect                               |
+=======================================================================+==================================================+===========================================+
| ``params.inputSections.generic.name``                                 | ``generic.name``                                 | Section name for an arbitrary ``&NAME``   |
|                                                                       |                                                  | block                                     |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.generic.directives``                           | ``generic.directives``                           | Keyword/value lines inside the generic    |
|                                                                       |                                                  | section                                   |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.atom.directives``                              | ``atom.directives``                              | Keyword/value lines inside ``&ATOM``      |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.atom.subsections``                             | ``atom.subsections``                             | Nested blocks inside ``&ATOM``            |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.symmetry``                              | ``system.symmetry``                              | ``&SYSTEM SYMMETRY``                      |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.angstrom``                              | ``system.angstrom``                              | ``&SYSTEM ANGSTROM``                      |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.cell``                                  | ``system.cell``                                  | ``&SYSTEM CELL``                          |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.cellAbsolute``                          | ``system.cellAbsolute``                          | Add ``ABSOLUTE`` to ``&SYSTEM CELL``      |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.cellDegree``                            | ``system.cellDegree``                            | Add ``DEGREE`` to ``&SYSTEM CELL``        |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.cellVectors``                           | ``system.cellVectors``                           | Add ``VECTORS`` to ``&SYSTEM CELL``       |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.referenceCell``                         | ``system.referenceCell``                         | ``&SYSTEM REFERENCE CELL``                |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.referenceCellAbsolute``                 | ``system.referenceCellAbsolute``                 | Add ``ABSOLUTE`` to ``REFERENCE CELL``    |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.referenceCellDegree``                   | ``system.referenceCellDegree``                   | Add ``DEGREE`` to ``REFERENCE CELL``      |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.referenceCellVectors``                  | ``system.referenceCellVectors``                  | Add ``VECTORS`` to ``REFERENCE CELL``     |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.classicalCell``                         | ``system.classicalCell``                         | ``&SYSTEM CLASSICAL CELL``                |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.classicalCellAbsolute``                 | ``system.classicalCellAbsolute``                 | Add ``ABSOLUTE`` to ``CLASSICAL CELL``    |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.classicalCellDegree``                   | ``system.classicalCellDegree``                   | Add ``DEGREE`` to ``CLASSICAL CELL``      |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.isotropicCell``                         | ``system.isotropicCell``                         | ``&SYSTEM ISOTROPIC CELL``                |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.zFlexibleCell``                         | ``system.zFlexibleCell``                         | ``&SYSTEM ZFLEXIBLE CELL``                |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.checkSymmetryPrecision``                | ``system.checkSymmetryPrecision``                | ``&SYSTEM CHECK SYMMETRY`` precision      |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.checkSymmetryOff``                      | ``system.checkSymmetryOff``                      | ``&SYSTEM CHECK SYMMETRY OFF``            |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.cutOffRy``                              | ``system.cutOffRy``                              | ``&SYSTEM CUTOFF``                        |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.cutoffShape``                           | ``system.cutoffShape``                           | ``&SYSTEM SPHERICAL/NOSPHERICAL CUTOFF``  |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.hfxCutoff``                             | ``system.hfxCutoff``                             | ``&SYSTEM HFX CUTOFF``                    |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.scale``                                 | ``system.scale``                                 | ``&SYSTEM SCALE S=<value>``               |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.scaleX``                                | ``system.scaleX``                                | ``&SYSTEM SCALE SX=<value>``              |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.scaleY``                                | ``system.scaleY``                                | ``&SYSTEM SCALE SY=<value>``              |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.scaleZ``                                | ``system.scaleZ``                                | ``&SYSTEM SCALE SZ=<value>``              |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.scaleCartesian``                        | ``system.scaleCartesian``                        | Add ``CARTESIAN`` to ``&SYSTEM SCALE``    |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.charge``                                | ``system.charge``                                | ``&SYSTEM CHARGE``                        |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.multiplicity``                          | ``system.multiplicity``                          | ``&SYSTEM MULTIPLICITY``                  |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.nSup``                                  | ``system.nSup``                                  | ``&SYSTEM NSUP``                          |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.states``                                | ``system.states``                                | ``&SYSTEM STATES``                        |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.occupation``                            | ``system.occupation``                            | ``&SYSTEM OCCUPATION`` values             |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.occupationFixed``                       | ``system.occupationFixed``                       | Add ``FIXED`` to ``&SYSTEM OCCUPATION``   |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.externalField``                         | ``system.externalField``                         | ``&SYSTEM EXTERNAL FIELD``                |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.wCut``                                  | ``system.wCut``                                  | ``&SYSTEM WCUT``                          |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.wGauss``                                | ``system.wGauss``                                | ``&SYSTEM WGAUSS`` sigma values           |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.couplingsFiniteDifference``             | ``system.couplingsFiniteDifference``             | Add ``FD`` to ``&SYSTEM COUPLINGS``       |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.couplingsFiniteDifferenceDisplacement`` | ``system.couplingsFiniteDifferenceDisplacement`` | ``&SYSTEM COUPLINGS FD=`` displacement    |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.couplingsProductDisplacement``          | ``system.couplingsProductDisplacement``          | ``&SYSTEM COUPLINGS PROD=`` displacement  |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.couplingsLinres``                       | ``system.couplingsLinres``                       | Add ``LINRES`` to ``&SYSTEM COUPLINGS``   |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.couplingsLinresTolerance``              | ``system.couplingsLinresTolerance``              | Add ``TOL=`` to ``COUPLINGS LINRES``      |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.couplingsLinresNvects``                 | ``system.couplingsLinresNvects``                 | Add ``NVECT=`` to ``COUPLINGS LINRES``    |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.couplingsLinresSpecify``                | ``system.couplingsLinresSpecify``                | Add ``SPECIFY`` to ``COUPLINGS LINRES``   |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.couplingsLinresThresholds``             | ``system.couplingsLinresThresholds``             | Add ``THRESHOLDS`` payload to             |
|                                                                       |                                                  | ``COUPLINGS LINRES``                      |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.couplingsLinresBruteForce``             | ``system.couplingsLinresBruteForce``             | Add ``BRUTE FORCE`` to                    |
|                                                                       |                                                  | ``COUPLINGS LINRES``                      |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.couplingsSurfaces``                     | ``system.couplingsSurfaces``                     | ``&SYSTEM COUPLINGS NSURF=`` surface      |
|                                                                       |                                                  | triples                                   |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.couplingsFiniteDifferenceAtoms``        | ``system.couplingsFiniteDifferenceAtoms``        | ``&SYSTEM COUPLINGS NAT=`` atom list      |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.cdftDonorAtoms``                        | ``system.cdftDonorAtoms``                        | ``&SYSTEM DONOR`` atom list               |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.cdftDonorWeights``                      | ``system.cdftDonorWeights``                      | Add ``WMULT`` weights to ``DONOR``        |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.cdftAcceptorAtoms``                     | ``system.cdftAcceptorAtoms``                     | ``&SYSTEM ACCEPTOR`` atom list            |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.cdftAcceptorHdasDonors``                | ``system.cdftAcceptorHdasDonors``                | Add ``HDAS`` donor list to ``ACCEPTOR``   |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.cdftAcceptorWeights``                   | ``system.cdftAcceptorWeights``                   | Add ``WMULT`` weights to ``ACCEPTOR``     |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.pointGroup``                            | ``system.pointGroup``                            | ``&SYSTEM POINT GROUP`` selector line     |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.pointGroupMolecule``                    | ``system.pointGroupMolecule``                    | Add ``MOLECULE`` to ``POINT GROUP``       |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.pointGroupDelta``                       | ``system.pointGroupDelta``                       | ``&SYSTEM POINT GROUP DELTA=`` accuracy   |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.lowSpinExcitation``                     | ``system.lowSpinExcitation``                     | ``&SYSTEM LOW SPIN EXCITATION`` options   |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.lowSpinExcitationLsets``                | ``system.lowSpinExcitationLsets``                | Add ``LSETS`` to ``LOW SPIN EXCITATION``  |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.lowSpinExcitationPenalty``              | ``system.lowSpinExcitationPenalty``              | ``&SYSTEM LOW SPIN EXCITATION PENALTY``   |
|                                                                       |                                                  | value                                     |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.lseParameters``                         | ``system.lseParameters``                         | ``&SYSTEM LSE PARAMETERS`` lsea/lseb      |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.modifiedGoedecker``                     | ``system.modifiedGoedecker``                     | ``&SYSTEM MODIFIED GOEDECKER``            |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.modifiedGoedeckerParameters``           | ``system.modifiedGoedeckerParameters``           | ``&SYSTEM MODIFIED GOEDECKER PARAMETERS`` |
|                                                                       |                                                  | lambda values                             |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.energyProfile``                         | ``system.energyProfile``                         | ``&SYSTEM ENERGY PROFILE``                |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.densityCutOffRy``                       | ``system.densityCutOffRy``                       | ``&SYSTEM DENSITY CUTOFF``                |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.densityCutoffNumber``                   | ``system.densityCutoffNumber``                   | ``&SYSTEM DENSITY CUTOFF NUMBER``         |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.dual``                                  | ``system.dual``                                  | ``&SYSTEM DUAL``                          |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.constantCutoff``                        | ``system.constantCutoff``                        | ``&SYSTEM CONSTANT CUTOFF``               |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.boxWalls``                              | ``system.boxWalls``                              | ``&SYSTEM BOX WALLS``                     |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.poissonSolver``                         | ``system.poissonSolver``                         | ``&SYSTEM POISSON SOLVER``                |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.poissonParameter``                      | ``system.poissonParameter``                      | Optional ``POISSON ... PARAMETER`` value  |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.mesh``                                  | ``system.mesh``                                  | ``&SYSTEM MESH``                          |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.kpoints``                               | ``system.kpoints``                               | ``&SYSTEM KPOINTS`` explicit weighted     |
|                                                                       |                                                  | points                                    |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.kpointsScaled``                         | ``system.kpointsScaled``                         | Add ``SCALED`` to ``&SYSTEM KPOINTS``     |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.kpointsOnlyDiagonal``                   | ``system.kpointsOnlyDiagonal``                   | Add ``ONLYDIAG`` to ``&SYSTEM KPOINTS``   |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.kpointsBlock``                          | ``system.kpointsBlock``                          | Add ``BLOCK=`` to ``&SYSTEM KPOINTS``     |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.kpointsBlockAll``                       | ``system.kpointsBlockAll``                       | Add ``ALL`` to ``KPOINTS BLOCK``          |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.kpointsBlockCalculated``                | ``system.kpointsBlockCalculated``                | Add ``CALCULATED`` to ``KPOINTS BLOCK``   |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.kpointsBlockNoSwap``                    | ``system.kpointsBlockNoSwap``                    | Add ``NOSWAP`` to ``KPOINTS BLOCK``       |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.kpointsMonkhorstPack``                  | ``system.kpointsMonkhorstPack``                  | ``&SYSTEM KPOINTS MONKHORST-PACK`` mesh   |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.kpointsMonkhorstSymmetrized``           | ``system.kpointsMonkhorstSymmetrized``           | Add ``SYMMETRIZED`` to                    |
|                                                                       |                                                  | ``KPOINTS MONKHORST-PACK``                |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.kpointsMonkhorstFull``                  | ``system.kpointsMonkhorstFull``                  | Add ``FULL`` to                           |
|                                                                       |                                                  | ``KPOINTS MONKHORST-PACK``                |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.kpointsMonkhorstKdp``                   | ``system.kpointsMonkhorstKdp``                   | Add ``KDP`` to ``KPOINTS MONKHORST-PACK`` |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.kpointsMonkhorstShift``                 | ``system.kpointsMonkhorstShift``                 | Add mesh-line ``SHIFT=`` to               |
|                                                                       |                                                  | ``KPOINTS MONKHORST-PACK``                |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.kpointBands``                           | ``system.kpointBands``                           | ``&SYSTEM KPOINTS BANDS`` segments        |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.doubleGrid``                            | ``system.doubleGrid``                            | ``&SYSTEM DOUBLE GRID`` ``ON`` or ``OFF`` |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.symmetrizeCoordinates``                 | ``system.symmetrizeCoordinates``                 | ``&SYSTEM SYMMETRIZE COORDINATES``        |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.tesr``                                  | ``system.tesr``                                  | ``&SYSTEM TESR``                          |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.surface``                               | ``system.surface``                               | ``&SYSTEM SURFACE`` direction             |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.polymer``                               | ``system.polymer``                               | ``&SYSTEM POLYMER``                       |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.cluster``                               | ``system.cluster``                               | ``&SYSTEM CLUSTER``                       |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.pressure``                              | ``system.pressure``                              | ``&SYSTEM PRESSURE``                      |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.stressTensor``                          | ``system.stressTensor``                          | ``&SYSTEM STRESS TENSOR``                 |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.shockVelocity``                         | ``system.shockVelocity``                         | ``&SYSTEM SHOCK VELOCITY``                |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.system.directives``                            | ``system.directives``                            | Additional ``&SYSTEM`` keyword/value      |
|                                                                       |                                                  | lines                                     |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.basis.directives``                             | ``basis.directives``                             | Keyword/value lines inside ``&BASIS``     |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.basis.subsections``                            | ``basis.subsections``                            | Nested blocks inside ``&BASIS``           |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.cpmd.directives``                              | ``cpmd.directives``                              | Additional ``&CPMD`` keyword/value lines  |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.dft.directives``                               | ``dft.directives``                               | Additional ``&DFT`` keyword/value lines   |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.clas.directives``                              | ``clas.directives``                              | Keyword/value lines inside ``&CLAS``      |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.clas.subsections``                             | ``clas.subsections``                             | Nested blocks inside ``&CLAS``            |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.eam.directives``                               | ``eam.directives``                               | Keyword/value lines inside ``&EAM``       |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.eam.subsections``                              | ``eam.subsections``                              | Nested blocks inside ``&EAM``             |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.exte.directives``                              | ``exte.directives``                              | Keyword/value lines inside ``&EXTE``      |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.exte.subsections``                             | ``exte.subsections``                             | Nested blocks inside ``&EXTE``            |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.hardness.directives``                          | ``hardness.directives``                          | Keyword/value lines inside ``&HARDNESS``  |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.hardness.subsections``                         | ``hardness.subsections``                         | Nested blocks inside ``&HARDNESS``        |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.info.directives``                              | ``info.directives``                              | Keyword/value lines inside ``&INFO``      |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.info.subsections``                             | ``info.subsections``                             | Nested blocks inside ``&INFO``            |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.linres.directives``                            | ``linres.directives``                            | Keyword/value lines inside ``&LINRES``    |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.linres.subsections``                           | ``linres.subsections``                           | Nested blocks inside ``&LINRES``          |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.molstates.directives``                         | ``molstates.directives``                         | Keyword/value lines inside ``&MOLSTATES`` |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.molstates.subsections``                        | ``molstates.subsections``                        | Nested blocks inside ``&MOLSTATES``       |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.mts.directives``                               | ``mts.directives``                               | Keyword/value lines inside ``&MTS``       |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.mts.subsections``                              | ``mts.subsections``                              | Nested blocks inside ``&MTS``             |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.nlcc.directives``                              | ``nlcc.directives``                              | Keyword/value lines inside ``&NLCC``      |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.nlcc.subsections``                             | ``nlcc.subsections``                             | Nested blocks inside ``&NLCC``            |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.path.directives``                              | ``path.directives``                              | Keyword/value lines inside ``&PATH``      |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.path.subsections``                             | ``path.subsections``                             | Nested blocks inside ``&PATH``            |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.pimd.directives``                              | ``pimd.directives``                              | Keyword/value lines inside ``&PIMD``      |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.pimd.subsections``                             | ``pimd.subsections``                             | Nested blocks inside ``&PIMD``            |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.potential.directives``                         | ``potential.directives``                         | Keyword/value lines inside ``&POTENTIAL`` |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.potential.subsections``                        | ``potential.subsections``                        | Nested blocks inside ``&POTENTIAL``       |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.prop.directives``                              | ``prop.directives``                              | Keyword/value lines inside ``&PROP``      |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.prop.subsections``                             | ``prop.subsections``                             | Nested blocks inside ``&PROP``            |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.ptddft.directives``                            | ``ptddft.directives``                            | Keyword/value lines inside ``&PTDDFT``    |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.ptddft.subsections``                           | ``ptddft.subsections``                           | Nested blocks inside ``&PTDDFT``          |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.resp.directives``                              | ``resp.directives``                              | Keyword/value lines inside ``&RESP``      |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.resp.subsections``                             | ``resp.subsections``                             | Nested blocks inside ``&RESP``            |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.tddft.directives``                             | ``tddft.directives``                             | Keyword/value lines inside ``&TDDFT``     |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.tddft.subsections``                            | ``tddft.subsections``                            | Nested blocks inside ``&TDDFT``           |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.vdw.directives``                               | ``vdw.directives``                               | Keyword/value lines inside ``&VDW``       |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.vdw.subsections``                              | ``vdw.subsections``                              | Nested blocks inside ``&VDW``, including  |
|                                                                       |                                                  | ``EMPIRICAL CORRECTION``                  |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.vectors.directives``                           | ``vectors.directives``                           | Keyword/value lines inside ``&VECTORS``   |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.vectors.subsections``                          | ``vectors.subsections``                          | Nested blocks inside ``&VECTORS``         |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.wavefunction.directives``                      | ``wavefunction.directives``                      | Keyword/value lines inside                |
|                                                                       |                                                  | ``&WAVEFUNCTION``                         |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.wavefunction.subsections``                     | ``wavefunction.subsections``                     | Nested blocks inside ``&WAVEFUNCTION``    |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.atoms.pseudopotentials``                       | ``atoms.pseudopotentials``                       | Pseudopotential entries grouped with      |
|                                                                       |                                                  | ``ForceInput`` coordinates                |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.atoms.directives``                             | ``atoms.directives``                             | Additional non-coordinate ``&ATOMS``      |
|                                                                       |                                                  | keyword/value lines                       |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.set.key``                                      | ``set.key``                                      | Dotted ``SECTION.KEYWORD`` merge target   |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.set.value``                                    | ``set.value``                                    | Optional value emitted under the dotted   |
|                                                                       |                                                  | key                                       |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+
| ``params.inputSections.raw``                                          | ``raw``                                          | Full section text inserted as-is          |
+-----------------------------------------------------------------------+--------------------------------------------------+-------------------------------------------+

Typed ``&CPMD`` controls
========================

``CPMDCpmdSection`` exposes these driver controls directly;
parser-accepted escape-hatch rows use ``directives``. The feature IDs
are the stable discovery keys returned by ``cpmdc_feature_table()``.

+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| Catalog feature ID                                       | Parameter feature ID                                              | Field                                   | Deck keyword                                |
+==========================================================+===================================================================+=========================================+=============================================+
| ``catalog.cpmd.OPTIMIZE_WAVEFUNCTION``                   | ``params.inputSections.cpmd.optimizeWavefunction``                | ``optimizeWavefunction``                | ``OPTIMIZE WAVEFUNCTION``                   |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.OPTIMIZE_GEOMETRY``                       | ``params.inputSections.cpmd.optimizeGeometry``                    | ``optimizeGeometry``                    | ``OPTIMIZE GEOMETRY``                       |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.OPTIMIZE_GEOMETRY_OPTIONS``               | ``params.inputSections.cpmd.optimizeGeometryOptions``             | ``optimizeGeometryOptions``             | ``OPTIMIZE GEOMETRY``                       |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.OPTIMIZE_GEOMETRY_SAMPLE``                | ``params.inputSections.cpmd.optimizeGeometrySample``              | ``optimizeGeometrySample``              | ``OPTIMIZE GEOMETRY SAMPLE``                |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.OPTIMIZE_COMBINED``                       | ``params.inputSections.cpmd.optimizeCombinedOptions``             | ``optimizeCombinedOptions``             | ``OPTIMIZE COMBINED``                       |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.OPTIMIZE_COMBINED_SAMPLE``                | ``params.inputSections.cpmd.optimizeCombinedSample``              | ``optimizeCombinedSample``              | ``OPTIMIZE COMBINED SAMPLE``                |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.MOLECULAR_DYNAMICS``                      | ``params.inputSections.cpmd.molecularDynamics``                   | ``molecularDynamics``                   | ``MOLECULAR DYNAMICS``                      |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.MOLECULAR_DYNAMICS_CP``                   | ``params.inputSections.cpmd.molecularDynamicsCp``                 | ``molecularDynamicsCp``                 | ``MOLECULAR DYNAMICS CP``                   |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.MOLECULAR_DYNAMICS_BO``                   | ``params.inputSections.cpmd.molecularDynamicsBo``                 | ``molecularDynamicsBo``                 | ``MOLECULAR DYNAMICS BO``                   |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.MOLECULAR_DYNAMICS_EH``                   | ``params.inputSections.cpmd.molecularDynamicsEh``                 | ``molecularDynamicsEh``                 | ``MOLECULAR DYNAMICS EH``                   |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.MOLECULAR_DYNAMICS_PT``                   | ``params.inputSections.cpmd.molecularDynamicsPt``                 | ``molecularDynamicsPt``                 | ``MOLECULAR DYNAMICS PT``                   |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.MOLECULAR_DYNAMICS_CLASSICAL``            | ``params.inputSections.cpmd.molecularDynamicsClassical``          | ``molecularDynamicsClassical``          | ``MOLECULAR DYNAMICS CLASSICAL``            |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.MOLECULAR_DYNAMICS_FILE``                 | ``params.inputSections.cpmd.molecularDynamicsFile``               | ``molecularDynamicsFile``               | ``MOLECULAR DYNAMICS FILE``                 |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.MOLECULAR_DYNAMICS_FILE_OPTIONS``         | ``params.inputSections.cpmd.molecularDynamicsFileOptions``        | ``molecularDynamicsFileOptions``        | ``MOLECULAR DYNAMICS FILE options``         |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.MOLECULAR_DYNAMICS_BD``                   | ``params.inputSections.cpmd.molecularDynamicsBdTrajectories``     | ``molecularDynamicsBdTrajectories``     | ``MOLECULAR DYNAMICS BD``                   |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.PARRINELLO_RAHMAN``                       | ``params.inputSections.cpmd.parrinelloRahmanOptions``             | ``parrinelloRahmanOptions``             | ``PARRINELLO-RAHMAN``                       |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.CHEBY``                                   | ``params.inputSections.cpmd.cheby``                               | ``cheby``                               | ``CHEBY``                                   |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.CAYLEY``                                  | ``params.inputSections.cpmd.cayley``                              | ``cayley``                              | ``CAYLEY``                                  |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.RUNGE_KUTTA``                             | ``params.inputSections.cpmd.rungeKutta``                          | ``rungeKutta``                          | ``RUNGE-KUTTA``                             |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.FORCEMATCH``                              | ``params.inputSections.cpmd.forceMatch``                          | ``forceMatch``                          | ``FORCEMATCH``                              |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.EMASS``                                   | ``params.inputSections.cpmd.electronMass``                        | ``electronMass``                        | ``EMASS``                                   |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.TIMESTEP``                                | ``params.inputSections.cpmd.timestep``                            | ``timestep``                            | ``TIMESTEP``                                |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.MAXRUNTIME``                              | ``params.inputSections.cpmd.maxRuntime``                          | ``maxRuntime``                          | ``MAXRUNTIME``                              |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.TIMESTEP_ELECTRONS``                      | ``params.inputSections.cpmd.timestepElectrons``                   | ``timestepElectrons``                   | ``TIMESTEP ELECTRONS``                      |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.TIMESTEP_IONS``                           | ``params.inputSections.cpmd.timestepIons``                        | ``timestepIons``                        | ``TIMESTEP IONS``                           |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.CMASS``                                   | ``params.inputSections.cpmd.cellMass``                            | ``cellMass``                            | ``CMASS``                                   |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.TEMPERATURE_ELECTRON``                    | ``params.inputSections.cpmd.temperatureElectron``                 | ``temperatureElectron``                 | ``TEMPERATURE ELECTRON``                    |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.TEMPERATURE``                             | ``params.inputSections.cpmd.temperature``                         | ``temperature``                         | ``TEMPERATURE``                             |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.TEMPERATURE_RAMP``                        | ``params.inputSections.cpmd.temperatureRamp``                     | ``temperatureRamp``                     | ``TEMPERATURE RAMP``                        |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.TEMPERATURE_RAMP_TIME``                   | ``params.inputSections.cpmd.temperatureRampTime``                 | ``temperatureRampTime``                 | ``TEMPERATURE RAMP payload``                |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.TEMPERATURE_RAMP_RATE``                   | ``params.inputSections.cpmd.temperatureRampRate``                 | ``temperatureRampRate``                 | ``TEMPERATURE RAMP payload``                |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.RESCALE_VELOCITIES``                      | ``params.inputSections.cpmd.rescaleOldVelocities``                | ``rescaleOldVelocities``                | ``RESCALE VELOCITIES``                      |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.RESCALE_OLD_VELOCITIES``                  | ``params.inputSections.cpmd.rescaleOldVelocities``                | ``rescaleOldVelocities``                | ``RESCALE OLD VELOCITIES``                  |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.REVERSE_VELOCITIES``                      | ``params.inputSections.cpmd.reverseVelocities``                   | ``reverseVelocities``                   | ``REVERSE VELOCITIES``                      |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.SUBTRACT_COMVEL``                         | ``params.inputSections.cpmd.subtractComVelocity``                 | ``subtractComVelocity``                 | ``SUBTRACT COMVEL``                         |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.SUBTRACT_ROTVEL``                         | ``params.inputSections.cpmd.subtractRotVelocity``                 | ``subtractRotVelocity``                 | ``SUBTRACT ROTVEL``                         |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.SUBTRACT``                                | ``params.inputSections.cpmd.subtractComVelocity``                 | ``subtractComVelocity``                 | ``SUBTRACT``                                |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.PRNGSEED``                                | ``params.inputSections.cpmd.prngSeed``                            | ``prngSeed``                            | ``PRNGSEED``                                |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.TEMPCONTROL_IONS``                        | ``params.inputSections.cpmd.tempControlIons``                     | ``tempControlIons``                     | ``TEMPCONTROL IONS``                        |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.TEMPCONTROL_ELECTRONS``                   | ``params.inputSections.cpmd.tempControlElectrons``                | ``tempControlElectrons``                | ``TEMPCONTROL ELECTRONS``                   |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.TEMPCONTROL_CELL``                        | ``params.inputSections.cpmd.tempControlCell``                     | ``tempControlCell``                     | ``TEMPCONTROL CELL``                        |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.TEMPCONTROL``                             | ``params.inputSections.cpmd.tempControlIons``                     | ``tempControlIons``                     | ``TEMPCONTROL``                             |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.BERENDSEN_IONS``                          | ``params.inputSections.cpmd.berendsenIons``                       | ``berendsenIons``                       | ``BERENDSEN IONS``                          |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.BERENDSEN_ELECTRONS``                     | ``params.inputSections.cpmd.berendsenElectrons``                  | ``berendsenElectrons``                  | ``BERENDSEN ELECTRONS``                     |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.BERENDSEN_CELL``                          | ``params.inputSections.cpmd.berendsenCell``                       | ``berendsenCell``                       | ``BERENDSEN CELL``                          |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.MAXSTEP``                                 | ``params.inputSections.cpmd.maxStep``                             | ``maxStep``                             | ``MAXSTEP``                                 |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.MAXITER``                                 | ``params.inputSections.cpmd.maxIter``                             | ``maxIter``                             | ``MAXITER``                                 |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.CONVERGENCE_ORBITALS``                    | ``params.inputSections.cpmd.convergenceOrbitals``                 | ``convergenceOrbitals``                 | ``CONVERGENCE ORBITALS``                    |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.CONVERGENCE_GEOMETRY``                    | ``params.inputSections.cpmd.convergenceGeometry``                 | ``convergenceGeometry``                 | ``CONVERGENCE GEOMETRY``                    |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.CONVERGENCE``                             | ``params.inputSections.cpmd.convergenceOrbitals``                 | ``convergenceOrbitals``                 | ``CONVERGENCE``                             |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.CONVERGENCE_CELL``                        | ``params.inputSections.cpmd.convergenceCell``                     | ``convergenceCell``                     | ``CONVERGENCE CELL``                        |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.CONVERGENCE_ADAPT``                       | ``params.inputSections.cpmd.convergenceAdapt``                    | ``convergenceAdapt``                    | ``CONVERGENCE ADAPT``                       |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.CONVERGENCE_ENERGY``                      | ``params.inputSections.cpmd.convergenceEnergy``                   | ``convergenceEnergy``                   | ``CONVERGENCE ENERGY``                      |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.CONVERGENCE_CALFOR``                      | ``params.inputSections.cpmd.convergenceCalfor``                   | ``convergenceCalfor``                   | ``CONVERGENCE CALFOR``                      |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.CONVERGENCE_RELAX``                       | ``params.inputSections.cpmd.convergenceRelax``                    | ``convergenceRelax``                    | ``CONVERGENCE RELAX``                       |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.CONVERGENCE_RHOFIX``                      | ``params.inputSections.cpmd.convergenceRhofix``                   | ``convergenceRhofix``                   | ``CONVERGENCE RHOFIX``                      |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.CONVERGENCE_INITIAL``                     | ``params.inputSections.cpmd.convergenceInitial``                  | ``convergenceInitial``                  | ``CONVERGENCE INITIAL``                     |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.CONVERGENCE_CONSTRAINT``                  | ``params.inputSections.cpmd.convergenceConstraint``               | ``convergenceConstraint``               | ``CONVERGENCE CONSTRAINT``                  |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.NOSE``                                    | ``params.inputSections.cpmd.nose``                                | ``nose``                                | ``NOSE``                                    |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.NOSE_IONS``                               | ``params.inputSections.cpmd.noseIons``                            | ``noseIons``                            | ``NOSE IONS``                               |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.NOSE_IONS_OPTIONS``                       | ``params.inputSections.cpmd.noseIonsOptions``                     | ``noseIonsOptions``                     | ``NOSE IONS options``                       |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.NOSE_ELECTRONS``                          | ``params.inputSections.cpmd.noseElectrons``                       | ``noseElectrons``                       | ``NOSE ELECTRONS``                          |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.NOSE_IONS_PAYLOAD``                       | ``params.inputSections.cpmd.noseIonsThermostat``                  | ``noseIonsThermostat``                  | ``NOSE IONS``                               |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.NOSE_IONS_LOCAL``                         | ``params.inputSections.cpmd.noseIonsLocalThermostatCount``        | ``noseIonsLocalThermostatCount``        | ``NOSE IONS LOCAL``                         |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.NOSE_IONS_LOCAL_T0``                      | ``params.inputSections.cpmd.noseIonsLocalT0``                     | ``noseIonsLocalT0``                     | ``NOSE IONS LOCAL T0``                      |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.NOSE_IONS_LOCAL_THERMOSTATS``             | ``params.inputSections.cpmd.noseIonsLocalThermostats``            | ``noseIonsLocalThermostats``            | ``NOSE IONS LOCAL thermostat payload``      |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.NOSE_IONS_LOCAL_RANGES``                  | ``params.inputSections.cpmd.noseIonsLocalRangeCount``             | ``noseIonsLocalRangeCount``             | ``NOSE IONS LOCAL range count``             |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.NOSE_IONS_LOCAL_RANGES``                  | ``params.inputSections.cpmd.noseIonsLocalRanges``                 | ``noseIonsLocalRanges``                 | ``NOSE IONS LOCAL range payload``           |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.NOSE_IONS_CAFES``                         | ``params.inputSections.cpmd.noseIonsCafesGroupCount``             | ``noseIonsCafesGroupCount``             | ``NOSE IONS CAFES``                         |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.NOSE_IONS_CAFES_GROUPS``                  | ``params.inputSections.cpmd.noseIonsCafesGroups``                 | ``noseIonsCafesGroups``                 | ``NOSE IONS CAFES group payload``           |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.NOSE_ELECTRONS_PAYLOAD``                  | ``params.inputSections.cpmd.noseElectronsThermostat``             | ``noseElectronsThermostat``             | ``NOSE ELECTRONS``                          |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.NOSE_CELL``                               | ``params.inputSections.cpmd.noseCellThermostat``                  | ``noseCellThermostat``                  | ``NOSE CELL``                               |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.NOSE_CELL_PAYLOAD``                       | ``params.inputSections.cpmd.noseCellThermostat``                  | ``noseCellThermostat``                  | ``NOSE CELL``                               |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.NOSE_PARAMETER``                          | ``params.inputSections.cpmd.noseParameters``                      | ``noseParameters``                      | ``NOSE PARAMETER``                          |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.NOSE_PARAMETERS``                         | ``params.inputSections.cpmd.noseParameters``                      | ``noseParameters``                      | ``NOSE PARAMETERS``                         |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.USE_IN_STREAM``                           | ``params.inputSections.cpmd.useInStream``                         | ``useInStream``                         | ``USE_IN_STREAM``                           |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.USE_OUT_STREAM``                          | ``params.inputSections.cpmd.useOutStream``                        | ``useOutStream``                        | ``USE_OUT_STREAM``                          |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.USE_CUBLAS``                              | ``params.inputSections.cpmd.useCublas``                           | ``useCublas``                           | ``USE_CUBLAS``                              |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.USE_CUFFT``                               | ``params.inputSections.cpmd.useCufft``                            | ``useCufft``                            | ``USE_CUFFT``                               |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.BLAS_N_STREAMS_PER_DEVICE``               | ``params.inputSections.cpmd.blasNStreamsPerDevice``               | ``blasNStreamsPerDevice``               | ``BLAS_N_STREAMS_PER_DEVICE``               |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.BLAS_N_DEVICES_PER_TASK``                 | ``params.inputSections.cpmd.blasNDevicesPerTask``                 | ``blasNDevicesPerTask``                 | ``BLAS_N_DEVICES_PER_TASK``                 |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.FFT_N_STREAMS_PER_DEVICE``                | ``params.inputSections.cpmd.fftNStreamsPerDevice``                | ``fftNStreamsPerDevice``                | ``FFT_N_STREAMS_PER_DEVICE``                |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.FFT_N_DEVICES_PER_TASK``                  | ``params.inputSections.cpmd.fftNDevicesPerTask``                  | ``fftNDevicesPerTask``                  | ``FFT_N_DEVICES_PER_TASK``                  |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.USE_MPI_IO``                              | ``params.inputSections.cpmd.useMpiIo``                            | ``useMpiIo``                            | ``USE_MPI_IO``                              |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.TRACE``                                   | ``params.inputSections.cpmd.traceOptions``                        | ``traceOptions``                        | ``TRACE``                                   |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.TRACE_PROCEDURE``                         | ``params.inputSections.cpmd.traceProcedure``                      | ``traceProcedure``                      | ``TRACE_PROCEDURE``                         |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.TRACE_MAX_DEPTH``                         | ``params.inputSections.cpmd.traceMaxDepth``                       | ``traceMaxDepth``                       | ``TRACE_MAX_DEPTH``                         |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.TRACE_MAX_CALLS``                         | ``params.inputSections.cpmd.traceMaxCalls``                       | ``traceMaxCalls``                       | ``TRACE_MAX_CALLS``                         |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.BERENDSEN``                               | ``params.inputSections.cpmd.berendsen``                           | ``berendsen``                           | ``BERENDSEN``                               |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.LANGEVIN``                                | ``params.inputSections.cpmd.langevin``                            | ``langevin``                            | ``LANGEVIN``                                |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.ANNEALING``                               | ``params.inputSections.cpmd.annealing``                           | ``annealing``                           | ``ANNEALING``                               |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.ANNEALING_IONS``                          | ``params.inputSections.cpmd.annealingIons``                       | ``annealingIons``                       | ``ANNEALING IONS``                          |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.ANNEALING_ELECTRONS``                     | ``params.inputSections.cpmd.annealingElectrons``                  | ``annealingElectrons``                  | ``ANNEALING ELECTRONS``                     |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.ANNEALING_CELL``                          | ``params.inputSections.cpmd.annealingCell``                       | ``annealingCell``                       | ``ANNEALING CELL``                          |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.DAMPING``                                 | ``params.inputSections.cpmd.damping``                             | ``damping``                             | ``DAMPING``                                 |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.DAMPING_IONS``                            | ``params.inputSections.cpmd.dampingIons``                         | ``dampingIons``                         | ``DAMPING IONS``                            |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.DAMPING_ELECTRONS``                       | ``params.inputSections.cpmd.dampingElectrons``                    | ``dampingElectrons``                    | ``DAMPING ELECTRONS``                       |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.DAMPING_CELL``                            | ``params.inputSections.cpmd.dampingCell``                         | ``dampingCell``                         | ``DAMPING CELL``                            |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.HESSIAN``                                 | ``params.inputSections.cpmd.hessian``                             | ``hessian``                             | ``HESSIAN``                                 |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.PROJECT``                                 | ``params.inputSections.cpmd.project``                             | ``project``                             | ``PROJECT``                                 |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.STRESS_TENSOR``                           | ``params.inputSections.cpmd.stressTensorSample``                  | ``stressTensorSample``                  | ``STRESS TENSOR``                           |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.STRESS_TENSOR_SAMPLE``                    | ``params.inputSections.cpmd.stressTensorSample``                  | ``stressTensorSample``                  | ``STRESS TENSOR``                           |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.STRESS_TENSOR_VIRIAL``                    | ``params.inputSections.cpmd.stressTensorVirial``                  | ``stressTensorVirial``                  | ``STRESS TENSOR VIRIAL``                    |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.CLASSTRESS``                              | ``params.inputSections.cpmd.classStressSample``                   | ``classStressSample``                   | ``CLASSTRESS``                              |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.QUENCH``                                  | ``params.inputSections.cpmd.quench``                              | ``quench``                              | ``QUENCH``                                  |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.RATTLE``                                  | ``params.inputSections.cpmd.rattle``                              | ``rattle``                              | ``RATTLE``                                  |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.SHAKE``                                   | ``params.inputSections.cpmd.shake``                               | ``shake``                               | ``SHAKE``                                   |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.CONSTRAINT``                              | ``params.inputSections.cpmd.constraint``                          | ``constraint``                          | ``CONSTRAINT``                              |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.TROTTER``                                 | ``params.inputSections.cpmd.trotter``                             | ``trotter``                             | ``TROTTER``                                 |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.RESTART``                                 | ``params.inputSections.cpmd.restart``                             | ``restart``                             | ``RESTART``                                 |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.RESTART_WAVEFUNCTION``                    | ``params.inputSections.cpmd.restartWavefunction``                 | ``restartWavefunction``                 | ``RESTART WAVEFUNCTION``                    |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.TRAJECTORY``                              | ``params.inputSections.cpmd.trajectory``                          | ``trajectory``                          | ``TRAJECTORY``                              |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.ISOLATED_MOLECULE``                       | ``params.inputSections.cpmd.isolatedMolecule``                    | ``isolatedMolecule``                    | ``ISOLATED MOLECULE``                       |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.PRINT``                                   | ``params.inputSections.cpmd.printOptions``                        | ``printOptions``                        | ``PRINT``                                   |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.STRUCTURE``                               | ``params.inputSections.cpmd.structureOptions``                    | ``structureOptions``                    | ``STRUCTURE``                               |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.STRUCTURE_SELECTION``                     | ``params.inputSections.cpmd.structureSelection``                  | ``structureSelection``                  | ``STRUCTURE SELECT payload``                |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.STORE``                                   | ``params.inputSections.cpmd.storeOptions``                        | ``storeOptions``                        | ``STORE``                                   |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.STORE_SELECTION``                         | ``params.inputSections.cpmd.storeSelection``                      | ``storeSelection``                      | ``STORE``                                   |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.STORE_INTERVAL``                          | ``params.inputSections.cpmd.storeInterval``                       | ``storeInterval``                       | ``STORE interval``                          |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.STORE_SC_INTERVAL``                       | ``params.inputSections.cpmd.storeSelfConsistentInterval``         | ``storeSelfConsistentInterval``         | ``STORE SC= interval``                      |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.STORE_OFF``                               | ``params.inputSections.cpmd.storeOffSelection``                   | ``storeOffSelection``                   | ``STORE OFF``                               |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.RESTFILE``                                | ``params.inputSections.cpmd.restFileCount``                       | ``restFileCount``                       | ``RESTFILE``                                |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.RESTFILE_SAMPLE``                         | ``params.inputSections.cpmd.restFileSample``                      | ``restFileSample``                      | ``RESTFILE SAMPLE``                         |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.TRAJECTORY_OPTIONS``                      | ``params.inputSections.cpmd.trajectoryOptions``                   | ``trajectoryOptions``                   | ``TRAJECTORY``                              |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.TRAJECTORY_SAMPLE``                       | ``params.inputSections.cpmd.trajectorySample``                    | ``trajectorySample``                    | ``TRAJECTORY SAMPLE``                       |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.TRAJECTORY_RANGE``                        | ``params.inputSections.cpmd.trajectoryRange``                     | ``trajectoryRange``                     | ``TRAJECTORY RANGE``                        |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.MOVIE_SAMPLE``                            | ``params.inputSections.cpmd.movieSample``                         | ``movieSample``                         | ``MOVIE SAMPLE``                            |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.MOVIE_OFF``                               | ``params.inputSections.cpmd.movieOff``                            | ``movieOff``                            | ``MOVIE OFF``                               |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.MOVIE``                                   | ``params.inputSections.cpmd.movieSample``                         | ``movieSample``                         | ``MOVIE``                                   |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.ENERGYBANDS``                             | ``params.inputSections.cpmd.energyBands``                         | ``energyBands``                         | ``ENERGYBANDS``                             |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.EXTERNAL_POTENTIAL``                      | ``params.inputSections.cpmd.externalPotential``                   | ``externalPotential``                   | ``EXTERNAL POTENTIAL``                      |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.EXTERNAL_POTENTIAL_ADD``                  | ``params.inputSections.cpmd.externalPotentialAdd``                | ``externalPotentialAdd``                | ``EXTERNAL POTENTIAL ADD``                  |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.ELECTROSTATIC_POTENTIAL``                 | ``params.inputSections.cpmd.electrostaticPotential``              | ``electrostaticPotential``              | ``ELECTROSTATIC POTENTIAL``                 |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.ELECTROSTATIC_POTENTIAL_SAMPLE``          | ``params.inputSections.cpmd.electrostaticPotentialSample``        | ``electrostaticPotentialSample``        | ``ELECTROSTATIC POTENTIAL SAMPLE``          |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.DIPOLE_DYNAMICS``                         | ``params.inputSections.cpmd.dipoleDynamics``                      | ``dipoleDynamics``                      | ``DIPOLE DYNAMICS``                         |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.DIPOLE_DYNAMICS_SAMPLE``                  | ``params.inputSections.cpmd.dipoleDynamicsSample``                | ``dipoleDynamicsSample``                | ``DIPOLE DYNAMICS SAMPLE``                  |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.DIPOLE_DYNAMICS_WANNIER``                 | ``params.inputSections.cpmd.dipoleDynamicsWannier``               | ``dipoleDynamicsWannier``               | ``DIPOLE DYNAMICS WANNIER``                 |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.RHOOUT``                                  | ``params.inputSections.cpmd.rhoOut``                              | ``rhoOut``                              | ``RHOOUT``                                  |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.RHOOUT_SAMPLE``                           | ``params.inputSections.cpmd.rhoOutSample``                        | ``rhoOutSample``                        | ``RHOOUT SAMPLE``                           |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.RHOOUT_BANDS_COUNT``                      | ``params.inputSections.cpmd.rhoOutBandsCount``                    | ``rhoOutBandsCount``                    | ``RHOOUT BANDS count``                      |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.RHOOUT_BANDS``                            | ``params.inputSections.cpmd.rhoOutBands``                         | ``rhoOutBands``                         | ``RHOOUT BANDS``                            |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.ELF``                                     | ``params.inputSections.cpmd.elf``                                 | ``elf``                                 | ``ELF``                                     |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.ELF_PARAMETERS``                          | ``params.inputSections.cpmd.elfParameters``                       | ``elfParameters``                       | ``ELF PARAMETER``                           |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.WANNIER_PARAMETER``                       | ``params.inputSections.cpmd.wannierParameters``                   | ``wannierParameters``                   | ``WANNIER PARAMETER``                       |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.WANNIER_OPTIMIZATION``                    | ``params.inputSections.cpmd.wannierOptimization``                 | ``wannierOptimization``                 | ``WANNIER OPTIMIZATION``                    |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.WANNIER_NPROC``                           | ``params.inputSections.cpmd.wannierNproc``                        | ``wannierNproc``                        | ``WANNIER NPROC``                           |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.WANNIER_RELOCALIZE_IN_SCF``               | ``params.inputSections.cpmd.wannierRelocalizeInScf``              | ``wannierRelocalizeInScf``              | ``WANNIER RELOCALIZE_IN_SCF``               |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.WANNIER_RECOMPUTE_DIPOLE_MATRICES_EVERY`` | ``params.inputSections.cpmd.wannierRecomputeDipoleMatricesEvery`` | ``wannierRecomputeDipoleMatricesEvery`` | ``WANNIER RECOMPUTE_DIPOLE_MATRICES_EVERY`` |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.WANNIER_RELOCALIZE_EVERY``                | ``params.inputSections.cpmd.wannierRelocalizeEvery``              | ``wannierRelocalizeEvery``              | ``WANNIER RELOCALIZE_EVERY``                |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.WANNIER_TYPE``                            | ``params.inputSections.cpmd.wannierType``                         | ``wannierType``                         | ``WANNIER TYPE``                            |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.WANNIER_REFERENCE``                       | ``params.inputSections.cpmd.wannierReference``                    | ``wannierReference``                    | ``WANNIER REFERENCE``                       |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.WANNIER_SERIAL``                          | ``params.inputSections.cpmd.wannierSerial``                       | ``wannierSerial``                       | ``WANNIER SERIAL``                          |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.WANNIER_DOS``                             | ``params.inputSections.cpmd.wannierDos``                          | ``wannierDos``                          | ``WANNIER DOS``                             |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.WANNIER_MOLECULAR``                       | ``params.inputSections.cpmd.wannierMolecular``                    | ``wannierMolecular``                    | ``WANNIER MOLECULAR``                       |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.WANNIER_WFNOUT_OPTIONS``                  | ``params.inputSections.cpmd.wannierWfnOutOptions``                | ``wannierWfnOutOptions``                | ``WANNIER WFNOUT``                          |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.WANNIER_WFNOUT_PAYLOAD``                  | ``params.inputSections.cpmd.wannierWfnOutPayload``                | ``wannierWfnOutPayload``                | ``WANNIER WFNOUT payload``                  |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.PARA_USE_MPI_IN_PLACE``                   | ``params.inputSections.cpmd.paraUseMpiInPlace``                   | ``paraUseMpiInPlace``                   | ``PARA_USE_MPI_IN_PLACE``                   |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.PARA_BUFF_SIZE``                          | ``params.inputSections.cpmd.paraBuffSize``                        | ``paraBuffSize``                        | ``PARA_BUFF_SIZE``                          |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.PARA_STACK_BUFF_SIZE``                    | ``params.inputSections.cpmd.paraStackBuffSize``                   | ``paraStackBuffSize``                   | ``PARA_STACK_BUFF_SIZE``                    |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.COMPRESS``                                | ``params.inputSections.cpmd.compress``                            | ``compress``                            | ``COMPRESS``                                |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.MEMORY``                                  | ``params.inputSections.cpmd.memory``                              | ``memory``                              | ``MEMORY``                                  |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.CHECK_MEMORY``                            | ``params.inputSections.cpmd.checkMemory``                         | ``checkMemory``                         | ``CHECK MEMORY``                            |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.REAL_SPACE_WFN``                          | ``params.inputSections.cpmd.realSpaceWfn``                        | ``realSpaceWfn``                        | ``REAL SPACE WFN``                          |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.REALSPACE_WFN``                           | ``params.inputSections.cpmd.realSpaceWfn``                        | ``realSpaceWfn``                        | ``REALSPACE WFN``                           |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.REAL_SPACE_WFN_KEEP``                     | ``params.inputSections.cpmd.realSpaceWfnKeep``                    | ``realSpaceWfnKeep``                    | ``REAL SPACE WFN KEEP``                     |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.REAL_SPACE_WFN_SIZE``                     | ``params.inputSections.cpmd.realSpaceWfnSize``                    | ``realSpaceWfnSize``                    | ``REAL SPACE WFN SIZE``                     |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.SPLINE``                                  | ``params.inputSections.cpmd.splineOptions``                       | ``splineOptions``                       | ``SPLINE``                                  |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.SPLINE_OPTIONS``                          | ``params.inputSections.cpmd.splineOptions``                       | ``splineOptions``                       | ``SPLINE``                                  |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.SPLINE_POINTS``                           | ``params.inputSections.cpmd.splinePoints``                        | ``splinePoints``                        | ``SPLINE POINTS``                           |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.SPLINE_RANGE``                            | ``params.inputSections.cpmd.splineRange``                         | ``splineRange``                         | ``SPLINE RANGE``                            |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.FINITE_DIFFERENCE``                       | ``params.inputSections.cpmd.finiteDifferences``                   | ``finiteDifferences``                   | ``FINITE DIFFERENCE``                       |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.FINITE_DIFFERENCES``                      | ``params.inputSections.cpmd.finiteDifferences``                   | ``finiteDifferences``                   | ``FINITE DIFFERENCES``                      |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.TASKGROUPS``                              | ``params.inputSections.cpmd.taskGroups``                          | ``taskGroups``                          | ``TASKGROUPS``                              |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.TASKGROUPS_COUNT``                        | ``params.inputSections.cpmd.taskGroupsCount``                     | ``taskGroupsCount``                     | ``TASKGROUPS count``                        |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.DISTRIBUTE_FNL``                          | ``params.inputSections.cpmd.distributeFnl``                       | ``distributeFnl``                       | ``DISTRIBUTE FNL``                          |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.CP_GROUPS``                               | ``params.inputSections.cpmd.cpGroups``                            | ``cpGroups``                            | ``CP_GROUPS``                               |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.FILE_FUSION``                             | ``params.inputSections.cpmd.fileFusionPayload``                   | ``fileFusionPayload``                   | ``FILE FUSION``                             |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.FILE_MERGE``                              | ``params.inputSections.cpmd.fileMergePayload``                    | ``fileMergePayload``                    | ``FILE MERGE``                              |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.FILE_SEPARATION``                         | ``params.inputSections.cpmd.fileSeparationPayload``               | ``fileSeparationPayload``               | ``FILE SEPARATION``                         |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.NO_RESET``                                | ``params.inputSections.cpmd.noReset``                             | ``noReset``                             | ``NO_RESET``                                |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.FILEPATH``                                | ``params.inputSections.cpmd.filePath``                            | ``filePath``                            | ``FILEPATH``                                |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.BENCHMARK``                               | ``params.inputSections.cpmd.benchmark``                           | ``benchmark``                           | ``BENCHMARK``                               |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.MIRROR``                                  | ``params.inputSections.cpmd.mirror``                              | ``mirror``                              | ``MIRROR``                                  |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.SHIFT_POTENTIAL``                         | ``params.inputSections.cpmd.shiftPotential``                      | ``shiftPotential``                      | ``SHIFT POTENTIAL``                         |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.GLOCALIZATION_PARAMETERS``                | ``params.inputSections.cpmd.glocalizationParameters``             | ``glocalizationParameters``             | ``GLOCALIZATION PARAMETERS``                |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.GLOCALIZATION_OPTIMIZATION``              | ``params.inputSections.cpmd.glocalizationOptimization``           | ``glocalizationOptimization``           | ``GLOCALIZATION OPTIMIZATION``              |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.GFUNCTIONAL_TYPE``                        | ``params.inputSections.cpmd.gfunctionalType``                     | ``gfunctionalType``                     | ``GFUNCTIONAL TYPE``                        |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.SPREAD_RSPACE``                           | ``params.inputSections.cpmd.spreadRspace``                        | ``spreadRspace``                        | ``SPREAD RSPACE=``                          |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.G_UNITARITY_OPTIONS``                     | ``params.inputSections.cpmd.gUnitarityOptions``                   | ``gUnitarityOptions``                   | ``PIPPO``                                   |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.STEP_TUNING``                             | ``params.inputSections.cpmd.stepTuning``                          | ``stepTuning``                          | ``STEP TUNING``                             |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.G_ANTISYM``                               | ``params.inputSections.cpmd.gAntisym``                            | ``gAntisym``                            | ``G_ANTISYM``                               |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.G_ANTISYM_PENALTY``                       | ``params.inputSections.cpmd.gAntisymPenalty``                     | ``gAntisymPenalty``                     | ``G_ANTISYM PENALTY``                       |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.G_KICK``                                  | ``params.inputSections.cpmd.gKick``                               | ``gKick``                               | ``G_KICK``                                  |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.G_COMPLEX``                               | ``params.inputSections.cpmd.gComplex``                            | ``gComplex``                            | ``G_COMPLEX``                               |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.G_REAL``                                  | ``params.inputSections.cpmd.gReal``                               | ``gReal``                               | ``G_REAL``                                  |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.READ_MATRIX``                             | ``params.inputSections.cpmd.readMatrix``                          | ``readMatrix``                          | ``READ MATRIX``                             |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.G_STEP_TUNE``                             | ``params.inputSections.cpmd.gStepTune``                           | ``gStepTune``                           | ``G_STEP TUNE``                             |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.GLOC_WFNOUT``                             | ``params.inputSections.cpmd.glocWfnOutOptions``                   | ``glocWfnOutOptions``                   | ``GLOC WFNOUT``                             |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.GLOC_WFNOUT_PAYLOAD``                     | ``params.inputSections.cpmd.glocWfnOutPayload``                   | ``glocWfnOutPayload``                   | ``GLOC WFNOUT payload``                     |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.NO_GEO_CHECK``                            | ``params.inputSections.cpmd.noGeoCheck``                          | ``noGeoCheck``                          | ``NO_GEO_CHECK``                            |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.BROKEN``                                  | ``params.inputSections.cpmd.brokenSymmetry``                      | ``brokenSymmetry``                      | ``BROKEN``                                  |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.DISTRIBUTED_LINALG``                      | ``params.inputSections.cpmd.distributedLinalg``                   | ``distributedLinalg``                   | ``DISTRIBUTED LINALG``                      |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.LINALG_NEWORTHO``                         | ``params.inputSections.cpmd.linalgNewOrtho``                      | ``linalgNewOrtho``                      | ``LINALG NEWORTHO``                         |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.DISORTHO_BSIZE``                          | ``params.inputSections.cpmd.disorthoBlockSize``                   | ``disorthoBlockSize``                   | ``DISORTHO_BSIZE``                          |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.BLOCKSIZE_STATES``                        | ``params.inputSections.cpmd.statesBlockSize``                     | ``statesBlockSize``                     | ``BLOCKSIZE STATES``                        |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.ALLTOALL``                                | ``params.inputSections.cpmd.allToAllPrecision``                   | ``allToAllPrecision``                   | ``ALLTOALL``                                |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.GSHELL``                                  | ``params.inputSections.cpmd.gshell``                              | ``gshell``                              | ``GSHELL``                                  |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.LOCAL_POTENTIAL``                         | ``params.inputSections.cpmd.localPotential``                      | ``localPotential``                      | ``LOCAL POTENTIAL``                         |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.CENTER_MOLECULE_OFF``                     | ``params.inputSections.cpmd.centerMoleculeOff``                   | ``centerMoleculeOff``                   | ``CENTER MOLECULE OFF``                     |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.CENTER_MOLECULE_ON``                      | ``params.inputSections.cpmd.centerMoleculeOn``                    | ``centerMoleculeOn``                    | ``CENTER MOLECULE ON``                      |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.CENTER_MOLECULE``                         | ``params.inputSections.cpmd.centerMoleculeOn``                    | ``centerMoleculeOn``                    | ``CENTER MOLECULE``                         |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.DIIS``                                    | ``params.inputSections.cpmd.diis``                                | ``diis``                                | ``DIIS``                                    |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.ODIIS``                                   | ``params.inputSections.cpmd.odiis``                               | ``odiis``                               | ``ODIIS``                                   |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.PCG``                                     | ``params.inputSections.cpmd.pcg``                                 | ``pcg``                                 | ``PCG``                                     |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.DIAGONALIZATION``                         | ``params.inputSections.cpmd.diagonalization``                     | ``diagonalization``                     | ``DIAGONALIZATION``                         |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.FREE_ENERGY``                             | ``params.inputSections.cpmd.freeEnergy``                          | ``freeEnergy``                          | ``FREE-ENERGY``                             |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.FREE_ENERGY_FUNCTIONAL``                  | ``params.inputSections.cpmd.freeEnergy``                          | ``freeEnergy``                          | ``FREE ENERGY FUNCTIONAL``                  |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.INTERFACE``                               | ``params.inputSections.cpmd.interface``                           | ``interface``                           | ``INTERFACE``                               |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.QMMM``                                    | ``params.inputSections.cpmd.qmmm``                                | ``qmmm``                                | ``QMMM``                                    |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.BICANONICAL_ENSEMBLE``                    | ``params.inputSections.cpmd.bicanonicalEnsemble``                 | ``bicanonicalEnsemble``                 | ``BICANONICAL ENSEMBLE``                    |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.CDFT``                                    | ``params.inputSections.cpmd.cdft``                                | ``cdft``                                | ``CDFT``                                    |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.CDFT_OPTIONS``                            | ``params.inputSections.cpmd.cdftOptions``                         | ``cdftOptions``                         | ``CDFT options``                            |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.CDFT_PAYLOAD``                            | ``params.inputSections.cpmd.cdftPayload``                         | ``cdftPayload``                         | ``CDFT payload``                            |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.CDFT_HDA_PAYLOAD``                        | ``params.inputSections.cpmd.cdftHdaPayload``                      | ``cdftHdaPayload``                      | ``CDFT HDA payload``                        |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.VGFACTOR``                                | ``params.inputSections.cpmd.vgfactor``                            | ``vgfactor``                            | ``VGFACTOR``                                |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.VMIRROR``                                 | ``params.inputSections.cpmd.vMirror``                             | ``vMirror``                             | ``VMIRROR``                                 |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.COMBINE_SYSTEMS``                         | ``params.inputSections.cpmd.combineSystemsOptions``               | ``combineSystemsOptions``               | ``COMBINE SYSTEMS``                         |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.COMBINE_SYSTEMS_PAYLOAD``                 | ``params.inputSections.cpmd.combineSystemsPayload``               | ``combineSystemsPayload``               | ``COMBINE SYSTEMS payload``                 |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.COMBINE_SYSTEMS_SAB_PAYLOAD``             | ``params.inputSections.cpmd.combineSystemsSabPayload``            | ``combineSystemsSabPayload``            | ``COMBINE SYSTEMS SAB payload``             |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.KSHAM``                                   | ``params.inputSections.cpmd.kshamOptions``                        | ``kshamOptions``                        | ``KSHAM``                                   |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.KSHAM_PAYLOAD``                           | ``params.inputSections.cpmd.kshamPayload``                        | ``kshamPayload``                        | ``KSHAM payload``                           |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.CZONES``                                  | ``params.inputSections.cpmd.czonesSet``                           | ``czonesSet``                           | ``CZONES SET``                              |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.WOUT``                                    | ``params.inputSections.cpmd.woutOptions``                         | ``woutOptions``                         | ``WOUT``                                    |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.WOUT_PAYLOAD``                            | ``params.inputSections.cpmd.woutPayload``                         | ``woutPayload``                         | ``WOUT payload``                            |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.XFMQC``                                   | ``params.inputSections.cpmd.xfmqcTrajectories``                   | ``xfmqcTrajectories``                   | ``XFMQC``                                   |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.DEBUG``                                   | ``params.inputSections.cpmd.debugOptions``                        | ``debugOptions``                        | ``DEBUG``                                   |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.DEBUG_CODE``                              | ``params.inputSections.cpmd.debugOptions``                        | ``debugOptions``                        | ``DEBUG CODE``                              |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.DEBUG_MEMORY``                            | ``params.inputSections.cpmd.debugOptions``                        | ``debugOptions``                        | ``DEBUG MEMORY``                            |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.DEBUG_FILEOPEN``                          | ``params.inputSections.cpmd.debugOptions``                        | ``debugOptions``                        | ``DEBUG FILEOPEN``                          |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.DEBUG_FORCES``                            | ``params.inputSections.cpmd.debugOptions``                        | ``debugOptions``                        | ``DEBUG FORCES``                            |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.DEBUG_NOACC``                             | ``params.inputSections.cpmd.debugOptions``                        | ``debugOptions``                        | ``DEBUG NOACC``                             |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.KOHN_SHAM_ENERGIES``                      | ``params.inputSections.cpmd.kohnShamEnergiesOptions``             | ``kohnShamEnergiesOptions``             | ``KOHN-SHAM ENERGIES``                      |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.KOHN_SHAM_ENERGIES_COUNT``                | ``params.inputSections.cpmd.kohnShamEnergiesCount``               | ``kohnShamEnergiesCount``               | ``KOHN-SHAM ENERGIES payload``              |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.SURFACE_HOPPING``                         | ``params.inputSections.cpmd.surfaceHoppingOptions``               | ``surfaceHoppingOptions``               | ``SURFACE HOPPING``                         |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.ROKS``                                    | ``params.inputSections.cpmd.roksOptions``                         | ``roksOptions``                         | ``ROKS``                                    |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.ROKS_EXPERT_PAYLOAD``                     | ``params.inputSections.cpmd.roksExpertPayload``                   | ``roksExpertPayload``                   | ``ROKS EXPERT payload``                     |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.PATH_SAMPLING``                           | ``params.inputSections.cpmd.pathSampling``                        | ``pathSampling``                        | ``PATH SAMPLING``                           |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.FIXRHO_UPWFN``                            | ``params.inputSections.cpmd.fixrhoUpwfnOptions``                  | ``fixrhoUpwfnOptions``                  | ``FIXRHO UPWFN``                            |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.FIXRHO_VECT``                             | ``params.inputSections.cpmd.fixrhoVectors``                       | ``fixrhoVectors``                       | ``FIXRHO VECT``                             |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.FIXRHO_LOOP``                             | ``params.inputSections.cpmd.fixrhoLoop``                          | ``fixrhoLoop``                          | ``FIXRHO LOOP``                             |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.FIXRHO_WFTOL``                            | ``params.inputSections.cpmd.fixrhoWftol``                         | ``fixrhoWftol``                         | ``FIXRHO WFTOL``                            |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.BOGOLIUBOV_CORRECTION``                   | ``params.inputSections.cpmd.bogoliubovCorrection``                | ``bogoliubovCorrection``                | ``BOGOLIUBOV CORRECTION``                   |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.VIBRATIONAL_ANALYSIS``                    | ``params.inputSections.cpmd.vibrationalAnalysisOptions``          | ``vibrationalAnalysisOptions``          | ``VIBRATIONAL ANALYSIS``                    |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.VIBRATIONAL_ANALYSIS_SAMPLE``             | ``params.inputSections.cpmd.vibrationalAnalysisSample``           | ``vibrationalAnalysisSample``           | ``VIBRATIONAL ANALYSIS SAMPLE``             |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.VIBRATIONAL_ANALYSIS_MODE``               | ``params.inputSections.cpmd.vibrationalAnalysisMode``             | ``vibrationalAnalysisMode``             | ``VIBRATIONAL ANALYSIS MODE=``              |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.ELECTRONIC_SPECTRA``                      | ``params.inputSections.cpmd.electronicSpectra``                   | ``electronicSpectra``                   | ``ELECTRONIC SPECTRA``                      |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.SPIN_ORBIT_COUPLING``                     | ``params.inputSections.cpmd.spinOrbitCouplingStates``             | ``spinOrbitCouplingStates``             | ``SPIN-ORBIT COUPLING``                     |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.SOC``                                     | ``params.inputSections.cpmd.spinOrbitCouplingStates``             | ``spinOrbitCouplingStates``             | ``SOC``                                     |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.PROPAGATION_SPECTRA``                     | ``params.inputSections.cpmd.propagationSpectra``                  | ``propagationSpectra``                  | ``PROPAGATION SPECTRA``                     |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.PROPAGATION_DISTRUB``                     | ``params.inputSections.cpmd.propagationDistrub``                  | ``propagationDistrub``                  | ``PROPAGATION DISTRUB``                     |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.GAUGEPULSE``                              | ``params.inputSections.cpmd.gaugePulse``                          | ``gaugePulse``                          | ``GAUGEPULSE``                              |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.GAUGEFIELD``                              | ``params.inputSections.cpmd.gaugeFieldFrequency``                 | ``gaugeFieldFrequency``                 | ``GAUGEFIELD``                              |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.NACV``                                    | ``params.inputSections.cpmd.nacv``                                | ``nacv``                                | ``NACV``                                    |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.ORBITAL_HARDNESS``                        | ``params.inputSections.cpmd.orbitalHardnessOptions``              | ``orbitalHardnessOptions``              | ``ORBITAL HARDNESS``                        |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.PATH_INTEGRAL``                           | ``params.inputSections.cpmd.pathIntegral``                        | ``pathIntegral``                        | ``PATH INTEGRAL``                           |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.PATH_MINIMIZATION``                       | ``params.inputSections.cpmd.pathMinimization``                    | ``pathMinimization``                    | ``PATH MINIMIZATION``                       |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.LANGEVIN_OPTIONS``                        | ``params.inputSections.cpmd.langevinOptions``                     | ``langevinOptions``                     | ``LANGEVIN``                                |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.LANGEVIN_PARAMETER``                      | ``params.inputSections.cpmd.langevinParameter``                   | ``langevinParameter``                   | ``LANGEVIN payload``                        |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.QMMMEASY``                                | ``params.inputSections.cpmd.qmmmEasy``                            | ``qmmmEasy``                            | ``QMMMEASY``                                |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.INTERFACE_OPTIONS``                       | ``params.inputSections.cpmd.interfaceOptions``                    | ``interfaceOptions``                    | ``INTERFACE``                               |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.TROTTER_FACTOR``                          | ``params.inputSections.cpmd.trotterFactorCount``                  | ``trotterFactorCount``                  | ``TROTTER FACTOR=``                         |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.TROTTER_FACTOR_PAYLOAD``                  | ``params.inputSections.cpmd.trotterFactorPayload``                | ``trotterFactorPayload``                | ``TROTTER FACTOR payload``                  |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.TROTTER_FACTORIZATION``                   | ``params.inputSections.cpmd.trotterFactorization``                | ``trotterFactorization``                | ``TROTTER FACTORIZATION``                   |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.TROTTER_FACTORIZATION_OFF``               | ``params.inputSections.cpmd.trotterFactorizationOff``             | ``trotterFactorizationOff``             | ``TROTTER FACTORIZATION OFF``               |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.LINEAR_RESPONSE``                         | ``params.inputSections.cpmd.linearResponse``                      | ``linearResponse``                      | ``LINEAR RESPONSE``                         |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.HARMONIC_REFERENCE``                      | ``params.inputSections.cpmd.harmonicReference``                   | ``harmonicReference``                   | ``HARMONIC REFERENCE``                      |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.HARMONIC_REFERENCE_SYSTEM``               | ``params.inputSections.cpmd.harmonicReference``                   | ``harmonicReference``                   | ``HARMONIC REFERENCE SYSTEM``               |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.SCALED_MASSES``                           | ``params.inputSections.cpmd.scaledMasses``                        | ``scaledMasses``                        | ``SCALED MASSES``                           |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.TDDFT``                                   | ``params.inputSections.cpmd.tddft``                               | ``tddft``                               | ``TDDFT``                                   |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.LSD``                                     | ``params.inputSections.cpmd.lsd``                                 | ``lsd``                                 | ``LSD``                                     |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.LOCAL_SPIN_DENSITY``                      | ``params.inputSections.cpmd.lsd``                                 | ``lsd``                                 | ``LOCAL SPIN DENSITY``                      |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.SSIC``                                    | ``params.inputSections.cpmd.ssic``                                | ``ssic``                                | ``SSIC``                                    |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.NONORTHOGONAL_ORBITALS``                  | ``params.inputSections.cpmd.nonorthogonalOrbitalsOptions``        | ``nonorthogonalOrbitalsOptions``        | ``NONORTHOGONAL ORBITALS``                  |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.NONORTHOGONAL_ORBITALS_LIMIT``            | ``params.inputSections.cpmd.nonorthogonalOrbitalsLimit``          | ``nonorthogonalOrbitalsLimit``          | ``NONORTHOGONAL ORBITALS payload``          |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.LANCZOS_DIAGONALIZATION``                 | ``params.inputSections.cpmd.lanczosDiagonalizationOptions``       | ``lanczosDiagonalizationOptions``       | ``LANCZOS DIAGONALIZATION``                 |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.LANCZOS_PARAMETERS``                      | ``params.inputSections.cpmd.lanczosParametersCount``              | ``lanczosParametersCount``              | ``LANCZOS PARAMETERS N=``                   |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.LANCZOS_PARAMETER``                       | ``params.inputSections.cpmd.lanczosParametersCount``              | ``lanczosParametersCount``              | ``LANCZOS PARAMETER``                       |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.LANCZOS_PARAMETERS_PAYLOAD``              | ``params.inputSections.cpmd.lanczosParametersPayload``            | ``lanczosParametersPayload``            | ``LANCZOS PARAMETERS payload``              |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.DAVIDSON_DIAGONALIZATION``                | ``params.inputSections.cpmd.davidsonDiagonalization``             | ``davidsonDiagonalization``             | ``DAVIDSON DIAGONALIZATION``                |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.DAVIDSON_PARAMETERS``                     | ``params.inputSections.cpmd.davidsonParameters``                  | ``davidsonParameters``                  | ``DAVIDSON PARAMETERS``                     |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.DAVIDSON_PARAMETER``                      | ``params.inputSections.cpmd.davidsonParameters``                  | ``davidsonParameters``                  | ``DAVIDSON PARAMETER``                      |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.ALEXANDER_MIXING``                        | ``params.inputSections.cpmd.alexanderMixing``                     | ``alexanderMixing``                     | ``ALEXANDER MIXING``                        |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.ANDERSON_MIXING``                         | ``params.inputSections.cpmd.andersonMixingGspace``                | ``andersonMixingGspace``                | ``ANDERSON MIXING G-SPACE``                 |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.ANDERSON_MIXING``                         | ``params.inputSections.cpmd.andersonMixingCount``                 | ``andersonMixingCount``                 | ``ANDERSON MIXING N=``                      |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.ANDERSON_MIXING_PAYLOAD``                 | ``params.inputSections.cpmd.andersonMixingPayload``               | ``andersonMixingPayload``               | ``ANDERSON MIXING payload``                 |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.BROYDEN_MIXING``                          | ``params.inputSections.cpmd.broydenMixingOptions``                | ``broydenMixingOptions``                | ``BROYDEN MIXING``                          |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.BROYDEN_MIXING_PAYLOAD``                  | ``params.inputSections.cpmd.broydenMixingPayload``                | ``broydenMixingPayload``                | ``BROYDEN MIXING payload``                  |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.DIIS_MIXING``                             | ``params.inputSections.cpmd.diisMixingCount``                     | ``diisMixingCount``                     | ``DIIS MIXING N=``                          |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.DIIS_MIXING_PAYLOAD``                     | ``params.inputSections.cpmd.diisMixingPayload``                   | ``diisMixingPayload``                   | ``DIIS MIXING payload``                     |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.MOVERHO``                                 | ``params.inputSections.cpmd.moverhoMixing``                       | ``moverhoMixing``                       | ``MOVERHO``                                 |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.EXTRAPOLATE_WFN``                         | ``params.inputSections.cpmd.extrapolateWfnOptions``               | ``extrapolateWfnOptions``               | ``EXTRAPOLATE WFN``                         |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.EXTRAPOLATE_WFN_ORDER``                   | ``params.inputSections.cpmd.extrapolateWfnOrder``                 | ``extrapolateWfnOrder``                 | ``EXTRAPOLATE WFN payload``                 |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.EXTRAPOLATE_CONSTRAINT``                  | ``params.inputSections.cpmd.extrapolateConstraintOrder``          | ``extrapolateConstraintOrder``          | ``EXTRAPOLATE CONSTRAINT``                  |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.TSDE``                                    | ``params.inputSections.cpmd.tsdeOptions``                         | ``tsdeOptions``                         | ``TSDE``                                    |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.TSDP``                                    | ``params.inputSections.cpmd.tsdpOptions``                         | ``tsdpOptions``                         | ``TSDP``                                    |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.TCGP``                                    | ``params.inputSections.cpmd.tcgp``                                | ``tcgp``                                | ``TCGP``                                    |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.TSDC``                                    | ``params.inputSections.cpmd.tsdc``                                | ``tsdc``                                | ``TSDC``                                    |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.STEEPEST_DESCENT``                        | ``params.inputSections.cpmd.steepestDescentOptions``              | ``steepestDescentOptions``              | ``STEEPEST DESCENT``                        |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.CONJUGATE_GRADIENT``                      | ``params.inputSections.cpmd.conjugateGradientOptions``            | ``conjugateGradientOptions``            | ``CONJUGATE GRADIENT``                      |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.CONJUGATE_GRADIENTS``                     | ``params.inputSections.cpmd.conjugateGradientOptions``            | ``conjugateGradientOptions``            | ``CONJUGATE GRADIENTS``                     |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.ODIIS_OPTIONS``                           | ``params.inputSections.cpmd.odiisOptions``                        | ``odiisOptions``                        | ``ODIIS``                                   |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.ODIIS_VECTORS``                           | ``params.inputSections.cpmd.odiisVectors``                        | ``odiisVectors``                        | ``ODIIS payload``                           |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.HAMILTONIAN_CUTOFF``                      | ``params.inputSections.cpmd.hamiltonianCutoff``                   | ``hamiltonianCutoff``                   | ``HAMILTONIAN CUTOFF``                      |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.GDIIS``                                   | ``params.inputSections.cpmd.gdiisVectors``                        | ``gdiisVectors``                        | ``GDIIS``                                   |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.LBFGS``                                   | ``params.inputSections.cpmd.lbfgsOptions``                        | ``lbfgsOptions``                        | ``LBFGS``                                   |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.LBFGS_PAYLOAD``                           | ``params.inputSections.cpmd.lbfgsPayload``                        | ``lbfgsPayload``                        | ``LBFGS payload``                           |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.LBFGS_NTRUST``                            | ``params.inputSections.cpmd.lbfgsNtrust``                         | ``lbfgsNtrust``                         | ``LBFGS NTRUST``                            |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.LBFGS_NRESTT``                            | ``params.inputSections.cpmd.lbfgsNrestt``                         | ``lbfgsNrestt``                         | ``LBFGS NRESTT``                            |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.LBFGS_NTRSTR``                            | ``params.inputSections.cpmd.lbfgsNtrstr``                         | ``lbfgsNtrstr``                         | ``LBFGS NTRSTR``                            |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.LBFGS_TRUSTR``                            | ``params.inputSections.cpmd.lbfgsTrustr``                         | ``lbfgsTrustr``                         | ``LBFGS TRUSTR``                            |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.PRFO``                                    | ``params.inputSections.cpmd.prfoOptions``                         | ``prfoOptions``                         | ``PRFO``                                    |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.PRFO_PAYLOAD``                            | ``params.inputSections.cpmd.prfoPayload``                         | ``prfoPayload``                         | ``PRFO payload``                            |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.PRFO_MODE``                               | ``params.inputSections.cpmd.prfoMode``                            | ``prfoMode``                            | ``PRFO MODE``                               |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.PRFO_MDLOCK``                             | ``params.inputSections.cpmd.prfoMdlock``                          | ``prfoMdlock``                          | ``PRFO MDLOCK``                             |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.PRFO_TOLENV``                             | ``params.inputSections.cpmd.prfoTolenv``                          | ``prfoTolenv``                          | ``PRFO TOLENV``                             |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.PRFO_TRUSTP``                             | ``params.inputSections.cpmd.prfoTrustp``                          | ``prfoTrustp``                          | ``PRFO TRUSTP``                             |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.PRFO_OMIN``                               | ``params.inputSections.cpmd.prfoOmin``                            | ``prfoOmin``                            | ``PRFO OMIN``                               |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.PRFO_NSVIB``                              | ``params.inputSections.cpmd.prfoNsvib``                           | ``prfoNsvib``                           | ``PRFO NSVIB``                              |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.PRFO_CORE``                               | ``params.inputSections.cpmd.prfoCoreAtoms``                       | ``prfoCoreAtoms``                       | ``PRFO CORE=``                              |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.PRFO_NSMAXP``                             | ``params.inputSections.cpmd.prfoNsmaxp``                          | ``prfoNsmaxp``                          | ``PRFO NSMAXP``                             |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.PRFO_PRJHES``                             | ``params.inputSections.cpmd.prfoProjectedHessian``                | ``prfoProjectedHessian``                | ``PRFO PRJHES``                             |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.PRFO_DISPL``                              | ``params.inputSections.cpmd.prfoDisplacement``                    | ``prfoDisplacement``                    | ``PRFO DISPL``                              |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.PRFO_HESSTYPE``                           | ``params.inputSections.cpmd.prfoHessianType``                     | ``prfoHessianType``                     | ``PRFO HESSTYPE``                           |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.HESSCORE``                                | ``params.inputSections.cpmd.hesscore``                            | ``hesscore``                            | ``HESSCORE``                                |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.BFGS``                                    | ``params.inputSections.cpmd.bfgs``                                | ``bfgs``                                | ``BFGS``                                    |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.RFO``                                     | ``params.inputSections.cpmd.rfoOrder``                            | ``rfoOrder``                            | ``RFO ORDER=``                              |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.RFO_ORDER``                               | ``params.inputSections.cpmd.rfoOrder``                            | ``rfoOrder``                            | ``RFO ORDER=nsorder``                       |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.INR_PARAMETER``                           | ``params.inputSections.cpmd.inrParametersCount``                  | ``inrParametersCount``                  | ``INR PARAMETER N=``                        |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.INR_PARAMETERS``                          | ``params.inputSections.cpmd.inrParametersCount``                  | ``inrParametersCount``                  | ``INR PARAMETERS N=``                       |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.TURBOCICCIO_PARAMETER``                   | ``params.inputSections.cpmd.inrParametersCount``                  | ``inrParametersCount``                  | ``TURBOCICCIO PARAMETER N=``                |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.TURBOCICCIO_PARAMETERS``                  | ``params.inputSections.cpmd.inrParametersCount``                  | ``inrParametersCount``                  | ``TURBOCICCIO PARAMETERS N=``               |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.INR_PARAMETERS_PAYLOAD``                  | ``params.inputSections.cpmd.inrParametersPayload``                | ``inrParametersPayload``                | ``INR PARAMETERS payload``                  |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.TURBOCICCIO``                             | ``params.inputSections.cpmd.implicitNewtonOptions``               | ``implicitNewtonOptions``               | ``TURBOCICCIO``                             |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.IMPLICIT_NEWTON``                         | ``params.inputSections.cpmd.implicitNewtonOptions``               | ``implicitNewtonOptions``               | ``IMPLICIT NEWTON``                         |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.IMPLICIT_NEWTON_MAXITER``                 | ``params.inputSections.cpmd.implicitNewtonMaxIter``               | ``implicitNewtonMaxIter``               | ``IMPLICIT NEWTON payload``                 |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.MIXSD``                                   | ``params.inputSections.cpmd.mixsd``                               | ``mixsd``                               | ``MIXSD``                                   |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.MIXDIIS``                                 | ``params.inputSections.cpmd.mixdiis``                             | ``mixdiis``                             | ``MIXDIIS``                                 |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.RESTART_OPTIONS``                         | ``params.inputSections.cpmd.restartOptions``                      | ``restartOptions``                      | ``RESTART``                                 |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.INTFILE``                                 | ``params.inputSections.cpmd.intFileOptions``                      | ``intFileOptions``                      | ``INTFILE``                                 |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.INTFILE_NAME``                            | ``params.inputSections.cpmd.intFileName``                         | ``intFileName``                         | ``INTFILE filename``                        |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.INITIALIZE_WAVEFUNCTION``                 | ``params.inputSections.cpmd.initializeWavefunctionOptions``       | ``initializeWavefunctionOptions``       | ``INITIALIZE WAVEFUNCTION``                 |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.RATTLE_PARAMETERS``                       | ``params.inputSections.cpmd.rattleParameters``                    | ``rattleParameters``                    | ``RATTLE payload``                          |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.ORTHOGONALIZATION``                       | ``params.inputSections.cpmd.orthogonalizationOptions``            | ``orthogonalizationOptions``            | ``ORTHOGONALIZATION``                       |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.QUENCH_OPTIONS``                          | ``params.inputSections.cpmd.quenchOptions``                       | ``quenchOptions``                       | ``QUENCH``                                  |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.RANDOMIZE``                               | ``params.inputSections.cpmd.randomizeOptions``                    | ``randomizeOptions``                    | ``RANDOMIZE``                               |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.RANDOMIZE_AMPLITUDE``                     | ``params.inputSections.cpmd.randomizeAmplitude``                  | ``randomizeAmplitude``                  | ``RANDOMIZE payload``                       |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.USE_MTS``                                 | ``params.inputSections.cpmd.useMts``                              | ``useMts``                              | ``USE_MTS``                                 |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.NABDY_ZMAX``                              | ``params.inputSections.cpmd.nabdyZmax``                           | ``nabdyZmax``                           | ``NABDY_ZMAX``                              |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.NABDY_SOFT``                              | ``params.inputSections.cpmd.nabdySoft``                           | ``nabdySoft``                           | ``NABDY_SOFT``                              |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.NABDY_REDISTR_AMPLI``                     | ``params.inputSections.cpmd.nabdyRedistributeAmplitude``          | ``nabdyRedistributeAmplitude``          | ``NABDY_REDISTR_AMPLI``                     |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.NABDY_SCALEP``                            | ``params.inputSections.cpmd.nabdyScaleP``                         | ``nabdyScaleP``                         | ``NABDY_SCALEP``                            |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.NABDY_THERMO``                            | ``params.inputSections.cpmd.nabdyThermo``                         | ``nabdyThermo``                         | ``NABDY_THERMO``                            |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.PROPERTIES``                              | ``params.inputSections.cpmd.properties``                          | ``properties``                          | ``PROPERTIES``                              |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.VDW_CORRECTION``                          | ``params.inputSections.cpmd.vdwCorrection``                       | ``vdwCorrection``                       | ``VDW CORRECTION``                          |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.VDW_WANNIER``                             | ``params.inputSections.cpmd.vdwWannier``                          | ``vdwWannier``                          | ``VDW WANNIER``                             |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+
| ``catalog.cpmd.DCACP``                                   | ``params.inputSections.cpmd.dcacp``                               | ``dcacp``                               | ``DCACP``                                   |
+----------------------------------------------------------+-------------------------------------------------------------------+-----------------------------------------+---------------------------------------------+

Typed ``&DFT`` controls
=======================

``CPMDDftSection`` exposes common DFT scalars directly. Parser-accepted
escape-hatch rows use ``directives`` so callers can still discover the
catalog ID without a dedicated scalar field.

+------------------------------------------+-------------------------------------------+------------------+------------------------------+
| Catalog feature ID                       | Parameter feature ID                      | Field            | Deck keyword                 |
+==========================================+===========================================+==================+==============================+
| ``catalog.dft.FUNCTIONAL``               | ``params.inputSections.dft.functional``   | ``functional``   | ``FUNCTIONAL``               |
+------------------------------------------+-------------------------------------------+------------------+------------------------------+
| ``catalog.dft.USE_CP``                   | ``params.inputSections.dft.directives``   | ``directives``   | ``USE_CP``                   |
+------------------------------------------+-------------------------------------------+------------------+------------------------------+
| ``catalog.dft.USE_LIBXC``                | ``params.inputSections.dft.directives``   | ``directives``   | ``USE_LIBXC``                |
+------------------------------------------+-------------------------------------------+------------------+------------------------------+
| ``catalog.dft.LSD``                      | ``params.inputSections.dft.lsd``          | ``lsd``          | ``LSD``                      |
+------------------------------------------+-------------------------------------------+------------------+------------------------------+
| ``catalog.dft.GC_CUTOFF``                | ``params.inputSections.dft.gcCutoff``     | ``gcCutoff``     | ``GC-CUTOFF``                |
+------------------------------------------+-------------------------------------------+------------------+------------------------------+
| ``catalog.dft.LIBRARY``                  | ``params.inputSections.dft.directives``   | ``directives``   | ``LIBRARY``                  |
+------------------------------------------+-------------------------------------------+------------------+------------------------------+
| ``catalog.dft.KERNEL_LIBRARY``           | ``params.inputSections.dft.directives``   | ``directives``   | ``KERNEL_LIBRARY``           |
+------------------------------------------+-------------------------------------------+------------------+------------------------------+
| ``catalog.dft.MTS_LOW_FUNC_LIBRARY``     | ``params.inputSections.dft.directives``   | ``directives``   | ``MTS_LOW_FUNC_LIBRARY``     |
+------------------------------------------+-------------------------------------------+------------------+------------------------------+
| ``catalog.dft.SCALES``                   | ``params.inputSections.dft.directives``   | ``directives``   | ``SCALES``                   |
+------------------------------------------+-------------------------------------------+------------------+------------------------------+
| ``catalog.dft.KERNEL_SCALES``            | ``params.inputSections.dft.directives``   | ``directives``   | ``KERNEL_SCALES``            |
+------------------------------------------+-------------------------------------------+------------------+------------------------------+
| ``catalog.dft.MTS_LOW_FUNC_SCALES``      | ``params.inputSections.dft.directives``   | ``directives``   | ``MTS_LOW_FUNC_SCALES``      |
+------------------------------------------+-------------------------------------------+------------------+------------------------------+
| ``catalog.dft.COMPATIBILITY``            | ``params.inputSections.dft.directives``   | ``directives``   | ``COMPATIBILITY``            |
+------------------------------------------+-------------------------------------------+------------------+------------------------------+
| ``catalog.dft.ANALYTICAL_DIV``           | ``params.inputSections.dft.directives``   | ``directives``   | ``ANALYTICAL_DIV``           |
+------------------------------------------+-------------------------------------------+------------------+------------------------------+
| ``catalog.dft.NUMERICAL_DIV``            | ``params.inputSections.dft.directives``   | ``directives``   | ``NUMERICAL_DIV``            |
+------------------------------------------+-------------------------------------------+------------------+------------------------------+
| ``catalog.dft.PBE_FLEX_KAPPA``           | ``params.inputSections.dft.directives``   | ``directives``   | ``PBE_FLEX_KAPPA``           |
+------------------------------------------+-------------------------------------------+------------------+------------------------------+
| ``catalog.dft.PBE_FLEX_MU``              | ``params.inputSections.dft.directives``   | ``directives``   | ``PBE_FLEX_MU``              |
+------------------------------------------+-------------------------------------------+------------------+------------------------------+
| ``catalog.dft.PBE_FLEX_BETA``            | ``params.inputSections.dft.directives``   | ``directives``   | ``PBE_FLEX_BETA``            |
+------------------------------------------+-------------------------------------------+------------------+------------------------------+
| ``catalog.dft.PBE_FLEX_GAMMA``           | ``params.inputSections.dft.directives``   | ``directives``   | ``PBE_FLEX_GAMMA``           |
+------------------------------------------+-------------------------------------------+------------------+------------------------------+
| ``catalog.dft.PBE_FLEX_UEG_CORRELATION`` | ``params.inputSections.dft.directives``   | ``directives``   | ``PBE_FLEX_UEG_CORRELATION`` |
+------------------------------------------+-------------------------------------------+------------------+------------------------------+
| ``catalog.dft.GRADIENT_CORRECTION``      | ``params.inputSections.dft.directives``   | ``directives``   | ``GRADIENT CORRECTION``      |
+------------------------------------------+-------------------------------------------+------------------+------------------------------+
| ``catalog.dft.HARTREE_FOCK``             | ``params.inputSections.dft.directives``   | ``directives``   | ``HARTREE-FOCK``             |
+------------------------------------------+-------------------------------------------+------------------+------------------------------+
| ``catalog.dft.HARTREE``                  | ``params.inputSections.dft.directives``   | ``directives``   | ``HARTREE``                  |
+------------------------------------------+-------------------------------------------+------------------+------------------------------+
| ``catalog.dft.SLATER``                   | ``params.inputSections.dft.directives``   | ``directives``   | ``SLATER``                   |
+------------------------------------------+-------------------------------------------+------------------+------------------------------+
| ``catalog.dft.LDA_CORRELATION``          | ``params.inputSections.dft.directives``   | ``directives``   | ``LDA CORRELATION``          |
+------------------------------------------+-------------------------------------------+------------------+------------------------------+
| ``catalog.dft.OLDCODE``                  | ``params.inputSections.dft.oldCode``      | ``oldCode``      | ``OLDCODE``                  |
+------------------------------------------+-------------------------------------------+------------------+------------------------------+
| ``catalog.dft.NEWCODE``                  | ``params.inputSections.dft.newCode``      | ``newCode``      | ``NEWCODE``                  |
+------------------------------------------+-------------------------------------------+------------------+------------------------------+
| ``catalog.dft.XC_DRIVER``                | ``params.inputSections.dft.xcDriver``     | ``xcDriver``     | ``XC_DRIVER``                |
+------------------------------------------+-------------------------------------------+------------------+------------------------------+
| ``catalog.dft.LIBXC``                    | ``params.inputSections.dft.libxc``        | ``libxc``        | ``LIBXC``                    |
+------------------------------------------+-------------------------------------------+------------------+------------------------------+
| ``catalog.dft.LR_KERNEL``                | ``params.inputSections.dft.lrKernel``     | ``lrKernel``     | ``LR KERNEL``                |
+------------------------------------------+-------------------------------------------+------------------+------------------------------+
| ``catalog.dft.REFUNCT``                  | ``params.inputSections.dft.refunct``      | ``refunct``      | ``REFUNCT``                  |
+------------------------------------------+-------------------------------------------+------------------+------------------------------+
| ``catalog.dft.MTS_HIGH_FUNC``            | ``params.inputSections.dft.mtsHighFunc``  | ``mtsHighFunc``  | ``MTS_HIGH_FUNC``            |
+------------------------------------------+-------------------------------------------+------------------+------------------------------+
| ``catalog.dft.MTS_LOW_FUNC``             | ``params.inputSections.dft.mtsLowFunc``   | ``mtsLowFunc``   | ``MTS_LOW_FUNC``             |
+------------------------------------------+-------------------------------------------+------------------+------------------------------+
| ``catalog.dft.HFX``                      | ``params.inputSections.dft.hfx``          | ``hfx``          | ``HFX``                      |
+------------------------------------------+-------------------------------------------+------------------+------------------------------+
| ``catalog.dft.HFX_BLOCK_SIZE``           | ``params.inputSections.dft.directives``   | ``directives``   | ``HFX_BLOCK_SIZE``           |
+------------------------------------------+-------------------------------------------+------------------+------------------------------+
| ``catalog.dft.SCALED_EXCHANGE``          | ``params.inputSections.dft.directives``   | ``directives``   | ``SCALED EXCHANGE``          |
+------------------------------------------+-------------------------------------------+------------------+------------------------------+
| ``catalog.dft.HFX_DISTRIBUTION``         | ``params.inputSections.dft.directives``   | ``directives``   | ``HFX_DISTRIBUTION``         |
+------------------------------------------+-------------------------------------------+------------------+------------------------------+
| ``catalog.dft.RANGE_SEPARATION``         | ``params.inputSections.dft.directives``   | ``directives``   | ``RANGE SEPARATION``         |
+------------------------------------------+-------------------------------------------+------------------+------------------------------+
| ``catalog.dft.SCREENED_EXCHANGE``        | ``params.inputSections.dft.directives``   | ``directives``   | ``SCREENED EXCHANGE``        |
+------------------------------------------+-------------------------------------------+------------------+------------------------------+
| ``catalog.dft.BECKE_BETA``               | ``params.inputSections.dft.directives``   | ``directives``   | ``BECKE BETA``               |
+------------------------------------------+-------------------------------------------+------------------+------------------------------+
| ``catalog.dft.EXCHANGE_TABLE``           | ``params.inputSections.dft.directives``   | ``directives``   | ``EXCHANGE TABLE``           |
+------------------------------------------+-------------------------------------------+------------------+------------------------------+
| ``catalog.dft.HFX_SCALE``                | ``params.inputSections.dft.directives``   | ``directives``   | ``HFX_SCALE``                |
+------------------------------------------+-------------------------------------------+------------------+------------------------------+
| ``catalog.dft.KERNEL_HFX_SCALE``         | ``params.inputSections.dft.directives``   | ``directives``   | ``KERNEL_HFX_SCALE``         |
+------------------------------------------+-------------------------------------------+------------------+------------------------------+
| ``catalog.dft.MTS_LOW_FUNC_HFX_SCALE``   | ``params.inputSections.dft.directives``   | ``directives``   | ``MTS_LOW_FUNC_HFX_SCALE``   |
+------------------------------------------+-------------------------------------------+------------------+------------------------------+
| ``catalog.dft.PHFX``                     | ``params.inputSections.dft.directives``   | ``directives``   | ``PHFX``                     |
+------------------------------------------+-------------------------------------------+------------------+------------------------------+
| ``catalog.dft.PXGC``                     | ``params.inputSections.dft.directives``   | ``directives``   | ``PXGC``                     |
+------------------------------------------+-------------------------------------------+------------------+------------------------------+
| ``catalog.dft.PCGC``                     | ``params.inputSections.dft.directives``   | ``directives``   | ``PCGC``                     |
+------------------------------------------+-------------------------------------------+------------------+------------------------------+
| ``catalog.dft.PXLDA``                    | ``params.inputSections.dft.directives``   | ``directives``   | ``PXLDA``                    |
+------------------------------------------+-------------------------------------------+------------------+------------------------------+
| ``catalog.dft.PCLDA``                    | ``params.inputSections.dft.directives``   | ``directives``   | ``PCLDA``                    |
+------------------------------------------+-------------------------------------------+------------------+------------------------------+
| ``catalog.dft.HFX_SCREENING``            | ``params.inputSections.dft.hfxScreening`` | ``hfxScreening`` | ``HFX-SCREENING``            |
+------------------------------------------+-------------------------------------------+------------------+------------------------------+
| ``catalog.dft.ACM0``                     | ``params.inputSections.dft.directives``   | ``directives``   | ``ACM0``                     |
+------------------------------------------+-------------------------------------------+------------------+------------------------------+
| ``catalog.dft.ACM1``                     | ``params.inputSections.dft.directives``   | ``directives``   | ``ACM1``                     |
+------------------------------------------+-------------------------------------------+------------------+------------------------------+
| ``catalog.dft.ACM3``                     | ``params.inputSections.dft.directives``   | ``directives``   | ``ACM3``                     |
+------------------------------------------+-------------------------------------------+------------------+------------------------------+
| ``catalog.dft.FORCE_NEW_HFX``            | ``params.inputSections.dft.directives``   | ``directives``   | ``FORCE_NEW_HFX``            |
+------------------------------------------+-------------------------------------------+------------------+------------------------------+
| ``catalog.dft.SMOOTHING``                | ``params.inputSections.dft.directives``   | ``directives``   | ``SMOOTHING``                |
+------------------------------------------+-------------------------------------------+------------------+------------------------------+
| ``catalog.dft.HUBBARD``                  | ``params.inputSections.dft.hubbard``      | ``hubbard``      | ``HUBBARD``                  |
+------------------------------------------+-------------------------------------------+------------------+------------------------------+
| ``catalog.dft.CORRELATION``              | ``params.inputSections.dft.correlation``  | ``correlation``  | ``CORRELATION``              |
+------------------------------------------+-------------------------------------------+------------------+------------------------------+
| ``catalog.dft.EXCHANGE``                 | ``params.inputSections.dft.exchange``     | ``exchange``     | ``EXCHANGE``                 |
+------------------------------------------+-------------------------------------------+------------------+------------------------------+
| ``catalog.dft.BECKE88``                  | ``params.inputSections.dft.becke88``      | ``becke88``      | ``BECKE88``                  |
+------------------------------------------+-------------------------------------------+------------------+------------------------------+
| ``catalog.dft.ALPHA``                    | ``params.inputSections.dft.alpha``        | ``alpha``        | ``ALPHA``                    |
+------------------------------------------+-------------------------------------------+------------------+------------------------------+
| ``catalog.dft.BETA``                     | ``params.inputSections.dft.beta``         | ``beta``         | ``BETA``                     |
+------------------------------------------+-------------------------------------------+------------------+------------------------------+

``catalog.dft.FUNCTIONAL_*`` IDs name accepted ``functional`` values:

======================================= ====================
Feature ID                              ``functional`` value
======================================= ====================
``catalog.dft.FUNCTIONAL_LDA``          ``LDA``
``catalog.dft.FUNCTIONAL_NONE``         ``NONE``
``catalog.dft.FUNCTIONAL_SONLY``        ``SONLY``
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
``catalog.dft.FUNCTIONAL_ZPBE0``        ``ZPBE0``
``catalog.dft.FUNCTIONAL_XPBE0``        ``XPBE0``
``catalog.dft.FUNCTIONAL_PBES``         ``PBES``
``catalog.dft.FUNCTIONAL_PBE1W``        ``PBE1W``
``catalog.dft.FUNCTIONAL_TPSS``         ``TPSS``
``catalog.dft.FUNCTIONAL_HARTREE``      ``HARTREE``
``catalog.dft.FUNCTIONAL_HARTREE_FOCK`` ``HARTREE-FOCK``
``catalog.dft.FUNCTIONAL_HSE06``        ``HSE06``
``catalog.dft.FUNCTIONAL_ACM0``         ``ACM0``
``catalog.dft.FUNCTIONAL_ACM1``         ``ACM1``
``catalog.dft.FUNCTIONAL_ACM3``         ``ACM3``
``catalog.dft.FUNCTIONAL_OLDX3LYP``     ``OLDX3LYP``
``catalog.dft.FUNCTIONAL_MODX3LYP``     ``MODX3LYP``
``catalog.dft.FUNCTIONAL_X3LYP``        ``X3LYP``
``catalog.dft.FUNCTIONAL_OLDB3LYP``     ``OLDB3LYP``
``catalog.dft.FUNCTIONAL_MODB3LYP``     ``MODB3LYP``
``catalog.dft.FUNCTIONAL_OLDB1LYP``     ``OLDB1LYP``
``catalog.dft.FUNCTIONAL_MODB1LYP``     ``MODB1LYP``
``catalog.dft.FUNCTIONAL_SAOP``         ``SAOP``
``catalog.dft.FUNCTIONAL_LB94``         ``LB94``
``catalog.dft.FUNCTIONAL_GLLB``         ``GLLB``
======================================= ====================

``PotentialConfig``
===================

For rgpot / multi-backend configure RPCs, ``PotentialConfig`` is a
tagged union with ``cpmd @2 :CPMDParams`` (ordinal aligned with the
shared schema; ``nwchem @1`` is reserved void in this package).
