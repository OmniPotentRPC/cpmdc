#include "cpmdc_features.h"
#include <setjmp.h>
#include <stdarg.h>
#include <stddef.h>
#include <cmocka.h>

static void test_feature_table_nonempty(void **state) {
  (void)state;
  assert_true(cpmdc_feature_count() > 50);
}

static void test_inscan_sections(void **state) {
  (void)state;
  static const char *const section_kinds[] = {
      "section.generic", "section.atom",      "section.atoms",
      "section.basis",   "section.clas",      "section.cpmd",
      "section.dft",     "section.eam",       "section.exte",
      "section.hardness", "section.info",     "section.linres",
      "section.molstates", "section.mts",     "section.nlcc",
      "section.path",    "section.pimd",      "section.potential",
      "section.prop",    "section.ptddft",    "section.resp",
      "section.system",  "section.tddft",     "section.vdw",
      "section.vectors", "section.wavefunction", "section.set",
      "section.raw",
  };

  for (size_t i = 0; i < sizeof(section_kinds) / sizeof(section_kinds[0]);
       i++) {
    const CPMDCFeatureEntry *entry = cpmdc_feature_find(section_kinds[i]);
    assert_non_null(entry);
    assert_string_equal(entry->feature_id, section_kinds[i]);
    assert_int_equal(entry->kind, CPMDC_FEATURE_SECTION);
    assert_int_equal(entry->stub_applicable, 1);
    assert_int_equal(entry->embed_applicable, 1);
  }

  assert_non_null(cpmdc_feature_find("catalog.section.EAM"));
  assert_non_null(cpmdc_feature_find("catalog.section.MOLSTATES"));
  assert_non_null(cpmdc_feature_find("catalog.section.NLCC"));
  assert_non_null(cpmdc_feature_find("catalog.section.VECTORS"));
  assert_non_null(cpmdc_feature_find("catalog.section.PIMD"));
  assert_non_null(cpmdc_feature_find("catalog.section.PROP"));
  /* keyword-only must not be catalog.section.* */
  assert_null(cpmdc_feature_find("catalog.section.QMMM"));
  assert_null(cpmdc_feature_find("catalog.section.PROPERTIES"));
  assert_null(cpmdc_feature_find("catalog.section.BICANONICAL"));
  assert_null(cpmdc_feature_find("catalog.section.CDFT"));
}

static void test_cp_keywords_not_sections(void **state) {
  (void)state;
  assert_non_null(cpmdc_feature_find("catalog.cpmd.QMMM"));
  assert_non_null(cpmdc_feature_find("catalog.cpmd.EMASS"));
  assert_non_null(cpmdc_feature_find("catalog.cpmd.VDW_CORRECTION"));
  assert_non_null(cpmdc_feature_find("catalog.cpmd.VDW_WANNIER"));
  assert_non_null(cpmdc_feature_find("catalog.cpmd.DCACP"));
  assert_non_null(cpmdc_feature_find("catalog.cpmd.ISOLATED_MOLECULE"));
  assert_non_null(cpmdc_feature_find("catalog.dft.FUNCTIONAL_BLYP"));
}

static void test_abi(void **state) {
  (void)state;
  static const char *const abi_features[] = {
      "abi.cpmdc_set_params",
      "abi.cpmdc_energy_gradient",
      "abi.cpmdc_energy",
      "abi.cpmdc_energy_forces",
      "abi.cpmdc_session_create",
      "abi.cpmdc_session_set_params",
      "abi.cpmdc_session_destroy",
      "abi.cpmdc_session_energy_gradient",
      "abi.cpmdc_session_energy",
      "abi.cpmdc_session_energy_forces",
      "abi.cpmdc_session_calculate_forces",
      "abi.cpmdc_session_calculate_result",
      "abi.cpmdc_calculate_result",
      "abi.cpmdc_potential_result_size_for_force_input",
      "abi.cpmdc_version",
      "abi.cpmdc_available",
      "abi.cpmdc_finalize",
      "abi.cpmdc_feature_count",
      "abi.cpmdc_feature_table",
      "abi.cpmdc_feature_find",
  };

  const CPMDCFeatureEntry *table = cpmdc_feature_table();
  assert_non_null(table);
  assert_true(cpmdc_feature_count() >=
              sizeof(abi_features) / sizeof(abi_features[0]));

  for (size_t i = 0; i < sizeof(abi_features) / sizeof(abi_features[0]); i++) {
    const CPMDCFeatureEntry *entry = cpmdc_feature_find(abi_features[i]);
    assert_non_null(entry);
    assert_string_equal(entry->feature_id, abi_features[i]);
    assert_int_equal(entry->kind, CPMDC_FEATURE_ABI);
    assert_int_equal(entry->stub_applicable, 1);
    assert_int_equal(entry->embed_applicable, 1);
  }

  assert_null(cpmdc_feature_find("abi.cpmdc_missing"));
  assert_non_null(cpmdc_feature_find("params.inputBlocks"));
}

static void test_param_applicability(void **state) {
  (void)state;
  const CPMDCFeatureEntry *memory = cpmdc_feature_find("params.memoryMb");
  const CPMDCFeatureEntry *engine = cpmdc_feature_find("params.enginePath");
  const CPMDCFeatureEntry *scratch = cpmdc_feature_find("params.scratchDir");
  const CPMDCFeatureEntry *permanent = cpmdc_feature_find("params.permanentDir");
  assert_non_null(memory);
  assert_non_null(engine);
  assert_non_null(scratch);
  assert_non_null(permanent);
  assert_int_equal(memory->embed_applicable, 0);
  assert_int_equal(engine->embed_applicable, 0);
  assert_int_equal(scratch->embed_applicable, 1);
  assert_int_equal(permanent->embed_applicable, 1);
}

static void test_structured_param_features(void **state) {
  (void)state;
  static const char *const param_features[] = {
      "params.inputSections.generic.name",
      "params.inputSections.generic.directives",
      "params.inputSections.atom.directives",
      "params.inputSections.atom.subsections",
      "params.inputSections.system.symmetry",
      "params.inputSections.system.angstrom",
      "params.inputSections.system.cell",
      "params.inputSections.system.cellAbsolute",
      "params.inputSections.system.cellDegree",
      "params.inputSections.system.cellVectors",
      "params.inputSections.system.referenceCell",
      "params.inputSections.system.referenceCellAbsolute",
      "params.inputSections.system.referenceCellDegree",
      "params.inputSections.system.referenceCellVectors",
      "params.inputSections.system.classicalCell",
      "params.inputSections.system.classicalCellAbsolute",
      "params.inputSections.system.classicalCellDegree",
      "params.inputSections.system.isotropicCell",
      "params.inputSections.system.zFlexibleCell",
      "params.inputSections.system.checkSymmetryPrecision",
      "params.inputSections.system.checkSymmetryOff",
      "params.inputSections.system.cutOffRy",
      "params.inputSections.system.cutoffShape",
      "params.inputSections.system.hfxCutoff",
      "params.inputSections.system.scale",
      "params.inputSections.system.scaleX",
      "params.inputSections.system.scaleY",
      "params.inputSections.system.scaleZ",
      "params.inputSections.system.scaleCartesian",
      "params.inputSections.system.charge",
      "params.inputSections.system.multiplicity",
      "params.inputSections.system.nSup",
      "params.inputSections.system.states",
      "params.inputSections.system.occupation",
      "params.inputSections.system.occupationFixed",
      "params.inputSections.system.externalField",
      "params.inputSections.system.wCut",
      "params.inputSections.system.wGauss",
      "params.inputSections.system.couplingsFiniteDifference",
      "params.inputSections.system.couplingsFiniteDifferenceDisplacement",
      "params.inputSections.system.couplingsProductDisplacement",
      "params.inputSections.system.couplingsLinres",
      "params.inputSections.system.couplingsLinresTolerance",
      "params.inputSections.system.couplingsLinresNvects",
      "params.inputSections.system.couplingsLinresSpecify",
      "params.inputSections.system.couplingsLinresThresholds",
      "params.inputSections.system.couplingsLinresBruteForce",
      "params.inputSections.system.couplingsSurfaces",
      "params.inputSections.system.couplingsFiniteDifferenceAtoms",
      "params.inputSections.system.cdftDonorAtoms",
      "params.inputSections.system.cdftDonorWeights",
      "params.inputSections.system.cdftAcceptorAtoms",
      "params.inputSections.system.cdftAcceptorHdasDonors",
      "params.inputSections.system.cdftAcceptorWeights",
      "params.inputSections.system.pointGroup",
      "params.inputSections.system.pointGroupMolecule",
      "params.inputSections.system.pointGroupDelta",
      "params.inputSections.system.lowSpinExcitation",
      "params.inputSections.system.lowSpinExcitationLsets",
      "params.inputSections.system.lowSpinExcitationPenalty",
      "params.inputSections.system.lseParameters",
      "params.inputSections.system.modifiedGoedecker",
      "params.inputSections.system.modifiedGoedeckerParameters",
      "params.inputSections.system.energyProfile",
      "params.inputSections.system.densityCutOffRy",
      "params.inputSections.system.densityCutoffNumber",
      "params.inputSections.system.dual",
      "params.inputSections.system.constantCutoff",
      "params.inputSections.system.boxWalls",
      "params.inputSections.system.poissonSolver",
      "params.inputSections.system.poissonParameter",
      "params.inputSections.system.mesh",
      "params.inputSections.system.kpoints",
      "params.inputSections.system.kpointsScaled",
      "params.inputSections.system.kpointsOnlyDiagonal",
      "params.inputSections.system.kpointsBlock",
      "params.inputSections.system.kpointsBlockAll",
      "params.inputSections.system.kpointsBlockCalculated",
      "params.inputSections.system.kpointsBlockNoSwap",
      "params.inputSections.system.kpointsMonkhorstPack",
      "params.inputSections.system.kpointsMonkhorstSymmetrized",
      "params.inputSections.system.kpointsMonkhorstFull",
      "params.inputSections.system.kpointsMonkhorstKdp",
      "params.inputSections.system.kpointsMonkhorstShift",
      "params.inputSections.system.kpointBands",
      "params.inputSections.system.doubleGrid",
      "params.inputSections.system.symmetrizeCoordinates",
      "params.inputSections.system.tesr",
      "params.inputSections.system.surface",
      "params.inputSections.system.polymer",
      "params.inputSections.system.cluster",
      "params.inputSections.system.pressure",
      "params.inputSections.system.stressTensor",
      "params.inputSections.system.shockVelocity",
      "params.inputSections.system.directives",
      "params.inputSections.basis.directives",
      "params.inputSections.basis.subsections",
      "params.inputSections.cpmd.optimizeWavefunction",
      "params.inputSections.cpmd.molecularDynamics",
      "params.inputSections.cpmd.convergenceOrbitals",
      "params.inputSections.cpmd.maxStep",
      "params.inputSections.cpmd.timestep",
      "params.inputSections.cpmd.maxRuntime",
      "params.inputSections.cpmd.timestepElectrons",
      "params.inputSections.cpmd.timestepIons",
      "params.inputSections.cpmd.cellMass",
      "params.inputSections.cpmd.convergenceCell",
      "params.inputSections.cpmd.convergenceAdapt",
      "params.inputSections.cpmd.convergenceEnergy",
      "params.inputSections.cpmd.convergenceCalfor",
      "params.inputSections.cpmd.convergenceRelax",
      "params.inputSections.cpmd.convergenceRhofix",
      "params.inputSections.cpmd.convergenceInitial",
      "params.inputSections.cpmd.convergenceConstraint",
      "params.inputSections.cpmd.temperatureElectron",
      "params.inputSections.cpmd.temperature",
      "params.inputSections.cpmd.temperatureRamp",
      "params.inputSections.cpmd.temperatureRampTime",
      "params.inputSections.cpmd.temperatureRampRate",
      "params.inputSections.cpmd.rescaleOldVelocities",
      "params.inputSections.cpmd.reverseVelocities",
      "params.inputSections.cpmd.subtractComVelocity",
      "params.inputSections.cpmd.subtractRotVelocity",
      "params.inputSections.cpmd.prngSeed",
      "params.inputSections.cpmd.tempControlIons",
      "params.inputSections.cpmd.tempControlElectrons",
      "params.inputSections.cpmd.tempControlCell",
      "params.inputSections.cpmd.berendsenIons",
      "params.inputSections.cpmd.berendsenElectrons",
      "params.inputSections.cpmd.berendsenCell",
      "params.inputSections.cpmd.noseIonsThermostat",
      "params.inputSections.cpmd.noseElectronsThermostat",
      "params.inputSections.cpmd.noseCellThermostat",
      "params.inputSections.cpmd.noseParameters",
      "params.inputSections.cpmd.restartWavefunction",
      "params.inputSections.cpmd.trajectory",
      "params.inputSections.cpmd.directives",
      "params.inputSections.cpmd.optimizeGeometry",
      "params.inputSections.cpmd.optimizeGeometryOptions",
      "params.inputSections.cpmd.optimizeGeometrySample",
      "params.inputSections.cpmd.optimizeCombinedOptions",
      "params.inputSections.cpmd.optimizeCombinedSample",
      "params.inputSections.cpmd.maxIter",
      "params.inputSections.cpmd.convergenceGeometry",
      "params.inputSections.cpmd.electronMass",
      "params.inputSections.cpmd.molecularDynamicsCp",
      "params.inputSections.cpmd.molecularDynamicsBo",
      "params.inputSections.cpmd.molecularDynamicsEh",
      "params.inputSections.cpmd.molecularDynamicsPt",
      "params.inputSections.cpmd.molecularDynamicsClassical",
      "params.inputSections.cpmd.molecularDynamicsFile",
      "params.inputSections.cpmd.molecularDynamicsFileOptions",
      "params.inputSections.cpmd.molecularDynamicsBdTrajectories",
      "params.inputSections.cpmd.parrinelloRahmanOptions",
      "params.inputSections.cpmd.cheby",
      "params.inputSections.cpmd.cayley",
      "params.inputSections.cpmd.rungeKutta",
      "params.inputSections.cpmd.forceMatch",
      "params.inputSections.cpmd.nose",
      "params.inputSections.cpmd.noseIons",
      "params.inputSections.cpmd.noseElectrons",
      "params.inputSections.cpmd.berendsen",
      "params.inputSections.cpmd.langevin",
      "params.inputSections.cpmd.annealing",
      "params.inputSections.cpmd.annealingIons",
      "params.inputSections.cpmd.annealingElectrons",
      "params.inputSections.cpmd.annealingCell",
      "params.inputSections.cpmd.dampingIons",
      "params.inputSections.cpmd.dampingElectrons",
      "params.inputSections.cpmd.dampingCell",
      "params.inputSections.cpmd.hessian",
      "params.inputSections.cpmd.project",
      "params.inputSections.cpmd.stressTensorSample",
      "params.inputSections.cpmd.stressTensorVirial",
      "params.inputSections.cpmd.classStressSample",
      "params.inputSections.cpmd.quench",
      "params.inputSections.cpmd.rattle",
      "params.inputSections.cpmd.shake",
      "params.inputSections.cpmd.constraint",
      "params.inputSections.cpmd.trotter",
      "params.inputSections.cpmd.restart",
      "params.inputSections.cpmd.printOptions",
      "params.inputSections.cpmd.storeOptions",
      "params.inputSections.cpmd.storeSelection",
      "params.inputSections.cpmd.storeInterval",
      "params.inputSections.cpmd.storeSelfConsistentInterval",
      "params.inputSections.cpmd.storeOffSelection",
      "params.inputSections.cpmd.restFileCount",
      "params.inputSections.cpmd.restFileSample",
      "params.inputSections.cpmd.trajectoryOptions",
      "params.inputSections.cpmd.trajectorySample",
      "params.inputSections.cpmd.trajectoryRange",
      "params.inputSections.cpmd.movieSample",
      "params.inputSections.cpmd.movieOff",
      "params.inputSections.cpmd.energyBands",
      "params.inputSections.cpmd.externalPotential",
      "params.inputSections.cpmd.externalPotentialAdd",
      "params.inputSections.cpmd.electrostaticPotential",
      "params.inputSections.cpmd.electrostaticPotentialSample",
      "params.inputSections.cpmd.dipoleDynamicsSample",
      "params.inputSections.cpmd.dipoleDynamicsWannier",
      "params.inputSections.cpmd.rhoOut",
      "params.inputSections.cpmd.rhoOutSample",
      "params.inputSections.cpmd.rhoOutBandsCount",
      "params.inputSections.cpmd.rhoOutBands",
      "params.inputSections.cpmd.elf",
      "params.inputSections.cpmd.elfParameters",
      "params.inputSections.cpmd.wannierParameters",
      "params.inputSections.cpmd.wannierOptimization",
      "params.inputSections.cpmd.wannierType",
      "params.inputSections.cpmd.wannierReference",
      "params.inputSections.cpmd.wannierSerial",
      "params.inputSections.cpmd.wannierDos",
      "params.inputSections.cpmd.wannierMolecular",
      "params.inputSections.cpmd.wannierWfnOutOptions",
      "params.inputSections.cpmd.wannierWfnOutPayload",
      "params.inputSections.cpmd.compress",
      "params.inputSections.cpmd.memory",
      "params.inputSections.cpmd.realSpaceWfnKeep",
      "params.inputSections.cpmd.realSpaceWfnSize",
      "params.inputSections.cpmd.splineOptions",
      "params.inputSections.cpmd.splinePoints",
      "params.inputSections.cpmd.splineRange",
      "params.inputSections.cpmd.finiteDifferences",
      "params.inputSections.cpmd.taskGroups",
      "params.inputSections.cpmd.taskGroupsCount",
      "params.inputSections.cpmd.distributeFnl",
      "params.inputSections.cpmd.filePath",
      "params.inputSections.cpmd.benchmark",
      "params.inputSections.cpmd.mirror",
      "params.inputSections.cpmd.shiftPotential",
      "params.inputSections.cpmd.glocalizationParameters",
      "params.inputSections.cpmd.glocalizationOptimization",
      "params.inputSections.cpmd.gfunctionalType",
      "params.inputSections.cpmd.spreadRspace",
      "params.inputSections.cpmd.gUnitarityOptions",
      "params.inputSections.cpmd.stepTuning",
      "params.inputSections.cpmd.gAntisym",
      "params.inputSections.cpmd.gAntisymPenalty",
      "params.inputSections.cpmd.gKick",
      "params.inputSections.cpmd.gComplex",
      "params.inputSections.cpmd.gReal",
      "params.inputSections.cpmd.readMatrix",
      "params.inputSections.cpmd.gStepTune",
      "params.inputSections.cpmd.glocWfnOutOptions",
      "params.inputSections.cpmd.glocWfnOutPayload",
      "params.inputSections.cpmd.noGeoCheck",
      "params.inputSections.cpmd.brokenSymmetry",
      "params.inputSections.cpmd.distributedLinalg",
      "params.inputSections.cpmd.linalgNewOrtho",
      "params.inputSections.cpmd.disorthoBlockSize",
      "params.inputSections.cpmd.statesBlockSize",
      "params.inputSections.cpmd.allToAllPrecision",
      "params.inputSections.cpmd.gshell",
      "params.inputSections.cpmd.localPotential",
      "params.inputSections.cpmd.centerMoleculeOff",
      "params.inputSections.cpmd.centerMoleculeOn",
      "params.inputSections.cpmd.diis",
      "params.inputSections.cpmd.odiis",
      "params.inputSections.cpmd.pcg",
      "params.inputSections.cpmd.diagonalization",
      "params.inputSections.cpmd.freeEnergy",
      "params.inputSections.cpmd.interface",
      "params.inputSections.cpmd.qmmm",
      "params.inputSections.cpmd.bicanonicalEnsemble",
      "params.inputSections.cpmd.cdft",
      "params.inputSections.cpmd.cdftOptions",
      "params.inputSections.cpmd.cdftPayload",
      "params.inputSections.cpmd.cdftHdaPayload",
      "params.inputSections.cpmd.vgfactor",
      "params.inputSections.cpmd.vMirror",
      "params.inputSections.cpmd.combineSystemsOptions",
      "params.inputSections.cpmd.combineSystemsPayload",
      "params.inputSections.cpmd.combineSystemsSabPayload",
      "params.inputSections.cpmd.kshamOptions",
      "params.inputSections.cpmd.kshamPayload",
      "params.inputSections.cpmd.czonesSet",
      "params.inputSections.cpmd.woutOptions",
      "params.inputSections.cpmd.woutPayload",
      "params.inputSections.cpmd.xfmqcTrajectories",
      "params.inputSections.cpmd.debugOptions",
      "params.inputSections.cpmd.kohnShamEnergiesOptions",
      "params.inputSections.cpmd.kohnShamEnergiesCount",
      "params.inputSections.cpmd.surfaceHoppingOptions",
      "params.inputSections.cpmd.roksOptions",
      "params.inputSections.cpmd.roksExpertPayload",
      "params.inputSections.cpmd.pathSampling",
      "params.inputSections.cpmd.fixrhoUpwfnOptions",
      "params.inputSections.cpmd.fixrhoVectors",
      "params.inputSections.cpmd.fixrhoLoop",
      "params.inputSections.cpmd.fixrhoWftol",
      "params.inputSections.cpmd.bogoliubovCorrection",
      "params.inputSections.cpmd.vibrationalAnalysisOptions",
      "params.inputSections.cpmd.vibrationalAnalysisSample",
      "params.inputSections.cpmd.vibrationalAnalysisMode",
      "params.inputSections.cpmd.electronicSpectra",
      "params.inputSections.cpmd.spinOrbitCouplingStates",
      "params.inputSections.cpmd.propagationSpectra",
      "params.inputSections.cpmd.propagationDistrub",
      "params.inputSections.cpmd.gaugePulse",
      "params.inputSections.cpmd.gaugeFieldFrequency",
      "params.inputSections.cpmd.nacv",
      "params.inputSections.cpmd.orbitalHardnessOptions",
      "params.inputSections.cpmd.pathIntegral",
      "params.inputSections.cpmd.pathMinimization",
      "params.inputSections.cpmd.langevinOptions",
      "params.inputSections.cpmd.langevinParameter",
      "params.inputSections.cpmd.qmmmEasy",
      "params.inputSections.cpmd.interfaceOptions",
      "params.inputSections.cpmd.trotterFactorCount",
      "params.inputSections.cpmd.trotterFactorPayload",
      "params.inputSections.cpmd.linearResponse",
      "params.inputSections.cpmd.harmonicReference",
      "params.inputSections.cpmd.scaledMasses",
      "params.inputSections.cpmd.tddft",
      "params.inputSections.cpmd.ssic",
      "params.inputSections.cpmd.nonorthogonalOrbitalsOptions",
      "params.inputSections.cpmd.nonorthogonalOrbitalsLimit",
      "params.inputSections.cpmd.lanczosDiagonalizationOptions",
      "params.inputSections.cpmd.lanczosParametersCount",
      "params.inputSections.cpmd.lanczosParametersPayload",
      "params.inputSections.cpmd.davidsonDiagonalization",
      "params.inputSections.cpmd.davidsonParameters",
      "params.inputSections.cpmd.alexanderMixing",
      "params.inputSections.cpmd.andersonMixingGspace",
      "params.inputSections.cpmd.andersonMixingCount",
      "params.inputSections.cpmd.andersonMixingPayload",
      "params.inputSections.cpmd.broydenMixingOptions",
      "params.inputSections.cpmd.broydenMixingPayload",
      "params.inputSections.cpmd.diisMixingCount",
      "params.inputSections.cpmd.diisMixingPayload",
      "params.inputSections.cpmd.moverhoMixing",
      "params.inputSections.cpmd.extrapolateWfnOptions",
      "params.inputSections.cpmd.extrapolateWfnOrder",
      "params.inputSections.cpmd.extrapolateConstraintOrder",
      "params.inputSections.cpmd.properties",
      "params.inputSections.cpmd.vdwCorrection",
      "params.inputSections.cpmd.vdwWannier",
      "params.inputSections.cpmd.dcacp",
      "params.inputSections.cpmd.isolatedMolecule",
      "params.inputSections.dft.functional",
      "params.inputSections.dft.lsd",
      "params.inputSections.dft.directives",
      "params.inputSections.dft.gcCutoff",
      "params.inputSections.dft.xcDriver",
      "params.inputSections.dft.libxc",
      "params.inputSections.dft.lrKernel",
      "params.inputSections.dft.refunct",
      "params.inputSections.dft.mtsHighFunc",
      "params.inputSections.dft.mtsLowFunc",
      "params.inputSections.dft.hfx",
      "params.inputSections.dft.hfxScreening",
      "params.inputSections.dft.hubbard",
      "params.inputSections.dft.alpha",
      "params.inputSections.dft.beta",
      "params.inputSections.dft.oldCode",
      "params.inputSections.dft.newCode",
      "params.inputSections.dft.correlation",
      "params.inputSections.dft.exchange",
      "params.inputSections.dft.becke88",
      "params.inputSections.clas.directives",
      "params.inputSections.clas.subsections",
      "params.inputSections.eam.directives",
      "params.inputSections.eam.subsections",
      "params.inputSections.exte.directives",
      "params.inputSections.exte.subsections",
      "params.inputSections.hardness.directives",
      "params.inputSections.hardness.subsections",
      "params.inputSections.info.directives",
      "params.inputSections.info.subsections",
      "params.inputSections.linres.directives",
      "params.inputSections.linres.subsections",
      "params.inputSections.molstates.directives",
      "params.inputSections.molstates.subsections",
      "params.inputSections.mts.directives",
      "params.inputSections.mts.subsections",
      "params.inputSections.nlcc.directives",
      "params.inputSections.nlcc.subsections",
      "params.inputSections.path.directives",
      "params.inputSections.path.subsections",
      "params.inputSections.pimd.directives",
      "params.inputSections.pimd.subsections",
      "params.inputSections.potential.directives",
      "params.inputSections.potential.subsections",
      "params.inputSections.prop.directives",
      "params.inputSections.prop.subsections",
      "params.inputSections.ptddft.directives",
      "params.inputSections.ptddft.subsections",
      "params.inputSections.resp.directives",
      "params.inputSections.resp.subsections",
      "params.inputSections.tddft.directives",
      "params.inputSections.tddft.subsections",
      "params.inputSections.vdw.directives",
      "params.inputSections.vdw.subsections",
      "params.inputSections.vectors.directives",
      "params.inputSections.vectors.subsections",
      "params.inputSections.wavefunction.directives",
      "params.inputSections.wavefunction.subsections",
      "params.inputSections.atoms.pseudopotentials",
      "params.inputSections.atoms.directives",
      "params.inputSections.set.key",
      "params.inputSections.set.value",
      "params.inputSections.raw",
  };

  for (size_t i = 0; i < sizeof(param_features) / sizeof(param_features[0]);
       i++) {
    const CPMDCFeatureEntry *entry = cpmdc_feature_find(param_features[i]);
    assert_non_null(entry);
    assert_string_equal(entry->feature_id, param_features[i]);
    assert_int_equal(entry->kind, CPMDC_FEATURE_PARAMS);
    assert_int_equal(entry->stub_applicable, 1);
    assert_int_equal(entry->embed_applicable, 1);
  }
}

int main(void) {
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(test_feature_table_nonempty),
      cmocka_unit_test(test_inscan_sections),
      cmocka_unit_test(test_cp_keywords_not_sections),
      cmocka_unit_test(test_abi),
      cmocka_unit_test(test_param_applicability),
      cmocka_unit_test(test_structured_param_features),
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
