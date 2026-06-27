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
      "params.inputSections.system.cutOffRy",
      "params.inputSections.system.scale",
      "params.inputSections.system.charge",
      "params.inputSections.system.multiplicity",
      "params.inputSections.system.directives",
      "params.inputSections.basis.directives",
      "params.inputSections.basis.subsections",
      "params.inputSections.cpmd.optimizeWavefunction",
      "params.inputSections.cpmd.molecularDynamics",
      "params.inputSections.cpmd.convergenceOrbitals",
      "params.inputSections.cpmd.maxStep",
      "params.inputSections.cpmd.timestep",
      "params.inputSections.cpmd.restartWavefunction",
      "params.inputSections.cpmd.trajectory",
      "params.inputSections.cpmd.directives",
      "params.inputSections.cpmd.optimizeGeometry",
      "params.inputSections.cpmd.maxIter",
      "params.inputSections.cpmd.convergenceGeometry",
      "params.inputSections.cpmd.electronMass",
      "params.inputSections.cpmd.molecularDynamicsCp",
      "params.inputSections.cpmd.molecularDynamicsBo",
      "params.inputSections.cpmd.molecularDynamicsEh",
      "params.inputSections.cpmd.molecularDynamicsPt",
      "params.inputSections.cpmd.molecularDynamicsClassical",
      "params.inputSections.cpmd.molecularDynamicsFile",
      "params.inputSections.cpmd.nose",
      "params.inputSections.cpmd.noseIons",
      "params.inputSections.cpmd.noseElectrons",
      "params.inputSections.cpmd.berendsen",
      "params.inputSections.cpmd.langevin",
      "params.inputSections.cpmd.annealing",
      "params.inputSections.cpmd.quench",
      "params.inputSections.cpmd.rattle",
      "params.inputSections.cpmd.shake",
      "params.inputSections.cpmd.constraint",
      "params.inputSections.cpmd.trotter",
      "params.inputSections.cpmd.restart",
      "params.inputSections.cpmd.printOptions",
      "params.inputSections.cpmd.storeOptions",
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
      "params.inputSections.cpmd.properties",
      "params.inputSections.cpmd.vdwCorrection",
      "params.inputSections.cpmd.vdwWannier",
      "params.inputSections.cpmd.dcacp",
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
