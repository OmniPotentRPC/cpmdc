#include "cpmdc_features.h"
#include <string.h>

static const CPMDCFeatureEntry g_features[] = {
    {"section.generic", CPMDC_FEATURE_SECTION, 1, 1},
    {"section.system", CPMDC_FEATURE_SECTION, 1, 1},
    {"section.cpmd", CPMDC_FEATURE_SECTION, 1, 1},
    {"section.dft", CPMDC_FEATURE_SECTION, 1, 1},
    {"section.atoms", CPMDC_FEATURE_SECTION, 1, 1},
    {"section.set", CPMDC_FEATURE_SECTION, 1, 0},
    {"section.raw", CPMDC_FEATURE_SECTION, 1, 1},
    {"params.functional", CPMDC_FEATURE_PARAMS, 1, 1},
    {"params.cutOffRy", CPMDC_FEATURE_PARAMS, 1, 1},
    {"params.charge", CPMDC_FEATURE_PARAMS, 1, 1},
    {"params.multiplicity", CPMDC_FEATURE_PARAMS, 1, 1},
    {"params.task", CPMDC_FEATURE_PARAMS, 1, 1},
    {"params.title", CPMDC_FEATURE_PARAMS, 1, 1},
    {"params.memoryMb", CPMDC_FEATURE_PARAMS, 1, 1},
    {"params.scratchDir", CPMDC_FEATURE_PARAMS, 1, 1},
    {"params.permanentDir", CPMDC_FEATURE_PARAMS, 1, 1},
    {"params.cpmdRoot", CPMDC_FEATURE_PARAMS, 1, 1},
    {"params.enginePath", CPMDC_FEATURE_PARAMS, 1, 1},
    {"params.inputBlocks", CPMDC_FEATURE_PARAMS, 1, 1},
    {"params.inputSections", CPMDC_FEATURE_PARAMS, 1, 1},
    {"abi.cpmdc_set_params", CPMDC_FEATURE_ABI, 1, 1},
    {"abi.cpmdc_energy_gradient", CPMDC_FEATURE_ABI, 0, 1},
    {"abi.cpmdc_energy", CPMDC_FEATURE_ABI, 0, 1},
    {"abi.cpmdc_energy_forces", CPMDC_FEATURE_ABI, 0, 1},
    {"abi.cpmdc_session_create", CPMDC_FEATURE_ABI, 1, 1},
    {"abi.cpmdc_session_set_params", CPMDC_FEATURE_ABI, 1, 1},
    {"abi.cpmdc_session_destroy", CPMDC_FEATURE_ABI, 1, 1},
    {"abi.cpmdc_session_energy_gradient", CPMDC_FEATURE_ABI, 0, 1},
    {"abi.cpmdc_session_energy", CPMDC_FEATURE_ABI, 0, 1},
    {"abi.cpmdc_session_energy_forces", CPMDC_FEATURE_ABI, 0, 1},
    {"abi.cpmdc_session_calculate_forces", CPMDC_FEATURE_ABI, 0, 1},
    {"abi.cpmdc_session_calculate_result", CPMDC_FEATURE_ABI, 0, 1},
    {"abi.cpmdc_calculate_result", CPMDC_FEATURE_ABI, 0, 1},
    {"abi.cpmdc_potential_result_size_for_force_input", CPMDC_FEATURE_ABI, 1, 1},
    {"abi.cpmdc_version", CPMDC_FEATURE_ABI, 1, 1},
    {"abi.cpmdc_available", CPMDC_FEATURE_ABI, 1, 1},
    {"abi.cpmdc_finalize", CPMDC_FEATURE_ABI, 1, 1},
};

size_t cpmdc_feature_count(void) {
  return sizeof(g_features) / sizeof(g_features[0]);
}

const CPMDCFeatureEntry *cpmdc_feature_table(void) { return g_features; }

const CPMDCFeatureEntry *cpmdc_feature_find(const char *feature_id) {
  if (!feature_id)
    return NULL;
  size_t n = cpmdc_feature_count();
  for (size_t i = 0; i < n; ++i) {
    if (strcmp(g_features[i].feature_id, feature_id) == 0)
      return &g_features[i];
  }
  return NULL;
}
