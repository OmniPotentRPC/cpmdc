#pragma once
#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif

typedef enum CPMDCFeatureKind {
  CPMDC_FEATURE_SECTION = 1,
  CPMDC_FEATURE_PARAMS = 2,
  CPMDC_FEATURE_ABI = 3
} CPMDCFeatureKind;

typedef struct CPMDCFeatureEntry {
  const char *feature_id; /* e.g. section.generic, params.cutOffRy, abi.cpmdc_energy */
  CPMDCFeatureKind kind;
  int stub_applicable;
  int embed_applicable;
} CPMDCFeatureEntry;

/** Number of interned features. */
size_t cpmdc_feature_count(void);

/** Pointer to intern table (length cpmdc_feature_count()). */
const CPMDCFeatureEntry *cpmdc_feature_table(void);

/** Find by feature_id; NULL if unknown. */
const CPMDCFeatureEntry *cpmdc_feature_find(const char *feature_id);

#ifdef __cplusplus
}
#endif
