#pragma once
#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif

typedef enum CPMDCFeatureKind {
  CPMDC_FEATURE_SECTION = 1,
  CPMDC_FEATURE_PARAMS = 2,
  CPMDC_FEATURE_ABI = 3,
  CPMDC_FEATURE_KEYWORD = 4
} CPMDCFeatureKind;

typedef struct CPMDCFeatureEntry {
  const char *feature_id;
  CPMDCFeatureKind kind;
  int stub_applicable;
  int embed_applicable;
} CPMDCFeatureEntry;

size_t cpmdc_feature_count(void);
const CPMDCFeatureEntry *cpmdc_feature_table(void);
const CPMDCFeatureEntry *cpmdc_feature_find(const char *feature_id);

#ifdef __cplusplus
}
#endif
