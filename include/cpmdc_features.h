#pragma once
#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif

/**
 * @file cpmdc_features.h
 * @brief Runtime feature discovery for the cpmdc C ABI.
 *
 * Feature IDs are stable strings used by embedders to discover available ABI
 * calls, CPMD catalog entries, and Cap'n Proto parameter fields. Examples:
 * `abi.cpmdc_session_calculate_result`, `catalog.section.PIMD`,
 * `catalog.cpmd.MAXITER`, and `params.inputSections.cpmd.maxIter`.
 */

/** Feature namespace represented by a `CPMDCFeatureEntry`. */
typedef enum CPMDCFeatureKind {
  /** CPMD input section support, such as `catalog.section.PIMD`. */
  CPMDC_FEATURE_SECTION = 1,
  /** Cap'n Proto `CPMDParams` field support. */
  CPMDC_FEATURE_PARAMS = 2,
  /** Exported C ABI symbol support. */
  CPMDC_FEATURE_ABI = 3,
  /** CPMD keyword support inside a catalog section. */
  CPMDC_FEATURE_KEYWORD = 4
} CPMDCFeatureKind;

/**
 * @brief One discoverable ABI, parameter, section, or keyword capability.
 */
typedef struct CPMDCFeatureEntry {
  /** Stable feature ID string. */
  const char *feature_id;
  /** Feature namespace. */
  CPMDCFeatureKind kind;
  /** Non-zero when the standalone stub build exposes the feature. */
  int stub_applicable;
  /** Non-zero when an embedded OpenCPMD build exposes the feature. */
  int embed_applicable;
} CPMDCFeatureEntry;

/** @brief Number of entries returned by `cpmdc_feature_table()`. */
size_t cpmdc_feature_count(void);

/**
 * @brief Return the contiguous feature table.
 *
 * The returned pointer has `cpmdc_feature_count()` entries and remains owned by
 * the library.
 */
const CPMDCFeatureEntry *cpmdc_feature_table(void);

/**
 * @brief Find one feature by stable ID.
 *
 * @return Pointer to the matching table entry, or `NULL` when the feature ID is
 * not present.
 */
const CPMDCFeatureEntry *cpmdc_feature_find(const char *feature_id);

#ifdef __cplusplus
}
#endif
