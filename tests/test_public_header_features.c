#include "cpmdc.h"

#include <stddef.h>

int main(void) {
  size_t count = cpmdc_feature_count();
  const CPMDCFeatureEntry *table = cpmdc_feature_table();
  const CPMDCFeatureEntry *entry = cpmdc_feature_find("abi.cpmdc_feature_find");
  return count > 0 && table != NULL && entry != NULL ? 0 : 1;
}
