#include "cpmdc.h"

#include <dlfcn.h>
#include <stdio.h>
#include <string.h>

typedef int (*cpmdc_available_fn)(void);
typedef const char *(*cpmdc_version_fn)(void);
typedef void (*cpmdc_finalize_fn)(void);
typedef size_t (*cpmdc_feature_count_fn)(void);
typedef const CPMDCFeatureEntry *(*cpmdc_feature_table_fn)(void);
typedef const CPMDCFeatureEntry *(*cpmdc_feature_find_fn)(const char *);

static void *required_symbol(void *handle, const char *name) {
  dlerror();
  void *sym = dlsym(handle, name);
  const char *err = dlerror();
  if (err != NULL || sym == NULL) {
    fprintf(stderr, "missing symbol %s: %s\n", name, err ? err : "null");
    return NULL;
  }
  return sym;
}

#define LOAD_REQUIRED_FUNCTION(handle, name, type)                             \
  ((type)required_symbol((handle), (name)))

static const char *const required_abi_symbols[] = {
    "cpmdc_set_params",
    "cpmdc_energy_gradient",
    "cpmdc_energy",
    "cpmdc_energy_forces",
    "cpmdc_session_create",
    "cpmdc_session_set_params",
    "cpmdc_session_destroy",
    "cpmdc_session_energy_gradient",
    "cpmdc_session_energy",
    "cpmdc_session_energy_forces",
    "cpmdc_session_calculate_forces",
    "cpmdc_session_calculate_result",
    "cpmdc_calculate_result",
    "cpmdc_potential_result_size_for_force_input",
    "cpmdc_version",
    "cpmdc_available",
    "cpmdc_finalize",
    "cpmdc_feature_count",
    "cpmdc_feature_table",
    "cpmdc_feature_find",
};

int main(int argc, char **argv) {
  if (argc != 2) {
    fprintf(stderr, "usage: %s /path/to/libcpmdc.so\n", argv[0]);
    return 1;
  }

  void *handle = dlopen(argv[1], RTLD_NOW | RTLD_LOCAL);
  if (handle == NULL) {
    fprintf(stderr, "dlopen failed for %s: %s\n", argv[1], dlerror());
    return 1;
  }

  for (size_t i = 0;
       i < sizeof(required_abi_symbols) / sizeof(required_abi_symbols[0]);
       ++i) {
    if (required_symbol(handle, required_abi_symbols[i]) == NULL) {
      dlclose(handle);
      return 1;
    }
  }

#if defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
#endif
  cpmdc_version_fn version = LOAD_REQUIRED_FUNCTION(
      handle, "cpmdc_version", cpmdc_version_fn);
  cpmdc_available_fn available = LOAD_REQUIRED_FUNCTION(
      handle, "cpmdc_available", cpmdc_available_fn);
  cpmdc_finalize_fn finalize = LOAD_REQUIRED_FUNCTION(
      handle, "cpmdc_finalize", cpmdc_finalize_fn);
  cpmdc_feature_count_fn feature_count = LOAD_REQUIRED_FUNCTION(
      handle, "cpmdc_feature_count", cpmdc_feature_count_fn);
  cpmdc_feature_table_fn feature_table = LOAD_REQUIRED_FUNCTION(
      handle, "cpmdc_feature_table", cpmdc_feature_table_fn);
  cpmdc_feature_find_fn feature_find = LOAD_REQUIRED_FUNCTION(
      handle, "cpmdc_feature_find", cpmdc_feature_find_fn);
#if defined(__GNUC__)
#pragma GCC diagnostic pop
#endif
  if (version == NULL || available == NULL || finalize == NULL ||
      feature_count == NULL || feature_table == NULL || feature_find == NULL) {
    dlclose(handle);
    return 1;
  }
  const CPMDCFeatureEntry *feature =
      feature_find("abi.cpmdc_session_calculate_result");
  if (feature_count() == 0 || feature_table() == NULL || feature == NULL ||
      feature->kind != CPMDC_FEATURE_ABI || feature->embed_applicable != 1) {
    fprintf(stderr, "shared cpmdc feature discovery is incomplete\n");
    dlclose(handle);
    return 1;
  }

  const char *ver = version();
  if (ver == NULL || strstr(ver, "cpmdc/0.1.0") == NULL) {
    fprintf(stderr, "unexpected cpmdc version: %s\n", ver ? ver : "(null)");
    dlclose(handle);
    return 1;
  }
  if (available() != 1) {
    fprintf(stderr, "shared cpmdc engine should report available\n");
    dlclose(handle);
    return 1;
  }
  finalize();
  if (available() != 0) {
    fprintf(stderr, "cpmdc_finalize did not mark runtime unavailable\n");
    dlclose(handle);
    return 1;
  }

  dlclose(handle);
  return 0;
}
