#define _GNU_SOURCE

#include "cpmdc.h"

#include <dlfcn.h>
#ifdef __linux__
#include <fcntl.h>
#include <limits.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef __linux__
#include <sys/mman.h>
#include <unistd.h>
#endif

typedef int (*cpmdc_available_fn)(void);
typedef const char *(*cpmdc_version_fn)(void);
typedef void (*cpmdc_finalize_fn)(void);
typedef size_t (*cpmdc_feature_count_fn)(void);
typedef const CPMDCFeatureEntry *(*cpmdc_feature_table_fn)(void);
typedef const CPMDCFeatureEntry *(*cpmdc_feature_find_fn)(const char *);

/* 1 only when this binary is the OpenCPMD link. The reference build is 0. */
#if defined(CPMDC_HAS_CPMD)
#define CPMDC_AVAILABLE_WHEN_LINKED 1
#else
#define CPMDC_AVAILABLE_WHEN_LINKED 0
#endif

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

#ifdef __linux__
/* dlclose unmaps libibverbs before LSan symbolizes its allocation stack.
 * Preserve just enough of the file mapping for the library-only suppression
 * to match without suppressing other allocations made by this test. */
static int preserve_lsan_symbol_map(const char *dso, void **base, size_t *size,
                                    char path[PATH_MAX]) {
  void *handle = dlopen(dso, RTLD_NOW | RTLD_LOCAL | RTLD_NOLOAD);
  void *symbol = handle ? dlsym(handle, "verbs_register_driver_59") : NULL;
  Dl_info info;
  long page_size = sysconf(_SC_PAGESIZE);
  if (symbol == NULL || dladdr(symbol, &info) == 0 || info.dli_fname == NULL ||
      strlen(info.dli_fname) >= PATH_MAX || page_size <= 0) {
    if (handle != NULL) {
      dlclose(handle);
    }
    return -1;
  }
  *base = info.dli_fbase;
  size_t symbol_end = (size_t)((char *)symbol - (char *)info.dli_fbase) + 1;
  *size = (symbol_end + (size_t)page_size - 1) & ~((size_t)page_size - 1);
  strcpy(path, info.dli_fname);
  dlclose(handle);
  return 0;
}

static int restore_lsan_symbol_map(void *base, size_t size, const char *path) {
  int fd = open(path, O_RDONLY);
  if (fd < 0) {
    return -1;
  }
  void *mapping = mmap(base, size, PROT_READ | PROT_EXEC,
                       MAP_PRIVATE | MAP_FIXED_NOREPLACE, fd, 0);
  close(fd);
  return mapping == base ? 0 : -1;
}
#endif

#define LOAD_REQUIRED_FUNCTION(handle, name, type)                             \
  ((type)required_symbol((handle), (name)))

static const char *const required_abi_symbols[] = {
    "cpmdc_set_params",
    "cpmdc_energy_gradient",
    "cpmdc_energy",
    "cpmdc_energy_forces",
    "cpmdc_bind_calculator",
    "cpmdc_adopt_calculator_comm",
    "cpmdc_adopted_comm",
    "cpmdc_session_create",
    "cpmdc_session_set_params",
    "cpmdc_session_set_constant_cutoff",
    "cpmdc_session_destroy",
    "cpmdc_session_select_orbitals",
    "cpmdc_session_copy_orbitals",
    "cpmdc_session_energy_gradient",
    "cpmdc_session_energy",
    "cpmdc_session_energy_forces",
    "cpmdc_last_energy_components",
    "cpmdc_last_property_snapshot",
    "cpmdc_last_gradient",
    "cpmdc_last_dipole",
    "cpmdc_last_stress",
    "cpmdc_last_rank_spread",
    "cpmdc_last_md_trajectory_row",
    "cpmdc_last_multi_state_energies",
    "cpmdc_last_charge_integrals",
    "cpmdc_session_last_energy_components",
    "cpmdc_session_last_charge_integrals",
    "cpmdc_session_last_multi_state_energies",
    "cpmdc_session_last_md_trajectory_row",
    "cpmdc_session_last_property_snapshot",
    "cpmdc_session_last_gradient",
    "cpmdc_session_last_dipole",
    "cpmdc_session_last_stress",
    "cpmdc_session_calculate_forces",
    "cpmdc_session_calculate_result",
    "cpmdc_session_calculate_au",
    "cpmdc_calculate_result",
    "cpmdc_calculate_result_from_config",
    "cpmdc_potential_result_size_for_force_input",
    "cpmdc_capabilities_result",
    "cpmdc_version",
    "cpmdc_abi_version",
    "cpmdc_configure",
    "cpmdc_session_create_from_config",
    "cpmdc_session_configure",
    "cpmdc_last_error",
    "cpmdc_available",
    "cpmdc_finalize",
    "cpmdc_feature_count",
    "cpmdc_feature_table",
    "cpmdc_feature_find",
};

/* Internal glue from each kind of source: C in libcpmdc, Fortran BIND(C) in
 * libcpmdc, and Fortran BIND(C) in a static archive linked into it. */
static const char *const hidden_symbols[] = {
    "cpmdc_params_render_input_deck",
    "cpmdc_stop_arm",
    "cpmdc_note_stop",
    "cpmdc_embed_init",
    "cpmdc_embed_energy_grad",
    "cpmdc_embed_apply_params",
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
  for (size_t i = 0; i < sizeof(hidden_symbols) / sizeof(hidden_symbols[0]);
       ++i) {
    if (dlsym(handle, hidden_symbols[i]) != NULL) {
      fprintf(stderr, "internal symbol %s is exported\n", hidden_symbols[i]);
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
  if (ver == NULL || strstr(ver, "cpmdc/") == NULL) {
    fprintf(stderr, "unexpected cpmdc version: %s\n", ver ? ver : "(null)");
    dlclose(handle);
    return 1;
  }
  int ready = available();
  if (ready != CPMDC_AVAILABLE_WHEN_LINKED) {
    fprintf(stderr, "cpmdc_available() returned %d\n", ready);
    dlclose(handle);
    return 1;
  }
  finalize();
  if (available() != 0) {
    fprintf(stderr, "cpmdc_finalize did not mark runtime unavailable\n");
    dlclose(handle);
    return 1;
  }

#ifdef __linux__
  const char *symbolize_dso = getenv("CPMDC_LSAN_SYMBOLIZE_DSO");
  void *symbolize_base = NULL;
  size_t symbolize_size = 0;
  char symbolize_path[PATH_MAX];
  int have_symbol_map =
      symbolize_dso != NULL && symbolize_dso[0] != '\0' &&
      preserve_lsan_symbol_map(symbolize_dso, &symbolize_base, &symbolize_size,
                               symbolize_path) == 0;
  dlclose(handle);
  if (have_symbol_map) {
    if (restore_lsan_symbol_map(symbolize_base, symbolize_size,
                                symbolize_path) != 0) {
      fprintf(stderr, "could not restore symbols for %s\n", symbolize_dso);
      return 1;
    }
  }
#else
  dlclose(handle);
#endif
  return 0;
}
