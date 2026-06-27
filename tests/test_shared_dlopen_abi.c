#include <dlfcn.h>
#include <stdio.h>
#include <string.h>

typedef int (*cpmdc_available_fn)(void);
typedef const char *(*cpmdc_version_fn)(void);
typedef void (*cpmdc_finalize_fn)(void);

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
#if defined(__GNUC__)
#pragma GCC diagnostic pop
#endif
  if (version == NULL || available == NULL || finalize == NULL) {
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
