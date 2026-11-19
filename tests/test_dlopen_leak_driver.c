#include <dlfcn.h>
#include <stdio.h>

int main(int argc, char **argv) {
  if (argc != 2) {
    fprintf(stderr, "usage: %s MODULE_PATH\n", argv[0]);
    return 2;
  }
  void *handle = dlopen(argv[1], RTLD_NOW | RTLD_LOCAL);
  if (handle == NULL) {
    fprintf(stderr, "dlopen failed for %s: %s\n", argv[1], dlerror());
    return 2;
  }
  dlclose(handle);
  return 0;
}
