#include <stddef.h>
#include <stdlib.h>

static void *volatile leak_sink;

__attribute__((constructor)) static void leak_on_dlopen(void) {
  void *ptr = malloc(77);
  leak_sink = ptr;
  leak_sink = NULL;
}
