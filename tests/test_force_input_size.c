#include "cpmdc.h"

#include <errno.h>
#include <setjmp.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cmocka.h>

static const char *g_step_a = NULL;
static const char *g_step_ev = NULL;

static unsigned char *read_file(const char *path, size_t *size) {
  FILE *fp = fopen(path, "rb");
  if (!fp)
    return NULL;
  if (fseek(fp, 0, SEEK_END) != 0) {
    fclose(fp);
    return NULL;
  }
  long n = ftell(fp);
  if (n <= 0) {
    fclose(fp);
    return NULL;
  }
  rewind(fp);
  unsigned char *buf = (unsigned char *)malloc((size_t)n);
  if (!buf) {
    fclose(fp);
    return NULL;
  }
  if (fread(buf, 1, (size_t)n, fp) != (size_t)n) {
    free(buf);
    fclose(fp);
    return NULL;
  }
  fclose(fp);
  *size = (size_t)n;
  return buf;
}

static void test_sizes_match_atom_count(void **state) {
  (void)state;
  size_t a_size = 0, ev_size = 0;
  unsigned char *a = read_file(g_step_a, &a_size);
  unsigned char *ev = read_file(g_step_ev, &ev_size);
  assert_non_null(a);
  assert_non_null(ev);

  size_t need_a = cpmdc_potential_result_size_for_force_input(a, a_size);
  size_t need_ev = cpmdc_potential_result_size_for_force_input(ev, ev_size);
  /* 2 atoms => 6 force components => 32 + 6*8 = 80 bytes flat layout estimate */
  assert_int_equal(need_a, 80);
  assert_int_equal(need_ev, 80);
  assert_int_equal(cpmdc_potential_result_size_for_force_input(NULL, 0), 0);

  free(a);
  free(ev);
}

int main(int argc, char **argv) {
  if (argc < 3) {
    fprintf(stderr, "usage: %s step_a.bin step_ev.bin\n", argv[0]);
    return 2;
  }
  g_step_a = argv[1];
  g_step_ev = argv[2];
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(test_sizes_match_atom_count),
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
