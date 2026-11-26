/**
 * CPMDCPropertySnapshot keeps the layout existing callers allocate.
 *
 * The offsets below are the LP64 layout of that struct: one int, then
 * size_t hessian_count, double hessian[4096], size_t dipole_count,
 * double dipole[3], size_t polarizability_count, and double
 * polarizability[9]. A caller buffer of that size, with a guard after
 * it, must stay intact when either snapshot accessor copies into it.
 */
#include "cpmdc.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <cmocka.h>

enum { k_snapshot_bytes = 32896, k_guard = 32, k_fill = 0xa5 };

typedef struct GuardedSnapshot {
  _Alignas(CPMDCPropertySnapshot) unsigned char raw[k_snapshot_bytes + k_guard];
} GuardedSnapshot;

_Static_assert(sizeof(CPMDCPropertySnapshot) == k_snapshot_bytes,
               "CPMDCPropertySnapshot size");
_Static_assert(offsetof(CPMDCPropertySnapshot, valid) == 0, "valid");
_Static_assert(offsetof(CPMDCPropertySnapshot, hessian_count) == 8,
               "hessian_count");
_Static_assert(offsetof(CPMDCPropertySnapshot, hessian) == 16, "hessian");
_Static_assert(offsetof(CPMDCPropertySnapshot, dipole_count) == 32784,
               "dipole_count");
_Static_assert(offsetof(CPMDCPropertySnapshot, dipole) == 32792, "dipole");
_Static_assert(offsetof(CPMDCPropertySnapshot, polarizability_count) == 32816,
               "polarizability_count");
_Static_assert(offsetof(CPMDCPropertySnapshot, polarizability) == 32824,
               "polarizability");
_Static_assert(offsetof(CPMDCPropertySnapshot, polarizability) + 9 * sizeof(double) ==
                   k_snapshot_bytes,
               "snapshot tail");
_Static_assert(_Alignof(CPMDCPropertySnapshot) == 8, "snapshot alignment");

static const char *g_params = NULL;

static unsigned char *read_file(const char *path, size_t *size) {
  FILE *fp = fopen(path, "rb");
  long n;
  unsigned char *buf;
  if (!fp)
    return NULL;
  if (fseek(fp, 0, SEEK_END) != 0) {
    fclose(fp);
    return NULL;
  }
  n = ftell(fp);
  if (n <= 0) {
    fclose(fp);
    return NULL;
  }
  rewind(fp);
  buf = (unsigned char *)malloc((size_t)n);
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

static void fill_guard(GuardedSnapshot *buf) {
  memset(buf->raw, k_fill, sizeof(buf->raw));
}

static void assert_guard(const GuardedSnapshot *buf) {
  size_t i;
  for (i = k_snapshot_bytes; i < sizeof(buf->raw); ++i)
    assert_int_equal(buf->raw[i], k_fill);
}

static void test_process_snapshot_stays_inside_old_buffer(void **state) {
  GuardedSnapshot buf;
  (void)state;
  fill_guard(&buf);
  (void)cpmdc_last_property_snapshot((CPMDCPropertySnapshot *)buf.raw);
  assert_guard(&buf);
}

static void test_session_snapshot_stays_inside_old_buffer(void **state) {
  size_t n = 0;
  unsigned char *params;
  CPMDCSession *session;
  GuardedSnapshot buf;
  (void)state;
  params = read_file(g_params, &n);
  assert_non_null(params);
  session = cpmdc_session_create(params, n);
  free(params);
  assert_non_null(session);
  fill_guard(&buf);
  (void)cpmdc_session_last_property_snapshot(session,
                                            (CPMDCPropertySnapshot *)buf.raw);
  assert_guard(&buf);
  fill_guard(&buf);
  (void)cpmdc_last_property_snapshot((CPMDCPropertySnapshot *)buf.raw);
  assert_guard(&buf);
  cpmdc_session_destroy(session);
}

int main(int argc, char **argv) {
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(test_process_snapshot_stays_inside_old_buffer),
      cmocka_unit_test(test_session_snapshot_stays_inside_old_buffer),
  };
  if (argc < 2)
    return 1;
  g_params = argv[1];
  return cmocka_run_group_tests(tests, NULL, NULL);
}
