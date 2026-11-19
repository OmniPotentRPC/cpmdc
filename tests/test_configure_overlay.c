/**
 * PotentialConfig common-overlay lowering: session creation, evaluation, and
 * loud rejection of unsound combinations.
 */
#define _GNU_SOURCE
#include "cpmdc.h"

#include <setjmp.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <cmocka.h>

static const char *g_overlay = NULL;
static const char *g_reject = NULL;
static const char *g_step = NULL;
static const char *g_bad_params = NULL;
static const char *g_arm_overlay = NULL;

static unsigned char *read_file(const char *path, size_t *size) {
  FILE *fp = fopen(path, "rb");
  if (!fp)
    return NULL;
  if (fseek(fp, 0, SEEK_END) != 0) {
    fclose(fp);
    return NULL;
  }
  long n = ftell(fp);
  if (n <= 0 || n > (16 << 20)) {
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

static void test_overlay_session_evaluates(void **state) {
  (void)state;
  size_t config_size = 0, step_size = 0;
  unsigned char *config = read_file(g_overlay, &config_size);
  unsigned char *step = read_file(g_step, &step_size);
  assert_non_null(config);
  assert_non_null(step);

  CPMDCSession *session =
      cpmdc_session_create_from_config(config, config_size);
  assert_non_null(session);

  size_t need = cpmdc_potential_result_size_for_force_input(step, step_size);
  assert_true(need > 0);
  unsigned char *out = (unsigned char *)malloc(need);
  assert_non_null(out);
  size_t wrote = 0;
  CPMDCResult r =
      cpmdc_session_calculate_result(session, step, step_size, out, need,
                                     &wrote);
  assert_int_equal(r.ok, 1);
  assert_true(wrote > 0);

  size_t wrote_one_shot = 0;
  CPMDCResult one_shot = cpmdc_calculate_result_from_config(
      config, config_size, step, step_size, out, need, &wrote_one_shot);
  assert_int_equal(one_shot.ok, 1);
  assert_int_equal(wrote_one_shot, wrote);

  cpmdc_session_destroy(session);
  free(out);
  free(step);
  free(config);
}

static void test_overlay_rejects_unlowerable_field(void **state) {
  (void)state;
  size_t config_size = 0;
  unsigned char *config = read_file(g_reject, &config_size);
  assert_non_null(config);
  assert_null(cpmdc_session_create_from_config(config, config_size));
  assert_non_null(strstr(cpmdc_last_error(), "relativityMethod"));
  assert_int_not_equal(cpmdc_configure(config, config_size), 0);
  free(config);
}

static void test_overlay_kmesh_reaches_deck(void **state) {
  (void)state;
  size_t config_size = 0;
  unsigned char *config = read_file(g_overlay, &config_size);
  assert_non_null(config);
  char deck_path[] = "overlay_kmesh_deck_XXXXXX";
  int fd = mkstemp(deck_path);
  assert_true(fd >= 0);
  close(fd);
  setenv("CPMDC_DECK_OUT", deck_path, 1);
  CPMDCSession *session =
      cpmdc_session_create_from_config(config, config_size);
  unsetenv("CPMDC_DECK_OUT");
  assert_non_null(session);
  size_t deck_size = 0;
  unsigned char *deck = read_file(deck_path, &deck_size);
  remove(deck_path);
  assert_non_null(deck);
  unsigned char *text = (unsigned char *)realloc(deck, deck_size + 1);
  assert_non_null(text);
  text[deck_size] = '\0';
  /* potential_config_common_cpmd.capnp.txt sets kMesh = [2, 2, 2]. */
  assert_non_null(strstr((const char *)text, "KPOINTS MONKHORST-PACK\n  2 2 2"));
  free(text);
  cpmdc_session_destroy(session);
  free(config);
}

/* Renders the deck of a constant-cutoff call on the session into a text
 * file, which shows the overlay the session currently renders with. */
static char *deck_after_constant_cutoff(CPMDCSession *session) {
  char deck_path[] = "overlay_atomic_deck_XXXXXX";
  int fd = mkstemp(deck_path);
  assert_true(fd >= 0);
  close(fd);
  setenv("CPMDC_DECK_OUT", deck_path, 1);
  int rc = cpmdc_session_set_constant_cutoff(session, 10.0, 0.0, 0.0);
  unsetenv("CPMDC_DECK_OUT");
  assert_int_equal(rc, 0);
  size_t size = 0;
  unsigned char *deck = read_file(deck_path, &size);
  remove(deck_path);
  assert_non_null(deck);
  char *text = (char *)realloc(deck, size + 1);
  assert_non_null(text);
  text[size] = '\0';
  return text;
}

static void test_failed_set_params_keeps_overlay(void **state) {
  (void)state;
  size_t config_size = 0, bad_size = 0;
  unsigned char *config = read_file(g_arm_overlay, &config_size);
  unsigned char *bad = read_file(g_bad_params, &bad_size);
  assert_non_null(config);
  assert_non_null(bad);

  CPMDCSession *control =
      cpmdc_session_create_from_config(config, config_size);
  assert_non_null(control);
  char *expected = deck_after_constant_cutoff(control);
  assert_non_null(strstr(expected, "FUNCTIONAL PBE"));
  assert_non_null(strstr(expected, "CHARGE\n  -1\n"));
  assert_non_null(strstr(expected, "MULTIPLICITY\n  3\n"));

  CPMDCSession *session =
      cpmdc_session_create_from_config(config, config_size);
  assert_non_null(session);

  const unsigned char garbage[] = {0xde, 0xad, 0xbe, 0xef, 0x01, 0x02,
                                   0x03, 0x04, 0x05, 0x06, 0x07, 0x08};
  assert_int_equal(cpmdc_session_set_params(session, garbage, sizeof(garbage)),
                   -1);
  assert_string_equal(cpmdc_last_error(), "CPMDParams parse failed");

  assert_int_equal(cpmdc_session_set_params(session, bad, bad_size), -1);
  assert_non_null(strstr(cpmdc_last_error(), "deck render failed"));

  char *actual = deck_after_constant_cutoff(session);
  assert_string_equal(actual, expected);

  free(actual);
  free(expected);
  cpmdc_session_destroy(session);
  cpmdc_session_destroy(control);
  free(bad);
  free(config);
}

int main(int argc, char **argv) {
  if (argc != 6) {
    fprintf(stderr,
            "usage: %s OVERLAY_BIN REJECT_BIN STEP_BIN BAD_PARAMS_BIN "
            "ARM_OVERLAY_BIN\n",
            argv[0]);
    return 2;
  }
  g_overlay = argv[1];
  g_reject = argv[2];
  g_step = argv[3];
  g_bad_params = argv[4];
  g_arm_overlay = argv[5];
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(test_overlay_session_evaluates),
      cmocka_unit_test(test_overlay_rejects_unlowerable_field),
      cmocka_unit_test(test_overlay_kmesh_reaches_deck),
      cmocka_unit_test(test_failed_set_params_keeps_overlay),
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
