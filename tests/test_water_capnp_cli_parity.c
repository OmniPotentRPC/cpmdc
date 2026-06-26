/**
 * Water BLYP via Cap'n Proto C ABI only (no caller INPUT). CLI reference energy.
 */
#include "cpmdc.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static const double CLI_WATER_E = -17.04926427;

static unsigned char *slurp(const char *path, size_t *n) {
  FILE *f = fopen(path, "rb");
  if (!f)
    return NULL;
  fseek(f, 0, SEEK_END);
  long L = ftell(f);
  rewind(f);
  unsigned char *b = malloc((size_t)L);
  *n = fread(b, 1, (size_t)L, f);
  fclose(f);
  return b;
}

static int run_once(const unsigned char *p, size_t ps, const unsigned char *fi,
                    size_t fs, double *e_out, double *f_out) {
  size_t need = cpmdc_potential_result_size_for_force_input(fi, fs);
  unsigned char *out = malloc(need);
  size_t outsz = 0;
  CPMDCResult r =
      cpmdc_calculate_result(p, ps, fi, fs, out, need, &outsz);
  if (!r.ok) {
    fprintf(stderr, "calculate_result failed: %s\n", r.message);
    free(out);
    return -1;
  }
  *e_out = r.energy_h;
  CPMDCSession *s = cpmdc_session_create(p, ps);
  if (!s) {
    free(out);
    return -1;
  }
  r = cpmdc_session_calculate_forces(s, fi, fs, f_out, 9);
  cpmdc_session_destroy(s);
  free(out);
  if (!r.ok) {
    fprintf(stderr, "session forces failed: %s\n", r.message);
    return -1;
  }
  return 0;
}

int main(int argc, char **argv) {
  const char *pb = argc > 1 ? argv[1] : "water_params.bin";
  const char *fb = argc > 2 ? argv[2] : "water_force.bin";
  if (!cpmdc_available()) {
    fprintf(stderr, "SKIP embed not available\n");
    return 77;
  }
  size_t ps = 0, fs = 0;
  unsigned char *p = slurp(pb, &ps);
  unsigned char *fi = slurp(fb, &fs);
  if (!p || !fi) {
    fprintf(stderr, "missing bins %s %s\n", pb, fb);
    return 2;
  }
  double e1 = 0, e2 = 0, f1[9] = {0}, f2[9] = {0};
  if (run_once(p, ps, fi, fs, &e1, f1) != 0)
    return 1;
  if (run_once(p, ps, fi, fs, &e2, f2) != 0)
    return 1;
  printf("run1 energy_h=%.12f\n", e1);
  printf("run2 energy_h=%.12f\n", e2);
  printf("cli_ref=%.12f delta1=%.3e delta2=%.3e\n", CLI_WATER_E, e1 - CLI_WATER_E,
         e2 - CLI_WATER_E);
  double maxf = 0;
  for (int i = 0; i < 9; ++i) {
    printf("F1[%d]=%.8e F2[%d]=%.8e\n", i, f1[i], i, f2[i]);
    if (fabs(f1[i]) > maxf)
      maxf = fabs(f1[i]);
  }
  printf("max|F|=%.8e\n", maxf);
  int ok = 1;
  if (fabs(e1 - CLI_WATER_E) > 1e-6 || fabs(e2 - CLI_WATER_E) > 1e-6)
    ok = 0;
  if (maxf < 1e-4)
    ok = 0; /* expect O z ~0.04 class */
  if (fabs(e1 - e2) > 1e-6)
    ok = 0;
  printf("PASS=%d\n", ok);
  free(p);
  free(fi);
  cpmdc_finalize();
  return ok ? 0 : 1;
}
