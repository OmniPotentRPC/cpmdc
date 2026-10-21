/**
 * The caller's box is the cell, also against a deck that states another
 * one. The Si3N4 Geo1 deck says CELL 14 A (SYMMETRY 8); the ForceInput
 * box is 14.5 A. The embed must return cpmd.x's energy for a 14.5 A deck,
 * not the 14 A one.
 *
 * References: cpmd.x (OpenCPMD 062582b, 4 ranks) on the same deck with
 * CELL 14.5 and with CELL 14, orbitals converged to 1e-7.
 */
#include "cpmdc.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static const double CPMD_X_E_BOX145 = -51.31218184;
static const double CPMD_X_E_BOX14 = -51.31196661;

static unsigned char *slurp(const char *path, size_t *n) {
  FILE *f = fopen(path, "rb");
  if (!f)
    return NULL;
  fseek(f, 0, SEEK_END);
  long len = ftell(f);
  rewind(f);
  unsigned char *b = malloc((size_t)len);
  *n = fread(b, 1, (size_t)len, f);
  fclose(f);
  return b;
}

int main(int argc, char **argv) {
  if (argc < 3) {
    fprintf(stderr, "usage: %s PARAMS_BIN FORCE_INPUT_BIN\n", argv[0]);
    return 2;
  }
  if (!cpmdc_available()) {
    fprintf(stderr, "SKIP embed not available\n");
    return 77;
  }
  size_t ps = 0, fs = 0;
  unsigned char *p = slurp(argv[1], &ps);
  unsigned char *fi = slurp(argv[2], &fs);
  if (!p || !fi)
    return 2;
  CPMDCSession *s = cpmdc_session_create(p, ps);
  if (!s)
    return 1;
  double f[21];
  CPMDCResult r = cpmdc_session_calculate_forces(s, fi, fs, f, 21);
  cpmdc_session_destroy(s);
  if (!r.ok) {
    fprintf(stderr, "calculate_forces failed: %s\n", r.message);
    return 1;
  }
  printf("energy_h=%.10f cpmd.x 14.5 A=%.8f delta=%.3e cpmd.x 14 A=%.8f\n",
         r.energy_h, CPMD_X_E_BOX145, r.energy_h - CPMD_X_E_BOX145,
         CPMD_X_E_BOX14);
  int ok = fabs(r.energy_h - CPMD_X_E_BOX145) < 1e-6 &&
           fabs(r.energy_h - CPMD_X_E_BOX14) > 1e-4;
  printf("PASS=%d\n", ok);
  free(p);
  free(fi);
  return ok ? 0 : 1;
}
