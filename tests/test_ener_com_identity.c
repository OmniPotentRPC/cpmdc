/**
 * OpenCPMD energy-table identity for water BLYP.
 *
 * wrener in src/wrener_utils.mod.F90 at OpenCPMD 062582b7cfd8 builds
 * ENERGY SUM(1) at lines 595-615 and prints it at lines 688-689 with F20.8.
 * Water BLYP is not SIC, not finite electronic temperature, not Hubbard,
 * and not the diagonalisation total, so
 *   elec1 = ehep + esr - eself
 *   etot1 = ekin + elec1 + epseu + enl + exc
 * Half a unit in the last printed place is 5e-9 Ha.
 *
 * Terms in ener_com that this sum does not use: eht (the stored
 * electrostatic energy; the sum uses elec1), ehee, ehii, vxc, egc, eeig,
 * eband, entropy, ecnstr, amu, ebogo, eext, etddft, ehsic, erestr, eefield.
 * Terms OpenCPMD prints from other modules, not from ener_com: hubbu%ehub,
 * vdwr%evdw, the CDFT constraint energy, and the QM/MM energies. Water
 * BLYP sets none of those.
 *
 * The printout is report-only across ranks: reductions live in OpenCPMD.
 */
#include "cpmdc.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static const char *g_params = NULL;
static const char *g_step = NULL;

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

static uint64_t fnv_bytes(const void *data, size_t n) {
  const unsigned char *p = (const unsigned char *)data;
  uint64_t h = 1469598103934665603ull;
  for (size_t i = 0; i < n; ++i) {
    h ^= p[i];
    h *= 1099511628211ull;
  }
  return h;
}

int main(int argc, char **argv) {
  if (argc < 3) {
    fprintf(stderr, "usage: %s params.bin force_input.bin\n", argv[0]);
    return 2;
  }
  g_params = argv[1];
  g_step = argv[2];
  size_t params_size = 0;
  size_t step_size = 0;
  unsigned char *params = read_file(g_params, &params_size);
  unsigned char *step = read_file(g_step, &step_size);
  if (!params || !step) {
    fprintf(stderr, "failed to read inputs\n");
    return 1;
  }
  CPMDCSession *session = cpmdc_session_create(params, params_size);
  if (!session) {
    fprintf(stderr, "session create failed: %s\n", cpmdc_last_error());
    return 1;
  }
  size_t need = cpmdc_potential_result_size_for_force_input(step, step_size);
  unsigned char *result_buf = (unsigned char *)malloc(need);
  if (!result_buf)
    return 1;
  size_t wrote = 0;
  CPMDCResult r = cpmdc_session_calculate_result(session, step, step_size,
                                                 result_buf, need, &wrote);
  if (!r.ok) {
    fprintf(stderr, "calculate failed: %s\n", r.message);
    return 1;
  }
  CPMDCEnergyComponents c;
  memset(&c, 0, sizeof(c));
  if (cpmdc_last_energy_components(&c) != 0 || !c.valid) {
    fprintf(stderr, "energy components missing\n");
    return 1;
  }
  const double *slot = &c.etot;
  static const char *names[24] = {
      "etot",    "ekin",   "epseu", "enl",    "eht",    "ehep",
      "ehee",    "ehii",   "exc",   "vxc",    "egc",    "esr",
      "eeig",    "eband",  "entropy", "eself", "ecnstr", "amu",
      "ebogo",   "eext",   "etddft", "ehsic",  "erestr", "eefield"};
  double max_abs = 0.0;
  for (int i = 0; i < 24; ++i) {
    if (!isfinite(slot[i])) {
      fprintf(stderr, "%s is not finite\n", names[i]);
      return 1;
    }
    if (fabs(slot[i]) > max_abs)
      max_abs = fabs(slot[i]);
    printf("ener_com %s %.16e\n", names[i], slot[i]);
  }
  printf("ener_com_fnv %016llx\n",
         (unsigned long long)fnv_bytes(slot, 24 * sizeof(double)));
  /* wrener lines 600-615, non-SIC, non-tfint, non-Hubbard grouping. */
  double elec1 = c.ehep + c.esr - c.eself;
  double etot1 = c.ekin + elec1 + c.epseu + c.enl + c.exc;
  double diff = etot1 - c.etot;
  printf("identity etot=%.12f sum1=%.12f diff=%.3e tol=5e-9\n", c.etot, etot1,
         diff);
  printf("identity_omitted eht ehee ehii vxc egc eeig eband entropy ecnstr "
         "amu ebogo eext etddft ehsic erestr eefield hubbu ehvdw cdft qmmm\n");
  if (fabs(diff) >= 5e-9) {
    fprintf(stderr, "ENERGY SUM(1) differs from etot by %.3e Ha\n", diff);
    return 1;
  }
  if (!(max_abs > 1e-8)) {
    fprintf(stderr, "components are all ~0\n");
    return 1;
  }
  cpmdc_session_destroy(session);
  free(result_buf);
  free(params);
  free(step);
  cpmdc_finalize();
  return 0;
}
