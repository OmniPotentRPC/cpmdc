/**
 * EXTRAPOLATE WFN on the in-process route, against cpmd.x BO-MD.
 *
 * tests/golden/water_aspc_bomd.txt is a cpmd.x run of the deck of
 * water_aspc_params.capnp.txt as MOLECULAR DYNAMICS BO (TIMESTEP 10,
 * EXTRAPOLATE WFN ASPC 4): the start and the positions and EKS of each
 * step, and the SCF iterations cpmd.x took over the steps with and
 * without the extrapolation. One session evaluates the same positions
 * in order. Each energy must match EKS, call k from the third on must
 * start from the prediction over min(k - 1, 4) earlier orbitals, and the
 * SCF iterations must come within 10 percent of cpmd.x with ASPC and
 * below cpmd.x without it.
 */
#include "cpmdc.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define NAT 3
#define MAXF 64

static const double ANGSTROM_PER_BOHR = 0.529177210859;

int main(int argc, char **argv) {
  if (argc < 3) {
    fprintf(stderr, "usage: %s PARAMS_BIN GOLDEN_TXT\n", argv[0]);
    return 2;
  }
  if (!cpmdc_available()) {
    fprintf(stderr, "SKIP embed not available\n");
    return 77;
  }
  FILE *g = fopen(argv[2], "r");
  FILE *pf = fopen(argv[1], "rb");
  if (!g || !pf)
    return 2;
  static unsigned char params[1 << 16];
  size_t ps = fread(params, 1, sizeof params, pf);
  fclose(pf);

  double pos[MAXF][3 * NAT], eks[MAXF];
  int nf = 0, scf_aspc = 0, scf_none = 0;
  char line[1024];
  while (fgets(line, sizeof line, g) && nf < MAXF) {
    double *x = pos[nf];
    if (sscanf(line, "scf_aspc %d", &scf_aspc) == 1 ||
        sscanf(line, "scf_none %d", &scf_none) == 1)
      continue;
    if (sscanf(line, "step %*d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
               &eks[nf], &x[0], &x[1], &x[2], &x[3], &x[4], &x[5], &x[6],
               &x[7], &x[8]) == 10)
      nf++;
  }
  fclose(g);
  if (nf < 6 || scf_aspc <= 0 || scf_none <= 0)
    return 2;

  setenv("CPMDC_SCF_STEPS", "1", 1);
  CPMDCSession *s = cpmdc_session_create(params, ps);
  if (!s)
    return 1;
  fflush(stderr);
  int saved = dup(2);
  FILE *cap = tmpfile();
  dup2(fileno(cap), 2);
  const int z[NAT] = {8, 1, 1};
  double e[MAXF], f[3 * NAT];
  int ok = 1;
  for (int k = 0; k < nf && ok; k++) {
    double ang[3 * NAT];
    for (int i = 0; i < 3 * NAT; i++)
      ang[i] = pos[k][i] * ANGSTROM_PER_BOHR;
    CPMDCResult r = cpmdc_session_energy_forces(s, NAT, ang, z, f);
    ok = r.ok;
    e[k] = r.energy_h;
  }
  fflush(stderr);
  dup2(saved, 2);
  close(saved);
  cpmdc_session_destroy(s);
  if (!ok) {
    fprintf(stderr, "evaluation failed: %s\n", cpmdc_last_error());
    return 1;
  }

  int call = -1, used[MAXF] = {0}, steps[MAXF] = {0}, v;
  rewind(cap);
  while (fgets(line, sizeof line, cap)) {
    if (strncmp(line, "cpmdc_warm_start ", 17) == 0 && call + 1 < MAXF)
      call++;
    else if (call >= 0 && sscanf(line, "cpmdc_wfn_extrapolated %d", &v) == 1)
      used[call] = v;
    else if (call >= 0 && sscanf(line, "cpmdc_scf_steps %d", &v) == 1)
      steps[call] += v;
  }
  fclose(cap);
  if (call != nf - 1) {
    fprintf(stderr, "%d calls reported, expected %d\n", call + 1, nf);
    return 1;
  }

  int total = 0;
  double worst = 0.0;
  for (int k = 0; k < nf; k++) {
    int want = k >= 2 ? (k < 4 ? k : 4) : 0;
    if (k > 0) {
      total += steps[k];
      worst = fmax(worst, fabs(e[k] - eks[k]));
    }
    printf("step %d energy_h=%.10f cpmd.x=%.10f scf=%d extrapolated=%d\n", k,
           e[k], eks[k], steps[k], used[k]);
    if (used[k] != want)
      ok = 0;
  }
  printf("max |dE|=%.3e scf %d, cpmd.x ASPC %d, cpmd.x none %d\n", worst,
         total, scf_aspc, scf_none);
  ok = ok && worst < 1e-8 && fabs(total - scf_aspc) <= 0.1 * scf_aspc &&
       total < scf_none;
  printf("PASS=%d\n", ok);
  return ok ? 0 : 1;
}
