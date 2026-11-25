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
#define _POSIX_C_SOURCE 200809L

#include "cpmdc.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define NAT 3
#define MAXF 64

static const double ANGSTROM_PER_BOHR = 0.529177210859;
static const double KEY_OFFSET_ANGSTROM = 0.001;

int main(int argc, char **argv) {
  if (argc < 4) {
    fprintf(stderr, "usage: %s PARAMS_BIN GOLDEN_TXT KEYED_GOLDEN_TXT\n", argv[0]);
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
  ok = ok && worst < 1e-8 && abs(total - scf_aspc) <= 0.1 * scf_aspc &&
       total < scf_none;

  FILE *kg = fopen(argv[3], "r");
  if (!kg)
    return 2;
  int keyed_want[2 * MAXF], keyed_want_n = 0;
  while (keyed_want_n < 2 * MAXF &&
         fscanf(kg, " extrapolated %d", &keyed_want[keyed_want_n]) == 1)
    keyed_want_n++;
  fclose(kg);
  if (keyed_want_n != 2 * nf)
    return 2;

  s = cpmdc_session_create(params, ps);
  if (!s)
    return 1;
  saved = dup(2);
  cap = tmpfile();
  dup2(fileno(cap), 2);
  double keyed_e[2 * MAXF];
  for (int k = 0; k < nf && ok; k++) {
    double ang[3 * NAT];
    for (int i = 0; i < 3 * NAT; i++)
      ang[i] = pos[k][i] * ANGSTROM_PER_BOHR;
    for (int bead = 0; bead < 2; bead++) {
      CPMDCResult r;
      double bead_ang[3 * NAT];
      memcpy(bead_ang, ang, sizeof(bead_ang));
      for (int atom = 0; atom < NAT; atom++)
        bead_ang[3 * atom] += bead * KEY_OFFSET_ANGSTROM;
      cpmdc_session_select_orbitals(s, 101 + bead);
      r = cpmdc_session_energy_forces(s, NAT, bead_ang, z, f);
      ok = r.ok;
      keyed_e[2 * k + bead] = r.energy_h;
      if (!ok)
        break;
    }
  }
  fflush(stderr);
  dup2(saved, 2);
  close(saved);
  cpmdc_session_destroy(s);
  if (!ok) {
    fprintf(stderr, "keyed evaluation failed: %s\n", cpmdc_last_error());
    return 1;
  }

  int keyed_call = -1, keyed_used[2 * MAXF] = {0};
  int keyed_steps[2 * MAXF] = {0};
  rewind(cap);
  while (fgets(line, sizeof line, cap)) {
    if (strncmp(line, "cpmdc_warm_start ", 17) == 0 &&
        keyed_call + 1 < 2 * MAXF)
      keyed_call++;
    else if (keyed_call >= 0 &&
             sscanf(line, "cpmdc_wfn_extrapolated %d", &v) == 1)
      keyed_used[keyed_call] = v;
    else if (keyed_call >= 0 && sscanf(line, "cpmdc_scf_steps %d", &v) == 1)
      keyed_steps[keyed_call] += v;
  }
  fclose(cap);
  if (keyed_call != 2 * nf - 1) {
    fprintf(stderr, "%d keyed calls reported, expected %d\n", keyed_call + 1,
            2 * nf);
    return 1;
  }
  for (int k = 0; k < nf; k++) {
    for (int bead = 0; bead < 2; bead++) {
      int i = 2 * k + bead;
      printf("key %d step %d energy_h=%.10f scf=%d extrapolated=%d\n",
             101 + bead, k, keyed_e[i], keyed_steps[i], keyed_used[i]);
      if (keyed_used[i] != keyed_want[i] ||
          fabs(keyed_e[i] - eks[k]) >= 1e-8)
        ok = 0;
    }
  }

  setenv("CPMDC_ORBITAL_SLOT_LIMIT", "9", 1);
  s = cpmdc_session_create(params, ps);
  if (!s)
    return 1;
  double last_ang[3 * NAT];
  for (int i = 0; i < 3 * NAT; i++)
    last_ang[i] = pos[nf - 1][i] * ANGSTROM_PER_BOHR;
  cpmdc_session_select_orbitals(s, 303);
  CPMDCResult limited =
      cpmdc_session_energy_forces(s, NAT, last_ang, z, f);
  if (limited.ok || strstr(limited.message, "CPMDC_ORBITAL_SLOT_LIMIT") == NULL ||
      strstr(limited.message, "unchanged") == NULL)
    ok = 0;
  setenv("CPMDC_ORBITAL_SLOT_LIMIT", "10", 1);
  limited = cpmdc_session_energy_forces(s, NAT, last_ang, z, f);
  if (!limited.ok)
    ok = 0;
  cpmdc_session_destroy(s);
  unsetenv("CPMDC_ORBITAL_SLOT_LIMIT");
  printf("PASS=%d\n", ok);
  return ok ? 0 : 1;
}
