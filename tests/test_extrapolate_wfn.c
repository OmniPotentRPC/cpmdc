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
 *
 * With a second golden file and "keys", the session alternates the two
 * trajectories under orbital keys 1 and 2 with CPMDC_EXTRAP_KEYS=1, and
 * each trajectory must meet the same bounds against its own cpmd.x run.
 * With "default" instead, CPMDC_EXTRAP_KEYS stays unset, and no call may
 * start from a prediction, since every change of key empties the ring.
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

typedef struct {
  double pos[MAXF][3 * NAT], eks[MAXF], e[MAXF];
  int nf, scf_aspc, scf_none, used[MAXF], steps[MAXF];
} Traj;

static int load(const char *path, Traj *t) {
  FILE *g = fopen(path, "r");
  if (!g)
    return -1;
  char line[1024];
  while (fgets(line, sizeof line, g) && t->nf < MAXF) {
    double *x = t->pos[t->nf];
    if (sscanf(line, "scf_aspc %d", &t->scf_aspc) == 1 ||
        sscanf(line, "scf_none %d", &t->scf_none) == 1)
      continue;
    if (sscanf(line, "step %*d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
               &t->eks[t->nf], &x[0], &x[1], &x[2], &x[3], &x[4], &x[5],
               &x[6], &x[7], &x[8]) == 10)
      t->nf++;
  }
  fclose(g);
  return t->nf < 6 || t->scf_aspc <= 0 || t->scf_none <= 0 ? -1 : 0;
}

int main(int argc, char **argv) {
  if (argc != 3 && argc != 5) {
    fprintf(stderr, "usage: %s PARAMS_BIN GOLDEN_TXT [GOLDEN_TXT keys|default]\n",
            argv[0]);
    return 2;
  }
  if (!cpmdc_available()) {
    fprintf(stderr, "SKIP embed not available\n");
    return 77;
  }
  FILE *pf = fopen(argv[1], "rb");
  if (!pf)
    return 2;
  static unsigned char params[1 << 16];
  size_t ps = fread(params, 1, sizeof params, pf);
  fclose(pf);

  static Traj tr[2];
  int nt = argc == 5 ? 2 : 1;
  int keys = nt == 2 && strcmp(argv[4], "keys") == 0;
  if (nt == 2 && !keys && strcmp(argv[4], "default") != 0)
    return 2;
  for (int t = 0; t < nt; t++)
    if (load(argv[2 + t], &tr[t]) != 0)
      return 2;
  int nf = nt == 2 && tr[1].nf < tr[0].nf ? tr[1].nf : tr[0].nf;

  setenv("CPMDC_SCF_STEPS", "1", 1);
  if (keys)
    setenv("CPMDC_EXTRAP_KEYS", "1", 1);
  else
    unsetenv("CPMDC_EXTRAP_KEYS");
  CPMDCSession *s = cpmdc_session_create(params, ps);
  if (!s)
    return 1;
  fflush(stderr);
  int saved = dup(2);
  FILE *cap = tmpfile();
  dup2(fileno(cap), 2);
  const int z[NAT] = {8, 1, 1};
  double f[3 * NAT];
  int ok = 1;
  for (int k = 0; k < nf && ok; k++)
    for (int t = 0; t < nt && ok; t++) {
      double ang[3 * NAT];
      for (int i = 0; i < 3 * NAT; i++)
        ang[i] = tr[t].pos[k][i] * ANGSTROM_PER_BOHR;
      if (nt == 2)
        cpmdc_session_select_orbitals(s, t + 1);
      CPMDCResult r = cpmdc_session_energy_forces(s, NAT, ang, z, f);
      ok = r.ok;
      tr[t].e[k] = r.energy_h;
    }
  fflush(stderr);
  dup2(saved, 2);
  close(saved);
  cpmdc_session_destroy(s);
  if (!ok) {
    fprintf(stderr, "evaluation failed: %s\n", cpmdc_last_error());
    return 1;
  }

  int call = -1, v;
  char line[1024];
  rewind(cap);
  while (fgets(line, sizeof line, cap)) {
    Traj *t = call >= 0 ? &tr[call % nt] : NULL;
    if (strncmp(line, "cpmdc_warm_start ", 17) == 0 && call + 1 < nt * MAXF)
      call++;
    else if (t && sscanf(line, "cpmdc_wfn_extrapolated %d", &v) == 1)
      t->used[call / nt] = v;
    else if (t && sscanf(line, "cpmdc_scf_steps %d", &v) == 1)
      t->steps[call / nt] += v;
  }
  fclose(cap);
  if (call != nt * nf - 1) {
    fprintf(stderr, "%d calls reported, expected %d\n", call + 1, nt * nf);
    return 1;
  }

  for (int t = 0; t < nt; t++) {
    Traj *x = &tr[t];
    int total = 0;
    double worst = 0.0;
    for (int k = 0; k < nf; k++) {
      int want = nt == 2 && !keys ? 0 : (k >= 2 ? (k < 4 ? k : 4) : 0);
      if (k > 0) {
        total += x->steps[k];
        worst = fmax(worst, fabs(x->e[k] - x->eks[k]));
      }
      printf("key %d step %d energy_h=%.10f cpmd.x=%.10f scf=%d extrapolated=%d\n",
             t + 1, k, x->e[k], x->eks[k], x->steps[k], x->used[k]);
      if (x->used[k] != want)
        ok = 0;
    }
    printf("key %d max |dE|=%.3e scf %d, cpmd.x ASPC %d, cpmd.x none %d\n", t + 1,
           worst, total, x->scf_aspc, x->scf_none);
    ok = ok && worst < 1e-8;
    if (nt == 1 || keys)
      ok = ok && abs(total - x->scf_aspc) <= 0.1 * x->scf_aspc &&
           total < x->scf_none;
  }
  printf("PASS=%d\n", ok);
  return ok ? 0 : 1;
}
