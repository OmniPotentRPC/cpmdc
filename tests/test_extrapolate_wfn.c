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
 * tests/golden/water_aspc_two_key.txt holds one line per keyed call of
 * the first KEYED_TEST_STEPS steps, two interleaved keys per step. A
 * key that an unkeyed call clears must not lend its history to that
 * call, and a destroyed session must not lend its history to the next.
 *
 * A build without the OpenCPMD embed checks the link stub and reference
 * session contract instead, fails when it breaks, and then skips the
 * cpmd.x comparison.
 */
#define _POSIX_C_SOURCE 200809L

#include "cpmdc.h"
#include "setup_child.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define NAT 3
#define MAXF 64
#define CLEARED_TEST_STEPS 6
#define CLEARED_KEYED_STEPS 4

static const double ANGSTROM_PER_BOHR = 0.529177210859;
static const double KEY_OFFSET_ANGSTROM = 0.05;
static const int KEYED_TEST_STEPS = 6;
static const long long CLEARED_KEY = 401;

/* Converged c0 in the selected extrapolation ring and the count the last
 * SCF combined. */
void cpmdc_embed_extrapolation_state(int *count, int *used);
void cpmdc_embed_clear_orbital_key(int drop_history);
int cpmdc_embed_copy_orbitals(long long from_key, long long to_key);

typedef struct {
  const unsigned char *params;
  size_t params_size;
  double prime_ang[3 * NAT];
  double first_ang[3 * NAT];
  double last_ang[3 * NAT];
} LimitedContext;

typedef struct {
  int ok;
  char message[256];
} LimitedOutput;

typedef struct {
  const unsigned char *params;
  size_t params_size;
  double ang[CLEARED_TEST_STEPS][3 * NAT];
} ClearedContext;

typedef struct {
  int ok;
  int keyed_used;
  int unkeyed_used;
  int unkeyed_ring;
  int destroyed_ring;
  char message[256];
} ClearedOutput;

/* Counts the calls in captured stderr and records, per call, the
 * extrapolation each reported and, when steps is not NULL, its SCF
 * iterations. */
static int read_capture(FILE *cap, int max_calls, int *used, int *steps) {
  char line[1024];
  int call = -1, v;
  rewind(cap);
  while (fgets(line, sizeof line, cap)) {
    if (strncmp(line, "cpmdc_warm_start ", 17) == 0 && call + 1 < max_calls)
      call++;
    else if (call >= 0 && sscanf(line, "cpmdc_wfn_extrapolated %d", &v) == 1)
      used[call] = v;
    else if (call >= 0 && steps &&
             sscanf(line, "cpmdc_scf_steps %d", &v) == 1)
      steps[call] += v;
  }
  return call + 1;
}

static int capture_stderr(FILE **cap, int *saved) {
  fflush(stderr);
  *cap = tmpfile();
  if (!*cap)
    return -1;
  *saved = dup(2);
  if (*saved < 0 || dup2(fileno(*cap), 2) < 0) {
    if (*saved >= 0)
      close(*saved);
    *saved = -1;
    fclose(*cap);
    *cap = NULL;
    return -1;
  }
  return 0;
}

static void restore_stderr(int *saved) {
  fflush(stderr);
  if (*saved < 0)
    return;
  dup2(*saved, 2);
  close(*saved);
  *saved = -1;
}

static void limited_session_child(void *opaque, void *result) {
  const LimitedContext *ctx = opaque;
  LimitedOutput *out = result;
  const int z[NAT] = {8, 1, 1};
  double f[3 * NAT];
  CPMDCSession *session = cpmdc_session_create(ctx->params, ctx->params_size);
  if (!session) {
    snprintf(out->message, sizeof(out->message), "create: %s",
             cpmdc_last_error());
    return;
  }
  if (cpmdc_session_set_constant_cutoff(session, 40.0, 0.0, 0.0) != 0) {
    snprintf(out->message, sizeof(out->message), "cutoff: %s",
             cpmdc_last_error());
    goto done;
  }
  CPMDCResult evaluated = cpmdc_session_energy_forces(
      session, NAT, ctx->prime_ang, z, f);
  if (!evaluated.ok) {
    snprintf(out->message, sizeof(out->message), "prime: %.240s",
             evaluated.message);
    goto done;
  }
  setenv("CPMDC_ORBITAL_SLOT_LIMIT", "1", 1);
  cpmdc_session_select_orbitals(session, 302);
  evaluated =
      cpmdc_session_energy_forces(session, NAT, ctx->first_ang, z, f);
  if (!evaluated.ok) {
    snprintf(out->message, sizeof(out->message), "key 302: %.238s",
             evaluated.message);
    goto done;
  }
  cpmdc_session_select_orbitals(session, 303);
  evaluated = cpmdc_session_energy_forces(session, NAT, ctx->last_ang, z, f);
  if (!evaluated.ok) {
    snprintf(out->message, sizeof(out->message), "key 303: %.238s",
             evaluated.message);
    goto done;
  }
  out->ok = 1;

done:
  cpmdc_session_destroy(session);
  unsetenv("CPMDC_ORBITAL_SLOT_LIMIT");
}

/* One keyed session fills its history along the trajectory, then a live
 * unkeyed session evaluates the next steps, clearing the key, and is
 * destroyed with a ring of its own. */
static void cleared_key_child(void *opaque, void *result) {
  const ClearedContext *ctx = opaque;
  ClearedOutput *out = result;
  const int z[NAT] = {8, 1, 1};
  double f[3 * NAT];
  int used[CLEARED_TEST_STEPS] = {0};
  int saved = -1, calls, last_used;
  FILE *cap = NULL;
  CPMDCSession *keyed = NULL, *unkeyed = NULL;
  CPMDCResult r;

  if (capture_stderr(&cap, &saved) != 0) {
    snprintf(out->message, sizeof(out->message), "stderr capture failed");
    return;
  }
  keyed = cpmdc_session_create(ctx->params, ctx->params_size);
  unkeyed = cpmdc_session_create(ctx->params, ctx->params_size);
  if (!keyed || !unkeyed) {
    snprintf(out->message, sizeof(out->message), "create: %.240s",
             cpmdc_last_error());
    goto done;
  }
  cpmdc_session_select_orbitals(keyed, CLEARED_KEY);
  for (int k = 0; k < CLEARED_TEST_STEPS; k++) {
    CPMDCSession *s = k < CLEARED_KEYED_STEPS ? keyed : unkeyed;
    r = cpmdc_session_energy_forces(s, NAT, ctx->ang[k], z, f);
    if (!r.ok) {
      snprintf(out->message, sizeof(out->message), "step %d: %.230s", k,
               r.message);
      goto done;
    }
  }
  cpmdc_embed_extrapolation_state(&out->unkeyed_ring, &last_used);
  cpmdc_session_destroy(unkeyed);
  unkeyed = NULL;
  cpmdc_embed_extrapolation_state(&out->destroyed_ring, &last_used);
  out->ok = 1;

done:
  restore_stderr(&saved);
  if (out->ok) {
    calls = read_capture(cap, CLEARED_TEST_STEPS, used, NULL);
    if (calls != CLEARED_TEST_STEPS) {
      snprintf(out->message, sizeof(out->message),
               "%d calls reported, expected %d", calls, CLEARED_TEST_STEPS);
      out->ok = 0;
    }
    out->keyed_used = used[CLEARED_KEYED_STEPS - 1];
    out->unkeyed_used = used[CLEARED_KEYED_STEPS];
  }
  fclose(cap);
  cpmdc_session_destroy(unkeyed);
  cpmdc_session_destroy(keyed);
}

/* Without OpenCPMD the ring entry points are the link stub: it reports an
 * empty ring and refuses every orbital copy, a copy of a key onto itself
 * included. The reference session has no orbital store, so a keyed
 * evaluation and a copy of its key onto itself both succeed and leave the
 * stub ring empty. */
static int reference_contract_holds(const unsigned char *params, size_t ps,
                                    const double *ang) {
  const int z[NAT] = {8, 1, 1};
  double f[3 * NAT];
  int count = -1, used = -1, ok;
  CPMDCSession *s;
  CPMDCResult r;

  cpmdc_embed_clear_orbital_key(1);
  cpmdc_embed_extrapolation_state(&count, &used);
  ok = count == 0 && used == 0 &&
       cpmdc_embed_copy_orbitals(CLEARED_KEY, CLEARED_KEY) == -1 &&
       cpmdc_embed_copy_orbitals(CLEARED_KEY, CLEARED_KEY + 1) == -1;
  printf("stub ring=%d used=%d ok=%d\n", count, used, ok);

  s = cpmdc_session_create(params, ps);
  if (!s)
    return 0;
  cpmdc_session_select_orbitals(s, CLEARED_KEY);
  r = cpmdc_session_energy_forces(s, NAT, ang, z, f);
  count = used = -1;
  cpmdc_embed_extrapolation_state(&count, &used);
  ok = ok && r.ok && isfinite(r.energy_h) &&
       cpmdc_session_copy_orbitals(s, CLEARED_KEY, CLEARED_KEY) == 0 &&
       count == 0 && used == 0;
  cpmdc_session_destroy(s);
  printf("reference keyed=%d self_copy ring=%d ok=%d\n", r.ok, count, ok);
  return ok;
}

int main(int argc, char **argv) {
  if (argc < 4) {
    fprintf(stderr, "usage: %s PARAMS_BIN GOLDEN_TXT KEYED_GOLDEN_TXT\n", argv[0]);
    return 2;
  }
  FILE *g = fopen(argv[2], "r");
  FILE *pf = fopen(argv[1], "rb");
  if (!g || !pf) {
    if (g)
      fclose(g);
    if (pf)
      fclose(pf);
    return 2;
  }
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
  if (!cpmdc_available()) {
    double ang[3 * NAT];
    for (int i = 0; i < 3 * NAT; i++)
      ang[i] = pos[0][i] * ANGSTROM_PER_BOHR;
    if (!reference_contract_holds(params, ps, ang)) {
      fprintf(stderr, "reference build broke the link stub contract: %s\n",
              cpmdc_last_error());
      return 1;
    }
    fprintf(stderr, "SKIP cpmd.x parity needs the OpenCPMD embed; "
                    "link stub contract checked\n");
    return 77;
  }

  int keyed_nf = nf < KEYED_TEST_STEPS ? nf : KEYED_TEST_STEPS;
  int keyed_want[2 * MAXF], keyed_want_n = 0;
  FILE *kg = fopen(argv[3], "r");
  if (!kg)
    return 2;
  while (fgets(line, sizeof line, kg)) {
    if (keyed_want_n == 2 * MAXF ||
        sscanf(line, "extrapolated %d", &keyed_want[keyed_want_n]) != 1) {
      keyed_want_n = -1;
      break;
    }
    keyed_want_n++;
  }
  fclose(kg);
  if (keyed_want_n != 2 * keyed_nf) {
    fprintf(stderr, "%s: %d lines, expected one per keyed call, %d\n",
            argv[3], keyed_want_n, 2 * keyed_nf);
    return 2;
  }

  setenv("CPMDC_SCF_STEPS", "1", 1);
  LimitedContext limited_ctx = {.params = params, .params_size = ps};
  for (int i = 0; i < 3 * NAT; i++) {
    limited_ctx.prime_ang[i] = pos[nf - 3][i] * ANGSTROM_PER_BOHR;
    limited_ctx.first_ang[i] = pos[nf - 2][i] * ANGSTROM_PER_BOHR;
    limited_ctx.last_ang[i] = pos[nf - 1][i] * ANGSTROM_PER_BOHR;
  }
  LimitedOutput isolated = {0};
  char isolated_why[256];
  int ok =
      run_setup_child(limited_session_child, &limited_ctx, &isolated,
                      sizeof(isolated), isolated_why, sizeof(isolated_why)) == 0;
  if (ok && !isolated.ok) {
    fprintf(stderr, "isolated limited_session failed: %s\n", isolated.message);
    ok = 0;
  } else if (!ok) {
    fprintf(stderr, "isolated limited_session child failed: %s\n",
            isolated_why);
  }
  printf("limited_session isolated=%d\n", ok);

  ClearedContext cleared_ctx = {.params = params, .params_size = ps};
  for (int k = 0; k < CLEARED_TEST_STEPS; k++)
    for (int i = 0; i < 3 * NAT; i++)
      cleared_ctx.ang[k][i] = pos[k][i] * ANGSTROM_PER_BOHR;
  ClearedOutput cleared = {0};
  char cleared_why[256];
  if (run_setup_child(cleared_key_child, &cleared_ctx, &cleared,
                      sizeof(cleared), cleared_why,
                      sizeof(cleared_why)) != 0) {
    fprintf(stderr, "cleared_key child failed: %s\n", cleared_why);
    ok = 0;
  } else if (!cleared.ok) {
    fprintf(stderr, "cleared_key failed: %s\n", cleared.message);
    ok = 0;
  } else if (cleared.keyed_used != CLEARED_KEYED_STEPS - 1 ||
             cleared.unkeyed_used != 0 || cleared.unkeyed_ring != 1 ||
             cleared.destroyed_ring != 0) {
    ok = 0;
  }
  printf("cleared_key keyed=%d unkeyed=%d unkeyed_ring=%d destroyed_ring=%d "
         "ok=%d\n",
         cleared.keyed_used, cleared.unkeyed_used, cleared.unkeyed_ring,
         cleared.destroyed_ring, ok);

  const int z[NAT] = {8, 1, 1};
  double e[MAXF], f[3 * NAT], keyed_e[2 * MAXF];
  int used[MAXF] = {0}, steps[MAXF] = {0}, calls;
  int keyed_used[2 * MAXF] = {0}, keyed_steps[2 * MAXF] = {0};
  int saved = -1, rc = 1;
  FILE *cap = NULL;
  CPMDCSession *s = NULL;
  LimitedOutput after_full = {0};

  s = cpmdc_session_create(params, ps);
  if (!s)
    goto out;
  if (capture_stderr(&cap, &saved) != 0)
    goto out;
  for (int k = 0; k < nf && ok; k++) {
    double ang[3 * NAT];
    for (int i = 0; i < 3 * NAT; i++)
      ang[i] = pos[k][i] * ANGSTROM_PER_BOHR;
    CPMDCResult r = cpmdc_session_energy_forces(s, NAT, ang, z, f);
    ok = r.ok;
    e[k] = r.energy_h;
  }
  restore_stderr(&saved);
  cpmdc_session_destroy(s);
  s = NULL;
  if (!ok) {
    fprintf(stderr, "evaluation failed: %s\n", cpmdc_last_error());
    goto out;
  }
  calls = read_capture(cap, MAXF, used, steps);
  fclose(cap);
  cap = NULL;
  if (calls != nf) {
    fprintf(stderr, "%d calls reported, expected %d\n", calls, nf);
    goto out;
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

  s = cpmdc_session_create(params, ps);
  if (!s)
    goto out;
  if (capture_stderr(&cap, &saved) != 0)
    goto out;
  for (int k = 0; k < keyed_nf && ok; k++) {
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
  restore_stderr(&saved);
  if (!ok) {
    fprintf(stderr, "keyed evaluation failed: %s\n", cpmdc_last_error());
    goto out;
  }
  calls = read_capture(cap, 2 * MAXF, keyed_used, keyed_steps);
  fclose(cap);
  cap = NULL;
  if (calls != 2 * keyed_nf) {
    fprintf(stderr, "%d keyed calls reported, expected %d\n", calls,
            2 * keyed_nf);
    goto out;
  }
  for (int k = 0; k < keyed_nf; k++) {
    for (int bead = 0; bead < 2; bead++) {
      int i = 2 * k + bead;
      printf("key %d step %d energy_h=%.10f scf=%d extrapolated=%d\n",
             101 + bead, k, keyed_e[i], keyed_steps[i], keyed_used[i]);
      if (keyed_used[i] != keyed_want[i] || !isfinite(keyed_e[i]))
        ok = 0;
    }
  }

  setenv("CPMDC_ORBITAL_SLOT_LIMIT", "9", 1);
  cpmdc_session_select_orbitals(s, 303);
  CPMDCResult limited = cpmdc_session_energy_forces(
      s, NAT, limited_ctx.last_ang, z, f);
  if (limited.ok || strstr(limited.message, "CPMDC_ORBITAL_SLOT_LIMIT") == NULL ||
      strstr(limited.message, "history") == NULL ||
      strstr(limited.message, "unchanged") == NULL)
    ok = 0;
  setenv("CPMDC_ORBITAL_SLOT_LIMIT", "10", 1);
  limited =
      cpmdc_session_energy_forces(s, NAT, limited_ctx.last_ang, z, f);
  if (!limited.ok)
    ok = 0;
  unsetenv("CPMDC_ORBITAL_SLOT_LIMIT");
  printf("history_limit required=10 reject=9 accept=10 ok=%d\n", ok);
  cpmdc_session_destroy(s);
  s = NULL;

  limited_session_child(&limited_ctx, &after_full);
  if (!after_full.ok) {
    fprintf(stderr, "full-run limited_session failed: %s\n",
            after_full.message);
    ok = 0;
  }
  printf("limited_session full_run=%d\n", after_full.ok);
  rc = ok ? 0 : 1;

out:
  restore_stderr(&saved);
  if (cap)
    fclose(cap);
  cpmdc_session_destroy(s);
  printf("PASS=%d\n", rc == 0);
  return rc;
}
