/**
 * Two ranks of one calculator under orbital-preflight divergence.
 *
 * Every rank must enter the same collectives during session evaluation and
 * agree on keyed state, stale-orbital cleanup, and orbital preparation before
 * committing the session topology. Each case forces one rank to diverge and
 * asserts that:
 *
 *   - both ranks finish without hanging;
 *   - both ranks fail with the exact same error message;
 *   - the session topology stays uncommitted on both ranks;
 *   - a subsequent valid evaluation with a different topology succeeds.
 *
 * Cases:
 *   (a) an orbital key is selected on rank 0 only;
 *   (b) both ranks select a key, but the values differ, so agreement fails
 *       before prepare;
 *   (c) cpmdc_embed_drop_stale_orbitals fails on rank 1 only;
 *   (d) an unusable CPMDC_ORBITAL_SLOT_LIMIT on one rank leaves the topology
 *       uncommitted.
 */
#define _POSIX_C_SOURCE 200809L
#include "cpmdc.h"

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

/* Fortran embed runtime controls linked through cpmdc_embed_shell. */
int cpmdc_embed_init(void);
void cpmdc_embed_finalize(void);

enum {
  kWaterAtoms = 3,
  kWaterCoords = 3 * kWaterAtoms,
  kH2Atoms = 2,
  kH2Coords = 3 * kH2Atoms,
  kRanks = 2,
  kMessage = 512,
#if defined(CPMDC_HAS_CPMD)
  kAlarmSeconds = 240
#else
  kAlarmSeconds = 45
#endif
};

static const int kWaterZ[kWaterAtoms] = {8, 1, 1};
static const double kWaterGeom[kWaterCoords] = {
    0.0, 0.0, 0.1173, 0.0, 0.7572, -0.4692, 0.0, -0.7572, -0.4692};
static const int kH2Z[kH2Atoms] = {1, 1};
static const double kH2Geom[kH2Coords] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.74};
static const char kLimitName[] = "CPMDC_ORBITAL_SLOT_LIMIT";

static int g_rank = 0;
static int g_fail = 0;
static int g_case_fail = 0;
static int g_step = 0;
static unsigned char *g_params = NULL;
static size_t g_params_size = 0;

struct outcome {
  int ok;
  char message[kMessage];
  char last_error[kMessage];
};

static void check(int cond, const char *what) {
  if (cond)
    return;
  fprintf(stderr, "rank %d FAIL: %s\n", g_rank, what);
  g_fail = 1;
  g_case_fail = 1;
}

static int read_params(const char *path) {
  FILE *file = path ? fopen(path, "rb") : NULL;
  long length;
  if (!file)
    return -1;
  if (fseek(file, 0, SEEK_END) != 0 || (length = ftell(file)) <= 0 ||
      fseek(file, 0, SEEK_SET) != 0 ||
      !(g_params = (unsigned char *)malloc((size_t)length)) ||
      fread(g_params, 1, (size_t)length, file) != (size_t)length) {
    fclose(file);
    return -1;
  }
  fclose(file);
  g_params_size = (size_t)length;
  return 0;
}

static void water_geometry(double *pos) {
  memcpy(pos, kWaterGeom, sizeof(kWaterGeom));
  g_step++;
  pos[0] += 0.02 * (double)g_step;
}

static void h2_geometry(double *pos) {
  memcpy(pos, kH2Geom, sizeof(kH2Geom));
  g_step++;
  pos[5] += 0.01 * (double)g_step;
}

static void share_outcome(const struct outcome *mine,
                          struct outcome all[kRanks]) {
  MPI_Allgather(mine, (int)sizeof(*mine), MPI_BYTE, all, (int)sizeof(*mine),
                MPI_BYTE, MPI_COMM_WORLD);
}

static struct outcome evaluate_system(CPMDCSession *session, int n_atoms,
                                      const int *z, const double *pos) {
  double grad[kWaterCoords];
  struct outcome out;
  CPMDCResult r = cpmdc_session_energy_gradient(session, n_atoms, pos, z, grad);
  out.ok = r.ok;
  snprintf(out.message, sizeof(out.message), "%s", r.message);
  snprintf(out.last_error, sizeof(out.last_error), "%s", cpmdc_last_error());
  return out;
}

static void report(const char *name) {
  int any = 0;
  MPI_Allreduce(&g_case_fail, &any, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  if (g_rank == 0)
    printf("%s %s\n", any ? "FAIL" : "ok", name);
  g_case_fail = 0;
}

static void assert_same_failure(const struct outcome all[kRanks],
                                const char *expected_substr,
                                const char *label) {
  check(!all[0].ok && !all[1].ok, label);
  check(all[0].message[0] != '\0' && all[1].message[0] != '\0',
        "both ranks report a non-empty error");
  check(strcmp(all[0].message, all[1].message) == 0,
        "both ranks report the exact same error message");
  check(strcmp(all[0].last_error, all[0].message) == 0 &&
            strcmp(all[1].last_error, all[1].message) == 0,
        "cpmdc_last_error matches CPMDCResult.message on both ranks");
  if (expected_substr) {
    check(strstr(all[0].message, expected_substr) != NULL &&
              strstr(all[1].message, expected_substr) != NULL,
          expected_substr);
  }
}

/* Proves topology_fixed is still 0 on both ranks: mutating session params and
 * constant cutoff succeed, and evaluating a different topology (H2 instead of
 * H2O) succeeds on both ranks. */
static void assert_topology_uncommitted(CPMDCSession *session,
                                        long long follow_key,
                                        const char *label) {
  int rc_params = cpmdc_session_set_params(session, g_params, g_params_size);
  int rc_cutoff = cpmdc_session_set_constant_cutoff(session, 0.0, 0.0, 0.0);
  int rcs_params[kRanks];
  int rcs_cutoff[kRanks];
  double pos[kH2Coords];
  struct outcome mine;
  struct outcome all[kRanks];
  MPI_Allgather(&rc_params, 1, MPI_INT, rcs_params, 1, MPI_INT, MPI_COMM_WORLD);
  MPI_Allgather(&rc_cutoff, 1, MPI_INT, rcs_cutoff, 1, MPI_INT, MPI_COMM_WORLD);
  check(rcs_params[0] == 0 && rcs_params[1] == 0,
        "cpmdc_session_set_params succeeds while topology is uncommitted");
  check(rcs_cutoff[0] == 0 && rcs_cutoff[1] == 0,
        "cpmdc_session_set_constant_cutoff succeeds while topology is uncommitted");
  check(cpmdc_session_select_orbitals(session, follow_key) == 0,
        "select follow-up key on both ranks");
  h2_geometry(pos);
  mine = evaluate_system(session, kH2Atoms, kH2Z, pos);
  share_outcome(&mine, all);
  check(all[0].ok && all[1].ok, label);
}

/* (a) Key present on rank 0 only. */
static void case_a_key_on_rank_0_only(void) {
  CPMDCSession *session;
  double pos[kWaterCoords];
  struct outcome mine;
  struct outcome all[kRanks];
  unsetenv(kLimitName);
  session = cpmdc_session_create(g_params, g_params_size);
  check(session != NULL, "case_a session_create");
  if (!session)
    return;
  if (g_rank == 0)
    check(cpmdc_session_select_orbitals(session, 10) == 0,
          "rank 0 selects orbital key");
  water_geometry(pos);
  mine = evaluate_system(session, kWaterAtoms, kWaterZ, pos);
  share_outcome(&mine, all);
  assert_same_failure(all, "orbital key",
                      "key on rank 0 only fails both ranks");
  assert_topology_uncommitted(
      session, 10, "case_a topology stays uncommitted and accepts H2");
  cpmdc_session_destroy(session);
  report("case_a: key present on rank 0 only fails both ranks cleanly");
}

/* (b) Both ranks have a key, but the values differ, after a saturated limit
 * of 1. Key agreement runs before the slot check, so both ranks report the
 * disagreement and neither hangs in prepare. */
static void case_b_prepare_fails_on_one_rank(void) {
  CPMDCSession *prime;
  CPMDCSession *session;
  double pos[kWaterCoords];
  struct outcome mine;
  struct outcome all[kRanks];
  setenv(kLimitName, "1", 1);
  prime = cpmdc_session_create(g_params, g_params_size);
  check(prime != NULL, "case_b prime session_create");
  if (!prime)
    return;
  cpmdc_session_select_orbitals(prime, 10);
  water_geometry(pos);
  mine = evaluate_system(prime, kWaterAtoms, kWaterZ, pos);
  share_outcome(&mine, all);
  check(all[0].ok && all[1].ok, "case_b prime key 10 succeeds");

  cpmdc_session_select_orbitals(prime, 21);
  water_geometry(pos);
  mine = evaluate_system(prime, kWaterAtoms, kWaterZ, pos);
  share_outcome(&mine, all);
  check(all[0].ok && all[1].ok, "case_b prime key 21 parks key 10");
  cpmdc_session_destroy(prime);

  session = cpmdc_session_create(g_params, g_params_size);
  check(session != NULL, "case_b fresh session_create");
  if (!session)
    return;
  cpmdc_session_select_orbitals(session, g_rank == 0 ? 22 : 21);
  water_geometry(pos);
  mine = evaluate_system(session, kWaterAtoms, kWaterZ, pos);
  share_outcome(&mine, all);
  assert_same_failure(all, "orbital key",
                      "differing keys fail both ranks before prepare");
  assert_topology_uncommitted(
      session, 21, "case_b topology stays uncommitted and accepts H2");
  cpmdc_session_destroy(session);
  unsetenv(kLimitName);
  report("case_b: differing keys fail both ranks cleanly");
}

/* (c) cpmdc_embed_drop_stale_orbitals fails on rank 1 only. Calling
 * cpmdc_embed_finalize on rank 1 after cpmdc_session_create leaves
 * session->embed_configured set so configure_embed_from_session is skipped
 * and cpmdc_embed_drop_stale_orbitals returns 0 on rank 1 and 1 on rank 0. */
static void case_c_drop_stale_fails_on_one_rank(void) {
  CPMDCSession *session;
  double pos[kWaterCoords];
  struct outcome mine;
  struct outcome all[kRanks];
  unsetenv(kLimitName);
  session = cpmdc_session_create(g_params, g_params_size);
  check(session != NULL, "case_c session_create");
  if (!session)
    return;
  check(cpmdc_session_select_orbitals(session, 30) == 0,
        "case_c select key 30 on both ranks");
  if (g_rank == 1)
    cpmdc_embed_finalize();
  water_geometry(pos);
  mine = evaluate_system(session, kWaterAtoms, kWaterZ, pos);
  if (g_rank == 1)
    check(cpmdc_embed_init() == 1, "rank 1 restores embed runtime");
  share_outcome(&mine, all);
  assert_same_failure(all, "stale orbitals",
                      "drop_stale_orbitals failure on rank 1 fails both ranks");
  assert_topology_uncommitted(
      session, 30, "case_c topology stays uncommitted and accepts H2");
  cpmdc_session_destroy(session);
  report("case_c: drop_stale_orbitals failing on one rank fails both cleanly");
}

/* (d) Unusable slot limit on one rank fails both ranks and leaves the fresh
 * session topology uncommitted. */
static void case_d_bad_limit_leaves_topology_uncommitted(void) {
  CPMDCSession *session;
  double pos[kWaterCoords];
  struct outcome mine;
  struct outcome all[kRanks];
  setenv(kLimitName, g_rank == 0 ? "1" : "abc", 1);
  session = cpmdc_session_create(g_params, g_params_size);
  check(session != NULL, "case_d session_create");
  if (!session)
    return;
  check(cpmdc_session_select_orbitals(session, 40) == 0,
        "case_d select key 40");
  water_geometry(pos);
  mine = evaluate_system(session, kWaterAtoms, kWaterZ, pos);
  share_outcome(&mine, all);
  assert_same_failure(all, "positive integer",
                      "bad limit on rank 1 fails both ranks");
  unsetenv(kLimitName);
  assert_topology_uncommitted(
      session, 40, "case_d topology stays uncommitted and accepts H2");
  cpmdc_session_destroy(session);
  report("case_d: bad limit on one rank leaves topology uncommitted");
}

int main(int argc, char **argv) {
  MPI_Comm world = MPI_COMM_WORLD;
  const char *only = argc > 2 ? argv[2] : NULL;
  int size = 0;
  int any = 0;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &g_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  alarm(kAlarmSeconds);
  if (size != kRanks || read_params(argc > 1 ? argv[1] : NULL) != 0) {
    if (g_rank == 0)
      fprintf(stderr, "needs 2 ranks and a CPMDParams file\n");
    MPI_Finalize();
    return 2;
  }
  check(cpmdc_adopt_calculator_comm(&world, sizeof(world), kRanks) == 0,
        "adopt the world as one calculator");
  unsetenv("CPMDC_RANK_CONSENSUS");

  if (!only || strcmp(only, "case_a") == 0)
    case_a_key_on_rank_0_only();
  if (!only || strcmp(only, "case_b") == 0)
    case_b_prepare_fails_on_one_rank();
  if (!only || strcmp(only, "case_c") == 0)
    case_c_drop_stale_fails_on_one_rank();
  if (!only || strcmp(only, "case_d") == 0)
    case_d_bad_limit_leaves_topology_uncommitted();

  free(g_params);
  MPI_Allreduce(&g_fail, &any, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Finalize();
  return any ? 1 : 0;
}
