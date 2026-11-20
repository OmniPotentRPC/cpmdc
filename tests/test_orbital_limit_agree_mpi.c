/**
 * Two ranks of one calculator and CPMDC_ORBITAL_SLOT_LIMIT.
 *
 * The ranks agree the limit in a collective that every evaluation and every
 * orbital copy runs after its own argument checks. This test drives each
 * way that agreement can go and checks that both ranks leave the call with
 * the answer the rules give, and that no rank is left waiting:
 *
 *   - a limit that is valid and equal on both ranks;
 *   - an unusable or different limit on one or both ranks;
 *   - a bad call that fails the same way on both ranks, which reports the
 *     bad call and not the limit;
 *   - a bad call that fails on one rank only, which fails the other rank
 *     too through the same collective;
 *   - the next valid call, which proves the ranks stayed in step.
 *
 * A hung rank ends the test through the alarm below.
 */
#define _POSIX_C_SOURCE 200809L
#include "cpmdc.h"

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

enum {
  kAtoms = 3,
  kCoords = 3 * kAtoms,
  kRanks = 2,
  kMessage = 256,
  kAlarmSeconds = 100
};

static const int kZ[kAtoms] = {8, 1, 1};
static const int kBadZ[kAtoms] = {8, 0, 1};
static const double kGeometry[kCoords] = {0.0, 0.0,     0.1173,  0.0,    0.7572,
                                          -0.4692, 0.0, -0.7572, -0.4692};
static const char kLimitName[] = "CPMDC_ORBITAL_SLOT_LIMIT";

static int g_rank = 0;
static int g_fail = 0;
static int g_step = 0;
static unsigned char *g_params = NULL;
static size_t g_params_size = 0;
static unsigned char *g_force_input = NULL;
static size_t g_force_input_size = 0;

struct outcome {
  int ok;
  char message[kMessage];
};

static void check(int cond, const char *what) {
  if (cond)
    return;
  fprintf(stderr, "rank %d FAIL: %s\n", g_rank, what);
  g_fail = 1;
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

static int read_force_input(const char *path) {
  FILE *file = path ? fopen(path, "rb") : NULL;
  long length;
  if (!file)
    return -1;
  if (fseek(file, 0, SEEK_END) != 0 || (length = ftell(file)) <= 0 ||
      fseek(file, 0, SEEK_SET) != 0 ||
      !(g_force_input = (unsigned char *)malloc((size_t)length)) ||
      fread(g_force_input, 1, (size_t)length, file) != (size_t)length) {
    fclose(file);
    return -1;
  }
  fclose(file);
  g_force_input_size = (size_t)length;
  return 0;
}

static void set_limit(const char *per_rank[kRanks]) {
  const char *mine = per_rank[g_rank];
  if (mine)
    setenv(kLimitName, mine, 1);
  else
    unsetenv(kLimitName);
}

static void geometry(double *pos) {
  memcpy(pos, kGeometry, sizeof(kGeometry));
  g_step++;
  pos[0] += 0.01 * (double)g_step;
}

static void record(struct outcome *out, int ok, const char *message) {
  out->ok = ok;
  snprintf(out->message, sizeof(out->message), "%s", message ? message : "");
}

/* Both ranks hold the same arrays after this, in rank order. */
static void share(const struct outcome *mine, struct outcome all[kRanks]) {
  MPI_Allgather(mine, (int)sizeof(*mine), MPI_BYTE, all, (int)sizeof(*mine),
                MPI_BYTE, MPI_COMM_WORLD);
}

static int mentions(const struct outcome *o, const char *text) {
  return strstr(o->message, text) != NULL;
}

static void report(const char *name) {
  int any = 0;
  MPI_Allreduce(&g_fail, &any, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  if (g_rank == 0)
    printf("%s %s\n", any ? "FAIL" : "ok", name);
}

static struct outcome evaluate(CPMDCSession *session, const int *z,
                               const double *pos) {
  double grad[kCoords];
  struct outcome out;
  CPMDCResult r = cpmdc_session_energy_gradient(session, kAtoms, pos, z, grad);
  record(&out, r.ok, r.message);
  return out;
}

static void case_agreed_limit(void) {
  const char *limits[kRanks] = {"2", "2"};
  struct outcome mine;
  struct outcome all[kRanks];
  double pos[kCoords];
  CPMDCSession *session;
  set_limit(limits);
  setenv("CPMDC_CELL_IN_PLACE", g_rank == 0 ? "0" : "1", 1);
  session = cpmdc_session_create(g_params, g_params_size);
  check(session != NULL, "session_create");
  if (!session)
    return;
  cpmdc_session_select_orbitals(session, 0);
  geometry(pos);
  mine = evaluate(session, kZ, pos);
  share(&mine, all);
  check(all[0].ok && all[1].ok, "an equal valid limit lets both ranks evaluate");
  cpmdc_session_destroy(session);
  unsetenv("CPMDC_CELL_IN_PLACE");
  report("root broadcasts CPMDC_CELL_IN_PLACE and an equal limit evaluates");
}

static void case_limit_is_refused_together(const char *label,
                                           const char *limits[kRanks],
                                           const char *expected) {
  struct outcome mine;
  struct outcome all[kRanks];
  double pos[kCoords];
  CPMDCSession *session;
  set_limit(limits);
  session = cpmdc_session_create(g_params, g_params_size);
  check(session != NULL, "session_create");
  if (!session)
    return;
  cpmdc_session_select_orbitals(session, 0);
  geometry(pos);
  mine = evaluate(session, kZ, pos);
  share(&mine, all);
  check(!all[0].ok && !all[1].ok, label);
  check(mentions(&all[0], kLimitName) && mentions(&all[1], kLimitName), label);
  check(mentions(&all[0], expected) && mentions(&all[1], expected), label);
  cpmdc_session_destroy(session);
  report(label);
}

static void case_null_session_on_both_ranks(void) {
  const char *limits[kRanks] = {"abc", "abc"};
  struct outcome mine;
  struct outcome all[kRanks];
  double pos[kCoords];
  double grad[kCoords];
  int rc;
  int rcs[kRanks];
  set_limit(limits);
  geometry(pos);
  {
    CPMDCResult r = cpmdc_session_energy_gradient(NULL, kAtoms, pos, kZ, grad);
    record(&mine, r.ok, r.message);
  }
  share(&mine, all);
  check(!all[0].ok && !all[1].ok, "a null session fails both ranks");
  check(!strcmp(all[0].message, "null session") &&
            !strcmp(all[1].message, "null session"),
        "a null session is reported, not the limit");
  rc = cpmdc_session_copy_orbitals(NULL, 1, 2);
  MPI_Allgather(&rc, 1, MPI_INT, rcs, 1, MPI_INT, MPI_COMM_WORLD);
  check(rcs[0] == -1 && rcs[1] == -1, "a null session copy fails both ranks");
  check(!strcmp(cpmdc_last_error(), "null session"),
        "a null session copy is reported, not the limit");
  report("a null session on both ranks reports the call, not the limit");
}

/* One rank passes a bad argument. The other rank, whose call is well
 * formed, must still fail instead of waiting in the collective. */
static void case_bad_call_on_one_rank(void) {
  const char *limits[kRanks] = {"1", "1"};
  struct outcome mine;
  struct outcome all[kRanks];
  struct outcome next[kRanks];
  struct outcome after;
  double pos[kCoords];
  CPMDCSession *session;
  int configure_rc;
  int configure_all[kRanks];
  set_limit(limits);
  session = cpmdc_session_create(g_params, g_params_size);
  check(session != NULL, "session_create");
  if (!session)
    return;
  cpmdc_session_select_orbitals(session, 0);
  geometry(pos);
  mine = evaluate(session, g_rank == 1 ? kBadZ : kZ, pos);
  share(&mine, all);
  check(!all[0].ok && !all[1].ok, "a bad call on one rank fails both ranks");
  check(mentions(&all[0], "atomic numbers") &&
            mentions(&all[1], "atomic numbers"),
        "both ranks report the agreed call error");
  report("a bad call on one rank fails both ranks");

  configure_rc = cpmdc_session_set_params(session, g_params, g_params_size);
  MPI_Allgather(&configure_rc, 1, MPI_INT, configure_all, 1, MPI_INT,
                MPI_COMM_WORLD);
  check(configure_all[0] == 0 && configure_all[1] == 0,
        "a rejected first call does not latch topology on either rank");
  report("a rejected first call leaves topology unfixed on every rank");

  geometry(pos);
  after = evaluate(session, kZ, pos);
  share(&after, next);
  check(next[0].ok && next[1].ok, "the ranks stay in step after the failure");
  cpmdc_session_destroy(session);
  report("the ranks stay in step after a failure on one rank");
}

/* Same, with a limit that is itself unusable: the rank with the bad call
 * keeps its own reason and the other rank still does not wait. */
static void case_bad_call_and_bad_limit(void) {
  const char *limits[kRanks] = {"abc", "abc"};
  struct outcome mine;
  struct outcome all[kRanks];
  double pos[kCoords];
  CPMDCSession *session;
  set_limit(limits);
  session = cpmdc_session_create(g_params, g_params_size);
  check(session != NULL, "session_create");
  if (!session)
    return;
  geometry(pos);
  mine = evaluate(session, g_rank == 1 ? kBadZ : kZ, pos);
  share(&mine, all);
  check(!all[0].ok && !all[1].ok, "a bad call and a bad limit fail both ranks");
  check(mentions(&all[1], "atomic numbers"),
        "the rank with the bad call keeps its own reason");
  check(!mentions(&all[1], kLimitName),
        "the rank with the bad call does not report the limit");
  cpmdc_session_destroy(session);
  report("a bad call keeps its own reason beside a bad limit");
}

static void case_orbital_fault(const char *variable, const char *fault,
                               const char *label) {
  const char *limits[kRanks] = {"2", "2"};
  struct outcome mine;
  struct outcome all[kRanks];
  double pos[kCoords];
  CPMDCSession *session;
  set_limit(limits);
  session = cpmdc_session_create(g_params, g_params_size);
  check(session != NULL, "session_create");
  if (!session)
    return;
  if (!strcmp(fault, "selector-allocation")) {
    cpmdc_session_select_orbitals(session, 70);
    geometry(pos);
    mine = evaluate(session, kZ, pos);
    share(&mine, all);
    check(all[0].ok && all[1].ok,
          "selector allocation fault has stored orbitals to park");
  }
  cpmdc_session_select_orbitals(session, 71);
  if (g_rank == 1)
    setenv(variable, fault, 1);
  else
    unsetenv(variable);
  geometry(pos);
  mine = evaluate(session, kZ, pos);
  share(&mine, all);
  check(!all[0].ok && !all[1].ok, label);
  unsetenv(variable);
  geometry(pos);
  mine = evaluate(session, kZ, pos);
  share(&mine, all);
  check(all[0].ok && all[1].ok,
        "ranks stay in step after an orbital preflight failure");
  cpmdc_session_destroy(session);
  report(label);
}

static void case_orbital_key_mismatch(void) {
  const char *limits[kRanks] = {"4", "4"};
  struct outcome mine;
  struct outcome all[kRanks];
  double pos[kCoords];
  CPMDCSession *session;
  int rc;
  int rcs[kRanks];
  set_limit(limits);
  session = cpmdc_session_create(g_params, g_params_size);
  check(session != NULL, "session_create");
  if (!session)
    return;

  cpmdc_session_select_orbitals(session, g_rank == 0 ? 81 : 82);
  geometry(pos);
  mine = evaluate(session, kZ, pos);
  share(&mine, all);
  check(!all[0].ok && !all[1].ok,
        "different selected keys fail both ranks");
  check(mentions(&all[0], "differs") && mentions(&all[1], "differs"),
        "different selected keys report the mismatch");

  cpmdc_session_destroy(session);
  session = cpmdc_session_create(g_params, g_params_size);
  check(session != NULL, "session_create");
  if (!session)
    return;
  if (g_rank == 1)
    cpmdc_session_select_orbitals(session, 83);
  geometry(pos);
  mine = evaluate(session, kZ, pos);
  share(&mine, all);
  check(!all[0].ok && !all[1].ok,
        "rank-different selected-key presence fails both ranks");

  cpmdc_session_select_orbitals(session, 84);
  geometry(pos);
  mine = evaluate(session, kZ, pos);
  share(&mine, all);
  check(all[0].ok && all[1].ok,
        "matching selected keys recover after a mismatch");

  rc = cpmdc_session_copy_orbitals(
      session, g_rank == 0 ? 84 : 85, g_rank == 0 ? 86 : 87);
  MPI_Allgather(&rc, 1, MPI_INT, rcs, 1, MPI_INT, MPI_COMM_WORLD);
  check(rcs[0] == -1 && rcs[1] == -1,
        "different copy key pairs fail both ranks");
  check(strstr(cpmdc_last_error(), "differs") != NULL,
        "different copy key pairs report the mismatch");
  cpmdc_session_destroy(session);
  report("orbital key mismatches fail collectively");
}

static void check_collective_failure(CPMDCResult result, const char *label) {
  struct outcome mine;
  struct outcome all[kRanks];
  record(&mine, result.ok, result.message);
  share(&mine, all);
  check(!all[0].ok && !all[1].ok, label);
  report(label);
}

static void case_wrapper_preflights(void) {
  const char *limits[kRanks] = {"2", "2"};
  double pos[kCoords];
  double energy = 0.0;
  double forces[kCoords];
  double virial[9];
  size_t result_size;
  unsigned char *result_buf;
  CPMDCSession *session;
  CPMDCResult result;
  int rc;
  int rcs[kRanks];
  set_limit(limits);
  geometry(pos);

  result = cpmdc_energy(g_rank == 1 ? 0 : kAtoms, pos, kZ, g_params,
                        g_params_size);
  check_collective_failure(result, "cpmdc_energy agrees a one-rank error");
  if (g_rank == 1)
    setenv("CPMDC_TEST_FAULT", "energy-allocation", 1);
  result = cpmdc_energy(kAtoms, pos, kZ, g_params, g_params_size);
  check_collective_failure(
      result, "cpmdc_energy agrees a one-rank allocation failure");
  unsetenv("CPMDC_TEST_FAULT");

  session = cpmdc_session_create(g_params, g_params_size);
  check(session != NULL, "session_create");
  if (!session)
    return;
  result = cpmdc_session_energy(session, g_rank == 1 ? 0 : kAtoms, pos, kZ);
  check_collective_failure(result,
                           "cpmdc_session_energy agrees a one-rank error");
  if (g_rank == 1)
    setenv("CPMDC_TEST_FAULT", "session-energy-allocation", 1);
  result = cpmdc_session_energy(session, kAtoms, pos, kZ);
  check_collective_failure(
      result, "cpmdc_session_energy agrees a one-rank allocation failure");
  unsetenv("CPMDC_TEST_FAULT");
  cpmdc_session_destroy(session);

  session = cpmdc_session_create(g_params, g_params_size);
  check(session != NULL, "session_create");
  if (!session)
    return;
  result = cpmdc_session_calculate_forces(
      session, g_force_input, g_force_input_size, forces,
      g_rank == 1 ? 0 : kCoords);
  check_collective_failure(
      result, "cpmdc_session_calculate_forces agrees a one-rank error");
  if (g_rank == 1)
    setenv("CPMDC_TEST_FAULT", "calculate-forces-reserve", 1);
  result = cpmdc_session_calculate_forces(
      session, g_force_input, g_force_input_size, forces, kCoords);
  check_collective_failure(
      result, "calculate_forces agrees a one-rank reserve failure");
  unsetenv("CPMDC_TEST_FAULT");
  cpmdc_session_destroy(session);

  result_size = cpmdc_potential_result_size_for_force_input(
      g_force_input, g_force_input_size);
  result_buf = (unsigned char *)malloc(result_size);
  check(result_size > 0 && result_buf != NULL, "allocate PotentialResult");
  session = cpmdc_session_create(g_params, g_params_size);
  check(session != NULL, "session_create");
  if (!session || !result_buf) {
    cpmdc_session_destroy(session);
    free(result_buf);
    return;
  }
  {
    size_t written = 0;
    result = cpmdc_session_calculate_result(
        session, g_force_input, g_force_input_size, result_buf,
        g_rank == 1 ? 0 : result_size, &written);
  }
  check_collective_failure(
      result, "cpmdc_session_calculate_result agrees a one-rank error");
  if (g_rank == 1)
    setenv("CPMDC_TEST_FAULT", "calculate-result-reserve", 1);
  {
    size_t written = 0;
    result = cpmdc_session_calculate_result(
        session, g_force_input, g_force_input_size, result_buf, result_size,
        &written);
  }
  check_collective_failure(
      result, "calculate_result agrees a one-rank reserve failure");
  unsetenv("CPMDC_TEST_FAULT");
  cpmdc_session_destroy(session);
  free(result_buf);

  result_buf = (unsigned char *)malloc(result_size);
  check(result_buf != NULL, "allocate PotentialResult for write failure");
  session = cpmdc_session_create(g_params, g_params_size);
  check(session != NULL, "session_create");
  if (!session || !result_buf) {
    cpmdc_session_destroy(session);
    free(result_buf);
    return;
  }
  if (g_rank == 1)
    setenv("CPMDC_TEST_FAULT", "potential-result-write", 1);
  else
    unsetenv("CPMDC_TEST_FAULT");
  {
    size_t written = 0;
    result = cpmdc_session_calculate_result(
        session, g_force_input, g_force_input_size, result_buf, result_size,
        &written);
  }
  check_collective_failure(
      result, "PotentialResult write failure is agreed across ranks");
  unsetenv("CPMDC_TEST_FAULT");
  cpmdc_session_destroy(session);
  free(result_buf);

  session = cpmdc_session_create(g_params, g_params_size);
  check(session != NULL, "session_create");
  if (!session)
    return;
  rc = cpmdc_session_calculate_au(
      session, kAtoms, g_rank == 1 ? NULL : pos, kZ, NULL, &energy, forces,
      virial);
  MPI_Allgather(&rc, 1, MPI_INT, rcs, 1, MPI_INT, MPI_COMM_WORLD);
  check(rcs[0] == -1 && rcs[1] == -1,
        "cpmdc_session_calculate_au agrees a one-rank error");
  if (g_rank == 1)
    setenv("CPMDC_TEST_FAULT", "calculate-au-reserve", 1);
  rc = cpmdc_session_calculate_au(session, kAtoms, pos, kZ, NULL, &energy,
                                  forces, virial);
  MPI_Allgather(&rc, 1, MPI_INT, rcs, 1, MPI_INT, MPI_COMM_WORLD);
  check(rcs[0] == -1 && rcs[1] == -1,
        "calculate_au agrees a one-rank reserve failure");
  unsetenv("CPMDC_TEST_FAULT");
  cpmdc_session_destroy(session);
  report("cpmdc_session_calculate_au agrees a one-rank error");
}

#if !defined(CPMDC_HAS_CPMD)
/* The reference evaluator holds the parked keys in the library, so the
 * third key refuses on both ranks and the first two stay usable. */
static void case_limit_reached_on_both_ranks(void) {
  const char *limits[kRanks] = {"1", "1"};
  struct outcome mine[3];
  struct outcome all[3][kRanks];
  double pos[kCoords];
  CPMDCSession *session;
  int i;
  set_limit(limits);
  session = cpmdc_session_create(g_params, g_params_size);
  check(session != NULL, "session_create");
  if (!session)
    return;
  for (i = 0; i < 3; ++i) {
    cpmdc_session_select_orbitals(session, i);
    geometry(pos);
    mine[i] = evaluate(session, kZ, pos);
    share(&mine[i], all[i]);
  }
  check(all[0][0].ok && all[0][1].ok, "the first key evaluates on both ranks");
  check(all[1][0].ok && all[1][1].ok, "the second key evaluates on both ranks");
  check(!all[2][0].ok && !all[2][1].ok,
        "a park past the limit fails both ranks");
  check(mentions(&all[2][0], "exceeded") && mentions(&all[2][1], "exceeded"),
        "both ranks name the exceeded limit");
  cpmdc_session_select_orbitals(session, 1);
  geometry(pos);
  mine[0] = evaluate(session, kZ, pos);
  share(&mine[0], all[0]);
  check(all[0][0].ok && all[0][1].ok, "a key already parked still evaluates");
  cpmdc_session_destroy(session);
  report("a park past the limit fails both ranks and keeps the stored keys");
}
#endif

int main(int argc, char **argv) {
  MPI_Comm world = MPI_COMM_WORLD;
  int size = 0;
  int any = 0;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &g_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  alarm(kAlarmSeconds);
  if (size != kRanks || read_params(argc > 1 ? argv[1] : NULL) != 0 ||
      read_force_input(argc > 2 ? argv[2] : NULL) != 0) {
    if (g_rank == 0)
      fprintf(stderr, "needs 2 ranks, CPMDParams and ForceInput files\n");
    MPI_Finalize();
    return 2;
  }
  check(cpmdc_adopt_calculator_comm(&world, sizeof(world), kRanks) == 0,
        "adopt the world as one calculator");
  unsetenv("CPMDC_RANK_CONSENSUS");

  case_agreed_limit();
  {
    const char *both_bad[kRanks] = {"abc", "abc"};
    const char *one_bad[kRanks] = {"1", "0"};
    const char *empty[kRanks] = {"", ""};
    const char *differ[kRanks] = {"1", "2"};
    const char *one_unset[kRanks] = {"1", NULL};
    case_limit_is_refused_together("an unusable limit on both ranks fails both",
                                   both_bad, "positive integer");
    case_limit_is_refused_together("an unusable limit on one rank fails both",
                                   one_bad, "positive integer");
    case_limit_is_refused_together("an empty limit fails both ranks", empty,
                                   "positive integer");
    case_limit_is_refused_together("different limits fail both ranks", differ,
                                   "differs");
    case_limit_is_refused_together("a limit on one rank only fails both",
                                   one_unset, "differs");
  }
  case_null_session_on_both_ranks();
  case_bad_call_on_one_rank();
  case_bad_call_and_bad_limit();
#if !defined(CPMDC_HAS_CPMD)
  case_limit_reached_on_both_ranks();
#endif
  case_orbital_fault("CPMDC_TEST_FAULT", "drop-stale",
                     "a stale-orbital failure is agreed");
  case_orbital_fault("CPMDC_TEST_ORBITAL_FAULT", "prepare",
                     "an orbital-prepare failure is agreed");
#if defined(CPMDC_HAS_CPMD)
  case_orbital_fault("CPMDC_TEST_ORBITAL_FAULT", "selector-allocation",
                     "a selector allocation failure is agreed");
#endif
  case_orbital_key_mismatch();
  case_wrapper_preflights();

  free(g_params);
  free(g_force_input);
  MPI_Allreduce(&g_fail, &any, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Finalize();
  return any ? 1 : 0;
}
