/**
 * Two calculator ranks with periodic cell and virial errors on one rank only.
 *
 * When a single rank receives a non-finite, rank-deficient, underflowing, or
 * overflowing periodic cell, or computes a non-finite virial after the step,
 * both ranks must return the agreed error without stranding a peer in an MPI
 * collective, leave topology unfixed on a rejected first call, and stay in
 * step for the next valid evaluation.
 */
#define _POSIX_C_SOURCE 200809L
#include "cpmdc.h"

#include "schema/Potentials.capnp.h"
#include <capn.h>
#include <float.h>
#include <math.h>
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
  kAlarmSeconds = 30
};

static const int kZ[kAtoms] = {8, 1, 1};
static const double kWaterBohr[kCoords] = {
    0.0, 0.0, 0.2216669, 0.0, 1.4308949, -0.8866678, 0.0, -1.4308949,
    -0.8866678};
static const double kWaterAng[kCoords] = {0.0,    0.0,     0.1173, 0.0,    0.7572,
                                          -0.4692, 0.0, -0.7572, -0.4692};
static const double kGoodCellBohr[9] = {18.9, 0.0, 0.0, 0.0, 18.9,
                                        0.0,  0.0, 0.0, 18.9};
static const double kGoodCellAng[9] = {10.0, 0.0, 0.0, 0.0, 10.0,
                                       0.0,  0.0, 0.0, 10.0};

static int g_rank = 0;
static int g_fail = 0;
static int g_step = 0;
static unsigned char *g_params = NULL;
static size_t g_params_size = 0;

struct outcome {
  int rc;
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

static void shifted_bohr(double *pos) {
  memcpy(pos, kWaterBohr, sizeof(kWaterBohr));
  g_step++;
  pos[0] += 0.01 * (double)g_step;
}

static void shifted_ang(double *pos) {
  memcpy(pos, kWaterAng, sizeof(kWaterAng));
  g_step++;
  pos[0] += 0.01 * (double)g_step;
}

static void record(struct outcome *out, int rc, const char *message) {
  out->rc = rc;
  snprintf(out->message, sizeof(out->message), "%s", message ? message : "");
}

static void share(const struct outcome *mine, struct outcome all[kRanks]) {
  MPI_Allgather(mine, (int)sizeof(*mine), MPI_BYTE, all, (int)sizeof(*mine),
                MPI_BYTE, MPI_COMM_WORLD);
}

static void report(const char *name) {
  int any = 0;
  MPI_Allreduce(&g_fail, &any, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  if (g_rank == 0)
    printf("%s %s\n", any ? "FAIL" : "ok", name);
}

static int encode_force_input(const double *pos_ang, const double *cell_ang,
                              unsigned char *out, size_t cap) {
  struct capn c;
  capn_init_malloc(&c);
  capn_ptr root = capn_root(&c);
  ForceInput_ptr fp = new_ForceInput(root.seg);
  struct ForceInput v;
  memset(&v, 0, sizeof(v));
  v.pos = capn_new_list64(root.seg, kCoords);
  v.atmnrs = capn_new_list32(root.seg, kAtoms);
  v.box = capn_new_list64(root.seg, 9);
  for (int i = 0; i < kCoords; ++i)
    capn_set64(v.pos, i, capn_from_f64(pos_ang[i]));
  for (int i = 0; i < kAtoms; ++i)
    capn_set32(v.atmnrs, i, (uint32_t)kZ[i]);
  for (int i = 0; i < 9; ++i)
    capn_set64(v.box, i, capn_from_f64(cell_ang[i]));
  v.lengthUnit = (capn_text){8, "angstrom", NULL};
  v.energyUnit = (capn_text){7, "hartree", NULL};
  write_ForceInput(&v, fp);
  int len =
      capn_setp(root, 0, fp.p) == 0 ? capn_write_mem(&c, out, cap, 0) : -1;
  capn_free(&c);
  return len;
}

static void check_au_bad_cell_one_rank(int bad_rank, const double bad_cell[9],
                                       int null_peer_virial,
                                       const char *expected,
                                       const char *label) {
  CPMDCSession *session = cpmdc_session_create(g_params, g_params_size);
  double pos[kCoords];
  double energy = -7.0;
  double forces[kCoords];
  double virial[9];
  struct outcome mine;
  struct outcome all[kRanks];
  int cfg_rc;
  int cfg_all[kRanks];
  check(session != NULL, "session_create");
  if (!session)
    return;

  shifted_bohr(pos);
  {
    const double *cell = g_rank == bad_rank ? bad_cell : kGoodCellBohr;
    double *vir =
        (null_peer_virial && g_rank != bad_rank) ? NULL : virial;
    int rc = cpmdc_session_calculate_au(session, kAtoms, pos, kZ, cell, &energy,
                                        forces, vir);
    record(&mine, rc, cpmdc_last_error());
  }
  share(&mine, all);
  check(all[0].rc == -1 && all[1].rc == -1, label);
  check(strstr(all[0].message, expected) != NULL &&
            strstr(all[1].message, expected) != NULL,
        label);

  cfg_rc = cpmdc_session_set_params(session, g_params, g_params_size);
  MPI_Allgather(&cfg_rc, 1, MPI_INT, cfg_all, 1, MPI_INT, MPI_COMM_WORLD);
  check(cfg_all[0] == 0 && cfg_all[1] == 0,
        "a rejected periodic preflight leaves session topology configurable");

  shifted_bohr(pos);
  {
    int rc = cpmdc_session_calculate_au(session, kAtoms, pos, kZ, kGoodCellBohr,
                                        &energy, forces, virial);
    record(&mine, rc, cpmdc_last_error());
  }
  share(&mine, all);
  check(all[0].rc == 0 && all[1].rc == 0,
        "ranks stay in step after a one-rank cell failure");
  check(isfinite(energy) && isfinite(virial[0]),
        "follow-up periodic evaluation writes finite energy and virial");

  cpmdc_session_destroy(session);
  report(label);
}

static void check_capnp_bad_cell_one_rank(int bad_rank,
                                          const double bad_cell_ang[9],
                                          const char *expected,
                                          const char *label) {
  CPMDCSession *session = cpmdc_session_create(g_params, g_params_size);
  unsigned char msg[4096];
  unsigned char result_buf[8192];
  double pos[kCoords];
  double forces[kCoords];
  struct outcome mine;
  struct outcome all[kRanks];
  int len;
  int cfg_rc;
  int cfg_all[kRanks];
  size_t written = 0;
  check(session != NULL, "session_create");
  if (!session)
    return;

  shifted_ang(pos);
  len = encode_force_input(
      pos, g_rank == bad_rank ? bad_cell_ang : kGoodCellAng, msg, sizeof(msg));
  check(len > 0, "encode ForceInput");
  {
    CPMDCResult r = cpmdc_session_calculate_forces(session, msg, (size_t)len,
                                                   forces, kCoords);
    record(&mine, r.ok ? 0 : -1, r.message);
  }
  share(&mine, all);
  check(all[0].rc == -1 && all[1].rc == -1, label);
  check(strstr(all[0].message, expected) != NULL &&
            strstr(all[1].message, expected) != NULL,
        label);

  {
    CPMDCResult r = cpmdc_session_calculate_result(
        session, msg, (size_t)len, result_buf, sizeof(result_buf), &written);
    record(&mine, r.ok ? 0 : -1, r.message);
  }
  share(&mine, all);
  check(all[0].rc == -1 && all[1].rc == -1, label);
  check(strstr(all[0].message, expected) != NULL &&
            strstr(all[1].message, expected) != NULL,
        label);

  cfg_rc = cpmdc_session_set_params(session, g_params, g_params_size);
  MPI_Allgather(&cfg_rc, 1, MPI_INT, cfg_all, 1, MPI_INT, MPI_COMM_WORLD);
  check(cfg_all[0] == 0 && cfg_all[1] == 0,
        "a rejected ForceInput cell leaves session topology configurable");

  shifted_ang(pos);
  len = encode_force_input(pos, kGoodCellAng, msg, sizeof(msg));
  check(len > 0, "encode valid ForceInput");
  {
    CPMDCResult r = cpmdc_session_calculate_result(
        session, msg, (size_t)len, result_buf, sizeof(result_buf), &written);
    record(&mine, r.ok ? 0 : -1, r.message);
  }
  share(&mine, all);
  check(all[0].rc == 0 && all[1].rc == 0,
        "ranks stay in step after a one-rank ForceInput cell failure");

  cpmdc_session_destroy(session);
  report(label);
}

static void case_fault_nonfinite_virial_one_rank(int null_peer_virial,
                                                 const char *label) {
  CPMDCSession *session = cpmdc_session_create(g_params, g_params_size);
  double pos[kCoords];
  double energy = -9.0;
  double forces[kCoords];
  double virial[9];
  struct outcome mine;
  struct outcome all[kRanks];
  check(session != NULL, "session_create");
  if (!session)
    return;

  for (int i = 0; i < 9; ++i)
    virial[i] = -3.0;
  shifted_bohr(pos);
  if (g_rank == 1)
    setenv("CPMDC_TEST_FAULT", "calculate-au-virial", 1);
  else
    unsetenv("CPMDC_TEST_FAULT");
  {
    double *vir = (null_peer_virial && g_rank == 0) ? NULL : virial;
    int rc = cpmdc_session_calculate_au(session, kAtoms, pos, kZ, kGoodCellBohr,
                                        &energy, forces, vir);
    record(&mine, rc, cpmdc_last_error());
  }
  unsetenv("CPMDC_TEST_FAULT");
  share(&mine, all);
  check(all[0].rc == -1 && all[1].rc == -1, label);
  check(!strcmp(all[0].message, "virial is not finite") &&
            !strcmp(all[1].message, "virial is not finite"),
        label);
  check(energy == -9.0 && virial[0] == -3.0,
        "failed virial step leaves energy and virial outputs untouched");

  shifted_bohr(pos);
  {
    int rc = cpmdc_session_calculate_au(session, kAtoms, pos, kZ, kGoodCellBohr,
                                        &energy, forces, virial);
    record(&mine, rc, cpmdc_last_error());
  }
  share(&mine, all);
  check(all[0].rc == 0 && all[1].rc == 0,
        "ranks stay in step after a one-rank virial failure");

  cpmdc_session_destroy(session);
  report(label);
}

static void case_fault_evaluation_one_rank(void) {
  CPMDCSession *session = cpmdc_session_create(g_params, g_params_size);
  double pos[kCoords];
  double energy = -9.0;
  double forces[kCoords];
  double virial[9];
  struct outcome mine;
  struct outcome all[kRanks];
  check(session != NULL, "session_create");
  if (!session)
    return;

  for (int i = 0; i < 9; ++i)
    virial[i] = -3.0;
  shifted_bohr(pos);
  if (g_rank == 1)
    setenv("CPMDC_TEST_FAULT", "embed-energy-gradient", 1);
  else
    unsetenv("CPMDC_TEST_FAULT");
  {
    int rc = cpmdc_session_calculate_au(session, kAtoms, pos, kZ, kGoodCellBohr,
                                        &energy, forces, virial);
    record(&mine, rc, cpmdc_last_error());
  }
  unsetenv("CPMDC_TEST_FAULT");
  share(&mine, all);
  check(all[0].rc == -1 && all[1].rc == -1,
        "calculate_au agrees a one-rank evaluation failure");
  check(!strcmp(all[0].message, "embedded energy gradient failed") &&
            !strcmp(all[1].message, "embedded energy gradient failed"),
        "both ranks report the failed evaluation");
  check(energy == -9.0 && virial[0] == -3.0,
        "failed evaluation leaves energy and virial outputs untouched");

  shifted_bohr(pos);
  {
    int rc = cpmdc_session_calculate_au(session, kAtoms, pos, kZ, kGoodCellBohr,
                                        &energy, forces, virial);
    record(&mine, rc, cpmdc_last_error());
  }
  share(&mine, all);
  check(all[0].rc == 0 && all[1].rc == 0,
        "ranks stay in step after a one-rank evaluation failure");

  cpmdc_session_destroy(session);
  report("calculate_au agrees a one-rank evaluation failure");
}

#if !defined(CPMDC_HAS_CPMD)
/* On the reference evaluator, a tiny full-rank cell paired with a huge finite
 * coordinate keeps the cell volume, energy, and gradient finite while the
 * stress tensor overflows, producing a non-finite virial on rank 1 alone. */
static void case_arithmetic_nonfinite_virial_one_rank(void) {
  CPMDCSession *session = cpmdc_session_create(g_params, g_params_size);
  double pos[kCoords];
  double cell[9] = {0};
  double energy = -5.0;
  double forces[kCoords];
  double virial[9];
  struct outcome mine;
  struct outcome all[kRanks];
  check(session != NULL, "session_create");
  if (!session)
    return;

  memcpy(pos, kWaterBohr, sizeof(pos));
  if (g_rank == 1) {
    pos[0] = 1.0e150;
    cell[0] = 1.0e-10;
    cell[4] = 1.0e-10;
    cell[8] = 1.0e-10;
  } else {
    memcpy(cell, kGoodCellBohr, sizeof(cell));
  }
  {
    int rc = cpmdc_session_calculate_au(session, kAtoms, pos, kZ, cell, &energy,
                                        forces, virial);
    record(&mine, rc, cpmdc_last_error());
  }
  share(&mine, all);
  check(all[0].rc == -1 && all[1].rc == -1,
        "a real non-finite virial on one rank fails both ranks");
  check(!strcmp(all[0].message, "virial is not finite") &&
            !strcmp(all[1].message, "virial is not finite"),
        "both ranks report virial is not finite");

  /* Rank 0 repeats the exact geometry it cached on the failed step while
   * rank 1 switches to it; rank_consensus_all keeps both ranks in step. */
  {
    int rc = cpmdc_session_calculate_au(session, kAtoms, kWaterBohr, kZ,
                                        kGoodCellBohr, &energy, forces, virial);
    record(&mine, rc, cpmdc_last_error());
  }
  share(&mine, all);
  check(all[0].rc == 0 && all[1].rc == 0,
        "ranks stay in step even when only one rank matches its cached step");
  check(isfinite(energy) && isfinite(virial[0]),
        "follow-up step produces finite energy and virial");

  cpmdc_session_destroy(session);
  report("a real non-finite virial on one rank fails both ranks cleanly");
}
#endif

int main(int argc, char **argv) {
  MPI_Comm world = MPI_COMM_WORLD;
  int size = 0;
  int any = 0;
  double nan_cell_bohr[9];
  double inf_cell_bohr[9];
  double coplanar_cell_bohr[9] = {18.9, 0.0,  0.0,  0.0, 18.9,
                                  0.0,  18.9, 18.9, 0.0};
  double underflow_cell_bohr[9] = {1.0e-110, 0.0, 0.0, 0.0,      1.0e-110,
                                   0.0,      0.0, 0.0, 1.0e-110};
  double overflow_cell_bohr[9] = {DBL_MAX / 2.0, 0.0, 0.0, 0.0,
                                  DBL_MAX / 2.0, 0.0, 0.0, 0.0,
                                  DBL_MAX / 2.0};
  double nan_cell_ang[9];
  double coplanar_cell_ang[9] = {10.0, 0.0,  0.0,  0.0, 10.0,
                                 0.0,  10.0, 10.0, 0.0};
  double underflow_cell_ang[9] = {1.0e-110, 0.0, 0.0, 0.0,      1.0e-110,
                                  0.0,      0.0, 0.0, 1.0e-110};

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &g_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  alarm(kAlarmSeconds);
  if (size != kRanks || read_params(argc > 1 ? argv[1] : NULL) != 0) {
    if (g_rank == 0)
      fprintf(stderr, "needs 2 ranks and a periodic CPMDParams file\n");
    MPI_Finalize();
    return 2;
  }
  check(cpmdc_adopt_calculator_comm(&world, sizeof(world), kRanks) == 0,
        "adopt the world as one calculator");
  unsetenv("CPMDC_RANK_CONSENSUS");
  unsetenv("CPMDC_ORBITAL_SLOT_LIMIT");

  memcpy(nan_cell_bohr, kGoodCellBohr, sizeof(nan_cell_bohr));
  nan_cell_bohr[4] = NAN;
  memcpy(inf_cell_bohr, kGoodCellBohr, sizeof(inf_cell_bohr));
  inf_cell_bohr[8] = INFINITY;
  memcpy(nan_cell_ang, kGoodCellAng, sizeof(nan_cell_ang));
  nan_cell_ang[0] = NAN;

  check_au_bad_cell_one_rank(1, nan_cell_bohr, 0,
                             "cell contains a non-finite component",
                             "calculate_au agrees a one-rank NaN cell");
  check_au_bad_cell_one_rank(0, inf_cell_bohr, 0,
                             "cell contains a non-finite component",
                             "calculate_au agrees a one-rank infinite cell");
  check_au_bad_cell_one_rank(
      0, coplanar_cell_bohr, 0, "cell does not have three independent vectors",
      "calculate_au agrees a one-rank rank-deficient cell");
  check_au_bad_cell_one_rank(1, underflow_cell_bohr, 0,
                             "cell volume is too small for binary64",
                             "calculate_au agrees a one-rank underflow cell");
  check_au_bad_cell_one_rank(1, overflow_cell_bohr, 0,
                             "cell volume is too large for binary64",
                             "calculate_au agrees a one-rank overflow cell");
  check_au_bad_cell_one_rank(
      1, overflow_cell_bohr, 1, "cell volume is too large for binary64",
      "calculate_au agrees an overflow cell when peer skips virial");

  check_capnp_bad_cell_one_rank(1, nan_cell_ang,
                                "cell contains a non-finite component",
                                "ForceInput agrees a one-rank NaN cell");
  check_capnp_bad_cell_one_rank(
      0, coplanar_cell_ang, "cell does not have three independent vectors",
      "ForceInput agrees a one-rank rank-deficient cell");
  check_capnp_bad_cell_one_rank(1, underflow_cell_ang,
                                "cell volume is too small for binary64",
                                "ForceInput agrees a one-rank underflow cell");

  case_fault_nonfinite_virial_one_rank(
      0, "calculate_au agrees a one-rank post-eval virial failure");
  case_fault_nonfinite_virial_one_rank(
      1, "calculate_au agrees a one-rank virial failure with peer virial NULL");
  case_fault_evaluation_one_rank();
#if !defined(CPMDC_HAS_CPMD)
  case_arithmetic_nonfinite_virial_one_rank();
#endif

  free(g_params);
  MPI_Allreduce(&g_fail, &any, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Finalize();
  return any ? 1 : 0;
}
