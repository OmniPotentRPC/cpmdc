#define _POSIX_C_SOURCE 200809L

#include "cpmdc.h"

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define NATOMS 3
#define N3 (NATOMS * 3)

static int g_world = 0;
static int g_fail = 0;
static const int k_z[NATOMS] = {8, 1, 1};

static void check(int cond, const char *msg) {
  if (cond)
    return;
  fprintf(stderr, "rank %d FAIL: %s\n", g_world, msg);
  g_fail = 1;
}

static void mark(const char *name) {
  int any = 0;
  MPI_Allreduce(&g_fail, &any, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  if (g_world == 0)
    printf("%s %s\n", any ? "FAIL" : "ok", name);
}

/* Water with the oxygen moved by a different amount on each rank, so the
 * ranks of one calculator disagree the way diverged host optimizers do. */
static void positions(double *pos, double shift) {
  static const double base[N3] = {0.0, 0.0,     0.1173,  0.0,    0.7572,
                                  -0.4692, 0.0, -0.7572, -0.4692};
  memcpy(pos, base, sizeof(base));
  pos[0] += 1.0e-3 * g_world + shift;
}

static unsigned char *g_params = NULL;
static size_t g_params_size = 0;

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

static CPMDCSession *open_session(void) {
  CPMDCSession *session = cpmdc_session_create(g_params, g_params_size);
  check(session != NULL, "session_create");
  return session;
}

static void eval(CPMDCSession *session, const double *pos, double *energy,
                 double *grad, const char *label) {
  CPMDCResult r = cpmdc_session_energy_gradient(session, NATOMS, pos, k_z, grad);
  check(r.ok, label);
  *energy = r.energy_h;
}

static int same(double e, const double *g, const double *e_ref,
                const double *g_ref) {
  return memcmp(&e, e_ref, sizeof(e)) == 0 &&
         memcmp(g, g_ref, N3 * sizeof(*g)) == 0;
}

int main(int argc, char **argv) {
  MPI_Comm world = MPI_COMM_WORLD;
  MPI_Comm pair = MPI_COMM_NULL;
  CPMDCSession *off = NULL;
  CPMDCSession *on = NULL;
  CPMDCSession *grouped = NULL;
  CPMDCRankSpread spread;
  double pos[N3], e, g[N3], e_all[4], g_all[4 * N3];
  double e_spread = 0.0;
  double g_spread = 0.0;
  int world_size = 0;
  int lead;
  int any = 0;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &g_world);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  if (world_size != 4 || read_params(argc > 1 ? argv[1] : NULL) != 0) {
    if (g_world == 0)
      fprintf(stderr, "needs 4 ranks and a CPMDParams file\n");
    MPI_Finalize();
    return 1;
  }
  check(cpmdc_adopt_calculator_comm(&world, sizeof(world), 4) == 0,
        "adopt the world as one calculator");

  unsetenv("CPMDC_RANK_CONSENSUS");
  off = open_session();
  positions(pos, 0.0);
  eval(off, pos, &e, g, "evaluation without consensus");
  MPI_Allgather(&e, 1, MPI_DOUBLE, e_all, 1, MPI_DOUBLE, world);
  MPI_Allgather(g, N3, MPI_DOUBLE, g_all, N3, MPI_DOUBLE, world);
  check(e_all[1] != e_all[0], "ranks disagree without consensus");
  check(cpmdc_last_rank_spread(&spread) == -1 && !spread.valid,
        "no spread without consensus");
  mark("each rank keeps its own result without consensus");

  setenv("CPMDC_RANK_CONSENSUS", "1", 1);
  on = open_session();
  eval(on, pos, &e, g, "evaluation with consensus");
  check(same(e, g, &e_all[0], g_all), "every rank holds rank 0's result");
  for (int r = 0; r < 4; ++r) {
    e_spread = fmax(e_spread, fabs(e_all[r] - e_all[0]));
    for (int i = 0; i < N3; ++i)
      g_spread = fmax(g_spread, fabs(g_all[r * N3 + i] - g_all[i]));
  }
  check(cpmdc_last_rank_spread(&spread) == 0 && spread.ranks == 4,
        "spread over 4 ranks");
  check(spread.energy_h == e_spread && spread.gradient_h_bohr == g_spread &&
            e_spread > 0.0 && g_spread > 0.0,
        "spread equals the gathered differences");
  mark("consensus gives every rank rank 0's result and the spread");

  /* Rank 0 repeats its stored positions while the others move. */
  if (g_world != 0)
    pos[4] += 1.0e-3;
  eval(on, pos, &e, g, "evaluation after a rank 0 cache hit");
  check(same(e, g, &e_all[0], g_all), "rank 0's stored positions win");
  mark("a cache hit on one rank does not split the calculator");

  positions(pos, 2.0e-3);
  if (g_world == 3)
    pos[0] = NAN;
  eval(on, pos, &e, g, "evaluation with a NaN on rank 3");
  {
    double root[N3 + 1];
    root[0] = e;
    memcpy(root + 1, g, sizeof(g));
    MPI_Bcast(root, N3 + 1, MPI_DOUBLE, 0, world);
    check(same(e, g, &root[0], root + 1) && isfinite(e),
          "rank 3 holds rank 0's finite result");
  }
  check(cpmdc_last_rank_spread(&spread) == 0 && isinf(spread.energy_h) &&
            isinf(spread.gradient_h_bohr),
        "a NaN difference is an infinite spread");
  mark("a NaN on one rank is reported as an infinite spread");

  cpmdc_session_destroy(off);
  cpmdc_session_destroy(on);
  MPI_Comm_split(world, g_world / 2, g_world % 2, &pair);
  check(cpmdc_adopt_calculator_comm(&pair, sizeof(pair), 2) == g_world / 2,
        "adopt two calculators of two ranks");
  grouped = open_session();
  positions(pos, 0.0);
  eval(grouped, pos, &e, g, "evaluation on two calculators");
  lead = (g_world / 2) * 2;
  check(same(e, g, &e_all[lead], g_all + lead * N3),
        "each calculator holds its own rank 0's result");
  check(cpmdc_last_rank_spread(&spread) == 0 && spread.ranks == 2 &&
            spread.energy_h == fabs(e_all[lead + 1] - e_all[lead]),
        "spread over the calculator, not the world");
  mark("consensus runs on the calculator communicator");

  unsetenv("CPMDC_RANK_CONSENSUS");
  positions(pos, 3.0e-3);
  eval(grouped, pos, &e, g, "evaluation after consensus is turned off");
  check(cpmdc_last_rank_spread(&spread) == -1,
        "an evaluation without consensus clears the spread");
  mark("turning consensus off clears the spread");

  cpmdc_session_destroy(grouped);
  MPI_Comm_free(&pair);
  free(g_params);
  MPI_Allreduce(&g_fail, &any, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Finalize();
  return any ? 1 : 0;
}
