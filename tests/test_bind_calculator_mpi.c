/*
 * cpmdc_bind_calculator under four MPI ranks.
 *
 * The program does not call MPI_Init before the first bind, because the
 * header says an OpenCPMD build starts MPI inside the call and no other
 * build does.
 *
 * OpenCPMD build (CPMDC_HAS_CPMD): the world size must be a multiple of
 * ranks_per_calc, a refused value splits nothing, the first accepted value
 * splits once, and a later call returns world_rank / ranks_per_calc without
 * splitting again.
 *
 * Reference build and link stub: every value returns -1, MPI is never
 * started by the library, and nothing is stored.
 */
#include "cpmdc.h"

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int g_world = 0;
static int g_fail = 0;

static void check(int cond, const char *msg) {
  if (cond)
    return;
  fprintf(stderr, "rank %d FAIL: %s\n", g_world, msg);
  g_fail = 1;
}

static int mpi_running(void) {
  int flag = 0;
  MPI_Initialized(&flag);
  return flag;
}

static int stored_size(void) {
  MPI_Comm comm = MPI_COMM_NULL;
  int size = -1;
  if (cpmdc_adopted_comm(&comm, sizeof(comm)) != 0 || comm == MPI_COMM_NULL)
    return -1;
  MPI_Comm_size(comm, &size);
  return size;
}

static void mark(const char *name) {
  int any = 0;
  MPI_Allreduce(&g_fail, &any, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  if (g_world == 0)
    printf("%s %s\n", any ? "FAIL" : "ok", name);
}

#if defined(CPMDC_HAS_CPMD)

static void checkf(int cond, const char *fmt, int a, int b) {
  char msg[160];
  if (cond)
    return;
  snprintf(msg, sizeof(msg), fmt, a, b);
  check(0, msg);
}

static int stored_rank(void) {
  MPI_Comm comm = MPI_COMM_NULL;
  int rank = -1;
  if (cpmdc_adopted_comm(&comm, sizeof(comm)) != 0 || comm == MPI_COMM_NULL)
    return -1;
  MPI_Comm_rank(comm, &rank);
  return rank;
}

static void run_opencpmd(void) {
  int world_size = 0;
  int index;
  const char *err;
  MPI_Comm world_split = MPI_COMM_NULL;
  MPI_Comm first = MPI_COMM_NULL;
  MPI_Comm now = MPI_COMM_NULL;
  int cmp = MPI_UNEQUAL;
  int gathered[4];
  int i;

  check(!mpi_running(), "MPI was running before the first bind");

  /* A value that does not divide the world size. The call has already
   * started MPI, then refuses the value and splits nothing.
   */
  index = cpmdc_bind_calculator(3);
  check(index == -1, "bind(3) with 4 ranks was not refused");
  check(mpi_running(), "bind did not start MPI");
  err = cpmdc_last_error();
  check(err && strstr(err, "refused"), "bind(3) last_error");
  MPI_Comm_rank(MPI_COMM_WORLD, &g_world);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  check(world_size == 4, "this test needs four ranks");
  check(stored_size() == -1, "a refused bind stored a communicator");
  mark("a value that does not divide the world size is refused");

  check(cpmdc_bind_calculator(8) == -1, "bind(8) above the world size");
  check(cpmdc_bind_calculator(5) == -1, "bind(5) above the world size");
  check(stored_size() == -1, "a too large bind stored a communicator");
  mark("a value above the world size is refused");

  /* After the refusals a valid value still splits. */
  index = cpmdc_bind_calculator(2);
  checkf(index == g_world / 2, "bind(2) index %d, expected rank/2 for rank %d",
         index, g_world);
  checkf(stored_size() == 2, "group size %d, expected 2 (%d)", stored_size(),
         0);
  checkf(stored_rank() == g_world % 2, "rank in group %d, expected %d",
         stored_rank(), g_world % 2);
  cpmdc_adopted_comm(&first, sizeof(first));
  {
    int sum = 0;
    int mine = g_world;
    MPI_Allreduce(&mine, &sum, 1, MPI_INT, MPI_SUM, first);
    checkf(sum == 2 * (g_world / 2) * 2 + 1, "group rank sum %d, expected %d",
           sum, 2 * (g_world / 2) * 2 + 1);
  }
  mark("bind(2) makes blocks of two consecutive ranks");

  /* Later calls check their own value and do not split again. */
  index = cpmdc_bind_calculator(4);
  checkf(index == g_world / 4, "bind(4) after bind(2): %d, expected %d", index,
         g_world / 4);
  cpmdc_adopted_comm(&now, sizeof(now));
  MPI_Comm_compare(first, now, &cmp);
  check(cmp == MPI_IDENT, "bind(4) replaced the communicator");
  checkf(stored_size() == 2, "size after bind(4) is %d, expected 2 (%d)",
         stored_size(), 0);

  index = cpmdc_bind_calculator(3);
  check(index == -1, "bind(3) after a split was not refused");
  index = cpmdc_bind_calculator(8);
  check(index == -1, "bind(8) after a split was not refused");
  checkf(stored_size() == 2, "size after refusals is %d, expected 2 (%d)",
         stored_size(), 0);

  index = cpmdc_bind_calculator(0);
  checkf(index == 0, "bind(0): %d, expected 0 (%d)", index, 0);
  index = cpmdc_bind_calculator(1);
  checkf(index == g_world, "bind(1): %d, expected %d", index, g_world);
  cpmdc_adopted_comm(&now, sizeof(now));
  MPI_Comm_compare(first, now, &cmp);
  check(cmp == MPI_IDENT, "a later bind replaced the communicator");
  mark("a later bind does not split again");

  /* A communicator the host split replaces the stored one, and a bind
   * after it does not split either.
   */
  MPI_Comm_split(MPI_COMM_WORLD, g_world / 4, g_world, &world_split);
  index = cpmdc_adopt_calculator_comm(&world_split, sizeof(world_split), 4);
  checkf(index == 0, "adopt of the whole world: %d, expected 0 (%d)", index, 0);
  checkf(stored_size() == 4, "size after adopt is %d, expected 4 (%d)",
         stored_size(), 0);
  index = cpmdc_bind_calculator(2);
  checkf(index == g_world / 2, "bind(2) after adopt: %d, expected %d", index,
         g_world / 2);
  checkf(stored_size() == 4, "bind after adopt split to size %d, expected 4 (%d)",
         stored_size(), 0);
  mark("a bind after an adopt does not split");

  index = cpmdc_bind_calculator(0);
  MPI_Allgather(&index, 1, MPI_INT, gathered, 1, MPI_INT, MPI_COMM_WORLD);
  for (i = 0; i < 4; i++)
    check(gathered[i] == 0, "bind(0) index differs between ranks");
  MPI_Comm_free(&world_split);
}

#else

static void run_without_backend(void) {
  int values[] = {-1, 0, 1, 2, 4, 8};
  size_t i;
  int started;
  const char *err;

  check(!mpi_running(), "MPI was running before the first bind");
  for (i = 0; i < sizeof(values) / sizeof(values[0]); i++) {
    char msg[96];
    snprintf(msg, sizeof(msg), "bind(%d) did not return -1", values[i]);
    check(cpmdc_bind_calculator(values[i]) == -1, msg);
  }
  err = cpmdc_last_error();
  check(err && err[0], "bind without a backend left no last_error");
  check(!mpi_running(), "bind started MPI in a build without OpenCPMD");
  started = mpi_running();

  MPI_Init(NULL, NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &g_world);
  mark(started ? "bind started MPI" : "bind returns -1 and does not start MPI");
  check(cpmdc_bind_calculator(2) == -1, "bind(2) with MPI running");
  check(stored_size() == -1, "bind stored a communicator without a backend");
  mark("bind with MPI running returns -1 and stores nothing");
}

#endif

int main(void) {
  int any = 0;

#if defined(CPMDC_HAS_CPMD)
  run_opencpmd();
#else
  run_without_backend();
#endif
  MPI_Allreduce(&g_fail, &any, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Finalize();
  return any ? 1 : 0;
}
