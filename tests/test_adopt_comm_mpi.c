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

static void mark(const char *name) {
  int any = 0;
  MPI_Allreduce(&g_fail, &any, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  if (g_world == 0)
    printf("%s %s\n", any ? "FAIL" : "ok", name);
}

static MPI_Comm split_color(int color, int key) {
  MPI_Comm comm = MPI_COMM_NULL;
  MPI_Comm_split(MPI_COMM_WORLD, color, key, &comm);
  return comm;
}

static int root_world_rank(MPI_Comm comm) {
  int local = 0;
  int root = -1;
  MPI_Comm_rank(comm, &local);
  root = (local == 0) ? g_world : -1;
  MPI_Allreduce(MPI_IN_PLACE, &root, 1, MPI_INT, MPI_MAX, comm);
  return root;
}

static void expect_distinct_roots(MPI_Comm comm, int index, const char *label);

static void expect_group(MPI_Comm comm, int ranks_per_calc, int blocked_index,
                         const char *label) {
  int index;
  int size = 0;
  int local = 0;
  int root;
  int expect;
  int gathered[4];
  int i;
  int cmp = MPI_UNEQUAL;
  MPI_Comm stored = MPI_COMM_NULL;
  char msg[160];

  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &local);
  root = root_world_rank(comm);
  expect = root / size;
  index = cpmdc_adopt_calculator_comm(&comm, sizeof(comm), ranks_per_calc);
  snprintf(msg, sizeof(msg),
           "%s: index %d, root world %d, size %d, expect %d, error %s", label,
           index, root, size, expect,
           cpmdc_last_error() ? cpmdc_last_error() : "");
  check(index == expect, msg);
  if (blocked_index >= 0) {
    snprintf(msg, sizeof(msg), "%s: blocked index %d, got %d", label,
             blocked_index, index);
    check(index == blocked_index, msg);
  }
  check(size <= 4, "group larger than the gather buffer");
  if (size > 0 && size <= 4) {
    MPI_Allgather(&index, 1, MPI_INT, gathered, 1, MPI_INT, comm);
    for (i = 0; i < size; i++) {
      snprintf(msg, sizeof(msg), "%s: rank %d of group returned %d", label, i,
               gathered[i]);
      check(gathered[i] == index, msg);
    }
  }
  check(cpmdc_adopted_comm(&stored, sizeof(stored)) == 0, "adopted_comm");
  if (stored != MPI_COMM_NULL) {
    int stored_size = -1;
    MPI_Comm_size(stored, &stored_size);
    snprintf(msg, sizeof(msg), "%s: stored size %d, group size %d", label,
             stored_size, size);
    check(stored_size == size, msg);
    MPI_Comm_compare(comm, stored, &cmp);
    check(cmp == MPI_IDENT || cmp == MPI_CONGRUENT, "stored comm compare");
  }
  expect_distinct_roots(comm, index, label);
}

/* One index per group root, gathered on MPI_COMM_WORLD. The roots' indices
 * are 0 .. ncalc-1 with no duplicates.
 */
static void expect_distinct_roots(MPI_Comm comm, int index, const char *label) {
  int local = 0;
  int size = 0;
  int world_size = 0;
  int from_root;
  int gathered[4];
  int roots[4];
  int nroots = 0;
  int ncalc;
  int i;
  int j;
  char msg[160];

  MPI_Comm_rank(comm, &local);
  MPI_Comm_size(comm, &size);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  from_root = (local == 0 && index >= 0) ? index : -1;
  check(world_size == 4 && world_size <= 4, "world size for the root gather");
  MPI_Allgather(&from_root, 1, MPI_INT, gathered, 1, MPI_INT, MPI_COMM_WORLD);
  for (i = 0; i < world_size && i < 4; i++) {
    if (gathered[i] >= 0)
      roots[nroots++] = gathered[i];
  }
  ncalc = size > 0 ? world_size / size : 0;
  snprintf(msg, sizeof(msg), "%s: %d group roots, expected %d", label, nroots,
           ncalc);
  check(nroots == ncalc, msg);
  for (i = 0; i < nroots; i++) {
    for (j = i + 1; j < nroots; j++) {
      if (roots[j] < roots[i]) {
        int tmp = roots[i];
        roots[i] = roots[j];
        roots[j] = tmp;
      }
    }
  }
  for (i = 0; i < nroots && i < ncalc; i++) {
    snprintf(msg, sizeof(msg), "%s: root index %d, expected %d", label,
             roots[i], i);
    check(roots[i] == i, msg);
  }
}

static void expect_mismatch(MPI_Comm comm, int ranks_per_calc) {
  int size = 0;
  int index;
  const char *err;
  char needle[64];
  MPI_Comm_size(comm, &size);
  index = cpmdc_adopt_calculator_comm(&comm, sizeof(comm), ranks_per_calc);
  check(index < 0, "size mismatch returned an index");
  err = cpmdc_last_error();
  check(err != NULL && strstr(err, "does not equal") != NULL,
        "size mismatch last_error");
  snprintf(needle, sizeof(needle), "size %d", size);
  check(err != NULL && strstr(err, needle) != NULL, "error names the size");
  snprintf(needle, sizeof(needle), "ranks_per_calc %d", ranks_per_calc);
  check(err != NULL && strstr(err, needle) != NULL,
        "error names ranks_per_calc");
}

static unsigned char *read_file(const char *path, size_t *nbytes) {
  FILE *file;
  long length;
  unsigned char *buf;
  file = fopen(path, "rb");
  if (!file)
    return NULL;
  if (fseek(file, 0, SEEK_END) != 0) {
    fclose(file);
    return NULL;
  }
  length = ftell(file);
  if (length < 0) {
    fclose(file);
    return NULL;
  }
  if (fseek(file, 0, SEEK_SET) != 0) {
    fclose(file);
    return NULL;
  }
  buf = (unsigned char *)malloc((size_t)length);
  if (!buf) {
    fclose(file);
    return NULL;
  }
  if (length > 0 && fread(buf, 1, (size_t)length, file) != (size_t)length) {
    free(buf);
    fclose(file);
    return NULL;
  }
  fclose(file);
  *nbytes = (size_t)length;
  return buf;
}

int main(int argc, char **argv) {
  int world_size = 0;
  int null_rc;
  const char *err;
  MPI_Comm pair;
  MPI_Comm reversed;
  MPI_Comm single;
  MPI_Comm striped;
  MPI_Comm uneven;
  MPI_Comm world = MPI_COMM_NULL;
  int stored_size = -1;
  MPI_Comm stored = MPI_COMM_NULL;
  int have_session = 0;
  int all_have = 0;
  CPMDCSession *session = NULL;
  const char *params_path = argc > 1 ? argv[1] : NULL;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &g_world);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  if (world_size != 4) {
    if (g_world == 0)
      fprintf(stderr, "adopt-comm test requires 4 ranks, found %d\n",
              world_size);
    MPI_Finalize();
    return 1;
  }

  null_rc = cpmdc_adopt_calculator_comm(NULL, 0, 2);
  check(null_rc < 0, "null communicator returned an index");
  err = cpmdc_last_error();
#if defined(CPMDC_TEST_STUB)
  check(err && strstr(err, "stub") && strstr(err, "cpmdc_adopt_calculator_comm"),
        "stub null error");
#else
  check(err && strstr(err, "calculator bind refused"), "null error text");
#endif
  mark("null communicator");

  pair = split_color(g_world / 2, g_world % 2);
  expect_mismatch(pair, 4);
  check(cpmdc_adopted_comm(&stored, sizeof(stored)) < 0,
        "mismatch stored a communicator");
  mark("2x2 size mismatch before a successful adopt");

  expect_group(pair, 2, g_world / 2, "2x2");
  mark("2x2 blocked");

  stored = MPI_COMM_NULL;
  check(cpmdc_adopted_comm(&stored, sizeof(stored)) == 0, "stored after 2x2");
  stored_size = -1;
  if (stored != MPI_COMM_NULL)
    MPI_Comm_size(stored, &stored_size);
  expect_mismatch(pair, 1);
  {
    MPI_Comm still = MPI_COMM_NULL;
    int still_size = -1;
    check(cpmdc_adopted_comm(&still, sizeof(still)) == 0,
          "mismatch cleared the stored communicator");
    if (still != MPI_COMM_NULL)
      MPI_Comm_size(still, &still_size);
    check(still_size == stored_size, "mismatch replaced the stored size");
  }
  mark("2x2 size mismatch leaves the previous adopt");

  expect_group(pair, 0, g_world / 2, "2x2 ranks_per_calc 0");
  mark("2x2 with ranks_per_calc <= 0");

  world = MPI_COMM_WORLD;
  expect_group(world, 4, 0, "1x4");
  mark("1x4");

  single = split_color(g_world, 0);
  expect_group(single, 1, g_world, "4x1");
  mark("4x1");

  reversed = split_color(g_world / 2, 100 - g_world);
  expect_group(reversed, 2, g_world / 2, "2x2 reversed key");
  mark("2x2 reversed key");

  /* World ranks 0 and 2 share a group. Both striped groups would receive
   * index 0, so the call is refused and the stored communicator stays.
   */
  striped = split_color(g_world % 2, g_world);
  {
    MPI_Comm before = MPI_COMM_NULL;
    MPI_Comm after = MPI_COMM_NULL;
    int before_size = -1;
    int after_size = -1;
    int cmp = MPI_UNEQUAL;
    int index;
    check(cpmdc_adopted_comm(&before, sizeof(before)) == 0,
          "comm before the striped adopt");
    if (before != MPI_COMM_NULL)
      MPI_Comm_size(before, &before_size);
    index = cpmdc_adopt_calculator_comm(&striped, sizeof(striped), 2);
    check(index < 0, "striped split returned an index");
    err = cpmdc_last_error();
    check(err && strstr(err, "aligned block"), "striped last_error");
    check(cpmdc_adopted_comm(&after, sizeof(after)) == 0,
          "striped cleared the stored communicator");
    if (after != MPI_COMM_NULL)
      MPI_Comm_size(after, &after_size);
    check(after_size == before_size, "striped replaced the stored size");
    if (before != MPI_COMM_NULL && after != MPI_COMM_NULL) {
      MPI_Comm_compare(before, after, &cmp);
      check(cmp == MPI_IDENT || cmp == MPI_CONGRUENT,
            "striped replaced the stored communicator");
    }
  }
  mark("striped placement");

  /* Colors 0, 1, 1, 2. Group {1, 2} has size 2 and starts at an odd rank.
   * The singletons are aligned blocks and receive their world ranks.
   */
  uneven = split_color(g_world == 0 ? 0 : (g_world == 3 ? 2 : 1), g_world);
  {
    int group_size = 0;
    int index;
    MPI_Comm_size(uneven, &group_size);
    index = cpmdc_adopt_calculator_comm(&uneven, sizeof(uneven), group_size);
    if (g_world == 1 || g_world == 2) {
      MPI_Comm after = MPI_COMM_NULL;
      int cmp = MPI_UNEQUAL;
      check(index < 0, "misaligned block returned an index");
      err = cpmdc_last_error();
      check(err && strstr(err, "aligned block"), "misaligned last_error");
      check(cpmdc_adopted_comm(&after, sizeof(after)) == 0,
            "misaligned cleared the stored communicator");
      if (after != MPI_COMM_NULL) {
        MPI_Comm_compare(reversed, after, &cmp);
        check(cmp == MPI_IDENT || cmp == MPI_CONGRUENT,
              "misaligned replaced the stored communicator");
      }
    } else if (g_world == 0) {
      check(index == 0, "singleton world rank 0 index");
    } else {
      check(index == 3, "singleton world rank 3 index");
    }
  }
  mark("misaligned block");

  if (params_path) {
    size_t nbytes = 0;
    unsigned char *bytes = read_file(params_path, &nbytes);
    if (bytes && nbytes > 0)
      session = cpmdc_session_create(bytes, nbytes);
    free(bytes);
    have_session = session != NULL;
  }
  MPI_Allreduce(&have_session, &all_have, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  if (all_have) {
    int before = -1;
    int refused;
    MPI_Comm before_comm = MPI_COMM_NULL;
    check(cpmdc_adopted_comm(&before_comm, sizeof(before_comm)) == 0,
          "comm before the live-session adopt");
    if (before_comm != MPI_COMM_NULL)
      MPI_Comm_size(before_comm, &before);
    refused = cpmdc_adopt_calculator_comm(&world, sizeof(world), 4);
    check(refused < 0, "live session allowed a rebind");
    err = cpmdc_last_error();
    check(err && strstr(err, "session is live"), "live session last_error");
    {
      MPI_Comm after = MPI_COMM_NULL;
      int after_size = -1;
      check(cpmdc_adopted_comm(&after, sizeof(after)) == 0,
            "live session cleared the communicator");
      if (after != MPI_COMM_NULL)
        MPI_Comm_size(after, &after_size);
      check(after_size == before, "live session replaced the communicator");
    }
  }
  if (session)
    cpmdc_session_destroy(session);
  if (all_have) {
    expect_group(world, 4, 0, "rebind after session destroy");
    mark("rebind with no live session");
  } else if (g_world == 0 && params_path) {
    printf("skip live-session rebind: session_create failed\n");
  }

  MPI_Comm_free(&pair);
  MPI_Comm_free(&reversed);
  MPI_Comm_free(&single);
  MPI_Comm_free(&striped);
  MPI_Comm_free(&uneven);
  {
    int any = 0;
    MPI_Allreduce(&g_fail, &any, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    MPI_Finalize();
    return any ? 1 : 0;
  }
}
