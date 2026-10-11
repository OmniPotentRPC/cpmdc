#include "cpmdc_adopt.h"

#include <stdio.h>
#include <string.h>

#if defined(CPMDC_HAVE_MPI)
#include <mpi.h>

static int g_have = 0;
static MPI_Comm g_comm;

static void set_err(char *err, size_t errcap, const char *msg) {
  if (!err || errcap == 0)
    return;
  snprintf(err, errcap, "%s", msg ? msg : "");
}

int cpmdc_adopt_resolve(const void *comm, size_t comm_bytes, int ranks_per_calc,
                        int *index_out, char *err, size_t errcap) {
  MPI_Comm ccomm;
  int inited = 0;
  int rank = 0;
  int size = 0;
  int world_rank = 0;
  int index = -1;

  if (!index_out) {
    set_err(err, errcap, "CPMD calculator bind refused");
    return -1;
  }
  if (!comm || comm_bytes != sizeof(MPI_Comm)) {
    set_err(err, errcap, "CPMD calculator bind refused");
    return -1;
  }
  memcpy(&ccomm, comm, sizeof(ccomm));
  if (ccomm == MPI_COMM_NULL) {
    set_err(err, errcap, "CPMD calculator bind refused");
    return -1;
  }
  if (MPI_Initialized(&inited) != MPI_SUCCESS || !inited) {
    set_err(err, errcap, "CPMD calculator bind refused");
    return -1;
  }
  if (MPI_Comm_rank(ccomm, &rank) != MPI_SUCCESS ||
      MPI_Comm_size(ccomm, &size) != MPI_SUCCESS ||
      MPI_Comm_rank(MPI_COMM_WORLD, &world_rank) != MPI_SUCCESS || size <= 0) {
    set_err(err, errcap, "CPMD calculator bind refused");
    return -1;
  }
  /* Every rank of comm sees the same size, so a mismatch returns before
   * the broadcast and does not leave the group waiting.
   */
  if (ranks_per_calc > 0 && size != ranks_per_calc) {
    if (err && errcap > 0) {
      snprintf(err, errcap,
               "adopted communicator size %d does not equal ranks_per_calc %d",
               size, ranks_per_calc);
    }
    return -1;
  }
  if (rank == 0)
    index = world_rank / size;
  if (MPI_Bcast(&index, 1, MPI_INT, 0, ccomm) != MPI_SUCCESS || index < 0) {
    set_err(err, errcap, "CPMD calculator bind refused");
    return -1;
  }
  *index_out = index;
  set_err(err, errcap, "");
  return 0;
}

int cpmdc_adopt_commit(const void *comm, size_t comm_bytes) {
  if (!comm || comm_bytes != sizeof(MPI_Comm))
    return -1;
  memcpy(&g_comm, comm, sizeof(g_comm));
  g_have = 1;
  return 0;
}

int cpmdc_adopt_copy(void *out, size_t nbytes) {
  if (!g_have || !out || nbytes != sizeof(MPI_Comm))
    return -1;
  memcpy(out, &g_comm, sizeof(g_comm));
  return 0;
}

#else

int cpmdc_adopt_resolve(const void *comm, size_t comm_bytes, int ranks_per_calc,
                        int *index_out, char *err, size_t errcap) {
  (void)comm;
  (void)comm_bytes;
  (void)ranks_per_calc;
  (void)index_out;
  if (err && errcap > 0)
    snprintf(err, errcap, "%s", "CPMD calculator bind refused");
  return -1;
}

int cpmdc_adopt_commit(const void *comm, size_t comm_bytes) {
  (void)comm;
  (void)comm_bytes;
  return -1;
}

int cpmdc_adopt_copy(void *out, size_t nbytes) {
  (void)out;
  (void)nbytes;
  return -1;
}

#endif
