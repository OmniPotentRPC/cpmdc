#define _POSIX_C_SOURCE 200809L
#include "cpmdc.h"

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

enum { kAtoms = 3, kCoords = 3 * kAtoms };

static const int kZ[kAtoms] = {8, 1, 1};
static const double kGeom[kCoords] = {0.0, 0.0,     0.1173,  0.0,    0.7572,
                                      -0.4692, 0.0, -0.7572, -0.4692};

static unsigned char *read_file(const char *path, size_t *size) {
  FILE *fp = path ? fopen(path, "rb") : NULL;
  long n;
  unsigned char *buf;
  if (!fp)
    return NULL;
  if (fseek(fp, 0, SEEK_END) != 0 || (n = ftell(fp)) <= 0 ||
      fseek(fp, 0, SEEK_SET) != 0) {
    fclose(fp);
    return NULL;
  }
  buf = (unsigned char *)malloc((size_t)n);
  if (!buf || fread(buf, 1, (size_t)n, fp) != (size_t)n) {
    free(buf);
    fclose(fp);
    return NULL;
  }
  fclose(fp);
  *size = (size_t)n;
  return buf;
}

int main(int argc, char **argv) {
  MPI_Comm world = MPI_COMM_WORLD;
  int rank = 0;
  int size = 0;
  int fail = 0;
  int any_fail = 0;
  size_t params_size = 0;
  unsigned char *params = NULL;
  CPMDCSession *session = NULL;
  double grad[kCoords];

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 2 || argc < 2) {
    if (rank == 0)
      fprintf(stderr, "usage: mpiexec -n 2 %s PARAMS_BIN\n", argv[0]);
    MPI_Finalize();
    return 2;
  }
  params = read_file(argv[1], &params_size);
  if (!params) {
    MPI_Finalize();
    return 2;
  }
  if (cpmdc_adopt_calculator_comm(&world, sizeof(world), 2) != 0) {
    fprintf(stderr, "rank %d: adopt_calculator_comm failed: %s\n", rank,
            cpmdc_last_error());
    free(params);
    MPI_Finalize();
    return 1;
  }

  setenv("CPMDC_RANK_CONSENSUS", "1", 1);
  session = cpmdc_session_create(params, params_size);
  if (!session) {
    fail = 1;
  } else {
    /* Prime consensus with a common explicit limit. */
    setenv("CPMDC_ORBITAL_SLOT_LIMIT", "1", 1);
    (void)cpmdc_session_select_orbitals(session, 0);
    CPMDCResult r_seed =
        cpmdc_session_energy_gradient(session, kAtoms, kGeom, kZ, grad);
    if (!r_seed.ok) {
      fprintf(stderr, "rank %d: seed call failed: %s\n", rank, r_seed.message);
      fail = 1;
    }

    /* Mixed set/unset must still enter a matched collective on both ranks. */
    if (rank == 0)
      setenv("CPMDC_ORBITAL_SLOT_LIMIT", "1", 1);
    else
      unsetenv("CPMDC_ORBITAL_SLOT_LIMIT");
    (void)cpmdc_session_select_orbitals(session, 1);
    CPMDCResult r_mixed =
        cpmdc_session_energy_gradient(session, kAtoms, kGeom, kZ, grad);
    if (!r_mixed.ok) {
      fprintf(stderr, "rank %d: mixed set/unset call failed: %s\n", rank,
              r_mixed.message);
      fail = 1;
    }
    cpmdc_session_destroy(session);
  }

  free(params);
  MPI_Allreduce(&fail, &any_fail, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  if (rank == 0)
    printf("%s orbital-slot-limit-mpi-sync\n", any_fail ? "FAIL" : "ok");
  MPI_Finalize();
  return any_fail ? 1 : 0;
}
