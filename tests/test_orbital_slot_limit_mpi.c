/**
 * Multi-rank MPI calculator consensus and bounded eviction for
 * CPMDC_ORBITAL_SLOT_LIMIT.
 */
#define _POSIX_C_SOURCE 200809L
#include "cpmdc.h"

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int cpmdc_embed_orbital_slots(void);

enum { kAtoms = 3, kCoords = 3 * kAtoms };

static const int kZ[kAtoms] = {8, 1, 1};
static const double kG0[kCoords] = {0.0, 0.0,     0.1173,  0.0,    0.7572,
                                    -0.4692, 0.0, -0.7572, -0.4692};
static const double kG1[kCoords] = {0.0, 0.0,     0.1173,  0.0,    0.9400,
                                    -0.5400, 0.0, -0.7572, -0.4692};
static const double kG2[kCoords] = {0.0, 0.0,     0.1173,  0.0,    1.1300,
                                    -0.6150, 0.0, -0.7572, -0.4692};

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
    /* 1) Invalid value on rank 1 fails every rank without deadlocking. */
    if (rank == 0)
      setenv("CPMDC_ORBITAL_SLOT_LIMIT", "1", 1);
    else
      setenv("CPMDC_ORBITAL_SLOT_LIMIT", "0", 1);
    (void)cpmdc_session_select_orbitals(session, 0);
    CPMDCResult r_bad =
        cpmdc_session_energy_gradient(session, kAtoms, kG0, kZ, grad);
    if (r_bad.ok != 0) {
      fprintf(stderr, "rank %d: expected failure on invalid rank limit\n",
              rank);
      fail = 1;
    }

    if (cpmdc_available()) {
      /* 2) Conflicting positive limits across ranks fail collectively. */
      setenv("CPMDC_ORBITAL_SLOT_LIMIT", rank == 0 ? "1" : "2", 1);
      CPMDCResult r_conflict =
          cpmdc_session_energy_gradient(session, kAtoms, kG0, kZ, grad);
      if (r_conflict.ok != 0 ||
          strstr(r_conflict.message, "CPMDC_ORBITAL_SLOT_LIMIT") == NULL) {
        fprintf(stderr,
                "rank %d: expected conflict failure, got ok=%d msg=%s\n", rank,
                r_conflict.ok, r_conflict.message);
        fail = 1;
      }

      /* 3) Rank 0 limit propagates to unset rank 1 so both cap slots at 1. */
      if (rank == 0)
        setenv("CPMDC_ORBITAL_SLOT_LIMIT", "1", 1);
      else
        unsetenv("CPMDC_ORBITAL_SLOT_LIMIT");
      const double *steps[3] = {kG0, kG1, kG2};
      for (int i = 0; i < 3; ++i) {
        (void)cpmdc_session_select_orbitals(session, (long long)i);
        CPMDCResult r =
            cpmdc_session_energy_gradient(session, kAtoms, steps[i], kZ, grad);
        if (!r.ok) {
          fprintf(stderr, "rank %d step %d failed: %s\n", rank, i, r.message);
          fail = 1;
          break;
        }
      }
      int slots = cpmdc_embed_orbital_slots();
      if (slots != 1) {
        fprintf(stderr, "rank %d: expected 1 parked slot, got %d\n", rank,
                slots);
        fail = 1;
      }
    } else {
      setenv("CPMDC_ORBITAL_SLOT_LIMIT", "1", 1);
      CPMDCResult r_ok =
          cpmdc_session_energy_gradient(session, kAtoms, kG0, kZ, grad);
      if (!r_ok.ok)
        fail = 1;
    }
    cpmdc_session_destroy(session);
  }

  free(params);
  MPI_Allreduce(&fail, &any_fail, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  if (rank == 0)
    printf("%s orbital-slot-limit-mpi\n", any_fail ? "FAIL" : "ok");
  MPI_Finalize();
  return any_fail ? 1 : 0;
}
