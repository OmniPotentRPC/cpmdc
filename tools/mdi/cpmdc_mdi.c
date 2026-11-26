/**
 * Command state for LAMMPS fix mdi/qm (stable 22 Jul 2025, src/MDI/fix_mdi_qm.cpp).
 *
 * LAMMPS sends, in Bohr: >NATOMS, then >CELL_DISPL and >CELL when the box is
 * fully periodic, then >ELEMENTS (atomic numbers, from the elements keyword)
 * or >TYPES, then >COORDS on every step. Coordinates are absolute Cartesian
 * values, unwrapped, and >CELL_DISPL is the box origin. This engine subtracts
 * that origin so the session sees a cell whose origin is zero.
 *
 * LAMMPS then reads <PE and <KE_ELEC when both exist, otherwise <ENERGY, then
 * <FORCES, then <STRESS when virial is requested. It adds <PE and <KE_ELEC.
 * The cpmdc energy already includes the electronic kinetic energy, so
 * <KE_ELEC is 0 and <PE equals <ENERGY. <STRESS is nine Hartree/Bohr^3
 * values, row-major xx xy xz yx yy yz zx zy zz, which is
 * cpmdc_session_last_stress(). LAMMPS symmetrizes that layout, multiplies by
 * the box volume, and stores the extensive virial. A step with no stress
 * sends zeros. One session call serves every read until the next '>' command.
 * >NATOMS starts a system: the previous cell, origin, species, and coordinates
 * are dropped, and LAMMPS sends the new ones before it reads.
 */
#include "cpmdc_mdi.h"

#include "cpmdc.h"

#include <limits.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct cpmdc_mdi {
  unsigned char *params;
  size_t params_bytes;
  CPMDCSession *session;
  int n;
  int *z;
  int *sorted;
  int *latched;
  int latched_n;
  double *coords;
  double *forces;
  double *shifted;
  double cell[9];
  double displ[3];
  double stress[9];
  double energy;
  int have_species;
  int have_coords;
  int have_cell;
  int fresh;
  int evaluations;
  int exited;
  char error[512];
};

static const char *const k_commands[] = {
    ">NATOMS",    ">CELL",    ">CELL_DISPL", ">ELEMENTS", ">TYPES", ">COORDS",
    "<ENERGY",    "<PE",      "<KE_ELEC",    "<FORCES",   "<STRESS", "<NATOMS",
    "<@",         "EXIT",     NULL};

static int fail(cpmdc_mdi *eng, const char *fmt, ...) {
  va_list ap;
  va_start(ap, fmt);
  vsnprintf(eng->error, sizeof eng->error, fmt, ap);
  va_end(ap);
  return -1;
}

static int cmp_int(const void *a, const void *b) {
  int ia = *(const int *)a;
  int ib = *(const int *)b;
  return (ia > ib) - (ia < ib);
}

static int finite_n(const double *v, int n) {
  for (int i = 0; i < n; i++) {
    if (!isfinite(v[i]))
      return 0;
  }
  return 1;
}

cpmdc_mdi *cpmdc_mdi_create(const void *params, size_t nbytes) {
  if (!params || nbytes == 0)
    return NULL;
  cpmdc_mdi *eng = calloc(1, sizeof *eng);
  if (!eng)
    return NULL;
  eng->params = malloc(nbytes);
  if (!eng->params) {
    free(eng);
    return NULL;
  }
  memcpy(eng->params, params, nbytes);
  eng->params_bytes = nbytes;
  return eng;
}

void cpmdc_mdi_destroy(cpmdc_mdi *eng) {
  if (!eng)
    return;
  cpmdc_session_destroy(eng->session);
  free(eng->params);
  free(eng->z);
  free(eng->sorted);
  free(eng->latched);
  free(eng->coords);
  free(eng->forces);
  free(eng->shifted);
  free(eng);
}

const char *const *cpmdc_mdi_commands(void) { return k_commands; }

const char *cpmdc_mdi_error(const cpmdc_mdi *eng) {
  return eng && eng->error[0] ? eng->error : "";
}

int cpmdc_mdi_evaluations(const cpmdc_mdi *eng) {
  return eng ? eng->evaluations : 0;
}

int cpmdc_mdi_exited(const cpmdc_mdi *eng) { return eng ? eng->exited : 0; }

static int set_count(cpmdc_mdi *eng, int n) {
  if (n < 1 || n > INT_MAX / 3)
    return fail(eng, "MDI >NATOMS is outside 1..INT_MAX/3");
  if (n != eng->n) {
    int *z = malloc((size_t)n * sizeof(int));
    int *sorted = malloc((size_t)n * sizeof(int));
    int *latched = malloc((size_t)n * sizeof(int));
    double *coords = malloc((size_t)n * 3 * sizeof(double));
    double *forces = malloc((size_t)n * 3 * sizeof(double));
    double *shifted = malloc((size_t)n * 3 * sizeof(double));
    if (!z || !sorted || !latched || !coords || !forces || !shifted) {
      free(z);
      free(sorted);
      free(latched);
      free(coords);
      free(forces);
      free(shifted);
      return fail(eng, "out of memory");
    }
    free(eng->z);
    free(eng->sorted);
    free(eng->latched);
    free(eng->coords);
    free(eng->forces);
    free(eng->shifted);
    eng->z = z;
    eng->sorted = sorted;
    eng->latched = latched;
    eng->coords = coords;
    eng->forces = forces;
    eng->shifted = shifted;
    eng->n = n;
    eng->latched_n = 0;
  }
  /* >NATOMS is a new system. LAMMPS sends the cell, species, and coordinates
   * again before it reads. Dropping them here is what makes a following
   * nonperiodic system stop using the previous cell. */
  eng->have_species = 0;
  eng->have_coords = 0;
  eng->have_cell = 0;
  eng->displ[0] = eng->displ[1] = eng->displ[2] = 0.0;
  eng->fresh = 0;
  return 0;
}

static int take_species(cpmdc_mdi *eng, const int *z) {
  if (eng->n < 1)
    return fail(eng, "MDI command needs >NATOMS first");
  if (!z)
    return fail(eng, "MDI command needs a buffer");
  memcpy(eng->z, z, (size_t)eng->n * sizeof(int));
  eng->have_species = 1;
  eng->fresh = 0;
  return 0;
}

static int ensure(cpmdc_mdi *eng) {
  const char *why;
  if (eng->fresh)
    return 0;
  if (eng->n < 1)
    return fail(eng, "MDI command needs >NATOMS first");
  if (!eng->have_species)
    return fail(eng,
                "MDI command needs >ELEMENTS or >TYPES before an evaluation");
  if (!eng->have_coords)
    return fail(eng, "MDI command needs >COORDS before an evaluation");

  memcpy(eng->sorted, eng->z, (size_t)eng->n * sizeof(int));
  qsort(eng->sorted, (size_t)eng->n, sizeof(int), cmp_int);
  if (!eng->session || eng->latched_n != eng->n ||
      memcmp(eng->sorted, eng->latched, (size_t)eng->n * sizeof(int)) != 0) {
    cpmdc_session_destroy(eng->session);
    eng->session = NULL;
    eng->latched_n = 0;
    eng->session = cpmdc_session_create(eng->params, eng->params_bytes);
    if (!eng->session) {
      why = cpmdc_last_error();
      return fail(eng, "cpmdc session create failed: %s",
                  why && why[0] ? why : "unknown error");
    }
  }

  for (int i = 0; i < eng->n; i++) {
    eng->shifted[3 * i] = eng->coords[3 * i] - eng->displ[0];
    eng->shifted[3 * i + 1] = eng->coords[3 * i + 1] - eng->displ[1];
    eng->shifted[3 * i + 2] = eng->coords[3 * i + 2] - eng->displ[2];
  }
  if (cpmdc_session_calculate_au(
          eng->session, eng->n, eng->shifted, eng->z,
          eng->have_cell ? eng->cell : NULL, &eng->energy, eng->forces,
          NULL) != 0) {
    why = cpmdc_last_error();
    return fail(eng, "cpmdc evaluation failed: %s",
                why && why[0] ? why : "unknown error");
  }
  if (!isfinite(eng->energy))
    return fail(eng, "the library returned a non-finite energy");
  if (!finite_n(eng->forces, 3 * eng->n))
    return fail(eng, "the library returned a non-finite force");

  CPMDCStressTensor tensor;
  memset(&tensor, 0, sizeof tensor);
  if (cpmdc_session_last_stress(eng->session, &tensor) == 0 && tensor.valid) {
    memcpy(eng->stress, tensor.values, sizeof eng->stress);
    if (!finite_n(eng->stress, 9))
      return fail(eng, "the library returned a non-finite stress");
  } else {
    memset(eng->stress, 0, sizeof eng->stress);
  }
  memcpy(eng->latched, eng->sorted, (size_t)eng->n * sizeof(int));
  eng->latched_n = eng->n;
  eng->fresh = 1;
  eng->evaluations++;
  return 0;
}

int cpmdc_mdi_payload(cpmdc_mdi *eng, const char *command, int *direction,
                      int *count) {
  if (!eng || !command || !direction || !count)
    return -1;
  eng->error[0] = '\0';
  *direction = CPMDC_MDI_NONE;
  *count = 0;
  if (strcmp(command, "EXIT") == 0)
    return 0;
  if (strcmp(command, ">NATOMS") == 0) {
    *direction = CPMDC_MDI_RECV_INT;
    *count = 1;
    return 0;
  }
  if (strcmp(command, ">CELL") == 0) {
    *direction = CPMDC_MDI_RECV_DOUBLE;
    *count = 9;
    return 0;
  }
  if (strcmp(command, ">CELL_DISPL") == 0) {
    *direction = CPMDC_MDI_RECV_DOUBLE;
    *count = 3;
    return 0;
  }
  if (strcmp(command, ">ELEMENTS") == 0 || strcmp(command, ">TYPES") == 0) {
    if (eng->n < 1)
      return fail(eng, "MDI command needs >NATOMS first");
    *direction = CPMDC_MDI_RECV_INT;
    *count = eng->n;
    return 0;
  }
  if (strcmp(command, ">COORDS") == 0) {
    if (eng->n < 1)
      return fail(eng, "MDI command needs >NATOMS first");
    *direction = CPMDC_MDI_RECV_DOUBLE;
    *count = 3 * eng->n;
    return 0;
  }
  if (strcmp(command, "<ENERGY") == 0 || strcmp(command, "<PE") == 0 ||
      strcmp(command, "<KE_ELEC") == 0) {
    *direction = CPMDC_MDI_SEND_DOUBLE;
    *count = 1;
    return 0;
  }
  if (strcmp(command, "<FORCES") == 0) {
    if (eng->n < 1)
      return fail(eng, "MDI command needs >NATOMS first");
    *direction = CPMDC_MDI_SEND_DOUBLE;
    *count = 3 * eng->n;
    return 0;
  }
  if (strcmp(command, "<STRESS") == 0) {
    *direction = CPMDC_MDI_SEND_DOUBLE;
    *count = 9;
    return 0;
  }
  if (strcmp(command, "<NATOMS") == 0) {
    *direction = CPMDC_MDI_SEND_INT;
    *count = 1;
    return 0;
  }
  if (strcmp(command, "<@") == 0) {
    *direction = CPMDC_MDI_SEND_NAME;
    *count = 9;
    return 0;
  }
  return fail(eng, "MDI command is not one fix mdi/qm sends: %s", command);
}

int cpmdc_mdi_command(cpmdc_mdi *eng, const char *command, const void *incoming,
                      void *outgoing) {
  if (!eng || !command)
    return -1;
  eng->error[0] = '\0';
  if (strcmp(command, "EXIT") == 0) {
    eng->exited = 1;
    return 0;
  }
  if (strcmp(command, ">NATOMS") == 0) {
    if (!incoming)
      return fail(eng, "MDI command needs a buffer");
    return set_count(eng, *(const int *)incoming);
  }
  if (strcmp(command, ">CELL") == 0) {
    if (!incoming)
      return fail(eng, "MDI command needs a buffer");
    memcpy(eng->cell, incoming, sizeof eng->cell);
    eng->have_cell = 1;
    eng->fresh = 0;
    return 0;
  }
  if (strcmp(command, ">CELL_DISPL") == 0) {
    if (!incoming)
      return fail(eng, "MDI command needs a buffer");
    memcpy(eng->displ, incoming, sizeof eng->displ);
    eng->fresh = 0;
    return 0;
  }
  if (strcmp(command, ">ELEMENTS") == 0 || strcmp(command, ">TYPES") == 0)
    return take_species(eng, incoming);
  if (strcmp(command, ">COORDS") == 0) {
    if (eng->n < 1)
      return fail(eng, "MDI command needs >NATOMS first");
    if (!incoming)
      return fail(eng, "MDI command needs a buffer");
    memcpy(eng->coords, incoming, (size_t)eng->n * 3 * sizeof(double));
    eng->have_coords = 1;
    eng->fresh = 0;
    return 0;
  }
  if (strcmp(command, "<KE_ELEC") == 0) {
    if (!outgoing)
      return fail(eng, "MDI command needs a buffer");
    /* Already inside the cpmdc total energy. LAMMPS adds this to <PE. */
    *(double *)outgoing = 0.0;
    return 0;
  }
  if (strcmp(command, "<ENERGY") == 0 || strcmp(command, "<PE") == 0) {
    if (!outgoing)
      return fail(eng, "MDI command needs a buffer");
    if (ensure(eng) != 0)
      return -1;
    *(double *)outgoing = eng->energy;
    return 0;
  }
  if (strcmp(command, "<FORCES") == 0) {
    if (!outgoing)
      return fail(eng, "MDI command needs a buffer");
    if (ensure(eng) != 0)
      return -1;
    memcpy(outgoing, eng->forces, (size_t)eng->n * 3 * sizeof(double));
    return 0;
  }
  if (strcmp(command, "<STRESS") == 0) {
    if (!outgoing)
      return fail(eng, "MDI command needs a buffer");
    if (ensure(eng) != 0)
      return -1;
    memcpy(outgoing, eng->stress, sizeof eng->stress);
    return 0;
  }
  if (strcmp(command, "<NATOMS") == 0) {
    if (!outgoing)
      return fail(eng, "MDI command needs a buffer");
    *(int *)outgoing = eng->n;
    return 0;
  }
  if (strcmp(command, "<@") == 0) {
    if (!outgoing)
      return fail(eng, "MDI command needs a buffer");
    memcpy(outgoing, "@DEFAULT", 9);
    return 0;
  }
  return fail(eng, "MDI command is not one fix mdi/qm sends: %s", command);
}
