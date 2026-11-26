/**
 * Command state for LAMMPS fix mdi/qm (stable 22 Jul 2025, src/MDI/fix_mdi_qm.cpp).
 *
 * LAMMPS sends, in Bohr: >NATOMS, then >CELL_DISPL and >CELL when the box is
 * fully periodic, then >ELEMENTS (atomic numbers, from the elements keyword)
 * or >TYPES, then >COORDS on every step. fix mdi/qm copies atom->x and calls
 * domain->remap before >COORDS, so each coordinate lies inside the box. The
 * image flags are not applied. >CELL_DISPL is the box origin. This engine
 * subtracts that origin. When a previous step and a cell are known, it then
 * shifts each atom by the lattice vector that brings it nearest to its
 * previous position. The session sees a continuous path, so a cached
 * geometry, a reused neighbor list, and wavefunction extrapolation follow
 * the atoms across a face. A singular cell on a later step cannot be
 * unwrapped and fails the evaluation.
 *
 * LAMMPS then reads <PE and <KE_ELEC when both exist, otherwise <ENERGY, then
 * <FORCES, then <STRESS when virial is requested. It adds <PE and <KE_ELEC.
 * The cpmdc energy already includes the electronic kinetic energy, so
 * <KE_ELEC is 0 and <PE equals <ENERGY. <STRESS is nine Hartree/Bohr^3
 * values, row-major xx xy xz yx yy yz zx zy zz, which is
 * cpmdc_session_last_stress(). LAMMPS symmetrizes that layout, multiplies by
 * the box volume, and stores the extensive virial. A step with no stress
 * sends zeros. One session call serves every read until the next '>' command.
 * >NATOMS starts a system: the previous cell, origin, species, coordinates,
 * and unwrapped path are dropped, and LAMMPS sends the new ones before it
 * reads.
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
  double *prev;
  double cell[9];
  double displ[3];
  double stress[9];
  double energy;
  int have_species;
  int have_coords;
  int have_cell;
  int have_prev;
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
  free(eng->prev);
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
    double *prev = malloc((size_t)n * 3 * sizeof(double));
    if (!z || !sorted || !latched || !coords || !forces || !shifted || !prev) {
      free(z);
      free(sorted);
      free(latched);
      free(coords);
      free(forces);
      free(shifted);
      free(prev);
      return fail(eng, "out of memory");
    }
    free(eng->z);
    free(eng->sorted);
    free(eng->latched);
    free(eng->coords);
    free(eng->forces);
    free(eng->shifted);
    free(eng->prev);
    eng->z = z;
    eng->sorted = sorted;
    eng->latched = latched;
    eng->coords = coords;
    eng->forces = forces;
    eng->shifted = shifted;
    eng->prev = prev;
    eng->n = n;
    eng->latched_n = 0;
  }
  /* >NATOMS is a new system. LAMMPS sends the cell, species, and coordinates
   * again before it reads. Dropping them here is what makes a following
   * nonperiodic system stop using the previous cell. The unwrapped path is
   * dropped with them, so the next frame is not shifted toward the old one. */
  eng->have_species = 0;
  eng->have_coords = 0;
  eng->have_cell = 0;
  eng->have_prev = 0;
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

static double det3(double m[3][3]) {
  return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
         m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
         m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
}

/* inv with x = s rows, s = x inv, for row vectors. */
static int invert3(double m[3][3], double inv[3][3]) {
  const double d = det3(m);
  if (d == 0.0 || !isfinite(d))
    return -1;
  inv[0][0] = (m[1][1] * m[2][2] - m[1][2] * m[2][1]) / d;
  inv[0][1] = (m[0][2] * m[2][1] - m[0][1] * m[2][2]) / d;
  inv[0][2] = (m[0][1] * m[1][2] - m[0][2] * m[1][1]) / d;
  inv[1][0] = (m[1][2] * m[2][0] - m[1][0] * m[2][2]) / d;
  inv[1][1] = (m[0][0] * m[2][2] - m[0][2] * m[2][0]) / d;
  inv[1][2] = (m[0][2] * m[1][0] - m[0][0] * m[1][2]) / d;
  inv[2][0] = (m[1][0] * m[2][1] - m[1][1] * m[2][0]) / d;
  inv[2][1] = (m[0][1] * m[2][0] - m[0][0] * m[2][1]) / d;
  inv[2][2] = (m[0][0] * m[1][1] - m[0][1] * m[1][0]) / d;
  return 0;
}

/* Shift x by the lattice vector that brings it nearest to ref.
 * nearbyint rounds each fractional coordinate on its own. That step
 * assumes a restricted (lower-triangular) triclinic cell, one lattice
 * vector per row. */
static void nearest_image(double rows[3][3], double inv[3][3],
                          const double *ref, double *x) {
  double d[3], s[3];
  for (int k = 0; k < 3; k++)
    d[k] = x[k] - ref[k];
  /* nearbyint on each fractional coordinate is the parallelepiped
   * (cell-wrapped) minimum image, not the Voronoi (true nearest-image)
   * one. It is exact when each single-step |delta*inv| < 0.5. */
  for (int a = 0; a < 3; a++)
    s[a] = nearbyint(d[0] * inv[0][a] + d[1] * inv[1][a] + d[2] * inv[2][a]);
  for (int k = 0; k < 3; k++)
    x[k] -= s[0] * rows[0][k] + s[1] * rows[1][k] + s[2] * rows[2][k];
}

/* >CELL is one lattice vector per row. A later step is folded onto the
 * previous origin-subtracted positions so a face crossing stays continuous. */
static int unwrap_coords(cpmdc_mdi *eng) {
  double rows[3][3], inv[3][3];
  if (!eng->have_cell || !eng->have_prev)
    return 0;
  for (int a = 0; a < 3; a++)
    for (int k = 0; k < 3; k++)
      rows[a][k] = eng->cell[3 * a + k];
  if (invert3(rows, inv) != 0)
    return fail(eng, "MDI >CELL is singular");
  for (int i = 0; i < eng->n; i++)
    nearest_image(rows, inv, &eng->prev[3 * i], &eng->shifted[3 * i]);
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
  if (unwrap_coords(eng) != 0)
    return -1;
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
  memcpy(eng->prev, eng->shifted, (size_t)eng->n * 3 * sizeof(double));
  eng->have_prev = 1;
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
