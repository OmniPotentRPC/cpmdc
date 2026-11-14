/**
 * libcpmdc_gromacs: cpmdc behind the libcp2k C API of GROMACS' QM/MM module.
 *
 * GROMACS built with -DGMX_CP2K=ON (checked: 2024.6, 2025.5, 2026.4) calls
 * libcp2k once per MD step from QMMMForceProvider
 * (src/gromacs/applied_forces/qmmm/qmmmforceprovider.cpp): set the
 * positions of every QM and MM atom and the box, evaluate, read the energy
 * and the forces. This library answers those calls with one
 * cpmdc_session_calculate_au() on the QM atoms, so mdrun integrates on
 * cpmdc forces in process. GROMACS adds the energy to its "Quantum En."
 * term (F_EQM), which enters the potential and the conserved energy.
 *
 * The input is the file of grompp -qmi with qmmm-cp2k-qmmethod = INPUT.
 * GROMACS checks it for CHARGE, MULTIPLICITY, RUN_TYPE ENERGY_FORCE,
 * COORD_FILE_FORMAT PDB, CHARGE_EXTENDED TRUE, and COORD_FILE_NAME, and at
 * mdrun writes it out with the name of the PDB file it generates. That PDB
 * gives every atom's element, whether it is QM or MM, and the classical
 * charge of each MM atom (QM atoms carry 0: grompp removes their charges).
 * The cpmdc settings sit in a &CPMDC section:
 *
 *   PARAMS path           CPMDParams message (required)
 *   CELL PERIODIC|NONE    send the GROMACS box to cpmdc, or the deck's cell
 *   EMBEDDING MECHANICAL|NONE   QM-MM Coulomb from fixed charges, or none
 *   QM_CHARGES q ...      classical charge of each QM atom, in index order
 *   EWALD_RCUT r          real-space cutoff of that Coulomb sum, nm
 *   EWALD_PRECISION p     erfc(alpha rcut) and the reciprocal-space cutoff
 *
 * Mechanical embedding: the QM-MM Coulomb energy of the QM_CHARGES and the
 * MM charges, an Ewald sum when the box is periodic and a direct sum when
 * it is not. GROMACS keeps the QM-MM Lennard-Jones and the bonded terms
 * that cross the boundary; the QM density is not polarized by the MM
 * charges. &LINK sections (QM-MM covalent bonds) are refused.
 *
 * Positions: GROMACS wraps every atom into the box. QM atoms are unwrapped
 * against their previous positions, so each follows a continuous path; with
 * CELL NONE the first step also makes the QM region whole about its first
 * atom.
 *
 * Errors cannot be returned through this API: the library prints
 * "cpmdc-gromacs: ..." to stderr and to the CP2K output file GROMACS names,
 * and exits with status 1 (MPI_Abort under MPI).
 */
#define _POSIX_C_SOURCE 200809L
#include "libcp2k.h"

#include "cpmdc.h"

#include <complex.h>
#include <ctype.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

#ifdef CPMDC_GROMACS_MPI
#include <mpi.h>
#endif

#define PROG "cpmdc-gromacs"
#define MAX_Z 118
#define MAX_ENVS 16
/* Same value as GROMACS' c_bohr2Nm and cpmdc.BOHR_TO_ANGSTROM / 10. */
#define BOHR_TO_NM 0.0529177210903
#define PI 3.14159265358979323846

static const char *const k_symbols[MAX_Z] = {
    "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na", "Mg",
    "Al", "Si", "P",  "S",  "Cl", "Ar", "K",  "Ca", "Sc", "Ti", "V",  "Cr",
    "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
    "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "In", "Sn", "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf",
    "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po",
    "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm",
    "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs",
    "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"};

enum { CELL_DEFAULT = -1, CELL_NONE = 0, CELL_PERIODIC = 1 };
enum { EMBED_UNSET = -1, EMBED_NONE = 0, EMBED_MECHANICAL = 1 };

typedef struct {
  int used;
  int rank, size;
#ifdef CPMDC_GROMACS_MPI
  MPI_Comm comm;
#endif
  FILE *log;
  /* rank 0 */
  char *params;
  size_t params_len;
  CPMDCSession *session;
  int adopted;
  int cell_mode, embedding;
  double ewald_rcut_nm, ewald_precision;
  int nqm, nmm, qm_charges_listed;
  int *qm_index, *qm_z, *mm_index;
  double *qm_charge, *mm_charge;
  double *qm_x, *qm_prev, *qm_f;
  int have_prev;
  /* every rank */
  int natoms;
  double *pos, *forces;
  double rows[3][3];
  int have_positions, have_cell;
  double energy;
  int evaluated;
} Env;

static Env g_envs[MAX_ENVS];

/* --- errors -------------------------------------------------------------- */

static void die(Env *env, const char *fmt, ...) {
  char text[2048];
  va_list ap;
  va_start(ap, fmt);
  vsnprintf(text, sizeof(text), fmt, ap);
  va_end(ap);
  fprintf(stderr, PROG ": %s\n", text);
  fflush(stderr);
  if (env && env->log) {
    fprintf(env->log, PROG ": error: %s\n", text);
    fflush(env->log);
  }
#ifdef CPMDC_GROMACS_MPI
  int up = 0;
  MPI_Initialized(&up);
  if (up && env && env->size > 1)
    MPI_Abort(env->comm, 1);
#endif
  exit(EXIT_FAILURE);
}

static Env *lookup(force_env_t id, const char *call) {
  if (id < 1 || id > MAX_ENVS || !g_envs[id - 1].used)
    die(NULL, "%s: no force environment %d", call, id);
  return &g_envs[id - 1];
}

static void *alloc(Env *env, size_t count, size_t size) {
  void *p = calloc(count ? count : 1, size);
  if (!p)
    die(env, "out of memory");
  return p;
}

/* --- input --------------------------------------------------------------- */

static int symbol_number(const char *s) {
  for (int i = 0; i < MAX_Z; ++i)
    if (strcasecmp(k_symbols[i], s) == 0)
      return i + 1;
  return 0;
}

static char *read_file(Env *env, const char *path, size_t *len, const char *what) {
  FILE *f = fopen(path, "rb");
  if (!f)
    die(env, "cannot open %s %s", what, path);
  size_t cap = 1 << 16, n = 0;
  char *buf = alloc(env, cap + 1, 1);
  for (;;) {
    size_t got = fread(buf + n, 1, cap - n, f);
    n += got;
    if (n < cap)
      break;
    cap *= 2;
    char *grown = realloc(buf, cap + 1);
    if (!grown)
      die(env, "out of memory");
    buf = grown;
  }
  fclose(f);
  buf[n] = '\0';
  if (len)
    *len = n;
  return buf;
}

/* Split a line into at most `max` words; comments start with # or !. */
static int words(char *line, char **out, int max) {
  int n = 0;
  for (char *p = line; *p; ++p)
    if (*p == '#' || *p == '!') {
      *p = '\0';
      break;
    }
  char *save = NULL;
  for (char *tok = strtok_r(line, " \t\r", &save); tok && n < max;
       tok = strtok_r(NULL, " \t\r", &save))
    out[n++] = tok;
  return n;
}

static double number(Env *env, const char *text, const char *key) {
  char *end = NULL;
  double v = strtod(text, &end);
  if (end == text || *end != '\0' || !isfinite(v))
    die(env, "&CPMDC %s: %s is not a number", key, text);
  return v;
}

static void parse_input(Env *env, const char *path, char *pdb_path, size_t pdb_cap) {
  char *text = read_file(env, path, NULL, "the QM/MM input");
  int in_cpmdc = 0, saw_cpmdc = 0, nq = 0, capq = 64;
  char *params_path = NULL;
  pdb_path[0] = '\0';
  env->cell_mode = CELL_DEFAULT;
  env->embedding = EMBED_UNSET;
  env->ewald_rcut_nm = 0.0;
  env->ewald_precision = 1e-8;
  env->qm_charge = alloc(env, (size_t)capq, sizeof(double));
  char *save = NULL;
  for (char *line = strtok_r(text, "\n", &save); line; line = strtok_r(NULL, "\n", &save)) {
    char *w[256];
    int n = words(line, w, 256);
    if (n == 0)
      continue;
    if (strcasecmp(w[0], "&LINK") == 0)
      die(env, "%s has a &LINK section: cpmdc-gromacs supports QM regions of whole "
               "molecules only, without QM-MM covalent bonds",
          path);
    if (strcasecmp(w[0], "COORD_FILE_NAME") == 0 && n >= 2) {
      snprintf(pdb_path, pdb_cap, "%s", w[1]);
      continue;
    }
    if (!in_cpmdc) {
      if (strcasecmp(w[0], "&CPMDC") == 0) {
        in_cpmdc = 1;
        saw_cpmdc = 1;
      }
      continue;
    }
    if (strcasecmp(w[0], "&END") == 0) {
      in_cpmdc = 0;
      continue;
    }
    if (w[0][0] == '&')
      die(env, "&CPMDC takes no subsections; found %s", w[0]);
    if (strcasecmp(w[0], "QM_CHARGES") == 0) {
      for (int i = 1; i < n; ++i) {
        if (nq == capq) {
          capq *= 2;
          double *grown = realloc(env->qm_charge, (size_t)capq * sizeof(double));
          if (!grown)
            die(env, "out of memory");
          env->qm_charge = grown;
        }
        env->qm_charge[nq++] = number(env, w[i], "QM_CHARGES");
      }
      continue;
    }
    if (n != 2)
      die(env, "&CPMDC %s takes one value", w[0]);
    if (strcasecmp(w[0], "PARAMS") == 0) {
      free(params_path);
      params_path = strdup(w[1]);
    } else if (strcasecmp(w[0], "CELL") == 0) {
      if (strcasecmp(w[1], "PERIODIC") == 0)
        env->cell_mode = CELL_PERIODIC;
      else if (strcasecmp(w[1], "NONE") == 0)
        env->cell_mode = CELL_NONE;
      else
        die(env, "&CPMDC CELL is PERIODIC or NONE, not %s", w[1]);
    } else if (strcasecmp(w[0], "EMBEDDING") == 0) {
      if (strcasecmp(w[1], "MECHANICAL") == 0)
        env->embedding = EMBED_MECHANICAL;
      else if (strcasecmp(w[1], "NONE") == 0)
        env->embedding = EMBED_NONE;
      else if (strcasecmp(w[1], "ELECTROSTATIC") == 0)
        die(env, "&CPMDC EMBEDDING ELECTROSTATIC is not implemented: cpmdc takes no "
                 "external point charges; use MECHANICAL or NONE");
      else
        die(env, "&CPMDC EMBEDDING is MECHANICAL or NONE, not %s", w[1]);
    } else if (strcasecmp(w[0], "EWALD_RCUT") == 0) {
      env->ewald_rcut_nm = number(env, w[1], "EWALD_RCUT");
      if (env->ewald_rcut_nm <= 0.0)
        die(env, "&CPMDC EWALD_RCUT must be positive");
    } else if (strcasecmp(w[0], "EWALD_PRECISION") == 0) {
      env->ewald_precision = number(env, w[1], "EWALD_PRECISION");
      if (!(env->ewald_precision > 0.0 && env->ewald_precision < 1.0))
        die(env, "&CPMDC EWALD_PRECISION must lie between 0 and 1");
    } else {
      die(env, "unknown &CPMDC keyword %s", w[0]);
    }
  }
  free(text);
  if (!saw_cpmdc)
    die(env, "%s has no &CPMDC section. Set qmmm-cp2k-qmmethod = INPUT and pass the "
             "cpmdc input to grompp with -qmi",
        path);
  if (!params_path)
    die(env, "&CPMDC needs PARAMS, the path of a CPMDParams message");
  if (pdb_path[0] == '\0')
    die(env, "%s names no COORD_FILE_NAME", path);
  env->params = read_file(env, params_path, &env->params_len, "the CPMDParams message");
  if (env->params_len == 0)
    die(env, "the CPMDParams message %s is empty", params_path);
  free(params_path);
  env->qm_charges_listed = nq;
}

static void parse_pdb(Env *env, const char *path) {
  char *text = read_file(env, path, NULL, "the QM/MM PDB file");
  int cap = 1024, n = 0, nqm = 0, nmm = 0;
  int *is_qm = alloc(env, (size_t)cap, sizeof(int));
  int *z = alloc(env, (size_t)cap, sizeof(int));
  double *q = alloc(env, (size_t)cap, sizeof(double));
  char *save = NULL;
  for (char *line = strtok_r(text, "\n", &save); line; line = strtok_r(NULL, "\n", &save)) {
    if (strncmp(line, "ATOM", 4) != 0 && strncmp(line, "HETATM", 6) != 0)
      continue;
    char *w[32];
    int nw = words(line, w, 32);
    if (nw < 6)
      die(env, "%s: cannot read the atom line %d", path, n + 1);
    if (n == cap) {
      cap *= 2;
      is_qm = realloc(is_qm, (size_t)cap * sizeof(int));
      z = realloc(z, (size_t)cap * sizeof(int));
      q = realloc(q, (size_t)cap * sizeof(double));
      if (!is_qm || !z || !q)
        die(env, "out of memory");
    }
    is_qm[n] = strcmp(w[3], "QM") == 0;
    if (!is_qm[n] && strcmp(w[3], "MM") != 0)
      die(env, "%s: atom %d is neither QM nor MM", path, n + 1);
    z[n] = symbol_number(w[nw - 2]);
    char *end = NULL;
    q[n] = strtod(w[nw - 1], &end);
    if (end == w[nw - 1] || *end != '\0')
      die(env, "%s: atom %d has no charge", path, n + 1);
    if (is_qm[n]) {
      if (z[n] == 0)
        die(env, "QM atom %d is %s, not an element (a virtual site?); keep virtual sites "
                 "out of the QM group",
            n + 1, w[nw - 2]);
      ++nqm;
    } else {
      ++nmm;
    }
    ++n;
  }
  free(text);
  if (nqm == 0)
    die(env, "%s has no QM atoms", path);
  const int nq = env->qm_charges_listed;
  env->natoms = n;
  env->nqm = nqm;
  env->nmm = nmm;
  env->qm_index = alloc(env, (size_t)nqm, sizeof(int));
  env->qm_z = alloc(env, (size_t)nqm, sizeof(int));
  env->mm_index = alloc(env, (size_t)nmm, sizeof(int));
  env->mm_charge = alloc(env, (size_t)nmm, sizeof(double));
  int iq = 0, im = 0;
  int charged_mm = 0;
  for (int i = 0; i < n; ++i) {
    if (is_qm[i]) {
      env->qm_index[iq] = i;
      env->qm_z[iq++] = z[i];
    } else {
      env->mm_index[im] = i;
      env->mm_charge[im++] = q[i];
      charged_mm |= q[i] != 0.0;
    }
  }
  free(is_qm);
  free(z);
  free(q);
  if (env->cell_mode == CELL_DEFAULT)
    env->cell_mode = nmm == 0 ? CELL_PERIODIC : CELL_NONE;
  if (env->embedding == EMBED_UNSET) {
    if (charged_mm)
      die(env, "the MM atoms carry charges: set &CPMDC EMBEDDING MECHANICAL (with "
               "QM_CHARGES) or EMBEDDING NONE");
    env->embedding = EMBED_NONE;
  }
  if (env->embedding == EMBED_MECHANICAL) {
    if (nq != nqm)
      die(env, "&CPMDC QM_CHARGES lists %d charges for %d QM atoms", nq, nqm);
  } else if (nq != 0) {
    die(env, "&CPMDC QM_CHARGES needs EMBEDDING MECHANICAL");
  }
  env->qm_x = alloc(env, 3 * (size_t)nqm, sizeof(double));
  env->qm_prev = alloc(env, 3 * (size_t)nqm, sizeof(double));
  env->qm_f = alloc(env, 3 * (size_t)nqm, sizeof(double));
}

/* --- geometry ------------------------------------------------------------ */

static double det3(double m[3][3]) {
  return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
         m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
         m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
}

/* inv with x = s rows, s = x inv (row vectors). */
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

/* Shift x by the lattice vector that brings it nearest to ref. */
static void nearest_image(double rows[3][3], double inv[3][3], const double *ref,
                          double *x) {
  double d[3], s[3];
  for (int k = 0; k < 3; ++k)
    d[k] = x[k] - ref[k];
  for (int a = 0; a < 3; ++a)
    s[a] = nearbyint(d[0] * inv[0][a] + d[1] * inv[1][a] + d[2] * inv[2][a]);
  for (int k = 0; k < 3; ++k)
    x[k] -= s[0] * rows[0][k] + s[1] * rows[1][k] + s[2] * rows[2][k];
}

/* Plane spacing along each lattice vector: volume over the opposite face. */
static void heights(double r[3][3], double h[3]) {
  const double v = fabs(det3(r));
  for (int a = 0; a < 3; ++a) {
    const double *b = r[(a + 1) % 3], *c = r[(a + 2) % 3];
    const double x = b[1] * c[2] - b[2] * c[1];
    const double y = b[2] * c[0] - b[0] * c[2];
    const double z = b[0] * c[1] - b[1] * c[0];
    h[a] = v / sqrt(x * x + y * y + z * z);
  }
}

/* --- QM-MM Coulomb ------------------------------------------------------- */

/* The splitting of GROMACS' calc_ewaldcoeff_q: erfc(alpha rc) = tol. */
static double ewald_alpha(double rc, double tol) {
  double beta = 5.0, low = 0.0, high;
  int i = 0;
  do {
    ++i;
    beta *= 2.0;
  } while (erfc(beta * rc) > tol);
  const int n = i + 60;
  high = beta;
  for (i = 0; i < n; ++i) {
    beta = 0.5 * (low + high);
    if (erfc(beta * rc) > tol)
      low = beta;
    else
      high = beta;
  }
  return beta;
}

/* Coulomb energy between the QM and the MM charges with every periodic
 * image (Ewald), Hartree with positions in Bohr; forces are added. The
 * real-space sum takes every image pair closer than rc, the reciprocal sum
 * every k with |k| <= kmax, and the neutralizing background its QM-MM
 * part, -pi Q_QM Q_MM / (V alpha^2). */
static double ewald_qm_mm(Env *env, double r[3][3], double rc, double tol,
                          const double *xq, double *fq, const double *xm, double *fm,
                          double *alpha_out, double *kmax_out) {
  const int nq = env->nqm, nm = env->nmm;
  const double *qq = env->qm_charge, *qm = env->mm_charge;
  double inv[3][3], h[3];
  if (invert3(r, inv) != 0)
    die(env, "the GROMACS box is singular");
  const double volume = fabs(det3(r));
  heights(r, h);
  const double alpha = ewald_alpha(rc, tol);
  const double kmax = 2.0 * alpha * sqrt(-log(tol));
  *alpha_out = alpha;
  *kmax_out = kmax;
  const double two_over_sqrt_pi = 2.0 / sqrt(PI);
  double energy = 0.0;

  int ni[3];
  for (int a = 0; a < 3; ++a)
    ni[a] = (int)ceil(rc / h[a] + 0.5);
  for (int i = 0; i < nq; ++i) {
    if (qq[i] == 0.0)
      continue;
    for (int j = 0; j < nm; ++j) {
      if (qm[j] == 0.0)
        continue;
      double d0[3];
      for (int k = 0; k < 3; ++k)
        d0[k] = xq[3 * i + k] - xm[3 * j + k];
      const double zero[3] = {0.0, 0.0, 0.0};
      nearest_image(r, inv, zero, d0);
      const double qij = qq[i] * qm[j];
      for (int n0 = -ni[0]; n0 <= ni[0]; ++n0)
        for (int n1 = -ni[1]; n1 <= ni[1]; ++n1)
          for (int n2 = -ni[2]; n2 <= ni[2]; ++n2) {
            double d[3], r2 = 0.0;
            for (int k = 0; k < 3; ++k) {
              d[k] = d0[k] + n0 * r[0][k] + n1 * r[1][k] + n2 * r[2][k];
              r2 += d[k] * d[k];
            }
            if (r2 >= rc * rc)
              continue;
            if (r2 == 0.0)
              die(env, "QM atom %d and MM atom %d coincide", env->qm_index[i] + 1,
                  env->mm_index[j] + 1);
            const double dist = sqrt(r2);
            const double e = erfc(alpha * dist) / dist;
            energy += qij * e;
            const double g =
                qij * (e + two_over_sqrt_pi * alpha * exp(-alpha * alpha * r2)) / r2;
            for (int k = 0; k < 3; ++k) {
              fq[3 * i + k] += g * d[k];
              fm[3 * j + k] -= g * d[k];
            }
          }
    }
  }

  int nk[3];
  for (int a = 0; a < 3; ++a) {
    const double len = sqrt(r[a][0] * r[a][0] + r[a][1] * r[a][1] + r[a][2] * r[a][2]);
    nk[a] = (int)floor(kmax * len / (2.0 * PI));
  }
  const int nall = nq + nm;
  double complex *e0 = alloc(env, (size_t)nall * (nk[0] + 1), sizeof(double complex));
  double complex *e1 = alloc(env, (size_t)nall * (2 * nk[1] + 1), sizeof(double complex));
  double complex *e2 = alloc(env, (size_t)nall * (2 * nk[2] + 1), sizeof(double complex));
  double complex *phase = alloc(env, (size_t)nall, sizeof(double complex));
  for (int p = 0; p < nall; ++p) {
    const double *x = p < nq ? &xq[3 * p] : &xm[3 * (p - nq)];
    double s[3];
    for (int a = 0; a < 3; ++a)
      s[a] = x[0] * inv[0][a] + x[1] * inv[1][a] + x[2] * inv[2][a];
    for (int n = 0; n <= nk[0]; ++n)
      e0[(size_t)p * (nk[0] + 1) + n] = cexp(I * 2.0 * PI * n * s[0]);
    for (int n = -nk[1]; n <= nk[1]; ++n)
      e1[(size_t)p * (2 * nk[1] + 1) + n + nk[1]] = cexp(I * 2.0 * PI * n * s[1]);
    for (int n = -nk[2]; n <= nk[2]; ++n)
      e2[(size_t)p * (2 * nk[2] + 1) + n + nk[2]] = cexp(I * 2.0 * PI * n * s[2]);
  }
  const double pref = 8.0 * PI / volume;
  const double inv4a2 = 1.0 / (4.0 * alpha * alpha);
  for (int n0 = 0; n0 <= nk[0]; ++n0)
    for (int n1 = -nk[1]; n1 <= nk[1]; ++n1)
      for (int n2 = -nk[2]; n2 <= nk[2]; ++n2) {
        if (n0 == 0 && (n1 < 0 || (n1 == 0 && n2 <= 0)))
          continue;
        double kv[3], k2 = 0.0;
        for (int c = 0; c < 3; ++c) {
          kv[c] = 2.0 * PI * (n0 * inv[c][0] + n1 * inv[c][1] + n2 * inv[c][2]);
          k2 += kv[c] * kv[c];
        }
        if (k2 > kmax * kmax)
          continue;
        const double f = exp(-k2 * inv4a2) / k2;
        double complex sq = 0.0, sm = 0.0;
        for (int p = 0; p < nall; ++p) {
          const double charge = p < nq ? qq[p] : qm[p - nq];
          phase[p] = e0[(size_t)p * (nk[0] + 1) + n0] *
                     e1[(size_t)p * (2 * nk[1] + 1) + n1 + nk[1]] *
                     e2[(size_t)p * (2 * nk[2] + 1) + n2 + nk[2]];
          if (p < nq)
            sq += charge * phase[p];
          else
            sm += charge * phase[p];
        }
        energy += pref * f * creal(sq * conj(sm));
        for (int p = 0; p < nall; ++p) {
          const double charge = p < nq ? qq[p] : qm[p - nq];
          if (charge == 0.0)
            continue;
          const double im = cimag(phase[p] * conj(p < nq ? sm : sq));
          double *out = p < nq ? &fq[3 * p] : &fm[3 * (p - nq)];
          for (int c = 0; c < 3; ++c)
            out[c] += pref * f * charge * im * kv[c];
        }
      }
  free(e0);
  free(e1);
  free(e2);
  free(phase);

  double sum_q = 0.0, sum_m = 0.0;
  for (int i = 0; i < nq; ++i)
    sum_q += qq[i];
  for (int j = 0; j < nm; ++j)
    sum_m += qm[j];
  energy -= PI * sum_q * sum_m / (volume * alpha * alpha);
  return energy;
}

static double direct_qm_mm(Env *env, const double *xq, double *fq, const double *xm,
                           double *fm) {
  double energy = 0.0;
  for (int i = 0; i < env->nqm; ++i)
    for (int j = 0; j < env->nmm; ++j) {
      const double qij = env->qm_charge[i] * env->mm_charge[j];
      if (qij == 0.0)
        continue;
      double d[3], r2 = 0.0;
      for (int k = 0; k < 3; ++k) {
        d[k] = xq[3 * i + k] - xm[3 * j + k];
        r2 += d[k] * d[k];
      }
      if (r2 == 0.0)
        die(env, "QM atom %d and MM atom %d coincide", env->qm_index[i] + 1,
            env->mm_index[j] + 1);
      const double dist = sqrt(r2);
      energy += qij / dist;
      const double g = qij / (r2 * dist);
      for (int k = 0; k < 3; ++k) {
        fq[3 * i + k] += g * d[k];
        fm[3 * j + k] -= g * d[k];
      }
    }
  return energy;
}

/* --- one evaluation ------------------------------------------------------ */

static void open_session(Env *env) {
#ifdef CPMDC_GROMACS_MPI
  if (env->size > 1 && !env->adopted) {
    MPI_Comm self = MPI_COMM_SELF;
    if (cpmdc_adopt_calculator_comm(&self, sizeof(self), 1) != 0 && cpmdc_available()) {
      const char *text = cpmdc_last_error();
      die(env, "%s", text && *text ? text : "cpmdc_adopt_calculator_comm failed");
    }
    env->adopted = 1;
  }
#endif
  env->session = cpmdc_session_create(env->params, env->params_len);
  if (!env->session) {
    const char *text = cpmdc_last_error();
    die(env, "%s", text && *text ? text : "cpmdc_session_create failed");
  }
}

/* Rank 0: energy and forces of every atom, Hartree and Hartree/Bohr. */
static int evaluate(Env *env, char *why, size_t why_cap) {
  const int nq = env->nqm;
  double inv[3][3];
  const int periodic = env->have_cell && invert3(env->rows, inv) == 0;
  if (env->have_cell && !periodic) {
    snprintf(why, why_cap, "the GROMACS box is singular");
    return -1;
  }
  if (env->cell_mode == CELL_PERIODIC && !periodic) {
    snprintf(why, why_cap, "&CPMDC CELL PERIODIC needs a GROMACS box");
    return -1;
  }
  for (int i = 0; i < nq; ++i)
    memcpy(&env->qm_x[3 * i], &env->pos[3 * (size_t)env->qm_index[i]], 3 * sizeof(double));
  if (periodic) {
    if (env->have_prev) {
      for (int i = 0; i < nq; ++i)
        nearest_image(env->rows, inv, &env->qm_prev[3 * i], &env->qm_x[3 * i]);
    } else if (env->cell_mode == CELL_NONE) {
      for (int i = 1; i < nq; ++i)
        nearest_image(env->rows, inv, &env->qm_x[0], &env->qm_x[3 * i]);
    }
  }
  if (!env->session)
    open_session(env);
  double qm_energy = 0.0;
  if (cpmdc_session_calculate_au(env->session, nq, env->qm_z, env->qm_x,
                                 env->cell_mode == CELL_PERIODIC ? &env->rows[0][0] : NULL,
                                 &qm_energy, env->qm_f, NULL) != 0) {
    const char *text = cpmdc_last_error();
    snprintf(why, why_cap, "%s", text && *text ? text : "cpmdc_session_calculate_au failed");
    return -1;
  }
  memcpy(env->qm_prev, env->qm_x, 3 * (size_t)nq * sizeof(double));
  env->have_prev = 1;

  memset(env->forces, 0, 3 * (size_t)env->natoms * sizeof(double));
  double coulomb = 0.0;
  if (env->embedding == EMBED_MECHANICAL && env->nmm > 0) {
    const int nm = env->nmm;
    double *xm = alloc(env, 3 * (size_t)nm, sizeof(double));
    double *fm = alloc(env, 3 * (size_t)nm, sizeof(double));
    double *fq = alloc(env, 3 * (size_t)nq, sizeof(double));
    for (int j = 0; j < nm; ++j)
      memcpy(&xm[3 * j], &env->pos[3 * (size_t)env->mm_index[j]], 3 * sizeof(double));
    if (periodic) {
      double h[3], rc, alpha, kmax;
      heights(env->rows, h);
      rc = env->ewald_rcut_nm > 0.0 ? env->ewald_rcut_nm / BOHR_TO_NM
                                    : fmin(1.0 / BOHR_TO_NM, 0.5 * fmin(h[0], fmin(h[1], h[2])));
      coulomb = ewald_qm_mm(env, env->rows, rc, env->ewald_precision, env->qm_x, fq, xm, fm,
                            &alpha, &kmax);
      if (env->log && env->evaluated == 0)
        fprintf(env->log, "ewald rcut_bohr %.17g alpha_per_bohr %.17g kmax_per_bohr %.17g\n",
                rc, alpha, kmax);
    } else {
      coulomb = direct_qm_mm(env, env->qm_x, fq, xm, fm);
    }
    for (int j = 0; j < nm; ++j)
      for (int k = 0; k < 3; ++k)
        env->forces[3 * (size_t)env->mm_index[j] + k] = fm[3 * j + k];
    for (int i = 0; i < nq; ++i)
      for (int k = 0; k < 3; ++k)
        env->forces[3 * (size_t)env->qm_index[i] + k] = fq[3 * i + k];
    free(xm);
    free(fm);
    free(fq);
  }
  for (int i = 0; i < nq; ++i)
    for (int k = 0; k < 3; ++k)
      env->forces[3 * (size_t)env->qm_index[i] + k] += env->qm_f[3 * i + k];
  env->energy = qm_energy + coulomb;
  if (env->log && env->evaluated == 0) {
    fprintf(env->log, "first step: qm_energy_hartree %.17g coulomb_qm_mm_hartree %.17g\n",
            qm_energy, coulomb);
    fflush(env->log);
  }
  return 0;
}

/* --- the libcp2k API ----------------------------------------------------- */

void cp2k_get_version(char *version_str, int str_length) {
  if (version_str && str_length > 0)
    snprintf(version_str, (size_t)str_length, "%s", cpmdc_version());
}

void cp2k_init(void) {}

void cp2k_init_without_mpi(void) {}

void cp2k_finalize(void) {}

void cp2k_finalize_without_mpi(void) {}

static void create(force_env_t *out, const char *input, const char *output, Env *proto) {
  int slot = -1;
  for (int i = 0; i < MAX_ENVS && slot < 0; ++i)
    if (!g_envs[i].used)
      slot = i;
  if (slot < 0)
    die(NULL, "more than %d force environments", MAX_ENVS);
  Env *env = &g_envs[slot];
  *env = *proto;
  env->used = 1;
  int natoms = 0;
  if (env->rank == 0) {
    if (output && *output) {
      env->log = fopen(output, "w");
      if (env->log)
        fprintf(env->log, PROG " on %s, input %s\n", cpmdc_version(), input);
    }
    char pdb[4096];
    parse_input(env, input, pdb, sizeof(pdb));
    parse_pdb(env, pdb);
    natoms = env->natoms;
    if (env->log) {
      fprintf(env->log,
              "atoms %d qm %d mm %d cell %s embedding %s ewald_precision %.3g ranks %d\n",
              env->natoms, env->nqm, env->nmm,
              env->cell_mode == CELL_PERIODIC ? "PERIODIC" : "NONE",
              env->embedding == EMBED_MECHANICAL ? "MECHANICAL" : "NONE",
              env->ewald_precision, env->size);
      fflush(env->log);
    }
  }
#ifdef CPMDC_GROMACS_MPI
  if (env->size > 1)
    MPI_Bcast(&natoms, 1, MPI_INT, 0, env->comm);
#endif
  env->natoms = natoms;
  env->pos = alloc(env, 3 * (size_t)natoms, sizeof(double));
  env->forces = alloc(env, 3 * (size_t)natoms, sizeof(double));
  *out = slot + 1;
}

void cp2k_create_force_env(force_env_t *new_force_env, const char *input_file_path,
                           const char *output_file_path) {
  Env proto = {0};
  proto.size = 1;
#ifdef CPMDC_GROMACS_MPI
  proto.comm = MPI_COMM_SELF;
#endif
  create(new_force_env, input_file_path, output_file_path, &proto);
}

void cp2k_create_force_env_comm(force_env_t *new_force_env, const char *input_file_path,
                                const char *output_file_path, int mpi_comm) {
  Env proto = {0};
#ifdef CPMDC_GROMACS_MPI
  proto.comm = MPI_Comm_f2c((MPI_Fint)mpi_comm);
  MPI_Comm_rank(proto.comm, &proto.rank);
  MPI_Comm_size(proto.comm, &proto.size);
#else
  (void)mpi_comm;
  (void)new_force_env;
  (void)input_file_path;
  (void)output_file_path;
  die(&proto, "libcpmdc_gromacs was built without MPI; use a GROMACS without an external "
              "MPI library (thread-MPI, -ntmpi 1) or rebuild with MPI");
#endif
  create(new_force_env, input_file_path, output_file_path, &proto);
}

void cp2k_destroy_force_env(force_env_t force_env) {
  Env *env = lookup(force_env, "cp2k_destroy_force_env");
  if (env->session)
    cpmdc_session_destroy(env->session);
  if (env->log)
    fclose(env->log);
  free(env->params);
  free(env->qm_index);
  free(env->qm_z);
  free(env->mm_index);
  free(env->qm_charge);
  free(env->mm_charge);
  free(env->qm_x);
  free(env->qm_prev);
  free(env->qm_f);
  free(env->pos);
  free(env->forces);
  memset(env, 0, sizeof(*env));
}

static void check_count(Env *env, int n_el, const char *call) {
  if (n_el != 3 * env->natoms)
    die(env, "%s got %d values; the PDB file has %d atoms (%d values)", call, n_el,
        env->natoms, 3 * env->natoms);
}

void cp2k_set_positions(force_env_t force_env, const double *new_pos, int n_el) {
  Env *env = lookup(force_env, "cp2k_set_positions");
  check_count(env, n_el, "cp2k_set_positions");
  memcpy(env->pos, new_pos, (size_t)n_el * sizeof(double));
  env->have_positions = 1;
}

void cp2k_set_cell(force_env_t force_env, const double *new_cell) {
  Env *env = lookup(force_env, "cp2k_set_cell");
  int any = 0;
  for (int a = 0; a < 3; ++a)
    for (int k = 0; k < 3; ++k) {
      env->rows[a][k] = new_cell[3 * k + a];
      any |= env->rows[a][k] != 0.0;
    }
  env->have_cell = any;
}

void cp2k_get_natom(force_env_t force_env, int *natom) {
  *natom = lookup(force_env, "cp2k_get_natom")->natoms;
}

void cp2k_get_positions(force_env_t force_env, double *pos, int n_el) {
  Env *env = lookup(force_env, "cp2k_get_positions");
  check_count(env, n_el, "cp2k_get_positions");
  memcpy(pos, env->pos, (size_t)n_el * sizeof(double));
}

void cp2k_calc_energy_force(force_env_t force_env) {
  Env *env = lookup(force_env, "cp2k_calc_energy_force");
  if (!env->have_positions)
    die(env, "cp2k_calc_energy_force before cp2k_set_positions");
  char why[1024] = {0};
  double head[2] = {0.0, 0.0};
  if (env->rank == 0) {
    if (evaluate(env, why, sizeof(why)) == 0)
      head[0] = 1.0;
    head[1] = env->energy;
  }
#ifdef CPMDC_GROMACS_MPI
  if (env->size > 1) {
    MPI_Bcast(head, 2, MPI_DOUBLE, 0, env->comm);
    if (head[0] == 1.0)
      MPI_Bcast(env->forces, 3 * env->natoms, MPI_DOUBLE, 0, env->comm);
  }
#endif
  if (head[0] != 1.0) {
    if (env->rank == 0)
      die(env, "step failed: %s", why);
    die(NULL, "step failed on rank 0");
  }
  env->energy = head[1];
  env->evaluated += 1;
}

void cp2k_calc_energy(force_env_t force_env) {
  cp2k_calc_energy_force(force_env);
}

void cp2k_get_potential_energy(force_env_t force_env, double *e_pot) {
  Env *env = lookup(force_env, "cp2k_get_potential_energy");
  if (!env->evaluated)
    die(env, "cp2k_get_potential_energy before cp2k_calc_energy_force");
  *e_pot = env->energy;
}

void cp2k_get_forces(force_env_t force_env, double *force, int n_el) {
  Env *env = lookup(force_env, "cp2k_get_forces");
  check_count(env, n_el, "cp2k_get_forces");
  if (!env->evaluated)
    die(env, "cp2k_get_forces before cp2k_calc_energy_force");
  memcpy(force, env->forces, (size_t)n_el * sizeof(double));
}
