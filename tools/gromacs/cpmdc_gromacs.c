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
 * The input is either the file of grompp -qmi with
 * qmmm-cp2k-qmmethod = INPUT, or GROMACS' generated PBE/BLYP input with the
 * CPMDParams message named by CPMDC_GROMACS_PARAMS. The generated
 * functional, charge, and multiplicity must match that message. On either
 * path the charge and multiplicity compared are the first CHARGE and
 * MULTIPLICITY lines of the file in any section, &CPMDC included, because
 * those are the values GROMACS reads. GROMACS
 * writes the input at mdrun with the name of the PDB file it generates.
 * That PDB gives every atom's element, whether it is QM or MM, and the
 * classical charge of each MM atom (QM atoms carry 0: grompp removes their
 * charges). An explicit &CPMDC section takes precedence and contains:
 *
 *   PARAMS path           CPMDParams message (required)
 *   CELL PERIODIC|NONE    send the GROMACS box to cpmdc, or the deck's cell
 *   EMBEDDING MECHANICAL|NONE   QM-MM Coulomb from fixed charges, or none
 *   QM_CHARGES q ...      classical charge of each QM atom, in index order
 *   EWALD_RCUT r          real-space cutoff of that Coulomb sum, nm
 *   EWALD_PRECISION p     erfc(alpha rcut) and the reciprocal-space cutoff
 *   EPSILON_R r           relative dielectric constant of that Coulomb sum
 *   QM_RANKS ALL|1        under an MPI GROMACS, run the QM calculation on
 *                         every rank of the communicator (default) or on
 *                         rank 0 alone
 *
 * Mechanical embedding: the QM-MM Coulomb energy of the QM_CHARGES and the
 * MM charges, an Ewald sum when the box is periodic and a direct sum when
 * it is not. GROMACS keeps the QM-MM Lennard-Jones and the bonded terms
 * that cross the boundary; the QM density is not polarized by the MM
 * charges.
 *
 * Link atoms: each &LINK section (QM_INDEX, MM_INDEX, optional ALPHA_IMOMM,
 * LINK_TYPE IMOMM) caps one QM-MM covalent bond with a hydrogen atom that is
 * part of the QM calculation only, at the IMOMM position of CP2K,
 *   x_H = x_QM + (x_MM - x_QM) / ALPHA_IMOMM,
 * ALPHA_IMOMM 1.38 unless given. GROMACS writes these sections in the input
 * it generates; on the explicit path they sit beside &CPMDC. The force on
 * each hydrogen goes to the two atoms of its bond by the chain rule, 1 - a
 * to the QM atom and a to the MM atom, a = 1 / ALPHA_IMOMM. With EMBEDDING
 * MECHANICAL the MM atom of each bond leaves the QM-MM Coulomb sum (its
 * charge is deleted), since it sits one bond from a QM charge.
 *
 * Virial: cp2k_get_stress_tensor gives the strain derivative of the energy
 * returned, for a symmetric strain that moves every atom and the box:
 * cpmdc's stress with CELL PERIODIC, the symmetric part of sum x (x) f over
 * the QM atoms with CELL NONE, and the strain derivative of the Ewald sum.
 * It is unavailable without a box, and with CELL PERIODIC when cpmdc
 * computes no stress (an isolated deck).
 *
 * Positions: GROMACS wraps every atom into the box. QM atoms are unwrapped
 * against their previous positions, so each follows a continuous path; with
 * CELL NONE the first step also makes the QM region whole about its first
 * atom.
 *
 * MPI: an MPI GROMACS calls cp2k_create_force_env_comm on each of its
 * particle-particle ranks with the simulation communicator, and gives every
 * rank the positions of all atoms. Rank 0 reads the files and broadcasts
 * the configuration, every QM rank evaluates rank 0's positions and cell on
 * one CPMD calculator (cpmdc_adopt_calculator_comm), and rank 0's energy
 * and forces go to every rank. GROMACS refuses separate PME ranks with
 * QM/MM, so every rank of that communicator makes each call.
 *
 * Errors cannot be returned through this API: the library prints
 * "cpmdc-gromacs: ..." to stderr and to the CP2K output file GROMACS names,
 * and exits with status 1 (MPI_Abort under MPI).
 */
#define _GNU_SOURCE
#include "libcp2k.h"

#include "cpmdc.h"

#include <complex.h>
#include <ctype.h>
#include <limits.h>
#ifdef __ELF__
#include <dlfcn.h>
#endif
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
/* CP2K's default ALPHA_IMOMM: the ratio of a C-C to a C-H bond length. */
#define LINK_ALPHA_DEFAULT 1.38

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
enum { QM_RANKS_ALL = 0, QM_RANKS_ONE = 1 };

typedef struct {
  int used;
  int rank, size;
#ifdef CPMDC_GROMACS_MPI
  MPI_Comm comm;
  /* GROMACS initialized MPI, so OpenCPMD does not set up a communicator of
   * its own and needs one handed over, even on one rank. */
  int host_mpi;
#endif
  FILE *log;
  /* every rank that runs the QM calculation */
  char *params;
  size_t params_len;
  CPMDCSession *session;
  int adopted;
  int cell_mode, qm_ranks;
  int generated_input;
  int nqm;
  /* QM-MM bonds capped with a link hydrogen; the QM calculation sees
   * ncalc = nqm + nlink atoms, the hydrogens last. */
  int nlink, ncalc;
  /* per link: the QM and MM atoms as GROMACS indices, the QM atom as a slot of
   * qm_index, and a = 1 / ALPHA_IMOMM */
  int *link_qm, *link_mm, *link_slot;
  double *link_scale;
  int *qm_index, *qm_z;
  double *qm_x, *qm_prev, *qm_f;
  int have_prev;
  /* rank 0 */
  int embedding, saw_epsilon_r;
  double ewald_rcut_nm, ewald_precision, epsilon_r;
  int nmm, qm_charges_listed;
  int *mm_index;
  double *qm_charge, *mm_charge;
  /* every rank */
  int natoms;
  double *pos, *forces;
  double rows[3][3];
  int have_positions, have_cell;
  double energy;
  /* -dE/d(strain) of that energy, Hartree, and whether it is known */
  double virial[3][3];
  int have_virial;
  int evaluated, stress_requests;
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
  if (up && env && env->host_mpi)
    MPI_Abort(env->comm, 1);
#endif
  exit(EXIT_FAILURE);
}

static Env *lookup(force_env_t id, const char *call) {
  if (id < 1 || id > MAX_ENVS || !g_envs[id - 1].used)
    die(NULL, "%s: no force environment %d", call, id);
  return &g_envs[id - 1];
}

static void *grow(Env *env, void *old, size_t count, size_t size) {
  void *p = realloc(old, count * size);
  if (!p)
    die(env, "out of memory");
  return p;
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
  for (char *tok = strtok_r(line, " \t\r", &save); tok;
       tok = strtok_r(NULL, " \t\r", &save)) {
    if (n >= max)
      return -1;
    out[n++] = tok;
  }
  return n;
}

static double number(Env *env, const char *text, const char *key) {
  char *end = NULL;
  double v = strtod(text, &end);
  if (end == text || *end != '\0' || !isfinite(v))
    die(env, "&CPMDC %s: %s is not a number", key, text);
  return v;
}

static int integer(Env *env, const char *text, const char *key) {
  char *end = NULL;
  long value = strtol(text, &end, 10);
  if (end == text || *end != '\0' || value < INT_MIN || value > INT_MAX)
    die(env, "%s value %s is not an integer", key, text);
  return (int)value;
}

static void copy_token_checked(Env *env, char *dst, size_t cap,
                               const char *src, const char *field) {
  size_t n = strlen(src);
  if (n >= cap)
    die(env, "%s token is too long (%zu bytes, max %zu)", field, n, cap - 1);
  memcpy(dst, src, n + 1);
}

static int is_section(const char *word, const char *name) {
  return word[0] == '&' && strcasecmp(word + 1, name) == 0;
}

static int in_section(char sections[][32], int depth, const char *name) {
  for (int i = 0; i < depth; ++i)
    if (strcasecmp(sections[i], name) == 0)
      return 1;
  return 0;
}

static void push_section(Env *env, char sections[][32], int *depth, const char *word) {
  size_t n = strlen(word + 1);
  if (n == 0)
    die(env, "empty section name %s", word);
  if (n >= sizeof(sections[0]))
    die(env, "section name %s is too long", word);
  if (*depth >= 32)
    die(env, "too many nested sections in generated input");
  memcpy(sections[*depth], word + 1, n + 1);
  *depth += 1;
}

static void pop_section(Env *env, char sections[][32], int *depth, char **w, int n) {
  if (*depth == 0)
    die(env, "unexpected &END without an open section");
  if (n > 2)
    die(env, "&END takes at most one section name");
  if (n == 2 && strcasecmp(sections[*depth - 1], w[1]) != 0)
    die(env, "&END %s closes %s", w[1], sections[*depth - 1]);
  *depth -= 1;
}

static void parse_input(Env *env, const char *path, char *pdb_path, size_t pdb_cap) {
  char *text = read_file(env, path, NULL, "the QM/MM input");
  int in_cpmdc = 0, saw_cpmdc = 0, nq = 0, capq = 64;
  char sections[32][32] = {{0}};
  int depth = 0;
  int saw_charge = 0, saw_multiplicity = 0, saw_functional = 0;
  int input_charge = 0, input_multiplicity = 0;
  char input_charge_text[64] = {0};
  char input_multiplicity_text[64] = {0};
  const char *input_charge_key = "&DFT CHARGE";
  const char *input_multiplicity_key = "&DFT MULTIPLICITY";
  char generated_functional[64] = {0};
  char xc_error[256] = {0};
  char *params_path = NULL;
  int link_cap = 0, link_open = 0, link_qm = 0, link_mm = 0;
  double link_alpha = LINK_ALPHA_DEFAULT;
  pdb_path[0] = '\0';
  env->cell_mode = CELL_DEFAULT;
  env->embedding = EMBED_UNSET;
  env->ewald_rcut_nm = 0.0;
  env->ewald_precision = 1e-8;
  env->epsilon_r = 1.0;
  env->saw_epsilon_r = 0;
  env->qm_ranks = QM_RANKS_ALL;
  env->qm_charge = alloc(env, (size_t)capq, sizeof(double));
  char *save = NULL;
  for (char *line = strtok_r(text, "\n", &save); line; line = strtok_r(NULL, "\n", &save)) {
    char *w[256];
    int n = words(line, w, 256);
    if (n < 0)
      die(env, "%s has a line with too many tokens", path);
    if (n == 0)
      continue;
    if (strcasecmp(w[0], "COORD_FILE_NAME") == 0) {
      if (n != 2)
        die(env, "COORD_FILE_NAME takes one value");
      copy_token_checked(env, pdb_path, pdb_cap, w[1], "COORD_FILE_NAME");
      continue;
    }
    if (!in_cpmdc && strcasecmp(w[0], "&END") == 0) {
      if (link_open && depth > 0 && strcasecmp(sections[depth - 1], "LINK") == 0) {
        if (link_qm == 0 || link_mm == 0)
          die(env, "%s has a &LINK section without %s", path,
              link_qm == 0 ? "QM_INDEX" : "MM_INDEX");
        if (env->nlink == link_cap) {
          link_cap = link_cap ? 2 * link_cap : 8;
          env->link_qm = grow(env, env->link_qm, (size_t)link_cap, sizeof(int));
          env->link_mm = grow(env, env->link_mm, (size_t)link_cap, sizeof(int));
          env->link_scale = grow(env, env->link_scale, (size_t)link_cap, sizeof(double));
        }
        env->link_qm[env->nlink] = link_qm - 1;
        env->link_mm[env->nlink] = link_mm - 1;
        env->link_scale[env->nlink++] = 1.0 / link_alpha;
        link_open = 0;
      }
      pop_section(env, sections, &depth, w, n);
      continue;
    }
    if (!in_cpmdc && link_open && w[0][0] != '&') {
      if (n != 2)
        die(env, "&LINK %s takes one value", w[0]);
      if (strcasecmp(w[0], "QM_INDEX") == 0 || strcasecmp(w[0], "MM_INDEX") == 0) {
        const int index = integer(env, w[1], w[0]);
        if (index < 1)
          die(env, "&LINK %s is a 1-based atom index, not %s", w[0], w[1]);
        if (toupper((unsigned char)w[0][0]) == 'Q')
          link_qm = index;
        else
          link_mm = index;
      } else if (strcasecmp(w[0], "ALPHA_IMOMM") == 0) {
        link_alpha = number(env, w[1], "&LINK ALPHA_IMOMM");
        if (!(link_alpha > 1.0))
          die(env, "&LINK ALPHA_IMOMM must exceed 1 so the hydrogen lies between the "
                   "bonded atoms, not %s",
              w[1]);
      } else if (strcasecmp(w[0], "LINK_TYPE") == 0) {
        if (strcasecmp(w[1], "IMOMM") != 0)
          die(env, "&LINK LINK_TYPE is IMOMM, the only scheme cpmdc-gromacs implements, "
                   "not %s",
              w[1]);
      } else {
        die(env, "unsupported &LINK keyword %s: QM_INDEX, MM_INDEX, ALPHA_IMOMM, and "
                 "LINK_TYPE IMOMM are read",
            w[0]);
      }
      continue;
    }
    if (!in_cpmdc && w[0][0] == '&') {
      if (is_section(w[0], "LINK")) {
        if (n != 1)
          die(env, "&LINK takes no value");
        if (link_open)
          die(env, "&LINK sections do not nest");
        link_open = 1;
        link_qm = link_mm = 0;
        link_alpha = LINK_ALPHA_DEFAULT;
      }
      if (is_section(w[0], "CPMDC")) {
        if (n != 1)
          die(env, "&CPMDC takes no value");
        in_cpmdc = 1;
        saw_cpmdc = 1;
        continue;
      }
      if (is_section(w[0], "XC_FUNCTIONAL") && in_section(sections, depth, "DFT")) {
        if (n != 2) {
          if (xc_error[0] == '\0')
            snprintf(xc_error, sizeof(xc_error),
                     "&XC_FUNCTIONAL in &DFT takes one value");
        } else {
          size_t flen = strlen(w[1]);
          if (flen >= sizeof(generated_functional)) {
            if (xc_error[0] == '\0')
              snprintf(xc_error, sizeof(xc_error),
                       "&XC_FUNCTIONAL token is too long (%zu bytes, max %zu)",
                       flen, sizeof(generated_functional) - 1);
          } else {
            memcpy(generated_functional, w[1], flen + 1);
            saw_functional = 1;
          }
        }
      }
      push_section(env, sections, &depth, w[0]);
      continue;
    }
    if (strcasecmp(w[0], "CHARGE") == 0) {
      if (n != 2)
        die(env, "CHARGE takes one value");
      if (!saw_charge) {
        copy_token_checked(env, input_charge_text, sizeof(input_charge_text),
                           w[1], "CHARGE");
        input_charge_key = in_section(sections, depth, "DFT") ? "&DFT CHARGE" : "input CHARGE";
      }
      saw_charge = 1;
      continue;
    }
    if (strcasecmp(w[0], "MULTIPLICITY") == 0) {
      if (n != 2)
        die(env, "MULTIPLICITY takes one value");
      if (!saw_multiplicity) {
        copy_token_checked(env, input_multiplicity_text,
                           sizeof(input_multiplicity_text), w[1],
                           "MULTIPLICITY");
        input_multiplicity_key = in_section(sections, depth, "DFT")
                                     ? "&DFT MULTIPLICITY"
                                     : "input MULTIPLICITY";
      }
      saw_multiplicity = 1;
      continue;
    }
    if (!in_cpmdc) {
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
    } else if (strcasecmp(w[0], "QM_RANKS") == 0) {
      if (strcasecmp(w[1], "ALL") == 0)
        env->qm_ranks = QM_RANKS_ALL;
      else if (strcmp(w[1], "1") == 0)
        env->qm_ranks = QM_RANKS_ONE;
      else
        die(env, "&CPMDC QM_RANKS is ALL or 1, not %s", w[1]);
    } else if (strcasecmp(w[0], "EWALD_RCUT") == 0) {
      env->ewald_rcut_nm = number(env, w[1], "EWALD_RCUT");
      if (env->ewald_rcut_nm <= 0.0)
        die(env, "&CPMDC EWALD_RCUT must be positive");
    } else if (strcasecmp(w[0], "EWALD_PRECISION") == 0) {
      env->ewald_precision = number(env, w[1], "EWALD_PRECISION");
      if (!(env->ewald_precision > 0.0 && env->ewald_precision < 1.0))
        die(env, "&CPMDC EWALD_PRECISION must lie between 0 and 1");
    } else if (strcasecmp(w[0], "EPSILON_R") == 0) {
      env->epsilon_r = number(env, w[1], "EPSILON_R");
      if (!(env->epsilon_r > 0.0))
        die(env, "&CPMDC EPSILON_R must be positive");
      env->saw_epsilon_r = 1;
    } else {
      die(env, "unknown &CPMDC keyword %s", w[0]);
    }
  }
  if (!in_cpmdc && depth != 0)
    die(env, "%s ended with an open section &%s", path, sections[depth - 1]);
  free(text);
  if (!saw_cpmdc) {
    const char *fallback = getenv("CPMDC_GROMACS_PARAMS");
    if (!fallback || !*fallback)
      die(env, "%s has no &CPMDC section and CPMDC_GROMACS_PARAMS is not set. "
               "Set it to a CPMDParams message, or use "
               "qmmm-cp2k-qmmethod = INPUT with an explicit &CPMDC section",
          path);
    if (xc_error[0] != '\0')
      die(env, "%s", xc_error);
    if (!saw_functional)
      die(env, "%s has no &XC_FUNCTIONAL value", path);
    if (!saw_charge)
      die(env, "%s has no CHARGE value", path);
    if (!saw_multiplicity)
      die(env, "%s has no MULTIPLICITY value", path);
    params_path = strdup(fallback);
    if (!params_path)
      die(env, "out of memory");
    env->generated_input = 1;
  } else if (!params_path) {
    die(env, "&CPMDC needs PARAMS, the path of a CPMDParams message");
  }
  if (pdb_path[0] == '\0')
    die(env, "%s names no COORD_FILE_NAME", path);
  env->params = read_file(env, params_path, &env->params_len, "the CPMDParams message");
  if (env->params_len == 0)
    die(env, "the CPMDParams message %s is empty", params_path);
  /* GROMACS scans the external input as a stream of tokens and keeps the first
   * CHARGE and MULTIPLICITY it sees, in generated input and with -qmi alike.
   * It removes the classical charge of the QM atoms by that CHARGE. The method
   * cpmdc runs is the CPMDParams message, so a value that the message does not
   * share would run another system than the one GROMACS set up. */
  if (env->generated_input || saw_charge || saw_multiplicity) {
    const char *input_name = env->generated_input ? "generated input" : "&DFT";
    const char *deck_name = env->generated_input ? "CPMDC_GROMACS_PARAMS" : "&CPMDC PARAMS";
    CPMDCMethodInfo method;
    if (cpmdc_params_method(env->params, env->params_len, &method) != 0) {
      const char *why = cpmdc_last_error();
      if (!why || !*why)
        die(env, "cannot read %s deck %s", deck_name, params_path);
      die(env, "cannot read %s deck %s: %s", deck_name, params_path, why);
    }
    if (env->generated_input &&
        strcasecmp(generated_functional, method.functional) != 0)
      die(env, "%s functional %s disagrees with %s functional %s", input_name,
          generated_functional, deck_name, method.functional);
    if (saw_charge) {
      input_charge = integer(env, input_charge_text,
                             env->generated_input ? "generated input CHARGE"
                                                  : input_charge_key);
      if (input_charge != method.charge)
        die(env, "%s CHARGE %d disagrees with %s charge %d", input_name,
            input_charge, deck_name, method.charge);
    }
    if (saw_multiplicity) {
      input_multiplicity =
          integer(env, input_multiplicity_text,
                  env->generated_input ? "generated input MULTIPLICITY"
                                       : input_multiplicity_key);
      if (input_multiplicity != method.multiplicity)
        die(env, "%s MULTIPLICITY %d disagrees with %s multiplicity %d",
            input_name, input_multiplicity, deck_name, method.multiplicity);
    }
  }
  free(params_path);
  env->qm_charges_listed = nq;
}

/* Match each &LINK to its QM atom and its MM atom, and append the link
 * hydrogens to the atoms the QM calculation sees. */
static void resolve_links(Env *env, const int *is_qm) {
  env->ncalc = env->nqm + env->nlink;
  env->link_slot = alloc(env, (size_t)env->nlink, sizeof(int));
  int *slot_of = alloc(env, (size_t)env->natoms, sizeof(int));
  for (int i = 0; i < env->nqm; ++i)
    slot_of[env->qm_index[i]] = i + 1;
  for (int k = 0; k < env->nlink; ++k) {
    const int q = env->link_qm[k], m = env->link_mm[k];
    if (q >= env->natoms || m >= env->natoms)
      die(env, "&LINK names atom %d, but the PDB file has %d atoms",
          (q >= env->natoms ? q : m) + 1, env->natoms);
    if (!is_qm[q])
      die(env, "&LINK QM_INDEX %d is an MM atom in the PDB file", q + 1);
    if (is_qm[m])
      die(env, "&LINK MM_INDEX %d is a QM atom in the PDB file", m + 1);
    for (int j = 0; j < k; ++j)
      if (env->link_qm[j] == q && env->link_mm[j] == m)
        die(env, "&LINK QM_INDEX %d MM_INDEX %d appears twice", q + 1, m + 1);
    env->link_slot[k] = slot_of[q] - 1;
  }
  free(slot_of);
  int *grown = realloc(env->qm_z, (size_t)env->ncalc * sizeof(int));
  if (!grown)
    die(env, "out of memory");
  env->qm_z = grown;
  for (int k = 0; k < env->nlink; ++k)
    env->qm_z[env->nqm + k] = 1;
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
    if (nw < 0)
      die(env, "%s: atom line %d has too many tokens", path, n + 1);
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
    if (end == w[nw - 1] || *end != '\0' || !isfinite(q[n]))
      die(env, "%s: atom %d has no finite charge", path, n + 1);
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
  resolve_links(env, is_qm);
  free(is_qm);
  free(z);
  free(q);
  if (env->cell_mode == CELL_DEFAULT)
    env->cell_mode = nmm == 0 ? CELL_PERIODIC : CELL_NONE;
  if (env->embedding == EMBED_UNSET) {
    if (charged_mm && env->generated_input)
      die(env,
          "generated GROMACS input has charged MM atoms but no QM_CHARGES; "
          "GROMACS writes zero instead of the QM topology charges. Use an "
          "explicit &CPMDC section with EMBEDDING MECHANICAL and QM_CHARGES, "
          "or select EMBEDDING NONE");
    if (charged_mm)
      die(env, "the MM atoms carry charges: set &CPMDC EMBEDDING MECHANICAL (with "
               "QM_CHARGES) or EMBEDDING NONE");
    env->embedding = EMBED_NONE;
  }
  if (env->embedding == EMBED_MECHANICAL) {
    if (nq != nqm)
      die(env, "&CPMDC QM_CHARGES lists %d charges for %d QM atoms", nq, nqm);
    /* The MM atom of a cut bond is one bond from a QM charge, where a point
     * charge model is not meaningful: delete its charge, as the Z1 scheme. */
    for (int k = 0; k < env->nlink; ++k)
      for (int j = 0; j < nmm; ++j)
        if (env->mm_index[j] == env->link_mm[k])
          env->mm_charge[j] = 0.0;
  } else {
    if (nq != 0)
      die(env, "&CPMDC QM_CHARGES needs EMBEDDING MECHANICAL");
    if (env->saw_epsilon_r)
      die(env, "&CPMDC EPSILON_R needs EMBEDDING MECHANICAL");
  }
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
 * image (Ewald), Hartree with positions in Bohr; forces are added, and so
 * is -dE/d(strain) to vir. The real-space sum takes every image pair closer
 * than rc, the reciprocal sum every k with |k| <= kmax, and the
 * neutralizing background its QM-MM part, -pi Q_QM Q_MM / (V alpha^2).
 * Under a strain the pair vectors and 1/k follow the cell and the structure
 * factors do not change: a pair term adds g d_a d_b, a k term
 * E_k (delta_ab - 2 k_a k_b (1/(4 alpha^2) + 1/k^2)), and the background,
 * which goes as 1/V, E_bg delta_ab. */
static double ewald_qm_mm(Env *env, double r[3][3], double rc, double tol,
                          const double *xq, double *fq, const double *xm, double *fm,
                          double vir[3][3], double *alpha_out, double *kmax_out) {
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
              for (int l = 0; l < 3; ++l)
                vir[k][l] += g * d[k] * d[l];
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
        const double ek = pref * f * creal(sq * conj(sm));
        energy += ek;
        const double c = 2.0 * (inv4a2 + 1.0 / k2);
        for (int a = 0; a < 3; ++a)
          for (int b = 0; b < 3; ++b)
            vir[a][b] += ek * ((a == b ? 1.0 : 0.0) - c * kv[a] * kv[b]);
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
  const double background = -PI * sum_q * sum_m / (volume * alpha * alpha);
  energy += background;
  for (int a = 0; a < 3; ++a)
    vir[a][a] += background;
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

/* Ranks that take part in the QM calculation: all of them, or rank 0. */
static int qm_rank(const Env *env) { return env->rank == 0 || env->qm_ranks == QM_RANKS_ALL; }

static void open_session(Env *env) {
#ifdef CPMDC_GROMACS_MPI
  if (env->host_mpi && !env->adopted) {
    MPI_Comm calc = env->qm_ranks == QM_RANKS_ALL ? env->comm : MPI_COMM_SELF;
    int ranks = env->qm_ranks == QM_RANKS_ALL ? env->size : 1;
    if (cpmdc_adopt_calculator_comm(&calc, sizeof(calc), ranks) != 0 && cpmdc_available()) {
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
#ifdef __ELF__
  /* OpenCPMD's BLAS calls bind to the first definition in the process. A
   * GROMACS built with its bundled BLAS and LAPACK exports dgemm_ and the
   * rest ahead of the BLAS libcpmdc links, and OpenCPMD's nonlocal
   * pseudopotential forces come out wrong. */
  Dl_info info;
  void *dgemm = dlsym(RTLD_DEFAULT, "dgemm_");
  if (cpmdc_available() && dgemm && dladdr(dgemm, &info) && info.dli_fname &&
      strstr(info.dli_fname, "libgromacs"))
    die(env,
        "dgemm_ resolves to %s: GROMACS was built with its own BLAS and LAPACK, which take "
        "the place of the BLAS OpenCPMD links and give wrong forces. Rebuild GROMACS with "
        "-DGMX_EXTERNAL_BLAS=ON -DGMX_EXTERNAL_LAPACK=ON against that BLAS",
        info.dli_fname);
#endif
}

/* Every QM rank: the QM energy and forces, from the same positions and cell
 * on each. Rank 0 then adds the embedding and fills the forces of every atom,
 * Hartree and Hartree/Bohr, and the virial -dE/d(strain) of the energy,
 * Hartree, for a strain that moves every atom and the GROMACS box with it.
 * With CELL PERIODIC the QM part is cpmdc's virial, which needs the stress
 * tensor of a periodic deck; with CELL NONE the QM energy does not depend on
 * the box, and its part is sum_i x_i (x) f_i over the QM atoms at the
 * positions cpmdc evaluated. */
static int evaluate(Env *env, char *why, size_t why_cap) {
  const int nq = env->nqm, nc = env->ncalc;
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
  for (int k = 0; k < env->nlink; ++k) {
    const double *xq = &env->qm_x[3 * env->link_slot[k]];
    double xm[3], *xh = &env->qm_x[3 * (nq + k)];
    memcpy(xm, &env->pos[3 * (size_t)env->link_mm[k]], sizeof(xm));
    if (periodic)
      nearest_image(env->rows, inv, xq, xm);
    for (int c = 0; c < 3; ++c)
      xh[c] = xq[c] + env->link_scale[k] * (xm[c] - xq[c]);
  }
  if (!env->session)
    open_session(env);
  double qm_energy = 0.0, qm_virial[9] = {0.0};
  if (cpmdc_session_calculate_au(env->session, nc, env->qm_x, env->qm_z,
                                 env->cell_mode == CELL_PERIODIC ? &env->rows[0][0] : NULL,
                                 &qm_energy, env->qm_f, qm_virial) != 0) {
    const char *text = cpmdc_last_error();
    snprintf(why, why_cap, "%s", text && *text ? text : "cpmdc_session_calculate_au failed");
    return -1;
  }
  memcpy(env->qm_prev, env->qm_x, 3 * (size_t)nq * sizeof(double));
  env->have_prev = 1;
  if (env->rank != 0) {
    env->energy = qm_energy;
    return 0;
  }

  memset(env->forces, 0, 3 * (size_t)env->natoms * sizeof(double));
  memset(env->virial, 0, sizeof(env->virial));
  env->have_virial = periodic;
  if (env->cell_mode == CELL_PERIODIC) {
    CPMDCStressTensor stress;
    if (cpmdc_session_last_stress(env->session, &stress) == 0 && stress.valid) {
      for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b)
          env->virial[a][b] = qm_virial[3 * a + b];
    } else {
      env->have_virial = 0;
    }
  } else {
    /* cpmdc removes the net force of an isolated deck, so this sum does not
     * depend on the origin of x unless CPMDC_RECENTER=off. The cluster grid
     * still leaves a torque, so sum x (x) f is not symmetric; a symmetric
     * strain sees only its symmetric part. */
    for (int i = 0; i < nc; ++i)
      for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b)
          env->virial[a][b] += 0.5 * (env->qm_x[3 * i + a] * env->qm_f[3 * i + b] +
                                      env->qm_x[3 * i + b] * env->qm_f[3 * i + a]);
  }
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
      double qmmm_vir[3][3] = {{0.0}};
      heights(env->rows, h);
      rc = env->ewald_rcut_nm > 0.0 ? env->ewald_rcut_nm / BOHR_TO_NM
                                    : fmin(1.0 / BOHR_TO_NM, 0.5 * fmin(h[0], fmin(h[1], h[2])));
      coulomb = ewald_qm_mm(env, env->rows, rc, env->ewald_precision, env->qm_x, fq, xm, fm,
                            qmmm_vir, &alpha, &kmax);
      if (env->epsilon_r != 1.0) {
        const double inv_eps = 1.0 / env->epsilon_r;
        coulomb *= inv_eps;
        for (int i = 0; i < 3 * nq; ++i)
          fq[i] *= inv_eps;
        for (int j = 0; j < 3 * nm; ++j)
          fm[j] *= inv_eps;
        for (int a = 0; a < 3; ++a)
          for (int b = 0; b < 3; ++b)
            qmmm_vir[a][b] *= inv_eps;
      }
      for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b)
          env->virial[a][b] += qmmm_vir[a][b];
      if (env->log && env->evaluated == 0)
        fprintf(env->log, "ewald rcut_bohr %.17g alpha_per_bohr %.17g kmax_per_bohr %.17g\n",
                rc, alpha, kmax);
    } else {
      coulomb = direct_qm_mm(env, env->qm_x, fq, xm, fm);
      if (env->epsilon_r != 1.0) {
        const double inv_eps = 1.0 / env->epsilon_r;
        coulomb *= inv_eps;
        for (int i = 0; i < 3 * nq; ++i)
          fq[i] *= inv_eps;
        for (int j = 0; j < 3 * nm; ++j)
          fm[j] *= inv_eps;
      }
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
  /* x_H = (1 - a) x_QM + a x_MM, so the force on a link hydrogen goes to the
   * two atoms of its bond in those shares. */
  for (int k = 0; k < env->nlink; ++k) {
    const double a = env->link_scale[k], *fh = &env->qm_f[3 * (nq + k)];
    double *fq = &env->forces[3 * (size_t)env->qm_index[env->link_slot[k]]];
    double *fm = &env->forces[3 * (size_t)env->link_mm[k]];
    for (int c = 0; c < 3; ++c) {
      fq[c] += (1.0 - a) * fh[c];
      fm[c] += a * fh[c];
    }
  }
  env->energy = qm_energy + coulomb;
  if (env->log && env->evaluated == 0) {
    fprintf(env->log,
            "first step: backend %s qm_energy_hartree %.17g coulomb_qm_mm_hartree %.17g\n",
            cpmdc_available() ? "opencpmd" : "reference", qm_energy, coulomb);
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

#ifdef CPMDC_GROMACS_MPI
static void bcast(Env *env, void *buf, int count, MPI_Datatype type) {
  if (count > 0 && MPI_Bcast(buf, count, type, 0, env->comm) != MPI_SUCCESS)
    die(env, "MPI_Bcast failed");
}

/* Rank 0 read the input and the PDB file; the other ranks get what the QM
 * calculation needs. GROMACS writes both files on its main rank just before
 * this call, so only rank 0 may read them. */
static void share_config(Env *env) {
  long head[6] = {env->natoms, env->nqm, env->cell_mode, env->qm_ranks,
                  (long)env->params_len, env->nlink};
  bcast(env, head, 6, MPI_LONG);
  if (env->rank != 0) {
    env->natoms = (int)head[0];
    env->nqm = (int)head[1];
    env->cell_mode = (int)head[2];
    env->qm_ranks = (int)head[3];
    env->params_len = (size_t)head[4];
    env->nlink = (int)head[5];
    env->ncalc = env->nqm + env->nlink;
  }
  if (env->qm_ranks != QM_RANKS_ALL)
    return;
  if (env->rank != 0) {
    env->qm_index = alloc(env, (size_t)env->nqm, sizeof(int));
    env->qm_z = alloc(env, (size_t)env->ncalc, sizeof(int));
    env->link_mm = alloc(env, (size_t)env->nlink, sizeof(int));
    env->link_slot = alloc(env, (size_t)env->nlink, sizeof(int));
    env->link_scale = alloc(env, (size_t)env->nlink, sizeof(double));
    env->params = alloc(env, env->params_len + 1, 1);
  }
  bcast(env, env->qm_index, env->nqm, MPI_INT);
  bcast(env, env->qm_z, env->ncalc, MPI_INT);
  bcast(env, env->link_mm, env->nlink, MPI_INT);
  bcast(env, env->link_slot, env->nlink, MPI_INT);
  bcast(env, env->link_scale, env->nlink, MPI_DOUBLE);
  for (size_t off = 0; off < env->params_len; off += 1 << 30) {
    size_t n = env->params_len - off < (size_t)1 << 30 ? env->params_len - off : (size_t)1 << 30;
    bcast(env, env->params + off, (int)n, MPI_BYTE);
  }
}
#endif

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
  if (env->rank == 0) {
    if (output && *output) {
      env->log = fopen(output, "w");
      if (env->log)
        fprintf(env->log, PROG " on %s, input %s\n", cpmdc_version(), input);
    }
    char pdb[4096];
    parse_input(env, input, pdb, sizeof(pdb));
    parse_pdb(env, pdb);
    if (env->log) {
      fprintf(env->log,
              "atoms %d qm %d mm %d cell %s embedding %s ewald_precision %.3g ranks %d "
              "qm_ranks %d epsilon_r %.17g links %d\n",
              env->natoms, env->nqm, env->nmm,
              env->cell_mode == CELL_PERIODIC ? "PERIODIC" : "NONE",
              env->embedding == EMBED_MECHANICAL ? "MECHANICAL" : "NONE",
              env->ewald_precision, env->size,
              env->qm_ranks == QM_RANKS_ALL ? env->size : 1,
              env->epsilon_r, env->nlink);
      fflush(env->log);
    }
  }
#ifdef CPMDC_GROMACS_MPI
  if (env->size > 1)
    share_config(env);
#endif
  if (qm_rank(env)) {
    env->qm_x = alloc(env, 3 * (size_t)env->ncalc, sizeof(double));
    env->qm_prev = alloc(env, 3 * (size_t)env->ncalc, sizeof(double));
    env->qm_f = alloc(env, 3 * (size_t)env->ncalc, sizeof(double));
  }
  env->pos = alloc(env, 3 * (size_t)env->natoms, sizeof(double));
  env->forces = alloc(env, 3 * (size_t)env->natoms, sizeof(double));
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
  MPI_Comm given = MPI_Comm_f2c((MPI_Fint)mpi_comm);
  MPI_Comm_rank(given, &proto.rank);
  MPI_Comm_size(given, &proto.size);
  proto.host_mpi = 1;
  if (MPI_Comm_dup(given, &proto.comm) != MPI_SUCCESS) {
    proto.comm = given;
    die(&proto, "MPI_Comm_dup failed");
  }
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
#ifdef CPMDC_GROMACS_MPI
  int finalized = 0;
  MPI_Finalized(&finalized);
  if (env->host_mpi && !finalized)
    MPI_Comm_free(&env->comm);
#endif
  free(env->params);
  free(env->qm_index);
  free(env->qm_z);
  free(env->mm_index);
  free(env->link_qm);
  free(env->link_mm);
  free(env->link_slot);
  free(env->link_scale);
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

/* GROMACS 2024 and 2025 send box_d[3 i + j] = box[j][i]; GROMACS 2026 sends
 * box_d[3 i + j] = box[i][j], one lattice vector after the other, which is
 * CP2K's layout. A GROMACS box is lower triangular (box[0][1], box[0][2]
 * and box[1][2] are zero), so the first layout has zeros in elements 3, 6
 * and 7 and the second in 1, 2 and 5; a rectangular box reads the same
 * either way. */
void cp2k_set_cell(force_env_t force_env, const double *new_cell) {
  Env *env = lookup(force_env, "cp2k_set_cell");
  const int transposed = new_cell[3] == 0.0 && new_cell[6] == 0.0 && new_cell[7] == 0.0;
  for (int a = 0; a < 3; ++a)
    for (int k = 0; k < 3; ++k)
      env->rows[a][k] = transposed ? new_cell[3 * k + a] : new_cell[3 * a + k];
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
#ifdef CPMDC_GROMACS_MPI
  /* Every QM rank evaluates rank 0's positions and cell, bit for bit. */
  if (env->size > 1 && env->qm_ranks == QM_RANKS_ALL) {
    bcast(env, env->pos, 3 * env->natoms, MPI_DOUBLE);
    bcast(env, &env->rows[0][0], 9, MPI_DOUBLE);
  }
#endif
  int any = 0;
  for (int a = 0; a < 3; ++a)
    for (int k = 0; k < 3; ++k)
      any |= env->rows[a][k] != 0.0;
  env->have_cell = any;
  char why[1024] = {0};
  /* ok, energy, have_virial, virial */
  double head[12] = {0.0};
  if (qm_rank(env)) {
    if (evaluate(env, why, sizeof(why)) == 0)
      head[0] = 1.0;
    head[1] = env->energy;
    head[2] = env->have_virial;
    memcpy(&head[3], &env->virial[0][0], 9 * sizeof(double));
  }
#ifdef CPMDC_GROMACS_MPI
  /* Rank 0's result is the result: CPMD's other ranks need not agree to the
   * last bit (see the MPI model in the cpmdc documentation). */
  if (env->size > 1) {
    bcast(env, head, 12, MPI_DOUBLE);
    if (head[0] == 1.0)
      bcast(env, env->forces, 3 * env->natoms, MPI_DOUBLE);
  }
#endif
  if (head[0] != 1.0) {
    if (env->rank == 0)
      die(env, "step failed: %s", why);
    die(env, "step failed on rank 0");
  }
  env->energy = head[1];
  env->have_virial = head[2] == 1.0;
  memcpy(&env->virial[0][0], &head[3], 9 * sizeof(double));
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

/* CP2K's contract: the configurational stress, pressure-positive, so that
 * stress times the cell volume is the virial -dE/d(strain); Hartree/Bohr^3,
 * column-major. available is 0, with a zero tensor, when it is not known. */
void cp2k_get_stress_tensor(force_env_t force_env, double *stress_tensor, int *available) {
  Env *env = lookup(force_env, "cp2k_get_stress_tensor");
  if (!env->evaluated)
    die(env, "cp2k_get_stress_tensor before cp2k_calc_energy_force");
  const double volume = env->have_cell ? fabs(det3(env->rows)) : 0.0;
  const int known = env->have_virial && volume > 0.0;
  for (int a = 0; a < 3; ++a)
    for (int b = 0; b < 3; ++b)
      stress_tensor[a + 3 * b] = known ? env->virial[a][b] / volume : 0.0;
  *available = known;
  if (env->log && env->stress_requests++ == 0) {
    if (known)
      fprintf(env->log, "first stress request: available 1 pressure_hartree_per_bohr3 %.17g\n",
              (stress_tensor[0] + stress_tensor[4] + stress_tensor[8]) / 3.0);
    else
      fprintf(env->log, "first stress request: available 0 (%s)\n",
              volume <= 0.0 ? "no GROMACS box"
                            : "cpmdc computed no stress for this cell: the deck is isolated, "
                              "or CPMDC_STRESS=0");
    fflush(env->log);
  }
}
