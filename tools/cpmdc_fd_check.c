/**
 * Central-difference check of the forces, and optionally the stress, that a
 * session returns, after CP2K RUN_TYPE DEBUG (DEBUG_FORCES,
 * DEBUG_STRESS_TENSOR). Every evaluation goes through
 * cpmdc_session_calculate_result() with a ForceInput in Bohr and Hartree.
 *
 * Exit status: 0 every checked component agrees, 1 a mismatch, 2 a usage,
 * input, or evaluation error.
 */
#include "cpmdc.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "schema/Potentials.capnp.h"
#include <capn.h>

static double g_dx = 1.0e-3, g_max_rel_pct = 0.2, g_eps = 1.0e-5;

struct fd {
  CPMDCSession *session;
  struct ForceInput in;
  size_t n;
  int has_box;
  unsigned char *msg, *out;
  size_t msg_cap, out_cap;
  double *x0, *x, *f0; /* Bohr, Bohr, Hartree/Bohr; 3n each */
  double h0[9], s0[9]; /* Bohr; Hartree/Bohr^3 */
  int s_valid;
};

struct fd_stats {
  double max_abs, max_rel;
  int checked, total, bad;
};

static unsigned char *read_file(const char *path, size_t *size) {
  FILE *fp = fopen(path, "rb");
  unsigned char *buf = NULL;
  long n = 0;
  if (fp && fseek(fp, 0, SEEK_END) == 0 && (n = ftell(fp)) > 0 &&
      fseek(fp, 0, SEEK_SET) == 0 && (buf = malloc((size_t)n)) != NULL &&
      fread(buf, 1, (size_t)n, fp) != (size_t)n) {
    free(buf);
    buf = NULL;
  }
  if (fp)
    fclose(fp);
  *size = buf ? (size_t)n : 0;
  return buf;
}

/* Energy for positions and cell in Bohr. f0, s0, and s_valid are filled
 * when want_all is set. */
static int evaluate(struct fd *st, const double *pos, const double *box,
                    double *energy, int want_all) {
  struct capn c;
  capn_init_malloc(&c);
  capn_ptr root = capn_root(&c);
  ForceInput_ptr fp = new_ForceInput(root.seg);
  struct ForceInput v = st->in;
  v.pos = capn_new_list64(root.seg, (int)(3 * st->n));
  v.atmnrs = capn_new_list32(root.seg, (int)st->n);
  v.box = capn_new_list64(root.seg, st->has_box ? 9 : 0);
  for (size_t i = 0; i < 3 * st->n; ++i)
    capn_set64(v.pos, (int)i, capn_from_f64(pos[i]));
  for (size_t i = 0; i < st->n; ++i)
    capn_set32(v.atmnrs, (int)i, capn_get32(st->in.atmnrs, (int)i));
  for (int i = 0; st->has_box && i < 9; ++i)
    capn_set64(v.box, i, capn_from_f64(box[i]));
  v.lengthUnit = (capn_text){4, "bohr", NULL};
  v.energyUnit = (capn_text){7, "hartree", NULL};
  write_ForceInput(&v, fp);
  int len = capn_setp(root, 0, fp.p) == 0
                ? capn_write_mem(&c, st->msg, st->msg_cap, 0)
                : -1;
  capn_free(&c);
  size_t wrote = 0;
  CPMDCResult r = {0, 0.0, "ForceInput encode failed"};
  if (len > 0)
    r = cpmdc_session_calculate_result(st->session, st->msg, (size_t)len,
                                       st->out, st->out_cap, &wrote);
  if (!r.ok || capn_init_mem(&c, st->out, wrote, 0) != 0) {
    fprintf(stderr, "cpmdc-fd-check: evaluation failed: %s\n", r.message);
    return -1;
  }
  PotentialResult_ptr pr;
  struct PotentialResult res;
  pr.p = capn_getp(capn_root(&c), 0, 1);
  read_PotentialResult(&res, pr);
  capn_resolve(&res.forces.p);
  capn_resolve(&res.stress.p);
  *energy = res.energy;
  int ok = !want_all || capn_len(res.forces) == (int)(3 * st->n);
  if (want_all && ok) {
    for (size_t i = 0; i < 3 * st->n; ++i)
      st->f0[i] = capn_to_f64(capn_get64(res.forces, (int)i));
    st->s_valid = res.stress.p.type != CAPN_NULL && capn_len(res.stress) == 9;
    for (int i = 0; st->s_valid && i < 9; ++i)
      st->s0[i] = capn_to_f64(capn_get64(res.stress, i));
  }
  capn_free(&c);
  return ok ? 0 : -1;
}

/* CP2K cp2k_debug.F: percent error against the analytic value, unchecked
 * below eps. A force at or below 1e-4 gets five times the limit. */
static void compare(const char *label, double numer, double analyt,
                    int is_force, struct fd_stats *s) {
  double diff = fabs(numer - analyt), rel = 0.0, limit = g_max_rel_pct;
  int checked = fabs(analyt) >= g_eps, bad = 0;
  if (is_force && fabs(analyt) <= 1.0e-4)
    limit *= 5.0;
  if (checked) {
    rel = 100.0 * diff / fabs(analyt);
    bad = rel > limit;
    s->checked++;
  }
  s->total++;
  s->bad |= bad;
  s->max_abs = fmax(s->max_abs, diff);
  s->max_rel = fmax(s->max_rel, rel);
  printf("%-10s %20.12e %20.12e %12.4e %s\n", label, numer, analyt, diff,
         checked ? (bad ? "MISMATCH" : "ok") : "-");
}

static void summary(const char *what, const char *unit,
                    const struct fd_stats *s) {
  printf("%s max_abs_error %.6e %s max_relative_error %.6e %% checked %d of "
         "%d %s\n",
         what, s->max_abs, unit, s->max_rel, s->checked, s->total,
         s->bad ? "MISMATCH" : "PASS");
}

static double volume(const double *h) {
  return fabs(h[0] * (h[4] * h[8] - h[5] * h[7]) -
              h[1] * (h[3] * h[8] - h[5] * h[6]) +
              h[2] * (h[3] * h[7] - h[4] * h[6]));
}

static void strain(const double *e, const double *x, double *y) {
  for (int i = 0; i < 3; ++i)
    y[i] = x[i] + e[3 * i] * x[0] + e[3 * i + 1] * x[1] + e[3 * i + 2] * x[2];
}

static int check(struct fd *st, int want_stress) {
  size_t n3 = 3 * st->n;
  double e0, ep, em, h[9], vol = volume(st->h0);
  struct fd_stats fst = {0}, sst = {0};
  if (evaluate(st, st->x0, st->h0, &e0, 1) != 0)
    return 2;
  printf("dx %.3e bohr  max_relative_error %.3g %%  eps_no_error_check %.3e\n",
         g_dx, g_max_rel_pct, g_eps);
  printf("%-10s %20s %20s %12s\n", "force", "numerical", "analytical",
         "abs diff");
  memcpy(st->x, st->x0, n3 * sizeof(double));
  for (size_t i = 0; i < n3; ++i) {
    char label[32];
    st->x[i] = st->x0[i] + g_dx;
    if (evaluate(st, st->x, st->h0, &ep, 0) != 0)
      return 2;
    st->x[i] = st->x0[i] - g_dx;
    if (evaluate(st, st->x, st->h0, &em, 0) != 0)
      return 2;
    st->x[i] = st->x0[i];
    snprintf(label, sizeof(label), "%zu%c", i / 3 + 1, (int)("xyz"[i % 3]));
    compare(label, -0.5 * (ep - em) / g_dx, st->f0[i], 1, &fst);
  }
  summary("forces", "ha/bohr", &fst);
  if (!want_stress)
    return fst.bad ? 1 : 0;
  if (!st->s_valid || !st->has_box || vol <= 0.0) {
    printf("stress not computed by this calculator for this cell\n");
    return 2;
  }
  printf("%-10s %20s %20s %12s\n", "stress", "numerical", "analytical",
         "abs diff");
  for (int a = 0; a < 3; ++a) {
    for (int b = a; b < 3; ++b) {
      double e[9] = {0}, en[2];
      for (int k = 0; k < 2; ++k) {
        e[3 * a + b] = e[3 * b + a] = (k ? -g_dx : g_dx) * (a == b ? 1.0 : 0.5);
        for (int v = 0; v < 3; ++v)
          strain(e, st->h0 + 3 * v, h + 3 * v);
        for (size_t i = 0; i < st->n; ++i)
          strain(e, st->x0 + 3 * i, st->x + 3 * i);
        if (evaluate(st, st->x, h, &en[k], 0) != 0)
          return 2;
      }
      double numer = -0.5 * (en[0] - en[1]) / (g_dx * vol);
      char label[8] = {'s', "xyz"[a], "xyz"[b], 0};
      compare(label, numer, st->s0[3 * a + b], 0, &sst);
      if (a != b) {
        char swapped[8] = {'s', "xyz"[b], "xyz"[a], 0};
        compare(swapped, numer, st->s0[3 * b + a], 0, &sst);
      }
    }
  }
  summary("stress", "ha/bohr^3", &sst);
  return fst.bad || sst.bad ? 1 : 0;
}

int main(int argc, char **argv) {
  int argi = 1, want_stress = 0, rc = 2, have_arena = 0;
  for (; argi < argc && strncmp(argv[argi], "--", 2) == 0; ++argi) {
    const char *val = argi + 1 < argc ? argv[argi + 1] : NULL;
    if (strcmp(argv[argi], "--stress") == 0)
      want_stress = 1;
    else if (val && strcmp(argv[argi], "--dx") == 0)
      g_dx = atof(argv[++argi]);
    else if (val && strcmp(argv[argi], "--max-relative-error") == 0)
      g_max_rel_pct = atof(argv[++argi]);
    else if (val && strcmp(argv[argi], "--eps-no-error-check") == 0)
      g_eps = atof(argv[++argi]);
    else
      break;
  }
  if (argc - argi != 2 || !(g_dx > 0.0)) {
    fprintf(stderr,
            "usage: cpmdc-fd-check [--dx BOHR] [--max-relative-error PCT]\n"
            "         [--eps-no-error-check A.U.] [--stress]\n"
            "         PARAMS.bin FORCE_INPUT.bin\n");
    return 2;
  }
  size_t ps = 0, fs = 0;
  unsigned char *params = read_file(argv[argi], &ps);
  unsigned char *step = read_file(argv[argi + 1], &fs);
  struct fd st = {0};
  struct capn in_arena;
  if (!params || !step || capn_init_mem(&in_arena, step, fs, 0) != 0) {
    fprintf(stderr, "cpmdc-fd-check: cannot read %s or %s\n", argv[argi],
            argv[argi + 1]);
    goto done;
  }
  have_arena = 1;
  ForceInput_ptr fin;
  fin.p = capn_getp(capn_root(&in_arena), 0, 1);
  read_ForceInput(&st.in, fin);
  const char *unit = st.in.lengthUnit.str ? st.in.lengthUnit.str : "angstrom";
  double to_bohr = strcmp(unit, "bohr") == 0       ? 1.0
                   : strcmp(unit, "angstrom") == 0 ? 1.0 / 0.529177210903
                                                   : 0.0;
  st.n = (size_t)capn_len(st.in.atmnrs);
  st.has_box = st.in.box.p.type != CAPN_NULL && capn_len(st.in.box) == 9;
  if (to_bohr == 0.0 || st.n == 0 || capn_len(st.in.pos) != (int)(3 * st.n)) {
    fprintf(stderr, "cpmdc-fd-check: need angstrom or bohr positions\n");
    goto done;
  }
  st.x0 = malloc(9 * st.n * sizeof(double));
  st.session = cpmdc_session_create(params, ps);
  st.out_cap = cpmdc_potential_result_size_for_force_input(step, fs);
  st.out = malloc(st.out_cap + 1);
  st.msg_cap = 1024 + 32 * st.n;
  st.msg = malloc(st.msg_cap);
  if (!st.x0 || !st.session || !st.out || st.out_cap == 0 || !st.msg) {
    fprintf(stderr, "cpmdc-fd-check: setup failed: %s\n", cpmdc_last_error());
    goto done;
  }
  st.x = st.x0 + 3 * st.n;
  st.f0 = st.x + 3 * st.n;
  for (size_t i = 0; i < 3 * st.n; ++i)
    st.x0[i] = capn_to_f64(capn_get64(st.in.pos, (int)i)) * to_bohr;
  for (int i = 0; st.has_box && i < 9; ++i)
    st.h0[i] = capn_to_f64(capn_get64(st.in.box, i)) * to_bohr;
  rc = check(&st, want_stress);
done:
  cpmdc_session_destroy(st.session);
  if (have_arena)
    capn_free(&in_arena);
  free(st.x0);
  free(st.out);
  free(st.msg);
  free(params);
  free(step);
  return rc;
}
