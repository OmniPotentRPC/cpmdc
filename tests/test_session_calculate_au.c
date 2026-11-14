/**
 * cpmdc_session_calculate_au gives the bits of cpmdc_session_calculate_result
 * for a ForceInput in bohr and hartree with the same geometry.
 *
 * Each path runs a fresh session in a setup child of its own over the same
 * steps, so both start from the same process state. Energies, forces, and
 * stresses are compared with memcmp; the virial against |det(cell)| times
 * the PotentialResult stress.
 */
#define _GNU_SOURCE
#include "cpmdc.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cmocka.h>

#include "setup_child.h"

#include "schema/Potentials.capnp.h"
#include <capn.h>

#define N_ATOMS 3
#define MAX_STEPS 4

static const char *g_cluster = NULL;
static const char *g_periodic = NULL;

static const int k_numbers[N_ATOMS] = {8, 1, 1};

/* Water in Bohr. The last step repeats the one before, which the session
 * answers from its stored result. */
static const double k_steps[MAX_STEPS][3 * N_ATOMS] = {
    {0.0, 0.0, 0.2216669, 0.0, 1.4308949, -0.8866678, 0.0, -1.4308949,
     -0.8866678},
    {0.03, -0.02, 0.25, 0.01, 1.47, -0.86, -0.02, -1.41, -0.93},
    {-0.01, 0.04, 0.19, 0.05, 1.39, -0.92, 0.0, -1.46, -0.84},
    {-0.01, 0.04, 0.19, 0.05, 1.39, -0.92, 0.0, -1.46, -0.84},
};

/* Lattice vectors as rows, Bohr: triclinic, cubic, compressed cubic. */
static const double k_cells[MAX_STEPS][9] = {
    {17.0, 0.0, 0.0, 2.3, 16.1, 0.0, -1.3, 1.7, 17.8},
    {18.9, 0.0, 0.0, 0.0, 18.9, 0.0, 0.0, 0.0, 18.9},
    {17.955, 0.0, 0.0, 0.0, 17.955, 0.0, 0.0, 0.0, 17.955},
    {17.955, 0.0, 0.0, 0.0, 17.955, 0.0, 0.0, 0.0, 17.955},
};

static const double k_zero_cell[9] = {0};

enum cell_mode { CELL_NONE, CELL_ZERO, CELL_STEPS };

struct au_case {
  const char *params_path;
  enum cell_mode cell;
  int steps;
};

struct au_step {
  int ok;
  double energy;
  double forces[3 * N_ATOMS];
  int stress_valid;
  double stress[9];
  double virial[9];
  char message[512];
};

struct au_out {
  struct au_step step[MAX_STEPS];
};

static unsigned char *read_file(const char *path, size_t *size) {
  FILE *fp = fopen(path, "rb");
  if (!fp)
    return NULL;
  if (fseek(fp, 0, SEEK_END) != 0) {
    fclose(fp);
    return NULL;
  }
  long n = ftell(fp);
  if (n <= 0) {
    fclose(fp);
    return NULL;
  }
  rewind(fp);
  unsigned char *buf = (unsigned char *)malloc((size_t)n);
  if (buf && fread(buf, 1, (size_t)n, fp) != (size_t)n) {
    free(buf);
    buf = NULL;
  }
  fclose(fp);
  *size = (size_t)n;
  return buf;
}

static CPMDCSession *open_session(const char *path) {
  size_t size = 0;
  unsigned char *params = read_file(path, &size);
  if (!params)
    return NULL;
  CPMDCSession *session = cpmdc_session_create(params, size);
  free(params);
  return session;
}

static const double *step_cell(const struct au_case *c, int i) {
  switch (c->cell) {
  case CELL_ZERO:
    return k_zero_cell;
  case CELL_STEPS:
    return k_cells[i];
  default:
    return NULL;
  }
}

static double det3(const double *m) {
  return m[0] * (m[4] * m[8] - m[5] * m[7]) -
         m[1] * (m[3] * m[8] - m[5] * m[6]) +
         m[2] * (m[3] * m[7] - m[4] * m[6]);
}

static int encode_force_input(const double *pos, const double *cell,
                              unsigned char *out, size_t cap) {
  struct capn c;
  capn_init_malloc(&c);
  capn_ptr root = capn_root(&c);
  ForceInput_ptr fp = new_ForceInput(root.seg);
  struct ForceInput v;
  memset(&v, 0, sizeof(v));
  v.pos = capn_new_list64(root.seg, 3 * N_ATOMS);
  v.atmnrs = capn_new_list32(root.seg, N_ATOMS);
  v.box = capn_new_list64(root.seg, cell ? 9 : 0);
  for (int i = 0; i < 3 * N_ATOMS; ++i)
    capn_set64(v.pos, i, capn_from_f64(pos[i]));
  for (int i = 0; i < N_ATOMS; ++i)
    capn_set32(v.atmnrs, i, (uint32_t)k_numbers[i]);
  for (int i = 0; cell && i < 9; ++i)
    capn_set64(v.box, i, capn_from_f64(cell[i]));
  v.lengthUnit = (capn_text){4, "bohr", NULL};
  v.energyUnit = (capn_text){7, "hartree", NULL};
  write_ForceInput(&v, fp);
  int len = capn_setp(root, 0, fp.p) == 0 ? capn_write_mem(&c, out, cap, 0) : -1;
  capn_free(&c);
  return len;
}

static void capnp_child(void *vctx, void *vout) {
  const struct au_case *c = (const struct au_case *)vctx;
  struct au_out *o = (struct au_out *)vout;
  CPMDCSession *session = open_session(c->params_path);
  if (!session)
    return;
  static unsigned char msg[4096];
  static unsigned char out[8192];
  for (int i = 0; i < c->steps; ++i) {
    struct au_step *s = &o->step[i];
    int len = encode_force_input(k_steps[i], step_cell(c, i), msg,
                                 sizeof(msg));
    if (len <= 0)
      return;
    size_t wrote = 0;
    CPMDCResult r = cpmdc_session_calculate_result(session, msg, (size_t)len,
                                                   out, sizeof(out), &wrote);
    snprintf(s->message, sizeof(s->message), "%s", r.message);
    if (!r.ok)
      return;
    struct capn arena;
    if (capn_init_mem(&arena, out, wrote, 0) != 0)
      return;
    PotentialResult_ptr pr;
    pr.p = capn_getp(capn_root(&arena), 0, 1);
    struct PotentialResult view;
    read_PotentialResult(&view, pr);
    capn_resolve(&view.forces.p);
    capn_resolve(&view.stress.p);
    s->energy = view.energy;
    for (int k = 0; k < 3 * N_ATOMS; ++k)
      s->forces[k] = capn_to_f64(capn_get64(view.forces, k));
    s->stress_valid =
        view.stress.p.type != CAPN_NULL && capn_len(view.stress) == 9;
    for (int k = 0; s->stress_valid && k < 9; ++k)
      s->stress[k] = capn_to_f64(capn_get64(view.stress, k));
    capn_free(&arena);
    s->ok = 1;
  }
  cpmdc_session_destroy(session);
}

static void au_child(void *vctx, void *vout) {
  const struct au_case *c = (const struct au_case *)vctx;
  struct au_out *o = (struct au_out *)vout;
  CPMDCSession *session = open_session(c->params_path);
  if (!session)
    return;
  for (int i = 0; i < c->steps; ++i) {
    struct au_step *s = &o->step[i];
    int rc = cpmdc_session_calculate_au(
        session, N_ATOMS, k_numbers, k_steps[i], step_cell(c, i),
        &s->energy, s->forces, s->virial);
    snprintf(s->message, sizeof(s->message), "%s", cpmdc_last_error());
    if (rc != 0)
      return;
    CPMDCStressTensor stress;
    memset(&stress, 0, sizeof(stress));
    s->stress_valid =
        cpmdc_session_last_stress(session, &stress) == 0 && stress.valid;
    memcpy(s->stress, stress.values, sizeof(s->stress));
    s->ok = 1;
  }
  cpmdc_session_destroy(session);
}

/* Returns how many steps carried a stress. */
static int compare_paths(const struct au_case *c) {
  static struct au_out via_capnp;
  static struct au_out via_au;
  RUN_SETUP_CHILD(capnp_child, (void *)c, &via_capnp);
  RUN_SETUP_CHILD(au_child, (void *)c, &via_au);
  int stressed = 0;
  for (int i = 0; i < c->steps; ++i) {
    const struct au_step *a = &via_capnp.step[i];
    const struct au_step *b = &via_au.step[i];
    if (!a->ok || !b->ok)
      fail_msg("step %d: capnp ok=%d (%s), au ok=%d (%s)", i, a->ok,
               a->message, b->ok, b->message);
    assert_memory_equal(&a->energy, &b->energy, sizeof(double));
    assert_memory_equal(a->forces, b->forces, sizeof(a->forces));
    assert_int_equal(a->stress_valid, b->stress_valid);
    double expected[9] = {0};
    const double *cell = step_cell(c, i);
    if (a->stress_valid) {
      assert_memory_equal(a->stress, b->stress, sizeof(a->stress));
      if (cell) {
        double volume = fabs(det3(cell));
        int nonzero = 0;
        for (int k = 0; k < 9; ++k) {
          expected[k] = volume * a->stress[k];
          nonzero |= expected[k] != 0.0;
        }
        assert_true(nonzero);
      }
      ++stressed;
    }
    assert_memory_equal(expected, b->virial, sizeof(expected));
  }
  return stressed;
}

static void test_cluster_matches_capnp(void **state) {
  (void)state;
  struct au_case c = {g_cluster, CELL_NONE, MAX_STEPS};
  assert_int_equal(compare_paths(&c), 0);
}

static void test_zero_cell_matches_capnp(void **state) {
  (void)state;
  struct au_case c = {g_cluster, CELL_ZERO, 2};
  assert_int_equal(compare_paths(&c), 0);
}

static void test_periodic_matches_capnp(void **state) {
  (void)state;
  struct au_case c = {g_periodic, CELL_STEPS, MAX_STEPS};
  assert_int_equal(compare_paths(&c), MAX_STEPS);
}

struct errors_out {
  int null_session;
  char null_session_text[128];
  int no_atoms;
  char no_atoms_text[128];
  int no_energy;
  int no_forces;
  int first;
  int virial_optional;
  int changed_count;
  char changed_count_text[128];
  int changed_species;
};

static void errors_child(void *vctx, void *vout) {
  const char *path = (const char *)vctx;
  struct errors_out *o = (struct errors_out *)vout;
  double energy = 0.0;
  double forces[3 * N_ATOMS];
  const double *pos = k_steps[0];
  o->null_session = cpmdc_session_calculate_au(NULL, N_ATOMS, k_numbers, pos,
                                               NULL, &energy, forces, NULL);
  snprintf(o->null_session_text, sizeof(o->null_session_text), "%s",
           cpmdc_last_error());
  CPMDCSession *session = open_session(path);
  if (!session)
    return;
  o->no_atoms = cpmdc_session_calculate_au(session, 0, k_numbers, pos, NULL,
                                           &energy, forces, NULL);
  snprintf(o->no_atoms_text, sizeof(o->no_atoms_text), "%s",
           cpmdc_last_error());
  o->no_energy = cpmdc_session_calculate_au(session, N_ATOMS, k_numbers, pos,
                                            NULL, NULL, forces, NULL);
  o->no_forces = cpmdc_session_calculate_au(session, N_ATOMS, k_numbers, pos,
                                            NULL, &energy, NULL, NULL);
  o->first = cpmdc_session_calculate_au(session, N_ATOMS, k_numbers, pos, NULL,
                                        &energy, forces, NULL);
  o->virial_optional = cpmdc_last_error()[0] == '\0';
  o->changed_count = cpmdc_session_calculate_au(
      session, N_ATOMS - 1, k_numbers, k_steps[1], NULL, &energy,
      forces, NULL);
  snprintf(o->changed_count_text, sizeof(o->changed_count_text), "%s",
           cpmdc_last_error());
  static const int other[N_ATOMS] = {8, 8, 1};
  o->changed_species = cpmdc_session_calculate_au(
      session, N_ATOMS, other, k_steps[1], NULL, &energy, forces,
      NULL);
  cpmdc_session_destroy(session);
}

static void test_errors_return_minus_one(void **state) {
  (void)state;
  static struct errors_out o;
  RUN_SETUP_CHILD(errors_child, (void *)g_cluster, &o);
  assert_int_equal(o.null_session, -1);
  assert_string_equal(o.null_session_text, "null session");
  assert_int_equal(o.no_atoms, -1);
  assert_string_equal(o.no_atoms_text, "invalid arguments");
  assert_int_equal(o.no_energy, -1);
  assert_int_equal(o.no_forces, -1);
  assert_int_equal(o.first, 0);
  assert_true(o.virial_optional);
  assert_int_equal(o.changed_count, -1);
  assert_string_equal(o.changed_count_text,
                      "topology change requires a new session");
  assert_int_equal(o.changed_species, -1);
}

int main(int argc, char **argv) {
  if (argc != 3) {
    fprintf(stderr, "usage: %s CLUSTER_PARAMS PERIODIC_PARAMS\n", argv[0]);
    return 2;
  }
  g_cluster = argv[1];
  g_periodic = argv[2];
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(test_cluster_matches_capnp),
      cmocka_unit_test(test_zero_cell_matches_capnp),
      cmocka_unit_test(test_periodic_matches_capnp),
      cmocka_unit_test(test_errors_return_minus_one),
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
