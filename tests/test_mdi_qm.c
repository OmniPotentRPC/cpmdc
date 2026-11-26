/**
 * cpmdc-mdi against the command sequence of LAMMPS fix mdi/qm.
 *
 * The command state is linked without libmdi. The test checks the
 * atomic-unit step, the cell origin, and the stress tensor LAMMPS reads
 * back as a virial. The reference evaluator supplies the numbers.
 */
#define _GNU_SOURCE
#include "cpmdc.h"
#include "cpmdc_mdi.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <setjmp.h>
#include <stdarg.h>
#include <cmocka.h>

#define N 3

static const char *g_isolated = NULL;
static const char *g_periodic = NULL;

/* Water in Bohr, and a restricted triclinic cell, one vector per row. */
static const int k_z[N] = {8, 1, 1};
static const double k_pos[3 * N] = {0.0, 0.0, 0.2216669, 0.0, 1.4308949,
                                    -0.8866678, 0.0, -1.4308949, -0.8866678};
static const double k_cell[9] = {17.0, 0.0, 0.0, 2.3, 16.1, 0.0, -1.3, 1.7, 17.8};
/* Powers of two larger than the coordinates, so subtracting them is exact. */
static const double k_displ[3] = {4.0, -4.0, 8.0};

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
  *size = buf ? (size_t)n : 0;
  return buf;
}

static cpmdc_mdi *open_engine(const char *path) {
  size_t n = 0;
  unsigned char *buf = read_file(path, &n);
  cpmdc_mdi *eng = buf ? cpmdc_mdi_create(buf, n) : NULL;
  free(buf);
  return eng;
}

static int direct(const char *path, const double *pos, const int *z, int n,
                  const double *cell, double *energy, double *forces,
                  CPMDCStressTensor *stress) {
  size_t nbytes = 0;
  unsigned char *buf = read_file(path, &nbytes);
  CPMDCSession *session = buf ? cpmdc_session_create(buf, nbytes) : NULL;
  free(buf);
  if (!session)
    return -1;
  int rc = cpmdc_session_calculate_au(session, n, pos, z, cell, energy, forces,
                                      NULL);
  if (rc == 0) {
    if (cpmdc_session_last_stress(session, stress) != 0)
      stress->valid = 0;
  }
  cpmdc_session_destroy(session);
  return rc;
}

static void add_displ(double *out, const double *pos, int n,
                      const double *displ) {
  for (int i = 0; i < n; i++) {
    out[3 * i] = pos[3 * i] + displ[0];
    out[3 * i + 1] = pos[3 * i + 1] + displ[1];
    out[3 * i + 2] = pos[3 * i + 2] + displ[2];
  }
}

/* LAMMPS fix mdi/qm reduces the nine MDI components to xx yy zz xy xz yz. */
static void lammps_symmetric(const double m[9], double out[6]) {
  out[0] = m[0];
  out[1] = m[4];
  out[2] = m[8];
  out[3] = 0.5 * (m[1] + m[3]);
  out[4] = 0.5 * (m[2] + m[6]);
  out[5] = 0.5 * (m[5] + m[7]);
}

static void require_ok(cpmdc_mdi *eng, const char *command, const void *in,
                       void *out) {
  int rc = cpmdc_mdi_command(eng, command, in, out);
  if (rc != 0)
    fprintf(stderr, "%s: %s\n", command, cpmdc_mdi_error(eng));
  assert_int_equal(rc, 0);
}

static void test_periodic_step_matches_and_serves_one_evaluation(void **state) {
  (void)state;
  setenv("CPMDC_STRESS", "1", 1);
  double energy = 0.0, forces[3 * N], wire_e = 0.0, wire_pe = 0.0, wire_ke = 1.0;
  double wire_f[3 * N], wire_s[9], sent[3 * N], recovered[3 * N];
  CPMDCStressTensor stress;
  memset(&stress, 0, sizeof stress);
  /* (r + d) - d is the position the engine must evaluate. For these
   * coordinates that is a few ulps away from r, so the reference call uses
   * the recovered position and the comparison stays bit-exact. */
  add_displ(sent, k_pos, N, k_displ);
  for (int i = 0; i < N; i++) {
    recovered[3 * i] = sent[3 * i] - k_displ[0];
    recovered[3 * i + 1] = sent[3 * i + 1] - k_displ[1];
    recovered[3 * i + 2] = sent[3 * i + 2] - k_displ[2];
    assert_double_equal(recovered[3 * i], k_pos[3 * i], 1e-9);
    assert_double_equal(recovered[3 * i + 1], k_pos[3 * i + 1], 1e-9);
    assert_double_equal(recovered[3 * i + 2], k_pos[3 * i + 2], 1e-9);
  }
  assert_int_equal(direct(g_periodic, recovered, k_z, N, k_cell, &energy,
                          forces, &stress),
                   0);
  assert_int_equal(stress.valid, 1);

  cpmdc_mdi *eng = open_engine(g_periodic);
  assert_non_null(eng);
  int n = N;
  require_ok(eng, ">NATOMS", &n, NULL);
  require_ok(eng, ">CELL_DISPL", k_displ, NULL);
  require_ok(eng, ">CELL", k_cell, NULL);
  require_ok(eng, ">ELEMENTS", k_z, NULL);
  require_ok(eng, ">COORDS", sent, NULL);
  assert_int_equal(cpmdc_mdi_evaluations(eng), 0);
  require_ok(eng, "<PE", NULL, &wire_pe);
  require_ok(eng, "<KE_ELEC", NULL, &wire_ke);
  require_ok(eng, "<ENERGY", NULL, &wire_e);
  require_ok(eng, "<FORCES", NULL, wire_f);
  require_ok(eng, "<STRESS", NULL, wire_s);
  assert_int_equal(cpmdc_mdi_evaluations(eng), 1);
  assert_double_equal(wire_ke, 0.0, 0.0);
  assert_memory_equal(&wire_pe, &energy, sizeof energy);
  assert_memory_equal(&wire_e, &energy, sizeof energy);
  assert_memory_equal(wire_f, forces, sizeof forces);
  assert_memory_equal(wire_s, stress.values, sizeof stress.values);

  double sym_wire[6], sym_ref[6];
  lammps_symmetric(wire_s, sym_wire);
  lammps_symmetric(stress.values, sym_ref);
  assert_memory_equal(sym_wire, sym_ref, sizeof sym_wire);

  sent[0] += 0.01;
  require_ok(eng, ">COORDS", sent, NULL);
  require_ok(eng, "<FORCES", NULL, wire_f);
  assert_int_equal(cpmdc_mdi_evaluations(eng), 2);
  cpmdc_mdi_destroy(eng);
}

static void test_isolated_stress_on_the_wire_is_zero(void **state) {
  (void)state;
  setenv("CPMDC_STRESS", "1", 1);
  double energy = 0.0, forces[3 * N], wire_e = 0.0, wire_f[3 * N], wire_s[9];
  CPMDCStressTensor stress;
  memset(&stress, 0, sizeof stress);
  memset(wire_s, 1, sizeof wire_s);
  assert_int_equal(
      direct(g_isolated, k_pos, k_z, N, k_cell, &energy, forces, &stress), 0);
  assert_int_equal(stress.valid, 0);

  cpmdc_mdi *eng = open_engine(g_isolated);
  assert_non_null(eng);
  int n = N;
  require_ok(eng, ">NATOMS", &n, NULL);
  require_ok(eng, ">CELL", k_cell, NULL);
  require_ok(eng, ">TYPES", k_z, NULL);
  require_ok(eng, ">COORDS", k_pos, NULL);
  require_ok(eng, "<ENERGY", NULL, &wire_e);
  require_ok(eng, "<FORCES", NULL, wire_f);
  require_ok(eng, "<STRESS", NULL, wire_s);
  assert_memory_equal(&wire_e, &energy, sizeof energy);
  assert_memory_equal(wire_f, forces, sizeof forces);
  for (int i = 0; i < 9; i++)
    assert_double_equal(wire_s[i], 0.0, 0.0);
  cpmdc_mdi_destroy(eng);
}

static void test_cluster_has_no_cell(void **state) {
  (void)state;
  double energy = 0.0, forces[3 * N], wire_e = 0.0, wire_f[3 * N];
  CPMDCStressTensor stress;
  memset(&stress, 0, sizeof stress);
  assert_int_equal(
      direct(g_isolated, k_pos, k_z, N, NULL, &energy, forces, &stress), 0);

  cpmdc_mdi *eng = open_engine(g_isolated);
  assert_non_null(eng);
  int n = N;
  require_ok(eng, ">NATOMS", &n, NULL);
  require_ok(eng, ">ELEMENTS", k_z, NULL);
  require_ok(eng, ">COORDS", k_pos, NULL);
  require_ok(eng, "<ENERGY", NULL, &wire_e);
  require_ok(eng, "<FORCES", NULL, wire_f);
  assert_memory_equal(&wire_e, &energy, sizeof energy);
  assert_memory_equal(wire_f, forces, sizeof forces);
  cpmdc_mdi_destroy(eng);
}

static void test_composition_change_reopens_the_session(void **state) {
  (void)state;
  cpmdc_mdi *eng = open_engine(g_periodic);
  assert_non_null(eng);
  int n = N;
  double first = 0.0, second = 0.0;
  require_ok(eng, ">NATOMS", &n, NULL);
  require_ok(eng, ">CELL", k_cell, NULL);
  require_ok(eng, ">ELEMENTS", k_z, NULL);
  require_ok(eng, ">COORDS", k_pos, NULL);
  require_ok(eng, "<ENERGY", NULL, &first);

  int hehh[N] = {2, 1, 1};
  require_ok(eng, ">NATOMS", &n, NULL);
  require_ok(eng, ">CELL", k_cell, NULL);
  require_ok(eng, ">ELEMENTS", hehh, NULL);
  require_ok(eng, ">COORDS", k_pos, NULL);
  require_ok(eng, "<ENERGY", NULL, &second);
  assert_int_equal(cpmdc_mdi_evaluations(eng), 2);
  /* Different elements, so the energies are not the same bits. */
  assert_true(memcmp(&first, &second, sizeof first) != 0);
  cpmdc_mdi_destroy(eng);
}

/* Dyadic triclinic cell, one vector per row, so a lattice shift is exact. */
static const double k_wrap_cell[9] = {16.0, 0.0, 0.0, 2.0, 16.0, 0.0,
                                     -2.0, 2.0, 16.0};
static const double k_wrap_pos[3 * N] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0,
                                        7.0, 8.0, 9.0};

static void open_wrap(cpmdc_mdi *eng, const double *pos) {
  int n = N;
  double sent[3 * N];
  require_ok(eng, ">NATOMS", &n, NULL);
  require_ok(eng, ">CELL_DISPL", k_displ, NULL);
  require_ok(eng, ">CELL", k_wrap_cell, NULL);
  require_ok(eng, ">ELEMENTS", k_z, NULL);
  add_displ(sent, pos, N, k_displ);
  require_ok(eng, ">COORDS", sent, NULL);
}

static void read_energy_forces(cpmdc_mdi *eng, double *energy, double *forces) {
  require_ok(eng, "<ENERGY", NULL, energy);
  require_ok(eng, "<FORCES", NULL, forces);
}

static void send_coords(cpmdc_mdi *eng, const double *pos) {
  double sent[3 * N];
  add_displ(sent, pos, N, k_displ);
  require_ok(eng, ">COORDS", sent, NULL);
}

static void copy_shift_atom0(double *dst, const double *src,
                             const double add[3]) {
  memcpy(dst, src, sizeof(double) * 3 * N);
  dst[0] += add[0];
  dst[1] += add[1];
  dst[2] += add[2];
}

/* Integer combination of the cell rows. The entries are dyadic, so the
 * shift is exact. */
static void lattice_image(double out[3], const int coeff[3]) {
  out[0] = out[1] = out[2] = 0.0;
  for (int row = 0; row < 3; row++) {
    out[0] += coeff[row] * k_wrap_cell[3 * row];
    out[1] += coeff[row] * k_wrap_cell[3 * row + 1];
    out[2] += coeff[row] * k_wrap_cell[3 * row + 2];
  }
}

/* fix mdi/qm remaps atom->x into the box. A frame that crosses a face is
 * an integer lattice image away from the continuous path. Each row is one
 * periodic direction, taken in both signs, and one corner uses all three
 * rows at once. The engine must recover that path, so the energy and
 * forces match an unwrapped run, and the next in-box step stays on it. */
static void test_wrapped_crossing_matches_unwrapped_run(void **state) {
  (void)state;
  setenv("CPMDC_STRESS", "0", 1);
  const double step[3] = {0.5, 0.0, 0.0};
  const int images[][3] = {
      {1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0},
      {0, 0, 1}, {0, 0, -1}, {1, 1, 1},
  };
  const int nimages = (int)(sizeof images / sizeof images[0]);

  for (int image = 0; image < nimages; image++) {
    double shift[3], step1[3 * N], step2[3 * N], wrapped[3 * N];
    double bare_e = 0.0, bare_f[3 * N], jump_e = 0.0, jump_f[3 * N];
    double plain_e = 0.0, plain_f[3 * N], wrap_e = 0.0, wrap_f[3 * N];
    double plain2 = 0.0, plain2_f[3 * N], wrap2 = 0.0, wrap2_f[3 * N];
    CPMDCStressTensor stress;
    memset(&stress, 0, sizeof stress);
    lattice_image(shift, images[image]);
    copy_shift_atom0(step1, k_wrap_pos, step);
    copy_shift_atom0(step2, step1, step);
    copy_shift_atom0(wrapped, step1, shift);

    assert_int_equal(direct(g_periodic, step1, k_z, N, k_wrap_cell, &bare_e,
                            bare_f, &stress),
                     0);
    assert_int_equal(direct(g_periodic, wrapped, k_z, N, k_wrap_cell, &jump_e,
                            jump_f, &stress),
                     0);
    if (memcmp(&bare_e, &jump_e, sizeof bare_e) == 0)
      fprintf(stderr, "image %d %d %d did not move the energy\n",
              images[image][0], images[image][1], images[image][2]);
    assert_true(memcmp(&bare_e, &jump_e, sizeof bare_e) != 0);

    cpmdc_mdi *plain = open_engine(g_periodic);
    cpmdc_mdi *crossed = open_engine(g_periodic);
    assert_non_null(plain);
    assert_non_null(crossed);
    open_wrap(plain, k_wrap_pos);
    open_wrap(crossed, k_wrap_pos);
    read_energy_forces(plain, &plain_e, plain_f);
    read_energy_forces(crossed, &wrap_e, wrap_f);
    assert_memory_equal(&plain_e, &wrap_e, sizeof plain_e);
    assert_memory_equal(plain_f, wrap_f, sizeof plain_f);

    send_coords(plain, step1);
    send_coords(crossed, wrapped);
    read_energy_forces(plain, &plain_e, plain_f);
    read_energy_forces(crossed, &wrap_e, wrap_f);
    if (memcmp(&wrap_e, &bare_e, sizeof bare_e) != 0)
      fprintf(stderr, "image %d %d %d was not unwrapped\n", images[image][0],
              images[image][1], images[image][2]);
    assert_memory_equal(&plain_e, &bare_e, sizeof bare_e);
    assert_memory_equal(plain_f, bare_f, sizeof bare_f);
    assert_memory_equal(&wrap_e, &bare_e, sizeof bare_e);
    assert_memory_equal(wrap_f, bare_f, sizeof bare_f);
    assert_int_equal(cpmdc_mdi_evaluations(plain), 2);
    assert_int_equal(cpmdc_mdi_evaluations(crossed), 2);

    assert_int_equal(direct(g_periodic, step2, k_z, N, k_wrap_cell, &bare_e,
                            bare_f, &stress),
                     0);
    send_coords(plain, step2);
    send_coords(crossed, step2);
    read_energy_forces(plain, &plain2, plain2_f);
    read_energy_forces(crossed, &wrap2, wrap2_f);
    assert_memory_equal(&plain2, &bare_e, sizeof bare_e);
    assert_memory_equal(plain2_f, bare_f, sizeof bare_f);
    assert_memory_equal(&wrap2, &bare_e, sizeof bare_e);
    assert_memory_equal(wrap2_f, bare_f, sizeof bare_f);
    cpmdc_mdi_destroy(plain);
    cpmdc_mdi_destroy(crossed);
  }
}

static void test_missing_inputs_and_unknown_commands_fail(void **state) {
  (void)state;
  cpmdc_mdi *eng = open_engine(g_isolated);
  assert_non_null(eng);
  int direction = -1, count = -1;
  assert_int_equal(cpmdc_mdi_payload(eng, ">COORDS", &direction, &count), -1);
  assert_non_null(strstr(cpmdc_mdi_error(eng), ">NATOMS"));
  double energy = 0.0;
  assert_int_equal(cpmdc_mdi_command(eng, "<ENERGY", NULL, &energy), -1);
  assert_non_null(strstr(cpmdc_mdi_error(eng), ">NATOMS"));

  int n = N;
  require_ok(eng, ">NATOMS", &n, NULL);
  require_ok(eng, ">ELEMENTS", k_z, NULL);
  assert_int_equal(cpmdc_mdi_command(eng, "<PE", NULL, &energy), -1);
  assert_non_null(strstr(cpmdc_mdi_error(eng), ">COORDS"));

  assert_int_equal(cpmdc_mdi_payload(eng, ">COORDS", &direction, &count), 0);
  assert_int_equal(direction, CPMDC_MDI_RECV_DOUBLE);
  assert_int_equal(count, 3 * N);
  assert_int_equal(cpmdc_mdi_payload(eng, "<STRESS", &direction, &count), 0);
  assert_int_equal(direction, CPMDC_MDI_SEND_DOUBLE);
  assert_int_equal(count, 9);
  assert_int_equal(cpmdc_mdi_payload(eng, "MD", &direction, &count), -1);

  char name[16];
  memset(name, 1, sizeof name);
  require_ok(eng, "<@", NULL, name);
  assert_string_equal(name, "@DEFAULT");
  require_ok(eng, "EXIT", NULL, NULL);
  assert_int_equal(cpmdc_mdi_exited(eng), 1);
  cpmdc_mdi_destroy(eng);
}

int main(int argc, char **argv) {
  if (argc != 3) {
    fprintf(stderr, "usage: %s ISOLATED_PARAMS PERIODIC_PARAMS\n", argv[0]);
    return 2;
  }
  g_isolated = argv[1];
  g_periodic = argv[2];
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(test_periodic_step_matches_and_serves_one_evaluation),
      cmocka_unit_test(test_isolated_stress_on_the_wire_is_zero),
      cmocka_unit_test(test_cluster_has_no_cell),
      cmocka_unit_test(test_composition_change_reopens_the_session),
      cmocka_unit_test(test_wrapped_crossing_matches_unwrapped_run),
      cmocka_unit_test(test_missing_inputs_and_unknown_commands_fail),
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
