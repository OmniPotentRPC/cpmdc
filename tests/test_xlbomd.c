/**
 * XL-BOMD auxiliary propagation, and a session evaluation on the reference
 * engine.
 *
 * The coefficient rows are Table I of Niklasson, Steneteg, Odell, Bock,
 * Challacombe, Tymczak, Holmstrom, Zheng, and Weber, J. Chem. Phys. 130,
 * 214109 (2009). The reversible map is Niklasson, Phys. Rev. Lett. 100,
 * 123004 (2008). The water drift is the nuclear Verlet energy of a harmonic
 * stand-in. The session test calls cpmdc_session_calculate_au().
 */
#include "cpmdc.h"
#include "cpmdc_params.h"
#include "cpmdc_restart.h"

#include <math.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdint.h>
#include <setjmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cmocka.h>

enum { ASPC_ORDER = 4, XL_ORDER = 5, WATER_STEPS = 12, REV_STEPS = 16 };

static const char *g_aspc_path;
static const char *g_xl_path;
static const char *g_bad_path;
static const char *g_golden_path;

static const double k_bohr0[9] = {
    9.44863062391448, 9.44863062391448, 9.44863062391448,
    9.44863062391448, 11.05489782997994, 8.40928125528389,
    9.44863062391448, 7.84236341784902, 8.40928125528389};

static const double k_r0 = 1.8088;
static const double k_theta0 = 104.52 * 3.14159265358979323846 / 180.0;
static const double k_kr = 0.40;
static const double k_kth = 0.12;
static const double k_mass_h = 1.00784 * 1822.888486;
static const double k_mass_o = 15.999 * 1822.888486;
static const double k_dt = 20.0;
static const double k_beta = 0.35;
static const double k_scf_tol = 1.0e-12;
static const int k_scf_max = 400;

static const int k_paper_c[7][10] = {
    {-2, 3, 0, -1, 0, 0, 0, 0, 0, 0},
    {-3, 6, -2, -2, 1, 0, 0, 0, 0, 0},
    {-6, 14, -8, -3, 4, -1, 0, 0, 0, 0},
    {-14, 36, -27, -2, 12, -6, 1, 0, 0, 0},
    {-36, 99, -88, 11, 32, -25, 8, -1, 0, 0},
    {-99, 286, -286, 78, 78, -90, 42, -10, 1, 0},
    {-286, 858, -936, 364, 168, -300, 184, -63, 12, -1},
};

static const double k_paper_kappa[] = {1.69, 1.75, 1.82, 1.84, 1.86, 1.88, 1.89};
static const double k_paper_milli[] = {150.0, 57.0, 18.0, 5.5, 1.6, 0.44, 0.12};

static int read_file(const char *path, unsigned char **buf, size_t *n) {
  FILE *fp;
  long sz;
  unsigned char *data;
  fp = fopen(path, "rb");
  if (!fp)
    return -1;
  if (fseek(fp, 0, SEEK_END) != 0) {
    fclose(fp);
    return -1;
  }
  sz = ftell(fp);
  if (sz < 0) {
    fclose(fp);
    return -1;
  }
  if (fseek(fp, 0, SEEK_SET) != 0) {
    fclose(fp);
    return -1;
  }
  data = (unsigned char *)malloc((size_t)sz + 1u);
  if (!data) {
    fclose(fp);
    return -1;
  }
  if (sz > 0 && fread(data, 1, (size_t)sz, fp) != (size_t)sz) {
    free(data);
    fclose(fp);
    return -1;
  }
  data[sz] = 0;
  fclose(fp);
  *buf = data;
  *n = (size_t)sz;
  return 0;
}

static CPMDCSession *open_session(const char *path) {
  unsigned char *buf = NULL;
  size_t n = 0;
  CPMDCSession *session;
  if (read_file(path, &buf, &n) != 0)
    return NULL;
  session = cpmdc_session_create(buf, n);
  free(buf);
  return session;
}

static void test_coefficients(void **state) {
  int order;
  (void)state;
  for (order = 3; order <= 9; ++order) {
    double kappa = 0.0;
    double alpha = 0.0;
    double c[10];
    int row = order - 3;
    int k;
    int sum = 0;
    assert_int_equal(cpmdc_xlbomd_coefficients(order, &kappa, &alpha, c, 10), 0);
    assert_true(kappa == k_paper_kappa[row]);
    assert_true(alpha == k_paper_milli[row] * 1.0e-3);
    for (k = 0; k <= order; ++k) {
      assert_true(c[k] == (double)k_paper_c[row][k]);
      sum += k_paper_c[row][k];
    }
    assert_int_equal(sum, 0);
  }
  {
    double kappa = 0.0;
    double alpha = 1.0;
    assert_int_equal(cpmdc_xlbomd_coefficients(0, &kappa, &alpha, NULL, 0), 0);
    assert_true(kappa == 2.0);
    assert_true(alpha == 0.0);
  }
  {
    double kappa = 0.0;
    assert_int_equal(cpmdc_xlbomd_coefficients(1, &kappa, NULL, NULL, 0), -1);
    assert_int_equal(cpmdc_xlbomd_coefficients(2, &kappa, NULL, NULL, 0), -1);
    assert_int_equal(cpmdc_xlbomd_coefficients(10, &kappa, NULL, NULL, 0), -1);
  }
  assert_int_equal(cpmdc_xlbomd_coefficients(5, NULL, NULL, NULL, 0), -1);
}

static uint32_t g_lcg = 1u;

static double lcg_unit(void) {
  g_lcg = g_lcg * 1664525u + 1013904223u;
  return (double)(g_lcg >> 8) * (1.0 / 16777216.0);
}

static void test_time_reversal(void **state) {
  enum { N = 8 };
  double p[REV_STEPS + 1][N];
  double back[REV_STEPS + 1][N];
  double density[REV_STEPS][N];
  double kappa = 0.0;
  int n;
  int i;
  int s;
  (void)state;
  g_lcg = 1u;
  assert_int_equal(cpmdc_xlbomd_coefficients(0, &kappa, NULL, NULL, 0), 0);
  for (i = 0; i < N; ++i) {
    p[0][i] = lcg_unit() - 0.5;
    p[1][i] = lcg_unit() - 0.5;
  }
  for (s = 1; s < REV_STEPS; ++s) {
    double hist[2 * N];
    for (i = 0; i < N; ++i) {
      density[s][i] = lcg_unit() - 0.5;
      hist[i] = p[s][i];
      hist[N + i] = p[s - 1][i];
    }
    assert_int_equal(
        cpmdc_xlbomd_propagate(0, density[s], N, hist, 2, p[s + 1]), 0);
  }
  memcpy(back[REV_STEPS], p[REV_STEPS], sizeof p[REV_STEPS]);
  memcpy(back[REV_STEPS - 1], p[REV_STEPS - 1], sizeof p[0]);
  for (s = REV_STEPS - 1; s >= 1; --s) {
    for (i = 0; i < N; ++i)
      back[s - 1][i] =
          (2.0 - kappa) * back[s][i] + kappa * density[s][i] - back[s + 1][i];
  }
  for (n = 0; n <= REV_STEPS; ++n) {
    for (i = 0; i < N; ++i)
      assert_true(fabs(back[n][i] - p[n][i]) < 1.0e-9);
  }
}

static void test_classify(void **state) {
  int order = -1;
  (void)state;
  assert_int_equal(cpmdc_xlbomd_classify("ASPC", 4, 4, &order),
                   CPMDC_WFN_EXTRAP_DECK);
  assert_int_equal(cpmdc_xlbomd_classify("POLY", 4, 3, &order),
                   CPMDC_WFN_EXTRAP_DECK);
  assert_int_equal(cpmdc_xlbomd_classify("", 0, 4, &order), CPMDC_WFN_EXTRAP_DECK);
  assert_int_equal(cpmdc_xlbomd_classify(NULL, 0, 0, &order),
                   CPMDC_WFN_EXTRAP_DECK);
  order = -1;
  assert_int_equal(cpmdc_xlbomd_classify("XLBOMD", 6, 0, &order),
                   CPMDC_WFN_EXTRAP_XLBOMD);
  assert_int_equal(order, 0);
  assert_int_equal(cpmdc_xlbomd_classify("XLBOMD", 6, 5, &order),
                   CPMDC_WFN_EXTRAP_XLBOMD);
  assert_int_equal(order, 5);
  assert_int_equal(cpmdc_xlbomd_classify("XLBOMD", 6, 9, &order),
                   CPMDC_WFN_EXTRAP_XLBOMD);
  assert_int_equal(cpmdc_xlbomd_classify("XLBOMD", 6, 2, &order), -1);
  assert_int_equal(cpmdc_xlbomd_classify("XLBOMD", 6, 10, &order), -1);
  assert_int_equal(cpmdc_xlbomd_classify("xlbomd", 6, 5, &order),
                   CPMDC_WFN_EXTRAP_DECK);
}

static void test_per_key(void **state) {
  double a[2] = {1.0, 0.0};
  double a2[2] = {1.1, 0.2};
  double b[2] = {0.0, 1.0};
  double b2[2] = {0.5, 0.25};
  double next[2];
  double bnext[2];
  double series[5][2];
  double got[2];
  double hist[4 * 2];
  double expect[2];
  CPMDCSession *session;
  CPMDCSession *other;
  int s;
  int i;
  (void)state;
  session = open_session(g_aspc_path);
  assert_non_null(session);
  other = open_session(g_aspc_path);
  assert_non_null(other);
  cpmdc_xlbomd_clear_key(session, 1, 1);
  cpmdc_xlbomd_clear_key(session, 2, 1);
  cpmdc_xlbomd_clear_key(session, 0, 0);
  cpmdc_xlbomd_clear_key(session, 0, 1);
  assert_int_equal(cpmdc_xlbomd_observe(session, 0, 1, 1, a, 2, next), 0);
  assert_int_equal(cpmdc_xlbomd_observe(session, 0, 2, 1, b, 2, bnext), 0);
  assert_int_equal(cpmdc_xlbomd_observe(session, 0, 1, 1, a2, 2, next), 1);
  assert_true(fabs(next[0] - 1.2) < 1.0e-12);
  assert_true(fabs(next[1] - 0.4) < 1.0e-12);
  assert_int_equal(cpmdc_xlbomd_observe(session, 0, 2, 1, b2, 2, bnext), 1);
  assert_true(fabs(bnext[0] - 1.0) < 1.0e-12);
  assert_true(fabs(bnext[1] + 0.5) < 1.0e-12);
  cpmdc_xlbomd_clear_key(session, 1, 1);
  assert_int_equal(cpmdc_xlbomd_observe(session, 0, 1, 1, a, 2, next), 0);
  assert_int_equal(cpmdc_xlbomd_observe(session, 0, 2, 1, b2, 2, bnext), 1);
  assert_int_equal(cpmdc_xlbomd_observe(other, 0, 2, 1, b2, 2, bnext), 0);

  for (s = 0; s < 5; ++s) {
    series[s][0] = 0.2 * (double)s;
    series[s][1] = 0.05 * (double)(s + 1);
  }
  cpmdc_xlbomd_clear_key(session, 11, 1);
  for (s = 0; s < 4; ++s)
    assert_int_equal(cpmdc_xlbomd_observe(session, 3, 11, 1, series[s], 2, got),
                     0);
  assert_int_equal(cpmdc_xlbomd_observe(session, 3, 11, 1, series[4], 2, got),
                   1);
  for (s = 0; s < 4; ++s) {
    hist[s * 2] = series[3 - s][0];
    hist[s * 2 + 1] = series[3 - s][1];
  }
  assert_int_equal(cpmdc_xlbomd_propagate(3, series[4], 2, hist, 4, expect), 0);
  for (i = 0; i < 2; ++i)
    assert_true(fabs(got[i] - expect[i]) < 1.0e-12);

  {
    double u = 3.0;
    double u2 = 3.5;
    double un = 0.0;
    double z = 9.0;
    double zn = 0.0;
    assert_int_equal(cpmdc_xlbomd_observe(session, 0, 0, 0, &u, 1, &un), 0);
    assert_int_equal(cpmdc_xlbomd_observe(session, 0, 0, 1, &z, 1, &zn), 0);
    assert_int_equal(cpmdc_xlbomd_observe(session, 0, 0, 0, &u2, 1, &un), 1);
    assert_true(fabs(un - 4.0) < 1.0e-12);
    assert_int_equal(cpmdc_xlbomd_observe(session, 0, 0, 1, &z, 1, &zn), 1);
    assert_true(fabs(zn - 9.0) < 1.0e-12);
  }
  cpmdc_session_destroy(other);
  cpmdc_session_destroy(session);
}

static char *render_message(const unsigned char *bytes, size_t n) {
  struct capn arena;
  CPMDParams_ptr root;
  char *deck = NULL;
  if (cpmdc_params_root(bytes, n, &arena, &root) != 0)
    return NULL;
  if (cpmdc_params_render_input_deck_alloc(root, NULL, &deck) != 0) {
    cpmdc_params_release(&arena);
    return NULL;
  }
  cpmdc_params_release(&arena);
  return deck;
}

static void test_parameter_surface(void **state) {
  unsigned char *aspc = NULL;
  unsigned char *xl = NULL;
  unsigned char *bad = NULL;
  size_t aspc_n = 0;
  size_t xl_n = 0;
  size_t bad_n = 0;
  char *aspc_deck = NULL;
  char *xl_deck = NULL;
  (void)state;
  assert_int_equal(read_file(g_aspc_path, &aspc, &aspc_n), 0);
  assert_int_equal(read_file(g_xl_path, &xl, &xl_n), 0);
  assert_int_equal(read_file(g_bad_path, &bad, &bad_n), 0);
  aspc_deck = render_message(aspc, aspc_n);
  xl_deck = render_message(xl, xl_n);
  assert_non_null(aspc_deck);
  assert_non_null(xl_deck);
  assert_non_null(strstr(aspc_deck, "EXTRAPOLATE WFN ASPC\n"));
  assert_non_null(strstr(aspc_deck, "\n  4\n"));
  assert_null(strstr(xl_deck, "EXTRAPOLATE WFN"));
  assert_non_null(strstr(xl_deck, "OPTIMIZE WAVEFUNCTION\n"));
  if (cpmdc_available()) {
    assert_int_equal(cpmdc_set_params(xl, xl_n), -1);
    assert_non_null(strstr(cpmdc_last_error(), "embed API"));
    assert_int_equal(cpmdc_xlbomd_configured_order(), -1);
    assert_int_equal(cpmdc_set_params(bad, bad_n), -1);
    assert_non_null(strstr(cpmdc_last_error(), "3 through 9"));
    assert_int_equal(cpmdc_xlbomd_configured_order(), -1);
    assert_null(cpmdc_session_create(xl, xl_n));
    assert_non_null(strstr(cpmdc_last_error(), "embed API"));
  } else {
    assert_int_equal(cpmdc_set_params(xl, xl_n), 0);
    assert_int_equal(cpmdc_xlbomd_configured_order(), 5);
    assert_int_equal(cpmdc_set_params(bad, bad_n), -1);
    assert_non_null(strstr(cpmdc_last_error(), "3 through 9"));
    assert_int_equal(cpmdc_xlbomd_configured_order(), 5);
  }
  assert_int_equal(cpmdc_set_params(aspc, aspc_n), 0);
  assert_int_equal(cpmdc_xlbomd_configured_order(), -1);
  free(aspc_deck);
  free(xl_deck);
  free(aspc);
  free(xl);
  free(bad);
}

static void test_copy_keeps_keys_apart(void **state) {
  unsigned char *aspc = NULL;
  size_t n = 0;
  CPMDCSession *session;
  double ang[9];
  double grad[9];
  int z[3] = {8, 1, 1};
  double src[2] = {1.0, 2.0};
  double src2[2] = {1.4, 2.2};
  double next[2];
  CPMDCResult result;
  int i;
  (void)state;
  assert_int_equal(read_file(g_aspc_path, &aspc, &n), 0);
  session = cpmdc_session_create(aspc, n);
  assert_non_null(session);
  for (i = 0; i < 9; ++i)
    ang[i] = k_bohr0[i] * CPMDC_BOHR_TO_ANGSTROM;
  result = cpmdc_session_energy_gradient(session, 3, ang, z, grad);
  assert_int_equal(result.ok, 1);
  cpmdc_xlbomd_clear_key(session, 7, 1);
  cpmdc_xlbomd_clear_key(session, 8, 1);
  assert_int_equal(cpmdc_xlbomd_observe(session, 0, 7, 1, src, 2, next), 0);
  assert_int_equal(cpmdc_xlbomd_observe(session, 0, 8, 1, src, 2, next), 0);
  assert_int_equal(cpmdc_session_copy_orbitals(session, 7, 7), 0);
  assert_int_equal(cpmdc_xlbomd_observe(session, 0, 7, 1, src2, 2, next), 1);
  assert_int_equal(cpmdc_session_copy_orbitals(session, 7, 8), 0);
  assert_int_equal(cpmdc_xlbomd_observe(session, 0, 8, 1, src2, 2, next), 0);
  assert_int_equal(cpmdc_xlbomd_observe(session, 0, 7, 1, src, 2, next), 1);
  cpmdc_session_destroy(session);
  free(aspc);
}

static double binom(int n, int k) {
  double c = 1.0;
  int kk;
  int i;
  if (k < 0 || k > n)
    return 0.0;
  if (k == 0 || k == n)
    return 1.0;
  kk = k > n - k ? n - k : k;
  c = 1.0;
  for (i = 1; i <= kk; ++i)
    c = c * (double)(n - kk + i) / (double)i;
  return c;
}

static void aspc_guess(const double *hist, int nhist, int n, double *out) {
  int k = nhist < ASPC_ORDER ? nhist : ASPC_ORDER;
  int j;
  int i;
  double denom;
  memset(out, 0, (size_t)n * sizeof(*out));
  if (k <= 0)
    return;
  denom = binom(2 * k - 2, k - 1);
  for (j = 1; j <= k; ++j) {
    double sign = (j % 2) ? 1.0 : -1.0;
    double b = sign * (double)j * binom(2 * k, k - j) / denom;
    for (i = 0; i < n; ++i)
      out[i] += b * hist[(j - 1) * n + i];
  }
}

static void invariants(const double *r, double *l1, double *l2, double *theta,
                       double *cosine, double *u1, double *u2) {
  double v1[3];
  double v2[3];
  double c;
  int k;
  for (k = 0; k < 3; ++k) {
    v1[k] = r[3 + k] - r[k];
    v2[k] = r[6 + k] - r[k];
  }
  *l1 = sqrt(v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]);
  *l2 = sqrt(v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2]);
  for (k = 0; k < 3; ++k) {
    u1[k] = v1[k] / *l1;
    u2[k] = v2[k] / *l2;
  }
  c = u1[0] * u2[0] + u1[1] * u2[1] + u1[2] * u2[2];
  if (c > 1.0)
    c = 1.0;
  if (c < -1.0)
    c = -1.0;
  *cosine = c;
  *theta = acos(c);
}

static double water_gradient(const double *r, double *grad) {
  double l1, l2, theta, cosine, u1[3], u2[3];
  double e;
  double dl1;
  double dl2;
  double sine;
  int k;
  invariants(r, &l1, &l2, &theta, &cosine, u1, u2);
  e = 0.5 * k_kr * (l1 - k_r0) * (l1 - k_r0);
  e += 0.5 * k_kr * (l2 - k_r0) * (l2 - k_r0);
  e += 0.5 * k_kth * (theta - k_theta0) * (theta - k_theta0);
  if (grad) {
    memset(grad, 0, 9 * sizeof(*grad));
    dl1 = k_kr * (l1 - k_r0);
    dl2 = k_kr * (l2 - k_r0);
    for (k = 0; k < 3; ++k) {
      grad[3 + k] += dl1 * u1[k];
      grad[k] -= dl1 * u1[k];
      grad[6 + k] += dl2 * u2[k];
      grad[k] -= dl2 * u2[k];
    }
    sine = sqrt(fmax(0.0, 1.0 - cosine * cosine));
    if (sine > 1.0e-12) {
      double dedc = k_kth * (theta - k_theta0) * (-1.0 / sine);
      double g1[3];
      double g2[3];
      for (k = 0; k < 3; ++k) {
        g1[k] = dedc * (u2[k] - cosine * u1[k]) / l1;
        g2[k] = dedc * (u1[k] - cosine * u2[k]) / l2;
        grad[3 + k] += g1[k];
        grad[6 + k] += g2[k];
        grad[k] -= g1[k] + g2[k];
      }
    }
  }
  return e;
}

static double mass_of(int atom) { return atom == 0 ? k_mass_o : k_mass_h; }

/* Same linear mix the reference session uses for an XL-BOMD step. */
static int mix_scf(const double *guess, const double *target, int n, double *out) {
  int steps = 0;
  int i;
  memcpy(out, guess, (size_t)n * sizeof(*out));
  for (;;) {
    double sum = 0.0;
    ++steps;
    for (i = 0; i < n; ++i) {
      double d = out[i] - target[i];
      sum += d * d;
    }
    if (sqrt(sum) < k_scf_tol || steps >= k_scf_max)
      break;
    for (i = 0; i < n; ++i)
      out[i] = (1.0 - k_beta) * out[i] + k_beta * target[i];
  }
  return steps;
}

static void aspc_push(double *hist, int *count, const double *d, int n) {
  int keep = *count < ASPC_ORDER ? *count : ASPC_ORDER - 1;
  if (keep > 0)
    memmove(hist + n, hist, (size_t)keep * (size_t)n * sizeof(*hist));
  memcpy(hist, d, (size_t)n * sizeof(*d));
  if (*count < ASPC_ORDER)
    *count += 1;
}

/* Geometry at the start of each nuclear step, before that step moves. */
static void water_frames(double geom[WATER_STEPS][9]) {
  double r[9];
  double v[9];
  double grad[9];
  int step;
  int i;
  memcpy(r, k_bohr0, sizeof r);
  memset(v, 0, sizeof v);
  water_gradient(r, grad);
  for (step = 0; step < WATER_STEPS; ++step) {
    memcpy(geom[step], r, sizeof r);
    for (i = 0; i < 9; ++i)
      v[i] += 0.5 * k_dt * (-grad[i] / mass_of(i / 3));
    for (i = 0; i < 9; ++i)
      r[i] += k_dt * v[i];
    water_gradient(r, grad);
    for (i = 0; i < 9; ++i)
      v[i] += 0.5 * k_dt * (-grad[i] / mass_of(i / 3));
  }
}

static void nuclear_drift(double *drift, double *max_abs) {
  double r[9];
  double v[9];
  double grad[9];
  int step;
  int i;
  double e0;
  memcpy(r, k_bohr0, sizeof r);
  memset(v, 0, sizeof v);
  e0 = water_gradient(r, grad);
  *drift = 0.0;
  *max_abs = 0.0;
  for (step = 0; step < WATER_STEPS; ++step) {
    double pot;
    double kin = 0.0;
    for (i = 0; i < 9; ++i)
      v[i] += 0.5 * k_dt * (-grad[i] / mass_of(i / 3));
    for (i = 0; i < 9; ++i)
      r[i] += k_dt * v[i];
    pot = water_gradient(r, grad);
    for (i = 0; i < 9; ++i)
      v[i] += 0.5 * k_dt * (-grad[i] / mass_of(i / 3));
    for (i = 0; i < 9; ++i)
      kin += 0.5 * mass_of(i / 3) * v[i] * v[i];
    {
      double delta = (pot + kin) - e0;
      if (fabs(delta) > *max_abs)
        *max_abs = fabs(delta);
      *drift = delta;
    }
  }
}

static void test_water_forces(void **state) {
  double grad[9];
  double r[9];
  double e0;
  double sum[3] = {0, 0, 0};
  int i;
  (void)state;
  memcpy(r, k_bohr0, sizeof r);
  e0 = water_gradient(r, grad);
  for (i = 0; i < 9; ++i)
    sum[i % 3] += grad[i];
  for (i = 0; i < 3; ++i)
    assert_true(fabs(sum[i]) < 1.0e-10);
  for (i = 0; i < 9; ++i) {
    double ep, em, fd;
    const double eps = 1.0e-6;
    r[i] += eps;
    ep = water_gradient(r, NULL);
    r[i] -= 2.0 * eps;
    em = water_gradient(r, NULL);
    r[i] += eps;
    fd = (ep - em) / (2.0 * eps);
    assert_true(fabs(fd - grad[i]) < 1.0e-6);
  }
  (void)e0;
  {
    double denom = binom(6, 3);
    double b1 = 1.0 * binom(8, 3) / denom;
    double b2 = -2.0 * binom(8, 2) / denom;
    double b3 = 3.0 * binom(8, 1) / denom;
    double b4 = -4.0 * binom(8, 0) / denom;
    double hist[4] = {1.0, 1.0, 1.0, 1.0};
    double out = 0.0;
    assert_true(fabs(b1 - 2.8) < 1.0e-12);
    assert_true(fabs(b2 + 2.8) < 1.0e-12);
    assert_true(fabs(b3 - 1.2) < 1.0e-12);
    assert_true(fabs(b4 + 0.2) < 1.0e-12);
    aspc_guess(hist, 4, 1, &out);
    assert_true(fabs(out - 1.0) < 1.0e-12);
  }
}

static int parse_double(const char *line, const char *key, double *out) {
  if (strncmp(line, key, strlen(key)) != 0)
    return -1;
  *out = strtod(line + strlen(key), NULL);
  return 0;
}

static void test_water_golden(void **state) {
  FILE *fp;
  char line[512];
  double drift = 0.0;
  double max_abs = 0.0;
  double want_drift = 0.0;
  double want_max = 0.0;
  int saw_drift = 0;
  int saw_max = 0;
  (void)state;
  nuclear_drift(&drift, &max_abs);
  fp = fopen(g_golden_path, "r");
  assert_non_null(fp);
  while (fgets(line, sizeof line, fp)) {
    if (line[0] == '#' || line[0] == '\n')
      continue;
    if (parse_double(line, "drift ", &want_drift) == 0)
      saw_drift = 1;
    else if (parse_double(line, "max_abs ", &want_max) == 0)
      saw_max = 1;
  }
  fclose(fp);
  assert_int_equal(saw_drift + saw_max, 2);
  assert_true(fabs(drift - want_drift) < 1.0e-9);
  assert_true(fabs(max_abs - want_max) < 1.0e-9);
}

/* Twelve nuclear steps. No extrapolation reuses the previous converged
 * vector. ASPC is order 4. XL-BOMD is the session propagation of order 5,
 * checked here with the same mixer and, on the reference build, with
 * cpmdc_session_calculate_au(). */
static void test_water_scf(void **state) {
  double geom[WATER_STEPS][9];
  double guess[9];
  double converged[9];
  double next[9];
  double hist[ASPC_ORDER * 9];
  int none_total = 0;
  int aspc_total = 0;
  int xl_total = 0;
  int nhist = 0;
  int step;
  int want_none = -1;
  int want_aspc = -1;
  int want_xl = -1;
  int saw = 0;
  FILE *fp;
  char line[512];
  CPMDCSession *session;
  (void)state;
  water_frames(geom);
  memset(guess, 0, sizeof guess);
  for (step = 0; step < WATER_STEPS; ++step) {
    none_total += mix_scf(guess, geom[step], 9, converged);
    memcpy(guess, converged, sizeof guess);
  }
  memset(guess, 0, sizeof guess);
  memset(hist, 0, sizeof hist);
  for (step = 0; step < WATER_STEPS; ++step) {
    aspc_guess(hist, nhist, 9, guess);
    aspc_total += mix_scf(guess, geom[step], 9, converged);
    aspc_push(hist, &nhist, converged, 9);
  }
  session = open_session(g_aspc_path);
  assert_non_null(session);
  memset(guess, 0, sizeof guess);
  for (step = 0; step < WATER_STEPS; ++step) {
    int rc;
    xl_total += mix_scf(guess, geom[step], 9, converged);
    rc = cpmdc_xlbomd_observe(session, XL_ORDER, 42, 1, converged, 9, next);
    assert_true(rc >= 0);
    memcpy(guess, next, sizeof guess);
  }
  cpmdc_session_destroy(session);
  if (!cpmdc_available()) {
    CPMDCSession *xl;
    int session_total = 0;
    int z[3] = {8, 1, 1};
    double energy = 0.0;
    double forces[9];
    int order = -1;
    int steps = -1;
    int propagated = -1;
    xl = open_session(g_xl_path);
    assert_non_null(xl);
    for (step = 0; step < WATER_STEPS; ++step) {
      assert_int_equal(cpmdc_session_calculate_au(xl, 3, geom[step], z, NULL,
                                                   &energy, forces, NULL),
                       0);
      assert_int_equal(
          cpmdc_session_xlbomd_last(xl, &order, &steps, &propagated), 0);
      session_total += steps;
    }
    cpmdc_session_destroy(xl);
    assert_int_equal(session_total, xl_total);
  }
  printf("scf_none %d\n", none_total);
  printf("scf_aspc %d\n", aspc_total);
  printf("scf_xlbomd %d\n", xl_total);
  fp = fopen(g_golden_path, "r");
  assert_non_null(fp);
  while (fgets(line, sizeof line, fp)) {
    if (sscanf(line, "scf_none %d", &want_none) == 1 ||
        sscanf(line, "scf_aspc %d", &want_aspc) == 1 ||
        sscanf(line, "scf_xlbomd %d", &want_xl) == 1)
      saw += 1;
  }
  fclose(fp);
  assert_int_equal(saw, 3);
  assert_int_equal(none_total, want_none);
  assert_int_equal(aspc_total, want_aspc);
  assert_int_equal(xl_total, want_xl);
  assert_true(xl_total > want_aspc);
}

static void geometry(int step, double *pos) {
  memcpy(pos, k_bohr0, 9 * sizeof(*pos));
  pos[0] += 0.02 * (double)step;
}

static int header_image(unsigned char **out, size_t *n) {
  char header[80];
  int32_t mark = 80;
  unsigned char *buf = (unsigned char *)malloc(88);
  if (!buf)
    return -1;
  memset(header, ' ', sizeof header);
  memcpy(header, " STUTTGART VERSION 3.0", 22);
  memcpy(buf, &mark, 4);
  memcpy(buf + 4, header, 80);
  memcpy(buf + 84, &mark, 4);
  *out = buf;
  *n = 88;
  return 0;
}

static void test_session_calculate(void **state) {
  CPMDCSession *session;
  CPMDCSession *other;
  CPMDCSession *restored;
  unsigned char *xl = NULL;
  size_t xl_n = 0;
  double pos[9];
  double forces[9];
  double energy = 0.0;
  double cold_energy = 0.0;
  int z[3] = {8, 1, 1};
  int step;
  int order = -2;
  int steps = -1;
  int propagated = -2;
  int held_steps = -1;
  int held_prop = -2;
  void *blob = NULL;
  size_t blob_n = 0;
  void *again = NULL;
  size_t again_n = 0;
  unsigned char *image = NULL;
  size_t image_n = 0;
  void *written = NULL;
  size_t written_n = 0;
  char err[128];
  CPMDCRestart *file;
  CPMDCRestart *back;
  void *got = NULL;
  size_t got_n = 0;
  double ang[9];
  double grad[9];
  CPMDCResult bare;
  int i;
  (void)state;
  assert_int_equal(read_file(g_xl_path, &xl, &xl_n), 0);
  if (cpmdc_available()) {
    assert_null(cpmdc_session_create(xl, xl_n));
    assert_non_null(strstr(cpmdc_last_error(), "embed API"));
    free(xl);
    return;
  }
  for (i = 0; i < 9; ++i)
    ang[i] = k_bohr0[i] * CPMDC_BOHR_TO_ANGSTROM;
  bare = cpmdc_energy_gradient(3, ang, z, xl, xl_n, grad);
  assert_int_equal(bare.ok, 0);
  assert_non_null(strstr(cpmdc_last_error(), "requires a session"));
  session = cpmdc_session_create(xl, xl_n);
  assert_non_null(session);
  for (step = 0; step < 7; ++step) {
    geometry(step, pos);
    assert_int_equal(
        cpmdc_session_calculate_au(session, 3, pos, z, NULL, &energy, forces,
                                   NULL),
        0);
    assert_int_equal(cpmdc_session_xlbomd_last(session, &order, &steps,
                                               &propagated),
                     0);
    assert_int_equal(order, XL_ORDER);
    assert_true(steps > 0 && steps < 400);
    assert_int_equal(propagated, step < 6 ? 0 : 1);
  }
  held_steps = steps;
  held_prop = propagated;
  assert_int_equal(cpmdc_session_xlbomd_export(session, &blob, &blob_n), 0);
  assert_non_null(blob);
  assert_true(blob_n > 0);
  geometry(6, pos);
  assert_int_equal(
      cpmdc_session_calculate_au(session, 3, pos, z, NULL, &energy, forces,
                                 NULL),
      0);
  assert_int_equal(
      cpmdc_session_xlbomd_last(session, &order, &steps, &propagated), 0);
  assert_int_equal(steps, held_steps);
  assert_int_equal(propagated, held_prop);
  assert_int_equal(cpmdc_session_xlbomd_export(session, &again, &again_n), 0);
  assert_int_equal(again_n, blob_n);
  assert_int_equal(memcmp(again, blob, blob_n), 0);
  free(again);
  other = cpmdc_session_create(xl, xl_n);
  assert_non_null(other);
  geometry(7, pos);
  assert_int_equal(
      cpmdc_session_calculate_au(other, 3, pos, z, NULL, &cold_energy, forces,
                                 NULL),
      0);
  assert_int_equal(
      cpmdc_session_xlbomd_last(other, &order, &steps, &propagated), 0);
  assert_int_equal(propagated, 0);
  assert_int_equal(header_image(&image, &image_n), 0);
  file = cpmdc_restart_read_mem(image, image_n, err, sizeof err);
  assert_non_null(file);
  assert_int_equal(cpmdc_restart_ncoords(file), 0);
  assert_int_equal(cpmdc_restart_set_xlbomd(file, blob, blob_n), 0);
  assert_int_equal(
      cpmdc_restart_write_mem(file, &written, &written_n, err, sizeof err), 0);
  cpmdc_restart_free(file);
  back = cpmdc_restart_read_mem(written, written_n, err, sizeof err);
  assert_non_null(back);
  assert_int_equal(cpmdc_restart_ncoords(back), 0);
  assert_int_equal(cpmdc_restart_xlbomd(back, &got, &got_n), 0);
  assert_int_equal(got_n, blob_n);
  assert_int_equal(memcmp(got, blob, blob_n), 0);
  restored = cpmdc_session_create(xl, xl_n);
  assert_non_null(restored);
  assert_int_equal(cpmdc_session_xlbomd_import(restored, got, got_n), 0);
  geometry(7, pos);
  assert_int_equal(
      cpmdc_session_calculate_au(session, 3, pos, z, NULL, &energy, forces,
                                 NULL),
      0);
  assert_int_equal(
      cpmdc_session_xlbomd_last(session, &order, &steps, &propagated), 0);
  held_steps = steps;
  held_prop = propagated;
  assert_int_equal(propagated, 1);
  assert_true(energy == cold_energy);
  assert_int_equal(
      cpmdc_session_calculate_au(restored, 3, pos, z, NULL, &cold_energy,
                                 forces, NULL),
      0);
  assert_int_equal(
      cpmdc_session_xlbomd_last(restored, &order, &steps, &propagated), 0);
  assert_int_equal(steps, held_steps);
  assert_int_equal(propagated, held_prop);
  assert_true(cold_energy == energy);
  cpmdc_restart_free(back);
  cpmdc_session_destroy(restored);
  cpmdc_session_destroy(other);
  cpmdc_session_destroy(session);
  free(got);
  free(written);
  free(image);
  free(blob);
  free(xl);
}

int main(int argc, char **argv) {
  const struct CMUnitTest tests[] = {
      cmocka_unit_test(test_coefficients),
      cmocka_unit_test(test_time_reversal),
      cmocka_unit_test(test_classify),
      cmocka_unit_test(test_per_key),
      cmocka_unit_test(test_parameter_surface),
      cmocka_unit_test(test_copy_keeps_keys_apart),
      cmocka_unit_test(test_water_forces),
      cmocka_unit_test(test_water_golden),
      cmocka_unit_test(test_water_scf),
      cmocka_unit_test(test_session_calculate),
  };
  if (argc != 5) {
    fprintf(stderr,
            "usage: test_xlbomd ASPC.bin XL.bin BAD.bin GOLDEN\n");
    return 2;
  }
  g_aspc_path = argv[1];
  g_xl_path = argv[2];
  g_bad_path = argv[3];
  g_golden_path = argv[4];
  return cmocka_run_group_tests(tests, NULL, NULL);
}
