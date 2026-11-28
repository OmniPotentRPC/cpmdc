/* SPDX-License-Identifier: MIT */
#include "cpmdc_damp.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

static int failures = 0;

static void check(int cond, const char *msg) {
  if (!cond) {
    fprintf(stderr, "fail: %s\n", msg);
    failures++;
  }
}

static int near(double got, double ref, double tol) {
  return isfinite(got) && fabs(got - ref) <= tol;
}

static void test_coefficients(void) {
  double fccc = 7.0, verl1 = 7.0, verl2 = 7.0, verl3 = 7.0;

  check(cpmdc_damp_coefficients(0.0, NULL, &verl1, &verl2, &verl3) == -1,
        "null fccc writes nothing");
  check(fccc == 7.0 && verl1 == 7.0 && verl2 == 7.0 && verl3 == 7.0,
        "null fccc leaves the other slots");
  check(cpmdc_damp_coefficients(-1.0, &fccc, &verl1, &verl2, &verl3) == -1,
        "friction of -1 is rejected");
  check(cpmdc_damp_coefficients(NAN, &fccc, &verl1, &verl2, &verl3) == -1,
        "non-finite friction is rejected");
  fccc = verl1 = verl2 = verl3 = 7.0;
  check(cpmdc_damp_coefficients(-2.0, &fccc, &verl1, &verl2, &verl3) == -1,
        "friction -2 is rejected");
  check(fccc == 7.0 && verl1 == 7.0 && verl2 == 7.0 && verl3 == 7.0,
        "friction -2 leaves weights");
  fccc = verl1 = verl2 = verl3 = 7.0;
  check(cpmdc_damp_coefficients(nextafter(-1.0, -INFINITY), &fccc, &verl1, &verl2,
                                &verl3) == -1,
        "friction just below -1 is rejected");
  check(fccc == 7.0 && verl1 == 7.0 && verl2 == 7.0 && verl3 == 7.0,
        "friction just below -1 leaves weights");
  check(cpmdc_damp_coefficients(0.0, &fccc, &verl1, &verl2, &verl3) == 0,
        "verlet weights");
  check(fccc == 1.0 && verl1 == 2.0 && verl2 == -1.0 && verl3 == 1.0,
        "friction 0 is undamped Verlet");
  check(cpmdc_damp_coefficients(1.0, &fccc, &verl1, &verl2, &verl3) == 0,
        "friction 1");
  check(fccc == 0.5 && verl1 == 1.0 && verl2 == 0.0 && verl3 == 0.5,
        "friction 1 weights");
}

static void test_one_step(void) {
  double c0_re[2] = {1.0, -0.25};
  double c0_im[2] = {0.5, 0.125};
  double cm_re[2] = {0.5, 0.0};
  double cm_im[2] = {-0.25, 0.25};
  double force_re[2] = {-1.0, 0.5};
  double force_im[2] = {0.25, -0.5};
  double ema[2] = {2.0, 1.0};
  double saved_re[2];
  double saved_im[2];
  double dt2 = 0.25;

  memcpy(saved_re, cm_re, sizeof saved_re);
  memcpy(saved_im, cm_im, sizeof saved_im);
  check(cpmdc_damp_step(0, c0_re, c0_im, cm_re, cm_im, force_re, force_im, NULL,
                        0.0, dt2) == -1,
        "n < 1");
  check(memcmp(cm_re, saved_re, sizeof saved_re) == 0 &&
            memcmp(cm_im, saved_im, sizeof saved_im) == 0,
        "rejected step leaves cm");
  check(cpmdc_damp_step(2, c0_re, c0_im, cm_re, cm_im, force_re, force_im, ema,
                        -1.0, dt2) == -1,
        "bad friction leaves cm");
  check(cpmdc_damp_step(2, c0_re, c0_im, cm_re, cm_im, force_re, force_im, ema,
                        -2.0, dt2) == -1,
        "friction -2 step");
  check(memcmp(cm_re, saved_re, sizeof saved_re) == 0 &&
            memcmp(cm_im, saved_im, sizeof saved_im) == 0,
        "friction -2 step leaves cm");
  check(cpmdc_damp_step(2, c0_re, c0_im, cm_re, cm_im, force_re, force_im, ema,
                        nextafter(-1.0, -INFINITY), dt2) == -1,
        "friction just below -1 step");
  check(memcmp(cm_re, saved_re, sizeof saved_re) == 0 &&
            memcmp(cm_im, saved_im, sizeof saved_im) == 0,
        "friction just below -1 step leaves cm");
  cm_re[0] = INFINITY;
  check(cpmdc_damp_step(2, c0_re, c0_im, cm_re, cm_im, force_re, force_im, NULL,
                        0.0, dt2) == -1,
        "non-finite cm is rejected");
  cm_re[0] = 0.5;

  check(cpmdc_damp_step(2, c0_re, c0_im, cm_re, cm_im, force_re, force_im, ema,
                        0.0, dt2) == 0,
        "verlet step");
  /* cm = 2*c0 - cm + dt2 * ema * force */
  check(near(cm_re[0], 2.0 * 1.0 - 0.5 + dt2 * 2.0 * -1.0, 0.0), "re0");
  check(near(cm_im[0], 2.0 * 0.5 - (-0.25) + dt2 * 2.0 * 0.25, 0.0), "im0");
  check(near(cm_re[1], 2.0 * -0.25 - 0.0 + dt2 * 1.0 * 0.5, 0.0), "re1");
  check(near(cm_im[1], 2.0 * 0.125 - 0.25 + dt2 * 1.0 * -0.5, 0.0), "im1");
  check(c0_re[0] == 1.0 && c0_im[1] == 0.125, "c0 is not written");

  c0_re[0] = 1.0;
  c0_im[0] = 0.5;
  cm_re[0] = 0.5;
  cm_im[0] = -0.25;
  check(cpmdc_damp_step(1, c0_re, c0_im, cm_re, cm_im, force_re, force_im, ema,
                        1.0, dt2) == 0,
        "damped step");
  /* verl1 = 1, verl2 = 0, scale = 0.5 * dt2 * ema */
  check(near(cm_re[0], 1.0 + 0.5 * dt2 * 2.0 * -1.0, 0.0), "damped re");
  check(near(cm_im[0], 0.5 + 0.5 * dt2 * 2.0 * 0.25, 0.0), "damped im");
}

static void test_null_ema_matches_ones(void) {
  double c0_re = 0.3, c0_im = -0.1, f_re = 0.2, f_im = -0.4;
  double cm_a_re = 0.1, cm_a_im = 0.2, cm_b_re = 0.1, cm_b_im = 0.2;
  double ones = 1.0;

  check(cpmdc_damp_step(1, &c0_re, &c0_im, &cm_a_re, &cm_a_im, &f_re, &f_im, NULL,
                        0.25, 0.04) == 0,
        "null ema");
  check(cpmdc_damp_step(1, &c0_re, &c0_im, &cm_b_re, &cm_b_im, &f_re, &f_im, &ones,
                        0.25, 0.04) == 0,
        "unit ema");
  check(cm_a_re == cm_b_re && cm_a_im == cm_b_im, "null ema is ones");
}

static void test_swap_and_advance(void) {
  double c0_re = 1.0, c0_im = 2.0, cm_re = 3.0, cm_im = 4.0;
  double f_re = -0.5, f_im = 0.25;
  double expect_re, expect_im;

  check(cpmdc_damp_swap(1, NULL, &c0_im, &cm_re, &cm_im) == -1, "null swap");
  check(c0_re == 1.0 && cm_re == 3.0, "rejected swap writes nothing");
  check(cpmdc_damp_swap(1, &c0_re, &c0_im, &cm_re, &cm_im) == 0, "swap");
  check(c0_re == 3.0 && c0_im == 4.0 && cm_re == 1.0 && cm_im == 2.0, "swapped");
  c0_re = NAN;
  check(cpmdc_damp_swap(1, &c0_re, &c0_im, &cm_re, &cm_im) == 0, "swap nan");
  check(isnan(cm_re) && c0_re == 1.0 && c0_im == 2.0 && cm_im == 4.0,
        "nan swapped");

  c0_re = 0.8;
  c0_im = -0.2;
  cm_re = 0.8;
  cm_im = -0.2;
  check(cpmdc_damp_step(1, &c0_re, &c0_im, &cm_re, &cm_im, &f_re, &f_im, NULL, 0.1,
                        0.01) == 0,
        "step before compare");
  expect_re = cm_re;
  expect_im = cm_im;
  c0_re = 0.8;
  c0_im = -0.2;
  cm_re = 0.8;
  cm_im = -0.2;
  check(cpmdc_damp_advance(1, &c0_re, &c0_im, &cm_re, &cm_im, &f_re, &f_im, NULL,
                           0.1, 0.01) == 0,
        "advance");
  check(c0_re == expect_re && c0_im == expect_im, "advance stores c(t+dt) in c0");
  check(cm_re == 0.8 && cm_im == -0.2, "advance stores c(t) in cm");

  c0_re = 0.8;
  cm_re = 0.2;
  check(cpmdc_damp_advance(1, &c0_re, &c0_im, &cm_re, &cm_im, &f_re, &f_im, NULL,
                           -1.0, 0.01) == -1,
        "failed advance");
  check(c0_re == 0.8 && cm_re == 0.2, "failed advance writes nothing");

  {
    double below = nextafter(-1.0, -INFINITY);
    double a_re = 0.8, a_im = -0.2, b_re = 0.2, b_im = 0.4;

    check(cpmdc_damp_advance(1, &a_re, &a_im, &b_re, &b_im, &f_re, &f_im, NULL, -2.0,
                             0.01) == -1,
          "friction -2 advance");
    check(a_re == 0.8 && a_im == -0.2 && b_re == 0.2 && b_im == 0.4,
          "friction -2 advance writes nothing");
    a_re = 0.8;
    a_im = -0.2;
    b_re = 0.2;
    b_im = 0.4;
    check(cpmdc_damp_advance(1, &a_re, &a_im, &b_re, &b_im, &f_re, &f_im, NULL, below,
                             0.01) == -1,
          "friction just below -1 advance");
    check(a_re == 0.8 && a_im == -0.2 && b_re == 0.2 && b_im == 0.4,
          "friction just below -1 advance writes nothing");
  }
}

static void test_grease(void) {
  double friction = 0.1;

  check(cpmdc_damp_grease(NULL, 1.0) == -1, "null grease");
  check(cpmdc_damp_grease(&friction, NAN) == -1, "non-finite grease");
  check(friction == 0.1, "rejected grease leaves friction");
  check(cpmdc_damp_grease(&friction, 1.0) == 0, "grease 1");
  check(friction == 0.1, "grease 1 is unchanged");
  check(cpmdc_damp_grease(&friction, 0.5) == 0, "grease one half");
  check(friction == 0.05, "friction halved");

  friction = -2.0;
  check(cpmdc_damp_grease(&friction, 1.0) == -1, "friction -2 grease 1");
  check(friction == -2.0, "friction -2 grease 1 unchanged");
  friction = -2.0;
  check(cpmdc_damp_grease(&friction, -0.1) == -1, "friction -2 grease -0.1");
  check(friction == -2.0, "friction -2 grease -0.1 unchanged");
  friction = -2.0;
  check(cpmdc_damp_grease(&friction, 0.0) == -1, "friction -2 grease 0");
  check(friction == -2.0, "friction -2 grease 0 unchanged");
  friction = -0.5;
  check(cpmdc_damp_grease(&friction, 4.0) == -1, "friction -0.5 grease 4");
  check(friction == -0.5, "friction -0.5 grease 4 unchanged");
  friction = -1.0;
  check(cpmdc_damp_grease(&friction, 1.0) == -1, "friction -1 grease 1");
  check(friction == -1.0, "friction -1 grease 1 unchanged");
  friction = 0.5;
  check(cpmdc_damp_grease(&friction, -2.0) == -1, "product -1");
  check(friction == 0.5, "product -1 unchanged");
  friction = -0.25;
  check(cpmdc_damp_grease(&friction, 2.0) == 0, "valid grease");
  check(friction == -0.5, "valid grease writes the product");
}

static void quench(double friction, int nstep, double *final_abs) {
  double ceq = 0.25;
  double dt = 0.05;
  double dt2 = dt * dt;
  double c0_re = 1.0, c0_im = 0.0, cm_re = 1.0, cm_im = 0.0;
  int i;

  for (i = 0; i < nstep; ++i) {
    double force = -(c0_re - ceq);
    double zero = 0.0;
    if (cpmdc_damp_advance(1, &c0_re, &c0_im, &cm_re, &cm_im, &force, &zero, NULL,
                           friction, dt2) != 0) {
      *final_abs = HUGE_VAL;
      return;
    }
  }
  *final_abs = fabs(c0_re - ceq);
}

static void test_quench(void) {
  double damped = 0.0;
  double undamped = 0.0;

  quench(0.1, 2000, &damped);
  quench(0.0, 2000, &undamped);
  check(damped <= 1e-9, "damped trajectory reaches the minimum");
  check(undamped > 1e-3, "undamped Verlet does not quench");
  fprintf(stdout, "quench friction 0.1 steps 2000 abs_error %.6e\n", damped);
  fprintf(stdout, "quench friction 0 steps 2000 abs_error %.6e\n", undamped);
}

static void test_absent_from_libcpmdc(const char *path) {
  FILE *fp;
  char buf[4096];
  size_t n;
  const char *needle = "cpmdc_damp_step";
  size_t needle_n = strlen(needle);
  char hold[64];
  size_t hold_n = 0;
  int found = 0;

  check(path != NULL && path[0] != '\0', "library path");
  fp = fopen(path, "rb");
  check(fp != NULL, "open libcpmdc");
  if (fp == NULL)
    return;
  while ((n = fread(buf, 1, sizeof buf, fp)) > 0) {
    size_t i;
    for (i = 0; i < n; ++i) {
      hold[hold_n++] = buf[i];
      if (hold_n == needle_n) {
        if (memcmp(hold, needle, needle_n) == 0)
          found = 1;
        memmove(hold, hold + 1, needle_n - 1);
        hold_n = needle_n - 1;
      }
    }
  }
  fclose(fp);
  check(!found, "libcpmdc does not contain cpmdc_damp_step");
}

int main(int argc, char **argv) {
  test_coefficients();
  test_one_step();
  test_null_ema_matches_ones();
  test_swap_and_advance();
  test_grease();
  test_quench();
  test_absent_from_libcpmdc(argc > 1 ? argv[1] : "");
  if (failures != 0) {
    fprintf(stderr, "%d damped-electron checks failed\n", failures);
    return 1;
  }
  return 0;
}
