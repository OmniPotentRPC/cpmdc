/* SPDX-License-Identifier: MIT */
#define _POSIX_C_SOURCE 200809L
#include "cpmdc_damp.h"

#include <math.h>
#include <stdio.h>
#include <time.h>

static double quench(double friction, int nstep) {
  double ceq = 0.25;
  double dt = 0.05;
  double dt2 = dt * dt;
  double c0_re = 1.0, c0_im = 0.0, cm_re = 1.0, cm_im = 0.0;
  int i;

  for (i = 0; i < nstep; ++i) {
    double force = -(c0_re - ceq);
    double zero = 0.0;
    if (cpmdc_damp_advance(1, &c0_re, &c0_im, &cm_re, &cm_im, &force, &zero, NULL,
                           friction, dt2) != 0)
      return HUGE_VAL;
  }
  return fabs(c0_re - ceq);
}

static double seconds_of(int repeats, double friction, int nstep) {
  struct timespec start, stop;
  int i;

  clock_gettime(CLOCK_MONOTONIC, &start);
  for (i = 0; i < repeats; ++i)
    (void)quench(friction, nstep);
  clock_gettime(CLOCK_MONOTONIC, &stop);
  return (double)(stop.tv_sec - start.tv_sec) +
         (double)(stop.tv_nsec - start.tv_nsec) * 1e-9;
}

int main(void) {
  const int nstep = 2000;
  const int repeats = 200;
  double damp_err = quench(0.1, nstep);
  double verlet_err = quench(0.0, nstep);
  double damp_s = seconds_of(repeats, 0.1, nstep);
  double verlet_s = seconds_of(repeats, 0.0, nstep);

  printf("damp friction 0.1 steps %d abs_error %.6e repeats %d seconds %.6f\n",
         nstep, damp_err, repeats, damp_s);
  printf("verlet friction 0 steps %d abs_error %.6e repeats %d seconds %.6f\n",
         nstep, verlet_err, repeats, verlet_s);
  if (!(damp_err <= 1e-9) || !(verlet_err > 1e-3))
    return 1;
  return 0;
}
