/* SPDX-License-Identifier: MIT */
#include "cpmdc_damp.h"

#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* CPV/src/cpr.f90 sets fccc = 1/(1+frice) unless the electrons are on a
 * Nose chain or on steepest descent. CPV/src/runcp.f90 then uses
 * verl1 = 2*fccc, verl2 = 1-verl1, verl3 = fccc, and
 * emaver = (dt^2/emass) * ema0bg * verl3.
 * wave_verlet writes cm = verl1*c0 + verl2*cm + emaver*grad.
 */

static int weights(double friction, double *fccc, double *verl1, double *verl2,
                   double *verl3) {
  double denom;
  double fc;
  double v1;
  double v2;
  double v3;

  if (!isfinite(friction) || !(friction > -1.0))
    return -1;
  denom = 1.0 + friction;
  fc = 1.0 / denom;
  v1 = 2.0 * fc;
  v2 = 1.0 - v1;
  v3 = fc;
  if (!isfinite(fc) || !isfinite(v1) || !isfinite(v2) || !isfinite(v3))
    return -1;
  *fccc = fc;
  *verl1 = v1;
  *verl2 = v2;
  *verl3 = v3;
  return 0;
}

static int finite_n(int n, const double *v) {
  int i;

  if (v == NULL)
    return 0;
  for (i = 0; i < n; ++i) {
    if (!isfinite(v[i]))
      return 0;
  }
  return 1;
}

int cpmdc_damp_coefficients(double friction, double *fccc, double *verl1,
                            double *verl2, double *verl3) {
  if (fccc == NULL || verl1 == NULL || verl2 == NULL || verl3 == NULL)
    return -1;
  return weights(friction, fccc, verl1, verl2, verl3);
}

int cpmdc_damp_step(int n, const double *c0_re, const double *c0_im,
                    double *cm_re, double *cm_im, const double *force_re,
                    const double *force_im, const double *ema0bg,
                    double friction, double dt2_over_emass) {
  double fccc;
  double verl1;
  double verl2;
  double verl3;
  double *out;
  int i;

  if (n < 1 || c0_re == NULL || c0_im == NULL || cm_re == NULL || cm_im == NULL)
    return -1;
  if (force_re == NULL || force_im == NULL)
    return -1;
  if ((size_t)n > SIZE_MAX / (2u * sizeof(double)))
    return -1;
  if (!isfinite(dt2_over_emass))
    return -1;
  if (weights(friction, &fccc, &verl1, &verl2, &verl3) != 0)
    return -1;
  if (!finite_n(n, c0_re) || !finite_n(n, c0_im) || !finite_n(n, cm_re) ||
      !finite_n(n, cm_im) || !finite_n(n, force_re) || !finite_n(n, force_im))
    return -1;
  if (ema0bg != NULL && !finite_n(n, ema0bg))
    return -1;

  out = malloc((size_t)n * 2u * sizeof(double));
  if (out == NULL)
    return -1;
  for (i = 0; i < n; ++i) {
    double ema = ema0bg == NULL ? 1.0 : ema0bg[i];
    double scale = fccc * dt2_over_emass * ema;
    double re = verl1 * c0_re[i] + verl2 * cm_re[i] + scale * force_re[i];
    double im = verl1 * c0_im[i] + verl2 * cm_im[i] + scale * force_im[i];
    if (!isfinite(re) || !isfinite(im)) {
      free(out);
      return -1;
    }
    out[(size_t)i] = re;
    out[(size_t)n + (size_t)i] = im;
  }
  memcpy(cm_re, out, (size_t)n * sizeof(double));
  memcpy(cm_im, out + n, (size_t)n * sizeof(double));
  free(out);
  return 0;
}

int cpmdc_damp_swap(int n, double *c0_re, double *c0_im, double *cm_re,
                    double *cm_im) {
  int i;

  if (n < 1 || c0_re == NULL || c0_im == NULL || cm_re == NULL || cm_im == NULL)
    return -1;
  for (i = 0; i < n; ++i) {
    double re = c0_re[i];
    double im = c0_im[i];
    c0_re[i] = cm_re[i];
    c0_im[i] = cm_im[i];
    cm_re[i] = re;
    cm_im[i] = im;
  }
  return 0;
}

int cpmdc_damp_grease(double *friction, double grease) {
  double next;

  if (friction == NULL || !isfinite(*friction) || !isfinite(grease))
    return -1;
  /* A friction at or below -1 makes the step denominator 1+friction
   * zero or negative. The stored product has to stay above that pole.
   */
  if (!(*friction > -1.0))
    return -1;
  next = *friction * grease;
  if (!isfinite(next) || !(next > -1.0))
    return -1;
  *friction = next;
  return 0;
}

int cpmdc_damp_advance(int n, double *c0_re, double *c0_im, double *cm_re,
                       double *cm_im, const double *force_re,
                       const double *force_im, const double *ema0bg,
                       double friction, double dt2_over_emass) {
  if (cpmdc_damp_step(n, c0_re, c0_im, cm_re, cm_im, force_re, force_im, ema0bg,
                      friction, dt2_over_emass) != 0)
    return -1;
  return cpmdc_damp_swap(n, c0_re, c0_im, cm_re, cm_im);
}
