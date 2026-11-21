#include "cpmdc_cell_rank.h"

#include <float.h>
#include <math.h>

/* python/cpmdc/_rank.py repeats every operation here in the same order.
 * Only +, -, *, /, sqrt, fabs, frexp and ldexp appear, each exact or
 * correctly rounded in binary64, and the build turns off contraction, so
 * both languages reach the same decision from the same bits. */

#define CELL_JACOBI_SWEEPS 64

static int scale_exponent(const double cell[9]) {
  double largest = 0.0;
  int exponent = 0;
  for (int i = 0; i < 9; ++i)
    if (fabs(cell[i]) > largest)
      largest = fabs(cell[i]);
  (void)frexp(largest, &exponent);
  return exponent;
}

static int rotate_columns(double a[9], int p, int q) {
  double alpha = 0.0;
  double beta = 0.0;
  double gamma = 0.0;
  for (int row = 0; row < 3; ++row) {
    double ap = a[3 * row + p];
    double aq = a[3 * row + q];
    alpha = alpha + ap * ap;
    beta = beta + aq * aq;
    gamma = gamma + ap * aq;
  }
  if (gamma == 0.0 || fabs(gamma) <= DBL_EPSILON * sqrt(alpha) * sqrt(beta))
    return 0;
  double zeta = (beta - alpha) / (2.0 * gamma);
  double t;
  if (fabs(zeta) > 1.0)
    t = 1.0 / (fabs(zeta) * (1.0 + sqrt(1.0 + 1.0 / (zeta * zeta))));
  else
    t = 1.0 / (fabs(zeta) + sqrt(1.0 + zeta * zeta));
  if (zeta < 0.0)
    t = -t;
  double cosine = 1.0 / sqrt(1.0 + t * t);
  double sine = cosine * t;
  for (int row = 0; row < 3; ++row) {
    double ap = a[3 * row + p];
    double aq = a[3 * row + q];
    a[3 * row + p] = cosine * ap - sine * aq;
    a[3 * row + q] = sine * ap + cosine * aq;
  }
  return 1;
}

/* |det| of a row-major 3 by 3 matrix by LU with partial pivoting. */
static double abs_det_lu(const double in[9]) {
  double m[9];
  for (int i = 0; i < 9; ++i)
    m[i] = in[i];
  for (int k = 0; k < 3; ++k) {
    int pivot = k;
    double best = fabs(m[3 * k + k]);
    for (int i = k + 1; i < 3; ++i)
      if (fabs(m[3 * i + k]) > best) {
        pivot = i;
        best = fabs(m[3 * i + k]);
      }
    if (best == 0.0)
      return 0.0;
    if (pivot != k)
      for (int j = 0; j < 3; ++j) {
        double swap = m[3 * k + j];
        m[3 * k + j] = m[3 * pivot + j];
        m[3 * pivot + j] = swap;
      }
    for (int i = k + 1; i < 3; ++i) {
      double factor = m[3 * i + k] / m[3 * k + k];
      for (int j = k + 1; j < 3; ++j)
        m[3 * i + j] = m[3 * i + j] - factor * m[3 * k + j];
    }
  }
  return fabs(m[0] * m[4] * m[8]);
}

int cpmdc_cell_rank3_state(const double cell[9]) {
  double a[9];
  int any = 0;
  for (int i = 0; i < 9; ++i) {
    if (!isfinite(cell[i]))
      return CPMDC_CELL_NONFINITE;
    if (cell[i] != 0.0)
      any = 1;
  }
  if (!any)
    return CPMDC_CELL_ZERO;
  int exponent = scale_exponent(cell);
  for (int i = 0; i < 9; ++i)
    a[i] = ldexp(cell[i], -exponent);
  double scaled[9];
  for (int i = 0; i < 9; ++i)
    scaled[i] = a[i];

  for (int sweep = 0; sweep < CELL_JACOBI_SWEEPS; ++sweep) {
    int rotated = rotate_columns(a, 0, 1);
    rotated |= rotate_columns(a, 0, 2);
    rotated |= rotate_columns(a, 1, 2);
    if (!rotated)
      break;
  }

  double singular[3];
  double largest = 0.0;
  for (int col = 0; col < 3; ++col) {
    double sum = 0.0;
    for (int row = 0; row < 3; ++row)
      sum = sum + a[3 * row + col] * a[3 * row + col];
    singular[col] = sqrt(sum);
    if (singular[col] > largest)
      largest = singular[col];
  }
  double cutoff = largest * 3.0 * DBL_EPSILON;
  for (int col = 0; col < 3; ++col)
    if (!(singular[col] > cutoff))
      return CPMDC_CELL_RANK_DEFICIENT;

  double determinant = abs_det_lu(scaled);
  if (determinant == 0.0)
    return CPMDC_CELL_RANK_DEFICIENT;
  int det_exponent = 0;
  (void)frexp(determinant, &det_exponent);
  /* The volume is in [2^(E-1), 2^E) with E below, and DBL_MIN is
   * 2^(DBL_MIN_EXP-1). */
  if (det_exponent + 3 * exponent < DBL_MIN_EXP)
    return CPMDC_CELL_VOLUME_UNDERFLOW;
  return CPMDC_CELL_FULL_RANK;
}

int cpmdc_cell_volume3(const double cell[9], double *volume) {
  double scaled[9];
  int any = 0;
  *volume = 0.0;
  for (int i = 0; i < 9; ++i) {
    if (!isfinite(cell[i]))
      return -1;
    if (cell[i] != 0.0)
      any = 1;
  }
  if (!any)
    return 0;
  int exponent = scale_exponent(cell);
  for (int i = 0; i < 9; ++i)
    scaled[i] = ldexp(cell[i], -exponent);
  double determinant = abs_det_lu(scaled);
  if (determinant == 0.0)
    return 0;
  int det_exponent = 0;
  double mantissa = frexp(determinant, &det_exponent);
  if (det_exponent + 3 * exponent > DBL_MAX_EXP) {
    *volume = HUGE_VAL;
    return -1;
  }
  *volume = ldexp(mantissa, det_exponent + 3 * exponent);
  return isfinite(*volume) ? 0 : -1;
}
