#include "cpmdc_cell_rank.h"

#include <float.h>
#include <math.h>

static void rotate_columns(double a[9], int p, int q) {
  double alpha = 0.0;
  double beta = 0.0;
  double gamma = 0.0;
  for (int row = 0; row < 3; ++row) {
    double ap = a[3 * row + p];
    double aq = a[3 * row + q];
    alpha += ap * ap;
    beta += aq * aq;
    gamma += ap * aq;
  }
  if (gamma == 0.0)
    return;
  double tau = (beta - alpha) / (2.0 * gamma);
  double t = copysign(1.0, tau) / (fabs(tau) + hypot(1.0, tau));
  double cosine = 1.0 / hypot(1.0, t);
  double sine = cosine * t;
  for (int row = 0; row < 3; ++row) {
    double ap = a[3 * row + p];
    double aq = a[3 * row + q];
    a[3 * row + p] = cosine * ap - sine * aq;
    a[3 * row + q] = sine * ap + cosine * aq;
  }
}

int cpmdc_cell_rank3_state(const double cell[9]) {
  double a[9];
  double scale = 0.0;
  int any = 0;
  for (int i = 0; i < 9; ++i) {
    if (!isfinite(cell[i]))
      return CPMDC_CELL_NONFINITE;
    if (cell[i] != 0.0)
      any = 1;
    scale = fmax(scale, fabs(cell[i]));
  }
  if (!any)
    return CPMDC_CELL_ZERO;
  for (int i = 0; i < 9; ++i)
    a[i] = cell[i] / scale;

  for (int sweep = 0; sweep < 24; ++sweep) {
    rotate_columns(a, 0, 1);
    rotate_columns(a, 0, 2);
    rotate_columns(a, 1, 2);
  }

  double smallest = HUGE_VAL;
  double largest = 0.0;
  for (int col = 0; col < 3; ++col) {
    double norm = 0.0;
    for (int row = 0; row < 3; ++row)
      norm = hypot(norm, a[3 * row + col]);
    smallest = fmin(smallest, norm);
    largest = fmax(largest, norm);
  }
  if (!(smallest > largest * DBL_EPSILON * 3.0))
    return CPMDC_CELL_RANK_DEFICIENT;

  double determinant =
      a[0] * (a[4] * a[8] - a[5] * a[7]) -
      a[1] * (a[3] * a[8] - a[5] * a[6]) +
      a[2] * (a[3] * a[7] - a[4] * a[6]);
  if (determinant == 0.0)
    return CPMDC_CELL_RANK_DEFICIENT;
  if (3.0 * log(scale) + log(fabs(determinant)) < log(DBL_MIN))
    return CPMDC_CELL_VOLUME_UNDERFLOW;
  return CPMDC_CELL_FULL_RANK;
}
