#ifndef CPMDC_CELL_RANK_H
#define CPMDC_CELL_RANK_H

#ifdef __cplusplus
extern "C" {
#endif

enum {
  CPMDC_CELL_VOLUME_UNDERFLOW = -3,
  CPMDC_CELL_NONFINITE = -2,
  CPMDC_CELL_RANK_DEFICIENT = -1,
  CPMDC_CELL_ZERO = 0,
  CPMDC_CELL_FULL_RANK = 1
};

/* Classify a row-major binary64 3 by 3 cell. Rank 3 means every singular
 * value from a one-sided Jacobi iteration on the cell, scaled by a power of
 * two, is above the largest times 3 times DBL_EPSILON: the
 * numpy.linalg.matrix_rank cutoff. cpmdc._rank repeats the iteration. */
int cpmdc_cell_rank3_state(const double cell[9]);

/* |det cell| by LU with partial pivoting on the cell scaled by a power of
 * two. Returns 0 with *volume set, or -1 for a non-finite entry or a
 * volume above DBL_MAX. An all-zero or singular cell gives 0. */
int cpmdc_cell_volume3(const double cell[9], double *volume);

#ifdef __cplusplus
}
#endif

#endif
