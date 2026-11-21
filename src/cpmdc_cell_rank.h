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


#ifdef __cplusplus
}
#endif

#endif
