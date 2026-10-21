/* SPDX-License-Identifier: MIT */
#include "cpmdc_geometry.h"

#include <float.h>
#include <math.h>
#include <minimage.h>

/* Positions and lattice rows use Angstrom; atom order is fixed by the caller. */
double cpmdc_geometry_rmsd(int n, const double *a, const double *b,
                          const double *rows, int periodic) {
    if (n <= 0 || !a || !b || (periodic && !rows)) return DBL_MAX;
    struct mi_cell cell;
    if (periodic && mi_cell_from_con(rows, &cell) != 0) return DBL_MAX;
    double sum = 0.0;
    for (int i = 0; i < n; ++i) {
        double dr[3];
        for (int k = 0; k < 3; ++k) {
            if (!isfinite(a[3*i+k]) || !isfinite(b[3*i+k])) return DBL_MAX;
            dr[k] = b[3*i+k] - a[3*i+k];
        }
        if (periodic && mi_displacement_euclidean(&cell, a+3*i, b+3*i, dr) != 0)
            return DBL_MAX;
        for (int k = 0; k < 3; ++k) sum += dr[k] * dr[k];
    }
    return isfinite(sum) ? sqrt(sum / n) : DBL_MAX;
}
