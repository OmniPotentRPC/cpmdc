/* SPDX-License-Identifier: MIT */
#include "cpmdc_geometry.h"
#include <float.h>
#include <math.h>
#include <stdio.h>

#define CHECK(condition) do { if (!(condition)) { \
  fprintf(stderr, "line %d: %s\n", __LINE__, #condition); return 1; \
} } while (0)

int main(void) {
    const double box[9] = {4,0,0, 0,5,0, 0,0,6};
    const double skew[9] = {4,0,0, .4,5,0, 0,0,6};
    const double zero[3] = {0,0,0};
    const double crossing[3] = {3.9,0,0};
    CHECK(fabs(cpmdc_geometry_rmsd(1, zero, crossing, box, 1) - .1) < 1e-12);
    CHECK(fabs(cpmdc_geometry_rmsd(1, zero, crossing, box, 0) - 3.9) < 1e-12);
    const double q[3] = {-5.375,-2.5,-3.125};
    double oracle = DBL_MAX;
    for (int i=-4; i<=4; ++i)
      for (int j=-4; j<=4; ++j)
        for (int k=-4; k<=4; ++k) {
          double d2=0;
          for (int c=0; c<3; ++c) {
            double d=q[c]+i*skew[c]+j*skew[3+c]+k*skew[6+c];
            d2+=d*d;
          }
          if (d2 < oracle) oracle=d2;
        }
    CHECK(fabs(cpmdc_geometry_rmsd(1, zero, q, skew, 1)-sqrt(oracle)) < 1e-12);
    CHECK(fabs(oracle-15.46625) < 1e-12);
    const double a[6] = {0,0,0, 0,0,0};
    const double b[6] = {3,0,0, 0,4,0};
    CHECK(fabs(cpmdc_geometry_rmsd(2,a,b,box,0)-sqrt(12.5)) < 1e-12);
    CHECK(cpmdc_geometry_rmsd(0,a,b,box,1) == DBL_MAX);
    const double bad[3]={NAN,0,0};
    CHECK(cpmdc_geometry_rmsd(1,zero,bad,box,1) == DBL_MAX);
    const double singular[9]={0};
    CHECK(cpmdc_geometry_rmsd(1,zero,q,singular,1) == DBL_MAX);
    return 0;
}
