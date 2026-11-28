#include "cpmdc_csvr.h"

#include <stdio.h>
#include <time.h>

int main(void) {
    cpmdc_csvr_rng rng;
    if (cpmdc_csvr_rng_seed(&rng, 1) != 0) {
        return 1;
    }
    double kinetic = 1.5;
    const int n = 200000;
    struct timespec start;
    struct timespec stop;
    if (timespec_get(&start, TIME_UTC) != TIME_UTC) {
        return 1;
    }
    for (int i = 0; i < n; ++i) {
        double scale = 0.0;
        double work = 0.0;
        if (cpmdc_csvr_scale(kinetic, 9, 1.0, 0.001, 0.1, &rng, &scale, &work) != 0) {
            return 1;
        }
        kinetic *= scale * scale;
        if (!(kinetic > 0.0)) {
            kinetic = 1.5;
        }
    }
    if (timespec_get(&stop, TIME_UTC) != TIME_UTC) {
        return 1;
    }
    const double seconds =
        (double)(stop.tv_sec - start.tv_sec) + (double)(stop.tv_nsec - start.tv_nsec) / 1e9;
    printf("calls %d seconds %.9f kinetic %.9f\n", n, seconds, kinetic);
    return 0;
}
