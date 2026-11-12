/* SPDX-License-Identifier: MIT */
#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 199309L
#endif

#include "cpmdc_cp.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static double wall_seconds(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

static void reset_state(double *R, double *V, double *C, double *U) {
    memset(V, 0, 9 * sizeof(double));
    memset(U, 0, 6 * sizeof(double));
    if (cpmdc_cp_water_positions(R) != 0 || cpmdc_cp_water_ground_state(R, C, NULL) != 0) {
        fprintf(stderr, "water model setup failed\n");
        exit(1);
    }
}

static int run_path(const char *name, int ortho, int steps, int repeats) {
    double M[3], R[9], V[9], C[6], U[6];
    if (cpmdc_cp_water_masses(M) != 0) return 1;
    reset_state(R, V, C, U);
    cpmdc_cp_system sys = {
        .n_basis = CPMDC_CP_WATER_BASIS,
        .n_state = CPMDC_CP_WATER_STATES,
        .n_atom = CPMDC_CP_WATER_ATOMS,
        .mu = 400.0,
        .mass = M,
    };
    cpmdc_cp_step_opts opt = {
        .dt = 0.25,
        .ortho = ortho,
        .shake_tol = 1e-12,
        .shake_max = 40,
    };
    cpmdc_cp_work *work = cpmdc_cp_work_new(sys.n_basis, sys.n_state, sys.n_atom);
    if (!work) return 1;
    cpmdc_cp_step_info info;
    if (cpmdc_cp_conserved(&sys, cpmdc_cp_water_force, NULL, R, V, C, U, &info) != 0) return 1;
    double e0 = info.econs;
    double drift = 0.0;
    long scf = 0;
    long ortho_iters = 0;
    for (int s = 0; s < steps; ++s) {
        if (cpmdc_cp_step(&sys, &opt, cpmdc_cp_water_force, NULL, work, R, V, C, U, &info) != 0)
            return 1;
        double d = info.econs - e0;
        if (d < 0.0) d = -d;
        if (d > drift) drift = d;
        scf += info.scf_iterations;
        ortho_iters += info.iterations;
    }
    double t0 = wall_seconds();
    for (int rep = 0; rep < repeats; ++rep) {
        reset_state(R, V, C, U);
        for (int s = 0; s < steps; ++s) {
            if (cpmdc_cp_step(&sys, &opt, cpmdc_cp_water_force, NULL, work, R, V, C, U, &info) !=
                0)
                return 1;
        }
    }
    double elapsed = wall_seconds() - t0;
    double timed = (double)repeats * (double)steps;
    printf("%s %d %d %.2f %.1f %.9f %.9e %ld %.6f %.9e\n", name, steps, repeats, opt.dt, sys.mu,
           elapsed, elapsed / timed, scf, (double)ortho_iters / (double)steps, drift);
    cpmdc_cp_work_free(work);
    return 0;
}

int main(int argc, char **argv) {
    int steps = 2000;
    int repeats = 40;
    if (argc > 1) steps = atoi(argv[1]);
    if (argc > 2) repeats = atoi(argv[2]);
    if (steps < 1 || repeats < 1) return 2;
    printf("path steps repeats dt_au mu time_s time_per_step_s scf_iterations "
           "ortho_iterations_mean econs_drift_ha\n");
    if (run_path("lowdin", CPMDC_CP_ORTHO_LOWDIN, steps, repeats) != 0) return 1;
    if (run_path("rattle", CPMDC_CP_ORTHO_RATTLE, steps, repeats) != 0) return 1;
    return 0;
}
