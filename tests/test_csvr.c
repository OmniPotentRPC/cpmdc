#include "cpmdc_csvr.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

static int fail(const char *msg) {
    fprintf(stderr, "csvr: %s\n", msg);
    return 1;
}

static int ideal_gas(int ndof, uint64_t seed) {
    const int n = 20000;
    const double kT = 1.0;
    const double expected_mean = 0.5 * (double)ndof * kT;
    const double expected_var = 0.5 * (double)ndof * kT * kT;
    const double excess_kurtosis = 12.0 / (double)ndof;
    cpmdc_csvr_rng rng;
    if (cpmdc_csvr_rng_seed(&rng, seed) != 0) {
        return fail("ideal seed");
    }
    const double k0 = expected_mean;
    double k = k0;
    double sum_work = 0.0;
    double mean = 0.0;
    double m2 = 0.0;
    for (int i = 1; i <= n; ++i) {
        double scale = 0.0;
        double work = 0.0;
        if (cpmdc_csvr_scale(k, ndof, kT, 1.0, 1e-8, &rng, &scale, &work) != 0) {
            fprintf(stderr, "csvr: ideal ndof %d step %d rejected\n", ndof, i);
            return 1;
        }
        if (!(scale >= 0.0) || !isfinite(scale) || !isfinite(work)) {
            return fail("ideal nonfinite");
        }
        const double knew = k * scale * scale;
        if (fabs((knew + work) - k) > 1e-9 * (1.0 + fabs(k))) {
            return fail("ideal work");
        }
        sum_work += work;
        k = knew;
        const double delta = k - mean;
        mean += delta / (double)i;
        m2 += delta * (k - mean);
    }
    const double var = m2 / (double)(n - 1);
    const double mean_limit = 8.0 * sqrt(expected_var / (double)n);
    const double var_limit = 10.0 * expected_var * sqrt((excess_kurtosis + 2.0) / (double)n);
    printf(
        "csvr_stat ideal ndof %d mean %.9f exp_mean %.9f mean_limit %.9f "
        "var %.9f exp_var %.9f var_limit %.9f\n",
        ndof, mean, expected_mean, mean_limit, var, expected_var, var_limit);
    if (fabs(mean - expected_mean) > mean_limit || fabs(var - expected_var) > var_limit) {
        fprintf(stderr, "csvr: ideal ndof %d outside the sampling bounds\n", ndof);
        return 1;
    }
    const double telescope = fabs((k + sum_work) - k0);
    printf("csvr_stat telescope ndof %d abs_err %.9e\n", ndof, telescope);
    if (telescope > 1e-9 * (1.0 + fabs(k0))) {
        return fail("ideal telescope");
    }
    return 0;
}

static int harmonic(void) {
    cpmdc_csvr_rng rng;
    if (cpmdc_csvr_rng_seed(&rng, 11) != 0) {
        return fail("harmonic seed");
    }
    const double dt = 0.01;
    const double kT = 1.0;
    const double tau = 0.2;
    double x = 0.0;
    double v = 1.0;
    double work_sum = 0.0;
    double sum_abs = 0.0;
    const double h0 = 0.5 * v * v;
    double max_drift = 0.0;
    for (int step = 0; step < 20000; ++step) {
        double a = -x;
        v += 0.5 * dt * a;
        x += dt * v;
        a = -x;
        v += 0.5 * dt * a;
        const double k = 0.5 * v * v;
        if (!(k > 0.0)) {
            return fail("harmonic kinetic vanished");
        }
        double scale = 0.0;
        double work = 0.0;
        if (cpmdc_csvr_scale(k, 1, kT, dt, tau, &rng, &scale, &work) != 0) {
            return fail("harmonic rejected");
        }
        const double knew = k * scale * scale;
        if (fabs((knew + work) - k) > 1e-9 * (1.0 + k)) {
            return fail("harmonic work");
        }
        v *= scale;
        work_sum += work;
        sum_abs += fabs(work);
        const double h = 0.5 * v * v + 0.5 * x * x + work_sum;
        const double drift = fabs(h - h0);
        if (drift > max_drift) {
            max_drift = drift;
        }
    }
    printf("csvr_stat harmonic_max_drift %.9e sum_abs_work %.9e\n", max_drift, sum_abs);
    if (!(max_drift < 1e-2) || !(sum_abs > max_drift)) {
        return fail("harmonic drift");
    }
    return 0;
}

static int same_trajectory(uint64_t seed_a, uint64_t seed_b, double *max_diff, int *differ) {
    cpmdc_csvr_rng a;
    cpmdc_csvr_rng b;
    if (cpmdc_csvr_rng_seed(&a, seed_a) != 0 || cpmdc_csvr_rng_seed(&b, seed_b) != 0) {
        return fail("repro seed");
    }
    double ka = 1.25;
    double kb = 1.25;
    *max_diff = 0.0;
    *differ = 0;
    for (int i = 0; i < 32; ++i) {
        double sa = 0.0;
        double sb = 0.0;
        double wa = 0.0;
        double wb = 0.0;
        if (cpmdc_csvr_scale(ka, 5, 0.8, 0.5, 2.0, &a, &sa, &wa) != 0 ||
            cpmdc_csvr_scale(kb, 5, 0.8, 0.5, 2.0, &b, &sb, &wb) != 0) {
            return fail("repro step");
        }
        const double diff = fabs(sa - sb);
        if (diff > *max_diff) {
            *max_diff = diff;
        }
        if (diff > 1e-9 || fabs(wa - wb) > 1e-9) {
            *differ = 1;
        }
        ka *= sa * sa;
        kb *= sb * sb;
    }
    return 0;
}

static int rejected(cpmdc_csvr_rng *rng, double kinetic, int ndof, double temperature,
                    double time_step, double time_constant) {
    cpmdc_csvr_rng copy = *rng;
    double scale = 11.0;
    double work = 13.0;
    if (cpmdc_csvr_scale(kinetic, ndof, temperature, time_step, time_constant, rng, &scale,
                         &work) != -1) {
        return 0;
    }
    if (scale != 11.0 || work != 13.0) {
        return 0;
    }
    return memcmp(rng, &copy, sizeof copy) == 0;
}

int main(void) {
    printf("csvr_stat rng_bytes %zu\n", sizeof(cpmdc_csvr_rng));

    cpmdc_csvr_rng rng;
    cpmdc_csvr_rng untouched;
    if (cpmdc_csvr_rng_seed(&rng, 99) != 0 || cpmdc_csvr_rng_seed(&untouched, 99) != 0) {
        return fail("nve seed");
    }
    double scale = 3.0;
    double work = 4.0;
    if (cpmdc_csvr_scale(1.25, 6, 0.8, 1.0, INFINITY, &rng, &scale, &work) != 0) {
        return fail("infinite tau");
    }
    printf("csvr_stat nve_inf_scale %.16e nve_inf_work %.16e\n", scale, work);
    if (fabs(scale - 1.0) > 1e-12 || fabs(work) > 1e-12 || memcmp(&rng, &untouched, sizeof rng) != 0) {
        return fail("infinite tau is not nve");
    }
    scale = 3.0;
    work = 4.0;
    if (cpmdc_csvr_scale(1.25, 6, 0.8, 1.0, 1e20, &rng, &scale, &work) != 0) {
        return fail("huge tau");
    }
    printf("csvr_stat nve_huge_scale %.16e nve_huge_work %.16e\n", scale, work);
    if (fabs(scale - 1.0) > 1e-12 || fabs(work) > 1e-12 || memcmp(&rng, &untouched, sizeof rng) != 0) {
        return fail("huge tau is not nve");
    }
    if (cpmdc_csvr_scale(1.25, 6, 0.8, 0.0, 1.0, &rng, &scale, &work) != 0 || fabs(scale - 1.0) > 1e-12) {
        return fail("zero time step");
    }

    if (cpmdc_csvr_rng_seed(&rng, 4) != 0) {
        return fail("cold seed");
    }
    if (cpmdc_csvr_scale(2.0, 4, 0.0, 2.0, 4.0, &rng, &scale, &work) != 0) {
        return fail("zero temperature");
    }
    const double cold_scale = exp(-0.25);
    const double cold_work = 2.0 * (1.0 - exp(-0.5));
    printf("csvr_stat cold_scale %.16e cold_work %.16e\n", scale, work);
    if (fabs(scale - cold_scale) > 1e-12 || fabs(work - cold_work) > 1e-12) {
        return fail("zero temperature decay");
    }

    double repro_diff = 0.0;
    int differ = 0;
    if (same_trajectory(5, 5, &repro_diff, &differ) != 0) {
        return 1;
    }
    printf("csvr_stat repro_max_abs_diff %.9e\n", repro_diff);
    if (differ || repro_diff > 1e-9) {
        return fail("same seed diverged");
    }
    double other_diff = 0.0;
    int other = 0;
    if (same_trajectory(5, 6, &other_diff, &other) != 0) {
        return 1;
    }
    printf("csvr_stat other_seed_max_abs_diff %.9e\n", other_diff);
    if (!other) {
        return fail("different seeds matched");
    }

    if (cpmdc_csvr_rng_seed(NULL, 1) != -1) {
        return fail("null seed");
    }
    scale = 11.0;
    work = 13.0;
    if (cpmdc_csvr_rng_seed(&rng, 1) != 0) {
        return fail("seed");
    }
    if (cpmdc_csvr_scale(1.0, 1, 1.0, 1.0, 1.0, NULL, &scale, &work) != -1 || scale != 11.0 ||
        work != 13.0) {
        return fail("null rng");
    }
    if (cpmdc_csvr_scale(1.0, 1, 1.0, 1.0, 1.0, &rng, NULL, &work) != -1 ||
        cpmdc_csvr_scale(1.0, 1, 1.0, 1.0, 1.0, &rng, &scale, NULL) != -1) {
        return fail("null output");
    }
    cpmdc_csvr_rng blank;
    memset(&blank, 0, sizeof blank);
    scale = 11.0;
    work = 13.0;
    if (cpmdc_csvr_scale(1.0, 1, 1.0, 1.0, 1.0, &blank, &scale, &work) != -1 || scale != 11.0) {
        return fail("unseeded rng");
    }
    if (cpmdc_csvr_rng_seed(&rng, 0) != 0 ||
        cpmdc_csvr_scale(1.0, 3, 1.0, 0.1, 1.0, &rng, &scale, &work) != 0 || !(scale > 0.0)) {
        return fail("seed zero");
    }

    if (cpmdc_csvr_rng_seed(&rng, 8) != 0) {
        return fail("reject seed");
    }
    const double bad_k[] = {0.0, -1.0, NAN, INFINITY};
    for (int i = 0; i < 4; ++i) {
        if (!rejected(&rng, bad_k[i], 3, 1.0, 1.0, 1.0)) {
            return fail("bad kinetic");
        }
    }
    if (!rejected(&rng, 1.0, 0, 1.0, 1.0, 1.0) || !rejected(&rng, 1.0, -3, 1.0, 1.0, 1.0)) {
        return fail("bad ndof");
    }
    if (!rejected(&rng, 1.0, 3, -0.1, 1.0, 1.0) || !rejected(&rng, 1.0, 3, NAN, 1.0, 1.0) ||
        !rejected(&rng, 1.0, 3, INFINITY, 1.0, 1.0)) {
        return fail("bad temperature");
    }
    if (!rejected(&rng, 1.0, 3, 1.0, -1.0, 1.0) || !rejected(&rng, 1.0, 3, 1.0, NAN, 1.0) ||
        !rejected(&rng, 1.0, 3, 1.0, INFINITY, 1.0)) {
        return fail("bad time step");
    }
    if (!rejected(&rng, 1.0, 3, 1.0, 1.0, 0.0) || !rejected(&rng, 1.0, 3, 1.0, 1.0, -1.0) ||
        !rejected(&rng, 1.0, 3, 1.0, 1.0, NAN) || !rejected(&rng, 1.0, 3, 1.0, 1.0, -INFINITY)) {
        return fail("bad time constant");
    }

    if (ideal_gas(1, 7) != 0 || ideal_gas(2, 9) != 0 || ideal_gas(3, 11) != 0 ||
        ideal_gas(6, 13) != 0 || harmonic() != 0) {
        return 1;
    }
    printf("csvr_stat ok\n");
    return 0;
}
