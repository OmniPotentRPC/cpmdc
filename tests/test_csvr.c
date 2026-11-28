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

static uint64_t replay_rotl(uint64_t x, int k) {
    return (x << k) | (x >> (64 - k));
}

static uint64_t replay_xoshiro(uint64_t s[4]) {
    const uint64_t result = replay_rotl(s[1] * 5, 7) * 9;
    const uint64_t t = s[1] << 17;
    s[2] ^= s[0];
    s[3] ^= s[1];
    s[1] ^= s[2];
    s[0] ^= s[3];
    s[2] ^= t;
    s[3] = replay_rotl(s[3], 45);
    return result;
}

static int replay_unit(cpmdc_csvr_rng *rng, double *out) {
    for (int attempt = 0; attempt < 64; ++attempt) {
        const uint64_t bits = replay_xoshiro(rng->state) >> 11;
        const double u = ldexp((double)bits, -53);
        if (u > 0.0) {
            *out = u;
            return 0;
        }
    }
    return -1;
}

static int replay_normal(cpmdc_csvr_rng *rng, double *out) {
    if (rng->spare_ready) {
        rng->spare_ready = 0;
        *out = rng->spare_normal;
        return 0;
    }
    double u1 = 0.0;
    double u2 = 0.0;
    if (replay_unit(rng, &u1) != 0 || replay_unit(rng, &u2) != 0) {
        return -1;
    }
    const double mag = sqrt(-2.0 * log(u1));
    const double ang = 2.0 * acos(-1.0) * u2;
    rng->spare_normal = mag * sin(ang);
    rng->spare_ready = 1;
    *out = mag * cos(ang);
    if (!isfinite(*out) || !isfinite(rng->spare_normal)) {
        return -1;
    }
    return 0;
}

/* Unit-scale gamma, copied so the test can replay the library draws. */
static int replay_gamma(cpmdc_csvr_rng *rng, double shape, double *out) {
    if (!(shape > 0.0) || !isfinite(shape)) {
        return -1;
    }
    if (shape < 1.0) {
        double boosted = 0.0;
        double u = 0.0;
        if (replay_gamma(rng, shape + 1.0, &boosted) != 0 || replay_unit(rng, &u) != 0) {
            return -1;
        }
        const double value = boosted * pow(u, 1.0 / shape);
        if (!isfinite(value) || value < 0.0) {
            return -1;
        }
        *out = value;
        return 0;
    }
    const double d = shape - (1.0 / 3.0);
    const double c = 1.0 / sqrt(9.0 * d);
    for (int draw = 0; draw < 1000000; ++draw) {
        double x = 0.0;
        double v = 0.0;
        int accepted_v = 0;
        for (int inner = 0; inner < 1000000; ++inner) {
            if (replay_normal(rng, &x) != 0) {
                return -1;
            }
            v = 1.0 + c * x;
            if (v > 0.0) {
                accepted_v = 1;
                break;
            }
        }
        if (!accepted_v) {
            return -1;
        }
        v = v * v * v;
        double u = 0.0;
        if (replay_unit(rng, &u) != 0 || !isfinite(v) || v <= 0.0) {
            return -1;
        }
        const double x2 = x * x;
        if (u < 1.0 - 0.0331 * x2 * x2 || log(u) < 0.5 * x2 + d * (1.0 - v + log(v))) {
            const double value = d * v;
            if (!isfinite(value) || value < 0.0) {
                return -1;
            }
            *out = value;
            return 0;
        }
    }
    return -1;
}

static int replay_sum(cpmdc_csvr_rng *rng, int count, double *out) {
    if (count <= 0) {
        *out = 0.0;
        return 0;
    }
    double unit = 0.0;
    if (replay_gamma(rng, 0.5 * (double)count, &unit) != 0) {
        return -1;
    }
    const double value = 2.0 * unit;
    if (!isfinite(value) || value < 0.0) {
        return -1;
    }
    *out = value;
    return 0;
}

/* Equation (22) of the arXiv HTML, (A7) in the arXiv PDF and the journal. */
static int paper_scale(const cpmdc_csvr_rng *seeded, double kinetic, int ndof, double temperature,
                        double time_step, double time_constant, double *scale, double *cross) {
    cpmdc_csvr_rng rng = *seeded;
    const double c = exp(-time_step / time_constant);
    const double ratio = (0.5 * temperature) / kinetic;
    double normal = 0.0;
    double rest = 0.0;
    if (replay_normal(&rng, &normal) != 0 || replay_sum(&rng, ndof - 1, &rest) != 0) {
        return -1;
    }
    const double one_c = 1.0 - c;
    *cross = 2.0 * exp(-0.5 * time_step / time_constant) * sqrt(ratio * one_c) * normal;
    const double alpha2 =
        c + ratio * one_c * (normal * normal + rest) + *cross;
    if (!(alpha2 >= 0.0) || !isfinite(alpha2)) {
        return -1;
    }
    *scale = sqrt(alpha2);
    return 0;
}

static int finite_tau(int ndof) {
    const double dt = 1.0;
    const double tau = 10.0;
    const double kinetic = 1.25;
    const double temperature = 0.8;
    const double c = exp(-dt / tau);
    if (!(c > 0.5)) {
        return fail("finite tau decayed");
    }
    for (uint64_t seed = 1; seed < 64; ++seed) {
        cpmdc_csvr_rng rng;
        if (cpmdc_csvr_rng_seed(&rng, seed) != 0) {
            return fail("finite seed");
        }
        double expected = 0.0;
        double cross = 0.0;
        if (paper_scale(&rng, kinetic, ndof, temperature, dt, tau, &expected, &cross) != 0) {
            return fail("finite paper");
        }
        if (!(fabs(cross) > 1e-3)) {
            continue;
        }
        double scale = 0.0;
        double work = 0.0;
        if (cpmdc_csvr_scale(kinetic, ndof, temperature, dt, tau, &rng, &scale, &work) != 0) {
            return fail("finite library");
        }
        const double diff = fabs(scale - expected);
        printf("csvr_stat finite_tau ndof %d diff %.9e cross %.9e c %.9e\n", ndof, diff, cross, c);
        if (diff > 1e-12) {
            return fail("finite tau mismatch");
        }
        return 0;
    }
    return fail("finite tau found no cross term");
}

static int copied_generator(void) {
    cpmdc_csvr_rng rng;
    if (cpmdc_csvr_rng_seed(&rng, 17) != 0) {
        return fail("copy seed");
    }
    double kinetic = 1.5;
    for (int step = 0; step < 4; ++step) {
        double scale = 0.0;
        double work = 0.0;
        if (cpmdc_csvr_scale(kinetic, 3, 1.0, 1.0, 10.0, &rng, &scale, &work) != 0) {
            return fail("copy warmup");
        }
        kinetic *= scale * scale;
    }
    cpmdc_csvr_rng copy;
    memcpy(&copy, &rng, sizeof copy);
    for (int step = 0; step < 6; ++step) {
        double left = 0.0;
        double right = 0.0;
        double left_work = 0.0;
        double right_work = 0.0;
        if (cpmdc_csvr_scale(kinetic, 3, 1.0, 1.0, 10.0, &rng, &left, &left_work) != 0 ||
            cpmdc_csvr_scale(kinetic, 3, 1.0, 1.0, 10.0, &copy, &right, &right_work) != 0) {
            return fail("copy step");
        }
        if (left != right || left_work != right_work) {
            return fail("copy diverged");
        }
        kinetic *= left * left;
    }
    printf("csvr_stat copy_replay ok\n");
    return 0;
}

static int positive_root(void) {
    const double dt = 1.0;
    const double tau = 10.0;
    const double kinetic = 0.05;
    const double temperature = 1.0;
    const double c = exp(-dt / tau);
    const double s = (1.0 - c) * (0.5 * temperature) / kinetic;
    for (uint64_t seed = 1; seed < 10000; ++seed) {
        cpmdc_csvr_rng rng;
        if (cpmdc_csvr_rng_seed(&rng, seed) != 0) {
            return fail("root seed");
        }
        cpmdc_csvr_rng draw = rng;
        double normal = 0.0;
        if (replay_normal(&draw, &normal) != 0) {
            return fail("root draw");
        }
        const double signed_term = sqrt(c) + normal * sqrt(s);
        if (!(signed_term < 0.0)) {
            continue;
        }
        double scale = 0.0;
        double work = 0.0;
        if (cpmdc_csvr_scale(kinetic, 1, temperature, dt, tau, &rng, &scale, &work) != 0) {
            return fail("root library");
        }
        const double positive = fabs(signed_term);
        printf("csvr_stat positive_root signed %.9e scale %.9e\n", signed_term, scale);
        if (!(scale > 0.0) || fabs(scale - positive) > 1e-12 || scale == signed_term) {
            return fail("root was not the positive square root");
        }
        return 0;
    }
    return fail("no negative root");
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

    if (finite_tau(1) != 0 || finite_tau(4) != 0 || copied_generator() != 0 || positive_root() != 0) {
        return 1;
    }
    if (ideal_gas(1, 7) != 0 || ideal_gas(2, 9) != 0 || ideal_gas(3, 11) != 0 ||
        ideal_gas(6, 13) != 0 || harmonic() != 0) {
        return 1;
    }
    printf("csvr_stat ok\n");
    return 0;
}
