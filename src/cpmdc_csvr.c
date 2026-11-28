#include "cpmdc_csvr.h"

#include <math.h>
#include <stddef.h>
#include <string.h>

_Static_assert(sizeof(cpmdc_csvr_rng) == 48, "cpmdc_csvr_rng size");
_Static_assert(offsetof(cpmdc_csvr_rng, spare_normal) == 32, "spare_normal offset");
_Static_assert(offsetof(cpmdc_csvr_rng, spare_ready) == 40, "spare_ready offset");

static uint64_t rotl(uint64_t x, int k) {
    return (x << k) | (x >> (64 - k));
}

/* Public-domain xoshiro256** of Blackman and Vigna. */
static uint64_t xoshiro256ss(uint64_t s[4]) {
    const uint64_t result = rotl(s[1] * 5, 7) * 9;
    const uint64_t t = s[1] << 17;
    s[2] ^= s[0];
    s[3] ^= s[1];
    s[1] ^= s[2];
    s[0] ^= s[3];
    s[2] ^= t;
    s[3] = rotl(s[3], 45);
    return result;
}

/* Public-domain SplitMix64 of Steele, Lea, and Flood. */
static uint64_t splitmix64(uint64_t *x) {
    uint64_t z = (*x += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

static int seeded(const cpmdc_csvr_rng *rng) {
    return (rng->state[0] | rng->state[1] | rng->state[2] | rng->state[3]) != 0;
}

static int unit_interval(cpmdc_csvr_rng *rng, double *out) {
    for (int attempt = 0; attempt < 64; ++attempt) {
        const uint64_t bits = xoshiro256ss(rng->state) >> 11;
        const double u = ldexp((double)bits, -53);
        if (u > 0.0) {
            *out = u;
            return 0;
        }
    }
    return -1;
}

static int standard_normal(cpmdc_csvr_rng *rng, double *out) {
    if (rng->spare_ready) {
        rng->spare_ready = 0;
        *out = rng->spare_normal;
        return 0;
    }
    double u1 = 0.0;
    double u2 = 0.0;
    if (unit_interval(rng, &u1) != 0 || unit_interval(rng, &u2) != 0) {
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

/* Unit-scale gamma (mean equals the shape) by Marsaglia and Tsang. */
static int gamma_unit(cpmdc_csvr_rng *rng, double shape, double *out) {
    if (!(shape > 0.0) || !isfinite(shape)) {
        return -1;
    }
    if (shape < 1.0) {
        double boosted = 0.0;
        double u = 0.0;
        if (gamma_unit(rng, shape + 1.0, &boosted) != 0 || unit_interval(rng, &u) != 0) {
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
            if (standard_normal(rng, &x) != 0) {
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
        if (unit_interval(rng, &u) != 0 || !isfinite(v) || v <= 0.0) {
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

/* Sum of `count` independent standard-normal squares. Zero when count is 0. */
static int sum_of_squares(cpmdc_csvr_rng *rng, int count, double *out) {
    if (count <= 0) {
        *out = 0.0;
        return 0;
    }
    double unit = 0.0;
    if (gamma_unit(rng, 0.5 * (double)count, &unit) != 0) {
        return -1;
    }
    const double value = 2.0 * unit;
    if (!isfinite(value) || value < 0.0) {
        return -1;
    }
    *out = value;
    return 0;
}

int cpmdc_csvr_rng_seed(cpmdc_csvr_rng *rng, uint64_t seed) {
    if (rng == NULL) {
        return -1;
    }
    memset(rng, 0, sizeof *rng);
    uint64_t mixed = seed;
    for (int i = 0; i < 4; ++i) {
        rng->state[i] = splitmix64(&mixed);
    }
    if (!seeded(rng)) {
        rng->state[0] = 0x9E3779B97F4A7C15ULL;
    }
    return 0;
}

int cpmdc_csvr_scale(double kinetic, int ndof, double temperature, double time_step,
                     double time_constant, cpmdc_csvr_rng *rng, double *scale, double *work) {
    if (rng == NULL || scale == NULL || work == NULL) {
        return -1;
    }
    if (ndof < 1 || !(kinetic > 0.0) || !isfinite(kinetic)) {
        return -1;
    }
    if (!(temperature >= 0.0) || !isfinite(temperature)) {
        return -1;
    }
    if (!(time_step >= 0.0) || !isfinite(time_step)) {
        return -1;
    }
    const int tau_inf = isinf(time_constant) && time_constant > 0.0;
    const int tau_pos = isfinite(time_constant) && time_constant > 0.0;
    if (!tau_inf && !tau_pos) {
        return -1;
    }
    if (!seeded(rng)) {
        return -1;
    }

    double c = 0.0;
    if (tau_inf || time_step == 0.0) {
        c = 1.0;
    } else {
        const double ratio = time_step / time_constant;
        c = isfinite(ratio) ? exp(-ratio) : 0.0;
        if (!isfinite(c)) {
            return -1;
        }
    }

    cpmdc_csvr_rng local = *rng;
    double alpha2 = 0.0;
    const int deterministic = (c == 1.0 || temperature == 0.0);
    if (deterministic) {
        alpha2 = c;
    } else {
        const double s = (1.0 - c) * (0.5 * temperature) / kinetic;
        if (!(s >= 0.0) || !isfinite(s)) {
            return -1;
        }
        double normal = 0.0;
        double rest = 0.0;
        if (standard_normal(&local, &normal) != 0 || sum_of_squares(&local, ndof - 1, &rest) != 0) {
            return -1;
        }
        if (!isfinite(normal) || !isfinite(rest)) {
            return -1;
        }
        const double mixed = sqrt(c) + normal * sqrt(s);
        alpha2 = mixed * mixed + s * rest;
        if (!(alpha2 >= 0.0) || !isfinite(alpha2)) {
            return -1;
        }
    }

    const double scale_out = sqrt(alpha2);
    const double work_out = kinetic * (1.0 - alpha2);
    if (!isfinite(scale_out) || !isfinite(work_out) || scale_out < 0.0) {
        return -1;
    }
    if (!deterministic) {
        *rng = local;
    }
    *scale = scale_out;
    *work = work_out;
    return 0;
}
