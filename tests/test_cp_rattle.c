/* SPDX-License-Identifier: MIT */
#include "cpmdc_cp.h"

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern void dsyev_(const char *jobz, const char *uplo, const int *n, double *a,
                   const int *lda, double *w, double *work, const int *lwork, int *info);

static int g_failed = 0;

static void expect(int cond, const char *msg) {
    if (cond) return;
    fprintf(stderr, "FAIL %s\n", msg);
    g_failed = 1;
}

static double max_abs_diff(int n, const double *a, const double *b) {
    double m = 0.0;
    for (int i = 0; i < n; ++i) {
        double d = fabs(a[i] - b[i]);
        if (d > m) m = d;
    }
    return m;
}

static int test_lowdin_scales_axes(void) {
    /* Columns (2,0,0) and (0,3,0) have overlap diag(4, 9). */
    double C[6] = {2.0, 0.0, 0.0, 0.0, 3.0, 0.0};
    double out[6];
    int sweeps = -1;
    double err = -1.0;
    expect(cpmdc_cp_lowdin(3, 2, C, out, &sweeps, &err) == 0, "lowdin returns");
    expect(sweeps >= 0, "lowdin reports sweeps");
    expect(err < 1e-12, "lowdin overlap");
    expect(fabs(fabs(out[0]) - 1.0) < 1e-12, "first column norm");
    expect(fabs(out[1]) < 1e-12 && fabs(out[2]) < 1e-12, "first column axis");
    expect(fabs(fabs(out[4]) - 1.0) < 1e-12, "second column norm");
    expect(fabs(out[3]) < 1e-12 && fabs(out[5]) < 1e-12, "second column axis");
    return 0;
}

static int test_lowdin_matches_dsyev(void) {
    double C[6] = {1.0, 0.2, 0.05, -0.1, 1.1, 0.3};
    double out[6];
    expect(cpmdc_cp_lowdin(3, 2, C, out, NULL, NULL) == 0, "lowdin of a tilted pair");
    double S[4];
    for (int j = 0; j < 2; ++j) {
        for (int i = 0; i < 2; ++i) {
            double s = 0.0;
            for (int p = 0; p < 3; ++p) s += C[p + i * 3] * C[p + j * 3];
            S[i + j * 2] = s;
        }
    }
    double a[4];
    memcpy(a, S, sizeof(S));
    double w[2], work[32];
    int n = 2, lda = 2, lwork = 32, info = 0;
    const char jobz = 'V', uplo = 'U';
    dsyev_(&jobz, &uplo, &n, a, &lda, w, work, &lwork, &info);
    expect(info == 0, "overlap dsyev");
    double inv[4] = {0.0, 0.0, 0.0, 0.0};
    for (int k = 0; k < 2; ++k) {
        double scale = 1.0 / sqrt(w[k]);
        for (int j = 0; j < 2; ++j)
            for (int i = 0; i < 2; ++i) inv[i + j * 2] += a[i + k * 2] * scale * a[j + k * 2];
    }
    double ref[6];
    for (int j = 0; j < 2; ++j) {
        for (int i = 0; i < 3; ++i) {
            double s = 0.0;
            for (int k = 0; k < 2; ++k) s += C[i + k * 3] * inv[k + j * 2];
            ref[i + j * 3] = s;
        }
    }
    expect(max_abs_diff(6, out, ref) < 1e-10, "Lowdin matches the dsyev inverse square root");
    return 0;
}

static int test_shake_stays_in_span(void) {
    double C[6] = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0};
    double pred[6] = {1.05, 0.02, 0.20, -0.03, 0.97, -0.15};
    double neu[6];
    int iters = -1;
    double err = -1.0;
    expect(cpmdc_cp_shake(3, 2, C, pred, neu, 1e-12, 40, &iters, &err) == 0,
           "shake returns");
    expect(iters >= 1, "shake iterates");
    expect(err < 1e-12, "shake overlap");
    /* The correction lies in the column space of C, so its third row is zero. */
    expect(fabs(neu[2] - pred[2]) < 1e-12, "shake leaves the virtual row of column 0");
    expect(fabs(neu[5] - pred[5]) < 1e-12, "shake leaves the virtual row of column 1");
    double low[6];
    expect(cpmdc_cp_lowdin(3, 2, pred, low, NULL, NULL) == 0, "lowdin of the same prediction");
    expect(max_abs_diff(6, neu, low) > 1e-6, "shake is not the Lowdin polar factor");
    return 0;
}

static int test_rattle_keeps_tangent(void) {
    double C[6] = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0};
    double U[6] = {0.4, 0.1, 0.3, -0.2, 0.5, -0.7};
    /* lambda = 1/2 (C^T U + U^T C) = [[0.4, -0.05], [-0.05, 0.5]],
       and the tangent velocity is U - C lambda. */
    double tangent[6] = {0.0, 0.15, 0.3, -0.15, 0.0, -0.7};
    expect(cpmdc_cp_rattle_velocities(3, 2, C, U) == 0, "rattle returns");
    expect(max_abs_diff(6, U, tangent) < 1e-12, "rattle removes the occupied part");
    double G00 = C[0] * U[0] + C[1] * U[1] + C[2] * U[2];
    double G11 = C[3] * U[3] + C[4] * U[4] + C[5] * U[5];
    double G01 = C[0] * U[3] + C[1] * U[4] + C[2] * U[5];
    double G10 = C[3] * U[0] + C[4] * U[1] + C[5] * U[2];
    expect(fabs(G00) < 1e-12 && fabs(G11) < 1e-12, "diagonal velocity constraints");
    expect(fabs(G01 + G10) < 1e-12, "symmetric velocity constraint");
    return 0;
}

static int read_fixture_positions(const char *path, double ang[9]) {
    FILE *fp = fopen(path, "r");
    if (!fp) return -1;
    char buf[4096];
    size_t n = fread(buf, 1, sizeof(buf) - 1, fp);
    fclose(fp);
    buf[n] = 0;
    const char *p = strstr(buf, "pos");
    if (!p) return -1;
    int got = 0;
    while (*p && got < 9) {
        if (isdigit((unsigned char)*p) || *p == '-' || *p == '.') {
            char *end = NULL;
            double v = strtod(p, &end);
            if (end == p) {
                ++p;
                continue;
            }
            ang[got++] = v;
            p = end;
        } else {
            ++p;
        }
    }
    return got == 9 ? 0 : -1;
}

static int test_water_matches_fixture_and_dsyev(const char *fixture) {
    double ang[9];
    expect(read_fixture_positions(fixture, ang) == 0, "read water fixture");
    double R[9];
    expect(cpmdc_cp_water_positions(R) == 0, "water positions");
    const double bohr = 0.529177210903;
    for (int i = 0; i < 9; ++i)
        expect(fabs(R[i] * bohr - ang[i]) < 1e-12, "positions match the water fixture");

    double H[9];
    expect(cpmdc_cp_water_hamiltonian(R, H) == 0, "hamiltonian");
    double a[9];
    memcpy(a, H, sizeof(a));
    double w[3];
    double work[64];
    int n = 3, lda = 3, lwork = 64, info = 0;
    const char jobz = 'N', uplo = 'U';
    dsyev_(&jobz, &uplo, &n, a, &lda, w, work, &lwork, &info);
    expect(info == 0, "dsyev");
    double ref = 2.0 * (w[0] + w[1]);
    double band = 0.0;
    expect(cpmdc_cp_water_band_energy(R, &band) == 0, "band energy");
    expect(fabs(band - ref) < 1e-10, "band energy matches dsyev");
    double C[6], gs = 0.0;
    expect(cpmdc_cp_water_ground_state(R, C, &gs) == 0, "ground state");
    expect(fabs(gs - ref) < 1e-10, "ground-state energy matches dsyev");
    printf("water_band_energy_ha=%.12f\n", ref);
    return 0;
}

static int test_water_forces(void) {
    double R[9], C[6], energy = 0.0, F[9], f[6];
    expect(cpmdc_cp_water_positions(R) == 0, "positions");
    expect(cpmdc_cp_water_ground_state(R, C, &energy) == 0, "ground state");
    expect(cpmdc_cp_water_force(NULL, R, C, &energy, F, f) == 0, "force");
    double fnorm = 0.0;
    for (int i = 0; i < 9; ++i) fnorm += F[i] * F[i];
    expect(sqrt(fnorm) > 1e-4, "the fixture is not a stationary point of the model");

    const double dx = 1e-5;
    for (int a = 0; a < 3; ++a) {
        for (int k = 0; k < 3; ++k) {
            double Rp[9], Rm[9], ep = 0.0, em = 0.0, Fp[9], Fm[9], fp[6], fm[6];
            memcpy(Rp, R, sizeof(R));
            memcpy(Rm, R, sizeof(R));
            Rp[3 * a + k] += dx;
            Rm[3 * a + k] -= dx;
            expect(cpmdc_cp_water_force(NULL, Rp, C, &ep, Fp, fp) == 0, "force plus");
            expect(cpmdc_cp_water_force(NULL, Rm, C, &em, Fm, fm) == 0, "force minus");
            double num = -(ep - em) / (2.0 * dx);
            expect(fabs(num - F[3 * a + k]) < 1e-7, "ionic force matches a finite difference");
        }
    }
    for (int i = 0; i < 6; ++i) {
        double Cp[6], Cm[6], ep = 0.0, em = 0.0, Fp[9], Fm[9], fp[6], fm[6];
        memcpy(Cp, C, sizeof(C));
        memcpy(Cm, C, sizeof(C));
        Cp[i] += dx;
        Cm[i] -= dx;
        expect(cpmdc_cp_water_force(NULL, R, Cp, &ep, Fp, fp) == 0, "elec plus");
        expect(cpmdc_cp_water_force(NULL, R, Cm, &em, Fm, fm) == 0, "elec minus");
        double num = -(ep - em) / (2.0 * dx);
        expect(fabs(num - f[i]) < 1e-7, "electronic force matches a finite difference");
    }
    return 0;
}

static int propagate(int ortho, int steps, double *R, double *V, double *C, double *U,
                     double *drift, double *rms_back, int *scf_sum, int *ortho_sum) {
    double M[3];
    expect(cpmdc_cp_water_masses(M) == 0, "masses");
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
    expect(work != NULL, "workspace");
    cpmdc_cp_step_info info;
    expect(cpmdc_cp_conserved(&sys, cpmdc_cp_water_force, NULL, R, V, C, U, &info) == 0,
           "initial conserved energy");
    double e0 = info.econs;
    double R0[9], C0[6];
    memcpy(R0, R, sizeof(R0));
    memcpy(C0, C, sizeof(C0));
    *drift = 0.0;
    *scf_sum = 0;
    *ortho_sum = 0;
    for (int s = 0; s < steps; ++s) {
        expect(cpmdc_cp_step(&sys, &opt, cpmdc_cp_water_force, NULL, work, R, V, C, U,
                             &info) == 0,
               "step");
        expect(info.scf_iterations == 0, "a CP step reports no SCF");
        expect(info.ortho_err < 1e-10, "orbitals stay orthonormal");
        double d = fabs(info.econs - e0);
        if (d > *drift) *drift = d;
        *scf_sum += info.scf_iterations;
        *ortho_sum += info.iterations;
    }
    for (int i = 0; i < 9; ++i) V[i] = -V[i];
    for (int i = 0; i < 6; ++i) U[i] = -U[i];
    for (int s = 0; s < steps; ++s) {
        expect(cpmdc_cp_step(&sys, &opt, cpmdc_cp_water_force, NULL, work, R, V, C, U,
                             &info) == 0,
               "return step");
    }
    double sum = 0.0;
    for (int i = 0; i < 9; ++i) {
        double d = R[i] - R0[i];
        sum += d * d;
    }
    *rms_back = sqrt(sum / 9.0);
    double csum = 0.0;
    for (int i = 0; i < 6; ++i) {
        double d = C[i] - C0[i];
        csum += d * d;
    }
    expect(sqrt(csum / 6.0) < 5e-8, "coefficients return");
    cpmdc_cp_work_free(work);
    return 0;
}

static int test_rattle_conserves_and_returns(void) {
    double R[9], V[9], C[6], U[6];
    memset(V, 0, sizeof(V));
    memset(U, 0, sizeof(U));
    expect(cpmdc_cp_water_positions(R) == 0, "positions");
    expect(cpmdc_cp_water_ground_state(R, C, NULL) == 0, "ground state");
    double drift = 0.0, rms = 0.0;
    int scf = 0, ortho = 0;
    const int steps = 80;
    propagate(CPMDC_CP_ORTHO_RATTLE, steps, R, V, C, U, &drift, &rms, &scf, &ortho);
    expect(drift < 2e-7, "RATTLE conserved-energy drift stays under 2e-7 Ha");
    expect(rms < 1e-10, "RATTLE returns the ions after the velocities are reversed");
    expect(scf == 0, "RATTLE runs no SCF");
    expect(ortho >= steps, "SHAKE ran on every forward step");
    printf("rattle_drift_ha=%.3e rattle_return_rms_bohr=%.3e shake_iterations=%d\n", drift,
           rms, ortho);
    return 0;
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "usage: %s tests/force_input_water.capnp.txt\n", argv[0]);
        return 2;
    }
    test_lowdin_scales_axes();
    test_lowdin_matches_dsyev();
    test_shake_stays_in_span();
    test_rattle_keeps_tangent();
    test_water_matches_fixture_and_dsyev(argv[1]);
    test_water_forces();
    test_rattle_conserves_and_returns();
    if (g_failed) {
        fprintf(stderr, "cp rattle tests failed\n");
        return 1;
    }
    printf("cp rattle tests passed\n");
    return 0;
}
