/* SPDX-License-Identifier: MIT */
#include "cpmdc_cp.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

/*
 * Velocity Verlet with orthonormality constraints.
 *
 * The constraint force on the coefficients lies in the column space of the
 * coefficients at the start of the step. SHAKE finds the symmetric multiplier
 * that puts the updated coefficients back on C^T C = I, and the same multiplier
 * is added to the half-step velocity, which is Andersen's RATTLE. The second
 * half-step velocity is then projected so that C^T U + U^T C = 0.
 *
 * Spectral Löwdin replaces that multiplier with the polar factor of the
 * unconstrained coefficients. The factor comes from a full diagonalisation of
 * the overlap, and the correction is not confined to the old column space.
 */

enum { CPMDC_CP_JACOBI_SWEEPS = 80 };

struct cpmdc_cp_work {
    int n_basis;
    int n_state;
    int n_atom;
    double *C_pred;
    double *C_old;
    double *lambda;
    double *S;
    double *Q;
    double *eval;
    double *evec;
    double *tmp;
    double *f;
    double *F;
    double *Acopy;
    double *block;
};

static int bad_shape(int n_basis, int n_state, int n_atom) {
    if (n_basis < 1 || n_state < 1 || n_atom < 1) return 1;
    if (n_state > n_basis) return 1;
    if (n_basis > 4096 || n_state > 1024 || n_atom > 4096) return 1;
    return 0;
}

cpmdc_cp_work *cpmdc_cp_work_new(int n_basis, int n_state, int n_atom) {
    if (bad_shape(n_basis, n_state, n_atom)) return NULL;
    cpmdc_cp_work *work = calloc(1, sizeof(*work));
    if (!work) return NULL;
    size_t nb = (size_t)n_basis;
    size_t ns = (size_t)n_state;
    size_t coeff = nb * ns;
    size_t sq = ns * ns;
    size_t ion = 3u * (size_t)n_atom;
    size_t n = coeff * 4u + sq * 5u + ns + ion;
    work->block = calloc(n, sizeof(double));
    if (!work->block) {
        free(work);
        return NULL;
    }
    double *p = work->block;
    work->C_pred = p; p += coeff;
    work->C_old = p; p += coeff;
    work->tmp = p; p += coeff;
    work->f = p; p += coeff;
    work->lambda = p; p += sq;
    work->S = p; p += sq;
    work->Q = p; p += sq;
    work->evec = p; p += sq;
    work->Acopy = p; p += sq;
    work->eval = p; p += ns;
    work->F = p;
    work->n_basis = n_basis;
    work->n_state = n_state;
    work->n_atom = n_atom;
    return work;
}

void cpmdc_cp_work_free(cpmdc_cp_work *work) {
    if (!work) return;
    free(work->block);
    free(work);
}

static void gemm(int m, int n, int k, const double *A, const double *B, double *C) {
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < m; ++i) {
            double s = 0.0;
            for (int p = 0; p < k; ++p) s += A[i + p * m] * B[p + j * k];
            C[i + j * m] = s;
        }
    }
}

static void gram(int n_basis, int n_state, const double *A, const double *B, double *G) {
    for (int j = 0; j < n_state; ++j) {
        for (int i = 0; i < n_state; ++i) {
            double s = 0.0;
            for (int p = 0; p < n_basis; ++p) s += A[p + i * n_basis] * B[p + j * n_basis];
            G[i + j * n_state] = s;
        }
    }
}

static double residual(int n, const double *Q) {
    double m = 0.0;
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i) {
            double d = Q[i + j * n] - (i == j ? 1.0 : 0.0);
            if (fabs(d) > m) m = fabs(d);
        }
    }
    return m;
}

/* Cyclic Jacobi. A is symmetric column-major and is overwritten. evec's
 * columns are the eigenvectors. Returns the number of sweeps, or -1. */
static int jacobi(int n, double *A, double *eval, double *evec) {
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i) evec[i + j * n] = (i == j) ? 1.0 : 0.0;
    int sweeps = 0;
    for (; sweeps < CPMDC_CP_JACOBI_SWEEPS; ++sweeps) {
        double off = 0.0;
        for (int p = 0; p < n; ++p)
            for (int q = p + 1; q < n; ++q) off += fabs(A[p + q * n]);
        if (off <= 1e-15 * (double)n) break;
        for (int p = 0; p < n; ++p) {
            for (int q = p + 1; q < n; ++q) {
                double apq = A[p + q * n];
                if (fabs(apq) <= 1e-18) continue;
                double app = A[p + p * n];
                double aqq = A[q + q * n];
                double tau = (aqq - app) / (2.0 * apq);
                double t = copysign(1.0 / (fabs(tau) + sqrt(1.0 + tau * tau)), tau);
                double c = 1.0 / sqrt(1.0 + t * t);
                double s = t * c;
                A[p + p * n] = app - t * apq;
                A[q + q * n] = aqq + t * apq;
                A[p + q * n] = 0.0;
                A[q + p * n] = 0.0;
                for (int r = 0; r < n; ++r) {
                    if (r == p || r == q) continue;
                    double arp = A[r + p * n];
                    double arq = A[r + q * n];
                    double np = c * arp - s * arq;
                    double nq = s * arp + c * arq;
                    A[r + p * n] = np;
                    A[p + r * n] = np;
                    A[r + q * n] = nq;
                    A[q + r * n] = nq;
                }
                for (int r = 0; r < n; ++r) {
                    double vp = evec[r + p * n];
                    double vq = evec[r + q * n];
                    evec[r + p * n] = c * vp - s * vq;
                    evec[r + q * n] = s * vp + c * vq;
                }
            }
        }
    }
    if (sweeps == CPMDC_CP_JACOBI_SWEEPS) return -1;
    for (int p = 0; p < n; ++p) eval[p] = A[p + p * n];
    return sweeps;
}

static int inverse_sqrt(cpmdc_cp_work *work, const double *S, double *out, int *sweeps) {
    int n = work->n_state;
    size_t sq = (size_t)n * (size_t)n;
    memcpy(work->Acopy, S, sq * sizeof(double));
    for (int j = 0; j < n; ++j) {
        for (int i = j + 1; i < n; ++i) {
            double v = 0.5 * (work->Acopy[i + j * n] + work->Acopy[j + i * n]);
            work->Acopy[i + j * n] = v;
            work->Acopy[j + i * n] = v;
        }
    }
    int sw = jacobi(n, work->Acopy, work->eval, work->evec);
    if (sw < 0) return -1;
    memset(out, 0, sq * sizeof(double));
    for (int k = 0; k < n; ++k) {
        if (!(work->eval[k] > 1e-18)) return -1;
        double scale = 1.0 / sqrt(work->eval[k]);
        for (int j = 0; j < n; ++j)
            for (int i = 0; i < n; ++i)
                out[i + j * n] += work->evec[i + k * n] * scale * work->evec[j + k * n];
    }
    if (sweeps) *sweeps = sw;
    return 0;
}

static int lowdin_w(cpmdc_cp_work *work, const double *C, double *C_out, int *sweeps,
                    double *err) {
    int nb = work->n_basis;
    int ns = work->n_state;
    gram(nb, ns, C, C, work->S);
    if (inverse_sqrt(work, work->S, work->Q, sweeps) != 0) return -1;
    double *dest = (C_out == C) ? work->tmp : C_out;
    gemm(nb, ns, ns, C, work->Q, dest);
    if (dest != C_out) memcpy(C_out, dest, (size_t)nb * (size_t)ns * sizeof(double));
    gram(nb, ns, C_out, C_out, work->S);
    double e = residual(ns, work->S);
    if (err) *err = e;
    return e <= 1e-8 ? 0 : -1;
}

static int shake_w(cpmdc_cp_work *work, const double *C_old, const double *C_pred,
                   double *C_new, double tol, int max_iter, int *iterations, double *err) {
    int nb = work->n_basis;
    int ns = work->n_state;
    size_t coeff = (size_t)nb * (size_t)ns;
    size_t sq = (size_t)ns * (size_t)ns;
    const double *pred = C_pred;
    if (C_new == C_pred) {
        memcpy(work->C_pred, C_pred, coeff * sizeof(double));
        pred = work->C_pred;
    }
    if (C_new == C_old) return -1;
    memset(work->lambda, 0, sq * sizeof(double));
    if (tol <= 0.0) tol = 1e-12;
    if (max_iter < 1) max_iter = 1;
    double e = 1.0;
    int iter = 0;
    for (; iter < max_iter; ++iter) {
        gemm(nb, ns, ns, C_old, work->lambda, work->tmp);
        for (size_t i = 0; i < coeff; ++i) C_new[i] = pred[i] + work->tmp[i];
        gram(nb, ns, C_new, C_new, work->Q);
        e = residual(ns, work->Q);
        if (e <= tol) {
            ++iter;
            break;
        }
        for (int j = 0; j < ns; ++j) {
            for (int i = 0; i < ns; ++i) {
                double q = 0.5 * (work->Q[i + j * ns] + work->Q[j + i * ns]);
                double step = 0.5 * ((i == j ? 1.0 : 0.0) - q);
                work->lambda[i + j * ns] += step;
            }
        }
    }
    if (iterations) *iterations = iter;
    if (err) *err = e;
    return e <= tol ? 0 : -1;
}

static int rattle_w(cpmdc_cp_work *work, const double *C, double *U) {
    int nb = work->n_basis;
    int ns = work->n_state;
    gram(nb, ns, C, U, work->Q);
    for (int j = 0; j < ns; ++j) {
        for (int i = 0; i < ns; ++i)
            work->lambda[i + j * ns] =
                0.5 * (work->Q[i + j * ns] + work->Q[j + i * ns]);
    }
    gemm(nb, ns, ns, C, work->lambda, work->tmp);
    size_t coeff = (size_t)nb * (size_t)ns;
    for (size_t i = 0; i < coeff; ++i) U[i] -= work->tmp[i];
    return 0;
}

int cpmdc_cp_lowdin(int n_basis, int n_state, const double *C, double *C_out, int *sweeps,
                    double *err) {
    if (!C || !C_out || bad_shape(n_basis, n_state, 1)) return -1;
    cpmdc_cp_work *work = cpmdc_cp_work_new(n_basis, n_state, 1);
    if (!work) return -1;
    int rc = lowdin_w(work, C, C_out, sweeps, err);
    cpmdc_cp_work_free(work);
    return rc;
}

int cpmdc_cp_shake(int n_basis, int n_state, const double *C_old, const double *C_pred,
                   double *C_new, double tol, int max_iter, int *iterations, double *err) {
    if (!C_old || !C_pred || !C_new || bad_shape(n_basis, n_state, 1)) return -1;
    cpmdc_cp_work *work = cpmdc_cp_work_new(n_basis, n_state, 1);
    if (!work) return -1;
    int rc = shake_w(work, C_old, C_pred, C_new, tol, max_iter, iterations, err);
    cpmdc_cp_work_free(work);
    return rc;
}

int cpmdc_cp_rattle_velocities(int n_basis, int n_state, const double *C, double *U) {
    if (!C || !U || bad_shape(n_basis, n_state, 1)) return -1;
    cpmdc_cp_work *work = cpmdc_cp_work_new(n_basis, n_state, 1);
    if (!work) return -1;
    int rc = rattle_w(work, C, U);
    cpmdc_cp_work_free(work);
    return rc;
}

static int kinetic(const cpmdc_cp_system *sys, const double *V, const double *U,
                   double *ekin_ion, double *ekinc) {
    double ion = 0.0;
    for (int a = 0; a < sys->n_atom; ++a) {
        if (!(sys->mass[a] > 0.0)) return -1;
        double v2 = 0.0;
        for (int k = 0; k < 3; ++k) {
            double v = V[3 * a + k];
            v2 += v * v;
        }
        ion += 0.5 * sys->mass[a] * v2;
    }
    double elec = 0.0;
    size_t coeff = (size_t)sys->n_basis * (size_t)sys->n_state;
    for (size_t i = 0; i < coeff; ++i) elec += U[i] * U[i];
    *ekin_ion = ion;
    *ekinc = 0.5 * sys->mu * elec;
    return 0;
}

static int system_ok(const cpmdc_cp_system *sys, cpmdc_cp_work *work) {
    if (!sys || !sys->mass || !(sys->mu > 0.0)) return 0;
    if (bad_shape(sys->n_basis, sys->n_state, sys->n_atom)) return 0;
    if (work && (work->n_basis != sys->n_basis || work->n_state != sys->n_state ||
                 work->n_atom != sys->n_atom))
        return 0;
    return 1;
}

int cpmdc_cp_conserved(const cpmdc_cp_system *sys, cpmdc_cp_force_fn force, void *ctx,
                       const double *R, const double *V, const double *C, const double *U,
                       cpmdc_cp_step_info *info) {
    if (!system_ok(sys, NULL) || !force || !R || !V || !C || !U || !info) return -1;
    size_t coeff = (size_t)sys->n_basis * (size_t)sys->n_state;
    double *F = calloc(3u * (size_t)sys->n_atom, sizeof(double));
    double *f = calloc(coeff, sizeof(double));
    if (!F || !f) {
        free(F);
        free(f);
        return -1;
    }
    double energy = 0.0;
    int rc = force(ctx, R, C, &energy, F, f);
    double ekin = 0.0, ekinc = 0.0;
    if (rc == 0) rc = kinetic(sys, V, U, &ekin, &ekinc);
    cpmdc_cp_work *work = rc == 0 ? cpmdc_cp_work_new(sys->n_basis, sys->n_state, sys->n_atom)
                                  : NULL;
    double ortho = 0.0;
    if (rc == 0 && work) {
        gram(sys->n_basis, sys->n_state, C, C, work->Q);
        ortho = residual(sys->n_state, work->Q);
    } else if (rc == 0) {
        rc = -1;
    }
    if (rc == 0) {
        info->iterations = 0;
        info->scf_iterations = 0;
        info->energy = energy;
        info->ekin_ion = ekin;
        info->ekinc = ekinc;
        info->econs = energy + ekin + ekinc;
        info->ortho_err = ortho;
    }
    cpmdc_cp_work_free(work);
    free(F);
    free(f);
    return rc;
}

static void kick(const cpmdc_cp_system *sys, double dt, const double *F, const double *f,
                 double *V, double *U) {
    for (int a = 0; a < sys->n_atom; ++a)
        for (int k = 0; k < 3; ++k)
            V[3 * a + k] += 0.5 * dt * F[3 * a + k] / sys->mass[a];
    size_t coeff = (size_t)sys->n_basis * (size_t)sys->n_state;
    for (size_t i = 0; i < coeff; ++i) U[i] += 0.5 * dt * f[i] / sys->mu;
}

int cpmdc_cp_step(const cpmdc_cp_system *sys, const cpmdc_cp_step_opts *opt,
                  cpmdc_cp_force_fn force, void *ctx, cpmdc_cp_work *work, double *R,
                  double *V, double *C, double *U, cpmdc_cp_step_info *info) {
    if (!system_ok(sys, work) || !opt || !force || !work || !R || !V || !C || !U) return -1;
    if (!(opt->dt > 0.0)) return -1;
    if (opt->ortho != CPMDC_CP_ORTHO_LOWDIN && opt->ortho != CPMDC_CP_ORTHO_RATTLE) return -1;
    double energy = 0.0;
    double ekin_chk = 0.0;
    double ekinc_chk = 0.0;
    if (force(ctx, R, C, &energy, work->F, work->f) != 0) return -1;
    if (kinetic(sys, V, U, &ekin_chk, &ekinc_chk) != 0) return -1;
    (void)ekin_chk;
    (void)ekinc_chk;

    kick(sys, opt->dt, work->F, work->f, V, U);
    size_t coeff = (size_t)sys->n_basis * (size_t)sys->n_state;
    memcpy(work->C_old, C, coeff * sizeof(double));
    for (size_t i = 0; i < coeff; ++i) work->C_pred[i] = work->C_old[i] + opt->dt * U[i];
    for (int a = 0; a < sys->n_atom; ++a)
        for (int k = 0; k < 3; ++k) R[3 * a + k] += opt->dt * V[3 * a + k];

    int iterations = 0;
    double oerr = 0.0;
    int rc;
    if (opt->ortho == CPMDC_CP_ORTHO_RATTLE)
        rc = shake_w(work, work->C_old, work->C_pred, C, opt->shake_tol, opt->shake_max,
                     &iterations, &oerr);
    else
        rc = lowdin_w(work, work->C_pred, C, &iterations, &oerr);
    if (rc != 0) return rc;
    for (size_t i = 0; i < coeff; ++i) U[i] += (C[i] - work->C_pred[i]) / opt->dt;

    if (force(ctx, R, C, &energy, work->F, work->f) != 0) return -1;
    kick(sys, opt->dt, work->F, work->f, V, U);
    if (rattle_w(work, C, U) != 0) return -1;

    double ekin = 0.0, ekinc = 0.0;
    if (kinetic(sys, V, U, &ekin, &ekinc) != 0) return -1;
    gram(sys->n_basis, sys->n_state, C, C, work->Q);
    if (info) {
        info->iterations = iterations;
        info->scf_iterations = 0;
        info->energy = energy;
        info->ekin_ion = ekin;
        info->ekinc = ekinc;
        info->econs = energy + ekin + ekinc;
        info->ortho_err = residual(sys->n_state, work->Q);
    }
    return 0;
}

/* CODATA 2018: one Bohr is 0.529177210903 Å, and one u is 1822.888486209 m_e. */
static const double CPMDC_CP_BOHR_PER_ANGSTROM = 0.529177210903;
static const double CPMDC_CP_ELECTRON_PER_U = 1822.888486209;
static const double CPMDC_CP_WATER_ANGSTROM[9] = {0.0, 0.0,  0.1173, 0.0, 0.7572,
                                                  -0.4692, 0.0, -0.7572, -0.4692};
static const double CPMDC_CP_WATER_EPS[3] = {-0.60, -0.20, -0.20};
static const double CPMDC_CP_WATER_T0 = -0.25;
static const double CPMDC_CP_WATER_ALPHA = 1.50;
static const double CPMDC_CP_WATER_K_OH = 0.50;
static const double CPMDC_CP_WATER_K_HH = 0.10;

int cpmdc_cp_water_positions(double R[9]) {
    if (!R) return -1;
    for (int i = 0; i < 9; ++i) R[i] = CPMDC_CP_WATER_ANGSTROM[i] / CPMDC_CP_BOHR_PER_ANGSTROM;
    return 0;
}

int cpmdc_cp_water_masses(double M[3]) {
    if (!M) return -1;
    M[0] = 16.0 * CPMDC_CP_ELECTRON_PER_U;
    M[1] = 1.0 * CPMDC_CP_ELECTRON_PER_U;
    M[2] = 1.0 * CPMDC_CP_ELECTRON_PER_U;
    return 0;
}

static double pair_distance(const double *R, int i, int j) {
    double s = 0.0;
    for (int k = 0; k < 3; ++k) {
        double d = R[3 * i + k] - R[3 * j + k];
        s += d * d;
    }
    return sqrt(s);
}

static double reference_bond(void) {
    double R[9];
    for (int i = 0; i < 9; ++i) R[i] = CPMDC_CP_WATER_ANGSTROM[i] / CPMDC_CP_BOHR_PER_ANGSTROM;
    return pair_distance(R, 0, 1);
}

static double hopping(double r, double r0) {
    return CPMDC_CP_WATER_T0 * exp(-CPMDC_CP_WATER_ALPHA * (r - r0));
}

int cpmdc_cp_water_hamiltonian(const double *R, double H[9]) {
    if (!R || !H) return -1;
    double r0 = reference_bond();
    memset(H, 0, 9 * sizeof(double));
    for (int i = 0; i < 3; ++i) H[i + i * 3] = CPMDC_CP_WATER_EPS[i];
    for (int i = 0; i < 3; ++i) {
        for (int j = i + 1; j < 3; ++j) {
            double r = pair_distance(R, i, j);
            if (!(r > 1e-12)) return -1;
            double t = hopping(r, r0);
            H[i + j * 3] = t;
            H[j + i * 3] = t;
        }
    }
    return 0;
}

static void sort3(const double eval[3], int order[3]) {
    order[0] = 0;
    order[1] = 1;
    order[2] = 2;
    for (int i = 0; i < 3; ++i) {
        for (int j = i + 1; j < 3; ++j) {
            if (eval[order[j]] < eval[order[i]]) {
                int tmp = order[i];
                order[i] = order[j];
                order[j] = tmp;
            }
        }
    }
}

int cpmdc_cp_water_ground_state(const double *R, double *C, double *energy) {
    if (!R || !C) return -1;
    double H[9], eval[3], evec[9];
    if (cpmdc_cp_water_hamiltonian(R, H) != 0) return -1;
    if (jacobi(3, H, eval, evec) < 0) return -1;
    int order[3];
    sort3(eval, order);
    for (int col = 0; col < 2; ++col) {
        int src = order[col];
        int pivot = 0;
        double best = fabs(evec[src * 3]);
        for (int i = 1; i < 3; ++i) {
            double a = fabs(evec[i + src * 3]);
            if (a > best) {
                best = a;
                pivot = i;
            }
        }
        double sign = evec[pivot + src * 3] >= 0.0 ? 1.0 : -1.0;
        for (int i = 0; i < 3; ++i) C[i + col * 3] = sign * evec[i + src * 3];
    }
    if (energy) *energy = 2.0 * (eval[order[0]] + eval[order[1]]);
    return 0;
}

int cpmdc_cp_water_band_energy(const double *R, double *energy) {
    if (!energy) return -1;
    double C[6];
    return cpmdc_cp_water_ground_state(R, C, energy);
}

int cpmdc_cp_water_force(void *ctx, const double *R, const double *C, double *energy,
                         double *F_ion, double *f_elec) {
    (void)ctx;
    if (!R || !C || !energy || !F_ion || !f_elec) return -1;
    double H[9];
    if (cpmdc_cp_water_hamiltonian(R, H) != 0) return -1;
    double P[9];
    for (int j = 0; j < 3; ++j) {
        for (int i = 0; i < 3; ++i) {
            double s = 0.0;
            for (int k = 0; k < 2; ++k) s += C[i + k * 3] * C[j + k * 3];
            P[i + j * 3] = s;
        }
    }
    double pot = 0.0;
    for (int j = 0; j < 3; ++j)
        for (int i = 0; i < 3; ++i) pot += H[i + j * 3] * P[i + j * 3];
    *energy = 2.0 * pot;

    memset(F_ion, 0, 9 * sizeof(double));
    double Rref[9];
    if (cpmdc_cp_water_positions(Rref) != 0) return -1;
    double r_oh = pair_distance(Rref, 0, 1);
    double r_hh = pair_distance(Rref, 1, 2);
    const int pairs[3][2] = {{0, 1}, {0, 2}, {1, 2}};
    const double kpair[3] = {CPMDC_CP_WATER_K_OH, CPMDC_CP_WATER_K_OH, CPMDC_CP_WATER_K_HH};
    const double rpair[3] = {r_oh, r_oh, r_hh};
    for (int p = 0; p < 3; ++p) {
        int i = pairs[p][0];
        int j = pairs[p][1];
        double r = pair_distance(R, i, j);
        if (!(r > 1e-12)) return -1;
        double ext = r - rpair[p];
        *energy += 0.5 * kpair[p] * ext * ext;
        for (int k = 0; k < 3; ++k) {
            double u = (R[3 * i + k] - R[3 * j + k]) / r;
            double fpair = kpair[p] * ext * u;
            F_ion[3 * i + k] -= fpair;
            F_ion[3 * j + k] += fpair;
        }
    }
    double r0 = reference_bond();
    for (int i = 0; i < 3; ++i) {
        for (int j = i + 1; j < 3; ++j) {
            double r = pair_distance(R, i, j);
            if (!(r > 1e-12)) return -1;
            double tp = -CPMDC_CP_WATER_ALPHA * hopping(r, r0);
            double pij = P[i + j * 3];
            for (int k = 0; k < 3; ++k) {
                double u = (R[3 * i + k] - R[3 * j + k]) / r;
                double dE = 4.0 * pij * tp * u;
                F_ion[3 * i + k] -= dE;
                F_ion[3 * j + k] += dE;
            }
        }
    }
    for (int col = 0; col < 2; ++col) {
        for (int row = 0; row < 3; ++row) {
            double s = 0.0;
            for (int p = 0; p < 3; ++p) s += H[row + p * 3] * C[p + col * 3];
            f_elec[row + col * 3] = -4.0 * s;
        }
    }
    return 0;
}
