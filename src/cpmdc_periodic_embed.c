/* Reference-build periodic electrostatic embedding.
 * Conducting-boundary Ewald field of the same Gaussian charges the isolated
 * kernel uses. A shell that fits in the Miller cap is the plain reciprocal
 * sum. A narrower Gaussian is split into a short-range real-space sum and a
 * screened reciprocal sum, both cut at the same decay, so the tail is not
 * dropped. Not a public symbol: libcpmdc hides it. */
#include "cpmdc_params.h"

#include <math.h>

/* Keep a reciprocal vector when sigma^2 |G|^2 / 2 is at most this. */
#define CPMDC_EWALD_HALF_DECAY 36.0
/* And when every Miller index lies in [-CAP, CAP]. */
#define CPMDC_EWALD_INDEX_CAP 48

static double dot3(const double a[3], const double b[3]) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static void cross3(const double a[3], const double b[3], double out[3]) {
  out[0] = a[1] * b[2] - a[2] * b[1];
  out[1] = a[2] * b[0] - a[0] * b[2];
  out[2] = a[0] * b[1] - a[1] * b[0];
}

/* Clamp in double before the cast. A narrow Gaussian makes the unclamped
 * ratio larger than INT_MAX, and converting that to int is undefined. */
static int miller_bound(double ratio) {
  int bound;
  if (!isfinite(ratio) || !(ratio > 0.0))
    return 0;
  if (ratio > (double)CPMDC_EWALD_INDEX_CAP)
    ratio = (double)CPMDC_EWALD_INDEX_CAP;
  bound = (int)ceil(ratio);
  if (bound < 0)
    return 0;
  if (bound > CPMDC_EWALD_INDEX_CAP)
    return CPMDC_EWALD_INDEX_CAP;
  return bound;
}

/* Field factor of [erf(alpha r) - erf(beta r)] / r, and (1/r) d(lambda)/dr.
 * E = q * lambda * r_vec. The series is the small-r form: the closed form
 * subtracts two terms of size O(1/r^2). */
static void split_lambda(double u, double alpha, double beta, double *lambda,
                         double *lambda_over_u) {
  const double inv_sqrt_pi = 1.0 / sqrt(acos(-1.0));
  double big = alpha > beta ? alpha : beta;
  if (!(u > 0.0) || u * big <= 0.125) {
    double u2 = u * u;
    double a2 = alpha * alpha;
    double b2 = beta * beta;
    double a_odd = alpha * a2;
    double b_odd = beta * b2;
    double lam = 0.0;
    double lov = 0.0;
    double pow_u = 1.0;
    double fact = 1.0;
    for (int m = 0; m <= 8; ++m) {
      double sign = (m & 1) ? -1.0 : 1.0;
      double sign_l = (m & 1) ? 1.0 : -1.0;
      lam += sign * (a_odd - b_odd) * pow_u / (fact * (double)(2 * m + 3));
      lov += sign_l * (a_odd * a2 - b_odd * b2) * pow_u /
             (fact * (double)(2 * m + 5));
      pow_u *= u2;
      fact *= (double)(m + 1);
      a_odd *= a2;
      b_odd *= b2;
    }
    *lambda = 4.0 * inv_sqrt_pi * lam;
    *lambda_over_u = 8.0 * inv_sqrt_pi * lov;
    return;
  }
  {
    double au = alpha * u;
    double bu = beta * u;
    double erf_diff = erf(au) - erf(bu);
    double exp_a = au > 20.0 ? 0.0 : exp(-au * au);
    double exp_b = bu > 20.0 ? 0.0 : exp(-bu * bu);
    double two = 2.0 * inv_sqrt_pi;
    double a_coef = two * (alpha * exp_a - beta * exp_b);
    double c_coef =
        two * (alpha * alpha * alpha * exp_a - beta * beta * beta * exp_b);
    double u2 = u * u;
    double u3 = u2 * u;
    double u4 = u2 * u2;
    *lambda = erf_diff / u3 - a_coef / u2;
    *lambda_over_u =
        3.0 * a_coef / u4 - 3.0 * erf_diff / (u4 * u) + 2.0 * c_coef / u2;
  }
}

/* Miller window of images with |r| <= rcut. Returns 2 when that window does
 * not fit in the cap, so the caller can fail instead of dropping images. */
static int real_window(const double recip[3][3], const double delta[3],
                       double rcut, int lo[3], int hi[3]) {
  double pi = acos(-1.0);
  for (int axis = 0; axis < 3; ++axis) {
    double blen = sqrt(dot3(recip[axis], recip[axis]));
    double center = dot3(recip[axis], delta) / (2.0 * pi);
    double span = blen * rcut / (2.0 * pi);
    if (!isfinite(center) || !isfinite(span) || fabs(center) > 2.0e9 ||
        span > (double)CPMDC_EWALD_INDEX_CAP)
      return 2;
    lo[axis] = (int)ceil(center - span - 1.0e-9);
    hi[axis] = (int)floor(center + span + 1.0e-9);
  }
  return 0;
}

static void add_split_reciprocal(const double recip[3][3], const int nbound[3],
                                 double volume, double beta, double q,
                                 const double delta[3], int deriv,
                                 double field[3], double dfield[3][3]) {
  double four_pi = 4.0 * acos(-1.0);
  for (int h = -nbound[0]; h <= nbound[0]; ++h) {
    for (int k = -nbound[1]; k <= nbound[1]; ++k) {
      for (int l = -nbound[2]; l <= nbound[2]; ++l) {
        double gvec[3];
        double g2;
        double decay;
        double phase;
        double pref;
        if (h == 0 && k == 0 && l == 0)
          continue;
        for (int c = 0; c < 3; ++c)
          gvec[c] = (double)h * recip[0][c] + (double)k * recip[1][c] +
                    (double)l * recip[2][c];
        g2 = dot3(gvec, gvec);
        if (!(g2 > 0.0) || !(beta > 0.0))
          continue;
        decay = g2 / (4.0 * beta * beta);
        if (decay > CPMDC_EWALD_HALF_DECAY)
          continue;
        phase = dot3(gvec, delta);
        pref = (four_pi * q / volume) * exp(-decay) / g2;
        if (!deriv) {
          double s = sin(phase);
          for (int c = 0; c < 3; ++c)
            field[c] += pref * gvec[c] * s;
        } else {
          double cosph = cos(phase);
          for (int b = 0; b < 3; ++b)
            for (int a = 0; a < 3; ++a)
              dfield[b][a] += pref * gvec[b] * cosph * gvec[a];
        }
      }
    }
  }
}

static int add_split_real(const double rows[3][3], const double recip[3][3],
                          double alpha, double beta, double q,
                          const double delta[3], int deriv, double field[3],
                          double dfield[3][3]) {
  double rcut;
  int lo[3];
  int hi[3];
  int rc;
  if (!(alpha > beta))
    return 0;
  rcut = sqrt(CPMDC_EWALD_HALF_DECAY) / beta;
  rc = real_window(recip, delta, rcut, lo, hi);
  if (rc != 0)
    return rc;
  for (int n1 = lo[0]; n1 <= hi[0]; ++n1) {
    for (int n2 = lo[1]; n2 <= hi[1]; ++n2) {
      for (int n3 = lo[2]; n3 <= hi[2]; ++n3) {
        double image[3];
        double disp[3];
        double radius;
        double lambda;
        double lambda_over_u;
        for (int c = 0; c < 3; ++c)
          image[c] = (double)n1 * rows[0][c] + (double)n2 * rows[1][c] +
                     (double)n3 * rows[2][c];
        for (int c = 0; c < 3; ++c)
          disp[c] = delta[c] - image[c];
        radius = sqrt(dot3(disp, disp));
        if (radius > rcut)
          continue;
        split_lambda(radius, alpha, beta, &lambda, &lambda_over_u);
        if (!deriv) {
          for (int c = 0; c < 3; ++c)
            field[c] += q * lambda * disp[c];
        } else {
          for (int b = 0; b < 3; ++b)
            for (int a = 0; a < 3; ++a)
              dfield[b][a] +=
                  q * (lambda * (a == b ? 1.0 : 0.0) +
                       lambda_over_u * disp[a] * disp[b]);
        }
      }
    }
  }
  return 0;
}

/* beta is the screening parameter. Charges narrower than beta use the
 * real-space complement; wider charges stay in the reciprocal sum. */
static int split_charge(const double rows[3][3], const double recip[3][3],
                        const int nbound[3], double volume, double beta_cap,
                        double q, double sigma, const double delta[3],
                        int deriv, double field[3], double dfield[3][3]) {
  double alpha = 1.0 / (sigma * sqrt(2.0));
  double beta = alpha < beta_cap ? alpha : beta_cap;
  int rc;
  add_split_reciprocal(recip, nbound, volume, beta, q, delta, deriv, field,
                       dfield);
  rc = add_split_real(rows, recip, alpha, beta, q, delta, deriv, field, dfield);
  return rc;
}

/* Add the induction energy and its derivatives. Positions, widths, and the
 * cell are Angstrom. grad receives dE/dR in Hartree/Bohr. ext_force receives
 * -dE/dR on the MM charges, in Hartree/Bohr. Returns 0, -1 when the cell
 * volume is not usable, or 2 when the real-space shell does not fit in the
 * Miller cap. */
int cpmdc_reference_periodic_induction(int n_atoms, const double *pos,
                                       const int *atomic_numbers, int n_ext,
                                       const double *ext_pos, const double *charge,
                                       const double *width, const double *cell,
                                       double *energy, double *grad,
                                       double *ext_force) {
  double cell_b[9];
  double rows[3][3];
  double recip[3][3];
  double det;
  double volume;
  double sigma_min;
  double gcut;
  double four_pi;
  double a_len[3];
  double a_max;
  double beta_cap;
  double pi;
  int nbound[3];
  int needs_split;
  int atom;
  int m;

  if (n_atoms <= 0 || n_ext < 0 || !pos || !atomic_numbers || !cell || !energy ||
      !grad)
    return -1;
  if (n_ext == 0)
    return 0;
  if (!ext_pos || !charge || !width || !ext_force)
    return -1;

  for (int i = 0; i < 9; ++i) {
    if (!isfinite(cell[i]))
      return -1;
    cell_b[i] = cell[i] / CPMDC_BOHR_TO_ANGSTROM;
  }
  for (int r = 0; r < 3; ++r)
    for (int c = 0; c < 3; ++c)
      rows[r][c] = cell_b[3 * r + c];

  det = rows[0][0] * (rows[1][1] * rows[2][2] - rows[1][2] * rows[2][1]) -
        rows[0][1] * (rows[1][0] * rows[2][2] - rows[1][2] * rows[2][0]) +
        rows[0][2] * (rows[1][0] * rows[2][1] - rows[1][1] * rows[2][0]);
  volume = fabs(det);
  if (!(volume > 0.0) || !isfinite(volume))
    return -1;

  {
    double bc[3], ca[3], ab[3];
    const double scale = 2.0 * acos(-1.0) / det;
    cross3(rows[1], rows[2], bc);
    cross3(rows[2], rows[0], ca);
    cross3(rows[0], rows[1], ab);
    for (int c = 0; c < 3; ++c) {
      recip[0][c] = scale * bc[c];
      recip[1][c] = scale * ca[c];
      recip[2][c] = scale * ab[c];
    }
  }

  sigma_min = INFINITY;
  for (m = 0; m < n_ext; ++m) {
    double sigma = width[m] / CPMDC_BOHR_TO_ANGSTROM;
    if (!(sigma > 0.0) || !isfinite(sigma) || !isfinite(charge[m]))
      return -1;
    for (int c = 0; c < 3; ++c)
      if (!isfinite(ext_pos[3 * m + c]))
        return -1;
    if (sigma < sigma_min)
      sigma_min = sigma;
  }

  pi = acos(-1.0);
  a_max = 0.0;
  for (int r = 0; r < 3; ++r) {
    double len = sqrt(dot3(recip[r], recip[r]));
    if (!(len > 0.0) || !isfinite(len))
      return -1;
    a_len[r] = sqrt(dot3(rows[r], rows[r]));
    if (!(a_len[r] > 0.0) || !isfinite(a_len[r]))
      return -1;
    if (a_len[r] > a_max)
      a_max = a_len[r];
  }

  /* |h_i| <= ceil(|G|_cut |a_i| / 2 pi). One bound from the shortest
   * reciprocal vector is enough only for an orthogonal cell. */
  gcut = sqrt(2.0 * CPMDC_EWALD_HALF_DECAY) / sigma_min;
  needs_split = 0;
  for (int r = 0; r < 3; ++r) {
    double ratio = gcut * a_len[r] / (2.0 * pi);
    if (!isfinite(ratio))
      return -1;
    if (ratio > (double)CPMDC_EWALD_INDEX_CAP)
      needs_split = 1;
    nbound[r] = miller_bound(ratio);
  }
  beta_cap = (double)CPMDC_EWALD_INDEX_CAP * pi / (6.0 * a_max);
  if (needs_split) {
    double gcut_split;
    if (!(beta_cap > 0.0) || !isfinite(beta_cap))
      return -1;
    gcut_split = 12.0 * beta_cap;
    for (int r = 0; r < 3; ++r) {
      double ratio = gcut_split * a_len[r] / (2.0 * pi);
      if (!isfinite(ratio))
        return -1;
      nbound[r] = miller_bound(ratio);
    }
    for (atom = 0; atom < n_atoms; ++atom) {
      for (m = 0; m < n_ext; ++m) {
        double delta[3];
        double sigma = width[m] / CPMDC_BOHR_TO_ANGSTROM;
        double alpha = 1.0 / (sigma * sqrt(2.0));
        double rcut;
        int lo[3];
        int hi[3];
        if (!(alpha > beta_cap))
          continue;
        rcut = sqrt(CPMDC_EWALD_HALF_DECAY) / beta_cap;
        for (int c = 0; c < 3; ++c) {
          if (!isfinite(pos[3 * atom + c]))
            return -1;
          delta[c] = pos[3 * atom + c] / CPMDC_BOHR_TO_ANGSTROM -
                     ext_pos[3 * m + c] / CPMDC_BOHR_TO_ANGSTROM;
        }
        if (real_window(recip, delta, rcut, lo, hi) != 0)
          return 2;
      }
    }
  }

  four_pi = 4.0 * acos(-1.0);

  for (atom = 0; atom < n_atoms; ++atom) {
    double centre[3];
    double field[3] = {0.0, 0.0, 0.0};
    double alpha;
    int z = atomic_numbers[atom];

    for (int c = 0; c < 3; ++c) {
      if (!isfinite(pos[3 * atom + c]))
        return -1;
      centre[c] = pos[3 * atom + c] / CPMDC_BOHR_TO_ANGSTROM;
    }
    alpha = 0.20 * (double)(z > 1 ? z : 1);

    for (m = 0; m < n_ext; ++m) {
      double delta[3];
      double sigma = width[m] / CPMDC_BOHR_TO_ANGSTROM;
      double q = charge[m];
      for (int c = 0; c < 3; ++c)
        delta[c] = centre[c] - ext_pos[3 * m + c] / CPMDC_BOHR_TO_ANGSTROM;
      if (needs_split) {
        int rc = split_charge(rows, recip, nbound, volume, beta_cap, q, sigma,
                              delta, 0, field, NULL);
        if (rc != 0)
          return rc;
        continue;
      }
      for (int h = -nbound[0]; h <= nbound[0]; ++h) {
        for (int k = -nbound[1]; k <= nbound[1]; ++k) {
          for (int l = -nbound[2]; l <= nbound[2]; ++l) {
            double gvec[3];
            double g2;
            double phase;
            double pref;
            double s;
            if (h == 0 && k == 0 && l == 0)
              continue;
            for (int c = 0; c < 3; ++c)
              gvec[c] = (double)h * recip[0][c] + (double)k * recip[1][c] +
                        (double)l * recip[2][c];
            g2 = dot3(gvec, gvec);
            if (!(g2 > 0.0))
              continue;
            if (0.5 * sigma_min * sigma_min * g2 > CPMDC_EWALD_HALF_DECAY)
              continue;
            phase = dot3(gvec, delta);
            pref = (four_pi * q / volume) * exp(-0.5 * sigma * sigma * g2) / g2;
            s = sin(phase);
            for (int c = 0; c < 3; ++c)
              field[c] += pref * gvec[c] * s;
          }
        }
      }
    }

    *energy += -0.5 * alpha * dot3(field, field);

    for (m = 0; m < n_ext; ++m) {
      double delta[3];
      double sigma = width[m] / CPMDC_BOHR_TO_ANGSTROM;
      double q = charge[m];
      double dfield[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
      for (int c = 0; c < 3; ++c)
        delta[c] = centre[c] - ext_pos[3 * m + c] / CPMDC_BOHR_TO_ANGSTROM;
      if (needs_split) {
        int rc = split_charge(rows, recip, nbound, volume, beta_cap, q, sigma,
                              delta, 1, NULL, dfield);
        if (rc != 0)
          return rc;
      } else {
        for (int h = -nbound[0]; h <= nbound[0]; ++h) {
          for (int k = -nbound[1]; k <= nbound[1]; ++k) {
            for (int l = -nbound[2]; l <= nbound[2]; ++l) {
              double gvec[3];
              double g2;
              double phase;
              double pref;
              double cosph;
              if (h == 0 && k == 0 && l == 0)
                continue;
              for (int c = 0; c < 3; ++c)
                gvec[c] = (double)h * recip[0][c] + (double)k * recip[1][c] +
                          (double)l * recip[2][c];
              g2 = dot3(gvec, gvec);
              if (!(g2 > 0.0))
                continue;
              if (0.5 * sigma_min * sigma_min * g2 > CPMDC_EWALD_HALF_DECAY)
                continue;
              phase = dot3(gvec, delta);
              pref = (four_pi * q / volume) * exp(-0.5 * sigma * sigma * g2) / g2;
              cosph = cos(phase);
              for (int b = 0; b < 3; ++b)
                for (int a = 0; a < 3; ++a)
                  dfield[b][a] += pref * gvec[b] * cosph * gvec[a];
            }
          }
        }
      }
      for (int a = 0; a < 3; ++a) {
        double sum = 0.0;
        for (int b = 0; b < 3; ++b)
          sum += field[b] * dfield[b][a];
        sum *= alpha;
        grad[3 * atom + a] -= sum;
        ext_force[3 * m + a] -= sum;
      }
    }
  }
  return 0;
}
