/* Reference-build periodic electrostatic embedding.
 * Conducting-boundary Ewald field of the same Gaussian charges the isolated
 * kernel uses. Not a public symbol: libcpmdc hides it. */
#include "cpmdc_params.h"

#include <math.h>

/* Keep a reciprocal vector when sigma_min^2 |G|^2 / 2 is at most this. */
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

/* Add the induction energy and its derivatives. Positions, widths, and the
 * cell are Angstrom. grad receives dE/dR in Hartree/Bohr. ext_force receives
 * -dE/dR on the MM charges, in Hartree/Bohr. Returns 0, or -1 when the cell
 * volume is not usable. */
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
  double bmin;
  double ratio;
  double four_pi;
  int nmax;
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

  gcut = sqrt(2.0 * CPMDC_EWALD_HALF_DECAY) / sigma_min;
  bmin = INFINITY;
  for (int r = 0; r < 3; ++r) {
    double len = sqrt(dot3(recip[r], recip[r]));
    if (!(len > 0.0) || !isfinite(len))
      return -1;
    if (len < bmin)
      bmin = len;
  }
  ratio = gcut / bmin;
  nmax = (int)ceil(ratio) + 1;
  if (nmax < 1)
    nmax = 1;
  if (nmax > CPMDC_EWALD_INDEX_CAP)
    nmax = CPMDC_EWALD_INDEX_CAP;

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
      for (int h = -nmax; h <= nmax; ++h) {
        for (int k = -nmax; k <= nmax; ++k) {
          for (int l = -nmax; l <= nmax; ++l) {
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
      for (int h = -nmax; h <= nmax; ++h) {
        for (int k = -nmax; k <= nmax; ++k) {
          for (int l = -nmax; l <= nmax; ++l) {
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
