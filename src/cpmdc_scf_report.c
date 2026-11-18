#include "cpmdc_scf_report.h"

#include <stdio.h>

const char *cpmdc_scf_minimizer_name(int minimizer) {
  switch (minimizer) {
  case CPMDC_SCF_MINIMIZER_ODIIS:
    return "ODIIS";
  case CPMDC_SCF_MINIMIZER_PCG:
    return "PCG";
  case CPMDC_SCF_MINIMIZER_PCG_MINIMIZE:
    return "PCG MINIMIZE";
  case CPMDC_SCF_MINIMIZER_STEEPEST_DESCENT:
    return "STEEPEST DESCENT";
  default:
    return "an unknown minimizer";
  }
}

int cpmdc_format_scf_unconverged(char *buf, size_t cap, int minimizer,
                                 int iterations, int maxiter, double gemax,
                                 double tolog, int first_minimizer,
                                 int first_iterations) {
  if (!buf || cap == 0)
    return 0;
  if (first_minimizer == CPMDC_SCF_MINIMIZER_UNKNOWN)
    return snprintf(buf, cap,
                    "CPMD SCF did not converge: %s stopped at iteration %d of "
                    "MAXITER %d with gemax %.6E above CONVERGENCE ORBITALS %.6E",
                    cpmdc_scf_minimizer_name(minimizer), iterations, maxiter,
                    gemax, tolog);
  return snprintf(buf, cap,
                  "CPMD SCF did not converge: %s stopped at iteration %d of "
                  "MAXITER %d with gemax %.6E above CONVERGENCE ORBITALS %.6E, "
                  "after %s stopped at iteration %d of MAXITER %d",
                  cpmdc_scf_minimizer_name(minimizer), iterations, maxiter, gemax,
                  tolog, cpmdc_scf_minimizer_name(first_minimizer),
                  first_iterations, maxiter);
}
