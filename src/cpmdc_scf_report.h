/* The message an embed SCF leaves when it stops short of the orbital
 * threshold. The Fortran embed reports the numbers; this file words them. */
#pragma once

#include <stddef.h>

/* The minimizer of one SCF pass, as the cntl flags of OpenCPMD name it. */
enum {
  CPMDC_SCF_MINIMIZER_UNKNOWN = 0,
  CPMDC_SCF_MINIMIZER_ODIIS = 1,
  CPMDC_SCF_MINIMIZER_PCG = 2,
  CPMDC_SCF_MINIMIZER_PCG_MINIMIZE = 3,
  CPMDC_SCF_MINIMIZER_STEEPEST_DESCENT = 4
};

/* The deck keyword of a minimizer code. */
const char *cpmdc_scf_minimizer_name(int minimizer);

/* Write the message of an SCF that did not converge to buf. The last pass
 * used `minimizer` and stopped after `iterations` of `maxiter` with the
 * largest orbital gradient `gemax` against the threshold `tolog`. When the
 * embed continued after a first pass, `first_minimizer` and
 * `first_iterations` describe it; first_minimizer 0 means there was one
 * pass. Returns what snprintf returns. */
int cpmdc_format_scf_unconverged(char *buf, size_t cap, int minimizer,
                                 int iterations, int maxiter, double gemax,
                                 double tolog, int first_minimizer,
                                 int first_iterations);

/* Format the message and record it as the embed failure of this call. */
void cpmdc_note_scf_unconverged(int minimizer, int iterations, int maxiter,
                                double gemax, double tolog,
                                int first_minimizer, int first_iterations);
