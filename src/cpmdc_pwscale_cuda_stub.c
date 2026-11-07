/* SPDX-License-Identifier: MIT */
#include "cpmdc_pwscale.h"

#include <string.h>

int cpmdc_pwscale_cuda_compiled(void) { return 0; }

int cpmdc_pwscale_cuda_device(char *name, int namelen) {
  if (name && namelen > 0)
    name[0] = '\0';
  return 0;
}

int cpmdc_pwscale_eval_cuda(const cpmdc_pwscale_problem *prob,
                            const double *positions_bohr,
                            const int *atomic_numbers, double *forces_h_bohr,
                            cpmdc_pwscale_stats *stats) {
  (void)prob;
  (void)positions_bohr;
  (void)atomic_numbers;
  (void)forces_h_bohr;
  if (stats)
    memset(stats, 0, sizeof(*stats));
  return CPMDC_PWSCALE_ERR_NOCUDA;
}
