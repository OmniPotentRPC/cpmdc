/* Caller-owned embed result. One calculator holds one of these.
 * Module variables do not. */
#pragma once

#include "cpmdc.h"

typedef struct CPMDCEmbedImage {
  CPMDCEnergyComponents energy;
  CPMDCChargeIntegrals charge;
  CPMDCMultiStateEnergies multi;
  CPMDCMDTrajectoryRow md;
  CPMDCPropertySnapshot prop;
  CPMDCStressTensor stress;
  double warm_cell[9];
  int warm_cell_set;
  int warm_has_cell;
  int cfg_warm_steps;
} CPMDCEmbedImage;
