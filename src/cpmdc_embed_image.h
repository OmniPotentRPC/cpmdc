/* Caller-owned embed result. One calculator holds one of these.
 * Module variables do not. */
#pragma once

#include "cpmdc.h"

/* Hirshfeld charges sit beside the public snapshot. Callers of
 * cpmdc_last_property_snapshot() still receive the historical
 * CPMDCPropertySnapshot and nothing past it. */
#define CPMDC_EMBED_CHARGE_CAP 1365

typedef struct CPMDCEmbedCharges {
  int valid;
  size_t count;
  double values[CPMDC_EMBED_CHARGE_CAP];
} CPMDCEmbedCharges;

typedef struct CPMDCEmbedImage {
  CPMDCEnergyComponents energy;
  CPMDCChargeIntegrals charge;
  CPMDCMultiStateEnergies multi;
  CPMDCMDTrajectoryRow md;
  CPMDCPropertySnapshot prop;
  CPMDCEmbedCharges charges;
  CPMDCStressTensor stress;
  double warm_cell[9];
  int warm_cell_set;
  int warm_has_cell;
  int cfg_warm_steps;
  int cfg_set;
  char functional[64];
  double cutoff_ry;
  int cfg_charge;
  int multiplicity;
  /* Owned method deck. NULL when unset. input_deck_len is the byte count
   * and does not include the trailing NUL; the allocation is one past it. */
  char *input_deck;
  int input_deck_len;
  char cpmd_root[1024];
  /* permanentDir, else scratchDir. Empty: host working directory. */
  char output_dir[1024];
} CPMDCEmbedImage;
