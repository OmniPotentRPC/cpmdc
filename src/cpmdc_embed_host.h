/* Host-side helpers for the OpenCPMD embed: pseudopotential directory,
 * species-block order, and the directory CPMD writes into. */
#pragma once

#include <stddef.h>

/* 0 and a directory path when CPMDC_PSEUDO_DIR, or else
 * CPMD_PP_LIBRARY_PATH, names an existing directory. -1 and an error
 * message otherwise. The message does not describe an SCF failure. */
int cpmdc_pseudopotential_directory(char *buf, size_t cap);

/* map_out[slot] is the 0-based ForceInput index that fills species block
 * `slot` (species_z / species_count in CPMD order). 0 on success, -1 when
 * the atomic numbers cannot fill those blocks. */
int cpmdc_species_order_map(int n_atoms, const int *atomic_numbers,
                            int n_species, const int *species_z,
                            const int *species_count, int *map_out);

/* species_grad is n_atoms * 3 doubles in species-block order. */
void cpmdc_scatter_species_gradient(int n_atoms, const int *map,
                                    const double *species_grad, double *grad);

/* Record an embed failure for the C result path. msg is a C string. */
void cpmdc_note_embed_failure(const char *msg);

int cpmdc_prepare_pp_cwd(const char *pseudo_dir);

/* After pseudopotentials are read: chdir to output_dir when it is non-empty,
 * otherwise back to the host directory saved by cpmdc_prepare_pp_cwd.
 * 0 on success, -1 when output_dir is set and is not a directory. */
int cpmdc_enter_output_cwd(const char *output_dir);

int cpmdc_restore_host_cwd(void);
