#pragma once

#include "Potentials.capnp.h"

#include <stddef.h>

#define CPMDC_BLOCKS 16384
#define CPMDC_BOHR_TO_ANGSTROM 0.529177210903

int cpmdc_params_root(const void *params_capnp, size_t params_capnp_size_bytes,
                      struct capn *arena, CPMDParams_ptr *params);

void cpmdc_params_release(struct capn *arena);

const char *cpmdc_params_text_or(capn_text text, const char *fallback);

int cpmdc_params_effective_config(CPMDParams_ptr params, char *functional,
                                  size_t functional_size, double *cutoff_ry,
                                  int *charge, int *multiplicity);

/** Render a full CPMD input deck (&SECTION ... &END) into dst. */
int cpmdc_params_render_input_deck(CPMDParams_ptr params, char *dst,
                                   size_t dst_size);

/**
 * Render deck with &ATOMS coordinates from Cap'n Proto geometry (Angstrom).
 * PP paths come from CPMDParams atoms.pseudopotentials; groups by Z via element
 * symbol (H=1, O=8, …) matching ForceInput atomic_numbers.
 */
int cpmdc_params_render_deck_with_geometry(
    CPMDParams_ptr params, int n_atoms, const double *positions_ang,
    const int *atomic_numbers, const double *cell_ang, int has_cell, char *dst,
    size_t dst_size);

int cpmdc_force_input_root(const void *force_input_capnp,
                           size_t force_input_capnp_size_bytes,
                           struct capn *arena, ForceInput_ptr *force_input);

int cpmdc_force_input_atom_count(ForceInput_ptr force_input, size_t *n_atoms,
                                 int *has_cell);

int cpmdc_force_input_copy_geometry(ForceInput_ptr force_input,
                                    double *positions_ang, int *atomic_numbers,
                                    size_t atom_capacity, double *cell_ang,
                                    int *has_cell);

int cpmdc_force_input_result_factors(ForceInput_ptr force_input,
                                     double *energy_factor,
                                     double *force_factor);

size_t cpmdc_potential_result_flat_size(size_t force_count);

int cpmdc_potential_result_write(double energy, const double *forces,
                                 size_t force_count,
                                 void *potential_result_capnp,
                                 size_t potential_result_capacity_bytes,
                                 size_t *potential_result_size_bytes);
