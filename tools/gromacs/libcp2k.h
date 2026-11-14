/* The libcp2k C API that GROMACS' QM/MM module calls, served by cpmdc.
 *
 * GROMACS 2022 and later (-DGMX_CP2K=ON) include "libcp2k.h" and link the
 * library that pkg-config names libcp2k. libcpmdc_gromacs implements these
 * declarations, with CP2K's signatures (cp2k/src/start/libcp2k.h), on top
 * of cpmdc_session_calculate_au(). Functions of CP2K's header that GROMACS
 * does not call are not provided.
 *
 * Units are CP2K's: positions and cell in Bohr, energy in Hartree, forces
 * in Hartree/Bohr. The cell is nine doubles: CP2K's layout, one lattice
 * vector after the other (a_x a_y a_z b_x ...), as GROMACS 2026 sends it, or
 * its transpose (a_x b_x c_x a_y ...), as GROMACS 2024 and 2025 send it.
 * The library tells them apart by the zeros of GROMACS' lower triangular box.
 */
#ifndef CPMDC_GROMACS_LIBCP2K_H
#define CPMDC_GROMACS_LIBCP2K_H

#ifdef __cplusplus
extern "C" {
#endif

typedef int force_env_t;

void cp2k_get_version(char *version_str, int str_length);
void cp2k_init(void);
void cp2k_init_without_mpi(void);
void cp2k_finalize(void);
void cp2k_finalize_without_mpi(void);
void cp2k_create_force_env(force_env_t *new_force_env, const char *input_file_path,
                           const char *output_file_path);
void cp2k_create_force_env_comm(force_env_t *new_force_env, const char *input_file_path,
                                const char *output_file_path, int mpi_comm);
void cp2k_destroy_force_env(force_env_t force_env);
void cp2k_set_positions(force_env_t force_env, const double *new_pos, int n_el);
void cp2k_set_cell(force_env_t force_env, const double *new_cell);
void cp2k_get_natom(force_env_t force_env, int *natom);
void cp2k_get_positions(force_env_t force_env, double *pos, int n_el);
void cp2k_get_forces(force_env_t force_env, double *force, int n_el);
/* As in CP2K (2026): the configurational stress after cp2k_calc_energy_force,
 * nine doubles, column-major, Hartree/Bohr^3, pressure-positive, so that
 * stress times the cell volume is the virial -dE/d(strain). available is 1
 * when it was computed and 0, with a zero tensor, when it was not. */
void cp2k_get_stress_tensor(force_env_t force_env, double *stress_tensor, int *available);
void cp2k_get_potential_energy(force_env_t force_env, double *e_pot);
void cp2k_calc_energy_force(force_env_t force_env);
void cp2k_calc_energy(force_env_t force_env);

#ifdef __cplusplus
}
#endif

#endif
