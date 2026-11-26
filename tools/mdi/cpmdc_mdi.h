/**
 * MDI engine state for LAMMPS fix mdi/qm, without libmdi.
 *
 * The executable in cpmdc_mdi_main.c moves bytes.
 * This object interprets them and calls cpmdc_session_calculate_au().
 * Wire units are Bohr and Hartree.
 * The nine stress components are row-major Hartree/Bohr^3, the tensor
 * cpmdc_session_last_stress() returns.
 */
#ifndef CPMDC_MDI_H
#define CPMDC_MDI_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

enum {
  CPMDC_MDI_NONE = 0,
  CPMDC_MDI_RECV_INT = 1,
  CPMDC_MDI_RECV_DOUBLE = 2,
  CPMDC_MDI_SEND_INT = 3,
  CPMDC_MDI_SEND_DOUBLE = 4,
  CPMDC_MDI_SEND_NAME = 5
};

typedef struct cpmdc_mdi cpmdc_mdi;

/* Copies the CPMDParams message. Returns NULL when params is missing or
 * the copy cannot be allocated. */
cpmdc_mdi *cpmdc_mdi_create(const void *params, size_t nbytes);
void cpmdc_mdi_destroy(cpmdc_mdi *eng);

/* @DEFAULT commands, NULL-terminated. Includes EXIT. */
const char *const *cpmdc_mdi_commands(void);

/* direction and count describe the payload of one command.
 * count is the number of ints, doubles, or, for a name, the bytes of the
 * C string including the NUL. Returns -1 when the command is unknown or
 * the atom count is not known yet. The reason is cpmdc_mdi_error(). */
int cpmdc_mdi_payload(cpmdc_mdi *eng, const char *command, int *direction,
                      int *count);

/* Applies one command. incoming is the driver payload for a '>' command.
 * outgoing receives a '<' command. Returns 0, or -1 with cpmdc_mdi_error(). */
int cpmdc_mdi_command(cpmdc_mdi *eng, const char *command, const void *incoming,
                      void *outgoing);

const char *cpmdc_mdi_error(const cpmdc_mdi *eng);
int cpmdc_mdi_evaluations(const cpmdc_mdi *eng);
int cpmdc_mdi_exited(const cpmdc_mdi *eng);

#ifdef __cplusplus
}
#endif

#endif
