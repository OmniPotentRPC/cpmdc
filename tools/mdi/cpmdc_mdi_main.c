/**
 * cpmdc-mdi: stand-alone MDI engine for LAMMPS fix mdi/qm.
 *
 * Launch beside LAMMPS, one engine process:
 *
 *   cpmdc-mdi params.bin -mdi
 *       "-role ENGINE -name cpmdc -method TCP -hostname localhost -port 8021"
 *   lmp -mdi
 *       "-role DRIVER -name LAMMPS -method TCP -hostname localhost -port 8021"
 *       -in in.mdi
 *
 * MDI_Init removes the -mdi argument. The remaining argument is the
 * CPMDParams message. The command meanings are in cpmdc_mdi.c.
 *
 * Exit status: 0 after EXIT, 1 on an MDI or evaluation error, 2 on usage.
 */
#include "cpmdc_mdi.h"

#include <mdi.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define PROG "cpmdc-mdi"

static unsigned char *read_file(const char *path, size_t *size) {
  FILE *fp = fopen(path, "rb");
  unsigned char *buf = NULL;
  long n = 0;
  if (fp && fseek(fp, 0, SEEK_END) == 0 && (n = ftell(fp)) > 0 &&
      fseek(fp, 0, SEEK_SET) == 0 && (buf = malloc((size_t)n)) != NULL &&
      fread(buf, 1, (size_t)n, fp) != (size_t)n) {
    free(buf);
    buf = NULL;
  }
  if (fp)
    fclose(fp);
  *size = buf ? (size_t)n : 0;
  return buf;
}

static int exchange(cpmdc_mdi *eng, MDI_Comm comm, const char *command) {
  int direction = 0, count = 0, rc;
  void *buf = NULL;
  int mdi_type = 0;
  size_t width = 0;

  if (strcmp(command, "EXIT") == 0) {
    if (cpmdc_mdi_command(eng, command, NULL, NULL) != 0) {
      fprintf(stderr, "%s: %s\n", PROG, cpmdc_mdi_error(eng));
      return -1;
    }
    return 1;
  }
  if (cpmdc_mdi_payload(eng, command, &direction, &count) != 0) {
    fprintf(stderr, "%s: %s\n", PROG, cpmdc_mdi_error(eng));
    return -1;
  }
  if (direction == CPMDC_MDI_RECV_INT || direction == CPMDC_MDI_SEND_INT) {
    width = sizeof(int);
    mdi_type = MDI_INT;
  } else if (direction == CPMDC_MDI_RECV_DOUBLE ||
             direction == CPMDC_MDI_SEND_DOUBLE) {
    width = sizeof(double);
    mdi_type = MDI_DOUBLE;
  } else if (direction == CPMDC_MDI_SEND_NAME) {
    buf = calloc((size_t)MDI_NAME_LENGTH, 1);
    if (!buf)
      return -1;
    rc = cpmdc_mdi_command(eng, command, NULL, buf);
    if (rc == 0)
      rc = MDI_Send(buf, MDI_NAME_LENGTH, MDI_CHAR, comm);
    free(buf);
    if (rc != 0) {
      fprintf(stderr, "%s: %s\n", PROG,
              cpmdc_mdi_error(eng)[0] ? cpmdc_mdi_error(eng) : "MDI send failed");
      return -1;
    }
    return 0;
  }

  if (width == 0)
    return cpmdc_mdi_command(eng, command, NULL, NULL) == 0 ? 0 : -1;

  buf = calloc((size_t)count, width);
  if (!buf)
    return -1;
  if (direction == CPMDC_MDI_RECV_INT || direction == CPMDC_MDI_RECV_DOUBLE) {
    rc = MDI_Recv(buf, count, mdi_type, comm);
    if (rc == 0)
      rc = cpmdc_mdi_command(eng, command, buf, NULL);
  } else {
    rc = cpmdc_mdi_command(eng, command, NULL, buf);
    if (rc == 0)
      rc = MDI_Send(buf, count, mdi_type, comm);
  }
  free(buf);
  if (rc != 0) {
    fprintf(stderr, "%s: %s\n", PROG,
            cpmdc_mdi_error(eng)[0] ? cpmdc_mdi_error(eng) : "MDI exchange failed");
    return -1;
  }
  return 0;
}

int main(int argc, char **argv) {
  size_t nbytes = 0;
  unsigned char *params;
  cpmdc_mdi *eng;
  const char *const *commands;
  MDI_Comm comm;
  char *command;
  int step;

  if (MDI_Init(&argc, &argv) != 0) {
    fprintf(stderr, "%s: MDI_Init failed\n", PROG);
    return 1;
  }
  if (argc != 2) {
    fprintf(stderr, "usage: %s PARAMS [-mdi \"MDI options\"]\n", PROG);
    return 2;
  }
  params = read_file(argv[1], &nbytes);
  eng = params ? cpmdc_mdi_create(params, nbytes) : NULL;
  free(params);
  if (!eng) {
    fprintf(stderr, "%s: cannot read %s\n", PROG, argv[1]);
    return 1;
  }
  if (MDI_Register_node("@DEFAULT") != 0) {
    fprintf(stderr, "%s: MDI_Register_node failed\n", PROG);
    cpmdc_mdi_destroy(eng);
    return 1;
  }
  commands = cpmdc_mdi_commands();
  for (int i = 0; commands[i]; i++) {
    if (MDI_Register_command("@DEFAULT", commands[i]) != 0) {
      fprintf(stderr, "%s: MDI_Register_command %s failed\n", PROG, commands[i]);
      cpmdc_mdi_destroy(eng);
      return 1;
    }
  }
  comm = MDI_COMM_NULL;
  if (MDI_Accept_communicator(&comm) != 0 || comm == MDI_COMM_NULL) {
    fprintf(stderr, "%s: MDI_Accept_communicator failed\n", PROG);
    cpmdc_mdi_destroy(eng);
    return 1;
  }
  command = calloc((size_t)MDI_COMMAND_LENGTH, 1);
  if (!command) {
    cpmdc_mdi_destroy(eng);
    return 1;
  }
  for (;;) {
    if (MDI_Recv_command(command, comm) != 0) {
      fprintf(stderr, "%s: MDI_Recv_command failed\n", PROG);
      free(command);
      cpmdc_mdi_destroy(eng);
      return 1;
    }
    step = exchange(eng, comm, command);
    if (step == 1) {
      free(command);
      cpmdc_mdi_destroy(eng);
      return 0;
    }
    if (step != 0) {
      free(command);
      cpmdc_mdi_destroy(eng);
      return 1;
    }
  }
}
