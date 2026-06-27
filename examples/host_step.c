#include "cpmdc.h"

#include <stdio.h>
#include <stdlib.h>

static unsigned char *read_file(const char *path, size_t *size_out) {
  FILE *fp = fopen(path, "rb");
  if (fp == NULL) {
    return NULL;
  }

  if (fseek(fp, 0, SEEK_END) != 0) {
    fclose(fp);
    return NULL;
  }

  long len = ftell(fp);
  if (len <= 0) {
    fclose(fp);
    return NULL;
  }
  rewind(fp);

  unsigned char *data = malloc((size_t)len);
  if (data == NULL) {
    fclose(fp);
    return NULL;
  }

  size_t got = fread(data, 1, (size_t)len, fp);
  fclose(fp);
  if (got != (size_t)len) {
    free(data);
    return NULL;
  }

  *size_out = got;
  return data;
}

int main(int argc, char **argv) {
  if (argc != 3) {
    fprintf(stderr, "usage: %s CPMDParams.bin ForceInput.bin\n", argv[0]);
    return 2;
  }

  size_t params_size = 0;
  size_t step_size = 0;
  unsigned char *params = read_file(argv[1], &params_size);
  unsigned char *step = read_file(argv[2], &step_size);
  if (params == NULL || step == NULL) {
    fprintf(stderr, "failed to read input messages\n");
    free(params);
    free(step);
    return 2;
  }

  CPMDCSession *session = cpmdc_session_create(params, params_size);
  if (session == NULL) {
    fprintf(stderr, "failed to create CPMDCSession\n");
    free(params);
    free(step);
    return 1;
  }

  size_t capacity =
      cpmdc_potential_result_size_for_force_input(step, step_size);
  if (capacity == 0) {
    fprintf(stderr, "failed to size PotentialResult\n");
    cpmdc_session_destroy(session);
    free(params);
    free(step);
    return 1;
  }

  unsigned char *potential_result = malloc(capacity);
  if (potential_result == NULL) {
    fprintf(stderr, "failed to allocate PotentialResult\n");
    cpmdc_session_destroy(session);
    free(params);
    free(step);
    return 1;
  }

  size_t potential_result_size = 0;
  CPMDCResult result = cpmdc_session_calculate_result(
      session, step, step_size, potential_result, capacity,
      &potential_result_size);
  if (!result.ok) {
    fprintf(stderr, "cpmdc_session_calculate_result failed: %s\n",
            result.message);
    free(potential_result);
    cpmdc_session_destroy(session);
    free(params);
    free(step);
    return 1;
  }

  printf("energy_h=%.12f\n", result.energy_h);
  printf("potential_result_size_bytes=%zu\n", potential_result_size);
  printf("message=%s\n", result.message[0] ? result.message : "ok");

  free(potential_result);
  cpmdc_session_destroy(session);
  free(params);
  free(step);
  cpmdc_finalize();
  return 0;
}
