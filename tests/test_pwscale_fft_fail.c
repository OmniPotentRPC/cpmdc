/* SPDX-License-Identifier: MIT */
/* Setup failures must return an error after releasing buffers already taken. */
#define _POSIX_C_SOURCE 200809L
#include "pwscale_fft.h"

#include "cpmdc_pwscale.h"

#include <stdio.h>
#include <stdlib.h>

static int fail_once(const char *fault, int want) {
  double seconds = 0, plan_s = 0, roundtrip = 0, versus = 0;
  long long alloc_bytes = 0, points = 0;
  int threads = 0, rc;
  setenv("CPMDC_PWSCALE_FFT_FAULT", fault, 1);
  rc = cpmdc_pwscale_fft(8, 2, 0, &seconds, &plan_s, &roundtrip, &versus,
                         &alloc_bytes, &threads, &points);
  if (rc != want) {
    fprintf(stderr, "fault %s rc %d want %d\n", fault, rc, want);
    return 1;
  }
  return 0;
}

int main(void) {
  double seconds = 0, plan_s = 0, roundtrip = 0, versus = 0;
  long long alloc_bytes = 0, points = 0;
  int threads = 0, rc;
  if (fail_once("serial-plan", CPMDC_PWSCALE_ERR_NOFFTW) ||
      fail_once("batch-buf", CPMDC_PWSCALE_ERR_ARG) ||
      fail_once("batch-plan", CPMDC_PWSCALE_ERR_NOFFTW))
    return 1;
  unsetenv("CPMDC_PWSCALE_FFT_FAULT");
  rc = cpmdc_pwscale_fft(8, 1, 0, &seconds, &plan_s, &roundtrip, &versus,
                         &alloc_bytes, &threads, &points);
  if (rc != CPMDC_PWSCALE_OK || points != 512) {
    fprintf(stderr, "follow-up fft rc %d points %lld\n", rc, points);
    return 1;
  }
  printf("pwscale fft fail cleanup ok\n");
  return 0;
}
