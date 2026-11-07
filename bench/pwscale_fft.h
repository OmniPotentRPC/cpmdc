/* SPDX-License-Identifier: MIT */
#ifndef CPMDC_PWSCALE_FFT_H
#define CPMDC_PWSCALE_FFT_H

/* FFTW round trip used only by the scaling harness.
 * plan_s is the untimed planning cost. seconds is execution.
 * A build without FFTW returns CPMDC_PWSCALE_ERR_NOFFTW.
 */
/* actual_threads is the OpenMP team inside a weak batch, or 0 when the
 * strong path's FFTW team is not visible. points_executed is the number of
 * real points this call transformed.
 */
int cpmdc_pwscale_fft(int n, int workers, int weak_batch, double *seconds,
                      double *plan_s, double *roundtrip_max_abs,
                      double *versus_serial_max_abs, long long *alloc_bytes,
                      int *actual_threads, long long *points_executed);

#endif
