/* SPDX-License-Identifier: MIT */
#include "pwscale_fft.h"

#include "cpmdc_pwscale.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#if defined(CPMDC_PWSCALE_HAVE_FFTW)
#include <fftw3.h>
#if defined(_OPENMP)
#include <omp.h>
#endif

static double fft_now(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

static void fill_grid(double *grid, ptrdiff_t nreal) {
  ptrdiff_t i;
  for (i = 0; i < nreal; ++i)
    grid[i] = sin(0.01 * (double)i) + 0.1;
}

static double max_abs_diff(const double *a, const double *b, ptrdiff_t n) {
  ptrdiff_t i;
  double peak = 0.0;
  for (i = 0; i < n; ++i) {
    double gap = fabs(a[i] - b[i]);
    if (gap > peak)
      peak = gap;
  }
  return peak;
}

static double roundtrip_peak(const double *got, const double *orig, ptrdiff_t nreal,
                             double scale) {
  ptrdiff_t i;
  double peak = 0.0;
  for (i = 0; i < nreal; ++i) {
    double gap = fabs(got[i] * scale - orig[i]);
    if (gap > peak)
      peak = gap;
  }
  return peak;
}

int cpmdc_pwscale_fft(int n, int workers, int weak_batch, double *seconds,
                      double *plan_s, double *roundtrip_max_abs,
                      double *versus_serial_max_abs, long long *alloc_bytes,
                      int *actual_threads, long long *points_executed) {
  const int nthreads = workers < 1 ? 1 : workers;
  const int batch = weak_batch ? nthreads : 1;
  const int plan_threads = weak_batch ? 1 : nthreads;
  ptrdiff_t nreal, ncomplex, i;
  double scale, t0, t1, plan0, plan1;
  double *serial_in = NULL, *serial_out = NULL;
  fftw_complex *serial_spec = NULL;
  double **in = NULL;
  double **out = NULL;
  fftw_complex **spec = NULL;
  fftw_plan *fwd = NULL;
  fftw_plan *rev = NULL;
  fftw_plan serial_fwd, serial_rev;
  int b;
  static int threads_ready = 0;

  if (n < 2 || (n % 2) != 0)
    return CPMDC_PWSCALE_ERR_ARG;
  if (!seconds || !plan_s || !roundtrip_max_abs || !versus_serial_max_abs ||
      !alloc_bytes || !actual_threads || !points_executed)
    return CPMDC_PWSCALE_ERR_ARG;
  *actual_threads = 0;
  *points_executed = 0;

  if (!threads_ready) {
    if (fftw_init_threads() == 0)
      return CPMDC_PWSCALE_ERR_NOFFTW;
    threads_ready = 1;
  }

  nreal = (ptrdiff_t)n * n * n;
  ncomplex = (ptrdiff_t)n * n * (n / 2 + 1);
  scale = 1.0 / (double)nreal;
  *alloc_bytes = 0;

  serial_in = (double *)fftw_malloc((size_t)nreal * sizeof(double));
  serial_out = (double *)fftw_malloc((size_t)nreal * sizeof(double));
  serial_spec = (fftw_complex *)fftw_malloc((size_t)ncomplex * sizeof(fftw_complex));
  if (!serial_in || !serial_out || !serial_spec)
    return CPMDC_PWSCALE_ERR_ARG;
  *alloc_bytes += (long long)nreal * 2 * (long long)sizeof(double);
  *alloc_bytes += (long long)ncomplex * (long long)sizeof(fftw_complex);
  fill_grid(serial_in, nreal);

  plan0 = fft_now();
  fftw_plan_with_nthreads(1);
  serial_fwd = fftw_plan_dft_r2c_3d(n, n, n, serial_in, serial_spec, FFTW_ESTIMATE);
  serial_rev = fftw_plan_dft_c2r_3d(n, n, n, serial_spec, serial_out, FFTW_ESTIMATE);
  if (!serial_fwd || !serial_rev)
    return CPMDC_PWSCALE_ERR_NOFFTW;

  in = (double **)calloc((size_t)batch, sizeof(double *));
  out = (double **)calloc((size_t)batch, sizeof(double *));
  spec = (fftw_complex **)calloc((size_t)batch, sizeof(fftw_complex *));
  fwd = (fftw_plan *)calloc((size_t)batch, sizeof(fftw_plan));
  rev = (fftw_plan *)calloc((size_t)batch, sizeof(fftw_plan));
  if (!in || !out || !spec || !fwd || !rev)
    return CPMDC_PWSCALE_ERR_ARG;

  fftw_plan_with_nthreads(plan_threads);
  for (b = 0; b < batch; ++b) {
    in[b] = (double *)fftw_malloc((size_t)nreal * sizeof(double));
    out[b] = (double *)fftw_malloc((size_t)nreal * sizeof(double));
    spec[b] = (fftw_complex *)fftw_malloc((size_t)ncomplex * sizeof(fftw_complex));
    if (!in[b] || !out[b] || !spec[b])
      return CPMDC_PWSCALE_ERR_ARG;
    *alloc_bytes += (long long)nreal * 2 * (long long)sizeof(double);
    *alloc_bytes += (long long)ncomplex * (long long)sizeof(fftw_complex);
    fill_grid(in[b], nreal);
    fwd[b] = fftw_plan_dft_r2c_3d(n, n, n, in[b], spec[b], FFTW_ESTIMATE);
    rev[b] = fftw_plan_dft_c2r_3d(n, n, n, spec[b], out[b], FFTW_ESTIMATE);
    if (!fwd[b] || !rev[b])
      return CPMDC_PWSCALE_ERR_NOFFTW;
  }
  plan1 = fft_now();
  *plan_s = plan1 - plan0;
  fftw_execute(serial_fwd);
  fftw_execute(serial_rev);

  if (weak_batch) {
    int observed = 1;
#if defined(_OPENMP)
    omp_set_dynamic(0);
#pragma omp parallel num_threads(batch)
    {
      if (omp_get_thread_num() == 0)
        observed = omp_get_num_threads();
    }
#endif
    *actual_threads = observed;
  }
  t0 = fft_now();
  if (weak_batch) {
#if defined(_OPENMP)
#pragma omp parallel for num_threads(batch) schedule(static)
#endif
    for (b = 0; b < batch; ++b) {
      fftw_execute(fwd[b]);
      fftw_execute(rev[b]);
    }
  } else {
    fftw_execute(fwd[0]);
    fftw_execute(rev[0]);
    *actual_threads = 0;
  }
  t1 = fft_now();
  *seconds = t1 - t0;
  *points_executed = (long long)batch * nreal;

  *roundtrip_max_abs = roundtrip_peak(out[0], in[0], nreal, scale);
  *versus_serial_max_abs = max_abs_diff(out[0], serial_out, nreal);
  for (i = 1; i < batch; ++i) {
    double rt = roundtrip_peak(out[i], in[i], nreal, scale);
    double vs = max_abs_diff(out[i], serial_out, nreal);
    if (rt > *roundtrip_max_abs)
      *roundtrip_max_abs = rt;
    if (vs > *versus_serial_max_abs)
      *versus_serial_max_abs = vs;
  }

  for (b = 0; b < batch; ++b) {
    if (fwd && fwd[b])
      fftw_destroy_plan(fwd[b]);
    if (rev && rev[b])
      fftw_destroy_plan(rev[b]);
    fftw_free(in[b]);
    fftw_free(out[b]);
    fftw_free(spec[b]);
  }
  fftw_destroy_plan(serial_fwd);
  fftw_destroy_plan(serial_rev);
  fftw_free(serial_in);
  fftw_free(serial_out);
  fftw_free(serial_spec);
  free(in);
  free(out);
  free(spec);
  free(fwd);
  free(rev);
  return CPMDC_PWSCALE_OK;
}

#else

int cpmdc_pwscale_fft(int n, int workers, int weak_batch, double *seconds,
                      double *plan_s, double *roundtrip_max_abs,
                      double *versus_serial_max_abs, long long *alloc_bytes,
                      int *actual_threads, long long *points_executed) {
  (void)n;
  (void)workers;
  (void)weak_batch;
  (void)seconds;
  (void)plan_s;
  (void)roundtrip_max_abs;
  (void)versus_serial_max_abs;
  (void)alloc_bytes;
  if (actual_threads)
    *actual_threads = 0;
  if (points_executed)
    *points_executed = 0;
  return CPMDC_PWSCALE_ERR_NOFFTW;
}

#endif
