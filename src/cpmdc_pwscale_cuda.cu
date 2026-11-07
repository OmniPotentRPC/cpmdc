/* SPDX-License-Identifier: MIT */
/* Optional CUDA backend. Meson builds this only with -Dwith_pwscale_cuda=true.
 * Total time includes the copies both ways and the device synchronize.
 * The kernel matches the CPU structure factor. Agreement uses the stated
 * tolerance. A machine with no device returns an error and no timing.
 */
#include "cpmdc_pwscale.h"

#include <cuda_runtime.h>

#include <math.h>
#include <stdio.h>
#include <string.h>

static const double PW_PI = 3.14159265358979323846;

static int grid_n_host(int npw) {
  int n = 1;
  while ((long long)n * n * n < (long long)npw)
    n += 2;
  return n;
}

__global__ void pwscale_forces(int npw, int natoms, int ngrid, double cell,
                               double alpha, double sigma, const double *pos,
                               const int *atomic_numbers, double *energy,
                               double *forces) {
  int k = (int)(blockIdx.x * blockDim.x + threadIdx.x);
  int half, x, y, z, ix, iy, iz, a;
  double scale, gx, gy, gz, g2, weight;
  if (k >= npw)
    return;
  half = ngrid / 2;
  x = k % ngrid;
  y = (k / ngrid) % ngrid;
  z = k / (ngrid * ngrid);
  ix = x - half;
  iy = y - half;
  iz = z - half;
  scale = (2.0 * PW_PI) / cell;
  gx = scale * (double)ix;
  gy = scale * (double)iy;
  gz = scale * (double)iz;
  g2 = (gx * gx + gy * gy) + gz * gz;
  if (g2 == 0.0)
    return;
  weight = ((4.0 * PW_PI) / g2) * exp(-g2 / (4.0 * alpha)) *
            exp(-g2 / (2.0 * sigma * sigma));
  for (a = 0; a < natoms; ++a) {
    double phase = (gx * pos[3 * a] + gy * pos[3 * a + 1]) + gz * pos[3 * a + 2];
    double c = cos(phase);
    double s = sin(phase);
    double wz = weight * (double)atomic_numbers[a];
    atomicAdd(energy, wz * c);
    atomicAdd(forces + 3 * a, wz * s * gx);
    atomicAdd(forces + 3 * a + 1, wz * s * gy);
    atomicAdd(forces + 3 * a + 2, wz * s * gz);
  }
}

extern "C" int cpmdc_pwscale_cuda_compiled(void) { return 1; }

extern "C" int cpmdc_pwscale_cuda_device(char *name, int namelen) {
  int count = 0;
  cudaDeviceProp prop;
  if (name && namelen > 0)
    name[0] = '\0';
  if (cudaGetDeviceCount(&count) != cudaSuccess || count < 1)
    return 0;
  if (cudaGetDeviceProperties(&prop, 0) != cudaSuccess)
    return 0;
  if (name && namelen > 0)
    snprintf(name, (size_t)namelen, "%s", prop.name);
  return 1;
}

extern "C" int cpmdc_pwscale_eval_cuda(const cpmdc_pwscale_problem *prob,
                                       const double *positions_bohr,
                                       const int *atomic_numbers,
                                       double *forces_h_bohr,
                                       cpmdc_pwscale_stats *stats) {
  int n3, ngrid, threads, blocks, k, nonzero;
  size_t pos_bytes, z_bytes, force_bytes;
  double *d_pos = NULL, *d_forces = NULL, *d_energy = NULL;
  int *d_z = NULL;
  double energy = 0.0;
  cudaEvent_t start, h2d_done, kernel_done, stop;
  float ms_total = 0.f, ms_h2d = 0.f, ms_d2h = 0.f;
  cudaError_t err;

  if (stats)
    memset(stats, 0, sizeof(*stats));
  if (!prob || !positions_bohr || !atomic_numbers || !forces_h_bohr)
    return CPMDC_PWSCALE_ERR_ARG;
  if (prob->natoms < 1 || prob->npw < 1 || !(prob->cell_bohr > 0.0) ||
      !(prob->alpha > 0.0) || !(prob->sigma > 0.0))
    return CPMDC_PWSCALE_ERR_ARG;
  if (!cpmdc_pwscale_cuda_device(NULL, 0))
    return CPMDC_PWSCALE_ERR_NOCUDA;

  n3 = prob->natoms * 3;
  ngrid = grid_n_host(prob->npw);
  pos_bytes = (size_t)n3 * sizeof(double);
  z_bytes = (size_t)prob->natoms * sizeof(int);
  force_bytes = pos_bytes;

  err = cudaMalloc(&d_pos, pos_bytes);
  if (err == cudaSuccess)
    err = cudaMalloc(&d_z, z_bytes);
  if (err == cudaSuccess)
    err = cudaMalloc(&d_forces, force_bytes);
  if (err == cudaSuccess)
    err = cudaMalloc(&d_energy, sizeof(double));
  if (err != cudaSuccess) {
    cudaFree(d_pos);
    cudaFree(d_z);
    cudaFree(d_forces);
    cudaFree(d_energy);
    return CPMDC_PWSCALE_ERR_CUDA;
  }
  cudaMemset(d_forces, 0, force_bytes);
  cudaMemset(d_energy, 0, sizeof(double));
  cudaEventCreate(&start);
  cudaEventCreate(&h2d_done);
  cudaEventCreate(&kernel_done);
  cudaEventCreate(&stop);

  cudaEventRecord(start);
  cudaMemcpy(d_pos, positions_bohr, pos_bytes, cudaMemcpyHostToDevice);
  cudaMemcpy(d_z, atomic_numbers, z_bytes, cudaMemcpyHostToDevice);
  cudaEventRecord(h2d_done);

  threads = 128;
  blocks = (prob->npw + threads - 1) / threads;
  pwscale_forces<<<blocks, threads>>>(prob->npw, prob->natoms, ngrid,
                                      prob->cell_bohr, prob->alpha, prob->sigma,
                                      d_pos, d_z, d_energy, d_forces);
  cudaEventRecord(kernel_done);
  cudaMemcpy(&energy, d_energy, sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(forces_h_bohr, d_forces, force_bytes, cudaMemcpyDeviceToHost);
  cudaEventRecord(stop);
  cudaEventSynchronize(stop);
  cudaDeviceSynchronize();

  cudaEventElapsedTime(&ms_total, start, stop);
  cudaEventElapsedTime(&ms_h2d, start, h2d_done);
  cudaEventElapsedTime(&ms_d2h, kernel_done, stop);

  cudaEventDestroy(start);
  cudaEventDestroy(h2d_done);
  cudaEventDestroy(kernel_done);
  cudaEventDestroy(stop);
  cudaFree(d_pos);
  cudaFree(d_z);
  cudaFree(d_forces);
  cudaFree(d_energy);

  err = cudaGetLastError();
  if (err != cudaSuccess)
    return CPMDC_PWSCALE_ERR_CUDA;
  if (stats) {
    stats->energy_h = energy;
    /* Event span covers both copies, the kernel, and the sync. */
    stats->seconds = (double)ms_total / 1e3;
    stats->h2d_seconds = (double)ms_h2d / 1e3;
    stats->d2h_seconds = (double)ms_d2h / 1e3;
    stats->g_visited = prob->npw;
    nonzero = 0;
    for (k = 0; k < prob->npw; ++k) {
      int half = ngrid / 2;
      int x = k % ngrid;
      int y = (k / ngrid) % ngrid;
      int z = k / (ngrid * ngrid);
      if ((x - half) != 0 || (y - half) != 0 || (z - half) != 0)
        nonzero += 1;
    }
    stats->atom_visits = (long long)nonzero * (long long)prob->natoms;
    stats->alloc_bytes = (long long)(pos_bytes + z_bytes + force_bytes + sizeof(double));
    stats->workers = 1;
  }
  return CPMDC_PWSCALE_OK;
}
