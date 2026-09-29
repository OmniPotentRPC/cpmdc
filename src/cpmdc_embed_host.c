#include "cpmdc_embed_host.h"

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

static int copy_trimmed_dir(const char *src, char *dst, size_t cap) {
  size_t n;
  if (!src || !dst || cap < 2)
    return -1;
  n = strnlen(src, cap);
  if (n == 0 || n >= cap)
    return -1;
  memcpy(dst, src, n);
  dst[n] = '\0';
  while (n > 1 && (dst[n - 1] == '/' || dst[n - 1] == ' '))
    dst[--n] = '\0';
  return n > 0 ? 0 : -1;
}

static int is_directory(const char *path) {
  struct stat st;
  return path && path[0] && stat(path, &st) == 0 && S_ISDIR(st.st_mode);
}

int cpmdc_pseudopotential_directory(char *buf, size_t cap) {
  const char *dir;
  char trimmed[1024];
  if (!buf || cap < 8)
    return -1;
  buf[0] = '\0';
  dir = getenv("CPMDC_PSEUDO_DIR");
  if (!dir || !dir[0])
    dir = getenv("CPMD_PP_LIBRARY_PATH");
  if (!dir || !dir[0]) {
    snprintf(buf, cap,
             "pseudopotential directory is not set: set CPMDC_PSEUDO_DIR or "
             "CPMD_PP_LIBRARY_PATH");
    return -1;
  }
  if (copy_trimmed_dir(dir, trimmed, sizeof(trimmed)) != 0) {
    snprintf(buf, cap, "pseudopotential directory path is unusable");
    return -1;
  }
  if (!is_directory(trimmed)) {
    snprintf(buf, cap,
             "pseudopotential directory is not a directory: set "
             "CPMDC_PSEUDO_DIR or CPMD_PP_LIBRARY_PATH");
    return -1;
  }
  snprintf(buf, cap, "%s", trimmed);
  return 0;
}

int cpmdc_species_order_map(int n_atoms, const int *atomic_numbers,
                            int n_species, const int *species_z,
                            const int *species_count, int *map_out) {
  int *used;
  int need = 0;
  int slot = 0;
  int s;
  if (n_atoms <= 0 || n_species <= 0 || !atomic_numbers || !species_z ||
      !species_count || !map_out)
    return -1;
  for (s = 0; s < n_species; ++s) {
    if (species_count[s] < 0 || need > INT_MAX - species_count[s])
      return -1;
    need += species_count[s];
  }
  if (need != n_atoms)
    return -1;
  used = (int *)calloc((size_t)n_atoms, sizeof(int));
  if (!used)
    return -1;
  for (s = 0; s < n_species; ++s) {
    int taken = 0;
    int j;
    for (j = 0; j < n_atoms && taken < species_count[s]; ++j) {
      if (used[j] || atomic_numbers[j] != species_z[s])
        continue;
      used[j] = 1;
      map_out[slot++] = j;
      taken++;
    }
    if (taken != species_count[s]) {
      free(used);
      return -1;
    }
  }
  free(used);
  return slot == n_atoms ? 0 : -1;
}

void cpmdc_scatter_species_gradient(int n_atoms, const int *map,
                                    const double *species_grad, double *grad) {
  int slot;
  if (n_atoms <= 0 || !map || !species_grad || !grad)
    return;
  for (slot = 0; slot < n_atoms; ++slot) {
    int atom = map[slot];
    int k;
    if (atom < 0 || atom >= n_atoms)
      continue;
    for (k = 0; k < 3; ++k)
      grad[atom * 3 + k] = species_grad[slot * 3 + k];
  }
}
