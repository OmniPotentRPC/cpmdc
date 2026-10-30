#include "cpmdc.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv) {
  if (argc != 2) return 2;
  FILE *file = fopen(argv[1], "rb");
  if (!file || fseek(file, 0, SEEK_END)) return 2;
  long size = ftell(file);
  if (size <= 0 || fseek(file, 0, SEEK_SET)) return 2;
  unsigned char *bytes = malloc((size_t)size);
  if (!bytes || fread(bytes, 1, (size_t)size, file) != (size_t)size) return 2;
  fclose(file);
  CPMDCSession *session = cpmdc_session_create(bytes, (size_t)size);
  free(bytes);
  if (!session || cpmdc_session_clear_result_cache(NULL) != -1 ||
      cpmdc_session_last_result_reused(NULL) != -1 ||
      cpmdc_session_last_result_reused(session) != 0) return 1;
  const int z[3] = {8, 1, 1};
  const double q[9] = {0,0,.1173,0,.7572,-.4692,0,-.7572,-.4692};
  double gradients[3][9], energies[3];
  for (int i = 0; i < 3; ++i) {
    if (i == 2 && cpmdc_session_clear_result_cache(session) != 0) return 1;
    if (i == 2 && cpmdc_session_last_result_reused(session) != 0) return 1;
    if (cpmdc_session_select_orbitals(session, 100) != 0) return 1;
    CPMDCResult result = cpmdc_session_energy_gradient(session, 3, q, z, gradients[i]);
    if (!result.ok || !isfinite(result.energy_h) ||
        cpmdc_session_last_result_reused(session) != (i == 1)) {
      fprintf(stderr, "call %d: %s\n", i, result.message); return 1;
    }
    energies[i] = result.energy_h;
    for (int k = 0; k < 9; ++k) if (!isfinite(gradients[i][k])) return 1;
  }
  if (energies[0] != energies[1] || fabs(energies[0]-energies[2]) > 1e-7) return 1;
  for (int k = 0; k < 9; ++k)
    if (gradients[0][k] != gradients[1][k] ||
        fabs(gradients[0][k]-gradients[2][k]) > 1e-4) return 1;
  cpmdc_session_destroy(session);
  puts("session_result_reuse=passed exact_geometry_refresh=passed orbital_key_preserved=passed");
  return 0;
}
