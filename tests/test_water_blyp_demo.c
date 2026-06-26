/* Water BLYP through shipped C ABI; requires -Dwith_cpmd + libcpmd.a. */
#include "cpmdc.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static unsigned char *read_all(const char *path, size_t *sz) {
  FILE *fp = fopen(path, "rb");
  if (!fp) return NULL;
  fseek(fp, 0, SEEK_END);
  long n = ftell(fp);
  rewind(fp);
  unsigned char *b = malloc((size_t)n);
  if (!b || fread(b, 1, (size_t)n, fp) != (size_t)n) { free(b); fclose(fp); return NULL; }
  fclose(fp);
  *sz = (size_t)n;
  return b;
}

int main(int argc, char **argv) {
  const char *out_path = argc > 3 ? argv[3] : "water_blyp_demo.txt";
  if (argc < 3) {
    fprintf(stderr, "usage: %s params.bin force_water.bin [out.txt]\n", argv[0]);
    return 2;
  }
  if (!cpmdc_available()) {
    fprintf(stderr, "cpmdc_available()==0 (need libcpmd embed build)\n");
    return 1;
  }
  size_t psz = 0, fsz = 0;
  unsigned char *params = read_all(argv[1], &psz);
  unsigned char *force = read_all(argv[2], &fsz);
  if (!params || !force) return 1;
  size_t need = cpmdc_potential_result_size_for_force_input(force, fsz);
  unsigned char *res = calloc(1, need ? need : 1);
  size_t rsz = 0;
  CPMDCResult r = cpmdc_calculate_result(params, psz, force, fsz, res, need, &rsz);
  FILE *o = fopen(out_path, "w");
  fprintf(o, "available=1\n");
  fprintf(o, "ok=%d\n", r.ok);
  fprintf(o, "energy_h=%.17g\n", r.energy_h);
  fprintf(o, "message=%s\n", r.message);
  fprintf(o, "result_bytes=%zu\n", rsz);
  /* Also session path (optimizer style, one step) */
  CPMDCSession *s = cpmdc_session_create(params, psz);
  double forces[9] = {0};
  CPMDCResult r2 = {0};
  if (s) {
    r2 = cpmdc_session_calculate_forces(s, force, fsz, forces, 9);
    fprintf(o, "session_ok=%d\n", r2.ok);
    fprintf(o, "session_energy_h=%.17g\n", r2.energy_h);
    for (int i = 0; i < 9; ++i)
      fprintf(o, "force[%d]=%.17g finite=%d\n", i, forces[i], isfinite(forces[i]));
    cpmdc_session_destroy(s);
  }
  fclose(o);
  printf("wrote %s ok=%d energy_h=%.10g\n", out_path, r.ok, r.energy_h);
  free(res); free(params); free(force);
  cpmdc_finalize();
  /* Success criteria for partial embed: available + session created; SCF may still fail */
  return (cpmdc_available() == 0 && r.ok == 0 && !s) ? 1 : 0;
}
