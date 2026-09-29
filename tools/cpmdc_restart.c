#include "cpmdc_restart.h"

#include <ctype.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void usage(void) {
  fprintf(stderr,
          "usage:\n"
          "  cpmdc-restart info FILE\n"
          "  cpmdc-restart positions FILE [--angstrom]\n"
          "  cpmdc-restart velocities FILE [--angstrom]\n"
          "  cpmdc-restart patch-positions FILE XYZ OUT [--angstrom]\n"
          "  cpmdc-restart patch-velocities FILE XYZ OUT [--angstrom]\n"
          "\n"
          "XYZ is one x y z triple per line, in Bohr unless --angstrom.\n"
          "Coordinates stay in CPMD species order. Other RESTART records,\n"
          "including the wavefunction, are copied unchanged.\n");
}

static int want_angstrom(int argc, char **argv) {
  for (int i = 0; i < argc; i++) {
    if (strcmp(argv[i], "--angstrom") == 0)
      return 1;
  }
  return 0;
}

static int print_triples(const cpmdc_restart *file, int which, int angstrom) {
  int n = cpmdc_restart_ncoords(file);
  double *xyz;
  int rc;
  if (which != 4)
    n = 0;
  if (which == 5 || which == 4) {
    int count = 0;
    if (cpmdc_restart_section(file, which, &count) != 0 || count < 0)
      return -1;
    n = count;
  }
  /* INT32-C: n*3 is a signed length. A count that does not fit is a different buffer. */
  if (n > INT_MAX / 3)
    return -1;
  xyz = calloc((size_t)(n > 0 ? n : 1) * 3, sizeof(double));
  if (!xyz)
    return -1;
  if (which == 4)
    rc = n ? cpmdc_restart_coordinates(file, xyz, n * 3) : 0;
  else
    rc = n ? cpmdc_restart_velocities(file, xyz, n * 3) : 0;
  if (rc != 0) {
    free(xyz);
    return -1;
  }
  printf("# unit=%s n=%d\n", angstrom ? "angstrom" : "bohr", n);
  for (int i = 0; i < n; i++) {
    double x = xyz[3 * i], y = xyz[3 * i + 1], z = xyz[3 * i + 2];
    if (angstrom) {
      x *= CPMDC_RESTART_ANGSTROM_PER_BOHR;
      y *= CPMDC_RESTART_ANGSTROM_PER_BOHR;
      z *= CPMDC_RESTART_ANGSTROM_PER_BOHR;
    }
    printf("%.17g %.17g %.17g\n", x, y, z);
  }
  free(xyz);
  return 0;
}

static int read_xyz(const char *path, double **out, int *n3, int angstrom, char *err,
                    size_t err_cap) {
  FILE *fp = fopen(path, "r");
  double *xyz = NULL;
  size_t n = 0, cap = 0;
  char line[512];
  if (!fp) {
    snprintf(err, err_cap, "cannot open %s", path);
    return -1;
  }
  while (fgets(line, sizeof line, fp)) {
    char *p = line;
    double v[3];
    int k = 0;
    while (*p && isspace((unsigned char)*p))
      p++;
    if (*p == '\0' || *p == '#')
      continue;
    while (k < 3 && *p) {
      char *end = NULL;
      v[k] = strtod(p, &end);
      if (end == p)
        break;
      p = end;
      k++;
    }
    if (k != 3) {
      snprintf(err, err_cap, "expected 3 floats per line in %s", path);
      free(xyz);
      fclose(fp);
      return -1;
    }
    if (n > SIZE_MAX - 3) {
      snprintf(err, err_cap, "out of memory");
      free(xyz);
      fclose(fp);
      return -1;
    }
    if (n + 3 > cap) {
      size_t ncap;
      /* INT30-C: a wrapped cap * 2 is a short realloc. */
      if (cap == 0)
        ncap = 48;
      else if (cap > SIZE_MAX / 2)
        ncap = 0;
      else
        ncap = cap * 2;
      if (ncap == 0 || ncap > SIZE_MAX / sizeof(double)) {
        snprintf(err, err_cap, "out of memory");
        free(xyz);
        fclose(fp);
        return -1;
      }
      double *grown = realloc(xyz, ncap * sizeof(double));
      if (!grown) {
        free(xyz);
        fclose(fp);
        snprintf(err, err_cap, "out of memory");
        return -1;
      }
      xyz = grown;
      cap = ncap;
    }
    if (angstrom) {
      v[0] /= CPMDC_RESTART_ANGSTROM_PER_BOHR;
      v[1] /= CPMDC_RESTART_ANGSTROM_PER_BOHR;
      v[2] /= CPMDC_RESTART_ANGSTROM_PER_BOHR;
    }
    xyz[n++] = v[0];
    xyz[n++] = v[1];
    xyz[n++] = v[2];
  }
  fclose(fp);
  if (n > (size_t)INT_MAX) {
    snprintf(err, err_cap, "coordinate count does not fit in int");
    free(xyz);
    return -1;
  }
  *out = xyz;
  *n3 = (int)n;
  return 0;
}

static int cmd_info(const cpmdc_restart *file) {
  double cell[6];
  int ibrav = 0, indpg = 0;
  const int *na = NULL;
  int nsp = 0;
  printf("header: %s\n", cpmdc_restart_header(file));
  printf("stream: %s\n", cpmdc_restart_is_stream(file) ? "yes" : "no");
  printf("coordinates: %d\n", cpmdc_restart_ncoords(file));
  if (cpmdc_restart_cell(file, &ibrav, &indpg, cell) == 0) {
    printf("cell: ibrav=%d indpg=%d\n", ibrav, indpg);
    printf("celldm: %.16g %.16g %.16g %.16g %.16g %.16g\n", cell[0], cell[1], cell[2],
           cell[3], cell[4], cell[5]);
  }
  if (cpmdc_restart_species(file, &nsp, &na) == 0) {
    printf("species: %d\n", nsp);
    for (int i = 0; i < nsp; i++)
      printf("  na[%d]=%d\n", i + 1, na[i]);
  }
  {
    double ecut, cdual;
    int dual_flag, nel, n1, n2, n3;
    if (cpmdc_restart_cutoff(file, &ecut, &cdual, &dual_flag, &nel, &n1, &n2, &n3) == 0)
      printf("cutoff: ecut=%.16g cdual=%.16g dual=%d nel=%d grid=%d %d %d\n", ecut, cdual,
             dual_flag, nel, n1, n2, n3);
  }
  {
    int n, nk, ngw, ngwl, nhg, nhgl;
    if (cpmdc_restart_states(file, &n, &nk, &ngw, &ngwl, &nhg, &nhgl) == 0)
      printf("states: n=%d nkpts=%d ngw=%d ngwl=%d nhg=%d nhgl=%d\n", n, nk, ngw, ngwl,
             nhg, nhgl);
  }
  for (int s = 1; s <= 99; s++) {
    int count = 0;
    if (cpmdc_restart_section(file, s, &count) == 0)
      printf("section %d: count=%d\n", s, count);
  }
  return 0;
}

int main(int argc, char **argv) {
  char err[256];
  cpmdc_restart *file;
  int angstrom;
  const char *cmd;
  if (argc < 3) {
    usage();
    return 2;
  }
  cmd = argv[1];
  angstrom = want_angstrom(argc, argv);
  err[0] = '\0';
  file = cpmdc_restart_read_path(argv[2], err, sizeof err);
  if (!file) {
    fprintf(stderr, "cpmdc-restart: %s\n", err[0] ? err : "read failed");
    return 1;
  }
  if (strcmp(cmd, "info") == 0) {
    cmd_info(file);
  } else if (strcmp(cmd, "positions") == 0) {
    if (print_triples(file, 4, angstrom) != 0) {
      fprintf(stderr, "cpmdc-restart: cannot read coordinates\n");
      cpmdc_restart_free(file);
      return 1;
    }
  } else if (strcmp(cmd, "velocities") == 0) {
    if (print_triples(file, 5, angstrom) != 0) {
      fprintf(stderr, "cpmdc-restart: cannot read velocities\n");
      cpmdc_restart_free(file);
      return 1;
    }
  } else if (strcmp(cmd, "patch-positions") == 0 || strcmp(cmd, "patch-velocities") == 0) {
    double *xyz = NULL;
    int n3 = 0;
    const char *out;
    if (argc < 5) {
      usage();
      cpmdc_restart_free(file);
      return 2;
    }
    out = argv[4];
    if (read_xyz(argv[3], &xyz, &n3, angstrom, err, sizeof err) != 0) {
      fprintf(stderr, "cpmdc-restart: %s\n", err);
      cpmdc_restart_free(file);
      return 1;
    }
    if (strcmp(cmd, "patch-positions") == 0) {
      if (cpmdc_restart_set_coordinates(file, xyz, n3) != 0) {
        fprintf(stderr, "cpmdc-restart: coordinate count does not match section 4\n");
        free(xyz);
        cpmdc_restart_free(file);
        return 1;
      }
    } else if (cpmdc_restart_set_velocities(file, xyz, n3) != 0) {
      fprintf(stderr, "cpmdc-restart: velocity count does not match section 5\n");
      free(xyz);
      cpmdc_restart_free(file);
      return 1;
    }
    free(xyz);
    if (cpmdc_restart_write_path(file, out, err, sizeof err) != 0) {
      fprintf(stderr, "cpmdc-restart: %s\n", err[0] ? err : "write failed");
      cpmdc_restart_free(file);
      return 1;
    }
  } else {
    usage();
    cpmdc_restart_free(file);
    return 2;
  }
  cpmdc_restart_free(file);
  return 0;
}
