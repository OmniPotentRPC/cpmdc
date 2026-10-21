#include "cpmdc_restart_con.h"

#include "readcon-core.h"

#include <ctype.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* One mass per symbol. Conventional atomic weights, shared by every atom of that symbol. */
static const struct {
  const char *symbol;
  double mass;
} k_masses[] = {
    {"H", 1.00784},  {"He", 4.002602}, {"Li", 6.94},      {"Be", 9.0121831},
    {"B", 10.81},    {"C", 12.011},    {"N", 14.007},     {"O", 15.999},
    {"F", 18.998403163}, {"Ne", 20.1797}, {"Na", 22.98976928}, {"Mg", 24.305},
    {"Al", 26.9815385}, {"Si", 28.085}, {"P", 30.973761998}, {"S", 32.06},
    {"Cl", 35.45},   {"Ar", 39.948},   {"K", 39.0983},    {"Ca", 40.078},
    {"Fe", 55.845},  {"Cu", 63.546},   {"Ge", 72.630},
};

static const char k_units_positions[] =
    "{\"length\":\"angstrom\",\"energy\":\"eV\",\"mass\":\"amu\"}";
static const char k_units_velocities[] =
    "{\"length\":\"angstrom\",\"energy\":\"eV\",\"mass\":\"amu\","
    "\"time\":\"(amu*angstrom^2/eV)^0.5\"}";

static void set_err(char *err, size_t cap, const char *msg) {
  if (!err || cap == 0)
    return;
  snprintf(err, cap, "%s", msg);
}

double cpmdc_restart_bohr_to_angstrom(double bohr) {
  return bohr * CPMDC_RESTART_ANGSTROM_PER_BOHR;
}

double cpmdc_restart_angstrom_to_bohr(double angstrom) {
  return angstrom / CPMDC_RESTART_ANGSTROM_PER_BOHR;
}

static double ratio_or_one(double ratio) {
  return ratio == 0.0 ? 1.0 : ratio;
}

static double cosine_to_degrees(double cosine) {
  if (cosine == 0.0)
    return 90.0;
  if (cosine == 1.0)
    return 0.0;
  if (cosine == -1.0)
    return 180.0;
  return acos(cosine) * (180.0 / acos(-1.0));
}

static double degrees_to_cosine(double degrees) {
  if (degrees == 90.0)
    return 0.0;
  if (degrees == 0.0)
    return 1.0;
  if (degrees == 180.0 || degrees == -180.0)
    return -1.0;
  return cos(degrees * (acos(-1.0) / 180.0));
}

int cpmdc_restart_celldm_to_con(const double celldm[6], double lengths[3],
                                double angles[3]) {
  double a, b, c;
  if (!celldm || !lengths || !angles)
    return -1;
  for (int k = 3; k < 6; k++) {
    if (celldm[k] < -1.0 || celldm[k] > 1.0)
      return -1;
  }
  a = celldm[0];
  b = a * ratio_or_one(celldm[1]);
  c = a * ratio_or_one(celldm[2]);
  lengths[0] = cpmdc_restart_bohr_to_angstrom(a);
  lengths[1] = cpmdc_restart_bohr_to_angstrom(b);
  lengths[2] = cpmdc_restart_bohr_to_angstrom(c);
  angles[0] = cosine_to_degrees(celldm[3]);
  angles[1] = cosine_to_degrees(celldm[4]);
  angles[2] = cosine_to_degrees(celldm[5]);
  return 0;
}

int cpmdc_restart_con_to_celldm(const double lengths[3], const double angles[3],
                                double celldm[6]) {
  double a, b, c;
  if (!lengths || !angles || !celldm)
    return -1;
  a = cpmdc_restart_angstrom_to_bohr(lengths[0]);
  b = cpmdc_restart_angstrom_to_bohr(lengths[1]);
  c = cpmdc_restart_angstrom_to_bohr(lengths[2]);
  celldm[0] = a;
  celldm[1] = a == 0.0 ? 0.0 : b / a;
  celldm[2] = a == 0.0 ? 0.0 : c / a;
  celldm[3] = degrees_to_cosine(angles[0]);
  celldm[4] = degrees_to_cosine(angles[1]);
  celldm[5] = degrees_to_cosine(angles[2]);
  return 0;
}

static int symbol_ok(const char *symbol) {
  size_t n;
  if (!symbol)
    return 0;
  n = strlen(symbol);
  if (n < 1 || n > 3)
    return 0;
  for (size_t i = 0; i < n; i++) {
    if (!isalpha((unsigned char)symbol[i]))
      return 0;
  }
  return 1;
}

static int mass_of(const char *symbol, double *mass) {
  if (!symbol_ok(symbol))
    return -1;
  for (size_t i = 0; i < sizeof k_masses / sizeof k_masses[0]; i++) {
    if (strcmp(k_masses[i].symbol, symbol) == 0) {
      *mass = k_masses[i].mass;
      return 0;
    }
  }
  return -1;
}

static void scale_bohr(double *xyz, int n3) {
  for (int i = 0; i < n3; i++)
    xyz[i] = cpmdc_restart_bohr_to_angstrom(xyz[i]);
}

int cpmdc_restart_con_export(const cpmdc_restart *file, const char *con_path,
                             const char *const *atom_symbols, int n_atoms, char *err,
                             size_t err_cap) {
  double *xyz = NULL;
  double *vel = NULL;
  double celldm[6];
  double lengths[3];
  double angles[3];
  int n;
  int n3;
  int vel_count = 0;
  int with_vel = 0;
  struct RKRConFrameBuilder *builder = NULL;
  struct RKRConFrame *frame = NULL;
  struct RKRConFrameWriter *writer = NULL;
  const struct RKRConFrame *frames[1];
  int rc = -1;

  if (!file || !con_path || !atom_symbols) {
    set_err(err, err_cap, "null export argument");
    return -1;
  }
  n = cpmdc_restart_ncoords(file);
  if (n <= 0 || n_atoms != n || n > INT_MAX / 3) {
    set_err(err, err_cap, "symbol count does not match section 4");
    return -1;
  }
  n3 = n * 3;
  if (cpmdc_restart_cell(file, NULL, NULL, celldm) != 0 ||
      cpmdc_restart_celldm_to_con(celldm, lengths, angles) != 0) {
    set_err(err, err_cap, "restart cell is missing or unusable");
    return -1;
  }
  xyz = calloc((size_t)n3, sizeof(double));
  if (!xyz) {
    set_err(err, err_cap, "out of memory");
    return -1;
  }
  if (cpmdc_restart_coordinates(file, xyz, n3) != 0) {
    set_err(err, err_cap, "cannot read coordinates");
    goto done;
  }
  scale_bohr(xyz, n3);
  if (cpmdc_restart_section(file, 5, &vel_count) == 0 && vel_count > 0) {
    if (vel_count != n) {
      set_err(err, err_cap, "velocity count does not match section 4");
      goto done;
    }
    vel = calloc((size_t)n3, sizeof(double));
    if (!vel) {
      set_err(err, err_cap, "out of memory");
      goto done;
    }
    if (cpmdc_restart_velocities(file, vel, n3) != 0) {
      set_err(err, err_cap, "cannot read velocities");
      goto done;
    }
    scale_bohr(vel, n3);
    with_vel = 1;
  }
  builder = rkr_frame_new(lengths, angles, "CPMDC RESTART", NULL, NULL, NULL);
  if (!builder) {
    set_err(err, err_cap, "cannot build a CON frame");
    goto done;
  }
  if (rkr_frame_builder_set_units(builder, with_vel ? k_units_velocities
                                                    : k_units_positions) !=
      RKR_STATUS_SUCCESS) {
    set_err(err, err_cap, "cannot set v3 units");
    goto done;
  }
  for (int i = 0; i < n; i++) {
    double mass = 0;
    const double *atom_vel = with_vel ? vel + 3 * i : NULL;
    if (mass_of(atom_symbols[i], &mass) != 0) {
      set_err(err, err_cap, "missing mass for a species symbol");
      goto done;
    }
    if (rkr_frame_add_atom_full(builder, atom_symbols[i], xyz[3 * i], xyz[3 * i + 1],
                                xyz[3 * i + 2], false, false, false, (uint64_t)i, mass,
                                atom_vel, NULL) != RKR_STATUS_SUCCESS) {
      set_err(err, err_cap, "cannot add an atom");
      goto done;
    }
  }
  frame = rkr_frame_builder_build(builder);
  builder = NULL;
  if (!frame) {
    set_err(err, err_cap, "CON build failed");
    goto done;
  }
  writer = create_writer_from_path_round_trip_c(con_path);
  if (!writer) {
    set_err(err, err_cap, "cannot open the CON writer");
    goto done;
  }
  if (rkr_writer_is_canonical(writer) != 0) {
    set_err(err, err_cap, "RoundTrip writer has canonical mode on");
    goto done;
  }
  frames[0] = frame;
  if (rkr_writer_extend(writer, frames, 1) != RKR_STATUS_SUCCESS ||
      rkr_writer_flush(writer) != RKR_STATUS_SUCCESS) {
    set_err(err, err_cap, "CON write failed");
    goto done;
  }
  rc = 0;
done:
  if (builder)
    free_rkr_frame_builder(builder);
  if (writer)
    free_rkr_writer(writer);
  if (frame)
    free_rkr_frame(frame);
  free(xyz);
  free(vel);
  return rc;
}

static int scatter_component(const double *grouped, const uint64_t *ids, uintptr_t n,
                             double *bohr) {
  unsigned char *seen;
  if (n > 100000000)
    return -1;
  seen = calloc(n, 1);
  if (!seen)
    return -1;
  for (uintptr_t i = 0; i < n; i++) {
    uint64_t id = ids[i];
    if (id >= n || seen[id]) {
      free(seen);
      return -1;
    }
    seen[id] = 1;
    for (int k = 0; k < 3; k++)
      bohr[3 * (int)id + k] = cpmdc_restart_angstrom_to_bohr(grouped[3 * i + k]);
  }
  free(seen);
  return 0;
}

int cpmdc_restart_con_import(cpmdc_restart *file, const char *con_path, char *err,
                             size_t err_cap) {
  struct RKRConFrame *frame = NULL;
  struct CFrame *view = NULL;
  const double *xyz = NULL;
  const double *vel = NULL;
  uint64_t *ids = NULL;
  double *bohr = NULL;
  double celldm[6];
  uintptr_t n = 0;
  uintptr_t nvel = 0;
  int ncoords;
  int rc = -1;

  if (!file || !con_path) {
    set_err(err, err_cap, "null import argument");
    return -1;
  }
  ncoords = cpmdc_restart_ncoords(file);
  if (ncoords <= 0 || ncoords > INT_MAX / 3) {
    set_err(err, err_cap, "restart has no coordinates");
    return -1;
  }
  frame = rkr_read_first_frame(con_path);
  if (!frame) {
    set_err(err, err_cap, "cannot read the CON file");
    return -1;
  }
  xyz = rkr_frame_xyz_f64(frame, &n);
  if (!xyz || n != (uintptr_t)ncoords) {
    set_err(err, err_cap, "CON atom count does not match section 4");
    goto done;
  }
  ids = calloc(n, sizeof(uint64_t));
  bohr = calloc(n * 3, sizeof(double));
  if (!ids || !bohr) {
    set_err(err, err_cap, "out of memory");
    goto done;
  }
  if (rkr_frame_copy_atom_ids(frame, ids, n) != RKR_STATUS_SUCCESS) {
    set_err(err, err_cap, "cannot read atom ids");
    goto done;
  }
  if (scatter_component(xyz, ids, n, bohr) != 0 ||
      cpmdc_restart_set_coordinates(file, bohr, ncoords * 3) != 0) {
    set_err(err, err_cap, "cannot scatter coordinates");
    goto done;
  }
  vel = rkr_frame_velocities_f64(frame, &nvel);
  if (vel) {
    if (nvel != n || scatter_component(vel, ids, n, bohr) != 0 ||
        cpmdc_restart_set_velocities(file, bohr, ncoords * 3) != 0) {
      set_err(err, err_cap, "cannot scatter velocities");
      goto done;
    }
  }
  view = rkr_frame_to_c_frame(frame);
  if (!view || cpmdc_restart_con_to_celldm(view->cell, view->angles, celldm) != 0 ||
      cpmdc_restart_set_cell(file, celldm) != 0) {
    set_err(err, err_cap, "cannot set the cell from the CON file");
    goto done;
  }
  rc = 0;
done:
  if (view)
    free_c_frame(view);
  free(ids);
  free(bohr);
  if (frame)
    free_rkr_frame(frame);
  return rc;
}
