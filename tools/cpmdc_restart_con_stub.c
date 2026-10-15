#include "cpmdc_restart_con.h"

#include <stdio.h>

static int unavailable(char *err, size_t err_cap) {
  if (err && err_cap)
    snprintf(err, err_cap, "cpmdc-restart was built without readcon");
  return -1;
}

int cpmdc_restart_con_export(const cpmdc_restart *file, const char *con_path,
                             const char *const *atom_symbols, int n_atoms, char *err,
                             size_t err_cap) {
  (void)file;
  (void)con_path;
  (void)atom_symbols;
  (void)n_atoms;
  return unavailable(err, err_cap);
}

int cpmdc_restart_con_import(cpmdc_restart *file, const char *con_path, char *err,
                             size_t err_cap) {
  (void)file;
  (void)con_path;
  return unavailable(err, err_cap);
}
