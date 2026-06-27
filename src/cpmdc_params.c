#include "cpmdc_params.h"

#include <ctype.h>
#include <limits.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>

static const capn_text empty_text = {0, "", NULL};

static int pointer_list_len(capn_ptr *ptr) {
  capn_resolve(ptr);
  if (ptr->type == CAPN_NULL)
    return 0;
  if (ptr->type != CAPN_PTR_LIST)
    return -1;
  return ptr->len;
}

static int struct_list_len(capn_ptr *ptr) {
  capn_resolve(ptr);
  if (ptr->type == CAPN_NULL)
    return 0;
  if (ptr->type != CAPN_LIST)
    return -1;
  return ptr->len;
}

static int list64_len(capn_list64 *list) {
  capn_resolve(&list->p);
  if (list->p.type == CAPN_NULL)
    return 0;
  if (list->p.type != CAPN_LIST || list->p.datasz != 8)
    return -1;
  return list->p.len;
}

static int list32_len(capn_list32 *list) {
  capn_resolve(&list->p);
  if (list->p.type == CAPN_NULL)
    return 0;
  if (list->p.type != CAPN_LIST || list->p.datasz != 4)
    return -1;
  return list->p.len;
}

static int text_eq_ci(capn_text text, const char *lit) {
  size_t n = strlen(lit);
  if (text.len != (int)n || !text.str)
    return 0;
  for (size_t i = 0; i < n; ++i) {
    char a = (char)tolower((unsigned char)text.str[i]);
    char b = (char)tolower((unsigned char)lit[i]);
    if (a != b)
      return 0;
  }
  return 1;
}

static int text_starts_ci(capn_text text, const char *prefix) {
  size_t n = strlen(prefix);
  if (text.len < (int)n || !text.str)
    return 0;
  for (size_t i = 0; i < n; ++i) {
    char a = (char)tolower((unsigned char)text.str[i]);
    char b = (char)tolower((unsigned char)prefix[i]);
    if (a != b)
      return 0;
  }
  return 1;
}

static int append_fmt(char *dst, size_t dst_size, size_t *used, const char *fmt,
                      ...) {
  if (!dst || !used || *used >= dst_size)
    return -1;
  va_list ap;
  va_start(ap, fmt);
  int n = vsnprintf(dst + *used, dst_size - *used, fmt, ap);
  va_end(ap);
  if (n < 0 || (size_t)n >= dst_size - *used)
    return -1;
  *used += (size_t)n;
  return 0;
}

static int append_text(char *dst, size_t dst_size, size_t *used, const char *s) {
  return append_fmt(dst, dst_size, used, "%s", s ? s : "");
}

static int append_capn_text(char *dst, size_t dst_size, size_t *used,
                            capn_text text) {
  if (!text.str || text.len <= 0)
    return 0;
  if (*used + (size_t)text.len >= dst_size)
    return -1;
  memcpy(dst + *used, text.str, (size_t)text.len);
  *used += (size_t)text.len;
  dst[*used] = '\0';
  return 0;
}

static int append_slice(char *dst, size_t dst_size, size_t *used,
                        const char *s, size_t len) {
  if (!s || len == 0)
    return 0;
  if (*used + len >= dst_size)
    return -1;
  memcpy(dst + *used, s, len);
  *used += len;
  dst[*used] = '\0';
  return 0;
}

static int append_directives(char *dst, size_t dst_size, size_t *used,
                             CPMDDirective_list directives) {
  int n = struct_list_len(&directives.p);
  if (n < 0)
    return -1;
  for (int i = 0; i < n; ++i) {
    struct CPMDDirective d;
    get_CPMDDirective(&d, directives, i);
    if (append_text(dst, dst_size, used, " ") != 0)
      return -1;
    if (append_capn_text(dst, dst_size, used, d.keyword) != 0)
      return -1;
    if (append_text(dst, dst_size, used, "\n") != 0)
      return -1;
    int na = pointer_list_len(&d.args);
    if (na < 0)
      return -1;
    for (int j = 0; j < na; ++j) {
      capn_text arg = capn_get_text(d.args, j, empty_text);
      if (append_text(dst, dst_size, used, "  ") != 0)
        return -1;
      if (append_capn_text(dst, dst_size, used, arg) != 0)
        return -1;
      if (append_text(dst, dst_size, used, "\n") != 0)
        return -1;
    }
  }
  return 0;
}

#define CPMDC_MAX_SET_DIRECTIVES 128

struct RenderSetDirective {
  const char *section;
  size_t section_len;
  const char *keyword;
  size_t keyword_len;
  capn_text value;
  int rendered;
};

struct RenderSetList {
  struct RenderSetDirective items[CPMDC_MAX_SET_DIRECTIVES];
  int len;
};

static int ascii_upper(int c) {
  if (c >= 'a' && c <= 'z')
    return c - 'a' + 'A';
  return c;
}

static int section_name_equals(const struct RenderSetDirective *set,
                               const char *name) {
  size_t name_len = strlen(name);
  if (set->section_len != name_len)
    return 0;
  for (size_t i = 0; i < name_len; ++i) {
    if (ascii_upper((unsigned char)set->section[i]) !=
        ascii_upper((unsigned char)name[i]))
      return 0;
  }
  return 1;
}

static int set_sections_equal(const struct RenderSetDirective *a,
                              const struct RenderSetDirective *b) {
  if (a->section_len != b->section_len)
    return 0;
  for (size_t i = 0; i < a->section_len; ++i) {
    if (ascii_upper((unsigned char)a->section[i]) !=
        ascii_upper((unsigned char)b->section[i]))
      return 0;
  }
  return 1;
}

static int slice_starts_ci(const char *text, size_t text_len,
                           const char *prefix) {
  size_t prefix_len = strlen(prefix);
  if (!text || text_len < prefix_len)
    return 0;
  for (size_t i = 0; i < prefix_len; ++i) {
    if (ascii_upper((unsigned char)text[i]) !=
        ascii_upper((unsigned char)prefix[i]))
      return 0;
  }
  return 1;
}

static int init_render_set(struct RenderSetDirective *out,
                           const struct CPMDSetDirective *set) {
  if (!set->key.str || set->key.len <= 2)
    return -1;
  const char *dot = memchr(set->key.str, '.', (size_t)set->key.len);
  if (!dot || dot == set->key.str ||
      dot == set->key.str + (size_t)set->key.len - 1)
    return -1;
  out->section = set->key.str;
  out->section_len = (size_t)(dot - set->key.str);
  out->keyword = dot + 1;
  out->keyword_len = (size_t)set->key.len - out->section_len - 1u;
  out->value = set->value;
  out->rendered = 0;
  return 0;
}

static int collect_set_directives(CPMDInputSection_list sections,
                                  struct RenderSetList *sets) {
  int nsec = struct_list_len(&sections.p);
  if (nsec < 0)
    return -1;
  sets->len = 0;
  for (int i = 0; i < nsec; ++i) {
    struct CPMDInputSection sec;
    get_CPMDInputSection(&sec, sections, i);
    if (sec.which != CPMDInputSection_set)
      continue;
    if (sets->len >= CPMDC_MAX_SET_DIRECTIVES)
      return -1;
    struct CPMDSetDirective body;
    read_CPMDSetDirective(&body, sec.set);
    if (init_render_set(&sets->items[sets->len], &body) != 0)
      return -1;
    ++sets->len;
  }
  return 0;
}

static int set_directives_have_prefix(struct RenderSetList *sets,
                                      const char *section,
                                      const char *prefix) {
  if (!sets)
    return 0;
  for (int i = 0; i < sets->len; ++i) {
    if (sets->items[i].rendered ||
        !section_name_equals(&sets->items[i], section))
      continue;
    if (slice_starts_ci(sets->items[i].keyword, sets->items[i].keyword_len,
                        prefix))
      return 1;
  }
  return 0;
}

static int append_set_payload(char *dst, size_t dst_size, size_t *used,
                              const struct RenderSetDirective *set) {
  if (append_text(dst, dst_size, used, " ") != 0)
    return -1;
  if (append_slice(dst, dst_size, used, set->keyword, set->keyword_len) != 0)
    return -1;
  if (append_text(dst, dst_size, used, "\n") != 0)
    return -1;
  if (set->value.str && set->value.len > 0) {
    if (append_text(dst, dst_size, used, "  ") != 0)
      return -1;
    if (append_capn_text(dst, dst_size, used, set->value) != 0)
      return -1;
    if (append_text(dst, dst_size, used, "\n") != 0)
      return -1;
  }
  return 0;
}

static int append_set_directives_for_section(char *dst, size_t dst_size,
                                             size_t *used,
                                             struct RenderSetList *sets,
                                             const char *section) {
  if (!sets)
    return 0;
  for (int i = 0; i < sets->len; ++i) {
    if (sets->items[i].rendered || !section_name_equals(&sets->items[i], section))
      continue;
    if (append_set_payload(dst, dst_size, used, &sets->items[i]) != 0)
      return -1;
    sets->items[i].rendered = 1;
  }
  return 0;
}

static capn_text selected_file_path(const struct CPMDParams *view) {
  if (view->permanentDir.str && view->permanentDir.len > 0)
    return view->permanentDir;
  if (view->scratchDir.str && view->scratchDir.len > 0)
    return view->scratchDir;
  return empty_text;
}

static int append_file_path_directive(char *dst, size_t dst_size, size_t *used,
                                      capn_text path) {
  if (!path.str || path.len <= 0)
    return 0;
  if (append_text(dst, dst_size, used, " FILEPATH\n  ") != 0)
    return -1;
  if (append_capn_text(dst, dst_size, used, path) != 0)
    return -1;
  if (path.str[path.len - 1] != '/') {
    if (append_text(dst, dst_size, used, "/") != 0)
      return -1;
  }
  return append_text(dst, dst_size, used, "\n");
}

static int render_remaining_set_sections(char *dst, size_t dst_size,
                                         size_t *used,
                                         struct RenderSetList *sets) {
  if (!sets)
    return 0;
  for (int i = 0; i < sets->len; ++i) {
    if (sets->items[i].rendered)
      continue;
    if (append_text(dst, dst_size, used, "&") != 0)
      return -1;
    if (append_slice(dst, dst_size, used, sets->items[i].section,
                     sets->items[i].section_len) != 0)
      return -1;
    if (append_text(dst, dst_size, used, "\n") != 0)
      return -1;
    for (int j = i; j < sets->len; ++j) {
      if (sets->items[j].rendered ||
          !set_sections_equal(&sets->items[i], &sets->items[j]))
        continue;
      if (append_set_payload(dst, dst_size, used, &sets->items[j]) != 0)
        return -1;
      sets->items[j].rendered = 1;
    }
    if (append_text(dst, dst_size, used, "&END\n\n") != 0)
      return -1;
  }
  return 0;
}

static int directives_have_prefix(CPMDDirective_list directives,
                                  const char *prefix) {
  int n = struct_list_len(&directives.p);
  if (n < 0)
    return -1;
  for (int i = 0; i < n; ++i) {
    struct CPMDDirective d;
    get_CPMDDirective(&d, directives, i);
    if (text_starts_ci(d.keyword, prefix))
      return 1;
  }
  return 0;
}

static int append_cpmd_cell(char *dst, size_t dst_size, size_t *used,
                            const double *cell, int ncell) {
  if (ncell != 6 && ncell != 9)
    return 0;
  if (append_text(dst, dst_size, used, " CELL\n ") != 0)
    return -1;
  if (ncell == 6) {
    double a = cell[0];
    double b = cell[1];
    double c = cell[2];
    double ab = cell[3];
    double ac = cell[4];
    double bc = cell[5];
    if (fabs(ab) > 1.0 || fabs(ac) > 1.0 || fabs(bc) > 1.0) {
      const double deg = 3.14159265358979323846 / 180.0;
      if (a != 0.0) {
        b /= a;
        c /= a;
      }
      ab = cos(ab * deg);
      ac = cos(ac * deg);
      bc = cos(bc * deg);
    }
    return append_fmt(dst, dst_size, used, " %.10g %.10g %.10g %.10g %.10g %.10g\n",
                      a, b, c, ab, ac, bc);
  }
  if (fabs(cell[1]) < 1e-12 && fabs(cell[2]) < 1e-12 &&
      fabs(cell[3]) < 1e-12 && fabs(cell[5]) < 1e-12 &&
      fabs(cell[6]) < 1e-12 && fabs(cell[7]) < 1e-12 && cell[0] != 0.0) {
    return append_fmt(dst, dst_size, used, " %.10g %.10g %.10g 0 0 0\n",
                      cell[0], cell[4] / cell[0], cell[8] / cell[0]);
  }
  for (int i = 0; i < ncell; ++i) {
    if (append_fmt(dst, dst_size, used, " %.10g", cell[i]) != 0)
      return -1;
  }
  return append_text(dst, dst_size, used, "\n");
}

static int render_system_section_with_cell(
    char *dst, size_t dst_size, size_t *used, struct CPMDSystemSection *sys,
    double default_cutoff, int default_charge, const double *cell_override,
    int override_ncell, struct RenderSetList *sets) {
  if (append_text(dst, dst_size, used, "&SYSTEM\n") != 0)
    return -1;
  int symmetry = sys->symmetry;
  int charge = sys->charge != 0 ? sys->charge : default_charge;
  double cutoff = sys->cutOffRy > 0.0 ? sys->cutOffRy : default_cutoff;
  if (append_fmt(dst, dst_size, used, " SYMMETRY\n  %d\n", symmetry) != 0)
    return -1;
  if (sys->angstrom) {
    if (append_text(dst, dst_size, used, " ANGSTROM\n") != 0)
      return -1;
  }
  int ncell = list64_len(&sys->cell);
  if (ncell < 0)
    return -1;
  if (override_ncell == 9 && cell_override) {
    if (append_cpmd_cell(dst, dst_size, used, cell_override, override_ncell) != 0)
      return -1;
  } else if (ncell == 6 || ncell == 9) {
    double cell[9] = {0};
    for (int i = 0; i < ncell; ++i)
      cell[i] = capn_to_f64(capn_get64(sys->cell, i));
    if (append_cpmd_cell(dst, dst_size, used, cell, ncell) != 0)
      return -1;
  }
  if (append_fmt(dst, dst_size, used, " CUTOFF\n  %.10g\n", cutoff) != 0)
    return -1;
  if (sys->scale != 0.0) {
    if (append_fmt(dst, dst_size, used, " SCALE\n  %.10g\n", sys->scale) != 0)
      return -1;
  }
  if (charge != 0) {
    if (append_fmt(dst, dst_size, used, " CHARGE\n  %d\n", charge) != 0)
      return -1;
  }
  int has_poisson = directives_have_prefix(sys->directives, "POISSON SOLVER");
  if (has_poisson < 0)
    return -1;
  int set_has_poisson =
      set_directives_have_prefix(sets, "SYSTEM", "POISSON SOLVER");
  if (set_has_poisson < 0)
    return -1;
  if (symmetry == 0 && !has_poisson && !set_has_poisson) {
    if (append_text(dst, dst_size, used, " POISSON SOLVER HOCKNEY\n") != 0)
      return -1;
  }
  if (append_directives(dst, dst_size, used, sys->directives) != 0)
    return -1;
  if (append_set_directives_for_section(dst, dst_size, used, sets, "SYSTEM") !=
      0)
    return -1;
  return append_text(dst, dst_size, used, "&END\n\n");
}

static int render_system_section(char *dst, size_t dst_size, size_t *used,
                                 struct CPMDSystemSection *sys,
                                 double default_cutoff, int default_charge,
                                 struct RenderSetList *sets) {
  return render_system_section_with_cell(dst, dst_size, used, sys,
                                         default_cutoff, default_charge, NULL, 0,
                                         sets);
}

static int render_cpmd_section(char *dst, size_t dst_size, size_t *used,
                               const struct CPMDCpmdSection *sec,
                               const char *task, capn_text file_path,
                               struct RenderSetList *sets) {
  if (append_text(dst, dst_size, used, "&CPMD\n") != 0)
    return -1;
  int opt_wf = sec->optimizeWavefunction;
  int do_md = sec->molecularDynamics;
  if (task && (strcmp(task, "md") == 0 || strcmp(task, "MD") == 0))
    do_md = 1;
  if (task && (strcmp(task, "energy") == 0 || strcmp(task, "gradient") == 0 ||
               strcmp(task, "optimize") == 0))
    opt_wf = 1;
  if (opt_wf && !do_md) {
    if (append_text(dst, dst_size, used, " OPTIMIZE WAVEFUNCTION\n") != 0)
      return -1;
  }
  if (do_md) {
    if (append_text(dst, dst_size, used, " MOLECULAR DYNAMICS\n") != 0)
      return -1;
  }
  if (sec->convergenceOrbitals > 0.0) {
    if (append_fmt(dst, dst_size, used, " CONVERGENCE ORBITALS\n  %.10g\n",
                   sec->convergenceOrbitals) != 0)
      return -1;
  }
  if (sec->maxStep > 0) {
    if (append_fmt(dst, dst_size, used, " MAXSTEP\n  %d\n", sec->maxStep) != 0)
      return -1;
  }
  if (sec->timestep > 0.0) {
    if (append_fmt(dst, dst_size, used, " TIMESTEP\n  %.10g\n", sec->timestep) !=
        0)
      return -1;
  }
  if (sec->restartWavefunction) {
    if (append_text(dst, dst_size, used, " RESTART WAVEFUNCTION\n") != 0)
      return -1;
  }
  if (sec->trajectory) {
    if (append_text(dst, dst_size, used, " TRAJECTORY\n") != 0)
      return -1;
  }
  if (append_directives(dst, dst_size, used, sec->directives) != 0)
    return -1;
  int has_filepath = directives_have_prefix(sec->directives, "FILEPATH");
  if (has_filepath < 0)
    return -1;
  int set_has_filepath = set_directives_have_prefix(sets, "CPMD", "FILEPATH");
  if (set_has_filepath < 0)
    return -1;
  if (!has_filepath && !set_has_filepath) {
    if (append_file_path_directive(dst, dst_size, used, file_path) != 0)
      return -1;
  }
  if (append_set_directives_for_section(dst, dst_size, used, sets, "CPMD") != 0)
    return -1;
  return append_text(dst, dst_size, used, "&END\n\n");
}

static int render_dft_section(char *dst, size_t dst_size, size_t *used,
                              const struct CPMDDftSection *dft,
                              const char *default_functional, int multiplicity,
                              struct RenderSetList *sets) {
  if (append_text(dst, dst_size, used, "&DFT\n") != 0)
    return -1;
  if (dft->functional.str && dft->functional.len > 0) {
    if (append_text(dst, dst_size, used, " FUNCTIONAL ") != 0)
      return -1;
    if (append_capn_text(dst, dst_size, used, dft->functional) != 0)
      return -1;
    if (append_text(dst, dst_size, used, "\n") != 0)
      return -1;
  } else {
    if (append_fmt(dst, dst_size, used, " FUNCTIONAL %s\n",
                   default_functional ? default_functional : "BLYP") != 0)
      return -1;
  }
  if (dft->lsd || multiplicity > 1) {
    if (append_text(dst, dst_size, used, " LSD\n") != 0)
      return -1;
  }
  if (append_directives(dst, dst_size, used, dft->directives) != 0)
    return -1;
  if (append_set_directives_for_section(dst, dst_size, used, sets, "DFT") != 0)
    return -1;
  return append_text(dst, dst_size, used, "&END\n\n");
}

/* CPMD &ATOMS: *filename.psp then LMAX=S|P|D and coordinates (Angstrom). */
static const char *lmax_letter(int lmax) {
  if (lmax <= 0)
    return "S";
  if (lmax == 1)
    return "P";
  if (lmax == 2)
    return "D";
  return "F";
}

static int render_atoms_section(char *dst, size_t dst_size, size_t *used,
                                struct CPMDAtomsSection *atoms,
                                struct RenderSetList *sets) {
  if (append_text(dst, dst_size, used, "&ATOMS\n") != 0)
    return -1;
  int npsp = struct_list_len(&atoms->pseudopotentials.p);
  if (npsp < 0)
    return -1;
  for (int i = 0; i < npsp; ++i) {
    struct CPMDAtomsPseudopotential psp;
    get_CPMDAtomsPseudopotential(&psp, atoms->pseudopotentials, i);
    /* path is the PP filename (e.g. O_MT_BLYP.psp); element is Z symbol only. */
    if (append_text(dst, dst_size, used, "*") != 0)
      return -1;
    if (psp.path.str && psp.path.len > 0) {
      if (append_capn_text(dst, dst_size, used, psp.path) != 0)
        return -1;
    } else if (psp.element.str && psp.element.len > 0) {
      if (append_capn_text(dst, dst_size, used, psp.element) != 0)
        return -1;
      if (append_text(dst, dst_size, used, ".psp") != 0)
        return -1;
    } else {
      return -1;
    }
    if (append_text(dst, dst_size, used, "\n") != 0)
      return -1;
    if (append_fmt(dst, dst_size, used, " LMAX=%s\n",
                   lmax_letter(psp.lmax >= 0 ? psp.lmax : 0)) != 0)
      return -1;
    /* Coordinate counts filled by geometry merge (see render_deck_with_geometry). */
    if (append_text(dst, dst_size, used, "   0\n") != 0)
      return -1;
  }
  if (append_directives(dst, dst_size, used, atoms->directives) != 0)
    return -1;
  if (append_set_directives_for_section(dst, dst_size, used, sets, "ATOMS") != 0)
    return -1;
  return append_text(dst, dst_size, used, "&END\n\n");
}

static int render_generic_section(char *dst, size_t dst_size, size_t *used,
                                  const struct CPMDGenericSection *sec,
                                  struct RenderSetList *sets) {
  if (append_text(dst, dst_size, used, "&") != 0)
    return -1;
  if (append_capn_text(dst, dst_size, used, sec->name) != 0)
    return -1;
  if (append_text(dst, dst_size, used, "\n") != 0)
    return -1;
  if (append_directives(dst, dst_size, used, sec->directives) != 0)
    return -1;
  if (sec->name.str && sec->name.len > 0) {
    char section[64];
    size_t n = (size_t)sec->name.len;
    if (n >= sizeof(section))
      return -1;
    memcpy(section, sec->name.str, n);
    section[n] = '\0';
    if (append_set_directives_for_section(dst, dst_size, used, sets, section) !=
        0)
      return -1;
  }
  return append_text(dst, dst_size, used, "&END\n\n");
}

int cpmdc_params_root(const void *params_capnp, size_t params_capnp_size_bytes,
                      struct capn *arena, CPMDParams_ptr *params) {
  if (!params_capnp || params_capnp_size_bytes == 0 || !arena || !params)
    return -1;
  memset(arena, 0, sizeof(*arena));
  memset(params, 0, sizeof(*params));
  if (capn_init_mem(arena, (const uint8_t *)params_capnp,
                    params_capnp_size_bytes, 0) != 0)
    return -1;
  params->p = capn_getp(capn_root(arena), 0, 1);
  if (params->p.type != CAPN_STRUCT) {
    cpmdc_params_release(arena);
    memset(params, 0, sizeof(*params));
    return -1;
  }
  return 0;
}

void cpmdc_params_release(struct capn *arena) {
  if (!arena)
    return;
  capn_free(arena);
  memset(arena, 0, sizeof(*arena));
}

const char *cpmdc_params_text_or(capn_text text, const char *fallback) {
  if (text.str && text.len > 0)
    return text.str;
  return fallback ? fallback : "";
}

int cpmdc_params_render_input_deck(CPMDParams_ptr params, char *dst,
                                   size_t dst_size) {
  if (params.p.type == CAPN_NULL || !dst || dst_size == 0)
    return -1;
  dst[0] = '\0';
  size_t used = 0;

  struct CPMDParams view;
  read_CPMDParams(&view, params);

  if (view.title.str && view.title.len > 0) {
    if (append_text(dst, dst_size, &used, "! ") != 0)
      return -1;
    if (append_capn_text(dst, dst_size, &used, view.title) != 0)
      return -1;
    if (append_text(dst, dst_size, &used, "\n\n") != 0)
      return -1;
  }

  int nblocks = pointer_list_len(&view.inputBlocks);
  if (nblocks < 0)
    return -1;
  for (int i = 0; i < nblocks; ++i) {
    capn_text block = capn_get_text(view.inputBlocks, i, empty_text);
    if (append_capn_text(dst, dst_size, &used, block) != 0)
      return -1;
    if (append_text(dst, dst_size, &used, "\n\n") != 0)
      return -1;
  }

  const char *task = cpmdc_params_text_or(view.task, "gradient");
  const char *functional = cpmdc_params_text_or(view.functional, "BLYP");
  capn_text file_path = selected_file_path(&view);
  double cutoff = view.cutOffRy > 0.0 ? view.cutOffRy : 70.0;
  int charge = view.charge;
  int mult = view.multiplicity > 0 ? view.multiplicity : 1;

  int has_system = 0, has_cpmd = 0, has_dft = 0, has_atoms = 0;
  int nsec = struct_list_len(&view.inputSections.p);
  if (nsec < 0)
    return -1;
  struct RenderSetList sets;
  if (collect_set_directives(view.inputSections, &sets) != 0)
    return -1;

  for (int i = 0; i < nsec; ++i) {
    struct CPMDInputSection sec;
    get_CPMDInputSection(&sec, view.inputSections, i);
    switch (sec.which) {
    case CPMDInputSection_generic: {
      struct CPMDGenericSection body;
      read_CPMDGenericSection(&body, sec.generic);
      if (render_generic_section(dst, dst_size, &used, &body, &sets) != 0)
        return -1;
      break;
    }
    case CPMDInputSection_system: {
      struct CPMDSystemSection body;
      read_CPMDSystemSection(&body, sec.system);
      has_system = 1;
      if (render_system_section(dst, dst_size, &used, &body, cutoff, charge,
                                &sets) != 0)
        return -1;
      break;
    }
    case CPMDInputSection_cpmd: {
      struct CPMDCpmdSection body;
      read_CPMDCpmdSection(&body, sec.cpmd);
      has_cpmd = 1;
      if (render_cpmd_section(dst, dst_size, &used, &body, task, file_path,
                              &sets) != 0)
        return -1;
      break;
    }
    case CPMDInputSection_dft: {
      struct CPMDDftSection body;
      read_CPMDDftSection(&body, sec.dft);
      has_dft = 1;
      if (render_dft_section(dst, dst_size, &used, &body, functional, mult,
                             &sets) != 0)
        return -1;
      break;
    }
    case CPMDInputSection_atoms: {
      struct CPMDAtomsSection body;
      read_CPMDAtomsSection(&body, sec.atoms);
      has_atoms = 1;
      if (render_atoms_section(dst, dst_size, &used, &body, &sets) != 0)
        return -1;
      break;
    }
    case CPMDInputSection_set: {
      break;
    }
    case CPMDInputSection_raw:
      if (append_capn_text(dst, dst_size, &used, sec.raw) != 0)
        return -1;
      if (append_text(dst, dst_size, &used, "\n\n") != 0)
        return -1;
      break;
    default:
      return -1;
    }
  }

  if (!has_cpmd) {
    struct CPMDCpmdSection def;
    memset(&def, 0, sizeof(def));
    def.optimizeWavefunction = 1;
    def.convergenceOrbitals = 1.0e-6;
    if (render_cpmd_section(dst, dst_size, &used, &def, task, file_path,
                            &sets) != 0)
      return -1;
  }
  if (!has_system) {
    struct CPMDSystemSection def;
    memset(&def, 0, sizeof(def));
    def.symmetry = 0;
    def.angstrom = 1;
    def.cutOffRy = cutoff;
    def.charge = charge;
    def.multiplicity = mult;
    if (render_system_section(dst, dst_size, &used, &def, cutoff, charge,
                              &sets) != 0)
      return -1;
  }
  if (!has_dft) {
    struct CPMDDftSection def;
    memset(&def, 0, sizeof(def));
    def.lsd = mult > 1;
    if (render_dft_section(dst, dst_size, &used, &def, functional, mult,
                           &sets) != 0)
      return -1;
  }
  if (!has_atoms) {
    struct CPMDAtomsSection def;
    memset(&def, 0, sizeof(def));
    if (render_atoms_section(dst, dst_size, &used, &def, &sets) != 0)
      return -1;
  }
  if (render_remaining_set_sections(dst, dst_size, &used, &sets) != 0)
    return -1;

  return 0;
}

static int element_z(const char *sym, size_t len) {
  if (!sym || len == 0)
    return -1;
  char a = sym[0];
  char b = len > 1 ? sym[1] : '\0';
  if (a >= 'a' && a <= 'z')
    a = (char)(a - 'a' + 'A');
  if (b >= 'a' && b <= 'z')
    b = (char)(b - 'a' + 'A');
  if (a == 'H' && (b == '\0' || b == ' '))
    return 1;
  if (a == 'O' && (b == '\0' || b == ' '))
    return 8;
  if (a == 'C' && (b == '\0' || b == ' '))
    return 6;
  if (a == 'N' && (b == '\0' || b == ' '))
    return 7;
  return -1;
}

static int render_atoms_with_geometry(char *dst, size_t dst_size, size_t *used,
                                      struct CPMDAtomsSection *atoms,
                                      int n_atoms, const double *pos,
                                      const int *z,
                                      struct RenderSetList *sets) {
  if (append_text(dst, dst_size, used, "&ATOMS\n") != 0)
    return -1;
  int npsp = struct_list_len(&atoms->pseudopotentials.p);
  if (npsp < 0)
    return -1;
  if (npsp == 0) {
    /* Default BLYP PPs for H/O when Cap'n Proto omits atoms list. */
    static const struct {
      int z;
      const char *path;
      const char *lmax;
    } defs[] = {{8, "O_MT_BLYP.psp", "P"}, {1, "H_CVB_BLYP.psp", "S"}};
    for (size_t d = 0; d < sizeof(defs) / sizeof(defs[0]); ++d) {
      int count = 0;
      for (int j = 0; j < n_atoms; ++j)
        if (z[j] == defs[d].z)
          ++count;
      if (count == 0)
        continue;
      if (append_fmt(dst, dst_size, used, "*%s\n LMAX=%s\n   %d\n", defs[d].path,
                     defs[d].lmax, count) != 0)
        return -1;
      for (int j = 0; j < n_atoms; ++j) {
        if (z[j] != defs[d].z)
          continue;
        if (append_fmt(dst, dst_size, used, "  %.6f  %.6f  %.6f\n", pos[3 * j],
                       pos[3 * j + 1], pos[3 * j + 2]) != 0)
          return -1;
      }
    }
    if (append_set_directives_for_section(dst, dst_size, used, sets, "ATOMS") !=
        0)
      return -1;
    return append_text(dst, dst_size, used, "&END\n\n");
  }
  for (int i = 0; i < npsp; ++i) {
    struct CPMDAtomsPseudopotential psp;
    get_CPMDAtomsPseudopotential(&psp, atoms->pseudopotentials, i);
    int zz = element_z(psp.element.str, (size_t)psp.element.len);
    if (zz < 0)
      return -1;
    int count = 0;
    for (int j = 0; j < n_atoms; ++j)
      if (z[j] == zz)
        ++count;
    if (append_text(dst, dst_size, used, "*") != 0)
      return -1;
    if (psp.path.str && psp.path.len > 0) {
      if (append_capn_text(dst, dst_size, used, psp.path) != 0)
        return -1;
    } else {
      return -1;
    }
    if (append_text(dst, dst_size, used, "\n") != 0)
      return -1;
    if (append_fmt(dst, dst_size, used, " LMAX=%s\n   %d\n",
                   lmax_letter(psp.lmax >= 0 ? psp.lmax : (zz == 8 ? 1 : 0)),
                   count) != 0)
      return -1;
    for (int j = 0; j < n_atoms; ++j) {
      if (z[j] != zz)
        continue;
      if (append_fmt(dst, dst_size, used, "  %.6f  %.6f  %.6f\n", pos[3 * j],
                     pos[3 * j + 1], pos[3 * j + 2]) != 0)
        return -1;
    }
  }
  if (append_directives(dst, dst_size, used, atoms->directives) != 0)
    return -1;
  if (append_set_directives_for_section(dst, dst_size, used, sets, "ATOMS") != 0)
    return -1;
  return append_text(dst, dst_size, used, "&END\n\n");
}

int cpmdc_params_render_deck_with_geometry(
    CPMDParams_ptr params, int n_atoms, const double *positions_ang,
    const int *atomic_numbers, const double *cell_ang, int has_cell, char *dst,
    size_t dst_size) {
  if (params.p.type == CAPN_NULL || !dst || dst_size == 0 || n_atoms <= 0 ||
      !positions_ang || !atomic_numbers)
    return -1;
  dst[0] = '\0';
  size_t used = 0;
  struct CPMDParams view;
  read_CPMDParams(&view, params);

  if (view.title.str && view.title.len > 0) {
    if (append_text(dst, dst_size, &used, "! ") != 0)
      return -1;
    if (append_capn_text(dst, dst_size, &used, view.title) != 0)
      return -1;
    if (append_text(dst, dst_size, &used, "\n\n") != 0)
      return -1;
  }

  int nblocks = pointer_list_len(&view.inputBlocks);
  if (nblocks < 0)
    return -1;
  for (int i = 0; i < nblocks; ++i) {
    capn_text block = capn_get_text(view.inputBlocks, i, empty_text);
    if (append_capn_text(dst, dst_size, &used, block) != 0)
      return -1;
    if (append_text(dst, dst_size, &used, "\n\n") != 0)
      return -1;
  }

  const char *task = cpmdc_params_text_or(view.task, "gradient");
  const char *functional = cpmdc_params_text_or(view.functional, "BLYP");
  capn_text file_path = selected_file_path(&view);
  double cutoff = view.cutOffRy > 0.0 ? view.cutOffRy : 70.0;
  int charge = view.charge;
  int mult = view.multiplicity > 0 ? view.multiplicity : 1;

  int has_system = 0, has_cpmd = 0, has_dft = 0, has_atoms = 0;
  int nsec = struct_list_len(&view.inputSections.p);
  if (nsec < 0)
    return -1;
  struct RenderSetList sets;
  if (collect_set_directives(view.inputSections, &sets) != 0)
    return -1;

  struct CPMDAtomsSection atoms_body;
  memset(&atoms_body, 0, sizeof(atoms_body));

  for (int i = 0; i < nsec; ++i) {
    struct CPMDInputSection sec;
    get_CPMDInputSection(&sec, view.inputSections, i);
    switch (sec.which) {
    case CPMDInputSection_generic: {
      struct CPMDGenericSection body;
      read_CPMDGenericSection(&body, sec.generic);
      if (render_generic_section(dst, dst_size, &used, &body, &sets) != 0)
        return -1;
      break;
    }
    case CPMDInputSection_system: {
      struct CPMDSystemSection body;
      read_CPMDSystemSection(&body, sec.system);
      has_system = 1;
      if (render_system_section_with_cell(
              dst, dst_size, &used, &body, cutoff, charge,
              has_cell ? cell_ang : NULL, has_cell ? 9 : 0, &sets) != 0)
        return -1;
      break;
    }
    case CPMDInputSection_cpmd: {
      struct CPMDCpmdSection body;
      read_CPMDCpmdSection(&body, sec.cpmd);
      has_cpmd = 1;
      if (render_cpmd_section(dst, dst_size, &used, &body, task, file_path,
                              &sets) != 0)
        return -1;
      break;
    }
    case CPMDInputSection_dft: {
      struct CPMDDftSection body;
      read_CPMDDftSection(&body, sec.dft);
      has_dft = 1;
      if (render_dft_section(dst, dst_size, &used, &body, functional, mult,
                             &sets) != 0)
        return -1;
      break;
    }
    case CPMDInputSection_atoms: {
      read_CPMDAtomsSection(&atoms_body, sec.atoms);
      has_atoms = 1;
      if (render_atoms_with_geometry(dst, dst_size, &used, &atoms_body, n_atoms,
                                     positions_ang, atomic_numbers, &sets) != 0)
        return -1;
      break;
    }
    case CPMDInputSection_set: {
      break;
    }
    case CPMDInputSection_raw:
      if (append_capn_text(dst, dst_size, &used, sec.raw) != 0)
        return -1;
      if (append_text(dst, dst_size, &used, "\n\n") != 0)
        return -1;
      break;
    default:
      return -1;
    }
  }

  if (!has_cpmd) {
    struct CPMDCpmdSection def;
    memset(&def, 0, sizeof(def));
    def.optimizeWavefunction = 1;
    def.convergenceOrbitals = 1.0e-5;
    if (render_cpmd_section(dst, dst_size, &used, &def, task, file_path,
                            &sets) != 0)
      return -1;
  }
  if (!has_system) {
    struct CPMDSystemSection def;
    memset(&def, 0, sizeof(def));
    def.symmetry = 0;
    def.angstrom = 1;
    def.cutOffRy = cutoff;
    def.charge = charge;
    def.multiplicity = mult;
    if (render_system_section_with_cell(dst, dst_size, &used, &def, cutoff,
                                        charge, has_cell ? cell_ang : NULL,
                                        has_cell ? 9 : 0, &sets) != 0)
      return -1;
  }
  if (!has_dft) {
    struct CPMDDftSection def;
    memset(&def, 0, sizeof(def));
    def.lsd = mult > 1;
    if (render_dft_section(dst, dst_size, &used, &def, functional, mult,
                           &sets) != 0)
      return -1;
  }
  if (!has_atoms) {
    if (render_atoms_with_geometry(dst, dst_size, &used, &atoms_body, n_atoms,
                                   positions_ang, atomic_numbers, &sets) != 0)
      return -1;
  }
  if (render_remaining_set_sections(dst, dst_size, &used, &sets) != 0)
    return -1;
  (void)has_cell;
  return 0;
}

static int force_input_length_factor(capn_text unit, double *factor) {
  if (!factor)
    return -1;
  if (!unit.str || unit.len <= 0 || text_eq_ci(unit, "angstrom") ||
      text_eq_ci(unit, "angstroms") || text_eq_ci(unit, "a") ||
      text_eq_ci(unit, "aa")) {
    *factor = 1.0;
    return 0;
  }
  if (text_eq_ci(unit, "bohr") || text_eq_ci(unit, "au") ||
      text_eq_ci(unit, "a.u.") || text_eq_ci(unit, "atomic")) {
    *factor = CPMDC_BOHR_TO_ANGSTROM;
    return 0;
  }
  if (text_eq_ci(unit, "nm") || text_eq_ci(unit, "nanometer")) {
    *factor = 10.0;
    return 0;
  }
  return -1;
}

static int force_input_energy_factor(capn_text unit, double *factor) {
  if (!factor)
    return -1;
  if (!unit.str || unit.len <= 0 || text_eq_ci(unit, "hartree") ||
      text_eq_ci(unit, "ha") || text_eq_ci(unit, "au") ||
      text_eq_ci(unit, "a.u.")) {
    *factor = 1.0;
    return 0;
  }
  if (text_eq_ci(unit, "ev") || text_eq_ci(unit, "electronvolt")) {
    *factor = 27.211386245988;
    return 0;
  }
  if (text_eq_ci(unit, "ry") || text_eq_ci(unit, "rydberg")) {
    *factor = 2.0;
    return 0;
  }
  if (text_eq_ci(unit, "kj/mol") || text_eq_ci(unit, "kjmol")) {
    *factor = 2625.499639479;
    return 0;
  }
  if (text_eq_ci(unit, "kcal/mol") || text_eq_ci(unit, "kcalmol")) {
    *factor = 627.5094740631;
    return 0;
  }
  return -1;
}

int cpmdc_force_input_root(const void *force_input_capnp,
                           size_t force_input_capnp_size_bytes,
                           struct capn *arena, ForceInput_ptr *force_input) {
  if (!force_input_capnp || force_input_capnp_size_bytes == 0 || !arena ||
      !force_input)
    return -1;
  memset(arena, 0, sizeof(*arena));
  memset(force_input, 0, sizeof(*force_input));
  if (capn_init_mem(arena, (const uint8_t *)force_input_capnp,
                    force_input_capnp_size_bytes, 0) != 0)
    return -1;
  force_input->p = capn_getp(capn_root(arena), 0, 1);
  if (force_input->p.type != CAPN_STRUCT) {
    cpmdc_params_release(arena);
    memset(force_input, 0, sizeof(*force_input));
    return -1;
  }
  return 0;
}

int cpmdc_force_input_atom_count(ForceInput_ptr force_input, size_t *n_atoms,
                                 int *has_cell) {
  if (force_input.p.type == CAPN_NULL || !n_atoms || !has_cell)
    return -1;
  struct ForceInput view;
  read_ForceInput(&view, force_input);
  int n_pos = list64_len(&view.pos);
  int n_z = list32_len(&view.atmnrs);
  int n_box = list64_len(&view.box);
  if (n_pos < 0 || n_z <= 0 || n_box < 0)
    return -1;
  if (n_pos != n_z * 3)
    return -1;
  if (n_box != 0 && n_box != 9)
    return -1;
  *n_atoms = (size_t)n_z;
  *has_cell = n_box == 9 ? 1 : 0;
  return 0;
}

int cpmdc_force_input_copy_geometry(ForceInput_ptr force_input,
                                    double *positions_ang, int *atomic_numbers,
                                    size_t atom_capacity, double *cell_ang,
                                    int *has_cell) {
  if (force_input.p.type == CAPN_NULL || !positions_ang || !atomic_numbers ||
      !cell_ang || !has_cell)
    return -1;
  struct ForceInput view;
  read_ForceInput(&view, force_input);
  size_t n_atoms = 0;
  int local_has_cell = 0;
  if (cpmdc_force_input_atom_count(force_input, &n_atoms, &local_has_cell) != 0)
    return -1;
  if (atom_capacity < n_atoms)
    return -1;
  double length_factor = 1.0;
  if (force_input_length_factor(view.lengthUnit, &length_factor) != 0)
    return -1;
  capn_resolve(&view.pos.p);
  capn_resolve(&view.atmnrs.p);
  capn_resolve(&view.box.p);
  for (size_t i = 0; i < n_atoms; ++i)
    atomic_numbers[i] = (int)(int32_t)capn_get32(view.atmnrs, (int)i);
  for (size_t i = 0; i < n_atoms * 3u; ++i)
    positions_ang[i] =
        capn_to_f64(capn_get64(view.pos, (int)i)) * length_factor;
  for (size_t i = 0; i < 9; ++i)
    cell_ang[i] = local_has_cell
                      ? capn_to_f64(capn_get64(view.box, (int)i)) * length_factor
                      : 0.0;
  *has_cell = local_has_cell;
  return 0;
}

int cpmdc_force_input_result_factors(ForceInput_ptr force_input,
                                     double *energy_factor,
                                     double *force_factor) {
  if (force_input.p.type == CAPN_NULL || !energy_factor || !force_factor)
    return -1;
  struct ForceInput view;
  read_ForceInput(&view, force_input);
  double length_factor = 1.0;
  double energy = 1.0;
  if (force_input_length_factor(view.lengthUnit, &length_factor) != 0 ||
      force_input_energy_factor(view.energyUnit, &energy) != 0)
    return -1;
  *energy_factor = energy;
  *force_factor = energy * length_factor / CPMDC_BOHR_TO_ANGSTROM;
  return 0;
}

size_t cpmdc_potential_result_flat_size(size_t force_count) {
  if (force_count > (SIZE_MAX - 32u) / 8u)
    return 0;
  return 32u + force_count * 8u;
}

int cpmdc_potential_result_write(double energy, const double *forces,
                                 size_t force_count,
                                 void *potential_result_capnp,
                                 size_t potential_result_capacity_bytes,
                                 size_t *potential_result_size_bytes) {
  if (!forces || !potential_result_capnp || !potential_result_size_bytes ||
      force_count > (size_t)INT_MAX)
    return -1;
  size_t required = cpmdc_potential_result_flat_size(force_count);
  *potential_result_size_bytes = required;
  if (required == 0 || potential_result_capacity_bytes < required)
    return -1;

  struct capn arena;
  capn_init_malloc(&arena);
  capn_ptr root = capn_root(&arena);
  if (root.type == CAPN_NULL) {
    capn_free(&arena);
    return -1;
  }
  PotentialResult_ptr result = new_PotentialResult(root.seg);
  capn_list64 force_list = capn_new_list64(root.seg, (int)force_count);
  if (result.p.type == CAPN_NULL ||
      (force_count > 0 && force_list.p.type == CAPN_NULL)) {
    capn_free(&arena);
    return -1;
  }
  for (size_t i = 0; i < force_count; ++i)
    capn_set64(force_list, (int)i, capn_from_f64(forces[i]));

  struct PotentialResult view;
  view.energy = energy;
  view.forces = force_list;
  write_PotentialResult(&view, result);
  if (capn_setp(capn_root(&arena), 0, result.p) != 0) {
    capn_free(&arena);
    return -1;
  }
  int wrote = capn_write_mem(&arena, (uint8_t *)potential_result_capnp,
                             (int)potential_result_capacity_bytes, 0);
  capn_free(&arena);
  if (wrote < 0 || (size_t)wrote > potential_result_capacity_bytes)
    return -1;
  *potential_result_size_bytes = (size_t)wrote;
  return 0;
}
