#define _GNU_SOURCE
#include "cpmdc_restart.h"

#include <errno.h>
#include <fcntl.h>
#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#ifdef _WIN32
#include <io.h>
#include <process.h>
#include <windows.h>
#else
#include <unistd.h>
#endif

/* gfortran's default maximum subrecord (libgfortran GFC_MAX_SUBRECORD_LENGTH). */
#define CPMDC_RESTART_MAX_SUB 2147483639u
#define CPMDC_RESTART_MAX_SECTIONS 99

struct rec {
  unsigned char *data;
  size_t len;
};

struct cpmdc_restart {
  struct rec *recs;
  size_t nrecs;
  int stream;
  char header[96];
  int sec_start[CPMDC_RESTART_MAX_SECTIONS + 1];
  int sec_nrec[CPMDC_RESTART_MAX_SECTIONS + 1];
  int *na;
  int nsp;
  int has_cell;
  int ibrav;
  int indpg;
  double celldm[6];
  int has_cut;
  double ecut;
  double cdual;
  int dual_flag;
  int nel;
  int nr1s;
  int nr2s;
  int nr3s;
  int has_states;
  int nstates;
  int nkpts;
  int ngw;
  int ngwl;
  int nhg;
  int nhgl;
};

static void set_err(char *err, size_t cap, const char *msg) {
  if (!err || cap == 0)
    return;
  snprintf(err, cap, "%s", msg);
}

/* INT30-C: off + width wraps when off is near SIZE_MAX, and the check then passes. */
static int field_fits(const struct rec *r, size_t off, size_t width) {
  if (!r || off > r->len || r->len - off < width)
    return 0;
  return 1;
}

static int read_i32(const struct rec *r, size_t off, int32_t *out) {
  if (!field_fits(r, off, 4))
    return -1;
  memcpy(out, r->data + off, 4);
  return 0;
}

static int read_f64(const struct rec *r, size_t off, double *out) {
  if (!field_fits(r, off, 8))
    return -1;
  memcpy(out, r->data + off, 8);
  return 0;
}

static int push_rec(struct rec **recs, size_t *n, size_t *cap, unsigned char *data,
                    size_t len) {
  if (*n == *cap) {
    size_t ncap;
    /* INT30-C: a wrapped cap * 2 or cap * sizeof is a short realloc. */
    if (*cap == 0)
      ncap = 64;
    else if (*cap > SIZE_MAX / 2)
      ncap = 0;
    else
      ncap = *cap * 2;
    if (ncap == 0 || ncap > SIZE_MAX / sizeof(*recs[0])) {
      free(data);
      return -1;
    }
    struct rec *grown = realloc(*recs, ncap * sizeof(*recs[0]));
    if (!grown) {
      free(data);
      return -1;
    }
    *recs = grown;
    *cap = ncap;
  }
  (*recs)[*n].data = data;
  (*recs)[*n].len = len;
  (*n)++;
  return 0;
}

static int append_bytes(unsigned char **buf, size_t *len, size_t *cap,
                        const unsigned char *src, size_t n) {
  if (*len + n < *len)
    return -1;
  if (*len + n > *cap) {
    size_t ncap = *cap ? *cap : 256;
    while (ncap < *len + n) {
      if (ncap > (SIZE_MAX / 2))
        return -1;
      ncap *= 2;
    }
    unsigned char *grown = realloc(*buf, ncap);
    if (!grown)
      return -1;
    *buf = grown;
    *cap = ncap;
  }
  if (n)
    memcpy(*buf + *len, src, n);
  *len += n;
  return 0;
}

static int parse_records(const unsigned char *buf, size_t n, struct rec **out,
                         size_t *nout, char *err, size_t err_cap) {
  size_t off = 0;
  struct rec *recs = NULL;
  size_t nrec = 0, cap = 0;
  *out = NULL;
  *nout = 0;
  while (off < n) {
    unsigned char *payload = NULL;
    size_t plen = 0, pcap = 0;
    if (n - off < 4) {
      set_err(err, err_cap, "truncated record marker");
      goto fail;
    }
    for (;;) {
      int32_t lead;
      uint32_t slen;
      int cont;
      if (n - off < 4) {
        set_err(err, err_cap, "truncated subrecord marker");
        free(payload);
        goto fail;
      }
      memcpy(&lead, buf + off, 4);
      off += 4;
      cont = lead < 0;
      if (lead == INT32_MIN) {
        set_err(err, err_cap, "record marker out of range");
        free(payload);
        goto fail;
      }
      slen = cont ? (uint32_t)(-lead) : (uint32_t)lead;
      if ((size_t)slen > n - off) {
        set_err(err, err_cap, "record extends past end of file");
        free(payload);
        goto fail;
      }
      if (append_bytes(&payload, &plen, &pcap, buf + off, slen) != 0) {
        set_err(err, err_cap, "out of memory");
        free(payload);
        goto fail;
      }
      off += slen;
      if (n - off < 4) {
        set_err(err, err_cap, "truncated trailing record marker");
        free(payload);
        goto fail;
      }
      {
        int32_t trail;
        memcpy(&trail, buf + off, 4);
        off += 4;
        if (trail != lead) {
          set_err(err, err_cap, "leading and trailing record markers differ");
          free(payload);
          goto fail;
        }
      }
      if (!cont)
        break;
    }
    if (push_rec(&recs, &nrec, &cap, payload, plen) != 0) {
      set_err(err, err_cap, "out of memory");
      goto fail;
    }
  }
  *out = recs;
  *nout = nrec;
  return 0;
fail:
  for (size_t i = 0; i < nrec; i++)
    free(recs[i].data);
  free(recs);
  return -1;
}

static void free_recs(struct cpmdc_restart *file) {
  if (!file)
    return;
  for (size_t i = 0; i < file->nrecs; i++)
    free(file->recs[i].data);
  free(file->recs);
  file->recs = NULL;
  file->nrecs = 0;
}

static int section_count_of(const struct rec *r, int32_t *count) {
  if (read_i32(r, 0, count) != 0)
    return -1;
  return 0;
}

static int index_sections(CPMDCRestart *file, char *err, size_t err_cap) {
  size_t i;
  const unsigned char *h;
  size_t hlen;
  if (file->nrecs < 1) {
    set_err(err, err_cap, "restart file has no records");
    return -1;
  }
  h = file->recs[0].data;
  hlen = file->recs[0].len;
  if (!memmem(h, hlen, "STUTTGART VERSION 3.0", 21)) {
    set_err(err, err_cap, "section 1 is not a STUTTGART VERSION 3.0 header");
    return -1;
  }
  if (hlen >= sizeof(file->header))
    hlen = sizeof(file->header) - 1;
  memcpy(file->header, h, hlen);
  file->header[hlen] = '\0';
  while (hlen > 0 && (file->header[hlen - 1] == ' ' || file->header[hlen - 1] == '\0')) {
    file->header[hlen - 1] = '\0';
    hlen--;
  }
  file->stream = strstr(file->header, "STREAM") != NULL;
  for (int s = 0; s <= CPMDC_RESTART_MAX_SECTIONS; s++) {
    file->sec_start[s] = -1;
    file->sec_nrec[s] = 0;
  }
  file->sec_start[1] = 0;
  file->sec_nrec[1] = 1;
  i = 1;
  for (int s = 2; s <= CPMDC_RESTART_MAX_SECTIONS && i < file->nrecs; s++) {
    int32_t count;
    int64_t nfollow;
    if (section_count_of(&file->recs[i], &count) != 0) {
      set_err(err, err_cap, "section header is shorter than one integer");
      return -1;
    }
    if (count == INT32_MIN) {
      set_err(err, err_cap, "section record count out of range");
      return -1;
    }
    nfollow = count < 0 ? -(int64_t)count : (int64_t)count;
    if (nfollow > (int64_t)(file->nrecs - i - 1)) {
      set_err(err, err_cap, "section asks for records past the end of the file");
      return -1;
    }
    /* INT31-C: sec_start and sec_nrec are int. A truncated index reads another record. */
    if (i > (size_t)INT_MAX || nfollow > (int64_t)INT_MAX - 1) {
      set_err(err, err_cap, "section index does not fit in int");
      return -1;
    }
    file->sec_start[s] = (int)i;
    file->sec_nrec[s] = (int)(1 + nfollow);
    i += (size_t)(1 + nfollow);
  }
  if (i != file->nrecs) {
    set_err(err, err_cap, "records remain after section 99");
    return -1;
  }
  return 0;
}

/* INT32-C: sec_start + offset is formed in size_t. A signed sum near INT_MAX is a different record. */
static struct rec *rec_at(CPMDCRestart *file, int section, int offset) {
  size_t idx;
  if (!file || section < 1 || section > CPMDC_RESTART_MAX_SECTIONS || offset < 0)
    return NULL;
  if (file->sec_start[section] < 0)
    return NULL;
  idx = (size_t)file->sec_start[section] + (size_t)offset;
  if (idx < (size_t)offset || idx >= file->nrecs)
    return NULL;
  return &file->recs[idx];
}

static int decode_typed(CPMDCRestart *file, char *err, size_t err_cap) {
  if (file->sec_nrec[2] == 3) {
    const struct rec *sym = rec_at(file, 2, 1);
    const struct rec *cell = rec_at(file, 2, 2);
    if (!sym || !cell) {
      set_err(err, err_cap, "section 2 symmetry or cell record is missing");
      return -1;
    }
    int32_t ibrav, indpg;
    if (sym->len < 8 || cell->len < 48) {
      set_err(err, err_cap, "section 2 symmetry or cell record has the wrong length");
      return -1;
    }
    if (read_i32(sym, 0, &ibrav) != 0 || read_i32(sym, 4, &indpg) != 0)
      return -1;
    file->ibrav = ibrav;
    file->indpg = indpg;
    for (int k = 0; k < 6; k++) {
      if (read_f64(cell, (size_t)k * 8, &file->celldm[k]) != 0)
        return -1;
    }
    file->has_cell = 1;
  }
  if (file->sec_nrec[3] == 3) {
    const struct rec *nsp_r = rec_at(file, 3, 1);
    const struct rec *na_r = rec_at(file, 3, 2);
    if (!nsp_r || !na_r) {
      set_err(err, err_cap, "section 3 species record is missing");
      return -1;
    }
    int32_t nsp;
    if (read_i32(nsp_r, 0, &nsp) != 0 || nsp < 0 || nsp > 100000) {
      set_err(err, err_cap, "section 3 species count is unusable");
      return -1;
    }
    if (na_r->len < (size_t)nsp * 4) {
      set_err(err, err_cap, "section 3 atom-count record is short");
      return -1;
    }
    file->na = calloc((size_t)nsp > 0 ? (size_t)nsp : 1, sizeof(int));
    if (!file->na) {
      set_err(err, err_cap, "out of memory");
      return -1;
    }
    file->nsp = nsp;
    for (int s = 0; s < nsp; s++) {
      int32_t na;
      if (read_i32(na_r, (size_t)s * 4, &na) != 0 || na < 0) {
        set_err(err, err_cap, "section 3 atom count is unusable");
        return -1;
      }
      file->na[s] = na;
    }
  }
  if (file->sec_nrec[7] == 2) {
    const struct rec *r = rec_at(file, 7, 1);
    if (!r) {
      set_err(err, err_cap, "section 7 cutoff record is missing");
      return -1;
    }
    int32_t dual_flag, nel, nr1, nr2, nr3;
    /* ecut, cdual (real*8), dual (logical*4), nel, nr1s, nr2s, nr3s. */
    if (r->len < 36) {
      set_err(err, err_cap, "section 7 cutoff record is short");
      return -1;
    }
    if (read_f64(r, 0, &file->ecut) != 0 || read_f64(r, 8, &file->cdual) != 0 ||
        read_i32(r, 16, &dual_flag) != 0 || read_i32(r, 20, &nel) != 0 ||
        read_i32(r, 24, &nr1) != 0 || read_i32(r, 28, &nr2) != 0 ||
        read_i32(r, 32, &nr3) != 0) {
      set_err(err, err_cap, "section 7 cutoff record is short");
      return -1;
    }
    file->dual_flag = dual_flag;
    file->nel = nel;
    file->nr1s = nr1;
    file->nr2s = nr2;
    file->nr3s = nr3;
    file->has_cut = 1;
  }
  if (file->sec_nrec[8] == 2) {
    const struct rec *r = rec_at(file, 8, 1);
    if (!r) {
      set_err(err, err_cap, "section 8 state record is missing");
      return -1;
    }
    int32_t n, nk, ngw, ngwl, nhg, nhgl;
    if (r->len < 24) {
      set_err(err, err_cap, "section 8 state record is short");
      return -1;
    }
    if (read_i32(r, 0, &n) != 0 || read_i32(r, 4, &nk) != 0 ||
        read_i32(r, 8, &ngw) != 0 || read_i32(r, 12, &ngwl) != 0 ||
        read_i32(r, 16, &nhg) != 0 || read_i32(r, 20, &nhgl) != 0) {
      set_err(err, err_cap, "section 8 state record is short");
      return -1;
    }
    file->nstates = n;
    file->nkpts = nk;
    file->ngw = ngw;
    file->ngwl = ngwl;
    file->nhg = nhg;
    file->nhgl = nhgl;
    file->has_states = 1;
  }
  return 0;
}

static CPMDCRestart *from_bytes(const void *bytes, size_t nbytes, char *err,
                                size_t err_cap) {
  CPMDCRestart *file;
  if (!bytes && nbytes) {
    set_err(err, err_cap, "null restart buffer");
    return NULL;
  }
  file = calloc(1, sizeof(*file));
  if (!file) {
    set_err(err, err_cap, "out of memory");
    return NULL;
  }
  if (parse_records(bytes, nbytes, &file->recs, &file->nrecs, err, err_cap) != 0) {
    free(file);
    return NULL;
  }
  if (index_sections(file, err, err_cap) != 0 || decode_typed(file, err, err_cap) != 0) {
    free(file->na);
    free_recs(file);
    free(file);
    return NULL;
  }
  return file;
}

CPMDCRestart *cpmdc_restart_read_mem(const void *bytes, size_t nbytes, char *err,
                                     size_t err_cap) {
  return from_bytes(bytes, nbytes, err, err_cap);
}

CPMDCRestart *cpmdc_restart_read_path(const char *path, char *err, size_t err_cap) {
  FILE *fp;
  unsigned char *buf = NULL;
  size_t n = 0, cap = 0;
  CPMDCRestart *file;
  unsigned char chunk[1 << 16];
  if (!path) {
    set_err(err, err_cap, "null path");
    return NULL;
  }
  fp = fopen(path, "rb");
  if (!fp) {
    set_err(err, err_cap, strerror(errno));
    return NULL;
  }
  for (;;) {
    size_t got = fread(chunk, 1, sizeof chunk, fp);
    if (got && append_bytes(&buf, &n, &cap, chunk, got) != 0) {
      set_err(err, err_cap, "out of memory");
      free(buf);
      fclose(fp);
      return NULL;
    }
    if (got < sizeof chunk) {
      if (ferror(fp)) {
        set_err(err, err_cap, "read failed");
        free(buf);
        fclose(fp);
        return NULL;
      }
      break;
    }
  }
  fclose(fp);
  file = from_bytes(buf, n, err, err_cap);
  free(buf);
  return file;
}

void cpmdc_restart_free(CPMDCRestart *file) {
  if (!file)
    return;
  free(file->na);
  free_recs(file);
  free(file);
}

int cpmdc_restart_is_stream(const CPMDCRestart *file) {
  return file && file->stream;
}

const char *cpmdc_restart_header(const CPMDCRestart *file) {
  return file ? file->header : NULL;
}

int cpmdc_restart_section(const CPMDCRestart *file, int section, int *count) {
  int32_t raw;
  if (!file || section < 1 || section > CPMDC_RESTART_MAX_SECTIONS ||
      file->sec_start[section] < 0)
    return -1;
  if (section == 1) {
    if (count)
      *count = 1;
    return 0;
  }
  if (section_count_of(&file->recs[file->sec_start[section]], &raw) != 0)
    return -1;
  if (count)
    *count = (int)raw;
  return 0;
}

int cpmdc_restart_cell(const CPMDCRestart *file, int *ibrav, int *indpg,
                       double celldm[6]) {
  if (!file || !file->has_cell || !celldm)
    return -1;
  if (ibrav)
    *ibrav = file->ibrav;
  if (indpg)
    *indpg = file->indpg;
  memcpy(celldm, file->celldm, sizeof(file->celldm));
  return 0;
}

int cpmdc_restart_species(const CPMDCRestart *file, int *nsp,
                          const int **na_per_species) {
  if (!file || file->nsp <= 0 || !file->na)
    return -1;
  if (nsp)
    *nsp = file->nsp;
  if (na_per_species)
    *na_per_species = file->na;
  return 0;
}

int cpmdc_restart_expand_species_symbols(const CPMDCRestart *file,
                                         const char *const *species_symbols,
                                         int nspecies, const char **atom_symbols,
                                         int n_atoms) {
  const int *na = NULL;
  int nsp = 0;
  int filled = 0;
  if (!file || !species_symbols || !atom_symbols || n_atoms < 0)
    return -1;
  if (cpmdc_restart_species(file, &nsp, &na) != 0 || nspecies != nsp)
    return -1;
  for (int s = 0; s < nsp; s++) {
    if (!species_symbols[s] || species_symbols[s][0] == '\0')
      return -1;
    if (na[s] < 0)
      return -1;
    for (int a = 0; a < na[s]; a++) {
      if (filled >= n_atoms)
        return -1;
      atom_symbols[filled++] = species_symbols[s];
    }
  }
  if (filled != n_atoms)
    return -1;
  return 0;
}

int cpmdc_restart_ncoords(const CPMDCRestart *file) {
  if (!file || file->sec_nrec[4] < 1)
    return 0;
  return file->sec_nrec[4] - 1;
}

/* INT32-C: n*3 and sec_start+a are signed. A section that does not fit is another record. */
static int triples_span(int sec_nrec, int sec_start, int *n_out, int *start_out) {
  int n;
  if (sec_nrec < 1)
    return -1;
  n = sec_nrec - 1;
  if (n > INT_MAX / 3)
    return -1;
  if (n > 0 && (sec_start < 0 || sec_start > INT_MAX - n))
    return -1;
  if (n_out)
    *n_out = n;
  if (start_out)
    *start_out = n > 0 ? sec_start + 1 : 0;
  return 0;
}

static int copy_triples(const CPMDCRestart *file, int section, double *xyz, int n3) {
  int n;
  int start;
  if (!file || !xyz)
    return -1;
  if (triples_span(file->sec_nrec[section], file->sec_start[section], &n, &start) != 0)
    return -1;
  if (n3 != n * 3)
    return -1;
  for (int a = 0; a < n; a++) {
    const struct rec *r = &file->recs[start + a];
    if (r->len != 24)
      return -1;
    for (int k = 0; k < 3; k++) {
      if (read_f64(r, (size_t)k * 8, &xyz[3 * a + k]) != 0)
        return -1;
    }
  }
  return 0;
}

static int store_triples(CPMDCRestart *file, int section, const double *xyz, int n3) {
  int n;
  int start;
  if (!file || !xyz)
    return -1;
  if (triples_span(file->sec_nrec[section], file->sec_start[section], &n, &start) != 0)
    return -1;
  if (n3 != n * 3)
    return -1;
  for (int a = 0; a < n; a++) {
    struct rec *r = &file->recs[start + a];
    if (r->len != 24)
      return -1;
    memcpy(r->data, xyz + 3 * a, 24);
  }
  return 0;
}

int cpmdc_restart_coordinates(const CPMDCRestart *file, double *xyz, int n3) {
  return copy_triples(file, 4, xyz, n3);
}

int cpmdc_restart_velocities(const CPMDCRestart *file, double *xyz, int n3) {
  return copy_triples(file, 5, xyz, n3);
}

int cpmdc_restart_initial_coordinates(const CPMDCRestart *file, double *xyz, int n3) {
  return copy_triples(file, 6, xyz, n3);
}

int cpmdc_restart_set_coordinates(CPMDCRestart *file, const double *xyz, int n3) {
  return store_triples(file, 4, xyz, n3);
}

int cpmdc_restart_set_velocities(CPMDCRestart *file, const double *xyz, int n3) {
  return store_triples(file, 5, xyz, n3);
}

int cpmdc_restart_set_cell(CPMDCRestart *file, const double celldm[6]) {
  struct rec *r;
  if (!file || !file->has_cell || !celldm || file->sec_nrec[2] != 3)
    return -1;
  r = rec_at(file, 2, 2);
  if (!r || r->len < 48)
    return -1;
  memcpy(r->data, celldm, 48);
  memcpy(file->celldm, celldm, sizeof(file->celldm));
  return 0;
}

int cpmdc_restart_cutoff(const CPMDCRestart *file, double *ecut, double *cdual,
                         int *dual_flag, int *nel, int *nr1s, int *nr2s, int *nr3s) {
  if (!file || !file->has_cut)
    return -1;
  if (ecut)
    *ecut = file->ecut;
  if (cdual)
    *cdual = file->cdual;
  if (dual_flag)
    *dual_flag = file->dual_flag;
  if (nel)
    *nel = file->nel;
  if (nr1s)
    *nr1s = file->nr1s;
  if (nr2s)
    *nr2s = file->nr2s;
  if (nr3s)
    *nr3s = file->nr3s;
  return 0;
}

int cpmdc_restart_states(const CPMDCRestart *file, int *n, int *nkpts, int *ngw,
                         int *ngwl, int *nhg, int *nhgl) {
  if (!file || !file->has_states)
    return -1;
  if (n)
    *n = file->nstates;
  if (nkpts)
    *nkpts = file->nkpts;
  if (ngw)
    *ngw = file->ngw;
  if (ngwl)
    *ngwl = file->ngwl;
  if (nhg)
    *nhg = file->nhg;
  if (nhgl)
    *nhgl = file->nhgl;
  return 0;
}

static int emit_record(unsigned char **buf, size_t *len, size_t *cap,
                       const unsigned char *data, size_t dlen) {
  size_t off = 0;
  if (dlen == 0) {
    int32_t z = 0;
    if (append_bytes(buf, len, cap, (unsigned char *)&z, 4) != 0)
      return -1;
    if (append_bytes(buf, len, cap, (unsigned char *)&z, 4) != 0)
      return -1;
    return 0;
  }
  while (off < dlen) {
    size_t chunk = dlen - off;
    int32_t marker;
    int last;
    if (chunk > CPMDC_RESTART_MAX_SUB)
      chunk = CPMDC_RESTART_MAX_SUB;
    /* INT30-C: off + chunk wraps when off is near SIZE_MAX. */
    last = (chunk == dlen - off);
    if (chunk > (size_t)INT32_MAX)
      return -1;
    marker = last ? (int32_t)chunk : -(int32_t)chunk;
    if (append_bytes(buf, len, cap, (unsigned char *)&marker, 4) != 0)
      return -1;
    if (append_bytes(buf, len, cap, data + off, chunk) != 0)
      return -1;
    if (append_bytes(buf, len, cap, (unsigned char *)&marker, 4) != 0)
      return -1;
    off += chunk;
  }
  return 0;
}

int cpmdc_restart_write_mem(const CPMDCRestart *file, void **bytes, size_t *nbytes,
                            char *err, size_t err_cap) {
  unsigned char *buf = NULL;
  size_t len = 0, cap = 0;
  if (!file || !bytes || !nbytes) {
    set_err(err, err_cap, "null write target");
    return -1;
  }
  *bytes = NULL;
  *nbytes = 0;
  for (size_t i = 0; i < file->nrecs; i++) {
    if (emit_record(&buf, &len, &cap, file->recs[i].data, file->recs[i].len) != 0) {
      set_err(err, err_cap, "out of memory");
      free(buf);
      return -1;
    }
  }
  *bytes = buf;
  *nbytes = len;
  return 0;
}

static void set_errf(char *err, size_t cap, const char *what, const char *why) {
  if (!err || cap == 0)
    return;
  snprintf(err, cap, "%s: %s", what, why);
}

#ifdef _WIN32
/* _O_NOINHERIT is the MSVC CRT counterpart of O_CLOEXEC on the POSIX openat
 * path, so a temporary is not inherited by child processes. */
#define publish_open(p)                                                                            \
  _open((p), _O_WRONLY | _O_CREAT | _O_EXCL | _O_BINARY | _O_NOINHERIT, _S_IREAD | _S_IWRITE)
#define publish_fdopen _fdopen
#define publish_close _close
#define publish_unlink _unlink
#define publish_strdup _strdup
#define publish_pid() ((unsigned long)_getpid())
#define CPMDC_IS_SEP(c) ((c) == '/' || (c) == '\\')
#else
#define publish_fdopen fdopen
#define publish_close close
#define publish_strdup strdup
#define publish_pid() ((unsigned long)getpid())
#define CPMDC_IS_SEP(c) ((c) == '/')
#endif

#ifdef _WIN32
/* Keep the temporary basename short even when the destination basename is at
 * NAME_MAX. The directory prefix is copied from the effective destination. */
static char *temporary_path(const char *path, unsigned long token) {
  static const char stem[] = ".cpmdc-tmp-";
  size_t dir = strlen(path);
  size_t cap;
  char *tmp;
  while (dir > 0 && !CPMDC_IS_SEP(path[dir - 1]))
    dir--;
  cap = dir + sizeof stem + 16;
  tmp = malloc(cap);
  if (!tmp)
    return NULL;
  if (dir)
    memcpy(tmp, path, dir);
  snprintf(tmp + dir, cap - dir, "%s%016lx", stem, token);
  return tmp;
}
#endif

#ifndef _WIN32
/* realpath resolves a chain in one call only when every hop exists. A
 * dangling chain makes it fail with ENOENT. The walk below is bounded so
 * a loop this code reaches on its own still ends in a clear error. */
#define CPMDC_RESTART_MAX_HOPS 40

/* One symlink hop recorded while resolving. dirfd is the directory that
 * contained the link at that moment, held open so a later check does not
 * look the path up again. */
struct publish_hop {
  int dirfd;
  char *name;
  char *text;
  dev_t dev;
  ino_t ino;
};

/* What publication resolved, and the directory inode the rename must use.
 * absent means the final component did not exist. mismatch means that
 * picture already disagreed with itself (realpath said the chain was
 * missing, but a walk found something else, or the reverse) and the
 * rename must refuse. */
struct publish_dest {
  int dirfd;
  char *final_name;
  int absent;
  int mismatch;
  int has_mode;
  mode_t mode;
  dev_t dev;
  ino_t ino;
  struct publish_hop hops[CPMDC_RESTART_MAX_HOPS];
  int nhops;
};

#ifdef CPMDC_RESTART_TEST_HOOK
/* Test seam, compiled out of the installed library. The restart tests set
 * this to run after realpath has reported a dangling chain and before the
 * fallback walk, which is the window before the publication rename. */
void (*cpmdc_restart_test_before_follow)(void);
void (*cpmdc_restart_test_after_temp_open)(int fd);
#endif

static void dest_init(struct publish_dest *d) {
  memset(d, 0, sizeof *d);
  d->dirfd = -1;
  for (int i = 0; i < CPMDC_RESTART_MAX_HOPS; i++)
    d->hops[i].dirfd = -1;
}

static void dest_clear(struct publish_dest *d) {
  if (d->dirfd >= 0)
    close(d->dirfd);
  free(d->final_name);
  for (int i = 0; i < d->nhops; i++) {
    if (d->hops[i].dirfd >= 0)
      close(d->hops[i].dirfd);
    free(d->hops[i].name);
    free(d->hops[i].text);
  }
  dest_init(d);
}

/* readlinkat(2) does not NUL-terminate, and the target can outgrow the
 * first buffer. Grow until the read fits. */
static char *read_linkat_text(int dirfd, const char *name, char *err, size_t err_cap) {
  size_t cap = 128;
  char *buf = NULL;
  for (;;) {
    char *grown = realloc(buf, cap);
    ssize_t n;
    if (!grown) {
      free(buf);
      set_err(err, err_cap, "out of memory");
      return NULL;
    }
    buf = grown;
    n = readlinkat(dirfd, name, buf, cap);
    if (n < 0) {
      set_errf(err, err_cap, "cannot read the RESTART symlink", strerror(errno));
      free(buf);
      return NULL;
    }
    if ((size_t)n < cap) {
      buf[n] = '\0';
      return buf;
    }
    if (cap > SIZE_MAX / 2) {
      set_err(err, err_cap, "the RESTART symlink target is too long");
      free(buf);
      return NULL;
    }
    cap *= 2;
  }
}

/* fopen follows a relative link from the directory that contains it. An
 * absolute target is already a path. */
static char *link_target_path(const char *path, const char *target, char *err,
                              size_t err_cap) {
  size_t dir, tlen;
  char *out;
  if (target[0] == '\0') {
    set_err(err, err_cap, "the RESTART symlink has an empty target");
    return NULL;
  }
  if (target[0] == '/') {
    out = publish_strdup(target);
    if (!out)
      set_err(err, err_cap, "out of memory");
    return out;
  }
  dir = strlen(path);
  while (dir > 0 && !CPMDC_IS_SEP(path[dir - 1]))
    dir--;
  tlen = strlen(target);
  if (tlen > SIZE_MAX - dir - 1) {
    set_err(err, err_cap, "the RESTART symlink target is too long");
    return NULL;
  }
  out = malloc(dir + tlen + 1);
  if (!out) {
    set_err(err, err_cap, "out of memory");
    return NULL;
  }
  if (dir)
    memcpy(out, path, dir);
  memcpy(out + dir, target, tlen + 1);
  return out;
}

/* Split `path` into a directory that can be opened and a final component.
 * A relative path with no slash is in the current directory. */
static int split_parent(const char *path, char **dir_out, char **name_out, char *err,
                        size_t err_cap) {
  size_t n, slash, dlen;
  char *dir, *name;
  *dir_out = NULL;
  *name_out = NULL;
  n = strlen(path);
  slash = n;
  while (slash > 0 && !CPMDC_IS_SEP(path[slash - 1]))
    slash--;
  if (slash == n) {
    set_err(err, err_cap, "the RESTART path has no final component");
    return -1;
  }
  if (slash == 0) {
    dir = publish_strdup(".");
  } else {
    dlen = slash;
    if (dlen > 1 && CPMDC_IS_SEP(path[dlen - 1]))
      dlen--;
    dir = malloc(dlen + 1);
    if (dir) {
      memcpy(dir, path, dlen);
      dir[dlen] = '\0';
    }
  }
  name = publish_strdup(path + slash);
  if (!dir || !name) {
    free(dir);
    free(name);
    set_err(err, err_cap, "out of memory");
    return -1;
  }
  *dir_out = dir;
  *name_out = name;
  return 0;
}

/* Open each directory component without following a symlink. The returned
 * descriptor pins the final directory inode, while rejecting directory
 * symlinks that are not represented in the recorded RESTART chain. */
static int open_directory_nofollow(const char *path) {
  const int flags = O_RDONLY | O_DIRECTORY | O_CLOEXEC | O_NOFOLLOW;
  char *copy = publish_strdup(path);
  char *part;
  int fd;
  if (!copy) {
    errno = ENOMEM;
    return -1;
  }
  fd = open(path[0] == '/' ? "/" : ".", flags);
  if (fd < 0) {
    free(copy);
    return -1;
  }
  part = copy;
  while (CPMDC_IS_SEP(*part))
    part++;
  while (*part) {
    char *end = part;
    char saved;
    int next;
    while (*end && !CPMDC_IS_SEP(*end))
      end++;
    saved = *end;
    *end = '\0';
    next = openat(fd, part, flags);
    if (next < 0) {
      int why = errno;
      close(fd);
      free(copy);
      errno = why;
      return -1;
    }
    close(fd);
    fd = next;
    if (!saved)
      break;
    part = end + 1;
    while (CPMDC_IS_SEP(*part))
      part++;
  }
  free(copy);
  return fd;
}

/* Walk `path` one symlink at a time. `expect_absent` is 1 when realpath
 * already reported the chain missing: an existing non symlink is then a
 * chain that changed, not a file this call may replace. The directory
 * descriptor of the destination, and one for every symlink followed, stay
 * open in `d` until publication finishes. */
static int walk_chain(const char *path, int expect_absent, struct publish_dest *d,
                      char *err, size_t err_cap) {
  char *cur = publish_strdup(path);
  if (!cur) {
    set_err(err, err_cap, "out of memory");
    return -1;
  }
  for (int hop = 0; hop < CPMDC_RESTART_MAX_HOPS; hop++) {
    char *dir = NULL, *name = NULL;
    struct stat st;
    int fd;
    if (split_parent(cur, &dir, &name, err, err_cap) != 0) {
      free(cur);
      return -1;
    }
    fd = open_directory_nofollow(dir);
    free(dir);
    if (fd < 0) {
      set_errf(err, err_cap, "cannot open the RESTART directory", strerror(errno));
      free(name);
      free(cur);
      return -1;
    }
    if (fstatat(fd, name, &st, AT_SYMLINK_NOFOLLOW) != 0) {
      if (errno != ENOENT) {
        set_errf(err, err_cap, "cannot stat the RESTART symlink target", strerror(errno));
        close(fd);
        free(name);
        free(cur);
        return -1;
      }
      /* The final component is missing. That is the dangling end when
       * realpath said so, and a chain that lost its file otherwise. */
      if (!expect_absent)
        d->mismatch = 1;
      d->absent = 1;
      d->has_mode = 0;
      d->dirfd = fd;
      d->final_name = name;
      free(cur);
      return 0;
    }
    if (S_ISLNK(st.st_mode)) {
      char *text = read_linkat_text(fd, name, err, err_cap);
      char *next;
      if (!text) {
        close(fd);
        free(name);
        free(cur);
        return -1;
      }
      d->hops[d->nhops].dirfd = fd;
      d->hops[d->nhops].name = name;
      d->hops[d->nhops].text = text;
      d->hops[d->nhops].dev = st.st_dev;
      d->hops[d->nhops].ino = st.st_ino;
      d->nhops++;
      next = link_target_path(cur, text, err, err_cap);
      free(cur);
      if (!next)
        return -1;
      cur = next;
      continue;
    }
    /* An existing non symlink. Replacing it is right only when the chain
     * was already known to end on a real file. After realpath reported
     * ENOENT this node is an intermediate that changed type, and the
     * rename must not adopt it. */
    d->absent = 0;
    d->dev = st.st_dev;
    d->ino = st.st_ino;
    d->mode = st.st_mode & 07777;
    d->has_mode = expect_absent ? 0 : 1;
    if (expect_absent)
      d->mismatch = 1;
    d->dirfd = fd;
    d->final_name = name;
    free(cur);
    return 0;
  }
  free(cur);
  set_err(err, err_cap, "the RESTART symlink chain has too many hops, likely a loop");
  return -1;
}

/* Probe the final component of `path` without lstat(path), whose path walk
 * follows a directory symlink in an intermediate component. The parent is
 * opened with open_directory_nofollow first, so a directory symlink that the
 * recorded RESTART chain does not hold is refused before any stat can pass
 * through it. */
static int probe_final(const char *path, struct stat *st, int *exists, char *err,
                       size_t err_cap) {
  char *dir = NULL, *name = NULL;
  int fd;
  if (split_parent(path, &dir, &name, err, err_cap) != 0)
    return -1;
  fd = open_directory_nofollow(dir);
  free(dir);
  if (fd < 0) {
    set_errf(err, err_cap, "cannot open the RESTART directory", strerror(errno));
    free(name);
    return -1;
  }
  if (fstatat(fd, name, st, AT_SYMLINK_NOFOLLOW) == 0) {
    *exists = 1;
  } else if (errno == ENOENT) {
    *exists = 0;
  } else {
    set_errf(err, err_cap, "cannot stat the RESTART file", strerror(errno));
    close(fd);
    free(name);
    return -1;
  }
  close(fd);
  free(name);
  return 0;
}

/* A path that is not itself a symlink. Hold its parent and remember
 * whether the final component exists, so the rename can refuse if that
 * component changes type. */
static int adopt_plain(const char *path, struct publish_dest *d, char *err,
                       size_t err_cap) {
  struct stat lst;
  char *dir = NULL, *name = NULL;
  int exists, fd;
  if (probe_final(path, &lst, &exists, err, err_cap) != 0)
    return -1;
  if (split_parent(path, &dir, &name, err, err_cap) != 0)
    return -1;
  fd = open_directory_nofollow(dir);
  free(dir);
  if (fd < 0) {
    set_errf(err, err_cap, "cannot open the RESTART directory", strerror(errno));
    free(name);
    return -1;
  }
  if (fstatat(fd, name, &lst, AT_SYMLINK_NOFOLLOW) != 0) {
    if (errno != ENOENT) {
      set_errf(err, err_cap, "cannot stat the RESTART file", strerror(errno));
      close(fd);
      free(name);
      return -1;
    }
    d->absent = 1;
    d->has_mode = 0;
    if (exists)
      d->mismatch = 1;
  } else if (S_ISLNK(lst.st_mode)) {
    d->mismatch = 1;
    d->absent = 0;
    d->has_mode = 0;
  } else {
    d->absent = 0;
    d->dev = lst.st_dev;
    d->ino = lst.st_ino;
    d->mode = lst.st_mode & 07777;
    d->has_mode = 1;
    if (!exists) {
      d->mismatch = 1;
      d->has_mode = 0;
    }
  }
  d->dirfd = fd;
  d->final_name = name;
  return 0;
}

static void set_chain_changed(char *err, size_t err_cap) {
  set_err(err, err_cap, "the RESTART symlink chain changed before publication");
}

/* Re-read every recorded symlink and the final component through the
 * directory descriptors taken during resolution. A symlink that is no
 * longer that link, or a final component that is not the absent name or
 * the same file, is a chain this call must not publish onto. */
static int chain_unchanged(const struct publish_dest *d, char *err, size_t err_cap) {
  struct stat st;
  if (d->mismatch || d->dirfd < 0 || !d->final_name) {
    set_chain_changed(err, err_cap);
    return -1;
  }
  for (int i = 0; i < d->nhops; i++) {
    char *text;
    const struct publish_hop *h = &d->hops[i];
    if (fstatat(h->dirfd, h->name, &st, AT_SYMLINK_NOFOLLOW) != 0 || !S_ISLNK(st.st_mode) ||
        st.st_dev != h->dev || st.st_ino != h->ino) {
      set_chain_changed(err, err_cap);
      return -1;
    }
    text = read_linkat_text(h->dirfd, h->name, err, err_cap);
    if (!text)
      return -1;
    if (strcmp(text, h->text) != 0) {
      free(text);
      set_chain_changed(err, err_cap);
      return -1;
    }
    free(text);
  }
  if (fstatat(d->dirfd, d->final_name, &st, AT_SYMLINK_NOFOLLOW) != 0) {
    if (errno != ENOENT) {
      set_errf(err, err_cap, "cannot stat the RESTART file", strerror(errno));
      return -1;
    }
    if (!d->absent) {
      set_chain_changed(err, err_cap);
      return -1;
    }
  } else if (d->absent || S_ISLNK(st.st_mode) || st.st_dev != d->dev || st.st_ino != d->ino) {
    set_chain_changed(err, err_cap);
    return -1;
  }
  return 0;
}

static int rename_onto_dest(const char *tmpname, const struct publish_dest *d, char *err,
                            size_t err_cap) {
#if defined(RENAME_NOREPLACE)
  if (d->absent) {
    if (renameat2(d->dirfd, tmpname, d->dirfd, d->final_name, RENAME_NOREPLACE) == 0)
      return 0;
    /* The final name appeared after the check. Do not replace it. A kernel
     * without this flag returns EINVAL or ENOSYS and the plain renameat
     * below is the POSIX fallback. */
    if (errno == EEXIST || errno == ENOTEMPTY || errno == EBUSY) {
      set_chain_changed(err, err_cap);
      return -1;
    }
    if (errno != EINVAL && errno != ENOSYS && errno != EPERM) {
      set_errf(err, err_cap, "cannot replace the RESTART file", strerror(errno));
      return -1;
    }
  }
#endif
  if (renameat(d->dirfd, tmpname, d->dirfd, d->final_name) != 0) {
    set_errf(err, err_cap, "cannot replace the RESTART file", strerror(errno));
    return -1;
  }
  return 0;
}

/* Check the resolved chain, then rename the temporary name onto the final
 * component. Both calls use the directory descriptor held from resolution. */
static int commit_dest(const char *tmpname, const struct publish_dest *d, char *err,
                       size_t err_cap) {
  if (chain_unchanged(d, err, err_cap) != 0)
    return -1;
  return rename_onto_dest(tmpname, d, err, err_cap);
}

/* rename(2) replaces a symlink itself, unlike fopen(3), which follows it.
 * Resolve an existing link first so CPMD's RESTART -> RESTART.1 convention
 * continues to update RESTART.1 while keeping RESTART as a link. realpath
 * fails with ENOENT when the target is absent; fopen created that target,
 * so a dangling link is published there as a new file. The walk that
 * follows records directory descriptors. The rename checks them again and
 * refuses if a symlink or the final component is no longer what was
 * resolved, instead of replacing a file that appeared in that window. */
static int publish_destination(const char *path, struct publish_dest *d, char *err,
                               size_t err_cap) {
  struct stat lst;
  char *resolved;
  int exists;
  dest_init(d);
  if (probe_final(path, &lst, &exists, err, err_cap) != 0)
    return -1;
  if (!(exists && S_ISLNK(lst.st_mode))) {
    if (adopt_plain(path, d, err, err_cap) != 0) {
      dest_clear(d);
      return -1;
    }
    return 0;
  }
  resolved = realpath(path, NULL);
  if (!resolved) {
    if (errno != ENOENT) {
      set_errf(err, err_cap, "cannot resolve the RESTART symlink", strerror(errno));
      return -1;
    }
#ifdef CPMDC_RESTART_TEST_HOOK
    /* The tests swap a link here. Production leaves the pointer unset,
     * so this is a null check and nothing else. */
    if (cpmdc_restart_test_before_follow)
      cpmdc_restart_test_before_follow();
#endif
    if (walk_chain(path, 1, d, err, err_cap) != 0) {
      dest_clear(d);
      return -1;
    }
    return 0;
  }
  free(resolved);
  if (walk_chain(path, 0, d, err, err_cap) != 0) {
    dest_clear(d);
    return -1;
  }
  return 0;
}
#endif

/* Flush the stream and the file it writes to stable storage. */
static int sync_stream(FILE *fp) {
  if (fflush(fp) != 0)
    return -1;
#ifdef _WIN32
  return _commit(_fileno(fp));
#else
  for (;;) {
    if (fsync(fileno(fp)) == 0)
      return 0;
    if (errno != EINTR)
      return -1;
  }
#endif
}

/* Make the rename durable. The descriptor is the directory inode held
 * from resolution, so this flush is not a second path lookup. Filesystems
 * that cannot sync a directory say EINVAL, and there is nothing more to
 * do on them. Windows has no directory handle; MOVEFILE_WRITE_THROUGH
 * covers that path. */
#ifndef _WIN32
static int sync_fd(int fd, char *err, size_t err_cap) {
  while (fsync(fd) != 0) {
    if (errno == EINTR)
      continue;
    if (errno != EINVAL) {
      set_errf(err, err_cap, "cannot sync the RESTART directory", strerror(errno));
      return -1;
    }
    break;
  }
  return 0;
}
#endif

#ifdef _WIN32
static int replace_path(const char *tmp, const char *path, char *err, size_t err_cap) {
  if (!MoveFileExA(tmp, path, MOVEFILE_REPLACE_EXISTING | MOVEFILE_WRITE_THROUGH)) {
    char why[64];
    snprintf(why, sizeof why, "MoveFileEx error %lu", (unsigned long)GetLastError());
    set_errf(err, err_cap, "cannot replace the RESTART file", why);
    return -1;
  }
  return 0;
}
#endif

/* Publish `path` by rename. A reader opens either the previous file or the
 * finished one, never a prefix of the new bytes. On POSIX the temporary
 * file is created in the directory descriptor held from resolution, and
 * renameat publishes it there, so both stay on one filesystem and on the
 * directory inode that was resolved. */
static int publish_bytes(const char *path, const void *bytes, size_t n, char *err,
                         size_t err_cap) {
  static _Thread_local unsigned long seq;
  char tmpname[sizeof ".cpmdc-tmp-" + 16];
  FILE *fp;
  int fd = -1;
  tmpname[0] = '\0';
#ifndef _WIN32
  struct publish_dest dest;
  if (publish_destination(path, &dest, err, err_cap) != 0)
    return -1;
#else
  char *effective = NULL;
  char *tmp = NULL;
  effective = publish_strdup(path);
  if (!effective) {
    set_err(err, err_cap, "out of memory");
    return -1;
  }
#endif
  /* The thread-local's address tells threads of one process apart. */
  for (int attempt = 0; attempt < 100 && fd < 0; attempt++) {
    unsigned long token =
        seq++ ^ publish_pid() ^ (unsigned long)(uintptr_t)&seq;
#ifndef _WIN32
    snprintf(tmpname, sizeof tmpname, ".cpmdc-tmp-%016lx", token);
    /* 0666 under the umask is what fopen(path, "wb") created before. */
    fd = openat(dest.dirfd, tmpname, O_WRONLY | O_CREAT | O_EXCL | O_CLOEXEC, 0666);
    if (fd < 0 && errno != EEXIST)
      break;
#else
    free(tmp);
    tmp = temporary_path(effective, token);
    if (!tmp) {
      set_err(err, err_cap, "out of memory");
      free(effective);
      return -1;
    }
    fd = publish_open(tmp);
    if (fd < 0 && errno != EEXIST)
      break;
#endif
  }
#ifdef CPMDC_RESTART_TEST_HOOK
  if (fd >= 0 && cpmdc_restart_test_after_temp_open)
    cpmdc_restart_test_after_temp_open(fd);
#endif
  if (fd < 0) {
    set_errf(err, err_cap, "cannot create a temporary RESTART file", strerror(errno));
#ifndef _WIN32
    dest_clear(&dest);
#else
    free(tmp);
    free(effective);
#endif
    return -1;
  }
  fp = publish_fdopen(fd, "wb");
  if (!fp) {
    set_errf(err, err_cap, "cannot open a temporary RESTART file", strerror(errno));
    publish_close(fd);
    goto discard;
  }
  if (n && fwrite(bytes, 1, n, fp) != n) {
    set_errf(err, err_cap, "cannot write the RESTART file", strerror(errno));
    fclose(fp);
    goto discard;
  }
#ifndef _WIN32
  if (dest.has_mode && fchmod(fileno(fp), dest.mode) != 0) {
    set_errf(err, err_cap, "cannot set the RESTART file mode", strerror(errno));
    fclose(fp);
    goto discard;
  }
#endif
  if (sync_stream(fp) != 0) {
    set_errf(err, err_cap, "cannot sync the RESTART file", strerror(errno));
    fclose(fp);
    goto discard;
  }
  if (fclose(fp) != 0) {
    set_errf(err, err_cap, "cannot close the RESTART file", strerror(errno));
    goto discard;
  }
#ifndef _WIN32
  if (commit_dest(tmpname, &dest, err, err_cap) != 0)
    goto discard;
  if (sync_fd(dest.dirfd, err, err_cap) != 0) {
    dest_clear(&dest);
    return -1;
  }
  dest_clear(&dest);
  return 0;
#else
  if (replace_path(tmp, effective, err, err_cap) != 0)
    goto discard;
  free(tmp);
  free(effective);
  return 0;
#endif
discard:
#ifndef _WIN32
  if (tmpname[0] != '\0')
    unlinkat(dest.dirfd, tmpname, 0);
  dest_clear(&dest);
#else
  publish_unlink(tmp);
  free(tmp);
  free(effective);
#endif
  return -1;
}

int cpmdc_restart_write_path(const CPMDCRestart *file, const char *path, char *err,
                             size_t err_cap) {
  void *bytes = NULL;
  size_t n = 0;
  if (!path) {
    set_err(err, err_cap, "null path");
    return -1;
  }
  if (cpmdc_restart_write_mem(file, &bytes, &n, err, err_cap) != 0)
    return -1;
  if (publish_bytes(path, bytes, n, err, err_cap) != 0) {
    free(bytes);
    return -1;
  }
  free(bytes);
  return 0;
}
