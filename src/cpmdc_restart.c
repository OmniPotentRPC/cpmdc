#define _GNU_SOURCE
#include "cpmdc_restart.h"

#include <errno.h>
#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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

static int read_i32(const struct rec *r, size_t off, int32_t *out) {
  if (!r || off + 4 > r->len)
    return -1;
  memcpy(out, r->data + off, 4);
  return 0;
}

static int read_f64(const struct rec *r, size_t off, double *out) {
  if (!r || off + 8 > r->len)
    return -1;
  memcpy(out, r->data + off, 8);
  return 0;
}

static int push_rec(struct rec **recs, size_t *n, size_t *cap, unsigned char *data,
                    size_t len) {
  if (*n == *cap) {
    size_t ncap = *cap ? *cap * 2 : 64;
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

static int index_sections(cpmdc_restart *file, char *err, size_t err_cap) {
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

static int decode_typed(cpmdc_restart *file, char *err, size_t err_cap) {
  if (file->sec_nrec[2] == 3) {
    const struct rec *sym = &file->recs[file->sec_start[2] + 1];
    const struct rec *cell = &file->recs[file->sec_start[2] + 2];
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
    const struct rec *nsp_r = &file->recs[file->sec_start[3] + 1];
    const struct rec *na_r = &file->recs[file->sec_start[3] + 2];
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
    const struct rec *r = &file->recs[file->sec_start[7] + 1];
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
    const struct rec *r = &file->recs[file->sec_start[8] + 1];
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

static cpmdc_restart *from_bytes(const void *bytes, size_t nbytes, char *err,
                                 size_t err_cap) {
  cpmdc_restart *file;
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

cpmdc_restart *cpmdc_restart_read_mem(const void *bytes, size_t nbytes, char *err,
                                      size_t err_cap) {
  return from_bytes(bytes, nbytes, err, err_cap);
}

cpmdc_restart *cpmdc_restart_read_path(const char *path, char *err, size_t err_cap) {
  FILE *fp;
  unsigned char *buf = NULL;
  size_t n = 0, cap = 0;
  cpmdc_restart *file;
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

void cpmdc_restart_free(cpmdc_restart *file) {
  if (!file)
    return;
  free(file->na);
  free_recs(file);
  free(file);
}

int cpmdc_restart_is_stream(const cpmdc_restart *file) {
  return file && file->stream;
}

const char *cpmdc_restart_header(const cpmdc_restart *file) {
  return file ? file->header : NULL;
}

int cpmdc_restart_section(const cpmdc_restart *file, int section, int *count) {
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

int cpmdc_restart_cell(const cpmdc_restart *file, int *ibrav, int *indpg,
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

int cpmdc_restart_species(const cpmdc_restart *file, int *nsp,
                          const int **na_per_species) {
  if (!file || file->nsp <= 0 || !file->na)
    return -1;
  if (nsp)
    *nsp = file->nsp;
  if (na_per_species)
    *na_per_species = file->na;
  return 0;
}

int cpmdc_restart_ncoords(const cpmdc_restart *file) {
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

static int copy_triples(const cpmdc_restart *file, int section, double *xyz, int n3) {
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

static int store_triples(cpmdc_restart *file, int section, const double *xyz, int n3) {
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

int cpmdc_restart_coordinates(const cpmdc_restart *file, double *xyz, int n3) {
  return copy_triples(file, 4, xyz, n3);
}

int cpmdc_restart_velocities(const cpmdc_restart *file, double *xyz, int n3) {
  return copy_triples(file, 5, xyz, n3);
}

int cpmdc_restart_initial_coordinates(const cpmdc_restart *file, double *xyz, int n3) {
  return copy_triples(file, 6, xyz, n3);
}

int cpmdc_restart_set_coordinates(cpmdc_restart *file, const double *xyz, int n3) {
  return store_triples(file, 4, xyz, n3);
}

int cpmdc_restart_set_velocities(cpmdc_restart *file, const double *xyz, int n3) {
  return store_triples(file, 5, xyz, n3);
}

int cpmdc_restart_set_cell(cpmdc_restart *file, const double celldm[6]) {
  struct rec *r;
  if (!file || !file->has_cell || !celldm || file->sec_nrec[2] != 3)
    return -1;
  r = &file->recs[file->sec_start[2] + 2];
  if (r->len < 48)
    return -1;
  memcpy(r->data, celldm, 48);
  memcpy(file->celldm, celldm, sizeof(file->celldm));
  return 0;
}

int cpmdc_restart_cutoff(const cpmdc_restart *file, double *ecut, double *cdual,
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

int cpmdc_restart_states(const cpmdc_restart *file, int *n, int *nkpts, int *ngw,
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
    last = (off + chunk == dlen);
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

int cpmdc_restart_write_mem(const cpmdc_restart *file, void **bytes, size_t *nbytes,
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

int cpmdc_restart_write_path(const cpmdc_restart *file, const char *path, char *err,
                             size_t err_cap) {
  void *bytes = NULL;
  size_t n = 0;
  FILE *fp;
  if (!path) {
    set_err(err, err_cap, "null path");
    return -1;
  }
  if (cpmdc_restart_write_mem(file, &bytes, &n, err, err_cap) != 0)
    return -1;
  fp = fopen(path, "wb");
  if (!fp) {
    set_err(err, err_cap, strerror(errno));
    free(bytes);
    return -1;
  }
  if (n && fwrite(bytes, 1, n, fp) != n) {
    set_err(err, err_cap, "write failed");
    fclose(fp);
    free(bytes);
    return -1;
  }
  if (fclose(fp) != 0) {
    set_err(err, err_cap, "close failed");
    free(bytes);
    return -1;
  }
  free(bytes);
  return 0;
}
