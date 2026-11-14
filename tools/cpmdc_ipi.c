/**
 * cpmdc-ipi-c: an i-PI socket client in C on cpmdc_session_calculate_au().
 *
 * It answers the servers that python/cpmdc/ipi.py (cpmdc-ipi) answers, the
 * same way: i-PI, ASE's SocketIOCalculator and SocketServer, and any server
 * with the i-PI wire format ("Distributed execution", "Communication
 * protocol", https://docs.ipi-code.org/distributed.html; i-PI's
 * ipi/interfaces/sockets.py; ASE's ase/calculators/socketio.py):
 *
 * - a header is 12 ASCII bytes padded with spaces;
 * - integers are int32 and reals float64, in host order;
 * - STATUS is answered with NEEDINIT, READY or HAVEDATA;
 * - INIT carries a bead index, a byte count, and that many bytes;
 * - POSDATA carries the cell h (lattice vectors as columns) and its inverse,
 *   nine reals each in row-major order, the atom count, and the positions;
 * - GETFORCE is answered with FORCEREADY, the energy, the atom count, the
 *   forces, the virial laid out as the cell, and an empty extra string;
 * - EXIT ends the run.
 *
 * Everything on the wire is atomic units, which the library takes and
 * returns, so nothing is converted. The virial is -dE/d(strain) in Hartree.
 * The protocol carries no species: they are fixed at start, and a step with
 * another atom count stops the client. One session serves every step.
 *
 * Exit status: 0 after EXIT or when the server closes the connection
 * between messages, 1 on an error, 2 on a usage error.
 */
#define _POSIX_C_SOURCE 200809L
#include "cpmdc.h"

#include <errno.h>
#include <math.h>
#include <netdb.h>
#include <netinet/in.h>
#include <netinet/tcp.h>
#include <signal.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/socket.h>
#include <sys/un.h>
#include <time.h>
#include <unistd.h>

#define PROG "cpmdc-ipi-c"
#define HEADER_BYTES 12
#define MAX_Z 118
/* Same value as CPMDC_BOHR_TO_ANGSTROM and cpmdc.BOHR_TO_ANGSTROM. */
#define BOHR_TO_ANGSTROM 0.529177210903

static const char *const k_symbols[MAX_Z] = {
    "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na", "Mg",
    "Al", "Si", "P",  "S",  "Cl", "Ar", "K",  "Ca", "Sc", "Ti", "V",  "Cr",
    "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
    "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "In", "Sn", "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf",
    "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po",
    "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm",
    "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs",
    "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"};

static char g_error[1024];
static int g_verbose = 0;

static int fail(const char *fmt, ...) {
  va_list ap;
  va_start(ap, fmt);
  vsnprintf(g_error, sizeof(g_error), fmt, ap);
  va_end(ap);
  return -1;
}

static void say(const char *fmt, ...) {
  if (!g_verbose)
    return;
  va_list ap;
  va_start(ap, fmt);
  fprintf(stderr, PROG ": ");
  vfprintf(stderr, fmt, ap);
  fputc('\n', stderr);
  va_end(ap);
}

/* --- species ------------------------------------------------------------ */

static int symbol_number(const char *s, size_t len) {
  for (int z = 0; z < MAX_Z; ++z)
    if (strlen(k_symbols[z]) == len && strncmp(k_symbols[z], s, len) == 0)
      return z + 1;
  return 0;
}

static int push_number(int **out, int *count, int *cap, int z) {
  if (*count == *cap) {
    int next = *cap ? 2 * *cap : 16;
    int *grown = (int *)realloc(*out, (size_t)next * sizeof(int));
    if (!grown)
      return fail("out of memory");
    *out = grown;
    *cap = next;
  }
  (*out)[(*count)++] = z;
  return 0;
}

/* One token: atomic numbers ("8"), or symbols with counts ("OH2"). */
static int species_token(const char *tok, size_t len, int **out, int *count,
                         int *cap) {
  size_t i = 0;
  int digits = 1;
  for (i = 0; i < len; ++i)
    digits &= tok[i] >= '0' && tok[i] <= '9';
  if (digits) {
    long z = 0;
    for (i = 0; i < len && z <= MAX_Z; ++i)
      z = 10 * z + (tok[i] - '0');
    if (z < 1 || z > MAX_Z)
      return fail("atomic number %.*s is outside 1..%d", (int)len, tok, MAX_Z);
    return push_number(out, count, cap, (int)z);
  }
  for (i = 0; i < len;) {
    if (tok[i] < 'A' || tok[i] > 'Z')
      return fail("cannot read species '%.*s'", (int)len, tok);
    ++i;
    if (i < len && tok[i] >= 'a' && tok[i] <= 'z')
      ++i;
    while (i < len && tok[i] >= '0' && tok[i] <= '9')
      ++i;
  }
  for (i = 0; i < len;) {
    size_t start = i++;
    if (i < len && tok[i] >= 'a' && tok[i] <= 'z')
      ++i;
    int z = symbol_number(tok + start, i - start);
    if (!z)
      return fail("unknown element '%.*s'", (int)(i - start), tok + start);
    long n = 1;
    if (i < len && tok[i] >= '0' && tok[i] <= '9') {
      n = 0;
      while (i < len && tok[i] >= '0' && tok[i] <= '9') {
        if (n > 1000000)
          return fail("species count in '%.*s' is too large", (int)len, tok);
        n = 10 * n + (tok[i++] - '0');
      }
    }
    for (long k = 0; k < n; ++k)
      if (push_number(out, count, cap, z) != 0)
        return -1;
  }
  return 0;
}

/* Symbols separated by spaces or commas ("O H H", "O,H,H"), run together
 * with counts ("OH2"), or atomic numbers ("8 1 1"), as cpmdc-ipi reads them. */
static int species_numbers(const char *text, int **out, int *count) {
  int cap = 0;
  *out = NULL;
  *count = 0;
  const char *p = text;
  while (*p) {
    while (*p == ',' || *p == ' ' || *p == '\t' || *p == '\n' || *p == '\r')
      ++p;
    const char *start = p;
    while (*p && *p != ',' && *p != ' ' && *p != '\t' && *p != '\n' &&
           *p != '\r')
      ++p;
    if (p > start && species_token(start, (size_t)(p - start), out, count,
                                   &cap) != 0)
      return -1;
  }
  if (*count == 0)
    return fail("species must name at least one atom");
  return 0;
}

/* --- socket ------------------------------------------------------------- */

/* i-PI's default sockets_prefix, which ASE's actualunixsocketname also
 * prepends: the system temporary directory under the root, then "ipi_".
 * Built from its parts as cpmdc-ipi builds it with os.path.join. */
static void default_socket_prefix(char *out, size_t cap) {
  const char sep = '/';
  snprintf(out, cap, "%c%s%c%s", sep, "tmp", sep, "ipi_");
}

static int unix_socket_path(const char *name, const char *prefix, char *out,
                            size_t cap) {
  char fallback[64];
  if (strchr(name, '/')) {
    if (strlen(name) >= cap)
      return fail("UNIX socket path %s is longer than %zu bytes", name,
                  cap - 1);
    snprintf(out, cap, "%s", name);
    return 0;
  }
  if (!prefix) {
    default_socket_prefix(fallback, sizeof(fallback));
    prefix = fallback;
  }
  if (strlen(prefix) + strlen(name) >= cap)
    return fail("UNIX socket path %s%s is longer than %zu bytes", prefix, name,
                cap - 1);
  snprintf(out, cap, "%s%s", prefix, name);
  return 0;
}

static double monotonic_seconds(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double)ts.tv_sec + 1.0e-9 * (double)ts.tv_nsec;
}

static int connect_unix(const char *path) {
  struct sockaddr_un addr;
  memset(&addr, 0, sizeof(addr));
  addr.sun_family = AF_UNIX;
  memcpy(addr.sun_path, path, strlen(path) + 1);
  int fd = socket(AF_UNIX, SOCK_STREAM, 0);
  if (fd < 0)
    return -1;
  if (connect(fd, (struct sockaddr *)&addr, sizeof(addr)) != 0) {
    int saved = errno;
    close(fd);
    errno = saved;
    return -1;
  }
  return fd;
}

static int connect_tcp(const char *host, int port, int *gai) {
  char service[16];
  struct addrinfo hints, *found = NULL;
  memset(&hints, 0, sizeof(hints));
  hints.ai_family = AF_UNSPEC;
  hints.ai_socktype = SOCK_STREAM;
  snprintf(service, sizeof(service), "%d", port);
  *gai = getaddrinfo(host, service, &hints, &found);
  if (*gai != 0)
    return -1;
  int fd = -1;
  int saved = ECONNREFUSED;
  for (struct addrinfo *ai = found; ai; ai = ai->ai_next) {
    fd = socket(ai->ai_family, ai->ai_socktype, ai->ai_protocol);
    if (fd < 0) {
      saved = errno;
      continue;
    }
    if (connect(fd, ai->ai_addr, ai->ai_addrlen) == 0)
      break;
    saved = errno;
    close(fd);
    fd = -1;
  }
  freeaddrinfo(found);
  if (fd < 0) {
    errno = saved;
    return -1;
  }
  int one = 1;
  setsockopt(fd, IPPROTO_TCP, TCP_NODELAY, &one, sizeof(one));
  return fd;
}

/* A server that is not listening yet is retried for wait seconds. */
static int open_socket(const char *unix_name, const char *prefix,
                       const char *host, int port, double wait) {
  char where[sizeof(((struct sockaddr_un *)0)->sun_path) + 64];
  char path[sizeof(((struct sockaddr_un *)0)->sun_path)];
  if (unix_name) {
    if (unix_socket_path(unix_name, prefix, path, sizeof(path)) != 0)
      return -1;
    snprintf(where, sizeof(where), "%s", path);
  } else {
    snprintf(where, sizeof(where), "%s:%d", host, port);
  }
  double deadline = monotonic_seconds() + (wait > 0.0 ? wait : 0.0);
  for (;;) {
    int gai = 0;
    int fd = unix_name ? connect_unix(path) : connect_tcp(host, port, &gai);
    if (fd >= 0)
      return fd;
    if (gai != 0)
      return fail("cannot connect to %s: %s", where, gai_strerror(gai));
    int saved = errno;
    if ((saved != ENOENT && saved != ECONNREFUSED) ||
        monotonic_seconds() >= deadline)
      return fail("cannot connect to %s: %s", where, strerror(saved));
    struct timespec pause = {0, 50000000L};
    nanosleep(&pause, NULL);
  }
}

/* 1 for every byte, 0 for a close before the first byte when at_boundary,
 * -1 otherwise. */
static int read_exact(int fd, void *buf, size_t n, const char *what,
                      int at_boundary) {
  unsigned char *p = (unsigned char *)buf;
  size_t got = 0;
  while (got < n) {
    ssize_t r = recv(fd, p + got, n - got, 0);
    if (r < 0 && errno == EINTR)
      continue;
    if (r < 0)
      return fail("reading %s: %s", what, strerror(errno));
    if (r == 0) {
      if (at_boundary && got == 0)
        return 0;
      return fail("the server closed the connection in the middle of %s",
                  what);
    }
    got += (size_t)r;
  }
  return 1;
}

static int send_all(int fd, const void *buf, size_t n) {
  const unsigned char *p = (const unsigned char *)buf;
  while (n > 0) {
    ssize_t w = send(fd, p, n, 0);
    if (w < 0 && errno == EINTR)
      continue;
    if (w < 0)
      return fail("sending to the server: %s", strerror(errno));
    p += w;
    n -= (size_t)w;
  }
  return 0;
}

static int send_header(int fd, const char *name) {
  char head[HEADER_BYTES];
  memset(head, ' ', sizeof(head));
  memcpy(head, name, strlen(name));
  if (send_all(fd, head, sizeof(head)) != 0)
    return -1;
  say("sent %s", name);
  return 0;
}

static int read_int32(int fd, int32_t *v, const char *what) {
  return read_exact(fd, v, sizeof(*v), what, 0) == 1 ? 0 : -1;
}

/* batch_size from an INIT string of key:value pairs separated by commas;
 * 1 when absent or not an integer. */
static long batch_size(const char *text) {
  const char *p = text;
  while (*p) {
    const char *end = strchr(p, ',');
    size_t len = end ? (size_t)(end - p) : strlen(p);
    const char *colon = memchr(p, ':', len);
    if (colon) {
      const char *k0 = p, *k1 = colon;
      while (k0 < k1 && (*k0 == ' ' || *k0 == '\t' || *k0 == '\n'))
        ++k0;
      while (k1 > k0 && (k1[-1] == ' ' || k1[-1] == '\t' || k1[-1] == '\n'))
        --k1;
      if ((size_t)(k1 - k0) == strlen("batch_size") &&
          strncmp(k0, "batch_size", (size_t)(k1 - k0)) == 0) {
        char value[64];
        size_t vlen = len - (size_t)(colon + 1 - p);
        if (vlen >= sizeof(value))
          return 1;
        memcpy(value, colon + 1, vlen);
        value[vlen] = '\0';
        char *stop = NULL;
        errno = 0;
        long n = strtol(value, &stop, 10);
        while (stop && (*stop == ' ' || *stop == '\t' || *stop == '\n'))
          ++stop;
        if (errno != 0 || stop == value || (stop && *stop != '\0'))
          return 1;
        return n;
      }
    }
    if (!end)
      break;
    p = end + 1;
  }
  return 1;
}

/* --- driver ------------------------------------------------------------- */

struct driver {
  CPMDCSession *session;
  int *numbers;
  int natoms;
  int cluster;
  int warned_virial;
  long steps;
  double *forces;
  double energy;
  double virial[9];
};

static double det3(const double *m) {
  return m[0] * (m[4] * m[8] - m[5] * m[7]) -
         m[1] * (m[3] * m[8] - m[5] * m[6]) +
         m[2] * (m[3] * m[7] - m[4] * m[6]);
}

static int check_count(const struct driver *d, int natoms) {
  if (natoms != d->natoms)
    return fail("the server sent %d atoms, but this driver was started for "
                "%d; the i-PI protocol carries no species, so a changed atom "
                "count or composition needs a new " PROG " process with "
                "matching --species",
                natoms, d->natoms);
  return 0;
}

/* The wire cell h has the lattice vectors as columns; the library takes
 * them as rows. Returns 0 with *periodic = 0 for a cluster. */
static int cell_rows(const struct driver *d, const double *h, double *rows,
                     int *periodic) {
  int any = 0;
  for (int i = 0; i < 9; ++i)
    any |= h[i] != 0.0;
  *periodic = 0;
  if (d->cluster || !any)
    return 0;
  double lengths[3];
  for (int k = 0; k < 3; ++k) {
    double sum = 0.0;
    for (int i = 0; i < 3; ++i) {
      rows[3 * k + i] = h[3 * i + k];
      double c = h[3 * i + k] * BOHR_TO_ANGSTROM;
      sum += c * c;
    }
    lengths[k] = sqrt(sum);
  }
  double smallest = fmin(lengths[0], fmin(lengths[1], lengths[2]));
  double b3 = BOHR_TO_ANGSTROM * BOHR_TO_ANGSTROM * BOHR_TO_ANGSTROM;
  if (!(smallest > 0.0) ||
      fabs(det3(h)) * b3 <= 1.0e-8 * (lengths[0] * lengths[1] * lengths[2]))
    return fail("the server sent a cell without three independent vectors; "
                "send an all-zero cell or start " PROG " with --cluster for "
                "an isolated system");
  *periodic = 1;
  return 0;
}

static int compute(struct driver *d, const double *h, const double *positions) {
  double rows[9];
  int periodic = 0;
  if (cell_rows(d, h, rows, &periodic) != 0)
    return -1;
  if (cpmdc_session_calculate_au(d->session, d->natoms, d->numbers, positions,
                                 periodic ? rows : NULL, &d->energy, d->forces,
                                 d->virial) != 0) {
    const char *why = cpmdc_last_error();
    return fail("%s", why && *why ? why : "evaluation failed");
  }
  if (periodic && !d->warned_virial) {
    CPMDCStressTensor stress;
    memset(&stress, 0, sizeof(stress));
    if (cpmdc_session_last_stress(d->session, &stress) != 0 || !stress.valid) {
      d->warned_virial = 1;
      fprintf(stderr, PROG ": warning: the library computed no stress for this "
                      "periodic cell (isolated deck or CPMDC_STRESS=0); the "
                      "virial sent is zero\n");
    }
  }
  d->steps += 1;
  return 0;
}

static int send_forces(int fd, const struct driver *d) {
  size_t n3 = 3u * (size_t)d->natoms;
  size_t size = HEADER_BYTES + sizeof(double) + sizeof(int32_t) +
                n3 * sizeof(double) + 9 * sizeof(double) + sizeof(int32_t);
  unsigned char *buf = (unsigned char *)malloc(size);
  if (!buf)
    return fail("out of memory");
  unsigned char *p = buf;
  memset(p, ' ', HEADER_BYTES);
  memcpy(p, "FORCEREADY", strlen("FORCEREADY"));
  p += HEADER_BYTES;
  memcpy(p, &d->energy, sizeof(double));
  p += sizeof(double);
  int32_t natoms = d->natoms;
  memcpy(p, &natoms, sizeof(natoms));
  p += sizeof(natoms);
  memcpy(p, d->forces, n3 * sizeof(double));
  p += n3 * sizeof(double);
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) {
      memcpy(p, &d->virial[3 * j + i], sizeof(double));
      p += sizeof(double);
    }
  int32_t extra = 0;
  memcpy(p, &extra, sizeof(extra));
  int rc = send_all(fd, buf, size);
  free(buf);
  if (rc == 0)
    say("sent FORCEREADY");
  return rc;
}

/* Answers the server until EXIT or a close between messages. */
static int serve(struct driver *d, int fd) {
  int initialised = 0;
  int pending = 0;
  double *positions = (double *)malloc(3u * (size_t)d->natoms * sizeof(double));
  if (!positions)
    return fail("out of memory");
  int rc = -1;
  for (;;) {
    char raw[HEADER_BYTES + 1];
    int got = read_exact(fd, raw, HEADER_BYTES, "a message header", 1);
    if (got < 0)
      break;
    if (got == 0) {
      say("the server closed the connection");
      rc = 0;
      break;
    }
    raw[HEADER_BYTES] = '\0';
    size_t a = 0, b = HEADER_BYTES;
    while (a < b && (raw[a] == ' ' || raw[a] == '\t' || raw[a] == '\0'))
      ++a;
    while (b > a && (raw[b - 1] == ' ' || raw[b - 1] == '\t' ||
                     raw[b - 1] == '\0' || raw[b - 1] == '\n'))
      --b;
    char head[HEADER_BYTES + 1];
    memcpy(head, raw + a, b - a);
    head[b - a] = '\0';
    for (size_t i = 0; i < b - a; ++i)
      if ((unsigned char)head[i] > 127) {
        fail("the server sent a header that is not ASCII");
        goto done;
      }
    say("received %s", head);
    if (strcmp(head, "EXIT") == 0) {
      rc = 0;
      break;
    } else if (strcmp(head, "STATUS") == 0) {
      const char *reply =
          pending ? "HAVEDATA" : (initialised ? "READY" : "NEEDINIT");
      if (send_header(fd, reply) != 0)
        break;
    } else if (strcmp(head, "INIT") == 0) {
      int32_t bead = 0, nbytes = 0;
      if (read_int32(fd, &bead, "INIT") != 0 ||
          read_int32(fd, &nbytes, "INIT") != 0)
        break;
      if (nbytes < 0) {
        fail("INIT announced %d bytes", (int)nbytes);
        break;
      }
      char *text = (char *)malloc((size_t)nbytes + 1);
      if (!text) {
        fail("out of memory");
        break;
      }
      if (nbytes > 0 && read_exact(fd, text, (size_t)nbytes, "INIT", 0) != 1) {
        free(text);
        break;
      }
      text[nbytes] = '\0';
      long batch = batch_size(text);
      free(text);
      if (batch > 1) {
        fail("the server asked for batch_size %ld; " PROG " evaluates one "
             "structure per request, so set batch_size to 1",
             batch);
        break;
      }
      say("bead %d, %d INIT bytes", (int)bead, (int)nbytes);
      initialised = 1;
    } else if (strcmp(head, "POSDATA") == 0) {
      double h[9], ih[9];
      int32_t natoms = 0;
      if (pending) {
        fail("POSDATA arrived before GETFORCE collected the last step");
        break;
      }
      if (read_exact(fd, h, sizeof(h), "POSDATA", 0) != 1 ||
          read_exact(fd, ih, sizeof(ih), "POSDATA", 0) != 1 ||
          read_int32(fd, &natoms, "POSDATA") != 0)
        break;
      if (natoms <= 0) {
        fail("POSDATA announced %d atoms", (int)natoms);
        break;
      }
      if (check_count(d, natoms) != 0)
        break;
      if (read_exact(fd, positions, 3u * (size_t)natoms * sizeof(double),
                     "POSDATA", 0) != 1)
        break;
      if (compute(d, h, positions) != 0)
        break;
      pending = 1;
      initialised = 1;
    } else if (strcmp(head, "GETFORCE") == 0) {
      if (!pending) {
        fail("GETFORCE arrived before POSDATA");
        break;
      }
      if (send_forces(fd, d) != 0)
        break;
      pending = 0;
    } else {
      fail("unknown message '%s' from the server", head);
      break;
    }
  }
done:
  free(positions);
  return rc;
}

/* --- command line ------------------------------------------------------- */

static void usage(FILE *to) {
  fprintf(to,
          "usage: " PROG " (--unix NAME | --port PORT) [--host HOST] "
          "--params FILE --species TEXT\n"
          "                   [--cluster] [--socket-prefix P] [--wait SECONDS] "
          "[--verbose]\n"
          "\n"
          "Serve cpmdc energies, forces, and virials to an i-PI protocol "
          "server (i-PI, ASE\n"
          "SocketIOCalculator or SocketServer) over a UNIX or TCP socket.\n"
          "\n"
          "  --unix NAME         UNIX socket address, as in i-PI and ASE\n"
          "  --port PORT         TCP port of the server\n"
          "  --host HOST         TCP host (default: localhost)\n"
          "  --params FILE       serialized CPMDParams file\n"
          "  --species TEXT      element symbols in server order: 'O H H', "
          "'O,H,H', or 'OH2'\n"
          "  --cluster           ignore the server's cell; the deck sets the "
          "box\n"
          "  --socket-prefix P   prefix for a UNIX socket name (default: "
          "i-PI's and ASE's)\n"
          "  --wait SECONDS      seconds to retry until the server listens\n"
          "  --verbose           log every protocol message\n");
}

static unsigned char *read_file(const char *path, size_t *size) {
  FILE *fp = fopen(path, "rb");
  unsigned char *buf = NULL;
  long n = 0;
  if (fp && fseek(fp, 0, SEEK_END) == 0 && (n = ftell(fp)) > 0 &&
      fseek(fp, 0, SEEK_SET) == 0 && (buf = malloc((size_t)n)) != NULL &&
      fread(buf, 1, (size_t)n, fp) != (size_t)n) {
    free(buf);
    buf = NULL;
  }
  if (fp)
    fclose(fp);
  *size = buf ? (size_t)n : 0;
  return buf;
}

struct options {
  const char *unix_name;
  const char *port_text;
  const char *host;
  const char *params;
  const char *species;
  const char *prefix;
  const char *wait_text;
  int cluster;
};

static int usage_error(const char *fmt, const char *arg) {
  usage(stderr);
  fprintf(stderr, PROG ": error: ");
  fprintf(stderr, fmt, arg);
  fputc('\n', stderr);
  return 2;
}

/* Returns -1 when the options parse, else the exit status. */
static int parse_options(int argc, char **argv, struct options *o) {
  static const struct {
    const char *name;
    size_t offset;
  } valued[] = {
      {"--unix", offsetof(struct options, unix_name)},
      {"--port", offsetof(struct options, port_text)},
      {"--host", offsetof(struct options, host)},
      {"--params", offsetof(struct options, params)},
      {"--species", offsetof(struct options, species)},
      {"--socket-prefix", offsetof(struct options, prefix)},
      {"--wait", offsetof(struct options, wait_text)},
  };
  for (int i = 1; i < argc; ++i) {
    const char *arg = argv[i];
    if (strcmp(arg, "-h") == 0 || strcmp(arg, "--help") == 0) {
      usage(stdout);
      return 0;
    }
    if (strcmp(arg, "--cluster") == 0) {
      o->cluster = 1;
      continue;
    }
    if (strcmp(arg, "--verbose") == 0) {
      g_verbose = 1;
      continue;
    }
    int matched = 0;
    for (size_t k = 0; k < sizeof(valued) / sizeof(valued[0]); ++k) {
      size_t len = strlen(valued[k].name);
      const char *value = NULL;
      if (strcmp(arg, valued[k].name) == 0) {
        if (i + 1 >= argc)
          return usage_error("%s needs a value", valued[k].name);
        value = argv[++i];
      } else if (strncmp(arg, valued[k].name, len) == 0 && arg[len] == '=') {
        value = arg + len + 1;
      } else {
        continue;
      }
      *(const char **)((char *)o + valued[k].offset) = value;
      matched = 1;
      break;
    }
    if (!matched)
      return usage_error("unrecognized argument %s", arg);
  }
  if ((o->unix_name == NULL) == (o->port_text == NULL))
    return usage_error("%s", "give exactly one of --unix and --port");
  if (!o->params)
    return usage_error("%s", "--params is required");
  if (!o->species)
    return usage_error("%s", "--species is required");
  return -1;
}

int main(int argc, char **argv) {
  struct options o;
  memset(&o, 0, sizeof(o));
  o.host = "localhost";
  int status = parse_options(argc, argv, &o);
  if (status >= 0)
    return status;

  long port = 0;
  double wait = 0.0;
  if (o.port_text) {
    char *stop = NULL;
    errno = 0;
    port = strtol(o.port_text, &stop, 10);
    if (errno != 0 || stop == o.port_text || *stop != '\0')
      return usage_error("--port: invalid int value: '%s'", o.port_text);
  }
  if (o.wait_text) {
    char *stop = NULL;
    wait = strtod(o.wait_text, &stop);
    if (stop == o.wait_text || *stop != '\0')
      return usage_error("--wait: invalid float value: '%s'", o.wait_text);
  }

  signal(SIGPIPE, SIG_IGN);
  struct driver d;
  memset(&d, 0, sizeof(d));
  d.cluster = o.cluster;
  int fd = -1;
  int rc = -1;
  unsigned char *params = NULL;
  size_t params_size = 0;
  if (o.port_text && (port <= 0 || port >= 65536)) {
    fail("port %ld is outside 1..65535", port);
    goto out;
  }
  if (species_numbers(o.species, &d.numbers, &d.natoms) != 0)
    goto out;
  if (d.natoms > INT32_MAX / 3) {
    fail("too many atoms in --species");
    goto out;
  }
  d.forces = (double *)malloc(3u * (size_t)d.natoms * sizeof(double));
  if (!d.forces) {
    fail("out of memory");
    goto out;
  }
  params = read_file(o.params, &params_size);
  if (!params) {
    fail("cannot read %s", o.params);
    goto out;
  }
  d.session = cpmdc_session_create(params, params_size);
  if (!d.session) {
    const char *why = cpmdc_last_error();
    fail("%s", why && *why ? why : "cpmdc_session_create failed");
    goto out;
  }
  fd = open_socket(o.unix_name, o.prefix, o.host, (int)port, wait);
  if (fd < 0)
    goto out;
  rc = serve(&d, fd);
  if (rc == 0)
    say("%ld steps", d.steps);

out:
  if (fd >= 0)
    close(fd);
  cpmdc_session_destroy(d.session);
  free(params);
  free(d.forces);
  free(d.numbers);
  if (rc != 0) {
    fprintf(stderr, PROG ": error: %s\n", g_error);
    return 1;
  }
  return 0;
}
