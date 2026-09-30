/* Anonymous memfd for one-time OpenCPMD cold parse (no disk INPUT). */
#define _GNU_SOURCE
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/stat.h>

/* Returns fd (>=0) and fills path as /proc/self/fd/N. Caller keeps fd open. */
int cpmdc_memfd_write(const char *bytes, int nbytes, char *path_out, int path_cap) {
  int fd;
  ssize_t w;
  if (!bytes || nbytes < 0 || !path_out || path_cap < 32)
    return -1;
  fd = memfd_create("cpmdc_cold_deck", MFD_CLOEXEC);
  if (fd < 0)
    return -1;
  w = write(fd, bytes, (size_t)nbytes);
  if (w != (ssize_t)nbytes) {
    close(fd);
    return -1;
  }
  if (lseek(fd, 0, SEEK_SET) != 0) {
    close(fd);
    return -1;
  }
  snprintf(path_out, (size_t)path_cap, "/proc/self/fd/%d", fd);
  return fd;
}

/* Host CWD saved before leaving it for the pseudopotential library.
 * Restored after the SCF. A force call writes no RESTART, LATEST, or GEOMETRY. */
static char g_host_cwd[1024];
static int g_host_cwd_saved = 0;

static int save_host_cwd(void) {
  if (g_host_cwd_saved)
    return 0;
  if (getcwd(g_host_cwd, sizeof(g_host_cwd)) == NULL)
    return -1;
  g_host_cwd_saved = 1;
  return 0;
}

/*
 * Point process CWD at the pseudopotential library for cold memfd decks.
 *
 * OpenCPMD recpnew/get_pplib uses argv[2] as the PP library whenever the host
 * process has argc>1, ignoring CPMD_PP_LIBRARY_PATH. Catch2/eOn filters always
 * set argc>1, so relative *PP basenames only resolve via the second-chance
 * CWD lookup (basename alone). chdir to the library directory first.
 *
 * Also exports CPMD_PP_LIBRARY_PATH with a trailing slash for hosts that do
 * honor the env (argc==1 CLI runs). The exported value is the absolute
 * library path: the next call reads it back through
 * cpmdc_pseudopotential_directory, and a relative value only resolves from
 * the directory the host stood in on the first call. A relative
 * CPMDC_PSEUDO_DIR is rewritten to the same absolute path for that reason.
 * Saves the prior CWD. cpmdc_enter_output_cwd leaves the library before
 * CPMD writes, and cpmdc_restore_host_cwd returns to the host directory.
 *
 * Returns 0 on success, -1 on failure (missing dir / chdir failed).
 */
int cpmdc_prepare_pp_cwd(const char *pseudo_dir) {
  char dir[1024];
  char abs_dir[1024];
  char libpath[1100];
  const char *caller_dir;
  size_t n;
  struct stat st;

  if (!pseudo_dir || !pseudo_dir[0])
    return -1;
  n = strnlen(pseudo_dir, sizeof(dir) - 1);
  if (n == 0 || n >= sizeof(dir) - 1)
    return -1;
  memcpy(dir, pseudo_dir, n);
  dir[n] = '\0';
  while (n > 1 && (dir[n - 1] == '/' || dir[n - 1] == ' '))
    dir[--n] = '\0';
  if (stat(dir, &st) != 0 || !S_ISDIR(st.st_mode))
    return -1;
  if (save_host_cwd() != 0)
    return -1;
  if (chdir(dir) != 0)
    return -1;
  if (getcwd(abs_dir, sizeof(abs_dir)) == NULL)
    return -1;
  n = strlen(abs_dir);
  /* Trailing slash required when CPMD_PP_LIBRARY_PATH is set (OpenCPMD get_pplib). */
  if (n + 2 < sizeof(libpath)) {
    memcpy(libpath, abs_dir, n);
    libpath[n] = '/';
    libpath[n + 1] = '\0';
    (void)setenv("CPMD_PP_LIBRARY_PATH", libpath, 1);
    (void)setenv("PP_LIBRARY_PATH", libpath, 1);
  }
  caller_dir = getenv("CPMDC_PSEUDO_DIR");
  if (caller_dir && caller_dir[0] && caller_dir[0] != '/')
    (void)setenv("CPMDC_PSEUDO_DIR", abs_dir, 1);
  return 0;
}

/* Move to the directory CPMD should write into. Empty means the host CWD.
 * A relative output_dir is resolved from the host CWD, not from the
 * pseudopotential directory. */
int cpmdc_enter_output_cwd(const char *output_dir) {
  char dest[1024];
  size_t n;
  struct stat st;

  if (save_host_cwd() != 0)
    return -1;
  if (!output_dir || !output_dir[0])
    return chdir(g_host_cwd) == 0 ? 0 : -1;
  n = strnlen(output_dir, sizeof(dest) - 1);
  if (n == 0 || n >= sizeof(dest) - 1)
    return -1;
  memcpy(dest, output_dir, n);
  dest[n] = '\0';
  while (n > 1 && (dest[n - 1] == '/' || dest[n - 1] == ' '))
    dest[--n] = '\0';
  if (dest[0] != '/') {
    if (chdir(g_host_cwd) != 0)
      return -1;
  }
  if (stat(dest, &st) != 0 || !S_ISDIR(st.st_mode))
    return -1;
  if (chdir(dest) != 0)
    return -1;
  return 0;
}

/* Restore CWD saved by cpmdc_prepare_pp_cwd. Idempotent; 0 on success. */
int cpmdc_restore_host_cwd(void) {
  if (!g_host_cwd_saved)
    return 0;
  if (chdir(g_host_cwd) != 0)
    return -1;
  g_host_cwd_saved = 0;
  return 0;
}
