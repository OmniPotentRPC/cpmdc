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

/*
 * Point process CWD at the pseudopotential library for cold memfd decks.
 *
 * OpenCPMD recpnew/get_pplib uses argv[2] as the PP library whenever the host
 * process has argc>1, ignoring CPMD_PP_LIBRARY_PATH. Catch2/eOn filters always
 * set argc>1, so relative *PP basenames only resolve via the second-chance
 * CWD lookup (basename alone). chdir to the library directory first.
 *
 * Also exports CPMD_PP_LIBRARY_PATH with a trailing slash for hosts that do
 * honor the env (argc==1 CLI runs).
 *
 * Returns 0 on success, -1 on failure (missing dir / chdir failed).
 */
int cpmdc_prepare_pp_cwd(const char *pseudo_dir) {
  char dir[1024];
  char libpath[1100];
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
  if (chdir(dir) != 0)
    return -1;
  /* Trailing slash required when CPMD_PP_LIBRARY_PATH is set (OpenCPMD get_pplib). */
  if (n + 2 < sizeof(libpath)) {
    memcpy(libpath, dir, n);
    libpath[n] = '/';
    libpath[n + 1] = '\0';
    (void)setenv("CPMD_PP_LIBRARY_PATH", libpath, 1);
    (void)setenv("PP_LIBRARY_PATH", libpath, 1);
  }
  return 0;
}
