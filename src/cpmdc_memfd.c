/* Anonymous memfd for one-time OpenCPMD cold parse (no disk INPUT). */
#define _GNU_SOURCE
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <sys/mman.h>

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
