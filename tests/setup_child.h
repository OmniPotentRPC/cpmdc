/**
 * One OpenCPMD setup per child process.
 *
 * OpenCPMD sets itself up once per process, and cpmdc refuses a second
 * setup (a new basis, or a call after a CPMD stop). A test that needs
 * several bases runs each one through run_setup_child: the child runs the
 * calls of that basis, fills a plain result struct, writes it to a pipe,
 * and exits. The parent reads the struct back and makes the assertions.
 */
#ifndef CPMDC_TESTS_SETUP_CHILD_H
#define CPMDC_TESTS_SETUP_CHILD_H

#include <errno.h>
#include <signal.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

typedef void (*setup_child_fn)(void *ctx, void *out);

/* Runs fn(ctx, out) in a child process and copies the child's *out (size
 * bytes) back into *out. Returns 0 when the child wrote every byte and
 * exited with status 0. Otherwise returns -1 and says why in why[]. A CPMD
 * STOP inside the child exits it before the write, which shows up here as
 * a short result. */
static int run_setup_child(setup_child_fn fn, void *ctx, void *out,
                           size_t size, char *why, size_t why_size) {
  int fds[2];
  why[0] = '\0';
  if (pipe(fds) != 0) {
    snprintf(why, why_size, "pipe: %s", strerror(errno));
    return -1;
  }
  fflush(NULL);
  pid_t pid = fork();
  if (pid < 0) {
    snprintf(why, why_size, "fork: %s", strerror(errno));
    close(fds[0]);
    close(fds[1]);
    return -1;
  }
  if (pid == 0) {
    close(fds[0]);
    /* The test runner's fault handlers would resume the runner inside the
     * child; a fault here ends the child and the parent reports it. */
    static const int faults[] = {SIGSEGV, SIGBUS, SIGILL, SIGFPE, SIGABRT};
    for (size_t i = 0; i < sizeof(faults) / sizeof(faults[0]); ++i)
      signal(faults[i], SIG_DFL);
    memset(out, 0, size);
    fn(ctx, out);
    const unsigned char *p = (const unsigned char *)out;
    size_t left = size;
    int rc = 0;
    while (left > 0) {
      ssize_t n = write(fds[1], p, left);
      if (n < 0) {
        if (errno == EINTR)
          continue;
        rc = 1;
        break;
      }
      p += n;
      left -= (size_t)n;
    }
    close(fds[1]);
    fflush(NULL);
    exit(rc);
  }
  close(fds[1]);
  unsigned char *q = (unsigned char *)out;
  size_t got = 0;
  while (got < size) {
    ssize_t n = read(fds[0], q + got, size - got);
    if (n < 0) {
      if (errno == EINTR)
        continue;
      break;
    }
    if (n == 0)
      break;
    got += (size_t)n;
  }
  close(fds[0]);
  int status = 0;
  while (waitpid(pid, &status, 0) < 0) {
    if (errno != EINTR) {
      snprintf(why, why_size, "waitpid: %s", strerror(errno));
      return -1;
    }
  }
  if (WIFSIGNALED(status)) {
    snprintf(why, why_size, "setup child killed by signal %d",
             WTERMSIG(status));
    return -1;
  }
  if (!WIFEXITED(status) || WEXITSTATUS(status) != 0) {
    snprintf(why, why_size, "setup child exited with status %d",
             WIFEXITED(status) ? WEXITSTATUS(status) : -1);
    return -1;
  }
  if (got != size) {
    snprintf(why, why_size, "setup child wrote %zu of %zu result bytes", got,
             size);
    return -1;
  }
  return 0;
}

/* Runs one setup child and fails the current cmocka test when it does not
 * hand back its whole result. Include cmocka.h before this header. */
#define RUN_SETUP_CHILD(fn, ctx, out)                                          \
  do {                                                                         \
    char setup_why_[256];                                                      \
    if (run_setup_child((fn), (ctx), (out), sizeof(*(out)), setup_why_,        \
                        sizeof(setup_why_)) != 0)                              \
      fail_msg("%s", setup_why_);                                              \
  } while (0)

#endif
