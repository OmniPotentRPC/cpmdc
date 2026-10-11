/* Print the faulting instruction when a rank dies with SIGILL.
 * The sanitize job preloads this. It is not part of the MPI launcher.
 */
#define _GNU_SOURCE

#include <dlfcn.h>
#include <execinfo.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <ucontext.h>
#include <unistd.h>

static void on_sigill(int sig, siginfo_t *info, void *ctx) {
  ucontext_t *uc = ctx;
  unsigned long ip = (unsigned long)uc->uc_mcontext.gregs[REG_RIP];
  void *site = info->si_addr != NULL ? info->si_addr : (void *)ip;

  fprintf(stderr, "trace_sigill: SIGILL rip=0x%lx si_addr=%p code=%d\n", ip, info->si_addr, info->si_code);
  Dl_info dli;
  if (dladdr(site, &dli) != 0 && dli.dli_fname != NULL) {
    fprintf(stderr, "trace_sigill: file=%s base=%p\n", dli.dli_fname, dli.dli_fbase);
    if (dli.dli_sname != NULL) {
      fprintf(stderr, "trace_sigill: symbol=%s + 0x%lx\n", dli.dli_sname,
              (unsigned long)site - (unsigned long)dli.dli_saddr);
    }
    fprintf(stderr, "trace_sigill: file_offset=0x%lx\n",
            (unsigned long)site - (unsigned long)dli.dli_fbase);
  }
  unsigned char *bytes = site;
  fprintf(stderr, "trace_sigill: bytes");
  for (int i = 0; i < 16; i++) {
    fprintf(stderr, " %02x", bytes[i]);
  }
  fprintf(stderr, "\n");
  void *frames[64];
  int n = backtrace(frames, 64);
  backtrace_symbols_fd(frames, n, STDERR_FILENO);
  _exit(128 + sig);
}

__attribute__((constructor(65535))) static void install(void) {
  struct sigaction sa;

  sa.sa_sigaction = on_sigill;
  sigemptyset(&sa.sa_mask);
  sa.sa_flags = SA_SIGINFO;
  sigaction(SIGILL, &sa, NULL);
  fprintf(stderr, "trace_sigill: handler installed pid=%d\n", (int)getpid());
}
