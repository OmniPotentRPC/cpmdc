/* Run on the sanitize runner. Prints whether two instructions that only
 * the failing runners advertise actually execute. A SIGILL is caught here
 * and reported. This is not linked into the MPI launcher.
 */
#define _GNU_SOURCE

#include <cpuid.h>
#include <setjmp.h>
#include <signal.h>
#include <stdio.h>
#include <string.h>

static sigjmp_buf jump_here;

static void on_sigill(int sig, siginfo_t *info, void *ctx) {
  (void)sig;
  (void)info;
  (void)ctx;
  siglongjmp(jump_here, 1);
}

static int cpuinfo_has(const char *flag) {
  FILE *fp = fopen("/proc/cpuinfo", "r");
  char line[4096];
  int found = 0;

  if (fp == NULL) {
    return 0;
  }
  while (fgets(line, sizeof line, fp) != NULL) {
    char *p;
    if (strncmp(line, "flags", 5) != 0) {
      continue;
    }
    p = line;
    while ((p = strstr(p, flag)) != NULL) {
      char after = p[strlen(flag)];
      if ((p == line || p[-1] == ' ' || p[-1] == '\t') && (after == ' ' || after == '\n' || after == '\0')) {
        found = 1;
      }
      p += strlen(flag);
    }
    break;
  }
  fclose(fp);
  return found;
}

static int bit_set(unsigned int leaf, unsigned int sub, int reg, unsigned int bit) {
  unsigned int a = 0, b = 0, c = 0, d = 0, value;
  if (__get_cpuid_count(leaf, sub, &a, &b, &c, &d) == 0) {
    return 0;
  }
  value = reg == 0 ? a : reg == 1 ? b : reg == 2 ? c : d;
  return (value & bit) != 0;
}

static void report(const char *name, int cpuinfo, int cpuid, void (*insn)(void)) {
  struct sigaction sa, old;
  int trapped;

  printf("probe_isa: %s cpuinfo=%d cpuid=%d", name, cpuinfo, cpuid);
  if (!cpuinfo && !cpuid) {
    printf(" skipped\n");
    return;
  }
  memset(&sa, 0, sizeof sa);
  sa.sa_sigaction = on_sigill;
  sigemptyset(&sa.sa_mask);
  sa.sa_flags = SA_SIGINFO;
  sigaction(SIGILL, &sa, &old);
  if (sigsetjmp(jump_here, 1) == 0) {
    insn();
    trapped = 0;
  } else {
    trapped = 1;
  }
  sigaction(SIGILL, &old, NULL);
  printf(trapped ? " SIGILL\n" : " executes\n");
}

static void insn_avx_vnni(void) {
  /* vpdpbusd xmm0, xmm0, xmm0. AVX-VNNI, not AVX-512 VNNI. */
  __asm__ volatile("vpxor %%xmm0, %%xmm0, %%xmm0\n\t"
                   ".byte 0xc4, 0xe2, 0x79, 0x50, 0xc0"
                   ::: "xmm0");
}

static void insn_vp2intersect(void) {
  /* vp2intersectd zmm0, zmm0, k0 */
  __asm__ volatile("vpxor %%xmm0, %%xmm0, %%xmm0\n\t"
                   ".byte 0x62, 0xf2, 0x7f, 0x48, 0x68, 0xc0\n\t"
                   "vzeroupper"
                   ::: "xmm0", "zmm0");
}

int main(void) {
  int vnni_cpuinfo = cpuinfo_has("avx_vnni");
  int vp2_cpuinfo = cpuinfo_has("avx512_vp2intersect");
  int vnni_cpuid = bit_set(7, 1, 0, bit_AVXVNNI);
  int vp2_cpuid = bit_set(7, 0, 3, bit_AVX512VP2INTERSECT);

  printf("probe_isa: avx512f_cpuinfo=%d\n", cpuinfo_has("avx512f"));
  report("vpdpbusd", vnni_cpuinfo, vnni_cpuid, insn_avx_vnni);
  report("vp2intersectd", vp2_cpuinfo, vp2_cpuid, insn_vp2intersect);
  return 0;
}
