/* Intercept OpenCPMD fatal path so embed can CALL library routines in-process. */
#include <setjmp.h>
#include <stdio.h>

static int g_code = 0;
static int g_armed = 0;
static jmp_buf g_jmp;

void cpmdc_stop_arm(void) { g_code = 0; g_armed = 1; }
void cpmdc_stop_disarm(void) { g_armed = 0; }
int cpmdc_stop_code(void) { return g_code; }
jmp_buf *cpmdc_stop_jmp(void) { return &g_jmp; }

void end_swap_(void) {}
void tistopgm_(int *file_unit) { (void)file_unit; }
void my_stopall_(int *code) {
  g_code = code ? *code : 1;
  if (g_armed)
    longjmp(g_jmp, 1);
}
