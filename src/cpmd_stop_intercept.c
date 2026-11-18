/* Arm cpmd_stopgm_hook for the duration of an embed call.
 * The patched archive exports that pointer. The catch records the stop
 * code and returns 1 while the call is armed, so stopgm returns to its
 * caller. Only an OpenCPMD build references the symbol. On ELF that
 * reference is weak, so a link without the archive leaves it absent.
 * ld64 rejects the same undefined weak reference, and the serial object
 * therefore does not mention the symbol at all.
 * cpmdc_embed_catch and cpmdc_note_stop stay exported. */
static int g_code = 0;
static int g_catch = 0;

typedef int (*cpmd_stopgm_hook_fn)(int code);
#if defined(CPMDC_HAS_CPMD)
extern cpmd_stopgm_hook_fn cpmd_stopgm_hook __attribute__((weak));

static void set_hook(cpmd_stopgm_hook_fn fn) {
  if (&cpmd_stopgm_hook != 0)
    cpmd_stopgm_hook = fn;
}
#else
static void set_hook(cpmd_stopgm_hook_fn fn) { (void)fn; }
#endif

int cpmdc_embed_catch(void) { return g_catch; }

void cpmdc_note_stop(int code) { g_code = code ? code : 1; }

static int cpmdc_stopgm_catch(int code) {
  cpmdc_note_stop(code);
  return g_catch ? 1 : 0;
}

void cpmdc_stop_arm(void) {
  g_code = 0;
  g_catch = 1;
  set_hook(cpmdc_stopgm_catch);
}

void cpmdc_stop_disarm(void) {
  g_catch = 0;
  set_hook(0);
}

int cpmdc_stop_code(void) { return g_code; }

void end_swap_(void) {}
void tistopgm_(int *file_unit) { (void)file_unit; }
void my_stopall_(int *code) { g_code = code ? *code : 1; }
