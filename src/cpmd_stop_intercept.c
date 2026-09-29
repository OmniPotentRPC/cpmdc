/* Record an OpenCPMD stopgm while an embed call is active.
 * stopgm returns instead of calling my_stopall, so this file does not
 * longjmp across the Fortran procedure. */
static int g_code = 0;
static int g_catch = 0;

void cpmdc_stop_arm(void) {
  g_code = 0;
  g_catch = 1;
}
void cpmdc_stop_disarm(void) { g_catch = 0; }
int cpmdc_stop_code(void) { return g_code; }

int cpmdc_embed_catch(void) { return g_catch; }
void cpmdc_note_stop(int code) { g_code = code ? code : 1; }

void end_swap_(void) {}
void tistopgm_(int *file_unit) { (void)file_unit; }
void my_stopall_(int *code) {
  g_code = code ? *code : 1;
}
