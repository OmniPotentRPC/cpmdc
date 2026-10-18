/* Ubuntu 24.04 cmocka.h uses jmp_buf and va_list and does not include their headers.
 * This file is passed with -include so those names exist.
 * _GNU_SOURCE has to be set first. A later define in the source is too late
 * once features.h has been included, and then MFD_CLOEXEC stays undeclared. */
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <setjmp.h>
#include <stdarg.h>
