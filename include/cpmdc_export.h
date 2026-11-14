#pragma once

/**
 * @file cpmdc_export.h
 * @brief Symbol visibility of the public cpmdc C ABI.
 *
 * `CPMDC_API` marks every function of `cpmdc.h` and `cpmdc_features.h`.
 * libcpmdc is compiled with hidden visibility, so a function without the
 * mark is not exported from the shared library.
 *
 * With GCC and Clang the mark is default visibility. On Windows it depends on
 * two macros:
 * - `CPMDC_BUILDING_LIBRARY`: defined while compiling libcpmdc itself;
 *   `dllexport`.
 * - `CPMDC_STATIC`: defined by code that links the cpmdc sources or a static
 *   archive instead of the DLL; empty.
 * - Neither: a consumer of the DLL; `dllimport`.
 */

#if defined(_WIN32) || defined(__CYGWIN__)
#if defined(CPMDC_BUILDING_LIBRARY)
#define CPMDC_API __declspec(dllexport)
#elif defined(CPMDC_STATIC)
#define CPMDC_API
#else
#define CPMDC_API __declspec(dllimport)
#endif
#elif defined(__GNUC__) && __GNUC__ >= 4
#define CPMDC_API __attribute__((visibility("default")))
#else
#define CPMDC_API
#endif
