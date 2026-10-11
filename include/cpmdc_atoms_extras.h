#pragma once

/* Host-owned &ATOMS directives. embed_method_deck_plus_atoms drops these
 * and rebuilds coordinates from the step. The array ends with NULL.
 * P04 (rejecting inert decks) reuses this list. */
extern const char *const cpmdc_atoms_dropped_directives[];

/* Copy non-coordinate &ATOMS directives from deck into out, in order.
 * deck_len is the byte count and deck is not necessarily NUL-terminated.
 * z is the caller atomic-number order (length n_atoms). Species grouping
 * matches embed_method_deck_plus_atoms: first-seen Z, then atoms of that Z
 * in caller order. Bare atom indices on CONSTRAINTS, DUMMY ATOMS, and
 * CHANGE BONDS lines are rewritten into that order. A line whose only
 * token is one integer is a count and is left alone. Floats and tokens
 * that start with '+' or '-' are left alone. ISOTOPE, CONFINEMENT
 * POTENTIAL, ATOMIC CHARGES, and MOVIE TYPE are per species and are
 * copied unchanged. VELOCITIES, METADYNAMICS, coordinate lines, and
 * GENERATE COORDINATES are omitted.
 * Returns 0. Returns -1 and writes a message that names the directive
 * when an atom index is outside 1..n_atoms. out is NUL-terminated when
 * out_cap > *out_len. */
int cpmdc_atoms_extras_keep(const char *deck, int deck_len, const int *z,
                            int n_atoms, char *out, int out_cap, int *out_len,
                            char *err, int err_cap);
