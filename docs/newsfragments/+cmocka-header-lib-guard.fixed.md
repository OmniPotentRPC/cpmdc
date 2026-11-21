Meson now performs a configure-time cmocka compile+link probe that catches header/library mismatches early and aborts with a clear diagnostic instead of deferring failure to the test link step.
