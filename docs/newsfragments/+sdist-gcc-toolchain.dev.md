The sdist job now sets `CC`, `CXX` and `FC` to the GCC toolchain it installs, so Meson and its cmake module no longer pick a clang `c++` that cannot link.
