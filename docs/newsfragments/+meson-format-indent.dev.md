Pin `meson.format` to two-space indentation so `meson format` matches the
in-tree Meson style.
`pixi run fmt-meson` formats the build files in place, and CI runs
`pixi run fmt-meson-check`, which fails when they would change.
