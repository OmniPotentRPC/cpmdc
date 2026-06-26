potentials_capnp_c = custom_target(
  'potentials-capnp-c',
  input: 'Potentials.capnp',
  output: ['Potentials.capnp.c', 'Potentials.capnp.h'],
  command: [
    python,
    meson.project_source_root() / 'tools/capnpc_c_compile.py',
    capnp,
    capnpc_c,
    meson.current_source_dir(),
    meson.current_build_dir(),
    '@INPUT@',
  ],
)

potentials_capnp_c_src = potentials_capnp_c[0]
