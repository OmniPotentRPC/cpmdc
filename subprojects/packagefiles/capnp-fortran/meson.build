# Consumer-facing Meson build for capnp-fortran (runtime + capnpc-fortran).
# Replaces the interop-only top-level meson.build when used as a wrap subproject.
project(
  'capnp-fortran',
  'fortran',
  version: '0.1.0',
  license: 'MIT',
  default_options: ['warning_level=1'],
)

# Core runtime + schema + RPC (generated Potentials.capnp may reference interfaces).
capnp_fsrc = files(
  'src/capnp_kinds.f90',
  'src/capnp_endian.f90',
  'src/capnp_pointer.f90',
  'src/capnp_arena.f90',
  'src/capnp_message.f90',
  'src/capnp_packed.f90',
  'src/capnp_union.f90',
  'src/capnp_serialize.f90',
  'src/capnp_stream.f90',
  'src/capnp_canonical.f90',
  'src/capnp.f90',
  'src/capnp_cabi.f90',
  'src/capnp_posix.f90',
  'src/capnp_rpc_transport.f90',
  'src/rpc_capnp.f90',
  'src/rpc_twoparty_capnp.f90',
  'src/capnp_rpc.f90',
  'src/capnp_schema.f90',
  'src/capnp_dynamic.f90',
  'src/stream_capnp.f90',
)

capnp_flib = static_library(
  'capnp_f',
  capnp_fsrc,
  install: false,
)

capnp_fortran_dep = declare_dependency(
  link_with: capnp_flib,
)

# Schema compiler plugin for `capnp compile -o...`.
capnpc_fortran = executable(
  'capnpc-fortran',
  files(
    'app/capnpc_emit.f90',
    'app/main.f90',
  ),
  link_with: capnp_flib,
  native: true,
  install: false,
)

meson.override_dependency('capnp-fortran', capnp_fortran_dep)
