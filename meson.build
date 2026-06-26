project(
  'cpmdc',
  ['c', 'fortran'],
  version: '0.1.0',
  default_options: ['warning_level=2', 'c_std=c11'],
)

cc = meson.get_compiler('c')
fc = meson.get_compiler('fortran')
python = find_program('python3')
math_dep = cc.find_library('m', required: false)
cmocka_dep = dependency('cmocka', required: get_option('with_tests'))

inc = include_directories('include', 'src')
cpmd_root = get_option('cpmd_root')
has_cpmd = get_option('with_cpmd') and cpmd_root != ''
need_capnp_c = get_option('with_tests') or has_cpmd

common_c_args = []
if cc.get_id() != 'msvc'
  common_c_args += ['-Wall', '-Wextra', '-Wpedantic']
endif

cpmdc_params_sources = []
cpmdc_params_deps = []
cpmdc_params_includes = []
if need_capnp_c
  capnp = find_program('capnp')
  capnp_c_project = subproject('c-capnproto')
  capnp_c_dep = capnp_c_project.get_variable('capnp_c_dep')
  capnpc_c = capnp_c_project.get_variable('capnpc_c')
  subdir('schema')
  potentials_capnp_c_dep = declare_dependency(sources: potentials_capnp_c)
  cpmdc_params_sources = ['src/cpmdc_params.c', potentials_capnp_c_src]
  cpmdc_params_deps = [capnp_c_dep, potentials_capnp_c_dep]
  cpmdc_params_includes = [include_directories('schema')]
endif

cpmdc_stub = static_library(
  'cpmdc_stub',
  'src/cpmdc_stub.c',
  include_directories: inc,
  c_args: common_c_args,
  install: false,
)

cpmdc_stub_dep = declare_dependency(
  link_with: cpmdc_stub,
  include_directories: inc,
)

if has_cpmd
  cpmd_root_check = run_command(
    python,
    [
      '-c',
      'import pathlib, sys; p = pathlib.Path(sys.argv[1]); sys.exit(0 if p.is_dir() else 1)',
      cpmd_root,
    ],
    check: false,
  )
  if cpmd_root_check.returncode() != 0
    error('cpmd_root is not a directory: ' + cpmd_root)
  endif

  # v0.1 installs the ISO_C embed shell + Cap'n Proto ABI. Full OpenCPMD
  # archive linking lands when embed evaluation is wired.
  cpmd_c_args = common_c_args + ['-DCPMDC_HAS_CPMD=1']
  cpmd_f90_args = ['-DCPMDC_HAS_CPMD=1']
  cpmdc = shared_library(
    'cpmdc',
    ['src/cpmdc.c', 'src/cpmd_embed_c_api.f90'] + cpmdc_params_sources,
    include_directories: [inc] + cpmdc_params_includes,
    dependencies: cpmdc_params_deps + [math_dep],
    c_args: cpmd_c_args,
    fortran_args: cpmd_f90_args,
    install: true,
  )
  install_headers('include/cpmdc.h')
endif

# Frontend embed shell: same C + Fortran ABI without CPMDC_HAS_CPMD so
# cpmdc_available() is 0 but Cap'n Proto socket helpers still work.
cpmdc_embed_shell = static_library(
  'cpmdc_embed_shell',
  ['src/cpmdc.c', 'src/cpmd_embed_c_api.f90'] + cpmdc_params_sources,
  include_directories: [inc] + cpmdc_params_includes,
  dependencies: cpmdc_params_deps + [math_dep],
  c_args: common_c_args,
  install: false,
)

cpmdc_embed_shell_dep = declare_dependency(
  link_with: cpmdc_embed_shell,
  include_directories: inc,
  dependencies: cpmdc_params_deps,
)

if get_option('with_tests')
  cpmd_parser_params_bin = custom_target(
    'cpmd_parser_params_bin',
    input: 'tests/cpmd_params_parser.capnp.txt',
    output: 'cpmd_params_parser.bin',
    command: [
      python,
      meson.current_source_dir() / 'tests/encode_capnp.py',
      capnp,
      meson.current_source_dir() / 'schema/Potentials.capnp',
      'CPMDParams',
      '@INPUT@',
      '@OUTPUT@',
    ],
  )
  force_input_step_a_bin = custom_target(
    'force_input_step_a_bin',
    input: 'tests/force_input_step_a.capnp.txt',
    output: 'force_input_step_a.bin',
    command: [
      python,
      meson.current_source_dir() / 'tests/encode_capnp.py',
      capnp,
      meson.current_source_dir() / 'schema/Potentials.capnp',
      'ForceInput',
      '@INPUT@',
      '@OUTPUT@',
    ],
  )
  force_input_step_ev_bin = custom_target(
    'force_input_step_ev_bin',
    input: 'tests/force_input_step_ev.capnp.txt',
    output: 'force_input_step_ev.bin',
    command: [
      python,
      meson.current_source_dir() / 'tests/encode_capnp.py',
      capnp,
      meson.current_source_dir() / 'schema/Potentials.capnp',
      'ForceInput',
      '@INPUT@',
      '@OUTPUT@',
    ],
  )

  test(
    'stub-abi',
    executable(
      'test_stub_abi',
      'tests/test_stub_abi.c',
      dependencies: [cpmdc_stub_dep, cmocka_dep],
      include_directories: inc,
      c_args: common_c_args,
    ),
    suite: ['stub'],
  )

  test(
    'capnp-parser',
    executable(
      'test_capnp_parser',
      'tests/test_capnp_parser.c',
      dependencies: [cmocka_dep] + cpmdc_params_deps,
      include_directories: [inc] + cpmdc_params_includes,
      sources: cpmdc_params_sources,
      c_args: common_c_args,
    ),
    args: [cpmd_parser_params_bin.full_path()],
    depends: [cpmd_parser_params_bin],
    suite: ['capnp'],
  )

  test(
    'force-input-size',
    executable(
      'test_force_input_size',
      'tests/test_force_input_size.c',
      dependencies: [cpmdc_embed_shell_dep, cmocka_dep],
      include_directories: [inc] + cpmdc_params_includes,
      c_args: common_c_args,
    ),
    args: [
      force_input_step_a_bin.full_path(),
      force_input_step_ev_bin.full_path(),
    ],
    depends: [force_input_step_a_bin, force_input_step_ev_bin],
    suite: ['capnp', 'socket'],
  )

  test(
    'session-socket',
    executable(
      'test_session_socket',
      'tests/test_session_socket.c',
      dependencies: [cpmdc_embed_shell_dep, cmocka_dep],
      include_directories: [inc] + cpmdc_params_includes,
      c_args: common_c_args,
    ),
    args: [
      cpmd_parser_params_bin.full_path(),
      force_input_step_a_bin.full_path(),
    ],
    depends: [cpmd_parser_params_bin, force_input_step_a_bin],
    suite: ['capnp', 'socket'],
  )

  test(
    'embed-shell-unavailable',
    executable(
      'test_embed_shell',
      'tests/test_embed_shell.c',
      dependencies: [cpmdc_embed_shell_dep, cmocka_dep],
      include_directories: inc,
      c_args: common_c_args,
    ),
    suite: ['embed'],
  )
endif
