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

# ISO_C embed without libcpmd: available()==0 (frontend / parser / sizing tests).
cpmd_embed_f_args = []
cpmd_embed_link = []
cpmd_embed_inc = [inc]
if has_cpmd
  cpmd_lib = cpmd_root / 'lib' / 'libcpmd.a'
  cpmd_obj = cpmd_root / 'obj'
  cpmd_src = cpmd_root / 'src'
  lib_check = run_command(
    python,
    [
      '-c',
      'import pathlib, sys; p = pathlib.Path(sys.argv[1]); sys.exit(0 if p.is_file() else 1)',
      cpmd_lib,
    ],
    check: false,
  )
  if lib_check.returncode() != 0
    error('with_cpmd requires ' + cpmd_lib + ' (build OpenCPMD in cpmd_root first)')
  endif
  cpmd_embed_f_args += [
    '-DCPMDC_HAS_CPMD=1',
    '-I' + cpmd_src,
    '-I' + cpmd_obj,
  ]
  cpmd_embed_link += [cpmd_lib]
  # LAPACK/BLAS as used by typical CPMD serial builds
  openblas = cc.find_library('openblas', required: false)
  if openblas.found()
    cpmd_embed_link += [openblas]
  else
    lapack = cc.find_library('lapack', required: false)
    blas = cc.find_library('blas', required: false)
    if lapack.found()
      cpmd_embed_link += [lapack]
    endif
    if blas.found()
      cpmd_embed_link += [blas]
    endif
  endif
  cpmd_embed_link += [math_dep, fc.find_library('gfortran', required: false)]
endif

cpmdc_embed_shell = static_library(
  'cpmdc_embed_shell',
  ['src/cpmdc.c', 'src/cpmd_embed_c_api.F90'] + cpmdc_params_sources,
  include_directories: cpmd_embed_inc + cpmdc_params_includes,
  dependencies: cpmdc_params_deps + [math_dep],
  c_args: common_c_args,
  fortran_args: cpmd_embed_f_args,
  link_args: has_cpmd ? [] : [],
  install: false,
)

# When HAS_CPMD, also produce installable shared library linked to libcpmd.a
if has_cpmd
  cpmdc = shared_library(
    'cpmdc',
    ['src/cpmdc.c', 'src/cpmd_embed_c_api.F90'] + cpmdc_params_sources,
    include_directories: cpmd_embed_inc + cpmdc_params_includes,
    dependencies: cpmdc_params_deps + [math_dep],
    c_args: common_c_args,
    fortran_args: cpmd_embed_f_args,
    link_whole: [],
    link_args: ['-Wl,--start-group', cpmd_root / 'lib' / 'libcpmd.a', '-Wl,--end-group']
      + (openblas.found() ? ['-lopenblas'] : ['-llapack', '-lblas'])
      + ['-lgfortran', '-lm'],
    install: true,
  )
  install_headers('include/cpmdc.h')
endif

cpmdc_embed_shell_dep = declare_dependency(
  link_with: cpmdc_embed_shell,
  include_directories: inc,
  dependencies: cpmdc_params_deps,
  link_args: has_cpmd ? ['-Wl,--start-group', cpmd_root / 'lib' / 'libcpmd.a', '-Wl,--end-group', '-lopenblas', '-lgfortran'] : [],
)

if get_option('with_tests')
  encode = [
    python,
    meson.current_source_dir() / 'tests/encode_capnp.py',
    capnp,
    meson.current_source_dir() / 'schema/Potentials.capnp',
  ]

  cpmd_parser_params_bin = custom_target(
    'cpmd_parser_params_bin',
    input: 'tests/cpmd_params_parser.capnp.txt',
    output: 'cpmd_params_parser.bin',
    command: encode + ['CPMDParams', '@INPUT@', '@OUTPUT@'],
  )
  force_input_step_a_bin = custom_target(
    'force_input_step_a_bin',
    input: 'tests/force_input_step_a.capnp.txt',
    output: 'force_input_step_a.bin',
    command: encode + ['ForceInput', '@INPUT@', '@OUTPUT@'],
  )
  force_input_step_ev_bin = custom_target(
    'force_input_step_ev_bin',
    input: 'tests/force_input_step_ev.capnp.txt',
    output: 'force_input_step_ev.bin',
    command: encode + ['ForceInput', '@INPUT@', '@OUTPUT@'],
  )
  force_input_step_species_bin = custom_target(
    'force_input_step_species_bin',
    input: 'tests/force_input_step_changed_species.capnp.txt',
    output: 'force_input_step_changed_species.bin',
    command: encode + ['ForceInput', '@INPUT@', '@OUTPUT@'],
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
    'params-render-cmocka',
    executable(
      'test_params_render_cmocka',
      'tests/cmocka/test_params_render_cmocka.c',
      dependencies: [cmocka_dep] + cpmdc_params_deps,
      include_directories: [inc] + cpmdc_params_includes,
      sources: cpmdc_params_sources,
      c_args: common_c_args,
    ),
    args: [cpmd_parser_params_bin.full_path()],
    depends: [cpmd_parser_params_bin],
    suite: ['cmocka', 'capnp'],
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
    'embed-shell-library',
    executable(
      'test_embed_shell',
      'tests/test_embed_shell.c',
      dependencies: [cpmdc_embed_shell_dep, cmocka_dep],
      include_directories: inc,
      c_args: common_c_args,
    ),
    suite: ['embed'],
  )

  # E2E suites need libcpmd.a (with_cpmd). Without it they skip at runtime.
  test(
    'e2e-single-point',
    executable(
      'test_e2e_single_point',
      'tests/test_e2e_single_point.c',
      dependencies: [cpmdc_embed_shell_dep, cmocka_dep, math_dep],
      include_directories: [inc] + cpmdc_params_includes,
      c_args: common_c_args,
    ),
    args: [
      cpmd_parser_params_bin.full_path(),
      force_input_step_a_bin.full_path(),
    ],
    depends: [cpmd_parser_params_bin, force_input_step_a_bin],
    suite: ['e2e', 'cmocka'],
  )

  test(
    'e2e-optimizer-session',
    executable(
      'test_e2e_optimizer_session',
      'tests/test_e2e_optimizer_session.c',
      dependencies: [cpmdc_embed_shell_dep, cmocka_dep, math_dep],
      include_directories: inc,
      c_args: common_c_args,
    ),
    args: [
      cpmd_parser_params_bin.full_path(),
      force_input_step_a_bin.full_path(),
      force_input_step_ev_bin.full_path(),
      force_input_step_species_bin.full_path(),
    ],
    depends: [
      cpmd_parser_params_bin,
      force_input_step_a_bin,
      force_input_step_ev_bin,
      force_input_step_species_bin,
    ],
    suite: ['e2e', 'cmocka'],
  )
endif
