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
dl_dep = dependency('dl', required: false)
if not dl_dep.found()
  dl_dep = cc.find_library('dl', required: false)
endif
cmocka_dep = dependency('cmocka', required: get_option('with_tests'))

inc = include_directories('include', 'src')
cpmd_root = get_option('cpmd_root')
has_cpmd = get_option('with_cpmd') and cpmd_root != ''
need_capnp_c = true

common_c_args = []
capnp_generated_c_args = []
if cc.get_id() != 'msvc'
  common_c_args += ['-Wall', '-Wextra', '-Wpedantic']
endif
if cc.get_id() in ['gcc', 'clang']
  capnp_generated_c_args += ['-Wno-missing-field-initializers']
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
  ['src/cpmdc_stub.c', 'src/cpmdc_features.c'],
  include_directories: inc,
  c_args: common_c_args + capnp_generated_c_args,
  install: false,
)

cpmdc_stub_dep = declare_dependency(
  link_with: cpmdc_stub,
  include_directories: inc,
)

# ISO_C embed without libcpmd: deterministic reference PEF for tests and dlopen users.
cpmd_embed_f_args = []
cpmd_embed_link = []
cpmd_embed_inc = [inc]
cpmd_extra_objects = []
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
  cpmd_timetag_obj = cpmd_obj / 'timetag.o'
  timetag_check = run_command(
    python,
    [
      '-c',
      'import pathlib, sys; p = pathlib.Path(sys.argv[1]); sys.exit(0 if p.is_file() else 1)',
      cpmd_timetag_obj,
    ],
    check: false,
  )
  if timetag_check.returncode() == 0
    # OpenCPMD links timetag outside libcpmd.a for executable build stamps.
    cpmd_extra_objects += [cpmd_timetag_obj]
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

cpmd_embed_sources = ['src/cpmdc.c', 'src/cpmd_embed_c_api.F90', 'src/cpmd_stop_intercept.c', 'src/cpmdc_features.c'] + cpmdc_params_sources
cpmd_link_multi = has_cpmd ? ['-Wl,--allow-multiple-definition'] : []
cpmdc_shared_link_args = []
if has_cpmd
  cpmdc_shared_link_args = cpmd_link_multi + ['-Wl,--start-group', cpmd_root / 'lib' / 'libcpmd.a', '-Wl,--end-group'] \
      + (openblas.found() ? ['-lopenblas'] : ['-llapack', '-lblas']) \
      + ['-lgfortran', '-lm']
endif

cpmdc_embed_shell = static_library(
  'cpmdc_embed_shell',
  cpmd_embed_sources,
  include_directories: cpmd_embed_inc + cpmdc_params_includes,
  dependencies: cpmdc_params_deps + [math_dep],
  c_args: common_c_args + capnp_generated_c_args,
  fortran_args: cpmd_embed_f_args,
  link_args: has_cpmd ? [] : [],
  install: false,
)

cpmdc_shared = shared_library(
  'cpmdc',
  cpmd_embed_sources,
  include_directories: cpmd_embed_inc + cpmdc_params_includes,
  dependencies: cpmdc_params_deps + [math_dep],
  c_args: common_c_args + capnp_generated_c_args,
  fortran_args: cpmd_embed_f_args,
  objects: cpmd_extra_objects,
  link_args: cpmdc_shared_link_args,
  install: true,
)

cpmdc_shared_dep = declare_dependency(
  link_with: cpmdc_shared,
  include_directories: inc,
  dependencies: cpmdc_params_deps,
)

install_headers('include/cpmdc.h', 'include/cpmdc_features.h')

cpmdc_embed_shell_dep = declare_dependency(
  link_with: cpmdc_embed_shell,
  include_directories: inc,
  dependencies: cpmdc_params_deps,
  link_args: has_cpmd ? (cpmd_extra_objects + cpmd_link_multi + ['-Wl,--start-group', cpmd_root / 'lib' / 'libcpmd.a', '-Wl,--end-group', '-lopenblas', '-lgfortran']) : [],
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
  cpmd_session_alt_params_bin = custom_target(
    'cpmd_session_alt_params_bin',
    input: 'tests/cpmd_params_session_alt.capnp.txt',
    output: 'cpmd_params_session_alt.bin',
    command: encode + ['CPMDParams', '@INPUT@', '@OUTPUT@'],
  )
  cpmd_set_params_bin = custom_target(
    'cpmd_set_params_bin',
    input: 'tests/params_set_directives.capnp.txt',
    output: 'params_set_directives.bin',
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
  water_params_bin = custom_target(
    'water_params_bin',
    input: 'tests/water_blyp_params.capnp.txt',
    output: 'water_blyp_params.bin',
    command: encode + ['CPMDParams', '@INPUT@', '@OUTPUT@'],
  )
  water_force_bin = custom_target(
    'water_force_bin',
    input: 'tests/force_input_water.capnp.txt',
    output: 'force_input_water.bin',
    command: encode + ['ForceInput', '@INPUT@', '@OUTPUT@'],
  )
  effective_top_params_bin = custom_target(
    'effective_top_params_bin',
    input: 'tests/params_effective_config_top.capnp.txt',
    output: 'params_effective_config_top.bin',
    command: encode + ['CPMDParams', '@INPUT@', '@OUTPUT@'],
  )
  effective_section_params_bin = custom_target(
    'effective_section_params_bin',
    input: 'tests/params_effective_config_sections.capnp.txt',
    output: 'params_effective_config_sections.bin',
    command: encode + ['CPMDParams', '@INPUT@', '@OUTPUT@'],
  )
  atoms_silicon_params_bin = custom_target(
    'atoms_silicon_params_bin',
    input: 'tests/params_atoms_silicon.capnp.txt',
    output: 'params_atoms_silicon.bin',
    command: encode + ['CPMDParams', '@INPUT@', '@OUTPUT@'],
  )
  atoms_default_directives_params_bin = custom_target(
    'atoms_default_directives_params_bin',
    input: 'tests/params_atoms_default_directives.capnp.txt',
    output: 'params_atoms_default_directives.bin',
    command: encode + ['CPMDParams', '@INPUT@', '@OUTPUT@'],
  )
  runtime_params_bin = has_cpmd ? water_params_bin : cpmd_parser_params_bin

  test(
    'stub-abi',
    executable(
      'test_stub_abi',
      'tests/test_stub_abi.c',
      dependencies: [cpmdc_stub_dep, cmocka_dep],
      include_directories: inc,
      c_args: common_c_args + capnp_generated_c_args,
    ),
    suite: ['stub'],
  )

  test(
    'public-header-features',
    executable(
      'test_public_header_features',
      'tests/test_public_header_features.c',
      dependencies: [cpmdc_embed_shell_dep],
      include_directories: inc,
      c_args: common_c_args,
    ),
    suite: ['abi'],
  )

  test(
    'params-render-cmocka',
    executable(
      'test_params_render_cmocka',
      'tests/cmocka/test_params_render_cmocka.c',
      dependencies: [cmocka_dep] + cpmdc_params_deps,
      include_directories: [inc] + cpmdc_params_includes,
      sources: cpmdc_params_sources,
      c_args: common_c_args + capnp_generated_c_args,
    ),
    args: [cpmd_parser_params_bin.full_path()],
    depends: [cpmd_parser_params_bin],
    suite: ['cmocka', 'capnp'],
  )

  test(
    'set-directives-render-cmocka',
    executable(
      'test_set_directives_render',
      'tests/cmocka/test_set_directives_render.c',
      dependencies: [cmocka_dep] + cpmdc_params_deps,
      include_directories: [inc] + cpmdc_params_includes,
      sources: cpmdc_params_sources + ['src/cpmdc_features.c'],
      c_args: common_c_args + capnp_generated_c_args,
    ),
    args: [cpmd_set_params_bin.full_path()],
    depends: [cpmd_set_params_bin],
    suite: ['cmocka', 'capnp', 'params'],
  )

  test(
    'force-input-size',
    executable(
      'test_force_input_size',
      'tests/test_force_input_size.c',
      dependencies: [cpmdc_embed_shell_dep, cmocka_dep],
      include_directories: [inc] + cpmdc_params_includes,
      c_args: common_c_args + capnp_generated_c_args,
    ),
    args: [
      force_input_step_a_bin.full_path(),
      force_input_step_ev_bin.full_path(),
    ],
    depends: [force_input_step_a_bin, force_input_step_ev_bin],
    suite: ['capnp', 'socket'],
  )

  test(
    'effective-config',
    executable(
      'test_effective_config',
      'tests/test_effective_config.c',
      dependencies: [cpmdc_embed_shell_dep, cmocka_dep, math_dep],
      include_directories: inc,
      c_args: common_c_args + capnp_generated_c_args,
    ),
    args: [
      effective_top_params_bin.full_path(),
      effective_section_params_bin.full_path(),
      force_input_step_a_bin.full_path(),
    ],
    depends: [
      effective_top_params_bin,
      effective_section_params_bin,
      force_input_step_a_bin,
    ],
    suite: ['e2e', 'cmocka', 'params'],
  )

  test(
    'atoms-symbols',
    executable(
      'test_atoms_symbols',
      'tests/test_atoms_symbols.c',
      dependencies: [cmocka_dep] + cpmdc_params_deps,
      include_directories: [inc] + cpmdc_params_includes,
      sources: cpmdc_params_sources,
      c_args: common_c_args + capnp_generated_c_args,
    ),
    args: [
      atoms_silicon_params_bin.full_path(),
      atoms_default_directives_params_bin.full_path(),
    ],
    depends: [
      atoms_silicon_params_bin,
      atoms_default_directives_params_bin,
    ],
    suite: ['cmocka', 'params'],
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
      runtime_params_bin.full_path(),
      force_input_step_a_bin.full_path(),
      force_input_step_ev_bin.full_path(),
      force_input_step_species_bin.full_path(),
      cpmd_session_alt_params_bin.full_path(),
    ],
    depends: [
      runtime_params_bin,
      force_input_step_a_bin,
      force_input_step_ev_bin,
      force_input_step_species_bin,
      cpmd_session_alt_params_bin,
    ],
    suite: ['capnp', 'socket'],
    timeout: 300,
  )
  example_host_step = executable(
    'example_host_step',
    'examples/host_step.c',
    dependencies: [cpmdc_embed_shell_dep],
    include_directories: inc,
    c_args: common_c_args + capnp_generated_c_args,
  )
  test(
    'example-host-step',
    example_host_step,
    args: [
      runtime_params_bin.full_path(),
      force_input_step_a_bin.full_path(),
    ],
    depends: [runtime_params_bin, force_input_step_a_bin],
    suite: ['examples', 'e2e'],
    timeout: 300,
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
  test(
    'shared-dlopen-abi',
    executable(
      'test_shared_dlopen_abi',
      'tests/test_shared_dlopen_abi.c',
      dependencies: [dl_dep],
      include_directories: inc,
      c_args: common_c_args,
    ),
    args: [cpmdc_shared.full_path()],
    depends: [cpmdc_shared],
    suite: ['embed', 'dlopen'],
  )

  # E2E suites use the deterministic reference PEF or the linked OpenCPMD runtime.
  e2e_single_point_exe = executable(
    'test_e2e_single_point',
    'tests/test_e2e_single_point.c',
    dependencies: [cpmdc_embed_shell_dep, cmocka_dep, math_dep],
    include_directories: [inc] + cpmdc_params_includes,
    c_args: common_c_args,
  )
  test(
    'e2e-single-point',
    e2e_single_point_exe,
    args: [
      runtime_params_bin.full_path(),
      force_input_step_a_bin.full_path(),
    ],
    depends: [runtime_params_bin, force_input_step_a_bin],
    suite: ['e2e', 'cmocka'],
    timeout: 300,
  )
  test(
    'e2e-set-directives-session',
    e2e_single_point_exe,
    args: [
      cpmd_set_params_bin.full_path(),
      force_input_step_a_bin.full_path(),
    ],
    depends: [cpmd_set_params_bin, force_input_step_a_bin],
    suite: ['e2e', 'cmocka', 'params'],
    timeout: 300,
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
      runtime_params_bin.full_path(),
      force_input_step_a_bin.full_path(),
      force_input_step_ev_bin.full_path(),
      force_input_step_species_bin.full_path(),
    ],
    depends: [
      runtime_params_bin,
      force_input_step_a_bin,
      force_input_step_ev_bin,
      force_input_step_species_bin,
    ],
    suite: ['e2e', 'cmocka'],
    timeout: 300,
  )

  params_escape_bin = custom_target(
    'params_escape_bin',
    input: 'tests/params_escape_generic.capnp.txt',
    output: 'params_escape_generic.bin',
    command: encode + ['CPMDParams', '@INPUT@', '@OUTPUT@'],
  )
  test(
    'feature-inventory',
    python,
    args: [meson.current_source_dir() / 'tools/check_feature_inventory.py'],
    suite: ['inventory'],
  )
  test(
    'examples-documented',
    python,
    args: [
      meson.current_source_dir() / 'tests/test_examples_documented.py',
      meson.current_source_dir(),
    ],
    suite: ['docs', 'examples'],
  )
  test(
    'params-escape-render',
    executable(
      'test_params_escape_render',
      'tests/test_params_escape_render.c',
      dependencies: [cpmdc_embed_shell_dep, cmocka_dep],
    ),
    args: [params_escape_bin.full_path()],
    depends: [params_escape_bin],
    suite: ['params'],
  )
  if has_cpmd
    test(
      'water-capnp-cli-parity',
      executable(
        'test_water_capnp_cli_parity',
        'tests/test_water_capnp_cli_parity.c',
        dependencies: [cpmdc_embed_shell_dep, math_dep],
      ),
      args: [water_params_bin.full_path(), water_force_bin.full_path()],
      depends: [water_params_bin, water_force_bin],
      suite: ['e2e', 'water'],
      timeout: 600,
    )
  endif


  test(
    'features-cmocka',
    executable(
      'test_features_cmocka',
      ['tests/test_features_cmocka.c', 'src/cpmdc_features.c'],
      dependencies: [cmocka_dep],
      include_directories: inc,
      c_args: common_c_args,
    ),
    suite: ['inventory', 'cmocka'],
  )
  water_params_70_txt = custom_target(
    'water_params_70_txt',
    output: 'water_blyp_params_70.capnp.txt',
    command: [python, '-c',
      'import pathlib,sys; t=pathlib.Path(sys.argv[1]).read_text().replace("50.0","70.0"); pathlib.Path(sys.argv[2]).write_text(t)',
      meson.current_source_dir() / 'tests/water_blyp_params.capnp.txt', '@OUTPUT@'],
  )
  water_params_70_bin = custom_target(
    'water_params_70_bin',
    input: water_params_70_txt,
    output: 'water_blyp_params_70.bin',
    command: encode + ['CPMDParams', '@INPUT@', '@OUTPUT@'],
  )
  test(
    'capnp-drives-config',
    executable(
      'test_capnp_drives_config',
      'tests/test_capnp_drives_config.c',
      dependencies: [cmocka_dep] + cpmdc_params_deps,
      include_directories: [inc] + cpmdc_params_includes,
      sources: cpmdc_params_sources,
      c_args: common_c_args,
    ),
    args: [water_params_bin.full_path(), water_params_70_bin.full_path()],
    depends: [water_params_bin, water_params_70_bin],
    suite: ['params', 'capnp'],
  )


  params_cp_md_bin = custom_target(
    'params_cp_md_bin',
    input: 'tests/params_cp_md.capnp.txt',
    output: 'params_cp_md.bin',
    command: encode + ['CPMDParams', '@INPUT@', '@OUTPUT@'],
  )
  params_dft_func_bin = custom_target(
    'params_dft_func_bin',
    input: 'tests/params_dft_functionals.capnp.txt',
    output: 'params_dft_functionals.bin',
    command: encode + ['CPMDParams', '@INPUT@', '@OUTPUT@'],
  )
  params_cpmd_geometry_bin = custom_target(
    'params_cpmd_geometry_bin',
    input: 'tests/params_cpmd_geometry.capnp.txt',
    output: 'params_cpmd_geometry.bin',
    command: encode + ['CPMDParams', '@INPUT@', '@OUTPUT@'],
  )
  params_dft_scalars_bin = custom_target(
    'params_dft_scalars_bin',
    input: 'tests/params_dft_scalars.capnp.txt',
    output: 'params_dft_scalars.bin',
    command: encode + ['CPMDParams', '@INPUT@', '@OUTPUT@'],
  )
  params_cpmd_dynamics_bin = custom_target(
    'params_cpmd_dynamics_bin',
    input: 'tests/params_cpmd_dynamics.capnp.txt',
    output: 'params_cpmd_dynamics.bin',
    command: encode + ['CPMDParams', '@INPUT@', '@OUTPUT@'],
  )
  params_cpmd_misc_bin = custom_target(
    'params_cpmd_misc_bin',
    input: 'tests/params_cpmd_misc.capnp.txt',
    output: 'params_cpmd_misc.bin',
    command: encode + ['CPMDParams', '@INPUT@', '@OUTPUT@'],
  )
  params_long_tail_sections_bin = custom_target(
    'params_long_tail_sections_bin',
    input: 'tests/params_long_tail_sections.capnp.txt',
    output: 'params_long_tail_sections.bin',
    command: encode + ['CPMDParams', '@INPUT@', '@OUTPUT@'],
  )
  params_vdw_controls_bin = custom_target(
    'params_vdw_controls_bin',
    input: 'tests/params_vdw_controls.capnp.txt',
    output: 'params_vdw_controls.bin',
    command: encode + ['CPMDParams', '@INPUT@', '@OUTPUT@'],
  )
  params_system_controls_bin = custom_target(
    'params_system_controls_bin',
    input: 'tests/params_system_controls.capnp.txt',
    output: 'params_system_controls.bin',
    command: encode + ['CPMDParams', '@INPUT@', '@OUTPUT@'],
  )
  params_system_kpoints_monkhorst_bin = custom_target(
    'params_system_kpoints_monkhorst_bin',
    input: 'tests/params_system_kpoints_monkhorst.capnp.txt',
    output: 'params_system_kpoints_monkhorst.bin',
    command: encode + ['CPMDParams', '@INPUT@', '@OUTPUT@'],
  )
  params_system_kpoints_bands_bin = custom_target(
    'params_system_kpoints_bands_bin',
    input: 'tests/params_system_kpoints_bands.capnp.txt',
    output: 'params_system_kpoints_bands.bin',
    command: encode + ['CPMDParams', '@INPUT@', '@OUTPUT@'],
  )
  params_system_cdft_acceptor_wmult_bin = custom_target(
    'params_system_cdft_acceptor_wmult_bin',
    input: 'tests/params_system_cdft_acceptor_wmult.capnp.txt',
    output: 'params_system_cdft_acceptor_wmult.bin',
    command: encode + ['CPMDParams', '@INPUT@', '@OUTPUT@'],
  )
  params_system_couplings_prod_bin = custom_target(
    'params_system_couplings_prod_bin',
    input: 'tests/params_system_couplings_prod.capnp.txt',
    output: 'params_system_couplings_prod.bin',
    command: encode + ['CPMDParams', '@INPUT@', '@OUTPUT@'],
  )
  params_system_couplings_linres_bin = custom_target(
    'params_system_couplings_linres_bin',
    input: 'tests/params_system_couplings_linres.capnp.txt',
    output: 'params_system_couplings_linres.bin',
    command: encode + ['CPMDParams', '@INPUT@', '@OUTPUT@'],
  )
  params_system_couplings_lists_bin = custom_target(
    'params_system_couplings_lists_bin',
    input: 'tests/params_system_couplings_lists.capnp.txt',
    output: 'params_system_couplings_lists.bin',
    command: encode + ['CPMDParams', '@INPUT@', '@OUTPUT@'],
  )
  params_system_cell_qualifiers_bin = custom_target(
    'params_system_cell_qualifiers_bin',
    input: 'tests/params_system_cell_qualifiers.capnp.txt',
    output: 'params_system_cell_qualifiers.bin',
    command: encode + ['CPMDParams', '@INPUT@', '@OUTPUT@'],
  )
  params_system_cell_vectors_bin = custom_target(
    'params_system_cell_vectors_bin',
    input: 'tests/params_system_cell_vectors.capnp.txt',
    output: 'params_system_cell_vectors.bin',
    command: encode + ['CPMDParams', '@INPUT@', '@OUTPUT@'],
  )
  params_system_bad_occupation_bin = custom_target(
    'params_system_bad_occupation_bin',
    input: 'tests/params_system_bad_occupation.capnp.txt',
    output: 'params_system_bad_occupation.bin',
    command: encode + ['CPMDParams', '@INPUT@', '@OUTPUT@'],
  )
  test(
    'params-cp-dft-render',
    executable(
      'test_params_cp_dft_render',
      ['tests/test_params_cp_dft_render.c', 'src/cpmdc_features.c'],
      dependencies: [cmocka_dep] + cpmdc_params_deps,
      include_directories: [inc] + cpmdc_params_includes,
      sources: cpmdc_params_sources,
      c_args: common_c_args,
    ),
    args: [
      params_cp_md_bin.full_path(),
      params_dft_func_bin.full_path(),
      params_cpmd_geometry_bin.full_path(),
      params_dft_scalars_bin.full_path(),
      params_cpmd_dynamics_bin.full_path(),
      params_cpmd_misc_bin.full_path(),
      params_long_tail_sections_bin.full_path(),
      params_vdw_controls_bin.full_path(),
      params_system_controls_bin.full_path(),
      params_system_kpoints_monkhorst_bin.full_path(),
      params_system_kpoints_bands_bin.full_path(),
      params_system_cdft_acceptor_wmult_bin.full_path(),
      params_system_couplings_prod_bin.full_path(),
      params_system_couplings_linres_bin.full_path(),
      params_system_couplings_lists_bin.full_path(),
      params_system_cell_qualifiers_bin.full_path(),
      params_system_cell_vectors_bin.full_path(),
      params_system_bad_occupation_bin.full_path(),
    ],
    depends: [
      params_cp_md_bin,
      params_dft_func_bin,
      params_cpmd_geometry_bin,
      params_dft_scalars_bin,
      params_cpmd_dynamics_bin,
      params_cpmd_misc_bin,
      params_long_tail_sections_bin,
      params_vdw_controls_bin,
      params_system_controls_bin,
      params_system_kpoints_monkhorst_bin,
      params_system_kpoints_bands_bin,
      params_system_cdft_acceptor_wmult_bin,
      params_system_couplings_prod_bin,
      params_system_couplings_linres_bin,
      params_system_couplings_lists_bin,
      params_system_cell_qualifiers_bin,
      params_system_cell_vectors_bin,
      params_system_bad_occupation_bin,
    ],
    suite: ['params', 'inventory'],
  )

endif
