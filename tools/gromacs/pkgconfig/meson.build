foreach where : [
  [
    'libcp2k.pc',
    meson.project_build_root(),
    meson.project_source_root() / 'tools' / 'gromacs',
  ],
  [
    'libcp2k-installed.pc',
    get_option('prefix') / get_option('libdir'),
    get_option('prefix') / get_option('includedir') / 'cpmdc' / 'gromacs',
  ],
]
  pc_conf = configuration_data()
  pc_conf.set('LIBDIR', where[1])
  pc_conf.set('INCLUDEDIR', where[2])
  pc_conf.set('VERSION', meson.project_version())
  pc_file = configure_file(
    input: 'libcp2k.pc.in',
    output: where[0],
    configuration: pc_conf,
  )
endforeach
install_data(
  pc_file,
  rename: 'libcp2k.pc',
  install_dir: get_option('libdir') / 'cpmdc' / 'gromacs' / 'pkgconfig',
)
