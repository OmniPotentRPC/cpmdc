# Installed next to libcpmdc_gromacs under <libdir>/cpmdc/gromacs/pkgconfig.
# The same relative_to rules as cpmdc.pc: libdir '.' and a libdir that
# contains '..' must not walk above the prefix, and an absolute includedir
# must not be copied verbatim.
gromacs_pc_dir = resolved_libdir / 'cpmdc' / 'gromacs' / 'pkgconfig'
gromacs_inc_dir = resolved_includedir / 'cpmdc' / 'gromacs'
cp2k_prefix_rel = fs.relative_to(install_prefix, gromacs_pc_dir).replace(
  '\\',
  '/',
)
cp2k_libdir_rel = fs.relative_to(resolved_libdir, gromacs_pc_dir).replace(
  '\\',
  '/',
)
cp2k_includedir_rel = fs.relative_to(gromacs_inc_dir, gromacs_pc_dir).replace(
  '\\',
  '/',
)
foreach rel : [cp2k_prefix_rel, cp2k_libdir_rel, cp2k_includedir_rel]
  if rel == '' or fs.is_absolute(rel) or rel.contains(':') or rel.contains('$')
    error(
      'cannot express a libcp2k pkg-config path from ' + gromacs_pc_dir + ' (' + rel + ')',
    )
  endif
endforeach
cp2k_prefix = cp2k_prefix_rel == '.' ? '${pcfiledir}' : '${pcfiledir}/' + cp2k_prefix_rel
cp2k_libdir = cp2k_libdir_rel == '.' ? '${pcfiledir}' : '${pcfiledir}/' + cp2k_libdir_rel
cp2k_includedir = cp2k_includedir_rel == '.' ? '${pcfiledir}' : '${pcfiledir}/' + cp2k_includedir_rel

foreach where : [
  [
    'libcp2k.pc',
    meson.project_build_root(),
    meson.project_source_root() / 'tools' / 'gromacs',
    meson.project_source_root(),
  ],
  [
    'libcp2k-installed.pc',
    cp2k_libdir,
    cp2k_includedir,
    cp2k_prefix,
  ],
]
  pc_conf = configuration_data()
  pc_conf.set('LIBDIR', where[1])
  pc_conf.set('INCLUDEDIR', where[2])
  pc_conf.set('PREFIX', where[3])
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
