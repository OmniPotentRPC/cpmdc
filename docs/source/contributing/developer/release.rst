Versioning
==========

``cpmdc`` uses semantic versioning on the public C ABI and Cap'n Proto
schema fields. Additive Cap'n Proto fields and new ABI entry points are
minor; breaking field renumbers or ABI signature changes are major.

Checklist
=========

#. Run ``meson test -C build --print-errorlogs`` on the default
   parser/session configuration.
#. Run ``pixi run -e docs sphinxbld`` after documentation or public
   header changes.
#. Run the OpenCPMD archive suite when the release changes the embed
   path:
   ``CPMDC_PSEUDO_DIR=/path/to/PP_LIBRARY meson test -C build-cpmd --print-errorlogs``.
#. Confirm ``cpmdc_feature_table()`` exposes any new ABI, schema,
   section, or keyword capability.
#. Update ``CHANGELOG.md`` for user-visible changes.
#. Tag ``vX.Y.Z`` when publishing.

0.1 covers the Cap'n Proto structure, C ABI, ISO\ :sub:`C` embed shell,
stub/reference backend, and archive-linked OpenCPMD energy/force path.
