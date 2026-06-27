Versioning
==========

``cpmdc`` uses semantic versioning on the public C ABI and Cap'n Proto
schema fields. Additive Cap'n Proto fields and new ABI entry points are
minor; breaking field renumbers or ABI signature changes are major.

Checklist
=========

#. ``meson test -C build --print-errorlogs`` on the stub/parser
   configuration.
#. Update ``CHANGELOG.md`` for user-visible changes.
#. Tag ``vX.Y.Z`` when publishing.

0.1 covers the Cap'n Proto structure, C ABI, ISO\ :sub:`C` embed shell,
stub/reference backend, and archive-linked OpenCPMD energy/force path.
