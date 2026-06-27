.. image:: /_static/logo-light.svg
   :class: light-only
   :width: 460
   :align: center

.. image:: /_static/logo-dark.svg
   :class: dark-only
   :width: 460
   :align: center

Overview
========

``cpmdc`` is the C ABI layer for embedding
`OpenCPMD <https://github.com/OpenCPMD/CPMD>`__ in OmniPotentRPC tools.
It keeps the public boundary small and language-neutral:

- method setup is an unpacked flat Cap'n Proto ``CPMDParams`` message
- per-step geometry is an unpacked flat Cap'n Proto ``ForceInput``
  message
- per-step output can be returned as native C values or as an unpacked
  flat ``PotentialResult`` message

The ABI does not expose C++ types, Rust types, or a second JSON/TOML
options language. All portable configuration lives in
``schema/Potentials.capnp``; the C library decodes those bytes with
generated ``capnp-c`` readers and renders the CPMD ``INPUT`` deck
internally.

Long-running drivers create one ``CPMDCSession`` from ``CPMDParams``,
then pass ``ForceInput`` per step. The result-carrier call,
``cpmdc_session_calculate_result()``, writes the same
``PotentialResult`` message used by ``rgpot`` RPC clients.

The first accepted session evaluation fixes atom count and ordered
atomic numbers. Later steps may change coordinates, units, and the 3x3
cell; topology changes require a new session.

.. code:: bash

   git clone https://github.com/OmniPotentRPC/cpmdc.git
   cd cpmdc
   meson setup build -Dwith_tests=true
   meson compile -C build
   meson test -C build --print-errorlogs

The default build includes a deterministic reference evaluator so the
ABI, parser, session, unit-conversion, and result-buffer contracts can
be tested without an OpenCPMD checkout. An OpenCPMD archive build links
the same C ABI surface against ``libcpmd.a``.

.. toctree::
   :maxdepth: 2
   :caption: Tutorials

   tutorials/quickstart

.. toctree::
   :maxdepth: 2
   :caption: How-To Guides

   howto/embedding

.. toctree::
   :maxdepth: 2
   :caption: Reference

   reference/architecture
   reference/cpmd-options

.. toctree::
   :maxdepth: 2
   :caption: Contributing

   contributing/index

.. toctree::
   :maxdepth: 1
   :caption: Project

   changelog
