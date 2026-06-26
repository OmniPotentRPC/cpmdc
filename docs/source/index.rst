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

``cpmdc`` provides a stable C ABI for in-process
`OpenCPMD <https://github.com/OpenCPMD/CPMD>`__ calculations under the
OmniPotentRPC org (same role as ``nwchemc`` for NWChem). Callers pass an
unpacked flat Cap'n Proto message whose root is ``CPMDParams``. The
runtime reads that message through generated ``capnp-c`` bindings and
applies the resulting method configuration through a Fortran
``iso_c_binding`` / ``iso_fortran_env`` embed shell.

Geometry for each evaluation is a separate Cap'n Proto ``ForceInput``
message. Long-running drivers create one ``CPMDCSession`` from
``CPMDParams``, then pass ``ForceInput`` per step and receive an
unpacked flat ``PotentialResult`` (direct-call "socket" style shared
with ``nwchemc`` and ``rgpot``). There is no parallel TOML/JSON options
language for the backend.

Wire types (``ForceInput``, ``PotentialResult``, ``CPMDParams``, section
unions) are defined only in ``schema/Potentials.capnp`` and materialized
as generated C structs in ``Potentials.capnp.h``. The public ABI in
``include/cpmdc.h`` is the language-neutral carrier around those byte
buffers (``CPMDCResult``, opaque ``CPMDCSession``).

.. code:: bash

   git clone https://github.com/OmniPotentRPC/cpmdc.git
   cd cpmdc
   meson setup build -Dwith_tests=true
   meson compile -C build
   meson test -C build --print-errorlogs

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
