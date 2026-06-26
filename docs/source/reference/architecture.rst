Layers
======

``cpmdc`` is deliberately small and mirrors ``nwchemc``:

-  ``include/cpmdc.h`` exposes the stable C ABI (``CPMDCResult``,
   ``CPMDCSession``, socket entry points).
-  ``schema/Potentials.capnp`` defines ``CPMDParams``, per-step
   ``ForceInput``, ``PotentialResult``, structured CPMD input sections,
   and ``PotentialConfig`` (``cpmd @2`` arm).
-  ``subprojects/c-capnproto`` provides the ``capnp-c`` runtime and
   compiler plugin.
-  ``src/cpmdc_params.c`` opens serialized Cap'n Proto messages, returns
   generated ``CPMDParams_ptr`` / ``ForceInput_ptr`` roots, renders
   ``&SECTION`` decks, and writes ``PotentialResult`` flat messages.
-  ``src/cpmdc.c`` owns sessions, topology checks, unit conversion for
   results, and calls into the Fortran embed shell.
-  ``src/cpmd_embed_c_api.f90`` owns the ``iso_c_binding`` /
   ``iso_fortran_env`` layer, saved method configuration, and a
   deterministic reference PEF for E2E tests until OpenCPMD archives are
   linked.
-  ``src/cpmdc_stub.c`` implements the same ABI symbols with
   ``cpmdc_available() == 0`` for frontend builds that do not link the
   embed shell.

Parameter Flow
==============

The ABI receives a byte pointer and byte count. ``capn_init_mem()``
opens the flat message, ``capn_getp(capn_root(...))`` obtains the
generated ``CPMDParams_ptr``, and generated ``read_CPMDParams()`` calls
are used as short-lived C views while the arena is alive. Session force
and result calls use the same pattern for generated ``ForceInput_ptr``
roots.

There is no second user-facing configuration struct in C. Structured
``CPMDInputSection`` unions (``cpmd``, ``system``, ``dft``, ``atoms``,
``generic``, ``raw``) are Cap'n Proto only; the renderer turns them into
CPMD-style ``&SECTION`` … ``&END`` text for embed configuration and
debugging.

Long-Running Sessions
=====================

``cpmdc_session_create()`` copies the ``CPMDParams`` bytes and renders
the input deck once. ``cpmdc_session_calculate_forces()`` and
``cpmdc_session_calculate_result()`` accept serialized ``ForceInput``
messages so MD, optimizers, and socket-style drivers can update
coordinates and the periodic cell without rebuilding the method carrier.
``cpmdc_potential_result_size_for_force_input()`` exposes result buffer
sizing without touching the embed runtime. ``cpmdc_calculate_result()``
creates a temporary session for one-shot callers.

The first accepted session evaluation records topology (atom count and
ordered atomic numbers). Session calls reject later steps with a
different topology. ``cpmdc_session_set_params`` may replace
configuration only before topology is accepted.

Branding
========

Wordmark and mark assets live under ``branding/logo/`` and
``docs/source/_static/`` (light/dark SVG, matching ``nwchemc`` docs
layout).
