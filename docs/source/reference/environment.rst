Read by libcpmdc at run time
============================

+--------------------------+-------------------------------+----------------------+
| Variable                 | Read in                       | Effect               |
+==========================+===============================+======================+
| ``CPMDC_PSEUDO_DIR``     | ``src/cpmd_embed_c_api.F90``, | directory holding    |
|                          | first call of an OpenCPMD     | the pseudopotential  |
|                          | session                       | files named in the   |
|                          |                               | deck; ``cpmdc``      |
|                          |                               | changes into it      |
|                          |                               | while CPMD reads     |
|                          |                               | them, then returns   |
|                          |                               | to the host's        |
|                          |                               | working directory    |
+--------------------------+-------------------------------+----------------------+
| ``CPMD_PP_LIBRARY_PATH`` | same place, when              | fallback for the     |
|                          | ``CPMDC_PSEUDO_DIR`` is unset | same directory       |
+--------------------------+-------------------------------+----------------------+
| ``CPMDC_DECK_OUT``       | ``src/cpmdc.c``, each time a  | file that receives   |
|                          | message is rendered           | the rendered method  |
|                          |                               | deck; each write     |
|                          |                               | replaces the last    |
+--------------------------+-------------------------------+----------------------+
| ``CPMDC_STRESS``         | ``src/cpmd_embed_c_api.F90``, | ``0`` skips the      |
|                          | each OpenCPMD evaluation      | stress tensor on     |
|                          |                               | periodic cells; any  |
|                          |                               | other value, or      |
|                          |                               | unset, computes it   |
+--------------------------+-------------------------------+----------------------+

The first call of an OpenCPMD session fails when neither
``CPMDC_PSEUDO_DIR`` nor ``CPMD_PP_LIBRARY_PATH`` names an existing
directory. While ``cpmdc`` is in the pseudopotential directory it also
sets ``CPMD_PP_LIBRARY_PATH`` and ``PP_LIBRARY_PATH`` to that directory
with a trailing slash, for OpenCPMD's own lookup. The change of
directory is needed because OpenCPMD's ``get_pplib`` takes ``argv[2]``
as the library path whenever the process has more than one argument,
which most hosts do.

``CPMDC_STRESS=0`` saves the stress calculation for callers that use
only energy and forces, and the stress snapshot stays invalid. Isolated
cells never compute stress, whatever the variable says, and the snapshot
stays invalid for them even when the box volume is positive.

Read by the tests and tools
===========================

+----------------------+----------------------------------------------+------------------------+
| Variable             | Read by                                      | Effect                 |
+======================+==============================================+========================+
| ``CPMD_ROOT``, then  | ``tools/check_feature_inventory.py``         | OpenCPMD source tree   |
| ``CPMDC_CPMD_ROOT``  | (``feature-inventory`` test)                 | to probe for           |
|                      |                                              | ``inscan('&SECTION')`` |
|                      |                                              | calls; without it the  |
|                      |                                              | live completeness      |
|                      |                                              | check is skipped with  |
|                      |                                              | a warning              |
+----------------------+----------------------------------------------+------------------------+
| ``CPMD_SRC``, then   | ``tests/test_cpmd_inscan_token_fidelity.py`` | OpenCPMD source tree   |
| ``OPENCPMD_SRC``     |                                              | for the parser token   |
|                      |                                              | check                  |
+----------------------+----------------------------------------------+------------------------+
| ``CPMDC_CPMD_EXE``,  | ``tools/cpmd_root_probe.py``                 | the ``cpmd.x`` to      |
| then ``CPMD_ROOT``   |                                              | check, or the tree     |
|                      |                                              | whose ``bin/cpmd.x``   |
|                      |                                              | it checks              |
+----------------------+----------------------------------------------+------------------------+
| ``MAKE``, ``AR``     | ``tools/rebuild_opencpmd_embed.sh``          | the Make and archiver  |
|                      |                                              | programs; default      |
|                      |                                              | ``make`` and ``ar``    |
+----------------------+----------------------------------------------+------------------------+

The live tests of an OpenCPMD build (``meson test -C build-cpmd``) run
``libcpmdc``, so they need ``CPMDC_PSEUDO_DIR`` like any other host.

Read by hosts that load libcpmdc
================================

These belong to the host, not to ``cpmdc``, and are listed here because
they decide which ``libcpmdc`` a host loads.

+-------------------------+----------------------+----------------------------+
| Variable                | Host                 | Effect                     |
+=========================+======================+============================+
| ``CPMDC_LIBRARY``,      | ``rgpot``            | path of the                |
| ``RGPOT_CPMDC_ENGINE``, | ``CPMDPot``          | ``libcpmdc.so`` to         |
| ``RGPOT_CPMD_ENGINE``   |                      | ``dlopen``, tried in that  |
|                         |                      | order after an explicit    |
|                         |                      | engine path                |
+-------------------------+----------------------+----------------------------+
| ``RGPOT_PARAMS_PATH``   | eOn ``RgpotPot``     | overrides                  |
|                         | builds that read     | ``[RgpotPot] params_path`` |
|                         | ``params_path``      |                            |
+-------------------------+----------------------+----------------------------+
| ``OMP_NUM_THREADS``     | OpenMP runtime       | threads per rank inside    |
|                         |                      | CPMD when OpenCPMD was     |
|                         |                      | built with OpenMP          |
+-------------------------+----------------------+----------------------------+
