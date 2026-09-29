.. index:: pair: struct; CPMDCEnergyComponents
.. _doxid-struct_c_p_m_d_c_energy_components:

struct CPMDCEnergyComponents
============================

.. toctree::
	:hidden:

Overview
~~~~~~~~

In-process snapshot of OpenCPMD ``ener_com`` scalars (Hartree a.u.). :ref:`More...<details-struct_c_p_m_d_c_energy_components>`


.. ref-code-block:: cpp
	:class: doxyrest-overview-code-block

	#include <cpmdc.h>
	
	struct CPMDCEnergyComponents {
		// fields
	
		int :ref:`valid<doxid-struct_c_p_m_d_c_energy_components_1ac6f6f27d3686bcdb5ffcbe470c33f089>`;
		double :target:`etot<doxid-struct_c_p_m_d_c_energy_components_1accbae2147ca6fb5aa149f2152b17d4c7>`;
		double :target:`ekin<doxid-struct_c_p_m_d_c_energy_components_1a28b4de9813ba1a996400d3667b8fcf2c>`;
		double :target:`epseu<doxid-struct_c_p_m_d_c_energy_components_1aa2a92a5dd822f90a86b650a2b89e7412>`;
		double :target:`enl<doxid-struct_c_p_m_d_c_energy_components_1a4c3da6c51c6c813a7baeed853bf40efe>`;
		double :target:`eht<doxid-struct_c_p_m_d_c_energy_components_1ac45c82cb8dccbf84940d23f791f017fe>`;
		double :target:`ehep<doxid-struct_c_p_m_d_c_energy_components_1ad8693b35625ee55d17cda05ea74e67b4>`;
		double :target:`ehee<doxid-struct_c_p_m_d_c_energy_components_1a3d2f2d2507e2fe0784b4d988efce489c>`;
		double :target:`ehii<doxid-struct_c_p_m_d_c_energy_components_1ac10cff885b025886efffb99775fc4cc2>`;
		double :target:`exc<doxid-struct_c_p_m_d_c_energy_components_1ae18cdd00a8bb4f77fea0d9ed23e4ae2e>`;
		double :target:`vxc<doxid-struct_c_p_m_d_c_energy_components_1a7ceb934333e4f5d1f4623c2b148707fa>`;
		double :target:`egc<doxid-struct_c_p_m_d_c_energy_components_1a9634fabd09a0ca6a2e9da9af466a17ad>`;
		double :target:`esr<doxid-struct_c_p_m_d_c_energy_components_1a6e29855f5e59dddd5855ae29bb9ab496>`;
		double :target:`eeig<doxid-struct_c_p_m_d_c_energy_components_1a7d3099aca2aad329b36f1e98337542fa>`;
		double :target:`eband<doxid-struct_c_p_m_d_c_energy_components_1aa6247b1a823db8e6ffc8805e162b65c8>`;
		double :target:`entropy<doxid-struct_c_p_m_d_c_energy_components_1a7b673a6dbcbbeac2cd8fc8b866da3706>`;
		double :target:`eself<doxid-struct_c_p_m_d_c_energy_components_1aa8653b2e42b6b411ad094905541da8cd>`;
		double :target:`ecnstr<doxid-struct_c_p_m_d_c_energy_components_1aba081817844e4757340fda2a092b8750>`;
		double :target:`amu<doxid-struct_c_p_m_d_c_energy_components_1ac1cd8bce5dbbc5d092ad21ae4f8da107>`;
		double :target:`ebogo<doxid-struct_c_p_m_d_c_energy_components_1ab6ec3c94b7e95daffc51106c0f11e7f6>`;
		double :target:`eext<doxid-struct_c_p_m_d_c_energy_components_1aa13195c90d973a2219d081a116d5dbed>`;
		double :target:`etddft<doxid-struct_c_p_m_d_c_energy_components_1a6ea61fdf5be3f31122b6cdd1c9c10b1a>`;
		double :target:`ehsic<doxid-struct_c_p_m_d_c_energy_components_1a01d7befc1b826aa857d14ec3f409284e>`;
		double :target:`erestr<doxid-struct_c_p_m_d_c_energy_components_1a2a3b5224cf5e67ca3d3db68a4bd6f3de>`;
		double :target:`eefield<doxid-struct_c_p_m_d_c_energy_components_1aa38bc6fd36ae04d4bd3fb5dadb9e6b49>`;
	};
.. _details-struct_c_p_m_d_c_energy_components:

Detailed Documentation
~~~~~~~~~~~~~~~~~~~~~~

In-process snapshot of OpenCPMD ``ener_com`` scalars (Hartree a.u.).

Filled after a successful embed SCF (``wfopts``) or reference PEF evaluation. Hosts read this via ``:ref:`cpmdc_last_energy_components() <doxid-cpmdc_8h_1a0d2884a4f1288d9f5e76ec4abb505246>``` without parsing CLI ENERGY files or opening a network socket. Field names mirror ``ener_com_t`` in OpenCPMD ``ener.mod.F90``. Zero fields are valid (not set for that run).

Fields
------

.. index:: pair: variable; valid
.. _doxid-struct_c_p_m_d_c_energy_components_1ac6f6f27d3686bcdb5ffcbe470c33f089:

.. ref-code-block:: cpp
	:class: doxyrest-title-code-block

	int valid

Non-zero when the snapshot was written by a successful evaluation.

