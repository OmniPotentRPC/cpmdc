.. index:: pair: struct; CPMDCChargeIntegrals
.. _doxid-struct_c_p_m_d_c_charge_integrals:

struct CPMDCChargeIntegrals
===========================

.. toctree::
	:hidden:

OpenCPMD ``chrg_t`` density integrals (post-SCF module state).


.. ref-code-block:: cpp
	:class: doxyrest-overview-code-block

	#include <cpmdc.h>
	
	struct CPMDCChargeIntegrals {
		// fields
	
		int :target:`valid<doxid-struct_c_p_m_d_c_charge_integrals_1a0d192c660b5e55121b4a0e001da82bc9>`;
		double :target:`csumg<doxid-struct_c_p_m_d_c_charge_integrals_1ac1647c1b8130835da878b36698acc333>`;
		double :target:`csumr<doxid-struct_c_p_m_d_c_charge_integrals_1afae80503e8a13771e129b991a900690a>`;
		double :target:`csums<doxid-struct_c_p_m_d_c_charge_integrals_1acb1e9588b5b90337f3fd1ac786bb6800>`;
		double :target:`csumsabs<doxid-struct_c_p_m_d_c_charge_integrals_1a7b2d031dcf7ede95440c7ad96cb80b17>`;
	};
