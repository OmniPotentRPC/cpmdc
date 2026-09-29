.. index:: pair: struct; CPMDCMultiStateEnergies
.. _doxid-struct_c_p_m_d_c_multi_state_energies:

struct CPMDCMultiStateEnergies
==============================

.. toctree::
	:hidden:

Overview
~~~~~~~~

Flattened CAS22-class multi-state catalog (``ener_c`` + ``ener_d``). :ref:`More...<details-struct_c_p_m_d_c_multi_state_energies>`


.. ref-code-block:: cpp
	:class: doxyrest-overview-code-block

	#include <cpmdc.h>
	
	struct CPMDCMultiStateEnergies {
		// fields
	
		int :target:`valid<doxid-struct_c_p_m_d_c_multi_state_energies_1a0d67163f5e59c14f9202d36570556d50>`;
		size_t :target:`count<doxid-struct_c_p_m_d_c_multi_state_energies_1aec6f9826107ef87accfbf42f6ea42e08>`;
		double :target:`values<doxid-struct_c_p_m_d_c_multi_state_energies_1aaff5c48a1ff24c43903c0b7dbfbbaae5>`[64];
	};
.. _details-struct_c_p_m_d_c_multi_state_energies:

Detailed Documentation
~~~~~~~~~~~~~~~~~~~~~~

Flattened CAS22-class multi-state catalog (``ener_c`` + ``ener_d``).

Layout is backend-defined; ``count`` is the number of doubles copied into ``values`` (caller provides capacity). Returns -1 when no snapshot.

