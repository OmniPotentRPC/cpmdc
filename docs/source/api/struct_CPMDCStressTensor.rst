.. index:: pair: struct; CPMDCStressTensor
.. _doxid-struct_c_p_m_d_c_stress_tensor:

struct CPMDCStressTensor
========================

.. toctree::
	:hidden:

Overview
~~~~~~~~

Cartesian stress tensor after a successful PEF evaluation. :ref:`More...<details-struct_c_p_m_d_c_stress_tensor>`


.. ref-code-block:: cpp
	:class: doxyrest-overview-code-block

	#include <cpmdc.h>
	
	struct CPMDCStressTensor {
		// fields
	
		int :target:`valid<doxid-struct_c_p_m_d_c_stress_tensor_1a040d79094dfec598b51035b87f1ceb1d>`;
		double :target:`values<doxid-struct_c_p_m_d_c_stress_tensor_1a7061f78949c193ef5285436a368780b3>`[9];
	};
.. _details-struct_c_p_m_d_c_stress_tensor:

Detailed Documentation
~~~~~~~~~~~~~~~~~~~~~~

Cartesian stress tensor after a successful PEF evaluation.

Layout is row-major [xx, xy, xz, yx, yy, yz, zx, zy, zz] in Hartree/Bohr^3 (OpenCPMD ``paiu/omega`` after ``totstr`` when ``cntltpres``). Returns 0 when ``out->valid`` is set; -1 when stress was not computed.

