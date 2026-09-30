.. index:: pair: struct; CPMDCPropertySnapshot
.. _doxid-struct_c_p_m_d_c_property_snapshot:

struct CPMDCPropertySnapshot
============================

.. toctree::
	:hidden:

Overview
~~~~~~~~

PROP-style property snapshot after a successful evaluation. :ref:`More...<details-struct_c_p_m_d_c_property_snapshot>`


.. ref-code-block:: cpp
	:class: doxyrest-overview-code-block

	#include <cpmdc.h>
	
	struct CPMDCPropertySnapshot {
		// fields
	
		int :target:`valid<doxid-struct_c_p_m_d_c_property_snapshot_1ab26d1b54dbe94a3381bd01ef061f2c87>`;
		size_t :target:`hessian_count<doxid-struct_c_p_m_d_c_property_snapshot_1a219308dd9232cb5340fc99466db536a5>`;
		double :target:`hessian<doxid-struct_c_p_m_d_c_property_snapshot_1af26b4dd2e4fbde6ebf3605d6eb778f2e>`[4096];
		size_t :target:`dipole_count<doxid-struct_c_p_m_d_c_property_snapshot_1af22f800017ca74f8ac3d58766c0c7d7e>`;
		double :target:`dipole<doxid-struct_c_p_m_d_c_property_snapshot_1a29b033740147bb25af1b397e20bc7aa9>`[3];
		size_t :target:`polarizability_count<doxid-struct_c_p_m_d_c_property_snapshot_1aa375aea8f4d2476b4d3e8857e1c1e128>`;
		double :target:`polarizability<doxid-struct_c_p_m_d_c_property_snapshot_1a0165ad7c4fe23bd81b91210082cf75cd>`[9];
	};
.. _details-struct_c_p_m_d_c_property_snapshot:

Detailed Documentation
~~~~~~~~~~~~~~~~~~~~~~

PROP-style property snapshot after a successful evaluation.

* dipole[3]: OpenCPMD ``ddippdipole`` (a.u.) when linked; PEF zeros

* polarizability[9]: filled when PROP/aoresponse available; else zeros with count 9

* hessian[]: nuclear gradient dE/dR packed as [natoms\*3] (full Hessian needs dedicated PROP/Hessian run; gradient is always available after force eval)

