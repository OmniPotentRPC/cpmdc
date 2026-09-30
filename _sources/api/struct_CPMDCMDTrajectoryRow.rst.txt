.. index:: pair: struct; CPMDCMDTrajectoryRow
.. _doxid-struct_c_p_m_d_c_m_d_trajectory_row:

struct CPMDCMDTrajectoryRow
===========================

.. toctree::
	:hidden:

Overview
~~~~~~~~

One ENERGY-file-equivalent trajectory row (Hartree). :ref:`More...<details-struct_c_p_m_d_c_m_d_trajectory_row>`


.. ref-code-block:: cpp
	:class: doxyrest-overview-code-block

	#include <cpmdc.h>
	
	struct CPMDCMDTrajectoryRow {
		// fields
	
		int :target:`valid<doxid-struct_c_p_m_d_c_m_d_trajectory_row_1ae8e42435497115257eb51b14aa4b7702>`;
		size_t :target:`count<doxid-struct_c_p_m_d_c_m_d_trajectory_row_1ad2d5cdd0881a4ab7988b448339eb6f8f>`;
		double :target:`values<doxid-struct_c_p_m_d_c_m_d_trajectory_row_1aba443c14d800f6c189cabec43a447708>`[32];
	};
.. _details-struct_c_p_m_d_c_m_d_trajectory_row:

Detailed Documentation
~~~~~~~~~~~~~~~~~~~~~~

One ENERGY-file-equivalent trajectory row (Hartree).

Layout (count >= 12 after a successful eval): [0] etot [1] ekin [2] epseu [3] enl [4] eht [5] exc [6] ehep [7] ehee [8] ehii [9] esr [10] eself [11] EKINC (fictitious electronic KE; 0 for BO/SCF-only wfopt, filled in MD)

