.. index:: pair: struct; CPMDCFeatureEntry
.. _doxid-struct_c_p_m_d_c_feature_entry:

struct CPMDCFeatureEntry
========================

.. toctree::
	:hidden:

Overview
~~~~~~~~

One discoverable ABI, parameter, section, or keyword capability. :ref:`More...<details-struct_c_p_m_d_c_feature_entry>`


.. ref-code-block:: cpp
	:class: doxyrest-overview-code-block

	#include <cpmdc_features.h>
	
	struct CPMDCFeatureEntry {
		// fields
	
		const char* :ref:`feature_id<doxid-struct_c_p_m_d_c_feature_entry_1aeea621b6e77fa5bb0ea64c79bf5788f1>`;
		:ref:`CPMDCFeatureKind<doxid-cpmdc__features_8h_1af730055345a93286960ba9ef1a3f05b2>` :ref:`kind<doxid-struct_c_p_m_d_c_feature_entry_1a25c110dc3efcf9c60c8e9cea5552b875>`;
		int :ref:`stub_applicable<doxid-struct_c_p_m_d_c_feature_entry_1aa6bb7b57531c1d4a83a7869d1f24b039>`;
		int :ref:`embed_applicable<doxid-struct_c_p_m_d_c_feature_entry_1aace20749543a7b4bbb6f5439c183bf41>`;
	};
.. _details-struct_c_p_m_d_c_feature_entry:

Detailed Documentation
~~~~~~~~~~~~~~~~~~~~~~

One discoverable ABI, parameter, section, or keyword capability.

Fields
------

.. index:: pair: variable; feature_id
.. _doxid-struct_c_p_m_d_c_feature_entry_1aeea621b6e77fa5bb0ea64c79bf5788f1:

.. ref-code-block:: cpp
	:class: doxyrest-title-code-block

	const char* feature_id

Stable feature ID string.

.. index:: pair: variable; kind
.. _doxid-struct_c_p_m_d_c_feature_entry_1a25c110dc3efcf9c60c8e9cea5552b875:

.. ref-code-block:: cpp
	:class: doxyrest-title-code-block

	:ref:`CPMDCFeatureKind<doxid-cpmdc__features_8h_1af730055345a93286960ba9ef1a3f05b2>` kind

Feature namespace.

.. index:: pair: variable; stub_applicable
.. _doxid-struct_c_p_m_d_c_feature_entry_1aa6bb7b57531c1d4a83a7869d1f24b039:

.. ref-code-block:: cpp
	:class: doxyrest-title-code-block

	int stub_applicable

Non-zero when the standalone stub build exposes the feature.

.. index:: pair: variable; embed_applicable
.. _doxid-struct_c_p_m_d_c_feature_entry_1aace20749543a7b4bbb6f5439c183bf41:

.. ref-code-block:: cpp
	:class: doxyrest-title-code-block

	int embed_applicable

Non-zero when an embedded OpenCPMD build exposes the feature.

