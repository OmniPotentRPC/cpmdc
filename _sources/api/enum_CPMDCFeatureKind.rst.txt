.. index:: pair: enum; CPMDCFeatureKind
.. _doxid-cpmdc__features_8h_1af730055345a93286960ba9ef1a3f05b2:

enum CPMDCFeatureKind
=====================

Overview
~~~~~~~~

Feature namespace represented by a ``:ref:`CPMDCFeatureEntry <doxid-struct_c_p_m_d_c_feature_entry>```. :ref:`More...<details-cpmdc__features_8h_1af730055345a93286960ba9ef1a3f05b2>`

.. ref-code-block:: cpp
	:class: doxyrest-overview-code-block

	#include <cpmdc_features.h>

	enum CPMDCFeatureKind {
	    :ref:`CPMDC_FEATURE_SECTION<doxid-cpmdc__features_8h_1af730055345a93286960ba9ef1a3f05b2a10b8c0f0846fd3083eaa32fd2ac815b9>` = 1,
	    :ref:`CPMDC_FEATURE_PARAMS<doxid-cpmdc__features_8h_1af730055345a93286960ba9ef1a3f05b2a51409d149d8bc82536d47a433093fd88>`  = 2,
	    :ref:`CPMDC_FEATURE_ABI<doxid-cpmdc__features_8h_1af730055345a93286960ba9ef1a3f05b2a64f399c7f6179d9eb652b50884d43bd8>`     = 3,
	    :ref:`CPMDC_FEATURE_KEYWORD<doxid-cpmdc__features_8h_1af730055345a93286960ba9ef1a3f05b2ab4937f616b6a49ac5102fdbf11372425>` = 4,
	};

.. _details-cpmdc__features_8h_1af730055345a93286960ba9ef1a3f05b2:

Detailed Documentation
~~~~~~~~~~~~~~~~~~~~~~

Feature namespace represented by a ``:ref:`CPMDCFeatureEntry <doxid-struct_c_p_m_d_c_feature_entry>```.

Enum Values
-----------

.. index:: pair: enumvalue; CPMDC_FEATURE_SECTION
.. _doxid-cpmdc__features_8h_1af730055345a93286960ba9ef1a3f05b2a10b8c0f0846fd3083eaa32fd2ac815b9:

.. ref-code-block:: cpp
	:class: doxyrest-title-code-block

	CPMDC_FEATURE_SECTION

CPMD input section support, such as ``catalog.section.PIMD``.

.. index:: pair: enumvalue; CPMDC_FEATURE_PARAMS
.. _doxid-cpmdc__features_8h_1af730055345a93286960ba9ef1a3f05b2a51409d149d8bc82536d47a433093fd88:

.. ref-code-block:: cpp
	:class: doxyrest-title-code-block

	CPMDC_FEATURE_PARAMS

Cap'n Proto ``CPMDParams`` field support.

.. index:: pair: enumvalue; CPMDC_FEATURE_ABI
.. _doxid-cpmdc__features_8h_1af730055345a93286960ba9ef1a3f05b2a64f399c7f6179d9eb652b50884d43bd8:

.. ref-code-block:: cpp
	:class: doxyrest-title-code-block

	CPMDC_FEATURE_ABI

Exported C ABI symbol support.

.. index:: pair: enumvalue; CPMDC_FEATURE_KEYWORD
.. _doxid-cpmdc__features_8h_1af730055345a93286960ba9ef1a3f05b2ab4937f616b6a49ac5102fdbf11372425:

.. ref-code-block:: cpp
	:class: doxyrest-title-code-block

	CPMDC_FEATURE_KEYWORD

CPMD keyword support inside a catalog section.

