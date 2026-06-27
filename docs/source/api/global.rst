.. _global:
.. index:: pair: namespace; global

Global Namespace
================

.. toctree::
	:hidden:

	enum_CPMDCFeatureKind.rst
	struct_CPMDCFeatureEntry.rst
	struct_CPMDCResult.rst

Overview
~~~~~~~~




.. ref-code-block:: cpp
	:class: doxyrest-overview-code-block

	
	// typedefs

	typedef struct CPMDCResult :ref:`CPMDCResult<doxid-cpmdc_8h_1a624f3d7d040e783aef1536f49e2a43ca>`;
	typedef struct :ref:`CPMDCSession<doxid-cpmdc_8h_1a49ff1835247d81d2ea1d2bea854ba239>` :ref:`CPMDCSession<doxid-cpmdc_8h_1a49ff1835247d81d2ea1d2bea854ba239>`;
	typedef enum :ref:`CPMDCFeatureKind<doxid-cpmdc__features_8h_1af730055345a93286960ba9ef1a3f05b2>` :ref:`CPMDCFeatureKind<doxid-cpmdc__features_8h_1ab8efa8b96095d9bcb86951f47cdf3ddd>`;
	typedef struct CPMDCFeatureEntry :ref:`CPMDCFeatureEntry<doxid-cpmdc__features_8h_1a78e933e2bdea5d5a752302a384a3137b>`;

	// enums

	enum :ref:`CPMDCFeatureKind<doxid-cpmdc__features_8h_1af730055345a93286960ba9ef1a3f05b2>`;

	// structs

	struct :ref:`CPMDCFeatureEntry<doxid-struct_c_p_m_d_c_feature_entry>`;
	struct :ref:`CPMDCResult<doxid-struct_c_p_m_d_c_result>`;

	// global functions

	int :ref:`cpmdc_set_params<doxid-cpmdc_8h_1a22da50d15419dadce4925c88100a1c23>`(const void* params_capnp, size_t params_capnp_size_bytes);
	:ref:`CPMDCResult<doxid-struct_c_p_m_d_c_result>` :ref:`cpmdc_energy_gradient<doxid-cpmdc_8h_1a525b75c1e47dd8595be42d2175f9f64f>`(int n_atoms, const double* positions_ang, const int* atomic_numbers, const void* params_capnp, size_t params_capnp_size_bytes, double* grad_h_bohr);
	:ref:`CPMDCResult<doxid-struct_c_p_m_d_c_result>` :ref:`cpmdc_energy<doxid-cpmdc_8h_1a6e17756faffc4d181352a5953991f584>`(int n_atoms, const double* positions_ang, const int* atomic_numbers, const void* params_capnp, size_t params_capnp_size_bytes);
	:ref:`CPMDCResult<doxid-struct_c_p_m_d_c_result>` :ref:`cpmdc_energy_forces<doxid-cpmdc_8h_1ad653481360bf399295aed47e14cef118>`(int n_atoms, const double* positions_ang, const int* atomic_numbers, const void* params_capnp, size_t params_capnp_size_bytes, double* forces_h_bohr);
	:ref:`CPMDCSession<doxid-cpmdc_8h_1a49ff1835247d81d2ea1d2bea854ba239>`* :ref:`cpmdc_session_create<doxid-cpmdc_8h_1a08b6b5b994b9754e6c853d3dd46c7745>`(const void* params_capnp, size_t params_capnp_size_bytes);
	int :ref:`cpmdc_session_set_params<doxid-cpmdc_8h_1a0ffb70752daab4a244e5ffd58eaa6db5>`(:ref:`CPMDCSession<doxid-cpmdc_8h_1a49ff1835247d81d2ea1d2bea854ba239>`* session, const void* params_capnp, size_t params_capnp_size_bytes);
	void :ref:`cpmdc_session_destroy<doxid-cpmdc_8h_1a6d07b0b1700b5352f1ce33f396a9f578>`(:ref:`CPMDCSession<doxid-cpmdc_8h_1a49ff1835247d81d2ea1d2bea854ba239>`* session);
	:ref:`CPMDCResult<doxid-struct_c_p_m_d_c_result>` :ref:`cpmdc_session_energy_gradient<doxid-cpmdc_8h_1a6d79318bc9f2eae253d4906376123570>`(:ref:`CPMDCSession<doxid-cpmdc_8h_1a49ff1835247d81d2ea1d2bea854ba239>`* session, int n_atoms, const double* positions_ang, const int* atomic_numbers, double* grad_h_bohr);
	:ref:`CPMDCResult<doxid-struct_c_p_m_d_c_result>` :ref:`cpmdc_session_energy<doxid-cpmdc_8h_1a3c36a0c86b25223e47eb139ec2d89a33>`(:ref:`CPMDCSession<doxid-cpmdc_8h_1a49ff1835247d81d2ea1d2bea854ba239>`* session, int n_atoms, const double* positions_ang, const int* atomic_numbers);
	:ref:`CPMDCResult<doxid-struct_c_p_m_d_c_result>` :ref:`cpmdc_session_energy_forces<doxid-cpmdc_8h_1ab94170b806559374915fb2510fd46d6b>`(:ref:`CPMDCSession<doxid-cpmdc_8h_1a49ff1835247d81d2ea1d2bea854ba239>`* session, int n_atoms, const double* positions_ang, const int* atomic_numbers, double* forces_h_bohr);
	:ref:`CPMDCResult<doxid-struct_c_p_m_d_c_result>` :ref:`cpmdc_session_calculate_forces<doxid-cpmdc_8h_1a9f5e4a955e689dcf3e06b52e68b6adb1>`(:ref:`CPMDCSession<doxid-cpmdc_8h_1a49ff1835247d81d2ea1d2bea854ba239>`* session, const void* force_input_capnp, size_t force_input_capnp_size_bytes, double* forces_h_bohr, size_t forces_len);
	:ref:`CPMDCResult<doxid-struct_c_p_m_d_c_result>` :ref:`cpmdc_session_calculate_result<doxid-cpmdc_8h_1ac63ef923626641269d88b5319f5eed8a>`(:ref:`CPMDCSession<doxid-cpmdc_8h_1a49ff1835247d81d2ea1d2bea854ba239>`* session, const void* force_input_capnp, size_t force_input_capnp_size_bytes, void* potential_result_capnp, size_t potential_result_capnp_capacity_bytes, size_t* potential_result_capnp_size_bytes);
	:ref:`CPMDCResult<doxid-struct_c_p_m_d_c_result>` :ref:`cpmdc_calculate_result<doxid-cpmdc_8h_1a9b206e2d2173eafb39c00977101da25b>`(const void* params_capnp, size_t params_capnp_size_bytes, const void* force_input_capnp, size_t force_input_capnp_size_bytes, void* potential_result_capnp, size_t potential_result_capnp_capacity_bytes, size_t* potential_result_capnp_size_bytes);
	size_t :ref:`cpmdc_potential_result_size_for_force_input<doxid-cpmdc_8h_1a50c58b775c95334c2842f3d07c247408>`(const void* force_input_capnp, size_t force_input_capnp_size_bytes);
	const char* :ref:`cpmdc_version<doxid-cpmdc_8h_1a4d5a0eade3deceabe9f7e48e7a45d430>`(void);
	int :ref:`cpmdc_available<doxid-cpmdc_8h_1a25f2bc651370d207f2049574e7073e0d>`(void);
	void :ref:`cpmdc_finalize<doxid-cpmdc_8h_1afae0c7a27c94036223f7a0cf60488f56>`(void);
	size_t :ref:`cpmdc_feature_count<doxid-cpmdc__features_8h_1a4cf7e2f669806606c1d1293533b3d38e>`(void);
	const :ref:`CPMDCFeatureEntry<doxid-struct_c_p_m_d_c_feature_entry>`* :ref:`cpmdc_feature_table<doxid-cpmdc__features_8h_1ae40eb0e96927dc5f3feacc20c3869315>`(void);
	const :ref:`CPMDCFeatureEntry<doxid-struct_c_p_m_d_c_feature_entry>`* :ref:`cpmdc_feature_find<doxid-cpmdc__features_8h_1a4d2414cdfa87d867cf1ab1e708dfb4ca>`(const char* feature_id);

.. _details-global:

Detailed Documentation
~~~~~~~~~~~~~~~~~~~~~~



Typedefs
--------

.. index:: pair: typedef; CPMDCResult
.. _doxid-cpmdc_8h_1a624f3d7d040e783aef1536f49e2a43ca:

.. ref-code-block:: cpp
	:class: doxyrest-title-code-block

	typedef struct CPMDCResult CPMDCResult

Result returned by energy / gradient / forces entry points.

.. index:: pair: typedef; CPMDCSession
.. _doxid-cpmdc_8h_1a49ff1835247d81d2ea1d2bea854ba239:

.. ref-code-block:: cpp
	:class: doxyrest-title-code-block

	typedef struct :ref:`CPMDCSession<doxid-cpmdc_8h_1a49ff1835247d81d2ea1d2bea854ba239>` CPMDCSession

Opaque handle for repeated evaluations with one Cap'n Proto parameter set.

.. index:: pair: typedef; CPMDCFeatureKind
.. _doxid-cpmdc__features_8h_1ab8efa8b96095d9bcb86951f47cdf3ddd:

.. ref-code-block:: cpp
	:class: doxyrest-title-code-block

	typedef enum :ref:`CPMDCFeatureKind<doxid-cpmdc__features_8h_1af730055345a93286960ba9ef1a3f05b2>` CPMDCFeatureKind

Feature namespace represented by a ``:ref:`CPMDCFeatureEntry <doxid-struct_c_p_m_d_c_feature_entry>```.

.. index:: pair: typedef; CPMDCFeatureEntry
.. _doxid-cpmdc__features_8h_1a78e933e2bdea5d5a752302a384a3137b:

.. ref-code-block:: cpp
	:class: doxyrest-title-code-block

	typedef struct CPMDCFeatureEntry CPMDCFeatureEntry

One discoverable ABI, parameter, section, or keyword capability.

Global Functions
----------------

.. index:: pair: function; cpmdc_set_params
.. _doxid-cpmdc_8h_1a22da50d15419dadce4925c88100a1c23:

.. ref-code-block:: cpp
	:class: doxyrest-title-code-block

	int cpmdc_set_params(const void* params_capnp, size_t params_capnp_size_bytes)

Apply CPMD method parameters from a Cap'n Proto message.



.. rubric:: Parameters:

.. list-table::
	:widths: 20 80

	*
		- params_capnp

		- Pointer to an unpacked flat ``CPMDParams`` message.

	*
		- params_capnp_size_bytes

		- Size of ``params_capnp`` in bytes.



.. rubric:: Returns:

0 on success, -1 on parse or configuration failure.

.. index:: pair: function; cpmdc_energy_gradient
.. _doxid-cpmdc_8h_1a525b75c1e47dd8595be42d2175f9f64f:

.. ref-code-block:: cpp
	:class: doxyrest-title-code-block

	:ref:`CPMDCResult<doxid-struct_c_p_m_d_c_result>` cpmdc_energy_gradient(int n_atoms, const double* positions_ang, const int* atomic_numbers, const void* params_capnp, size_t params_capnp_size_bytes, double* grad_h_bohr)

Compute energy and nuclear gradient for an atomic configuration.

Positions are Angstrom; gradient is Hartree/Bohr (CPMD ionic forces are negated into a nuclear gradient for API symmetry with nwchemc).

.. index:: pair: function; cpmdc_energy
.. _doxid-cpmdc_8h_1a6e17756faffc4d181352a5953991f584:

.. ref-code-block:: cpp
	:class: doxyrest-title-code-block

	:ref:`CPMDCResult<doxid-struct_c_p_m_d_c_result>` cpmdc_energy(int n_atoms, const double* positions_ang, const int* atomic_numbers, const void* params_capnp, size_t params_capnp_size_bytes)

Compute total energy only (no gradient allocation).

.. index:: pair: function; cpmdc_energy_forces
.. _doxid-cpmdc_8h_1ad653481360bf399295aed47e14cef118:

.. ref-code-block:: cpp
	:class: doxyrest-title-code-block

	:ref:`CPMDCResult<doxid-struct_c_p_m_d_c_result>` cpmdc_energy_forces(int n_atoms, const double* positions_ang, const int* atomic_numbers, const void* params_capnp, size_t params_capnp_size_bytes, double* forces_h_bohr)

Compute energy and nuclear forces (negative gradient, Hartree/Bohr).

.. index:: pair: function; cpmdc_session_create
.. _doxid-cpmdc_8h_1a08b6b5b994b9754e6c853d3dd46c7745:

.. ref-code-block:: cpp
	:class: doxyrest-title-code-block

	:ref:`CPMDCSession<doxid-cpmdc_8h_1a49ff1835247d81d2ea1d2bea854ba239>`* cpmdc_session_create(const void* params_capnp, size_t params_capnp_size_bytes)

Create a persistent evaluation session from a Cap'n Proto message.

The session owns a copy of the serialized message so callers may release the input buffer after this call returns.

.. index:: pair: function; cpmdc_session_set_params
.. _doxid-cpmdc_8h_1a0ffb70752daab4a244e5ffd58eaa6db5:

.. ref-code-block:: cpp
	:class: doxyrest-title-code-block

	int cpmdc_session_set_params(:ref:`CPMDCSession<doxid-cpmdc_8h_1a49ff1835247d81d2ea1d2bea854ba239>`* session, const void* params_capnp, size_t params_capnp_size_bytes)

Replace Cap'n Proto parameters before the session accepts topology.

.. index:: pair: function; cpmdc_session_destroy
.. _doxid-cpmdc_8h_1a6d07b0b1700b5352f1ce33f396a9f578:

.. ref-code-block:: cpp
	:class: doxyrest-title-code-block

	void cpmdc_session_destroy(:ref:`CPMDCSession<doxid-cpmdc_8h_1a49ff1835247d81d2ea1d2bea854ba239>`* session)

Release a persistent evaluation session.

.. index:: pair: function; cpmdc_session_energy_gradient
.. _doxid-cpmdc_8h_1a6d79318bc9f2eae253d4906376123570:

.. ref-code-block:: cpp
	:class: doxyrest-title-code-block

	:ref:`CPMDCResult<doxid-struct_c_p_m_d_c_result>` cpmdc_session_energy_gradient(:ref:`CPMDCSession<doxid-cpmdc_8h_1a49ff1835247d81d2ea1d2bea854ba239>`* session, int n_atoms, const double* positions_ang, const int* atomic_numbers, double* grad_h_bohr)

Compute energy and nuclear gradient with session-owned parameters.

Positions are Angstrom. The gradient buffer must have ``n_atoms * 3`` doubles and is filled in Hartree/Bohr.

.. index:: pair: function; cpmdc_session_energy
.. _doxid-cpmdc_8h_1a3c36a0c86b25223e47eb139ec2d89a33:

.. ref-code-block:: cpp
	:class: doxyrest-title-code-block

	:ref:`CPMDCResult<doxid-struct_c_p_m_d_c_result>` cpmdc_session_energy(:ref:`CPMDCSession<doxid-cpmdc_8h_1a49ff1835247d81d2ea1d2bea854ba239>`* session, int n_atoms, const double* positions_ang, const int* atomic_numbers)

Compute total energy with session-owned parameters.

Positions are Angstrom and the returned energy is Hartree.

.. index:: pair: function; cpmdc_session_energy_forces
.. _doxid-cpmdc_8h_1ab94170b806559374915fb2510fd46d6b:

.. ref-code-block:: cpp
	:class: doxyrest-title-code-block

	:ref:`CPMDCResult<doxid-struct_c_p_m_d_c_result>` cpmdc_session_energy_forces(:ref:`CPMDCSession<doxid-cpmdc_8h_1a49ff1835247d81d2ea1d2bea854ba239>`* session, int n_atoms, const double* positions_ang, const int* atomic_numbers, double* forces_h_bohr)

Compute energy and nuclear forces with session-owned parameters.

Positions are Angstrom. The forces buffer must have ``n_atoms * 3`` doubles and is filled in Hartree/Bohr.

.. index:: pair: function; cpmdc_session_calculate_forces
.. _doxid-cpmdc_8h_1a9f5e4a955e689dcf3e06b52e68b6adb1:

.. ref-code-block:: cpp
	:class: doxyrest-title-code-block

	:ref:`CPMDCResult<doxid-struct_c_p_m_d_c_result>` cpmdc_session_calculate_forces(:ref:`CPMDCSession<doxid-cpmdc_8h_1a49ff1835247d81d2ea1d2bea854ba239>`* session, const void* force_input_capnp, size_t force_input_capnp_size_bytes, double* forces_h_bohr, size_t forces_len)

Compute energy and forces for one Cap'n Proto ``ForceInput`` step.

Session keeps persistent ``CPMDParams``; each call supplies geometry. Returned energy/forces use CPMD native units: Hartree and Hartree/Bohr.

.. index:: pair: function; cpmdc_session_calculate_result
.. _doxid-cpmdc_8h_1ac63ef923626641269d88b5319f5eed8a:

.. ref-code-block:: cpp
	:class: doxyrest-title-code-block

	:ref:`CPMDCResult<doxid-struct_c_p_m_d_c_result>` cpmdc_session_calculate_result(:ref:`CPMDCSession<doxid-cpmdc_8h_1a49ff1835247d81d2ea1d2bea854ba239>`* session, const void* force_input_capnp, size_t force_input_capnp_size_bytes, void* potential_result_capnp, size_t potential_result_capnp_capacity_bytes, size_t* potential_result_capnp_size_bytes)

Compute forces for one ``ForceInput`` step and write ``PotentialResult``.

Direct-call socket entry point: method state in the session, geometry in ``ForceInput``, output energy/forces converted to ``ForceInput.energyUnit`` and ``energyUnit / lengthUnit``.

When ``potential_result_capnp_capacity_bytes`` is too small, returns ``ok == 0``, writes the required byte count to ``potential_result_capnp_size_bytes``, and does not evaluate CPMD.

.. index:: pair: function; cpmdc_calculate_result
.. _doxid-cpmdc_8h_1a9b206e2d2173eafb39c00977101da25b:

.. ref-code-block:: cpp
	:class: doxyrest-title-code-block

	:ref:`CPMDCResult<doxid-struct_c_p_m_d_c_result>` cpmdc_calculate_result(const void* params_capnp, size_t params_capnp_size_bytes, const void* force_input_capnp, size_t force_input_capnp_size_bytes, void* potential_result_capnp, size_t potential_result_capnp_capacity_bytes, size_t* potential_result_capnp_size_bytes)

One-shot Cap'n Proto entry point (params + ForceInput -> PotentialResult).

Multi-step callers should create one session and call ``:ref:`cpmdc_session_calculate_result() <doxid-cpmdc_8h_1ac63ef923626641269d88b5319f5eed8a>``` per step.

.. index:: pair: function; cpmdc_potential_result_size_for_force_input
.. _doxid-cpmdc_8h_1a50c58b775c95334c2842f3d07c247408:

.. ref-code-block:: cpp
	:class: doxyrest-title-code-block

	size_t cpmdc_potential_result_size_for_force_input(const void* force_input_capnp, size_t force_input_capnp_size_bytes)

Byte count needed for a ``PotentialResult`` for the given ``ForceInput``.

Parses geometry only; does not initialize or evaluate CPMD. Returns 0 when the message is invalid or too large for the C ABI.

.. index:: pair: function; cpmdc_version
.. _doxid-cpmdc_8h_1a4d5a0eade3deceabe9f7e48e7a45d430:

.. ref-code-block:: cpp
	:class: doxyrest-title-code-block

	const char* cpmdc_version(void)

Compiled library version string.

.. index:: pair: function; cpmdc_available
.. _doxid-cpmdc_8h_1a25f2bc651370d207f2049574e7073e0d:

.. ref-code-block:: cpp
	:class: doxyrest-title-code-block

	int cpmdc_available(void)

1 when the embedded OpenCPMD runtime is available.

.. index:: pair: function; cpmdc_finalize
.. _doxid-cpmdc_8h_1afae0c7a27c94036223f7a0cf60488f56:

.. ref-code-block:: cpp
	:class: doxyrest-title-code-block

	void cpmdc_finalize(void)

Finalize an owned embedded CPMD runtime.

.. index:: pair: function; cpmdc_feature_count
.. _doxid-cpmdc__features_8h_1a4cf7e2f669806606c1d1293533b3d38e:

.. ref-code-block:: cpp
	:class: doxyrest-title-code-block

	size_t cpmdc_feature_count(void)

Number of entries returned by ``:ref:`cpmdc_feature_table() <doxid-cpmdc__features_8h_1ae40eb0e96927dc5f3feacc20c3869315>```.

.. index:: pair: function; cpmdc_feature_table
.. _doxid-cpmdc__features_8h_1ae40eb0e96927dc5f3feacc20c3869315:

.. ref-code-block:: cpp
	:class: doxyrest-title-code-block

	const :ref:`CPMDCFeatureEntry<doxid-struct_c_p_m_d_c_feature_entry>`* cpmdc_feature_table(void)

Return the contiguous feature table.

The returned pointer has ``:ref:`cpmdc_feature_count() <doxid-cpmdc__features_8h_1a4cf7e2f669806606c1d1293533b3d38e>``` entries and remains owned by the library.

.. index:: pair: function; cpmdc_feature_find
.. _doxid-cpmdc__features_8h_1a4d2414cdfa87d867cf1ab1e708dfb4ca:

.. ref-code-block:: cpp
	:class: doxyrest-title-code-block

	const :ref:`CPMDCFeatureEntry<doxid-struct_c_p_m_d_c_feature_entry>`* cpmdc_feature_find(const char* feature_id)

Find one feature by stable ID.



.. rubric:: Returns:

Pointer to the matching table entry, or ``NULL`` when the feature ID is not present.

