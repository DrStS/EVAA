#pragma once
#include <mkl.h>
#include "ReadXML.h"
#include "MathLibrary.h"

template <class T>
class Car {
private:
	void construct_11DOF_mass(T* Global_mass, T* Global_momemnt_Inertia, T* Mass_11DOF) {
		temp_linear[0] = Global_mass[0];
		temp_linear[1] = Global_momemnt_Inertia[0];
		temp_linear[2] = Global_momemnt_Inertia[4];
		cblas_dcopy(vec_DIM - 1, Global_mass + 1, 1, temp_linear + 3, 1);
		MathLibrary::allocate_to_diagonal(M_linear, temp_linear, DOF);
	}
	void construct_ALE_vectors(T* Global_vector, T* local_vector) {
		T* start_pointer, *current_ptr;
		start_pointer = Global_vector; // copy x and y and move next
		current_ptr = local_vector;
		for (size_t i = 0; i < vec_DIM; ++i) {
			cblas_dcopy(DIM - 1, start_pointer, 1, current_ptr, 1);
			start_pointer += DIM;
			current_ptr += DIM - 1;
		}
	}
	void populate_results() {
		/*
		Not implemented
		assign the ALE components to x and y direction and 11 dof components to z direction
		*/
	}
public:
	/*
	Using Stefan's convention
	W1 = front left wheel
	W2 = front right wheel
	W3 = rear left wheel
	W4 = rear right wheel
	T1 = front left tyre
	T2 = front right tyre
	T3 = rear left tyre
	T4 = rear right tyre
	*/

	// consider replacing double with floatEVAA
	// alignment and g are defined globally === REPLACE THEM!!!!
	const int alignment = 64;
	// MKL / vector constants:
	const int DIM = 3, vec_DIM = 9, incx = 1; // consider DIM = 4 for efficiency!!! // 10 dimension because of torque of the body
	const size_t num_wheels = 4;
	const int NUM_LEGS = 4, num_tyre = 4;
	T* Position_vec; // [CG, W1, T1, W2, T2, W3, T3, W4, T4] 9 x 3 !!! Consider alignment (3+1),(3+1),... 
	T* Velocity_vec; // [CG, W1, T1, W2, T2, W3, T3, W4, T4] 9 x 3
	T* Mass_vec; // [CG,  W1, T1, W2, T2, W3, T3, W4, T4]     9
	T* angle_CG; // [x, y, z]
	T* w_CG; // [x, y, z]
	T* I_CG; // [Ixx, Ixy, Ixz, Iyx, Iyy, Iyz, Izx, Izy, Izz]
	

	// Initial Conditions of the car
	T* initial_position; // [CG, W1, T1, W2, T2, W3, T3, W4, T4]
	T* initial_velocity_vec; // [CG, W1, T1, W2, T2, W3, T3, W4, T4]
	T* initial_angle; // [x, y, z]
	T* initial_angular_velocity; // [x, y, z]


	/////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////// Members from 11 DOF system //////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////
	int DOF;
	EVAAComputeStiffness* lookupStiffness;
	T k_body_fl;
	T k_tyre_fl;
	T k_body_fr;
	T k_tyre_fr;
	T k_body_rl;
	T k_tyre_rl;
	T k_body_rr;
	T k_tyre_rr;
	
	int i;
	T *quad_angle_init;
	T* M_linear, *temp_linear, *K, *K_trans, *D;
	T *spring_length, *current_spring_length;
	T* u_prev_linear, *u_current_linear;
	T* k_vec, *l_lat, *l_long;
	T* velocity_current_linear;

	////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////// ALE Vectors ///////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////
	T *Position_vec_xy, *Angle_z, *Velocity_vec_xy, *w_z;


	
	Car(const Simulation_Parameters &params, EVAAComputeStiffness* interpolator) {
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////// Generte Lookup Table /////////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		lookupStiffness = interpolator;
		DOF = params.DOF;
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////// System Mono ///////////////////////////////////////////////////////
		///////////////////////////////// Memory Allocation and matrix formulation ///////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////// Params for Global coordinate ////////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		Position_vec = (T*)mkl_malloc(DIM*vec_DIM * sizeof(T), alignment);
		Velocity_vec = (T*)mkl_malloc(DIM*vec_DIM * sizeof(T), alignment);
		Mass_vec = (T*)mkl_malloc(vec_DIM * sizeof(T), alignment);
		angle_CG = (T*)mkl_malloc(DIM * sizeof(T), alignment);
		w_CG = (T*)mkl_malloc(DIM * sizeof(T), alignment);
		I_CG = (T*)mkl_malloc(DIM*DIM * sizeof(T), alignment);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////// Initial Params for Global coordinate ////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		initial_position = (T*)mkl_malloc(DIM*vec_DIM * sizeof(T), alignment);
		initial_velocity_vec = (T*)mkl_malloc(DIM*vec_DIM * sizeof(T), alignment);
		quad_angle_init = (T*)mkl_calloc(4, sizeof(T), alignment);
		initial_angle = (T*)mkl_malloc(DIM * sizeof(T), alignment);
		initial_angular_velocity = (T*)mkl_malloc(DIM * sizeof(T), alignment);
		


		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////// Params for 11 DOF system ////////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		M_linear = (T*)mkl_calloc(DOF * DOF, sizeof(T), alignment);
		temp_linear = (T*)mkl_calloc(DOF , sizeof(T), alignment);
		K = (T*)mkl_calloc(DOF * DOF, sizeof(T), alignment);
		K_trans = (T*)mkl_calloc(DOF * DOF, sizeof(T), alignment);
		D = (T*)mkl_calloc(DOF * DOF, sizeof(T), alignment);
		u_prev_linear = (T*)mkl_malloc(DOF*sizeof(T), alignment);
		u_current_linear = (T*)mkl_malloc(DOF*sizeof(T), alignment);
		velocity_current_linear = (T*)mkl_malloc(DOF * sizeof(T), alignment);
		k_vec = (T*)mkl_malloc(num_wheels * 2 * sizeof(T), alignment);
		l_lat = (T*)mkl_malloc(num_wheels * sizeof(T), alignment);
		l_long = (T*)mkl_malloc(num_wheels * sizeof(T), alignment);
		spring_length = (T*)mkl_malloc(2 * num_wheels * sizeof(T), alignment);
		current_spring_length = (T*)mkl_malloc(2 * num_wheels * sizeof(T), alignment);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////// Extract Data from parser /////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		
		
		l_long[0] = params.l_long[2];
		l_long[1] = params.l_long[3];
		l_long[2] = params.l_long[1];
		l_long[3] = params.l_long[0];
		
		l_lat[0] = params.l_lat[2];
		l_lat[1] = params.l_lat[3];
		l_lat[2] = params.l_lat[1];
		l_lat[3] = params.l_lat[0];


		I_CG[0] = params.I_body[0];
		I_CG[1] = params.I_body[3];
		I_CG[2] = params.I_body[4];
		I_CG[3] = params.I_body[3];
		I_CG[4] = params.I_body[1];
		I_CG[5] = params.I_body[5];
		I_CG[6] = params.I_body[4];
		I_CG[7] = params.I_body[5];
		I_CG[8] = params.I_body[2];
		
		Mass_vec[0] = params.mass_body;
		Mass_vecMass_vec[1] = params.mass_wheel[2];
		Mass_vec[2] = params.mass_tyre[2];
		Mass_vec[3] = params.mass_wheel[3];
		Mass_vec[4] = params.mass_tyre[3];
		Mass_vec[5] = params.mass_wheel[1];
		Mass_vec[6] = params.mass_tyre[1];
		Mass_vec[7] = params.mass_wheel[0];
		Mass_vec[8] = params.mass_tyre[0];
		
		spring_length[0] = params.upper_spring_length[2];
		spring_length[1] = params.lower_spring_length[2];
		spring_length[2] = params.upper_spring_length[3];
		spring_length[3] = params.lower_spring_length[3];
		spring_length[4] = params.lower_spring_length[1];
		spring_length[5] = params.upper_spring_length[1];
		spring_length[6] = params.upper_spring_length[0];
		spring_length[7] = params.lower_spring_length[0];


		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////// Initial Iteration vector ////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		// Initial Angles
		quad_angle_init[0] = params.initial_angle[0];
		quad_angle_init[1] = params.initial_angle[1];
		quad_angle_init[2] = params.initial_angle[2];
		quad_angle_init[3] = params.initial_angle[3];
		MathLibrary::ToEulerAngles<T>(quad_angle_init, initial_angle);
		cblas_dcopy(DIM, initial_angle, 1, angle_CG, 1);
		
		// Spring lengths
		current_spring_length[0] = params.initial_upper_spring_length[2];
		current_spring_length[1] = params.initial_lower_spring_length[2];
		current_spring_length[2] = params.initial_upper_spring_length[3];
		current_spring_length[3] = params.initial_lower_spring_length[3];
		current_spring_length[4] = params.initial_upper_spring_length[1];
		current_spring_length[5] = params.initial_lower_spring_length[1];
		current_spring_length[6] = params.initial_upper_spring_length[0];
		current_spring_length[7] = params.initial_lower_spring_length[0];

		// Filling the position vector with initial condition
		// CG
		cblas_dcopy(DIM, params.initial_pos_body, 1, initial_position, 1); // copy the center of mass position
		

		const T* xml_start;
		T *position_start;
		if (params.initial_leg_flag) {
			// if prescribed initial position (add a check for consistency with spring lengths)
			// W1 = W_fl
			xml_start = params.initial_pos_wheel + 2 * 3;
			position_start = initial_position + 3;
			cblas_dcopy(DIM, xml_start, 1, position_start, 1);
			// W2 = W_fr
			xml_start = params.initial_pos_wheel + 3 * 3;
			position_start += 6; // skip 3 for tyre
			cblas_dcopy(DIM, xml_start, 1, position_start, 1);
			// W3 = W_rl
			xml_start = params.initial_pos_wheel + 1 * 3;
			position_start += 6; // skip 3 for tyre
			cblas_dcopy(DIM, xml_start, 1, position_start, 1);
			// W2 = W_rr
			xml_start = params.initial_pos_wheel + 0 * 3;
			position_start += 6; // skip 3 for tyre
			cblas_dcopy(DIM, xml_start, 1, position_start, 1);

			// T1 = T_fl
			xml_start = params.initial_pos_tyre + 2 * 3;
			position_start = initial_position + 6; // skip 3 for center of mass and 3 for the wheel
			cblas_dcopy(DIM, xml_start, 1, position_start, 1);
			// T2 = T_fr
			xml_start = params.initial_pos_tyre + 3 * 3;
			position_start += 6; // skip 3 for the wheel
			cblas_dcopy(DIM, xml_start, 1, position_start, 1);
			// T3 = T_rl
			xml_start = params.initial_pos_tyre + 1 * 3;
			position_start += 6; // skip 3 for the wheel
			cblas_dcopy(DIM, xml_start, 1, position_start, 1);
			// T4 = T_rr
			xml_start = params.initial_pos_tyre + 0 * 3;
			position_start += 6; // skip 3 for the wheel
			cblas_dcopy(DIM, xml_start, 1, position_start, 1);
		}
		else {
			/*Not implemented
			compute positions of the wheel and tyre
			*/
			
		}
		// copy the initial position to the position vector
		cblas_dcopy(DIM*vec_DIM, initial_position, 1, Position_vec, 1);

		

		// Initial Velocity (Reuse the pointers)
		cblas_dcopy(DIM, params.initial_vel_body, 1, initial_velocity_vec, 1);
		// W1 = W_fl
		xml_start = params.initial_vel_body + 2 * 3;
		position_start = initial_velocity_vec + 3;
		cblas_dcopy(DIM, xml_start, 1, position_start, 1);
		// W2 = W_fr
		xml_start = params.initial_vel_body + 3 * 3;
		position_start += 6; // skip 3 for tyre
		cblas_dcopy(DIM, xml_start, 1, position_start, 1);
		// W3 = W_rl
		xml_start = params.initial_vel_body + 1 * 3;
		position_start += 6; // skip 3 for tyre
		cblas_dcopy(DIM, xml_start, 1, position_start, 1);
		// W2 = W_rr
		xml_start = params.initial_vel_body + 0 * 3;
		position_start += 6; // skip 3 for tyre
		cblas_dcopy(DIM, xml_start, 1, position_start, 1);

		// T1 = T_fl
		xml_start = params.initial_vel_body + 2 * 3;
		position_start = initial_velocity_vec + 6; // skip 3 for center of mass and 3 for the wheel
		cblas_dcopy(DIM, xml_start, 1, position_start, 1);
		// T2 = T_fr
		xml_start = params.initial_vel_body + 3 * 3;
		position_start += 6; // skip 3 for the wheel
		cblas_dcopy(DIM, xml_start, 1, position_start, 1);
		// T3 = T_rl
		xml_start = params.initial_vel_body + 1 * 3;
		position_start += 6; // skip 3 for the wheel
		cblas_dcopy(DIM, xml_start, 1, position_start, 1);
		// T4 = T_rr
		xml_start = params.initial_vel_body + 0 * 3;
		position_start += 6; // skip 3 for the wheel
		cblas_dcopy(DIM, xml_start, 1, position_start, 1);

		// copy the initial position to the position vector
		cblas_dcopy(DIM*vec_DIM, initial_velocity_vec, 1, Velocity_vec, 1);

		// Initial Angular velocity
		cblas_dcopy(DIM, params.initial_ang_vel_body, 1, initial_angular_velocity, 1);
		cblas_dcopy(DIM, initial_angular_velocity, 1, w_CG, 1);
		
		/*
		Global assignments Done 
		11 DOF assignments
		*/
		construct_11DOF_mass(Mass_vec, I_CG, M_linear);
		construct_11DOF_vector(initial_position, initial_angle, u_prev_linear);
		construct_11DOF_vector(initial_velocity_vec, initial_angular_velocity, velocity_current_linear);
		
		/* This stays in the 11 DOF
		cblas_dscal(DOF, -h_, u_n_m_1, 1);
		cblas_daxpy(DOF, 1, u_n, 1, u_n_m_1, 1);
		*/
		if (params.interpolation) {
			lookupStiffness->getStiffness(current_spring_length, k_vec);
		}
		else {
			k_body_fl = params.k_body[2];
			k_tyre_fl = params.k_tyre[2];
			k_body_fr = params.k_body[3];
			k_tyre_fr = params.k_tyre[3];
			k_body_rl = params.k_body[1];
			k_tyre_rl = params.k_tyre[1];
			k_body_rr = params.k_body[0];
			k_tyre_rr = params.k_tyre[0];
			populate_K(k_vec, k_body_fl, k_tyre_fl, k_body_fr, k_tyre_fr, k_body_rl, k_tyre_rl, k_body_rr, k_tyre_rr);
		}
		construct_K(K, k_vec, l_lat, l_long);


		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////// ALE Buffer Initialization /////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		Position_vec_xy = (T*)mkl_malloc((DIM - 1)*vec_DIM, alignment);
		Angle_z = new T;
		Velocity_vec_xy = (T*)mkl_malloc((DIM - 1)*vec_DIM, alignment);
		w_z = new T;
		construct_ALE_vectors(Position_vec, Position_vec_xy);
		construct_ALE_vectors(Velocity_vec, Velocity_vec_xy);
		*Angle_z = angle_CG[2];
		*w_z = w_CG[2];


	}

	void apply_normal_force(T* force, T* u, size_t* index, size_t n) {
		#pragma loop( ivdep )
		for (int i = 0; i < n; ++i) {
			force[index[i]] = force[index[i]] > 0 ? force[index[i]] : 0;
		}
		#pragma loop( ivdep )
		for (int i = 0; i < n; ++i) {
			u[index[i]] = u[index[i]] > 0 ? u[index[i]] : 0;
		}
	}
	void compute_normal_force(T* K, T* u, T* force, size_t* index, size_t dim, size_t n) {
		#pragma loop( ivdep )
		for (int i = 0; i < n; ++i) {
			force[index[i]] = -K[index[i] * dim + index[i]] * u[index[i]];
		}
	}
	void construct_K(T* K, T *k_vect, T *l_lat, T* l_long) {
		int i;
		i = 0;
		T k_body_fl_, k_body_fr_, k_body_rl_, k_body_rr_, k_tyre_fl_, k_tyre_fr_, k_tyre_rl_, k_tyre_rr_;
		k_body_fl_ = k_vect[0];
		k_tyre_fl_ = k_vect[1];
		k_body_fr_ = k_vect[2];
		k_tyre_fr_ = k_vect[3];
		k_body_rl_ = k_vect[4];
		k_tyre_rl_ = k_vect[5];
		k_body_rr_ = k_vect[6];
		k_tyre_rr_ = k_vect[7];
		T l_lat_fl_, l_lat_fr_, l_lat_rl_, l_lat_rr_, l_long_fl_, l_long_fr_, l_long_rl_, l_long_rr_;
		l_lat_fl_ = l_lat[0];
		l_lat_fr_ = l_lat[1];
		l_lat_rl_ = l_lat[2];
		l_lat_rr_ = l_lat[3];
		l_long_fl_ = l_long[0];
		l_long_fr_ = l_long[1];
		l_long_rl_ = l_long[2];
		l_long_rr_ = l_long[3];
		temp_linear[0] = k_body_fl_ + k_body_fr_ + k_body_rl_ + k_body_rr_; // K[i*DOF + 0] to be set later
		K[i * DOF + 1] = k_body_fl_ * l_lat_fl_ - k_body_fr_ * l_lat_fr_ + k_body_rl_ * l_lat_rl_ - k_body_rr_ * l_lat_rr_;
		K[i * DOF + 2] = -k_body_fl_ * l_long_fl_ - k_body_fr_ * l_long_fr_ + k_body_rl_ * l_long_rl_ + k_body_rr_ * l_long_rr_;
		K[i * DOF + 3] = -k_body_fl_;
		K[i * DOF + 4] = 0;
		K[i * DOF + 5] = -k_body_fr_;
		K[i * DOF + 6] = 0;
		K[i * DOF + 7] = -k_body_rl_;
		K[i * DOF + 8] = 0;
		K[i * DOF + 9] = -k_body_rr_;
		K[i * DOF + 10] = 0;

		i = 1;
		temp_linear[1] = l_lat_fl_ * l_lat_fl_ * k_body_fl_ + l_lat_fr_ * l_lat_fr_ * k_body_fr_ + l_lat_rl_ * l_lat_rl_ * k_body_rl_ + l_lat_rr_ * l_lat_rr_ * k_body_rr_; // K[i*DOF + 1] to be set later
		K[i * DOF + 2] = -l_long_fl_ * l_lat_fl_ * k_body_fl_ + l_lat_fr_ * l_long_fr_ * k_body_fr_ + l_long_rl_ * l_lat_rl_ * k_body_rl_ - l_long_rr_ * l_lat_rr_ * k_body_rr_;
		K[i * DOF + 3] = -l_lat_fl_ * k_body_fl_;
		K[i * DOF + 4] = 0;
		K[i * DOF + 5] = l_lat_fr_ * k_body_fr_;
		K[i * DOF + 6] = 0;
		K[i * DOF + 7] = -l_lat_rl_ * k_body_rl_;
		K[i * DOF + 8] = 0;
		K[i * DOF + 9] = l_lat_rr_ * k_body_rr_;
		K[i * DOF + 10] = 0;

		i = 2;
		temp_linear[2] = l_long_fl_ * l_long_fl_ * k_body_fl_ + l_long_fr_ * l_long_fr_ * k_body_fr_ + l_long_rl_ * l_long_rl_ * k_body_rl_ + l_long_rr_ * l_long_rr_ * k_body_rr_; // K[i*DOF + 2] to be set later
		K[i * DOF + 3] = l_long_fl_ * k_body_fl_;
		K[i * DOF + 4] = 0;
		K[i * DOF + 5] = l_long_fr_ * k_body_fr_;
		K[i * DOF + 6] = 0;
		K[i * DOF + 7] = -l_long_rl_ * k_body_rl_;
		K[i * DOF + 8] = 0;
		K[i * DOF + 9] = -l_long_rr_ * k_body_rr_;
		K[i * DOF + 10] = 0;

		i = 3;

		temp_linear[3] = k_body_fl_ + k_tyre_fl_; // K[i*DOF + 3]
		K[i * DOF + 4] = -k_tyre_fl_;
		K[i * DOF + 5] = 0;
		K[i * DOF + 6] = 0;
		K[i * DOF + 7] = 0;
		K[i * DOF + 8] = 0;
		K[i * DOF + 9] = 0;
		K[i * DOF + 10] = 0;
		// all others are zero

		i = 4;
		temp_linear[4] = k_tyre_fl_; //K[i*DOF + 4]
		K[i * DOF + 5] = 0;
		K[i * DOF + 6] = 0;
		K[i * DOF + 7] = 0;
		K[i * DOF + 8] = 0;
		K[i * DOF + 9] = 0;
		K[i * DOF + 10] = 0;

		i = 5;
		temp_linear[5] = k_body_fr_ + k_tyre_fr_; // K[i*DOF + 5]
		K[i * DOF + 6] = -k_tyre_fr_;
		K[i * DOF + 7] = 0;
		K[i * DOF + 8] = 0;
		K[i * DOF + 9] = 0;
		K[i * DOF + 10] = 0;

		i = 6;
		temp_linear[6] = k_tyre_fr_; // K[i*DOF + 6]
		K[i * DOF + 7] = 0;
		K[i * DOF + 8] = 0;
		K[i * DOF + 9] = 0;
		K[i * DOF + 10] = 0;

		i = 7;
		temp_linear[7] = k_body_rl_ + k_tyre_rl_; // K[i*DOF + 7]
		K[i * DOF + 8] = -k_tyre_rl_;
		K[i * DOF + 9] = 0;
		K[i * DOF + 10] = 0;


		i = 8;
		temp_linear[8] = k_tyre_rl_; // K[i*DOF + 8]
		K[i * DOF + 9] = 0;
		K[i * DOF + 10] = 0;

		i = 9;
		temp_linear[9] = k_body_rr_ + k_tyre_rr_; // K[i*DOF + 9]
		K[i * DOF + 10] = -k_tyre_rr_;

		i = 10;
		temp_linear[10] = k_tyre_rr_; // K[i*DOF + 10]

		// K=K+K'-diag(diag(K));
		cblas_dcopy(DOF * DOF, K, 1, K_trans, 1);
		mkl_dimatcopy('R', 'T', DOF, DOF, 1.0, K_trans, DOF, DOF); // get transpose of matrix
		cblas_daxpy(DOF * DOF, 1.0, K_trans, 1, K, 1); // K = K + K'
		MathLibrary::allocate_to_diagonal(K, temp_linear, DOF); // K = K + K'+ diag(K)
	}

	void construct_11DOF_vector(T* Global_position, T* Global_angle, T* Position_11dof) {
		Position_11dof[0] = Global_position[2]; // y coordinate of CG
		Position_11dof[1] = Global_angle[0]; // x angle of the CG
		Position_11dof[2] = Global_angle[1]; // y angle of the CG
		// copy y coordinate in order wheel, tyre, wheel, tyre, wheel, tyre, ...
		cblas_dcopy(vec_DIM - 1, Global_position + 5, 1, Position_11dof, 1);
	}

	void populate_K(T* k_vect, T k_body_fl, T k_tyre_fl, T k_body_fr,
		T k_tyre_fr,
		T k_body_rl,
		T k_tyre_rl,
		T k_body_rr,
		T k_tyre_rr) {

		/*
		This is needed when interpolation is turned off
		*/
		k_vect[0] = k_body_fl;
		k_vect[1] = k_tyre_fl;
		k_vect[2] = k_body_fr;
		k_vect[3] = k_tyre_fr;
		k_vect[4] = k_body_rl;
		k_vect[5] = k_tyre_rl;
		k_vect[6] = k_body_rr;
		k_vect[7] = k_tyre_rr;
	}
	void get_length(
		T* initial_orientation_,
		const T* r1_,
		const T* r2_,
		const T* r3_,
		const T* r4_,
		const T* pcc_,
		const T* initial_upper_spring_length_,
		const T* initial_lower_spring_length_,
		T* wheel_coordinate1_,
		T* wheel_coordinate2_,
		T* wheel_coordinate3_,
		T* wheel_coordinate4_,
		T* tyre_coordinate1_,
		T* tyre_coordinate2_,
		T* tyre_coordinate3_,
		T* tyre_coordinate4_)

	{
		/*
		To reduce memory trace and better use cache this function is implemented in following fashion:
		Original steps for computation of one component:
														1.	qc = qc/norm(qc);
														2.	C_Nc = get_basis(qc);
														3.	global_y = C_Nc(:,2);
														4.	global_y = -global_y / norm(global_y);
														5.	global_r1 = pcc + C_Nc*r1;
														6.	upper_global_spring_1 = upper_length(1)*global_y;
														7.	lower_global_spring_1 = lower_length(1)*global_y;
														8.	pw1 = global_r1 + upper_global_spring_1;
														9.	pt1 = pw1 + lower_global_spring_1;
		Modified steps for computation of one component:
														1.	qc = qc/norm(qc);
														2.	C_Nc = get_basis(qc);
														3.	global_y = C_Nc(:,2);
														4.	global_y = -global_y / norm(global_y);
														5.	pw1 = pcc;
														6.	pw1 = pw1 + C_Nc*r1;
														7.	pw1 = pw1 + upper_length(1)*global_y;
														8.	pt1 = pw1
														8.	pt1	= pt1 + lower_length(1)*global_y;
		*/

		T* global_y = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		T* C_Nc = (T*)mkl_calloc((this->DIM) * (this->DIM), sizeof(T), this->alignment);

		//	1. qc = qc/norm(qc); This is in quaternions 
		T nrm = cblas_dnrm2(this->NUM_LEGS, initial_orientation_, 1);
		cblas_dscal(this->NUM_LEGS, 1.0 / nrm, initial_orientation_, 1);

		// 2.	C_Nc = get_basis(qc);
		MathLibrary::get_basis<T>(initial_orientation_, C_Nc);
		// 3.	global_y = C_Nc(:,2);
		cblas_dcopy(this->DIM, C_Nc + 1, this->DIM, global_y, 1);

		// 4.	global_y = -global_y / norm(global_y);
		nrm = cblas_dnrm2(this->DIM, global_y, 1);
		cblas_dscal(this->DIM, -1.0 / nrm, global_y, 1);

		/////////////////////////////////////////// Leg 1 ////////////////////////////////////////////////////////
		// 5.	pw1 = pcc;
		cblas_dcopy(this->DIM, pcc_, 1, wheel_coordinate1_, 1);

		// 6.	pw1 = pw1 + C_Nc*r1;
		cblas_dgemv(CblasRowMajor, CblasNoTrans, this->DIM, this->DIM, 1, C_Nc, this->DIM, r1_, 1, 1, wheel_coordinate1_, 1);

		// 7.	pw1 = pw1 + upper_length(1)*global_y;
		cblas_daxpy(this->DIM, initial_upper_spring_length_[0], global_y, 1, wheel_coordinate1_, 1);

		// 8.	pt1 = pw1
		cblas_dcopy(this->DIM, wheel_coordinate1_, 1, tyre_coordinate1_, 1);

		// 9.	pt1 = pw1 + lower_length(1)*global_y;
		cblas_daxpy(this->DIM, initial_lower_spring_length_[0], global_y, 1, tyre_coordinate1_, 1);

		/////////////////////////////////////////// Leg 2 ////////////////////////////////////////////////////////
		// 5.	pw2 = pcc;
		cblas_dcopy(this->DIM, pcc_, 1, wheel_coordinate2_, 1);

		// 6.	pw2 = pw2 + C_Nc*r2;
		cblas_dgemv(CblasRowMajor, CblasNoTrans, this->DIM, this->DIM, 1, C_Nc, this->DIM, r2_, 1, 1, wheel_coordinate2_, 1);

		// 7.	pw2 = pw2 + upper_length(2)*global_y;
		cblas_daxpy(this->DIM, initial_upper_spring_length_[1], global_y, 1, wheel_coordinate2_, 1);

		// 8.	pt2 = pw2
		cblas_dcopy(this->DIM, wheel_coordinate2_, 1, tyre_coordinate2_, 1);

		// 9.	pt2 = pw2 + lower_length(2)*global_y;
		cblas_daxpy(this->DIM, initial_lower_spring_length_[1], global_y, 1, tyre_coordinate2_, 1);

		/////////////////////////////////////////// Leg 3 ////////////////////////////////////////////////////////
		// 5.	pw3 = pcc;
		cblas_dcopy(this->DIM, pcc_, 1, wheel_coordinate3_, 1);

		// 6.	pw3 = pw3 + C_Nc*r3;
		cblas_dgemv(CblasRowMajor, CblasNoTrans, this->DIM, this->DIM, 1, C_Nc, this->DIM, r3_, 1, 1, wheel_coordinate3_, 1);

		// 7.	pw3 = pw3 + upper_length(3)*global_y;
		cblas_daxpy(this->DIM, initial_upper_spring_length_[2], global_y, 1, wheel_coordinate3_, 1);

		// 8.	pt3 = pw3
		cblas_dcopy(this->DIM, wheel_coordinate3_, 1, tyre_coordinate3_, 1);

		// 9.	pt3 = pw3 + lower_length(3)*global_y;
		cblas_daxpy(this->DIM, initial_lower_spring_length_[2], global_y, 1, tyre_coordinate3_, 1);

		/////////////////////////////////////////// Leg 4 ////////////////////////////////////////////////////////
		// 5.	pw4 = pcc;
		cblas_dcopy(this->DIM, pcc_, 1, wheel_coordinate4_, 1);

		// 6.	pw4 = pw4 + C_Nc*r4;
		cblas_dgemv(CblasRowMajor, CblasNoTrans, this->DIM, this->DIM, 1, C_Nc, this->DIM, r4_, 1, 1, wheel_coordinate4_, 1);

		// 7.	pw4 = pw4 + upper_length(4)*global_y;
		cblas_daxpy(this->DIM, initial_upper_spring_length_[3], global_y, 1, wheel_coordinate4_, 1);

		// 8.	pt4 = pw4
		cblas_dcopy(this->DIM, wheel_coordinate4_, 1, tyre_coordinate4_, 1);

		// 9.	pt4 = pw4 + lower_length(4)*global_y;
		cblas_daxpy(this->DIM, initial_lower_spring_length_[3], global_y, 1, tyre_coordinate4_, 1);


		mkl_free(global_y);
		mkl_free(C_Nc);
	}

	inline void compute_dx(const T* current_length, T* dx) {
		/*
		the dx follows the order
		[w1, t1, w2, t2, w3, t3, w4, t4]
		current_length has input of type
		[w1, t1, w2, t2, w3, t3, w4, t4]

		*/

		// vdSub(n, a, b, y);  <---> y = a - b: 
		cblas_dcopy(2 * (this->num_tyre), spring_length, current_length, dx);
	}

	inline void compute_dx(T* dx) {
		compute_dx(current_spring_length, dx);
	}

	void get_Position_vec(T* Pos) {
		if (Pos != NULL) {
			cblas_dcopy(DIM * vec_DIM, Position_vec, 1, Pos, 1);
		}
	}
	void set_Position_vec(const T* Pos) {
		if (Pos != NULL) {
			cblas_dcopy(DIM * vec_DIM, Pos, 1, Position_vec, 1);
		}
	}
	void get_Position_vec_CG(T* Pos_CG) {
		if (Pos_CG != NULL) {
			cblas_dcopy(DIM, Position_vec, 1, Pos_CG, 1);
		}
	}
	void set_Position_vec_CG(const T* Pos_CG) {
		if (Pos_CG != NULL) {
			cblas_dcopy(DIM, Pos_CG, 1, Position_vec, 1);
		}
	}

	void get_Velocity_vec(T* Vel) {
		if (Vel != NULL) {
			// b=a, cblas_dcopy(n,a,inc,b,inc)
			cblas_dcopy(DIM * vec_DIM, Velocity_vec, 1, Vel, 1);
		}
	}
	void set_Velocity_vec(const T* Vel) {
		if (Vel != NULL) {
			cblas_dcopy(DIM * vec_DIM, Vel, 1, Velocity_vec, 1);
		}
	}
	void get_Velocity_vec_CG(T* Vel_CG) {
		if (Vel_CG != NULL) {
			// b=a, cblas_dcopy(n,a,inc,b,inc)
			cblas_dcopy(DIM, Velocity_vec, 1, Vel_CG, 1);
		}
	}
	void set_Velocity_vec_CG(const T* Vel_CG) {
		if (Vel_CG != NULL) {
			cblas_dcopy(DIM, Vel_CG, 1, Velocity_vec, 1);
		}
	}

	void get_k_vec(T* k) {
		if (k != NULL) {
			// b=a, cblas_dcopy(n,a,inc,b,inc)
			cblas_dcopy(DIM * vec_DIM, k_vec, 1, k, 1);
		}
	}
	void set_k_vec(const T* k) {
		if (k != NULL) {
			cblas_dcopy(DIM * vec_DIM, k, 1, k_vec, 1);
		}
	}

	void get_Mass_vec(T* M) {
		if (M != NULL) {
			// b=a, cblas_dcopy(n,a,inc,b,inc)
			cblas_dcopy(DIM * vec_DIM, Mass_vec, 1, M, 1);
		}
	}
	double get_Mass_vec_CG() const {
		return Mass_vec[0];
	}
	void set_Mass_vec(const T* M) {
		if (M != NULL) {
			cblas_dcopy(DIM * vec_DIM, M, 1, Mass_vec, 1);
		}
	}

	void get_dist_vector(T* Point_P, T* dist_vector) {
		// get distance vector from each important Point of the car (9: CG, 4*W_i, 4*T_i)
	// source: Point_P, dest: each entry from Position_vec
		if (Point_P != NULL && dist_vector != NULL) {
			for (auto i = 0; i < vec_DIM; ++i) {
				cblas_dcopy(DIM, Point_P, incx, &dist_vector[DIM * i], incx);
			}
			// y=a-b, vdSub(n,a,b,y)
			vdSub(DIM * vec_DIM, Position_vec, dist_vector, dist_vector);
		}
	
	} // 9 * 3 - from each important point to a fixed Point_P
	void get_dist_vector_CG(T* Point_P, T* dist_vector) {
		// get distance vector from Center of Gravity of the car to a Point P 
		// source: Point_P, dest: CG
		if (Point_P != NULL && dist_vector != NULL) {
			// y=a-b, vdSub(n,a,b,y)
			vdSub(DIM, Position_vec, Point_P, dist_vector);
		}
	} // 3 - from Center of Gravity of a fixed Point_P

	double get_I_body_xx() const {
		return I_CG[0];
	}
	void set_I_body_xx(const T& I_body_xx_val) {
		I_CG[0] = I_body_xx_val;
	}
	double get_I_body_yy() const {
		return I_CG[4];
	}
	void set_I_body_yy(const T& I_body_yy_val) {
		I_CG[4] = I_body_yy_val;
	}
	~Car() {
		mkl_free(Position_vec); 
		mkl_free(Velocity_vec);
		mkl_free(Mass_vec);
		mkl_free(angle_CG);
		mkl_free(w_CG);
		mkl_free(I_CG);


		// Initial Conditions of the car
		mkl_free(initial_position);
		mkl_free(initial_velocity_vec); 
		mkl_free(initial_angle);
		mkl_free(initial_angular_velocity);

		
		
		mkl_free(quad_angle_init);
		mkl_free(M_linear); 
		mkl_free(temp_linear); 
		mkl_free(K);
		mkl_free(K_trans);
		mkl_free(D);
		mkl_free(spring_length); 
		mkl_free(current_spring_length);
		mkl_free(u_prev_linear); 
		mkl_free(u_current_linear);
		mkl_free(k_vec); 
		mkl_free(l_lat); 
		mkl_free(l_long);
		mkl_free(velocity_current_linear);

		////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////// ALE Vectors ///////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////////////////
		mkl_free(Position_vec_xy);
		mkl_free(Velocity_vec_xy);
		delete Angle_z;
		delete w_z;
	}
};