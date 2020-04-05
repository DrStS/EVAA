#pragma once
#include <mkl.h>
#include "ReadXML.h"
#include "MathLibrary.h"
#include <cmath>

template <class T>
class Car {
private:


	/*
	Create the diagonal mass matrix M_linear to solve the 11DOF system
	\param Global_mass contains the masses of the 9 bodies (CG, 4 * W, 4 * T)
	\param Global_moment_Inertia contains the tensor of inertia of the body
	\param Mass_11DOF unused
	*/
	void construct_11DOF_mass(T* Global_mass, T* Global_momemnt_Inertia) {
		temp_linear[0] = Global_mass[0];
		temp_linear[1] = Global_momemnt_Inertia[0];
		temp_linear[2] = Global_momemnt_Inertia[4];
		cblas_dcopy(vec_DIM - 1, Global_mass + 1, 1, temp_linear + 3, 1);
		MathLibrary::allocate_to_diagonal(M_linear, temp_linear, DOF);
	}
	/*
	Copy all X and Y coordinates of the global vector to the local vector
	\param Global_vector vector with coordinates [X,Y,Z,X,Y,Z,...]
	\param local_vector vector with coordinates [X,Y,X,Y,...]
	*/
	void construct_ALE_vectors(T* Global_vector, T* local_vector) {
		T* start_pointer, *current_ptr;
		start_pointer = Global_vector; // copy x and y and move next
		current_ptr = local_vector;
		for (size_t i = 0; i < vec_DIM - 1; ++i) {
			cblas_dcopy(DIM - 1, start_pointer, 1, current_ptr, 1);
			start_pointer += DIM;
			current_ptr += DIM - 1;
		}
		cblas_dcopy(DIM - 1, start_pointer, 1, current_ptr, 1);
	}
	/*
	Construct corner initilizer
	*/
	void construct_corner(T* pos_CG, T* corners) {
		T c, s;
		c = std::cos(angle_CG[2]);
		s = std::sin(angle_CG[2]);
		corners[0] = pos_CG[0] + l_long[0] * c - l_lat[0] * s; // fl
		corners[4] = pos_CG[1] + l_lat[0] * c + l_long[0]*s ; // fl
		corners[8] = pos_CG[2];
		corners[1] = pos_CG[0] + l_long[1] * c + l_lat[1] * s; // fr
		corners[5] = pos_CG[1] - l_lat[1] * c + l_long[1] * s; // fr
		corners[9] = pos_CG[2];
		corners[2] = pos_CG[0] - l_long[2]*c - l_lat[2]*s; // rl
		corners[6] = pos_CG[1] + l_lat[2]*c - l_long[2]*s; // rl
		corners[8] = pos_CG[2];
		corners[3] = pos_CG[0] - l_long[3]*c + l_lat[3]*s; // rr
		corners[7] = pos_CG[1] - l_lat[3]*c - l_long[3]*s; // rr
		corners[11] = pos_CG[2];
	}
	/*
	Calculates the values of Corners for general angles
	*/
	void update_corners_11DOF(T* angles, T* rotation_mat_buffer, T* initial_corners, T* updated_corners)
	{
		// zz, yy, xx
		MathLibrary::get_rotation_matrix(angles[2], angles[1], angles[0], rotation_mat_buffer);

		// do rotation: rotationMat * r
		//void cblas_dgemm(const CBLAS_LAYOUT Layout, const CBLAS_TRANSPOSE transa, const CBLAS_TRANSPOSE transb, const MKL_INT m, const MKL_INT n, const MKL_INT k, const double alpha, const double* a, const MKL_INT lda, const double* b, const MKL_INT ldb, const double beta, double* c, const MKL_INT ldc);
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, DIM, num_wheels, DIM, 1, rotation_mat_buffer, DIM, initial_corners, num_wheels, 0, updated_corners, num_wheels);
	}
	/*
	Calculates the values of Corners_current according to the current orientation
	*/
	void update_corners_11DOF() {
		pos_buffer[0] = Position_vec_xy[0];
		pos_buffer[1] = Position_vec_xy[1];
		pos_buffer[2] = u_current_linear[0];
		angle_buffer[0] = u_current_linear[1];
		angle_buffer[1] = u_current_linear[2];
		angle_buffer[3] = 0;
		construct_corner(pos_buffer, Corners_init);
		update_corners_11DOF(angle_buffer, Corners_rot, Corners_init, Corners_current);
	}
	void set_ALE2global(T* vect, T* global_vect) {
#pragma loop(ivdep)
		for (size_t i = 0; i < vec_DIM; ++i) {
			global_vect[DIM*i] = vect[(DIM - 1)*i];
			global_vect[DIM*i + 1] = vect[(DIM - 1)*i + 1];
		}
	}
	void set_11DOF2global(T* vect, T* global_vect) {
		global_vect[2] = vect[0];
#pragma loop(ivdep)
		for (size_t i = 1; i < vec_DIM; ++i) {
			global_vect[DIM*i + 2] += vect[(DIM - 1) + i];
		}
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
	const size_t malloc_factor = 1;
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
	///////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////// Interpolator Members ///////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////
	T *Corners_rot, *Corners_current, *Corners_init;
	T* angle_buffer, *pos_buffer;

	////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////// ALE Vectors ///////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////
	T *Position_vec_xy, *Angle_z, *Velocity_vec_xy, *w_z;


	/*
	Constructor
	*/
	Car(const Simulation_Parameters &params, EVAAComputeStiffness* interpolator) {
		std::cout << "I am the flying car" << std::endl;
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
		Position_vec = (T*)mkl_malloc(malloc_factor*DIM*vec_DIM * sizeof(T), alignment); //27 dim
		Velocity_vec = (T*)mkl_malloc(malloc_factor*DIM*vec_DIM * sizeof(T), alignment); //27 dim
		Mass_vec = (T*)mkl_malloc(malloc_factor*vec_DIM * sizeof(T), alignment); //9 dim
		angle_CG = (T*)mkl_malloc(malloc_factor*DIM * sizeof(T), alignment); //3 dim
		w_CG = (T*)mkl_malloc(malloc_factor*DIM * sizeof(T), alignment); //3 dim
		I_CG = (T*)mkl_malloc(malloc_factor*DIM*DIM * sizeof(T), alignment); //9 dim

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////// Initial Params for Global coordinate ////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


		initial_position = (T*)mkl_malloc(malloc_factor*DIM*vec_DIM * sizeof(T), alignment); // 27 dim
		initial_velocity_vec = (T*)mkl_malloc(malloc_factor*DIM*vec_DIM * sizeof(T), alignment);// 27 dim
		quad_angle_init = (T*)mkl_calloc(malloc_factor * 4, sizeof(T), alignment); // 4 dim
		initial_angle = (T*)mkl_malloc(malloc_factor*DIM * sizeof(T), alignment); // 3 dim
		initial_angular_velocity = (T*)mkl_malloc(malloc_factor*DIM * sizeof(T), alignment); // 3 dim

		


		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////// Params for 11 DOF system ////////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		M_linear = (T*)mkl_calloc(malloc_factor*DOF * DOF, sizeof(T), alignment); // 121 dim
		temp_linear = (T*)mkl_calloc(malloc_factor*DOF , sizeof(T), alignment);  // 11 dim
		K = (T*)mkl_calloc(malloc_factor*DOF * DOF, sizeof(T), alignment); // 121 dim
		K_trans = (T*)mkl_calloc(malloc_factor*DOF * DOF, sizeof(T), alignment); // 121 dim
		D = (T*)mkl_calloc(malloc_factor*DOF * DOF, sizeof(T), alignment); // 121 dim
		u_prev_linear = (T*)mkl_malloc(malloc_factor*DOF*sizeof(T), alignment); // 11 dim
		u_current_linear = (T*)mkl_malloc(malloc_factor*DOF*sizeof(T), alignment); // 11 dim
		velocity_current_linear = (T*)mkl_malloc(malloc_factor*DOF * sizeof(T), alignment); // 3 dim
		k_vec = (T*)mkl_malloc(malloc_factor*num_wheels * 2 * sizeof(T), alignment); // 4 dim
		l_lat = (T*)mkl_malloc(malloc_factor*num_wheels * sizeof(T), alignment); // 4 dim
		l_long = (T*)mkl_malloc(malloc_factor*num_wheels * sizeof(T), alignment); // 4 dim
		spring_length = (T*)mkl_malloc(malloc_factor * 2 * num_wheels * sizeof(T), alignment); // 8 dim
		current_spring_length = (T*)mkl_malloc(malloc_factor * 2 * num_wheels * sizeof(T), alignment); // 8 dim

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////// Memory allocation for interpolator /////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		Corners_current = (T*)mkl_malloc(malloc_factor*num_tyre*DIM * sizeof(T), alignment); // 12 dim
		Corners_rot = (T*)mkl_malloc(malloc_factor*DIM*DIM * sizeof(T), alignment); // 9 dim
		Corners_init = (T*)mkl_malloc(malloc_factor*num_tyre*DIM * sizeof(T), alignment); // 12 dim
		angle_buffer = (T*)mkl_malloc(malloc_factor*DIM * sizeof(T), alignment); // 3 dim
		pos_buffer = (T*)mkl_malloc(malloc_factor*DIM * sizeof(T), alignment); // 3 dim

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////// ALE Buffer Allocation /////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		Position_vec_xy = (T*)mkl_malloc((DIM - 1) * vec_DIM * sizeof(T), alignment);
		Angle_z = new T;
		Velocity_vec_xy = (T*)mkl_malloc((DIM - 1) * vec_DIM * sizeof(T), alignment);
		w_z = new T;

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


		/// Good indices???
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
		Mass_vec[1] = params.mass_wheel[2];
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
		
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////// Interpolator Initialization///////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// read init corners vectors into matrix

		construct_corner(initial_position, Corners_init); // only CG position is used to construct corners
		cblas_dcopy(DIM, angle_CG, 1, angle_buffer, 1);
		angle_buffer[2] = 0;
		update_corners_11DOF(angle_buffer, Corners_rot, Corners_init, Corners_current);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////// Remaining position initialization ////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		const T* xml_start;
		T *position_start;
		if (params.initial_leg_flag) {
			// if prescribed initial position (add a check for consistency with spring lengths)
			// W1 = W_fl
			xml_start = params.initial_pos_wheel + 2 * 3;
			position_start = initial_position + 3; //(end at 5)
			cblas_dcopy(DIM, xml_start, 1, position_start, 1);
			// W2 = W_fr
			xml_start = params.initial_pos_wheel + 3 * 3;
			position_start += 6; // skip 3 for tyre (end at 11)
			cblas_dcopy(DIM, xml_start, 1, position_start, 1);
			// W3 = W_rl
			xml_start = params.initial_pos_wheel + 1 * 3;
			position_start += 6; // skip 3 for tyre (end at 17)
			cblas_dcopy(DIM, xml_start, 1, position_start, 1);
			// W2 = W_rr
			xml_start = params.initial_pos_wheel + 0 * 3;
			position_start += 6; // skip 3 for tyre (end at 23)
			cblas_dcopy(DIM, xml_start, 1, position_start, 1);

			// T1 = T_fl
			xml_start = params.initial_pos_tyre + 2 * 3;
			position_start = initial_position + 6; // skip 3 for center of mass and 3 for the wheel
			cblas_dcopy(DIM, xml_start, 1, position_start, 1); // (end at 8)
			// T2 = T_fr
			xml_start = params.initial_pos_tyre + 3 * 3;
			position_start += 6; // skip 3 for the wheel
			cblas_dcopy(DIM, xml_start, 1, position_start, 1); // (end at 14)
			// T3 = T_rl
			xml_start = params.initial_pos_tyre + 1 * 3;
			position_start += 6; // skip 3 for the wheel
			cblas_dcopy(DIM, xml_start, 1, position_start, 1);// (end at 20)
			// T4 = T_rr
			xml_start = params.initial_pos_tyre + 0 * 3;
			position_start += 6; // skip 3 for the wheel
			cblas_dcopy(DIM, xml_start, 1, position_start, 1); // (end at 26)
		} 
		else {
			T* W_fl = initial_position + 3;
			T* W_fr = initial_position + 9;
			T* W_rl = initial_position + 15;
			T* W_rr = initial_position + 21;
			T* T_fl = initial_position + 6;
			T* T_fr = initial_position + 12;
			T* T_rl = initial_position + 18;
			T* T_rr = initial_position + 24;
			get_length(Corners_current, current_spring_length, W_fl, T_fl, W_fr, T_fr, W_rl, T_rl, W_rr, T_rr);
		}
		// copy the initial position to the position vector
		cblas_dcopy(DIM*vec_DIM, initial_position, 1, Position_vec, 1);
		
		

		// Initial Velocity (Reuse the pointers)
		cblas_dcopy(DIM, params.initial_vel_body, 1, initial_velocity_vec, 1);
		// W1 = W_fl
		xml_start = params.initial_vel_wheel + 2 * 3;
		position_start = initial_velocity_vec + 3;
		cblas_dcopy(DIM, xml_start, 1, position_start, 1); // (end at 5)
		// W2 = W_fr
		xml_start = params.initial_vel_wheel + 3 * 3;
		position_start += 6; // skip 3 for tyre
		cblas_dcopy(DIM, xml_start, 1, position_start, 1);// (end at 11)
		// W3 = W_rl
		xml_start = params.initial_vel_wheel + 1 * 3;
		position_start += 6; // skip 3 for tyre
		cblas_dcopy(DIM, xml_start, 1, position_start, 1); // (end at 17)
		// W2 = W_rr
		xml_start = params.initial_vel_wheel + 0 * 3;
		position_start += 6; // skip 3 for tyre
		cblas_dcopy(DIM, xml_start, 1, position_start, 1); // (end at 23)

		// T1 = T_fl
		xml_start = params.initial_vel_tyre + 2 * 3;
		position_start = initial_velocity_vec + 6; // skip 3 for center of mass and 3 for the wheel
		cblas_dcopy(DIM, xml_start, 1, position_start, 1); // (end at 8)
		// T2 = T_fr
		xml_start = params.initial_vel_tyre + 3 * 3;
		position_start += 6; // skip 3 for the wheel
		cblas_dcopy(DIM, xml_start, 1, position_start, 1);// (end at 14)
		// T3 = T_rl
		xml_start = params.initial_vel_tyre + 1 * 3;
		position_start += 6; // skip 3 for the wheel
		cblas_dcopy(DIM, xml_start, 1, position_start, 1); // (end at 20)
		// T4 = T_rr
		xml_start = params.initial_vel_tyre + 0 * 3;
		position_start += 6; // skip 3 for the wheel
		cblas_dcopy(DIM, xml_start, 1, position_start, 1); // (end at 26)

		// copy the initial position to the position vector
		cblas_dcopy(DIM*vec_DIM, initial_velocity_vec, 1, Velocity_vec, 1);

		// Initial Angular velocity
		cblas_dcopy(DIM, params.initial_ang_vel_body, 1, initial_angular_velocity, 1); 
		cblas_dcopy(DIM, initial_angular_velocity, 1, w_CG, 1);
		
		/*
		Global assignments Done 
		11 DOF assignments
		*/
		construct_11DOF_mass(Mass_vec, I_CG);
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
		update_K(k_vec);

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////// ALE Buffer Initialization /////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		construct_ALE_vectors(Position_vec, Position_vec_xy);
		construct_ALE_vectors(Velocity_vec, Velocity_vec_xy);
		*Angle_z = angle_CG[2];
		*w_z = w_CG[2];
		


	}
	/*
	If forces and positiosn are negative, set them to zero, elsewise, keep them
	*/
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

	/*
	compute a reaction which is opposite to the internal force acting on the tyre
	*/
	void compute_normal_force(T* K, T* u, T* force, size_t* index, size_t dim, size_t n) {
		#pragma loop( ivdep )
		for (int i = 0; i < n; ++i) {
			force[index[i]] = -K[index[i] * dim + index[i]] * u[index[i]];
		}
	}
	/*
	Calculate the entries of the stiffness matrix
	\param k_vect vector with all spring stiffnesses (in Stefans order)
	*/
	void update_K(T* k_vect) {
		cblas_dscal(DOF * DOF, 0.0, K, 1);

		temp_linear[0] = k_vect[0] + k_vect[2] + k_vect[4] + k_vect[6];
		K[1] = k_vect[0] * l_lat[0] - k_vect[2] * l_lat[1] + k_vect[4] * l_lat[2] - k_vect[6] * l_lat[3];
		K[2] = -k_vect[0] * l_long[0] - k_vect[2] * l_long[1] + k_vect[4] * l_long[2] + k_vect[6] * l_long[3];
		K[3] = -k_vect[0];
		K[4] = 0.0;
		K[5] = -k_vect[2];
		K[6] = 0.0;
		K[7] = -k_vect[4];
		K[8] = 0.0;
		K[9] = -k_vect[6];
		K[10] = 0.0;

		temp_linear[1] = l_lat[0] * l_lat[0] * k_vect[0] + l_lat[1] * l_lat[1] * k_vect[2] + l_lat[2] * l_lat[2] * k_vect[4] + l_lat[3] * l_lat[3] * k_vect[6];
		K[DOF + 2] = -l_long[0] * l_lat[0] * k_vect[0] + l_lat[1] * l_long[1] * k_vect[2] + l_long[2] * l_lat[2] * k_vect[4] - l_long[3] * l_lat[3] * k_vect[6];
		K[DOF + 3] = -l_lat[0] * k_vect[0];
		K[DOF + 4] = 0;
		K[DOF + 5] = l_lat[1] * k_vect[2];
		K[DOF + 6] = 0;
		K[DOF + 7] = -l_lat[2] * k_vect[4];
		K[DOF + 8] = 0;
		K[DOF + 9] = l_lat[3] * k_vect[6];
		K[DOF + 10] = 0;

		temp_linear[2] = l_long[0] * l_long[0] * k_vect[0] + l_long[1] * l_long[1] * k_vect[2] + l_long[2] * l_long[2] * k_vect[4] + l_long[3] * l_long[3] * k_vect[6];
		K[2 * DOF + 3] = l_long[0] * k_vect[0];
		K[2 * DOF + 4] = 0;
		K[2 * DOF + 5] = l_long[1] * k_vect[2];
		K[2 * DOF + 6] = 0;
		K[2 * DOF + 7] = -l_long[2] * k_vect[4];
		K[2 * DOF + 8] = 0;
		K[2 * DOF + 9] = -l_long[3] * k_vect[6];
		K[2 * DOF + 10] = 0;

		temp_linear[3] = k_vect[0] + k_vect[1];
		K[3 * DOF + 4] = -k_vect[1];
		K[3 * DOF + 5] = 0;
		K[3 * DOF + 6] = 0;
		K[3 * DOF + 7] = 0;
		K[3 * DOF + 8] = 0;
		K[3 * DOF + 9] = 0;
		K[3 * DOF + 10] = 0;
		// all others are zero

		temp_linear[4] = k_vect[1];
		K[4 * DOF + 5] = 0;
		K[4 * DOF + 6] = 0;
		K[4 * DOF + 7] = 0;
		K[4 * DOF + 8] = 0;
		K[4 * DOF + 9] = 0;
		K[4 * DOF + 10] = 0;

		temp_linear[5] = k_vect[2] + k_vect[3];
		K[5 * DOF + 6] = -k_vect[3];
		K[5 * DOF + 7] = 0;
		K[5 * DOF + 8] = 0;
		K[5 * DOF + 9] = 0;
		K[5 * DOF + 10] = 0;

		temp_linear[6] = k_vect[3];
		K[6 * DOF + 7] = 0;
		K[6 * DOF + 8] = 0;
		K[6 * DOF + 9] = 0;
		K[6 * DOF + 10] = 0;

		temp_linear[7] = k_vect[4] + k_vect[5];
		K[7 * DOF + 8] = -k_vect[5];
		K[7 * DOF + 9] = 0;
		K[7 * DOF + 10] = 0;


		temp_linear[8] = k_vect[5];
		K[8 * DOF + 9] = 0;
		K[8 * DOF + 10] = 0;


		temp_linear[9] = k_vect[6] + k_vect[7];
		K[9 * DOF + 10] = -k_vect[7];

		temp_linear[10] = k_vect[7];
		
		cblas_dcopy(DOF * DOF, K, 1, K_trans, 1);
		mkl_dimatcopy('R', 'T', DOF, DOF, 1.0, K_trans, DOF, DOF); // get transpose of matrix
		cblas_daxpy(DOF * DOF, 1.0, K_trans, 1, K, 1); // K = K + K'
		MathLibrary::allocate_to_diagonal(K, temp_linear, DOF); // K = K + K'+ diag(K)
	}

	/*
	Get the solution vector as required for the 11DOF system
	\param Global_position in the format [GC:XYZ,W1:XYZ,T1:XYZ,W2:XYZ,T2:XYZ,...]
	\param Global_angle with three angles [X,Y,Z]
	\return Position_11dof in the format [GC:Y,angle:XY,W1:Y,T1:Y,W2:Y,T2:Y,...]
	*/
	void construct_11DOF_vector(T* Global_position, T* Global_angle, T* Position_11dof) {
		Position_11dof[0] = Global_position[2]; // z coordinate of CG
		Position_11dof[1] = Global_angle[0]; // x angle of the CG
		Position_11dof[2] = Global_angle[1]; // y angle of the CG

		// copy y coordinate in order wheel, tyre, wheel, tyre, wheel, tyre, ...
		cblas_dcopy(vec_DIM - 1, Global_position + 5, 3, Position_11dof + 3, 1);
	}
	/*
	fill the vector with all stiffness with the constant values from the XML (if the lookup table is not used)
	\param k_tyre_** stiffnesses of the lower springs
	\param k_body_** stiffnesses of the upper springs
	\return k_vect with all stiffnesses in Stefan's ordering
	*/
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
	void get_length(T* Corners, T* curr_spring_len, T* W_fl, T* T_fl, T* W_fr, T* T_fr, T* W_rl, T* T_rl, T* W_rr, T* T_rr) {
		// W_fl & T_fl
		cblas_dcopy(DIM, Corners, 4, W_fl, 1);
		W_fl[2] -= curr_spring_len[0];
		cblas_dcopy(DIM, W_fl, 1, T_fl, 1);
		T_fl[2] -= curr_spring_len[1];

		// W_fr & T_fr
		cblas_dcopy(DIM, Corners + 1, 4, W_fr, 1);
		W_fr[2] -= curr_spring_len[2];
		cblas_dcopy(DIM, W_fr, 1, T_fr, 1);
		T_fr[2] -= curr_spring_len[3];

		// W_rl & T_rl
		cblas_dcopy(DIM, Corners + 2, 4, W_rl, 1);
		W_rl[2] -= curr_spring_len[4];
		cblas_dcopy(DIM, W_rl, 1, T_rl, 1);
		T_rl[2] -= curr_spring_len[5];

		// W_rr & T_rr
		cblas_dcopy(DIM, Corners + 3, 4, W_rr, 1);
		W_rr[2] -= curr_spring_len[6];
		cblas_dcopy(DIM, W_rr, 1, T_rr, 1);
		T_rr[2] -= curr_spring_len[7];
	}


	inline void compute_dx(const T* current_length, T* dx) {
		/*
		the dx follows the order
		[w1, t1, w2, t2, w3, t3, w4, t4]
		current_length has input of type
		[w1, t1, w2, t2, w3, t3, w4, t4]

		*/

		// vdSub(n, a, b, y);  <---> y = a - b: 
		vdSub(2 * (this->num_tyre), spring_length, current_length, dx);
	}
	/*
	From the current elongations, calucale the difference to the rest position
	\return length differences
	*/
	inline void compute_dx(T* dx) {
		compute_dx(current_spring_length, dx);
	}
	/*
	First updates the corner and afterwards compute the lengths of the springs
	*/
	void update_lengths_11DOF() {
		update_corners_11DOF();
		current_spring_length[0] = std::abs(Corners_current[8] - u_current_linear[3]);
		current_spring_length[1] = std::abs(u_current_linear[3] - u_current_linear[4]);
		current_spring_length[2] = std::abs(Corners_current[9] - u_current_linear[5]);
		current_spring_length[3] = std::abs(u_current_linear[5] - u_current_linear[6]);
		current_spring_length[4] = std::abs(Corners_current[10] - u_current_linear[7]);
		current_spring_length[5] = std::abs(u_current_linear[7] - u_current_linear[8]);
		current_spring_length[6] = std::abs(Corners_current[11] - u_current_linear[9]);
		current_spring_length[7] = std::abs(u_current_linear[9] - u_current_linear[10]);
	}

	/* Fills the global vector with all entries
	\param ALE_vectors contains X and Y components [GC:XY,W1:XY,T1:XY,W2:XY,T2:XY,...]
	\param vector 11DOF contains Z components [GC:Z,W1:Z,T1:Z,W2:Z,T2:Z,...]
	\return global_vector [GC:XYZ,W1:XYZ,T1:XYZ,W2:XYZ,T2:XYZ,...]
	*/
	void populate_results(T* ALE_vector, T * vector_11DOF, T* global_vector) {
		set_ALE2global(ALE_vector, global_vector);
		set_11DOF2global(vector_11DOF, global_vector);

	}
	

	void combine_results() {
		set_ALE2global(Position_vec_xy, Position_vec);
		set_11DOF2global(u_current_linear, Position_vec);
		// Angles manually
		angle_CG[0] = u_current_linear[1];
		angle_CG[1] = u_current_linear[2];
		angle_CG[2] = *Angle_z;
		w_CG[2] = *w_z;
		set_ALE2global(Velocity_vec_xy, Velocity_vec);
	}

	void get_Position_vec(T* Pos) {
		cblas_dcopy(DIM * vec_DIM, Position_vec, 1, Pos, 1);
	}
	void set_Position_vec(const T* Pos) {		
		cblas_dcopy(DIM * vec_DIM, Pos, 1, Position_vec, 1);
	}
	void get_Position_vec_CG(T* Pos_CG) {
		cblas_dcopy(DIM, Position_vec, 1, Pos_CG, 1);
	}

	// Sums up all the 9 masses
	inline double get_global_mass() {
		return (Mass_vec[0] + // CG
			Mass_vec[1] + Mass_vec[2] + Mass_vec[3] + Mass_vec[4] +
			Mass_vec[5] + Mass_vec[6] + Mass_vec[7] + Mass_vec[8]);
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

	void get_Velocity_vec_xy(T* Vel) {
		// b=a, cblas_dcopy(n,a,inc,b,inc)
		cblas_dcopy((DIM - 1) * vec_DIM, Velocity_vec_xy, 1, Vel, 1);
	}

	void get_Position_vec_xy(T* Vel) {
		// b=a, cblas_dcopy(n,a,inc,b,inc)
		cblas_dcopy((DIM - 1) * vec_DIM, Position_vec_xy, 1, Vel, 1);
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
			cblas_dcopy( 2 * num_wheels, k_vec, 1, k, 1);
		}
	}
	void set_k_vec(const T* k) {
		if (k != NULL) {
			cblas_dcopy(2* num_wheels, k, 1, k_vec, 1);
		}
	}

	void get_Mass_vec(T* M) {
		if (M != NULL) {
			// b=a, cblas_dcopy(n,a,inc,b,inc)
			cblas_dcopy(vec_DIM, Mass_vec, 1, M, 1);
		}
	}
	double get_Mass_vec_CG() const {
		return Mass_vec[0];
	}
	void set_Mass_vec(const T* M) {
		if (M != NULL) {
			cblas_dcopy(vec_DIM, M, 1, Mass_vec, 1);
		}
	}
	/*
	 get distance vector from each important Point of the car (9: CG, 4*W_i, 4*T_i)
	 \param Point_P, 
	 \return each entry from Position_vec
	*/

	void get_dist_vector_xy(T* Point_P, T* dist_vector) {
		for (auto i = 0; i < vec_DIM; ++i) {
			cblas_dcopy(DIM-1, Point_P, incx, &dist_vector[(DIM-1) * i], incx);
		}
		// y=a-b, vdSub(n,a,b,y)
		vdSub((DIM-1) * vec_DIM, Position_vec_xy, dist_vector, dist_vector);
	}

	void get_dist_vector(T* Point_P, T* dist_vector) {	
		for (auto i = 0; i < vec_DIM; ++i) {
			cblas_dcopy(DIM, Point_P, incx, &dist_vector[DIM * i], incx);
		}
		// y=a-b, vdSub(n,a,b,y)
		vdSub(DIM * vec_DIM, Position_vec, dist_vector, dist_vector);
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
	void do_ALE_update(T* change, T* global_vect, size_t dim, size_t incx) {
#pragma loop(ivdep)
		for (size_t i = 0; i < dim; ++i) {
			global_vect[i*incx] += change[0];
			global_vect[i*incx + 1] += change[1];
		}
	}
	void get_ALE_change(T* current_ALE_vect, T* global_vect, T* change_vect) {
		change_vect[0] = current_ALE_vect[0] - global_vect[0];
		change_vect[1] = current_ALE_vect[1] - global_vect[1];
	}
	void get_vel_pos_change(T* velocity_change, T* position_change, T* angle_change) {
		get_ALE_change(Position_vec_xy, Position_vec_prev_xy, position_change);
		get_ALE_change(Velocity_vec_xy, Velocity_vec_prev_xy, velocity_change);
	}
	void apply_ALE_change() {
		/*Now both vector are at current state. swap pointer and CG location in new previous will be updated and following will be obselete which */

		T c, s;
		c = std::cos(*Angle_z);
		s = std::sin(*Angle_z);
		Position_vec_xy[2] = Position_vec_xy[0] + l_long[0] * c - l_lat[0] * s; // fl
		Position_vec_xy[3] = Position_vec_xy[1] + l_lat[0] * c + l_long[0] * s; // fl
		Position_vec_xy[4] = Position_vec_xy[2]; // fl
		Position_vec_xy[5] = Position_vec_xy[3]; // fl
		Position_vec_xy[6] = Position_vec_xy[0] + l_long[1] * c + l_lat[1] * s; // fr
		Position_vec_xy[7] = Position_vec_xy[1] - l_lat[1] * c + l_long[1] * s; // fr
		Position_vec_xy[8] = Position_vec_xy[6]; // fl
		Position_vec_xy[9] = Position_vec_xy[7]; // fl
		Position_vec_xy[10] = Position_vec_xy[0] - l_long[2] * c - l_lat[2] * s; // rl
		Position_vec_xy[11] = Position_vec_xy[1] + l_lat[2] * c - l_long[2] * s; // rl
		Position_vec_xy[12] = Position_vec_xy[10]; // fl
		Position_vec_xy[13] = Position_vec_xy[11]; // fl
		Position_vec_xy[14] = Position_vec_xy[0] - l_long[3] * c + l_lat[3] * s; // rr
		Position_vec_xy[15] = Position_vec_xy[1] - l_lat[3] * c - l_long[3] * s; // rr
		Position_vec_xy[16] = Position_vec_xy[14]; // fl
		Position_vec_xy[17] = Position_vec_xy[15]; // fl
	}

	void get_final_vel_pos_change(T* velocity_change, T* position_change) {
		get_ALE_change(Position_vec_xy, Position_vec, position_change);
		get_ALE_change(Velocity_vec_xy, Velocity_vec, velocity_change);
	}
	void apply_final_ALE_change(T* velocity_change, T* position_change) {
		do_ALE_update(position_change, Position_vec);
		do_ALE_update(velocity_change, Velocity_vec);
	}
	
	~Car() {
		mkl_free_buffers();
		mkl_free(Position_vec); 
		Position_vec = nullptr;
		mkl_free(Velocity_vec);
		Velocity_vec = nullptr;
		mkl_free(Mass_vec);
		Mass_vec = nullptr;
		mkl_free(angle_CG);
		angle_CG = nullptr;
		mkl_free(w_CG);
		w_CG = nullptr;
		mkl_free(I_CG);
		I_CG = nullptr;


		// Initial Conditions of the car
		mkl_free(initial_position);
		initial_position = nullptr;
		mkl_free(initial_velocity_vec);
		initial_velocity_vec = nullptr;
		mkl_free(initial_angle);
		initial_angle = nullptr;
		mkl_free(initial_angular_velocity);
		initial_angular_velocity = nullptr;

		
		
		mkl_free(quad_angle_init);
		quad_angle_init = nullptr;
		mkl_free(M_linear);
		M_linear = nullptr;
		mkl_free(temp_linear); 
		temp_linear = nullptr;
		mkl_free(K);
		K = nullptr;
		mkl_free(K_trans);
		K_trans = nullptr;
		mkl_free(D);
		D = nullptr;
		mkl_free(spring_length); 
		spring_length = nullptr;
		mkl_free(current_spring_length);
		current_spring_length = nullptr;
		mkl_free(u_prev_linear); 
		u_prev_linear = nullptr;
		mkl_free(u_current_linear);
		u_current_linear = nullptr;
		mkl_free(k_vec); 
		k_vec = nullptr;
		mkl_free(l_lat); 
		l_lat = nullptr;
		mkl_free(l_long);
		l_long = nullptr;
		mkl_free(velocity_current_linear);
		velocity_current_linear = nullptr;

		////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////// ALE Vectors ///////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////////////////
		mkl_free(Position_vec_xy);

		Position_vec_xy = nullptr;
		mkl_free(Velocity_vec_xy);

		Velocity_vec_xy = nullptr;
		delete Angle_z;
		Angle_z = nullptr;
		delete w_z;
		w_z = nullptr;
		///////////////////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////// Interpolator objects ///////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////
		mkl_free(Corners_current);
		Corners_current = nullptr;
		mkl_free(Corners_rot);
		Corners_rot = nullptr;
		mkl_free(Corners_init);
		Corners_init = nullptr;
		mkl_free(angle_buffer);
		angle_buffer = nullptr;
		mkl_free(pos_buffer);
		pos_buffer = nullptr;
	}
};
