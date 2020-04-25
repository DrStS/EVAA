#pragma once
#include <cmath>
#include <mkl.h>
#include "BLAS.h"

#include "MetaDataBase.h"
#include "MathLibrary.h"
#include "Constants.h"

template <typename T>
class Car {
private:
	/*
	Create the diagonal mass matrix M_linear to solve the 11DOF system
	\param Global_mass contains the masses of the 9 bodies (CG, 4 * W, 4 * T)
	\param Global_moment_Inertia contains the tensor of inertia of the body
	\param Mass_11DOF unused
	*/
	void Construct11DOFMass() {
		M_linear[0] = Mass_vec[0];
		M_linear[Constants::DOF + 1] = I_CG[0];
		M_linear[2 * Constants::DOF + 2] = I_CG[4];
		mkl<T>::copy(Constants::DOF - 3, Mass_vec + 1, 1, M_linear + 3 * Constants::DOF + 3, Constants::DOF + 1); // M_linear = diagonal matrix
	}

	/*
	Copy all X and Y coordinates of the global vector to the local vector
	\param Global_vector vector with coordinates [X,Y,Z,X,Y,Z,...]
	\param local_vector vector with coordinates [X,Y,X,Y,...]
	*/
	void ConstructALEVectors(T* globalVector, T* localVector) {
		for (size_t i = 0; i < Constants::VEC_DIM; ++i) {
			mkl<T>::copy(Constants::DIM - 1, globalVector + i * Constants::DIM, 1, localVector + i * (Constants::DIM - 1), 1);
		}
	}

	/*
	Construct corner initilizer
	*/
	void ConstructCorner(T* posCG, T* corners) {
		T c = std::cos(angle_CG[2]);
		T s = std::sin(angle_CG[2]);

		for (auto i = 0; i < Constants::NUM_LEGS; ++i) {
			mkl<T>::copy(Constants::DIM, posCG, Constants::INCX, corners + i, Constants::NUM_LEGS);
		}
		corners[0] += + l_long[0] * c - l_lat[0] * s; // fl
		corners[4] += + l_lat[0] * c + l_long[0] * s; // fl
	
		corners[1] += + l_long[1] * c + l_lat[1] * s; // fr
		corners[5] += - l_lat[1] * c + l_long[1] * s; // fr
	
		corners[2] += - l_long[2] * c - l_lat[2] * s; // rl
		corners[6] += + l_lat[2] * c - l_long[2] * s; // rl
	
		corners[3] += - l_long[3] * c + l_lat[3] * s; // rr
		corners[7] += - l_lat[3] * c - l_long[3] * s; // rr
	}

	/*
	Calculates the values of Corners for general angles.
	TODO: Consider moving to MathLibrary or anonymous namespace.
	*/
	void UpdateCorners11DOF(T* angles, T* rotation_mat_buffer, T* initial_corners, T* updated_corners)
	{
		// zz, yy, xx
		MathLibrary::get_rotation_matrix(angles[2], angles[1], angles[0], rotation_mat_buffer);

		// do rotation: rotationMat * r
		//void cblas_dgemm(const CBLAS_LAYOUT Layout, const CBLAS_TRANSPOSE transa, const CBLAS_TRANSPOSE transb, const MKL_INT m, const MKL_INT n, const MKL_INT k, const double alpha, const double* a, const MKL_INT lda, const double* b, const MKL_INT ldb, const double beta, double* c, const MKL_INT ldc);
		mkl<T>::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, Constants::DIM, Constants::NUM_LEGS, Constants::DIM, 1, rotation_mat_buffer, Constants::DIM, initial_corners, Constants::NUM_LEGS, 0, updated_corners, Constants::NUM_LEGS);
	}

	/*
	Calculates the values of Corners_current according to the current orientation
	*/
	void UpdateCorners11DOF() {
		pos_buffer[0] = Position_vec_xy[0];
		pos_buffer[1] = Position_vec_xy[1];
		pos_buffer[2] = Position_vec[2] + u_current_linear[0];
		angle_buffer[0] = angle_CG[1] + u_current_linear[1];
		angle_buffer[1] = angle_CG[2] + u_current_linear[2];
		angle_buffer[2] = 0;
		
		ConstructCorner(pos_buffer, Corners_init);
		UpdateCorners11DOF(angle_buffer, Corners_rot, Corners_init, Corners_current);
	}

	/*
	Was ist das? Maybe to MathLibrary
	*/
	void ConvertALEToGlobal(T* vect, T* global_vect) {
#pragma loop(ivdep)
		for (size_t i = 0; i < Constants::VEC_DIM; ++i) {
			global_vect[Constants::DIM * i] = vect[(Constants::DIM - 1) * i];
			global_vect[Constants::DIM * i + 1] = vect[(Constants::DIM - 1) * i + 1];
		}
	}

	/*
	Was ist das? Maybe to MathLibrary
	*/
	void Convert11DOFToGlobal(T* vect, T* global_vect) {
		global_vect[2] += vect[0];
#pragma loop(ivdep)
		for (size_t i = 1; i < Constants::VEC_DIM; ++i) {
			global_vect[Constants::DIM*i + 2] += vect[(Constants::DIM - 1) + i];
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

	T* Position_vec; // [CG, W1, T1, W2, T2, W3, T3, W4, T4] 9 x 3 !!! Consider Constants::ALIGNMENT (3+1),(3+1),... 
	T* Velocity_vec; // [CG, W1, T1, W2, T2, W3, T3, W4, T4] 9 x 3
	T* Mass_vec; // [CG,  W1, T1, W2, T2, W3, T3, W4, T4]    9 x 1
	T* global_mass;
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
	int DOF; // comes from xml (Consider extracting as constant, move to Constants.h)
	EVAALookup<Constants::floatEVAA>* lookupStiffness;
	T k_body_fl;
	T k_tyre_fl;
	T k_body_fr;
	T k_tyre_fr;
	T k_body_rl;
	T k_tyre_rl;
	T k_body_rr;
	T k_tyre_rr;
	
	T *quad_angle_init;
	T* M_linear, *temp_linear, *K, *K_trans, *D;
	T *spring_length, *current_spring_length;
	T* pos_nickpol, *distance_nickpol;
	T* u_prev_linear, *u_current_linear;
	T* k_vec, *l_lat, *l_long;
	T* velocity_current_linear;
	size_t* tyre_index_set;
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
	Car(EVAALookup<Constants::floatEVAA> * interpolator) {
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////// Generte Lookup Table /////////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		lookupStiffness = interpolator;
		DOF = MetaDataBase::DataBase()->getDOF();

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////// System Mono ///////////////////////////////////////////////////////
		///////////////////////////////// Memory Allocation and matrix formulation ///////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		const int positionAllocSize = Constants::DIM * Constants::VEC_DIM * sizeof(T); // 27 dimensions
		const int pointAllocSize = Constants::VEC_DIM * sizeof(T); // 9 dimensions
		const int dimAllocSize = Constants::DIM * sizeof(T); // 3 dimensions

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////// Params for Global coordinate ////////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		Position_vec = (T*)mkl_malloc(positionAllocSize, Constants::ALIGNMENT); 
		Velocity_vec = (T*)mkl_malloc(positionAllocSize, Constants::ALIGNMENT);
		Mass_vec = (T*)mkl_malloc(pointAllocSize, Constants::ALIGNMENT);
		angle_CG = (T*)mkl_malloc(dimAllocSize, Constants::ALIGNMENT);
		w_CG = (T*)mkl_malloc(dimAllocSize, Constants::ALIGNMENT);
		I_CG = (T*)mkl_malloc(Constants::DIM * Constants::DIM * sizeof(T), Constants::ALIGNMENT); // 9 Dimensions


		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////// Initial Params for Global coordinate ////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		initial_position = (T*)mkl_malloc(positionAllocSize, Constants::ALIGNMENT);
		initial_velocity_vec = (T*)mkl_malloc(positionAllocSize, Constants::ALIGNMENT);
		quad_angle_init = (T*)mkl_calloc(4, sizeof(T), Constants::ALIGNMENT); // 4 Constants::DIM
		initial_angle = (T*)mkl_malloc(dimAllocSize, Constants::ALIGNMENT); // 3 Constants::DIM
		initial_angular_velocity = (T*)mkl_malloc(dimAllocSize, Constants::ALIGNMENT); // 3 Constants::DIM


		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////// Params for 11 DOF system ////////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		M_linear = (T*)mkl_calloc(Constants::DOF * Constants::DOF, sizeof(T), Constants::ALIGNMENT); // 121 Constants::DIM
		temp_linear = (T*)mkl_malloc(Constants::DOF * sizeof(T), Constants::ALIGNMENT);  // 11 Constants::DIM
		K = (T*)mkl_calloc(Constants::DOF * Constants::DOF, sizeof(T), Constants::ALIGNMENT); // 121 Constants::DIM
		K_trans = (T*)mkl_malloc(Constants::DOF * Constants::DOF * sizeof(T), Constants::ALIGNMENT); // 121 Constants::DIM
		D = (T*)mkl_calloc(Constants::DOF * Constants::DOF, sizeof(T), Constants::ALIGNMENT); // 121 Constants::DIM
		u_prev_linear = (T*)mkl_malloc(Constants::DOF * sizeof(T), Constants::ALIGNMENT); // 11 Constants::DIM
		u_current_linear = (T*)mkl_malloc(Constants::DOF * sizeof(T), Constants::ALIGNMENT); // 11 Constants::DIM
		velocity_current_linear = (T*)mkl_malloc(Constants::DOF * sizeof(T), Constants::ALIGNMENT); // 3 Constants::DIM
		k_vec = (T*)mkl_malloc(Constants::NUM_LEGS * 2 * sizeof(T), Constants::ALIGNMENT); // 4 Constants::DIM
		l_lat = (T*)mkl_malloc(Constants::NUM_LEGS * sizeof(T), Constants::ALIGNMENT); // 4 Constants::DIM
		l_long = (T*)mkl_malloc(Constants::NUM_LEGS * sizeof(T), Constants::ALIGNMENT); // 4 Constants::DIM
		spring_length = (T*)mkl_malloc(2 * Constants::NUM_LEGS * sizeof(T), Constants::ALIGNMENT); // 8 Constants::DIM
		pos_nickpol = (T*)mkl_malloc(Constants::DIM * sizeof(T), Constants::ALIGNMENT); // 3 Constants::DIM
		distance_nickpol = (T*)mkl_malloc((2*Constants::NUM_LEGS + 1) * sizeof(T), Constants::ALIGNMENT); // 9 
		current_spring_length = (T*)mkl_malloc(2 * Constants::NUM_LEGS * sizeof(T), Constants::ALIGNMENT); // 8 Constants::DIM
		tyre_index_set = (size_t*)mkl_malloc(Constants::NUM_LEGS * sizeof(size_t), Constants::ALIGNMENT);

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////// Memory allocation for interpolator /////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		Corners_current = (T*)mkl_malloc(Constants::NUM_LEGS * dimAllocSize, Constants::ALIGNMENT); // 12 Constants::DIM
		Corners_rot = (T*)mkl_malloc(Constants::DIM * dimAllocSize, Constants::ALIGNMENT); // 9 Constants::DIM
		Corners_init = (T*)mkl_malloc(Constants::NUM_LEGS * dimAllocSize, Constants::ALIGNMENT); // 12 Constants::DIM
		angle_buffer = (T*)mkl_malloc(dimAllocSize, Constants::ALIGNMENT); // 3 Constants::DIM
		pos_buffer = (T*)mkl_malloc(dimAllocSize, Constants::ALIGNMENT); // 3 Constants::DIM


		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////// ALE Buffer Allocation /////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		Position_vec_xy = (T*)mkl_malloc((Constants::DIM - 1) * pointAllocSize, Constants::ALIGNMENT);
		Angle_z = new T;
		Velocity_vec_xy = (T*)mkl_malloc((Constants::DIM - 1) * pointAllocSize, Constants::ALIGNMENT);
		w_z = new T;
		global_mass = new T;


		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////// Extract Data from parser /////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		mkl<T>::copy(Constants::NUM_LEGS, MetaDataBase::DataBase()->getLongitudalLegPositionVector(), 1, l_long, 1);

		mkl<T>::copy(Constants::NUM_LEGS, MetaDataBase::DataBase()->getLatidudalLegPositionVector(), 1, l_lat, 1);

		mkl<T>::copy(Constants::DIM * Constants::DIM, MetaDataBase::DataBase()->getMomentOfInertiaVector(), 1, I_CG, 1);

		Mass_vec[0] =  MetaDataBase::DataBase()->getBodyMass();
		Mass_vec[1] =  MetaDataBase::DataBase()->getWheelMassFrontLeft();
		Mass_vec[2] =  MetaDataBase::DataBase()->getTyreMassFrontLeft();
		Mass_vec[3] = MetaDataBase::DataBase()->getWheelMassFrontRight();
		Mass_vec[4] = MetaDataBase::DataBase()->getTyreMassFrontRight();
		Mass_vec[5] = MetaDataBase::DataBase()->getWheelMassRearLeft();
		Mass_vec[6] = MetaDataBase::DataBase()->getTyreMassRearLeft();
		Mass_vec[7] = MetaDataBase::DataBase()->getWheelMassRearRight();
		Mass_vec[8] = MetaDataBase::DataBase()->getTyreMassRearRight();
		*global_mass = get_global_mass();

		spring_length[0] = MetaDataBase::DataBase()->getBodySpringLengthFrontLeft();
		spring_length[1] = MetaDataBase::DataBase()->getTyreSpringLengthFrontLeft();
		spring_length[2] = MetaDataBase::DataBase()->getBodySpringLengthFrontRight();
		spring_length[3] = MetaDataBase::DataBase()->getTyreSpringLengthFrontRight();
		spring_length[4] = MetaDataBase::DataBase()->getBodySpringLengthRearLeft();
		spring_length[5] = MetaDataBase::DataBase()->getTyreSpringLengthRearLeft();
		spring_length[6] = MetaDataBase::DataBase()->getBodySpringLengthRearRight();
		spring_length[7] = MetaDataBase::DataBase()->getTyreSpringLengthRearRight();

		mkl<T>::copy(Constants::DIM, MetaDataBase::DataBase()->getPositionCenterOfInstantaneousRotation(), 1, pos_nickpol, 1);


		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////// Initial Iteration vector ////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Initial Angles
		/*quad_angle_init[0] = params.initial_angle[0];
		quad_angle_init[1] = params.initial_angle[1];
		quad_angle_init[2] = params.initial_angle[2];
		quad_angle_init[3] = params.initial_angle[3];*/
		mkl<T>::copy(4, MetaDataBase::DataBase()->getBodyInitialOrientation(), 1, quad_angle_init, 1);
		MathLibrary::ToEulerAngles<T>(quad_angle_init, initial_angle);
		mkl<T>::copy(Constants::DIM, initial_angle, 1, angle_CG, 1);

		// Spring lengths
		current_spring_length[0] = MetaDataBase::DataBase()->getBodySpringInitialLengthFrontLeft();
		current_spring_length[1] = MetaDataBase::DataBase()->getTyreSpringInitialLengthFrontLeft();
		current_spring_length[2] = MetaDataBase::DataBase()->getBodySpringInitialLengthFrontRight();
		current_spring_length[3] = MetaDataBase::DataBase()->getTyreSpringInitialLengthFrontRight();
		current_spring_length[4] = MetaDataBase::DataBase()->getBodySpringInitialLengthRearLeft();
		current_spring_length[5] = MetaDataBase::DataBase()->getTyreSpringInitialLengthRearLeft();
		current_spring_length[6] = MetaDataBase::DataBase()->getBodySpringInitialLengthRearRight();
		current_spring_length[7] = MetaDataBase::DataBase()->getTyreSpringInitialLengthRearRight();


		// Filling the position vector with initial condition
		// CG
		mkl<T>::copy(Constants::DIM, MetaDataBase::DataBase()->getBodyInitialPosition(), 1, initial_position, 1); // copy the center of mass position

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////// Interpolator Initialization///////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// read init corners vectors into matrix

		ConstructCorner(initial_position, Corners_init); // only CG position is used to construct corners
		mkl<T>::copy(Constants::DIM, angle_CG, 1, angle_buffer, 1);
		angle_buffer[2] = 0;
		UpdateCorners11DOF(angle_buffer, Corners_rot, Corners_init, Corners_current);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////// Remaining position initialization ////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		const T* xml_start;
		T* position_start;
		if (MetaDataBase::DataBase()->getFlagInitialLeg()) {
			// if prescribed initial position (add a check for consistency with spring lengths)
			// W1 = W_fl
			xml_start = MetaDataBase::DataBase()->getWheelInitialPositionFrontLeft();
			position_start = initial_position + 3; //(end at 5)
			mkl<T>::copy(Constants::DIM, xml_start, 1, position_start, 1);
			// W2 = W_fr
			xml_start = MetaDataBase::DataBase()->getWheelInitialPositionFrontRight();
			position_start += 6; // skip 3 for tyre (end at 11)
			mkl<T>::copy(Constants::DIM, xml_start, 1, position_start, 1);
			// W3 = W_rl
			xml_start = MetaDataBase::DataBase()->getWheelInitialPositionRearLeft();
			position_start += 6; // skip 3 for tyre (end at 17)
			mkl<T>::copy(Constants::DIM, xml_start, 1, position_start, 1);
			// W2 = W_rr
			xml_start = MetaDataBase::DataBase()->getWheelInitialPositionRearRight();
			position_start += 6; // skip 3 for tyre (end at 23)
			mkl<T>::copy(Constants::DIM, xml_start, 1, position_start, 1);

			// T1 = T_fl
			xml_start = MetaDataBase::DataBase()->getTyreInitialPositionFrontLeft();
			position_start = initial_position + 6; // skip 3 for center of mass and 3 for the wheel
			mkl<T>::copy(Constants::DIM, xml_start, 1, position_start, 1); // (end at 8)
			// T2 = T_fr
			xml_start = MetaDataBase::DataBase()->getTyreInitialPositionFrontRight();
			position_start += 6; // skip 3 for the wheel
			mkl<T>::copy(Constants::DIM, xml_start, 1, position_start, 1); // (end at 14)
			// T3 = T_rl
			xml_start = MetaDataBase::DataBase()->getTyreInitialPositionRearLeft();
			position_start += 6; // skip 3 for the wheel
			mkl<T>::copy(Constants::DIM, xml_start, 1, position_start, 1);// (end at 20)
			// T4 = T_rr
			xml_start = MetaDataBase::DataBase()->getTyreInitialPositionRearRight();
			position_start += 6; // skip 3 for the wheel
			mkl<T>::copy(Constants::DIM, xml_start, 1, position_start, 1); // (end at 26)
			mkl<T>::copy(Constants::DIM * Constants::VEC_DIM, initial_position, 1, Position_vec, 1);
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
			/* Update the mean position where changes are to be added*/
			mkl<T>::copy(Constants::DIM, initial_position, 1, Position_vec, 1);
			W_fl = Position_vec + 3;
			W_fr = Position_vec + 9;
			W_rl = Position_vec + 15;
			W_rr = Position_vec + 21;
			T_fl = Position_vec + 6;
			T_fr = Position_vec + 12;
			T_rl = Position_vec + 18;
			T_rr = Position_vec + 24;
			get_length(Corners_current, spring_length, W_fl, T_fl, W_fr, T_fr, W_rl, T_rl, W_rr, T_rr);
		}
		
		// we updated the initial global position, so we 
		UpdateNickpolRadius();

		// Initial Velocity (Reuse the pointers)
		mkl<T>::copy(Constants::DIM, MetaDataBase::DataBase()->getBodyInitialVelocity(), 1, initial_velocity_vec, 1);
		// W1 = W_fl
		xml_start = MetaDataBase::DataBase()->getWheelInitialVelocityFrontLeft();
		position_start = initial_velocity_vec + 3;
		mkl<T>::copy(Constants::DIM, xml_start, 1, position_start, 1); // (end at 5)
		// W2 = W_fr
		xml_start = MetaDataBase::DataBase()->getWheelInitialVelocityFrontRight();
		position_start += 6; // skip 3 for tyre
		mkl<T>::copy(Constants::DIM, xml_start, 1, position_start, 1);// (end at 11)
		// W3 = W_rl
		xml_start = MetaDataBase::DataBase()->getWheelInitialVelocityRearLeft();
		position_start += 6; // skip 3 for tyre
		mkl<T>::copy(Constants::DIM, xml_start, 1, position_start, 1); // (end at 17)
		// W2 = W_rr
		xml_start = MetaDataBase::DataBase()->getWheelInitialVelocityRearRight();
		position_start += 6; // skip 3 for tyre
		mkl<T>::copy(Constants::DIM, xml_start, 1, position_start, 1); // (end at 23)

		// T1 = T_fl
		xml_start = MetaDataBase::DataBase()->getTyreInitialVelocityFrontLeft();
		position_start = initial_velocity_vec + 6; // skip 3 for center of mass and 3 for the wheel
		mkl<T>::copy(Constants::DIM, xml_start, 1, position_start, 1); // (end at 8)
		
		// T2 = T_fr
		xml_start = MetaDataBase::DataBase()->getTyreInitialVelocityFrontRight();
		position_start += 6; // skip 3 for the Tyre
		mkl<T>::copy(Constants::DIM, xml_start, 1, position_start, 1);// (end at 14)
		// T3 = T_rl
		xml_start = MetaDataBase::DataBase()->getTyreInitialVelocityRearLeft();
		position_start += 6; // skip 3 for the wheel
		mkl<T>::copy(Constants::DIM, xml_start, 1, position_start, 1); // (end at 20)
		
		// T4 = T_rr
		xml_start = MetaDataBase::DataBase()->getTyreInitialVelocityRearRight();
		position_start += 6; // skip 3 for the wheel
		mkl<T>::copy(Constants::DIM, xml_start, 1, position_start, 1); // (end at 26)

		// copy the initial position to the position vector
		mkl<T>::copy(Constants::DIM * Constants::VEC_DIM, initial_velocity_vec, 1, Velocity_vec, 1);

		// Initial Angular velocity
		mkl<T>::copy(Constants::DIM, MetaDataBase::DataBase()->getBodyInitialAngularVelocity(), 1, initial_angular_velocity, 1);
		mkl<T>::copy(Constants::DIM, initial_angular_velocity, 1, w_CG, 1);


		/*
		Global assignments Done */
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////// 11 DOF Buffer Initialization //////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		Construct11DOFMass();
		//compute_dx(u_prev_linear + 3);
		//construct_11DOF_vector(initial_position, initial_angle, u_prev_linear);
		compute_dx(u_prev_linear + 3);
		mkl<T>::copy(Constants::DOF, u_prev_linear, 1, u_current_linear, 1);
		tyre_index_set[0] = 2;
		tyre_index_set[1] = 4;
		tyre_index_set[2] = 6;
		tyre_index_set[3] = 8;

		construct_11DOF_vector(initial_velocity_vec, initial_angular_velocity, velocity_current_linear);
		/* This stays in the 11 DOF
		cblas_dscal(DOF, -h_, u_n_m_1, 1);
		cblas_daxpy(DOF, 1, u_n, 1, u_n_m_1, 1);
		*/
		if (MetaDataBase::DataBase()->getUseInterpolation()) {
			lookupStiffness->getInterpolation(current_spring_length, k_vec);
		}
		else {
			k_body_fl = MetaDataBase::DataBase()->getBodyStiffnessFrontLeft();
			k_tyre_fl = MetaDataBase::DataBase()->getTyreStiffnessFrontLeft();
			k_body_fr = MetaDataBase::DataBase()->getBodyStiffnessFrontRight();
			k_tyre_fr = MetaDataBase::DataBase()->getTyreStiffnessFrontRight();
			k_body_rl = MetaDataBase::DataBase()->getBodyStiffnessRearLeft();
			k_tyre_rl = MetaDataBase::DataBase()->getTyreStiffnessRearLeft();
			k_body_rr = MetaDataBase::DataBase()->getBodyStiffnessRearRight();
			k_tyre_rr = MetaDataBase::DataBase()->getTyreStiffnessRearRight();
			populate_K(k_vec, k_body_fl, k_tyre_fl, k_body_fr, k_tyre_fr, k_body_rl, k_tyre_rl, k_body_rr, k_tyre_rr);
		}
		update_K(k_vec);

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////// ALE Buffer Initialization /////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		ConstructALEVectors(Position_vec, Position_vec_xy);
		ConstructALEVectors(Velocity_vec, Velocity_vec_xy);
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
	void compute_normal_force(T* K, T* u, T* force, size_t* index, size_t n) {
		#pragma loop( ivdep )
		for (int i = 0; i < n; ++i) {
			force[index[i]] = -K[index[i] * Constants::DIM + index[i]] * u[index[i]];
		}
	}

	/*
	Calculate the entries of the stiffness matrix
	\param k_vect vector with all spring stiffnesses (in Stefans order)
	*/
	void update_K(T* k_vect) {
		// K = 0
		//cblas_dscal(DOF * DOF, 0.0, K, 1); // by copying upper to lower, we avoid reinitialization with 0

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
		K[Constants::DOF + 2] = -l_long[0] * l_lat[0] * k_vect[0] + l_lat[1] * l_long[1] * k_vect[2] + l_long[2] * l_lat[2] * k_vect[4] - l_long[3] * l_lat[3] * k_vect[6];
		K[Constants::DOF + 3] = -l_lat[0] * k_vect[0];
		K[Constants::DOF + 4] = 0;
		K[Constants::DOF + 5] = l_lat[1] * k_vect[2];
		K[Constants::DOF + 6] = 0;
		K[Constants::DOF + 7] = -l_lat[2] * k_vect[4];
		K[Constants::DOF + 8] = 0;
		K[Constants::DOF + 9] = l_lat[3] * k_vect[6];
		K[Constants::DOF + 10] = 0;

		temp_linear[2] = l_long[0] * l_long[0] * k_vect[0] + l_long[1] * l_long[1] * k_vect[2] + l_long[2] * l_long[2] * k_vect[4] + l_long[3] * l_long[3] * k_vect[6];
		K[2 * Constants::DOF + 3] = l_long[0] * k_vect[0];
		K[2 * Constants::DOF + 4] = 0;
		K[2 * Constants::DOF + 5] = l_long[1] * k_vect[2];
		K[2 * Constants::DOF + 6] = 0;
		K[2 * Constants::DOF + 7] = -l_long[2] * k_vect[4];
		K[2 * Constants::DOF + 8] = 0;
		K[2 * Constants::DOF + 9] = -l_long[3] * k_vect[6];
		K[2 * Constants::DOF + 10] = 0;

		temp_linear[3] = k_vect[0] + k_vect[1];
		K[3 * Constants::DOF + 4] = -k_vect[1];
		K[3 * Constants::DOF + 5] = 0;
		K[3 * Constants::DOF + 6] = 0;
		K[3 * Constants::DOF + 7] = 0;
		K[3 * Constants::DOF + 8] = 0;
		K[3 * Constants::DOF + 9] = 0;
		K[3 * Constants::DOF + 10] = 0;
		// all others are zero

		temp_linear[4] = k_vect[1];
		K[4 * Constants::DOF + 5] = 0;
		K[4 * Constants::DOF + 6] = 0;
		K[4 * Constants::DOF + 7] = 0;
		K[4 * Constants::DOF + 8] = 0;
		K[4 * Constants::DOF + 9] = 0;
		K[4 * Constants::DOF + 10] = 0;

		temp_linear[5] = k_vect[2] + k_vect[3];
		K[5 * Constants::DOF + 6] = -k_vect[3];
		K[5 * Constants::DOF + 7] = 0;
		K[5 * Constants::DOF + 8] = 0;
		K[5 * Constants::DOF + 9] = 0;
		K[5 * Constants::DOF + 10] = 0;

		temp_linear[6] = k_vect[3];
		K[6 * Constants::DOF + 7] = 0;
		K[6 * Constants::DOF + 8] = 0;
		K[6 * Constants::DOF + 9] = 0;
		K[6 * Constants::DOF + 10] = 0;

		temp_linear[7] = k_vect[4] + k_vect[5];
		K[7 * Constants::DOF + 8] = -k_vect[5];
		K[7 * Constants::DOF + 9] = 0;
		K[7 * Constants::DOF + 10] = 0;


		temp_linear[8] = k_vect[5];
		K[8 * Constants::DOF + 9] = 0;
		K[8 * Constants::DOF + 10] = 0;


		temp_linear[9] = k_vect[6] + k_vect[7];
		K[9 * Constants::DOF + 10] = -k_vect[7];

		temp_linear[10] = k_vect[7];
		

		// symmetrize K
		//cblas_dcopy(DOF * DOF, K, 1, K_trans, 1);
		mkl<T>::lacpy(LAPACK_ROW_MAJOR, 'U', Constants::DOF, Constants::DOF, K, Constants::DOF, K_trans, Constants::DOF);

		mkl<T>::imatcopy('R', 'T', Constants::DOF, Constants::DOF, 1.0, K_trans, Constants::DOF, Constants::DOF); // get transpose of matrix
		
		mkl<T>::lacpy(LAPACK_ROW_MAJOR, 'L', Constants::DOF, Constants::DOF, K_trans, Constants::DOF, K, Constants::DOF); // copy lower triangular in the orig matrix
		//cblas_daxpy(DOF * DOF, 1.0, K_trans, 1, K, 1); // K = K + K'

		// add the diagonal to K
		MathLibrary::allocate_to_diagonal(K, temp_linear, Constants::DOF); // K = K + K'+ diag(K)


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
		mkl<T>::copy(Constants::VEC_DIM - 1, Global_position + 5, 3, Position_11dof + 3, 1);
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
		mkl<T>::copy(Constants::DIM, Corners, 4, W_fl, 1);
		W_fl[2] -= curr_spring_len[0];
		mkl<T>::copy(Constants::DIM, W_fl, 1, T_fl, 1);
		T_fl[2] -= curr_spring_len[1];

		// W_fr & T_fr
		mkl<T>::copy(Constants::DIM, Corners + 1, 4, W_fr, 1);
		W_fr[2] -= curr_spring_len[2];
		mkl<T>::copy(Constants::DIM, W_fr, 1, T_fr, 1);
		T_fr[2] -= curr_spring_len[3];

		// W_rl & T_rl
		mkl<T>::copy(Constants::DIM, Corners + 2, 4, W_rl, 1);
		W_rl[2] -= curr_spring_len[4];
		mkl<T>::copy(Constants::DIM, W_rl, 1, T_rl, 1);
		T_rl[2] -= curr_spring_len[5];

		// W_rr & T_rr
		mkl<T>::copy(Constants::DIM, Corners + 3, 4, W_rr, 1);
		W_rr[2] -= curr_spring_len[6];
		mkl<T>::copy(Constants::DIM, W_rr, 1, T_rr, 1);
		T_rr[2] -= curr_spring_len[7];
	}


	inline void compute_dx(const T* current_length, T* dx) {
		/*
		the dx follows the order
		[w1, t1, w2, t2, w3, t3, w4, t4]
		current_length has input of type
		[w1, t1, w2, t2, w3, t3, w4, t4]

		*/

		mkl<T>::vSub(2 * (Constants::NUM_LEGS), spring_length, current_length, dx);
	}
	/*
	From the current elongations, calculate the differences to the rest positions
	\return length differences
	*/
	inline void compute_dx(T* dx) {
		compute_dx(current_spring_length, dx);
	}
	inline void compute_dx_tyre(T* dx) {
		mkl<T>::copy(Constants::NUM_LEGS, u_current_linear + 4, 2, dx, 1);
	}

	/*
	First updates the corner and afterwards compute the lengths of the springs
	*/
	void update_lengths_11DOF() {
		UpdateCorners11DOF();
		
		// after establishing the order of tyres & wheels, fix this.
		// With the current formulation, improvement:
		//	- current_spring_length[0,2,4,6] = Corners_current[8,9,10,11] + (u_current_linear[3,5,7,9] - Position_vec[1*Constants::DIM+2,3*Constants::DIM+2,5*Constants::DIM+2,7*Constants::DIM+2])
		//	- current_spring_length[1,3,5,7] = (Position_vec[1*Constants::DIM+2,3*Constants::DIM+2,5*Constants::DIM+2,7*Constants::DIM+2] - u_current_linear[3,5,7,9]) + (u_current_linear[4,6,8,10] - Position_vec[2*Constants::DIM+2,4*Constants::DIM+2,6*Constants::DIM+2,8*Constants::DIM+2])
		//  - current_spring_length = abs(current_spring_length)

		// current_spring_length[:] = u_current_linear[3:10]
		mkl<T>::copy(8, u_current_linear + 3, 1, current_spring_length, 1);

		// current_spring_length[:] -= Position_vec[..3..] -> in 3rd dimension =>
		//		=> current_spring_length[:] = u_current_linear[3:10] - Position_vec[5,8,11,...,26]
		mkl<T>::axpy(8, -1, Position_vec + (1 * Constants::DIM + 2), Constants::DIM, current_spring_length, 1);

		// subtract even indices from odd indices in current_spring_length: current_spring_length[1,3,5,7] -= current_spring_length[0,2,4,6]
		mkl<T>::axpy(4, -1, current_spring_length, 2, current_spring_length + 1, 2);
		
		// add Corners_current to even indices: current_spring_length[0,2,4,6] += Corners_current[8,9,10,11]
		mkl<T>::axpy(4, 1, Corners_current + 8, 1, current_spring_length, 2);

		// get absolute value for each element of the vector
		mkl<T>::vAbs(8, current_spring_length, current_spring_length);

		//std::cout << " Current spring length blas: \n\n";
		//for (auto i = 0; i < 8; ++i)
		//	std::cout << current_spring_length[i] << "\t ";

	/*	current_spring_length[0] = std::abs(Corners_current[8]									- (-u_current_linear[3] + Position_vec[1 * Constants::DIM + 2]));
		current_spring_length[1] = std::abs((-u_current_linear[3] + Position_vec[1 * Constants::DIM + 2])  - (-u_current_linear[4] + Position_vec[2 * Constants::DIM + 2]));
		current_spring_length[2] = std::abs(Corners_current[9]									- (-u_current_linear[5] + Position_vec[3 * Constants::DIM + 2]));
		current_spring_length[3] = std::abs((-u_current_linear[5] + Position_vec[3 * Constants::DIM + 2])  - (-u_current_linear[6] + Position_vec[4 * Constants::DIM + 2]));
		current_spring_length[4] = std::abs(Corners_current[10]									- (-u_current_linear[7] + Position_vec[5 * Constants::DIM + 2]));
		current_spring_length[5] = std::abs((-u_current_linear[7] + Position_vec[5 * Constants::DIM + 2])  - (-u_current_linear[8] + Position_vec[6 * Constants::DIM + 2]));
		current_spring_length[6] = std::abs(Corners_current[11]			 					    - (-u_current_linear[9] + Position_vec[7 * Constants::DIM + 2]));
		current_spring_length[7] = std::abs((-u_current_linear[9] + Position_vec[7 * Constants::DIM + 2])  - (-u_current_linear[10] + Position_vec[8 * Constants::DIM + 2]));*/

		//std::cout << " \n\n Current spring length classic: \n\n";
		//for (auto i = 0; i < 8; ++i)
		//	std::cout << current_spring_length[i] << "\t ";

		//exit(11);
		UpdateNickpolRadius();
	}

	void UpdateNickpolRadius() {
		T globalNickpolPosition = Position_vec[2] + pos_nickpol[2];		// get global Z position of the Nickpol
		for (int i = 0; i < 2 * Constants::NUM_LEGS + 1; ++i) {
			distance_nickpol[i] = Position_vec[2 + i * Constants::DIM] - globalNickpolPosition;
		}
	}

	/* Fills the global vector with all entries
	\param ALE_vectors contains X and Y components [GC:XY,W1:XY,T1:XY,W2:XY,T2:XY,...]
	\param vector 11DOF contains Z components [GC:Z,W1:Z,T1:Z,W2:Z,T2:Z,...]
	\return global_vector [GC:XYZ,W1:XYZ,T1:XYZ,W2:XYZ,T2:XYZ,...]
	*/
	void populate_results(T* ALE_vector, T * vector_11DOF, T* global_vector) {
		ConvertALEToGlobal(ALE_vector, global_vector);
		Convert11DOFToGlobal(vector_11DOF, global_vector);
	}
	

	void combine_results() {
		ConvertALEToGlobal(Position_vec_xy, Position_vec);
		Convert11DOFToGlobal(u_current_linear, Position_vec);

		// Angles manually
		angle_CG[0] = u_current_linear[1];
		angle_CG[1] = u_current_linear[2];
		angle_CG[2] = *Angle_z;
		w_CG[2] = *w_z;
		ConvertALEToGlobal(Velocity_vec_xy, Velocity_vec);
	}

	void get_Position_vec(T* Pos) {
		mkl<T>::copy(Constants::DIM * Constants::VEC_DIM, Position_vec, 1, Pos, 1);
	}
	void set_Position_vec(const T* Pos) {		
		mkl<T>::copy(Constants::DIM * Constants::VEC_DIM, Pos, 1, Position_vec, 1);
	}
	void get_Position_vec_CG(T* Pos_CG) {
		mkl<T>::copy(Constants::DIM, Position_vec, 1, Pos_CG, 1);
	}

	// Sums up all the 9 masses
	inline T get_global_mass() {
		return (Mass_vec[0] + // CG
			Mass_vec[1] + Mass_vec[2] + Mass_vec[3] + Mass_vec[4] +
			Mass_vec[5] + Mass_vec[6] + Mass_vec[7] + Mass_vec[8]);
	}

	void set_Position_vec_CG(const T* Pos_CG) {
		mkl<T>::copy(Constants::DIM, Pos_CG, 1, Position_vec, 1);
	}

	void get_Velocity_vec(T* Vel) {
		mkl<T>::copy(Constants::DIM * Constants::VEC_DIM, Velocity_vec, 1, Vel, 1);
	}

	void get_Velocity_vec_xy(T* Vel) {
		mkl<T>::copy((Constants::DIM - 1) * Constants::VEC_DIM, Velocity_vec_xy, 1, Vel, 1);
	}

	void get_Position_vec_xy(T* Vel) {
		mkl<T>::copy((Constants::DIM - 1) * Constants::VEC_DIM, Position_vec_xy, 1, Vel, 1);
	}

	void set_Velocity_vec(const T* Vel) {
		mkl<T>::copy(Constants::DIM * Constants::VEC_DIM, Vel, 1, Velocity_vec, 1);
	}
	void get_Velocity_vec_CG(T* Vel_CG) {
		mkl<T>::copy(Constants::DIM, Velocity_vec, 1, Vel_CG, 1);
	}
	void set_Velocity_vec_CG(const T* Vel_CG) {
		mkl<T>::copy(Constants::DIM, Vel_CG, 1, Velocity_vec, 1);
	}

	void get_k_vec(T* k) {
		mkl<T>::copy(2 * Constants::NUM_LEGS, k_vec, 1, k, 1);
	}
	void get_k_vec_tyre(T* k) {
		mkl<T>::copy(Constants::NUM_LEGS, k_vec + 1, 2, k, 1);
	}
	void get_k_vec_wheel(T* k) {
		mkl<T>::copy(Constants::NUM_LEGS, k_vec, 2, k, 1);
	}

	void set_k_vec(const T* k) {
		mkl<T>::copy(2 * Constants::NUM_LEGS, k, 1, k_vec, 1);
	}

	void get_Mass_vec(T* M) {
		mkl<T>::copy(Constants::VEC_DIM, Mass_vec, 1, M, 1);
	}
	T get_Mass_vec_CG() const {
		return Mass_vec[0];
	}
	void set_Mass_vec(const T* M) {
		mkl<T>::copy(Constants::VEC_DIM, M, 1, Mass_vec, 1);
	}

	/*
	 get distance vector from each important Point of the car (9: CG, 4*W_i, 4*T_i)
	 \param Point_P, 
	 \return each entry from Position_vec
	*/
	void get_dist_vector_xy(T* Point_P, T* dist_vector) {
		for (auto i = 0; i < Constants::VEC_DIM; ++i) {
			mkl<T>::copy(Constants::DIM-1, Point_P, Constants::INCX, &dist_vector[(Constants::DIM-1) * i], Constants::INCX);
		}
		mkl<T>::vSub((Constants::DIM-1) * Constants::VEC_DIM, Position_vec_xy, dist_vector, dist_vector);
	}

	void get_dist_vector(T* Point_P, T* dist_vector) {	
		for (auto i = 0; i < Constants::VEC_DIM; ++i) {
			mkl<T>::copy(Constants::DIM, Point_P, Constants::INCX, &dist_vector[Constants::DIM * i], Constants::INCX);
		}
		mkl<T>::vSub(Constants::DIM * Constants::VEC_DIM, Position_vec, dist_vector, dist_vector);
	} // 9 * 3 - from each important point to a fixed Point_P

	void get_dist_vector_CG(T* Point_P, T* dist_vector) {
		// get distance vector from Center of Gravity of the car to a Point P 
		// source: Point_P, dest: CG
		mkl<T>::vSub(Constants::DIM, Position_vec, Point_P, dist_vector);
	} // 3 - from Center of Gravity of a fixed Point_P

	T get_I_body_xx() const {
		return I_CG[0];
	}
	void set_I_body_xx(const T& I_body_xx_val) {
		I_CG[0] = I_body_xx_val;
	}
	T get_I_body_yy() const {
		return I_CG[4];
	}
	void set_I_body_yy(const T& I_body_yy_val) {
		I_CG[4] = I_body_yy_val;
	}
	void do_ALE_update(T* change, T* global_vect) {
#pragma loop(ivdep)
		for (size_t i = 0; i < Constants::DIM; ++i) {
			global_vect[i * Constants::INCX] += change[0];
			global_vect[i * Constants::INCX + 1] += change[1];
		}
	}
	void get_ALE_change(T* current_ALE_vect, T* global_vect, T* change_vect) {
		change_vect[0] = current_ALE_vect[0] - global_vect[0];
		change_vect[1] = current_ALE_vect[1] - global_vect[1];
	}

	void apply_ALE_change() {
		/*Now both vector are at current state. swap pointer and CG location in new previous will be updated and following will be obselete which */

		T c = std::cos(*Angle_z);
		T s = std::sin(*Angle_z);
		
		Position_vec_xy[2] = Position_vec_xy[0] + l_long[0] * c - l_lat[0] * s; // fl
		Position_vec_xy[3] = Position_vec_xy[1] + l_lat[0] * c + l_long[0] * s; // fl
		
		Position_vec_xy[6] = Position_vec_xy[0] + l_long[1] * c + l_lat[1] * s; // fr
		Position_vec_xy[7] = Position_vec_xy[1] - l_lat[1] * c + l_long[1] * s; // fr
		
		Position_vec_xy[10] = Position_vec_xy[0] - l_long[2] * c - l_lat[2] * s; // rl
		Position_vec_xy[11] = Position_vec_xy[1] + l_lat[2] * c - l_long[2] * s; // rl
		
		Position_vec_xy[14] = Position_vec_xy[0] - l_long[3] * c + l_lat[3] * s; // rr
		Position_vec_xy[15] = Position_vec_xy[1] - l_lat[3] * c - l_long[3] * s; // rr

		// copy values
		//Position_vec_xy[4] = Position_vec_xy[2]; // fl
		//Position_vec_xy[5] = Position_vec_xy[3]; // fl
		//
		//Position_vec_xy[8] = Position_vec_xy[6]; // fl
		//Position_vec_xy[9] = Position_vec_xy[7]; // fl
		//
		//Position_vec_xy[12] = Position_vec_xy[10]; // fl
		//Position_vec_xy[13] = Position_vec_xy[11]; // fl
		//
		//Position_vec_xy[16] = Position_vec_xy[14]; // fl
		//Position_vec_xy[17] = Position_vec_xy[15]; // fl

		//std::cout << "\n\n orig: " << Position_vec_xy[4] << "\t" << Position_vec_xy[5] << "\t"
		//	<< Position_vec_xy[8] << "\t" << Position_vec_xy[9] << "\t"
		//	<< Position_vec_xy[12] << "\t" << Position_vec_xy[13] << "\t"
		//	<< Position_vec_xy[16] << "\t" << Position_vec_xy[17] << "\t";

	    // Position_vec_xy[4,5,8,9,12,13,16,17] = Position_vec_xy[2,3,6,7,10,11,14,15]
		mkl<T>::copy(Constants::NUM_LEGS, Position_vec_xy + 2, 4, Position_vec_xy + 4, 4);
		mkl<T>::copy(Constants::NUM_LEGS, Position_vec_xy + 3, 4, Position_vec_xy + 5, 4);

	}

	void get_final_vel_pos_change(T* velocity_change, T* position_change) {
		get_ALE_change(Position_vec_xy, Position_vec, position_change);
		get_ALE_change(Velocity_vec_xy, Velocity_vec, velocity_change);
	}
	//void apply_final_ALE_change(T* velocity_change, T* position_change) {
	//	do_ALE_update(position_change, Position_vec);
	//	do_ALE_update(velocity_change, Velocity_vec);
	//}

	
	~Car() {
		mkl_free_buffers();
		
		/*
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
		mkl_free(pos_nickpol);
		pos_nickpol = nullptr;
		mkl_free(distance_nickpol);
		distance_nickpol = nullptr;
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
		mkl_free(tyre_index_set);
		tyre_index_set = nullptr;

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
		delete global_mass;
		global_mass = nullptr;

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
		*/
	}
};
