#pragma once

#include <chrono>
#include <mkl.h>
#include <string>
#include "BLAS.h"
#include "Constants.h"
#include "MathLibrary.h"
#include "MetaDataBase.h"

/*
Handles the whole MBD simulation (from the Matlab code)
*/
template <class T>
class MBD_method
{
private:
	////////////////////////////// Simulation Parameters ///////////////////////////////////////////////////////////////
	T h;
	size_t num_iter;
	int max_iter;
	T tol;
	size_t solution_dim; /// this is by the formulation
	std::string solver_name;
	int used_solver;
	
	////////////////////////////// Environment conditions //////////////////////////////////////////////////////////////////
	int boundary_conditions;
	T radius_circular_path;
	T* center_of_circle;

	////////////////////////////// Car Definition ///////////////////////////////////////////////////////////////////////
	bool use_interpolation;
	T k_body_fl;
	T k_tyre_fl;
	T k_body_fr;
	T k_tyre_fr;
	T k_body_rl;
	T k_tyre_rl;
	T k_body_rr;
	T k_tyre_rr;
	T k_body_rot_fl;
	T k_body_rot_fr;
	T k_body_rot_rl;
	T k_body_rot_rr;
	T k_tyre_rot_fl;
	T k_tyre_rot_fr;
	T k_tyre_rot_rl;
	T k_tyre_rot_rr;
	T c_body_fl;
	T c_tyre_fl;
	T c_body_fr;
	T c_tyre_fr;
	T c_body_rl;
	T c_tyre_rl;
	T c_body_rr;
	T c_tyre_rr;
	T l_long_fl;
	T l_long_fr;
	T l_long_rl;
	T l_long_rr;
	T l_lat_fl;
	T l_lat_fr;
	T l_lat_rl;
	T l_lat_rr;
	T mass;
	T I_body_xx;
	T I_body_yy;
	T I_body_zz;
	T mass_wheel_fl;
	T mass_tyre_fl;
	T mass_wheel_fr;
	T mass_tyre_fr;
	T mass_wheel_rl;
	T mass_tyre_rl;
	T mass_wheel_rr;
	T mass_tyre_rr;
	T upper_spring_length_rr;
	T upper_spring_length_rl;
	T upper_spring_length_fl;
	T upper_spring_length_fr;
	T lower_spring_length_rr;
	T lower_spring_length_rl;
	T lower_spring_length_fl;
	T lower_spring_length_fr;
	T* Ic;
	T* mass_wheel, * mass_tyre;
	T* upper_spring_length, * lower_spring_length;
	T* upper_spring_stiffness, * lower_spring_stiffness, * upper_spring_damping, * lower_spring_damping;
	T* upper_rotational_stiffness, * lower_rotational_stiffness;
	T* vc, * vw_fl, * vw_fr, * vw_rl, * vw_rr, * vt_fl, * vt_fr, * vt_rl, * vt_rr; // velocity of the center of mass, wheel and tyre.

	//////////////////////////	External Forces terms //////////////////////////////////////////////////////
	T g;
	T* FC;
	T FT_fl, FT_fr, FT_rl, FT_rr;
	T FW_fl, FW_fr, FW_rl, FW_rr;
	T FR_fl, FR_fr, FR_rl, FR_rr;

	//////////////////////////// Initial condition params //////////////////////////////////////////////////
	T* initial_upper_spring_length, * initial_lower_spring_length, * initial_orientation, * initial_angular_velocity;


	//////////////////////////// Arrays for intermediate steps /////////////////////////////////////////////
	T* r_fl, * r_fr, * r_rl, * r_rr;
	T* pcc; // position of center of mass
	T* x_vector;


	////////////////////////// Auxillary Parameters for Compute_f function //////////////////////////////////
	T* r_fl_tilda, * r_fr_tilda, * r_rl_tilda, * r_rr_tilda, ** FW, * FT, * A_Ic, * A_rem;


	/////////////////////////	Variables needed in compute_f function  /////////////////////////////////////
	T* cf_C_cN, * cf_r_up_fl, * cf_r_up_fr, * cf_r_up_rl, * cf_r_up_rr, * cf_r_low_fl, * cf_r_low_fr, * cf_r_low_rl, * cf_r_low_rr;
	T* cf_upper_normal_fl, * cf_upper_normal_fr, * cf_upper_normal_rl, * cf_upper_normal_rr, * cf_lower_normal_fl, * cf_lower_normal_fr, * cf_lower_normal_rl, * cf_lower_normal_rr, * cf_col_dat;
	T* cf_upper_force_fl, * cf_upper_force_fr, * cf_upper_force_rl, * cf_upper_force_rr, * cf_lower_force_fl, * cf_lower_force_fr, * cf_lower_force_rl, * cf_lower_force_rr;
	T* cf_upper_dampf_fl, * cf_upper_dampf_fr, * cf_upper_dampf_rl, * cf_upper_dampf_rr, * cf_lower_dampf_fl, * cf_lower_dampf_fr, * cf_lower_dampf_rl, * cf_lower_dampf_rr, * cf_temp;
	T* cf_upper_angle_fl, * cf_upper_angle_fr, * cf_upper_angle_rl, * cf_upper_angle_rr, * cf_lower_angle_fl, * cf_lower_angle_fr, * cf_lower_angle_rl, * cf_lower_angle_rr;
	T* cf_upper_S_fl, * cf_upper_S_fr, * cf_upper_S_rl, * cf_upper_S_rr, * cf_lower_S_fl, * cf_lower_S_fr, * cf_lower_S_rl, * cf_lower_S_rr;
	T* cf_lower_rot_force_fl, * cf_lower_rot_force_fr, * cf_lower_rot_force_rl, * cf_lower_rot_force_rr, * cf_upper_rot_force_fl, * cf_upper_rot_force_fr, * cf_upper_rot_force_rl, * cf_upper_rot_force_rr;
	T* cf_car_rot_force_fl, * cf_car_rot_force_fr, * cf_car_rot_force_rl, * cf_car_rot_force_rr, * cf_sum_car_force_fl, * cf_sum_car_force_fr, * cf_sum_car_force_rl, * cf_sum_car_force_rr;
	T* cf_local_FR_fl, * cf_local_FR_fr, * cf_local_FR_rl, * cf_local_FR_rr;
	T* cf_Hc, * cf_sum_torque_spring_car, * cf_Tc, * cf_wc_tilda;
	T* cf_b_rem, * cf_Qc, * cf_qc_dot;
	T* current_spring_lengths, * stiffness_vector;


	////////////////////////// Lookup table ///////////////////////////////////////////
	EVAALookup<Constants::floatEVAA>* lookupStiffness;

	/*
	Calculate the positions of the tyres and wheels according to the initial orientation of the car
	The legs always form a 90° angle to the car body, such that the rotational springs are at rest
	*/
	void get_initial_length(
		T* initial_orientation_,
		const T* r_fl_,
		const T* r_fr_,
		const T* r_rl_,
		const T* r_rr_,
		const T* pcc_,
		const T* initial_upper_spring_length_,
		const T* initial_lower_spring_length_,
		T* wheel_coordinate_fl_,
		T* wheel_coordinate_fr_,
		T* wheel_coordinate_rl_,
		T* wheel_coordinate_rr_,
		T* tyre_coordinate_fl_,
		T* tyre_coordinate_fr_,
		T* tyre_coordinate_rl_,
		T* tyre_coordinate_rr_)

	{
		/*
		To reduce memory trace and better use cache this function is implemented in following fashion:
		Original steps for computation of one component:
														1.	qc = qc/norm(qc);
														2.	C_Nc = get_basis(qc);
														3.	global_z = C_Nc(:,2);
														4.	global_z = -global_z / norm(global_z);
														5.	global_r_fl = pcc + C_Nc*r_fl;
														6.	upper_global_spring__fl = upper_length(_fl)*global_z;
														7.	lower_global_spring__fl = lower_length(_fl)*global_z;
														8.	pw_fl = global_r_fl + upper_global_spring__fl;
														9.	pt_fl = pw_fl + lower_global_spring__fl;
		Modified steps for computation of one component:
														_fl.	qc = qc/norm(qc);
														2.	C_Nc = get_basis(qc);
														3.	global_z = C_Nc(:,2);
														4.	global_z = -global_z / norm(global_z);
														5.	pw_fl = pcc;
														6.	pw_fl = pw_fl + C_Nc*r_fl;
														7.	pw_fl = pw_fl + upper_length(1)*global_z;
														8.	pt_fl = pw_fl
														8.	pt_fl	= pt_fl + lower_length(_fl)*global_z;
		*/

		T* global_z = (T*)mkl_calloc((Constants::DIM), sizeof(T), Constants::ALIGNMENT);
		T* C_Nc = (T*)mkl_calloc((Constants::DIM) * (Constants::DIM), sizeof(T), Constants::ALIGNMENT);

		//	1. qc = qc/norm(qc); This is in quaternions 
		T nrm = mkl<T>::nrm2(Constants::NUM_LEGS, initial_orientation_, 1);
		mkl<T>::scal(Constants::NUM_LEGS, 1.0 / nrm, initial_orientation_, 1);

		// 2.	C_Nc = get_basis(qc);
		MathLibrary::get_basis<T>(initial_orientation_, C_Nc);

		// 3.	global_z = C_Nc(:,2);
		mkl<T>::copy(Constants::DIM, C_Nc + 1, Constants::DIM, global_z, 1);

		// 4.	global_z = -global_z / norm(global_z);
		nrm = mkl<T>::nrm2(Constants::DIM, global_z, 1);
		mkl<T>::scal(Constants::DIM, -1.0 / nrm, global_z, 1);

		/////////////////////////////////////////// Leg 1 ////////////////////////////////////////////////////////
		// 5.	pw1 = pcc;
		mkl<T>::copy(Constants::DIM, pcc_, 1, wheel_coordinate_fl_, 1);

		// 6.	pw1 = pw1 + C_Nc*r1;
		mkl<T>::gemv(CblasRowMajor, CblasNoTrans, Constants::DIM, Constants::DIM, 1, C_Nc, Constants::DIM, r_fl_, 1, 1, wheel_coordinate_fl_, 1);

		// 7.	pw1 = pw1 + upper_length(1)*global_z;
		mkl<T>::axpy(Constants::DIM, initial_upper_spring_length_[0], global_z, 1, wheel_coordinate_fl_, 1);

		// 8.	pt1 = pw1
		mkl<T>::copy(Constants::DIM, wheel_coordinate_fl_, 1, tyre_coordinate_fl_, 1);

		// 9.	pt1 = pw1 + lower_length(1)*global_z;
		mkl<T>::axpy(Constants::DIM, initial_lower_spring_length_[0], global_z, 1, tyre_coordinate_fl_, 1);

		/////////////////////////////////////////// Leg 2 ////////////////////////////////////////////////////////
		// 5.	pw2 = pcc;
		mkl<T>::copy(Constants::DIM, pcc_, 1, wheel_coordinate_fr_, 1);

		// 6.	pw2 = pw2 + C_Nc*r2;
		mkl<T>::gemv(CblasRowMajor, CblasNoTrans, Constants::DIM, Constants::DIM, 1, C_Nc, Constants::DIM, r_fr_, 1, 1, wheel_coordinate_fr_, 1);

		// 7.	pw2 = pw2 + upper_length(2)*global_z;
		mkl<T>::axpy(Constants::DIM, initial_upper_spring_length_[1], global_z, 1, wheel_coordinate_fr_, 1);

		// 8.	pt2 = pw2
		mkl<T>::copy(Constants::DIM, wheel_coordinate_fr_, 1, tyre_coordinate_fr_, 1);

		// 9.	pt2 = pw2 + lower_length(2)*global_z;
		mkl<T>::axpy(Constants::DIM, initial_lower_spring_length_[1], global_z, 1, tyre_coordinate_fr_, 1);

		/////////////////////////////////////////// Leg 3 ////////////////////////////////////////////////////////
		// 5.	pw3 = pcc;
		mkl<T>::copy(Constants::DIM, pcc_, 1, wheel_coordinate_rl_, 1);

		// 6.	pw3 = pw3 + C_Nc*r3;
		mkl<T>::gemv(CblasRowMajor, CblasNoTrans, Constants::DIM, Constants::DIM, 1, C_Nc, Constants::DIM, r_rl_, 1, 1, wheel_coordinate_rl_, 1);

		// 7.	pw3 = pw3 + upper_length(3)*global_z;
		mkl<T>::axpy(Constants::DIM, initial_upper_spring_length_[2], global_z, 1, wheel_coordinate_rl_, 1);

		// 8.	pt3 = pw3
		mkl<T>::copy(Constants::DIM, wheel_coordinate_rl_, 1, tyre_coordinate_rl_, 1);

		// 9.	pt3 = pw3 + lower_length(3)*global_z;
		mkl<T>::axpy(Constants::DIM, initial_lower_spring_length_[2], global_z, 1, tyre_coordinate_rl_, 1);

		/////////////////////////////////////////// Leg 4 ////////////////////////////////////////////////////////
		// 5.	pw4 = pcc;
		mkl<T>::copy(Constants::DIM, pcc_, 1, wheel_coordinate_rr_, 1);

		// 6.	pw4 = pw4 + C_Nc*r4;
		mkl<T>::gemv(CblasRowMajor, CblasNoTrans, Constants::DIM, Constants::DIM, 1, C_Nc, Constants::DIM, r_rr_, 1, 1, wheel_coordinate_rr_, 1);

		// 7.	pw4 = pw4 + upper_length(4)*global_z;
		mkl<T>::axpy(Constants::DIM, initial_upper_spring_length_[3], global_z, 1, wheel_coordinate_rr_, 1);

		// 8.	pt4 = pw4
		mkl<T>::copy(Constants::DIM, wheel_coordinate_rr_, 1, tyre_coordinate_rr_, 1);

		// 9.	pt4 = pw4 + lower_length(4)*global_z;
		mkl<T>::axpy(Constants::DIM, initial_lower_spring_length_[3], global_z, 1, tyre_coordinate_rr_, 1);


		mkl_free(global_z);
		mkl_free(C_Nc);
	}

	/*
	For Debug purposes
	*/
	void write_matrix(T* vect, int count) {
		std::cout << "Debug mode print" << std::endl;
		for (size_t i = 0; i < count; ++i) {
			//std::cout << vect[i] << std::endl;
			std::cout.precision(5);
			for (size_t j = 0; j < count; ++j) {
				std::cout << std::scientific << vect[i * count + j] << "  ";
			}
			std::cout << "\n" << std::endl;
			//printf("%1.15f\n", vect[i]);
		}
		//exit(5);
	}

	/*
	For Debug purposes
	*/
	void write_vector(T* vect, int count) {
		std::cout << "Debug mode print" << std::endl;
		for (size_t i = 0; i < count; ++i) {
			//std::cout << vect[i] << std::endl;
			std::cout.precision(15);
			std::cout << std::scientific << vect[i] << std::endl;
			//printf("%1.15f\n", vect[i]);
		}
		//exit(5);
	}


public:
	/*
	Constructor
	*/
	MBD_method(EVAALookup<Constants::floatEVAA> * interpolator) {

		////////////////////////////// Simulation Parameters ///////////////////////////////////////////////////////////////
		h = MetaDataBase::DataBase()->getTimeStepSize();
		num_iter = MetaDataBase::DataBase()->getNumberOfTimeIterations();
		max_iter = MetaDataBase::DataBase()->getMaxNumberOfBroydenIterationForMBD();
		tol = MetaDataBase::DataBase()->getToleranceBroydenIterationForMBD();
		solution_dim = MetaDataBase::DataBase()->getSolutionVectorSize(); /// this is by the formulation
		used_solver = MetaDataBase::DataBase()->getUsedSolverForMBD();
		boundary_conditions = MetaDataBase::DataBase()->getRoadConditions();
		radius_circular_path = MetaDataBase::DataBase()->getCircularRoadRadius();
		use_interpolation = MetaDataBase::DataBase()->getUseInterpolation();
		lookupStiffness = interpolator;

		////////////////////////////// Car Definition ///////////////////////////////////////////////////////////////////////
		k_body_fl = MetaDataBase::DataBase()->getBodyStiffnessFrontLeft();
		k_tyre_fl = MetaDataBase::DataBase()->getTyreStiffnessFrontLeft();
		k_body_fr = MetaDataBase::DataBase()->getBodyStiffnessFrontRight();
		k_tyre_fr = MetaDataBase::DataBase()->getTyreStiffnessFrontRight();
		k_body_rl = MetaDataBase::DataBase()->getBodyStiffnessRearLeft();
		k_tyre_rl = MetaDataBase::DataBase()->getTyreStiffnessRearLeft();
		k_body_rr = MetaDataBase::DataBase()->getBodyStiffnessRearRight();
		k_tyre_rr = MetaDataBase::DataBase()->getTyreStiffnessRearRight();

		k_body_rot_fl = 1e5;
		k_body_rot_fr = 1e5;
		k_body_rot_rl = 1e5;
		k_body_rot_rr = 1e5;
		k_tyre_rot_fl = 1e5;
		k_tyre_rot_fr = 1e5;
		k_tyre_rot_rl = 1e5;
		k_tyre_rot_rr = 1e5;

		c_body_fl = MetaDataBase::DataBase()->getBodyDampingFrontLeft();
		c_tyre_fl = MetaDataBase::DataBase()->getTyreDampingFrontLeft();
		c_body_fr = MetaDataBase::DataBase()->getBodyDampingFrontRight();
		c_tyre_fr = MetaDataBase::DataBase()->getTyreDampingFrontRight();
		c_body_rl = MetaDataBase::DataBase()->getBodyDampingRearLeft();
		c_tyre_rl = MetaDataBase::DataBase()->getTyreDampingRearLeft();
		c_body_rr = MetaDataBase::DataBase()->getBodyDampingRearRight();
		c_tyre_rr = MetaDataBase::DataBase()->getTyreDampingRearRight();

		l_long_fl = MetaDataBase::DataBase()->getLongitudalLegPositionFrontLeft();
		l_long_fr = MetaDataBase::DataBase()->getLongitudalLegPositionFrontRight();
		l_long_rl = MetaDataBase::DataBase()->getLongitudalLegPositionRearLeft();
		l_long_rr = MetaDataBase::DataBase()->getLongitudalLegPositionRearRight();

		l_lat_fl = MetaDataBase::DataBase()->getLatidudalLegPositionFrontLeft();
		l_lat_fr = MetaDataBase::DataBase()->getLatidudalLegPositionFrontRight();
		l_lat_rl = MetaDataBase::DataBase()->getLatidudalLegPositionRearLeft();
		l_lat_rr = MetaDataBase::DataBase()->getLatidudalLegPositionRearRight();

		mass = MetaDataBase::DataBase()->getBodyMass();

		I_body_xx = MetaDataBase::DataBase()->getMomentOfInertiaXX();
		I_body_yy = MetaDataBase::DataBase()->getMomentOfInertiaYY();
		I_body_zz = MetaDataBase::DataBase()->getMomentOfInertiaZZ();

		mass_wheel_fl = MetaDataBase::DataBase()->getWheelMassFrontLeft();
		mass_tyre_fl = MetaDataBase::DataBase()->getTyreMassFrontLeft();
		mass_wheel_fr = MetaDataBase::DataBase()->getWheelMassFrontRight();
		mass_tyre_fr = MetaDataBase::DataBase()->getTyreMassFrontRight();
		mass_wheel_rl = MetaDataBase::DataBase()->getWheelMassRearLeft();
		mass_tyre_rl = MetaDataBase::DataBase()->getTyreMassRearLeft();
		mass_wheel_rr = MetaDataBase::DataBase()->getWheelMassRearRight();
		mass_tyre_rr = MetaDataBase::DataBase()->getTyreMassRearRight();

		upper_spring_length_fl = MetaDataBase::DataBase()->getBodySpringLengthFrontLeft();
		upper_spring_length_fr = MetaDataBase::DataBase()->getBodySpringLengthFrontRight();
		upper_spring_length_rl = MetaDataBase::DataBase()->getBodySpringLengthRearLeft();
		upper_spring_length_rr = MetaDataBase::DataBase()->getBodySpringLengthRearRight();

		lower_spring_length_fl = MetaDataBase::DataBase()->getTyreSpringLengthFrontLeft();
		lower_spring_length_fr = MetaDataBase::DataBase()->getTyreSpringLengthFrontRight();
		lower_spring_length_rl = MetaDataBase::DataBase()->getTyreSpringLengthRearLeft();
		lower_spring_length_rr = MetaDataBase::DataBase()->getTyreSpringLengthRearRight();

		g = MetaDataBase::DataBase()->getGravityField()[2];

		int i;

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////// Memory Allocation and matrix formulation ///////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		r_fl = (T*)mkl_calloc(Constants::DIM, sizeof(T), Constants::ALIGNMENT);
		r_fr = (T*)mkl_calloc(Constants::DIM, sizeof(T), Constants::ALIGNMENT);
		r_rl = (T*)mkl_calloc(Constants::DIM, sizeof(T), Constants::ALIGNMENT);
		r_rr = (T*)mkl_calloc(Constants::DIM, sizeof(T), Constants::ALIGNMENT);
		Ic = (T*)mkl_calloc((Constants::DIM) * (Constants::DIM), sizeof(T), Constants::ALIGNMENT);
		mass_wheel = (T*)mkl_calloc(Constants::NUM_LEGS, sizeof(T), Constants::ALIGNMENT);
		mass_tyre = (T*)mkl_calloc(Constants::NUM_LEGS, sizeof(T), Constants::ALIGNMENT);
		upper_spring_length = (T*)mkl_calloc(Constants::NUM_LEGS, sizeof(T), Constants::ALIGNMENT);
		lower_spring_length = (T*)mkl_calloc(Constants::NUM_LEGS, sizeof(T), Constants::ALIGNMENT);
		upper_spring_stiffness = (T*)mkl_calloc(Constants::NUM_LEGS, sizeof(T), Constants::ALIGNMENT);
		lower_spring_stiffness = (T*)mkl_calloc(Constants::NUM_LEGS, sizeof(T), Constants::ALIGNMENT);
		upper_rotational_stiffness = (T*)mkl_calloc(Constants::NUM_LEGS, sizeof(T), Constants::ALIGNMENT);
		lower_rotational_stiffness = (T*)mkl_calloc(Constants::NUM_LEGS, sizeof(T), Constants::ALIGNMENT);
		initial_upper_spring_length = (T*)mkl_calloc(Constants::NUM_LEGS, sizeof(T), Constants::ALIGNMENT);
		initial_lower_spring_length = (T*)mkl_calloc(Constants::NUM_LEGS, sizeof(T), Constants::ALIGNMENT);
		initial_orientation = (T*)mkl_calloc(Constants::NUM_LEGS, sizeof(T), Constants::ALIGNMENT);
		initial_angular_velocity = (T*)mkl_calloc(Constants::DIM, sizeof(T), Constants::ALIGNMENT);
		upper_spring_damping = (T*)mkl_calloc(Constants::NUM_LEGS, sizeof(T), Constants::ALIGNMENT);
		lower_spring_damping = (T*)mkl_calloc(Constants::NUM_LEGS, sizeof(T), Constants::ALIGNMENT);
		vc = (T*)mkl_calloc(Constants::DIM, sizeof(T), Constants::ALIGNMENT);
		vw_fl = (T*)mkl_calloc(Constants::DIM, sizeof(T), Constants::ALIGNMENT);
		vw_fr = (T*)mkl_calloc(Constants::DIM, sizeof(T), Constants::ALIGNMENT);
		vw_rl = (T*)mkl_calloc(Constants::DIM, sizeof(T), Constants::ALIGNMENT);
		vw_rr = (T*)mkl_calloc(Constants::DIM, sizeof(T), Constants::ALIGNMENT);
		vt_fl = (T*)mkl_calloc(Constants::DIM, sizeof(T), Constants::ALIGNMENT);
		vt_fr = (T*)mkl_calloc(Constants::DIM, sizeof(T), Constants::ALIGNMENT);
		vt_rl = (T*)mkl_calloc(Constants::DIM, sizeof(T), Constants::ALIGNMENT);
		vt_rr = (T*)mkl_calloc(Constants::DIM, sizeof(T), Constants::ALIGNMENT);
		pcc = (T*)mkl_calloc(Constants::DIM, sizeof(T), Constants::ALIGNMENT);
		FC = (T*)mkl_calloc(Constants::DIM, sizeof(T), Constants::ALIGNMENT);
		r_fl_tilda = (T*)mkl_calloc((Constants::DIM) * (Constants::DIM), sizeof(T), Constants::ALIGNMENT);
		r_fr_tilda = (T*)mkl_calloc((Constants::DIM) * (Constants::DIM), sizeof(T), Constants::ALIGNMENT);
		r_rl_tilda = (T*)mkl_calloc((Constants::DIM) * (Constants::DIM), sizeof(T), Constants::ALIGNMENT);
		r_rr_tilda = (T*)mkl_calloc((Constants::DIM) * (Constants::DIM), sizeof(T), Constants::ALIGNMENT);
		FW = (T**)mkl_calloc((Constants::NUM_LEGS), sizeof(T*), Constants::ALIGNMENT);
		FT = (T*)mkl_calloc((Constants::NUM_LEGS), sizeof(T*), Constants::ALIGNMENT);
		A_Ic = (T*)mkl_calloc((Constants::DIM) * (Constants::DIM), sizeof(T), Constants::ALIGNMENT);
		A_rem = (T*)mkl_calloc(9 * Constants::DIM, sizeof(T), Constants::ALIGNMENT);
		center_of_circle = (T*)mkl_calloc(Constants::DIM, sizeof(T), Constants::ALIGNMENT);

		i = 0;
		r_fl[i] = -l_long_fl; r_fr[i] = -l_long_fr; r_rl[i] = l_long_fl; r_rr[i] = l_long_rr;
		Ic[i * Constants::DIM + i] = I_body_xx;
		mass_wheel[i] = mass_wheel_fl;
		mass_tyre[i] = mass_tyre_fl;
		upper_spring_length[i] = upper_spring_length_fl;
		lower_spring_length[i] = lower_spring_length_fl;
		upper_spring_stiffness[i] = k_body_fl;
		lower_spring_stiffness[i] = k_tyre_fl;
		upper_rotational_stiffness[i] = k_body_rot_fl;
		lower_rotational_stiffness[i] = k_tyre_rot_fl;
		FT_fl = -mass_tyre[i] * g;
		FW_fl = -mass_wheel[i] * g;
		upper_spring_damping[i] = c_body_fl;
		lower_spring_damping[i] = c_tyre_fl;

		initial_upper_spring_length[i] = MetaDataBase::DataBase()->getBodySpringInitialLengthFrontLeft();
		initial_lower_spring_length[i] = MetaDataBase::DataBase()->getTyreSpringInitialLengthFrontLeft();
		initial_orientation[i] = MetaDataBase::DataBase()->getBodyInitialOrientation()[i];
		vc[i] = MetaDataBase::DataBase()->getBodyInitialVelocity()[i];
		vw_fl[i] = MetaDataBase::DataBase()->getWheelInitialVelocityFrontLeft()[i];
		vw_fr[i] = MetaDataBase::DataBase()->getWheelInitialVelocityFrontRight()[i];
		vw_rl[i] = MetaDataBase::DataBase()->getWheelInitialVelocityRearLeft()[i];
		vw_rr[i] = MetaDataBase::DataBase()->getWheelInitialVelocityRearRight()[i];
		vt_fl[i] = MetaDataBase::DataBase()->getTyreInitialVelocityFrontLeft()[i];
		vt_fr[i] = MetaDataBase::DataBase()->getTyreInitialVelocityFrontRight()[i];
		vt_rl[i] = MetaDataBase::DataBase()->getTyreInitialVelocityRearLeft()[i];
		vt_rr[i] = MetaDataBase::DataBase()->getTyreInitialVelocityRearRight()[i];
		pcc[i] = MetaDataBase::DataBase()->getBodyInitialPosition()[i];
		FC[i] = 0;
		center_of_circle[i] = MetaDataBase::DataBase()->getCircularRoadCenter()[i];

		i = 1;
		r_fl[i] = l_lat_fl; r_fr[i] = -l_lat_fr; r_rl[i] = -l_lat_fl; r_rr[i] = l_lat_rr;
		Ic[i * Constants::DIM + i] = I_body_yy;
		mass_wheel[i] = mass_wheel_fr;
		mass_tyre[i] = mass_tyre_fr;
		upper_spring_length[i] = upper_spring_length_fr;
		lower_spring_length[i] = lower_spring_length_fr;
		upper_spring_stiffness[i] = k_body_fr;
		lower_spring_stiffness[i] = k_tyre_fr;
		upper_rotational_stiffness[i] = k_body_rot_fr;
		lower_rotational_stiffness[i] = k_tyre_rot_fr;
		FT_fr = -mass_tyre[i] * g;
		FW_fr = -mass_wheel[i] * g;
		upper_spring_damping[i] = c_body_fr;
		lower_spring_damping[i] = c_tyre_fr;
		initial_upper_spring_length[i] = MetaDataBase::DataBase()->getBodySpringInitialLengthFrontRight();
		initial_lower_spring_length[i] = MetaDataBase::DataBase()->getTyreSpringInitialLengthFrontRight();
		initial_orientation[i] = MetaDataBase::DataBase()->getBodyInitialOrientation()[i];
		vc[i] = MetaDataBase::DataBase()->getBodyInitialVelocity()[i];
		vw_fl[i] = MetaDataBase::DataBase()->getWheelInitialVelocityFrontLeft()[i];
		vw_fr[i] = MetaDataBase::DataBase()->getWheelInitialVelocityFrontRight()[i];
		vw_rl[i] = MetaDataBase::DataBase()->getWheelInitialVelocityRearLeft()[i];
		vw_rr[i] = MetaDataBase::DataBase()->getWheelInitialVelocityRearRight()[i];
		vt_fl[i] = MetaDataBase::DataBase()->getTyreInitialVelocityFrontLeft()[i];
		vt_fr[i] = MetaDataBase::DataBase()->getTyreInitialVelocityFrontRight()[i];
		vt_rl[i] = MetaDataBase::DataBase()->getTyreInitialVelocityRearLeft()[i];
		vt_rr[i] = MetaDataBase::DataBase()->getTyreInitialVelocityRearRight()[i];
		pcc[i] = MetaDataBase::DataBase()->getBodyInitialPosition()[i];
		FC[i] = -mass*g;
		center_of_circle[i] = MetaDataBase::DataBase()->getCircularRoadCenter()[i];

		i = 2;
		r_fl[i] = 0; r_fr[i] = 0; r_rl[i] = 0; r_rr[i] = 0;
		Ic[i * Constants::DIM + i] = I_body_zz;
		mass_wheel[i] = mass_wheel_rl;
		mass_tyre[i] = mass_tyre_rl;
		upper_spring_length[i] = upper_spring_length_rl;
		lower_spring_length[i] = lower_spring_length_rl;
		upper_spring_stiffness[i] = k_body_rl;
		lower_spring_stiffness[i] = k_tyre_rl;
		upper_rotational_stiffness[i] = k_body_rot_rl;
		lower_rotational_stiffness[i] = k_tyre_rot_rl;
		FT_rl = -mass_tyre[i] * g;
		FW_rl = -mass_wheel[i] * g;
		upper_spring_damping[i] = c_body_rl;
		lower_spring_damping[i] = c_tyre_rl;
		initial_upper_spring_length[i] = MetaDataBase::DataBase()->getBodySpringInitialLengthRearLeft();
		initial_lower_spring_length[i] = MetaDataBase::DataBase()->getTyreSpringInitialLengthRearLeft();
		initial_orientation[i] = MetaDataBase::DataBase()->getBodyInitialOrientation()[i];
		vc[i] = MetaDataBase::DataBase()->getBodyInitialVelocity()[i];
		vw_fl[i] = MetaDataBase::DataBase()->getWheelInitialVelocityFrontLeft()[i];
		vw_fr[i] = MetaDataBase::DataBase()->getWheelInitialVelocityFrontRight()[i];
		vw_rl[i] = MetaDataBase::DataBase()->getWheelInitialVelocityRearLeft()[i];
		vw_rr[i] = MetaDataBase::DataBase()->getWheelInitialVelocityRearRight()[i];
		vt_fl[i] = MetaDataBase::DataBase()->getTyreInitialVelocityFrontLeft()[i];
		vt_fr[i] = MetaDataBase::DataBase()->getTyreInitialVelocityFrontRight()[i];
		vt_rl[i] = MetaDataBase::DataBase()->getTyreInitialVelocityRearLeft()[i];
		vt_rr[i] = MetaDataBase::DataBase()->getTyreInitialVelocityRearRight()[i];
		pcc[i] = MetaDataBase::DataBase()->getBodyInitialPosition()[i];
		FC[i] = 0;
		center_of_circle[i] = MetaDataBase::DataBase()->getCircularRoadCenter()[i];

		i = 3;
		mass_wheel[i] = mass_wheel_rr;
		mass_tyre[i] = mass_tyre_rr;
		upper_spring_length[i] = upper_spring_length_rr;
		lower_spring_length[i] = lower_spring_length_rr;
		upper_spring_stiffness[i] = k_body_rr;
		lower_spring_stiffness[i] = k_tyre_rr;
		upper_rotational_stiffness[i] = k_body_rot_rr;
		lower_rotational_stiffness[i] = k_tyre_rot_rr;
		FT_rr = -mass_tyre[i] * g;
		FW_rr = -mass_wheel[i] * g;
		upper_spring_damping[i] = c_body_rr;
		lower_spring_damping[i] = c_tyre_rr;
		initial_upper_spring_length[i] = MetaDataBase::DataBase()->getBodySpringInitialLengthRearRight();
		initial_lower_spring_length[i] = MetaDataBase::DataBase()->getTyreSpringInitialLengthRearRight();
		initial_orientation[i] = MetaDataBase::DataBase()->getBodyInitialOrientation()[i];

		// A_Ic has cholesky factorization of Ic
		mkl<T>::copy(Constants::DIM * Constants::DIM, Ic, 1, A_Ic, 1);
		mkl<T>::potrf(LAPACK_ROW_MAJOR, 'L', Constants::DIM, A_Ic, Constants::DIM);
	}

	/*
	Updates the velocity of the 4 wheels and tyres as well as the angular velocity, such that the car is already in the trajectory of the circle
	overwrites the velocity in the tyres and wheels and the angular velocities
	only keeps the tangential component of the velocity of the car body 
	*/
	void circular_path_initialization(T* vc, T* vw_fl, T* vw_fr, T* vw_rl, T* vw_rr, 
		T* vt_fl, T* vt_fr, T* vt_rl, T* vt_rr, T* omega, T* pcc, 
		T* pt_fl , T* pt_fr, T* pt_rl, T* pt_rr, T &radius_param) {

		const MKL_INT dim = Constants::DIM;
		const MKL_INT incx = 1;
 
		vc[2] = 0;

		T* perpendicular_dir = (T*)mkl_calloc(Constants::DIM, sizeof(T), Constants::ALIGNMENT);
		T* tangential_dir = (T*)mkl_calloc(Constants::DIM, sizeof(T), Constants::ALIGNMENT);
		T* radial_vector = (T*)mkl_calloc(Constants::DIM, sizeof(T), Constants::ALIGNMENT);

		radial_vector[0] = pcc[0] - center_of_circle[0];
		radial_vector[1] = pcc[1] - center_of_circle[1];
		radial_vector[2] = 0;

		perpendicular_dir[2] = 1;

		T radius = mkl<T>::nrm2(Constants::DIM, radial_vector, 1);

		if (abs(radius - radius_param) > 0.01)
			std::cout << "Warning! the initial position of the car is not on the trajectory provided in the circular path. \n The expected radius is " << radius_circular_path << ", but the car is at an initial distance of " << radius << " from the center of the circle.\n The execution procedes with the current spatial configuration and with the current distance to the center of the circle." << std::endl;

		T inv_radius = 1. / radius;

		mkl<T>::scal(dim, inv_radius, radial_vector, incx);

		MathLibrary::crossProduct(radial_vector, perpendicular_dir, tangential_dir);

		T magnitude = mkl<T>::dot(dim, vc, incx, tangential_dir, incx);

		mkl<T>::copy(Constants::DIM, tangential_dir, 1, vc, 1);

		mkl<T>::scal(dim, magnitude, vc, incx);

		mkl<T>::scal(dim, radius, radial_vector, incx);

		MathLibrary::crossProduct(radial_vector, vc, omega);

		mkl<T>::scal(Constants::DIM, inv_radius * inv_radius, omega, 1);

		MathLibrary::crossProduct(omega, pt_fl, vt_fl);
		MathLibrary::crossProduct(omega, pt_fr, vt_fr);
		MathLibrary::crossProduct(omega, pt_rl, vt_rl);
		MathLibrary::crossProduct(omega, pt_rr, vt_rr);

		mkl<T>::copy(Constants::DIM, vt_fl, 1, vw_fl, 1);
		mkl<T>::copy(Constants::DIM, vt_fr, 1, vw_fr, 1);
		mkl<T>::copy(Constants::DIM, vt_rl, 1, vw_rl, 1);
		mkl<T>::copy(Constants::DIM, vt_rr, 1, vw_rr, 1);
	}


	/*
	Fixes the tyres to their initial position
	The forces acting on the tyres are now always zero
	\param v is the velocity of the tyre
	\param m is the mass of the tyre
	\param p is the global position of the tyre
	\result Fr The only force acting on the tyre
	*/
	void get_fixed_road_force(T* Fr_fl, T* Fr_fr, T* Fr_rl, T* Fr_rr) {
		Fr_fl[0] = 0.0;
		Fr_fl[1] = 0.0;
		Fr_fl[2] = 0.0;

		Fr_fr[0] = 0.0;
		Fr_fr[1] = 0.0;
		Fr_fr[2] = 0.0;

		Fr_rl[0] = 0.0;
		Fr_rl[1] = 0.0;
		Fr_rl[2] = 0.0;

		Fr_rr[0] = 0.0;
		Fr_rr[1] = 0.0;
		Fr_rr[2] = 0.0;
	}


	/*
	No interaction with the road, no additional forces on the tyres
	\param v is the velocity of the tyre
	\param m is the mass of the tyre
	\param p is the global position of the tyre
	\result Fr The only force acting on the tyre
	*/
	void get_nonfixed_road_force(T* Fr_fl, T* Fr_fr, T* Fr_rl, T* Fr_rr) {
	}

	/* calculates the force in the tyre only with respect to its velocity, mass and position
	\param v is the velocity of the tyre
	\param m is the mass of the tyre
	\param p is the global position of the tyre
	\result Fr The only force acting on the tyre
	*/
	void get_circular_road_force(T* Fr, T* v, T &m, T* p) {

		T *unit_z_vector, *velocity_direction_tyre;
		T velocity_magnitude_tyre, inv_radius, force_magnitude_tyre;

		unit_z_vector = (T*)mkl_calloc(Constants::DIM, sizeof(T), Constants::ALIGNMENT);
		velocity_direction_tyre = (T*)mkl_calloc(Constants::DIM, sizeof(T), Constants::ALIGNMENT);


		const MKL_INT mkl_DIM = Constants::DIM;
		const MKL_INT mkl_incx = 1;
		const MKL_INT mkl_incy = 1;

		mkl<T>::copy(Constants::DIM, p, 1, Fr, 1);

		Fr[2] = 0;		// path only in XZ-plane

		inv_radius = 1.0 / mkl<T>::nrm2(Constants::DIM, p, 1);		// corresponds to the (inverse) radius of the trajectory at the considered tyre

		mkl<T>::scal(Constants::DIM, -inv_radius, Fr, 1);

		unit_z_vector[0] = 0.;
		unit_z_vector[1] = 0.;
		unit_z_vector[2] = 1.;

		MathLibrary::crossProduct(Fr, unit_z_vector, velocity_direction_tyre);

		velocity_magnitude_tyre = mkl<T>::dot(mkl_DIM, v, mkl_incx, velocity_direction_tyre, mkl_incy);

		force_magnitude_tyre = m * velocity_magnitude_tyre * velocity_magnitude_tyre * inv_radius;

		mkl<T>::scal(Constants::DIM, force_magnitude_tyre, Fr, 1);

		MKL_free(unit_z_vector);
		MKL_free(velocity_direction_tyre);
	}


	/*
	Memory allocation of all the variables required in the solve function
	To increase performance by removing repeting memory allocations
	The same locations are overwritten at each timestep
	*/
	void compute_f_mem_alloc() {
		// add to members


		cf_C_cN = (T*)mkl_calloc((Constants::DIM) * (Constants::DIM), sizeof(T), Constants::ALIGNMENT);
		cf_r_up_fl = (T*)mkl_calloc(Constants::DIM, sizeof(T), Constants::ALIGNMENT);
		cf_r_up_fr = (T*)mkl_calloc(Constants::DIM, sizeof(T), Constants::ALIGNMENT);
		cf_r_up_rl = (T*)mkl_calloc(Constants::DIM, sizeof(T), Constants::ALIGNMENT);	
		cf_r_up_rr = (T*)mkl_calloc(Constants::DIM, sizeof(T), Constants::ALIGNMENT);
		cf_r_low_fl = (T*)mkl_calloc(Constants::DIM, sizeof(T), Constants::ALIGNMENT);
		cf_r_low_fr = (T*)mkl_calloc(Constants::DIM, sizeof(T), Constants::ALIGNMENT);
		cf_r_low_rl = (T*)mkl_calloc(Constants::DIM, sizeof(T), Constants::ALIGNMENT);
		cf_r_low_rr = (T*)mkl_calloc(Constants::DIM, sizeof(T), Constants::ALIGNMENT);

		cf_upper_normal_fl = (T*)mkl_calloc((Constants::DIM), sizeof(T), Constants::ALIGNMENT);
		cf_upper_normal_fr = (T*)mkl_calloc((Constants::DIM), sizeof(T), Constants::ALIGNMENT);
		cf_upper_normal_rl = (T*)mkl_calloc((Constants::DIM), sizeof(T), Constants::ALIGNMENT);
		cf_upper_normal_rr = (T*)mkl_calloc((Constants::DIM), sizeof(T), Constants::ALIGNMENT);
		cf_lower_normal_fl = (T*)mkl_calloc((Constants::DIM), sizeof(T), Constants::ALIGNMENT);
		cf_lower_normal_fr = (T*)mkl_calloc((Constants::DIM), sizeof(T), Constants::ALIGNMENT);
		cf_lower_normal_rl = (T*)mkl_calloc((Constants::DIM), sizeof(T), Constants::ALIGNMENT);
		cf_lower_normal_rr = (T*)mkl_calloc((Constants::DIM), sizeof(T), Constants::ALIGNMENT);
		cf_col_dat = (T*)mkl_calloc((Constants::DIM), sizeof(T), Constants::ALIGNMENT);

		cf_upper_force_fl = (T*)mkl_calloc((Constants::DIM), sizeof(T), Constants::ALIGNMENT);
		cf_upper_force_fr = (T*)mkl_calloc((Constants::DIM), sizeof(T), Constants::ALIGNMENT);
		cf_upper_force_rl = (T*)mkl_calloc((Constants::DIM), sizeof(T), Constants::ALIGNMENT);
		cf_upper_force_rr = (T*)mkl_calloc((Constants::DIM), sizeof(T), Constants::ALIGNMENT);
		cf_lower_force_fl = (T*)mkl_calloc((Constants::DIM), sizeof(T), Constants::ALIGNMENT);
		cf_lower_force_fr = (T*)mkl_calloc((Constants::DIM), sizeof(T), Constants::ALIGNMENT);
		cf_lower_force_rl = (T*)mkl_calloc((Constants::DIM), sizeof(T), Constants::ALIGNMENT);
		cf_lower_force_rr = (T*)mkl_calloc((Constants::DIM), sizeof(T), Constants::ALIGNMENT);

		cf_upper_dampf_fl = (T*)mkl_calloc((Constants::DIM), sizeof(T), Constants::ALIGNMENT);
		cf_upper_dampf_fr = (T*)mkl_calloc((Constants::DIM), sizeof(T), Constants::ALIGNMENT);
		cf_upper_dampf_rl = (T*)mkl_calloc((Constants::DIM), sizeof(T), Constants::ALIGNMENT);
		cf_upper_dampf_rr = (T*)mkl_calloc((Constants::DIM), sizeof(T), Constants::ALIGNMENT);
		cf_lower_dampf_fl = (T*)mkl_calloc((Constants::DIM), sizeof(T), Constants::ALIGNMENT);
		cf_lower_dampf_fr = (T*)mkl_calloc((Constants::DIM), sizeof(T), Constants::ALIGNMENT);
		cf_lower_dampf_rl = (T*)mkl_calloc((Constants::DIM), sizeof(T), Constants::ALIGNMENT);
		cf_lower_dampf_rr = (T*)mkl_calloc((Constants::DIM), sizeof(T), Constants::ALIGNMENT);
		cf_temp = (T*)mkl_calloc((Constants::DIM), sizeof(T), Constants::ALIGNMENT);

		cf_upper_angle_fl = new T;
		cf_upper_angle_fr = new T;
		cf_upper_angle_rl = new T;
		cf_upper_angle_rr = new T;
		cf_lower_angle_fl = new T;
		cf_lower_angle_fr = new T;
		cf_lower_angle_rl = new T;
		cf_lower_angle_rr = new T;

		cf_upper_S_fl = (T*)mkl_calloc((Constants::DIM), sizeof(T), Constants::ALIGNMENT);
		cf_upper_S_fr = (T*)mkl_calloc((Constants::DIM), sizeof(T), Constants::ALIGNMENT);
		cf_upper_S_rl = (T*)mkl_calloc((Constants::DIM), sizeof(T), Constants::ALIGNMENT);
		cf_upper_S_rr = (T*)mkl_calloc((Constants::DIM), sizeof(T), Constants::ALIGNMENT);
		cf_lower_S_fl = (T*)mkl_calloc((Constants::DIM), sizeof(T), Constants::ALIGNMENT);
		cf_lower_S_fr = (T*)mkl_calloc((Constants::DIM), sizeof(T), Constants::ALIGNMENT);
		cf_lower_S_rl = (T*)mkl_calloc((Constants::DIM), sizeof(T), Constants::ALIGNMENT);
		cf_lower_S_rr = (T*)mkl_calloc((Constants::DIM), sizeof(T), Constants::ALIGNMENT);

		cf_lower_rot_force_fl = (T*)mkl_calloc((Constants::DIM), sizeof(T), Constants::ALIGNMENT);
		cf_lower_rot_force_fr = (T*)mkl_calloc((Constants::DIM), sizeof(T), Constants::ALIGNMENT);
		cf_lower_rot_force_rl = (T*)mkl_calloc((Constants::DIM), sizeof(T), Constants::ALIGNMENT);
		cf_lower_rot_force_rr = (T*)mkl_calloc((Constants::DIM), sizeof(T), Constants::ALIGNMENT);
		cf_upper_rot_force_fl = (T*)mkl_calloc((Constants::DIM), sizeof(T), Constants::ALIGNMENT);
		cf_upper_rot_force_fr = (T*)mkl_calloc((Constants::DIM), sizeof(T), Constants::ALIGNMENT);
		cf_upper_rot_force_rl = (T*)mkl_calloc((Constants::DIM), sizeof(T), Constants::ALIGNMENT);
		cf_upper_rot_force_rr = (T*)mkl_calloc((Constants::DIM), sizeof(T), Constants::ALIGNMENT);
		cf_car_rot_force_fl = (T*)mkl_calloc((Constants::DIM), sizeof(T), Constants::ALIGNMENT);
		cf_car_rot_force_fr = (T*)mkl_calloc((Constants::DIM), sizeof(T), Constants::ALIGNMENT);
		cf_car_rot_force_rl = (T*)mkl_calloc((Constants::DIM), sizeof(T), Constants::ALIGNMENT);
		cf_car_rot_force_rr = (T*)mkl_calloc((Constants::DIM), sizeof(T), Constants::ALIGNMENT);
		cf_sum_car_force_fl = (T*)mkl_calloc((Constants::DIM), sizeof(T), Constants::ALIGNMENT);
		cf_sum_car_force_fr = (T*)mkl_calloc((Constants::DIM), sizeof(T), Constants::ALIGNMENT);
		cf_sum_car_force_rl = (T*)mkl_calloc((Constants::DIM), sizeof(T), Constants::ALIGNMENT);
		cf_sum_car_force_rr = (T*)mkl_calloc((Constants::DIM), sizeof(T), Constants::ALIGNMENT);

		cf_local_FR_fl = (T*)mkl_calloc((Constants::DIM), sizeof(T), Constants::ALIGNMENT);
		cf_local_FR_fr = (T*)mkl_calloc((Constants::DIM), sizeof(T), Constants::ALIGNMENT);
		cf_local_FR_rl = (T*)mkl_calloc((Constants::DIM), sizeof(T), Constants::ALIGNMENT);
		cf_local_FR_rr = (T*)mkl_calloc((Constants::DIM), sizeof(T), Constants::ALIGNMENT);

		cf_Hc = (T*)mkl_calloc(Constants::DIM, sizeof(T), Constants::ALIGNMENT);
		cf_sum_torque_spring_car = (T*)mkl_calloc(Constants::DIM, sizeof(T), Constants::ALIGNMENT);
		cf_Tc = (T*)mkl_calloc(Constants::DIM, sizeof(T), Constants::ALIGNMENT);
		cf_wc_tilda = (T*)mkl_calloc((Constants::DIM) * (Constants::DIM), sizeof(T), Constants::ALIGNMENT);

		cf_Tc[0] = 0.0;
		cf_Tc[1] = 0.0;
		cf_Tc[2] = 0.0;

		cf_b_rem = (T*)mkl_calloc(9 * Constants::DIM, sizeof(T), Constants::ALIGNMENT);
		cf_Qc = (T*)mkl_calloc((Constants::NUM_LEGS) * (Constants::DIM), sizeof(T), Constants::ALIGNMENT);
		cf_qc_dot = (T*)mkl_calloc((Constants::NUM_LEGS), sizeof(T), Constants::ALIGNMENT);

		// required for the look up table
		current_spring_lengths = (T*)mkl_calloc(8 * Constants::DIM, sizeof(T), Constants::ALIGNMENT);
		stiffness_vector = (T*)mkl_calloc(8 * Constants::DIM, sizeof(T), Constants::ALIGNMENT);

	}

	/*
	Clean the memory allocated for the main solver
	*/
	void compute_f_clean() {
		mkl_free_buffers();
		mkl_free(cf_C_cN);
		mkl_free(cf_r_up_fl);
		mkl_free(cf_r_up_fr);
		mkl_free(cf_r_up_rl);
		mkl_free(cf_r_up_rr);
		mkl_free(cf_r_low_fl);
		mkl_free(cf_r_low_fr);
		mkl_free(cf_r_low_rl);
		mkl_free(cf_r_low_rr);
		mkl_free(cf_upper_normal_fl);
		mkl_free(cf_upper_normal_fr);
		mkl_free(cf_upper_normal_rl);
		mkl_free(cf_upper_normal_rr);
		mkl_free(cf_lower_normal_fl);
		mkl_free(cf_lower_normal_fr);
		mkl_free(cf_lower_normal_rl);
		mkl_free(cf_lower_normal_rr);
		mkl_free(cf_col_dat);
		mkl_free(cf_upper_force_fl);
		mkl_free(cf_upper_force_fr);
		mkl_free(cf_upper_force_rl);
		mkl_free(cf_upper_force_rr);
		mkl_free(cf_lower_force_fl);
		mkl_free(cf_lower_force_fr);
		mkl_free(cf_lower_force_rl);
		mkl_free(cf_lower_force_rr);
		mkl_free(cf_upper_dampf_fl);
		mkl_free(cf_upper_dampf_fr);
		mkl_free(cf_upper_dampf_rl);
		mkl_free(cf_upper_dampf_rr);
		mkl_free(cf_lower_dampf_fl);
		mkl_free(cf_lower_dampf_fr);
		mkl_free(cf_lower_dampf_rl);
		mkl_free(cf_lower_dampf_rr);
		mkl_free(cf_temp);
		delete cf_upper_angle_fl;
		delete cf_upper_angle_fr;
		delete cf_upper_angle_rl;
		delete cf_upper_angle_rr;
		delete cf_lower_angle_fl;
		delete cf_lower_angle_fr;
		delete cf_lower_angle_rl;
		delete cf_lower_angle_rr;
		mkl_free(cf_upper_S_fl);
		mkl_free(cf_upper_S_fr);
		mkl_free(cf_upper_S_rl);
		mkl_free(cf_upper_S_rr);
		mkl_free(cf_lower_S_fl);
		mkl_free(cf_lower_S_fr);
		mkl_free(cf_lower_S_rl);
		mkl_free(cf_lower_S_rr);
		mkl_free(cf_lower_rot_force_fl);
		mkl_free(cf_lower_rot_force_fr);
		mkl_free(cf_lower_rot_force_rl);
		mkl_free(cf_lower_rot_force_rr);
		mkl_free(cf_upper_rot_force_fl);
		mkl_free(cf_upper_rot_force_fr);
		mkl_free(cf_upper_rot_force_rl);
		mkl_free(cf_upper_rot_force_rr);
		mkl_free(cf_car_rot_force_fl);
		mkl_free(cf_car_rot_force_fr);
		mkl_free(cf_car_rot_force_rl);
		mkl_free(cf_car_rot_force_rr);
		mkl_free(cf_sum_car_force_fl);
		mkl_free(cf_sum_car_force_fr);
		mkl_free(cf_sum_car_force_rl);
		mkl_free(cf_sum_car_force_rr);
		mkl_free(cf_local_FR_fl);
		mkl_free(cf_local_FR_fr);
		mkl_free(cf_local_FR_rl);
		mkl_free(cf_local_FR_rr);
		mkl_free(cf_Hc);
		mkl_free(cf_sum_torque_spring_car);
		mkl_free(cf_Tc);
		mkl_free(cf_wc_tilda);
		mkl_free(cf_b_rem);
		mkl_free(cf_Qc);
		mkl_free(cf_qc_dot);
		mkl_free(current_spring_lengths);
		mkl_free(stiffness_vector);
	}

	/*
	Solver which is called at each time step
	Computes the forces and torques on each point mass and computes the right hand side of the ODE
	\param x_ current solution of the system
	\param t_ current simulation time
	\return f_ the load vector
	*/
	void compute_f3D_reduced(T* x_, T t_, T* f_) {
		/*
		Small performance gain might be possible by transforming C_cN to column major
		Note: corresponding MKL function call have to be changed too
		*/
		// Aliasing for readability
		T* wc_ = x_;
		T* vc_ = x_ + Constants::DIM;
		T* vw_fl_ = x_ + 2 * Constants::DIM;
		T* vw_fr_ = x_ + 3 * Constants::DIM;
		T* vw_rl_ = x_ + 4 * Constants::DIM;
		T* vw_rr_ = x_ + 5 * Constants::DIM;
		T* vt_fl_ = x_ + 6 * Constants::DIM;
		T* vt_fr_ = x_ + 7 * Constants::DIM;
		T* vt_rl_ = x_ + 8 * Constants::DIM;
		T* vt_rr_ = x_ + 9 * Constants::DIM;
		T* qc_ = x_ + 10 * Constants::DIM;
		T* pcc_ = x_ + 10 * Constants::DIM + Constants::NUM_LEGS;
		T* pw_fl_ = x_ + 11 * Constants::DIM + Constants::NUM_LEGS;
		T* pw_fr_ = x_ + 12 * Constants::DIM + Constants::NUM_LEGS;
		T* pw_rl_ = x_ + 13 * Constants::DIM + Constants::NUM_LEGS;
		T* pw_rr_ = x_ + 14 * Constants::DIM + Constants::NUM_LEGS;
		T* pt_fl_ = x_ + 15 * Constants::DIM + Constants::NUM_LEGS;
		T* pt_fr_ = x_ + 16 * Constants::DIM + Constants::NUM_LEGS;
		T* pt_rl_ = x_ + 17 * Constants::DIM + Constants::NUM_LEGS;
		T* pt_rr_ = x_ + 18 * Constants::DIM + Constants::NUM_LEGS;
		T inv_norm_r_up_fl, inv_norm_r_up_fr, inv_norm_r_up_rl, inv_norm_r_up_rr;
		T inv_norm_r_low_fl, inv_norm_r_low_fr, inv_norm_r_low_rl, inv_norm_r_low_rr;
		T norm_r_up_fl, norm_r_up_fr, norm_r_up_rl, norm_r_up_rr;
		T norm_r_low_fl, norm_r_low_fr, norm_r_low_rl, norm_r_low_rr;
		//write_vector(r_fl, Constants::DIM);
		//write_matrix(r_fl_tilda, Constants::DIM);
		/*
		/////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////                 Basis        //////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////
		get cosine transforms (C_Nc means r_N = C_Nc * r_c)
		compute local base vectors
		basis_c = get_basis(qc);
		*/

		MathLibrary::get_basis<T>(qc_, cf_C_cN);
		const MKL_INT mkl_DIM = Constants::DIM;
		const MKL_INT mkl_incx = 1;
		const MKL_INT mkl_incy = 1;

		/*
		% global positions of the tyre connections
			pc1 = C_cN' * r1 + pcc;
			pc2 = C_cN' * r2 + pcc;
			pc3 = C_cN' * r3 + pcc;
			pc4 = C_cN' * r4 + pcc;

			% currrent length of the legs
			r_up1 = pc1 - pw1;
			r_up2 = pc2 - pw2;
			r_up3 = pc3 - pw3;
			r_up4 = pc4 - pw4;

			r_low1 = pw1 - pt1;
			r_low2 = pw2 - pt2;
			r_low3 = pw3 - pt3;
			r_low4 = pw4 - pt4;
		*/
		mkl<T>::copy(Constants::DIM, pcc_, 1, cf_r_up_fl, 1);
		mkl<T>::gemv(CblasRowMajor, CblasNoTrans, Constants::DIM, Constants::DIM, 1, cf_C_cN, Constants::DIM, this->r_fl, 1, 1, cf_r_up_fl, 1);
		mkl<T>::axpy(Constants::DIM, -1.0, pw_fl_, 1, cf_r_up_fl, 1);

		// r_up_fr
		mkl<T>::copy(Constants::DIM, pcc_, 1, cf_r_up_fr, 1);
		mkl<T>::gemv(CblasRowMajor, CblasNoTrans, Constants::DIM, Constants::DIM, 1, cf_C_cN, Constants::DIM, this->r_fr, 1, 1, cf_r_up_fr, 1);
		mkl<T>::axpy(Constants::DIM, -1.0, pw_fr_, 1, cf_r_up_fr, 1);
		// r_up_rl
		mkl<T>::copy(Constants::DIM, pcc_, 1, cf_r_up_rl, 1);
		mkl<T>::gemv(CblasRowMajor, CblasNoTrans, Constants::DIM, Constants::DIM, 1, cf_C_cN, Constants::DIM, this->r_rl, 1, 1, cf_r_up_rl, 1);
		mkl<T>::axpy(Constants::DIM, -1.0, pw_rl_, 1, cf_r_up_rl, 1);
		// r_up4
		mkl<T>::copy(Constants::DIM, pcc_, 1, cf_r_up_rr, 1);
		mkl<T>::gemv(CblasRowMajor, CblasNoTrans, Constants::DIM, Constants::DIM, 1, cf_C_cN, Constants::DIM, this->r_rr, 1, 1, cf_r_up_rr, 1);
		mkl<T>::axpy(Constants::DIM, -1.0, pw_rr_, 1, cf_r_up_rr, 1);
		
		// r_low1
		mkl<T>::copy(Constants::DIM, pw_fl_, 1, cf_r_low_fl, 1);
		mkl<T>::axpy(Constants::DIM, -1.0, pt_fl_, 1, cf_r_low_fl, 1);

		// r_low2
		mkl<T>::copy(Constants::DIM, pw_fr_, 1, cf_r_low_fr, 1);
		mkl<T>::axpy(Constants::DIM, -1.0, pt_fr_, 1, cf_r_low_fr, 1);

		// r_low3
		mkl<T>::copy(Constants::DIM, pw_rl_, 1, cf_r_low_rl, 1);
		mkl<T>::axpy(Constants::DIM, -1.0, pt_rl_, 1, cf_r_low_rl, 1);
		// r_low4
		mkl<T>::copy(Constants::DIM, pw_rr_, 1, cf_r_low_rr, 1);
		mkl<T>::axpy(Constants::DIM, -1.0, pt_rr_, 1, cf_r_low_rr, 1);

		/* Compute spring lengths and their inverses */
		norm_r_up_fl = mkl<T>::nrm2(Constants::DIM, cf_r_up_fl, 1);
		norm_r_up_fr = mkl<T>::nrm2(Constants::DIM, cf_r_up_fr, 1);
		norm_r_up_rl = mkl<T>::nrm2(Constants::DIM, cf_r_up_rl, 1);
		norm_r_up_rr = mkl<T>::nrm2(Constants::DIM, cf_r_up_rr, 1);

		inv_norm_r_up_fl = 1.0 / norm_r_up_fl;
		inv_norm_r_up_fr = 1.0 / norm_r_up_fr;
		inv_norm_r_up_rl = 1.0 / norm_r_up_rl;
		inv_norm_r_up_rr = 1.0 / norm_r_up_rr;

		norm_r_low_fl = mkl<T>::nrm2(Constants::DIM, cf_r_low_fl, 1);
		norm_r_low_fr = mkl<T>::nrm2(Constants::DIM, cf_r_low_fr, 1);
		norm_r_low_rl = mkl<T>::nrm2(Constants::DIM, cf_r_low_rl, 1);
		norm_r_low_rr = mkl<T>::nrm2(Constants::DIM, cf_r_low_rr, 1);

		inv_norm_r_low_fl = 1.0 / norm_r_low_fl;
		inv_norm_r_low_fr = 1.0 / norm_r_low_fr;
		inv_norm_r_low_rl = 1.0 / norm_r_low_rl;
		inv_norm_r_low_rr = 1.0 / norm_r_low_rr;


		/* Compute stiffness from lookup table*/
		if (use_interpolation) {
			//populate the lenght_vector
			current_spring_lengths[0] = norm_r_up_fl;
			current_spring_lengths[1] = norm_r_low_fl;
			current_spring_lengths[2] = norm_r_up_fr;
			current_spring_lengths[3] = norm_r_low_fr;
			current_spring_lengths[4] = norm_r_up_rl;
			current_spring_lengths[5] = norm_r_low_rl;
			current_spring_lengths[6] = norm_r_up_rr;
			current_spring_lengths[7] = norm_r_low_rr;
			// calculate the new stiffnesses
			lookupStiffness->getInterpolation(current_spring_lengths, stiffness_vector);

			// overwrite stiffness values
			upper_spring_stiffness[0] = stiffness_vector[0];
			lower_spring_stiffness[0] = stiffness_vector[1];
			upper_spring_stiffness[1] = stiffness_vector[2];
			lower_spring_stiffness[1] = stiffness_vector[3];
			upper_spring_stiffness[2] = stiffness_vector[4];
			lower_spring_stiffness[2] = stiffness_vector[5];
			upper_spring_stiffness[3] = stiffness_vector[6];
			lower_spring_stiffness[3] = stiffness_vector[7];
		}

		/*
		get angle and normal vectors at the legs

			[~, upper_angle1, upper_normal1] = get_quaternion(r_up1, C_cN(2,:)');
			[~, upper_angle2, upper_normal2] = get_quaternion(r_up2, C_cN(2,:)');
			[~, upper_angle3, upper_normal3] = get_quaternion(r_up3, C_cN(2,:)');
			[~, upper_angle4, upper_normal4] = get_quaternion(r_up4, C_cN(2,:)');

			[~, lower_angle1, lower_normal1] = get_quaternion(r_low1, r_up1);
			[~, lower_angle2, lower_normal2] = get_quaternion(r_low2, r_up2);
			[~, lower_angle3, lower_normal3] = get_quaternion(r_low3, r_up3);
			[~, lower_angle4, lower_normal4] = get_quaternion(r_low4, r_up4);

		*/

		mkl<T>::copy(Constants::DIM, cf_C_cN + 1, Constants::DIM, cf_col_dat, 1);

		MathLibrary::get_quaternion<T>(cf_r_up_fl, cf_col_dat, cf_upper_angle_fl, cf_upper_normal_fl, Constants::DIM);
		MathLibrary::get_quaternion<T>(cf_r_up_fr, cf_col_dat, cf_upper_angle_fr, cf_upper_normal_fr, Constants::DIM);
		MathLibrary::get_quaternion<T>(cf_r_up_rl, cf_col_dat, cf_upper_angle_rl, cf_upper_normal_rl, Constants::DIM);
		MathLibrary::get_quaternion<T>(cf_r_up_rr, cf_col_dat, cf_upper_angle_rr, cf_upper_normal_rr, Constants::DIM);

		MathLibrary::get_quaternion(cf_r_low_fl, cf_r_up_fl, cf_lower_angle_fl, cf_lower_normal_fl, Constants::DIM);
		MathLibrary::get_quaternion(cf_r_low_fr, cf_r_up_fr, cf_lower_angle_fr, cf_lower_normal_fr, Constants::DIM);
		MathLibrary::get_quaternion(cf_r_low_rl, cf_r_up_rl, cf_lower_angle_rl, cf_lower_normal_rl, Constants::DIM);
		MathLibrary::get_quaternion(cf_r_low_rr, cf_r_up_rr, cf_lower_angle_rr, cf_lower_normal_rr, Constants::DIM);


		/*
		///////////////////////////////////////////////////////////////
		////////////           Forces  and Torques          ///////////
		///////////////////////////////////////////////////////////////
		calculate the elongational spring forces (in global basis)

		upper_force1 = upper_spring_stiffness(1) * (r_up1) * (1 - upper_spring_length(1) * inv_norm_r_up1);
		upper_force2 = upper_spring_stiffness(2) * (r_up2) * (1 - upper_spring_length(2) * inv_norm_r_up2);
		upper_force3 = upper_spring_stiffness(3) * (r_up3) * (1 - upper_spring_length(3) * inv_norm_r_up3);
		upper_force4 = upper_spring_stiffness(4) * (r_up4) * (1 - upper_spring_length(4) * inv_norm_r_up4);

		lower_force1 = lower_spring_stiffness(1) * (r_low1) * (1 - lower_spring_length(1) * inv_norm_r_low1);
		lower_force2 = lower_spring_stiffness(2) * (r_low2) * (1 - lower_spring_length(2) * inv_norm_r_low2);
		lower_force3 = lower_spring_stiffness(3) * (r_low3) * (1 - lower_spring_length(3) * inv_norm_r_low3);
		lower_force4 = lower_spring_stiffness(4) * (r_low4) * (1 - lower_spring_length(4) * inv_norm_r_low4);

		*/

		T scale;
		mkl<T>::copy(Constants::DIM, cf_r_up_fl, 1, cf_upper_force_fl, 1);
		scale = this->upper_spring_stiffness[0] * (1.0 - this->upper_spring_length[0] * inv_norm_r_up_fl);
		mkl<T>::scal(Constants::DIM, scale, cf_upper_force_fl, 1);
		mkl<T>::copy(Constants::DIM, cf_r_up_fr, 1, cf_upper_force_fr, 1);
		scale = this->upper_spring_stiffness[1] * (1.0 - this->upper_spring_length[1] * inv_norm_r_up_fr);
		mkl<T>::scal(Constants::DIM, scale, cf_upper_force_fr, 1);
		mkl<T>::copy(Constants::DIM, cf_r_up_rl, 1, cf_upper_force_rl, 1);
		scale = this->upper_spring_stiffness[2] * (1.0 - this->upper_spring_length[2] * inv_norm_r_up_rl);
		mkl<T>::scal(Constants::DIM, scale, cf_upper_force_rl, 1);
		mkl<T>::copy(Constants::DIM, cf_r_up_rr, 1, cf_upper_force_rr, 1);
		scale = this->upper_spring_stiffness[3] * (1.0 - this->upper_spring_length[3] * inv_norm_r_up_rr);
		mkl<T>::scal(Constants::DIM, scale, cf_upper_force_rr, 1);

		mkl<T>::copy(Constants::DIM, cf_r_low_fl, 1, cf_lower_force_fl, 1);
		scale = this->lower_spring_stiffness[0] * (1.0 - this->lower_spring_length[0] * inv_norm_r_low_fl);
		mkl<T>::scal(Constants::DIM, scale, cf_lower_force_fl, 1);
		mkl<T>::copy(Constants::DIM, cf_r_low_fr, 1, cf_lower_force_fr, 1);
		scale = this->lower_spring_stiffness[1] * (1.0 - this->lower_spring_length[1] * inv_norm_r_low_fr);
		mkl<T>::scal(Constants::DIM, scale, cf_lower_force_fr, 1);
		mkl<T>::copy(Constants::DIM, cf_r_low_rl, 1, cf_lower_force_rl, 1);
		scale = this->lower_spring_stiffness[2] * (1.0 - this->lower_spring_length[2] * inv_norm_r_low_rl);
		mkl<T>::scal(Constants::DIM, scale, cf_lower_force_rl, 1);
		mkl<T>::copy(Constants::DIM, cf_r_low_rr, 1, cf_lower_force_rr, 1);
		scale = this->lower_spring_stiffness[3] * (1.0 - this->lower_spring_length[3] * inv_norm_r_low_rr);
		mkl<T>::scal(Constants::DIM, scale, cf_lower_force_rr, 1);

		/*
		calculate forces from damping effects
			upper_vdiff1 = (dot((vc - C_cN' * (r1_tilda * wc)), r_up1) - dot(vw1, r_up1)) * r_up1 * inv_norm_r_up1 * inv_norm_r_up1;
			upper_vdiff2 = (dot((vc - C_cN' * (r2_tilda * wc)), r_up2) - dot(vw2, r_up2)) * r_up2 * inv_norm_r_up2 * inv_norm_r_up2;
			upper_vdiff3 = (dot((vc - C_cN' * (r3_tilda * wc)), r_up3) - dot(vw3, r_up3)) * r_up3 * inv_norm_r_up3 * inv_norm_r_up3;
			upper_vdiff4 = (dot((vc - C_cN' * (r4_tilda * wc)), r_up4) - dot(vw4, r_up4)) * r_up4 * inv_norm_r_up4 * inv_norm_r_up4;

			lower_vdiff1 = (dot(vw1, r_low1) - dot(vt1, r_low1)) * r_low1 * inv_norm_r_low1 * inv_norm_r_low1;
			lower_vdiff2 = (dot(vw2, r_low2) - dot(vt2, r_low2)) * r_low2 * inv_norm_r_low2 * inv_norm_r_low2;
			lower_vdiff3 = (dot(vw3, r_low3) - dot(vt3, r_low3)) * r_low3 * inv_norm_r_low3 * inv_norm_r_low3;
			lower_vdiff4 = (dot(vw4, r_low4) - dot(vt4, r_low4)) * r_low4 * inv_norm_r_low4 * inv_norm_r_low4;

			upper_dampf1 = upper_spring_damping1(upper_vdiff1);
			upper_dampf2 = upper_spring_damping2(upper_vdiff2);
			upper_dampf3 = upper_spring_damping3(upper_vdiff3);
			upper_dampf4 = upper_spring_damping4(upper_vdiff4);

			lower_dampf1 = lower_spring_damping1(lower_vdiff1);
			lower_dampf2 = lower_spring_damping2(lower_vdiff2);
			lower_dampf3 = lower_spring_damping3(lower_vdiff3);
			lower_dampf4 = lower_spring_damping4(lower_vdiff4);

		*/


		//// upper_dampf1
		// compute: vc - C_cN' * (r1_tilda * wc))
		mkl<T>::gemv(CblasRowMajor, CblasNoTrans, Constants::DIM, Constants::DIM, 1.0, this->r_fl_tilda, Constants::DIM, wc_, 1, 0.0, cf_temp, 1);
		mkl<T>::copy(Constants::DIM, vc_, 1, cf_upper_dampf_fl, 1);
		mkl<T>::gemv(CblasRowMajor, CblasNoTrans, Constants::DIM, Constants::DIM, -1.0, cf_C_cN, Constants::DIM, cf_temp, 1, 1.0, cf_upper_dampf_fl, 1);
		// dot((vc - C_cN' * (r1_tilda * wc)), r_up1)
		//std::cout << mkl<T>::dot(Constants::DIM, upper_dampf1, 1, r_up1, 1) << std::endl;
		//std::cout << MathLibrary::dot_product<T>(upper_dampf1, r_up1, Constants::DIM) << std::endl;
		scale = mkl<T>::dot(mkl_DIM, cf_upper_dampf_fl, mkl_incx, cf_r_up_fl, mkl_incy);
		//scale = MathLibrary::dot_product<T>(upper_dampf1, r_up1, Constants::DIM);

		// dot((vc - C_cN' * (r1_tilda * wc)), r_up1) - dot(vw1, r_up1)
		scale -= mkl<T>::dot(mkl_DIM, vw_fl_, mkl_incx, cf_r_up_fl, mkl_incx);
		//scale -= MathLibrary::dot_product<T>(vw1_, r_up1, Constants::DIM);
		// (dot((vc - C_cN' * (r1_tilda * wc)), r_up1) - dot(vw1, r_up1))* inv_norm_r_up1 * inv_norm_r_up1
		scale = scale * inv_norm_r_up_fl * inv_norm_r_up_fl * this->upper_spring_damping[0];
		mkl<T>::copy(Constants::DIM, cf_r_up_fl, 1, cf_upper_dampf_fl, 1);
		mkl<T>::scal(Constants::DIM, scale, cf_upper_dampf_fl, 1);


		//// upper_dampf_fr
		// compute: vc - C_cN' * (r_fr_tilda * wc))
		mkl<T>::gemv(CblasRowMajor, CblasNoTrans, Constants::DIM, Constants::DIM, 1.0, this->r_fr_tilda, Constants::DIM, wc_, 1, 0.0, cf_temp, 1);
		mkl<T>::copy(Constants::DIM, vc_, 1, cf_upper_dampf_fr, 1);
		mkl<T>::gemv(CblasRowMajor, CblasNoTrans, Constants::DIM, Constants::DIM, -1.0, cf_C_cN, Constants::DIM, cf_temp, 1, 1.0, cf_upper_dampf_fr, 1);
		// dot((vc - C_cN' * (r_fr_tilda * wc)), r_up_fr)
		scale = mkl<T>::dot(mkl_DIM, cf_upper_dampf_fr, mkl_incx, cf_r_up_fr, mkl_incy);
		//scale = MathLibrary::dot_product<T>(upper_dampf_fr, r_up_fr, Constants::DIM);
		// dot((vc - C_cN' * (r_fr_tilda * wc)), r_up_fr) - dot(vw_fr, r_up_fr)
		scale -= mkl<T>::dot(mkl_DIM, vw_fr_, mkl_incx, cf_r_up_fr, mkl_incy);
		//scale -= MathLibrary::dot_product<T>(vw_fr_, r_up_fr, Constants::DIM);
		// (dot((vc - C_cN' * (r_fr_tilda * wc)), r_up_fr) - dot(vw_fr, r_up_fr))* inv_norm_r_up_fr * inv_norm_r_up_fr
		scale = scale * inv_norm_r_up_fr * inv_norm_r_up_fr * this->upper_spring_damping[1];
		mkl<T>::copy(Constants::DIM, cf_r_up_fr, 1, cf_upper_dampf_fr, 1);
		mkl<T>::scal(Constants::DIM, scale, cf_upper_dampf_fr, 1);

		//// upper_dampf3
		// compute: vc - C_cN' * (r3_tilda * wc))
		mkl<T>::gemv(CblasRowMajor, CblasNoTrans, Constants::DIM, Constants::DIM, 1.0, this->r_rl_tilda, Constants::DIM, wc_, 1, 0.0, cf_temp, 1);
		mkl<T>::copy(Constants::DIM, vc_, 1, cf_upper_dampf_rl, 1);
		mkl<T>::gemv(CblasRowMajor, CblasNoTrans, Constants::DIM, Constants::DIM, -1.0, cf_C_cN, Constants::DIM, cf_temp, 1, 1.0, cf_upper_dampf_rl, 1);
		// dot((vc - C_cN' * (r_rl_tilda * wc)), r_up_rl)
		scale = mkl<T>::dot(mkl_DIM, cf_upper_dampf_rl, mkl_incx, cf_r_up_rl, mkl_incy);
		//scale = MathLibrary::dot_product<T>(upper_dampf_rl, r_up_rl, Constants::DIM);
		// dot((vc - C_cN' * (r_rl_tilda * wc)), r_up_rl) - dot(vw_rl, r_up_rl)
		scale -= mkl<T>::dot(mkl_DIM, vw_rl_, mkl_incx, cf_r_up_rl, mkl_incy);
		//scale -= MathLibrary::dot_product<T>(vw_rl_, r_up_rl, Constants::DIM);
		// (dot((vc - C_cN' * (r_rl_tilda * wc)), r_up_rl) - dot(vw_rl, r_up_rl))* inv_norm_r_up_rl * inv_norm_r_up_rl
		scale = scale * inv_norm_r_up_rl * inv_norm_r_up_rl * this->upper_spring_damping[2];
		mkl<T>::copy(Constants::DIM, cf_r_up_rl, 1, cf_upper_dampf_rl, 1);
		mkl<T>::scal(Constants::DIM, scale, cf_upper_dampf_rl, 1);

		//// upper_dampf4
		// compute: vc - C_cN' * (r4_tilda * wc))
		mkl<T>::gemv(CblasRowMajor, CblasNoTrans, Constants::DIM, Constants::DIM, 1.0, this->r_rr_tilda, Constants::DIM, wc_, 1, 0.0, cf_temp, 1);
		mkl<T>::copy(Constants::DIM, vc_, 1, cf_upper_dampf_rr, 1);
		mkl<T>::gemv(CblasRowMajor, CblasNoTrans, Constants::DIM, Constants::DIM, -1.0, cf_C_cN, Constants::DIM, cf_temp, 1, 1.0, cf_upper_dampf_rr, 1);
		// dot((vc - C_cN' * (r4_tilda * wc)), r_up4)
		scale = mkl<T>::dot(mkl_DIM, cf_upper_dampf_rr, mkl_incx, cf_r_up_rr, mkl_incy);
		//scale = MathLibrary::dot_product<T>(upper_dampf4, r_up4, Constants::DIM);
		// dot((vc - C_cN' * (r4_tilda * wc)), r_up4) - dot(vw4, r_up4)
		scale -= mkl<T>::dot(mkl_DIM, vw_rr_, mkl_incx, cf_r_up_rr, mkl_incy);
		//scale -= MathLibrary::dot_product<T>(r_up4, vw4_, Constants::DIM);
		// (dot((vc - C_cN' * (r4_tilda * wc)), r_up4) - dot(vw4, r_up4))* inv_norm_r_up4 * inv_norm_r_up4
		scale = scale * inv_norm_r_up_rr * inv_norm_r_up_rr * this->upper_spring_damping[3];
		mkl<T>::copy(Constants::DIM, cf_r_up_rr, 1, cf_upper_dampf_rr, 1);
		mkl<T>::scal(Constants::DIM, scale, cf_upper_dampf_rr, 1);

		//// lower_dampf1
		// dot(vw1, r_low1)
		scale = mkl<T>::dot(mkl_DIM, vw_fl_, mkl_incx, cf_r_low_fl, mkl_incy);
		/*scale = MathLibrary::dot_product<T>(vw1_, r_low1, Constants::DIM);*/
		// (dot(vw1, r_low1) - dot(vt1, r_low1))
		scale -= mkl<T>::dot(mkl_DIM, vt_fl_, mkl_incx, cf_r_low_fl, mkl_incx);
		//scale -= MathLibrary::dot_product<T>(vt1_, r_low1, Constants::DIM);
		//(dot(vw1, r_low1) - dot(vt1, r_low1)) * inv_norm_r_low1 * inv_norm_r_low1
		scale = scale * inv_norm_r_low_fl * inv_norm_r_low_fl * this->lower_spring_damping[0];
		mkl<T>::copy(Constants::DIM, cf_r_low_fl, 1, cf_lower_dampf_fl, 1);
		mkl<T>::scal(Constants::DIM, scale, cf_lower_dampf_fl, 1);

		//// lower_dampf2
		// dot(vw2, r_low2)
		scale = mkl<T>::dot(mkl_DIM, vw_fr_, mkl_incx, cf_r_low_fr, mkl_incy);
		//scale = MathLibrary::dot_product<T>(vw2_, r_low2, Constants::DIM);
		// (dot(vw2, r_low2) - dot(vt2, r_low2))
		scale -= mkl<T>::dot(mkl_DIM, vt_fr_, mkl_incx, cf_r_low_fr, mkl_incy);
		//scale -= MathLibrary::dot_product<T>(vt2_, r_low2, Constants::DIM);
		//(dot(vw2, r_low2) - dot(vt2, r_low2)) * inv_norm_r_low2 * inv_norm_r_low2
		scale = scale * inv_norm_r_low_fr * inv_norm_r_low_fr * this->lower_spring_damping[1];
		mkl<T>::copy(Constants::DIM, cf_r_low_fr, 1, cf_lower_dampf_fr, 1);
		mkl<T>::scal(Constants::DIM, scale, cf_lower_dampf_fr, 1);

		//// lower_dampf3
		// dot(vw3, r_low3)
		scale = mkl<T>::dot(mkl_DIM, vw_rl_, mkl_incx, cf_r_low_rl, mkl_incy);
		//scale = MathLibrary::dot_product<T>(vw_rl_, r_low_rl, Constants::DIM);
		// (dot(vw_rl, r_low_rl) - dot(vt_rl, r_low_rl))
		scale -= mkl<T>::dot(mkl_DIM, vt_rl_, mkl_incx, cf_r_low_rl, mkl_incy);
		//scale -= MathLibrary::dot_product<T>(vt_rl_, r_low_rl, Constants::DIM);
		//(dot(vw_rl, r_low_rl) - dot(vt_rl, r_low_rl)) * inv_norm_r_low_rl * inv_norm_r_low_rl
		scale = scale * inv_norm_r_low_rl * inv_norm_r_low_rl * this->lower_spring_damping[2];
		mkl<T>::copy(Constants::DIM, cf_r_low_rl, 1, cf_lower_dampf_rl, 1);
		mkl<T>::scal(Constants::DIM, scale, cf_lower_dampf_rl, 1);

		//// lower_dampf4
		// dot(vw4, r_low4)
		scale = mkl<T>::dot(mkl_DIM, vw_rr_, mkl_incx, cf_r_low_rr, mkl_incy);
		//scale = MathLibrary::dot_product<T>(vw4_, r_low4, Constants::DIM);
		// (dot(vw4, r_low4) - dot(vt4, r_low4))
		scale -= mkl<T>::dot(mkl_DIM, vt_rr_, mkl_incx, cf_r_low_rr, mkl_incy);
		//scale -= MathLibrary::dot_product<T>(vt4_, r_low4, Constants::DIM);
		//(dot(vw4, r_low4) - dot(vt4, r_low4)) * inv_norm_r_low4 * inv_norm_r_low4
		scale = scale * inv_norm_r_low_rr * inv_norm_r_low_rr * this->lower_spring_damping[3];
		mkl<T>::copy(Constants::DIM, cf_r_low_rr, 1, cf_lower_dampf_rr, 1);
		mkl<T>::scal(Constants::DIM, scale, cf_lower_dampf_rr, 1);

		/*
		torque from the rotational spring
		upper_S1 = upper_rotational_stiffness(1) * upper_angle1 * upper_normal1;        % in global basis
		upper_S2 = upper_rotational_stiffness(2) * upper_angle2 * upper_normal2;
		upper_S3 = upper_rotational_stiffness(3) * upper_angle3 * upper_normal3;
		upper_S4 = upper_rotational_stiffness(4) * upper_angle4 * upper_normal4;

		lower_S1 = lower_rotational_stiffness(1) * lower_angle1 * lower_normal1;        % in global basis
		lower_S2 = lower_rotational_stiffness(2) * lower_angle2 * lower_normal2;
		lower_S3 = lower_rotational_stiffness(3) * lower_angle3 * lower_normal3;
		lower_S4 = lower_rotational_stiffness(4) * lower_angle4 * lower_normal4;
		*/



		// upper_S1 = upper_rotational_stiffness(1) * upper_angle1 * upper_normal1;
		scale = (this->upper_rotational_stiffness[0]) * (*cf_upper_angle_fl);
		mkl<T>::copy(Constants::DIM, cf_upper_normal_fl, 1, cf_upper_S_fl, 1);
		mkl<T>::scal(Constants::DIM, scale, cf_upper_S_fl, 1);


		// upper_S2 = upper_rotational_stiffness(2) * upper_angle2 * upper_normal2;
		scale = (this->upper_rotational_stiffness[1]) * (*cf_upper_angle_fr);
		mkl<T>::copy(Constants::DIM, cf_upper_normal_fr, 1, cf_upper_S_fr, 1);
		mkl<T>::scal(Constants::DIM, scale, cf_upper_S_fr, 1);

		// upper_S3 = upper_rotational_stiffness(3) * upper_angle3 * upper_normal3;
		scale = (this->upper_rotational_stiffness[2]) * (*cf_upper_angle_rl);
		mkl<T>::copy(Constants::DIM, cf_upper_normal_rl, 1, cf_upper_S_rl, 1);
		mkl<T>::scal(Constants::DIM, scale, cf_upper_S_rl, 1);

		// upper_S4 = upper_rotational_stiffness(4) * upper_angle4 * upper_normal4;
		scale = (this->upper_rotational_stiffness[3]) * (*cf_upper_angle_rr);
		mkl<T>::copy(Constants::DIM, cf_upper_normal_rr, 1, cf_upper_S_rr, 1);
		mkl<T>::scal(Constants::DIM, scale, cf_upper_S_rr, 1);

		// lower_S1 = lower_rotational_stiffness(1) * lower_angle1 * lower_normal1
		scale = (this->lower_rotational_stiffness[0]) * (*cf_lower_angle_fl);
		mkl<T>::copy(Constants::DIM, cf_lower_normal_fl, 1, cf_lower_S_fl, 1);
		mkl<T>::scal(Constants::DIM, scale, cf_lower_S_fl, 1);

		// lower_S2 = lower_rotational_stiffness(2) * lower_angle2 * lower_normal2
		scale = (this->lower_rotational_stiffness[1]) * (*cf_lower_angle_fr);
		mkl<T>::copy(Constants::DIM, cf_lower_normal_fr, 1, cf_lower_S_fr, 1);
		mkl<T>::scal(Constants::DIM, scale, cf_lower_S_fr, 1);

		// lower_S3 = lower_rotational_stiffness(3) * lower_angle3 * lower_normal3
		scale = (this->lower_rotational_stiffness[2]) * (*cf_lower_angle_rl);
		mkl<T>::copy(Constants::DIM, cf_lower_normal_rl, 1, cf_lower_S_rl, 1);
		mkl<T>::scal(Constants::DIM, scale, cf_lower_S_rl, 1);

		// lower_S_rr = lower_rotational_stiffness(4) * lower_angle4 * lower_normal4
		scale = (this->lower_rotational_stiffness[3]) * (*cf_lower_angle_rr);
		mkl<T>::copy(Constants::DIM, cf_lower_normal_rr, 1, cf_lower_S_rr, 1);
		mkl<T>::scal(Constants::DIM, scale, cf_lower_S_rr, 1);

		/*
		calculate the effect of the rotational springs (in global basis)
		lower_rot_force1 = -cross( lower_S1, r_low1) / (r_low1'*r_low1);
		lower_rot_force2 = -cross( lower_S2, r_low2) / (r_low2'*r_low2);
		lower_rot_force3 = -cross( lower_S3, r_low3) / (r_low3'*r_low3);
		lower_rot_force4 = -cross( lower_S4, r_low4) / (r_low4'*r_low4);

		upper_rot_force1 = -cross( upper_S1, r_up1) / (r_up1'*r_up1);
		upper_rot_force2 = -cross( upper_S2, r_up2) / (r_up2'*r_up2);
		upper_rot_force3 = -cross( upper_S3, r_up3) / (r_up3'*r_up3);
		upper_rot_force4 = -cross( upper_S4, r_up4) / (r_up4'*r_up4);

		car_rot_force1 = -cross( lower_S1, r_up1) / (r_up1'*r_up1);
		car_rot_force2 = -cross( lower_S2, r_up2) / (r_up2'*r_up2);
		car_rot_force3 = -cross( lower_S3, r_up3) / (r_up3'*r_up3);
		car_rot_force4 = -cross( lower_S4, r_up4) / (r_up4'*r_up4);

		sum_car_force1 = car_rot_force1 - upper_force1 - upper_dampf1 - upper_rot_force1;
		sum_car_force2 = car_rot_force2 - upper_force2 - upper_dampf2 - upper_rot_force2;
		sum_car_force3 = car_rot_force3 - upper_force3 - upper_dampf3 - upper_rot_force3;
		sum_car_force4 = car_rot_force4 - upper_force4 - upper_dampf4 - upper_rot_force4;
		*/

		T scale_u_fl, scale_u_fr, scale_u_rl, scale_u_rr;
		// lower_rot_force1 = -cross( lower_S1, r_low1) / (r_low1'*r_low1);
		scale = -1.0 / mkl<T>::dot(mkl_DIM, cf_r_low_fl, mkl_incx, cf_r_low_fl, mkl_incy);
		//scale = -1.0 / MathLibrary::dot_product<T>(r_low1, r_low1, Constants::DIM);
		MathLibrary::crossProduct(cf_lower_S_fl, cf_r_low_fl, cf_lower_rot_force_fl);
		mkl<T>::scal(Constants::DIM, scale, cf_lower_rot_force_fl, 1);

		// lower_rot_force2 = -cross( lower_S2, r_low2) / (r_low2'*r_low2);
		scale = -1.0 / mkl<T>::dot(mkl_DIM, cf_r_low_fr, mkl_incx, cf_r_low_fr, mkl_incy);
		//scale = -1.0 / MathLibrary::dot_product<T>(r_low_fr, r_low_fr, Constants::DIM);
		MathLibrary::crossProduct(cf_lower_S_fr, cf_r_low_fr, cf_lower_rot_force_fr);
		mkl<T>::scal(Constants::DIM, scale, cf_lower_rot_force_fr, 1);

		// lower_rot_force3 = -cross( lower_S3, r_low3) / (r_low3'*r_low3);
		scale = -1.0 / mkl<T>::dot(mkl_DIM, cf_r_low_rl, mkl_incx, cf_r_low_rl, mkl_incy);
		//scale = -1.0 / MathLibrary::dot_product<T>(r_low3, r_low3, Constants::DIM);
		MathLibrary::crossProduct(cf_lower_S_rl, cf_r_low_rl, cf_lower_rot_force_rl);
		mkl<T>::scal(Constants::DIM, scale, cf_lower_rot_force_rl, 1);

		// lower_rot_force4 = -cross( lower_S4, r_low4) / (r_low4'*r_low4);
		scale = -1.0 / mkl<T>::dot(mkl_DIM, cf_r_low_rr, mkl_incx, cf_r_low_rr, mkl_incx);
		//scale = -1.0 / MathLibrary::dot_product<T>(r_low4, r_low4, Constants::DIM);
		MathLibrary::crossProduct(cf_lower_S_rr, cf_r_low_rr, cf_lower_rot_force_rr);
		mkl<T>::scal(Constants::DIM, scale, cf_lower_rot_force_rr, 1);

		// upper_rot_force1 = -cross( upper_S1, r_up1) / (r_up1'*r_up1);
		scale_u_fl = -1.0 / mkl<T>::dot(mkl_DIM, cf_r_up_fl, mkl_incx, cf_r_up_fl, mkl_incy);
		//scale_u1 = -1.0 / MathLibrary::dot_product<T>(r_up1, r_up1, Constants::DIM);
		MathLibrary::crossProduct(cf_upper_S_fl, cf_r_up_fl, cf_upper_rot_force_fl);
		mkl<T>::scal(Constants::DIM, scale_u_fl, cf_upper_rot_force_fl, 1);

		// upper_rot_force_fr = -cross( upper_S_fr, r_up_fr) / (r_up_fr'*r_up_fr);
		scale_u_fr = -1.0 / mkl<T>::dot(mkl_DIM, cf_r_up_fr, mkl_incx, cf_r_up_fr, mkl_incy);
		//scale_u_fr = -1.0 / MathLibrary::dot_product<T>(r_up_fr, r_up_fr, Constants::DIM);
		MathLibrary::crossProduct(cf_upper_S_fr, cf_r_up_fr, cf_upper_rot_force_fr);
		mkl<T>::scal(Constants::DIM, scale_u_fr, cf_upper_rot_force_fr, 1);

		// upper_rot_force3 = -cross( upper_S3, r_up3) / (r_up3'*r_up3);
		scale_u_rl = -1.0 / mkl<T>::dot(mkl_DIM, cf_r_up_rl, mkl_incx, cf_r_up_rl, mkl_incy);
		//scale_u3 = -1.0 / MathLibrary::dot_product<T>(r_up3, r_up3, Constants::DIM);
		MathLibrary::crossProduct(cf_upper_S_rl, cf_r_up_rl, cf_upper_rot_force_rl);
		mkl<T>::scal(Constants::DIM, scale_u_rl, cf_upper_rot_force_rl, 1);

		// upper_rot_force4 = -cross( upper_S4, r_up4) / (r_up4'*r_up4);
		scale_u_rr = -1.0 / mkl<T>::dot(mkl_DIM, cf_r_up_rr, mkl_incx, cf_r_up_rr, mkl_incy);
		//scale_u4 = -1.0 / MathLibrary::dot_product<T>(r_up4, r_up4, Constants::DIM);
		MathLibrary::crossProduct(cf_upper_S_rr, cf_r_up_rr, cf_upper_rot_force_rr);
		mkl<T>::scal(Constants::DIM, scale_u_rr, cf_upper_rot_force_rr, 1);

		// car_rot_force1 = -cross( lower_S1, r_up1) / (r_up1'*r_up1);
		MathLibrary::crossProduct(cf_lower_S_fl, cf_r_up_fl, cf_car_rot_force_fl);
		mkl<T>::scal(Constants::DIM, scale_u_fl, cf_car_rot_force_fl, 1);

		// car_rot_force_fr = -cross( lower_S2, r_up2) / (r_up2'*r_up2);
		MathLibrary::crossProduct(cf_lower_S_fr, cf_r_up_fr, cf_car_rot_force_fr);
		mkl<T>::scal(Constants::DIM, scale_u_fr, cf_car_rot_force_fr, 1);

		// car_rot_force3 = -cross( lower_S3, r_up3) / (r_up3'*r_up3);
		MathLibrary::crossProduct(cf_lower_S_rl, cf_r_up_rl, cf_car_rot_force_rl);
		mkl<T>::scal(Constants::DIM, scale_u_rl, cf_car_rot_force_rl, 1);

		// car_rot_force4 = -cross( lower_S4, r_up4) / (r_up4'*r_up4);
		MathLibrary::crossProduct(cf_lower_S_rr, cf_r_up_rr, cf_car_rot_force_rr);
		mkl<T>::scal(Constants::DIM, scale_u_rr, cf_car_rot_force_rr, 1);

		// sum_car_force1 = car_rot_force1 - upper_force1 - upper_dampf1 - upper_rot_force1;
		mkl<T>::copy(Constants::DIM, cf_car_rot_force_fl, 1, cf_sum_car_force_fl, 1);
		mkl<T>::axpy(Constants::DIM, -1.0, cf_upper_force_fl, 1, cf_sum_car_force_fl, 1);
		mkl<T>::axpy(Constants::DIM, -1.0, cf_upper_dampf_fl, 1, cf_sum_car_force_fl, 1);
		mkl<T>::axpy(Constants::DIM, -1.0, cf_upper_rot_force_fl, 1, cf_sum_car_force_fl, 1);

		// sum_car_force_fr = car_rot_force_fr - upper_force_fr - upper_dampf_fr - upper_rot_force_fr;
		mkl<T>::copy(Constants::DIM, cf_car_rot_force_fr, 1, cf_sum_car_force_fr, 1);
		mkl<T>::axpy(Constants::DIM, -1.0, cf_upper_force_fr, 1, cf_sum_car_force_fr, 1);
		mkl<T>::axpy(Constants::DIM, -1.0, cf_upper_dampf_fr, 1, cf_sum_car_force_fr, 1);
		mkl<T>::axpy(Constants::DIM, -1.0, cf_upper_rot_force_fr, 1, cf_sum_car_force_fr, 1);

		// sum_car_force3 = car_rot_force3 - upper_force3 - upper_dampf3 - upper_rot_force3;
		mkl<T>::copy(Constants::DIM, cf_car_rot_force_rl, 1, cf_sum_car_force_rl, 1);
		mkl<T>::axpy(Constants::DIM, -1.0, cf_upper_force_rl, 1, cf_sum_car_force_rl, 1);
		mkl<T>::axpy(Constants::DIM, -1.0, cf_upper_dampf_rl, 1, cf_sum_car_force_rl, 1);
		mkl<T>::axpy(Constants::DIM, -1.0, cf_upper_rot_force_rl, 1, cf_sum_car_force_rl, 1);

		// sum_car_force4 = car_rot_force4 - upper_force4 - upper_dampf4 - upper_rot_force4;
		mkl<T>::copy(Constants::DIM, cf_car_rot_force_rr, 1, cf_sum_car_force_rr, 1);
		mkl<T>::axpy(Constants::DIM, -1.0, cf_upper_force_rr, 1, cf_sum_car_force_rr, 1);
		mkl<T>::axpy(Constants::DIM, -1.0, cf_upper_dampf_rr, 1, cf_sum_car_force_rr, 1);
		mkl<T>::axpy(Constants::DIM, -1.0, cf_upper_rot_force_rr, 1, cf_sum_car_force_rr, 1);

		/*
		get external forces
		for the legs
		local_FW1 = [0; FW(1); 0];      %in global basis
		local_FW_fr = [0; FW(_fr); 0];
		local_FW3 = [0; FW(3); 0];
		local_FW4 = [0; FW(4); 0];

		local_FT1 = [0; FT(1); 0];      %in global basis
		local_FT_fr = [0; FT(2); 0];
		local_FT3 = [0; FT(3); 0];
		local_FT4 = [0; FT(4); 0];
		*/

		// Calculate the sum of all forces in the tyre (and keep it in local_FR)

		// local_FR1 = lower_force1 + lower_dampf1 + local_FT1 + local_FR1 + lower_rot_force1; ...          %vt1_dot
		mkl<T>::copy(Constants::DIM, cf_lower_force_fl, 1, cf_local_FR_fl, 1);
		mkl<T>::axpy(Constants::DIM, 1.0, cf_lower_dampf_fl, 1, cf_local_FR_fl, 1);
		//	mkl<T>::axpy(Constants::DIM, 1.0, local_FT1, 1, cf_local_FR1, 1); not implemented since not used now
		mkl<T>::axpy(Constants::DIM, 1.0, cf_local_FR_fl, 1, cf_local_FR_fl, 1);
		mkl<T>::axpy(Constants::DIM, 1.0, cf_lower_rot_force_fl, 1, cf_local_FR_fl, 1);

		// local_FR2 = lower_force2 + lower_dampf2 + local_FT2 + local_FR2 + lower_rot_force2; ...          %vt2_dot
		mkl<T>::copy(Constants::DIM, cf_lower_force_fr, 1, cf_local_FR_fr, 1);
		mkl<T>::axpy(Constants::DIM, 1.0, cf_lower_dampf_fr, 1, cf_local_FR_fr, 1);
		//	mkl<T>::axpy(Constants::DIM, 1.0, local_FT2, 1, cf_local_FR2, 1); not implemented since not used now
		mkl<T>::axpy(Constants::DIM, 1.0, cf_local_FR_fr, 1, cf_local_FR_fr, 1);
		mkl<T>::axpy(Constants::DIM, 1.0, cf_lower_rot_force_fr, 1, cf_local_FR_fr, 1);

		// local_FR3 = lower_force3 + lower_dampf3 + local_FT3 + local_FR3 + lower_rot_force3; ...          %vt3_dot
		mkl<T>::copy(Constants::DIM, cf_lower_force_rl, 1, cf_local_FR_rl, 1);
		mkl<T>::axpy(Constants::DIM, 1.0, cf_lower_dampf_rl, 1, cf_local_FR_rl, 1);
		//	mkl<T>::axpy(Constants::DIM, 1.0, local_FT3, 1, cf_local_FR3, 1); not implemented since not used now
		mkl<T>::axpy(Constants::DIM, 1.0, cf_local_FR_rl, 1, cf_local_FR_rl, 1);
		mkl<T>::axpy(Constants::DIM, 1.0, cf_lower_rot_force_rl, 1, cf_local_FR_rl, 1);

		// local_FR4 = lower_force4 + lower_dampf4 + local_FT4 + local_FR4 + lower_rot_force4];             %vt4_dot
		mkl<T>::copy(Constants::DIM, cf_lower_force_rr, 1, cf_local_FR_rr, 1);
		mkl<T>::axpy(Constants::DIM, 1.0, cf_lower_dampf_rr, 1, cf_local_FR_rr, 1);
		//	mkl<T>::axpy(Constants::DIM, 1.0, local_FT4, 1, cf_local_FR4, 1); not implemented since not used now
		mkl<T>::axpy(Constants::DIM, 1.0, cf_local_FR_rr, 1, cf_local_FR_rr, 1);
		mkl<T>::axpy(Constants::DIM, 1.0, cf_lower_rot_force_rr, 1, cf_local_FR_rr, 1);


		if (boundary_conditions == FIXED) {
			get_fixed_road_force(cf_local_FR_fl, cf_local_FR_fr, cf_local_FR_rl, cf_local_FR_rr);
		}
		else if (boundary_conditions == NONFIXED) {
			get_nonfixed_road_force(cf_local_FR_fl, cf_local_FR_fr, cf_local_FR_rl, cf_local_FR_rr);
		}
		else if (boundary_conditions == CIRCULAR) {
			get_circular_road_force(cf_local_FR_fl, vt_fl_, mass_tyre_rr, pt_fl_);
			get_circular_road_force(cf_local_FR_fr, vt_fr_, mass_tyre_rl, pt_fr_);
			get_circular_road_force(cf_local_FR_rl, vt_rl_, mass_tyre_fl, pt_rl_);
			get_circular_road_force(cf_local_FR_rr, vt_rr_, mass_tyre_fr, pt_rr_);
		}
		
		/*
		get H=I*w
		Hc = A(1:3, 1:3) * wc;           %in local car basis

		external torque on the car body (use later for rotational damping)
		Tc = zeros(3,1);     %in local car basis

		sum of all torques induced by forces
		for the car body
		sum_torque_spring_car = r1_tilda * (C_cN * sum_car_force1) + ...     % from the elongational springs
								r2_tilda * (C_cN * sum_car_force2) + ...
								r3_tilda * (C_cN * sum_car_force3) + ...
								r_rr_tilda * (C_cN * sum_car_force4) + ...
							   -C_cN * upper_S1 - C_cN * upper_S2 - C_cN * upper_S3 - C_cN * upper_S4 + ...     % ??from the rotational spring
							   -get_tilda(wc) * Hc + Tc;                       % from angular momentum and external torques

		*/


		// Hc = A(1:3, 1:3) * wc;
		mkl<T>::gemv(CblasRowMajor, CblasNoTrans, Constants::DIM, Constants::DIM, 1.0, this->Ic, Constants::DIM, wc_, 1, 0.0, cf_Hc, 1);
		MathLibrary::get_tilda<T>(wc_, cf_wc_tilda);
		mkl<T>::copy(Constants::DIM, cf_Hc, 1, cf_temp, 1);
		mkl<T>::gemv(CblasRowMajor, CblasNoTrans, Constants::DIM, Constants::DIM, -1.0, cf_wc_tilda, Constants::DIM, cf_temp, 1, 0.0, cf_Hc, 1);
		mkl<T>::copy(Constants::DIM, cf_Hc, 1, cf_sum_torque_spring_car, 1);
		mkl<T>::axpy(Constants::DIM, 1.0, cf_Tc, 1, cf_sum_torque_spring_car, 1);
		mkl<T>::gemv(CblasRowMajor, CblasTrans, Constants::DIM, Constants::DIM, -1.0, cf_C_cN, Constants::DIM, cf_upper_S_rr, 1, 1.0, cf_sum_torque_spring_car, 1);
		mkl<T>::gemv(CblasRowMajor, CblasTrans, Constants::DIM, Constants::DIM, -1.0, cf_C_cN, Constants::DIM, cf_upper_S_rl, 1, 1.0, cf_sum_torque_spring_car, 1);
		mkl<T>::gemv(CblasRowMajor, CblasTrans, Constants::DIM, Constants::DIM, -1.0, cf_C_cN, Constants::DIM, cf_upper_S_fr, 1, 1.0, cf_sum_torque_spring_car, 1);
		mkl<T>::gemv(CblasRowMajor, CblasTrans, Constants::DIM, Constants::DIM, -1.0, cf_C_cN, Constants::DIM, cf_upper_S_fl, 1, 1.0, cf_sum_torque_spring_car, 1);
		mkl<T>::gemv(CblasRowMajor, CblasTrans, Constants::DIM, Constants::DIM, 1.0, cf_C_cN, Constants::DIM, cf_sum_car_force_fl, 1, 0.0, cf_temp, 1);
		mkl<T>::gemv(CblasRowMajor, CblasNoTrans, Constants::DIM, Constants::DIM, 1.0, this->r_fl_tilda, Constants::DIM, cf_temp, 1, 1.0, cf_sum_torque_spring_car, 1);
		mkl<T>::gemv(CblasRowMajor, CblasTrans, Constants::DIM, Constants::DIM, 1.0, cf_C_cN, Constants::DIM, cf_sum_car_force_fr, 1, 0.0, cf_temp, 1);
		mkl<T>::gemv(CblasRowMajor, CblasNoTrans, Constants::DIM, Constants::DIM, 1.0, this->r_fr_tilda, Constants::DIM, cf_temp, 1, 1.0, cf_sum_torque_spring_car, 1);
		mkl<T>::gemv(CblasRowMajor, CblasTrans, Constants::DIM, Constants::DIM, 1.0, cf_C_cN, Constants::DIM, cf_sum_car_force_rl, 1, 0.0, cf_temp, 1);
		mkl<T>::gemv(CblasRowMajor, CblasNoTrans, Constants::DIM, Constants::DIM, 1.0, this->r_rl_tilda, Constants::DIM, cf_temp, 1, 1.0, cf_sum_torque_spring_car, 1);
		mkl<T>::gemv(CblasRowMajor, CblasTrans, Constants::DIM, Constants::DIM, 1.0, cf_C_cN, Constants::DIM, cf_sum_car_force_rr, 1, 0.0, cf_temp, 1);
		mkl<T>::gemv(CblasRowMajor, CblasNoTrans, Constants::DIM, Constants::DIM, 1.0, this->r_rr_tilda, Constants::DIM, cf_temp, 1, 1.0, cf_sum_torque_spring_car, 1);

		/*
			///////////////////////////////////////////////////////////////
			///////////                 Solve                  ////////////
			///////////////////////////////////////////////////////////////
			b =  [sum_torque_spring_car;...                                                           %w_dot_c
				 FC + sum_car_force1 + sum_car_force2 + sum_car_force3 + sum_car_force4; ...          %vc_dot
				 upper_force1 - lower_force1 + upper_dampf1 - lower_dampf1 + local_FW1 + ...
						upper_rot_force1 - car_rot_force1 - lower_rot_force1; ...                     %vw1_dot
				 upper_force2 - lower_force2 + upper_dampf2 - lower_dampf2 + local_FW2 + ...
						upper_rot_force2 - car_rot_force2 - lower_rot_force2; ...                     %vw2_dot
				 upper_force3 - lower_force3 + upper_dampf3 - lower_dampf3 + local_FW3 + ...
						upper_rot_force3 - car_rot_force3 - lower_rot_force3; ...                     %vw3_dot
				 upper_force4 - lower_force4 + upper_dampf4 - lower_dampf4 + local_FW4 + ...
						upper_rot_force4 - car_rot_force4 - lower_rot_force4; ...                     %vw4_dot
				 lower_force1 + lower_dampf1 + local_FT1 + local_FR1 + lower_rot_force1; ...          %vt1_dot
				 lower_force2 + lower_dampf2 + local_FT2 + local_FR2 + lower_rot_force2; ...          %vt2_dot
				 lower_force3 + lower_dampf3 + local_FT3 + local_FR3 + lower_rot_force3; ...          %vt3_dot
				 lower_force4 + lower_dampf4 + local_FT4 + local_FR4 + lower_rot_force4];             %vt4_dot

		*/
		// FC + sum_car_force1 + sum_car_force2 + sum_car_force3 + sum_car_force4; ...          %vc_dot  
		T* brem_start = cf_b_rem;
		mkl<T>::copy(Constants::DIM, this->FC, 1, brem_start, 1);
		mkl<T>::axpy(Constants::DIM, 1.0, cf_sum_car_force_fl, 1, brem_start, 1);
		mkl<T>::axpy(Constants::DIM, 1.0, cf_sum_car_force_fr, 1, brem_start, 1);
		mkl<T>::axpy(Constants::DIM, 1.0, cf_sum_car_force_rl, 1, brem_start, 1);
		mkl<T>::axpy(Constants::DIM, 1.0, cf_sum_car_force_rr, 1, brem_start, 1);

		//  upper_force1 - lower_force1 + upper_dampf1 - lower_dampf1 + local_FW1 + upper_rot_force1 - car_rot_force1 - lower_rot_force1;
		brem_start += Constants::DIM;
		mkl<T>::copy(Constants::DIM, cf_upper_force_fl, 1, brem_start, 1);
		mkl<T>::axpy(Constants::DIM, -1.0, cf_lower_force_fl, 1, brem_start, 1);
		mkl<T>::axpy(Constants::DIM, 1.0, cf_upper_dampf_fl, 1, brem_start, 1);
		mkl<T>::axpy(Constants::DIM, -1.0, cf_lower_dampf_fl, 1, brem_start, 1);
		//mkl<T>::axpy(Constants::DIM, 1.0, local_FW1, 1, brem_start, 1);
		mkl<T>::axpy(Constants::DIM, 1.0, cf_upper_rot_force_fl, 1, brem_start, 1);
		mkl<T>::axpy(Constants::DIM, -1.0, cf_car_rot_force_fl, 1, brem_start, 1);
		mkl<T>::axpy(Constants::DIM, -1.0, cf_lower_rot_force_fl, 1, brem_start, 1);

		// upper_force2 - lower_force2 + upper_dampf2 - lower_dampf2 + local_FW2 + upper_rot_force2 - car_rot_force2 - lower_rot_force2;
		brem_start += Constants::DIM;
		mkl<T>::copy(Constants::DIM, cf_upper_force_fr, 1, brem_start, 1);
		mkl<T>::axpy(Constants::DIM, -1.0, cf_lower_force_fr, 1, brem_start, 1);
		mkl<T>::axpy(Constants::DIM, 1.0, cf_upper_dampf_fr, 1, brem_start, 1);
		mkl<T>::axpy(Constants::DIM, -1.0, cf_lower_dampf_fr, 1, brem_start, 1);
		//mkl<T>::axpy(Constants::DIM, 1.0, local_FW2, 1, brem_start, 1);
		mkl<T>::axpy(Constants::DIM, 1.0, cf_upper_rot_force_fr, 1, brem_start, 1);
		mkl<T>::axpy(Constants::DIM, -1.0, cf_car_rot_force_fr, 1, brem_start, 1);
		mkl<T>::axpy(Constants::DIM, -1.0, cf_lower_rot_force_fr, 1, brem_start, 1);

		// upper_force3 - lower_force3 + upper_dampf3 - lower_dampf3 + local_FW3 + upper_rot_force3 - car_rot_force3 - lower_rot_force3;
		brem_start += Constants::DIM;
		mkl<T>::copy(Constants::DIM, cf_upper_force_rl, 1, brem_start, 1);
		mkl<T>::axpy(Constants::DIM, -1.0, cf_lower_force_rl, 1, brem_start, 1);
		mkl<T>::axpy(Constants::DIM, 1.0, cf_upper_dampf_rl, 1, brem_start, 1);
		mkl<T>::axpy(Constants::DIM, -1.0, cf_lower_dampf_rl, 1, brem_start, 1);
		//mkl<T>::axpy(Constants::DIM, 1.0, local_FW_rl, 1, brem_start, 1);
		mkl<T>::axpy(Constants::DIM, 1.0, cf_upper_rot_force_rl, 1, brem_start, 1);
		mkl<T>::axpy(Constants::DIM, -1.0, cf_car_rot_force_rl, 1, brem_start, 1);
		mkl<T>::axpy(Constants::DIM, -1.0, cf_lower_rot_force_rl, 1, brem_start, 1);

		//  upper_force4 - lower_force4 + upper_dampf4 - lower_dampf4 + local_FW4 + upper_rot_force4 - car_rot_force4 - lower_rot_force4;
		brem_start += Constants::DIM;
		mkl<T>::copy(Constants::DIM, cf_upper_force_rr, 1, brem_start, 1);
		mkl<T>::axpy(Constants::DIM, -1.0, cf_lower_force_rr, 1, brem_start, 1);
		mkl<T>::axpy(Constants::DIM, 1.0, cf_upper_dampf_rr, 1, brem_start, 1);
		mkl<T>::axpy(Constants::DIM, -1.0, cf_lower_dampf_rr, 1, brem_start, 1);
		//mkl<T>::axpy(Constants::DIM, 1.0, local_FW4, 1, brem_start, 1);
		mkl<T>::axpy(Constants::DIM, 1.0, cf_upper_rot_force_rr, 1, brem_start, 1);
		mkl<T>::axpy(Constants::DIM, -1.0, cf_car_rot_force_rr, 1, brem_start, 1);
		mkl<T>::axpy(Constants::DIM, -1.0, cf_lower_rot_force_rr, 1, brem_start, 1);

		// local_FR1; ...					  %vt1_dot
		brem_start += Constants::DIM;
		mkl<T>::copy(Constants::DIM, cf_local_FR_fl, 1, brem_start, 1);

		// local_FR2; ...			 	      %vt2_dot
		brem_start += Constants::DIM;
		mkl<T>::copy(Constants::DIM, cf_local_FR_fr, 1, brem_start, 1);

		// local_FR3; ...			          %vt3_dot
		brem_start += Constants::DIM;
		mkl<T>::copy(Constants::DIM, cf_local_FR_rl, 1, brem_start, 1);
	
		// local_FR4];			              %vt4_dot
		brem_start += Constants::DIM;
		mkl<T>::copy(Constants::DIM, cf_local_FR_rr, 1, brem_start, 1);

		mkl<T>::potrs(LAPACK_ROW_MAJOR, 'L', Constants::DIM, 1, A_Ic, Constants::DIM, cf_sum_torque_spring_car, 1);
		mkl<T>::copy(Constants::DIM, cf_sum_torque_spring_car, 1, f_, 1);
		T* start_next = f_ + Constants::DIM;

		MathLibrary::vector_elem_wise_product<T>(A_rem, cf_b_rem, start_next, 9 * Constants::DIM);
		start_next += 9 * Constants::DIM;

		/*
		get the derivative of the altitude (expressed in quaternions) from the angular velocities
		Qc = 0.5 * [qc(4) -qc(3) qc(2); qc(3) qc(4) -qc(1); -qc(2) qc(1) qc(4); -qc(1) -qc(2) -qc(3)];
		qc_dot = Qc * wc;
		*/

		cf_Qc[0] = 0.5 * qc_[3];
		cf_Qc[1] = -0.5 * qc_[2];
		cf_Qc[2] = 0.5 * qc_[1];
		cf_Qc[3] = 0.5 * qc_[2];
		cf_Qc[4] = 0.5 * qc_[3];
		cf_Qc[5] = -0.5 * qc_[0];
		cf_Qc[6] = -0.5 * qc_[1];
		cf_Qc[7] = 0.5 * qc_[0];
		cf_Qc[8] = 0.5 * qc_[3];
		cf_Qc[9] = -0.5 * qc_[0];
		cf_Qc[10] = -0.5 * qc_[1];
		cf_Qc[11] = -0.5 * qc_[2];
		mkl<T>::gemv(CblasRowMajor, CblasNoTrans, Constants::NUM_LEGS, Constants::DIM, 1.0, cf_Qc, Constants::DIM, wc_, 1, 0.0, cf_qc_dot, 1);

		mkl<T>::copy(Constants::NUM_LEGS, cf_qc_dot, 1, start_next, 1);
		start_next += Constants::NUM_LEGS;
		// add vc to f vector
		mkl<T>::copy(Constants::DIM, vc_, 1, start_next, 1);
		start_next += Constants::DIM;

		// add vw1 to f vector

		mkl<T>::copy(Constants::DIM, vw_fl_, 1, start_next, 1);
		start_next += Constants::DIM;

		// add vw2 to f vector
		mkl<T>::copy(Constants::DIM, vw_fr_, 1, start_next, 1);
		start_next += Constants::DIM;

		// add vw3 to f vector
		mkl<T>::copy(Constants::DIM, vw_rl_, 1, start_next, 1);
		start_next += Constants::DIM;

		// add vw_rr to f vector
		mkl<T>::copy(Constants::DIM, vw_rr_, 1, start_next, 1);
		start_next += Constants::DIM;

		// add vt1 to f vector
		mkl<T>::copy(Constants::DIM, vt_fl_, 1, start_next, 1);
		start_next += Constants::DIM;

		// add vt2 to f vector
		mkl<T>::copy(Constants::DIM, vt_fr_, 1, start_next, 1);
		start_next += Constants::DIM;

		// add vt3 to f vector
		mkl<T>::copy(Constants::DIM, vt_rl_, 1, start_next, 1);
		start_next += Constants::DIM;

		// add vt_rr to f vector
		mkl<T>::copy(Constants::DIM, vt_rr_, 1, start_next, 1);
		start_next += Constants::DIM;


	}

	/*
	Initializes the time iteration and handles the numerical scheme
	*/
	void solve(T* solution_vector) {
		// From the formulation we have 61 dimensions in the solution vector
		size_t solution_size = (this->num_iter + 1) * this->solution_dim;
		T* complete_vector = (T*)mkl_calloc(solution_size, sizeof(T), Constants::ALIGNMENT);
		x_vector = (T*)mkl_calloc(solution_dim, sizeof(T), Constants::ALIGNMENT);
		MathLibrary::get_tilda<T>(r_fl, r_fl_tilda);
		MathLibrary::get_tilda<T>(r_fr, r_fr_tilda);
		MathLibrary::get_tilda<T>(r_rl, r_rl_tilda);
		MathLibrary::get_tilda<T>(r_rr, r_rr_tilda);
		/*
		Preparing x_vector in the form of
		x_vector = [wc; ...     % 3 1:3
				vc; ...     % 3 4:6
				vw1; ...    % 3 7:9
				vw2; ...    % 3 10:12
				vw3; ...    % 3 13:15
				vw4; ...    % 3 16:18
				vt1; ...    % 3 19:21
				vt2; ...    % 3 22:24
				vt3; ...    % 3 25:27
				vt4; ...    % 3 28:30
				qc; ...     % 4 31:34
				pcc; ...    % 3 35:37
				pw1; ...    % 3 38:40
				pw2; ...    % 3 41:43
				pw3; ...    % 3 44:46
				pw4; ...    % 3 47:49
				pt1; ...    % 3 50:52
				pt2; ...    % 3 53:55
				pt3; ...    % 3 56:58
				pt4];       % 3 59:61
		*/
		size_t i, j;
		i = 10;
		j = 0;
		// qc
		T* start_orient = x_vector + i * (Constants::DIM) + j * (Constants::NUM_LEGS);
		mkl<T>::copy(Constants::NUM_LEGS, initial_orientation, 1, x_vector + i * (Constants::DIM) + j * (Constants::NUM_LEGS), 1);
		j++;

		// pcc
		mkl<T>::copy(Constants::DIM, pcc, 1, x_vector + i * (Constants::DIM) + j * (Constants::NUM_LEGS), 1);
		i++;
		// pw1
		T* pw_fl = x_vector + i * (Constants::DIM) + j * (Constants::NUM_LEGS);
		i++;
		// pw2
		T* pw_fr = x_vector + i * (Constants::DIM) + j * (Constants::NUM_LEGS);
		i++;
		// pw3
		T* pw_rl = x_vector + i * (Constants::DIM) + j * (Constants::NUM_LEGS);
		i++;
		// pw_rr
		T* pw_rr = x_vector + i * (Constants::DIM) + j * (Constants::NUM_LEGS);
		i++;
		// pt1
		T* pt_fl = x_vector + i * (Constants::DIM) + j * (Constants::NUM_LEGS);
		i++;
		// pt2
		T* pt_fr = x_vector + i * (Constants::DIM) + j * (Constants::NUM_LEGS);
		i++;
		// pt3
		T* pt_rl = x_vector + i * (Constants::DIM) + j * (Constants::NUM_LEGS);
		i++;
		// pt_rr
		T* pt_rr = x_vector + i * (Constants::DIM) + j * (Constants::NUM_LEGS);


		get_initial_length(start_orient, r_fl, r_fr, r_rl, r_rr, pcc, initial_upper_spring_length, initial_lower_spring_length, pw_fl, pw_fr, pw_rl, pw_rr, pt_fl, pt_fr, pt_rl, pt_rr);

		// overwrites the initial velocity values
		if (boundary_conditions == CIRCULAR) 
			circular_path_initialization(vc, vw_fl, vw_fr, vw_rl, vw_rr, vt_fl, vt_fr, vt_rl, vt_rr, initial_angular_velocity, pcc, pt_fl, pt_fr, pt_rl, pt_rr, radius_circular_path);

		i = 0;
		j = 0;
		// wc
		mkl<T>::copy(Constants::DIM, initial_angular_velocity, 1, x_vector + i * (Constants::DIM) + j * (Constants::NUM_LEGS), 1);
		i++;
		// vc
		mkl<T>::copy(Constants::DIM, vc, 1, x_vector + i * (Constants::DIM) + j * (Constants::NUM_LEGS), 1);
		i++;
		// vw1
		mkl<T>::copy(Constants::DIM, vw_fl, 1, x_vector + i * (Constants::DIM) + j * (Constants::NUM_LEGS), 1);
		i++;
		// vw2
		mkl<T>::copy(Constants::DIM, vw_fr, 1, x_vector + i * (Constants::DIM) + j * (Constants::NUM_LEGS), 1);
		i++;
		// vw_rl
		mkl<T>::copy(Constants::DIM, vw_rl, 1, x_vector + i * (Constants::DIM) + j * (Constants::NUM_LEGS), 1);
		i++;
		// vw4
		mkl<T>::copy(Constants::DIM, vw_rr, 1, x_vector + i * (Constants::DIM) + j * (Constants::NUM_LEGS), 1);
		i++;
		// vt1
		mkl<T>::copy(Constants::DIM, vt_fl, 1, x_vector + i * (Constants::DIM) + j * (Constants::NUM_LEGS), 1);
		i++;
		// vt2
		mkl<T>::copy(Constants::DIM, vt_fr, 1, x_vector + i * (Constants::DIM) + j * (Constants::NUM_LEGS), 1);
		i++;
		// vt_rl
		mkl<T>::copy(Constants::DIM, vt_rl, 1, x_vector + i * (Constants::DIM) + j * (Constants::NUM_LEGS), 1);
		i++;
		// vt4
		mkl<T>::copy(Constants::DIM, vt_rr, 1, x_vector + i * (Constants::DIM) + j * (Constants::NUM_LEGS), 1);
		i++;

		j = 0;
		i = 0;
		while (i < Constants::DIM) {
			A_rem[j] = 1.0 / this->mass;
			i++;
			j++;
		}
		i = 0;
		while (i < Constants::DIM) {
			A_rem[j] = 1.0 / mass_wheel[0];
			i++;
			j++;
		}
		i = 0;
		while (i < Constants::DIM) {
			A_rem[j] = 1.0 / mass_wheel[1];
			i++;
			j++;
		}
		i = 0;
		while (i < Constants::DIM) {
			A_rem[j] = 1.0 / mass_wheel[2];
			i++;
			j++;
		}
		i = 0;
		while (i < Constants::DIM) {
			A_rem[j] = 1.0 / mass_wheel[3];
			i++;
			j++;
		}
		i = 0;
		while (i < Constants::DIM) {
			A_rem[j] = 1.0 / mass_tyre[0];
			i++;
			j++;
		}
		i = 0;
		while (i < Constants::DIM) {
			A_rem[j] = 1.0 / mass_tyre[1];
			i++;
			j++;
		}
		i = 0;
		while (i < Constants::DIM) {
			A_rem[j] = 1.0 / mass_tyre[2];
			i++;
			j++;
		}
		i = 0;
		while (i < Constants::DIM) {
			A_rem[j] = 1.0 / mass_tyre[3];
			i++;
			j++;
		}

		compute_f_mem_alloc();

		if (used_solver == BROYDEN_CN) {
			MathLibrary::Solvers<T, MBD_method>::Broyden_CN(this, x_vector, complete_vector, this->h, this->num_iter, this->tol, this->max_iter);
		}
		else if (used_solver == RUNGE_KUTTA_4) {
			MathLibrary::Solvers<T, MBD_method>::RK4(this, x_vector, complete_vector, this->h, this->num_iter, this->tol, this->max_iter);
		}
		else if (used_solver == BROYDEN_BDF2) {
			MathLibrary::Solvers<T, MBD_method>::Broyden_PDF2(this, x_vector, complete_vector, this->h, this->num_iter, this->tol, this->max_iter);
		}
		else if (used_solver == BROYDEN_EULER) {
			MathLibrary::Solvers<T, MBD_method>::Broyden_Euler(this, x_vector, complete_vector, this->h, this->num_iter, this->tol, this->max_iter);
		}
		else if (used_solver == EXPLICIT_EULER) {
			std::cout << "Explicit solver hasn't been implemented, you don't want to use it" << std::endl;
		}
		else {
			std::cout << "sorry man, the solver you picked for MBD is weird and hasn't been implemented yet" << std::endl;
		}

		compute_f_clean();

		T* start = complete_vector + (this->num_iter) * this->solution_dim;
		mkl<T>::copy(this->solution_dim, start, 1, solution_vector, 1);
		//	std::cout << "Solution copied!\n" << std::endl;
		mkl_free(complete_vector);
		mkl_free(x_vector);
	}

	/*
	Return the vector alignment (of the system)
	*/
	size_t get_alignment() {
		return Constants::ALIGNMENT;
	}

	/*
	Returns the dimension of the solution vector (=61)
	*/
	size_t get_solution_dimension() {
		return this->solution_dim;
	}

	/*
	Beatiful output of the result
	\param sln solution vector
	*/
	void print_final_result(T* sln) {
		std::cout << "MBD: angular velocity w=\n\t["<<sln[0] << "\n\t " << sln[1] << "\n\t " << sln[2] << "]" << std::endl;
		std::cout << "MBD: car body velocity vc=\n\t[" << sln[3] << "\n\t " << sln[4] << "\n\t " << sln[5] << "]" << std::endl;
		std::cout << "MBD: rear-right wheel velocity vw1=\n\t[" << sln[6] << "\n\t " << sln[7] << "\n\t " << sln[8] << "]" << std::endl;
		std::cout << "MBD: rear-left wheel velocity vw2=\n\t[" << sln[9] << "\n\t " << sln[10] << "\n\t " << sln[11] << "]" << std::endl;
		std::cout << "MBD: front-left wheel velocity vw3=\n\t[" << sln[12] << "\n\t " << sln[13] << "\n\t " << sln[14] << "]" << std::endl;
		std::cout << "MBD: front-right wheel velocity vw4=\n\t[" << sln[15] << "\n\t " << sln[16] << "\n\t " << sln[17] << "]" << std::endl;
		std::cout << "MBD: rear-right tyre velocity vt1=\n\t[" << sln[18] << "\n\t " << sln[19] << "\n\t " << sln[20] << "]" << std::endl;
		std::cout << "MBD: rear-left tyre velocity vt2=\n\t[" << sln[21] << "\n\t " << sln[22] << "\n\t " << sln[23] << "]" << std::endl;
		std::cout << "MBD: front-left tyre velocity vt3=\n\t[" << sln[24] << "\n\t " << sln[25] << "\n\t " << sln[26] << "]" << std::endl;
		std::cout << "MBD: front-right tyre velocity vt4=\n\t[" << sln[27] << "\n\t " << sln[28] << "\n\t " << sln[29] << "]" << std::endl;
		std::cout << "MBD: orientation q=\n\t[" << sln[30] << "\n\t " << sln[31] << "\n\t " << sln[32] << "\n\t " << sln[33] << "]" << std::endl;
		std::cout << "MBD: car body position pc=\n\t[" << sln[34] << "\n\t " << sln[35] << "\n\t " << sln[36] << "]" << std::endl;
		std::cout << "MBD: rear-right wheel position pw1=\n\t[" << sln[37] << "\n\t " << sln[38] << "\n\t " << sln[39] << "]" << std::endl;
		std::cout << "MBD: rear-left wheel position pw2=\n\t[" << sln[40] << "\n\t " << sln[41] << "\n\t " << sln[42] << "]" << std::endl;
		std::cout << "MBD: front-left wheel position pw3=\n\t[" << sln[43] << "\n\t " << sln[44] << "\n\t " << sln[47] << "]" << std::endl;
		std::cout << "MBD: front-right wheel position pw4=\n\t[" << sln[46] << "\n\t " << sln[47] << "\n\t " << sln[48] << "]" << std::endl;
		std::cout << "MBD: rear-right tyre position pt1=\n\t[" << sln[49] << "\n\t " << sln[50] << "\n\t " << sln[51] << "]" << std::endl;
		std::cout << "MBD: rear-left tyre position pt2=\n\t[" << sln[52] << "\n\t " << sln[53] << "\n\t " << sln[54] << "]" << std::endl;
		std::cout << "MBD: front-left tyre position pt3=\n\t[" << sln[55] << "\n\t " << sln[56] << "\n\t " << sln[57] << "]" << std::endl;
		std::cout << "MBD: front-right tyre position pt4=\n\t[" << sln[58] << "\n\t " << sln[59] << "\n\t " << sln[60] << "]" << std::endl;

	}

	/*
	Destructor
	*/
	~MBD_method() {
		mkl_free_buffers();
		mkl_free(r_fl);
		mkl_free(r_fr);
		mkl_free(r_rl);
		mkl_free(r_rr);
		mkl_free(Ic);
		mkl_free(mass_wheel);
		mkl_free(mass_tyre);
		mkl_free(upper_spring_length);
		mkl_free(lower_spring_length);
		mkl_free(upper_spring_stiffness);
		mkl_free(lower_spring_stiffness);
		mkl_free(upper_rotational_stiffness);
		mkl_free(lower_rotational_stiffness);
		mkl_free(initial_upper_spring_length);
		mkl_free(initial_lower_spring_length);
		mkl_free(initial_orientation);
		mkl_free(initial_angular_velocity);
		mkl_free(upper_spring_damping);
		mkl_free(lower_spring_damping);
		mkl_free(vc);
		mkl_free(vw_fl);
		mkl_free(vw_fr);
		mkl_free(vw_rl);
		mkl_free(vw_rr);
		mkl_free(vt_fl);
		mkl_free(vt_fr);
		mkl_free(vt_rl);
		mkl_free(vt_rr);
		mkl_free(pcc);
		mkl_free(FC);
		mkl_free(r_fl_tilda);
		mkl_free(r_fr_tilda);
		mkl_free(r_rl_tilda);
		mkl_free(r_rr_tilda);
		mkl_free(FW);
		mkl_free(FT);
		mkl_free(A_Ic);
		mkl_free(A_rem);
	}

};
