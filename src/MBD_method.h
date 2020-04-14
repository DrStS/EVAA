#pragma once
#include <string>
#include <mkl.h>
#include <chrono>
#include "MathLibrary.h"
#include "ReadXML.h"

/*
Handles the whole MBD simulation (from the Matlab code)
*/
template <class T>
class MBD_method
{
private:
	////////////////////////////// Simulation Parameters ///////////////////////////////////////////////////////////////
	const int DIM = 3;
	const int NUM_LEGS = 4;
	const int alignment = 64;
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
	T* vc, * vw1, * vw2, * vw3, * vw4, * vt1, * vt2, * vt3, * vt4; // velocity of the center of mass, wheel and tyre.

	//////////////////////////	External Forces terms //////////////////////////////////////////////////////
	T g;
	T* FC;
	T FT1, FT2, FT3, FT4;
	T FW1, FW2, FW3, FW4;
	T FR1, FR2, FR3, FR4;

	//////////////////////////// Initial condition params //////////////////////////////////////////////////
	T* initial_upper_spring_length, * initial_lower_spring_length, * initial_orientation, * initial_angular_velocity;


	//////////////////////////// Arrays for intermediate steps /////////////////////////////////////////////
	T* r1, * r2, * r3, * r4;
	T* pcc; // position of center of mass
	T* x_vector;


	////////////////////////// Auxillary Parameters for Compute_f function //////////////////////////////////
	T* r1_tilda, * r2_tilda, * r3_tilda, * r4_tilda, ** FW, * FT, * A_Ic, * A_rem;


	/////////////////////////	Variables needed in compute_f function  /////////////////////////////////////
	T* cf_C_cN, * cf_r_up1, * cf_r_up2, * cf_r_up3, * cf_r_up4, * cf_r_low1, * cf_r_low2, * cf_r_low3, * cf_r_low4;
	T* cf_upper_normal1, * cf_upper_normal2, * cf_upper_normal3, * cf_upper_normal4, * cf_lower_normal1, * cf_lower_normal2, * cf_lower_normal3, * cf_lower_normal4, * cf_col_dat;
	T* cf_upper_force1, * cf_upper_force2, * cf_upper_force3, * cf_upper_force4, * cf_lower_force1, * cf_lower_force2, * cf_lower_force3, * cf_lower_force4;
	T* cf_upper_dampf1, * cf_upper_dampf2, * cf_upper_dampf3, * cf_upper_dampf4, * cf_lower_dampf1, * cf_lower_dampf2, * cf_lower_dampf3, * cf_lower_dampf4, * cf_temp;
	T* cf_upper_angle1, * cf_upper_angle2, * cf_upper_angle3, * cf_upper_angle4, * cf_lower_angle1, * cf_lower_angle2, * cf_lower_angle3, * cf_lower_angle4;
	T* cf_upper_S1, * cf_upper_S2, * cf_upper_S3, * cf_upper_S4, * cf_lower_S1, * cf_lower_S2, * cf_lower_S3, * cf_lower_S4;
	T* cf_lower_rot_force1, * cf_lower_rot_force2, * cf_lower_rot_force3, * cf_lower_rot_force4, * cf_upper_rot_force1, * cf_upper_rot_force2, * cf_upper_rot_force3, * cf_upper_rot_force4;
	T* cf_car_rot_force1, * cf_car_rot_force2, * cf_car_rot_force3, * cf_car_rot_force4, * cf_sum_car_force1, * cf_sum_car_force2, * cf_sum_car_force3, * cf_sum_car_force4;
	T* cf_local_FR1, * cf_local_FR2, * cf_local_FR3, * cf_local_FR4;
	T* cf_Hc, * cf_sum_torque_spring_car, * cf_Tc, * cf_wc_tilda;
	T* cf_b_rem, * cf_Qc, * cf_qc_dot;
	T* current_spring_lengths, * stiffness_vector;


	////////////////////////// Lookup table ///////////////////////////////////////////
	EVAAComputeStiffness* lookupStiffness;

	/*
	Calculate the positions of the tyres and wheels according to the initial orientation of the car
	The legs always form a 90° angle to the car body, such that the rotational springs are at rest
	*/
	void get_initial_length(
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
	MBD_method(const Simulation_Parameters& params, const Load_Params& load_params, EVAAComputeStiffness* interpolator) {

		////////////////////////////// Simulation Parameters ///////////////////////////////////////////////////////////////
		h = params.timestep;
		num_iter = params.num_time_iter;
		max_iter = params.max_num_iter;
		tol = params.tolerance;
		solution_dim = params.solution_dim; /// this is by the formulation
		used_solver = params.solver;
		boundary_conditions = load_params.boundary_condition_road;
		radius_circular_path = load_params.profile_radius;
		use_interpolation = params.interpolation;
		lookupStiffness = interpolator;

		////////////////////////////// Car Definition ///////////////////////////////////////////////////////////////////////
		k_body_fl = params.k_body[2];
		k_tyre_fl = params.k_tyre[2];
		k_body_fr = params.k_body[3];
		k_tyre_fr = params.k_tyre[3];
		k_body_rl = params.k_body[1];
		k_tyre_rl = params.k_tyre[1];
		k_body_rr = params.k_body[0];
		k_tyre_rr = params.k_tyre[0];

		k_body_rot_fl = 1e4;
		k_body_rot_fr = 1e4;
		k_body_rot_rl = 1e4;
		k_body_rot_rr = 1e4;
		k_tyre_rot_fl = 1e4;
		k_tyre_rot_fr = 1e4;
		k_tyre_rot_rl = 1e4;
		k_tyre_rot_rr = 1e4;

		c_body_fl = params.c_body[2];
		c_tyre_fl = params.c_tyre[2];
		c_body_fr = params.c_body[3];
		c_tyre_fr = params.c_tyre[3];
		c_body_rl = params.c_body[1];
		c_tyre_rl = params.c_tyre[1];
		c_body_rr = params.c_body[0];
		c_tyre_rr = params.c_tyre[0];

		l_long_fl = params.l_long[2];
		l_long_fr = params.l_long[3];
		l_long_rl = params.l_long[1];
		l_long_rr = params.l_long[0];

		l_lat_fl = params.l_lat[2];
		l_lat_fr = params.l_lat[3];
		l_lat_rl = params.l_lat[1];
		l_lat_rr = params.l_lat[0];

		mass = params.mass_body;

		I_body_xx = params.I_body[0];
		I_body_yy = params.I_body[1];
		I_body_zz = params.I_body[2];

		mass_wheel_fl = params.mass_wheel[2];
		mass_tyre_fl = params.mass_tyre[2];
		mass_wheel_fr = params.mass_wheel[3];
		mass_tyre_fr = params.mass_tyre[3];
		mass_wheel_rl = params.mass_wheel[1];
		mass_tyre_rl = params.mass_tyre[1];
		mass_wheel_rr = params.mass_wheel[0];
		mass_tyre_rr = params.mass_tyre[0];

		upper_spring_length_rr = params.upper_spring_length[0];
		upper_spring_length_rl = params.upper_spring_length[1];
		upper_spring_length_fl = params.upper_spring_length[2];
		upper_spring_length_fr = params.upper_spring_length[3];

		lower_spring_length_rr = params.lower_spring_length[0];
		lower_spring_length_rl = params.lower_spring_length[1];
		lower_spring_length_fl = params.lower_spring_length[2];
		lower_spring_length_fr = params.lower_spring_length[3];

		g = params.gravity[2];

		int i;

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////// Memory Allocation and matrix formulation ///////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		r1 = (T*)mkl_calloc(this->DIM, sizeof(T), this->alignment);
		r2 = (T*)mkl_calloc(this->DIM, sizeof(T), this->alignment);
		r3 = (T*)mkl_calloc(this->DIM, sizeof(T), this->alignment);
		r4 = (T*)mkl_calloc(this->DIM, sizeof(T), this->alignment);
		Ic = (T*)mkl_calloc((this->DIM) * (this->DIM), sizeof(T), this->alignment);
		mass_wheel = (T*)mkl_calloc(this->NUM_LEGS, sizeof(T), this->alignment);
		mass_tyre = (T*)mkl_calloc(this->NUM_LEGS, sizeof(T), this->alignment);
		upper_spring_length = (T*)mkl_calloc(this->NUM_LEGS, sizeof(T), this->alignment);
		lower_spring_length = (T*)mkl_calloc(this->NUM_LEGS, sizeof(T), this->alignment);
		upper_spring_stiffness = (T*)mkl_calloc(this->NUM_LEGS, sizeof(T), this->alignment);
		lower_spring_stiffness = (T*)mkl_calloc(this->NUM_LEGS, sizeof(T), this->alignment);
		upper_rotational_stiffness = (T*)mkl_calloc(this->NUM_LEGS, sizeof(T), this->alignment);
		lower_rotational_stiffness = (T*)mkl_calloc(this->NUM_LEGS, sizeof(T), this->alignment);
		initial_upper_spring_length = (T*)mkl_calloc(this->NUM_LEGS, sizeof(T), this->alignment);
		initial_lower_spring_length = (T*)mkl_calloc(this->NUM_LEGS, sizeof(T), this->alignment);
		initial_orientation = (T*)mkl_calloc(this->NUM_LEGS, sizeof(T), this->alignment);
		initial_angular_velocity = (T*)mkl_calloc(this->DIM, sizeof(T), this->alignment);
		upper_spring_damping = (T*)mkl_calloc(this->NUM_LEGS, sizeof(T), this->alignment);
		lower_spring_damping = (T*)mkl_calloc(this->NUM_LEGS, sizeof(T), this->alignment);
		vc = (T*)mkl_calloc(this->DIM, sizeof(T), this->alignment);
		vw1 = (T*)mkl_calloc(this->DIM, sizeof(T), this->alignment);
		vw2 = (T*)mkl_calloc(this->DIM, sizeof(T), this->alignment);
		vw3 = (T*)mkl_calloc(this->DIM, sizeof(T), this->alignment);
		vw4 = (T*)mkl_calloc(this->DIM, sizeof(T), this->alignment);
		vt1 = (T*)mkl_calloc(this->DIM, sizeof(T), this->alignment);
		vt2 = (T*)mkl_calloc(this->DIM, sizeof(T), this->alignment);
		vt3 = (T*)mkl_calloc(this->DIM, sizeof(T), this->alignment);
		vt4 = (T*)mkl_calloc(this->DIM, sizeof(T), this->alignment);
		pcc = (T*)mkl_calloc(this->DIM, sizeof(T), this->alignment);
		FC = (T*)mkl_calloc(this->DIM, sizeof(T), this->alignment);
		r1_tilda = (T*)mkl_calloc((this->DIM) * (this->DIM), sizeof(T), this->alignment);
		r2_tilda = (T*)mkl_calloc((this->DIM) * (this->DIM), sizeof(T), this->alignment);
		r3_tilda = (T*)mkl_calloc((this->DIM) * (this->DIM), sizeof(T), this->alignment);
		r4_tilda = (T*)mkl_calloc((this->DIM) * (this->DIM), sizeof(T), this->alignment);
		FW = (T**)mkl_calloc((this->NUM_LEGS), sizeof(T*), this->alignment);
		FT = (T*)mkl_calloc((this->NUM_LEGS), sizeof(T*), this->alignment);
		A_Ic = (T*)mkl_calloc((this->DIM) * (this->DIM), sizeof(T), this->alignment);
		A_rem = (T*)mkl_calloc(9 * this->DIM, sizeof(T), this->alignment);
		center_of_circle = (T*)mkl_calloc(this->DIM, sizeof(T), this->alignment);


		i = 0;
		r1[i] = -l_long_rr; r2[i] = -l_long_rl; r3[i] = l_long_fl; r4[i] = l_long_fr;
		Ic[i * DIM + i] = I_body_xx;
		mass_wheel[i] = mass_wheel_rr;
		mass_tyre[i] = mass_tyre_rr;
		upper_spring_length[i] = upper_spring_length_rr;
		lower_spring_length[i] = lower_spring_length_rr;
		upper_spring_stiffness[i] = k_body_rr;
		lower_spring_stiffness[i] = k_tyre_rr;
		upper_rotational_stiffness[i] = k_body_rot_rr;
		lower_rotational_stiffness[i] = k_tyre_rot_rr;
		FT1 = -mass_tyre[i] * g;
		FW1 = -mass_wheel[i] * g;
		upper_spring_damping[i] = c_body_rr;
		lower_spring_damping[i] = c_tyre_rr;
		initial_upper_spring_length[i] = params.initial_upper_spring_length[i];
		initial_lower_spring_length[i] = params.initial_lower_spring_length[i];
		initial_orientation[i] = params.initial_angle[i];
		vc[i] = params.initial_vel_body[i];
		vw1[i] = params.initial_vel_wheel[i];
		vw2[i] = params.initial_vel_wheel[i+3];
		vw3[i] = params.initial_vel_wheel[i+6];
		vw4[i] = params.initial_vel_wheel[i+9];
		vt1[i] = params.initial_vel_tyre[i];
		vt2[i] = params.initial_vel_tyre[i+3];
		vt3[i] = params.initial_vel_tyre[i+6];
		vt4[i] = params.initial_vel_tyre[i+9];
		pcc[i] = params.initial_pos_body[i]; 
		FC[i] = 0;
		center_of_circle[i] = load_params.profile_center[i];

		i = 1;
		r1[i] = 0; r2[i] = 0; r3[i] = 0; r4[i] = 0;
		Ic[i * DIM + i] = I_body_zz;
		mass_wheel[i] = mass_wheel_rl;
		mass_tyre[i] = mass_tyre_rl;
		upper_spring_length[i] = upper_spring_length_rl;
		lower_spring_length[i] = lower_spring_length_rl;
		upper_spring_stiffness[i] = k_body_rl;
		lower_spring_stiffness[i] = k_tyre_rl;
		upper_rotational_stiffness[i] = k_body_rot_rl;
		lower_rotational_stiffness[i] = k_tyre_rot_rl;
		FT2 = -mass_tyre[i] * g;
		FW2 = -mass_wheel[i] * g;
		upper_spring_damping[i] = c_body_rl;
		lower_spring_damping[i] = c_tyre_rl;
		initial_upper_spring_length[i] = params.initial_upper_spring_length[i];
		initial_lower_spring_length[i] = params.initial_lower_spring_length[i];
		initial_orientation[i] = params.initial_angle[i];
		vc[i] = params.initial_vel_body[i+1];
		vw1[i] = params.initial_vel_wheel[i+1];
		vw2[i] = params.initial_vel_wheel[i+1 + 3];
		vw3[i] = params.initial_vel_wheel[i+1 + 6];
		vw4[i] = params.initial_vel_wheel[i+1 + 9];
		vt1[i] = params.initial_vel_tyre[i+1];
		vt2[i] = params.initial_vel_tyre[i+1 + 3];
		vt3[i] = params.initial_vel_tyre[i+1 + 6];
		vt4[i] = params.initial_vel_tyre[i+1 + 9];
		pcc[i] = params.initial_pos_body[i+1];
		FC[i] = -mass * g;
		center_of_circle[i] = load_params.profile_center[i+1];

		i = 2;
		r1[i] = l_lat_rr; r2[i] = -l_lat_rl; r3[i] = -l_lat_fl; r4[i] = l_lat_fr;
		Ic[i * DIM + i] = I_body_yy;
		mass_wheel[i] = mass_wheel_fl;
		mass_tyre[i] = mass_tyre_fl;
		upper_spring_length[i] = upper_spring_length_fl;
		lower_spring_length[i] = lower_spring_length_fl;
		upper_spring_stiffness[i] = k_body_fl;
		lower_spring_stiffness[i] = k_tyre_fl;
		upper_rotational_stiffness[i] = k_body_rot_fl;
		lower_rotational_stiffness[i] = k_tyre_rot_fl;
		FT3 = -mass_tyre[i] * g;
		FW3 = -mass_wheel[i] * g;
		upper_spring_damping[i] = c_body_fl;
		lower_spring_damping[i] = c_tyre_fl;
		initial_upper_spring_length[i] = params.initial_upper_spring_length[i];
		initial_lower_spring_length[i] = params.initial_lower_spring_length[i];
		initial_orientation[i] = params.initial_angle[i];
		vc[i] = params.initial_vel_body[i-1];
		vw1[i] = params.initial_vel_wheel[i-1];
		vw2[i] = params.initial_vel_wheel[i-1 + 3];
		vw3[i] = params.initial_vel_wheel[i-1 + 6];
		vw4[i] = params.initial_vel_wheel[i-1 + 9];
		vt1[i] = params.initial_vel_tyre[i-1];
		vt2[i] = params.initial_vel_tyre[i-1 + 3];
		vt3[i] = params.initial_vel_tyre[i-1 + 6];
		vt4[i] = params.initial_vel_tyre[i-1 + 9];
		pcc[i] = params.initial_pos_body[i-1];
		FC[i] = 0;
		center_of_circle[i] = load_params.profile_center[i-1];

		i = 3;
		mass_wheel[i] = mass_wheel_fr;
		mass_tyre[i] = mass_tyre_fr;
		upper_spring_length[i] = upper_spring_length_fr;
		lower_spring_length[i] = lower_spring_length_fr;
		upper_spring_stiffness[i] = k_body_fr;
		lower_spring_stiffness[i] = k_tyre_fr;
		upper_rotational_stiffness[i] = k_body_rot_fr;
		lower_rotational_stiffness[i] = k_tyre_rot_fr;
		FT4 = -mass_tyre[i] * g;
		FW4 = -mass_wheel[i] * g;
		upper_spring_damping[i] = c_body_fr;
		lower_spring_damping[i] = c_tyre_fr;
		initial_upper_spring_length[i] = params.initial_upper_spring_length[i];
		initial_lower_spring_length[i] = params.initial_lower_spring_length[i];
		initial_orientation[i] = params.initial_angle[i];

		// A_Ic has cholesky factorization of Ic
		cblas_dcopy(this->DIM * this->DIM, Ic, 1, A_Ic, 1);
		LAPACKE_dpotrf(LAPACK_ROW_MAJOR, 'L', this->DIM, A_Ic, this->DIM);
	}

	/*
	Updates the velocity of the 4 wheels and tyres as well as the angular velocity, such that the car is already in the trajectory of the circle
	overwrites the velocity in the tyres and wheels and the angular velocities
	only keeps the tangential component of the velocity of the car body 
	*/
	void circular_path_initialization(T* vc, T* vw1, T* vw2, T* vw3, T* vw4, 
		T* vt1, T* vt2, T* vt3, T* vt4, T* omega, T* pcc, 
		T* pt1 , T* pt2, T* pt3, T* pt4, T &radius_param) {
		const MKL_INT dim = this->DIM;
		const MKL_INT incx = 1;
 
		vc[1] = 0;

		T* perpendicular_dir = (T*)mkl_calloc(this->DIM, sizeof(T), this->alignment);
		T* tangential_dir = (T*)mkl_calloc(this->DIM, sizeof(T), this->alignment);
		T* radial_vector = (T*)mkl_calloc(this->DIM, sizeof(T), this->alignment);

		radial_vector[0] = pcc[0] - center_of_circle[0];
		radial_vector[1] = 0;
		radial_vector[2] = pcc[2] - center_of_circle[2];

		perpendicular_dir[1] = 1;

		T radius = cblas_dnrm2(this->DIM, radial_vector, 1);

		if (abs(radius - radius_param) > 0.01)
			std::cout << "Warning! the initial position of the car is not on the trajectory provided in the circular path. \n The expected radius is " << radius_circular_path << ", but the car is at an initial distance of " << radius << " from the center of the circle.\n The execution procedes with the current spatial configuration and with the current distance to the center of the circle." << std::endl;

		T inv_radius = 1. / radius;

		cblas_dscal(dim, inv_radius, radial_vector, incx);

		MathLibrary::crossProduct(radial_vector, perpendicular_dir, tangential_dir);

		T magnitude = cblas_ddot(dim, vc, incx, tangential_dir, incx);

		cblas_dcopy(this->DIM, tangential_dir, 1, vc, 1);

		cblas_dscal(dim, magnitude, vc, incx);

		cblas_dscal(dim, radius, radial_vector, incx);

		MathLibrary::crossProduct(radial_vector, vc, omega);

		cblas_dscal(this->DIM, inv_radius * inv_radius, omega, 1);

		MathLibrary::crossProduct(omega, pt1, vt1);
		MathLibrary::crossProduct(omega, pt2, vt2);
		MathLibrary::crossProduct(omega, pt3, vt3);
		MathLibrary::crossProduct(omega, pt4, vt4);

		cblas_dcopy(this->DIM, vt1, 1, vw1, 1);
		cblas_dcopy(this->DIM, vt2, 1, vw2, 1);
		cblas_dcopy(this->DIM, vt3, 1, vw3, 1);
		cblas_dcopy(this->DIM, vt4, 1, vw4, 1);

	}


	/*
	Fixes the tyres to their initial position
	The forces acting on the tyres are now always zero
	\param v is the velocity of the tyre
	\param m is the mass of the tyre
	\param p is the global position of the tyre
	\result Fr The only force acting on the tyre
	*/
	void get_fixed_road_force(T* Fr1, T* Fr2, T* Fr3, T* Fr4) {
		Fr1[0] = 0.0;
		Fr1[1] = 0.0;
		Fr1[2] = 0.0;

		Fr2[0] = 0.0;
		Fr2[1] = 0.0;
		Fr2[2] = 0.0;

		Fr3[0] = 0.0;
		Fr3[1] = 0.0;
		Fr3[2] = 0.0;

		Fr4[0] = 0.0;
		Fr4[1] = 0.0;
		Fr4[2] = 0.0;
	}


	/*
	No interaction with the road, no additional forces on the tyres
	\param v is the velocity of the tyre
	\param m is the mass of the tyre
	\param p is the global position of the tyre
	\result Fr The only force acting on the tyre
	*/
	void get_nonfixed_road_force(T* Fr1, T* Fr2, T* Fr3, T* Fr4) {
	}

	/* calculates the force in the tyre only with respect to its velocity, mass and position
	\param v is the velocity of the tyre
	\param m is the mass of the tyre
	\param p is the global position of the tyre
	\result Fr The only force acting on the tyre
	*/
	void get_circular_road_force(T* Fr, T* v, T &m, T* p) {

		T *unit_y_vector, *velocity_direction_tyre;
		T velocity_magnitude_tyre, inv_radius, force_magnitude_tyre;

		unit_y_vector = (T*)mkl_calloc(this->DIM, sizeof(T), this->alignment);
		velocity_direction_tyre = (T*)mkl_calloc(this->DIM, sizeof(T), this->alignment);


		const MKL_INT mkl_DIM = this->DIM;
		const MKL_INT mkl_incx = 1;
		const MKL_INT mkl_incy = 1;

		cblas_dcopy(this->DIM, p, 1, Fr, 1);

		Fr[1] = 0;		// path only in XZ-plane

		inv_radius = 1.0 / cblas_dnrm2(this->DIM, p, 1);		// corresponds to the (inverse) radius of the trajectory at the considered tyre

		cblas_dscal(this->DIM, -inv_radius, Fr, 1);

		unit_y_vector[0] = 0.;
		unit_y_vector[1] = 1.;
		unit_y_vector[2] = 0.;

		MathLibrary::crossProduct(Fr, unit_y_vector, velocity_direction_tyre);

		velocity_magnitude_tyre = cblas_ddot(mkl_DIM, v, mkl_incx, velocity_direction_tyre, mkl_incy);

		force_magnitude_tyre = m * velocity_magnitude_tyre * velocity_magnitude_tyre * inv_radius;

		cblas_dscal(this->DIM, force_magnitude_tyre, Fr, 1);

		MKL_free(unit_y_vector);
		MKL_free(velocity_direction_tyre);
	}


	/*
	Memory allocation of all the variables required in the solve function
	To increase performance by removing repeting memory allocations
	The same locations are overwritten at each timestep
	*/
	void compute_f_mem_alloc() {
		// add to members


		cf_C_cN = (T*)mkl_calloc((this->DIM) * (this->DIM), sizeof(T), this->alignment);
		cf_r_up1 = (T*)mkl_calloc(this->DIM, sizeof(T), this->alignment);
		cf_r_up2 = (T*)mkl_calloc(this->DIM, sizeof(T), this->alignment);
		cf_r_up3 = (T*)mkl_calloc(this->DIM, sizeof(T), this->alignment);	
		cf_r_up4 = (T*)mkl_calloc(this->DIM, sizeof(T), this->alignment);
		cf_r_low1 = (T*)mkl_calloc(this->DIM, sizeof(T), this->alignment);
		cf_r_low2 = (T*)mkl_calloc(this->DIM, sizeof(T), this->alignment);
		cf_r_low3 = (T*)mkl_calloc(this->DIM, sizeof(T), this->alignment);
		cf_r_low4 = (T*)mkl_calloc(this->DIM, sizeof(T), this->alignment);

		cf_upper_normal1 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		cf_upper_normal2 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		cf_upper_normal3 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		cf_upper_normal4 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		cf_lower_normal1 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		cf_lower_normal2 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		cf_lower_normal3 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		cf_lower_normal4 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		cf_col_dat = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);

		cf_upper_force1 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		cf_upper_force2 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		cf_upper_force3 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		cf_upper_force4 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		cf_lower_force1 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		cf_lower_force2 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		cf_lower_force3 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		cf_lower_force4 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);

		cf_upper_dampf1 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		cf_upper_dampf2 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		cf_upper_dampf3 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		cf_upper_dampf4 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		cf_lower_dampf1 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		cf_lower_dampf2 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		cf_lower_dampf3 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		cf_lower_dampf4 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		cf_temp = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);

		cf_upper_angle1 = new T;
		cf_upper_angle2 = new T;
		cf_upper_angle3 = new T;
		cf_upper_angle4 = new T;
		cf_lower_angle1 = new T;
		cf_lower_angle2 = new T;
		cf_lower_angle3 = new T;
		cf_lower_angle4 = new T;

		cf_upper_S1 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		cf_upper_S2 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		cf_upper_S3 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		cf_upper_S4 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		cf_lower_S1 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		cf_lower_S2 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		cf_lower_S3 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		cf_lower_S4 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);

		cf_lower_rot_force1 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		cf_lower_rot_force2 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		cf_lower_rot_force3 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		cf_lower_rot_force4 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		cf_upper_rot_force1 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		cf_upper_rot_force2 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		cf_upper_rot_force3 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		cf_upper_rot_force4 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		cf_car_rot_force1 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		cf_car_rot_force2 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		cf_car_rot_force3 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		cf_car_rot_force4 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		cf_sum_car_force1 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		cf_sum_car_force2 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		cf_sum_car_force3 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		cf_sum_car_force4 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);

		cf_local_FR1 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		cf_local_FR2 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		cf_local_FR3 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		cf_local_FR4 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);

		cf_Hc = (T*)mkl_calloc(this->DIM, sizeof(T), this->alignment);
		cf_sum_torque_spring_car = (T*)mkl_calloc(this->DIM, sizeof(T), this->alignment);
		cf_Tc = (T*)mkl_calloc(this->DIM, sizeof(T), this->alignment);
		cf_wc_tilda = (T*)mkl_calloc((this->DIM) * (this->DIM), sizeof(T), this->alignment);
		cf_Tc[0] = 0.0;
		cf_Tc[1] = 0.0;
		cf_Tc[2] = 0.0;

		cf_b_rem = (T*)mkl_calloc(9 * this->DIM, sizeof(T), this->alignment);
		cf_Qc = (T*)mkl_calloc((this->NUM_LEGS) * (this->DIM), sizeof(T), this->alignment);
		cf_qc_dot = (T*)mkl_calloc((this->NUM_LEGS), sizeof(T), this->alignment);

		// required for the look up table
		current_spring_lengths = (T*)mkl_calloc(8 * this->DIM, sizeof(T), this->alignment);
		stiffness_vector = (T*)mkl_calloc(8 * this->DIM, sizeof(T), this->alignment);

	}

	/*
	Clean the memory allocated for the main solver
	*/
	void compute_f_clean() {
		mkl_free_buffers();
		mkl_free(cf_C_cN);
		mkl_free(cf_r_up1);
		mkl_free(cf_r_up2);
		mkl_free(cf_r_up3);
		mkl_free(cf_r_up4);
		mkl_free(cf_r_low1);
		mkl_free(cf_r_low2);
		mkl_free(cf_r_low3);
		mkl_free(cf_r_low4);
		mkl_free(cf_upper_normal1);
		mkl_free(cf_upper_normal2);
		mkl_free(cf_upper_normal3);
		mkl_free(cf_upper_normal4);
		mkl_free(cf_lower_normal1);
		mkl_free(cf_lower_normal2);
		mkl_free(cf_lower_normal3);
		mkl_free(cf_lower_normal4);
		mkl_free(cf_col_dat);
		mkl_free(cf_upper_force1);
		mkl_free(cf_upper_force2);
		mkl_free(cf_upper_force3);
		mkl_free(cf_upper_force4);
		mkl_free(cf_lower_force1);
		mkl_free(cf_lower_force2);
		mkl_free(cf_lower_force3);
		mkl_free(cf_lower_force4);
		mkl_free(cf_upper_dampf1);
		mkl_free(cf_upper_dampf2);
		mkl_free(cf_upper_dampf3);
		mkl_free(cf_upper_dampf4);
		mkl_free(cf_lower_dampf1);
		mkl_free(cf_lower_dampf2);
		mkl_free(cf_lower_dampf3);
		mkl_free(cf_lower_dampf4);
		mkl_free(cf_temp);
		delete cf_upper_angle1;
		delete cf_upper_angle2;
		delete cf_upper_angle3;
		delete cf_upper_angle4;
		delete cf_lower_angle1;
		delete cf_lower_angle2;
		delete cf_lower_angle3;
		delete cf_lower_angle4;
		mkl_free(cf_upper_S1);
		mkl_free(cf_upper_S2);
		mkl_free(cf_upper_S3);
		mkl_free(cf_upper_S4);
		mkl_free(cf_lower_S1);
		mkl_free(cf_lower_S2);
		mkl_free(cf_lower_S3);
		mkl_free(cf_lower_S4);
		mkl_free(cf_lower_rot_force1);
		mkl_free(cf_lower_rot_force2);
		mkl_free(cf_lower_rot_force3);
		mkl_free(cf_lower_rot_force4);
		mkl_free(cf_upper_rot_force1);
		mkl_free(cf_upper_rot_force2);
		mkl_free(cf_upper_rot_force3);
		mkl_free(cf_upper_rot_force4);
		mkl_free(cf_car_rot_force1);
		mkl_free(cf_car_rot_force2);
		mkl_free(cf_car_rot_force3);
		mkl_free(cf_car_rot_force4);
		mkl_free(cf_sum_car_force1);
		mkl_free(cf_sum_car_force2);
		mkl_free(cf_sum_car_force3);
		mkl_free(cf_sum_car_force4);
		mkl_free(cf_local_FR1);
		mkl_free(cf_local_FR2);
		mkl_free(cf_local_FR3);
		mkl_free(cf_local_FR4);
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
		T* vc_ = x_ + this->DIM;
		T* vw1_ = x_ + 2 * this->DIM;
		T* vw2_ = x_ + 3 * this->DIM;
		T* vw3_ = x_ + 4 * this->DIM;
		T* vw4_ = x_ + 5 * this->DIM;
		T* vt1_ = x_ + 6 * this->DIM;
		T* vt2_ = x_ + 7 * this->DIM;
		T* vt3_ = x_ + 8 * this->DIM;
		T* vt4_ = x_ + 9 * this->DIM;
		T* qc_ = x_ + 10 * this->DIM;
		T* pcc_ = x_ + 10 * this->DIM + this->NUM_LEGS;
		T* pw1_ = x_ + 11 * this->DIM + this->NUM_LEGS;
		T* pw2_ = x_ + 12 * this->DIM + this->NUM_LEGS;
		T* pw3_ = x_ + 13 * this->DIM + this->NUM_LEGS;
		T* pw4_ = x_ + 14 * this->DIM + this->NUM_LEGS;
		T* pt1_ = x_ + 15 * this->DIM + this->NUM_LEGS;
		T* pt2_ = x_ + 16 * this->DIM + this->NUM_LEGS;
		T* pt3_ = x_ + 17 * this->DIM + this->NUM_LEGS;
		T* pt4_ = x_ + 18 * this->DIM + this->NUM_LEGS;
		T inv_norm_r_up1, inv_norm_r_up2, inv_norm_r_up3, inv_norm_r_up4;
		T inv_norm_r_low1, inv_norm_r_low2, inv_norm_r_low3, inv_norm_r_low4;
		T norm_r_up1, norm_r_up2, norm_r_up3, norm_r_up4;
		T norm_r_low1, norm_r_low2, norm_r_low3, norm_r_low4;
		//write_vector(r1, this->DIM);
		//write_matrix(r1_tilda, this->DIM);
		/*
		/////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////                 Basis        //////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////
		get cosine transforms (C_Nc means r_N = C_Nc * r_c)
		compute local base vectors
		basis_c = get_basis(qc);
		*/

		MathLibrary::get_basis<T>(qc_, cf_C_cN);
		const MKL_INT mkl_DIM = this->DIM;
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
		// r_up1
		cblas_dcopy(this->DIM, pcc_, 1, cf_r_up1, 1);
		cblas_dgemv(CblasRowMajor, CblasNoTrans, this->DIM, this->DIM, 1, cf_C_cN, this->DIM, this->r1, 1, 1, cf_r_up1, 1);
		cblas_daxpy(this->DIM, -1.0, pw1_, 1, cf_r_up1, 1);

		// r_up2
		cblas_dcopy(this->DIM, pcc_, 1, cf_r_up2, 1);
		cblas_dgemv(CblasRowMajor, CblasNoTrans, this->DIM, this->DIM, 1, cf_C_cN, this->DIM, this->r2, 1, 1, cf_r_up2, 1);
		cblas_daxpy(this->DIM, -1.0, pw2_, 1, cf_r_up2, 1);
		// r_up3
		cblas_dcopy(this->DIM, pcc_, 1, cf_r_up3, 1);
		cblas_dgemv(CblasRowMajor, CblasNoTrans, this->DIM, this->DIM, 1, cf_C_cN, this->DIM, this->r3, 1, 1, cf_r_up3, 1);
		cblas_daxpy(this->DIM, -1.0, pw3_, 1, cf_r_up3, 1);
		// r_up4
		cblas_dcopy(this->DIM, pcc_, 1, cf_r_up4, 1);
		cblas_dgemv(CblasRowMajor, CblasNoTrans, this->DIM, this->DIM, 1, cf_C_cN, this->DIM, this->r4, 1, 1, cf_r_up4, 1);
		cblas_daxpy(this->DIM, -1.0, pw4_, 1, cf_r_up4, 1);
		
		// r_low1
		cblas_dcopy(this->DIM, pw1_, 1, cf_r_low1, 1);
		cblas_daxpy(this->DIM, -1.0, pt1_, 1, cf_r_low1, 1);

		// r_low2
		cblas_dcopy(this->DIM, pw2_, 1, cf_r_low2, 1);
		cblas_daxpy(this->DIM, -1.0, pt2_, 1, cf_r_low2, 1);

		// r_low3
		cblas_dcopy(this->DIM, pw3_, 1, cf_r_low3, 1);
		cblas_daxpy(this->DIM, -1.0, pt3_, 1, cf_r_low3, 1);
		// r_low4
		cblas_dcopy(this->DIM, pw4_, 1, cf_r_low4, 1);
		cblas_daxpy(this->DIM, -1.0, pt4_, 1, cf_r_low4, 1);

		/* Compute spring lengths and their inverses */
		norm_r_up1 = cblas_dnrm2(this->DIM, cf_r_up1, 1);
		norm_r_up2 = cblas_dnrm2(this->DIM, cf_r_up2, 1);
		norm_r_up3 = cblas_dnrm2(this->DIM, cf_r_up3, 1);
		norm_r_up4 = cblas_dnrm2(this->DIM, cf_r_up4, 1);

		inv_norm_r_up1 = 1.0 / norm_r_up1;
		inv_norm_r_up2 = 1.0 / norm_r_up2;
		inv_norm_r_up3 = 1.0 / norm_r_up3;
		inv_norm_r_up4 = 1.0 / norm_r_up4;

		norm_r_low1 = cblas_dnrm2(this->DIM, cf_r_low1, 1);
		norm_r_low2 = cblas_dnrm2(this->DIM, cf_r_low2, 1);
		norm_r_low3 = cblas_dnrm2(this->DIM, cf_r_low3, 1);
		norm_r_low4 = cblas_dnrm2(this->DIM, cf_r_low4, 1);

		inv_norm_r_low1 = 1.0 / norm_r_low1;
		inv_norm_r_low2 = 1.0 / norm_r_low2;
		inv_norm_r_low3 = 1.0 / norm_r_low3;
		inv_norm_r_low4 = 1.0 / norm_r_low4;


		/* Compute stiffness from lookup table*/
		if (use_interpolation) {
			//populate the lenght_vector
			current_spring_lengths[0] = norm_r_up3;
			current_spring_lengths[1] = norm_r_low3;
			current_spring_lengths[2] = norm_r_up4;
			current_spring_lengths[3] = norm_r_low4;
			current_spring_lengths[4] = norm_r_up2;
			current_spring_lengths[5] = norm_r_low2;
			current_spring_lengths[6] = norm_r_up1;
			current_spring_lengths[7] = norm_r_low1;
			// calculate the new stiffnesses
			lookupStiffness->getStiffness(current_spring_lengths, stiffness_vector);

			// overwrite stiffness values
			upper_spring_stiffness[2] = stiffness_vector[0];
			lower_spring_stiffness[2] = stiffness_vector[1];
			upper_spring_stiffness[3] = stiffness_vector[2];
			lower_spring_stiffness[3] = stiffness_vector[3];
			upper_spring_stiffness[1] = stiffness_vector[4];
			lower_spring_stiffness[1] = stiffness_vector[5];
			upper_spring_stiffness[0] = stiffness_vector[6];
			lower_spring_stiffness[0] = stiffness_vector[7];
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

		cblas_dcopy(this->DIM, cf_C_cN + 1, this->DIM, cf_col_dat, 1);

		MathLibrary::get_quaternion<T>(cf_r_up1, cf_col_dat, cf_upper_angle1, cf_upper_normal1, this->DIM);
		MathLibrary::get_quaternion<T>(cf_r_up2, cf_col_dat, cf_upper_angle2, cf_upper_normal2, this->DIM);
		MathLibrary::get_quaternion<T>(cf_r_up3, cf_col_dat, cf_upper_angle3, cf_upper_normal3, this->DIM);
		MathLibrary::get_quaternion<T>(cf_r_up4, cf_col_dat, cf_upper_angle4, cf_upper_normal4, this->DIM);

		MathLibrary::get_quaternion(cf_r_low1, cf_r_up1, cf_lower_angle1, cf_lower_normal1, this->DIM);
		MathLibrary::get_quaternion(cf_r_low2, cf_r_up2, cf_lower_angle2, cf_lower_normal2, this->DIM);
		MathLibrary::get_quaternion(cf_r_low3, cf_r_up3, cf_lower_angle3, cf_lower_normal3, this->DIM);
		MathLibrary::get_quaternion(cf_r_low4, cf_r_up4, cf_lower_angle4, cf_lower_normal4, this->DIM);


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

		cblas_dcopy(this->DIM, cf_r_up1, 1, cf_upper_force1, 1);
		scale = this->upper_spring_stiffness[0] * (1.0 - this->upper_spring_length[0] * inv_norm_r_up1);
		cblas_dscal(this->DIM, scale, cf_upper_force1, 1);
		cblas_dcopy(this->DIM, cf_r_up2, 1, cf_upper_force2, 1);
		scale = this->upper_spring_stiffness[1] * (1.0 - this->upper_spring_length[1] * inv_norm_r_up2);
		cblas_dscal(this->DIM, scale, cf_upper_force2, 1);
		cblas_dcopy(this->DIM, cf_r_up3, 1, cf_upper_force3, 1);
		scale = this->upper_spring_stiffness[2] * (1.0 - this->upper_spring_length[2] * inv_norm_r_up3);
		cblas_dscal(this->DIM, scale, cf_upper_force3, 1);
		cblas_dcopy(this->DIM, cf_r_up4, 1, cf_upper_force4, 1);
		scale = this->upper_spring_stiffness[3] * (1.0 - this->upper_spring_length[3] * inv_norm_r_up4);
		cblas_dscal(this->DIM, scale, cf_upper_force4, 1);

		cblas_dcopy(this->DIM, cf_r_low1, 1, cf_lower_force1, 1);
		scale = this->lower_spring_stiffness[0] * (1.0 - this->lower_spring_length[0] * inv_norm_r_low1);
		cblas_dscal(this->DIM, scale, cf_lower_force1, 1);
		cblas_dcopy(this->DIM, cf_r_low2, 1, cf_lower_force2, 1);
		scale = this->lower_spring_stiffness[1] * (1.0 - this->lower_spring_length[1] * inv_norm_r_low2);
		cblas_dscal(this->DIM, scale, cf_lower_force2, 1);
		cblas_dcopy(this->DIM, cf_r_low3, 1, cf_lower_force3, 1);
		scale = this->lower_spring_stiffness[2] * (1.0 - this->lower_spring_length[2] * inv_norm_r_low3);
		cblas_dscal(this->DIM, scale, cf_lower_force3, 1);
		cblas_dcopy(this->DIM, cf_r_low4, 1, cf_lower_force4, 1);
		scale = this->lower_spring_stiffness[3] * (1.0 - this->lower_spring_length[3] * inv_norm_r_low4);
		cblas_dscal(this->DIM, scale, cf_lower_force4, 1);

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
		cblas_dgemv(CblasRowMajor, CblasNoTrans, this->DIM, this->DIM, 1.0, this->r1_tilda, this->DIM, wc_, 1, 0.0, cf_temp, 1);
		cblas_dcopy(this->DIM, vc_, 1, cf_upper_dampf1, 1);
		cblas_dgemv(CblasRowMajor, CblasNoTrans, this->DIM, this->DIM, -1.0, cf_C_cN, this->DIM, cf_temp, 1, 1.0, cf_upper_dampf1, 1);
		// dot((vc - C_cN' * (r1_tilda * wc)), r_up1)
		//std::cout << cblas_ddot(this->DIM, upper_dampf1, 1, r_up1, 1) << std::endl;
		//std::cout << MathLibrary::dot_product<T>(upper_dampf1, r_up1, this->DIM) << std::endl;
		scale = cblas_ddot(mkl_DIM, cf_upper_dampf1, mkl_incx, cf_r_up1, mkl_incy);
		//scale = MathLibrary::dot_product<T>(upper_dampf1, r_up1, this->DIM);

		// dot((vc - C_cN' * (r1_tilda * wc)), r_up1) - dot(vw1, r_up1)
		scale -= cblas_ddot(mkl_DIM, vw1_, mkl_incx, cf_r_up1, mkl_incx);
		//scale -= MathLibrary::dot_product<T>(vw1_, r_up1, this->DIM);
		// (dot((vc - C_cN' * (r1_tilda * wc)), r_up1) - dot(vw1, r_up1))* inv_norm_r_up1 * inv_norm_r_up1
		scale = scale * inv_norm_r_up1 * inv_norm_r_up1 * this->upper_spring_damping[0];
		cblas_dcopy(this->DIM, cf_r_up1, 1, cf_upper_dampf1, 1);
		cblas_dscal(this->DIM, scale, cf_upper_dampf1, 1);


		//// upper_dampf2
		// compute: vc - C_cN' * (r2_tilda * wc))
		cblas_dgemv(CblasRowMajor, CblasNoTrans, this->DIM, this->DIM, 1.0, this->r2_tilda, this->DIM, wc_, 1, 0.0, cf_temp, 1);
		cblas_dcopy(this->DIM, vc_, 1, cf_upper_dampf2, 1);
		cblas_dgemv(CblasRowMajor, CblasNoTrans, this->DIM, this->DIM, -1.0, cf_C_cN, this->DIM, cf_temp, 1, 1.0, cf_upper_dampf2, 1);
		// dot((vc - C_cN' * (r2_tilda * wc)), r_up2)
		scale = cblas_ddot(mkl_DIM, cf_upper_dampf2, mkl_incx, cf_r_up2, mkl_incy);
		//scale = MathLibrary::dot_product<T>(upper_dampf2, r_up2, this->DIM);
		// dot((vc - C_cN' * (r2_tilda * wc)), r_up2) - dot(vw2, r_up2)
		scale -= cblas_ddot(mkl_DIM, vw2_, mkl_incx, cf_r_up2, mkl_incy);
		//scale -= MathLibrary::dot_product<T>(vw2_, r_up2, this->DIM);
		// (dot((vc - C_cN' * (r2_tilda * wc)), r_up2) - dot(vw2, r_up2))* inv_norm_r_up2 * inv_norm_r_up2
		scale = scale * inv_norm_r_up2 * inv_norm_r_up2 * this->upper_spring_damping[1];
		cblas_dcopy(this->DIM, cf_r_up2, 1, cf_upper_dampf2, 1);
		cblas_dscal(this->DIM, scale, cf_upper_dampf2, 1);

		//// upper_dampf3
		// compute: vc - C_cN' * (r3_tilda * wc))
		cblas_dgemv(CblasRowMajor, CblasNoTrans, this->DIM, this->DIM, 1.0, this->r3_tilda, this->DIM, wc_, 1, 0.0, cf_temp, 1);
		cblas_dcopy(this->DIM, vc_, 1, cf_upper_dampf3, 1);
		cblas_dgemv(CblasRowMajor, CblasNoTrans, this->DIM, this->DIM, -1.0, cf_C_cN, this->DIM, cf_temp, 1, 1.0, cf_upper_dampf3, 1);
		// dot((vc - C_cN' * (r3_tilda * wc)), r_up3)
		scale = cblas_ddot(mkl_DIM, cf_upper_dampf3, mkl_incx, cf_r_up3, mkl_incy);
		//scale = MathLibrary::dot_product<T>(upper_dampf3, r_up3, this->DIM);
		// dot((vc - C_cN' * (r3_tilda * wc)), r_up3) - dot(vw3, r_up3)
		scale -= cblas_ddot(mkl_DIM, vw3_, mkl_incx, cf_r_up3, mkl_incy);
		//scale -= MathLibrary::dot_product<T>(vw3_, r_up3, this->DIM);
		// (dot((vc - C_cN' * (r3_tilda * wc)), r_up3) - dot(vw3, r_up3))* inv_norm_r_up3 * inv_norm_r_up3
		scale = scale * inv_norm_r_up3 * inv_norm_r_up3 * this->upper_spring_damping[2];
		cblas_dcopy(this->DIM, cf_r_up3, 1, cf_upper_dampf3, 1);
		cblas_dscal(this->DIM, scale, cf_upper_dampf3, 1);

		//// upper_dampf4
		// compute: vc - C_cN' * (r4_tilda * wc))
		cblas_dgemv(CblasRowMajor, CblasNoTrans, this->DIM, this->DIM, 1.0, this->r4_tilda, this->DIM, wc_, 1, 0.0, cf_temp, 1);
		cblas_dcopy(this->DIM, vc_, 1, cf_upper_dampf4, 1);
		cblas_dgemv(CblasRowMajor, CblasNoTrans, this->DIM, this->DIM, -1.0, cf_C_cN, this->DIM, cf_temp, 1, 1.0, cf_upper_dampf4, 1);
		// dot((vc - C_cN' * (r4_tilda * wc)), r_up4)
		scale = cblas_ddot(mkl_DIM, cf_upper_dampf4, mkl_incx, cf_r_up4, mkl_incy);
		//scale = MathLibrary::dot_product<T>(upper_dampf4, r_up4, this->DIM);
		// dot((vc - C_cN' * (r4_tilda * wc)), r_up4) - dot(vw4, r_up4)
		scale -= cblas_ddot(mkl_DIM, vw4_, mkl_incx, cf_r_up4, mkl_incy);
		//scale -= MathLibrary::dot_product<T>(r_up4, vw4_, this->DIM);
		// (dot((vc - C_cN' * (r4_tilda * wc)), r_up4) - dot(vw4, r_up4))* inv_norm_r_up4 * inv_norm_r_up4
		scale = scale * inv_norm_r_up4 * inv_norm_r_up4 * this->upper_spring_damping[3];
		cblas_dcopy(this->DIM, cf_r_up4, 1, cf_upper_dampf4, 1);
		cblas_dscal(this->DIM, scale, cf_upper_dampf4, 1);

		//// lower_dampf1
		// dot(vw1, r_low1)
		scale = cblas_ddot(mkl_DIM, vw1_, mkl_incx, cf_r_low1, mkl_incy);
		/*scale = MathLibrary::dot_product<T>(vw1_, r_low1, this->DIM);*/
		// (dot(vw1, r_low1) - dot(vt1, r_low1))
		scale -= cblas_ddot(mkl_DIM, vt1_, mkl_incx, cf_r_low1, mkl_incx);
		//scale -= MathLibrary::dot_product<T>(vt1_, r_low1, this->DIM);
		//(dot(vw1, r_low1) - dot(vt1, r_low1)) * inv_norm_r_low1 * inv_norm_r_low1
		scale = scale * inv_norm_r_low1 * inv_norm_r_low1 * this->lower_spring_damping[0];
		cblas_dcopy(this->DIM, cf_r_low1, 1, cf_lower_dampf1, 1);
		cblas_dscal(this->DIM, scale, cf_lower_dampf1, 1);

		//// lower_dampf2
		// dot(vw2, r_low2)
		scale = cblas_ddot(mkl_DIM, vw2_, mkl_incx, cf_r_low2, mkl_incy);
		//scale = MathLibrary::dot_product<T>(vw2_, r_low2, this->DIM);
		// (dot(vw2, r_low2) - dot(vt2, r_low2))
		scale -= cblas_ddot(mkl_DIM, vt2_, mkl_incx, cf_r_low2, mkl_incy);
		//scale -= MathLibrary::dot_product<T>(vt2_, r_low2, this->DIM);
		//(dot(vw2, r_low2) - dot(vt2, r_low2)) * inv_norm_r_low2 * inv_norm_r_low2
		scale = scale * inv_norm_r_low2 * inv_norm_r_low2 * this->lower_spring_damping[1];
		cblas_dcopy(this->DIM, cf_r_low2, 1, cf_lower_dampf2, 1);
		cblas_dscal(this->DIM, scale, cf_lower_dampf2, 1);

		//// lower_dampf3
		// dot(vw3, r_low3)
		scale = cblas_ddot(mkl_DIM, vw3_, mkl_incx, cf_r_low3, mkl_incy);
		//scale = MathLibrary::dot_product<T>(vw3_, r_low3, this->DIM);
		// (dot(vw3, r_low3) - dot(vt3, r_low3))
		scale -= cblas_ddot(mkl_DIM, vt3_, mkl_incx, cf_r_low3, mkl_incy);
		//scale -= MathLibrary::dot_product<T>(vt3_, r_low3, this->DIM);
		//(dot(vw3, r_low3) - dot(vt3, r_low3)) * inv_norm_r_low3 * inv_norm_r_low3
		scale = scale * inv_norm_r_low3 * inv_norm_r_low3 * this->lower_spring_damping[2];
		cblas_dcopy(this->DIM, cf_r_low3, 1, cf_lower_dampf3, 1);
		cblas_dscal(this->DIM, scale, cf_lower_dampf3, 1);

		//// lower_dampf4
		// dot(vw4, r_low4)
		scale = cblas_ddot(mkl_DIM, vw4_, mkl_incx, cf_r_low4, mkl_incy);
		//scale = MathLibrary::dot_product<T>(vw4_, r_low4, this->DIM);
		// (dot(vw4, r_low4) - dot(vt4, r_low4))
		scale -= cblas_ddot(mkl_DIM, vt4_, mkl_incx, cf_r_low4, mkl_incy);
		//scale -= MathLibrary::dot_product<T>(vt4_, r_low4, this->DIM);
		//(dot(vw4, r_low4) - dot(vt4, r_low4)) * inv_norm_r_low4 * inv_norm_r_low4
		scale = scale * inv_norm_r_low4 * inv_norm_r_low4 * this->lower_spring_damping[3];
		cblas_dcopy(this->DIM, cf_r_low4, 1, cf_lower_dampf4, 1);
		cblas_dscal(this->DIM, scale, cf_lower_dampf4, 1);

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
		scale = (this->upper_rotational_stiffness[0]) * (*cf_upper_angle1);
		cblas_dcopy(this->DIM, cf_upper_normal1, 1, cf_upper_S1, 1);
		cblas_dscal(this->DIM, scale, cf_upper_S1, 1);


		// upper_S2 = upper_rotational_stiffness(2) * upper_angle2 * upper_normal2;
		scale = (this->upper_rotational_stiffness[1]) * (*cf_upper_angle2);
		cblas_dcopy(this->DIM, cf_upper_normal2, 1, cf_upper_S2, 1);
		cblas_dscal(this->DIM, scale, cf_upper_S2, 1);

		// upper_S3 = upper_rotational_stiffness(3) * upper_angle3 * upper_normal3;
		scale = (this->upper_rotational_stiffness[2]) * (*cf_upper_angle3);
		cblas_dcopy(this->DIM, cf_upper_normal3, 1, cf_upper_S3, 1);
		cblas_dscal(this->DIM, scale, cf_upper_S3, 1);

		// upper_S4 = upper_rotational_stiffness(4) * upper_angle4 * upper_normal4;
		scale = (this->upper_rotational_stiffness[3]) * (*cf_upper_angle4);
		cblas_dcopy(this->DIM, cf_upper_normal4, 1, cf_upper_S4, 1);
		cblas_dscal(this->DIM, scale, cf_upper_S4, 1);

		// lower_S1 = lower_rotational_stiffness(1) * lower_angle1 * lower_normal1
		scale = (this->lower_rotational_stiffness[0]) * (*cf_lower_angle1);
		cblas_dcopy(this->DIM, cf_lower_normal1, 1, cf_lower_S1, 1);
		cblas_dscal(this->DIM, scale, cf_lower_S1, 1);

		// lower_S2 = lower_rotational_stiffness(2) * lower_angle2 * lower_normal2
		scale = (this->lower_rotational_stiffness[1]) * (*cf_lower_angle2);
		cblas_dcopy(this->DIM, cf_lower_normal2, 1, cf_lower_S2, 1);
		cblas_dscal(this->DIM, scale, cf_lower_S2, 1);

		// lower_S3 = lower_rotational_stiffness(3) * lower_angle3 * lower_normal3
		scale = (this->lower_rotational_stiffness[2]) * (*cf_lower_angle3);
		cblas_dcopy(this->DIM, cf_lower_normal3, 1, cf_lower_S3, 1);
		cblas_dscal(this->DIM, scale, cf_lower_S3, 1);

		// lower_S4 = lower_rotational_stiffness(4) * lower_angle4 * lower_normal4
		scale = (this->lower_rotational_stiffness[3]) * (*cf_lower_angle4);
		cblas_dcopy(this->DIM, cf_lower_normal4, 1, cf_lower_S4, 1);
		cblas_dscal(this->DIM, scale, cf_lower_S4, 1);

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

		T scale_u1, scale_u2, scale_u3, scale_u4;
		// lower_rot_force1 = -cross( lower_S1, r_low1) / (r_low1'*r_low1);
		scale = -1.0 / cblas_ddot(mkl_DIM, cf_r_low1, mkl_incx, cf_r_low1, mkl_incy);
		//scale = -1.0 / MathLibrary::dot_product<T>(r_low1, r_low1, this->DIM);
		MathLibrary::crossProduct(cf_lower_S1, cf_r_low1, cf_lower_rot_force1);
		cblas_dscal(this->DIM, scale, cf_lower_rot_force1, 1);

		// lower_rot_force2 = -cross( lower_S2, r_low2) / (r_low2'*r_low2);
		scale = -1.0 / cblas_ddot(mkl_DIM, cf_r_low2, mkl_incx, cf_r_low2, mkl_incy);
		//scale = -1.0 / MathLibrary::dot_product<T>(r_low2, r_low2, this->DIM);
		MathLibrary::crossProduct(cf_lower_S2, cf_r_low2, cf_lower_rot_force2);
		cblas_dscal(this->DIM, scale, cf_lower_rot_force2, 1);

		// lower_rot_force3 = -cross( lower_S3, r_low3) / (r_low3'*r_low3);
		scale = -1.0 / cblas_ddot(mkl_DIM, cf_r_low3, mkl_incx, cf_r_low3, mkl_incy);
		//scale = -1.0 / MathLibrary::dot_product<T>(r_low3, r_low3, this->DIM);
		MathLibrary::crossProduct(cf_lower_S3, cf_r_low3, cf_lower_rot_force3);
		cblas_dscal(this->DIM, scale, cf_lower_rot_force3, 1);

		// lower_rot_force4 = -cross( lower_S4, r_low4) / (r_low4'*r_low4);
		scale = -1.0 / cblas_ddot(mkl_DIM, cf_r_low4, mkl_incx, cf_r_low4, mkl_incx);
		//scale = -1.0 / MathLibrary::dot_product<T>(r_low4, r_low4, this->DIM);
		MathLibrary::crossProduct(cf_lower_S4, cf_r_low4, cf_lower_rot_force4);
		cblas_dscal(this->DIM, scale, cf_lower_rot_force4, 1);

		// upper_rot_force1 = -cross( upper_S1, r_up1) / (r_up1'*r_up1);
		scale_u1 = -1.0 / cblas_ddot(mkl_DIM, cf_r_up1, mkl_incx, cf_r_up1, mkl_incy);
		//scale_u1 = -1.0 / MathLibrary::dot_product<T>(r_up1, r_up1, this->DIM);
		MathLibrary::crossProduct(cf_upper_S1, cf_r_up1, cf_upper_rot_force1);
		cblas_dscal(this->DIM, scale_u1, cf_upper_rot_force1, 1);

		// upper_rot_force2 = -cross( upper_S2, r_up2) / (r_up2'*r_up2);
		scale_u2 = -1.0 / cblas_ddot(mkl_DIM, cf_r_up2, mkl_incx, cf_r_up2, mkl_incy);
		//scale_u2 = -1.0 / MathLibrary::dot_product<T>(r_up2, r_up2, this->DIM);
		MathLibrary::crossProduct(cf_upper_S2, cf_r_up2, cf_upper_rot_force2);
		cblas_dscal(this->DIM, scale_u2, cf_upper_rot_force2, 1);

		// upper_rot_force3 = -cross( upper_S3, r_up3) / (r_up3'*r_up3);
		scale_u3 = -1.0 / cblas_ddot(mkl_DIM, cf_r_up3, mkl_incx, cf_r_up3, mkl_incy);
		//scale_u3 = -1.0 / MathLibrary::dot_product<T>(r_up3, r_up3, this->DIM);
		MathLibrary::crossProduct(cf_upper_S3, cf_r_up3, cf_upper_rot_force3);
		cblas_dscal(this->DIM, scale_u3, cf_upper_rot_force3, 1);

		// upper_rot_force4 = -cross( upper_S4, r_up4) / (r_up4'*r_up4);
		scale_u4 = -1.0 / cblas_ddot(mkl_DIM, cf_r_up4, mkl_incx, cf_r_up4, mkl_incy);
		//scale_u4 = -1.0 / MathLibrary::dot_product<T>(r_up4, r_up4, this->DIM);
		MathLibrary::crossProduct(cf_upper_S4, cf_r_up4, cf_upper_rot_force4);
		cblas_dscal(this->DIM, scale_u4, cf_upper_rot_force4, 1);

		// car_rot_force1 = -cross( lower_S1, r_up1) / (r_up1'*r_up1);
		MathLibrary::crossProduct(cf_lower_S1, cf_r_up1, cf_car_rot_force1);
		cblas_dscal(this->DIM, scale_u1, cf_car_rot_force1, 1);

		// car_rot_force2 = -cross( lower_S2, r_up2) / (r_up2'*r_up2);
		MathLibrary::crossProduct(cf_lower_S2, cf_r_up2, cf_car_rot_force2);
		cblas_dscal(this->DIM, scale_u2, cf_car_rot_force2, 1);

		// car_rot_force3 = -cross( lower_S3, r_up3) / (r_up3'*r_up3);
		MathLibrary::crossProduct(cf_lower_S3, cf_r_up3, cf_car_rot_force3);
		cblas_dscal(this->DIM, scale_u3, cf_car_rot_force3, 1);

		// car_rot_force4 = -cross( lower_S4, r_up4) / (r_up4'*r_up4);
		MathLibrary::crossProduct(cf_lower_S4, cf_r_up4, cf_car_rot_force4);
		cblas_dscal(this->DIM, scale_u4, cf_car_rot_force4, 1);

		// sum_car_force1 = car_rot_force1 - upper_force1 - upper_dampf1 - upper_rot_force1;
		cblas_dcopy(this->DIM, cf_car_rot_force1, 1, cf_sum_car_force1, 1);
		cblas_daxpy(this->DIM, -1.0, cf_upper_force1, 1, cf_sum_car_force1, 1);
		cblas_daxpy(this->DIM, -1.0, cf_upper_dampf1, 1, cf_sum_car_force1, 1);
		cblas_daxpy(this->DIM, -1.0, cf_upper_rot_force1, 1, cf_sum_car_force1, 1);

		// sum_car_force2 = car_rot_force2 - upper_force2 - upper_dampf2 - upper_rot_force2;
		cblas_dcopy(this->DIM, cf_car_rot_force2, 1, cf_sum_car_force2, 1);
		cblas_daxpy(this->DIM, -1.0, cf_upper_force2, 1, cf_sum_car_force2, 1);
		cblas_daxpy(this->DIM, -1.0, cf_upper_dampf2, 1, cf_sum_car_force2, 1);
		cblas_daxpy(this->DIM, -1.0, cf_upper_rot_force2, 1, cf_sum_car_force2, 1);

		// sum_car_force3 = car_rot_force3 - upper_force3 - upper_dampf3 - upper_rot_force3;
		cblas_dcopy(this->DIM, cf_car_rot_force3, 1, cf_sum_car_force3, 1);
		cblas_daxpy(this->DIM, -1.0, cf_upper_force3, 1, cf_sum_car_force3, 1);
		cblas_daxpy(this->DIM, -1.0, cf_upper_dampf3, 1, cf_sum_car_force3, 1);
		cblas_daxpy(this->DIM, -1.0, cf_upper_rot_force3, 1, cf_sum_car_force3, 1);

		// sum_car_force4 = car_rot_force4 - upper_force4 - upper_dampf4 - upper_rot_force4;
		cblas_dcopy(this->DIM, cf_car_rot_force4, 1, cf_sum_car_force4, 1);
		cblas_daxpy(this->DIM, -1.0, cf_upper_force4, 1, cf_sum_car_force4, 1);
		cblas_daxpy(this->DIM, -1.0, cf_upper_dampf4, 1, cf_sum_car_force4, 1);
		cblas_daxpy(this->DIM, -1.0, cf_upper_rot_force4, 1, cf_sum_car_force4, 1);

		/*
		get external forces
		for the legs
		local_FW1 = [0; FW(1); 0];      %in global basis
		local_FW2 = [0; FW(2); 0];
		local_FW3 = [0; FW(3); 0];
		local_FW4 = [0; FW(4); 0];

		local_FT1 = [0; FT(1); 0];      %in global basis
		local_FT2 = [0; FT(2); 0];
		local_FT3 = [0; FT(3); 0];
		local_FT4 = [0; FT(4); 0];
		*/

		// Calculate the sum of all forces in the tyre (and keep it in local_FR)

		// local_FR1 = lower_force1 + lower_dampf1 + local_FT1 + local_FR1 + lower_rot_force1; ...          %vt1_dot
		cblas_dcopy(this->DIM, cf_lower_force1, 1, cf_local_FR1, 1);
		cblas_daxpy(this->DIM, 1.0, cf_lower_dampf1, 1, cf_local_FR1, 1);
		//	cblas_daxpy(this->DIM, 1.0, local_FT1, 1, cf_local_FR1, 1); not implemented since not used now
		cblas_daxpy(this->DIM, 1.0, cf_local_FR1, 1, cf_local_FR1, 1);
		cblas_daxpy(this->DIM, 1.0, cf_lower_rot_force1, 1, cf_local_FR1, 1);

		// local_FR2 = lower_force2 + lower_dampf2 + local_FT2 + local_FR2 + lower_rot_force2; ...          %vt2_dot
		cblas_dcopy(this->DIM, cf_lower_force2, 1, cf_local_FR2, 1);
		cblas_daxpy(this->DIM, 1.0, cf_lower_dampf2, 1, cf_local_FR2, 1);
		//	cblas_daxpy(this->DIM, 1.0, local_FT2, 1, cf_local_FR2, 1); not implemented since not used now
		cblas_daxpy(this->DIM, 1.0, cf_local_FR2, 1, cf_local_FR2, 1);
		cblas_daxpy(this->DIM, 1.0, cf_lower_rot_force2, 1, cf_local_FR2, 1);

		// local_FR3 = lower_force3 + lower_dampf3 + local_FT3 + local_FR3 + lower_rot_force3; ...          %vt3_dot
		cblas_dcopy(this->DIM, cf_lower_force3, 1, cf_local_FR3, 1);
		cblas_daxpy(this->DIM, 1.0, cf_lower_dampf3, 1, cf_local_FR3, 1);
		//	cblas_daxpy(this->DIM, 1.0, local_FT3, 1, cf_local_FR3, 1); not implemented since not used now
		cblas_daxpy(this->DIM, 1.0, cf_local_FR3, 1, cf_local_FR3, 1);
		cblas_daxpy(this->DIM, 1.0, cf_lower_rot_force3, 1, cf_local_FR3, 1);

		// local_FR4 = lower_force4 + lower_dampf4 + local_FT4 + local_FR4 + lower_rot_force4];             %vt4_dot
		cblas_dcopy(this->DIM, cf_lower_force4, 1, cf_local_FR4, 1);
		cblas_daxpy(this->DIM, 1.0, cf_lower_dampf4, 1, cf_local_FR4, 1);
		//	cblas_daxpy(this->DIM, 1.0, local_FT4, 1, cf_local_FR4, 1); not implemented since not used now
		cblas_daxpy(this->DIM, 1.0, cf_local_FR4, 1, cf_local_FR4, 1);
		cblas_daxpy(this->DIM, 1.0, cf_lower_rot_force4, 1, cf_local_FR4, 1);


		if (boundary_conditions == FIXED) {
			get_fixed_road_force(cf_local_FR1, cf_local_FR2, cf_local_FR3, cf_local_FR4);
		}
		else if (boundary_conditions == NONFIXED) {
			get_nonfixed_road_force(cf_local_FR1, cf_local_FR2, cf_local_FR3, cf_local_FR4);
		}
		else if (boundary_conditions == CIRCULAR) {
			get_circular_road_force(cf_local_FR1, vt1_, mass_tyre_rr, pt1_);
			get_circular_road_force(cf_local_FR2, vt2_, mass_tyre_rl, pt2_);
			get_circular_road_force(cf_local_FR3, vt3_, mass_tyre_fl, pt3_);
			get_circular_road_force(cf_local_FR4, vt4_, mass_tyre_fr, pt4_);
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
								r4_tilda * (C_cN * sum_car_force4) + ...
							   -C_cN * upper_S1 - C_cN * upper_S2 - C_cN * upper_S3 - C_cN * upper_S4 + ...     % ??from the rotational spring
							   -get_tilda(wc) * Hc + Tc;                       % from angular momentum and external torques

		*/


		// Hc = A(1:3, 1:3) * wc;
		cblas_dgemv(CblasRowMajor, CblasNoTrans, this->DIM, this->DIM, 1.0, this->Ic, this->DIM, wc_, 1, 0.0, cf_Hc, 1);
		MathLibrary::get_tilda<T>(wc_, cf_wc_tilda);
		cblas_dcopy(this->DIM, cf_Hc, 1, cf_temp, 1);
		cblas_dgemv(CblasRowMajor, CblasNoTrans, this->DIM, this->DIM, -1.0, cf_wc_tilda, this->DIM, cf_temp, 1, 0.0, cf_Hc, 1);
		cblas_dcopy(this->DIM, cf_Hc, 1, cf_sum_torque_spring_car, 1);
		cblas_daxpy(this->DIM, 1.0, cf_Tc, 1, cf_sum_torque_spring_car, 1);
		cblas_dgemv(CblasRowMajor, CblasTrans, this->DIM, this->DIM, -1.0, cf_C_cN, this->DIM, cf_upper_S4, 1, 1.0, cf_sum_torque_spring_car, 1);
		cblas_dgemv(CblasRowMajor, CblasTrans, this->DIM, this->DIM, -1.0, cf_C_cN, this->DIM, cf_upper_S3, 1, 1.0, cf_sum_torque_spring_car, 1);
		cblas_dgemv(CblasRowMajor, CblasTrans, this->DIM, this->DIM, -1.0, cf_C_cN, this->DIM, cf_upper_S2, 1, 1.0, cf_sum_torque_spring_car, 1);
		cblas_dgemv(CblasRowMajor, CblasTrans, this->DIM, this->DIM, -1.0, cf_C_cN, this->DIM, cf_upper_S1, 1, 1.0, cf_sum_torque_spring_car, 1);
		cblas_dgemv(CblasRowMajor, CblasTrans, this->DIM, this->DIM, 1.0, cf_C_cN, this->DIM, cf_sum_car_force1, 1, 0.0, cf_temp, 1);
		cblas_dgemv(CblasRowMajor, CblasNoTrans, this->DIM, this->DIM, 1.0, this->r1_tilda, this->DIM, cf_temp, 1, 1.0, cf_sum_torque_spring_car, 1);
		cblas_dgemv(CblasRowMajor, CblasTrans, this->DIM, this->DIM, 1.0, cf_C_cN, this->DIM, cf_sum_car_force2, 1, 0.0, cf_temp, 1);
		cblas_dgemv(CblasRowMajor, CblasNoTrans, this->DIM, this->DIM, 1.0, this->r2_tilda, this->DIM, cf_temp, 1, 1.0, cf_sum_torque_spring_car, 1);
		cblas_dgemv(CblasRowMajor, CblasTrans, this->DIM, this->DIM, 1.0, cf_C_cN, this->DIM, cf_sum_car_force3, 1, 0.0, cf_temp, 1);
		cblas_dgemv(CblasRowMajor, CblasNoTrans, this->DIM, this->DIM, 1.0, this->r3_tilda, this->DIM, cf_temp, 1, 1.0, cf_sum_torque_spring_car, 1);
		cblas_dgemv(CblasRowMajor, CblasTrans, this->DIM, this->DIM, 1.0, cf_C_cN, this->DIM, cf_sum_car_force4, 1, 0.0, cf_temp, 1);
		cblas_dgemv(CblasRowMajor, CblasNoTrans, this->DIM, this->DIM, 1.0, this->r4_tilda, this->DIM, cf_temp, 1, 1.0, cf_sum_torque_spring_car, 1);

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
		cblas_dcopy(this->DIM, this->FC, 1, brem_start, 1);
		cblas_daxpy(this->DIM, 1.0, cf_sum_car_force1, 1, brem_start, 1);
		cblas_daxpy(this->DIM, 1.0, cf_sum_car_force2, 1, brem_start, 1);
		cblas_daxpy(this->DIM, 1.0, cf_sum_car_force3, 1, brem_start, 1);
		cblas_daxpy(this->DIM, 1.0, cf_sum_car_force4, 1, brem_start, 1);

		//  upper_force1 - lower_force1 + upper_dampf1 - lower_dampf1 + local_FW1 + upper_rot_force1 - car_rot_force1 - lower_rot_force1;
		brem_start += this->DIM;
		cblas_dcopy(this->DIM, cf_upper_force1, 1, brem_start, 1);
		cblas_daxpy(this->DIM, -1.0, cf_lower_force1, 1, brem_start, 1);
		cblas_daxpy(this->DIM, 1.0, cf_upper_dampf1, 1, brem_start, 1);
		cblas_daxpy(this->DIM, -1.0, cf_lower_dampf1, 1, brem_start, 1);
		//cblas_daxpy(this->DIM, 1.0, local_FW1, 1, brem_start, 1);
		cblas_daxpy(this->DIM, 1.0, cf_upper_rot_force1, 1, brem_start, 1);
		cblas_daxpy(this->DIM, -1.0, cf_car_rot_force1, 1, brem_start, 1);
		cblas_daxpy(this->DIM, -1.0, cf_lower_rot_force1, 1, brem_start, 1);

		// upper_force2 - lower_force2 + upper_dampf2 - lower_dampf2 + local_FW2 + upper_rot_force2 - car_rot_force2 - lower_rot_force2;
		brem_start += this->DIM;
		cblas_dcopy(this->DIM, cf_upper_force2, 1, brem_start, 1);
		cblas_daxpy(this->DIM, -1.0, cf_lower_force2, 1, brem_start, 1);
		cblas_daxpy(this->DIM, 1.0, cf_upper_dampf2, 1, brem_start, 1);
		cblas_daxpy(this->DIM, -1.0, cf_lower_dampf2, 1, brem_start, 1);
		//cblas_daxpy(this->DIM, 1.0, local_FW2, 1, brem_start, 1);
		cblas_daxpy(this->DIM, 1.0, cf_upper_rot_force2, 1, brem_start, 1);
		cblas_daxpy(this->DIM, -1.0, cf_car_rot_force2, 1, brem_start, 1);
		cblas_daxpy(this->DIM, -1.0, cf_lower_rot_force2, 1, brem_start, 1);

		// upper_force3 - lower_force3 + upper_dampf3 - lower_dampf3 + local_FW3 + upper_rot_force3 - car_rot_force3 - lower_rot_force3;
		brem_start += this->DIM;
		cblas_dcopy(this->DIM, cf_upper_force3, 1, brem_start, 1);
		cblas_daxpy(this->DIM, -1.0, cf_lower_force3, 1, brem_start, 1);
		cblas_daxpy(this->DIM, 1.0, cf_upper_dampf3, 1, brem_start, 1);
		cblas_daxpy(this->DIM, -1.0, cf_lower_dampf3, 1, brem_start, 1);
		//cblas_daxpy(this->DIM, 1.0, local_FW3, 1, brem_start, 1);
		cblas_daxpy(this->DIM, 1.0, cf_upper_rot_force3, 1, brem_start, 1);
		cblas_daxpy(this->DIM, -1.0, cf_car_rot_force3, 1, brem_start, 1);
		cblas_daxpy(this->DIM, -1.0, cf_lower_rot_force3, 1, brem_start, 1);

		//  upper_force4 - lower_force4 + upper_dampf4 - lower_dampf4 + local_FW4 + upper_rot_force4 - car_rot_force4 - lower_rot_force4;
		brem_start += this->DIM;
		cblas_dcopy(this->DIM, cf_upper_force4, 1, brem_start, 1);
		cblas_daxpy(this->DIM, -1.0, cf_lower_force4, 1, brem_start, 1);
		cblas_daxpy(this->DIM, 1.0, cf_upper_dampf4, 1, brem_start, 1);
		cblas_daxpy(this->DIM, -1.0, cf_lower_dampf4, 1, brem_start, 1);
		//cblas_daxpy(this->DIM, 1.0, local_FW4, 1, brem_start, 1);
		cblas_daxpy(this->DIM, 1.0, cf_upper_rot_force4, 1, brem_start, 1);
		cblas_daxpy(this->DIM, -1.0, cf_car_rot_force4, 1, brem_start, 1);
		cblas_daxpy(this->DIM, -1.0, cf_lower_rot_force4, 1, brem_start, 1);

		// local_FR1; ...					  %vt1_dot
		brem_start += this->DIM;
		cblas_dcopy(this->DIM, cf_local_FR1, 1, brem_start, 1);

		// local_FR2; ...			 	      %vt2_dot
		brem_start += this->DIM;
		cblas_dcopy(this->DIM, cf_local_FR2, 1, brem_start, 1);

		// local_FR3; ...			          %vt3_dot
		brem_start += this->DIM;
		cblas_dcopy(this->DIM, cf_local_FR3, 1, brem_start, 1);
	
		// local_FR4];			              %vt4_dot
		brem_start += this->DIM;
		cblas_dcopy(this->DIM, cf_local_FR4, 1, brem_start, 1);

		LAPACKE_dpotrs(LAPACK_ROW_MAJOR, 'L', this->DIM, 1, A_Ic, this->DIM, cf_sum_torque_spring_car, 1);
		cblas_dcopy(this->DIM, cf_sum_torque_spring_car, 1, f_, 1);
		T* start_next = f_ + this->DIM;

		MathLibrary::vector_elem_wise_product<T>(A_rem, cf_b_rem, start_next, 9 * this->DIM);
		start_next += 9 * this->DIM;

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
		cblas_dgemv(CblasRowMajor, CblasNoTrans, this->NUM_LEGS, this->DIM, 1.0, cf_Qc, this->DIM, wc_, 1, 0.0, cf_qc_dot, 1);

		cblas_dcopy(this->NUM_LEGS, cf_qc_dot, 1, start_next, 1);
		start_next += this->NUM_LEGS;
		// add vc to f vector
		cblas_dcopy(this->DIM, vc_, 1, start_next, 1);
		start_next += this->DIM;

		// add vw1 to f vector
		cblas_dcopy(this->DIM, vw1_, 1, start_next, 1);
		start_next += this->DIM;

		// add vw2 to f vector
		cblas_dcopy(this->DIM, vw2_, 1, start_next, 1);
		start_next += this->DIM;

		// add vw3 to f vector
		cblas_dcopy(this->DIM, vw3_, 1, start_next, 1);
		start_next += this->DIM;

		// add vw4 to f vector
		cblas_dcopy(this->DIM, vw4_, 1, start_next, 1);
		start_next += this->DIM;

		// add vt1 to f vector
		cblas_dcopy(this->DIM, vt1_, 1, start_next, 1);
		start_next += this->DIM;

		// add vt2 to f vector
		cblas_dcopy(this->DIM, vt2_, 1, start_next, 1);
		start_next += this->DIM;

		// add vt3 to f vector
		cblas_dcopy(this->DIM, vt3_, 1, start_next, 1);
		start_next += this->DIM;

		// add vt4 to f vector
		cblas_dcopy(this->DIM, vt4_, 1, start_next, 1);
		start_next += this->DIM;


	}

	void solve(T* solution_vector, T* complete_vector) {
		x_vector = (T*)mkl_calloc(solution_dim, sizeof(T), this->alignment);
		MathLibrary::get_tilda<T>(r1, r1_tilda);
		MathLibrary::get_tilda<T>(r2, r2_tilda);
		MathLibrary::get_tilda<T>(r3, r3_tilda);
		MathLibrary::get_tilda<T>(r4, r4_tilda);
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
		T* start_orient = x_vector + i * (this->DIM) + j * (this->NUM_LEGS);
		cblas_dcopy(this->NUM_LEGS, initial_orientation, 1, x_vector + i * (this->DIM) + j * (this->NUM_LEGS), 1);
		j++;

		// pcc
		cblas_dcopy(this->DIM, pcc, 1, x_vector + i * (this->DIM) + j * (this->NUM_LEGS), 1);
		i++;
		// pw1
		T* pw1 = x_vector + i * (this->DIM) + j * (this->NUM_LEGS);
		i++;
		// pw2
		T* pw2 = x_vector + i * (this->DIM) + j * (this->NUM_LEGS);
		i++;
		// pw3
		T* pw3 = x_vector + i * (this->DIM) + j * (this->NUM_LEGS);
		i++;
		// pw4
		T* pw4 = x_vector + i * (this->DIM) + j * (this->NUM_LEGS);
		i++;
		// pt1
		T* pt1 = x_vector + i * (this->DIM) + j * (this->NUM_LEGS);
		i++;
		// pt2
		T* pt2 = x_vector + i * (this->DIM) + j * (this->NUM_LEGS);
		i++;
		// pt3
		T* pt3 = x_vector + i * (this->DIM) + j * (this->NUM_LEGS);
		i++;
		// pt4
		T* pt4 = x_vector + i * (this->DIM) + j * (this->NUM_LEGS);


		get_initial_length(start_orient, r1, r2, r3, r4, pcc, initial_upper_spring_length, initial_lower_spring_length, pw1, pw2, pw3, pw4, pt1, pt2, pt3, pt4);

		// overwrites the initial velocity values
		if (boundary_conditions == CIRCULAR)
			circular_path_initialization(vc, vw1, vw2, vw3, vw4, vt1, vt2, vt3, vt4, initial_angular_velocity, pcc, pt1, pt2, pt3, pt4, radius_circular_path);

		i = 0;
		j = 0;
		// wc
		cblas_dcopy(this->DIM, initial_angular_velocity, 1, x_vector + i * (this->DIM) + j * (this->NUM_LEGS), 1);
		i++;
		// vc
		cblas_dcopy(this->DIM, vc, 1, x_vector + i * (this->DIM) + j * (this->NUM_LEGS), 1);
		i++;
		// vw1
		cblas_dcopy(this->DIM, vw1, 1, x_vector + i * (this->DIM) + j * (this->NUM_LEGS), 1);
		i++;
		// vw2
		cblas_dcopy(this->DIM, vw2, 1, x_vector + i * (this->DIM) + j * (this->NUM_LEGS), 1);
		i++;
		// vw3
		cblas_dcopy(this->DIM, vw3, 1, x_vector + i * (this->DIM) + j * (this->NUM_LEGS), 1);
		i++;
		// vw4
		cblas_dcopy(this->DIM, vw4, 1, x_vector + i * (this->DIM) + j * (this->NUM_LEGS), 1);
		i++;
		// vt1
		cblas_dcopy(this->DIM, vt1, 1, x_vector + i * (this->DIM) + j * (this->NUM_LEGS), 1);
		i++;
		// vt2
		cblas_dcopy(this->DIM, vt2, 1, x_vector + i * (this->DIM) + j * (this->NUM_LEGS), 1);
		i++;
		// vt3
		cblas_dcopy(this->DIM, vt3, 1, x_vector + i * (this->DIM) + j * (this->NUM_LEGS), 1);
		i++;
		// vt4
		cblas_dcopy(this->DIM, vt4, 1, x_vector + i * (this->DIM) + j * (this->NUM_LEGS), 1);
		i++;

		j = 0;
		i = 0;
		while (i < this->DIM) {
			A_rem[j] = 1.0 / this->mass;
			i++;
			j++;
		}
		i = 0;
		while (i < this->DIM) {
			A_rem[j] = 1.0 / mass_wheel[0];
			i++;
			j++;
		}
		i = 0;
		while (i < this->DIM) {
			A_rem[j] = 1.0 / mass_wheel[1];
			i++;
			j++;
		}
		i = 0;
		while (i < this->DIM) {
			A_rem[j] = 1.0 / mass_wheel[2];
			i++;
			j++;
		}
		i = 0;
		while (i < this->DIM) {
			A_rem[j] = 1.0 / mass_wheel[3];
			i++;
			j++;
		}
		i = 0;
		while (i < this->DIM) {
			A_rem[j] = 1.0 / mass_tyre[0];
			i++;
			j++;
		}
		i = 0;
		while (i < this->DIM) {
			A_rem[j] = 1.0 / mass_tyre[1];
			i++;
			j++;
		}
		i = 0;
		while (i < this->DIM) {
			A_rem[j] = 1.0 / mass_tyre[2];
			i++;
			j++;
		}
		i = 0;
		while (i < this->DIM) {
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
		cblas_dcopy(this->solution_dim, start, 1, solution_vector, 1);
		//	std::cout << "Solution copied!\n" << std::endl;

		mkl_free(x_vector);
	}

	/*
	Initializes the time iteration and handles the numerical scheme
	*/
	void solve(T* solution_vector) {
		// From the formulation we have 61 dimensions in the solution vector
		size_t solution_size = (this->num_iter + 1) * this->solution_dim;
		T* complete_vector = (T*)mkl_calloc(solution_size, sizeof(T), this->alignment);
		solve(solution_vector, complete_vector);
		mkl_free(complete_vector);
		
	}

	/*
	Return the vector alignment (of the system)
	*/
	size_t get_alignment() {
		return this->alignment;
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
		std::cout << "MBD: angular velocity w=\n\t["<<sln[0] << "\n\t " << sln[2] << "\n\t " << sln[1] << "]" << std::endl;
		std::cout << "MBD: car body velocity vc=\n\t[" << sln[3] << "\n\t " << sln[5] << "\n\t " << sln[4] << "]" << std::endl;
		std::cout << "MBD: rear-right wheel velocity vw1=\n\t[" << sln[6] << "\n\t " << sln[8] << "\n\t " << sln[7] << "]" << std::endl;
		std::cout << "MBD: rear-left wheel velocity vw2=\n\t[" << sln[9] << "\n\t " << sln[11] << "\n\t " << sln[10] << "]" << std::endl;
		std::cout << "MBD: front-left wheel velocity vw3=\n\t[" << sln[12] << "\n\t " << sln[14] << "\n\t " << sln[13] << "]" << std::endl;
		std::cout << "MBD: front-right wheel velocity vw4=\n\t[" << sln[15] << "\n\t " << sln[17] << "\n\t " << sln[16] << "]" << std::endl;
		std::cout << "MBD: rear-right tyre velocity vt1=\n\t[" << sln[18] << "\n\t " << sln[20] << "\n\t " << sln[19] << "]" << std::endl;
		std::cout << "MBD: rear-left tyre velocity vt2=\n\t[" << sln[21] << "\n\t " << sln[23] << "\n\t " << sln[22] << "]" << std::endl;
		std::cout << "MBD: front-left tyre velocity vt3=\n\t[" << sln[24] << "\n\t " << sln[26] << "\n\t " << sln[25] << "]" << std::endl;
		std::cout << "MBD: front-right tyre velocity vt4=\n\t[" << sln[27] << "\n\t " << sln[29] << "\n\t " << sln[28] << "]" << std::endl;
		std::cout << "MBD: orientation q=\n\t[" << sln[30] << "\n\t " << sln[31] << "\n\t " << sln[32] << "\n\t " << sln[33] << "]" << std::endl;
		std::cout << "MBD: car body position pc=\n\t[" << sln[34] << "\n\t " << sln[36] << "\n\t " << sln[35] << "]" << std::endl;
		std::cout << "MBD: rear-right wheel position pw1=\n\t[" << sln[37] << "\n\t " << sln[39] << "\n\t " << sln[38] << "]" << std::endl;
		std::cout << "MBD: rear-left wheel position pw2=\n\t[" << sln[40] << "\n\t " << sln[42] << "\n\t " << sln[41] << "]" << std::endl;
		std::cout << "MBD: front-left wheel position pw3=\n\t[" << sln[43] << "\n\t " << sln[45] << "\n\t " << sln[44] << "]" << std::endl;
		std::cout << "MBD: front-right wheel position pw4=\n\t[" << sln[46] << "\n\t " << sln[48] << "\n\t " << sln[47] << "]" << std::endl;
		std::cout << "MBD: rear-right tyre position pt1=\n\t[" << sln[49] << "\n\t " << sln[51] << "\n\t " << sln[50] << "]" << std::endl;
		std::cout << "MBD: rear-left tyre position pt2=\n\t[" << sln[52] << "\n\t " << sln[54] << "\n\t " << sln[53] << "]" << std::endl;
		std::cout << "MBD: front-left tyre position pt3=\n\t[" << sln[55] << "\n\t " << sln[57] << "\n\t " << sln[56] << "]" << std::endl;
		std::cout << "MBD: front-right tyre position pt4=\n\t[" << sln[58] << "\n\t " << sln[60] << "\n\t " << sln[59] << "]" << std::endl;

	}

	/*
	Destructor
	*/
	~MBD_method() {
		mkl_free_buffers();
		mkl_free(r1);
		mkl_free(r2);
		mkl_free(r3);
		mkl_free(r4);
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
		mkl_free(vw1);
		mkl_free(vw2);
		mkl_free(vw3);
		mkl_free(vw4);
		mkl_free(vt1);
		mkl_free(vt2);
		mkl_free(vt3);
		mkl_free(vt4);
		mkl_free(pcc);
		mkl_free(FC);
		mkl_free(r1_tilda);
		mkl_free(r2_tilda);
		mkl_free(r3_tilda);
		mkl_free(r4_tilda);
		mkl_free(FW);
		mkl_free(FT);
		mkl_free(A_Ic);
		mkl_free(A_rem);
	}

};
