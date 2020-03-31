/*  Copyright &copy; 2019, Dr. Stefan Sicklinger, Munich \n
*
*  All rights reserved.
*
*  This file is part of EVAA.
*
*  EVAA is free software: you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation, either version 3 of the License, or
*  (at your option) any later version.
*
*  EVAA is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with EVAA.  If not, see http://www.gnu.org/licenses/.
*/
/***********************************************************************************************//**
* \file EVAAComputeEngine.h
* This file holds the class of ComputeEngine.
* \date 6/13/2019
**************************************************************************************************/
#pragma once
#include <string>
#include <mkl.h>
#include <chrono>
#include "MathLibrary.h"
#include "ReadXML.h"


/********//**
* \brief Class ComputeEngine the core of STACCATO
***********/
class EVAAComputeEngine {

public:
	/***********************************************************************************************
	* \brief Constructor
	* \param[in] file name of xml file to call singelton constructor of metadatabase
	* \author Stefan Sicklinger
	***********/
	EVAAComputeEngine(std::string _xmlFileName);
	/***********************************************************************************************
	* \brief Destructor
	* \author Stefan Sicklinger
	***********/
	~EVAAComputeEngine();
	/***********************************************************************************************
	* \brief prepare compute engine
	* \author Stefan Sicklinger
	***********/
	void prepare(void);
	/***********************************************************************************************
	* \brief compute engine
	* \author Stefan Sicklinger
	***********/
	// EIGEN
	void computeEigen11DOF(void);
	/***********************************************************************************************
	* \brief compute engine
	* \author Stefan Sicklinger
	***********/
	// BLAZE
	void computeBlaze11DOF(void);
	/***********************************************************************************************
	* \brief compute engine
	* \author Stefan Sicklinger
	***********/
	// MKL
	void computeMKL11DOF(void);
	/***********************************************************************************************
	* \brief compute engine
	* \author Stefan Sicklinger
	***********/
	// MKL Linear solver
	void computeMKLlinear11dof(void);
	void computeMKLlinear11dof_reduced(void);
	/***********************************************************************************************
	* \brief compute engine
	* \author Stefan Sicklinger
	***********/

	// MKL MBD solver
	void computeMBD(void);

	// MKL
/*	void computeMKLNasa(void);
	void computeMKLNasa_example(void);
*/
	/***********************************************************************************************
	* \brief clean compute engine free memory
	* \author Stefan Sicklinger
	***********/
	void clean(void);
private:
	std::string _xmlFileName;
	Simulation_Parameters _parameters;
};

template <class T>
class linear11dof
{
private:
	// define constants
	const int alignment = 64;
	const int dim = 3;
	const int num_wheels = 4;
	const int dim_x_dim = dim * dim;
	const int num_wheels_x_dim = num_wheels * dim;
	const int DOF_diag = 9; // the diagonal elements from A
	int DOF;

	T k_body_fl;
	T k_tyre_fl;
	T k_body_fr;
	T k_tyre_fr;
	T k_body_rl;
	T k_tyre_rl;
	T k_body_rr;
	T k_tyre_rr;
	T l_long_fl;
	T l_long_fr;
	T l_long_rl;
	T l_long_rr;
	T l_lat_fl;
	T l_lat_fr;
	T l_lat_rl;
	T l_lat_rr;
	T mass_Body;
	T I_body_xx;
	T I_body_zz;
	T mass_wheel_fl;
	T mass_tyre_fl;
	T mass_wheel_fr;
	T mass_tyre_fr;
	T mass_wheel_rl;
	T mass_tyre_rl;
	T mass_wheel_rr;
	T mass_tyre_rr;
	T u_init_body;
	T theta_x_init_body;
	T theta_z_init_body;

	//// Solver type selection based on type of boundary condition
	std::string condition_type;
	std::string cond1 = "fixed_to_road";
	std::string cond2 = "road_force";
	////////////////////////////////////////////////////////////////////////////////////////////
	T tend_;
	T u_init_;
	T du_init_;
	T h_;
	int i;
	T *quad_angle_init, *euler_angle_init;
	T* M, * temp, * K, * K_trans, * D, * M_red, * D_red, * K_red;
	T* u_sol, * u_sol_red, * u_n_p_1, * u_n_p_1_red, * u_n_m_1, * u_n, * u_n_red, * u_n_m_1_red, * A, * Ared, * B, * Bred, * f_n_p_1, * f_n_p_1_red;
	size_t* tyre_index_set;
	size_t num_tyre = 4;
	T* time; // this is not necessary

	void check_status(lapack_int status) {
		if (status == 1) {
			std::cout<<"Matrix non Positive Definite"<<std::endl;
			exit(5);
		}
		else if (status == -1) {
			std::cout << "Matrix contain illegal value" << std::endl;
			exit(5);
		}
	}
	void write_vector(T* vect, int count) {
		std::cout << "Debug mode print" << std::endl;
		for (size_t i = 0; i < count; ++i) {
			//std::cout << vect[i] << std::endl;
			std::cout.precision(15);
			std::cout << std::scientific << vect[i] << std::endl;
			//printf("%1.15f\n", vect[i]);
		}
	}

	void write_matrix(T* vect, int count) {
		std::cout << "Debug mode print" << std::endl;
		for (size_t i = 0; i < count; ++i) {
			//std::cout << vect[i] << std::endl;
			std::cout.precision(5);
			for (size_t j = 0; j < count; ++j) {
				std::cout << std::scientific << vect[i*count+j] << "  ";
			}
			std::cout <<"\n"<< std::endl;
			//printf("%1.15f\n", vect[i]);
		}
	}

	void apply_normal_force(T* force, T* u, size_t* index, size_t n) {
		#pragma loop( ivdep )
		for (int i = 0; i < n; ++i) {
			force[index[i]] = force[index[i]] > 0 ? force[index[i]] : 0;
		}
		#pragma loop( ivdep )
		for (int i = 0; i < n; ++i) {
			u[index[i]] = u[index[i]]>0 ? u[index[i]] : 0;
		}
	}
	void compute_normal_force(T* K, T* u, T* force, size_t* index, size_t dim,size_t n) {
		#pragma loop( ivdep )
		for (int i = 0; i < n; ++i) {
			force[index[i]] = -K[index[i]*dim + index[i]]*u[index[i]];
		}
	}



public:
	linear11dof(const Simulation_Parameters &params){

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////// Extract Data from parser /////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		DOF = params.DOF;
		k_body_fl = params.k_body[2];
		k_tyre_fl = params.k_tyre[2];
		k_body_fr = params.k_body[3];
		k_tyre_fr = params.k_tyre[3];
		k_body_rl = params.k_body[1];
		k_tyre_rl = params.k_tyre[1];
		k_body_rr = params.k_body[0];
		k_tyre_rr = params.k_tyre[0];
		l_long_fl = params.l_long[2];
		l_long_fr = params.l_long[3];
		l_long_rl = params.l_long[1];
		l_long_rr = params.l_long[0];
		l_lat_fl = params.l_lat[2];
		l_lat_fr = params.l_lat[3];
		l_lat_rl = params.l_lat[1];
		l_lat_rr = params.l_lat[0];
		mass_Body = params.mass_body;
		I_body_xx = params.I_body[0];
		I_body_zz = params.I_body[2];
		mass_wheel_fl = params.mass_wheel[2];
		mass_tyre_fl = params.mass_tyre[2];
		mass_wheel_fr = params.mass_wheel[3];
		mass_tyre_fr = params.mass_tyre[3];
		mass_wheel_rl = params.mass_wheel[1];
		mass_tyre_rl = params.mass_tyre[1];
		mass_wheel_rr = params.mass_wheel[0];
		mass_tyre_rr = params.mass_tyre[0];

		h_ = params.timestep;
		tend_ = params.num_time_iter*h_;
		quad_angle_init = (T*)mkl_calloc(4, sizeof(T), alignment);
		euler_angle_init = (T*)mkl_calloc(3, sizeof(T), alignment);
		quad_angle_init[0] = params.initial_angle[0];
		quad_angle_init[1] = params.initial_angle[1];
		quad_angle_init[2] = params.initial_angle[2];
		quad_angle_init[3] = params.initial_angle[3];
		MathLibrary::ToEulerAngles(quad_angle_init, euler_angle_init);
		u_init_body = params.initial_pos_body[1];
		theta_x_init_body = euler_angle_init[0];
		theta_z_init_body = euler_angle_init[2];

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////// System Mono ///////////////////////////////////////////////////////
		///////////////////////////////// Memory Allocation and matrix formulation ///////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		M = (T*)mkl_calloc(DOF * DOF, sizeof(T), alignment);
		temp = (T*)mkl_malloc(DOF * sizeof(T), alignment);
		K = (T*)mkl_calloc(DOF * DOF, sizeof(T), alignment);
		K_trans = (T*)mkl_calloc(DOF * DOF, sizeof(T), alignment);
		D = (T*)mkl_calloc(DOF * DOF, sizeof(T), alignment);
		u_n_m_1 = (T*)mkl_calloc(DOF, sizeof(T), alignment);
		u_n = (T*)mkl_calloc(DOF, sizeof(T), alignment);
		f_n_p_1 = (T*)mkl_calloc(DOF, sizeof(T), alignment);

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////// Initial Iteration vector ////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		u_n[0] = u_init_body;
		u_n[1] = theta_x_init_body;
		u_n[2] = theta_z_init_body;
		u_n[3] = params.initial_pos_wheel[2 * 3 + 1];
		u_n[4] = params.initial_pos_tyre[2 * 3 + 1];
		u_n[5] = params.initial_pos_wheel[3 * 3 + 1];
		u_n[6] = params.initial_pos_tyre[3 * 3 + 1];
		u_n[7] = params.initial_pos_wheel[1 * 3 + 1];
		u_n[8] = params.initial_pos_tyre[1 * 3 + 1];
		u_n[9] = params.initial_pos_wheel[0 * 3 + 1];
		u_n[10] = params.initial_pos_tyre[0 * 3 + 1];

		u_n_m_1[0] = params.initial_vel_body[1];
		u_n_m_1[1] = params.initial_ang_vel_body[0];
		u_n_m_1[2] = params.initial_ang_vel_body[2];
		u_n_m_1[3] = params.initial_vel_wheel[2 * 3 + 1];
		u_n_m_1[4] = params.initial_vel_tyre[2 * 3 + 1];
		u_n_m_1[5] = params.initial_vel_wheel[3 * 3 + 1];
		u_n_m_1[6] = params.initial_vel_tyre[3 * 3 + 1];
		u_n_m_1[7] = params.initial_vel_wheel[1 * 3 + 1];
		u_n_m_1[8] = params.initial_vel_tyre[1 * 3 + 1];
		u_n_m_1[9] = params.initial_vel_wheel[0 * 3 + 1];
		u_n_m_1[10] = params.initial_vel_tyre[0 * 3 + 1];

		f_n_p_1[0] = params.external_force_body[1];
		f_n_p_1[3] = params.external_force_wheel[2 * 3 + 1];
		f_n_p_1[4] = params.external_force_tyre[2 * 3 + 1];
		f_n_p_1[5] = params.external_force_wheel[3 * 3 + 1];
		f_n_p_1[6] = params.external_force_tyre[3 * 3 + 1];
		f_n_p_1[7] = params.external_force_wheel[1 * 3 + 1];
		f_n_p_1[8] = params.external_force_tyre[1 * 3 + 1];
		f_n_p_1[9] = params.external_force_wheel[0 * 3 + 1];
		f_n_p_1[10] = params.external_force_tyre[0 * 3 + 1];

		

		cblas_dscal(DOF, -h_, u_n_m_1, 1);
		cblas_daxpy(DOF, 1, u_n, 1, u_n_m_1, 1);

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////////// Mass //////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		temp[0] = mass_Body;
		temp[1] = I_body_xx;
		temp[2] = I_body_zz;
		temp[3] = mass_wheel_fl;
		temp[4] = mass_tyre_fl;
		temp[5] = mass_wheel_fr;
		temp[6] = mass_tyre_fr;
		temp[7] = mass_wheel_rl;
		temp[8] = mass_tyre_rl;
		temp[9] = mass_wheel_rr;
		temp[10] = mass_tyre_rr;

		MathLibrary::allocate_to_diagonal(M, temp, DOF);
		//M = diag([mass_Body, I_body_xx, I_body_yy, mass_wheel_fl, mass_tyre_fl, mass_wheel_fr, mass_tyre_fr, mass_wheel_rl, mass_tyre_rl, mass_wheel_rr, mass_tyre_rr]);
		i = 0;
		temp[0] = k_body_fl + k_body_fr + k_body_rl + k_body_rr; // K[i*DOF + 0] to be set later
		K[i * DOF + 1] = k_body_fl * l_lat_fl - k_body_fr * l_lat_fr + k_body_rl * l_lat_rl - k_body_rr * l_lat_rr;
		K[i * DOF + 2] = -k_body_fl * l_long_fl - k_body_fr * l_long_fr + k_body_rl * l_long_rl + k_body_rr * l_long_rr;
		K[i * DOF + 3] = -k_body_fl;
		K[i * DOF + 4] = 0;
		K[i * DOF + 5] = -k_body_fr;
		K[i * DOF + 6] = 0;
		K[i * DOF + 7] = -k_body_rl;
		K[i * DOF + 8] = 0;
		K[i * DOF + 9] = -k_body_rr;
		K[i * DOF + 10] = 0;

		i = 1;
		temp[1] = l_lat_fl * l_lat_fl * k_body_fl + l_lat_fr * l_lat_fr * k_body_fr + l_lat_rl * l_lat_rl * k_body_rl + l_lat_rr * l_lat_rr * k_body_rr; // K[i*DOF + 1] to be set later
		K[i * DOF + 2] = -l_long_fl * l_lat_fl * k_body_fl + l_lat_fr * l_long_fr * k_body_fr + l_long_rl * l_lat_rl * k_body_rl - l_long_rr * l_lat_rr * k_body_rr;
		K[i * DOF + 3] = -l_lat_fl * k_body_fl;
		K[i * DOF + 4] = 0;
		K[i * DOF + 5] = l_lat_fr * k_body_fr;
		K[i * DOF + 6] = 0;
		K[i * DOF + 7] = -l_lat_rl * k_body_rl;
		K[i * DOF + 8] = 0;
		K[i * DOF + 9] = l_lat_rr * k_body_rr;
		K[i * DOF + 10] = 0;

		i = 2;
		temp[2] = l_long_fl * l_long_fl * k_body_fl + l_long_fr * l_long_fr * k_body_fr + l_long_rl * l_long_rl * k_body_rl + l_long_rr * l_long_rr * k_body_rr; // K[i*DOF + 2] to be set later
		K[i * DOF + 3] = l_long_fl * k_body_fl;
		K[i * DOF + 4] = 0;
		K[i * DOF + 5] = l_long_fr * k_body_fr;
		K[i * DOF + 6] = 0;
		K[i * DOF + 7] = -l_long_rl * k_body_rl;
		K[i * DOF + 8] = 0;
		K[i * DOF + 9] = -l_long_rr * k_body_rr;
		K[i * DOF + 10] = 0;

		i = 3;

		temp[3] = k_body_fl + k_tyre_fl; // K[i*DOF + 3]
		K[i * DOF + 4] = -k_tyre_fl;
		K[i * DOF + 5] = 0;
		K[i * DOF + 6] = 0;
		K[i * DOF + 7] = 0;
		K[i * DOF + 8] = 0;
		K[i * DOF + 9] = 0;
		K[i * DOF + 10] = 0;
		// all others are zero

		i = 4;
		temp[4] = k_tyre_fl; //K[i*DOF + 4]
		K[i * DOF + 5] = 0;
		K[i * DOF + 6] = 0;
		K[i * DOF + 7] = 0;
		K[i * DOF + 8] = 0;
		K[i * DOF + 9] = 0;
		K[i * DOF + 10] = 0;

		i = 5;
		temp[5] = k_body_fr + k_tyre_fr; // K[i*DOF + 5]
		K[i * DOF + 6] = -k_tyre_fr;
		K[i * DOF + 7] = 0;
		K[i * DOF + 8] = 0;
		K[i * DOF + 9] = 0;
		K[i * DOF + 10] = 0;

		i = 6;
		temp[6] = k_tyre_fr; // K[i*DOF + 6]
		K[i * DOF + 7] = 0;
		K[i * DOF + 8] = 0;
		K[i * DOF + 9] = 0;
		K[i * DOF + 10] = 0;

		i = 7;
		temp[7] = k_body_rl + k_tyre_rl; // K[i*DOF + 7]
		K[i * DOF + 8] = -k_tyre_rl;
		K[i * DOF + 9] = 0;
		K[i * DOF + 10] = 0;


		i = 8;
		temp[8] = k_tyre_rl; // K[i*DOF + 8]
		K[i * DOF + 9] = 0;
		K[i * DOF + 10] = 0;

		i = 9;
		temp[9] = k_body_rr + k_tyre_rr; // K[i*DOF + 9]
		K[i * DOF + 10] = -k_tyre_rr;

		i = 10;
		temp[10] = k_tyre_rr; // K[i*DOF + 10]

		// K=K+K'-diag(diag(K));
		cblas_dcopy(DOF * DOF, K, 1, K_trans, 1);
		mkl_dimatcopy('R', 'T', DOF, DOF, 1.0, K_trans, DOF, DOF); // get transpose of matrix
		cblas_daxpy(DOF * DOF, 1.0, K_trans, 1, K, 1); // K = K + K'
		MathLibrary::allocate_to_diagonal(K, temp, DOF); // K = K + K'+ diag(K)

		// D = K *0;
		// default D value is 0
	}

	void apply_boundary_condition(std::string s) {
		
		condition_type = s;

		if (s == cond1) {
			T* start_loc_k_trans, * start_loc_k, * start_loc_D_red, * start_loc_D;
			// using K_trans as temperory variable copy all element of K to k_trans
			// cblas_dcopy(DOF*DOF, K, 1, K_trans, 1);
			// change the DOF for further reference
			DOF = DOF - 4; // DOF 7 now
			K_red = (T*)mkl_calloc(DOF * DOF, sizeof(T), alignment);
			M_red = (T*)mkl_calloc(DOF * DOF, sizeof(T), alignment);
			D_red = (T*)mkl_calloc(DOF * DOF, sizeof(T), alignment);
			u_n_red = (T*)mkl_calloc((DOF), sizeof(T), alignment);
			u_n_m_1_red = (T*)mkl_calloc((DOF), sizeof(T), alignment);
			// memory allocation for the reduced force field
			f_n_p_1_red = (T*)mkl_calloc(DOF, sizeof(T), alignment);
			

			// reassign the values of K in the new K for first 4 DOF 
			for (int i = 0; i < DOF - 3; ++i) {
				start_loc_k_trans = K + i * (DOF + 4);
				start_loc_k = K_red + i * DOF;
				start_loc_D_red = D_red + i * DOF;
				start_loc_D = D + i * (DOF + 4);

				cblas_dcopy(DOF - 3, start_loc_k_trans, 1, start_loc_k, 1);

				cblas_dcopy(DOF - 3, start_loc_D, 1, start_loc_D_red, 1);

				start_loc_k[DOF - 3] = start_loc_k_trans[DOF - 3 + 1];
				start_loc_k[DOF - 3 + 1] = start_loc_k_trans[DOF - 3 + 3];
				start_loc_k[DOF - 3 + 2] = start_loc_k_trans[DOF - 3 + 5];

				// needs rechecking
				start_loc_D_red[DOF - 3] = start_loc_D[DOF - 3 + 1];
				start_loc_D_red[DOF - 3 + 1] = start_loc_D[DOF - 3 + 3];
				start_loc_D_red[DOF - 3 + 2] = start_loc_D[DOF - 3 + 5];

				M_red[i * DOF + i] = M[i * (DOF + 4) + i];
				f_n_p_1_red[i] = f_n_p_1[i];
				u_n_red[i] = u_n[i];
				u_n_m_1_red[i] = u_n_m_1[i];
			}
			int j = 4;
			for (int i = 5; i < (DOF + 4); i = i + 2) {
				start_loc_k_trans = K + (i) * (DOF + 4);
				start_loc_k = K_red + j * DOF;
				start_loc_D_red = D_red + j * DOF;
				start_loc_D = D + i * (DOF + 4);

				cblas_dcopy(DOF - 3, start_loc_k_trans, 1, start_loc_k, 1);
				cblas_dcopy(DOF - 3, start_loc_D, 1, start_loc_D_red, 1);

				start_loc_k[DOF - 3] = start_loc_k_trans[DOF - 3 + 1];
				start_loc_k[DOF - 3 + 1] = start_loc_k_trans[DOF - 3 + 3];
				start_loc_k[DOF - 3 + 2] = start_loc_k_trans[DOF - 3 + 5];

				// needs rechecking
				start_loc_D_red[DOF - 3] = start_loc_D[DOF - 3 + 1];
				start_loc_D_red[DOF - 3 + 1] = start_loc_D[DOF - 3 + 3];
				start_loc_D_red[DOF - 3 + 2] = start_loc_D[DOF - 3 + 5];

				M_red[j * DOF + j] = M[i * (DOF + 4) + i];
				f_n_p_1_red[j] = f_n_p_1[i];
				u_n_red[j] = u_n[i];
				u_n_m_1_red[j] = u_n_m_1[i];
				j++;
			}
		}
		else if (s == cond2) {
			// memory allocation for the force field
			tyre_index_set = (size_t*)mkl_calloc(num_tyre, sizeof(size_t), alignment);
			for (int i = 4, j=0; i < DOF && j<num_tyre; i = i + 2, j++) {
				tyre_index_set[j] = i;
			}
		}
		else {
			throw "Incorrect boundary condition";
		}
	}
	void solve(T* sol_vect) {
		if (condition_type == cond1) {
			solve_reduced(sol_vect);
		}
		else if (condition_type == cond2) {
			solve_full(sol_vect);
		}
		else {
			throw "Inappropriate boundary condition!";
		}
	}

	void solve_full(T* sol_vect) {
		///////////////////////////////////// For 7 DOF System //////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////// Memory Allocation for Intermediate Solution Vectors //////////////////
		////////////////////////////////////////////////////////////////////////////////////////////////////////////
		int sol_size = (floor(tend_ / h_) + 1);
		T factor_h2 = (1 / (h_ * h_));
		T factor_h = (1 / (h_));
		int mat_len = (DOF) * (DOF);
		u_sol = (T*)mkl_calloc(sol_size * (DOF), sizeof(T), alignment);
		u_n_p_1 = (T*)mkl_calloc((DOF), sizeof(T), alignment);
		A = (T*)mkl_calloc(mat_len, sizeof(T), alignment);
		B = (T*)mkl_calloc(mat_len, sizeof(T), alignment);
		time = (T*)mkl_calloc(sol_size, sizeof(T), alignment);
	

		// A=((1/(h*h))*M+(1/h)*D+K);
		cblas_daxpy(mat_len, factor_h2, M, 1, A, 1);
		cblas_daxpy(mat_len, factor_h, D, 1, A, 1);
		cblas_daxpy(mat_len, 1, K, 1, A, 1);
		// Cholesky factorization of A
		lapack_int status;
		//lapack_int* piv = (lapack_int*)mkl_calloc((DOF + 4), sizeof(lapack_int), alignment);
		//write_matrix(A, DOF);
		status = LAPACKE_dpotrf(LAPACK_ROW_MAJOR, 'L', DOF, A, DOF);
		//status = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, DOF + 4, DOF + 4, A, DOF + 4, piv);
		check_status(status);
		
		// B=((2/(h*h))*M+(1/h)*D);
		cblas_daxpy(mat_len, 2 * factor_h2, M, 1, B, 1);
		cblas_daxpy(mat_len, factor_h, D, 1, B, 1);
		int iter = 1;
		T t = h_;
		double eps = h_/100;
		/*auto start = std::chrono::steady_clock::now();*/
		while (std::abs(t-(tend_+h_)) > eps) {
			// u_n_p_1=A\(B*u_n-((1/(h*h))*M)*u_n_m_1+f_n_p_1);
			// u_n_p_1 = B*u_n

			cblas_dgemv(CblasRowMajor, CblasNoTrans, DOF, DOF, 1, B, DOF, u_n, 1, 0, u_n_p_1, 1);
			// u_n_p_1 = -((1/(h*h))*M)*u_n_m_1 + u_n_p_1
			cblas_dgemv(CblasRowMajor, CblasNoTrans, DOF, DOF, -factor_h2, M, DOF, u_n_m_1, 1, 1, u_n_p_1, 1);
			// u_n_p_1 = f_n_p_1 + u_n_p_1
			cblas_daxpy(DOF, 1, f_n_p_1, 1, u_n_p_1, 1);
			// u_n_p_1=A\u_n_p_1
			LAPACKE_dpotrs(LAPACK_ROW_MAJOR, 'L', DOF, 1, A, DOF, u_n_p_1, 1);
			//LAPACKE_dgetrs(LAPACK_ROW_MAJOR, 'N', DOF + 4, 1, A, DOF + 4, piv, u_n_p_1, 1);
			///////////////// Debug print statements ///////////////////////////////
			/*if (iter==2){
				write_vector(u_n_p_1, DOF + 4);
			}*/
			time[iter] = t;
			/*
			f_n_p_1(idx) = -K(idx,idx)*u_n_p_1(idx);
			f_n_p_1(idx) = f_update(f_n_p_1, idx);
			u_n_p_1(idx) = 0;
			*/
			////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//////////////////////////////// Normal force computation here /////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////////////////////////////
			// compute_normal_force(K, u_n_p_1, f_n_p_1, tyre_index_set, DOF, num_tyre);
			
			
			///////////////// Debug print statements ///////////////////////////////
			/*if (iter==10){
				write_vector(u_n_p_1, DOF);
			}*/
			
			
			
			// apply_normal_force(f_n_p_1, u_n_p_1, tyre_index_set, num_tyre);
			
			////////////////////////////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////////////////////////////

			// u_sol(j,:)=u_n_p_1;
			cblas_dcopy(DOF, u_n_p_1, 1, u_sol + iter * (DOF), 1);
			// u_n_m_1=u_n;
			cblas_dcopy(DOF, u_n, 1, u_n_m_1, 1);
			// u_n    =u_n_p_1;
			cblas_dcopy(DOF, u_n_p_1, 1, u_n, 1);
			iter++;
			t += h_;
		}
		/*auto end = std::chrono::steady_clock::now();
		std::cout << "Elapsed time in nanoseconds : "
			<< std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
			<< " ms" << std::endl;*/
		std::cout << "iter = " << iter << " sol_size = "<< sol_size <<"\n\n" << std::endl;
		cblas_dcopy(DOF, u_sol + (iter-1)*(DOF), 1, sol_vect, 1);
		clean("full");

	}
	void solve_reduced(T* sol_vect) {
		///////////////////////////////////// For 7 DOF System //////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////// Memory Allocation for Intermediate Solution Vectors //////////////////
		////////////////////////////////////////////////////////////////////////////////////////////////////////////
		int sol_size = (floor(tend_ / h_) + 1);
		T factor_h2 = (1 / (h_ * h_));
		T factor_h = (1 / (h_));
		int mat_len = (DOF) * (DOF);
		//u_sol = (T*)mkl_calloc(sol_size * (DOF + 4), sizeof(T), alignment);
		u_sol_red = (T*)mkl_calloc(sol_size * (DOF), sizeof(T), alignment);
		//u_n_p_1 = (T*)mkl_calloc((DOF + 4), sizeof(T), alignment);
		u_n_p_1_red = (T*)mkl_calloc((DOF), sizeof(T), alignment);
		//u_n_m_1 = (T*)mkl_calloc((DOF + 4), sizeof(T), alignment);
		//u_n = (T*)mkl_calloc((DOF + 4), sizeof(T), alignment);
		
		//A = (T*)mkl_calloc(mat_len, sizeof(T), alignment);
		Ared = (T*)mkl_calloc((DOF) * (DOF), sizeof(T), alignment);
		//B = (T*)mkl_calloc(mat_len, sizeof(T), alignment);
		Bred = (T*)mkl_calloc((DOF) * (DOF), sizeof(T), alignment);
		time = (T*)mkl_calloc(sol_size, sizeof(T), alignment);
		

		// A=((1/(h*h))*M+(1/h)*D+K);
		cblas_daxpy(mat_len, factor_h2, M_red, 1, Ared, 1);
		cblas_daxpy(mat_len, factor_h, D_red, 1, Ared, 1);
		cblas_daxpy(mat_len, 1, K_red, 1, Ared, 1);
		// Cholesky factorization of A
		lapack_int status;
		//lapack_int* piv = (lapack_int*)mkl_calloc((DOF + 4), sizeof(lapack_int), alignment);
		status = LAPACKE_dpotrf(LAPACK_ROW_MAJOR, 'L', DOF, Ared, DOF);
		check_status(status);
		// B=((2/(h*h))*M+(1/h)*D);
		cblas_daxpy(mat_len, 2 * factor_h2, M_red, 1, Bred, 1);
		cblas_daxpy(mat_len, factor_h, D_red, 1, Bred, 1);
		int iter = 1;
		T t = h_;
		/*auto start = std::chrono::steady_clock::now();*/
		while (t < tend_) {
			// u_n_p_1=A\(B*u_n-((1/(h*h))*M)*u_n_m_1+f_n_p_1);
			// u_n_p_1 = B*u_n
			cblas_dgemv(CblasRowMajor, CblasNoTrans, DOF, DOF, 1, Bred, DOF, u_n_red, 1, 1, u_n_p_1_red, 1);
			// u_n_p_1 = -((1/(h*h))*M)*u_n_m_1 + u_n_p_1
			cblas_dgemv(CblasRowMajor, CblasNoTrans, DOF, DOF, -factor_h2, M_red, DOF, u_n_m_1_red, 1, 1, u_n_p_1_red, 1);
			// u_n_p_1 = f_n_p_1 + u_n_p_1
			cblas_daxpy(DOF, 1, f_n_p_1_red, 1, u_n_p_1_red, 1);
			// u_n_p_1=A\u_n_p_1
			LAPACKE_dpotrs(LAPACK_ROW_MAJOR, 'L', DOF, 1, Ared, DOF, u_n_p_1_red, 1);
			//LAPACKE_dgetrs(LAPACK_ROW_MAJOR, 'N', DOF + 4, 1, A, DOF + 4, piv, u_n_p_1, 1);
			///////////////// Debug print statements ///////////////////////////////
			/*if (iter==2){
				write_vector(u_n_p_1_red, DOF);
			}*/
			time[iter] = t;
			// u_sol(j,:)=u_n_p_1;
			cblas_dcopy(DOF, u_n_p_1_red, 1, u_sol_red + iter * (DOF), 1);
			// u_n_m_1=u_n;
			cblas_dcopy(DOF, u_n_red, 1, u_n_m_1_red, 1);
			// u_n    =u_n_p_1;
			cblas_dcopy(DOF, u_n_p_1_red, 1, u_n_red, 1);
			iter++;
			t += h_;
		}
		/*auto end = std::chrono::steady_clock::now();
		std::cout << "Elapsed time in nanoseconds : "
			<< std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
			<< " ms" << std::endl;*/
			//std::cout << "iter = " << iter << " sol_size = "<< sol_size <<"\n\n" << std::endl;
		cblas_dcopy(DOF, u_sol_red + (iter - 1)*(DOF), 1, sol_vect, 1);
		clean("reduced");

	}

	void clean(std::string s){
		std::string full = "full";
		std::string reduced = "reduced";
		if (s == full) {
			mkl_free(u_sol);
			mkl_free(u_n_p_1);
			mkl_free(tyre_index_set);
			mkl_free(A);
			mkl_free(B);
			mkl_free(time);
		}
		else if (s == reduced) {
			mkl_free(u_sol_red);
			mkl_free(u_n_red);
			mkl_free(u_n_m_1_red);
			mkl_free(u_n_p_1_red);
			mkl_free(f_n_p_1_red);
			mkl_free(Ared);
			mkl_free(Bred);
			mkl_free(time);
			mkl_free(K_red);
			mkl_free(M_red);
			mkl_free(D_red);
		}
	}

	~linear11dof() {
		mkl_free(M);
		mkl_free(temp);
		mkl_free(K);
		mkl_free(K_trans);
		mkl_free(D);
		mkl_free(u_n_m_1);
		mkl_free(u_n);
		mkl_free(quad_angle_init);
		mkl_free(euler_angle_init);
		mkl_free(f_n_p_1);
	}
};

template <class T>
class MBD_method
{
private:
	////////////////////////////// Simulation Parameters ///////////////////////////////////////////////////////////////
	int DIM = 3;
	int NUM_LEGS = 4;
	int alignment = 64;
	T h = 1e-3;
	size_t num_iter = 1e3;
	int max_iter = 1e3;
	T tol = 1e-7;
	size_t solution_dim = 61; /// this is by the formulation
	////////////////////////////// Car Definition ///////////////////////////////////////////////////////////////////////
	T k_body_fl = 28e3*0.69;
	T k_tyre_fl = 260e3;
	T k_body_fr = 28e3*0.69;
	T k_tyre_fr = 260e3;
	T k_body_rl = 16e3*0.82;
	T k_tyre_rl = 260e3;
	T k_body_rr = 16e3*0.82;
	T k_tyre_rr = 260e3;
	T k_body_rot_fl = 1e4;
	T k_body_rot_fr = 1e4;
	T k_body_rot_rl = 1e4;
	T k_body_rot_rr = 1e4;
	T k_tyre_rot_fl = 1e4;
	T k_tyre_rot_fr = 1e4;
	T k_tyre_rot_rl = 1e4;
	T k_tyre_rot_rr = 1e4;
	T c_body_fl = 0.0;    
	T c_tyre_fl = 0.0;
	T c_body_fr = 0.0;
	T c_tyre_fr = 0.0;
	T c_body_rl = 0.0;
	T c_tyre_rl = 0.0;
	T c_body_rr = 0.0;
	T c_tyre_rr = 0.0;
	T l_long_fl = 1.395;
	T l_long_fr = 1.395;
	T l_long_rl = 1.596;
	T l_long_rr = 1.596;
	T l_lat_fl = 2.0 * 0.84;
	T l_lat_fr = 2.0 * 0.84;
	T l_lat_rl = 2.0 * 0.84;
	T l_lat_rr = 2.0 * 0.84;
	T mass = 1936.0;
	T I_body_xx = 6400.0;
	T I_body_yy = 4800.0;
	T I_body_zz = 5000.0;
	T mass_wheel_fl = 145.0 / 2.0;
	T mass_tyre_fl = 30.0;
	T mass_wheel_fr = 145.0 / 2.0;
	T mass_tyre_fr = 30.0;
	T mass_wheel_rl = 135.0 / 2.0;
	T mass_tyre_rl = 30.0;
	T mass_wheel_rr = 135.0 / 2.0;
	T mass_tyre_rr = 30.0;
	T upper_spring_length_rr = 0.2;
	T upper_spring_length_rl = 0.2;
	T upper_spring_length_fl = 0.2;
	T upper_spring_length_fr = 0.2;
	T lower_spring_length_rr = 0.2;
	T lower_spring_length_rl = 0.2;
	T lower_spring_length_fl = 0.2;
	T lower_spring_length_fr = 0.2;
	T *Ic;
	T *mass_wheel, *mass_tyre;
	T *upper_spring_length, *lower_spring_length;
	T *upper_spring_stiffness, *lower_spring_stiffness, *upper_spring_damping, *lower_spring_damping;
	T *upper_rotational_stiffness, *lower_rotational_stiffness;
	T *vc, *vw1, *vw2, *vw3, *vw4, *vt1, *vt2, *vt3, *vt4; // velocity of the center of mass, wheel and tyre.
	//////////////////////////	External Forces terms //////////////////////////////////////////////////////
	T g = 0;
	T *FC;
	T FT1, FT2, FT3, FT4;
	T FW1, FW2, FW3, FW4;
	T FR1, FR2, FR3, FR4;
	//////////////////////////// Initial condition params //////////////////////////////////////////////////
	T *initial_upper_spring_length, *initial_lower_spring_length, *initial_orientation, *initial_angular_velocity;
	//////////////////////////// Arrays for intermediate steps /////////////////////////////////////////////
	T *r1, *r2, *r3, *r4;
	T *pcc; // position of center of mass
	T* x_vector;
	////////////////////////// Auxillary Parameters for Compute_f function //////////////////////////////////
	T *r1_tilda, *r2_tilda, *r3_tilda, *r4_tilda, **FW, *FT, *A_Ic, *A_rem;

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
		T* C_Nc = (T*)mkl_calloc((this->DIM)*(this->DIM), sizeof(T), this->alignment);

		//	1. qc = qc/norm(qc); This is in quaternions 
		T nrm = cblas_dnrm2(this->NUM_LEGS, initial_orientation_, 1);
		cblas_dscal(this->NUM_LEGS, 1.0 / nrm, initial_orientation_, 1);

		// 2.	C_Nc = get_basis(qc);
		get_basis(initial_orientation_, C_Nc);

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

	void get_basis(const T* initial_orientation, T* transfrormed_basis) {
		// quaternion initial_orientation yields basis(calculates basically the matrix rotation
		// the euclidian basis to the local basis)

		// get local basis vectors using unit quaternion rotation
		// s = 1 / norm(q) ^ 2;      %normalizer, only to remove numerical stuffy stuff
		T nrm = cblas_dnrm2(this->NUM_LEGS, initial_orientation, 1);
		T s = 1.0 / (nrm * nrm);
		cblas_dscal((this->DIM) * (this->DIM), 0.0, transfrormed_basis, 1);
		size_t i, j;

		T quad_sum = 0.0;
		for (i = 0; i < this->DIM; ++i) {
			quad_sum += initial_orientation[i] * initial_orientation[i];
		}
		i = 0;
		j = 0;
		transfrormed_basis[i * this->DIM + j] = 1 - 2 * s * (quad_sum - initial_orientation[i] * initial_orientation[i]);
		j = 1;
		transfrormed_basis[i * this->DIM + j] = 2 * s * (initial_orientation[i] * initial_orientation[j] - initial_orientation[i + 2] * initial_orientation[j + 2]);
		j = 2;
		transfrormed_basis[i * this->DIM + j] = 2 * s * (initial_orientation[i] * initial_orientation[j] + initial_orientation[i + 1] * initial_orientation[j + 1]);
		i = 1;
		j = 0;
		transfrormed_basis[i * this->DIM + j] = 2 * s * (initial_orientation[i] * initial_orientation[j] + initial_orientation[i + 2] * initial_orientation[j + 2]);
		j = 1;
		transfrormed_basis[i * this->DIM + j] = 1 - 2 * s * (quad_sum - initial_orientation[i] * initial_orientation[i]);
		j = 2;
		transfrormed_basis[i * this->DIM + j] = 2 * s * (initial_orientation[i] * initial_orientation[j] - initial_orientation[i - 1] * initial_orientation[j + 1]);
		i = 2;
		j = 0;
		transfrormed_basis[i * this->DIM + j] = 2 * s * (initial_orientation[i] * initial_orientation[j] - initial_orientation[i + 1] * initial_orientation[j + 1]);
		j = 1;
		transfrormed_basis[i * this->DIM + j] = 2 * s * (initial_orientation[i] * initial_orientation[j] + initial_orientation[i + 1] * initial_orientation[j - 1]);
		j = 2;
		transfrormed_basis[i * this->DIM + j] = 1 - 2 * s * (quad_sum - initial_orientation[i] * initial_orientation[i]);

	}

	void get_tilda(const T* input_vector, T* tilda_output) {
		/*
		This function is only suitable for a 3 dimensional system and renders unusable, might throw exceptions when used with other dimensions.
		given y: x_tilda*y = cross(x,y)
		[Stoneking, page 3 bottom]
		x_tilda = [	0, -x(3), x(2)
					x(3), 0, -x(1)
					-x(2), x(1), 0	]

		*/
		tilda_output[0] = 0;
		tilda_output[1] = -input_vector[2];
		tilda_output[2] = input_vector[1];
		tilda_output[3] = input_vector[2];
		tilda_output[4] = 0;
		tilda_output[5] = -input_vector[0];
		tilda_output[6] = -input_vector[1];
		tilda_output[7] = input_vector[0];
		tilda_output[8] = 0;
	}
	void write_matrix(T* vect, int count) {
		std::cout << "Debug mode print" << std::endl;
		for (size_t i = 0; i < count; ++i) {
			//std::cout << vect[i] << std::endl;
			std::cout.precision(5);
			for (size_t j = 0; j < count; ++j) {
				std::cout << std::scientific << vect[i*count + j] << "  ";
			}
			std::cout << "\n" << std::endl;
			//printf("%1.15f\n", vect[i]);
		}
		//exit(5);
	}
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
	MBD_method() {
		
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
		A_rem = (T*)mkl_calloc(9*this->DIM, sizeof(T), this->alignment);
		
		i = 0;
		r1[i] = -l_long_rr; r2[i] = -l_long_rl; r3[i] = l_long_fl; r4[i] = l_long_fr;
		Ic[i*DIM + i] = I_body_xx;
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
		initial_upper_spring_length[i] = 0.24;
		initial_lower_spring_length[i] = 0.24;
		initial_orientation[i] = 1;
		vc[i] = 0;
		vw1[i] = 0;
		vw2[i] = 0;
		vw3[i] = 0;
		vw4[i] = 0;
		vt1[i] = 0;
		vt2[i] = 0;
		vt3[i] = 0;
		vt4[i] = 0;
		pcc[i] = 0; // initial position of the center of mass is at 0
		FC[i] = 0;

		i = 1;
		r1[i] = 0; r2[i] = 0; r3[i] = 0; r4[i] = 0;
		Ic[i*DIM + i] = I_body_zz;
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
		initial_upper_spring_length[i] = 0.18;
		initial_lower_spring_length[i] = 0.20;
		initial_orientation[i] = 1;
		vc[i] = 0;
		vw1[i] = 0;
		vw2[i] = 0;
		vw3[i] = 0;
		vw4[i] = 0;
		vt1[i] = 0;
		vt2[i] = 0;
		vt3[i] = 0;
		vt4[i] = 0;
		pcc[i] = 0; // initial position of the center of mass is at 0
		FC[i] = -mass * g;

		i = 2;
		r1[i] = l_lat_rr; r2[i] = -l_lat_rl; r3[i] = -l_lat_fl; r4[i] = l_lat_fr;
		Ic[i*DIM + i] = I_body_yy;
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
		initial_upper_spring_length[i] = 0.21;
		initial_lower_spring_length[i] = 0.16;
		initial_orientation[i] = 1;
		vc[i] = 0;
		vw1[i] = 0;
		vw2[i] = 0;
		vw3[i] = 0;
		vw4[i] = 0;
		vt1[i] = 0;
		vt2[i] = 0;
		vt3[i] = 0;
		vt4[i] = 0;
		pcc[i] = 0; // initial position of the center of mass is at 0
		FC[i] = 0;

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
		initial_upper_spring_length[i] = 0.14;
		initial_lower_spring_length[i] = 0.18;
		initial_orientation[i] = 1;

		// A_Ic has cholesky factorization of Ic
		cblas_dcopy(this->DIM * this->DIM, Ic, 1, A_Ic, 1);
		LAPACKE_dpotrf(LAPACK_ROW_MAJOR, 'L', this->DIM, A_Ic, this->DIM);
		
	}
	
	void apply_boundary_condition() {
		// from line 119 in matlab
	}

	void get_road_force(T t, T y, T v, T m, T F, T* Fr) {
		*Fr = 0.0;
	}
	

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
		T* C_cN = (T*)mkl_calloc((this->DIM) * (this->DIM), sizeof(T), this->alignment);
		T* r_up1 = (T*)mkl_calloc(this->DIM, sizeof(T), this->alignment);
		T* r_up2 = (T*)mkl_calloc(this->DIM, sizeof(T), this->alignment);
		T* r_up3 = (T*)mkl_calloc(this->DIM, sizeof(T), this->alignment);
		T* r_up4 = (T*)mkl_calloc(this->DIM, sizeof(T), this->alignment);
		T* r_low1 = (T*)mkl_calloc(this->DIM, sizeof(T), this->alignment);
		T* r_low2 = (T*)mkl_calloc(this->DIM, sizeof(T), this->alignment);
		T* r_low3 = (T*)mkl_calloc(this->DIM, sizeof(T), this->alignment);
		T* r_low4 = (T*)mkl_calloc(this->DIM, sizeof(T), this->alignment);
		get_basis(qc_, C_cN);
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
		cblas_dcopy(this->DIM, pcc_, 1, r_up1, 1);
		cblas_dgemv(CblasRowMajor, CblasNoTrans, this->DIM, this->DIM, 1, C_cN, this->DIM, this->r1, 1, 1, r_up1, 1);
		cblas_daxpy(this->DIM, -1.0, pw1_, 1, r_up1, 1);
		
		// r_up2
		cblas_dcopy(this->DIM, pcc_, 1, r_up2, 1);
		cblas_dgemv(CblasRowMajor, CblasNoTrans, this->DIM, this->DIM, 1, C_cN, this->DIM, this->r2, 1, 1, r_up2, 1);
		cblas_daxpy(this->DIM, -1.0, pw2_, 1, r_up2, 1);
		// r_up3
		cblas_dcopy(this->DIM, pcc_, 1, r_up3, 1);
		cblas_dgemv(CblasRowMajor, CblasNoTrans, this->DIM, this->DIM, 1, C_cN, this->DIM, this->r3, 1, 1, r_up3, 1);
		cblas_daxpy(this->DIM, -1.0, pw3_, 1, r_up3, 1);
		// r_up4
		cblas_dcopy(this->DIM, pcc_, 1, r_up4, 1);
		cblas_dgemv(CblasRowMajor, CblasNoTrans, this->DIM, this->DIM, 1, C_cN, this->DIM, this->r4, 1, 1, r_up4, 1);
		cblas_daxpy(this->DIM, -1.0, pw4_, 1, r_up4, 1);
		// r_low1
		cblas_dcopy(this->DIM, pw1_, 1, r_low1, 1);
		cblas_daxpy(this->DIM, -1.0, pt1_, 1, r_low1, 1);
		
		// r_low2
		cblas_dcopy(this->DIM, pw2_, 1, r_low2, 1);
		cblas_daxpy(this->DIM, -1.0, pt2_, 1, r_low2, 1);
		
		// r_low3
		cblas_dcopy(this->DIM, pw3_, 1, r_low3, 1);
		cblas_daxpy(this->DIM, -1.0, pt3_, 1, r_low3, 1);
		// r_low4
		cblas_dcopy(this->DIM, pw4_, 1, r_low4, 1);
		cblas_daxpy(this->DIM, -1.0, pt4_, 1, r_low4, 1);

		inv_norm_r_up1 = 1.0 / cblas_dnrm2(this->DIM, r_up1, 1);
		inv_norm_r_up2 = 1.0 / cblas_dnrm2(this->DIM, r_up2, 1);
		inv_norm_r_up3 = 1.0 / cblas_dnrm2(this->DIM, r_up3, 1);
		inv_norm_r_up4 = 1.0 / cblas_dnrm2(this->DIM, r_up4, 1);

		inv_norm_r_low1 = 1.0 / cblas_dnrm2(this->DIM, r_low1, 1);
		inv_norm_r_low2 = 1.0 / cblas_dnrm2(this->DIM, r_low2, 1);
		inv_norm_r_low3 = 1.0 / cblas_dnrm2(this->DIM, r_low3, 1);
		inv_norm_r_low4 = 1.0 / cblas_dnrm2(this->DIM, r_low4, 1);

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
		T* upper_normal1 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		T* upper_normal2 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		T* upper_normal3 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		T* upper_normal4 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		T* lower_normal1 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		T* lower_normal2 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		T* lower_normal3 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		T* lower_normal4 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		T* col_dat = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		T *upper_angle1, *upper_angle2, *upper_angle3, *upper_angle4, *lower_angle1, *lower_angle2, *lower_angle3, *lower_angle4;
		upper_angle1 = new T;
		upper_angle2 = new T;
		upper_angle3 = new T;
		upper_angle4 = new T;
		lower_angle1 = new T;
		lower_angle2 = new T;
		lower_angle3 = new T;
		lower_angle4 = new T;
		cblas_dcopy(this->DIM, C_cN + 1, this->DIM, col_dat, 1);
		
		MathLibrary::get_quaternion<T>(r_up1, col_dat, upper_angle1, upper_normal1, this->DIM);
		MathLibrary::get_quaternion<T>(r_up2, col_dat, upper_angle2, upper_normal2, this->DIM);
		MathLibrary::get_quaternion<T>(r_up3, col_dat, upper_angle3, upper_normal3, this->DIM);
		MathLibrary::get_quaternion<T>(r_up4, col_dat, upper_angle4, upper_normal4, this->DIM);
		
		MathLibrary::get_quaternion(r_low1, r_up1, lower_angle1, lower_normal1, this->DIM);
		MathLibrary::get_quaternion(r_low2, r_up2, lower_angle2, lower_normal2, this->DIM);
		MathLibrary::get_quaternion(r_low3, r_up3, lower_angle3, lower_normal3, this->DIM);
		MathLibrary::get_quaternion(r_low4, r_up4, lower_angle4, lower_normal4, this->DIM);
		
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
		T* upper_force1 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		T* upper_force2 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		T* upper_force3 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		T* upper_force4 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		T* lower_force1 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		T* lower_force2 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		T* lower_force3 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		T* lower_force4 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		T scale;
		
		cblas_dcopy(this->DIM, r_up1, 1, upper_force1, 1); 
		scale = this->upper_spring_stiffness[0] * (1.0 - this->upper_spring_length[0] * inv_norm_r_up1);
		cblas_dscal(this->DIM, scale, upper_force1, 1);
		cblas_dcopy(this->DIM, r_up2, 1, upper_force2, 1);
		scale = this->upper_spring_stiffness[1] * (1.0 - this->upper_spring_length[1] * inv_norm_r_up2);
		cblas_dscal(this->DIM, scale, upper_force2, 1);
		cblas_dcopy(this->DIM, r_up3, 1, upper_force3, 1);
		scale = this->upper_spring_stiffness[2] * (1.0 - this->upper_spring_length[2] * inv_norm_r_up3);
		cblas_dscal(this->DIM, scale, upper_force3, 1);
		cblas_dcopy(this->DIM, r_up4, 1, upper_force4, 1);
		scale = this->upper_spring_stiffness[3] * (1.0 - this->upper_spring_length[3] * inv_norm_r_up4);
		cblas_dscal(this->DIM, scale, upper_force4, 1);

		cblas_dcopy(this->DIM, r_low1, 1, lower_force1, 1);
		scale = this->lower_spring_stiffness[0] * (1.0 - this->lower_spring_length[0] * inv_norm_r_low1);
		cblas_dscal(this->DIM, scale, lower_force1, 1);
		cblas_dcopy(this->DIM, r_low2, 1, lower_force2, 1);
		scale = this->lower_spring_stiffness[1] * (1.0 - this->lower_spring_length[1] * inv_norm_r_low2);
		cblas_dscal(this->DIM, scale, lower_force2, 1);
		cblas_dcopy(this->DIM, r_low3, 1, lower_force3, 1);
		scale = this->lower_spring_stiffness[2] * (1.0 - this->lower_spring_length[2] * inv_norm_r_low3);
		cblas_dscal(this->DIM, scale, lower_force3, 1);
		cblas_dcopy(this->DIM, r_low4, 1, lower_force4, 1);
		scale = this->lower_spring_stiffness[3] * (1.0 - this->lower_spring_length[3] * inv_norm_r_low4);
		cblas_dscal(this->DIM, scale, lower_force4, 1);

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
		T* upper_dampf1 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		T* upper_dampf2 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		T* upper_dampf3 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		T* upper_dampf4 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		T* lower_dampf1 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		T* lower_dampf2 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		T* lower_dampf3 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		T* lower_dampf4 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		T* temp = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);

		//// upper_dampf1
		// compute: vc - C_cN' * (r1_tilda * wc))
		cblas_dgemv(CblasRowMajor, CblasNoTrans, this->DIM, this->DIM, 1.0, this->r1_tilda, this->DIM, wc_, 1, 0.0, temp, 1);
		cblas_dcopy(this->DIM, vc_, 1, upper_dampf1, 1);
		cblas_dgemv(CblasRowMajor, CblasNoTrans, this->DIM, this->DIM, -1.0, C_cN, this->DIM, temp, 1, 1.0, upper_dampf1, 1);
		// dot((vc - C_cN' * (r1_tilda * wc)), r_up1)
		//std::cout << cblas_ddot(this->DIM, upper_dampf1, 1, r_up1, 1) << std::endl;
		//std::cout << MathLibrary::dot_product<T>(upper_dampf1, r_up1, this->DIM) << std::endl;
		scale = cblas_ddot(mkl_DIM, upper_dampf1, mkl_incx, r_up1, mkl_incy);
		//scale = MathLibrary::dot_product<T>(upper_dampf1, r_up1, this->DIM);
		
		// dot((vc - C_cN' * (r1_tilda * wc)), r_up1) - dot(vw1, r_up1)
		scale -= cblas_ddot(mkl_DIM, vw1_, mkl_incx, r_up1, mkl_incx);
		//scale -= MathLibrary::dot_product<T>(vw1_, r_up1, this->DIM);
		// (dot((vc - C_cN' * (r1_tilda * wc)), r_up1) - dot(vw1, r_up1))* inv_norm_r_up1 * inv_norm_r_up1
		scale = scale * inv_norm_r_up1*inv_norm_r_up1 * this->upper_spring_damping[0];
		cblas_dcopy(this->DIM, r_up1, 1, upper_dampf1, 1);
		cblas_dscal(this->DIM, scale, upper_dampf1, 1);
		

		//// upper_dampf2
		// compute: vc - C_cN' * (r2_tilda * wc))
		cblas_dgemv(CblasRowMajor, CblasNoTrans, this->DIM, this->DIM, 1.0, this->r2_tilda, this->DIM, wc_, 1, 0.0, temp, 1);
		cblas_dcopy(this->DIM, vc_, 1, upper_dampf2, 1);
		cblas_dgemv(CblasRowMajor, CblasNoTrans, this->DIM, this->DIM, -1.0, C_cN, this->DIM, temp, 1, 1.0, upper_dampf2, 1);
		// dot((vc - C_cN' * (r2_tilda * wc)), r_up2)
		scale = cblas_ddot(mkl_DIM, upper_dampf2, mkl_incx, r_up2, mkl_incy);
		//scale = MathLibrary::dot_product<T>(upper_dampf2, r_up2, this->DIM);
		// dot((vc - C_cN' * (r2_tilda * wc)), r_up2) - dot(vw2, r_up2)
		scale -= cblas_ddot(mkl_DIM, vw2_, mkl_incx, r_up2, mkl_incy);
		//scale -= MathLibrary::dot_product<T>(vw2_, r_up2, this->DIM);
		// (dot((vc - C_cN' * (r2_tilda * wc)), r_up2) - dot(vw2, r_up2))* inv_norm_r_up2 * inv_norm_r_up2
		scale = scale * inv_norm_r_up2*inv_norm_r_up2* this->upper_spring_damping[1];
		cblas_dcopy(this->DIM, r_up2, 1, upper_dampf2, 1);
		cblas_dscal(this->DIM, scale, upper_dampf2, 1);

		//// upper_dampf3
		// compute: vc - C_cN' * (r3_tilda * wc))
		cblas_dgemv(CblasRowMajor, CblasNoTrans, this->DIM, this->DIM, 1.0, this->r3_tilda, this->DIM, wc_, 1, 0.0, temp, 1);
		cblas_dcopy(this->DIM, vc_, 1, upper_dampf3, 1);
		cblas_dgemv(CblasRowMajor, CblasNoTrans, this->DIM, this->DIM, -1.0, C_cN, this->DIM, temp, 1, 1.0, upper_dampf3, 1);
		// dot((vc - C_cN' * (r3_tilda * wc)), r_up3)
		scale = cblas_ddot(mkl_DIM, upper_dampf3, mkl_incx, r_up3, mkl_incy);
		//scale = MathLibrary::dot_product<T>(upper_dampf3, r_up3, this->DIM);
		// dot((vc - C_cN' * (r3_tilda * wc)), r_up3) - dot(vw3, r_up3)
		scale -= cblas_ddot(mkl_DIM, vw3_, mkl_incx, r_up3, mkl_incy);
		//scale -= MathLibrary::dot_product<T>(vw3_, r_up3, this->DIM);
		// (dot((vc - C_cN' * (r3_tilda * wc)), r_up3) - dot(vw3, r_up3))* inv_norm_r_up3 * inv_norm_r_up3
		scale = scale * inv_norm_r_up3*inv_norm_r_up3* this->upper_spring_damping[2];
		cblas_dcopy(this->DIM, r_up3, 1, upper_dampf3, 1);
		cblas_dscal(this->DIM, scale, upper_dampf3, 1);

		//// upper_dampf4
		// compute: vc - C_cN' * (r4_tilda * wc))
		cblas_dgemv(CblasRowMajor, CblasNoTrans, this->DIM, this->DIM, 1.0, this->r4_tilda, this->DIM, wc_, 1, 0.0, temp, 1);
		cblas_dcopy(this->DIM, vc_, 1, upper_dampf4, 1);
		cblas_dgemv(CblasRowMajor, CblasNoTrans, this->DIM, this->DIM, -1.0, C_cN, this->DIM, temp, 1, 1.0, upper_dampf4, 1);
		// dot((vc - C_cN' * (r4_tilda * wc)), r_up4)
		scale = cblas_ddot(mkl_DIM, upper_dampf4, mkl_incx, r_up4, mkl_incy);
		//scale = MathLibrary::dot_product<T>(upper_dampf4, r_up4, this->DIM);
		// dot((vc - C_cN' * (r4_tilda * wc)), r_up4) - dot(vw4, r_up4)
		scale -= cblas_ddot(mkl_DIM, vw4_, mkl_incx, r_up4, mkl_incy);
		//scale -= MathLibrary::dot_product<T>(r_up4, vw4_, this->DIM);
		// (dot((vc - C_cN' * (r4_tilda * wc)), r_up4) - dot(vw4, r_up4))* inv_norm_r_up4 * inv_norm_r_up4
		scale = scale * inv_norm_r_up4*inv_norm_r_up4* this->upper_spring_damping[3];
		cblas_dcopy(this->DIM, r_up4, 1, upper_dampf4, 1);
		cblas_dscal(this->DIM, scale, upper_dampf4, 1);

		//// lower_dampf1
		// dot(vw1, r_low1)
		scale = cblas_ddot(mkl_DIM, vw1_, mkl_incx, r_low1, mkl_incy);
		/*scale = MathLibrary::dot_product<T>(vw1_, r_low1, this->DIM);*/
		// (dot(vw1, r_low1) - dot(vt1, r_low1))
		scale -= cblas_ddot(mkl_DIM, vt1_, mkl_incx, r_low1, mkl_incx);
		//scale -= MathLibrary::dot_product<T>(vt1_, r_low1, this->DIM);
		//(dot(vw1, r_low1) - dot(vt1, r_low1)) * inv_norm_r_low1 * inv_norm_r_low1
		scale = scale * inv_norm_r_low1*inv_norm_r_low1* this->lower_spring_damping[0];
		cblas_dcopy(this->DIM, r_low1, 1, lower_dampf1, 1);
		cblas_dscal(this->DIM, scale, lower_dampf1, 1);

		//// lower_dampf2
		// dot(vw2, r_low2)
		scale = cblas_ddot(mkl_DIM, vw2_, mkl_incx, r_low2, mkl_incy);
		//scale = MathLibrary::dot_product<T>(vw2_, r_low2, this->DIM);
		// (dot(vw2, r_low2) - dot(vt2, r_low2))
		scale -= cblas_ddot(mkl_DIM, vt2_, mkl_incx, r_low2, mkl_incy);
		//scale -= MathLibrary::dot_product<T>(vt2_, r_low2, this->DIM);
		//(dot(vw2, r_low2) - dot(vt2, r_low2)) * inv_norm_r_low2 * inv_norm_r_low2
		scale = scale * inv_norm_r_low2*inv_norm_r_low2* this->lower_spring_damping[1];
		cblas_dcopy(this->DIM, r_low2, 1, lower_dampf2, 1);
		cblas_dscal(this->DIM, scale, lower_dampf2, 1);

		//// lower_dampf3
		// dot(vw3, r_low3)
		scale = cblas_ddot(mkl_DIM, vw3_, mkl_incx, r_low3, mkl_incy);
		//scale = MathLibrary::dot_product<T>(vw3_, r_low3, this->DIM);
		// (dot(vw3, r_low3) - dot(vt3, r_low3))
		scale -= cblas_ddot(mkl_DIM, vt3_, mkl_incx, r_low3, mkl_incy);
		//scale -= MathLibrary::dot_product<T>(vt3_, r_low3, this->DIM);
		//(dot(vw3, r_low3) - dot(vt3, r_low3)) * inv_norm_r_low3 * inv_norm_r_low3
		scale = scale * inv_norm_r_low3*inv_norm_r_low3* this->lower_spring_damping[2];
		cblas_dcopy(this->DIM, r_low3, 1, lower_dampf3, 1);
		cblas_dscal(this->DIM, scale, lower_dampf3, 1);

		//// lower_dampf4
		// dot(vw4, r_low4)
		scale = cblas_ddot(mkl_DIM, vw4_, mkl_incx, r_low4, mkl_incy);
		//scale = MathLibrary::dot_product<T>(vw4_, r_low4, this->DIM);
		// (dot(vw4, r_low4) - dot(vt4, r_low4))
		scale -= cblas_ddot(mkl_DIM, vt4_, mkl_incx, r_low4, mkl_incy);
		//scale -= MathLibrary::dot_product<T>(vt4_, r_low4, this->DIM);
		//(dot(vw4, r_low4) - dot(vt4, r_low4)) * inv_norm_r_low4 * inv_norm_r_low4
		scale = scale * inv_norm_r_low4*inv_norm_r_low4* this->lower_spring_damping[3];
		cblas_dcopy(this->DIM, r_low4, 1, lower_dampf4, 1);
		cblas_dscal(this->DIM, scale, lower_dampf4, 1);

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

		T* upper_S1 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		T* upper_S2 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		T* upper_S3 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		T* upper_S4 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		T* lower_S1 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		T* lower_S2 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		T* lower_S3 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		T* lower_S4 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);

		// upper_S1 = upper_rotational_stiffness(1) * upper_angle1 * upper_normal1;
		scale = (this->upper_rotational_stiffness[0]) * (*upper_angle1);
		cblas_dcopy(this->DIM, upper_normal1, 1, upper_S1, 1);
		cblas_dscal(this->DIM, scale, upper_S1, 1);
		

		// upper_S2 = upper_rotational_stiffness(2) * upper_angle2 * upper_normal2;
		scale = (this->upper_rotational_stiffness[1]) * (*upper_angle2);
		cblas_dcopy(this->DIM, upper_normal2, 1, upper_S2, 1);
		cblas_dscal(this->DIM, scale, upper_S2, 1);

		// upper_S3 = upper_rotational_stiffness(3) * upper_angle3 * upper_normal3;
		scale = (this->upper_rotational_stiffness[2]) * (*upper_angle3);
		cblas_dcopy(this->DIM, upper_normal3, 1, upper_S3, 1);
		cblas_dscal(this->DIM, scale, upper_S3, 1);

		// upper_S4 = upper_rotational_stiffness(4) * upper_angle4 * upper_normal4;
		scale = (this->upper_rotational_stiffness[3]) * (*upper_angle4);
		cblas_dcopy(this->DIM, upper_normal4, 1, upper_S4, 1);
		cblas_dscal(this->DIM, scale, upper_S4, 1);

		// lower_S1 = lower_rotational_stiffness(1) * lower_angle1 * lower_normal1
		scale = (this->lower_rotational_stiffness[0]) * (*lower_angle1);
		cblas_dcopy(this->DIM, lower_normal1, 1, lower_S1, 1);
		cblas_dscal(this->DIM, scale, lower_S1, 1);

		// lower_S2 = lower_rotational_stiffness(2) * lower_angle2 * lower_normal2
		scale = (this->lower_rotational_stiffness[1]) * (*lower_angle2);
		cblas_dcopy(this->DIM, lower_normal2, 1, lower_S2, 1);
		cblas_dscal(this->DIM, scale, lower_S2, 1);

		// lower_S3 = lower_rotational_stiffness(3) * lower_angle3 * lower_normal3
		scale = (this->lower_rotational_stiffness[2]) * (*lower_angle3);
		cblas_dcopy(this->DIM, lower_normal3, 1, lower_S3, 1);
		cblas_dscal(this->DIM, scale, lower_S3, 1);

		// lower_S4 = lower_rotational_stiffness(4) * lower_angle4 * lower_normal4
		scale = (this->lower_rotational_stiffness[3]) * (*lower_angle4);
		cblas_dcopy(this->DIM, lower_normal4, 1, lower_S4, 1);
		cblas_dscal(this->DIM, scale, lower_S4, 1);

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
		T* lower_rot_force1 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		T* lower_rot_force2 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		T* lower_rot_force3 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		T* lower_rot_force4 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		T* upper_rot_force1 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		T* upper_rot_force2 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		T* upper_rot_force3 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		T* upper_rot_force4 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		T* car_rot_force1 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		T* car_rot_force2 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		T* car_rot_force3 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		T* car_rot_force4 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		T* sum_car_force1 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		T* sum_car_force2 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		T* sum_car_force3 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		T* sum_car_force4 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		T scale_u1, scale_u2, scale_u3, scale_u4;
		// lower_rot_force1 = -cross( lower_S1, r_low1) / (r_low1'*r_low1);
		scale = -1.0 / cblas_ddot(mkl_DIM, r_low1, mkl_incx, r_low1, mkl_incy);
		//scale = -1.0 / MathLibrary::dot_product<T>(r_low1, r_low1, this->DIM);
		MathLibrary::crossProduct(lower_S1, r_low1, lower_rot_force1);
		cblas_dscal(this->DIM, scale, lower_rot_force1, 1);

		// lower_rot_force2 = -cross( lower_S2, r_low2) / (r_low2'*r_low2);
		scale = -1.0 / cblas_ddot(mkl_DIM, r_low2, mkl_incx, r_low2, mkl_incy);
		//scale = -1.0 / MathLibrary::dot_product<T>(r_low2, r_low2, this->DIM);
		MathLibrary::crossProduct(lower_S2, r_low2, lower_rot_force2);
		cblas_dscal(this->DIM, scale, lower_rot_force2, 1);

		// lower_rot_force3 = -cross( lower_S3, r_low3) / (r_low3'*r_low3);
		scale = -1.0 / cblas_ddot(mkl_DIM, r_low3, mkl_incx, r_low3, mkl_incy);
		//scale = -1.0 / MathLibrary::dot_product<T>(r_low3, r_low3, this->DIM);
		MathLibrary::crossProduct(lower_S3, r_low3, lower_rot_force3);
		cblas_dscal(this->DIM, scale, lower_rot_force3, 1);

		// lower_rot_force4 = -cross( lower_S4, r_low4) / (r_low4'*r_low4);
		scale = -1.0 / cblas_ddot(mkl_DIM, r_low4, mkl_incx, r_low4, mkl_incx);
		//scale = -1.0 / MathLibrary::dot_product<T>(r_low4, r_low4, this->DIM);
		MathLibrary::crossProduct(lower_S4, r_low4, lower_rot_force4);
		cblas_dscal(this->DIM, scale, lower_rot_force4, 1);

		// upper_rot_force1 = -cross( upper_S1, r_up1) / (r_up1'*r_up1);
		scale_u1 = -1.0 / cblas_ddot(mkl_DIM, r_up1, mkl_incx, r_up1, mkl_incy);
		//scale_u1 = -1.0 / MathLibrary::dot_product<T>(r_up1, r_up1, this->DIM);
		MathLibrary::crossProduct(upper_S1, r_up1, upper_rot_force1);
		cblas_dscal(this->DIM, scale_u1, upper_rot_force1, 1);

		// upper_rot_force2 = -cross( upper_S2, r_up2) / (r_up2'*r_up2);
		scale_u2 = -1.0 / cblas_ddot(mkl_DIM, r_up2, mkl_incx, r_up2, mkl_incy);
		//scale_u2 = -1.0 / MathLibrary::dot_product<T>(r_up2, r_up2, this->DIM);
		MathLibrary::crossProduct(upper_S2, r_up2, upper_rot_force2);
		cblas_dscal(this->DIM, scale_u2, upper_rot_force2, 1);
		
		// upper_rot_force3 = -cross( upper_S3, r_up3) / (r_up3'*r_up3);
		scale_u3 = -1.0 / cblas_ddot(mkl_DIM, r_up3, mkl_incx, r_up3, mkl_incy);
		//scale_u3 = -1.0 / MathLibrary::dot_product<T>(r_up3, r_up3, this->DIM);
		MathLibrary::crossProduct(upper_S3, r_up3, upper_rot_force3);
		cblas_dscal(this->DIM, scale_u3, upper_rot_force3, 1);

		// upper_rot_force4 = -cross( upper_S4, r_up4) / (r_up4'*r_up4);
		scale_u4 = -1.0 / cblas_ddot(mkl_DIM, r_up4, mkl_incx, r_up4, mkl_incy);
		//scale_u4 = -1.0 / MathLibrary::dot_product<T>(r_up4, r_up4, this->DIM);
		MathLibrary::crossProduct(upper_S4, r_up4, upper_rot_force4);
		cblas_dscal(this->DIM, scale_u4, upper_rot_force4, 1);
		
		// car_rot_force1 = -cross( lower_S1, r_up1) / (r_up1'*r_up1);
		MathLibrary::crossProduct(lower_S1, r_up1, car_rot_force1);
		cblas_dscal(this->DIM, scale_u1, car_rot_force1, 1);

		// car_rot_force2 = -cross( lower_S2, r_up2) / (r_up2'*r_up2);
		MathLibrary::crossProduct(lower_S2, r_up2, car_rot_force2);
		cblas_dscal(this->DIM, scale_u2, car_rot_force2, 1);

		// car_rot_force3 = -cross( lower_S3, r_up3) / (r_up3'*r_up3);
		MathLibrary::crossProduct(lower_S3, r_up3, car_rot_force3);
		cblas_dscal(this->DIM, scale_u3, car_rot_force3, 1);

		// car_rot_force4 = -cross( lower_S4, r_up4) / (r_up4'*r_up4);
		MathLibrary::crossProduct(lower_S4, r_up4, car_rot_force4);
		cblas_dscal(this->DIM, scale_u4, car_rot_force4, 1);

		// sum_car_force1 = car_rot_force1 - upper_force1 - upper_dampf1 - upper_rot_force1;
		cblas_dcopy(this->DIM, car_rot_force1, 1, sum_car_force1, 1);
		cblas_daxpy(this->DIM, -1.0, upper_force1, 1, sum_car_force1, 1);
		cblas_daxpy(this->DIM, -1.0, upper_dampf1, 1, sum_car_force1, 1);
		cblas_daxpy(this->DIM, -1.0, upper_rot_force1, 1, sum_car_force1, 1);

		// sum_car_force2 = car_rot_force2 - upper_force2 - upper_dampf2 - upper_rot_force2;
		cblas_dcopy(this->DIM, car_rot_force2, 1, sum_car_force2, 1);
		cblas_daxpy(this->DIM, -1.0, upper_force2, 1, sum_car_force2, 1);
		cblas_daxpy(this->DIM, -1.0, upper_dampf2, 1, sum_car_force2, 1);
		cblas_daxpy(this->DIM, -1.0, upper_rot_force2, 1, sum_car_force2, 1);

		// sum_car_force3 = car_rot_force3 - upper_force3 - upper_dampf3 - upper_rot_force3;
		cblas_dcopy(this->DIM, car_rot_force3, 1, sum_car_force3, 1);
		cblas_daxpy(this->DIM, -1.0, upper_force3, 1, sum_car_force3, 1);
		cblas_daxpy(this->DIM, -1.0, upper_dampf3, 1, sum_car_force3, 1);
		cblas_daxpy(this->DIM, -1.0, upper_rot_force3, 1, sum_car_force3, 1);

		// sum_car_force4 = car_rot_force4 - upper_force4 - upper_dampf4 - upper_rot_force4;
		cblas_dcopy(this->DIM, car_rot_force4, 1, sum_car_force4, 1);
		cblas_daxpy(this->DIM, -1.0, upper_force4, 1, sum_car_force4, 1);
		cblas_daxpy(this->DIM, -1.0, upper_dampf4, 1, sum_car_force4, 1);
		cblas_daxpy(this->DIM, -1.0, upper_rot_force4, 1, sum_car_force4, 1);

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

		road forces    
		local_FR1 = [0;  aux_vals.FR1(t, pt1(2), vt1(2), A(19,19), lower_force1(2) + FT(1) + lower_rot_force1(2)); 0];      %in global basis
		local_FR2 = [0;  aux_vals.FR2(t, pt2(2), vt2(2), A(22,22), lower_force2(2) + FT(2) + lower_rot_force2(2)); 0];   
		local_FR3 = [0;  aux_vals.FR3(t, pt3(2), vt3(2), A(25,25), lower_force3(2) + FT(3) + lower_rot_force3(2)); 0]; 
		local_FR4 = [0;  aux_vals.FR4(t, pt4(2), vt4(2), A(28,28), lower_force4(2) + FT(4) + lower_rot_force4(2)); 0]; 
		*/
		T* local_FR1 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		T* local_FR2 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		T* local_FR3 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		T* local_FR4 = (T*)mkl_calloc((this->DIM), sizeof(T), this->alignment);
		
		// local_FR1
		local_FR1[0] = 0;
		get_road_force(t_, pt1_[1], vt1_[1], this->mass_tyre[0], lower_force1[1] + FT[0] + lower_rot_force1[1], local_FR1 + 1);
		local_FR1[2] = 0;

		// local_FR2
		local_FR2[0] = 0;
		get_road_force(t_, pt2_[1], vt2_[1], this->mass_tyre[1], lower_force2[1] + FT[1] + lower_rot_force2[1], local_FR2 + 1);
		local_FR2[2] = 0;

		// local_FR3
		local_FR3[0] = 0;
		get_road_force(t_, pt3_[1], vt3_[1], this->mass_tyre[2], lower_force3[1] + FT[2] + lower_rot_force3[1], local_FR3 + 1);
		local_FR3[2] = 0;

		// local_FR4
		local_FR4[0] = 0;
		get_road_force(t_, pt4_[1], vt4_[1], this->mass_tyre[3], lower_force4[1] + FT[3] + lower_rot_force4[1], local_FR4 + 1);
		local_FR4[2] = 0;

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
		T* Hc = (T*)mkl_calloc(this->DIM, sizeof(T), this->alignment);
		T* sum_torque_spring_car = (T*)mkl_calloc(this->DIM, sizeof(T), this->alignment);
		T* Tc = (T*)mkl_calloc(this->DIM, sizeof(T), this->alignment);
		T* wc_tilda = (T*)mkl_calloc((this->DIM) * (this->DIM), sizeof(T), this->alignment);
		
		// Hc = A(1:3, 1:3) * wc;
		cblas_dgemv(CblasRowMajor, CblasNoTrans, this->DIM, this->DIM, 1.0, this->Ic, this->DIM, wc_, 1, 0.0, Hc, 1);
		get_tilda(wc_, wc_tilda);
		cblas_dcopy(this->DIM, Hc, 1, temp, 1);
		cblas_dgemv(CblasRowMajor, CblasNoTrans, this->DIM, this->DIM, -1.0, wc_tilda, this->DIM, temp, 1, 0.0, Hc, 1);
		cblas_dcopy(this->DIM, Hc, 1, sum_torque_spring_car, 1);
		cblas_daxpy(this->DIM, 1.0, Tc, 1, sum_torque_spring_car, 1);
		cblas_dgemv(CblasRowMajor, CblasTrans, this->DIM, this->DIM, -1.0, C_cN, this->DIM, upper_S4, 1, 1.0, sum_torque_spring_car, 1);
		cblas_dgemv(CblasRowMajor, CblasTrans, this->DIM, this->DIM, -1.0, C_cN, this->DIM, upper_S3, 1, 1.0, sum_torque_spring_car, 1);
		cblas_dgemv(CblasRowMajor, CblasTrans, this->DIM, this->DIM, -1.0, C_cN, this->DIM, upper_S2, 1, 1.0, sum_torque_spring_car, 1);
		cblas_dgemv(CblasRowMajor, CblasTrans, this->DIM, this->DIM, -1.0, C_cN, this->DIM, upper_S1, 1, 1.0, sum_torque_spring_car, 1);
		cblas_dgemv(CblasRowMajor, CblasTrans, this->DIM, this->DIM, 1.0, C_cN, this->DIM, sum_car_force1, 1, 0.0, temp, 1);
		cblas_dgemv(CblasRowMajor, CblasNoTrans, this->DIM, this->DIM, 1.0, this->r1_tilda, this->DIM, temp, 1, 1.0, sum_torque_spring_car, 1);
		cblas_dgemv(CblasRowMajor, CblasTrans, this->DIM, this->DIM, 1.0, C_cN, this->DIM, sum_car_force2, 1, 0.0, temp, 1);
		cblas_dgemv(CblasRowMajor, CblasNoTrans, this->DIM, this->DIM, 1.0, this->r2_tilda, this->DIM, temp, 1, 1.0, sum_torque_spring_car, 1);
		cblas_dgemv(CblasRowMajor, CblasTrans, this->DIM, this->DIM, 1.0, C_cN, this->DIM, sum_car_force3, 1, 0.0, temp, 1);
		cblas_dgemv(CblasRowMajor, CblasNoTrans, this->DIM, this->DIM, 1.0, this->r3_tilda, this->DIM, temp, 1, 1.0, sum_torque_spring_car, 1);
		cblas_dgemv(CblasRowMajor, CblasTrans, this->DIM, this->DIM, 1.0, C_cN, this->DIM, sum_car_force4, 1, 0.0, temp, 1);
		cblas_dgemv(CblasRowMajor, CblasNoTrans, this->DIM, this->DIM, 1.0, this->r4_tilda, this->DIM, temp, 1, 1.0, sum_torque_spring_car, 1);
		
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
		T* b_rem = (T*)mkl_calloc(this->DIM, sizeof(T), this->alignment);
		// FC + sum_car_force1 + sum_car_force2 + sum_car_force3 + sum_car_force4; ...          %vc_dot  
		T* brem_start = b_rem;
		cblas_dcopy(this->DIM, this->FC, 1, brem_start, 1);
		cblas_daxpy(this->DIM, 1.0, sum_car_force1, 1, brem_start, 1);
		cblas_daxpy(this->DIM, 1.0, sum_car_force2, 1, brem_start, 1);
		cblas_daxpy(this->DIM, 1.0, sum_car_force3, 1, brem_start, 1);
		cblas_daxpy(this->DIM, 1.0, sum_car_force4, 1, brem_start, 1);

		//  upper_force1 - lower_force1 + upper_dampf1 - lower_dampf1 + local_FW1 + upper_rot_force1 - car_rot_force1 - lower_rot_force1;
		brem_start += this->DIM;
		cblas_dcopy(this->DIM, upper_force1, 1, brem_start, 1);
		cblas_daxpy(this->DIM, -1.0, lower_force1, 1, brem_start, 1);
		cblas_daxpy(this->DIM, 1.0, upper_dampf1, 1, brem_start, 1);
		cblas_daxpy(this->DIM, -1.0, lower_dampf1, 1, brem_start, 1);
		//cblas_daxpy(this->DIM, 1.0, local_FW1, 1, brem_start, 1);
		cblas_daxpy(this->DIM, 1.0, upper_rot_force1, 1, brem_start, 1);
		cblas_daxpy(this->DIM, -1.0, car_rot_force1, 1, brem_start, 1);
		cblas_daxpy(this->DIM, -1.0, lower_rot_force1, 1, brem_start, 1);

		// upper_force2 - lower_force2 + upper_dampf2 - lower_dampf2 + local_FW2 + upper_rot_force2 - car_rot_force2 - lower_rot_force2;
		brem_start += this->DIM;
		cblas_dcopy(this->DIM, upper_force2, 1, brem_start, 1);
		cblas_daxpy(this->DIM, -1.0, lower_force2, 1, brem_start, 1);
		cblas_daxpy(this->DIM, 1.0, upper_dampf2, 1, brem_start, 1);
		cblas_daxpy(this->DIM, -1.0, lower_dampf2, 1, brem_start, 1);
		//cblas_daxpy(this->DIM, 1.0, local_FW2, 1, brem_start, 1);
		cblas_daxpy(this->DIM, 1.0, upper_rot_force2, 1, brem_start, 1);
		cblas_daxpy(this->DIM, -1.0, car_rot_force2, 1, brem_start, 1);
		cblas_daxpy(this->DIM, -1.0, lower_rot_force2, 1, brem_start, 1);

		// upper_force3 - lower_force3 + upper_dampf3 - lower_dampf3 + local_FW3 + upper_rot_force3 - car_rot_force3 - lower_rot_force3;
		brem_start += this->DIM;
		cblas_dcopy(this->DIM, upper_force3, 1, brem_start, 1);
		cblas_daxpy(this->DIM, -1.0, lower_force3, 1, brem_start, 1);
		cblas_daxpy(this->DIM, 1.0, upper_dampf3, 1, brem_start, 1);
		cblas_daxpy(this->DIM, -1.0, lower_dampf3, 1, brem_start, 1);
		//cblas_daxpy(this->DIM, 1.0, local_FW3, 1, brem_start, 1);
		cblas_daxpy(this->DIM, 1.0, upper_rot_force3, 1, brem_start, 1);
		cblas_daxpy(this->DIM, -1.0, car_rot_force3, 1, brem_start, 1);
		cblas_daxpy(this->DIM, -1.0, lower_rot_force3, 1, brem_start, 1);

		//  upper_force4 - lower_force4 + upper_dampf4 - lower_dampf4 + local_FW4 + upper_rot_force4 - car_rot_force4 - lower_rot_force4;
		brem_start += this->DIM;
		cblas_dcopy(this->DIM, upper_force4, 1, brem_start, 1);
		cblas_daxpy(this->DIM, -1.0, lower_force4, 1, brem_start, 1);
		cblas_daxpy(this->DIM, 1.0, upper_dampf4, 1, brem_start, 1);
		cblas_daxpy(this->DIM, -1.0, lower_dampf4, 1, brem_start, 1);
		//cblas_daxpy(this->DIM, 1.0, local_FW4, 1, brem_start, 1);
		cblas_daxpy(this->DIM, 1.0, upper_rot_force4, 1, brem_start, 1);
		cblas_daxpy(this->DIM, -1.0, car_rot_force4, 1, brem_start, 1);
		cblas_daxpy(this->DIM, -1.0, lower_rot_force4, 1, brem_start, 1);

		// lower_force1 + lower_dampf1 + local_FT1 + local_FR1 + lower_rot_force1; ...          %vt1_dot
		brem_start += this->DIM;
		cblas_dcopy(this->DIM, lower_force1, 1, brem_start, 1);
		cblas_daxpy(this->DIM, 1.0, lower_dampf1, 1, brem_start, 1);
	//	cblas_daxpy(this->DIM, 1.0, local_FT1, 1, brem_start, 1); not implemented since not used now
		cblas_daxpy(this->DIM, 1.0, local_FR1, 1, brem_start, 1);
		cblas_daxpy(this->DIM, 1.0, lower_rot_force1, 1, brem_start, 1);

		// lower_force2 + lower_dampf2 + local_FT2 + local_FR2 + lower_rot_force2; ...          %vt2_dot
		brem_start += this->DIM;
		cblas_dcopy(this->DIM, lower_force2, 1, brem_start, 1);
		cblas_daxpy(this->DIM, 1.0, lower_dampf2, 1, brem_start, 1);
		//	cblas_daxpy(this->DIM, 1.0, local_FT2, 1, brem_start, 1); not implemented since not used now
		cblas_daxpy(this->DIM, 1.0, local_FR2, 1, brem_start, 1);
		cblas_daxpy(this->DIM, 1.0, lower_rot_force2, 1, brem_start, 1);

		// lower_force3 + lower_dampf3 + local_FT3 + local_FR3 + lower_rot_force3; ...          %vt3_dot
		brem_start += this->DIM;
		cblas_dcopy(this->DIM, lower_force3, 1, brem_start, 1);
		cblas_daxpy(this->DIM, 1.0, lower_dampf3, 1, brem_start, 1);
		//	cblas_daxpy(this->DIM, 1.0, local_FT3, 1, brem_start, 1); not implemented since not used now
		cblas_daxpy(this->DIM, 1.0, local_FR3, 1, brem_start, 1);
		cblas_daxpy(this->DIM, 1.0, lower_rot_force3, 1, brem_start, 1);

		// lower_force4 + lower_dampf4 + local_FT4 + local_FR4 + lower_rot_force4];             %vt4_dot
		brem_start += this->DIM;
		cblas_dcopy(this->DIM, lower_force4, 1, brem_start, 1);
		cblas_daxpy(this->DIM, 1.0, lower_dampf4, 1, brem_start, 1);
		//	cblas_daxpy(this->DIM, 1.0, local_FT4, 1, brem_start, 1); not implemented since not used now
		cblas_daxpy(this->DIM, 1.0, local_FR4, 1, brem_start, 1);
		cblas_daxpy(this->DIM, 1.0, lower_rot_force4, 1, brem_start, 1);

		LAPACKE_dpotrs(LAPACK_ROW_MAJOR, 'L', this->DIM, 1, A_Ic, this->DIM, sum_torque_spring_car, 1);
		cblas_dcopy(this->DIM, sum_torque_spring_car, 1, f_, 1);
		T* start_next = f_ + this->DIM;

		MathLibrary::vector_elem_wise_product<T>(A_rem, b_rem, start_next, 9 * this->DIM);
		start_next += 9 * this->DIM;
		
		/*
		get the derivative of the altitude (expressed in quaternions) from the angular velocities
		Qc = 0.5 * [qc(4) -qc(3) qc(2); qc(3) qc(4) -qc(1); -qc(2) qc(1) qc(4); -qc(1) -qc(2) -qc(3)];
		qc_dot = Qc * wc;
		*/
		T* Qc = (T*)mkl_calloc((this->NUM_LEGS)*(this->DIM), sizeof(T), this->alignment);
		T* qc_dot = (T*)mkl_calloc((this->NUM_LEGS), sizeof(T), this->alignment);
		Qc[0] = 0.5 * qc_[3];
		Qc[1] = -0.5 * qc_[2];
		Qc[2] = 0.5 * qc_[1];
		Qc[3] = 0.5 * qc_[2];
		Qc[4] = 0.5 * qc_[3];
		Qc[5] = -0.5 * qc_[0];
		Qc[6] = -0.5 * qc_[1];
		Qc[7] = 0.5 * qc_[0];
		Qc[8] = 0.5 * qc_[3];
		Qc[9] = -0.5 * qc_[0];
		Qc[10] = -0.5 * qc_[1];
		Qc[11] = -0.5 * qc_[2];
		cblas_dgemv(CblasRowMajor, CblasNoTrans, this->NUM_LEGS, this->DIM, 1.0, Qc, this->DIM, wc_, 1, 0.0, qc_dot, 1);
		
		cblas_dcopy(this->NUM_LEGS, qc_dot, 1, start_next, 1);
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
		
		// Clean allocated memory
		mkl_free(C_cN);
		mkl_free(r_up1);
		mkl_free(r_up2);
		mkl_free(r_up3);
		mkl_free(r_up4);
		mkl_free(r_low1);
		mkl_free(r_low2);
		mkl_free(r_low3);
		mkl_free(r_low4);
		mkl_free(upper_normal1);
		mkl_free(upper_normal2);
		mkl_free(upper_normal3);
		mkl_free(upper_normal4);
		mkl_free(lower_normal1);
		mkl_free(lower_normal2);
		mkl_free(lower_normal3);
		mkl_free(lower_normal4);
		mkl_free(col_dat);
		mkl_free(upper_force1);
		mkl_free(upper_force2);
		mkl_free(upper_force3);
		mkl_free(upper_force4);
		mkl_free(lower_force1);
		mkl_free(lower_force2);
		mkl_free(lower_force3);
		mkl_free(lower_force4);
		mkl_free(upper_dampf1);
		mkl_free(upper_dampf2);
		mkl_free(upper_dampf3);
		mkl_free(upper_dampf4);
		mkl_free(lower_dampf1);
		mkl_free(lower_dampf2);
		mkl_free(lower_dampf3);
		mkl_free(lower_dampf4);
		mkl_free(temp);
		mkl_free(upper_S1);
		mkl_free(upper_S2);
		mkl_free(upper_S3);
		mkl_free(upper_S4);
		mkl_free(lower_S1);
		mkl_free(lower_S2);
		mkl_free(lower_S3);
		mkl_free(lower_S4);
		mkl_free(lower_rot_force1);
		mkl_free(lower_rot_force2);
		mkl_free(lower_rot_force3);
		mkl_free(lower_rot_force4);
		mkl_free(upper_rot_force2);
		mkl_free(upper_rot_force1);
		mkl_free(upper_rot_force3);
		mkl_free(upper_rot_force4);
		mkl_free(car_rot_force1);
		mkl_free(car_rot_force2);
		mkl_free(car_rot_force3);
		mkl_free(car_rot_force4);
		mkl_free(sum_car_force1);
		mkl_free(sum_car_force2);
		mkl_free(sum_car_force3);
		mkl_free(sum_car_force4);
		mkl_free(local_FR1);
		mkl_free(local_FR2);
		mkl_free(local_FR3);
		mkl_free(local_FR4);
		mkl_free(Hc);
		mkl_free(sum_torque_spring_car);
		mkl_free(Tc);
		mkl_free(wc_tilda);
		mkl_free(b_rem);
		mkl_free(Qc);
		mkl_free(qc_dot);
	}

	void solve(T* solution_vector) {
		// From the formulation we have 61 dimensions in the solution vector
		std::cout << "Solver Triggered!" << std::endl;
		size_t solution_size = (this->num_iter+1) * this->solution_dim;
		T* complete_vector = (T*)mkl_calloc(solution_size, sizeof(T), this->alignment);
		x_vector = (T*)mkl_calloc(solution_dim, sizeof(T), this->alignment);
		get_tilda(r1, r1_tilda);
		get_tilda(r2, r2_tilda);
		get_tilda(r3, r3_tilda);
		get_tilda(r4, r4_tilda);
		
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
		// qc
		T *start_orient = x_vector + i * (this->DIM) + j * (this->NUM_LEGS);
		cblas_dcopy(this->NUM_LEGS, initial_orientation, 1, x_vector + i * (this->DIM) + j * (this->NUM_LEGS), 1);
		j++;
		
		// pcc
		cblas_dcopy(this->DIM, pcc, 1, x_vector + i * (this->DIM) + j * (this->NUM_LEGS), 1);
		i++;
		// pw1
		T *pw1 = x_vector + i * (this->DIM) + j * (this->NUM_LEGS);
		i++;
		// pw2
		T *pw2 = x_vector + i * (this->DIM) + j * (this->NUM_LEGS);
		i++;
		// pw3
		T *pw3 = x_vector + i * (this->DIM) + j * (this->NUM_LEGS);
		i++;
		// pw4
		T *pw4 = x_vector + i * (this->DIM) + j * (this->NUM_LEGS);
		i++;
		// pt1
		T *pt1 = x_vector + i * (this->DIM) + j * (this->NUM_LEGS);
		i++;
		// pt2
		T *pt2 = x_vector + i * (this->DIM) + j * (this->NUM_LEGS);
		i++;
		// pt3
		T *pt3 = x_vector + i * (this->DIM) + j * (this->NUM_LEGS);
		i++;
		// pt4
		T *pt4 = x_vector + i * (this->DIM) + j * (this->NUM_LEGS);
		
		
		get_initial_length(start_orient, r1, r2, r3, r4, pcc, initial_upper_spring_length, initial_lower_spring_length, pw1, pw2, pw3, pw4, pt1, pt2, pt3, pt4);
		
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
			std::cout << "mass wheel = " << mass_wheel[0] << std::endl;
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

		MathLibrary::Solvers<T, MBD_method>::Broyden_CN(this, x_vector, complete_vector, this->h, this->num_iter, this->tol, this->max_iter);
		T* start = complete_vector + (this->num_iter)*this->solution_dim;
		cblas_dcopy(this->solution_dim, start, 1, solution_vector, 1);
		std::cout << "Solution copied!\n" << std::endl;
		mkl_free(complete_vector);
	}

	size_t get_alignment() {
		return this->alignment;
	}

	size_t get_solution_dimension() {
		return this->solution_dim;
	}

	~MBD_method() {

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


//#ifndef ComputeNasa
//#define ComputeNasa
//template <typename T>
//class ComputeNasa {
//private:
//	// define constants
//	const int dim = 3;
//	const int alignment = 64;
//	const int num_wheels = 4;
//	const int num_wheels_x_dim = num_wheels * dim;
//	const int DOF_diag = 9; // the diagonal elements from A
//	const int DOF = 11;
//	const int matrixElements = dim * dim;
//	// const int matrixElements = DOF * DOF;
//	T FC;
//	T* r, r_tilda, FW, FT, FR, lower_spring_length, upper_spring_length, lower_spring_stiffness, upper_spring_stiffness;
//	// The system matrix A. 0:3 - 2x2 upper-left symmetric block; 4:12 . We will store its inverse on diagonal blocks
//	// b vector, 11 DOF (from the 12 components we prohibit rotation along y-axis, b[2]=0 is ignored)
//	// x - the solution vector
//	T* A, rhs, x; 
//	T* Ic, Tc, basis_car, C_Nc, r_global, car_corners, ones4;
//	T* upper_length, lower_length, diff_vector, upper_force, lower_force, upper_F;
//	T* C_Nc_extended, sum_torque, w_tilda, Qc;
//
//public:
//	ComputeNasa(); // constructor without parameters (default, testing reasons); initialization
//	
//	ComputeNasa(T* r, T* r_tilda, T* FW, T* FT, T* FR, T* lower_spring_length, T* upper_spring_length, T* lower_spring_stiffness, T* upper_spring_stiffness, T* A, T FC, T* IC) :
//		r(r), r_tilda(r_tilda), FW(FW), FT(FT), FR(FR), lower_spring_length(lower_spring_length), upper_spring_length(upper_spring_length), lower_spring_stiffness(lower_spring_stiffness), upper_spring_stiffness(upper_spring_stiffness), A(A), FC(FC), IC(IC) {
//		// memset(ones4, 1, num_wheels);
//	}
//
//	~ComputeNasa();
//
//	void compute_f3D_reduced(T time, T* f);
//};
//#endif

