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
			throw "Matrix non Positive Definite";
		}
		else if (status == -1) {
			throw "Matrix contain illegal value";
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
		////////////////////////////////////// Mass ////////////////////////////////////////////////////////////////////////
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
		try {
			//write_matrix(A, DOF);
			status = LAPACKE_dpotrf(LAPACK_ROW_MAJOR, 'L', DOF, A, DOF);
			//status = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, DOF + 4, DOF + 4, A, DOF + 4, piv);
			check_status(status);
			//write_matrix(A, DOF);
		}
		catch (const char* msg) {
			std::cerr << msg << std::endl;
		}
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
		try {
			status = LAPACKE_dpotrf(LAPACK_ROW_MAJOR, 'L', DOF, Ared, DOF);
			//status = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, DOF + 4, DOF + 4, A, DOF + 4, piv);
			check_status(status);
		}
		catch (const char* msg) {
			std::cerr << msg << std::endl;
		}
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
			MKL_free(u_sol);
			MKL_free(u_n_p_1);
			MKL_free(tyre_index_set);
			MKL_free(A); 
			MKL_free(B);
			MKL_free(time);
		}
		else if (s == reduced) {
			MKL_free(u_sol_red);
			MKL_free(u_n_red);
			MKL_free(u_n_m_1_red);
			MKL_free(u_n_p_1_red);
			MKL_free(f_n_p_1_red);
			MKL_free(Ared);
			MKL_free(Bred);
			MKL_free(time);
			MKL_free(K_red);
			MKL_free(M_red);
			MKL_free(D_red);
		}
	}

	~linear11dof() {
		MKL_free(M);
		MKL_free(temp);
		MKL_free(K);
		MKL_free(K_trans);
		MKL_free(D);
		MKL_free(u_n_m_1);
		MKL_free(u_n);
		MKL_free(quad_angle_init);
		MKL_free(euler_angle_init);
		MKL_free(f_n_p_1);
	}
};

template <class T>
class MBD_method
{
public:
	MBD_method(T* initial_orientation, 
		T* initial_upper_spring_length, 
		T* initial_lower_spring_length,
		T* vc,
		T* vw1,
		T* vw2,
		T* vw3, 
		T* vw4,
		T* vt1, 
		T* vt2,
		T* vt3,
		T* vt4,
		T* wc) {
		int i;
		r1 = (T*)mkl_calloc(DIM, sizeof(T), alignment);
		r2 = (T*)mkl_calloc(DIM, sizeof(T), alignment);
		r3 = (T*)mkl_calloc(DIM, sizeof(T), alignment);
		r4 = (T*)mkl_calloc(DIM, sizeof(T), alignment);
		Ic = (T*)mkl_calloc(DIM * DIM, sizeof(T), alignment);
		mass_wheel = (T*)mkl_calloc(NUM_LEGS, sizeof(T), alignment);
		mass_tyre = (T*)mkl_calloc(NUM_LEGS, sizeof(T), alignment);
		upper_spring_length = (T*)mkl_calloc(NUM_LEGS, sizeof(T), alignment);
		lower_spring_length = (T*)mkl_calloc(NUM_LEGS, sizeof(T), alignment);
		upper_spring_stiffness = (T*)mkl_calloc(NUM_LEGS, sizeof(T), alignment);
		lower_spring_stiffness = (T*)mkl_calloc(NUM_LEGS, sizeof(T), alignment);
		upper_rotational_stiffness = (T*)mkl_calloc(NUM_LEGS, sizeof(T), alignment);
		lower_rotational_stiffness = (T*)mkl_calloc(NUM_LEGS, sizeof(T), alignment);

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

	}
	void apply_boundary_condition() {
		// from line 119 in matlab
	}
	~MBD_method() {
		MKL_free(r1);
		MKL_free(r2);
		MKL_free(r3);
		MKL_free(r4);
		MKL_free(Ic);
		MKL_free(mass_wheel);
		MKL_free(mass_tyre);
		MKL_free(upper_spring_length);
		MKL_free(lower_spring_length);
		MKL_free(upper_spring_stiffness);
		MKL_free(lower_spring_stiffness);
		MKL_free(upper_rotational_stiffness);
		MKL_free(lower_rotational_stiffness);
		
	}

private:
	const int DIM = 3;
	const int NUM_LEGS = 4;
	T k_body_fl = 28e3*0.69;
	T k_tyre_fl = 260e3;
	T k_body_fr = 28e3*0.69;
	T k_tyre_fr = 260e3;
	T k_body_rl = 28e3*0.69;
	T k_tyre_rl = 260e3;
	T k_body_rr = 28e3*0.69;
	T k_tyre_rr = 260e3;
	T k_body_rot_fl = 1e4;
	T k_body_rot_fr = 1e4;
	T k_body_rot_rl = 1e4;
	T k_body_rot_rr = 1e4;
	T k_tyre_rot_fl = 1e4;
	T k_tyre_rot_fr = 1e4;
	T k_tyre_rot_rl = 1e4;
	T k_tyre_rot_rr = 1e4;
	T l_long_fl = 1.395;
	T l_long_fr = 1.395;
	T l_long_rl = 1.596;
	T l_long_rr = 1.596;
	T l_lat_fl = 2 * 0.84;
	T l_lat_fr = 2 * 0.84;
	T l_lat_rl = 2 * 0.84;
	T l_lat_rr = 2 * 0.84;
	T mass = 1936;
	T I_body_xx = 640;
	T I_body_yy = 4800;
	T I_body_zz = 10000;
	T mass_wheel_fl = 145 / 2;
	T mass_tyre_fl = 30;
	T mass_wheel_fr = 145 / 2;
	T mass_tyre_fr = 30;
	T mass_wheel_rl = 135 / 2;
	T mass_tyre_rl = 30;
	T mass_wheel_rr = 135 / 2;
	T mass_tyre_rr = 30;
	T upper_spring_length_rr = 0.2;
	T upper_spring_length_rl = 0.2;
	T upper_spring_length_fl = 0.2;
	T upper_spring_length_fr = 0.2;
	T lower_spring_length_rr = 0.2;
	T lower_spring_length_rl = 0.2;
	T lower_spring_length_fl = 0.2;
	T lower_spring_length_fr = 0.2;
	T g = 0;
	T FC = -mass * g;
	T FT1, FT2, FT3, FT4;
	T FW1, FW2, FW3, FW4;
	T h = 1e-3;
	T num_iter;
	int max_iter;
	T tol = 1e-10;
	T *r1, *r2, *r3, *r4;
	T *Ic;
	T *mass_wheel;
	T *upper_spring_length, *lower_spring_length;
	T *upper_spring_stiffness, *lower_spring_stiffness;
	T *upper_rotational_stiffness, *lower_rotational_stiffness;
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

