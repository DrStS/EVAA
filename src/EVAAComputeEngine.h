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
#include "MBD_method.h"


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

	EVAAComputeEngine(std::string xmlFileName, std::string loadxml);
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
	std::string _xmlLoadFileName;
	Simulation_Parameters _parameters;
	Load_Params _load_module_parameter;
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
	T* k_vect, *l_lat, *l_long;
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
		temp[0] = k_body_fl_ + k_body_fr_ + k_body_rl_ + k_body_rr_; // K[i*DOF + 0] to be set later
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
		temp[1] = l_lat_fl_ * l_lat_fl_ * k_body_fl_ + l_lat_fr_ * l_lat_fr_ * k_body_fr_ + l_lat_rl_ * l_lat_rl_ * k_body_rl_ + l_lat_rr_ * l_lat_rr_ * k_body_rr_; // K[i*DOF + 1] to be set later
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
		temp[2] = l_long_fl_ * l_long_fl_ * k_body_fl_ + l_long_fr_ * l_long_fr_ * k_body_fr_ + l_long_rl_ * l_long_rl_ * k_body_rl_ + l_long_rr_ * l_long_rr_ * k_body_rr_; // K[i*DOF + 2] to be set later
		K[i * DOF + 3] = l_long_fl_ * k_body_fl_;
		K[i * DOF + 4] = 0;
		K[i * DOF + 5] = l_long_fr_ * k_body_fr_;
		K[i * DOF + 6] = 0;
		K[i * DOF + 7] = -l_long_rl_ * k_body_rl_;
		K[i * DOF + 8] = 0;
		K[i * DOF + 9] = -l_long_rr_ * k_body_rr_;
		K[i * DOF + 10] = 0;

		i = 3;

		temp[3] = k_body_fl_ + k_tyre_fl_; // K[i*DOF + 3]
		K[i * DOF + 4] = -k_tyre_fl_;
		K[i * DOF + 5] = 0;
		K[i * DOF + 6] = 0;
		K[i * DOF + 7] = 0;
		K[i * DOF + 8] = 0;
		K[i * DOF + 9] = 0;
		K[i * DOF + 10] = 0;
		// all others are zero

		i = 4;
		temp[4] = k_tyre_fl_; //K[i*DOF + 4]
		K[i * DOF + 5] = 0;
		K[i * DOF + 6] = 0;
		K[i * DOF + 7] = 0;
		K[i * DOF + 8] = 0;
		K[i * DOF + 9] = 0;
		K[i * DOF + 10] = 0;

		i = 5;
		temp[5] = k_body_fr_ + k_tyre_fr_; // K[i*DOF + 5]
		K[i * DOF + 6] = -k_tyre_fr_;
		K[i * DOF + 7] = 0;
		K[i * DOF + 8] = 0;
		K[i * DOF + 9] = 0;
		K[i * DOF + 10] = 0;

		i = 6;
		temp[6] = k_tyre_fr_; // K[i*DOF + 6]
		K[i * DOF + 7] = 0;
		K[i * DOF + 8] = 0;
		K[i * DOF + 9] = 0;
		K[i * DOF + 10] = 0;

		i = 7;
		temp[7] = k_body_rl_ + k_tyre_rl_; // K[i*DOF + 7]
		K[i * DOF + 8] = -k_tyre_rl_;
		K[i * DOF + 9] = 0;
		K[i * DOF + 10] = 0;


		i = 8;
		temp[8] = k_tyre_rl_; // K[i*DOF + 8]
		K[i * DOF + 9] = 0;
		K[i * DOF + 10] = 0;

		i = 9;
		temp[9] = k_body_rr_ + k_tyre_rr_; // K[i*DOF + 9]
		K[i * DOF + 10] = -k_tyre_rr_;

		i = 10;
		temp[10] = k_tyre_rr_; // K[i*DOF + 10]

		// K=K+K'-diag(diag(K));
		cblas_dcopy(DOF * DOF, K, 1, K_trans, 1);
		mkl_dimatcopy('R', 'T', DOF, DOF, 1.0, K_trans, DOF, DOF); // get transpose of matrix
		cblas_daxpy(DOF * DOF, 1.0, K_trans, 1, K, 1); // K = K + K'
		MathLibrary::allocate_to_diagonal(K, temp, DOF); // K = K + K'+ diag(K)
	}
	void populate_K(T* k_vect, T k_body_fl, T k_tyre_fl, T k_body_fr,
		T k_tyre_fr,
		T k_body_rl,
		T k_tyre_rl,
		T k_body_rr,
		T k_tyre_rr) {
		k_vect[0] = k_body_fl;
		k_vect[1] = k_tyre_fl;
		k_vect[2] = k_body_fr;
		k_vect[3] = k_tyre_fr;
		k_vect[4] = k_body_rl;
		k_vect[5] = k_tyre_rl;
		k_vect[6] = k_body_rr;
		k_vect[7] = k_tyre_rr;
	}


public:
	linear11dof(const Simulation_Parameters &params, const Load_Params &load_param){

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
		k_vect = (T*)mkl_calloc(num_wheels*2, sizeof(T), alignment);
		l_lat = (T*)mkl_calloc(num_wheels, sizeof(T), alignment);
		l_long = (T*)mkl_calloc(num_wheels, sizeof(T), alignment);

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////// Initial Iteration vector ////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		u_n[0] = u_init_body;
		u_n[1] = theta_x_init_body;
		u_n[2] = theta_z_init_body;
		u_n[3] = params.initial_pos_wheel[2 * 3 + 1];		// TODO: calculate this from initial orientation
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

		f_n_p_1[0] = load_param.external_force_body[1];
		f_n_p_1[3] = load_param.external_force_wheel[2 * 3 + 1];
		f_n_p_1[4] = load_param.external_force_tyre[2 * 3 + 1];
		f_n_p_1[5] = load_param.external_force_wheel[3 * 3 + 1];
		f_n_p_1[6] = load_param.external_force_tyre[3 * 3 + 1];
		f_n_p_1[7] = load_param.external_force_wheel[1 * 3 + 1];
		f_n_p_1[8] = load_param.external_force_tyre[1 * 3 + 1];
		f_n_p_1[9] = load_param.external_force_wheel[0 * 3 + 1];
		f_n_p_1[10] = load_param.external_force_tyre[0 * 3 + 1];

		

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
		populate_K(k_vect, k_body_fl, k_tyre_fl, k_body_fr, k_tyre_fr, k_body_rl, k_tyre_rl, k_body_rr, k_tyre_rr);
		l_lat[0] = l_lat_fl;
		l_lat[1] = l_lat_fr;
		l_lat[2] = l_lat_rl;
		l_lat[3] = l_lat_rr;
		l_long[0] = l_long_fl;
		l_long[1] = l_long_fr;
		l_long[2] = l_long_rl;
		l_long[3] = l_long_rr;
		construct_K(K, k_vect, l_lat, l_long);
		//write_matrix(K, DOF);
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
	

		
		
		// B=((2/(h*h))*M+(1/h)*D);
		cblas_daxpy(mat_len, 2 * factor_h2, M, 1, B, 1);
		cblas_daxpy(mat_len, factor_h, D, 1, B, 1);
		int iter = 1;
		T t = h_;
		double eps = h_/100;
		/*auto start = std::chrono::steady_clock::now();*/
		while (std::abs(t-(tend_+h_)) > eps) {

			// K update here


			// A=((1/(h*h))*M+(1/h)*D+K);
			cblas_dcopy(mat_len, M, 1, A, 1);
			cblas_daxpy(mat_len, factor_h, D, 1, A, 1);
			cblas_daxpy(mat_len, 1, K, 1, A, 1);
			// Cholesky factorization of A
			lapack_int status;
			//lapack_int* piv = (lapack_int*)mkl_calloc((DOF + 4), sizeof(lapack_int), alignment);
			//write_matrix(A, DOF);
			status = LAPACKE_dpotrf(LAPACK_ROW_MAJOR, 'L', DOF, A, DOF);
			//status = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, DOF + 4, DOF + 4, A, DOF + 4, piv);
			check_status(status);
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
		//std::cout << "iter = " << iter << " sol_size = "<< sol_size <<"\n\n" << std::endl;
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
