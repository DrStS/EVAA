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
	int DOF = 11;

	T k_body_fl = 28e3 * 0.69;
	T k_tyre_fl = 260e3;
	T k_body_fr = 28e3 * 0.69;
	T k_tyre_fr = 260e3;
	T k_body_rl = 16e3 * 0.82;
	T k_tyre_rl = 260e3;
	T k_body_rr = 16e3 * 0.82;
	T k_tyre_rr = 260e3;
	T l_long_fl = 1.395;
	T l_long_fr = 1.395;
	T l_long_rl = 1.596;
	T l_long_rr = 1.596;
	T l_lat_fl = 2 * 0.8458;
	T l_lat_fr = 2 * 0.8458;
	T l_lat_rl = 2 * 0.84;
	T l_lat_rr = 2 * 0.84;
	T mass_Body = 1936.0;
	T I_body_xx = 640.0;
	T I_body_yy = 4800.0;
	T mass_wheel_fl = 145.0 / 2.0;
	T mass_tyre_fl = 0.0;
	T mass_wheel_fr = 145.0 / 2.0;
	T mass_tyre_fr = 0.0;
	T mass_wheel_rl = 135.0 / 2.0;
	T mass_tyre_rl = 0.0;
	T mass_wheel_rr = 135.0 / 2.0;
	T mass_tyre_rr = 0.0;
	////////////////////////////////////////////////////////////////////////////////////////////
	T tend_;
	T u_init_;
	T du_init_;
	T h_;
	int i;
	T* M, * temp, * K, * K_trans, * D, * M_red, * D_red, * Kred;
	T* u_sol, * u_sol_red, * u_n_p_1, * u_n_p_1_red, * u_n_m_1, * u_n, * u_n_red, * u_n_m_1_red, * A, * Ared, * B, * Bred, * f_n_p_1, * f_n_p_1_red;
	T* time; // this is not necessary

	void check_status(lapack_int status) {
		if (status == 1) {
			throw "Matrix non Positive Definite";
		}
		else if (status == -1) {
			throw "Matrix contain illegal value";
		}
	}



public:
	linear11dof(T tend, T h, T u_init, T du_init) :tend_(tend), h_(h), u_init_(u_init), du_init_(du_init) {

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////// System Mono ///////////////////////////////////////////////////////
		///////////////////////////////// Memory Allocation and matrix formulation ///////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		M = (T*)mkl_calloc(DOF * DOF, sizeof(T), alignment);
		temp = (T*)mkl_malloc(DOF * sizeof(T), alignment);
		K = (T*)mkl_calloc(DOF * DOF, sizeof(T), alignment);
		K_trans = (T*)mkl_calloc(DOF * DOF, sizeof(T), alignment);
		D = (T*)mkl_calloc(DOF * DOF, sizeof(T), alignment);

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////// Mass ////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		temp[0] = mass_Body;
		temp[1] = I_body_xx;
		temp[2] = I_body_yy;
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

	void apply_boundary_condition(std::string s, T* force) {
		std::string cond1 = "fixed_to_road";
		std::string cond2 = "road_force";

		if (s == cond1) {
			T* start_loc_k_trans, * start_loc_k, * start_loc_D_red, * start_loc_D;
			// using K_trans as temperory variable copy all element of K to k_trans
			// cblas_dcopy(DOF*DOF, K, 1, K_trans, 1);
			// change the DOF for further reference
			DOF = DOF - 4; // DOF 7 now
			Kred = (T*)mkl_calloc(DOF * DOF, sizeof(T), alignment);
			M_red = (T*)mkl_calloc(DOF * DOF, sizeof(T), alignment);
			D_red = (T*)mkl_calloc(DOF * DOF, sizeof(T), alignment);
			// memory allocation for the reduced force field
			f_n_p_1 = (T*)mkl_calloc(DOF + 4, sizeof(T), alignment);
			f_n_p_1_red = (T*)mkl_calloc(DOF, sizeof(T), alignment);
			if (*force != NULL) {
				cblas_dcopy(DOF + 4, force, 1, f_n_p_1, 1);
				cblas_dcopy(DOF - 3, f_n_p_1, 1, f_n_p_1_red, 1);
			}

			// reassign the values of K in the new K for first 4 DOF 
			for (int i = 0; i < DOF - 3; ++i) {
				start_loc_k_trans = K + i * (DOF + 4);
				start_loc_k = Kred + i * DOF;
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
			}
			int j = 4;
			for (int i = 5; i < (DOF + 4); i = i + 2) {
				start_loc_k_trans = K + (i) * (DOF + 4);
				start_loc_k = Kred + j * DOF;
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
				j++;
			}
		}

	}

	void solve(T* sol_vect) {
		///////////////////////////////////// For 7 DOF System //////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////// Memory Allocation for Intermediate Solution Vectors //////////////////
		////////////////////////////////////////////////////////////////////////////////////////////////////////////
		int sol_size = (floor(tend_ / h_) + 1);
		T factor_h2 = (1 / (h_ * h_));
		T factor_h = (1 / (h_));
		int mat_len = (DOF + 4) * (DOF + 4);
		u_sol = (T*)mkl_calloc(sol_size * (DOF + 4), sizeof(T), alignment);
		u_sol_red = (T*)mkl_calloc(sol_size * (DOF), sizeof(T), alignment);
		u_n_p_1 = (T*)mkl_calloc((DOF + 4), sizeof(T), alignment);
		u_n_p_1_red = (T*)mkl_calloc((DOF), sizeof(T), alignment);
		u_n_m_1 = (T*)mkl_calloc((DOF + 4), sizeof(T), alignment);
		u_n = (T*)mkl_calloc((DOF + 4), sizeof(T), alignment);
		u_n_red = (T*)mkl_calloc((DOF), sizeof(T), alignment);
		u_n_m_1_red = (T*)mkl_calloc((DOF), sizeof(T), alignment);;
		A = (T*)mkl_calloc(mat_len, sizeof(T), alignment);
		Ared = (T*)mkl_calloc((DOF) * (DOF), sizeof(T), alignment);;
		B = (T*)mkl_calloc(mat_len, sizeof(T), alignment);
		Bred = (T*)mkl_calloc((DOF) * (DOF), sizeof(T), alignment);
		time = (T*)mkl_calloc(sol_size, sizeof(T), alignment);
		u_n[0] = u_init_;
		cblas_dcopy(DOF + 4, u_n, 1, u_n_m_1, 1);
		u_n_m_1[0] = u_init_ - h_ * du_init_;

		// A=((1/(h*h))*M+(1/h)*D+K);
		cblas_daxpy(mat_len, factor_h2, M, 1, A, 1);
		cblas_daxpy(mat_len, factor_h, D, 1, A, 1);
		cblas_daxpy(mat_len, 1, K, 1, A, 1);
		// Cholesky factorization of A
		lapack_int status;
		try {
			status = LAPACKE_dpotrf(LAPACK_ROW_MAJOR, 'L', DOF + 4, A, DOF + 4);
			check_status(status);
		}
		catch (const char* msg) {
			std::cerr << msg << std::endl;
		}
		// B=((2/(h*h))*M+(1/h)*D);
		cblas_daxpy(mat_len, 2 * factor_h2, M, 1, B, 1);
		cblas_daxpy(mat_len, factor_h, D, 1, B, 1);
		int iter = 1;
		T t = h_;
		/*auto start = std::chrono::steady_clock::now();*/
		while (t < tend_) {
			// u_n_p_1=A\(B*u_n-((1/(h*h))*M)*u_n_m_1+f_n_p_1);
			// u_n_p_1 = B*u_n
			cblas_dgemv(CblasRowMajor, CblasNoTrans, DOF + 4, DOF + 4, 1, B, DOF + 4, u_n, 1, 1, u_n_p_1, 1);
			// u_n_p_1 = -((1/(h*h))*M)*u_n_m_1 + u_n_p_1
			cblas_dgemv(CblasRowMajor, CblasNoTrans, DOF + 4, DOF + 4, -factor_h2, M, DOF + 4, u_n_m_1, 1, 1, u_n_p_1, 1);
			// u_n_p_1 = f_n_p_1 + u_n_p_1
			cblas_daxpy(DOF + 4, 1, f_n_p_1, 1, u_n_p_1, 1);
			// u_n_p_1=A\u_n_p_1
			LAPACKE_dpotrs(LAPACK_ROW_MAJOR, 'L', DOF + 4, 1, A, DOF + 4, u_n_p_1, 1);
			time[iter] = t;
			// u_sol(j,:)=u_n_p_1;
			cblas_dcopy(DOF + 4, u_n_p_1, 1, u_sol + iter * (DOF + 4), 1);
			// u_n_m_1=u_n;
			cblas_dcopy(DOF + 4, u_n, 1, u_n_m_1, 1);
			// u_n    =u_n_p_1;
			cblas_dcopy(DOF + 4, u_n_p_1, 1, u_n, 1);
			iter++;
			t += h_;
		}
		/*auto end = std::chrono::steady_clock::now();
		std::cout << "Elapsed time in nanoseconds : "
			<< std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
			<< " ms" << std::endl;*/
		//std::cout << "iter = " << iter << " sol_size = "<< sol_size <<"\n\n" << std::endl;
		cblas_dcopy(DOF + 4, u_sol + (iter - 1)*(DOF + 4), 1, sol_vect, 1);

	}

	~linear11dof() {
		MKL_free(M);
		MKL_free(temp);
		MKL_free(K);
		MKL_free(K_trans);
		MKL_free(D);
		MKL_free(Kred);
		MKL_free(M_red);
		MKL_free(D_red);
		// memory allocation for the reduced force field
		MKL_free(f_n_p_1);
		MKL_free(f_n_p_1_red);
		MKL_free(u_sol);
		MKL_free(u_sol_red);
		MKL_free(u_n_p_1);
		MKL_free(u_n_p_1_red);
		MKL_free(u_n_m_1);
		MKL_free(u_n);
		MKL_free(u_n_red);
		MKL_free(u_n_m_1_red);
		MKL_free(A);
		MKL_free(Ared);
		MKL_free(B);
		MKL_free(Bred);
		MKL_free(time);
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