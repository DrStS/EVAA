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
#include <iostream>
#include <vector>
#include "EVAAComputeEngine.h"
#include "MathLibrary.h"
#include <limits>
#ifdef USE_INTEL_MKL
#include <mkl.h>
#define USE_GEMM
#endif

#ifdef USE_EIGEN
#include <Eigen/Dense>
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Matrix3d;
using Eigen::IOFormat;
#endif

#ifdef USE_BLAZE
#include <blaze/Math.h>
#endif

typedef std::numeric_limits< double > dbl;
typedef double floatEVAA;

EVAAComputeEngine::EVAAComputeEngine(std::string xmlFileName) {
	// Intialize XML metadatabase singelton
	_xmlFileName = xmlFileName;
	std::cout << "Read XML file: " << _xmlFileName << std::endl;
	ReadXML reader(_xmlFileName);
	reader.ReadParameters(_parameters);
}


EVAAComputeEngine::~EVAAComputeEngine() {
}

void EVAAComputeEngine::prepare(void) {
	MathLibrary::printMKLInfo();
}

void EVAAComputeEngine::computeEigen11DOF(void) {

	ReadXML reader(_xmlFileName);

	int DOF = _parameters.DOF;

	// K stiffness
	double k_11 = _parameters.k_body[0];
	double k_12 = _parameters.k_tyre[0];
	double k_21 = _parameters.k_body[1];
	double k_22 = _parameters.k_tyre[1];
	double k_31 = _parameters.k_body[2];
	double k_32 = _parameters.k_tyre[2];
	double k_41 = _parameters.k_body[3];
	double k_42 = _parameters.k_tyre[3];
	double l_1 = _parameters.l_lat[0];
	double l_2 = _parameters.l_lat[2];
	double l_3 = _parameters.l_long[0];
	double l_4 = _parameters.l_long[2];


	// D - damping matrix
	double d_11 = _parameters.c_body[0];
	double d_12 = _parameters.c_tyre[0];
	double d_21 = _parameters.c_body[1];
	double d_22 = _parameters.c_tyre[1];
	double d_31 = _parameters.c_body[2];
	double d_32 = _parameters.c_tyre[2];
	double d_41 = _parameters.c_body[3];
	double d_42 = _parameters.c_tyre[3];
	double d_l1 = 0.1*_parameters.l_lat[0];
	double d_l2 = 0.1*_parameters.l_lat[2];
	double d_l3 = 0.1*_parameters.l_long[0];
	double d_l4 = 0.1*_parameters.l_long[2];

	// M - mass matrix
	double m_1 = _parameters.mass_body;
	double m_2 = _parameters.I_body[0];
	double m_3 = _parameters.I_body[2];
	double m_4 = _parameters.mass_wheel[2];
	double m_5 = _parameters.mass_tyre[2];
	double m_6 = _parameters.mass_wheel[3];
	double m_7 = _parameters.mass_tyre[3];
	double m_8 = _parameters.mass_wheel[1];
	double m_9 = _parameters.mass_tyre[1];
	double m_10 = _parameters.mass_wheel[0];
	double m_11 = _parameters.mass_tyre[0];

	MatrixXd A(DOF, DOF);
	MatrixXd B(DOF, DOF);
	MatrixXd M(DOF, DOF);
	MatrixXd D(DOF, DOF);
	MatrixXd K(DOF, DOF);
	MatrixXd K_aux(DOF, DOF);
	MatrixXd D_aux(DOF, DOF);

	VectorXd u_n_p_1(DOF);
	VectorXd u_n(DOF);
	VectorXd u_n_m_1(DOF);

	// Define the mass matrix
	M = MatrixXd::Zero(DOF, DOF);
	M(0, 0) = m_1;
	M(1, 1) = m_2;
	M(2, 2) = m_3;
	M(3, 3) = m_4;
	M(4, 4) = m_5;
	M(5, 5) = m_6;
	M(6, 6) = m_7;
	M(7, 7) = m_8;
	M(8, 8) = m_9;
	M(9, 9) = m_10;
	M(10, 10) = m_11;

	// Define the stiffness matrix
	K << k_11 + k_21 + k_31 + k_41, -k_11 * l_1 - k_41 * l_1 + k_21 * l_2 + k_31 * l_2, -k_11 * l_3 - k_21 * l_3 + k_31 * l_4 + k_41 * l_4, -k_11, 0.0, -k_21, 0.0, -k_31, 0.0, -k_41, 0.0,
		0.0, l_1* l_1* k_11 + l_2 * l_2 * k_21 + l_2 * l_2 * k_31 + l_1 * l_1 * k_41, l_1* l_3* k_11 - l_3 * l_2 * k_21 + l_2 * l_4 * k_31 - l_1 * l_4 * k_41, l_1* k_11, 0.0, -l_2 * k_21, 0.0, -l_2 * k_31, 0.0, l_1* k_41, 0.0,
		0.0, 0.0, l_3* l_3* k_11 + l_3 * l_3 * k_21 + l_4 * l_4 * k_31 + l_4 * l_4 * k_41, l_3* k_11, 0.0, l_3* k_21, 0.0, -l_4 * k_31, 0.0, -l_4 * k_41, 0.0,
		0.0, 0.0, 0.0, k_11 + k_12, -k_12, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
		0.0, 0.0, 0.0, 0.0, k_12, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
		0.0, 0.0, 0.0, 0.0, 0.0, k_21 + k_22, -k_22, 0.0, 0.0, 0.0, 0.0,
		0.0, 0.0, 0.0, 0.0, 0.0, 0.0, k_22, 0.0, 0.0, 0.0, 0.0,
		0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, k_31 + k_32, -k_32, 0.0, 0.0,
		0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, k_32, 0.0, 0.0,
		0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, k_41 + k_42, -k_42,
		0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, k_42;
	MatrixXd K_diag(DOF, DOF);
	K_aux = K + K.transpose();
	K_aux.diagonal() -= K.diagonal();
	K = K_aux;

	// Define the damping matrix
	D << d_11 + d_21 + d_31 + d_41, -d_11 * l_1 - d_41 * l_1 + d_21 * l_2 + d_31 * l_2, -d_11 * l_3 - d_21 * l_3 + d_31 * l_4 + d_41 * l_4, -d_11, 0.0, -d_21, 0.0, -d_31, 0.0, -d_41, 0.0,
		0.0, l_1* l_1* d_11 + l_2 * l_2 * d_21 + l_2 * l_2 * d_31 + l_1 * l_1 * d_41, l_1* l_3* d_11 - l_3 * l_2 * d_21 + l_2 * l_4 * d_31 - l_1 * l_4 * d_41, l_1* d_11, 0.0, -l_2 * d_21, 0.0, -l_2 * d_31, 0.0, l_1* d_41, 0.0,
		0.0, 0.0, l_3* l_3* d_11 + l_3 * l_3 * d_21 + l_4 * l_4 * d_31 + l_4 * l_4 * d_41, l_3* d_11, 0.0, l_3* d_21, 0.0, -l_4 * d_31, 0.0, -l_4 * d_41, 0.0,
		0.0, 0.0, 0.0, d_11 + d_12, -d_12, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
		0.0, 0.0, 0.0, 0.0, d_12, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
		0.0, 0.0, 0.0, 0.0, 0.0, d_21 + d_22, -d_22, 0.0, 0.0, 0.0, 0.0,
		0.0, 0.0, 0.0, 0.0, 0.0, 0.0, d_22, 0.0, 0.0, 0.0, 0.0,
		0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, d_31 + d_32, -d_32, 0.0, 0.0,
		0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, d_32, 0.0, 0.0,
		0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, d_41 + d_42, -d_42,
		0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, d_42;
	MatrixXd D_diag(DOF, DOF);
	D_aux = D + D.transpose();
	D_aux.diagonal() -= D.diagonal();
	D = D_aux;
	A = MatrixXd::Zero(DOF, DOF);
	B = MatrixXd::Zero(DOF, DOF);

	u_n_p_1 = VectorXd::Zero(DOF);
	u_n = VectorXd::Zero(DOF);
	u_n(0) = 1;
	u_n_m_1 = u_n;

	//	int nRefinement = 10;
	//	int numTimeSteps = pow(2, nRefinement);
	int numTimeSteps = _parameters.num_time_iter;
	//time step size 
	//double h =  / (numTimeSteps);
	double h = _parameters.timestep;
	std::cout << "Time step h is: " << h << std::scientific << std::endl;

	/*IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
	std::cout << "M: \n\n" << M.format(CleanFmt) << std::endl << std::endl;
	std::cout << "D: \n\n" << D.format(CleanFmt) << std::endl << std::endl;
	std::cout << "K: \n\n" << K.format(CleanFmt) << std::endl << std::endl;*/

	/// Build dynamic stiffness matrix
	// A = (1.0/(h*h))*M + (1.0/h)*D + K
	A = (1.0 / (h * h)) * M + 1.0 / h * D + K;

	///Build rhs for BE integrator
	//B = (2.0 / (h*h))*M + 1.0 / h * D + B
	B = 2.0 / (h * h) * M + 1.0 / h * D;

	/*std::cout << "A: \n\n" << A.format(CleanFmt) << std::endl << std::endl;
	std::cout << "B: \n\n" << B.format(CleanFmt) << std::endl << std::endl;*/
	// LU Decomposition
	Eigen::PartialPivLU<MatrixXd> lu(A);

	M *= (-1.0 / (h * h));
	for (int iTime = 0; iTime < numTimeSteps; iTime++) {
		// B*u_n + (-1.0 / (h*h))*M
		// Solve system
		u_n_p_1 = lu.solve(B * u_n + M * u_n_m_1);

		u_n_m_1 = u_n;
		u_n = u_n_p_1;
		reader.ReadVariableParameters(_parameters);
	}
	std::cout << "We ran #" << numTimeSteps << " time steps!" << std::endl;
	std::cout << u_n_p_1 << std::scientific << std::endl;
}

void EVAAComputeEngine::computeMKL11DOF(void) {

	ReadXML reader(_xmlFileName);

	int DOF = _parameters.DOF;

	// K stiffness
	double k_11 = _parameters.k_body[0];
	double k_12 = _parameters.k_tyre[0];
	double k_21 = _parameters.k_body[1];
	double k_22 = _parameters.k_tyre[1];
	double k_31 = _parameters.k_body[2];
	double k_32 = _parameters.k_tyre[2];
	double k_41 = _parameters.k_body[3];
	double k_42 = _parameters.k_tyre[3];
	double l_1 = _parameters.l_lat[0];
	double l_2 = _parameters.l_lat[2];
	double l_3 = _parameters.l_long[0];
	double l_4 = _parameters.l_long[2];


	// D - damping matrix
	double d_11 = _parameters.c_body[0];
	double d_12 = _parameters.c_tyre[0];
	double d_21 = _parameters.c_body[1];
	double d_22 = _parameters.c_tyre[1];
	double d_31 = _parameters.c_body[2];
	double d_32 = _parameters.c_tyre[2];
	double d_41 = _parameters.c_body[3];
	double d_42 = _parameters.c_tyre[3];
	double d_l1 = 0.1 * _parameters.l_lat[0];
	double d_l2 = 0.1 * _parameters.l_lat[2];
	double d_l3 = 0.1 * _parameters.l_long[0];
	double d_l4 = 0.1 * _parameters.l_long[2];

	// M - mass matrix
	double m_1 = _parameters.mass_body;
	double m_2 = _parameters.I_body[0];
	double m_3 = _parameters.I_body[2];
	double m_4 = _parameters.mass_wheel[2];
	double m_5 = _parameters.mass_tyre[2];
	double m_6 = _parameters.mass_wheel[3];
	double m_7 = _parameters.mass_tyre[3];
	double m_8 = _parameters.mass_wheel[1];
	double m_9 = _parameters.mass_tyre[1];
	double m_10 = _parameters.mass_wheel[0];
	double m_11 = _parameters.mass_tyre[0];

	int alignment = 64;
	int matrixElements = DOF * DOF;

	// allocate matrices of zeros
	double* B = (double*)mkl_calloc(matrixElements, sizeof(double), alignment);
	double* M = (double*)mkl_calloc(matrixElements, sizeof(double), alignment);
	double* D = (double*)mkl_calloc(matrixElements, sizeof(double), alignment);
	double* K = (double*)mkl_calloc(matrixElements, sizeof(double), alignment);

	//std::vector<double> B(121);
	//std::vector<double> M(121);
	//std::vector<double> D(121);
	//std::vector<double> K(121);

	double* u_n_p_1 = (double*)mkl_calloc(DOF, sizeof(double), alignment);
	double* u_n = (double*)mkl_calloc(DOF, sizeof(double), alignment);
	double* u_n_m_1 = (double*)mkl_calloc(DOF, sizeof(double), alignment);
	double* tmp = (double*)mkl_calloc(DOF, sizeof(double), alignment);


	//std::vector<double> f_n_p_1(11);
	//std::vector<double> u_n_p_1(11);
	//std::vector<double> u_n(11);
	//std::vector<double> u_n_m_1(11);
	//std::vector<double> tmp(11);

	// Mass matrix initialization
	M[0] = m_1;
	M[12] = m_2;
	M[24] = m_3;
	M[36] = m_4;
	M[48] = m_5;
	M[60] = m_6;
	M[72] = m_7;
	M[84] = m_8;
	M[96] = m_9;
	M[108] = m_10;
	M[120] = m_11;

	/*std::cout << "\nM:\n";
	for (int i = 0; i < 11; ++i) {
		for (int j = 0; j < 11; ++j) {
			std::cout << M[i * 11 + j] << "\t";
		}
		std::cout << std::endl;
	}*/

	// Stiffness matrix
	K[0] = k_11 + k_21 + k_31 + k_41;
	K[1] = -k_11 * l_1 - k_41 * l_1 + k_21 * l_2 + k_31 * l_2;
	K[2] = -k_11 * l_3 - k_21 * l_3 + k_31 * l_4 + k_41 * l_4;
	K[3] = -k_11;
	K[5] = -k_21;
	K[7] = -k_31;
	K[9] = -k_41;
	K[12] = l_1 * l_1 * k_11 + l_2 * l_2 * k_21 + l_2 * l_2 * k_31 + l_1 * l_1 * k_41;
	K[13] = l_1 * l_3 * k_11 - l_3 * l_2 * k_21 + l_2 * l_4 * k_31 - l_1 * l_4 * k_41;
	K[14] = l_1 * k_11;
	K[16] = -l_2 * k_21;
	K[18] = -l_2 * k_31;
	K[20] = l_1 * k_41;
	K[24] = l_3 * l_3 * k_11 + l_3 * l_3 * k_21 + l_4 * l_4 * k_31 + l_4 * l_4 * k_41;
	K[25] = l_3 * k_11;
	K[27] = l_3 * k_21;
	K[29] = -l_4 * k_31;
	K[31] = -l_4 * k_41;
	K[36] = k_11 + k_12;
	K[37] = -k_12;
	K[48] = k_12;
	K[60] = k_21 + k_22;
	K[61] = -k_22;
	K[72] = k_22;
	K[84] = k_31 + k_32;
	K[85] = -k_32;
	K[96] = k_32;
	K[108] = k_41 + k_42;
	K[109] = -k_42;
	K[120] = k_42;

	// Symmetry
	for (int i = 1; i < DOF; ++i)
		for (int j = 0; j < i; ++j)
			K[i * DOF + j] = K[j * DOF + i];

	/*std::cout << "\nK:\n";
	for (int i = 0; i < 11; ++i) {
		for (int j = 0; j < 11; ++j) {
			std::cout << K[i * 11 + j] << "\t";
		}
		std::cout << std::endl;
	}*/

	// Damping matrix
	D[0] = d_11 + d_21 + d_31 + d_41;
	D[1] = -d_11 * l_1 - d_41 * l_1 + d_21 * l_2 + d_31 * l_2;
	D[2] = -d_11 * l_3 - d_21 * l_3 + d_31 * l_4 + d_41 * l_4;
	D[3] = -d_11;
	D[5] = -d_21;
	D[7] = -d_31;
	D[9] = -d_41;
	D[12] = l_1 * l_1 * d_11 + l_2 * l_2 * d_21 + l_2 * l_2 * d_31 + l_1 * l_1 * d_41;
	D[13] = l_1 * l_3 * d_11 - l_3 * l_2 * d_21 + l_2 * l_4 * d_31 - l_1 * l_4 * d_41;
	D[14] = l_1 * d_11;
	D[16] = -l_2 * d_21;
	D[18] = -l_2 * d_31;
	D[20] = l_1 * d_41;
	D[24] = l_3 * l_3 * d_11 + l_3 * l_3 * d_21 + l_4 * l_4 * d_31 + l_4 * l_4 * d_41;
	D[25] = l_3 * d_11;
	D[27] = l_3 * d_21;
	D[29] = -l_4 * d_31;
	D[31] = -l_4 * d_41;
	D[36] = d_11 + d_12;
	D[37] = -d_12;
	D[48] = d_12;
	D[60] = d_21 + d_22;
	D[61] = -d_22;
	D[72] = d_22;
	D[84] = d_31 + d_32;
	D[85] = -d_32;
	D[96] = d_32;
	D[108] = d_41 + d_42;
	D[109] = -d_42;
	D[120] = d_42;

	// Symmetry
	for (int i = 1; i < DOF; ++i)
		for (int j = 0; j < i; ++j)
			D[i * DOF + j] = D[j * DOF + i];

	/*std::cout << "\nD:\n";
	for (int i = 0; i < 11; ++i) {
		for (int j = 0; j < 11; ++j) {
			std::cout << D[i * 11 + j] << "\t";
		}
		std::cout << std::endl;
	}*/

	//Initial conditions
	u_n[0] = 1.;
	u_n_m_1[0] = 1.;
	/*std::cout << "u_n:\n";
	for (int i = 0; i < 11; ++i) {
		std::cout << u_n[i] << ", ";
	}*/


	mkl_set_num_threads(8);
	//	int nRefinement = 10;
	//	int numTimeSteps = pow(2, nRefinement);
	int numTimeSteps = _parameters.num_time_iter;
	//time step size 
	//double h =  / (numTimeSteps);
	double h = _parameters.timestep;
	std::cout << "Time step h is: " << h << std::scientific << std::endl;
	/// Build dynamic stiffness matrix
	// gets written in rear vector
	// K' <- (1.0/(h*h))*M + K
//	MathLibrary::computeDenseVectorAddition(M.data(), K.data(), (1.0 / (h*h)), 9);
	cblas_daxpy(121, (1.0 / (h * h)), M, 1, K, 1);
	// K <- (1.0/h)*D + K'
//	MathLibrary::computeDenseVectorAddition(D.data(), K.data(), (1.0 / h), 121);
	cblas_daxpy(121, (1.0 / h), D, 1, K, 1);
	/// K holds now dynamic stiffness matrix  for BE integrator
	///Build rhs for BE integrator
	//B' <-(2.0 / (h*h))*M + B
//	MathLibrary::computeDenseVectorAddition(M.data(), B.data(), (2.0 / (h*h)), 121);
	cblas_daxpy(matrixElements, (2.0 / (h * h)), M, 1, B, 1);

	//B <-(1.0 / (h))*D + B'
//	MathLibrary::computeDenseVectorAddition(D.data(), B.data(), (1.0 / h), 121);
	cblas_daxpy(matrixElements, (1.0 / h), D, 1, B, 1);
	//A*u_n_p_1=B*u_n+C*u_n_m_1+f_n_p_1 <== BE	

	std::vector<int> pivot(DOF);
	// LU Decomposition
//	MathLibrary::computeDenseSymLUFactorisation(11, K, pivot);
	LAPACKE_dgetrf(LAPACK_ROW_MAJOR, DOF, DOF, K, DOF, pivot.data());

	// Time loop
	double tmpScalar = (-1.0 / (h * h));
	for (int i = 0; i < DOF; ++i)
		M[(DOF + 1) * i] *= tmpScalar;

	void* jitter;
	// adds matrix matrix product to scalar matrix product
	std::cout << mkl_jit_create_dgemm(&jitter, MKL_COL_MAJOR, MKL_NOTRANS, MKL_NOTRANS, DOF, 1, DOF, 1.0, DOF, DOF, 0.0, DOF) << std::endl;
	dgemm_jit_kernel_t myDGEMMKernel = mkl_jit_get_dgemm_ptr(jitter);

	for (int iTime = 0; iTime < numTimeSteps; iTime++) {
		//timeVec[iTime] = iTime * h;
		// y: = alpha * A*x + beta * y
		// u_n_p_1 = B*u_n
#ifndef USE_GEMM
/* void cblas_dgemv(const CBLAS_LAYOUT Layout, const CBLAS_TRANSPOSE trans, const MKL_INT m,
	const MKL_INT n, const double alpha, const double* a, const MKL_INT lda, const double* x,
	const MKL_INT incx, const double beta, double* y, const MKL_INT incy); */
		cblas_dgemv(CblasColMajor, CblasNoTrans, 11, 11, 1.0, B, 11, u_n, 1, 0.0, u_n_p_1, 1);
#endif
#ifdef USE_GEMM 
		//cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 11, 1, 11, 1.0, B, 11, u_n, 11, 0.0, u_n_p_1, 11);
		myDGEMMKernel(jitter, B, u_n, u_n_p_1);
#endif
		// y: = alpha*A*x + beta*y
		// tmp = ((-1.0 / (h*h))*M) * u_n_m_1
#ifndef USE_GEMM
		cblas_dgemv(CblasColMajor, CblasNoTrans, DOF, DOF, (-1.0 / (h * h)), M, DOF, u_n_m_1, 1, 0.0, tmp, 1);
#endif
#ifdef USE_GEMM 
		//cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 11, 1, 11, 1.0, M, 11, u_n_m_1, 11, 0.0, tmp, 11);
		myDGEMMKernel(jitter, M, u_n_m_1, tmp);
#endif
		// u_n_p_1 <- 1.0 tmp + u_n_p_1
//		MathLibrary::computeDenseVectorAddition(tmp.data(), u_n_p_1.data(), 1.0, 11);
		//void cblas_daxpy (const MKL_INT n, const double a, const double *x, const MKL_INT incx, double *y, const MKL_INT incy);
		cblas_daxpy(DOF, 1.0, tmp, 1, u_n_p_1, 1);
		// Solve system
//		MathLibrary::computeDenseSymSolution(11, K, pivot, u_n_p_1);
// lapack_int LAPACKE_dgetrs (int matrix_layout , char trans , lapack_int n , lapack_int nrhs , const double * a , lapack_int lda , const lapack_int * ipiv , double * b , lapack_int ldb );
		LAPACKE_dgetrs(LAPACK_ROW_MAJOR, 'N', DOF, 1, K, DOF, pivot.data(), u_n_p_1, 1);

		for (int i = 0; i < DOF; ++i) {
			u_n_m_1[i] = u_n[i];
			u_n[i] = u_n_p_1[i];
		}
		reader.ReadVariableParameters(_parameters);
	}
	std::cout << "We ran #" << numTimeSteps << " time steps!" << std::endl;
	for (int i = 0; i < DOF; ++i) {
		std::cout << u_n_p_1[i] << std::scientific << std::endl;
	}

	// free the memory
	mkl_free(B);
	mkl_free(M);
	mkl_free(D);
	mkl_free(K);

	mkl_free(u_n_p_1);
	mkl_free(u_n);
	mkl_free(u_n_m_1);
	mkl_free(tmp);
}

/*
template <typename T>
class ComputeNasa {
public:
	// define constants
	const int alignment = 64;
	const int dim = 3;
	const int num_wheels = 4;
	const int dim_x_dim = dim * dim;
	const int num_wheels_x_dim = num_wheels * dim;
	const int DOF_diag = 9; // the diagonal elements from A
	const int DOF = 11;

	// simulation specifications
	int num_iter;
	T delta_t;
	T tol;
	int max_iter;
	T FC;
	T* r;
	T* r_tilda;
	T* FW;
	T* FT;
	T* FR;
	T* lower_spring_length;
	T* upper_spring_length;
	T* lower_spring_stiffness;
	T* upper_spring_stiffness;
	T* A; // The system matrix A. 0:3 - 2x2 upper-left symmetric block; 4:12 . We will store its inverse on diagonal blocks
	T* Ic;
	T* x; // the solution vector
	//
	T* Tc; // 4 * C_Nc; external torque on the car body (use later for rotational damping)
	T* basis_car;
	T* C_Nc;
	T* r_global;
	T* car_corners;
	// T* ones4;
	T* upper_length;
	T* lower_length;
	T* diff_vector;
	T* upper_force;
	T* lower_force;
	T* upper_F;
	T* C_Nc_extended;
	T* sum_torque;
	T* rhs; // b vector, 11 DOF (from the 12 components we prohibit rotation along y-axis, b[2]=0 is ignored)
	T* wc_tilda; // 9 elements of \tilde{wc}
	T* Qc; // the derivative of the quaternion
public:
	ComputeNasa() { // constructor without parameters (default, testing reasons)
	// initialization
		Tc = (T*)mkl_calloc(num_wheels, sizeof(T), alignment); // 4 * C_Nc; external torque on the car body (use later for rotational damping)
		basis_car = (T*)mkl_calloc(dim_x_dim, sizeof(T), alignment); // 9 elements
		C_Nc = (T*)mkl_calloc(dim_x_dim, sizeof(T), alignment); // 9 elements
		r_global = (T*)mkl_calloc(num_wheels_x_dim, sizeof(T), alignment); // 12 elements
		car_corners = (T*)mkl_calloc(num_wheels, sizeof(T), alignment); // 4 elements
		// ones4 = (T*)mkl_calloc(num_wheels, sizeof(T), alignment);
		upper_length = (T*)mkl_calloc(num_wheels, sizeof(T), alignment);  // 4 elements
		lower_length = (T*)mkl_calloc(num_wheels, sizeof(T), alignment);  // 4 elements
		diff_vector = (T*)mkl_calloc(num_wheels, sizeof(T), alignment);  // 4 elements
		upper_force = (T*)mkl_calloc(num_wheels, sizeof(T), alignment);  // 4 elements
		lower_force = (T*)mkl_calloc(num_wheels, sizeof(T), alignment);  // 4 elements
		upper_F = (T*)mkl_calloc(num_wheels_x_dim, sizeof(T), alignment);  // 12 elements
		C_Nc_extended = (T*)mkl_calloc(dim * dim * num_wheels, sizeof(T), alignment); // 4 * C_Nc;   // 36 elements
		sum_torque = (T*)mkl_calloc(dim, sizeof(T), alignment);   // 3 elements
		rhs = (T*)mkl_calloc(11, sizeof(T), alignment); // b vector, 11 DOF (from the 12 components we prohibit rotation along y-axis, b[2]=0 is ignored)
		wc_tilda = (T*)mkl_calloc(dim_x_dim, sizeof(T), alignment); // 9 elements of \tilde{wc}
		Qc = (T*)mkl_calloc(12, sizeof(T), alignment); // the derivative of the quaternion
	// simulation specifications
		num_iter = 1e3;
		delta_t = 1e-3;
		tol = 1e-10;
		max_iter = 200;
		mkl_set_num_threads(8);
		// force parameters
		const T g = 0; // actually 9.81
		// initial conditions - read XML data
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
		T mass_Body = 1936;
		T I_body_xx = 640;
		T I_body_yy = 4800;
		T mass_wheel_fl = 145 / 2;
		T mass_tyre_fl = 0;
		T mass_wheel_fr = 145 / 2;
		T mass_tyre_fr = 0;
		T mass_wheel_rl = 135 / 2;
		T mass_tyre_rl = 0;
		T mass_wheel_rr = 135 / 2;
		T mass_tyre_rr = 0;
		// Dimensions of the main car body (the center of rotation is at the origin)
		r = (T*)mkl_calloc(num_wheels * dim, sizeof(T), alignment);
		r[0] = -l_long_rr; // r1
		r[8] = l_lat_rr; // r1
		r[1] = -l_long_rl; // r2
		r[9] = l_lat_rl; // r2
		r[2] = -l_long_fl; // r3
		r[10] = l_lat_fl; // r3
		r[3] = -l_long_fr; // r4
		r[11] = l_lat_fr; // r4
		// moment of inertia of the car
		Ic = (T*)mkl_calloc(dim_x_dim, sizeof(T), alignment);
		Ic[0] = I_body_xx;
		Ic[4] = 1;
		Ic[8] = I_body_yy;
		// initial orientation of the car body as quaternion
		T* initial_orientation = (T*)mkl_calloc(4, sizeof(T), alignment);
		initial_orientation[3] = 1;
		// wheel parameters (provided as vectors [right-back, left-back, left-front, right_front])
		T* mass_wheel = (T*)mkl_calloc(4, sizeof(T), alignment);
		mass_wheel[0] = mass_wheel_rr;
		mass_wheel[1] = mass_wheel_rl;
		mass_wheel[2] = mass_wheel_fl;
		mass_wheel[3] = mass_wheel_fr;
		// tyre parameters (provided as vectors [right-back, left-back, left-front, right_front])
		T* mass_tyre = (T*)mkl_calloc(4, sizeof(T), alignment);
		// do not put to zero to avoid singularities, Dirichlet condition are enforced via a force vector
		mass_tyre[0] = 3;
		mass_tyre[1] = 3;
		mass_tyre[2] = 3;
		mass_tyre[3] = 3;
		// lengths of the springs around the tyre
		lower_spring_length = (T*)mkl_calloc(4, sizeof(T), alignment);
		lower_spring_length[0] = 0.2;
		lower_spring_length[1] = 0.2;
		lower_spring_length[2] = 0.2;
		lower_spring_length[3] = 0.2;
		upper_spring_length = (T*)mkl_calloc(4, sizeof(T), alignment);
		upper_spring_length[0] = 0.2;
		upper_spring_length[1] = 0.2;
		upper_spring_length[2] = 0.2;
		upper_spring_length[3] = 0.2;
		T* initial_lower_spring_length = (T*)mkl_calloc(4, sizeof(T), alignment);
		initial_lower_spring_length[0] = -0.2;
		initial_lower_spring_length[1] = -0.2;
		initial_lower_spring_length[2] = -0.2;
		initial_lower_spring_length[3] = -0.2;
		// initially the spring lengths with a negative sign as we substract by adding negative numbers afterwards
		T* initial_upper_spring_length = (T*)mkl_calloc(4, sizeof(T), alignment);
		initial_upper_spring_length[0] = -0.2;
		initial_upper_spring_length[1] = -0.2;
		initial_upper_spring_length[2] = -0.2;
		initial_upper_spring_length[3] = -0.2;
		lower_spring_stiffness = (T*)mkl_calloc(4, sizeof(T), alignment);
		lower_spring_stiffness[0] = k_body_rr;
		lower_spring_stiffness[1] = k_body_rl;
		lower_spring_stiffness[2] = k_body_fl;
		lower_spring_stiffness[3] = k_body_fr;
		upper_spring_stiffness = (T*)mkl_calloc(4, sizeof(T), alignment);
		upper_spring_stiffness[0] = k_tyre_rr;
		upper_spring_stiffness[1] = k_tyre_rl;
		upper_spring_stiffness[2] = k_tyre_fl;
		upper_spring_stiffness[3] = k_tyre_fr;
		// initial velocities (only y-components)
			// can be used only x such that we don't allocate and deallocate not used small matrices
		T vc = 0;
		// init to 0
		T* vw = (T*)mkl_calloc(4, sizeof(T), alignment);
		// init ot 0
		T* vt = (T*)mkl_calloc(4, sizeof(T), alignment);
		// initial angular velocity
		T* wc = (T*)mkl_calloc(dim, sizeof(T), alignment);
		FC = -1.1e3; // actually - mass_Body * g
		FT = (T*)mkl_calloc(4, sizeof(T), alignment);
		FW = (T*)mkl_calloc(4, sizeof(T), alignment);
		FR = (T*)mkl_calloc(4, sizeof(T), alignment);
		for (int i = 0; i < 4; ++i) {
			FT[i] = -mass_tyre[i] * g;
			FW[i] = -mass_wheel[i] * g;
		}
		// ===================================== FINISHED Data Initialization =========================================
		// ======================== main nasa car.m ===============================
		// quaternion & normalization
		T* qc = (T*)mkl_calloc(4, sizeof(T), alignment); // quaternions
		qc[3] = 1; // initial orientation
		create_basis_car(qc, basis_car);
		//T* basis_N = (T*)mkl_calloc(dim_x_dim, sizeof(T), alignment);
		//basis_N[0] = 1;	basis_N[4] = 1; basis_N[8] = 1; // basis_N = eye(3)
		// change of basis matrices
		C_cos_transf(basis_car, C_Nc);
		// initial position of the car if it is different from 0 it has to be added to pc_1!!!
		double pcc = 0.;
		// positions of the lower corners of the car body
		T* pc = (T*)mkl_calloc(num_wheels * dim, sizeof(T), alignment);
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dim, num_wheels, dim, 1., C_Nc, dim, r, num_wheels, 0., pc, num_wheels);
		// position of the wheels
		T* pw = (T*)mkl_calloc(num_wheels * dim, sizeof(T), alignment);
		T* pt = (T*)mkl_calloc(num_wheels * dim, sizeof(T), alignment);
		cblas_dcopy(num_wheels * dim, pc, 1, pw, 1);
		vdAdd(num_wheels, pw + 4, initial_upper_spring_length, pw + 4);
		cblas_dcopy(num_wheels * dim, pw, 1, pt, 1);
		vdAdd(num_wheels, pt + 4, initial_lower_spring_length, pt + 4);
		// get r_tilda for each component - r1, r2, r3, r4; r_tilda is a (3x3)x4 matrix, r_tilda=[r1_tilda, r2_tilda, r3_tilda, r4_tilda]
		r_tilda = (T*)mkl_calloc(num_wheels * dim_x_dim, sizeof(T), alignment);
		get_tilda_r(r, dim, num_wheels, r_tilda);
		// initialize x vector (the solution vector)
		x = (T*)mkl_calloc(25, sizeof(T), alignment);
		// x = [wc (3x1), vc(1x1), vw(4x1), vt(4x1), pt(4x1), pcc(1x1), pw(4x1) (only y component), pt(4x1) (only y component)]
		cblas_dcopy(dim, wc, 1, x, 1);
		x[3] = vc;
		cblas_dcopy(num_wheels, vw, 1, x + 4, 1);
		cblas_dcopy(num_wheels, vt, 1, x + 8, 1);
		cblas_dcopy(dim + 1, qc, 1, x + 12, 1);
		x[16] = pcc;
		cblas_dcopy(num_wheels, pw + num_wheels, 1, x + 17, 1);
		cblas_dcopy(num_wheels, pt + num_wheels, 1, x + 21, 1);
		// create A
		printf("Ic[0]=%lf, Ic[2]=%lf, Ic[8]=%lf\n\n", Ic[0], Ic[2], Ic[8]);
		A = (T*)mkl_calloc(13, sizeof(T), alignment); // 0:3 - 2x2 upper-left symmetric block; 4:12 - diag elements
		// compute the explicit inverse of the symmetric 2x2 upper-left block from A
			// A_upper_left = [Ixx, Ixz; Izx=Ixz, Izz] = [Ic[0], Ic[2]; Ic[2], Ic[8]]
			// A_inverse_2x2 = 1/(Ic[0]*Ic[8]-Ic[2]*Ic[2])*[Ic[8], -Ic[2]; -Ic[2], Ic[0]]
		double det_2x2_block_inverse = 1. / (Ic[0] * Ic[8] - Ic[2] * Ic[2]);
		A[0] = det_2x2_block_inverse * Ic[8]; // A_upper_left[0]
		A[2] = A[1] = -det_2x2_block_inverse * Ic[2]; // A_upper_left[1]=A_upper_left[2]
		A[3] = det_2x2_block_inverse * Ic[0]; // A_upper_left[3]
		// Store the inverse of the diagonal elements
			// A[4:12] is diagonal: = diag([mass_Body, mass_wheel, mass_tyre])
		A[4] = mass_Body;
		cblas_dcopy(num_wheels, mass_wheel, 1, A + 5, 1);
		cblas_dcopy(num_wheels, mass_tyre, 1,  A + 9, 1);
		// We store the inverse of A[4:12]:= diag([1./mass_Body, 1./mass_wheel, 1./mass_tyre])
		MathLibrary::elementwise_inversion<T>(A + 4, 1 + num_wheels + num_wheels);
		// free the allocated elements; mostly 3 or 4 elements
		MKL_free(pw);
		MKL_free(pt);
		MKL_free(pc);
		MKL_free(qc);
		MKL_free(vw);
		MKL_free(vt);
		MKL_free(wc);
		MKL_free(mass_tyre);
		MKL_free(mass_wheel);
		MKL_free(initial_lower_spring_length);
		MKL_free(initial_upper_spring_length);
	}
	ComputeNasa(T* r, T* r_tilda, T* FW, T* FT, T* FR, T* lower_spring_length, T* upper_spring_length, T* lower_spring_stiffness, T* upper_spring_stiffness, T* A, T FC, T* IC, T* x) :
		r(r), r_tilda(r_tilda), FW(FW), FT(FT), FR(FR), lower_spring_length(lower_spring_length), upper_spring_length(upper_spring_length), lower_spring_stiffness(lower_spring_stiffness), upper_spring_stiffness(upper_spring_stiffness), A(A), FC(FC), IC(IC), x(x) {
		// memset(ones4, 1, num_wheels);
		num_iter = 1e3;
		delta_t = 1e-3;
		tol = 1e-10;
		max_iter = 200;
	}

	~ComputeNasa() {
		MKL_free(r);
		MKL_free(r_tilda);
		MKL_free(FW);
		MKL_free(FT);
		MKL_free(FR);
		MKL_free(lower_spring_length);
		MKL_free(upper_spring_length);
		MKL_free(lower_spring_stiffness);
		MKL_free(upper_spring_stiffness);
		MKL_free(A);
		MKL_free(Ic);
		MKL_free(x);
		MKL_free(Tc);
		////
		MKL_free(basis_car);
		MKL_free(C_Nc);
		MKL_free(r_global);
		MKL_free(car_corners);
		// MKL_free(ones4)
		MKL_free(upper_length);
		MKL_free(lower_length);
		MKL_free(diff_vector);
		MKL_free(upper_force);
		MKL_free(lower_force);
		MKL_free(upper_F);
		MKL_free(C_Nc_extended);
		MKL_free(sum_torque);
		MKL_free(rhs);
		MKL_free(wc_tilda);
		MKL_free(Qc);
	}
	void compute_f3D_reduced(T time, T* x, T* f) {

		// extract data from x_vector
		// wc = x_vector(1:3);
		// vc = x_vector(4);
		// vw = x_vector(5:8);
		// vt = x_vector(9:12);
		// qc = x_vector(13:16);
		// pcc = x_vector(17);

		// =================================== UPDATE rhs vector ==============================================
		// get the new local basis vectors using quaternion rotation
		create_basis_car(x + 12, basis_car);
		// compute the cosine transformation
		C_cos_transf(basis_car, C_Nc);
		// compute r_global = C_Nc * r; (matrix-vector multiplication)
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dim, num_wheels, dim, 1., C_Nc, dim, r, num_wheels, 0., r_global, num_wheels);
		// get positions of the corners of the car body
		// car_corners = pcc + [r1_global(2); r2_global(2); r3_global(2); r4_global(2)];
		cblas_dcopy(num_wheels, r_global + 4, 1, car_corners, 1);
		for (auto i = 0; i < num_wheels; ++i) {
			car_corners[i] += x[16];
		}
		// calculate the length of the springs
		// upper_length = car_corners - x_vector(18:21);
		cblas_dcopy(num_wheels, car_corners, 1, upper_length, 1);
		cblas_daxpy(num_wheels, -1, x + 17, 1, upper_length, 1);
		// lower_length = x_vector(18:21) - x_vector(22:25);
		cblas_dcopy(num_wheels, x + 17, 1, lower_length, 1);
		cblas_daxpy(num_wheels, -1, x + 21, 1, lower_length, 1);
		// calculate the spring forces (take opposite if top node is considered)
		// upper_force = upper_spring_stiffness .* (upper_length - upper_spring_length);
		cblas_dcopy(num_wheels, upper_length, 1, diff_vector, 1);
		cblas_daxpy(num_wheels, -1, upper_spring_length, 1, diff_vector, 1);
		vdmul(&num_wheels, upper_spring_stiffness, diff_vector, upper_force);

		// lower_force = lower_spring_stiffness .* (lower_length - lower_spring_length);
		cblas_dcopy(num_wheels, lower_length, 1, diff_vector, 1);
		cblas_daxpy(num_wheels, -1, lower_spring_length, 1, diff_vector, 1);
		vdmul(&num_wheels, lower_spring_stiffness, diff_vector, lower_force);
		// get road forces: FR := -(lower_force + FT); for the moment, we don't have the minus (-)
		vdAdd(num_wheels, lower_force, FT, FR);
		// convert spring forces to local basis
			// upper_F = C_Nc' * [0; -upper_force(1); 0; 0;-upper_force(2); 0; 0; -upper_force(3); 0; 0; -upper_force(4); 0];
			// upper_F = [C_Nc(2,:) * upper_force(1); C_Nc(2,:) * upper_force(2); C_Nc(2,:) * upper_force(3); C_Nc(2,:) * upper_force(4)]
			//					only the second line from C_Nc counts
			// upper_F is a (3x4) vector, 3 components for each wheel
		for (auto i = 0; i < num_wheels; ++i) {
			// upper_F = [upper_force[0], upper_force[0], upper_force[0], upper_force[1], upper_force[1], upper_force[1],
			//			  upper_force[2], upper_force[2], upper_force[2], upper_force[3], upper_force[3], upper_force[3]]
			upper_F[dim * i + 2] = upper_F[dim * i + 1] = upper_F[dim * i] = upper_force[i];
			// C_Nc_extended = [C_Nc, C_Nc, C_Nc, C_Nc]
			cblas_dcopy(dim, C_Nc, 1, C_Nc_extended + dim*i, 1);
		}
		// perform element-wise vector multiplication: upper_F <- C_Nc_extended .* upper_F
		vdmul(&num_wheels_x_dim, upper_F, C_Nc_extended, upper_F);
		// perform matrix-vector multiplication: sum_torque = r_tilda * upper_F = ((3x3)x4) x (3x4)x1 = (3x1)
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dim, 1, dim*num_wheels, 1., r_tilda, dim, upper_F, dim*num_wheels, 0., sum_torque, dim);
		// get Hc = Ic * wc; however, we store Hc directly in the first 3 components of rhs (b); wc = [x[0], x[1], x[2]]
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dim, 1, dim, 1., Ic, dim, x, dim, 0., rhs, dim);
		// get first "3" (actually 2) components of rhs (b)
			// rhs(1:3) = -get_tilda(wc)*Hc
		get_tilda(x, wc_tilda);
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dim, 1, dim, -1., wc_tilda, dim, rhs, dim, 0., rhs, dim);
			// rhs(1:3) += sum_torque
		vdAdd(dim, rhs, sum_torque, rhs);
			// rhs(1:3) =+ sum(Tc) not implemented!!!
		rhs[1] = rhs[2]; // suppresing the y-component (b(2)=0)

		// rhs[2] = FC-sum(upper_force)
			// may be optimized: see ippsSum_64f((Ipp64f *) upper_force, num_wheels, (Ipp64f *) (rhs +2));
		rhs[2] = FC;
		for (auto i = 0; i < num_wheels; ++i) {
			rhs[2] -= upper_force[i];
		}
		// rhs[3:6] = upper_force - lower_force + FW = (upper_force + FW) - lower_force
			// rhs[3:6] = upper_force + FW
		vdAdd(num_wheels, upper_force, FW, rhs + 3);
			// rhs[3:6] -= lower_force
		cblas_daxpy(num_wheels, -1., lower_force, 1, rhs + 3, 1);

		// rhs[7:10] = lower_force + FR + FT = (lower_force + FT) - (FR) - we defined FR without the minus
			// rhs[7:10] = lower_force + FT
		vdAdd(num_wheels, lower_force, FT, rhs + 7);
			// rhs[7:10] -= FR
		cblas_daxpy(num_wheels, -1., FR, 1, rhs + 7, 1);
		// ====================== FINISHED Computing rhs (b) ==================================================
		// ===================== SOLVING The Linear System: x = A (inverse) * b ===============================
		// f = A_inv * rhs -> A_inv is our A that contains the inverse on blocks
			// solve the upper_left system with the 2x2 symmetric block from A
			// f[0;1] = A[0,1;2,3] * rhs[0;1]
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 2, 1, 2, 1., A, 2, rhs, 2, 0., f, 2);
		// impose the 0 in the y element of wc_dot block from f -> just for the moment!
		f[2] = f[1];
		f[1] = 0;
			// f[3;...;11] = A[4;5;...;12] * rhs[2;3;...;10]
		vdMul(9, A + 4, rhs + 2, f + 3);
		// ===================== FINISHED SOLVING The Linear System ===========================================
		// ================== Prepair the second part of f - accelerations become velocities ==================
		get_quaternion_derivative(x + 12, Qc);
		// (f[12:15] =) qc_dot = Qc * wc; (wc=[x[0],x[1],x[2])
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 4, 1, dim, 1., Qc, 4, x, dim, 0., f + 12, dim);
		// copy from x_vector the components [vc, vw, vt] (= accelerations) to f (=> velocities)
		cblas_dcopy(num_wheels + num_wheels + 1, x + 3, 1, f + 13, 1);
		// ==================================== END OF FUNCTION ================================================
	}
	void print_solution() {
		std::cout << "Print x: \n\n";
		for (auto i = 0; i < 25; ++i) {
			std::cout << x[i] << "\n";
		}
		std::cout << "\nEnd of printing!";
	}
};
*/
/*
void EVAAComputeEngine::computeMKLNasa_example(void) {
	ComputeNasa<floatEVAA> obj_example;
	obj_example.print_solution();
	floatEVAA *f = (floatEVAA*)mkl_calloc(25, sizeof(floatEVAA), 64);
	obj_example.compute_f3D_reduced(1, obj_example.x, f);
	std::cout << "\n\nFunction after first computation, f = \n\n\n";
	for (auto i = 0; i < 25; ++i) {
		std::cout << f[i] << "\n";
	}
	MKL_free(f);
}
*/
void EVAAComputeEngine::computeMKLlinear11dof(void) {
	floatEVAA h = _parameters.timestep;
	floatEVAA tend = _parameters.num_time_iter*h;
	floatEVAA u_init = 0.0;
	floatEVAA du_init = 0.0;
	int DOF = _parameters.DOF;
	const int alignment = 64;
	floatEVAA* soln = (floatEVAA*)mkl_calloc(DOF, sizeof(floatEVAA), alignment);
	linear11dof<floatEVAA> solver(_parameters);
	solver.apply_boundary_condition("road_force");
	solver.solve(soln);
	size_t steps = floor(tend / h);
	std::cout << "Solution after " << steps << " timesteps, f =" << std::endl;
	for (auto i = 0; i < DOF; ++i) {
		std::cout << soln[i] << std::endl;
	}
	mkl_free(soln);
}

void EVAAComputeEngine::computeMKLlinear11dof_reduced(void) {
	floatEVAA tend = 1.0;
	floatEVAA h = 1.0 / 1000.0;
	floatEVAA u_init = 0.0;
	floatEVAA du_init = 0.0;
	const int alignment = 64;
	floatEVAA* soln = (floatEVAA*)mkl_calloc(7, sizeof(floatEVAA), alignment);
	linear11dof<floatEVAA> solver(_parameters);
	solver.apply_boundary_condition("fixed_to_road");
	solver.solve(soln);
	std::cout << "Solution after 1000 timesteps, f =" << std::endl;
	for (auto i = 0; i < 7; ++i) {
		std::cout << soln[i] << std::endl;
	}
	std::cout << std::endl;
	mkl_free(soln);
}

void EVAAComputeEngine::computeBlaze11DOF(void) {
#ifdef USE_BLAZE

	ReadXML reader(_xmlFileName);

	using blaze::CompressedMatrix;
	using blaze::DynamicMatrix;
	using blaze::SymmetricMatrix;
	using blaze::DynamicVector;
	using blaze::columnVector;
	using blaze::rowMajor;

	int DOF = _parameters.DOF;

	// K stiffness
	double k_11 = _parameters.k_body[0];
	double k_12 = _parameters.k_tyre[0];
	double k_21 = _parameters.k_body[1];
	double k_22 = _parameters.k_tyre[1];
	double k_31 = _parameters.k_body[2];
	double k_32 = _parameters.k_tyre[2];
	double k_41 = _parameters.k_body[3];
	double k_42 = _parameters.k_tyre[3];
	double l_1 = _parameters.l_lat[0];
	double l_2 = _parameters.l_lat[2];
	double l_3 = _parameters.l_long[0];
	double l_4 = _parameters.l_long[2];


	// D - damping matrix
	double d_11 = _parameters.c_body[0];
	double d_12 = _parameters.c_tyre[0];
	double d_21 = _parameters.c_body[1];
	double d_22 = _parameters.c_tyre[1];
	double d_31 = _parameters.c_body[2];
	double d_32 = _parameters.c_tyre[2];
	double d_41 = _parameters.c_body[3];
	double d_42 = _parameters.c_tyre[3];
	double d_l1 = 0.1 * _parameters.l_lat[0];
	double d_l2 = 0.1 * _parameters.l_lat[2];
	double d_l3 = 0.1 * _parameters.l_long[0];
	double d_l4 = 0.1 * _parameters.l_long[2];

	// M - mass matrix
	double m_1 = _parameters.mass_body;
	double m_2 = _parameters.I_body[0];
	double m_3 = _parameters.I_body[2];
	double m_4 = _parameters.mass_wheel[2];
	double m_5 = _parameters.mass_tyre[2];
	double m_6 = _parameters.mass_wheel[3];
	double m_7 = _parameters.mass_tyre[3];
	double m_8 = _parameters.mass_wheel[1];
	double m_9 = _parameters.mass_tyre[1];
	double m_10 = _parameters.mass_wheel[0];
	double m_11 = _parameters.mass_tyre[0];

	// Using symmetric matrix enforce the symmetry and is safer
	SymmetricMatrix <DynamicMatrix<double>, rowMajor> K(DOF);

	K(0, 0) = k_11 + k_21 + k_31 + k_41;
	K(0, 1) = -k_11 * l_1 - k_41 * l_1 + k_21 * l_2 + k_31 * l_2;
	K(0, 2) = -k_11 * l_3 - k_21 * l_3 + k_31 * l_4 + k_41 * l_4;
	K(0, 3) = -k_11;
	K(0, 5) = -k_21;
	K(0, 7) = -k_31;
	K(0, 9) = -k_41;
	K(1, 1) = l_1 * l_1 * k_11 + l_2 * l_2 * k_21 + l_2 * l_2 * k_31 + l_1 * l_1 * k_41;
	K(1, 2) = l_1 * l_3 * k_11 - l_3 * l_2 * k_21 + l_2 * l_4 * k_31 - l_1 * l_4 * k_41;
	K(1, 3) = l_1 * k_11;
	K(1, 5) = -l_2 * k_21;
	K(1, 7) = -l_2 * k_31;
	K(1, 9) = l_1 * k_41;
	K(2, 2) = l_3 * l_3 * k_11 + l_3 * l_3 * k_21 + l_4 * l_4 * k_31 + l_4 * l_4 * k_41;
	K(2, 3) = l_3 * k_11;
	K(2, 5) = l_3 * k_21;
	K(2, 7) = -l_4 * k_31;
	K(2, 9) = -l_4 * k_41;
	K(3, 3) = k_11 + k_12;
	K(3, 4) = -k_12;
	K(4, 4) = k_12;
	K(5, 5) = k_21 + k_22;
	K(5, 6) = -k_22;
	K(6, 6) = k_22;
	K(7, 7) = k_31 + k_32;
	K(7, 8) = -k_32;
	K(8, 8) = k_32;
	K(9, 9) = k_41 + k_42;
	K(9, 10) = -k_42;
	K(10, 10) = k_42;

	SymmetricMatrix <DynamicMatrix<double>, rowMajor> D(DOF);
	D(0, 0) = d_11 + d_21 + d_31 + d_41;
	D(0, 1) = -d_11 * l_1 - d_41 * l_1 + d_21 * l_2 + d_31 * l_2;
	D(0, 2) = -d_11 * l_3 - d_21 * l_3 + d_31 * l_4 + d_41 * l_4;
	D(0, 3) = -d_11;
	D(0, 5) = -d_21;
	D(0, 7) = -d_31;
	D(0, 9) = -d_41;
	D(1, 1) = l_1 * l_1 * d_11 + l_2 * l_2 * d_21 + l_2 * l_2 * d_31 + l_1 * l_1 * d_41;
	D(1, 2) = l_1 * l_3 * d_11 - l_3 * l_2 * d_21 + l_2 * l_4 * d_31 - l_1 * l_4 * d_41;
	D(1, 3) = l_1 * d_11;
	D(1, 5) = -l_2 * d_21;
	D(1, 7) = -l_2 * d_31;
	D(1, 9) = l_1 * d_41;
	D(2, 2) = l_3 * l_3 * d_11 + l_3 * l_3 * d_21 + l_4 * l_4 * d_31 + l_4 * l_4 * d_41;
	D(2, 3) = l_3 * d_11;
	D(2, 5) = l_3 * d_21;
	D(2, 7) = -l_4 * d_31;
	D(2, 9) = -l_4 * d_41;
	D(3, 3) = d_11 + d_12;
	D(3, 4) = -d_12;
	D(4, 4) = d_12;
	D(5, 5) = d_21 + d_22;
	D(5, 6) = -d_22;
	D(6, 6) = d_22;
	D(7, 7) = d_31 + d_32;
	D(7, 8) = -d_32;
	D(8, 8) = d_32;
	D(9, 9) = d_41 + d_42;
	D(9, 10) = -d_42;
	D(10, 10) = d_42;

	SymmetricMatrix <DynamicMatrix<double>, rowMajor> M(DOF);
	M(0, 0) = m_1;
	M(1, 1) = m_2;
	M(2, 2) = m_3;
	M(3, 3) = m_4;
	M(4, 4) = m_5;
	M(5, 5) = m_6;
	M(6, 6) = m_7;
	M(7, 7) = m_8;
	M(8, 8) = m_9;
	M(9, 9) = m_10;
	M(10, 10) = m_11;

	// Define the solution vectors
	DynamicVector<double, columnVector> u_n_p_1(DOF, (double)0.0); // initialize vector of dimension 11 and null elements
	DynamicVector<double, columnVector> u_n(DOF, (double)0.0); // initialize vector of dimension 11 and null elements
	DynamicVector<double, columnVector> u_n_m_1(DOF, (double)0.0); // initialize vector of dimension 11 and null elements

	// Perform the iterations
	//	int nRefinement = 10;
	//	int numTimeSteps = pow(2, nRefinement);
	int numTimeSteps = _parameters.num_time_iter;

	//time step size 
	//double h =  / (numTimeSteps);
	double h = _parameters.timestep;
	std::cout << "Time step h is: " << h << std::scientific << std::endl;

	// Initial conditions
	//const double u_init = 1;
	//const double du_init = 0;
	u_n[0] = 1;
	u_n_m_1[0] = 1;

	/// Build dynamic stiffness matrix
	// A = (1.0/(h*h))*M + (1.0/h)*D + K
	DynamicMatrix<double, rowMajor> A(DOF, DOF);
	A = 1. / (h * h) * M + 1. / h * D + K;

	///Build rhs for BE integrator
	//B = (2.0 / (h*h))*M + 1.0 / h * D + B
	DynamicMatrix<double, rowMajor> B(DOF, DOF);
	B = 2. / (h * h) * M + 1. / h * D;

	// LU Decomposition
	DynamicVector<int, columnVector> ipiv(DOF);   // Pivoting indices
	DynamicMatrix<double, rowMajor>  A_LU(A);  // Temporary matrix to be decomposed

	getrf(A_LU, ipiv.data());
	M *= -1. / (h * h);
	for (int iTime = 0; iTime < numTimeSteps; iTime++) {

		// Solve system: A*u_n_p_1 = B*u_n - M*u_n_m_1
		// rhs = B*u_n + (-1.0 / (h*h))*M
		u_n_p_1 = B * u_n + M * u_n_m_1; // rhs

		getrs(A_LU, u_n_p_1, 'T', ipiv.data());
		u_n_m_1 = u_n;
		u_n = u_n_p_1;
		reader.ReadVariableParameters(_parameters);
	}
	std::cout << "We ran #" << numTimeSteps << " time steps!" << std::endl;
	std::cout << u_n_p_1 << std::scientific << std::endl;

#endif
}

void EVAAComputeEngine::computeMBD(void) {
	size_t num_iter = _parameters.num_time_iter;
	const int alignment = 64;
	size_t solution_dim = _parameters.solution_dim;
	floatEVAA* soln = (floatEVAA*)mkl_calloc(solution_dim, sizeof(floatEVAA), alignment);
	MBD_method<floatEVAA> solver(_parameters);
	solver.solve(soln);
	std::cout << "Solution after " << num_iter << " timesteps, f =" << std::endl;
	for (auto i = 0; i < solution_dim; ++i) {
		std::cout << soln[i] << std::endl;
	}
	std::cout << std::endl;
	mkl_free(soln);
}

void create_basis_car(floatEVAA* qc, floatEVAA* basis_c) {
	const MKL_INT n = 4, incx = 1;
	floatEVAA q_norm_inv = 2 / (floatEVAA)cblas_ddot(n, qc, incx, qc, incx); // 2 / (||q||_2)^2
	basis_c[0] = 1 - q_norm_inv * (qc[1] * qc[1] + qc[2] * qc[2]);
	basis_c[1] = q_norm_inv * (qc[0] * qc[1] + qc[2] * qc[3]);
	basis_c[2] = q_norm_inv * (qc[0] * qc[2] - qc[1] * qc[3]);
	basis_c[3] = q_norm_inv * (qc[0] * qc[1] - qc[2] * qc[3]);
	basis_c[4] = 1 - q_norm_inv * (qc[0] * qc[0] + qc[2] * qc[2]);
	basis_c[5] = q_norm_inv * (qc[1] * qc[2] + qc[0] * qc[3]);
	basis_c[6] = q_norm_inv * (qc[0] * qc[2] + qc[1] * qc[3]);
	basis_c[7] = q_norm_inv * (qc[1] * qc[2] - qc[0] * qc[3]);
	basis_c[8] = 1 - q_norm_inv * (qc[0] * qc[0] + qc[1] * qc[1]);
}

void get_tilda(floatEVAA* x, floatEVAA* x_tilda) {
	// x - vector 3 elems (floatEVAA* x); x - tilda matrix 3x3 (floatEVAA* x_tilda)
	// given y: x_tilda*y=cross(x,y) (Stoneking, page 3 bottom)
	x_tilda[0] = 0;
	x_tilda[1] = -x[2];
	x_tilda[2] = x[1];
	x_tilda[3] = x[2];
	x_tilda[4] = 0;
	x_tilda[5] = -x[0];
	x_tilda[6] = -x[1];
	x_tilda[7] = x[0];
	x_tilda[8] = 0;
}

void get_tilda_r(floatEVAA* r, const int dim, const int num_wheels, floatEVAA* r_tilda) {
	// r - matrix 3x4 elems (floatEVAA* r); r - tilda matrix 3x12 (floatEVAA* x_tilda)
	// given y: x_tilda*y=cross(x,y) (Stoneking, page 3 bottom)
	for (int i = 0; i < num_wheels; i++) {
		r_tilda[0 + i * dim] = 0;
		r_tilda[1 + i * dim] = -r[2 * num_wheels + i];
		r_tilda[2 + i * dim] = r[num_wheels + i];
		r_tilda[dim * num_wheels + i * dim] = r[2 * num_wheels + i];
		r_tilda[dim * num_wheels + 1 + i * dim] = 0;
		r_tilda[dim * num_wheels + 2 + i * dim] = -r[i];
		r_tilda[2 * dim * num_wheels + i * dim] = -r[num_wheels + i];
		r_tilda[2 * dim * num_wheels + 1 + i * dim] = r[i];
		r_tilda[2 * dim * num_wheels + 2 + i * dim] = 0;
	}
}

void C_cos_transf(floatEVAA* Y, floatEVAA* C_transf) {
	// This function defines the cosine transformation matrix, used in changing coordinates from a frame to another. The X and Y matrices represent basis in the X and, respectively, Y frame. x = C * y
	// Cosinus transformation to N base (eye(3)). X basis is the identity implicitly; Y basis is the basis_car
	// C = [ Y_1 /||Y_1||, Y_2 /||Y_2||, Y_3 /||Y_3||]; Y_i - columns of Y
	floatEVAA Y_norms_on_columns;
	for (auto i = 0; i < 3; ++i) {
		Y_norms_on_columns = 1. / (floatEVAA)(cblas_dnrm2(3, Y + i, 3));
		Y[i] /= Y_norms_on_columns;
		Y[i + 3] /= Y_norms_on_columns;
		Y[i + 6] /= Y_norms_on_columns;
	}

}

void get_quaternion_derivative(floatEVAA* qc, floatEVAA* Qc) {
	// input: qc - the quaternion, 4x1 vector
	// output: Qc - the matrix of quaternion's derivative; 4x3 matrix in vector representation with row major
	Qc[0] = 0.5 * qc[3];
	Qc[1] = -0.5 * qc[2];
	Qc[2] = 0.5 * qc[1];

	Qc[3] = 0.5 * qc[2];
	Qc[4] = 0.5 * qc[3];
	Qc[5] = -0.5 * qc[0];

	Qc[6] = -0.5 * qc[1];
	Qc[7] = 0.5 * qc[0];
	Qc[8] = 0.5 * qc[3];

	Qc[9] = -0.5 * qc[0];
	Qc[10] = -0.5 * qc[1];
	Qc[11] = -0.5 * qc[2];
}


void EVAAComputeEngine::clean(void) {

}