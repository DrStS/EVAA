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


EVAAComputeEngine::EVAAComputeEngine(std::string _xmlFileName) {
	// Intialize XML metadatabase singelton 
}


EVAAComputeEngine::~EVAAComputeEngine() {
}

void EVAAComputeEngine::prepare(void) {
	MathLibrary::printMKLInfo();
}

// EIGEN
void EVAAComputeEngine::computeEigen(void) {

	int DOF = 3;

	typedef double floatEVAA;
	floatEVAA k_1 = 1.;
	floatEVAA k_2 = 2.;
	floatEVAA k_3 = 3.;
	floatEVAA d_1 = 1. / 10.;
	floatEVAA d_2 = 1. / 2.;
	floatEVAA m_1 = 1e-1;
	floatEVAA m_2 = 2e-1;
	floatEVAA m_3 = 3e-1;

	MatrixXd B(DOF, DOF);
	MatrixXd M(DOF, DOF);
	MatrixXd D(DOF, DOF);
	MatrixXd K(DOF, DOF);

	VectorXd u_n_p_1(DOF);
	VectorXd u_n(DOF);
	VectorXd u_n_m_1(DOF);

	M << m_1, 0.0, 0.0,
		0.0, m_2, 0.0,
		0.0, 0.0, m_3;
	D << d_1, 0.0, 0.0,
		0.0, d_2, -d_2,
		0.0, -d_2, d_2;
	K << k_1 + k_2, -k_2, 0.0,
		-k_2, k_2, 0.0,
		0.0, 0.0, k_3;
	B << 0.0, 0.0, 0.0,
		0.0, 0.0, 0.0,
		0.0, 0.0, 0.0;

	u_n_p_1 << 0.0, 0.0, 0.0;
	u_n << 1.0, 0.0, 0.0;
	u_n_m_1 << 1.0, 0.0, 0.0;

	int nRefinement = 10;
	int numTimeSteps = pow(2, nRefinement);
	//time step size 
	floatEVAA h = 1.0 / (numTimeSteps);
	std::cout << "Time step h is: " << h << std::scientific << std::endl;

	/// Build dynamic stiffness matrix
	// K = (1.0/(h*h))*M + (1.0/h)*D + K
	K += (1.0 / (h * h)) * M + 1.0 / h * D;

	///Build rhs for BE integrator
	//B = (2.0 / (h*h))*M + 1.0 / h * D + B
	B += 2.0 / (h * h) * M + 1.0 / h * D;


	// LU Decomposition
	Eigen::PartialPivLU<Matrix3d> lu(K);

	M *= (-1.0 / (h * h));
	for (int iTime = 0; iTime < numTimeSteps; iTime++) {
		// B*u_n + (-1.0 / (h*h))*M
		// Solve system
		u_n_p_1 = lu.solve(B * u_n + M * u_n_m_1);

		u_n_m_1 = u_n;
		u_n = u_n_p_1;
	}
	std::cout << "We ran #" << numTimeSteps << " time steps!" << std::endl;
	std::cout << u_n_p_1 << std::scientific << std::endl;
}

void EVAAComputeEngine::computeEigen11DOF(void) {

	int DOF = 11;

	typedef double floatEVAA;
	// K stiffness
	floatEVAA k = 10.;
	floatEVAA k_11 = 1.1;
	floatEVAA k_12 = k;
	floatEVAA k_21 = 1.2;
	floatEVAA k_22 = k;
	floatEVAA k_31 = 1.3;
	floatEVAA k_32 = k;
	floatEVAA k_41 = 1.4;
	floatEVAA k_42 = k;
	floatEVAA k_l1 = 2.;
	floatEVAA k_l2 = 1.;
	floatEVAA k_l3 = 0.8;
	floatEVAA k_l4 = 1.25;
	floatEVAA l_1 = 2;
	floatEVAA l_2 = 1;
	floatEVAA l_3 = 0.8;
	floatEVAA l_4 = 1.25;

	// D - damping matrix
	floatEVAA d_11 = 1.1 * 0.1;
	floatEVAA d_12 = k * 0.1;
	floatEVAA d_21 = 1.2 * 0.1;
	floatEVAA d_22 = k * 0.1;
	floatEVAA d_31 = 1.3 * 0.1;
	floatEVAA d_32 = k * 0.1;
	floatEVAA d_41 = 1.4 * 0.1;
	floatEVAA d_42 = k * 0.1;
	floatEVAA d_l1 = 2. * 0.1;
	floatEVAA d_l2 = 1. * 0.1;
	floatEVAA d_l3 = 0.8 * 0.1;
	floatEVAA d_l4 = 1.25 * 0.1;

	// M - mass matrix
	floatEVAA m_1 = 1e-1;
	floatEVAA m_2 = 2e-1;
	floatEVAA m_3 = 3e-1;
	floatEVAA m_4 = 4e-1;
	floatEVAA m_5 = 5e-1;
	floatEVAA m_6 = 6e-1;
	floatEVAA m_7 = 7e-1;
	floatEVAA m_8 = 8e-1;
	floatEVAA m_9 = 9e-1;
	floatEVAA m_10 = 10e-1;
	floatEVAA m_11 = 11e-1;

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
	/*K_diag = MatrixXd::Zero(DOF, DOF);
	K_diag.diagonal = K.diagonal();*/
	//std::cout << "K is: " << K_diag << std::endl;
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
	/*D_diag = MatrixXd::Zero(DOF, DOF);
	D_diag.diagonal = D.diagonal();*/
	D_aux = D + D.transpose();
	D_aux.diagonal() -= D.diagonal();
	D = D_aux;
	A = MatrixXd::Zero(DOF, DOF);
	B = MatrixXd::Zero(DOF, DOF);

	u_n_p_1 = VectorXd::Zero(DOF);
	u_n = VectorXd::Zero(DOF);
	u_n(0) = 1;
	u_n_m_1 = u_n;

	int nRefinement = 10;
	int numTimeSteps = pow(2, nRefinement);
	//time step size 
	floatEVAA h = 1.0 / (numTimeSteps);
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
	}
	std::cout << "We ran #" << numTimeSteps << " time steps!" << std::endl;
	std::cout << u_n_p_1 << std::scientific << std::endl;
}

// INTEL MKL
void EVAAComputeEngine::computeMKL(void) {

	typedef double floatEVAA;
	floatEVAA k_1 = 1.;
	floatEVAA k_2 = 2.;
	floatEVAA k_3 = 3.;
	floatEVAA d_1 = 1. / 10.;
	floatEVAA d_2 = 1. / 2.;
	floatEVAA m_1 = 1e-1;
	floatEVAA m_2 = 2e-1;
	floatEVAA m_3 = 3e-1;

	int alignment = 64;
	int DOF = 3;
	int matrixElements = DOF * DOF;

	floatEVAA* B = (floatEVAA*)mkl_malloc(sizeof(floatEVAA) * matrixElements, alignment);
	floatEVAA* M = (floatEVAA*)mkl_malloc(sizeof(floatEVAA) * matrixElements, alignment);
	floatEVAA* D = (floatEVAA*)mkl_malloc(sizeof(floatEVAA) * matrixElements, alignment);
	floatEVAA* K = (floatEVAA*)mkl_malloc(sizeof(floatEVAA) * matrixElements, alignment);

	//std::vector<floatEVAA> B(9);
	//std::vector<floatEVAA> M(9);
	//std::vector<floatEVAA> D(9);
	//std::vector<floatEVAA> K(9);

	floatEVAA* u_n_p_1 = (floatEVAA*)mkl_malloc(sizeof(floatEVAA) * DOF, alignment);
	floatEVAA* u_n = (floatEVAA*)mkl_malloc(sizeof(floatEVAA) * DOF, alignment);
	floatEVAA* u_n_m_1 = (floatEVAA*)mkl_malloc(sizeof(floatEVAA) * DOF, alignment);
	floatEVAA* tmp = (floatEVAA*)mkl_malloc(sizeof(floatEVAA) * DOF, alignment);

	//std::vector<floatEVAA> f_n_p_1(DOF);
	//std::vector<floatEVAA> u_n_p_1(DOF);
	//std::vector<floatEVAA> u_n(DOF);
	//std::vector<floatEVAA> u_n_m_1(DOF);
	//std::vector<floatEVAA> tmp(3);

	M[0] = m_1;
	M[1] = 0.0;
	M[2] = 0.0;
	M[3] = 0.0;
	M[4] = m_2;
	M[5] = 0.0;
	M[6] = 0.0;
	M[7] = 0.0;
	M[8] = m_3;

	D[0] = d_1;
	D[1] = 0.0;
	D[2] = 0.0;
	D[3] = 0.0;
	D[4] = d_2;
	D[5] = -d_2;
	D[6] = 0.0;
	D[7] = -d_2;
	D[8] = d_2;

	K[0] = k_1 + k_2;
	K[1] = -k_2;
	K[2] = 0.0;
	K[3] = -k_2;
	K[4] = k_2;
	K[5] = 0.0;
	K[6] = 0.0;
	K[7] = 0.0;
	K[8] = k_3;

	B[0] = 0.0;
	B[1] = 0.0;
	B[2] = 0.0;
	B[3] = 0.0;
	B[4] = 0.0;
	B[5] = 0.0;
	B[6] = 0.0;
	B[7] = 0.0;
	B[8] = 0.0;

	//Initial conditions
	u_n_p_1[0] = 0.;
	u_n_p_1[1] = 0.;
	u_n_p_1[2] = 0.;

	u_n[0] = 1.;
	u_n[1] = 0.;
	u_n[2] = 0.;
	//
	u_n_m_1[0] = 1.;
	u_n_m_1[1] = 0.;
	u_n_m_1[2] = 0.;

	mkl_set_num_threads(1);
	int nRefinement = 10;
	int numTimeSteps = pow(2, nRefinement);
	//time step size 
	floatEVAA h = 1.0 / (numTimeSteps);
	std::cout << "Time step h is: " << h << std::scientific << std::endl;
	/// Build dynamic stiffness matrix
	// K' <- (1.0/(h*h))*M + K
//	MathLibrary::computeDenseVectorAddition(M.data(), K.data(), (1.0 / (h*h)), 9);
	cblas_daxpy(matrixElements, (1.0 / (h * h)), M, 1, K, 1);
	// K <- (1.0/h)*D + K'
//	MathLibrary::computeDenseVectorAddition(D.data(), K.data(), (1.0 / h), 9);
	cblas_daxpy(matrixElements, (1.0 / h), D, 1, K, 1);
	/// K holds now dynamic stiffness matrix  for BE integrator

	///Build rhs for BE integrator
	//B' <-(2.0 / (h*h))*M + B
//	MathLibrary::computeDenseVectorAddition(M.data(), B.data(), (2.0 / (h*h)), 9);
	cblas_daxpy(matrixElements, (2.0 / (h * h)), M, 1, B, 1);
	//B <-(1.0 / (h))*D + B'
//	MathLibrary::computeDenseVectorAddition(D.data(), B.data(), (1.0 / h), 9);
	cblas_daxpy(matrixElements, (1.0 / h), D, 1, B, 1);
	//A*u_n_p_1=B*u_n+C*u_n_m_1+f_n_p_1 <== BE

	std::vector<int> pivot(3);
	// LU Decomposition
//	MathLibrary::computeDenseSymLUFactorisation(DOF, K, pivot);
	LAPACKE_dgetrf(LAPACK_COL_MAJOR, DOF, DOF, K, DOF, pivot.data());

	std::vector<double> timeVec;
	timeVec.resize(numTimeSteps);
	///Time loop

	double tmpScalar = (-1.0 / (h * h));
	M[0] = M[0] * tmpScalar;
	M[1] = M[1] * tmpScalar;
	M[2] = M[2] * tmpScalar;
	M[3] = M[3] * tmpScalar;
	M[4] = M[4] * tmpScalar;
	M[5] = M[5] * tmpScalar;
	M[6] = M[6] * tmpScalar;
	M[7] = M[7] * tmpScalar;
	M[8] = M[8] * tmpScalar;

	void* jitter;
	std::cout << mkl_jit_create_dgemm(&jitter, MKL_COL_MAJOR, MKL_NOTRANS, MKL_NOTRANS, DOF, 1, DOF, 1.0, DOF, DOF, 0.0, DOF) << std::endl;
	dgemm_jit_kernel_t myDGEMMKernel = mkl_jit_get_dgemm_ptr(jitter);
	LAPACKE_set_nancheck(0);
	for (int iTime = 0; iTime < numTimeSteps; iTime++) {
		//timeVec[iTime] = iTime * h;
		// y: = alpha * A*x + beta * y
		// u_n_p_1 = B*u_n
#ifndef USE_GEMM
		cblas_dgemv(CblasColMajor, CblasNoTrans, DOF, DOF, 1.0, B, DOF, u_n, 1, 0.0, u_n_p_1, 1);
#endif
#ifdef USE_GEMM 
		//cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, DOF, 1, DOF, 1.0, B, DOF, u_n, DOF, 0.0, u_n_p_1, DOF);
		myDGEMMKernel(jitter, B, u_n, u_n_p_1);
#endif
		// y: = alpha*A*x + beta*y
		// tmp = ((-1.0 / (h*h))*M) * u_n_m_1
#ifndef USE_GEMM
		cblas_dgemv(CblasColMajor, CblasNoTrans, DOF, DOF, (-1.0 / (h * h)), M, DOF, u_n_m_1, 1, 0.0, tmp, 1);
#endif
#ifdef USE_GEMM 
		//cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, DOF, 1, DOF, 1.0, M, DOF, u_n_m_1, DOF, 0.0, tmp, DOF);
		myDGEMMKernel(jitter, M, u_n_m_1, tmp);
#endif
		// u_n_p_1 <- 1.0 tmp + u_n_p_1
//		MathLibrary::computeDenseVectorAddition(tmp.data(), u_n_p_1.data(), 1.0, DOF);
		cblas_daxpy(DOF, 1.0, tmp, 1, u_n_p_1, 1);
		// Solve system
//		MathLibrary::computeDenseSymSolution(DOF, K, pivot, u_n_p_1);
		LAPACKE_dgetrs(LAPACK_COL_MAJOR, 'N', DOF, 1, K, DOF, pivot.data(), u_n_p_1, DOF);

		u_n_m_1[0] = u_n[0];
		u_n_m_1[1] = u_n[1];
		u_n_m_1[2] = u_n[2];

		u_n[0] = u_n_p_1[0];
		u_n[1] = u_n_p_1[1];
		u_n[2] = u_n_p_1[2];
	}
	std::cout << "We ran #" << numTimeSteps << " time steps!" << std::endl;
	std::cout << u_n_p_1[0] << "\n" << u_n_p_1[1] << "\n" << (double)u_n_p_1[2] << std::scientific << std::endl;

}

void EVAAComputeEngine::computeMKL11DOF(void) {

	typedef double floatEVAA;
	// K stiffness
	floatEVAA k = 10.;
	floatEVAA k_11 = 1.1;
	floatEVAA k_12 = k;
	floatEVAA k_21 = 1.2;
	floatEVAA k_22 = k;
	floatEVAA k_31 = 1.3;
	floatEVAA k_32 = k;
	floatEVAA k_41 = 1.4;
	floatEVAA k_42 = k;
	floatEVAA l_1 = 2;
	floatEVAA l_2 = 1;
	floatEVAA l_3 = 0.8;
	floatEVAA l_4 = 1.25;

	// D - damping matrix
	floatEVAA d_11 = 1.1 * 0.1;
	floatEVAA d_12 = k * 0.1;
	floatEVAA d_21 = 1.2 * 0.1;
	floatEVAA d_22 = k * 0.1;
	floatEVAA d_31 = 1.3 * 0.1;
	floatEVAA d_32 = k * 0.1;
	floatEVAA d_41 = 1.4 * 0.1;
	floatEVAA d_42 = k * 0.1;

	// M - mass matrix
	floatEVAA m_1 = 1e-1;
	floatEVAA m_2 = 2e-1;
	floatEVAA m_3 = 3e-1;
	floatEVAA m_4 = 4e-1;
	floatEVAA m_5 = 5e-1;
	floatEVAA m_6 = 6e-1;
	floatEVAA m_7 = 7e-1;
	floatEVAA m_8 = 8e-1;
	floatEVAA m_9 = 9e-1;
	floatEVAA m_10 = 10e-1;
	floatEVAA m_11 = 11e-1;

	int alignment = 64;
	int DOF = 11;
	int matrixElements = DOF * DOF;

	// allocate matrices of zeros
	floatEVAA* B = (floatEVAA*)mkl_calloc(matrixElements, sizeof(floatEVAA), alignment);
	floatEVAA* M = (floatEVAA*)mkl_calloc(matrixElements, sizeof(floatEVAA), alignment);
	floatEVAA* D = (floatEVAA*)mkl_calloc(matrixElements, sizeof(floatEVAA), alignment);
	floatEVAA* K = (floatEVAA*)mkl_calloc(matrixElements, sizeof(floatEVAA), alignment);

	//std::vector<floatEVAA> B(121);
	//std::vector<floatEVAA> M(121);
	//std::vector<floatEVAA> D(121);
	//std::vector<floatEVAA> K(121);

	floatEVAA* u_n_p_1 = (floatEVAA*)mkl_calloc(DOF, sizeof(floatEVAA), alignment);
	floatEVAA* u_n = (floatEVAA*)mkl_calloc(DOF, sizeof(floatEVAA), alignment);
	floatEVAA* u_n_m_1 = (floatEVAA*)mkl_calloc(DOF, sizeof(floatEVAA), alignment);
	floatEVAA* tmp = (floatEVAA*)mkl_calloc(DOF, sizeof(floatEVAA), alignment);

	//std::vector<floatEVAA> f_n_p_1(11);
	//std::vector<floatEVAA> u_n_p_1(11);
	//std::vector<floatEVAA> u_n(11);
	//std::vector<floatEVAA> u_n_m_1(11);
	//std::vector<floatEVAA> tmp(11);

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
	int nRefinement = 10;
	int numTimeSteps = pow(2, nRefinement);
	//time step size 
	floatEVAA h = 1.0 / (numTimeSteps);
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
	LAPACKE_dgetrf(LAPACK_COL_MAJOR, DOF, DOF, K, DOF, pivot.data());

	// Time loop
	double tmpScalar = (-1.0 / (h * h));
	for (int i = 0; i < DOF; ++i)
		M[(DOF + 1) * i] *= tmpScalar;

		void* jitter;
	// adds matrix matrix product to scalar matrix product
	std::cout << mkl_jit_create_dgemm(&jitter, MKL_COL_MAJOR, MKL_NOTRANS, MKL_NOTRANS, DOF, 1, DOF, 1.0, DOF, DOF, 0.0, DOF) << std::endl;
	dgemm_jit_kernel_t myDGEMMKernel = mkl_jit_get_dgemm_ptr(jitter);
	LAPACKE_set_nancheck(0);



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
		LAPACKE_dgetrs(LAPACK_COL_MAJOR, 'N', DOF, 1, K, DOF, pivot.data(), u_n_p_1, DOF);

		for (int i = 0; i < DOF; ++i) {
			u_n_m_1[i] = u_n[i];
			u_n[i] = u_n_p_1[i];
		}
	}
	std::cout << "We ran #" << numTimeSteps << " time steps!" << std::endl;
	for (int i = 0; i < DOF; ++i) {
		std::cout << u_n_p_1[i] << std::scientific << std::endl;
	}

}


// BLAZE
void EVAAComputeEngine::computeBlazetoy(void) {
#ifdef USE_BLAZE
	using blaze::StaticVector;
	using blaze::DynamicVector;
	using blaze::DynamicMatrix;
	using blaze::rowMajor;
	using blaze::columnVector;

	int DOF = 3;

	DynamicMatrix<double, rowMajor> A(DOF, DOF);  // The system matrix A
	DynamicVector<double, columnVector> b(DOF);   // The right-hand side vector b
	DynamicVector<int, columnVector> ipiv(DOF);   // Pivoting indices
	DynamicMatrix<double, rowMajor> L(DOF, DOF), U(DOF, DOF), P(DOF, DOF);  // The system matrix A
	// ... Initialization
	A(0, 0) = 5;
	A(0, 1) = 2;
	A(0, 2) = 1;
	A(1, 0) = 10;
	A(1, 1) = 9;
	A(1, 2) = 5;
	A(2, 0) = 15;
	A(2, 1) = 26;
	A(2, 2) = 21;

	b[0] = 12;
	b[1] = 43;
	b[2] = 130;

	DynamicMatrix<double, rowMajor>     D(A);  // Temporary matrix to be decomposed
	DynamicVector<double, columnVector> x(b);  // Temporary vector for the solution

	std::cout << "A: \n" << A << " \n";
	std::cout << "b: \n" << b << " \n";
	getrf(D, ipiv.data());
	std::cout << "\tD before: D=LU: \n" << D << " \n"; // A^T=P*L*U
	std::cout << "\tPivot: \n" << ipiv << " \n";
	getrs(D, x, 'T', ipiv.data());
	std::cout << "\tx: x=D^{-1}b: \n" << x << " \n";

	assert(trans(A) * x == b);

	//D = trans(A);
	lu(A, L, U, P); // A = P*L*U
	std::cout << " L = : \n" << L << "\n";
	std::cout << " U = : \n" << U << "\n";
	std::cout << " P = : \n" << P << "\n";


#endif
}

void EVAAComputeEngine::computeBlaze(void) {
#ifdef USE_BLAZE
	using blaze::CompressedMatrix;
	using blaze::DynamicMatrix;
	using blaze::SymmetricMatrix;
	using blaze::DynamicVector;
	using blaze::columnVector;
	using blaze::rowMajor;

	int DOF = 3;

	typedef double floatEVAA;
	// K stiffness
	floatEVAA k_1 = 1.;
	floatEVAA k_2 = 2.;
	floatEVAA k_3 = 3.;

	// Using symmetric matrix enforce the symmetry and is safer
	SymmetricMatrix <DynamicMatrix<floatEVAA>, rowMajor> K(DOF);
	K(0, 0) = k_1 + k_2;
	K(0, 1) = -k_2;
	K(1, 0) = -k_2;
	K(1, 1) = k_2;
	K(2, 2) = k_3;

	// D - damping matrix
	floatEVAA d_1 = 0.1;
	floatEVAA d_2 = 0.5;

	SymmetricMatrix <DynamicMatrix<floatEVAA>, rowMajor> D(DOF);
	D(0, 0) = d_1;
	D(1, 1) = d_2;
	D(1, 2) = -d_2;
	D(2, 1) = -d_2;
	D(2, 2) = d_2;

	// M - mass matrix
	floatEVAA m_1 = 0.1;
	floatEVAA m_2 = 0.2;
	floatEVAA m_3 = 0.3;

	SymmetricMatrix <DynamicMatrix<floatEVAA>, rowMajor> M(DOF);
	M(0, 0) = m_1;
	M(1, 1) = m_2;
	M(2, 2) = m_3;

	// Define the solution vectors
	DynamicVector<floatEVAA, columnVector> u_n_p_1(DOF, (floatEVAA)0.0); // initialize vector of dimension 11 and null elements
	DynamicVector<floatEVAA, columnVector> u_n(DOF, (floatEVAA)0.0); // initialize vector of dimension 11 and null elements
	DynamicVector<floatEVAA, columnVector> u_n_m_1(DOF, (floatEVAA)0.0); // initialize vector of dimension 11 and null elements

	// Perform the iterations
	int nRefinement = 10;
	int numTimeSteps = pow(2, nRefinement);

	//time step size 
	floatEVAA h = 1.0 / ((floatEVAA)numTimeSteps);
	std::cout << "Time step h is: " << h << std::scientific << std::endl;

	// Initial conditions
	const floatEVAA u_init = 1;
	const floatEVAA du_init = 0;
	u_n[0] = u_init;
	u_n_m_1[0] = u_init - h * du_init;

	/// Build dynamic stiffness matrix
	// A = (1.0/(h*h))*M + (1.0/h)*D + K
	DynamicMatrix<floatEVAA, rowMajor> A(DOF, DOF);
	A = 1. / (h * h) * M + 1. / h * D + K;

	///Build rhs for BE integrator
	//B = (2.0 / (h*h))*M + 1.0 / h * D + B
	DynamicMatrix<floatEVAA, rowMajor> B(DOF, DOF);
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
	}
	std::cout << "We ran #" << numTimeSteps << " time steps!" << std::endl;
	std::cout << u_n_p_1 << std::scientific << std::endl;

#endif
}

void EVAAComputeEngine::computeBlaze11DOF(void) {
#ifdef USE_BLAZE
	using blaze::CompressedMatrix;
	using blaze::DynamicMatrix;
	using blaze::SymmetricMatrix;
	using blaze::DynamicVector;
	using blaze::columnVector;
	using blaze::rowMajor;

	typedef double floatEVAA;

	int DOF = 11;

	// K stiffness
	floatEVAA k = 10.;
	floatEVAA k_11 = 1.1;
	floatEVAA k_12 = k;
	floatEVAA k_21 = 1.2;
	floatEVAA k_22 = k;
	floatEVAA k_31 = 1.3;
	floatEVAA k_32 = k;
	floatEVAA k_41 = 1.4;
	floatEVAA k_42 = k;
	floatEVAA l_1 = 2;
	floatEVAA l_2 = 1;
	floatEVAA l_3 = 0.8;
	floatEVAA l_4 = 1.25;

	// Using symmetric matrix enforce the symmetry and is safer
	SymmetricMatrix <DynamicMatrix<floatEVAA>, rowMajor> K(DOF);

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

	// D - damping matrix
	floatEVAA d_11 = 1.1 * 0.1;
	floatEVAA d_12 = k * 0.1;
	floatEVAA d_21 = 1.2 * 0.1;
	floatEVAA d_22 = k * 0.1;
	floatEVAA d_31 = 1.3 * 0.1;
	floatEVAA d_32 = k * 0.1;
	floatEVAA d_41 = 1.4 * 0.1;
	floatEVAA d_42 = k * 0.1;

	SymmetricMatrix <DynamicMatrix<floatEVAA>, rowMajor> D(DOF);
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

	// M - mass matrix
	floatEVAA m_1 = 1e-1;
	floatEVAA m_2 = 2e-1;
	floatEVAA m_3 = 3e-1;
	floatEVAA m_4 = 4e-1;
	floatEVAA m_5 = 5e-1;
	floatEVAA m_6 = 6e-1;
	floatEVAA m_7 = 7e-1;
	floatEVAA m_8 = 8e-1;
	floatEVAA m_9 = 9e-1;
	floatEVAA m_10 = 10e-1;
	floatEVAA m_11 = 11e-1;

	SymmetricMatrix <DynamicMatrix<floatEVAA>, rowMajor> M(DOF);
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
	DynamicVector<floatEVAA, columnVector> u_n_p_1(DOF, (floatEVAA)0.0); // initialize vector of dimension 11 and null elements
	DynamicVector<floatEVAA, columnVector> u_n(DOF, (floatEVAA)0.0); // initialize vector of dimension 11 and null elements
	DynamicVector<floatEVAA, columnVector> u_n_m_1(DOF, (floatEVAA)0.0); // initialize vector of dimension 11 and null elements

	// Perform the iterations
	int nRefinement = 10;
	int numTimeSteps = pow(2, nRefinement);

	//time step size 
	floatEVAA h = 1.0 / ((floatEVAA)numTimeSteps);
	std::cout << "Time step h is: " << h << std::scientific << std::endl;

	// Initial conditions
	//const floatEVAA u_init = 1;
	//const floatEVAA du_init = 0;
	u_n[0] = 1;
	u_n_m_1[0] = 1;

	/// Build dynamic stiffness matrix
	// A = (1.0/(h*h))*M + (1.0/h)*D + K
	DynamicMatrix<floatEVAA, rowMajor> A(DOF, DOF);
	A = 1. / (h * h) * M + 1. / h * D + K;

	///Build rhs for BE integrator
	//B = (2.0 / (h*h))*M + 1.0 / h * D + B
	DynamicMatrix<floatEVAA, rowMajor> B(DOF, DOF);
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
	}
	std::cout << "We ran #" << numTimeSteps << " time steps!" << std::endl;
	std::cout << u_n_p_1 << std::scientific << std::endl;

#endif
}


void EVAAComputeEngine::clean(void) {

}


