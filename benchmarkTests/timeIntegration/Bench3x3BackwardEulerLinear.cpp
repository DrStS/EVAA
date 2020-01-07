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
#include <benchmark/benchmark.h>
#include <mkl.h>
#include <iostream>

#ifdef USE_INTEL_MKL
#include <mkl.h>
#define USE_GEMM
#endif
#ifdef USE_EIGEN
#include <Eigen/Dense>
using namespace Eigen;
#endif
#ifdef USE_BLAZE
#include <blaze/Math.h>
#endif


#ifdef USE_INTEL_MKL
static void BM_MKL3x3BackwardEulerLinear(benchmark::State& state) {
	for (auto _ : state) {
		double k_1 = 1.;
		double k_2 = 2.;
		double k_3 = 3.;
		double d_1 = 1. / 10.;
		double d_2 = 1. / 2.;
		double m_1 = 1e-1;
		double m_2 = 2e-1;
		double m_3 = 3e-1;

		int alignment = 64;

		double* B = (double*)mkl_malloc(sizeof(double) * 9, alignment);
		double* M = (double*)mkl_malloc(sizeof(double) * 9, alignment);
		double* D = (double*)mkl_malloc(sizeof(double) * 9, alignment);
		double* K = (double*)mkl_malloc(sizeof(double) * 9, alignment);

		//std::vector<double> B(9);
		//std::vector<double> M(9);
		//std::vector<double> D(9);
		//std::vector<double> K(9);

		double* u_n_p_1 = (double*)mkl_malloc(sizeof(double) * 3, alignment);
		double* u_n = (double*)mkl_malloc(sizeof(double) * 3, alignment);
		double* u_n_m_1 = (double*)mkl_malloc(sizeof(double) * 3, alignment);
		double* tmp = (double*)mkl_malloc(sizeof(double) * 3, alignment);

		//std::vector<double> f_n_p_1(3);
		//std::vector<double> u_n_p_1(3);
		//std::vector<double> u_n(3);
		//std::vector<double> u_n_m_1(3);
		//std::vector<double> tmp(3);

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
		double h = 1.0 / (numTimeSteps);
		/// Build dynamic stiffness matrix
		// K' <- (1.0/(h*h))*M + K
	//	MathLibrary::computeDenseVectorAddition(M.data(), K.data(), (1.0 / (h*h)), 9);
		cblas_daxpy(9, (1.0 / (h * h)), M, 1, K, 1);
		// K <- (1.0/h)*D + K'
	//	MathLibrary::computeDenseVectorAddition(D.data(), K.data(), (1.0 / h), 9);
		cblas_daxpy(9, (1.0 / h), D, 1, K, 1);
		/// K holds now dynamic stiffness matrix  for BE integrator

		///Build rhs for BE integrator
		//B' <-(2.0 / (h*h))*M + B
	//	MathLibrary::computeDenseVectorAddition(M.data(), B.data(), (2.0 / (h*h)), 9);
		cblas_daxpy(9, (2.0 / (h * h)), M, 1, B, 1);
		//B <-(1.0 / (h))*D + B'
	//	MathLibrary::computeDenseVectorAddition(D.data(), B.data(), (1.0 / h), 9);
		cblas_daxpy(9, (1.0 / h), D, 1, B, 1);
		//A*u_n_p_1=B*u_n+C*u_n_m_1+f_n_p_1 <== BE

		std::vector<int> pivot(3);
		// LU Decomposition
	//	MathLibrary::computeDenseSymLUFactorisation(3, K, pivot);
		LAPACKE_dgetrf(LAPACK_COL_MAJOR, 3, 3, K, 3, pivot.data());

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
		mkl_jit_create_dgemm(&jitter, MKL_COL_MAJOR, MKL_NOTRANS, MKL_NOTRANS, 3, 1, 3, 1.0, 3, 3, 0.0, 3) ;
		dgemm_jit_kernel_t myDGEMMKernel = mkl_jit_get_dgemm_ptr(jitter);
		LAPACKE_set_nancheck(0);
		for (int iTime = 0; iTime < numTimeSteps; iTime++) {
			//timeVec[iTime] = iTime * h;
			// y: = alpha * A*x + beta * y
			// u_n_p_1 = B*u_n
#ifndef USE_GEMM
			cblas_dgemv(CblasColMajor, CblasNoTrans, 3, 3, 1.0, B, 3, u_n, 1, 0.0, u_n_p_1, 1);
#endif
#ifdef USE_GEMM 
			//cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 3, 1, 3, 1.0, B, 3, u_n, 3, 0.0, u_n_p_1, 3);
			myDGEMMKernel(jitter, B, u_n, u_n_p_1);
#endif
			// y: = alpha*A*x + beta*y
			// tmp = ((-1.0 / (h*h))*M) * u_n_m_1
#ifndef USE_GEMM
			cblas_dgemv(CblasColMajor, CblasNoTrans, 3, 3, (-1.0 / (h * h)), M, 3, u_n_m_1, 1, 0.0, tmp, 1);
#endif
#ifdef USE_GEMM 
			//cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 3, 1, 3, 1.0, M, 3, u_n_m_1, 3, 0.0, tmp, 3);
			myDGEMMKernel(jitter, M, u_n_m_1, tmp);
#endif
			// u_n_p_1 <- 1.0 tmp + u_n_p_1
	//		MathLibrary::computeDenseVectorAddition(tmp.data(), u_n_p_1.data(), 1.0, 3);
			cblas_daxpy(3, 1.0, tmp, 1, u_n_p_1, 1);
			// Solve system
	//		MathLibrary::computeDenseSymSolution(3, K, pivot, u_n_p_1);
			LAPACKE_dgetrs(LAPACK_COL_MAJOR, 'N', 3, 1, K, 3, pivot.data(), u_n_p_1, 3);

			u_n_m_1[0] = u_n[0];
			u_n_m_1[1] = u_n[1];
			u_n_m_1[2] = u_n[2];

			u_n[0] = u_n_p_1[0];
			u_n[1] = u_n_p_1[1];
			u_n[2] = u_n_p_1[2];
		}
		//std::cout << "We ran #" << numTimeSteps << " time steps!" << std::endl;
		//std::cout << u_n_p_1[0] << " " << u_n_p_1[1] << " " << (double)u_n_p_1[2] << std::scientific << std::endl;
	}
}
BENCHMARK(BM_MKL3x3BackwardEulerLinear)->UseRealTime();
#endif

#ifdef USE_EIGEN
static void BM_Eigen3x3BackwardEulerLinear(benchmark::State& state) {
	for (auto _ : state) {
		double k_1 = 1.;
		double k_2 = 2.;
		double k_3 = 3.;
		double d_1 = 1. / 10.;
		double d_2 = 1. / 2.;
		double m_1 = 1e-1;
		double m_2 = 2e-1;
		double m_3 = 3e-1;

		MatrixXd B(3, 3);
		MatrixXd M(3, 3);
		MatrixXd D(3, 3);
		MatrixXd K(3, 3);

		VectorXd u_n_p_1(3);
		VectorXd u_n(3);
		VectorXd u_n_m_1(3);

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
		double h = 1.0 / (numTimeSteps);
		//std::cout << "Time step h is: " << h << std::scientific << std::endl;

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
		//std::cout << "We ran #" << numTimeSteps << " time steps!" << std::endl;
		//std::cout << u_n_p_1 << std::scientific << std::endl;
	}
}
BENCHMARK(BM_Eigen3x3BackwardEulerLinear)->UseRealTime();
#endif

#ifdef USE_BLAZE
static void BM_Blaze3x3BackwardEulerLinear(benchmark::State& state) {
	for (auto _ : state) {
		using blaze::CompressedMatrix;
		using blaze::DynamicMatrix;
		using blaze::SymmetricMatrix;
		using blaze::DynamicVector;
		using blaze::columnVector;
		using blaze::rowMajor;

		// K stiffness
		double k_1 = 1.;
		double k_2 = 2.;
		double k_3 = 3.;

		// Using symmetric matrix enforce the symmetry and is safer
		SymmetricMatrix <DynamicMatrix<double>, rowMajor> K(3);
		K(0, 0) = k_1 + k_2;
		K(0, 1) = -k_2;
		K(1, 0) = -k_2;
		K(1, 1) = k_2;
		K(2, 2) = k_3;

		// D - damping matrix
		double d_1 = 0.1;
		double d_2 = 0.5;

		SymmetricMatrix <DynamicMatrix<double>, rowMajor> D(3);
		D(0, 0) = d_1;
		D(1, 1) = d_2;
		D(1, 2) = -d_2;
		D(2, 1) = -d_2;
		D(2, 2) = d_2;

		// M - mass matrix
		double m_1 = 0.1;
		double m_2 = 0.2;
		double m_3 = 0.3;

		SymmetricMatrix <DynamicMatrix<double>, rowMajor> M(3);
		M(0, 0) = m_1;
		M(1, 1) = m_2;
		M(2, 2) = m_3;

		// Define the solution vectors
		DynamicVector<double, columnVector> u_n_p_1(3, (double)0.0); // initialize vector of dimension 11 and null elements
		DynamicVector<double, columnVector> u_n(3, (double)0.0); // initialize vector of dimension 11 and null elements
		DynamicVector<double, columnVector> u_n_m_1(3, (double)0.0); // initialize vector of dimension 11 and null elements

		// Perform the iterations
		int nRefinement = 10;
		int numTimeSteps = pow(2, nRefinement);

		//time step size 
		double h = 1.0 / ((double)numTimeSteps);
		std::cout << "Time step h is: " << h << std::scientific << std::endl;

		// Initial conditions
		const double u_init = 1;
		const double du_init = 0;
		u_n[0] = u_init;
		u_n_m_1[0] = u_init - h * du_init;

		/// Build dynamic stiffness matrix
		// A = (1.0/(h*h))*M + (1.0/h)*D + K
		DynamicMatrix<double, rowMajor> A(3, 3);
		A = 1. / (h * h) * M + 1. / h * D + K;

		///Build rhs for BE integrator
		//B = (2.0 / (h*h))*M + 1.0 / h * D + B
		DynamicMatrix<double, rowMajor> B(3, 3);
		B = 2. / (h * h) * M + 1. / h * D;

		// LU Decomposition
		DynamicVector<int, columnVector> ipiv(3);   // Pivoting indices
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
	}
}
BENCHMARK(BM_Eigen3x3BackwardEulerLinear)->UseRealTime();
#endif

