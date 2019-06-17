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
#include "Mathlibrary.h"
#ifdef USE_INTEL_MKL
#include <mkl.h>
#endif

EVAAComputeEngine::EVAAComputeEngine(std::string _xmlFileName){
	// Intialize XML metadatabase singelton 
}


EVAAComputeEngine::~EVAAComputeEngine(){
}

void EVAAComputeEngine::prepare(void) {
   
}

void EVAAComputeEngine::compute(void) {

	double k_1 = 1.;
	double k_2 = 2.;
	double k_3 = 3.;
	double d_1 = 1. / 10.;
	double d_2 = 1. / 2.;
	double m_1 = 1e-1;
	double m_2 = 2e-1;
	double m_3 = 3e-1;
	std::vector<double> B(9);
	std::vector<double> M(9);
	std::vector<double> D(9);
	std::vector<double> K(9);

	std::vector<double> f_n_p_1(3);
	std::vector<double> u_n_p_1(3);
	std::vector<double> u_n(3);
	std::vector<double> u_n_m_1(3);
	std::vector<double> tmp(3);
	
	M[0] = m_1;
	M[1] = 0.0;
	M[2] = 0.0;
	M[3] = 0.0;
	M[4] = m_2;
	M[5] = 0.0;
	M[6] = 0.0;
	M[7] = 0.0;
	M[8] = m_3;

	D[0] =  d_1;
	D[1] =  0.0;
	D[2] =  0.0;
	D[3] =  0.0;
	D[4] =  d_2;
	D[5] = -d_2;
	D[6] =  0.0;
	D[7] = -d_2;
	D[8] =  d_2;

	K[0] =  k_1 + k_2;
	K[1] = -k_2;
	K[2] =  0.0;
	K[3] = -k_2;
	K[4] =  k_2;
	K[5] =  0.0;
	K[6] =  0.0;
	K[7] =  0.0;
	K[8] =  k_3;

	//Initial conditions
	u_n[0] = 1.;
	u_n[1] = 0.;
	u_n[2] = 0.;
	//
	u_n_m_1[0] = 1.;
	u_n_m_1[1] = 0.;
	u_n_m_1[2] = 0.;


	int nRefinement = 20;
	int numTimeSteps=pow(2, nRefinement);
	//time step size 
	double h = 1.0 / (numTimeSteps);
	std::cout << "Time step h is: " << h << std::scientific << std::endl;
	/// Build dynamic stiffness matrix
	// K' <- (1.0/(h*h))*M + K
	MathLibrary::computeDenseVectorAddition(M.data(), K.data(), (1.0 / (h*h)), 9);
	// K <- (1.0/h)*D + K'
	MathLibrary::computeDenseVectorAddition(D.data(), K.data(), (1.0 / h), 9);
	/// K holds now dynamic stiffness matrix  for BE integrator

	///Build rhs for BE integrator
	//B' <-(2.0 / (h*h))*M + B
	MathLibrary::computeDenseVectorAddition(M.data(), B.data(), (2.0 / (h*h)), 9);
	//B <-(1.0 / (h))*D + B'
	MathLibrary::computeDenseVectorAddition(D.data(), B.data(), (1.0 / h), 9);
	//A*u_n_p_1=B*u_n+C*u_n_m_1+f_n_p_1 <== BE
    
	std::vector<int> pivot(3);
	// LU Decomposition
	MathLibrary::computeDenseSymLUFactorisation(3, K, pivot);

	std::vector<double> timeVec;
	timeVec.resize(numTimeSteps);
	///Time loop
	for (int iTime = 0; iTime < numTimeSteps; iTime++) {
		timeVec[iTime] = iTime * h;
		// y: = alpha * A*x + beta * y
		// tmp1 = B*u_n
		cblas_dgemv(CblasColMajor, CblasNoTrans, 3, 3, 1.0, B.data(), 3, u_n.data(), 1, 0.0, u_n_p_1.data(), 1);
		// y: = alpha * A*x + beta * y
		// tmp = ((-1.0 / (h*h))*M) * u_n_m_1
		cblas_dgemv(CblasColMajor, CblasNoTrans, 3, 3, (-1.0 / (h*h)), M.data(), 3, u_n_m_1.data(), 1, 0.0, tmp.data(), 1);
		//rhs = tmp1 + tmp2
		// tmp1 <- 1.0 tmp2 + tmp1
		MathLibrary::computeDenseVectorAddition(tmp.data(), u_n_p_1.data(), 1.0, 3);
		//tmp1 <- rhs

		// Solve system
		MathLibrary::computeDenseSymSolution(3, K, pivot, u_n_p_1);


		u_n_m_1[0] = u_n[0];
		u_n_m_1[1] = u_n[1];
		u_n_m_1[2] = u_n[2];

		u_n[0] = u_n_p_1[0];
		u_n[1] = u_n_p_1[1];
		u_n[2] = u_n_p_1[2];


	}

	std::cout << u_n_p_1[0] << " " << u_n_p_1[1] << " " << u_n_p_1[2] << std::endl;

}


void EVAAComputeEngine::clean(void) {

}


