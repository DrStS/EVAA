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
#define MKL_DIRECT_CALL 1
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

	double k_1 = 1;
	double k_2 = 2;
	double k_3 = 3;
	std::vector<double> K(9);
	std::vector<double> rhs(3);
	
	K[0] = k_1 + k_2;
	K[1] = -k_2;
	K[2] = 0.0;
	K[3] = -k_2;
	K[4] = k_2;
	K[5] = 0.0;
	K[6] = 0.0;
	K[7] = 0.0;
	K[8] = k_3;
	rhs[0] = 1.0;
	rhs[1] = 1.0;
	rhs[2] = 1.0;

	std::vector<int> pivot(3);
	// LU Decomposition
	MathLibrary::computeDenseSymLUFactorisation(3, K, pivot);
	// Solve system
	MathLibrary::computeDenseSymSolution(3, K, pivot, rhs);

	std::vector<double> resultReference = {2.,5./2.,1./3.};
	double absError = fabs(resultReference[0] - rhs[0]) + fabs(resultReference[1] - rhs[1]) + fabs(resultReference[2] - rhs[2]);
	std::cout << "Done! Absolut error is: " << absError << std::scientific << std::endl;
}


void EVAAComputeEngine::clean(void) {

}


