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
#ifdef USE_INTEL_MKL
#define USE_INTEL_MKL_BLAS
#endif

#include "MathLibrary.h"

namespace MathLibrary {

double computeDenseDotProduct(const double *vec1, const double *vec2, const int elements) {
#ifdef USE_INTEL_MKL
	mkl_set_num_threads(EVAA::AuxiliaryParameters::denseVectorMatrixThreads);
	return cblas_ddot(elements, vec1, 1, vec2, 1);
#endif
#ifndef USE_INTEL_MKL
	return 0;
#endif
}

double computeDenseDotProduct(const std::vector<double> &vec1, const std::vector<double> &vec2) {
#ifdef USE_INTEL_MKL
	return cblas_ddot(vec1.size(), vec1.data(), 1, vec2.data(), 1);
#endif
#ifndef USE_INTEL_MKL
	return 0;
#endif
}

void computeDenseSymLUFactorisation(const int nElements, std::vector<double> &A, std::vector<int> &pivots){
#ifdef USE_INTEL_MKL 
   LAPACKE_dgetrf(LAPACK_COL_MAJOR, nElements, nElements, A.data(), nElements, pivots.data());
#endif
}

void computeDenseSymSolution(const int nElements, std::vector<double> &A, std::vector<int> &pivots, std::vector<double> &rhs){
#ifdef USE_INTEL_MKL 
	LAPACKE_dgetrs(LAPACK_COL_MAJOR, 'N', nElements, 1, A.data(), nElements, pivots.data(), rhs.data(), nElements);
#endif
}

void computeDenseVectorAddition(double *vec1, double *vec2, const double a, const int elements) {
#ifdef USE_INTEL_MKL
	mkl_set_num_threads(EVAA::AuxiliaryParameters::denseVectorMatrixThreads);
	cblas_daxpy(elements, a, vec1, 1, vec2, 1);
#endif
}


} /* namespace Math */


	