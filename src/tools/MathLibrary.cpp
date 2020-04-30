/*
 * Copyright &copy; 2019, Dr. Stefan Sicklinger, Munich \n
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

#include "MathLibrary.h"

#include <mkl.h>
#include <mkl_cblas.h>

#include <iostream>

namespace EVAA {
namespace MathLibrary {

void printMKLInfo()
{
#ifdef USE_INTEL_MKL
    int len = 198;
    char buf[198];
    mkl_get_version_string(buf, len);
    std::cout << buf << std::endl;
#endif
}

void computeDenseSymLUFactorisation(const int nElements, std::vector<double> &A,
                                    std::vector<int> &pivots)
{
#ifdef USE_INTEL_MKL
    LAPACKE_dgetrf(LAPACK_COL_MAJOR, nElements, nElements, A.data(), nElements, pivots.data());
#endif
}

void computeDenseSymSolution(const int nElements, std::vector<double> &A, std::vector<int> &pivots,
                             std::vector<double> &rhs)
{
#ifdef USE_INTEL_MKL
    mkl_set_num_threads(EVAA::AuxiliaryParameters::denseVectorMatrixThreads);
    LAPACKE_dgetrs(LAPACK_COL_MAJOR, 'N', nElements, 1, A.data(), nElements, pivots.data(),
                   rhs.data(), nElements);
#endif
}

void computeDenseVectorAddition(double *vec1, double *vec2, const double a, const int elements)
{
#ifdef USE_INTEL_MKL
    mkl_set_num_threads(EVAA::AuxiliaryParameters::denseVectorMatrixThreads);
    cblas_daxpy(elements, a, vec1, 1, vec2, 1);
#endif
}

}  // namespace MathLibrary
}  // namespace EVAA
