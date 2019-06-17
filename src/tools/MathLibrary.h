/*  Copyright &copy; 2019, Stefan Sicklinger, Munich
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
/*************************************************************************************************
* \file MathLibrary.h
* This file holds the math funcktion of EVAA.
* \date 6/13/2019
**************************************************************************************************/
#pragma once

#include "AuxiliaryParameters.h"
#include <vector>

#ifdef USE_INTEL_MKL
#define MKL_DIRECT_CALL 1
#include <mkl.h>
#endif

namespace MathLibrary {
	/***********************************************************************************************
	* \brief Compute the dot product of two dense vectors
	* \param[in] vec1 the 1st vector
	* \param[in] vec2 the 2nd vector
	* \return dot product
	* \author Stefan Sicklinger
	***********/
	double computeDenseDotProduct(const std::vector<double> &vec1, const std::vector<double> &vec2);
	/***********************************************************************************************
	* \brief Compute the dot product of two dense vectors
	* \param[in] vec1 the 1st vector
	* \param[in] vec2 the 2nd vector
	* \param[in] elements number of elements in vec1 (vec2)
	* \return dot product
	* \author Stefan Sicklinger
	***********/
	double computeDenseDotProduct(const double *vec1, const double *vec2, const int elements);
   /***********************************************************************************************
	* \brief Compute dense symmetrix matrix LU factorisation
	* \param[in] nElements number of rows = number of columns
	* \param[in] A matrix 
	* \param[in] pivot elements
	* \author Stefan Sicklinger
	***********/
	void computeDenseSymLUFactorisation(const int nElements, std::vector<double> &A, std::vector<int> &pivots);
	/***********************************************************************************************
	* \brief Compute backward/forward substitution 
	* \param[in] nElements number of rows = number of columns
	* \param[in] A matrix
	 * \param[in] pivot elements
	 * \author Stefan Sicklinger
	 ***********/
	void computeDenseSymSolution(const int nElements, std::vector<double> &A, std::vector<int> &pivots, std::vector<double> &rhs);
   
} /* namespace Math */
