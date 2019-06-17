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
	* \param[in] _vec1 the 1st vector
	* \param[in] _vec2 the 2nd vector
	* \return dot product
	* \author Stefan Sicklinger
	***********/
	double computeDenseDotProduct(const std::vector<double> &_vec1, const std::vector<double> &_vec2);
	/***********************************************************************************************
	* \brief Compute the dot product of two dense vectors
	* \param[in] _vec1 the 1st vector
	* \param[in] _vec2 the 2nd vector
	* \param[in] _nElements number of elements in vec1 (vec2)
	* \return dot product
	* \author Stefan Sicklinger
	***********/
	double computeDenseDotProduct(const double *_vec1, const double *_vec2, const int _nElements);
   /***********************************************************************************************
	* \brief Compute dense symmetrix matrix LU factorisation
	* \param[in] _nElements number of rows = number of columns
	* \param[in] _A matrix 
	* \param[in] _pivot elements
	* \author Stefan Sicklinger
	***********/
	void computeDenseSymLUFactorisation(const int _nElements, std::vector<double> &_A, std::vector<int> &_pivots);
	/***********************************************************************************************
	* \brief Compute backward/forward substitution 
	* \param[in] _nElements number of rows = number of columns
	* \param[in] _A matrix
	 * \param[in] _pivot elements
	 * \author Stefan Sicklinger
	 ***********/
	void computeDenseSymSolution(const int _nElements, std::vector<double> &_A, std::vector<int> &_pivots, std::vector<double> &_rhs);
	/***********************************************************************************************
	* \brief Computes a vector-scalar product and adds the result to a vector. vec2 <- a*vec1 + vec2
	* \param[in] _vec1 the 1st vector
	* \param[in] _vec2 the 2nd vector
	* \param[in] _alpha   scalar
	* \param[in] _nElements number of elements in vec1
	* \author Stefan Sicklinger
	***********/
	void computeDenseVectorAddition(double *vec1, double *vec2, const double _alpha, const int _nElements);
	
   
} /* namespace Math */
