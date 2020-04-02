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
* \file EVAAComputeStiffness.h
* This file holds the class of ComputeStiffness.
* \date 03/18/2020
**************************************************************************************************/

#include <vector>
#include "mkl.h"

/*
* this class is only used for testing
*/
class EVAAComputeGrid {
private:
	/*
	* \brief function which the grid point values correspond to
	*/
	static double responseFunction(const double a, const double b, const double c, const double length);
public:
	/*
	* \brief for linear interpolation we need a linear grid
	*/
	static void buildLinearGrid(double* grid, double* axis, int size, double l_min, double l_max, double a, double b, double c, int k);
	/*
	* \brief for spline interpolation we need a Chebyshev grid
	*/
	static void buildChebyshevGrid(double* grid, double* axis, int size, double l_min, double l_max, double a, double b, double c, int k);
};

class EVAAComputeStiffness {
private:
	DFTaskPtr* task;									// Data Fitting task descriptor
	MKL_INT nx;											// number of break points
	MKL_INT xhint = DF_NON_UNIFORM_PARTITION;			// additional info about break points
	MKL_INT ny;											// number of functions
	MKL_INT yhint = DF_NO_HINT;							// additional info about function
	MKL_INT scoeffhint = DF_NO_HINT;					// additional info about spline coefficients
	MKL_INT bc_type = DF_NO_BC;							// boundary conditions type
	MKL_INT ic_type = DF_NO_IC;							// internal conditions type
	MKL_INT rhint = DF_NO_HINT;							// interpolation results storage format
	double* axis;										// array of break points
	double* grid;										// function values
	double* ic = 0;										// internal conditions
	double* bc = 0;										// boundary conditions
	double* scoeff;										// array of spline coefficients
	double* datahint = 0;								// additional info about the structure
	MKL_INT stype;										// linear = DF_PP_DEFAULT, spline = DF_PP_NATURAL
	MKL_INT sorder;										// linear = DF_PP_LINEAR, spline = DF_PP_CUBIC
public:
	/*
	* \brief Constructor
	*/
	EVAAComputeStiffness(int size, double a, double b, double c, double l_min, double l_max, int k, int type, int order);
	/*
	* \brief free all the allocated storage
	*/
	~EVAAComputeStiffness();
	/*
	* \brief get stiffness value via linear interpolation
	*/
	void getStiffness(double* length, double* stiffness); // or double* x and calculate length
	/*
	* \brief get the derivative ob k after l
	*/
	void getDerivative(double* length, double* deriv);
};