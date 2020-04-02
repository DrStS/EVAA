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

class EVAAComputeStiffness_Linear {
private:
	DFTaskPtr task;                     // Data Fitting task descriptor
	MKL_INT nx;                         // number of break points
	MKL_INT xhint = DF_NON_UNIFORM_PARTITION;                      // additional info about break points
	MKL_INT ny;                         // number of functions
	MKL_INT yhint = DF_NO_HINT;                      // additional info about function

	double* axis;                        // array of break points
	double* grid;                    // function values

	MKL_INT nsite;                      // total number of interpolation sites
	MKL_INT sitehint = DF_NON_UNIFORM_PARTITION;                   // additional info about interpolation
	double* site; ;                 // array of interpolation sites

	MKL_INT ndorder = 1;                    // size of array describing derivative
										// orders
	MKL_INT dorder[1] = { 1 };           // only value to calculate
										// will be computed

	MKL_INT rhint = DF_MATRIX_STORAGE_ROWS;                      // interpolation results storage format
	double* result;                    // spline evaluation results

	MKL_INT stype = DF_LOOKUP_INTERPOLANT;
	MKL_INT sorder = DF_PP_STD;

	double* datahint = 0;                   // additional info about structure
										// of arrays x and y
public:
	/*
	* \brief Constructor
	*/
	EVAAComputeStiffness_Linear(int size, double a, double b, double c, double l_min, double l_max, int k);
	/*
	* \brief free all the allocated storage
	*/
	~EVAAComputeStiffness_Linear();
	/*
	* \brief get stiffness value via linear interpolation
	*/
	void getStiffness(double* stiffness, double* x); // or double* x and calculate length
	/*
	* \brief get derivative of linear interpolation
	*/
	double getDerivative(double* x);  // or double* x and calculate length
};


class EVAAComputeStiffness_Spline {
private:
	/*
	* \brief boundary value for spring length (later for all params)
	*/
	const double _l_min, _l_max;
	/*
	* \brief function evaluations of all param combinations
	*/
	double* grid;
	/*
	* \brief _l_min to _l_max in _density steps
	*/
	double* axis;
public:
	/*
	* \brief Constructor
	*/
	EVAAComputeStiffness_Spline(double size, double a, double b, double c, double l_max, double l_min);
	/*
	* \brief get stiffness value via linear interpolation
	*/
	void getStiffness(double* stiffness, double* x); // or double* x and calculate length
	/*
	* \brief get derivative of linear interpolation
	*/
	double getDerivative(double* x);  // or double* x and calculate length
};