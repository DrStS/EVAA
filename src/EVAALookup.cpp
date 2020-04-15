/***********************************************************************************************//**
* \file EVAALookup.cpp
* This file holds the function definitions of EVAALookup and EVAAComputeGrid.
* \date 04/14/2020
**************************************************************************************************/

#pragma once
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <stdlib.h>
#include <vector>
#ifndef U_Lookup
#define U_Lookup
#include "EVAALookup.h"
#endif
#include "Constants.h"

/**
* function returns a readable output for mkl errors
*/
void CheckDfError(int num);

/*************************************************/
/*			Functions for building a grid        */
/*************************************************/
/**
* \brief Function to feed values to the grid, used for the lookup table
*
* \return double y
*/
double EVAAComputeGrid::responseFunction(
	const double a /**< [in] coefficient for const part */,
	const double b /**< [in] coefficient for linear part*/,
	const double c /**< [in] coefficient for quadratic part */,
	const double length /**< [in] docs for input parameter v. */
) {
	return c * length * length + b * length + a;
}

/**
* \brief build equidistant grid and the corresponding axis
*/
void EVAAComputeGrid::buildLinearGrid(
	double* grid /**< [out] pointer to grid of size size*k to store y values */,
	double* axis /**< [out] pointer to axis of size size to write x value */,
	int size /**< [in] size of one grid */,
	double l_min /**< [in] min length of spring in lookup */,
	double l_max /**< [in] max length of spring in lookup */,
	double* a /**< [in] pointer to array of size k with constant coefficient for each spring */,
	double b /**< [in] coefficient for linear part */,
	double c /**< [in] coefficient for quadratic part */,
	int k /**< [in] number of springs */
) {
	double density = (l_max - l_min) / (size - 1);
	for (auto i = 0; i < size; i++) {
		axis[i] = l_min + i * density;
	}
	for (auto j = 0; j < k; j++) {
		for (auto i = 0; i < size; i++) {
			grid[i + j * size] = responseFunction(a[j], b, c, axis[i]);
		}
	}
}

/**
* \brief build Chebyshev grid and the corresponding axis
*/
void EVAAComputeGrid::buildChebyshevGrid(
	double* grid /**< [out] pointer to grid of size size*k to store y values */,
	double* axis /**< [out] pointer to axis of size size to write x value */,
	int size /**< [in] size of one grid */,
	double l_min /**< [in] min length of spring in lookup */,
	double l_max /**< [in] max length of spring in lookup */,
	double* a /**< [in] pointer to array of size k with constant coefficient for each spring */,
	double b /**< [in] coefficient for linear part */,
	double c /**< [in] coefficient for quadratic part */,
	int k /**< [in] number of springs */
) {
	for (auto i = 0; i < size; i++) {
		axis[i] = (1 + cos((2 * i + 1) / (2 * size) * M_PI)) / 2 * (l_max - l_min) + l_min;
	}
	for (auto j = 0; j < k; j++) {
		for (auto i = 0; i < size; i++) {
			grid[i + j * size] = responseFunction(a[j], b, c, axis[i]);
		}
	}
}


/*************************************************/
/*	Functions for EVAALookup    */
/*************************************************/
/**
* \brief interpolates the ny = k grids ob the lookuptable
* 
* The lookup table has been generated in the initialisation of the object
* this function uses the calculated coefficients to interpolate certain values
*/
void EVAALookup::getInterpolation(
	double* length /**< [in] pointer to array of size k with lenght values of springs*/,
	double* inter /**< [out] pointer to array of size k to store interpolation values*/
) {
	const MKL_INT ndorder = 1;						// size of array describing derivative (dorder), which is definde two lines below
	const MKL_INT dorder[1] = { 1 };				// only the values are computed

	for (auto i = 0; i < ny; i++) {
		dfdInterpolate1D(task[i], DF_INTERP, DF_METHOD_PP,
			1, &length[i], DF_NO_HINT, ndorder,
			dorder, datahint, &inter[i], rhint, 0);
	}


}

/*
* \brief derivative from stifftness k after length
*
* The lookup table has been generated in the initialisation of the object
* this function uses the calculated coefficients to interpolate certain values
*/
void EVAALookup::getDerivative(
	double* length /**< [in] pointer to array of size k with lenght values of springs*/,
	double* deriv /**< [out] pointer to array of size k to store values of the derivative*/
) {
	const MKL_INT ndorder = 2;						// size of array describing derivative (dorder), which is definde two lines below
	const MKL_INT dorder[2] = { 0, 1 };				// only the derivative values are computed
	for (auto i = 0; i < ny; i++) {
		dfdInterpolate1D(task[i], DF_INTERP, DF_METHOD_PP,
			1, &length[i], DF_NO_HINT, ndorder,
			dorder, datahint, &deriv[i], rhint, 0);
	}
}

/**
* \brief Constructor of lookup class
*
*	for linear interpolation:
* 	type = DF_PP_DEFAULT, order = DF_PP_LINEAR
*	for spline interpolation:
*	type = DF_PP_NATURAL, order = DF_PP_CUBIC
*/
EVAALookup::EVAALookup(
	int size /**< [in] size of one grid */,
	double* a /**< [in] pointer to array of size k with constant coefficient for each spring */,
	double b /**< [in] coefficient for linear part of grid function */,
	double c /**< [in] coefficient for quadratic part of grid function */,
	double l_min /**< [in] min length of spring in lookup */,
	double l_max  /**< [in] max length of spring in lookup */,
	int k /**< [in] number of springs */,
	int type /**< [in] int which corersponds to a certain type of interpolation */,
	int order /**< [in] order of the interpolation: depends on type */
) : nx(size), ny(k), sorder(order), stype(type) {
	if (sorder == DF_PP_CUBIC) {
		bc_type = DF_BC_FREE_END;
	}
	// we will get the grid aftwards from a file. That why I do not directly write size, l_min, l_max into the variables
	task = (DFTaskPtr*)mkl_malloc(ny * sizeof(DFTaskPtr), Constants::ALIGNMENT);
	grid = (double*)mkl_calloc(nx * ny, sizeof(double), Constants::ALIGNMENT);
	axis = (double*)mkl_calloc(nx, sizeof(double), Constants::ALIGNMENT);
	scoeff = (double*)mkl_calloc(ny * (nx - 1) * sorder, sizeof(double), Constants::ALIGNMENT);
	
	/* create grid */
	EVAAComputeGrid::buildLinearGrid(grid, axis, nx, l_min, l_max, a, b, c, ny);

	int err = 0;
	for (auto i = 0; i < ny; i++) {
		/***** Create Data Fitting task *****/
		err = dfdNewTask1D(&task[i], nx, axis, xhint, 1, &grid[i * nx], yhint);

		CheckDfError(err);

		/***** Edit task parameters for look up interpolant *****/
		err = dfdEditPPSpline1D(task[i], sorder, stype, bc_type, bc,
			ic_type, ic, &scoeff[i * (nx - 1) * sorder], scoeffhint);

		CheckDfError(err);

		/***** Construct linear spline using STD method *****/
		err = dfdConstruct1D(task[i], DF_PP_SPLINE, DF_METHOD_STD);
		CheckDfError(err);
	}
	mkl_free(grid);
	grid = nullptr;
	mkl_free(axis);
	axis = nullptr;
	delete ic;
	ic = nullptr;
	delete bc;
	bc = nullptr;
}

/**
* \brief Destructor of lookup class
*
* free all the allocated space
*/
EVAALookup::~EVAALookup() {
	/***** Delete Data Fitting task *****/
	for (auto i = 0; i < ny; i++) {
		dfDeleteTask(&task[i]);
	}
	mkl_free(scoeff);
	scoeff = nullptr;
	delete datahint;
	datahint = nullptr;
}


// the following is just to know what exactly the error code of a function means
/*******************************************************************************
* Copyright 2010-2019 Intel Corporation.
*
* This software and the related documents are Intel copyrighted  materials,  and
* your use of  them is  governed by the  express license  under which  they were
* provided to you (License).  Unless the License provides otherwise, you may not
* use, modify, copy, publish, distribute,  disclose or transmit this software or
* the related documents without Intel's prior written permission.
*
* This software and the related documents  are provided as  is,  with no express
* or implied  warranties,  other  than those  that are  expressly stated  in the
* License.
*******************************************************************************/


void CheckDfError(int num)
{
	switch (num)
	{
	case DF_ERROR_NULL_TASK_DESCRIPTOR:
	{
		printf("Error: null task descriptor (code %d).\n", num);
		break;
	}
	case DF_ERROR_MEM_FAILURE:
	{
		printf("Error: memory allocation failure in DF functionality (code %d).\n", num);
		break;
	}
	case DF_ERROR_BAD_NX:
	{
		printf("Error: the number of breakpoints is invalid (code %d).\n", num);
		break;
	}
	case DF_ERROR_BAD_X:
	{
		printf("Error: the array which contains the breakpoints is not defined (code %d).\n", num);
		break;
	}
	case DF_ERROR_BAD_X_HINT:
	{
		printf("Error: invalid flag describing structure of partition (code %d).\n", num);
		break;
	}
	case DF_ERROR_BAD_NY:
	{
		printf("Error: invalid dimension of vector-valued function y (code %d).\n", num);
		break;
	}
	case DF_ERROR_BAD_Y:
	{
		printf("Error: the array which contains function values is invalid (code %d).\n", num);
		break;
	}
	case DF_ERROR_BAD_Y_HINT:
	{
		printf("Error: invalid flag describing structure of function y (code %d).\n", num);
		break;
	}
	case DF_ERROR_BAD_SPLINE_ORDER:
	{
		printf("Error: invalid spline order (code %d).\n", num);
		break;
	}
	case DF_ERROR_BAD_SPLINE_TYPE:
	{
		printf("Error: invalid type of the spline (code %d).\n", num);
		break;
	}
	case DF_ERROR_BAD_IC_TYPE:
	{
		printf("Error: invalid type of internal conditions used in the spline construction (code %d).\n", num);
		break;
	}
	case DF_ERROR_BAD_IC:
	{
		printf("Error: array of internal conditions for spline construction is not defined (code %d).\n", num);
		break;
	}
	case DF_ERROR_BAD_BC_TYPE:
	{
		printf("Error: invalid type of boundary conditions used in the spline construction (code %d).\n", num);
		break;
	}
	case DF_ERROR_BAD_BC:
	{
		printf("Error: array which presents boundary conditions for spline construction is not defined (code %d).\n", num);
		break;
	}
	case DF_ERROR_BAD_PP_COEFF:
	{
		printf("Error: array of piece-wise polynomial spline coefficients is not defined (code %d).\n", num);
		break;
	}
	case DF_ERROR_BAD_PP_COEFF_HINT:
	{
		printf("Error: invalid flag describing structure of the piece-wise polynomial spline coefficients (code %d).\n", num);
		break;
	}
	case DF_ERROR_BAD_PERIODIC_VAL:
	{
		printf("Error: function values at the end points of the interpolation interval are not equal as required in periodic boundary conditions (code %d).\n", num);
		break;
	}
	case DF_ERROR_BAD_DATA_ATTR:
	{
		printf("Error: invalid attribute of the pointer to be set or modified in Data Fitting task descriptor with EditIdxPtr editor (code %d).\n", num);
		break;
	}
	case DF_ERROR_BAD_DATA_IDX:
	{
		printf("Error: index of pointer to be set or modified in Data Fitting task descriptor with EditIdxPtr editor is out of range (code %d).\n", num);
		break;
	}
	case DF_ERROR_BAD_NSITE:
	{
		printf("Error: invalid number of interpolation sites (code %d).\n", num);
		break;
	}
	case DF_ERROR_BAD_SITE:
	{
		printf("Error: array of interpolation sites is not defined (code %d).\n", num);
		break;
	}
	case DF_ERROR_BAD_SITE_HINT:
	{
		printf("Error: invalid flag describing structure of interpolation sites (code %d).\n", num);
		break;
	}
	case DF_ERROR_BAD_NDORDER:
	{
		printf("Error: invalid size of array that defines order of the derivatives to be computed at the interpolation sites (code %d).\n", num);
		break;
	}
	case DF_ERROR_BAD_DORDER:
	{
		printf("Error: array defining derivative orders to be computed at interpolation sites is not defined (code %d).\n", num);
		break;
	}
	case DF_ERROR_BAD_DATA_HINT:
	{
		printf("Error: invalid flag providing a-priori information about partition and/or interpolation sites (code %d).\n", num);
		break;
	}
	case DF_ERROR_BAD_INTERP:
	{
		printf("Error: array of spline based interpolation results is not defined (code %d).\n", num);
		break;
	}
	case DF_ERROR_BAD_INTERP_HINT:
	{
		printf("Error: invalid flag defining structure of spline based interpolation results (code %d).\n", num);
		break;
	}
	case DF_ERROR_BAD_CELL_IDX:
	{
		printf("Error: array of indices of partition cells containing interpolation sites is not defined (code %d).\n", num);
		break;
	}
	case DF_ERROR_BAD_NLIM:
	{
		printf("Error: invalid size of arrays containing integration limits (code %d).\n", num);
		break;
	}
	case DF_ERROR_BAD_LLIM:
	{
		printf("Error: array of left integration limits is not defined (code %d).\n", num);
		break;
	}
	case DF_ERROR_BAD_RLIM:
	{
		printf("Error: array of right integration limits is not defined (code %d).\n", num);
		break;
	}
	case DF_ERROR_BAD_INTEGR:
	{
		printf("Error: array of spline based integration results is not defined (code %d).\n", num);
		break;
	}
	case DF_ERROR_BAD_INTEGR_HINT:
	{
		printf("Error: invalid flag defining structure of spline based integration results (code %d).\n", num);
		break;
	}
	case DF_ERROR_BAD_LOOKUP_INTERP_SITE:
	{
		printf("Error: bad site provided for interpolation with look-up interpolator (code %d).\n", num);
		break;
	}
	case DF_ERROR_NULL_PTR:
	{
		printf("Error: bad pointer provided in DF function (code %d).\n", num);
		break;
	}
	default: break;
	}

	if (num < 0) {
		exit(1);
	}
}
