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
#pragma once
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <stdlib.h>
#include <vector>
#ifndef U_COMPSTIFF
#define U_COMPSTIFF
#include "EVAAComputeStiffness.h"
#endif

void CheckDfError(int num);

/*************************************************/
/*			Functions for building a grid        */
/*************************************************/
double EVAAComputeGrid::responseFunction(const double a, const double b, const double c, const double length) {
	return c * length * length + b * length + a;
	// return a;
}

/*
* build equidistant grid and the corresponding axis
*/
void EVAAComputeGrid::buildLinearGrid(double* grid, double* axis, int size, double l_min, double l_max, double* a, double b, double c, int k) {
	double density = (l_max - l_min) / (size-1);
	for (auto i = 0; i < size; i++) {
		axis[i] = l_min + i * density;
	}
	for (auto j = 0; j < k; j++) {
		for (auto i = 0; i < size; i++) {
			grid[i + j * size] = responseFunction(a[j], b, c, axis[i]);
		}
	}
}

/*
* build Chebyshev grid and the corresponding axis
*/
void EVAAComputeGrid::buildChebyshevGrid(double* grid, double* axis, int size, double l_min, double l_max, double* a, double b, double c, int k) {
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
/*	Functions for EVAAComputeStiffness_Linear    */
/*************************************************/
void EVAAComputeStiffness::getStiffness(double* length, double* stiffness) {
	int err = 0;
	MKL_INT ndorder = 1;                    // size of array describing derivative (dorder), which is definde two lines below
	MKL_INT dorder[1] = { 1 }; // only the values are computed
	for (auto i = 0; i < ny; i++) {
		err = dfdInterpolate1D(task[i], DF_INTERP, DF_METHOD_PP,
			1, &length[i], DF_NO_HINT, ndorder,
			dorder, datahint, &stiffness[i], rhint, 0);

		CheckDfError(err);
		CheckDfError(err);
	}
}

/*
* \brief derivative from stifftness k after length
*/
void EVAAComputeStiffness::getDerivative(double* length, double* deriv) {
	int err = 0;
	MKL_INT ndorder = 2;                    // size of array describing derivative (dorder), which is definde two lines below
	MKL_INT dorder[2] = { 0, 1 }; // only the derivative values are computed
	for (auto i = 0; i < ny; i++) {
		err = dfdInterpolate1D(task[i], DF_INTERP, DF_METHOD_PP,
			1, &length[i], DF_NO_HINT, ndorder,
			dorder, datahint, &deriv[i], rhint, 0);
		CheckDfError(err);
	}
}
/*
*	for linear interpolation:
* 	type = DF_PP_DEFAULT, order = DF_PP_LINEAR
*
*	for spline interpolation:
*	type = DF_PP_NATURAL, order = DF_PP_CUBIC
*/
EVAAComputeStiffness::EVAAComputeStiffness(int size, double* a, double b, double c, double l_min, double l_max, int k, int type, int order) : nx(size), ny(k), sorder(order), stype(type) {
	if (sorder == DF_PP_CUBIC) {
		bc_type = DF_BC_FREE_END;
	}
	// we will get the grid aftwards from a file. That why I do not directly write size, l_min, l_max into the variables
	task = (DFTaskPtr*)mkl_malloc(ny * sizeof(DFTaskPtr), alignment);
	grid = (double*)mkl_calloc(nx * ny, sizeof(double), alignment);
	axis = (double*)mkl_calloc(nx, sizeof(double), alignment);
	scoeff = (double*)mkl_calloc(ny * (nx - 1) * sorder, sizeof(double), alignment);
	
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

EVAAComputeStiffness::~EVAAComputeStiffness() {
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
