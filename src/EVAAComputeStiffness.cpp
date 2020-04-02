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
#include <vector>
#include "EVAAComputeStiffness.h"

#define ALIGNMENT 64

/*************************************************/
/*			Functions for building a grid        */
/*************************************************/
double EVAAComputeGrid::responseFunction(const double a, const double b, const double c, const double length) {
	return c * length * length + b * length + a;
}

void EVAAComputeGrid::buildLinearGrid(double* grid, double* axis, int size, double l_min, double l_max, double a, double b, double c, int k) {
	double density = (l_max - l_min) / size;
#pragma simd
	for (auto i = 0; i < size; i++) {
		axis[i] = l_min + i * density;
	}
#pragma simd
	for (auto j = 0; j < k; j++) {
		for (auto i = 0; i < size; i++) {
			grid[i] = responseFunction(a, b, c * j, axis[i]);
		}
	}
}

void EVAAComputeGrid::buildChebyshevGrid(double* grid, double* axis, int size, double l_min, double l_max, double a, double b, double c, int k) {
#pragma simd
	for (auto i = 0; i < size; i++) {
		axis[i] = (1 + cos((2 * i + 1) / (2 * size) * M_PI)) / 2 * (l_max - l_min) + l_min;
	}
#pragma simd
	for (auto j = 0; j < k; j++) {
		for (auto i = 0; i < size; i++) {
			grid[i + j * size] = responseFunction(a, b, c * j, axis[i]);
		}
	}
}


/*************************************************/
/*	Functions for EVAAComputeStiffness_Linear    */
/*************************************************/
void EVAAComputeStiffness_Linear::getStiffness(double* stiffness, double* x) {
	/*
	int ind_low = floor((length - _l_min));
	if (length > _l_max) { // over max
		stiffness = _grid[_size - 1];
	}
	else {
		if (length < _l_min) { // below min
			stiffness = grid[0];
		}
		else { // in boundaries
			double val_low = grid[ind_low];
			double val_up = grid[ceil((length - _l_min) / _density)];
			double weight = (length - _l_min - _density * ind_low) / _density;
			stiffness = (1 - weight) * val_low + weight * val_up;
		}
	}*/
	/***** Interpolate using lookup method *****/
	dfdInterpolate1D(task, DF_INTERP, DF_METHOD_PP,
		nsite, site, sitehint, ndorder,
		dorder, datahint, result, rhint, 0);
}

/*
double EVAAComputeStiffness_Linear::getDerivative(double length) {
	double val_low, val_up, val_der;
	if (length > _l_max) { // over max
		val_low = _grid[_grid.size() - 2];
		val_up = _grid[_grid.size() - 1];
		val_der =  (val_up - val_low) / _density;
	}
	else {
		if (length < _l_min) { // below min
			val_low = _grid[0];
			val_up = _grid[1];
			val_der = (val_up - val_low) / _density;
		}
		else { // in boundaries
			val_low = _grid[floor((length - _l_min) / _density)];
			val_up = _grid[ceil((length - _l_min) / _density)];
			val_der = (val_up - val_low) / _density;
		}
	}
	return val_der;
}*/

EVAAComputeStiffness_Linear::EVAAComputeStiffness_Linear(MKL_INT size, double a, double b, double c, double l_min, double l_max, MKL_INT k) {
	// we will get the grid aftwards from a file. That why I do not directly write size, l_min, l_max into the variables
	nx = size;
	ny = k;
	nsite = nx;
	grid = (double*)mkl_calloc(nx * ny, sizeof(double), ALIGNMENT);
	axis = (double*)mkl_calloc(nx, sizeof(double), ALIGNMENT);
	site = (double*)mkl_calloc(nx, sizeof(double), ALIGNMENT);
	result = (double*)mkl_calloc(nx, sizeof(double), ALIGNMENT);
	EVAAComputeGrid::buildLinearGrid(grid, axis, nx, l_min, l_max, a, b, c, ny);

	/***** Generate interpolation sites *****/
	for (auto i = 0; i < nx; i++)
	{
		site[i] = axis[i];
	}

	/***** Create Data Fitting task *****/
	dfdNewTask1D(&task, nx, axis, xhint, ny, grid, yhint);
	
	/***** Edit task parameters for look up interpolant *****/
	dfdEditPPSpline1D(task, sorder, stype, 0, 0, 0, 0, 0, 0);
	
}

EVAAComputeStiffness_Linear::~EVAAComputeStiffness_Linear() {
	/***** Delete Data Fitting task *****/
	dfDeleteTask(&task);

	int errnums = 0;
	/***** Check interpolation results *****/
	for (auto j = 0; j < nx; j++)
	{
		if (abs(result[j] - grid[j]) > std::numeric_limits<double>::epsilon())
			errnums++;
	}

	/***** Print summary of the test *****/
	if (errnums != 0)
	{
		printf("\n\nError: Computed interpolation results");
		printf(" are incorrect\n");
	}
	else
	{
		printf("\n\nComputed interpolation results");
		printf(" are correct\n");
	}

	/* free used space */
	MKL_free(grid);
	MKL_free(axis);
	MKL_free(site);
	MKL_free(result);
}
