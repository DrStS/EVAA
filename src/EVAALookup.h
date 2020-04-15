/***********************************************************************************************//**
* \file EVAALookup.h
* This file holds the class of EVAALookup and EVAAComputeGrid.
* \date 04/14/2020
**************************************************************************************************/
#pragma once
#include <vector>
#include "mkl.h"

/*
* \brief class to generate a grid for interpolation
* 
* this class is only used for testing purposes
*/
class EVAAComputeGrid {
private:
	static double responseFunction(const double a, const double b, const double c, const double length);
public:
	static void buildLinearGrid(double* grid, double* axis, int size, double l_min, double l_max, double* a, double b, double c, int k);
	static void buildChebyshevGrid(double* grid, double* axis, int size, double l_min, double l_max, double* a, double b, double c, int k);
};

/**
* \brief Class in case there are values which are not available in analytical form
*
* generate an instance of this class before the simulation and call the interpolate function during the simulation
*/
class EVAALookup {
private:
	DFTaskPtr* task = nullptr;								/**< Data Fitting task descriptor */
	MKL_INT nx;												/**< number of break points (data points in the grid)*/
	MKL_INT xhint = DF_NON_UNIFORM_PARTITION;				/**< additional info about break points*/
	MKL_INT ny;												/**< number of functions ( in our case 1 per task but 8 in total)*/
	MKL_INT yhint = DF_NO_HINT;								/**< additional info about function*/
	MKL_INT scoeffhint = DF_NO_HINT;						/**< additional info about spline coefficients*/
	MKL_INT bc_type = DF_NO_BC;								/**< boundary conditions type*/
	MKL_INT ic_type = DF_NO_IC;								/**< internal conditions type*/
	MKL_INT rhint = DF_NO_HINT;								/**< interpolation results storage format*/
	double* axis;											/**< array of break points (axis with length values of the grid points)*/
	double* grid;											/**< function values */
	double* ic = nullptr;									/**< internal conditions*/
	double* bc = nullptr;									/**< boundary conditions*/
	double* scoeff = nullptr;								/**< array of spline coefficients*/
	double* datahint = nullptr;								/**< additional info about the structure*/
	MKL_INT stype;											/**< spline type: linear = DF_PP_DEFAULT, spline = DF_PP_NATURAL*/
	MKL_INT sorder;											/**< spline order: linear = DF_PP_LINEAR, spline = DF_PP_CUBIC*/
	int alignment = 64;										/**< alignment for mkl */
public:
	EVAALookup(int size, double* a, double b, double c, double l_min, double l_max, int k, int type, int order);
	~EVAALookup();
	void getInterpolation(double* length, double* inter);
	void getDerivative(double* length, double* deriv);
};