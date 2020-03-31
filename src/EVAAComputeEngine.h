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
* \file EVAAComputeEngine.h
* This file holds the class of ComputeEngine.
* \date 6/13/2019
**************************************************************************************************/
#pragma once
#include <string>
/********//**
* \brief Class ComputeEngine the core of STACCATO
***********/
class EVAAComputeEngine {

public:
	/***********************************************************************************************
	* \brief Constructor
	* \param[in] file name of xml file to call singelton constructor of metadatabase
	* \author Stefan Sicklinger
	***********/
	EVAAComputeEngine(std::string _xmlFileName);
	/***********************************************************************************************
	* \brief Destructor
	* \author Stefan Sicklinger
	***********/
	~EVAAComputeEngine();
	/***********************************************************************************************
	* \brief prepare compute engine
	* \author Stefan Sicklinger
	***********/
	void prepare(void);
	/***********************************************************************************************
	* \brief compute engine
	* \author Stefan Sicklinger
	***********/
	// EIGEN
	void computeEigen11DOF(void);
	/***********************************************************************************************
	* \brief compute engine
	* \author Stefan Sicklinger
	***********/
	// BLAZE
	void computeBlaze11DOF(void);
	/***********************************************************************************************
	* \brief compute engine
	* \author Stefan Sicklinger
	***********/
	// MKL
	void computeMKL11DOF(void);
	/***********************************************************************************************
	* \brief compute engine
	* \author Stefan Sicklinger
	***********/
	// MKL
	void computeMKLNasa(void);
	void computeMKLNasa_example(void);
	/***********************************************************************************************
	* \brief clean compute engine free memory
	* \author Stefan Sicklinger
	***********/
	void clean(void);
private:
};

//#ifndef ComputeNasa
//#define ComputeNasa
//template <typename T>
//class ComputeNasa {
//private:
//	// define constants
//	const int dim = 3;
//	const int alignment = 64;
//	const int num_wheels = 4;
//	const int num_wheels_x_dim = num_wheels * dim;
//	const int DOF_diag = 9; // the diagonal elements from A
//	const int DOF = 11;
//	const int matrixElements = dim * dim;
//	// const int matrixElements = DOF * DOF;
//	T FC;
//	T* r, r_tilda, FW, FT, FR, lower_spring_length, upper_spring_length, lower_spring_stiffness, upper_spring_stiffness;
//	// The system matrix A. 0:3 - 2x2 upper-left symmetric block; 4:12 . We will store its inverse on diagonal blocks
//	// b vector, 11 DOF (from the 12 components we prohibit rotation along y-axis, b[2]=0 is ignored)
//	// x - the solution vector
//	T* A, rhs, x; 
//	T* Ic, Tc, basis_car, C_Nc, r_global, car_corners, ones4;
//	T* upper_length, lower_length, diff_vector, upper_force, lower_force, upper_F;
//	T* C_Nc_extended, sum_torque, w_tilda, Qc;
//
//public:
//	ComputeNasa(); // constructor without parameters (default, testing reasons); initialization
//	
//	ComputeNasa(T* r, T* r_tilda, T* FW, T* FT, T* FR, T* lower_spring_length, T* upper_spring_length, T* lower_spring_stiffness, T* upper_spring_stiffness, T* A, T FC, T* IC) :
//		r(r), r_tilda(r_tilda), FW(FW), FT(FT), FR(FR), lower_spring_length(lower_spring_length), upper_spring_length(upper_spring_length), lower_spring_stiffness(lower_spring_stiffness), upper_spring_stiffness(upper_spring_stiffness), A(A), FC(FC), IC(IC) {
//		// memset(ones4, 1, num_wheels);
//	}
//
//	~ComputeNasa();
//
//	void compute_f3D_reduced(T time, T* f);
//};
//#endif