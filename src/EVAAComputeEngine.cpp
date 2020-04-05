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
#include <iostream>
#include <vector>
#include "EVAAComputeEngine.h"
#include "MathLibrary.h"
#include <limits>
#include <fstream>
#include "car.h"
#include "11DOF.h"
#ifdef USE_INTEL_MKL
#include <mkl.h>
#define USE_GEMM
#endif

#ifdef USE_EIGEN
#include <Eigen/Dense>
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Matrix3d;
using Eigen::IOFormat;
#endif

#ifdef USE_BLAZE
#include <blaze/Math.h>
#endif

typedef std::numeric_limits< double > dbl;
typedef double floatEVAA;

EVAAComputeEngine::EVAAComputeEngine(std::string xmlFileName) {
	// Intialize XML metadatabase singelton
	_xmlFileName = xmlFileName;
	std::ifstream f(xmlFileName.c_str());
	if (f.good()) {
		std::cout << "Read general simulation parameters and car input data at " << _xmlFileName << std::endl;
	}
	else {
		std::cout << "XML file at " << _xmlFileName << " does not exist!" << std::endl;
		exit(2);
	}
	ReadXML reader(_xmlFileName);
	reader.ReadParameters(_parameters);
}

EVAAComputeEngine::EVAAComputeEngine(std::string xmlFileName, std::string loadxml) {
	// Intialize XML metadatabase singelton
	_xmlFileName = xmlFileName;
	_xmlLoadFileName = loadxml;

	std::ifstream f(xmlFileName.c_str());
	if (f.good()) {
		std::cout << "Read general simulation parameters and car input data at " << _xmlFileName << std::endl;
	}
	else {
		std::cout << "XML file at " << _xmlFileName << " for general and car settings does not exist!" << std::endl;
		exit(2);
	}

	std::ifstream ff(loadxml.c_str());
	if (ff.good()) {
		std::cout << "Read load parameters in file: " << _xmlLoadFileName << std::endl;
	}
	else {
		std::cout << "XML file at " << _xmlLoadFileName << " for load parameters does not exist!" << std::endl;
		exit(2);
	}

	ReadXML reader(_xmlFileName, _xmlLoadFileName);
	reader.ReadParameters(_parameters);
	reader.ReadLoadParameters(_load_module_parameter);
	reader.ReadLookupParameters(&lookupStiffness, _parameters);
	std::cout << "interpolation = " << _parameters.interpolation << std::endl;
}


EVAAComputeEngine::~EVAAComputeEngine() {
	if (_parameters.interpolation) {
		delete lookupStiffness;
	}
}

void EVAAComputeEngine::prepare(void) {
	MathLibrary::printMKLInfo();
}

void EVAAComputeEngine::computeEigen11DOF(void) {

	ReadXML reader(_xmlLoadFileName);

	int DOF = _parameters.DOF;

	// K stiffness
	double k_11 = _parameters.k_body[0];
	double k_12 = _parameters.k_tyre[0];
	double k_21 = _parameters.k_body[1];
	double k_22 = _parameters.k_tyre[1];
	double k_31 = _parameters.k_body[2];
	double k_32 = _parameters.k_tyre[2];
	double k_41 = _parameters.k_body[3];
	double k_42 = _parameters.k_tyre[3];
	double l_1 = _parameters.l_lat[0];
	double l_2 = _parameters.l_lat[2];
	double l_3 = _parameters.l_long[0];
	double l_4 = _parameters.l_long[2];


	// D - damping matrix
	double d_11 = _parameters.c_body[0];
	double d_12 = _parameters.c_tyre[0];
	double d_21 = _parameters.c_body[1];
	double d_22 = _parameters.c_tyre[1];
	double d_31 = _parameters.c_body[2];
	double d_32 = _parameters.c_tyre[2];
	double d_41 = _parameters.c_body[3];
	double d_42 = _parameters.c_tyre[3];
	double d_l1 = 0.1*_parameters.l_lat[0];
	double d_l2 = 0.1*_parameters.l_lat[2];
	double d_l3 = 0.1*_parameters.l_long[0];
	double d_l4 = 0.1*_parameters.l_long[2];

	// M - mass matrix
	double m_1 = _parameters.mass_body;
	double m_2 = _parameters.I_body[0];
	double m_3 = _parameters.I_body[2];
	double m_4 = _parameters.mass_wheel[2];
	double m_5 = _parameters.mass_tyre[2];
	double m_6 = _parameters.mass_wheel[3];
	double m_7 = _parameters.mass_tyre[3];
	double m_8 = _parameters.mass_wheel[1];
	double m_9 = _parameters.mass_tyre[1];
	double m_10 = _parameters.mass_wheel[0];
	double m_11 = _parameters.mass_tyre[0];

	MatrixXd A(DOF, DOF);
	MatrixXd B(DOF, DOF);
	MatrixXd M(DOF, DOF);
	MatrixXd D(DOF, DOF);
	MatrixXd K(DOF, DOF);
	MatrixXd K_aux(DOF, DOF);
	MatrixXd D_aux(DOF, DOF);

	VectorXd u_n_p_1(DOF);
	VectorXd u_n(DOF);
	VectorXd u_n_m_1(DOF);

	// Define the mass matrix
	M = MatrixXd::Zero(DOF, DOF);
	M(0, 0) = m_1;
	M(1, 1) = m_2;
	M(2, 2) = m_3;
	M(3, 3) = m_4;
	M(4, 4) = m_5;
	M(5, 5) = m_6;
	M(6, 6) = m_7;
	M(7, 7) = m_8;
	M(8, 8) = m_9;
	M(9, 9) = m_10;
	M(10, 10) = m_11;

	// Define the stiffness matrix
	K << k_11 + k_21 + k_31 + k_41, -k_11 * l_1 - k_41 * l_1 + k_21 * l_2 + k_31 * l_2, -k_11 * l_3 - k_21 * l_3 + k_31 * l_4 + k_41 * l_4, -k_11, 0.0, -k_21, 0.0, -k_31, 0.0, -k_41, 0.0,
		0.0, l_1* l_1* k_11 + l_2 * l_2 * k_21 + l_2 * l_2 * k_31 + l_1 * l_1 * k_41, l_1* l_3* k_11 - l_3 * l_2 * k_21 + l_2 * l_4 * k_31 - l_1 * l_4 * k_41, l_1* k_11, 0.0, -l_2 * k_21, 0.0, -l_2 * k_31, 0.0, l_1* k_41, 0.0,
		0.0, 0.0, l_3* l_3* k_11 + l_3 * l_3 * k_21 + l_4 * l_4 * k_31 + l_4 * l_4 * k_41, l_3* k_11, 0.0, l_3* k_21, 0.0, -l_4 * k_31, 0.0, -l_4 * k_41, 0.0,
		0.0, 0.0, 0.0, k_11 + k_12, -k_12, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
		0.0, 0.0, 0.0, 0.0, k_12, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
		0.0, 0.0, 0.0, 0.0, 0.0, k_21 + k_22, -k_22, 0.0, 0.0, 0.0, 0.0,
		0.0, 0.0, 0.0, 0.0, 0.0, 0.0, k_22, 0.0, 0.0, 0.0, 0.0,
		0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, k_31 + k_32, -k_32, 0.0, 0.0,
		0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, k_32, 0.0, 0.0,
		0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, k_41 + k_42, -k_42,
		0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, k_42;
	MatrixXd K_diag(DOF, DOF);
	K_aux = K + K.transpose();
	K_aux.diagonal() -= K.diagonal();
	K = K_aux;

	// Define the damping matrix
	D << d_11 + d_21 + d_31 + d_41, -d_11 * l_1 - d_41 * l_1 + d_21 * l_2 + d_31 * l_2, -d_11 * l_3 - d_21 * l_3 + d_31 * l_4 + d_41 * l_4, -d_11, 0.0, -d_21, 0.0, -d_31, 0.0, -d_41, 0.0,
		0.0, l_1* l_1* d_11 + l_2 * l_2 * d_21 + l_2 * l_2 * d_31 + l_1 * l_1 * d_41, l_1* l_3* d_11 - l_3 * l_2 * d_21 + l_2 * l_4 * d_31 - l_1 * l_4 * d_41, l_1* d_11, 0.0, -l_2 * d_21, 0.0, -l_2 * d_31, 0.0, l_1* d_41, 0.0,
		0.0, 0.0, l_3* l_3* d_11 + l_3 * l_3 * d_21 + l_4 * l_4 * d_31 + l_4 * l_4 * d_41, l_3* d_11, 0.0, l_3* d_21, 0.0, -l_4 * d_31, 0.0, -l_4 * d_41, 0.0,
		0.0, 0.0, 0.0, d_11 + d_12, -d_12, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
		0.0, 0.0, 0.0, 0.0, d_12, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
		0.0, 0.0, 0.0, 0.0, 0.0, d_21 + d_22, -d_22, 0.0, 0.0, 0.0, 0.0,
		0.0, 0.0, 0.0, 0.0, 0.0, 0.0, d_22, 0.0, 0.0, 0.0, 0.0,
		0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, d_31 + d_32, -d_32, 0.0, 0.0,
		0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, d_32, 0.0, 0.0,
		0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, d_41 + d_42, -d_42,
		0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, d_42;
	MatrixXd D_diag(DOF, DOF);
	D_aux = D + D.transpose();
	D_aux.diagonal() -= D.diagonal();
	D = D_aux;
	A = MatrixXd::Zero(DOF, DOF);
	B = MatrixXd::Zero(DOF, DOF);

	u_n_p_1 = VectorXd::Zero(DOF);
	u_n = VectorXd::Zero(DOF);
	u_n(0) = 1;
	u_n_m_1 = u_n;

	//	int nRefinement = 10;
	//	int numTimeSteps = pow(2, nRefinement);
	int numTimeSteps = _parameters.num_time_iter;
	//time step size 
	//double h =  / (numTimeSteps);
	double h = _parameters.timestep;
	std::cout << "Time step h is: " << h << std::scientific << std::endl;

	/*IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
	std::cout << "M: \n\n" << M.format(CleanFmt) << std::endl << std::endl;
	std::cout << "D: \n\n" << D.format(CleanFmt) << std::endl << std::endl;
	std::cout << "K: \n\n" << K.format(CleanFmt) << std::endl << std::endl;*/

	/// Build dynamic stiffness matrix
	// A = (1.0/(h*h))*M + (1.0/h)*D + K
	A = (1.0 / (h * h)) * M + 1.0 / h * D + K;

	///Build rhs for BE integrator
	//B = (2.0 / (h*h))*M + 1.0 / h * D + B
	B = 2.0 / (h * h) * M + 1.0 / h * D;

	/*std::cout << "A: \n\n" << A.format(CleanFmt) << std::endl << std::endl;
	std::cout << "B: \n\n" << B.format(CleanFmt) << std::endl << std::endl;*/
	// LU Decomposition
	Eigen::PartialPivLU<MatrixXd> lu(A);

	M *= (-1.0 / (h * h));
	for (int iTime = 0; iTime < numTimeSteps; iTime++) {
		// B*u_n + (-1.0 / (h*h))*M
		// Solve system
		u_n_p_1 = lu.solve(B * u_n + M * u_n_m_1);

		u_n_m_1 = u_n;
		u_n = u_n_p_1;
		//reader.ReadLoadParameters(_load_module_parameter);
	}
	std::cout << "We ran #" << numTimeSteps << " time steps!" << std::endl;
	std::cout << u_n_p_1 << std::scientific << std::endl;
}

void EVAAComputeEngine::computeMKL11DOF(void) {

	ReadXML reader(_xmlLoadFileName);

	int DOF = _parameters.DOF;

	// K stiffness
	double k_11 = _parameters.k_body[0];
	double k_12 = _parameters.k_tyre[0];
	double k_21 = _parameters.k_body[1];
	double k_22 = _parameters.k_tyre[1];
	double k_31 = _parameters.k_body[2];
	double k_32 = _parameters.k_tyre[2];
	double k_41 = _parameters.k_body[3];
	double k_42 = _parameters.k_tyre[3];
	double l_1 = _parameters.l_lat[0];
	double l_2 = _parameters.l_lat[2];
	double l_3 = _parameters.l_long[0];
	double l_4 = _parameters.l_long[2];


	// D - damping matrix
	double d_11 = _parameters.c_body[0];
	double d_12 = _parameters.c_tyre[0];
	double d_21 = _parameters.c_body[1];
	double d_22 = _parameters.c_tyre[1];
	double d_31 = _parameters.c_body[2];
	double d_32 = _parameters.c_tyre[2];
	double d_41 = _parameters.c_body[3];
	double d_42 = _parameters.c_tyre[3];
	double d_l1 = 0.1 * _parameters.l_lat[0];
	double d_l2 = 0.1 * _parameters.l_lat[2];
	double d_l3 = 0.1 * _parameters.l_long[0];
	double d_l4 = 0.1 * _parameters.l_long[2];

	// M - mass matrix
	double m_1 = _parameters.mass_body;
	double m_2 = _parameters.I_body[0];
	double m_3 = _parameters.I_body[2];
	double m_4 = _parameters.mass_wheel[2];
	double m_5 = _parameters.mass_tyre[2];
	double m_6 = _parameters.mass_wheel[3];
	double m_7 = _parameters.mass_tyre[3];
	double m_8 = _parameters.mass_wheel[1];
	double m_9 = _parameters.mass_tyre[1];
	double m_10 = _parameters.mass_wheel[0];
	double m_11 = _parameters.mass_tyre[0];

	int alignment = 64;
	int matrixElements = DOF * DOF;

	// allocate matrices of zeros
	double* B = (double*)mkl_calloc(matrixElements, sizeof(double), alignment);
	double* M = (double*)mkl_calloc(matrixElements, sizeof(double), alignment);
	double* D = (double*)mkl_calloc(matrixElements, sizeof(double), alignment);
	double* K = (double*)mkl_calloc(matrixElements, sizeof(double), alignment);

	//std::vector<double> B(121);
	//std::vector<double> M(121);
	//std::vector<double> D(121);
	//std::vector<double> K(121);

	double* u_n_p_1 = (double*)mkl_calloc(DOF, sizeof(double), alignment);
	double* u_n = (double*)mkl_calloc(DOF, sizeof(double), alignment);
	double* u_n_m_1 = (double*)mkl_calloc(DOF, sizeof(double), alignment);
	double* tmp = (double*)mkl_calloc(DOF, sizeof(double), alignment);


	//std::vector<double> f_n_p_1(11);
	//std::vector<double> u_n_p_1(11);
	//std::vector<double> u_n(11);
	//std::vector<double> u_n_m_1(11);
	//std::vector<double> tmp(11);

	// Mass matrix initialization
	M[0] = m_1;
	M[12] = m_2;
	M[24] = m_3;
	M[36] = m_4;
	M[48] = m_5;
	M[60] = m_6;
	M[72] = m_7;
	M[84] = m_8;
	M[96] = m_9;
	M[108] = m_10;
	M[120] = m_11;

	/*std::cout << "\nM:\n";
	for (int i = 0; i < 11; ++i) {
		for (int j = 0; j < 11; ++j) {
			std::cout << M[i * 11 + j] << "\t";
		}
		std::cout << std::endl;
	}*/

	// Stiffness matrix
	K[0] = k_11 + k_21 + k_31 + k_41;
	K[1] = -k_11 * l_1 - k_41 * l_1 + k_21 * l_2 + k_31 * l_2;
	K[2] = -k_11 * l_3 - k_21 * l_3 + k_31 * l_4 + k_41 * l_4;
	K[3] = -k_11;
	K[5] = -k_21;
	K[7] = -k_31;
	K[9] = -k_41;
	K[12] = l_1 * l_1 * k_11 + l_2 * l_2 * k_21 + l_2 * l_2 * k_31 + l_1 * l_1 * k_41;
	K[13] = l_1 * l_3 * k_11 - l_3 * l_2 * k_21 + l_2 * l_4 * k_31 - l_1 * l_4 * k_41;
	K[14] = l_1 * k_11;
	K[16] = -l_2 * k_21;
	K[18] = -l_2 * k_31;
	K[20] = l_1 * k_41;
	K[24] = l_3 * l_3 * k_11 + l_3 * l_3 * k_21 + l_4 * l_4 * k_31 + l_4 * l_4 * k_41;
	K[25] = l_3 * k_11;
	K[27] = l_3 * k_21;
	K[29] = -l_4 * k_31;
	K[31] = -l_4 * k_41;
	K[36] = k_11 + k_12;
	K[37] = -k_12;
	K[48] = k_12;
	K[60] = k_21 + k_22;
	K[61] = -k_22;
	K[72] = k_22;
	K[84] = k_31 + k_32;
	K[85] = -k_32;
	K[96] = k_32;
	K[108] = k_41 + k_42;
	K[109] = -k_42;
	K[120] = k_42;

	// Symmetry
	for (int i = 1; i < DOF; ++i)
		for (int j = 0; j < i; ++j)
			K[i * DOF + j] = K[j * DOF + i];


	/*std::cout << "\nK:\n";
	for (int i = 0; i < 11; ++i) {
		for (int j = 0; j < 11; ++j) {
			std::cout << K[i * 11 + j] << "\t";
		}
		std::cout << std::endl;
	}*/

	// Damping matrix
	D[0] = d_11 + d_21 + d_31 + d_41;
	D[1] = -d_11 * l_1 - d_41 * l_1 + d_21 * l_2 + d_31 * l_2;
	D[2] = -d_11 * l_3 - d_21 * l_3 + d_31 * l_4 + d_41 * l_4;
	D[3] = -d_11;
	D[5] = -d_21;
	D[7] = -d_31;
	D[9] = -d_41;
	D[12] = l_1 * l_1 * d_11 + l_2 * l_2 * d_21 + l_2 * l_2 * d_31 + l_1 * l_1 * d_41;
	D[13] = l_1 * l_3 * d_11 - l_3 * l_2 * d_21 + l_2 * l_4 * d_31 - l_1 * l_4 * d_41;
	D[14] = l_1 * d_11;
	D[16] = -l_2 * d_21;
	D[18] = -l_2 * d_31;
	D[20] = l_1 * d_41;
	D[24] = l_3 * l_3 * d_11 + l_3 * l_3 * d_21 + l_4 * l_4 * d_31 + l_4 * l_4 * d_41;
	D[25] = l_3 * d_11;
	D[27] = l_3 * d_21;
	D[29] = -l_4 * d_31;
	D[31] = -l_4 * d_41;
	D[36] = d_11 + d_12;
	D[37] = -d_12;
	D[48] = d_12;
	D[60] = d_21 + d_22;
	D[61] = -d_22;
	D[72] = d_22;
	D[84] = d_31 + d_32;
	D[85] = -d_32;
	D[96] = d_32;
	D[108] = d_41 + d_42;
	D[109] = -d_42;
	D[120] = d_42;

	// Symmetry
	for (int i = 1; i < DOF; ++i)
		for (int j = 0; j < i; ++j)
			D[i * DOF + j] = D[j * DOF + i];

	/*std::cout << "\nD:\n";
	for (int i = 0; i < 11; ++i) {
		for (int j = 0; j < 11; ++j) {
			std::cout << D[i * 11 + j] << "\t";
		}
		std::cout << std::endl;
	}*/

	//Initial conditions
	u_n[0] = 1.;
	u_n_m_1[0] = 1.;
	/*std::cout << "u_n:\n";
	for (int i = 0; i < 11; ++i) {
		std::cout << u_n[i] << ", ";
	}*/


	mkl_set_num_threads(8);
	//	int nRefinement = 10;
	//	int numTimeSteps = pow(2, nRefinement);
	int numTimeSteps = _parameters.num_time_iter;
	//time step size 
	//double h =  / (numTimeSteps);
	double h = _parameters.timestep;
	std::cout << "Time step h is: " << h << std::scientific << std::endl;
	/// Build dynamic stiffness matrix
	// gets written in rear vector
	// K' <- (1.0/(h*h))*M + K
//	MathLibrary::computeDenseVectorAddition(M.data(), K.data(), (1.0 / (h*h)), 9);
	cblas_daxpy(121, (1.0 / (h * h)), M, 1, K, 1);
	// K <- (1.0/h)*D + K'
//	MathLibrary::computeDenseVectorAddition(D.data(), K.data(), (1.0 / h), 121);
	cblas_daxpy(121, (1.0 / h), D, 1, K, 1);
	/// K holds now dynamic stiffness matrix  for BE integrator
	///Build rhs for BE integrator
	//B' <-(2.0 / (h*h))*M + B
//	MathLibrary::computeDenseVectorAddition(M.data(), B.data(), (2.0 / (h*h)), 121);
	cblas_daxpy(matrixElements, (2.0 / (h * h)), M, 1, B, 1);

	//B <-(1.0 / (h))*D + B'
//	MathLibrary::computeDenseVectorAddition(D.data(), B.data(), (1.0 / h), 121);
	cblas_daxpy(matrixElements, (1.0 / h), D, 1, B, 1);
	//A*u_n_p_1=B*u_n+C*u_n_m_1+f_n_p_1 <== BE	

	std::vector<int> pivot(DOF);
	// LU Decomposition
//	MathLibrary::computeDenseSymLUFactorisation(11, K, pivot);
	LAPACKE_dgetrf(LAPACK_ROW_MAJOR, DOF, DOF, K, DOF, pivot.data());

	// Time loop
	double tmpScalar = (-1.0 / (h * h));
	for (int i = 0; i < DOF; ++i)
		M[(DOF + 1) * i] *= tmpScalar;

	void* jitter;
	// adds matrix matrix product to scalar matrix product
	std::cout << mkl_jit_create_dgemm(&jitter, MKL_COL_MAJOR, MKL_NOTRANS, MKL_NOTRANS, DOF, 1, DOF, 1.0, DOF, DOF, 0.0, DOF) << std::endl;
	dgemm_jit_kernel_t myDGEMMKernel = mkl_jit_get_dgemm_ptr(jitter);

	for (int iTime = 0; iTime < numTimeSteps; iTime++) {
		//timeVec[iTime] = iTime * h;
		// y: = alpha * A*x + beta * y
		// u_n_p_1 = B*u_n
#ifndef USE_GEMM
/* void cblas_dgemv(const CBLAS_LAYOUT Layout, const CBLAS_TRANSPOSE trans, const MKL_INT m,
	const MKL_INT n, const double alpha, const double* a, const MKL_INT lda, const double* x,
	const MKL_INT incx, const double beta, double* y, const MKL_INT incy); */
		cblas_dgemv(CblasColMajor, CblasNoTrans, 11, 11, 1.0, B, 11, u_n, 1, 0.0, u_n_p_1, 1);
#endif
#ifdef USE_GEMM 
		//cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 11, 1, 11, 1.0, B, 11, u_n, 11, 0.0, u_n_p_1, 11);
		myDGEMMKernel(jitter, B, u_n, u_n_p_1);
#endif
		// y: = alpha*A*x + beta*y
		// tmp = ((-1.0 / (h*h))*M) * u_n_m_1
#ifndef USE_GEMM
		cblas_dgemv(CblasColMajor, CblasNoTrans, DOF, DOF, (-1.0 / (h * h)), M, DOF, u_n_m_1, 1, 0.0, tmp, 1);
#endif
#ifdef USE_GEMM 
		//cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 11, 1, 11, 1.0, M, 11, u_n_m_1, 11, 0.0, tmp, 11);
		myDGEMMKernel(jitter, M, u_n_m_1, tmp);
#endif
		// u_n_p_1 <- 1.0 tmp + u_n_p_1
//		MathLibrary::computeDenseVectorAddition(tmp.data(), u_n_p_1.data(), 1.0, 11);
		//void cblas_daxpy (const MKL_INT n, const double a, const double *x, const MKL_INT incx, double *y, const MKL_INT incy);
		cblas_daxpy(DOF, 1.0, tmp, 1, u_n_p_1, 1);
		// Solve system
//		MathLibrary::computeDenseSymSolution(11, K, pivot, u_n_p_1);
// lapack_int LAPACKE_dgetrs (int matrix_layout , char trans , lapack_int n , lapack_int nrhs , const double * a , lapack_int lda , const lapack_int * ipiv , double * b , lapack_int ldb );
		LAPACKE_dgetrs(LAPACK_ROW_MAJOR, 'N', DOF, 1, K, DOF, pivot.data(), u_n_p_1, 1);

		for (int i = 0; i < DOF; ++i) {
			u_n_m_1[i] = u_n[i];
			u_n[i] = u_n_p_1[i];
		}
		//reader.ReadLoadParameters(_load_module_parameter);
	}
	std::cout << "We ran #" << numTimeSteps << " time steps!" << std::endl;
	for (int i = 0; i < DOF; ++i) {
		std::cout << u_n_p_1[i] << std::scientific << std::endl;
	}

	// free the memory
	mkl_free(B);
	mkl_free(M);
	mkl_free(D);
	mkl_free(K);

	mkl_free(u_n_p_1);
	mkl_free(u_n);
	mkl_free(u_n_m_1);
	mkl_free(tmp);
}


void EVAAComputeEngine::computeMKLlinear11dof() {
	if (_load_module_parameter.boundary_condition_road == NONFIXED) {
		floatEVAA h = _parameters.timestep;
		floatEVAA tend = _parameters.num_time_iter * h;
		int DOF = _parameters.DOF;
		const int alignment = 64;
		floatEVAA* soln = (floatEVAA*)mkl_calloc(DOF, sizeof(floatEVAA), alignment);
		Car<floatEVAA>* Car1 = new Car<floatEVAA>(_parameters, lookupStiffness);
		linear11dof_full<floatEVAA> solver(_parameters, _load_module_parameter, Car1);
		solver.apply_boundary_condition(_load_module_parameter.boundary_condition_road);
		solver.solve(soln);
		size_t steps = floor(tend / h);
		std::cout << "Solution after " << steps << " timesteps, f =" << std::endl;
		for (auto i = 0; i < DOF; ++i) {
			std::cout << soln[i] << std::endl;
		}
		mkl_free(soln);
		/*floatEVAA* vect1 = (floatEVAA*)mkl_calloc(27, sizeof(floatEVAA), alignment);
		floatEVAA* vect2 = (floatEVAA*)mkl_calloc(11, sizeof(floatEVAA), alignment);
		floatEVAA* vect3 = (floatEVAA*)mkl_calloc(18, sizeof(floatEVAA), alignment);
		for (int i = 0; i < 27; ++i) {
			vect2[i] = rand();
			vect3[i] = rand()/2.0;
		}
		Car1->populate_results(vect3, vect2, vect1);
		MathLibrary::write_vector(vect1, 27);
		MathLibrary::write_vector(vect2, 11);
		MathLibrary::write_vector(vect3, 18);
		exit(5);*/
		delete Car1;
	}
	else {
		std::cout << "Linear11dof will only work with NONFIXED, computation skipped" << std::endl;
	}
}


void EVAAComputeEngine::computeBlaze11DOF(void) {
#ifdef USE_BLAZE

	ReadXML reader(_xmlLoadFileName);

	using blaze::CompressedMatrix;
	using blaze::DynamicMatrix;
	using blaze::SymmetricMatrix;
	using blaze::DynamicVector;
	using blaze::columnVector;
	using blaze::rowMajor;

	int DOF = _parameters.DOF;

	// K stiffness
	double k_11 = _parameters.k_body[0];
	double k_12 = _parameters.k_tyre[0];
	double k_21 = _parameters.k_body[1];
	double k_22 = _parameters.k_tyre[1];
	double k_31 = _parameters.k_body[2];
	double k_32 = _parameters.k_tyre[2];
	double k_41 = _parameters.k_body[3];
	double k_42 = _parameters.k_tyre[3];
	double l_1 = _parameters.l_lat[0];
	double l_2 = _parameters.l_lat[2];
	double l_3 = _parameters.l_long[0];
	double l_4 = _parameters.l_long[2];


	// D - damping matrix
	double d_11 = _parameters.c_body[0];
	double d_12 = _parameters.c_tyre[0];
	double d_21 = _parameters.c_body[1];
	double d_22 = _parameters.c_tyre[1];
	double d_31 = _parameters.c_body[2];
	double d_32 = _parameters.c_tyre[2];
	double d_41 = _parameters.c_body[3];
	double d_42 = _parameters.c_tyre[3];
	double d_l1 = 0.1 * _parameters.l_lat[0];
	double d_l2 = 0.1 * _parameters.l_lat[2];
	double d_l3 = 0.1 * _parameters.l_long[0];
	double d_l4 = 0.1 * _parameters.l_long[2];

	// M - mass matrix
	double m_1 = _parameters.mass_body;
	double m_2 = _parameters.I_body[0];
	double m_3 = _parameters.I_body[2];
	double m_4 = _parameters.mass_wheel[2];
	double m_5 = _parameters.mass_tyre[2];
	double m_6 = _parameters.mass_wheel[3];
	double m_7 = _parameters.mass_tyre[3];
	double m_8 = _parameters.mass_wheel[1];
	double m_9 = _parameters.mass_tyre[1];
	double m_10 = _parameters.mass_wheel[0];
	double m_11 = _parameters.mass_tyre[0];

	// Using symmetric matrix enforce the symmetry and is safer
	SymmetricMatrix <DynamicMatrix<double>, rowMajor> K(DOF);

	K(0, 0) = k_11 + k_21 + k_31 + k_41;
	K(0, 1) = -k_11 * l_1 - k_41 * l_1 + k_21 * l_2 + k_31 * l_2;
	K(0, 2) = -k_11 * l_3 - k_21 * l_3 + k_31 * l_4 + k_41 * l_4;
	K(0, 3) = -k_11;
	K(0, 5) = -k_21;
	K(0, 7) = -k_31;
	K(0, 9) = -k_41;
	K(1, 1) = l_1 * l_1 * k_11 + l_2 * l_2 * k_21 + l_2 * l_2 * k_31 + l_1 * l_1 * k_41;
	K(1, 2) = l_1 * l_3 * k_11 - l_3 * l_2 * k_21 + l_2 * l_4 * k_31 - l_1 * l_4 * k_41;
	K(1, 3) = l_1 * k_11;
	K(1, 5) = -l_2 * k_21;
	K(1, 7) = -l_2 * k_31;
	K(1, 9) = l_1 * k_41;
	K(2, 2) = l_3 * l_3 * k_11 + l_3 * l_3 * k_21 + l_4 * l_4 * k_31 + l_4 * l_4 * k_41;
	K(2, 3) = l_3 * k_11;
	K(2, 5) = l_3 * k_21;
	K(2, 7) = -l_4 * k_31;
	K(2, 9) = -l_4 * k_41;
	K(3, 3) = k_11 + k_12;
	K(3, 4) = -k_12;
	K(4, 4) = k_12;
	K(5, 5) = k_21 + k_22;
	K(5, 6) = -k_22;
	K(6, 6) = k_22;
	K(7, 7) = k_31 + k_32;
	K(7, 8) = -k_32;
	K(8, 8) = k_32;
	K(9, 9) = k_41 + k_42;
	K(9, 10) = -k_42;
	K(10, 10) = k_42;

	SymmetricMatrix <DynamicMatrix<double>, rowMajor> D(DOF);
	D(0, 0) = d_11 + d_21 + d_31 + d_41;
	D(0, 1) = -d_11 * l_1 - d_41 * l_1 + d_21 * l_2 + d_31 * l_2;
	D(0, 2) = -d_11 * l_3 - d_21 * l_3 + d_31 * l_4 + d_41 * l_4;
	D(0, 3) = -d_11;
	D(0, 5) = -d_21;
	D(0, 7) = -d_31;
	D(0, 9) = -d_41;
	D(1, 1) = l_1 * l_1 * d_11 + l_2 * l_2 * d_21 + l_2 * l_2 * d_31 + l_1 * l_1 * d_41;
	D(1, 2) = l_1 * l_3 * d_11 - l_3 * l_2 * d_21 + l_2 * l_4 * d_31 - l_1 * l_4 * d_41;
	D(1, 3) = l_1 * d_11;
	D(1, 5) = -l_2 * d_21;
	D(1, 7) = -l_2 * d_31;
	D(1, 9) = l_1 * d_41;
	D(2, 2) = l_3 * l_3 * d_11 + l_3 * l_3 * d_21 + l_4 * l_4 * d_31 + l_4 * l_4 * d_41;
	D(2, 3) = l_3 * d_11;
	D(2, 5) = l_3 * d_21;
	D(2, 7) = -l_4 * d_31;
	D(2, 9) = -l_4 * d_41;
	D(3, 3) = d_11 + d_12;
	D(3, 4) = -d_12;
	D(4, 4) = d_12;
	D(5, 5) = d_21 + d_22;
	D(5, 6) = -d_22;
	D(6, 6) = d_22;
	D(7, 7) = d_31 + d_32;
	D(7, 8) = -d_32;
	D(8, 8) = d_32;
	D(9, 9) = d_41 + d_42;
	D(9, 10) = -d_42;
	D(10, 10) = d_42;

	SymmetricMatrix <DynamicMatrix<double>, rowMajor> M(DOF);
	M(0, 0) = m_1;
	M(1, 1) = m_2;
	M(2, 2) = m_3;
	M(3, 3) = m_4;
	M(4, 4) = m_5;
	M(5, 5) = m_6;
	M(6, 6) = m_7;
	M(7, 7) = m_8;
	M(8, 8) = m_9;
	M(9, 9) = m_10;
	M(10, 10) = m_11;

	// Define the solution vectors
	DynamicVector<double, columnVector> u_n_p_1(DOF, (double)0.0); // initialize vector of dimension 11 and null elements
	DynamicVector<double, columnVector> u_n(DOF, (double)0.0); // initialize vector of dimension 11 and null elements
	DynamicVector<double, columnVector> u_n_m_1(DOF, (double)0.0); // initialize vector of dimension 11 and null elements

	// Perform the iterations
	//	int nRefinement = 10;
	//	int numTimeSteps = pow(2, nRefinement);
	int numTimeSteps = _parameters.num_time_iter;

	//time step size 
	//double h =  / (numTimeSteps);
	double h = _parameters.timestep;
	std::cout << "Time step h is: " << h << std::scientific << std::endl;

	// Initial conditions
	//const double u_init = 1;
	//const double du_init = 0;
	u_n[0] = 1;
	u_n_m_1[0] = 1;

	/// Build dynamic stiffness matrix
	// A = (1.0/(h*h))*M + (1.0/h)*D + K
	DynamicMatrix<double, rowMajor> A(DOF, DOF);
	A = 1. / (h * h) * M + 1. / h * D + K;

	///Build rhs for BE integrator
	//B = (2.0 / (h*h))*M + 1.0 / h * D + B
	DynamicMatrix<double, rowMajor> B(DOF, DOF);
	B = 2. / (h * h) * M + 1. / h * D;

	// LU Decomposition
	DynamicVector<int, columnVector> ipiv(DOF);   // Pivoting indices
	DynamicMatrix<double, rowMajor>  A_LU(A);  // Temporary matrix to be decomposed

	getrf(A_LU, ipiv.data());
	M *= -1. / (h * h);
	for (int iTime = 0; iTime < numTimeSteps; iTime++) {

		// Solve system: A*u_n_p_1 = B*u_n - M*u_n_m_1
		// rhs = B*u_n + (-1.0 / (h*h))*M
		u_n_p_1 = B * u_n + M * u_n_m_1; // rhs

		getrs(A_LU, u_n_p_1, 'T', ipiv.data());
		u_n_m_1 = u_n;
		u_n = u_n_p_1;
		//reader.ReadLoadParameters(_load_module_parameter);
	}
	std::cout << "We ran #" << numTimeSteps << " time steps!" << std::endl;
	std::cout << u_n_p_1 << std::scientific << std::endl;

#endif
}



void EVAAComputeEngine::computeMBD(void) {
	size_t num_iter = _parameters.num_time_iter;
	const int alignment = 64;
	size_t solution_dim = _parameters.solution_dim;
	floatEVAA* soln = (floatEVAA*)mkl_calloc(solution_dim, sizeof(floatEVAA), alignment);
	MBD_method<floatEVAA> solver(_parameters, _load_module_parameter, lookupStiffness);
	solver.solve(soln);
	std::cout << "mbd: solution after " << num_iter << " timesteps" << std::endl;
	solver.print_final_result(soln);
	mkl_free(soln);
}


void EVAAComputeEngine::clean(void) {

}

void EVAAComputeEngine::computeALE(void) {
	size_t num_iter = _parameters.num_time_iter;
	Car<floatEVAA>* Car1 = new Car<floatEVAA>(_parameters, lookupStiffness);
	size_t solution_dim = Car1->DIM * (size_t) Car1->vec_DIM;
	Profile* Road_Profile;
	if (_load_module_parameter.boundary_condition_road == CIRCULAR) {
		Road_Profile = new Circular(_load_module_parameter.profile_center,
			_load_module_parameter.profile_radius);
	}
	else if (_load_module_parameter.boundary_condition_road == NONFIXED) {
		Road_Profile = new Nonfixed(_load_module_parameter.profile_center,
			_load_module_parameter.profile_radius);
	}
	else {
		std::cout << "ALE will only work with a circular pathor nonfixed boundaries, computation skipped" << std::endl;
		exit(5);
	}

	Road_Profile->update_initial_condition(Car1);

	Load_module* Load_module1 = new Load_module(Road_Profile, Car1, _load_module_parameter);
	linear11dof<floatEVAA>* linear11dof_sys = new linear11dof<floatEVAA>(Car1);
	ALE<floatEVAA>* Ale_sys = new ALE<floatEVAA>(Car1, Load_module1, linear11dof_sys, lookupStiffness, _parameters);

	floatEVAA* soln = (floatEVAA*)mkl_malloc(solution_dim * sizeof(floatEVAA), Car1->alignment);

	Ale_sys->solve(soln);
	std::cout << "ALE: Solution after " << num_iter << " timesteps, f =" << std::endl;
	for (auto i = 0; i < solution_dim; ++i) {
		std::cout << soln[i] << std::endl;
	}
	std::cout << std::endl;
	Ale_sys->print_final_results();


	delete Car1;
	delete Load_module1;

	mkl_free(soln);
}

void EVAAComputeEngine::computeALEtest(void) {

	size_t num_iter = _parameters.num_time_iter;
	size_t solution_dim = _parameters.solution_dim;
	Car<floatEVAA>* Car1 = new Car<floatEVAA>(_parameters, lookupStiffness);
	Profile* Road_Profile = new Circular(_load_module_parameter.profile_center,
		_load_module_parameter.profile_radius);
	Road_Profile->update_initial_condition(Car1);
	Load_module* Load_module1 = new Load_module(Road_Profile, Car1, _load_module_parameter);
	std::cout << "Load module initialized!\n";
	delete Load_module1;
	delete Road_Profile;
	delete Car1;

}


void EVAAComputeEngine::Car_test(void) {

	Car<floatEVAA> Car1(_parameters, lookupStiffness);

}