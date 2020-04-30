/*
 * Copyright &copy; 2019, Dr. Stefan Sicklinger, Munich \n
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

#include "EVAAComputeEngine.h"

#include <fstream>
#include <iostream>
#include <limits>
#include <vector>

#include "11DOF.h"
#include "Car.h"
#include "Constants.h"
#include "LoadModule.h"
#include "MathLibrary.h"
#include "MetaDataBase.h"
#include "Output.h"

#ifdef USE_EIGEN
#include <Eigen/Dense>
using Eigen::IOFormat;
using Eigen::Matrix3d;
using Eigen::MatrixXd;
using Eigen::VectorXd;
#endif

#ifdef USE_BLAZE
#include <blaze/Math.h>
#endif

namespace EVAA {
EVAAComputeEngine::EVAAComputeEngine(std::string xmlCarFileName, std::string xmlLoadFileName) {
    std::ifstream f(xmlCarFileName.c_str());
    if (f.good()) {
        std::cout << "Read general simulation parameters and car input data at " << xmlCarFileName
                  << std::endl;
        // remove cout!
    }
    else {
        std::cout << "XML file at " << xmlCarFileName
                  << " for general and car settings does not exist!" << std::endl;
        exit(2);
    }

    std::ifstream ff(xmlLoadFileName.c_str());
    if (ff.good()) {
        std::cout << "Read load parameters in file: " << xmlLoadFileName << std::endl;
    }
    else {
        std::cout << "XML file at " << xmlLoadFileName << " for load parameters does not exist!"
                  << std::endl;
        exit(2);
    }
    MetaDataBase::getDataBase().readParameters(xmlCarFileName);
    MetaDataBase::getDataBase().readLoadParameters(xmlLoadFileName);
}

void EVAAComputeEngine::printInfo(void) {
    Math::printMKLInfo();

    auto& db = MetaDataBase::getDataBase();
    std::cout << "\n\nCalculate the solution after "
              << db.getNumberOfTimeIterations() * db.getTimeStepSize()
              << "s with dt = " << db.getTimeStepSize() << " for " << db.getNumberOfTimeIterations()
              << " iterations\n\n\n";
}

void EVAAComputeEngine::clean(void) {}

void EVAAComputeEngine::computeEigen11DOF(void) {
    auto& db = MetaDataBase::getDataBase();
    // K stiffness
    double k_11 = db.getBodyStiffnessFrontLeft();
    double k_12 = db.getTyreStiffnessFrontLeft();
    double k_21 = db.getBodyStiffnessFrontRight();
    double k_22 = db.getTyreStiffnessFrontRight();
    double k_31 = db.getBodyStiffnessRearLeft();
    double k_32 = db.getTyreStiffnessRearLeft();
    double k_41 = db.getBodyStiffnessRearRight();
    double k_42 = db.getTyreStiffnessRearRight();
    double l_1 = db.getLatidudalLegPositionFrontLeft();
    double l_2 = db.getLatidudalLegPositionRearLeft();
    double l_3 = db.getLongitudalLegPositionFrontLeft();
    double l_4 = db.getLongitudalLegPositionRearLeft();

    // D - damping matrix
    double d_11 = db.getBodyDampingFrontLeft();
    double d_12 = db.getTyreDampingFrontLeft();
    double d_21 = db.getBodyDampingFrontRight();
    double d_22 = db.getTyreDampingFrontRight();
    double d_31 = db.getBodyDampingRearLeft();
    double d_32 = db.getTyreDampingRearLeft();
    double d_41 = db.getBodyDampingRearRight();
    double d_42 = db.getTyreDampingRearRight();
    double d_l1 = 0.1 * l_1;
    double d_l2 = 0.1 * l_2;
    double d_l3 = 0.1 * l_3;
    double d_l4 = 0.1 * l_4;

    // M - mass matrix
    double m_1 = db.getBodyMass();
    double m_2 = db.getMomentOfInertiaXX();
    double m_3 = db.getMomentOfInertiaYY();
    double m_4 = db.getWheelMassFrontLeft();
    double m_5 = db.getTyreMassFrontLeft();
    double m_6 = db.getWheelMassFrontRight();
    double m_7 = db.getTyreMassFrontRight();
    double m_8 = db.getWheelMassRearLeft();
    double m_9 = db.getTyreMassRearLeft();
    double m_10 = db.getWheelMassRearRight();
    double m_11 = db.getTyreMassRearRight();

    MatrixXd A(Constants::DOF, Constants::DOF);
    MatrixXd B(Constants::DOF, Constants::DOF);
    MatrixXd M(Constants::DOF, Constants::DOF);
    MatrixXd D(Constants::DOF, Constants::DOF);
    MatrixXd K(Constants::DOF, Constants::DOF);
    MatrixXd K_aux(Constants::DOF, Constants::DOF);
    MatrixXd D_aux(Constants::DOF, Constants::DOF);

    VectorXd u_n_p_1(Constants::DOF);
    VectorXd u_n(Constants::DOF);
    VectorXd u_n_m_1(Constants::DOF);

    // Define the mass matrix
    M = MatrixXd::Zero(Constants::DOF, Constants::DOF);
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
    K << k_11 + k_21 + k_31 + k_41, -k_11 * l_1 - k_41 * l_1 + k_21 * l_2 + k_31 * l_2,
        -k_11 * l_3 - k_21 * l_3 + k_31 * l_4 + k_41 * l_4, -k_11, 0., -k_21, 0., -k_31, 0., -k_41,
        0., 0., l_1 * l_1 * k_11 + l_2 * l_2 * k_21 + l_2 * l_2 * k_31 + l_1 * l_1 * k_41,
        l_1 * l_3 * k_11 - l_3 * l_2 * k_21 + l_2 * l_4 * k_31 - l_1 * l_4 * k_41, l_1 * k_11, 0.,
        -l_2 * k_21, 0., -l_2 * k_31, 0., l_1 * k_41, 0., 0., 0.,
        l_3 * l_3 * k_11 + l_3 * l_3 * k_21 + l_4 * l_4 * k_31 + l_4 * l_4 * k_41, l_3 * k_11, 0.,
        l_3 * k_21, 0., -l_4 * k_31, 0., -l_4 * k_41, 0., 0., 0., 0., k_11 + k_12, -k_12, 0., 0.,
        0., 0., 0., 0., 0., 0., 0., 0., k_12, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        k_21 + k_22, -k_22, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., k_22, 0., 0., 0., 0., 0., 0.,
        0., 0., 0., 0., 0., k_31 + k_32, -k_32, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., k_32, 0.,
        0., 0., 0., 0., 0., 0., 0., 0., 0., 0., k_41 + k_42, -k_42, 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., k_42;
    MatrixXd K_diag(Constants::DOF, Constants::DOF);
    K_aux = K + K.transpose();
    K_aux.diagonal() -= K.diagonal();
    K = K_aux;

    // Define the damping matrix
    D << d_11 + d_21 + d_31 + d_41, -d_11 * l_1 - d_41 * l_1 + d_21 * l_2 + d_31 * l_2,
        -d_11 * l_3 - d_21 * l_3 + d_31 * l_4 + d_41 * l_4, -d_11, 0., -d_21, 0., -d_31, 0., -d_41,
        0., 0., l_1 * l_1 * d_11 + l_2 * l_2 * d_21 + l_2 * l_2 * d_31 + l_1 * l_1 * d_41,
        l_1 * l_3 * d_11 - l_3 * l_2 * d_21 + l_2 * l_4 * d_31 - l_1 * l_4 * d_41, l_1 * d_11, 0.,
        -l_2 * d_21, 0., -l_2 * d_31, 0., l_1 * d_41, 0., 0., 0.,
        l_3 * l_3 * d_11 + l_3 * l_3 * d_21 + l_4 * l_4 * d_31 + l_4 * l_4 * d_41, l_3 * d_11, 0.,
        l_3 * d_21, 0., -l_4 * d_31, 0., -l_4 * d_41, 0., 0., 0., 0., d_11 + d_12, -d_12, 0., 0.,
        0., 0., 0., 0., 0., 0., 0., 0., d_12, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        d_21 + d_22, -d_22, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., d_22, 0., 0., 0., 0., 0., 0.,
        0., 0., 0., 0., 0., d_31 + d_32, -d_32, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., d_32, 0.,
        0., 0., 0., 0., 0., 0., 0., 0., 0., 0., d_41 + d_42, -d_42, 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., d_42;
    MatrixXd D_diag(Constants::DOF, Constants::DOF);
    D_aux = D + D.transpose();
    D_aux.diagonal() -= D.diagonal();
    D = D_aux;
    A = MatrixXd::Zero(Constants::DOF, Constants::DOF);
    B = MatrixXd::Zero(Constants::DOF, Constants::DOF);

    u_n_p_1 = VectorXd::Zero(Constants::DOF);
    u_n = VectorXd::Zero(Constants::DOF);
    u_n(0) = 1;
    u_n_m_1 = u_n;

    int numTimeSteps = db.getNumberOfTimeIterations();
    // time step size
    double h = db.getTimeStepSize();
    std::cout << "Time step h is: " << h << std::scientific << std::endl;

    /// Build dynamic stiffness matrix
    // A = (1./(h*h))*M + (1./h)*D + K
    A = (1. / (h * h)) * M + 1. / h * D + K;

    /// Build rhs for BE integrator
    // B = (2.0 / (h*h))*M + 1. / h * D + B
    B = 2.0 / (h * h) * M + 1. / h * D;

    // LU Decomposition
    Eigen::PartialPivLU<MatrixXd> lu(A);

    M *= (-1. / (h * h));
    for (int iTime = 0; iTime < numTimeSteps; iTime++) {
        // B*u_n + (-1. / (h*h))*M
        // Solve system
        u_n_p_1 = lu.solve(B * u_n + M * u_n_m_1);

        u_n_m_1 = u_n;
        u_n = u_n_p_1;
    }
    std::cout << "We ran #" << numTimeSteps << " time steps!" << std::endl;
    std::cout << u_n_p_1 << std::scientific << std::endl;
}

void EVAAComputeEngine::computeBlaze11DOF(void) {
#ifdef USE_BLAZE

    using blaze::columnVector;
    using blaze::CompressedMatrix;
    using blaze::DynamicMatrix;
    using blaze::DynamicVector;
    using blaze::rowMajor;
    using blaze::SymmetricMatrix;

    auto& db = MetaDataBase::getDataBase();

    // K stiffness
    double k_11 = db.getBodyStiffnessFrontLeft();
    double k_12 = db.getTyreStiffnessFrontLeft();
    double k_21 = db.getBodyStiffnessFrontRight();
    double k_22 = db.getTyreStiffnessFrontRight();
    double k_31 = db.getBodyStiffnessRearLeft();
    double k_32 = db.getTyreStiffnessRearLeft();
    double k_41 = db.getBodyStiffnessRearRight();
    double k_42 = db.getTyreStiffnessRearRight();
    double l_1 = db.getLatidudalLegPositionFrontLeft();
    double l_2 = db.getLatidudalLegPositionRearLeft();
    double l_3 = db.getLongitudalLegPositionFrontLeft();
    double l_4 = db.getLongitudalLegPositionRearLeft();

    // D - damping matrix
    double d_11 = db.getBodyDampingFrontLeft();
    double d_12 = db.getTyreDampingFrontLeft();
    double d_21 = db.getBodyDampingFrontRight();
    double d_22 = db.getTyreDampingFrontRight();
    double d_31 = db.getBodyDampingRearLeft();
    double d_32 = db.getTyreDampingRearLeft();
    double d_41 = db.getBodyDampingRearRight();
    double d_42 = db.getTyreDampingRearRight();
    double d_l1 = 0.1 * l_1;
    double d_l2 = 0.1 * l_2;
    double d_l3 = 0.1 * l_3;
    double d_l4 = 0.1 * l_4;

    // M - mass matrix
    double m_1 = db.getBodyMass();
    double m_2 = db.getMomentOfInertiaXX();
    double m_3 = db.getMomentOfInertiaYY();
    double m_4 = db.getWheelMassFrontLeft();
    double m_5 = db.getTyreMassFrontLeft();
    double m_6 = db.getWheelMassFrontRight();
    double m_7 = db.getTyreMassFrontRight();
    double m_8 = db.getWheelMassRearLeft();
    double m_9 = db.getTyreMassRearLeft();
    double m_10 = db.getWheelMassRearRight();
    double m_11 = db.getTyreMassRearRight();

    // Using symmetric matrix enforce the symmetry and is safer
    SymmetricMatrix<DynamicMatrix<double>, rowMajor> K(Constants::DOF);

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

    SymmetricMatrix<DynamicMatrix<double>, rowMajor> D(Constants::DOF);
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

    SymmetricMatrix<DynamicMatrix<double>, rowMajor> M(Constants::DOF);
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
    DynamicVector<double, columnVector> u_n_p_1(
        Constants::DOF, (double)0.);  // initialize vector of dimension 11 and null elements
    DynamicVector<double, columnVector> u_n(
        Constants::DOF, (double)0.);  // initialize vector of dimension 11 and null elements
    DynamicVector<double, columnVector> u_n_m_1(
        Constants::DOF, (double)0.);  // initialize vector of dimension 11 and null elements

    // Perform the iterations
    int numTimeSteps = db.getNumberOfTimeIterations();
    // time step size
    double h = db.getTimeStepSize();

    std::cout << "Time step h is: " << h << std::scientific << std::endl;

    // Initial conditions
    u_n[0] = 1;
    u_n_m_1[0] = 1;

    /// Build dynamic stiffness matrix
    // A = (1./(h*h))*M + (1./h)*D + K
    DynamicMatrix<double, rowMajor> A(Constants::DOF, Constants::DOF);
    A = 1. / (h * h) * M + 1. / h * D + K;

    /// Build rhs for BE integrator
    // B = (2.0 / (h*h))*M + 1. / h * D + B
    DynamicMatrix<double, rowMajor> B(Constants::DOF, Constants::DOF);
    B = 2. / (h * h) * M + 1. / h * D;

    // LU Decomposition
    DynamicVector<int, columnVector> ipiv(Constants::DOF);  // Pivoting indices
    DynamicMatrix<double, rowMajor> A_LU(A);                // Temporary matrix to be decomposed

    getrf(A_LU, ipiv.data());
    M *= -1. / (h * h);
    for (int iTime = 0; iTime < numTimeSteps; iTime++) {
        // Solve system: A*u_n_p_1 = B*u_n - M*u_n_m_1
        // rhs = B*u_n + (-1. / (h*h))*M
        u_n_p_1 = B * u_n + M * u_n_m_1;  // rhs

        getrs(A_LU, u_n_p_1, 'T', ipiv.data());
        u_n_m_1 = u_n;
        u_n = u_n_p_1;
    }
    std::cout << "We ran #" << numTimeSteps << " time steps!" << std::endl;
    std::cout << u_n_p_1 << std::scientific << std::endl;

#endif
}

void EVAAComputeEngine::computeMKL11DOF(void) {
    auto& db = MetaDataBase::getDataBase();

    // TODO: Templated floating type for all doubles below.

    // K stiffness
    double k_11 = db.getBodyStiffnessFrontLeft();
    double k_12 = db.getTyreStiffnessFrontLeft();
    double k_21 = db.getBodyStiffnessFrontRight();
    double k_22 = db.getTyreStiffnessFrontRight();
    double k_31 = db.getBodyStiffnessRearLeft();
    double k_32 = db.getTyreStiffnessRearLeft();
    double k_41 = db.getBodyStiffnessRearRight();
    double k_42 = db.getTyreStiffnessRearRight();
    double l_1 = db.getLatidudalLegPositionFrontLeft();
    double l_2 = db.getLatidudalLegPositionRearLeft();
    double l_3 = db.getLongitudalLegPositionFrontLeft();
    double l_4 = db.getLongitudalLegPositionRearLeft();

    // D - damping matrix
    double d_11 = db.getBodyDampingFrontLeft();
    double d_12 = db.getTyreDampingFrontLeft();
    double d_21 = db.getBodyDampingFrontRight();
    double d_22 = db.getTyreDampingFrontRight();
    double d_31 = db.getBodyDampingRearLeft();
    double d_32 = db.getTyreDampingRearLeft();
    double d_41 = db.getBodyDampingRearRight();
    double d_42 = db.getTyreDampingRearRight();
    double d_l1 = 0.1 * l_1;
    double d_l2 = 0.1 * l_2;
    double d_l3 = 0.1 * l_3;
    double d_l4 = 0.1 * l_4;

    // M - mass matrix
    double m_1 = db.getBodyMass();
    double m_2 = db.getMomentOfInertiaXX();
    double m_3 = db.getMomentOfInertiaYY();
    double m_4 = db.getWheelMassFrontLeft();
    double m_5 = db.getTyreMassFrontLeft();
    double m_6 = db.getWheelMassFrontRight();
    double m_7 = db.getTyreMassFrontRight();
    double m_8 = db.getWheelMassRearLeft();
    double m_9 = db.getTyreMassRearLeft();
    double m_10 = db.getWheelMassRearRight();
    double m_11 = db.getTyreMassRearRight();

    // allocate matrices of zeros
    double* B = Math::calloc<double>(Constants::DOFDOF);
    double* M = Math::calloc<double>(Constants::DOFDOF);
    double* D = Math::calloc<double>(Constants::DOFDOF);
    double* K = Math::calloc<double>(Constants::DOFDOF);

    double* u_n_p_1 = Math::calloc<double>(Constants::DOF);
    double* u_n = Math::calloc<double>(Constants::DOF);
    double* u_n_m_1 = Math::calloc<double>(Constants::DOF);
    double* tmp = Math::calloc<double>(Constants::DOF);

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
    for (int i = 1; i < Constants::DOF; ++i)
        for (int j = 0; j < i; ++j) K[i * Constants::DOF + j] = K[j * Constants::DOF + i];

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
    for (int i = 1; i < Constants::DOF; ++i)
        for (int j = 0; j < i; ++j) D[i * Constants::DOF + j] = D[j * Constants::DOF + i];

    // Initial conditions
    u_n[0] = 1.;
    u_n_m_1[0] = 1.;

    mkl_set_num_threads(8);
    int numTimeSteps = db.getNumberOfTimeIterations();
    // time step size
    double h = db.getTimeStepSize();
    std::cout << "Time step h is: " << h << std::scientific << std::endl;
    /// Build dynamic stiffness matrix
    // gets written in rear vector
    // K' <- (1./(h*h))*M + K
    //	Math::computeDenseVectorAddition(M.data(), K.data(), (1. / (h*h)), 9);
    Math::axpy(121, (1. / (h * h)), M, 1, K, 1);
    // K <- (1./h)*D + K'
    //	Math::computeDenseVectorAddition(D.data(), K.data(), (1. / h), 121);
    Math::axpy(121, (1. / h), D, 1, K, 1);
    /// K holds now dynamic stiffness matrix  for BE integrator
    /// Build rhs for BE integrator
    // B' <-(2.0 / (h*h))*M + B
    //	Math::computeDenseVectorAddition(M.data(), B.data(), (2.0 / (h*h)), 121);
    Math::axpy(Constants::DOFDOF, (2.0 / (h * h)), M, 1, B, 1);

    // B <-(1. / (h))*D + B'
    //	Math::computeDenseVectorAddition(D.data(), B.data(), (1. / h), 121);
    Math::axpy(Constants::DOFDOF, (1. / h), D, 1, B, 1);
    // A*u_n_p_1=B*u_n+C*u_n_m_1+f_n_p_1 <== BE

    std::vector<int> pivot(Constants::DOF);
    // LU Decomposition
    //	Math::computeDenseSymLUFactorisation(11, K, pivot);
    LAPACKE_dgetrf(LAPACK_ROW_MAJOR, Constants::DOF, Constants::DOF, K, Constants::DOF,
                   pivot.data());

    // Time loop
    double tmpScalar = (-1. / (h * h));
    for (int i = 0; i < Constants::DOF; ++i) M[(Constants::DOF + 1) * i] *= tmpScalar;

    void* jitter;
    // adds matrix matrix product to scalar matrix product
    std::cout << mkl_jit_create_dgemm(&jitter, MKL_COL_MAJOR, MKL_NOTRANS, MKL_NOTRANS,
                                      Constants::DOF, 1, Constants::DOF, 1., Constants::DOF,
                                      Constants::DOF, 0., Constants::DOF)
              << std::endl;
    dgemm_jit_kernel_t myDGEMMKernel = mkl_jit_get_dgemm_ptr(jitter);

    for (int iTime = 0; iTime < numTimeSteps; iTime++) {
        // timeVec[iTime] = iTime * h;
        // y: = alpha * A*x + beta * y
        // u_n_p_1 = B*u_n
#ifndef USE_GEMM
        /* void cblas_dgemv(const CBLAS_LAYOUT Layout, const CBLAS_TRANSPOSE trans, const MKL_INT m,
                const MKL_INT n, const double alpha, const double* a, const MKL_INT lda, const
           double* x, const MKL_INT incx, const double beta, double* y, const MKL_INT incy); */
        Math::gemv(CblasColMajor, CblasNoTrans, 11, 11, 1., B, 11, u_n, 1, 0., u_n_p_1, 1);
#endif
#ifdef USE_GEMM
        // cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 11, 1, 11, 1., B, 11, u_n, 11,
        // 0., u_n_p_1, 11);
        myDGEMMKernel(jitter, B, u_n, u_n_p_1);
#endif
        // y: = alpha*A*x + beta*y
        // tmp = ((-1. / (h*h))*M) * u_n_m_1
#ifndef USE_GEMM
        Math::gemv(CblasColMajor, CblasNoTrans, Constants::DOF, Constants::DOF, (-1. / (h * h)), M,
                   Constants::DOF, u_n_m_1, 1, 0., tmp, 1);
#endif
#ifdef USE_GEMM
        // cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 11, 1, 11, 1., M, 11, u_n_m_1,
        // 11, 0., tmp, 11);
        myDGEMMKernel(jitter, M, u_n_m_1, tmp);
#endif
        // u_n_p_1 <- 1. tmp + u_n_p_1
        //		Math::computeDenseVectorAddition(tmp.data(), u_n_p_1.data(), 1.,
        // 11); void cblas_daxpy (const MKL_INT n, const double a, const double *x, const MKL_INT
        // incx, double *y, const MKL_INT incy);
        Math::axpy(Constants::DOF, 1., tmp, 1, u_n_p_1, 1);
        // Solve system
        //		Math::computeDenseSymSolution(11, K, pivot, u_n_p_1);
        // lapack_int LAPACKE_dgetrs (int matrix_layout , char trans , lapack_int n , lapack_int
        // nrhs , const double * a , lapack_int lda , const lapack_int * ipiv , double * b ,
        // lapack_int ldb );
        LAPACKE_dgetrs(LAPACK_ROW_MAJOR, 'N', Constants::DOF, 1, K, Constants::DOF, pivot.data(),
                       u_n_p_1, 1);

        for (int i = 0; i < Constants::DOF; ++i) {
            u_n_m_1[i] = u_n[i];
            u_n[i] = u_n_p_1[i];
        }
    }
    std::cout << "We ran #" << numTimeSteps << " time steps!" << std::endl;
    for (int i = 0; i < Constants::DOF; ++i) {
        std::cout << u_n_p_1[i] << std::scientific << std::endl;
    }

    // free the memory
    Math::free(B);
    Math::free(M);
    Math::free(D);
    Math::free(K);

    Math::free(u_n_p_1);
    Math::free(u_n);
    Math::free(u_n_m_1);
    Math::free(tmp);
}

void EVAAComputeEngine::computeMKLTwoTrackModelBE(void) {
    if (MetaDataBase::getDataBase().getRoadConditions() == NONFIXED) {
        Constants::floatEVAA* sol = Math::malloc<Constants::floatEVAA>(Constants::DOF);
        Car<Constants::floatEVAA>* car = new Car<Constants::floatEVAA>();
        TwoTrackModelFull<Constants::floatEVAA> solver(car);
        solver.apply_boundary_condition(MetaDataBase::getDataBase().getRoadConditions());
        solver.solve(sol);

        solver.print_final_results(sol);

        Math::free(sol);
        delete car;
    }
    else {
        std::cout << "Linear11dof solver will only work with NONFIXED boundary conditions, "
                     "computation skipped"
                  << std::endl;
    }
}

#if MIGHT_BE_USEFUL
void EVAAComputeEngine::computeMKLTwoTrackModel() {
    if (MetaDataBase::getDataBase().getRoadConditions() == NONFIXED) {
        Constants::floatEVAA* sol = Math::calloc<Constants::floatEVAA>(Constants::DOF);
        Car<Constants::floatEVAA>* car = new Car<Constants::floatEVAA>();
        TwoTrackModelFull<Constants::floatEVAA, TwoTrackModelBDF2<Constants::floatEVAA>> solver(
            car);
        solver.apply_boundary_condition(MetaDataBase::getDataBase().getRoadConditions());
        solver.solve(sol);

        solver.print_final_results(sol);

        Math::free(sol);
        delete car;
    }
    else {
        std::cout << "TwoTrackModel solver will only work with NONFIXED boundary conditions, "
                     "computation skipped "
                  << std::endl;
    }
}
#endif

void EVAAComputeEngine::computeMBD(void) {
    size_t num_iter = MetaDataBase::getDataBase().getNumberOfTimeIterations();
    size_t solution_dim = MetaDataBase::getDataBase().getSolutionVectorSize();
    Constants::floatEVAA* sol = Math::calloc<Constants::floatEVAA>(solution_dim);
    MBD_method<Constants::floatEVAA> solver;
    solver.solve(sol);
    solver.print_final_result(sol);
    Math::free(sol);
}

void EVAAComputeEngine::computeALE(void) {
    Profile<Constants::floatEVAA>* roadProfile;

    Car<Constants::floatEVAA>* car = new Car<Constants::floatEVAA>();

    auto& db = MetaDataBase::getDataBase();
    switch (db.getRoadConditions()) {
    case CIRCULAR:
        roadProfile = new Circular<Constants::floatEVAA>(db.getCircularRoadCenter(),
                                                         db.getCircularRoadRadius());
        break;
    case NONFIXED:
        roadProfile = new Nonfixed<Constants::floatEVAA>(db.getCircularRoadCenter(),
                                                         db.getCircularRoadRadius());
        break;
    case FIXED:
        roadProfile = new Fixed<Constants::floatEVAA>(db.getGravityField()[1]);
        roadProfile->set_fixed_index(car->tyre_index_set);
        break;
    default:
        std::cout << "ALE will only work with a circular path, fixed or nonfixed boundaries, "
                     "computation skipped!"
                  << std::endl;
        delete car;
        exit(5);
        break;
    }

    roadProfile->update_initial_condition(car);
    LoadModule<Constants::floatEVAA>* loadModule =
        new LoadModule<Constants::floatEVAA>(roadProfile, car);
    TwoTrackModelParent<Constants::floatEVAA>* TwoTrackModel_obj =
        new TwoTrackModelBE<Constants::floatEVAA>(car);
    ALE<Constants::floatEVAA>* ale =
        new ALE<Constants::floatEVAA>(car, loadModule, TwoTrackModel_obj);

    size_t solutionDim = Constants::DIM * (size_t)Constants::VEC_DIM;
    Constants::floatEVAA* sol = Math::malloc<Constants::floatEVAA>(solutionDim);

    ale->solve(sol);

    ale->print_final_results();

    delete car;
    delete loadModule;
    delete roadProfile;
    delete TwoTrackModel_obj;
    delete ale;

    Math::free(sol);
}

void EVAAComputeEngine::computeALEtest(void) {
    auto& db = MetaDataBase::getDataBase();
    size_t num_iter = db.getNumberOfTimeIterations();
    size_t solution_dim = db.getSolutionVectorSize();
    Car<Constants::floatEVAA>* car = new Car<Constants::floatEVAA>();
    Profile<Constants::floatEVAA>* roadProfile =
        new Circular<Constants::floatEVAA>(db.getCircularRoadCenter(), db.getCircularRoadRadius());
    roadProfile->update_initial_condition(car);
    LoadModule<Constants::floatEVAA>* loadModule =
        new LoadModule<Constants::floatEVAA>(roadProfile, car);
    std::cout << "Load module initialized!\n";
    delete loadModule;
    delete roadProfile;
    delete car;
}

#if MIGHT_BE_USEFUL
void EVAAComputeEngine::compare_ALE_MBD(void) {
    // MBD Call
    size_t num_iter = _parameters.num_time_iter;
    size_t solution_dim = _parameters.solution_dim;
    Constants::floatEVAA* soln = Math::calloc<Constants::floatEVAA>(solution_dim);
    MBD_method<Constants::floatEVAA> solver(_parameters);
    size_t solution_size = (num_iter + 1) * solution_dim;
    Constants::floatEVAA* complete_soln = Math::malloc<Constants::floatEVAA>(solution_size);
    solver.solve(soln, complete_soln);
    solver.print_final_result(soln);
    std::cout << "(num_iter + 1) = " << (num_iter + 1) << "solution_dim = " << solution_dim
              << std::endl;
#ifdef IO
    IO::write_matrix(complete_soln, "MBD_result.dat", (num_iter + 1), solution_dim);
#endif
    // IO
    Math::scal(solution_dim, 0., soln, 1);
    Math::free(complete_soln);
    // ALE call

    Profile* Road_Profile;

    Car<Constants::floatEVAA>* Car1 = new Car<Constants::floatEVAA>(_parameters, _lookupStiffness);

    if (_loadModuleParameter.boundary_condition_road == CIRCULAR) {
        Road_Profile =
            new Circular(_loadModuleParameter.profile_center, _loadModuleParameter.profile_radius);
    }
    else if (_loadModuleParameter.boundary_condition_road == NONFIXED) {
        Road_Profile =
            new Nonfixed(_loadModuleParameter.profile_center, _loadModuleParameter.profile_radius);
    }
    else if (_loadModuleParameter.boundary_condition_road == FIXED) {
        Road_Profile = new Fixed(_parameters.gravity[2], _loadModuleParameter);
        Road_Profile->set_fixed_index(Car1->tyre_index_set);
    }
    else {
        std::cout << "ALE will only work with a circular path, fixed or nonfixed boundaries, "
                     "computation skipped "
                  << std::endl;
        delete Car1;
        exit(5);
    }

    solution_dim = Constants::DIM * Constants::VEC_DIM;
    solution_size = (num_iter + 1) * solution_dim;
    Constants::floatEVAA* complete_soln2 = Math::calloc<Constants::floatEVAA>(solution_size);
#ifdef IO
    IO::write_matrix(Car1->Position_vec, "initial_car_pos_vec.dat", 1, solution_dim);
#endif  // IO
    Road_Profile->update_initial_condition(Car1);

    Load_module* Load_module1 = new Load_module(Road_Profile, Car1, _loadModuleParameter);
    TwoTrackModel<Constants::floatEVAA>* TwoTrackModel_sys =
        new TwoTrackModel<Constants::floatEVAA>(Car1);
    ALE<Constants::floatEVAA>* Ale_sys = new ALE<Constants::floatEVAA>(
        Car1, Load_module1, TwoTrackModel_sys, _lookupStiffness, _parameters);

    Ale_sys->solve(soln, complete_soln2);

    Ale_sys->print_final_results();
#ifdef IO
    IO::write_matrix(complete_soln2, "ALE_result.dat", (num_iter + 1), solution_dim);
#endif  // IO
    delete Car1;
    delete Load_module1;
    delete Road_Profile;
    delete TwoTrackModel_sys;
    delete Ale_sys;

    Math::free(soln);
    Math::free(complete_soln2);
}
#endif

}  // namespace EVAA
