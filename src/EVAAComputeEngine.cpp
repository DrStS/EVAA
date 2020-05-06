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
#include "MetaDatabase.h"
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
    IO::checkFileExists(xmlCarFileName);
    IO::checkFileExists(xmlLoadFileName);

    MetaDatabase<Constants::floatEVAA>::getDatabase().readParameters(xmlCarFileName);
    MetaDatabase<Constants::floatEVAA>::getDatabase().readLoadParameters(xmlLoadFileName);
}

void EVAAComputeEngine::printInfo(void) {
    Math::PrintMKLInfo();
    auto& db = MetaDatabase<Constants::floatEVAA>::getDatabase();
    std::cout << "\n\nCalculate the solution after " << db.getNumberOfTimeIterations() * db.getTimeStepSize() << "s with dt = " << db.getTimeStepSize() << " for " << db.getNumberOfTimeIterations() << " iterations\n\n\n";
}

void EVAAComputeEngine::computeEigen11DOF(void) {
    auto& db = MetaDatabase<Constants::floatEVAA>::getDatabase();

    // K stiffness
    Constants::floatEVAA k_11 = db.getBodyStiffnessFrontLeft();
    Constants::floatEVAA k_12 = db.getTyreStiffnessFrontLeft();
    Constants::floatEVAA k_21 = db.getBodyStiffnessFrontRight();
    Constants::floatEVAA k_22 = db.getTyreStiffnessFrontRight();
    Constants::floatEVAA k_31 = db.getBodyStiffnessRearLeft();
    Constants::floatEVAA k_32 = db.getTyreStiffnessRearLeft();
    Constants::floatEVAA k_41 = db.getBodyStiffnessRearRight();
    Constants::floatEVAA k_42 = db.getTyreStiffnessRearRight();
    Constants::floatEVAA l_1 = db.getLatidudalLegPositionFrontLeft();
    Constants::floatEVAA l_2 = db.getLatidudalLegPositionRearLeft();
    Constants::floatEVAA l_3 = db.getLongitudalLegPositionFrontLeft();
    Constants::floatEVAA l_4 = db.getLongitudalLegPositionRearLeft();

    // D - damping matrix
    Constants::floatEVAA d_11 = db.getBodyDampingFrontLeft();
    Constants::floatEVAA d_12 = db.getTyreDampingFrontLeft();
    Constants::floatEVAA d_21 = db.getBodyDampingFrontRight();
    Constants::floatEVAA d_22 = db.getTyreDampingFrontRight();
    Constants::floatEVAA d_31 = db.getBodyDampingRearLeft();
    Constants::floatEVAA d_32 = db.getTyreDampingRearLeft();
    Constants::floatEVAA d_41 = db.getBodyDampingRearRight();
    Constants::floatEVAA d_42 = db.getTyreDampingRearRight();
    Constants::floatEVAA d_l1 = 0.1 * l_1;
    Constants::floatEVAA d_l2 = 0.1 * l_2;
    Constants::floatEVAA d_l3 = 0.1 * l_3;
    Constants::floatEVAA d_l4 = 0.1 * l_4;

    // M - mass matrix
    Constants::floatEVAA m_1 = db.getBodyMass();
    Constants::floatEVAA m_2 = db.getMomentOfInertiaXX();
    Constants::floatEVAA m_3 = db.getMomentOfInertiaYY();
    Constants::floatEVAA m_4 = db.getWheelMassFrontLeft();
    Constants::floatEVAA m_5 = db.getTyreMassFrontLeft();
    Constants::floatEVAA m_6 = db.getWheelMassFrontRight();
    Constants::floatEVAA m_7 = db.getTyreMassFrontRight();
    Constants::floatEVAA m_8 = db.getWheelMassRearLeft();
    Constants::floatEVAA m_9 = db.getTyreMassRearLeft();
    Constants::floatEVAA m_10 = db.getWheelMassRearRight();
    Constants::floatEVAA m_11 = db.getTyreMassRearRight();

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
    K << k_11 + k_21 + k_31 + k_41, -k_11 * l_1 - k_41 * l_1 + k_21 * l_2 + k_31 * l_2, -k_11 * l_3 - k_21 * l_3 + k_31 * l_4 + k_41 * l_4, -k_11, 0, -k_21, 0, -k_31, 0, -k_41, 0, 0, l_1 * l_1 * k_11 + l_2 * l_2 * k_21 + l_2 * l_2 * k_31 + l_1 * l_1 * k_41, l_1 * l_3 * k_11 - l_3 * l_2 * k_21 + l_2 * l_4 * k_31 - l_1 * l_4 * k_41, l_1 * k_11, 0, -l_2 * k_21, 0, -l_2 * k_31, 0, l_1 * k_41, 0, 0, 0, l_3 * l_3 * k_11 + l_3 * l_3 * k_21 + l_4 * l_4 * k_31 + l_4 * l_4 * k_41, l_3 * k_11, 0, l_3 * k_21, 0, -l_4 * k_31, 0, -l_4 * k_41, 0, 0, 0, 0, k_11 + k_12, -k_12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, k_12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, k_21 + k_22, -k_22, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, k_22, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, k_31 + k_32, -k_32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, k_32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, k_41 + k_42, -k_42, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, k_42;
    MatrixXd K_diag(Constants::DOF, Constants::DOF);
    K_aux = K + K.transpose();
    K_aux.diagonal() -= K.diagonal();
    K = K_aux;

    // Define the damping matrix
    D << d_11 + d_21 + d_31 + d_41, -d_11 * l_1 - d_41 * l_1 + d_21 * l_2 + d_31 * l_2, -d_11 * l_3 - d_21 * l_3 + d_31 * l_4 + d_41 * l_4, -d_11, 0, -d_21, 0, -d_31, 0, -d_41, 0, 0, l_1 * l_1 * d_11 + l_2 * l_2 * d_21 + l_2 * l_2 * d_31 + l_1 * l_1 * d_41, l_1 * l_3 * d_11 - l_3 * l_2 * d_21 + l_2 * l_4 * d_31 - l_1 * l_4 * d_41, l_1 * d_11, 0, -l_2 * d_21, 0, -l_2 * d_31, 0, l_1 * d_41, 0, 0, 0, l_3 * l_3 * d_11 + l_3 * l_3 * d_21 + l_4 * l_4 * d_31 + l_4 * l_4 * d_41, l_3 * d_11, 0, l_3 * d_21, 0, -l_4 * d_31, 0, -l_4 * d_41, 0, 0, 0, 0, d_11 + d_12, -d_12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, d_12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, d_21 + d_22, -d_22, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, d_22, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, d_31 + d_32, -d_32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, d_32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, d_41 + d_42, -d_42, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, d_42;
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
    Constants::floatEVAA h = db.getTimeStepSize();
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

    auto& db = MetaDatabase<Constants::floatEVAA>::getDatabase();

    // K stiffness
    Constants::floatEVAA k_11 = db.getBodyStiffnessFrontLeft();
    Constants::floatEVAA k_12 = db.getTyreStiffnessFrontLeft();
    Constants::floatEVAA k_21 = db.getBodyStiffnessFrontRight();
    Constants::floatEVAA k_22 = db.getTyreStiffnessFrontRight();
    Constants::floatEVAA k_31 = db.getBodyStiffnessRearLeft();
    Constants::floatEVAA k_32 = db.getTyreStiffnessRearLeft();
    Constants::floatEVAA k_41 = db.getBodyStiffnessRearRight();
    Constants::floatEVAA k_42 = db.getTyreStiffnessRearRight();
    Constants::floatEVAA l_1 = db.getLatidudalLegPositionFrontLeft();
    Constants::floatEVAA l_2 = db.getLatidudalLegPositionRearLeft();
    Constants::floatEVAA l_3 = db.getLongitudalLegPositionFrontLeft();
    Constants::floatEVAA l_4 = db.getLongitudalLegPositionRearLeft();

    // D - damping matrix
    Constants::floatEVAA d_11 = db.getBodyDampingFrontLeft();
    Constants::floatEVAA d_12 = db.getTyreDampingFrontLeft();
    Constants::floatEVAA d_21 = db.getBodyDampingFrontRight();
    Constants::floatEVAA d_22 = db.getTyreDampingFrontRight();
    Constants::floatEVAA d_31 = db.getBodyDampingRearLeft();
    Constants::floatEVAA d_32 = db.getTyreDampingRearLeft();
    Constants::floatEVAA d_41 = db.getBodyDampingRearRight();
    Constants::floatEVAA d_42 = db.getTyreDampingRearRight();
    Constants::floatEVAA d_l1 = 0.1 * l_1;
    Constants::floatEVAA d_l2 = 0.1 * l_2;
    Constants::floatEVAA d_l3 = 0.1 * l_3;
    Constants::floatEVAA d_l4 = 0.1 * l_4;

    // M - mass matrix
    Constants::floatEVAA m_1 = db.getBodyMass();
    Constants::floatEVAA m_2 = db.getMomentOfInertiaXX();
    Constants::floatEVAA m_3 = db.getMomentOfInertiaYY();
    Constants::floatEVAA m_4 = db.getWheelMassFrontLeft();
    Constants::floatEVAA m_5 = db.getTyreMassFrontLeft();
    Constants::floatEVAA m_6 = db.getWheelMassFrontRight();
    Constants::floatEVAA m_7 = db.getTyreMassFrontRight();
    Constants::floatEVAA m_8 = db.getWheelMassRearLeft();
    Constants::floatEVAA m_9 = db.getTyreMassRearLeft();
    Constants::floatEVAA m_10 = db.getWheelMassRearRight();
    Constants::floatEVAA m_11 = db.getTyreMassRearRight();

    // Using symmetric matrix enforce the symmetry and is safer
    SymmetricMatrix<DynamicMatrix<Constants::floatEVAA>, rowMajor> K(Constants::DOF);

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

    SymmetricMatrix<DynamicMatrix<Constants::floatEVAA>, rowMajor> D(Constants::DOF);
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

    SymmetricMatrix<DynamicMatrix<Constants::floatEVAA>, rowMajor> M(Constants::DOF);
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
    DynamicVector<Constants::floatEVAA, columnVector> u_n_p_1(Constants::DOF,
                                                              0);  // initialize vector of dimension 11 and null elements
    DynamicVector<Constants::floatEVAA, columnVector> u_n(Constants::DOF,
                                                          0);  // initialize vector of dimension 11 and null elements
    DynamicVector<Constants::floatEVAA, columnVector> u_n_m_1(Constants::DOF,
                                                              0);  // initialize vector of dimension 11 and null elements

    // Perform the iterations
    int numTimeSteps = db.getNumberOfTimeIterations();
    // time step size
    Constants::floatEVAA h = db.getTimeStepSize();

    std::cout << "Time step h is: " << h << std::scientific << std::endl;

    // Initial conditions
    u_n[0] = 1;
    u_n_m_1[0] = 1;

    /// Build dynamic stiffness matrix
    // A = (1./(h*h))*M + (1./h)*D + K
    DynamicMatrix<Constants::floatEVAA, rowMajor> A(Constants::DOF, Constants::DOF);
    A = 1. / (h * h) * M + 1. / h * D + K;

    /// Build rhs for BE integrator
    // B = (2.0 / (h*h))*M + 1. / h * D + B
    DynamicMatrix<Constants::floatEVAA, rowMajor> B(Constants::DOF, Constants::DOF);
    B = 2. / (h * h) * M + 1. / h * D;

    // LU Decomposition
    DynamicVector<int, columnVector> ipiv(Constants::DOF);  // Pivoting indices
    DynamicMatrix<Constants::floatEVAA, rowMajor> A_LU(A);  // Temporary matrix to be decomposed

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
#if TODO_TEMPLATE_JIT_DONE
    auto& db = MetaDatabase<Constants::floatEVAA>::getDatabase();

    // K stiffness
    Constants::floatEVAA k_11 = db.getBodyStiffnessFrontLeft();
    Constants::floatEVAA k_12 = db.getTyreStiffnessFrontLeft();
    Constants::floatEVAA k_21 = db.getBodyStiffnessFrontRight();
    Constants::floatEVAA k_22 = db.getTyreStiffnessFrontRight();
    Constants::floatEVAA k_31 = db.getBodyStiffnessRearLeft();
    Constants::floatEVAA k_32 = db.getTyreStiffnessRearLeft();
    Constants::floatEVAA k_41 = db.getBodyStiffnessRearRight();
    Constants::floatEVAA k_42 = db.getTyreStiffnessRearRight();
    Constants::floatEVAA l_1 = db.getLatidudalLegPositionFrontLeft();
    Constants::floatEVAA l_2 = db.getLatidudalLegPositionRearLeft();
    Constants::floatEVAA l_3 = db.getLongitudalLegPositionFrontLeft();
    Constants::floatEVAA l_4 = db.getLongitudalLegPositionRearLeft();

    // D - damping matrix
    Constants::floatEVAA d_11 = db.getBodyDampingFrontLeft();
    Constants::floatEVAA d_12 = db.getTyreDampingFrontLeft();
    Constants::floatEVAA d_21 = db.getBodyDampingFrontRight();
    Constants::floatEVAA d_22 = db.getTyreDampingFrontRight();
    Constants::floatEVAA d_31 = db.getBodyDampingRearLeft();
    Constants::floatEVAA d_32 = db.getTyreDampingRearLeft();
    Constants::floatEVAA d_41 = db.getBodyDampingRearRight();
    Constants::floatEVAA d_42 = db.getTyreDampingRearRight();
    Constants::floatEVAA d_l1 = 0.1 * l_1;
    Constants::floatEVAA d_l2 = 0.1 * l_2;
    Constants::floatEVAA d_l3 = 0.1 * l_3;
    Constants::floatEVAA d_l4 = 0.1 * l_4;

    // M - mass matrix
    Constants::floatEVAA m_1 = db.getBodyMass();
    Constants::floatEVAA m_2 = db.getMomentOfInertiaXX();
    Constants::floatEVAA m_3 = db.getMomentOfInertiaYY();
    Constants::floatEVAA m_4 = db.getWheelMassFrontLeft();
    Constants::floatEVAA m_5 = db.getTyreMassFrontLeft();
    Constants::floatEVAA m_6 = db.getWheelMassFrontRight();
    Constants::floatEVAA m_7 = db.getTyreMassFrontRight();
    Constants::floatEVAA m_8 = db.getWheelMassRearLeft();
    Constants::floatEVAA m_9 = db.getTyreMassRearLeft();
    Constants::floatEVAA m_10 = db.getWheelMassRearRight();
    Constants::floatEVAA m_11 = db.getTyreMassRearRight();

    // allocate matrices of zeros
    Constants::floatEVAA* B = Math::calloc<Constants::floatEVAA>(Constants::DOFDOF);
    Constants::floatEVAA* M = Math::calloc<Constants::floatEVAA>(Constants::DOFDOF);
    Constants::floatEVAA* D = Math::calloc<Constants::floatEVAA>(Constants::DOFDOF);
    Constants::floatEVAA* K = Math::calloc<Constants::floatEVAA>(Constants::DOFDOF);

    Constants::floatEVAA* u_n_p_1 = Math::calloc<Constants::floatEVAA>(Constants::DOF);
    Constants::floatEVAA* u_n = Math::calloc<Constants::floatEVAA>(Constants::DOF);
    Constants::floatEVAA* u_n_m_1 = Math::calloc<Constants::floatEVAA>(Constants::DOF);
    Constants::floatEVAA* tmp = Math::calloc<Constants::floatEVAA>(Constants::DOF);

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
    u_n[0] = 1;
    u_n_m_1[0] = 1;

    mkl_set_num_threads(8);
    int numTimeSteps = db.getNumberOfTimeIterations();
    // time step size
    Constants::floatEVAA h = db.getTimeStepSize();
    std::cout << "Time step h is: " << h << std::scientific << std::endl;
    /// Build dynamic stiffness matrix
    // gets written in rear vector
    // K' <- (1./(h*h))*M + K
    //	Math::ComputeDenseVectorAddition<Constants::floatEVAA>(M.data(),
    // K.data(), (1. / (h*h)), 9);
    Math::axpy<Constants::floatEVAA>(121, (1. / (h * h)), M, 1, K, 1);
    // K <- (1./h)*D + K'
    //	Math::ComputeDenseVectorAddition<Constants::floatEVAA>(D.data(),
    // K.data(), (1. / h), 121);
    Math::axpy<Constants::floatEVAA>(121, (1. / h), D, 1, K, 1);
    /// K holds now dynamic stiffness matrix  for BE integrator
    /// Build rhs for BE integrator
    // B' <-(2.0 / (h*h))*M + B
    //	Math::ComputeDenseVectorAddition<Constants::floatEVAA>(M.data(),
    // B.data(), (2.0 / (h*h)),
    // 121);
    Math::axpy<Constants::floatEVAA>(Constants::DOFDOF, (2.0 / (h * h)), M, 1, B, 1);

    // B <-(1. / (h))*D + B'
    //	Math::ComputeDenseVectorAddition<Constants::floatEVAA>(D.data(),
    // B.data(), (1. / h), 121);
    Math::axpy<Constants::floatEVAA>(Constants::DOFDOF, (1. / h), D, 1, B, 1);
    // A*u_n_p_1=B*u_n+C*u_n_m_1+f_n_p_1 <== BE

    std::vector<int> pivot(Constants::DOF);
    // LU Decomposition
    //	Math::ComputeDenseSymLUFactorisation<Constants::floatEVAA>(11, K,
    // pivot);
    Math::getrf(LAPACK_ROW_MAJOR, Constants::DOF, Constants::DOF, K, Constants::DOF, pivot.data());

    // Time loop
    Constants::floatEVAA tmpScalar = (-1. / (h * h));
    for (int i = 0; i < Constants::DOF; ++i) M[(Constants::DOF + 1) * i] *= tmpScalar;

    void* jitter;
    // adds matrix matrix product to scalar matrix product
    std::cout << mkl_jit_create_dgemm(&jitter, MKL_COL_MAJOR, MKL_NOTRANS, MKL_NOTRANS, Constants::DOF, 1, Constants::DOF, 1, Constants::DOF, Constants::DOF, 0, Constants::DOF) << std::endl;
    dgemm_jit_kernel_t myDGEMMKernel = mkl_jit_get_dgemm_ptr(jitter);

    for (int iTime = 0; iTime < numTimeSteps; iTime++) {
        // timeVec[iTime] = iTime * h;
        // y: = alpha * A*x + beta * y
        // u_n_p_1 = B*u_n
#ifndef USE_GEMM
        /* void cblas_dgemv(const CBLAS_LAYOUT Layout, const CBLAS_TRANSPOSE
           trans, const MKL_INT m, const MKL_INT n, const double alpha, const
           double* a, const MKL_INT lda, const double* x, const MKL_INT incx,
           const double beta, double* y, const MKL_INT incy); */
        Math::gemv<Constants::floatEVAA>(CblasColMajor, CblasNoTrans, 11, 11, 1, B, 11, u_n, 1, 0, u_n_p_1, 1);
#endif
#ifdef USE_GEMM
        // cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 11, 1, 11,1,
        // B, 11, u_n, 11, 0, u_n_p_1, 11);
        myDGEMMKernel(jitter, B, u_n, u_n_p_1);
#endif
        // y: = alpha*A*x + beta*y
        // tmp = ((-1. / (h*h))*M) * u_n_m_1
#ifndef USE_GEMM
        Math::gemv<Constants::floatEVAA>(CblasColMajor, CblasNoTrans, Constants::DOF, Constants::DOF,  //
                                         (-1. / (h * h)), M, Constants::DOF, u_n_m_1, 1, 0, tmp, 1);   //
#endif
#ifdef USE_GEMM
        // cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 11, 1, 11,1,
        // M, 11, u_n_m_1, 11,0, tmp, 11);
        myDGEMMKernel(jitter, M, u_n_m_1, tmp);
#endif
        // u_n_p_1 <- 1. tmp + u_n_p_1
        //		Math::ComputeDenseVectorAddition<Constants::floatEVAA>(tmp.data(),
        // u_n_p_1.data(),1,
        // 11); void cblas_daxpy (const MKL_INT n, const double a, const double
        // *x, const MKL_INT incx, double *y, const MKL_INT incy);
        Math::axpy<Constants::floatEVAA>(Constants::DOF, 1, tmp, 1, u_n_p_1, 1);
        // Solve system
        //		Math::ComputeDenseSymSolution<Constants::floatEVAA>(11,
        // K, pivot, u_n_p_1);
        // lapack_int LAPACKE_dgetrs (int matrix_layout , char trans ,
        // lapack_int n , lapack_int nrhs , const double * a , lapack_int lda ,
        // const lapack_int * ipiv , double * b , lapack_int ldb );
        LAPACKE_dgetrs(LAPACK_ROW_MAJOR, 'N', Constants::DOF, 1, K, Constants::DOF, pivot.data(), u_n_p_1, 1);

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
    Math::free<Constants::floatEVAA>(B);
    Math::free<Constants::floatEVAA>(M);
    Math::free<Constants::floatEVAA>(D);
    Math::free<Constants::floatEVAA>(K);

    Math::free<Constants::floatEVAA>(u_n_p_1);
    Math::free<Constants::floatEVAA>(u_n);
    Math::free<Constants::floatEVAA>(u_n_m_1);
    Math::free<Constants::floatEVAA>(tmp);
#endif  // TODO_TEMPLATE_JIT_DONE
}

void EVAAComputeEngine::computeMKLTwoTrackModelBE(void) {
    auto& db = MetaDatabase<Constants::floatEVAA>::getDatabase();
    if (true) { // TODO remove this
        Constants::floatEVAA* sol = Math::malloc<Constants::floatEVAA>(Constants::DOF);
        Car<Constants::floatEVAA>* car = new Car<Constants::floatEVAA>();
		Lagrange<Constants::floatEVAA>* lagrange = new Straight<Constants::floatEVAA>();
		Euler<Constants::floatEVAA>* euler = new Nonfixed<Constants::floatEVAA>(db.getGravityField()[2]);
		LoadModule<Constants::floatEVAA>* load = new LoadModule<Constants::floatEVAA>(lagrange, euler, car);
        TwoTrackModelFull<Constants::floatEVAA> solver(car, load);
		euler->ApplyProfileInitialCondition(car);
        solver.Solve(sol);

        solver.PrintFinalResults(sol);

        Math::free<Constants::floatEVAA>(sol);
        delete car;
        delete lagrange;
        delete euler;
        delete load;
    }
    else {
        std::cout << "Linear11dof solver will only work with NONFIXED boundary "
                     "conditions, "
                     "computation skipped"
                  << std::endl;
    }
}

void EVAAComputeEngine::computeMBD(void) {
    size_t num_iter = MetaDatabase<Constants::floatEVAA>::getDatabase().getNumberOfTimeIterations();
    size_t solution_dim = MetaDatabase<Constants::floatEVAA>::getDatabase().getSolutionVectorSize();
    Constants::floatEVAA* sol = Math::calloc<Constants::floatEVAA>(solution_dim);
    MBDMethod<Constants::floatEVAA> solver;

    solver.Solve(sol);
    solver.PrintFinalResult(sol);
    Math::free<Constants::floatEVAA>(sol);
}

void EVAAComputeEngine::computeALE(void) {
    Lagrange<Constants::floatEVAA>* lagrangeProfile;
	Euler<Constants::floatEVAA>* eulerProfile;
	Car<Constants::floatEVAA>* car = new Car<Constants::floatEVAA>();

    auto& db = MetaDatabase<Constants::floatEVAA>::getDatabase();
	lagrangeProfile = new Circular<Constants::floatEVAA>(db.getCircularRoadCenter());
    //lagrangeProfile = new Straight<Constants::floatEVAA>();
	eulerProfile = new Fixed<Constants::floatEVAA>(db.getGravityField()[2]);
	lagrangeProfile->ApplyProfileInitialCondition(car);
	eulerProfile->ApplyProfileInitialCondition(car);
	LoadModule<Constants::floatEVAA>* loadModule = new LoadModule<Constants::floatEVAA>(lagrangeProfile, eulerProfile, car);
    TwoTrackModelParent<Constants::floatEVAA>* TwoTrackModel_obj = new TwoTrackModelBE<Constants::floatEVAA>(car, loadModule);
    ALE<Constants::floatEVAA>* ale = new ALE<Constants::floatEVAA>(car, loadModule, TwoTrackModel_obj);

    size_t solutionDim = Constants::DIM * (size_t)Constants::VEC_DIM;
    Constants::floatEVAA* sol = Math::malloc<Constants::floatEVAA>(solutionDim);

    ale->Solve(sol);

    ale->PrintFinalResults();

	delete TwoTrackModel_obj;
    delete car;
    delete loadModule;
    delete lagrangeProfile;
	delete eulerProfile;
    delete ale;

    Math::free<Constants::floatEVAA>(sol);
}
}  // namespace EVAA
