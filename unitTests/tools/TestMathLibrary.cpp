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
#include "MathLibrary.h"
#include "MKLRoutines.h"
#include "gtest/gtest.h"


class TestMathLibrary: public ::testing::Test {
private:

protected: 
virtual void SetUp(){ 
} 

virtual void TearDown(){ 
} 
};
TEST_F(TestMathLibrary, solverLU3x3a){
	double k_1 = 1;
	double k_2 = 2;
	double k_3 = 3;
	std::vector<double> myK3x3(9);
	std::vector<double> myRhs3(3);
	myK3x3[0] = k_1 + k_2;
	myK3x3[1] = -k_2;
	myK3x3[2] = 0.0;
	myK3x3[3] = -k_2;
	myK3x3[4] = k_2;
	myK3x3[5] = 0.0;
	myK3x3[6] = 0.0;
	myK3x3[7] = 0.0;
	myK3x3[8] = k_3;
	myRhs3[0] = 1.0;
	myRhs3[1] = 1.0;
	myRhs3[2] = 1.0;
	std::vector<int> pivot(3);
	// LU Decomposition
	EVAA::Math::ComputeDenseSymLUFactorisation(3, myK3x3, pivot);
	// Solve system
    EVAA::Math::ComputeDenseSymSolution(3, myK3x3, pivot, myRhs3);

	std::vector<double> resultReference = { 2.,5. / 2.,1. / 3. };
	double absError = fabs(resultReference[0] - myRhs3[0]) + fabs(resultReference[1] - myRhs3[1]) + fabs(resultReference[2] - myRhs3[2]);
	EXPECT_NEAR(0.0, absError, EVAA::AuxiliaryParameters::machineEpsilon*3 );
	std::cout << "Absolut error is: " << absError << " versus machineEpsilon: "<< EVAA::AuxiliaryParameters::machineEpsilon << std::scientific << std::endl;
}

TEST_F(TestMathLibrary, solverLU3x3b) {
	double k_1 = 1;
	double k_2 = 2;
	double k_3 = 3;
	std::vector<double> myK3x3(9);
	std::vector<double> myRhs3(3);
	myK3x3[0] = k_1 + k_2;
	myK3x3[1] = -k_2;
	myK3x3[2] = 0.0;
	myK3x3[3] = -k_2;
	myK3x3[4] = k_2;
	myK3x3[5] = 0.0;
	myK3x3[6] = 0.0;
	myK3x3[7] = 0.0;
	myK3x3[8] = k_3;
	myRhs3[0] = 2.0;
	myRhs3[1] = 2.0;
	myRhs3[2] = 2.0;
	std::vector<int> pivot(3);
	// LU Decomposition
    EVAA::Math::ComputeDenseSymLUFactorisation(3, myK3x3, pivot);
	// Solve system
    EVAA::Math::ComputeDenseSymSolution(3, myK3x3, pivot, myRhs3);

	std::vector<double> resultReference = { 2.,5. / 2.,1. / 3. };
	double absError = fabs(2.0*resultReference[0] - myRhs3[0]) + fabs(2.0*resultReference[1] - myRhs3[1]) + fabs(2.0*resultReference[2] - myRhs3[2]);
	EXPECT_NEAR(0.0, absError, EVAA::AuxiliaryParameters::machineEpsilon * 8);
	std::cout << "Absolut error is: " << absError << " versus machineEpsilon: " << EVAA::AuxiliaryParameters::machineEpsilon << std::scientific << std::endl;
}



	