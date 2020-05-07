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
#include "EVAAComputeEngine.h"
#include "gtest/gtest.h"
#include "../helper/CSVReader.h"

#ifdef INTERPOLATION
#ifdef WRITECSV
class TestTwoTrackModel: public ::testing::Test {
public:
EVAA::EVAAComputeEngine* engine;
 
virtual void SetUp(){
	engine = new EVAA::EVAAComputeEngine("C:\\software\\repos\\EVAA\\unitTests\\testXML\\car_TwoTrackTest.xml", "C:\\software\\repos\\EVAA\\unitTests\\testXML\\load_TwoTrackTest.xml");
} 

virtual void TearDown(){ 
	delete engine;
} 
};

TEST_F(TestTwoTrackModel, TwoTrackModelNonFixed){
	engine->computeMKLTwoTrackModelBE();
	std::ifstream newFile("C:\\software\\repos\\EVAA\\output\\newtonOutput.txt");
	std::ifstream refFile("C:\\software\\repos\\EVAA\\unitTests\\refCSV\\Newton11DOF_lookup_noDamping_dt_1e-3_tol_1e-8.txt");

	CSVRow rowNew, rowRef;
    while(refFile >> rowRef)
    {
		newFile >> rowNew;
		EXPECT_EQ(rowRef[0], rowNew[0]);
		EXPECT_NEAR(::atof(rowRef[1].c_str()), ::atof(rowNew[1].c_str()), EVAA::AuxiliaryParameters::machineEpsilon*3 );
    }
}
#endif // WRITECSV
#endif // INTERPOLATION