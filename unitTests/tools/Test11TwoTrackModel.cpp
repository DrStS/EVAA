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
#include "../helper/CSVReader.h"
#include "EVAAComputeEngine.h"
#include "gtest/gtest.h"

#ifdef INTERPOLATION

namespace EVAA {

class TestTwoTrackModel : public ::testing::Test {
public:
    /** Hardcoded path for simulation parameters*/
    std::string simulationParametersFileNameXML =
        "C:\\software\\repos\\EVAA\\inputFiles\\SimulationParameters.xml";
    std::string loadProfileFileNameXML;

    /** Hardcoded paths for car settings */
    const std::string carWithInterpolationFileNameXML =  //
        "C:\\software\\repos\\EVAA\\inputFiles\\CarWithInterpolation.xml";
    const std::string carConstantStiffnessFileNameXML =  //
        "C:\\software\\repos\\EVAA\\inputFiles\\CarConstantStiffness.xml";
#ifdef INTERPOLATION
    std::string carSettingsFileNameXML = carWithInterpolationFileNameXML;
#else
    std::string carSettingsFileNameXML = carConstantStiffnessFileNameXML;
#endif
    /** Hardcoded paths for load profiles*/
    const std::string loadArbitraryFileNameXML =
        "C:\\software\\repos\\EVAA\\inputFiles\\LoadArbitraryCar.xml";
    const std::string loadCircularFileNameXML =
        "C:\\software\\repos\\EVAA\\inputFiles\\LoadCircularCar.xml";
    const std::string loadStraightFileNameXML =
        "C:\\software\\repos\\EVAA\\inputFiles\\LoadStraightCar.xml";

    EVAAComputeEngine* engine = nullptr;

public:    

    virtual void SetUp() {}
    virtual void TearDown() { delete engine; }
};

class TestTwoTrackModelArbitrary : public TestTwoTrackModel {
public:
    virtual void SetUp() {
        loadProfileFileNameXML = loadArbitraryFileNameXML;
        engine = new EVAAComputeEngine(simulationParametersFileNameXML, carSettingsFileNameXML,
                                       loadProfileFileNameXML);
    }
};

class TestTwoTrackModelCircular : public TestTwoTrackModel {
public:
    virtual void SetUp() {
        loadProfileFileNameXML = loadCircularFileNameXML;
        engine = new EVAAComputeEngine(simulationParametersFileNameXML, carSettingsFileNameXML,
                                       loadProfileFileNameXML);
    }
};

class TestTwoTrackModelStraight : public TestTwoTrackModel {
public:
    virtual void SetUp() {
        loadProfileFileNameXML = loadStraightFileNameXML;
        engine = new EVAAComputeEngine(simulationParametersFileNameXML, carSettingsFileNameXML,
                                       loadProfileFileNameXML);
    }
};

//TEST_F(TestTwoTrackModelStraight, Straight_Fixed_10000it) {    
//    auto& db = MetaDatabase<Constants::floatEVAA>::getDatabase();
//
//    Car<Constants::floatEVAA>* car = new Car<Constants::floatEVAA>();      
//    Lagrange<Constants::floatEVAA>* lagrange = new Straight<Constants::floatEVAA>();
//    Euler<Constants::floatEVAA>*  euler = new Fixed<Constants::floatEVAA>(db.getGravityField()[2]);
//    LoadModule<Constants::floatEVAA>* load =
//        new LoadModule<Constants::floatEVAA>(lagrange, euler, car);
//    // TwoTrackModelFull<Constants::floatEVAA> solver(car, load);
//
//    //euler->ApplyProfileInitialCondition(car);
//
//    Constants::floatEVAA* sol = Math::malloc<Constants::floatEVAA>(Constants::DOF);
//    //solver.Solve(sol);
//
//    //// solution for 10,000 iterations, dt = 1e-3 => t = 10sec
//    //static const Constants::floatEVAA refSol[11]{
//    //    -2.857594795257409e-02, 3.103609667343786e-11, -1.781939846572105e-02,
//    //    2.606202178817649e-01,  2.830000001221397e-01, 2.606202177281093e-01,
//    //    2.829999999648761e-01,  2.642462719927013e-01, 2.830000001833303e-01,
//    //    2.642462719576703e-01,  2.830000001517805e-01};
//
//    //for (auto i = 0; i < 11; ++i) {
//    //    EXPECT_NEAR(sol[i], refSol[i], 1e-13);
//    //}
//
//    delete car;
//    delete lagrange;
//    delete euler;
//    delete load;
//    free(sol);
//}


#ifdef WRITECSV
TEST_F(TestTwoTrackModel, TwoTrackModelNonFixed) {
    engine->computeMKLTwoTrackModelBE();
    std::ifstream newFile("C:\\software\\repos\\EVAA\\output\\newtonOutput.txt");
    std::ifstream refFile(
        "C:\\software\\repos\\EVAA\\unitTests\\refCSV\\Newton11DOF_lookup_noDamping_dt_1e-3_tol_1e-"
        "8.txt");

    CSVRow rowNew, rowRef;
    while (refFile >> rowRef) {
        newFile >> rowNew;
        EXPECT_EQ(rowRef[0], rowNew[0]);
        EXPECT_NEAR(::atof(rowRef[1].c_str()), ::atof(rowNew[1].c_str()), 1e-12);
    }
}
#endif  // WRITECSV

#endif  // INTERPOLATION
}
