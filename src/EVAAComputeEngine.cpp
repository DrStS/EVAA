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

namespace EVAA {
EVAAComputeEngine::EVAAComputeEngine(std::string xmlSimulationFileName, std::string xmlCarFileName, std::string xmlLoadFileName) {
    IO::checkFileExists(xmlSimulationFileName);
    IO::checkFileExists(xmlCarFileName);
    IO::checkFileExists(xmlLoadFileName);

    MetaDatabase<Constants::floatEVAA>::getDatabase().readSimulationParameters(xmlSimulationFileName);
    MetaDatabase<Constants::floatEVAA>::getDatabase().readVehicleParameters(xmlCarFileName);
    MetaDatabase<Constants::floatEVAA>::getDatabase().readLoadParameters(xmlLoadFileName);
}

void EVAAComputeEngine::printInfo(void) {
    Math::PrintMKLInfo();
    auto& db = MetaDatabase<Constants::floatEVAA>::getDatabase();
    std::cout << "\n\nCalculate the solution after " << db.getNumberOfTimeIterations() * db.getTimeStepSize() << "s with dt = " << db.getTimeStepSize() << " for " << db.getNumberOfTimeIterations() << " iterations\n\n\n";
}

void EVAAComputeEngine::computeMKLTwoTrackModelBE(void) {
    auto& db = MetaDatabase<Constants::floatEVAA>::getDatabase();
    if (true) {  // TODO remove this
        Constants::floatEVAA* sol = Math::malloc<Constants::floatEVAA>(Constants::DOF);
        Car<Constants::floatEVAA>* car = new Car<Constants::floatEVAA>();
        Lagrange<Constants::floatEVAA>* lagrange = new Straight<Constants::floatEVAA>();
        Euler<Constants::floatEVAA>* euler = new Fixed<Constants::floatEVAA>(db.getGravityField()[2]);
        LoadModule<Constants::floatEVAA>* load = new LoadModule<Constants::floatEVAA>(lagrange, euler, car);
        TwoTrackModelFull<Constants::floatEVAA> solver(car, load);
        euler->ApplyProfileInitialCondition(car);

        solver.Solve(sol);

        solver.PrintFinalResults(sol);
        
        // solution for 10,000 iterations, dt = 1e-3 => t = 10sec
        static const Constants::floatEVAA refSolStraightFixed[11]{
            -2.857594795257409e-02, 3.103609667343786e-11, -1.781939846572105e-02, 2.606202178817649e-01, 2.830000001221397e-01, 2.606202177281093e-01, 2.829999999648761e-01, 2.642462719927013e-01, 2.830000001833303e-01, 2.642462719576703e-01, 2.830000001517805e-01
        };

        static size_t tmp = 0;
        size_t count = 0;
        for (auto i = 0; i < 11; ++i) {
            if (std::abs(refSolStraightFixed[i] - sol[i]) > 1e-13) {
                count++;
                tmp++;
                std::cout << i << ": " << std::abs(sol[i] - refSolStraightFixed[i]) << "\n";
            }
        }

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
    Constants::floatEVAA* sol = Math::calloc<Constants::floatEVAA>(Constants::MBD_SOLUTION_SIZE);
    MBDMethod<Constants::floatEVAA> solver;

    solver.Solve(sol);
    solver.PrintFinalResult(sol);

#ifdef USE_HDF5
    solver.WriteFinalResult(sol);
    solver.WriteFinalResultFormatted(sol);
#endif  // USE_HDF5

    Math::free<Constants::floatEVAA>(sol);
}


void EVAAComputeEngine::computeALE(void) {    
    
    Lagrange<Constants::floatEVAA>* lagrangeProfile;
    Euler<Constants::floatEVAA>* eulerProfile;
    Car<Constants::floatEVAA>* car = new Car<Constants::floatEVAA>();

    auto& db = MetaDatabase<Constants::floatEVAA>::getDatabase();
    if (db.getLagrangianRoadConditions() == LagrangianBoundaryConditionRoad::CIRCULAR) {
        lagrangeProfile = new Circular<Constants::floatEVAA>(db.getCircularRoadCenter());
    }
    else if (db.getLagrangianRoadConditions() == LagrangianBoundaryConditionRoad::STRAIGHT) {
        lagrangeProfile = new Straight<Constants::floatEVAA>();
    }
    else if (db.getLagrangianRoadConditions() == LagrangianBoundaryConditionRoad::ARBITRARY) {
        lagrangeProfile = new Arbitrary<Constants::floatEVAA>(db.getArbitraryTrajectory());
    }
    else {
        throw "wrong lagrangian road profile condition";
    }

    if (db.getEulerianRoadConditions() == EulerianBoundaryConditionRoad::FIXED) {
        eulerProfile = new Fixed<Constants::floatEVAA>(db.getGravityField()[2]);
    }
    else if (db.getEulerianRoadConditions() == EulerianBoundaryConditionRoad::NONFIXED) {
        eulerProfile = new Nonfixed<Constants::floatEVAA>(db.getGravityField()[2]);
    }
    else if (db.getEulerianRoadConditions() == EulerianBoundaryConditionRoad::SINUSOIDAL) {
        eulerProfile = new Sinusoidal<Constants::floatEVAA>(db.getArbitraryTrajectory(), db.getGravityField()[2]);
    }
    else {
        throw "wrong eulerian road profile condition";
    }
    lagrangeProfile->ApplyProfileInitialCondition(car);
    eulerProfile->ApplyProfileInitialCondition(car);
    LoadModule<Constants::floatEVAA>* loadModule = new LoadModule<Constants::floatEVAA>(lagrangeProfile, eulerProfile, car);
    TwoTrackModelParent<Constants::floatEVAA>* TwoTrackModel_obj = nullptr;
    if (db.getALESolver() == ALESolver::IMPLICIT_EULER) {
        TwoTrackModel_obj = new TwoTrackModelBDF2<Constants::floatEVAA>(car, loadModule);
    }
    else if (db.getALESolver() == ALESolver::BDF2) {
		std::cout << "BDF2 type eulerain solver" << std::endl;
        TwoTrackModel_obj = new TwoTrackModelBDF2<Constants::floatEVAA>(car, loadModule);
    }

    else {
        throw "Only Backward Euler (IMPLICIT_EULER) and BDF2 (BDF2) implemented!";
    }
    
#ifndef USE_HDF5
    ALE<Constants::floatEVAA>* ale = new ALE<Constants::floatEVAA>(car, loadModule, TwoTrackModel_obj);
#else
    std::string fileName = "ALE_Checkpoints.hdf5";
    std::string filePath = "";
    ALE<Constants::floatEVAA>* ale = new ALE<Constants::floatEVAA>(car, loadModule, TwoTrackModel_obj, filePath, fileName);
#endif  // ! USE_HDF5

    size_t solutionDim = Constants::DIM * (size_t)Constants::VEC_DIM;
    Constants::floatEVAA* sol = Math::malloc<Constants::floatEVAA>(solutionDim);

    ale->Solve(sol);

    ale->PrintFinalResults();

#ifndef DAMPING

    static const Constants::floatEVAA referenceSol10000[27]{
        -2.080734004508660e+01, 4.546487545598364e+01, -3.139735789643002e-01,
        -2.277937596041519e+01, 4.642269570993442e+01, -7.291467887135016e-01,
        -2.277937596041519e+01, 4.642269570993442e+01, -9.069000002374381e-01,
        -2.137224395060193e+01, 4.334600552806688e+01, -7.272225898769764e-01,
        -2.137224395060193e+01, 4.334600552806688e+01, -9.069000014894851e-01,
        -2.005484271404052e+01, 4.765684183026740e+01, -7.274584751808845e-01,
        -2.005484271404052e+01, 4.765684183026740e+01, -9.068999977574306e-01,
        -1.865735998970482e+01, 4.460124978395771e+01, -7.261445267890266e-01,
        -1.865735998970482e+01, 4.460124978395771e+01, -9.068999987423678e-01};
#else
    static const Constants::floatEVAA referenceSol10000[27]{-2.080734004508660e+01, 4.546487545598364e+01, -3.115758548353433e-01, -2.278023469533989e+01, 4.642230270340894e+01, -7.272668308078089e-01, -2.278023469533989e+01, 4.642230270340894e+01, -9.068999999925623e-01, -2.137138521567723e+01, 4.334639853459237e+01, -7.312925657635747e-01, -2.137138521567723e+01, 4.334639853459237e+01, -9.069000000016965e-01, -2.005569556026510e+01, 4.765645151875001e+01, -7.242654770287074e-01, -2.005569556026510e+01, 4.765645151875001e+01, -9.068999999177775e-01, -1.865650714348024e+01, 4.460164009547510e+01, -7.270421696854118e-01, -1.865650714348024e+01, 4.460164009547510e+01, -9.069000000166272e-01};
#endif  // ! DAMPING

    static size_t tmp = 0;
    size_t count = 0;
    for (auto i = 0; i < 27; ++i) {
        if (std::abs(referenceSol10000[i] - sol[i]) > 1e-14) {
            count++;
            std::cout << i << ": " << std::abs(sol[i] - referenceSol10000[i]) << "\n";
        }        
    }   
    if (count) tmp++;    
    if (tmp) std::cout<< "\n\n \t # of Wrong solutions" << tmp << " \n";

    delete TwoTrackModel_obj;
    delete car;
    delete loadModule;
    delete lagrangeProfile;
    delete eulerProfile;
    delete ale;

    Math::free<Constants::floatEVAA>(sol);
}
}  // namespace EVAA
