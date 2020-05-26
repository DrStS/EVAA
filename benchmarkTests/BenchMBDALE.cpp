#pragma once 

#include <benchmark/benchmark.h>

//#include "../../src/EVAAComputeEngine.h"
//
//using namespace EVAA;
//
//// car
//std::string carWithInterpolationFileNameXML =
//    "C:\\software\\repos\\EVAA\\inputFiles\\CarWithInterpolation.xml";
//std::string carConstantStiffnessFileNameXML =
//    "C:\\software\\repos\\EVAA\\inputFiles\\CarConstantStiffness.xml";
//
///** Hardcoded paths for load profiles*/
//std::string loadArbitraryFileNameXML =
//    "C:\\software\\repos\\EVAA\\inputFiles\\LoadArbitraryCar.xml";
//std::string loadCircularFileNameXML =
//    "C:\\software\\repos\\EVAA\\inputFiles\\LoadCircularCar.xml";
//std::string loadStraightFileNameXML =
//    "C:\\software\\repos\\EVAA\\inputFiles\\LoadStraightCar.xml";
//
///** Hardcoded path for simulation parameters*/
//std::string simulationParametersFileNameXML =
//    "C:\\software\\repos\\EVAA\\inputFiles\\SimulationParameters.xml";
//
//static void BM_MBD(benchmark::State& state) {
//    std::string loadProfileFileNameXML = loadCircularFileNameXML;
//    std::string simulationParametersFileNameXML = "SimulationParameters.xml";
//#ifdef INTERPOLATION
//    std::string carSettingsFileNameXML = carWithInterpolationFileNameXML;
//#else
//    std::string carSettingsFileNameXML = carConstantStiffnessFileNameXML;
//#endif
//    EVAA::EVAAComputeEngine* myComputeEngine = new EVAA::EVAAComputeEngine(
//        simulationParametersFileNameXML, carSettingsFileNameXML, loadProfileFileNameXML);
//
//    auto& db = MetaDatabase<Constants::floatEVAA>::getDatabase();
//    db.readSimulationParameters(simulationParametersFileNameXML);
//    
//    for (auto _ : state) {
//        // ...-
//    }
//    delete myComputeEngine;
//}
//BENCHMARK(BM_MBD)->UseRealTime()->Arg(10)->Arg(15)->Arg(20);