/**
 * \mainpage
 * \section LICENSE
 *  Copyright &copy; 2019, Dr. Stefan Sicklinger, Munich \n
 *  All rights reserved. \n
 *
 *  This file is part of EVAA.
 *
 *  EVAA is free software: you can redistribute it and/or modify \n
 *  it under the terms of the GNU General Public  License as published by \n
 *  the Free Software Foundation, either version 3 of the License, or \n
 *  (at your option) any later version. \n
 *
 *  EVAA is distributed in the hope that it will be useful, \n
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of \n
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the \n
 *  GNU General Public License for more details. \n
 *
 *  You should have received a copy of the GNU General Public License \n
 *  along with EVAA.  If not, see <a
 * href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses/</a>.
 *
 *  Additional permission under GNU GPL version 3 section 7
 *
 * If you modify this Program, or any covered work, by linking or combining it
 * with Intel Math Kernel Libraries(MKL) (or a modified version of that
 * library), containing parts covered by the terms of the license of the MKL,
 * the licensors of this Program grant you additional permission to convey the
 * resulting work.
 *
 * \section DESCRIPTION
 *  This is the main file of EVAA
 *  EVAA: Efficient Vehicle dynAmics simulAtor
 *
 *
 *
 * \section HOWTO
 * Please find all further information on
 * <a href="https://github.com/DrStS/EVAA">EVAA Project</a>
 *
 *
 * <EM> Note: The Makefile suppresses per default all compile and linking
 * command output to the terminal. You may enable this information by make
 * VEREBOSE=1</EM>
 */

/**
 * \file main.cpp
 * This file holds the main function of EVAA.
 * \author Stefan Sicklinger
 * \date 6/11/2019
 * \version alpha
 */

#ifdef EVAA_COMMANDLINE_ON
#include <iostream>
#include <string>
#include <vector>

#include "AuxiliaryParameters.h"
#include "EVAAComputeEngine.h"
#include "Timer.h"
#endif  // EVAA_COMMANDLINE_ON

#ifndef EVAA_COMMANDLINE_ON
#include "EVAAMainWindow.h"
#endif  // EVAA_COMMANDLINE_ON

int main(int argc, char **argv) {
#ifdef EVAA_COMMANDLINE_ON
    std::cout << "Hello EVAA is fired up!" << std::endl;
    std::cout << "GIT: " << EVAA::AuxiliaryParameters::gitSHA1 << std::endl;
    std::vector<std::string> allArgs(argv, argv + argc);
    for (std::vector<std::string>::iterator it = allArgs.begin(); it != allArgs.end(); ++it) {
        std::cout << *it << std::endl;
    }

    /** Define the paths to xml configuration files */
    std::string simulationParametersFileNameXML;
    std::string carSettingsFileNameXML;
    std::string loadProfileFileNameXML;

    if (allArgs.size() > 3) {
        /** Take arguments from command line */
        simulationParametersFileNameXML = allArgs[1];
        carSettingsFileNameXML = allArgs[2];
        loadProfileFileNameXML = allArgs[3];
    }
    else {
        /** Hardcoded paths*/

        /** Hardcoded paths for car settings */
        const std::string carWithInterpolationFileNameXML =
            "C:\\software\\repos\\EVAA\\inputFiles\\CarWithInterpolation.xml";
        const std::string carConstantStiffnessFileNameXML =
            "C:\\software\\repos\\EVAA\\inputFiles\\CarConstantStiffness.xml";
#ifdef INTERPOLATION
        carSettingsFileNameXML = carWithInterpolationFileNameXML;
#else
        carSettingsFileNameXML = carConstantStiffnessFileNameXML;
#endif
        /** Hardcoded path for simulation parameters*/
        simulationParametersFileNameXML =
            "C:\\software\\repos\\EVAA\\inputFiles\\SimulationParameters.xml";
        /** Hardcoded paths for load profiles*/
        const std::string loadArbitraryFileNameXML =
            "C:\\software\\repos\\EVAA\\inputFiles\\LoadArbitraryCar.xml";
        const std::string loadCircularFileNameXML =
            "C:\\software\\repos\\EVAA\\inputFiles\\LoadCircularCar.xml";
        const std::string loadStraightFileNameXML =
            "C:\\software\\repos\\EVAA\\inputFiles\\LoadStraightCar.xml";
        loadProfileFileNameXML = loadStraightFileNameXML;
    }

    /** Construct the car*/
    EVAA::EVAAComputeEngine *myComputeEngine = new EVAA::EVAAComputeEngine(
        simulationParametersFileNameXML, carSettingsFileNameXML, loadProfileFileNameXML);
    myComputeEngine->printInfo();

    auto &timer1 = EVAA::anaysisTimer01;
    constexpr size_t numIterations = 1;

    std::cout << std::defaultfloat;
    unsigned long timeMBD = 0.;    
    for (auto i = 0; i < numIterations; ++i) {
        timer1.start();
        myComputeEngine->computeMBD();
        timer1.stop();
        timeMBD += timer1.getDurationMilliSec();
    }
    std::cout << "It took " << std::defaultfloat << timeMBD / numIterations << " ms to run the solver(computeMBD).\n\n\n" << std::endl;
    
    unsigned long time11DOF = 0.;
    for (auto i = 0; i < numIterations; ++i) {
        timer1.start();
        myComputeEngine->computeMKLTwoTrackModelBE();
        timer1.stop();
        time11DOF += timer1.getDurationMilliSec();
    }
    std::cout << "It took " << std::defaultfloat << time11DOF / numIterations << " ms to run the solver 11dofBE.\n\n\n" << std::endl;

    unsigned long timeALE = 0.;
    for (auto i = 0; i < numIterations; ++i) {
        timer1.start();
        myComputeEngine->computeALE();
        timer1.stop();
        timeALE += timer1.getDurationMilliSec();
    }
    std::cout << "It took " << std::defaultfloat << timeALE / numIterations << " ms to run the solver(computeALE).\n\n\n" << std::endl;

    delete myComputeEngine;

    std::cout << "\nWe did a great job! Awesome!" << std::endl;

#endif  // EVAA_COMMANDLINE_ON

#ifndef EVAA_COMMANDLINE_ON
#ifdef __linux
    putenv((char *)"MESA_GL_VERSION_OVERRIDE=3.2");
    // Fixes decimal point issue in vtkSTLReader
    putenv((char *)"LC_NUMERIC=C");
#endif  // LINUX

    EVAAMainWindow(argc, argv);
#endif  // EVAA_COMMANDLINE_ON

    return 0;
}
