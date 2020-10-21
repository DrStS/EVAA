/**
 * \mainpage
 * \section LICENSE
 *  Copyright &copy; 2020, Dr. Stefan Sicklinger, Munich \n
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
 * \section HOWTO
 * Please find all further information on
 * <a href="https://github.com/DrStS/EVAA">EVAA Project</a>
 */

 /**
  * \file main.cpp
  * This file holds the main function of EVAA.
  * \author Stefan Sicklinger
  * \date 09/08/2020
  * \version 1.0
 **/
#include <iostream>
#include <vector>

#include "AuxiliaryParameters.h"
#include "EVAAComputeEngine.h"
#include "Message.h"
#include "InputSchemaEVAA.h"



static void showHelp(void) {
	LOG_INFO << "Valid command line options:" << std::endl;
	LOG_INFO << "\t-h,--help\t\t\tShow this help message." << std::endl;
	LOG_INFO << "\t-i,--inputfile INPUT_FILE\tSpecify the file name to the input file for the run." << std::endl;
}

int main(int argc, char** argv) {
	EVAA::Message myMessage;
	myMessage.initLogging();
	LOG_INFO << "EVAA version hash: " << EVAA::AuxiliaryParameters::gitSHA1 << std::endl;
	std::vector<std::string> allArgs(argv, argv + argc);
	std::string inputFile;
	for (std::vector<std::string>::iterator itArg = allArgs.begin() + 1; itArg != allArgs.end(); ++itArg) {
		if ((*itArg == "-h") || (*itArg == "--help")) {
			showHelp();
			return 0;
		}
		else if ((*itArg == "-i") || (*itArg == "--inputfile")) {
			if ((itArg - allArgs.begin() + 1) < allArgs.size()) {
				inputFile = *(itArg + 1);
			}
			else {
				LOG_ERROR << "No iput file selected: use --help for more details." << std::endl;
				return 1;
			}
		}
		else if (!((*(itArg - 1) == "-i") || (*(itArg - 1) == "--inputfile"))){
				LOG_WARNING << "Unknown option " << *itArg << ": use --help for more details." << std::endl;
		}
	}
	if (!inputFile.empty()) {
		LOG_INFO << "Hello EVAA is fired up with input file: " << inputFile  << std::endl;
	}
	else {
		LOG_ERROR << "No iput file selected: use --help for more details." << std::endl;
		return 1;
	}
	EVAA::EVAAComputeEngine myComputeEngine(inputFile);
	myMessage.writeASCIIArt();
	LOG_DEBUG << "Start HDF5 test" << std::endl;
	#ifdef USE_HDF5
	myComputeEngine.testHDF5("testOutput.h5");
	#endif

	try
	{
		const auto parsedInputFile = InputFileEVAA("inputFileLinear.xml");
	}
	catch (const xml_schema::exception& e)
	{
		std::cerr << e << std::endl;
	}

	return 0;
}
