/*           This file has been prepared for Doxygen automatic documentation generation.          */
/***********************************************************************************************//**
 * \mainpage
 * \section LICENSE
 *  Copyright &copy; 2019, Dr. Stefan Sicklinger, Munich \n
 *  All rights reserved. \n
 *
 *  This file is part of EVAA.
 *
 *  EVAA is free software: you can redistribute it and/or modify \n
 *  it under the terms of the GNU General Public License as published by \n
 *  the Free Software Foundation, either version 3 of the License, or \n
 *  (at your option) any later version. \n
 *
 *  EVAA is distributed in the hope that it will be useful, \n
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of \n
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the \n
 *  GNU General Public License for more details. \n
 *
 *  You should have received a copy of the GNU General Public License \n
 *  along with EVAA.  If not, see <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses/</a>.
 *
 *  Additional permission under GNU GPL version 3 section 7
 *
 * If you modify this Program, or any covered work, by linking or combining it with Intel Math Kernel Libraries(MKL) 
 * (or a modified version of that library), containing parts covered by the terms of the license of the MKL, 
 * the licensors of this Program grant you additional permission to convey the resulting work. 
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
 * <EM> Note: The Makefile suppresses per default all compile and linking command output to the terminal.
 *       You may enable this information by make VEREBOSE=1</EM>
 *
 *
 *
 **************************************************************************************************/
/***********************************************************************************************//**
 * \file main.cpp
 * This file holds the main function of EVAA.
 * \author Stefan Sicklinger
 * \date 6/11/2019
 * \version alpha
 **************************************************************************************************/
#ifdef EVAA_COMMANDLINE_ON
   #include <iostream>
   #include <string>
   #include <vector>
   #include "AuxiliaryParameters.h"
   #include "EVAAComputeEngine.h"
   #include "Timer.h"
#endif // EVAA_COMMANDLINE_ON

#ifndef EVAA_COMMANDLINE_ON
#include "EVAAMainWindow.h"
#endif // EVAA_COMMANDLINE_ON

int main(int argc, char **argv) {
#ifdef EVAA_COMMANDLINE_ON
	std::cout << "Hello EVAA is fired up!" << std::endl;
	std::cout << "GIT: " << EVAA::AuxiliaryParameters::gitSHA1 << std::endl;
	std::vector<std::string> allArgs(argv, argv + argc);
	for (std::vector<std::string>::iterator it = allArgs.begin(); it != allArgs.end(); ++it) {
		 std::cout << *it << std::endl;
	}
	//allArgs[0] = "EVAA.exe"
	if (allArgs.size()>1) {
		EVAAComputeEngine* myComputeEngine = new EVAAComputeEngine(allArgs[1]);
	}
	EVAAComputeEngine* myComputeEngine = new EVAAComputeEngine("default");
		myComputeEngine->prepare();
		/*anaysisTimer01.start();
		myComputeEngine->compute();
		anaysisTimer01.stop();
		std::cout << "It took " << anaysisTimer01.getDurationMilliSec() << " ms to run the solver(MKL)." << std::endl;
		myComputeEngine->clean();*/
		/*anaysisTimer01.start();
		myComputeEngine->computeEigen();
		anaysisTimer01.stop();
		std::cout << "It took " << anaysisTimer01.getDurationMilliSec() << " ms to run the solver(Eigen)." << std::endl;*/
		anaysisTimer01.start();
		myComputeEngine->computeBlaze11DOF();
		anaysisTimer01.stop();
		std::cout << "It took " << anaysisTimer01.getDurationMilliSec() << " ms to run the solver(Blaze)." << std::endl;
	std::cout << "We did a great job! Awesome!" << std::endl;
#endif // EVAA_COMMANDLINE_ON
#ifndef EVAA_COMMANDLINE_ON
#ifdef __linux
	putenv((char *)"MESA_GL_VERSION_OVERRIDE=3.2");
	// Fixes decimal point issue in vtkSTLReader
	putenv((char *)"LC_NUMERIC=C");
#endif //LINUX

	EVAAMainWindow(argc, argv);
#endif // EVAA_COMMANDLINE_ON


}