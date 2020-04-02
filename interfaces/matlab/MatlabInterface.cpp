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
#include "Timer.h"
#include "AuxiliaryParameters.h"
#include "mex.hpp"
#include "mexAdapter.hpp"

using namespace matlab::data;
using matlab::mex::ArgumentList;

class MexFunction : public matlab::mex::Function {
public:
	void operator()(ArgumentList outputs, ArgumentList inputs) {
		checkArguments(outputs, inputs);

		std::cout << "Hello EVAA is fired up!" << std::endl;
		std::cout << "GIT: " << EVAA::AuxiliaryParameters::gitSHA1 << std::endl;
		matlab::data::CharArray inputFile = inputs[0];
		std::cout << "EVAA configuration file is: " << inputFile.toAscii() << std::endl;
		EVAAComputeEngine* myComputeEngine = new EVAAComputeEngine(inputFile.toAscii());
		myComputeEngine->prepare();
		anaysisTimer01.start();
		myComputeEngine->compute();
		anaysisTimer01.stop();
		std::cout << "It took " << anaysisTimer01.getDurationMilliSec() << " ms to run the solver." << std::endl;
		myComputeEngine->clean();
		std::cout << "We did a great job! Awesome!" << std::endl;
	}

	void checkArguments(ArgumentList outputs, ArgumentList inputs) {
		// Get pointer to engine
		std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();

		// Get array factory
		ArrayFactory factory;

		// Check first input argument
		if (inputs[0].getType() != ArrayType::CHAR )
		{

			matlabPtr->feval(u"error",
				0,
				std::vector<Array>({ factory.createScalar("First input must be string") }));
		}

		// Check number of outputs
		if (outputs.size() > 0) {
			matlabPtr->feval(u"error",
				0,
				std::vector<Array>({ factory.createScalar("No return value") }));
		}
	}
};