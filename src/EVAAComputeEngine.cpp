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
#include "H5Cpp.h"
namespace EVAA {
/** Constructor to instantiate the singletone and create the engine.*/
	EVAAComputeEngine::EVAAComputeEngine(std::string inputFileName) {
	}

	void EVAAComputeEngine::testHDF5(std::string outputFileName) {
		try
		{

			// Turn off the auto-printing when failure occurs so that we can
			// handle the errors appropriately.

			H5::Exception::dontPrint();

			// Create a new file using default properties.

			H5::H5File file(outputFileName, H5F_ACC_TRUNC);

			// Create group "MyGroup" in the root group using an absolute name.

			H5::Group group1(file.createGroup("/MyGroup"));

			// Create group "Group_A" in group "MyGroup" using an
			// absolute name.

			H5::Group group2(file.createGroup("/MyGroup/Group_A"));

			// Create group "Group_B" in group "MyGroup" using a
			// relative name.

			H5::Group group3(group1.createGroup("Group_B"));

			// Close the groups and file.

			group1.close();
			group2.close();
			group3.close();
			file.close();

		} // end of try block

		// catch failure caused by the File operations
		catch (H5::FileIException error)
		{
			error.printErrorStack();
		}

		// catch failure caused by the Group operations
		catch (H5::GroupIException error)
		{
			error.printErrorStack();
		}
	}
}  // namespace EVAA




