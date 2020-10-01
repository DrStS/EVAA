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

 /**
  * \file EVAAComputeEngine.h
  * This file holds the class of ComputeEngine.
  * \date 6/13/2019
  */

#pragma once
#include <string>
namespace EVAA {
	/**
	 * \brief Class ComputeEngine the core of EVAA
	**/
	class EVAAComputeEngine {
	public:
		/**
		 * \brief Reads parameters from XML files.
		 * \param[in] file name of xml file to call singelton constructor of
		 * \author Stefan Sicklinger
		**/
		EVAAComputeEngine(std::string t_inputFileName);
		/**
		 * \brief Test hdf5 output demo write
		 * \param[in] file name of hdf5 file
		 * \author Stefan Sicklinger
		**/
		void testHDF5(std::string t_outputFileName);
	};

}  // namespace EVAA
