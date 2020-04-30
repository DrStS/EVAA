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
/***********************************************************************************************/ /**
                                                                                                   * \file EVAAComputeEngine.h
                                                                                                   * This file holds the class of ComputeEngine.
                                                                                                   * \date 6/13/2019
                                                                                                   **************************************************************************************************/
#pragma once
#include <mkl.h>

#include <chrono>
#include <string>

#include "11DOF.h"
#include "ALE.h"
#ifndef U_Lookup
#define U_Lookup
#include "EVAALookup.h"
#endif
#include "MBDMethod.h"
#include "MathLibrary.h"

/********/ /**
            * \brief Class ComputeEngine the core of STACCATO
            ***********/
class EVAAComputeEngine {
public:
    /***********************************************************************************************
     * \brief Reads parameters from XML files.
     * \param[in] file name of xml file to call singelton constructor of metadatabase
     * \author Stefan Sicklinger
     ***********/
    EVAAComputeEngine(std::string xmlCarFileName, std::string xmlLoadFileName);
    /***********************************************************************************************
     * \brief Destructor
     * \author Stefan Sicklinger
     ***********/
    ~EVAAComputeEngine();
    /***********************************************************************************************
     * \brief print info
     * \author Stefan Sicklinger
     ***********/
    void printInfo(void);
    /***********************************************************************************************
     * \brief compute engine
     * \author Stefan Sicklinger
     ***********/
    // EIGEN
    void computeEigen11DOF(void);
    /***********************************************************************************************
     * \brief compute engine
     * \author Stefan Sicklinger
     ***********/
    // BLAZE
    void computeBlaze11DOF(void);
    /***********************************************************************************************
     * \brief compute engine
     * \author Stefan Sicklinger
     ***********/
    // MKL
    void computeMKL11DOF(void);
    /***********************************************************************************************
     * \brief compute engine
     * \author Stefan Sicklinger
     ***********/
    // MKL Linear solver
    void computeMKLTwoTrackModel();
    /***********************************************************************************************
     * \brief compute engine
     * \author Stefan Sicklinger
     ***********/

    /**
     * \brief compute 11Dof system with BE to prove correctness
     */
    void computeMKLTwoTrackModelBE();

    /*
    Solve the MBD system
    */
    void computeMBD(void);

    /*
    Fancy solver ALE
    */
    void computeALE(void);

    // just for testing purposes
    void computeALEtest(void);
    void compare_ALE_MBD(void){};

    // MKL
    /*	void computeMKLNasa(void);
            void computeMKLNasa_example(void);
    */
    /***********************************************************************************************
     * \brief clean compute engine free memory
     * \author Stefan Sicklinger
     ***********/
    void clean(void);

private:
    std::string _xmlCarFileName;
    std::string _xmlLoadFileName;
    EVAALookup<Constants::floatEVAA>* _lookupStiffness;
    EVAALookup<Constants::floatEVAA>* _lookupDamping;
};
