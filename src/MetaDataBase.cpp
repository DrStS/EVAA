/***********************************************************************************************//**
* \file MetaDataBase.cpp
* This file holds the function definitions of MetaDataBase.
* \date 04/14/2020
**************************************************************************************************/
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>

#include "MetaDataBase.h"


MetaDataBase* MetaDataBase::_database = NULL;

MetaDataBase* MetaDataBase::DataBase() {
    if (!_database)
        _database = new MetaDataBase;
    return _database;
}


/**
* \brief blank private Constructor
*/
MetaDataBase::MetaDataBase() {
    _filename = "";
    _load_filename = "";
}


/**
* \brief setter for the Filename
*/
void MetaDataBase::setFileName(
    const std::string& filename /**< [in] reference to the filename*/
) {
    _filename = filename;
    settings = EVAA_settings(filename, xml_schema::flags::dont_validate);
}

/**
* \brief setter for the load Filename
*/
void MetaDataBase::setloadFileName(
    const std::string& filename /**< [in] reference to the load filename*/
) {
    _load_filename = filename;
    load_data = EVAA_load_module(filename, xml_schema::flags::dont_validate);

}


/**
* \brief reads the vectors of the positions of the legs relative to the center of mass
*/
template<typename T> void MetaDataBase::readVectorLegs(double* storage, T vec){
       readVector(storage, vec.FrontLeft());
       readVector(storage+3,vec.FrontRight());
       readVector(storage+6,vec.ReerLeft());
       readVector(storage+9,vec.ReerRight());
 }

/**
* \brief reads a 3 dim vector
*/
template<typename T> void MetaDataBase::readVector(double* storage, T vec) {
    storage[0] = vec.x();
    storage[1] = vec.y();
    storage[2] = vec.z();
}

/**
* \brief reads paramaters which are given for all the legs 
*/
template<typename T> void MetaDataBase::readLegs(double* storage, T vec) {
    storage[0] = vec.FrontLeft();
    storage[1] = vec.FrontRight();
    storage[2] = vec.ReerLeft();
    storage[3] = vec.ReerRight();
}

/**
* \brief reads quaternion
*/
template<typename T> void MetaDataBase::readangles(double* storage, T vec) {
	storage[0] = vec.x();
	storage[1] = vec.y();
	storage[2] = vec.z();
	storage[3] = vec.w();
}



/**
* \brief read car, initial and simulaiton parameter
*/
void MetaDataBase::ReadParameters(){

    //--------------------------------------------------
    // Load car parameters
    //--------------------------------------------------
 
    mass_body = settings->Vehicle().TwoTrackModel().Mass().Body();
    readLegs(mass_wheel, settings->Vehicle().TwoTrackModel().Mass().UnsprungMass());
    readLegs(mass_tyre, settings->Vehicle().TwoTrackModel().Mass().Tyre());
    I_body[0] = settings->Vehicle().TwoTrackModel().Inertia().XX();
    I_body[1] = settings->Vehicle().TwoTrackModel().Inertia().XY();
    I_body[2] = settings->Vehicle().TwoTrackModel().Inertia().XZ();
    I_body[3] = settings->Vehicle().TwoTrackModel().Inertia().YX();
    I_body[4] = settings->Vehicle().TwoTrackModel().Inertia().YY();
    I_body[5] = settings->Vehicle().TwoTrackModel().Inertia().YZ();
    I_body[6] = settings->Vehicle().TwoTrackModel().Inertia().ZX();
    I_body[7] = settings->Vehicle().TwoTrackModel().Inertia().ZY();
    I_body[8] = settings->Vehicle().TwoTrackModel().Inertia().ZZ();
    if (settings->Vehicle().TwoTrackModel().Stiffness().Constant().present()) {
        std::cout << "Take constant stiffness without lookup table" << std::endl;
        readLegs(k_tyre, settings->Vehicle().TwoTrackModel().Stiffness().Constant().get().Tyre());
        readLegs(k_body, settings->Vehicle().TwoTrackModel().Stiffness().Constant().get().Body());
        _lookup_filename = "NO_FILE_SPECIFIED";
    }
    else {
        _lookup_filename = settings->Vehicle().TwoTrackModel().Stiffness().LookupTable().get().FilePath();
        std::ifstream f(_lookup_filename.c_str());
        if (f.good()) {
            std::cout << "Read lookup table from " << _lookup_filename << std::endl;
        }
        else {
            std::cout << "Lookup table at " << _lookup_filename << " does not exist!" << std::endl;
            exit(2);
        }
    }
    readLegs(c_tyre, settings->Vehicle().TwoTrackModel().DampingCoefficients().Tyre());
    readLegs(c_body, settings->Vehicle().TwoTrackModel().DampingCoefficients().Body());
    readLegs(l_long, settings->Vehicle().TwoTrackModel().Geometry().LongitudinalReferenceToWheel());
    readLegs(l_lat, settings->Vehicle().TwoTrackModel().Geometry().LateralReferenceToWheel());
    readLegs(lower_spring_length, settings->Vehicle().TwoTrackModel().Geometry().SuspensionSprings());
    readLegs(upper_spring_length, settings->Vehicle().TwoTrackModel().Geometry().TyreSprings());


//--------------------------------------------------
// Load initial parameters
//--------------------------------------------------
    readVector(initial_vel_body, settings->InitialConditions().Velocities().Body());
    readVector(initial_ang_vel_body, settings->InitialConditions().Velocities().angularBody());

    readVectorLegs(initial_vel_wheel, settings->InitialConditions().Velocities().UnsprungMass());
    readVectorLegs(initial_vel_tyre, settings->InitialConditions().Velocities().Tyre());

    readLegs(initial_lower_spring_length, settings->InitialConditions().SpringElongation().Tyre());
    readLegs(initial_upper_spring_length, settings->InitialConditions().SpringElongation().Body());

	readVector(initial_pos_body, settings->InitialConditions().Position().Body());
	if (settings->InitialConditions().Position().UnsprungMass().present()) {
		initial_leg_flag = 1;
		readVectorLegs(initial_pos_wheel, settings->InitialConditions().Position().UnsprungMass().get());
		readVectorLegs(initial_pos_tyre, settings->InitialConditions().Position().Tyre().get());
	}

	readangles(initial_angle, settings->InitialConditions().Orientation());


//--------------------------------------------------
// Load simulation parameters
//--------------------------------------------------
    readVector(gravity, settings->SimulationParameters().GeneralSettings().Gravity());
    num_time_iter = settings->SimulationParameters().GeneralSettings().NumberOfIterations();
    timestep = settings->SimulationParameters().GeneralSettings().TimestepSize();
    

    DOF = settings->SimulationParameters().LinearALE().DOF();

    max_num_iter = settings->SimulationParameters().MultyBodyDynamics().MaximalIterationNumber();
    tolerance = settings->SimulationParameters().MultyBodyDynamics().Tolerance();
    solution_dim = settings->SimulationParameters().MultyBodyDynamics().SolutionDimension();

    std::string solver = settings->SimulationParameters().MultyBodyDynamics().Solver();


    if (solver=="explicit_Euler"){
        MBD_solver = EXPLICIT_EULER;
    } else if (solver=="RK4"){
        MBD_solver = RUNGE_KUTTA_4;
    } else if (solver=="Broyden_Euler"){
        MBD_solver = BROYDEN_EULER;
    } else if (solver=="Broyden_CN"){
        MBD_solver = BROYDEN_CN;
    } else if (solver=="Broyden_BDF2"){
        MBD_solver = BROYDEN_BDF2;
    } else {
        std::cerr<<"Wrong MBD-solver specified in the XML! Please type: \n   -explicit_Euler \n   -RK4\n   -Broyden_Euler\n   -Broyden_CN\n   -Broyden_BDF2"<<std::endl;
        exit(2);
    }
}

/**
* \brief read boundary condition and forces on tyres, wheel and body
*/
void MetaDataBase::ReadLoadParameters() {

    //--------------------------------------------------
    // Load external parameters
    //--------------------------------------------------
	std::string boundary_conditions = load_data->boundary_description().BoundaryConditions();
	if (boundary_conditions == "fixed") {
		boundary_condition_road = FIXED;
        std::cout << "Run the simulation with fixed tyres" << std::endl;
	}
	else if (boundary_conditions == "nonfixed") {
		boundary_condition_road = NONFIXED;
        std::cout << "Run the simulation without any tyre constraints" << std::endl;
    }
	else if (boundary_conditions == "circle") {
		boundary_condition_road = CIRCULAR;
        std::cout << "Run the simulation on a circular road" << std::endl;
    }
	else {
		std::cerr << "Wrong boundary conditions! Only circle, fixed and nonfixed implemented so far" << std::endl;
		exit(2);
	}
	if (load_data->boundary_description().circular().present()) {
		profile_radius = load_data->boundary_description().circular()->radius();
		readVector(profile_center, load_data->boundary_description().circular()->center());
	}
	readVectorLegs(external_force_tyre, load_data->forces().force_tyre());
	readVectorLegs(external_force_wheel, load_data->forces().force_wheel());
    readVector(external_force_body, load_data->forces().force_body());
}



/**
* \brief read and order params used for the loookup tables and generate the stiffness and damping loookup table
*/
void MetaDataBase::ReadLookupParameters(
    EVAALookup<Constants::floatEVAA>** lookupStiffness /**< [out] pointer to the pointer to the stiffness lookup from the compute engine*/,
    EVAALookup<Constants::floatEVAA>** lookupDamping /**< [out] pointer to the pointer to the damping lookup from the compute engine*/
    ) {
    if (_lookup_filename == "NO_FILE_SPECIFIED") return;
    lookup_table = LookupHandler(_lookup_filename, xml_schema::flags::dont_validate);
    if (lookup_table->LookupTableGenerator().present()) {
		double* a, * k_body, * k_tyre;
		double b, c, l_min, l_max;
		int size, k, type, order;

		a = new(double[8]);
		k_body = new(double[4]);
		k_tyre = new(double[4]);

		std::cout << "Generate look up table from parameters." << std::endl;

		size = lookup_table->LookupTableGenerator().get().Size();
		b = lookup_table->LookupTableGenerator().get().TableParameters().b();
		c = lookup_table->LookupTableGenerator().get().TableParameters().c();
		l_min = lookup_table->LookupTableGenerator().get().Range().l_min();
		l_max = lookup_table->LookupTableGenerator().get().Range().l_max();
		k = lookup_table->LookupTableGenerator().get().InterpolationMethod().k();
		type = lookup_table->LookupTableGenerator().get().InterpolationMethod().type();
		order = lookup_table->LookupTableGenerator().get().InterpolationMethod().order();

		readLegs(k_body, lookup_table->LookupTableGenerator().get().Magnitude().Body());
		readLegs(k_tyre, lookup_table->LookupTableGenerator().get().Magnitude().Tyre());
		a[0] = k_body[0];
		a[1] = k_tyre[0];
		a[2] = k_body[1];
		a[3] = k_tyre[1];
		a[4] = k_body[2];
		a[5] = k_tyre[2];
		a[6] = k_body[3];
		a[7] = k_tyre[3];

        *lookupStiffness = new EVAALookup<Constants::floatEVAA>(size, a, b, c, l_min, l_max, k, type, order);
		interpolation = 1;  // to switch from constant to interpolation type

        // damping is /100 from the stiffness for the start
        for (auto j = 0; j < k; j++) {
            a[j] /= 100;
        }
        *lookupDamping = new EVAALookup<Constants::floatEVAA>(size, a, b, c, l_min, l_max, k, type, order);

        delete[] a;
		delete[] k_body;
		delete[] k_tyre;
    }
}