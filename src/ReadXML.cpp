/***********************************************************************************************//**
* \file ReadXML.cpp
* This file holds the function definitions of ReadXML.
* \date 04/14/2020
**************************************************************************************************/
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>

#include "ReadXML.h"


/**
* \brief reads the vectors of the positions of the legs relative to the center of mass
*/
template<typename T> void ReadXML::readVectorLegs(double* storage, T vec){
       readVector(storage, vec.ReerRight());
       readVector(storage+3,vec.ReerLeft());
       readVector(storage+6,vec.FrontLeft());
       readVector(storage+9,vec.FrontRight());
 }

/**
* \brief reads a 3 dim vector
*/
template<typename T> void ReadXML::readVector(double* storage, T vec) {
    storage[0] = vec.x();
    storage[1] = vec.y();
    storage[2] = vec.z();
}

/**
* \brief reads paramaters which are given for all the legs 
*/
template<typename T> void ReadXML::readLegs(double* storage, T vec) {
    storage[0] = vec.ReerRight();
    storage[1] = vec.ReerLeft();
    storage[2] = vec.FrontLeft();
    storage[3] = vec.FrontRight();
}

/**
* \brief reads quaternion?
*/
template<typename T> void ReadXML::readangles(double* storage, T vec) {
	storage[0] = vec.x();
	storage[1] = vec.y();
	storage[2] = vec.z();
	storage[3] = vec.w();
}


/**
* \brief blank Constructor
*/
ReadXML::ReadXML(){
    _filename = "";
	_load_filename = "";
}

/**
* \brief Constructor
*/
ReadXML::ReadXML(
    const std::string& load_filename /**< [in] reference to the filename*/
) :
    _load_filename(load_filename),
    load_data(EVAA_load_module(load_filename, xml_schema::flags::dont_validate)) {
}

/**
* \brief Constructor
*/
ReadXML::ReadXML(
    const std::string & filename /**< [in] reference to the filename*/,
    const std::string & load_filename /**< [in] reference to the load filename*/
) :_filename(filename), _load_filename(load_filename),
	settings(EVAA_settings(filename, xml_schema::flags::dont_validate)), 
    load_data(EVAA_load_module(load_filename, xml_schema::flags::dont_validate)) {
}

/**
* \brief setter for the Filename
*/
void ReadXML::setFileName(
    const std::string & filename /**< [in] reference to the filename*/
){
    _filename = filename;
}

/**
* \brief setter for the load Filename
*/
void ReadXML::setloadFileName(
    const std::string & filename /**< [in] reference to the load filename*/
) {
	_load_filename = filename;
}

/**
* \brief read car, initial and simulaiton parameter
*/
void ReadXML::ReadParameters(Simulation_Parameters & parameters /**< [out] reference to parameter to where everything gets stored*/){

    //--------------------------------------------------
    // Load car parameters
    //--------------------------------------------------
 
    parameters.mass_body = settings->Vehicle().TwoTrackModel().Mass().Body();
    readLegs(parameters.mass_wheel, settings->Vehicle().TwoTrackModel().Mass().UnsprungMass());
    readLegs(parameters.mass_tyre, settings->Vehicle().TwoTrackModel().Mass().Tyre());
    parameters.I_body[0] = settings->Vehicle().TwoTrackModel().Inertia().XX();
    parameters.I_body[1] = settings->Vehicle().TwoTrackModel().Inertia().YY();
    parameters.I_body[2] = settings->Vehicle().TwoTrackModel().Inertia().ZZ();
    parameters.I_body[3] = settings->Vehicle().TwoTrackModel().Inertia().XY();
    parameters.I_body[4] = settings->Vehicle().TwoTrackModel().Inertia().XZ();
    parameters.I_body[5] = settings->Vehicle().TwoTrackModel().Inertia().YZ();
    if (settings->Vehicle().TwoTrackModel().Stiffness().Constant().present()) {
        std::cout << "Take constant stiffness without lookup table" << std::endl;
        readLegs(parameters.k_tyre, settings->Vehicle().TwoTrackModel().Stiffness().Constant().get().Tyre());
        readLegs(parameters.k_body, settings->Vehicle().TwoTrackModel().Stiffness().Constant().get().Body());
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
    readLegs(parameters.c_tyre, settings->Vehicle().TwoTrackModel().DampingCoefficients().Tyre());
    readLegs(parameters.c_body, settings->Vehicle().TwoTrackModel().DampingCoefficients().Body());
    readLegs(parameters.l_long, settings->Vehicle().TwoTrackModel().Geometry().LongitudinalReferenceToWheel());
    readLegs(parameters.l_lat, settings->Vehicle().TwoTrackModel().Geometry().LateralReferenceToWheel());
    readLegs(parameters.lower_spring_length, settings->Vehicle().TwoTrackModel().Geometry().SuspensionSprings());
    readLegs(parameters.upper_spring_length, settings->Vehicle().TwoTrackModel().Geometry().TyreSprings());


//--------------------------------------------------
// Load initial parameters
//--------------------------------------------------
    readVector(parameters.initial_vel_body, settings->InitialConditions().Velocities().Body());
    readVector(parameters.initial_ang_vel_body, settings->InitialConditions().Velocities().angularBody());

    readVectorLegs(parameters.initial_vel_wheel, settings->InitialConditions().Velocities().UnsprungMass());
    readVectorLegs(parameters.initial_vel_tyre, settings->InitialConditions().Velocities().Tyre());

    readLegs(parameters.initial_lower_spring_length, settings->InitialConditions().SpringElongation().Tyre());
    readLegs(parameters.initial_upper_spring_length, settings->InitialConditions().SpringElongation().Body());

	readVector(parameters.initial_pos_body, settings->InitialConditions().Position().Body());
	if (settings->InitialConditions().Position().UnsprungMass().present()) {
		parameters.initial_leg_flag = 1;
		readVectorLegs(parameters.initial_pos_wheel, settings->InitialConditions().Position().UnsprungMass().get());
		readVectorLegs(parameters.initial_pos_tyre, settings->InitialConditions().Position().Tyre().get());
	}

	readangles(parameters.initial_angle, settings->InitialConditions().Orientation());


//--------------------------------------------------
// Load simulation parameters
//--------------------------------------------------
   readVector(parameters.gravity, settings->SimulationParameters().GeneralSettings().Gravity());
    parameters.num_time_iter = settings->SimulationParameters().GeneralSettings().NumberOfIterations();
    parameters.timestep = settings->SimulationParameters().GeneralSettings().TimestepSize();
    

    parameters.DOF = settings->SimulationParameters().LinearALE().DOF();

    parameters.max_num_iter = settings->SimulationParameters().MultyBodyDynamics().MaximalIterationNumber();
    parameters.tolerance = settings->SimulationParameters().MultyBodyDynamics().Tolerance();
    parameters.solution_dim = settings->SimulationParameters().MultyBodyDynamics().SolutionDimension();

    std::string solver = settings->SimulationParameters().MultyBodyDynamics().Solver();


    if (solver=="explicit_Euler"){
        parameters.solver = EXPLICIT_EULER;
    } else if (solver=="RK4"){
        parameters.solver = RUNGE_KUTTA_4;
    } else if (solver=="Broyden_Euler"){
        parameters.solver = BROYDEN_EULER;
    } else if (solver=="Broyden_CN"){
        parameters.solver = BROYDEN_CN;
    } else if (solver=="Broyden_BDF2"){
        parameters.solver = BROYDEN_BDF2;
    } else {
        std::cerr<<"Wrong MBD-solver specified in the XML! Please type: \n   -explicit_Euler \n   -RK4\n   -Broyden_Euler\n   -Broyden_CN\n   -Broyden_BDF2"<<std::endl;
        exit(2);
    }
}

/**
* \brief read boundary condition and forces on tyres, wheel and body
*/
void ReadXML::ReadLoadParameters(Load_Params& parameters /**< [out] reference to parameter*/) {

    //--------------------------------------------------
    // Load external parameters
    //--------------------------------------------------
	std::string boundary_conditions = load_data->boundary_description().BoundaryConditions();
	if (boundary_conditions == "fixed") {
		parameters.boundary_condition_road = FIXED;
        std::cout << "Run the simulation with fixed tyres" << std::endl;
	}
	else if (boundary_conditions == "nonfixed") {
		parameters.boundary_condition_road = NONFIXED;
        std::cout << "Run the simulation without any tyre constraints" << std::endl;
    }
	else if (boundary_conditions == "circle") {
		parameters.boundary_condition_road = CIRCULAR;
        std::cout << "Run the simulation on a circular road" << std::endl;
    }
	else {
		std::cerr << "Wrong boundary conditions! Only circle, fixed and nonfixed implemented so far" << std::endl;
		exit(2);
	}
	if (load_data->boundary_description().circular().present()) {
		parameters.profile_radius = load_data->boundary_description().circular()->radius();
		readVector(parameters.profile_center, load_data->boundary_description().circular()->center());
	}
	readVectorLegs(parameters.external_force_tyre, load_data->forces().force_tyre());
	readVectorLegs(parameters.external_force_wheel, load_data->forces().force_wheel());
    readVector(parameters.external_force_body, load_data->forces().force_body());
}



/**
* \brief read and order params used for the loookup tables and generate the stiffness and damping loookup table
*/
void ReadXML::ReadLookupParameters(
    EVAALookup** lookupStiffness /**< [out] pointer to the pointer to the stiffness lookup from the compute engine*/,
    EVAALookup** lookupDamping /**< [out] pointer to the pointer to the damping lookup from the compute engine*/,
    Simulation_Parameters & parameters /**< [in] reference to the parameters and to set set the interpolation tag*/
) {
    if (_lookup_filename == "NO_FILE_SPECIFIED") return;
    lookup_table = LookupHandler(_lookup_filename, xml_schema::flags::dont_validate);
    if (lookup_table->LookupTableGenerator().present()) {
		double* a, * k_body, * k_tyre;
		double b, c, l_min, l_max;
		int size, k, type, order;

		a = new(double[8]);
		k_body = new(double[8]);
		k_tyre = new(double[8]);

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
		a[0] = k_body[2];
		a[1] = k_tyre[0];
		a[2] = k_body[3];
		a[3] = k_tyre[1];
		a[4] = k_body[0];
		a[5] = k_tyre[2];
		a[6] = k_body[1];
		a[7] = k_tyre[3];

        *lookupStiffness = new EVAALookup(size, a, b, c, l_min, l_max, k, type, order);
		parameters.interpolation = 1;  // to switch from constant to interpolation type

        // damping is /100 from the stiffness for the start
        for (auto j = 0; j < k; j++) {
            a[j] /= 100;
        }
        *lookupDamping = new EVAALookup(size, a, b, c, l_min, l_max, k, type, order);

        delete[] a;
		delete[] k_body;
		delete[] k_tyre;
    }
}