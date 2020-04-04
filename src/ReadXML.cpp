#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>

#include "ReadXML.h"


template<typename T> void ReadXML::readVectorLegs(double* storage, T vec){
       readVector(storage, vec.ReerRight());
       readVector(storage+3,vec.ReerLeft());
       readVector(storage+6,vec.FrontLeft());
       readVector(storage+9,vec.FrontRight());
 }
template<typename T> void ReadXML::readVector(double* storage, T vec) {
    storage[0] = vec.x();
    storage[1] = vec.y();
    storage[2] = vec.z();
}
template<typename T> void ReadXML::readLegs(double* storage, T vec) {
    storage[0] = vec.ReerRight();
    storage[1] = vec.ReerLeft();
    storage[2] = vec.FrontLeft();
    storage[3] = vec.FrontRight();
}

template<typename T> void ReadXML::readangles(double* storage, T vec) {
	storage[0] = vec.x();
	storage[1] = vec.y();
	storage[2] = vec.z();
	storage[3] = vec.w();
}


ReadXML::ReadXML(){
    _filename = "";
	_load_filename = "";
}

ReadXML::ReadXML(const std::string& load_filename) :
    _load_filename(load_filename),
    load_data(EVAA_load_module(load_filename, xml_schema::flags::dont_validate)) {
}

ReadXML::ReadXML(const std::string & filename, const std::string & load_filename) :
	_filename(filename), _load_filename(load_filename),
	settings(EVAA_settings(filename, xml_schema::flags::dont_validate)), load_data(EVAA_load_module(load_filename, xml_schema::flags::dont_validate)) {
}

void ReadXML::setFileName(const std::string & filename){
    _filename = filename;
}

void ReadXML::setloadFileName(const std::string & filename) {
	_load_filename = filename;
}

void ReadXML::ReadParameters(Simulation_Parameters & parameters){

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

void ReadXML::ReadLoadParameters(Load_Params& parameters) {

    //--------------------------------------------------
    // Load external parameters
    //--------------------------------------------------
	std::string boundary_conditions = load_data->boundary_description().BoundaryConditions();
	if (boundary_conditions == "fixed") {
		parameters.boundary_condition_road = FIXED;
	}
	else if (boundary_conditions == "notfixed") {
		parameters.boundary_condition_road = NONFIXED;
	}
	else if (boundary_conditions == "circle") {
		parameters.boundary_condition_road = CIRCULAR;
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

void ReadXML::ReadLookupParameters(EVAAComputeStiffness* lookupStiffness, Simulation_Parameters & parameters) {
    if (_lookup_filename == "NO_FILE_SPECIFIED") return;
    lookup_table = LookupHandler(_lookup_filename, xml_schema::flags::dont_validate);
    if (lookup_table->LookupTableGenerator().present()) {
        double* a;
        double b, c, l_min, l_max;
        int size, k, type, order;

        a = new(double[8]);

        std::cout << "Generate look up table from parameters." << std::endl;

        size = lookup_table->LookupTableGenerator().get().Size();
        b = lookup_table->LookupTableGenerator().get().TableParameters().b();
        c = lookup_table->LookupTableGenerator().get().TableParameters().c();
        l_min = lookup_table->LookupTableGenerator().get().Range().l_min();
        l_max = lookup_table->LookupTableGenerator().get().Range().l_max();
        k = lookup_table->LookupTableGenerator().get().InterpolationMethod().k();
        type = lookup_table->LookupTableGenerator().get().InterpolationMethod().type();
        order = lookup_table->LookupTableGenerator().get().InterpolationMethod().order();

        readLegs(a, lookup_table->LookupTableGenerator().get().Magnitude().Body());
        readLegs(a+4, lookup_table->LookupTableGenerator().get().Magnitude().Tyre());

        //@Felix, initialize your lookup table here, add and change function signatures, private variables as you need it
		// EVAAComputeStiffness(int size, double a, double b, double c, double l_min, double l_max, int k, int type, int order);
		lookupStiffness = new EVAAComputeStiffness(size, a, b, c, l_min, l_max, k, type, order);
		parameters.interpolation = 1;  // to switch from constant to interpolation type
        delete[] a;
    }
}