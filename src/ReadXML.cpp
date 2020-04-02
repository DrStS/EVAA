#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>

#include "ReadXML.h"


void ReadXML::readVectorLegs(double* storage, legs_vector_t vec){
       readVector(storage, vec.ReerRight());
       readVector(storage+3,vec.ReerLeft());
       readVector(storage+6,vec.FrontLeft());
       readVector(storage+9,vec.FrontRight());
 }
void ReadXML::readVector(double* storage, vector_t vec) {
    storage[0] = vec.x();
    storage[1] = vec.y();
    storage[2] = vec.z();
}
void ReadXML::readLegs(double* storage, legs_t vec) {
    storage[0] = vec.ReerRight();
    storage[1] = vec.ReerLeft();
    storage[2] = vec.FrontLeft();
    storage[3] = vec.FrontRight();
}

void ReadXML::readangles(double* storage, quad vec) {
	storage[0] = vec.x();
	storage[1] = vec.y();
	storage[2] = vec.z();
	storage[3] = vec.w();
}


ReadXML::ReadXML(){
    _filename = "";
}

ReadXML::ReadXML(const std::string & filename):
            _filename(filename), 
            settings(EVAA_settings(filename, xml_schema::flags::dont_validate)){}

void ReadXML::setFileName(const std::string & filename){
    _filename = filename;
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
    readLegs(parameters.k_tyre, settings->Vehicle().TwoTrackModel().Stiffness().Tyre());
    readLegs(parameters.k_body, settings->Vehicle().TwoTrackModel().Stiffness().Body());
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

/*	readVector(parameters.initial_pos_body, settings->initial().pos_body());
	readVectorLegs(parameters.initial_pos_wheel, settings->initial().pos_wheel());
	readVectorLegs(parameters.initial_pos_tyre, settings->initial().pos_tyre());
*/

	readangles(parameters.initial_angle, settings->InitialConditions().Orientation());

//--------------------------------------------------
// Load external parameters
//--------------------------------------------------
/*    readVector(parameters.external_force_body, settings->external().force_body());
    readVector(parameters.gravity, settings->external().gravity());
    readVectorLegs(parameters.external_force_tyre, settings->external().force_tyre());
    readVectorLegs(parameters.external_force_wheel, settings->external().force_wheel());

    */

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
	std::string boundary_conditions = load_data->boundary_description.BoundaryConditions();
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
		std::cerr << "Wrong boundary conditions! Only fixed and nonfixed implemented so far" << std::endl;
		exit(2);
	}
	if (!load_data->boundary_description().circular().empty()) {
		for (auto iter = load_data->boundary_description().circular().begin(); iter < load_data->boundary_description().circular().end(); ++iter) {
			parameters.profile_radius = iter->radius();
			readVector(parameters.profile_center, iter->center());
		}
	}
	readVectorLegs(parameters.external_force_tyre, load_data->forces().force_tyre());
	readVectorLegs(parameters.external_force_wheel, load_data->forces().force_wheel());
    readVector(parameters.external_force_body, load_data->forces().force_body());
}
