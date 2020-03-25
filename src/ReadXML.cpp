#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>

#include "ReadXML.h"


void ReadXML::readVectorLegs(double* storage, legs_vector_t vec){
       readVector(storage, vec.rr());
       readVector(storage+3,vec.rl());
       readVector(storage+6,vec.fl());
       readVector(storage+9,vec.fr());
 }
void ReadXML::readVector(double* storage, vector_t vec) {
    storage[0] = vec.z();
    storage[1] = vec.y();
    storage[2] = vec.z();
}
void ReadXML::readLegs(double* storage, legs_t vec) {
    storage[0] = vec.rr();
    storage[1] = vec.rl();
    storage[2] = vec.fl();
    storage[3] = vec.fr();
}

void ReadXML::readangles(double* storage, quad vec) {
	storage[0] = vec.w();
	storage[1] = vec.x();
	storage[2] = vec.y();
	storage[3] = vec.z();
}


ReadXML::ReadXML(){
    _filename = "";
}

ReadXML::ReadXML(const std::string & filename):
            _filename(filename), 
            settings(car_settings(filename, xml_schema::flags::dont_validate)){}

void ReadXML::setFileName(const std::string & filename){
    _filename = filename;
}

void ReadXML::ReadParameters(Simulation_Parameters & parameters){

    //--------------------------------------------------
    // Load car parameters
    //--------------------------------------------------
 
    parameters.mass_body = settings->car().mass_body();
    parameters.I_body[0] = settings->car().I_body_xx();
    parameters.I_body[1] = settings->car().I_body_yy();
    parameters.I_body[2] = settings->car().I_body_zz();
    parameters.I_body[3] = settings->car().I_body_xy();
    parameters.I_body[4] = settings->car().I_body_xz();
    parameters.I_body[5] = settings->car().I_body_yz();
    readLegs(parameters.k_tyre, settings->car().k_tyre());
    readLegs(parameters.k_body, settings->car().k_body());
    readLegs(parameters.c_tyre, settings->car().c_tyre());
    readLegs(parameters.c_body, settings->car().c_body());
    readLegs(parameters.l_long, settings->car().l_long());
    readLegs(parameters.l_lat, settings->car().l_lat());
    readLegs(parameters.mass_wheel, settings->car().mass_wheel());
    readLegs(parameters.mass_tyre, settings->car().mass_tyre());
    readLegs(parameters.lower_spring_length, settings->car().lower_spring_length());
    readLegs(parameters.upper_spring_length, settings->car().lower_spring_length());


//--------------------------------------------------
// Load initial parameters
//--------------------------------------------------

    readVector(parameters.initial_vel_body, settings->initial().vel_body());
    readVector(parameters.initial_ang_vel_body, settings->initial().ang_vel_body());
    
    readLegs(parameters.initial_lower_spring_length, settings->initial().lower_spring_length());
    readLegs(parameters.initial_upper_spring_length, settings->initial().upper_spring_length());
    readVectorLegs(parameters.initial_vel_wheel, settings->initial().vel_wheel());
    readVectorLegs(parameters.initial_vel_tyre, settings->initial().vel_tyre());

	readVector(parameters.initial_pos_body, settings->initial().pos_body());
	readVectorLegs(parameters.initial_pos_wheel, settings->initial().pos_wheel());
	readVectorLegs(parameters.initial_pos_tyre, settings->initial().pos_tyre());

	readangles(parameters.initial_angle, settings->initial().angle_body());

//--------------------------------------------------
// Load external parameters
//--------------------------------------------------
    readVector(parameters.external_force_body, settings->external().force_body());
    readVector(parameters.gravity, settings->external().gravity());
    readVectorLegs(parameters.external_force_tyre, settings->external().force_tyre());
    readVectorLegs(parameters.external_force_wheel, settings->external().force_wheel());
    std::string boundary_conditions=settings->external().boundary_conditions();

    if (boundary_conditions=="fixed"){
        parameters.boundary_condition_road = FIXED;
    } else if (boundary_conditions=="notfixed"){
        parameters.boundary_condition_road = NONFIXED;
    } else {
        std::cerr<<"Wrong boundary conditions! Only fixed and nonfixed implemented so far"<<std::endl;
        exit(2);
    }


//--------------------------------------------------
// Load simulation parameters
//--------------------------------------------------
    
    parameters.DOF = settings->simulation().DOF();
    parameters.max_num_iter = settings->simulation().max_num_iter();
    parameters.num_time_iter = settings->simulation().num_time_iter();
    parameters.timestep = settings->simulation().timestep();
    parameters.tolerance = settings->simulation().tolerance();
    
    std::string solver= settings->simulation().solver();
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
        std::cerr<<"Wrong solver! Only fixed and nonfixed implemented so far"<<std::endl;
        exit(2);
    }
}


void ReadXML::ReadVariableParameters(Simulation_Parameters& parameters) {

    //--------------------------------------------------
    // Load external parameters
    //--------------------------------------------------
    readVector(parameters.external_force_body, settings->external().force_body());
    readVector(parameters.gravity, settings->external().gravity());
    readVectorLegs(parameters.external_force_tyre, settings->external().force_tyre());
    readVectorLegs(parameters.external_force_wheel, settings->external().force_wheel());
    std::string boundary_conditions = settings->external().boundary_conditions();

    if (boundary_conditions == "fixed") {
        parameters.boundary_condition_road = FIXED;
    }
    else if (boundary_conditions == "notfixed") {
        parameters.boundary_condition_road = NONFIXED;
    }
    else {
        std::cerr << "Wrong boundary conditions! Only fixed and nonfixed implemented so far" << std::endl;
        exit(2);
    }


}
