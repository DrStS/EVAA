#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>

#include "ReadXML.h"


void ReadXML::readVectorLegs(double* storage, vector_legs_t vec){
       storage[0]=vec.rr();
       storage[1]=vec.rl();
       storage[2]=vec.fl();
       storage[3]=vec.fr();
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
    parameters.I_body[1] = settings->car().I_body_zz();
    parameters.I_body[2] = settings->car().I_body_xz();
    parameters.I_body[3] = settings->car().I_body_zx();
     readVectorLegs(parameters.k_tyre, settings->car().k_tyre());
     readVectorLegs(parameters.k_body, settings->car().k_body());
     readVectorLegs(parameters.l_long, settings->car().l_long());
     readVectorLegs(parameters.l_lat, settings->car().l_lat());
     readVectorLegs(parameters.mass_wheel, settings->car().mass_wheel());
     readVectorLegs(parameters.mass_tyre, settings->car().mass_tyre());
     readVectorLegs(parameters.lower_spring_length, settings->car().lower_spring_length());
     readVectorLegs(parameters.upper_spring_length, settings->car().lower_spring_length());


//--------------------------------------------------
// Load initial parameters
//--------------------------------------------------
    parameters.initial_orientation[0] = settings->initial().qi();
    parameters.initial_orientation[1] = settings->initial().qj();
    parameters.initial_orientation[2] = settings->initial().qk();
    parameters.initial_orientation[3] = settings->initial().qr();
    parameters.initial_vel_body = settings->initial().vel_body();
    parameters.initial_ang_vel_body[0] = settings->initial().ang_vel_body_x();
    parameters.initial_ang_vel_body[1] = settings->initial().ang_vel_body_z();
    
     readVectorLegs(parameters.initial_lower_spring_length, settings->initial().lower_spring_length());
     readVectorLegs(parameters.initial_upper_spring_length, settings->initial().upper_spring_length());
     readVectorLegs(parameters.initial_vel_wheel, settings->initial().vel_wheel());
     readVectorLegs(parameters.initial_vel_tyre, settings->initial().vel_tyre());

//--------------------------------------------------
// Load external parameters
//--------------------------------------------------
    parameters.external_force_body = settings->external().force_body();
    parameters.gravity = settings->external().gravity();
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
