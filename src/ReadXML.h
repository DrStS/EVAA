#pragma once
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include "IP_EVAA_XML.h"
#include "LOAD_EVAA_XML.h"
#include "LOOKUP_EVAA_XML.h"
#ifndef U_COMPSTIFF
#define U_COMPSTIFF
#include "EVAAComputeStiffness.h"
#endif


struct Simulation_Parameters {int DOF; 
            double k_tyre[4]; double k_body[4]; 
            double c_tyre[4]; double c_body[4];
            double l_long[4]; double l_lat[4];
            double mass_body; double I_body[6];
            double mass_tyre[4]; double mass_wheel[4];
            double lower_spring_length[4]; double upper_spring_length[4];
            double initial_lower_spring_length[4]; double initial_upper_spring_length[4];
            double initial_vel_body[3]; double initial_vel_wheel[12]; double initial_vel_tyre[12];
            double initial_ang_vel_body[3]; double gravity[3];
			double initial_pos_body[3]; double initial_angle[4];
			double initial_pos_wheel[12]; double initial_pos_tyre[12]; // this has to be removed or used only if it is prescribed
			bool initial_leg_flag = 0;
			bool interpolation = 0;
            int solver; 
            int max_num_iter; double tolerance;
            double timestep; int num_time_iter;
			int solution_dim;
            };


struct Load_Params {
	double external_force_body[3]; double external_force_wheel[12]; double external_force_tyre[12];
	// circular profile params
	double profile_radius; double profile_center[3];
	// boundary condition type
	int boundary_condition_road;

};

enum boundary_condition_road{FIXED, NONFIXED, CIRCULAR};

enum solver {EXPLICIT_EULER, RUNGE_KUTTA_4, BROYDEN_EULER, BROYDEN_CN, BROYDEN_BDF2}; 



class ReadXML{
    private:
        std::string _filename;
		std::string _load_filename;
        std::string _lookup_filename;
        std::auto_ptr<car_settings_t> settings;
		std::auto_ptr<load_t> load_data;
        std::auto_ptr<lookup_handler_t> lookup_table;
        template<typename T> void readVectorLegs(double* storage, T vec);
        template<typename T> void readLegs(double* storage, T vec);
        template<typename T> void readVector(double* storage, T vec);
        template<typename T> void readangles(double* storage, T vec);

    public:
        ReadXML();
        ReadXML(const std::string & load_filename);
		ReadXML(const std::string & filename, const std::string & load_filename);
        void setFileName (const std::string & filename);
		void setloadFileName(const std::string & filename);
        void ReadParameters(Simulation_Parameters& parameters);
        void ReadLoadParameters(Load_Params & parameters);
        void ReadLookupParameters(EVAAComputeStiffness** lookupStiffness, Simulation_Parameters & parameters);
};

