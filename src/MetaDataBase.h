#pragma once
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include "IP_EVAA_XML.h"
#include "LOAD_EVAA_XML.h"
#include "LOOKUP_EVAA_XML.h"
#ifndef U_Lookup
#define U_Lookup
#include "EVAALookup.h"
#endif


enum boundary_condition_road{FIXED, NONFIXED, CIRCULAR};

/*
Solver for the MBD system
*/
enum solver {EXPLICIT_EULER, RUNGE_KUTTA_4, BROYDEN_EULER, BROYDEN_CN, BROYDEN_BDF2}; 


/*
Handles the XML parsers
*/
class MetaDataBase{
    private:
        MetaDataBase();
        MetaDataBase(const std::string& load_filename);
        MetaDataBase(const std::string& filename, const std::string& load_filename);

        // singleton instance
        static MetaDataBase* DataBase;

        // filenames
        std::string _filename;
		std::string _load_filename;
        std::string _lookup_filename;

        // pointer to the file parsers
        std::auto_ptr<car_settings_t> settings;
		std::auto_ptr<load_t> load_data;
        std::auto_ptr<lookup_handler_t> lookup_table;

        /*
        Read 4 legs which contain each 3 vectors
        \param vec XML parser
        \return storage all components in one vector with 12 elements [rr:XYZ,rl:XYZ,fl:XYZ,rl:XYZ]
        */
        template<typename T> void readVectorLegs(double* storage, T vec);

        /*
        Read 4 legs which contain each 1 double
        \param vec XML parser
        \return storage all components in one vector with 4 elements [rr,rl,fl,rl]
        */
        template<typename T> void readLegs(double* storage, T vec);

        /*
        Read a vector with 3 doubles
        \param vec XML parser
        \return storage all components in one vector with 3 elements [XYZ]
        */
        template<typename T> void readVector(double* storage, T vec);


        /*
        Read a quaternion with 4 doubles
        \param vec XML parser
        \return storage all components in one vector with 4 elements [XYZW]
        */
        template<typename T> void readangles(double* storage, T vec);

        /*
        Holds all general simulation parameters and car specific parameters (such as geometry and initial conditions)
        */
        int DOF;
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

        /*
        Environment parameters (road conditions and external force fields)
        */
        double external_force_body[3]; double external_force_wheel[12]; double external_force_tyre[12];
        // circular profile params
        double profile_radius; double profile_center[3];
        // boundary condition type
        int boundary_condition_road;


    public:
        /*
        \return static pointer to Database
        */
        static MetaDataBase* DataBase();

        /*
        Initializes the parser for global and car parameters
        \param filename
        */
        void setFileName (const std::string & filename); 

        /*
        Initializes the parser for environment parameters
        \param filename
        */
        void setloadFileName(const std::string & filename); 

        /*
        Read global and car parameters from the XML
        */
        void ReadParameters(); 

        /*
        Read environment parameters
        */
        void ReadLoadParameters();

        /*
        Read lookup table parameters
        \param parameters for global simulation parameters
        \return lookupStiffness class instance of the lookup handler
        */
        void ReadLookupParameters(EVAALookup** lookupStiffness, EVAALookup** lookupDamping);
};

