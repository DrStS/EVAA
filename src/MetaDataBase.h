#pragma once
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include "Constants.h"
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

        // singleton instance
        static MetaDataBase* _database;

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
        double mass_body; double I_body[9];
        double mass_tyre[4]; double mass_wheel[4];
        double lower_spring_length[4]; double upper_spring_length[4];
        double initial_lower_spring_length[4]; double initial_upper_spring_length[4];
        double initial_vel_body[3]; double initial_vel_wheel[12]; double initial_vel_tyre[12];
        double initial_ang_vel_body[3]; double gravity[3];
        double initial_pos_body[3]; double initial_angle[4];
        double initial_pos_wheel[12]; double initial_pos_tyre[12]; // this has to be removed or used only if it is prescribed
        bool initial_leg_flag = 0;
        bool interpolation = 0;
        int MBD_solver;
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
        void ReadLookupParameters(EVAALookup<Constants::floatEVAA>** lookupStiffness, EVAALookup<Constants::floatEVAA>** lookupDamping);

        /*
        A bunch of getter functions
        */
        int getDOF() {
            return DOF;
        }
        double getTyreStiffnessFrontLeft() {
            return k_tyre[0];
        }
        double getTyreStiffnessFrontRight() {
            return k_tyre[1];
        }
        double getTyreStiffnessRearLeft() {
            return k_tyre[2];
        }
        double getTyreStiffnessRearRight() {
            return k_tyre[3];
        }
        // vector in the format fl fr rl rr
        double* getTyreStiffnessVector() {
            return k_tyre;
        }


        double getBodyStiffnessFrontLeft() {
            return k_body[0];
        }
        double getBodyStiffnessFrontRight() {
            return k_body[1];
        }
        double getBodyStiffnessRearLeft() {
            return k_body[2];
        }
        double getBodyStiffnessRearRight() {
            return k_body[3];
        }
        // vector in the format fl fr rl rr
        double* getBodyStiffnessVector() {
            return k_body;
        }


        double getTyreDampingFrontLeft() {
            return c_tyre[0];
        }
        double getTyreDampingFrontRight() {
            return c_tyre[1];
        }
        double getTyreDampingRearLeft() {
            return c_tyre[2];
        }
        double getTyreDampingRearRight() {
            return c_tyre[3];
        }
        // vector in the format fl fr rl rr
        double* getTyreDampingVector() {
            return c_tyre;
        }


        double getBodyDampingFrontLeft() {
            return c_body[0];
        }
        double getBodyDampingFrontRight() {
            return c_body[1];
        }
        double getBodyDampingRearLeft() {
            return c_body[2];
        }
        double getBodyDampingRearRight() {
            return c_body[3];
        }
        // vector in the format fl fr rl rr
        double* getBodyDampingVector() {
            return c_body;
        }


        double getLongitudalLegPositionFrontLeft() {
            return l_long[0];
        }
        double getLongitudalLegPositionFrontRight() {
            return l_long[1];
        }
        double getLongitudalLegPositionRearLeft() {
            return l_long[2];
        }
        double getLongitudalLegPositionRearRight() {
            return l_long[3];
        }
        // vector in the format fl fr rl rr
        double* getLongitudalLegPositionVector() {
            return l_long;
        }


        double getLatidudalLegPositionFrontLeft() {
            return l_lat[0];
        }
        double getLatidudalLegPositionFrontRight() {
            return l_lat[1];
        }
        double getLatidudalLegPositionRearLeft() {
            return l_lat[2];
        }
        double getLatidudalLegPositionRearRight() {
            return l_lat[3];
        }
        // vector in the format fl fr rl rr
        double* getLatidudalLegPositionVector() {
            return l_lat;
        }


        double getTyreMassFrontLeft() {
            return mass_tyre[0];
        }
        double getTyreMassFrontRight() {
            return mass_tyre[1];
        }
        double getTyreMassRearLeft() {
            return mass_tyre[2];
        }
        double getTyreMassRearRight() {
            return mass_tyre[3];
        }
        // vector in the format fl fr rl rr
        double* getTyreMassVector() {
            return mass_tyre;
        }


        double getWheelMassFrontLeft() {
            return mass_wheel[0];
        }
        double getWheelMassFrontRight() {
            return mass_wheel[1];
        }
        double getWheelMassRearLeft() {
            return mass_wheel[2];
        }
        double getWheelMassRearRight() {
            return mass_wheel[3];
        }
        // vector in the format fl fr rl rr
        double* getWheelMassVector() {
            return mass_wheel;
        }


        double getBodyMass() {
            return mass_body;
        }


        double getMomentOfInertiaXX() {
            return I_body[0];
        }
        double getMomentOfInertiaXY() {
            return I_body[1];
        }
        double getMomentOfInertiaXZ() {
            return I_body[2];
        }
        double getMomentOfInertiaYX() {
            return I_body[3];
        }
        double getMomentOfInertiaYY() {
            return I_body[4];
        }
        double getMomentOfInertiaYZ() {
            return I_body[5];
        }
        double getMomentOfInertiaZX() {
            return I_body[6];
        }
        double getMomentOfInertiaZY() {
            return I_body[7];
        }
        double getMomentOfInertiaZZ() {
            return I_body[8];
        }
        // contains all elements of the InertiaMatrix in the format [XX,XY,XZ,YX,YY,YZ,ZX,ZY,ZX] TODO: change it to Eigen/MKL output
        double* getMomentOfInertiaVector() {
            return I_body;
        }


        double getTyreSpringLengthFrontLeft() {
            return lower_spring_length[0];
        }
        double getTyreSpringLengthFrontRight() {
            return lower_spring_length[1];
        }
        double getTyreSpringLengthRearLeft() {
            return lower_spring_length[2];
        }
        double getTyreSpringLengthRearRight() {
            return lower_spring_length[3];
        }
        // vector in the format fl fr rl rr
        double* getTyreSpringLengthVector() {
            return lower_spring_length;
        }


        double getBodySpringLengthFrontLeft() {
            return upper_spring_length[0];
        }
        double getBodySpringLengthFrontRight() {
            return upper_spring_length[1];
        }
        double getBodySpringLengthRearLeft() {
            return upper_spring_length[2];
        }
        double getBodySpringLengthRearRight() {
            return upper_spring_length[3];
        }
        // vector in the format fl fr rl rr
        double* getBodySpringLengthVector() {
            return upper_spring_length;
        }


        double getTyreSpringInitialLengthFrontLeft() {
            return initial_lower_spring_length[0];
        }
        double getTyreSpringInitialLengthFrontRight() {
            return initial_lower_spring_length[1];
        }
        double getTyreSpringInitialLengthRearLeft() {
            return initial_lower_spring_length[2];
        }
        double getTyreSpringInitialLengthRearRight() {
            return initial_lower_spring_length[3];
        }
        // vector in the format fl fr rl rr
        double* getTyreSpringInitialLengthVector() {
            return initial_lower_spring_length;
        }


        double getBodySpringInitialLengthFrontLeft() {
            return initial_upper_spring_length[0];
        }
        double getBodySpringInitialLengthFrontRight() {
            return initial_upper_spring_length[1];
        }
        double getBodySpringInitialLengthRearLeft() {
            return initial_upper_spring_length[2];
        }
        double getBodySpringInitialLengthRearRight() {
            return initial_upper_spring_length[3];
        }
        // vector in the format fl fr rl rr
        double* getBodySpringInitialLengthVector() {
            return initial_upper_spring_length;
        }


        double* getBodyInitialVelocity() {
            return initial_vel_body;
        }


        double* getWheelInitialVelocityFrontLeft() {
            return initial_vel_wheel;
        }
        double* getWheelInitialVelocityFrontRight() {
            return initial_vel_wheel + 3;
        }
        double* getWheelInitialVelocityRearLeft() {
            return initial_vel_wheel + 6;
        }
        double* getWheelInitialVelocityRearRight() {
            return initial_vel_wheel + 9;
        }



        double* getTyreInitialVelocityFrontLeft() {
            return initial_vel_tyre;
        }
        double* getTyreInitialVelocityFrontRight() {
            return initial_vel_tyre + 3;
        }
        double* getTyreInitialVelocityRearLeft() {
            return initial_vel_tyre + 6;
        }
        double* getTyreInitialVelocityRearRight() {
            return initial_vel_tyre + 9;
        }


        double* getBodyInitialAngularVelocity() {
            return initial_ang_vel_body;
        }


        double* getGravityField() {
            return gravity;
        }

        double* getBodyInitialPosition() {
            return initial_pos_body;
        }


        double* getWheelInitialPositionFrontLeft() {
            return initial_pos_wheel;
        }
        double* getWheelInitialPositionFrontRight() {
            return initial_pos_wheel + 3;
        }
        double* getWheelInitialPositionRearLeft() {
            return initial_pos_wheel + 6;
        }
        double* getWheelInitialPositionRearRight() {
            return initial_pos_wheel + 9;
        }



        double* getTyreInitialPositionFrontLeft() {
            return initial_pos_tyre;
        }
        double* getTyreInitialPositionFrontRight() {
            return initial_pos_tyre + 3;
        }
        double* getTyreInitialPositionRearLeft() {
            return initial_pos_tyre + 6;
        }
        double* getTyreInitialPositionRearRight() {
            return initial_pos_tyre + 9;
        }


        double* getBodyInitialOrientation() {
            return initial_angle;
        }
        bool getFlagInitialLeg() {
            return initial_leg_flag;
        }

        bool getUseInterpolation() {
            return interpolation;
        }

        int getUsedSolverForMBD() {
            return MBD_solver;
        }

        int getMaxNumberOfBroydenIterationForMBD() {
            return max_num_iter;
        }

        int getToleranceBroydenIterationForMBD() {
            return tolerance;
        }

        double getTimeStepSize() {
            return timestep;
        }

        double getNumberOfTimeIterations() {
            return num_time_iter;
        }

        double getSolutionVectorSize() {
            return solution_dim;
        }


        double* getBodyExternalForce() {
            return external_force_body;
        }


        double* getWheelExternalForceFrontLeft() {
            return external_force_wheel;
        }
        double* getWheelExternalForceFrontRight() {
            return external_force_wheel + 3;
        }
        double* getWheelExternalForceRearLeft() {
            return external_force_wheel + 6;
        }
        double* getWheelExternalForceRearRight() {
            return external_force_wheel + 9;
        }



        double* getTyreExternalForceFrontLeft() {
            return external_force_tyre;
        }
        double* getTyreExternalForceFrontRight() {
            return external_force_tyre + 3;
        }
        double* getTyreExternalForceRearLeft() {
            return external_force_tyre + 6;
        }
        double* getTyreExternalForceRearRight() {
            return external_force_tyre + 9;
        }


        double getCircularRoadRadius() {
            return profile_radius;
        }

        double* getCircularRoadCenter() {
            return profile_center;
        }


        int getRoadConditions() {
            return boundary_condition_road;
        }
};

