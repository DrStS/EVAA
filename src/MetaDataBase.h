// TODO: copyright header

#pragma once

// TODO: remove cstdio
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>

#include "arbitraryTrajectory.h"

#include "Constants.h"

// TODO: U_Lookup feels artificial, find a way to not use it.
#ifndef U_Lookup
#define U_Lookup
#include "EVAALookup.h"
#endif

namespace EVAA {

/** Road condition. */
enum class BoundaryConditionRoad { FIXED, NONFIXED, CIRCULAR, ARBITRARY };

/** Solver for the MBD system. */
enum class MBDSolver { EXPLICIT_EULER, RUNGE_KUTTA_4, BROYDEN_EULER, BROYDEN_CN, BROYDEN_BDF2 };

/** Singleton handling input data parsing from XML files. */
class MetaDataBase {
public:
    /**
     * \return The singleton instance.
     */
    static MetaDataBase& DataBase();

    /** Deleted copy constructor. */
    MetaDataBase(MetaDataBase const&) = delete;

    /** Deleted copy operator. */
    void operator=(MetaDataBase const&) = delete;

    /**
     * Reads the car, initial and simulation parameters from an XML file.
     * \param[in] filename The XML file.
     */
    void readParameters(const std::string& filename);

    /**
     * Reads the load parameters from an XML file.
     * \param[in] loadFilename The XML file.
     */
    void readLoadParameters(const std::string& loadFilename);

    /* Getters. TODO: make all "double *get<X>" return const ref, mark the getter as const. */
    double getTyreStiffnessFrontLeft() const { return _k_tyre[Constants::FRONT_LEFT]; }
    double getTyreStiffnessFrontRight() const { return _k_tyre[Constants::FRONT_RIGHT]; }
    double getTyreStiffnessRearLeft() const { return _k_tyre[Constants::REAR_LEFT]; }
    double getTyreStiffnessRearRight() const { return _k_tyre[Constants::REAR_RIGHT]; }
    // vector in the format fl fr rl rr
    double* getTyreStiffnessVector() { return _k_tyre; }

    double getBodyStiffnessFrontLeft() const { return _k_body[Constants::FRONT_LEFT]; }
    double getBodyStiffnessFrontRight() const { return _k_body[Constants::FRONT_RIGHT]; }
    double getBodyStiffnessRearLeft() const { return _k_body[Constants::REAR_LEFT]; }
    double getBodyStiffnessRearRight() const { return _k_body[Constants::REAR_RIGHT]; }
    // vector in the format fl fr rl rr
    double* getBodyStiffnessVector() { return _k_body; }

    double getTyreDampingFrontLeft() const { return _c_tyre[Constants::FRONT_LEFT]; }
    double getTyreDampingFrontRight() const { return _c_tyre[Constants::FRONT_RIGHT]; }
    double getTyreDampingRearLeft() const { return _c_tyre[Constants::REAR_LEFT]; }
    double getTyreDampingRearRight() const { return _c_tyre[Constants::REAR_RIGHT]; }
    // vector in the format fl fr rl rr
    double* getTyreDampingVector() { return _c_tyre; }

    double getBodyDampingFrontLeft() const { return _c_body[Constants::FRONT_LEFT]; }
    double getBodyDampingFrontRight() const { return _c_body[Constants::FRONT_RIGHT]; }
    double getBodyDampingRearLeft() const { return _c_body[Constants::REAR_LEFT]; }
    double getBodyDampingRearRight() const { return _c_body[Constants::REAR_RIGHT]; }
    // vector in the format fl fr rl rr
    double* getBodyDampingVector() { return _c_body; }

    double getLongitudalLegPositionFrontLeft() { return _l_long[Constants::FRONT_LEFT]; }
    double getLongitudalLegPositionFrontRight() { return _l_long[Constants::FRONT_RIGHT]; }
    double getLongitudalLegPositionRearLeft() { return _l_long[Constants::REAR_LEFT]; }
    double getLongitudalLegPositionRearRight() { return _l_long[Constants::REAR_RIGHT]; }
    // vector in the format fl fr rl rr
    double* getLongitudalLegPositionVector() { return _l_long; }

    double getLatidudalLegPositionFrontLeft() { return _l_lat[Constants::FRONT_LEFT]; }
    double getLatidudalLegPositionFrontRight() { return _l_lat[Constants::FRONT_RIGHT]; }
    double getLatidudalLegPositionRearLeft() { return _l_lat[Constants::REAR_LEFT]; }
    double getLatidudalLegPositionRearRight() { return _l_lat[Constants::REAR_RIGHT]; }
    // vector in the format fl fr rl rr
    double* getLatidudalLegPositionVector() { return _l_lat; }

    double getTyreMassFrontLeft() { return _mass_tyre[Constants::FRONT_LEFT]; }
    double getTyreMassFrontRight() { return _mass_tyre[Constants::FRONT_RIGHT]; }
    double getTyreMassRearLeft() { return _mass_tyre[Constants::REAR_LEFT]; }
    double getTyreMassRearRight() { return _mass_tyre[Constants::REAR_RIGHT]; }
    // vector in the format fl fr rl rr
    double* getTyreMassVector() { return _mass_tyre; }

    double getWheelMassFrontLeft() { return _mass_wheel[Constants::FRONT_LEFT]; }
    double getWheelMassFrontRight() { return _mass_wheel[Constants::FRONT_RIGHT]; }
    double getWheelMassRearLeft() { return _mass_wheel[Constants::REAR_LEFT]; }
    double getWheelMassRearRight() { return _mass_wheel[Constants::REAR_RIGHT]; }
    // vector in the format fl fr rl rr
    double* getWheelMassVector() { return _mass_wheel; }

    double getBodyMass() { return _mass_body; }

    double getMomentOfInertiaXX() { return _I_body[0]; }
    double getMomentOfInertiaXY() { return _I_body[1]; }
    double getMomentOfInertiaXZ() { return _I_body[2]; }
    double getMomentOfInertiaYX() { return _I_body[3]; }
    double getMomentOfInertiaYY() { return _I_body[4]; }
    double getMomentOfInertiaYZ() { return _I_body[5]; }
    double getMomentOfInertiaZX() { return _I_body[6]; }
    double getMomentOfInertiaZY() { return _I_body[7]; }
    double getMomentOfInertiaZZ() { return _I_body[8]; }
    // contains all elements of the InertiaMatrix in the format [XX,XY,XZ,YX,YY,YZ,ZX,ZY,ZZ]
    double* getMomentOfInertiaVector() { return _I_body; }

    double getTyreSpringLengthFrontLeft() { return _lower_spring_length[Constants::FRONT_LEFT]; }
    double getTyreSpringLengthFrontRight() { return _lower_spring_length[Constants::FRONT_RIGHT]; }
    double getTyreSpringLengthRearLeft() { return _lower_spring_length[Constants::REAR_LEFT]; }
    double getTyreSpringLengthRearRight() { return _lower_spring_length[Constants::REAR_RIGHT]; }
    // vector in the format fl fr rl rr
    double* getTyreSpringLengthVector() { return _lower_spring_length; }

    double getBodySpringLengthFrontLeft() { return _upper_spring_length[Constants::FRONT_LEFT]; }
    double getBodySpringLengthFrontRight() { return _upper_spring_length[Constants::FRONT_RIGHT]; }
    double getBodySpringLengthRearLeft() { return _upper_spring_length[Constants::REAR_LEFT]; }
    double getBodySpringLengthRearRight() { return _upper_spring_length[Constants::REAR_RIGHT]; }
    // vector in the format fl fr rl rr
    double* getBodySpringLengthVector() { return _upper_spring_length; }

    double getTyreSpringInitialLengthFrontLeft() {
        return _initial_lower_spring_length[Constants::FRONT_LEFT];
    }
    double getTyreSpringInitialLengthFrontRight() {
        return _initial_lower_spring_length[Constants::FRONT_RIGHT];
    }
    double getTyreSpringInitialLengthRearLeft() {
        return _initial_lower_spring_length[Constants::REAR_LEFT];
    }
    double getTyreSpringInitialLengthRearRight() {
        return _initial_lower_spring_length[Constants::REAR_RIGHT];
    }
    // vector in the format fl fr rl rr
    double* getTyreSpringInitialLengthVector() { return _initial_lower_spring_length; }

    double getBodySpringInitialLengthFrontLeft() {
        return _initial_upper_spring_length[Constants::FRONT_LEFT];
    }
    double getBodySpringInitialLengthFrontRight() {
        return _initial_upper_spring_length[Constants::FRONT_RIGHT];
    }
    double getBodySpringInitialLengthRearLeft() {
        return _initial_upper_spring_length[Constants::REAR_LEFT];
    }
    double getBodySpringInitialLengthRearRight() {
        return _initial_upper_spring_length[Constants::REAR_RIGHT];
    }
    // vector in the format fl fr rl rr
    double* getBodySpringInitialLengthVector() { return _initial_upper_spring_length; }

    double* getBodyInitialVelocity() { return _initial_vel_body; }

    double* getWheelInitialVelocityFrontLeft() {
        return _initial_vel_wheel + Constants::DIM * Constants::FRONT_LEFT;
    }
    double* getWheelInitialVelocityFrontRight() {
        return _initial_vel_wheel + Constants::DIM * Constants::FRONT_RIGHT;
    }
    double* getWheelInitialVelocityRearLeft() {
        return _initial_vel_wheel + Constants::DIM * Constants::REAR_LEFT;
    }
    double* getWheelInitialVelocityRearRight() {
        return _initial_vel_wheel + Constants::DIM * Constants::REAR_RIGHT;
    }

    double* getTyreInitialVelocityFrontLeft() {
        return _initial_vel_tyre + Constants::DIM * Constants::FRONT_LEFT;
    }
    double* getTyreInitialVelocityFrontRight() {
        return _initial_vel_tyre + Constants::DIM * Constants::FRONT_RIGHT;
    }
    double* getTyreInitialVelocityRearLeft() {
        return _initial_vel_tyre + Constants::DIM * Constants::REAR_LEFT;
    }
    double* getTyreInitialVelocityRearRight() {
        return _initial_vel_tyre + Constants::DIM * Constants::REAR_RIGHT;
    }

    double* getBodyInitialAngularVelocity() { return _initial_ang_vel_body; }

    double* getGravityField() { return _gravity; }

    double* getBodyInitialPosition() { return _initial_pos_body; }

    double* getWheelInitialPositionFrontLeft() {
        return _initial_pos_wheel + Constants::DIM * Constants::FRONT_LEFT;
    }
    double* getWheelInitialPositionFrontRight() {
        return _initial_pos_wheel + Constants::DIM * Constants::FRONT_RIGHT;
    }
    double* getWheelInitialPositionRearLeft() {
        return _initial_pos_wheel + Constants::DIM * Constants::REAR_LEFT;
    }
    double* getWheelInitialPositionRearRight() {
        return _initial_pos_wheel + Constants::DIM * Constants::REAR_RIGHT;
    }

    double* getTyreInitialPositionFrontLeft() {
        return _initial_pos_tyre + Constants::DIM * Constants::FRONT_LEFT;
    }
    double* getTyreInitialPositionFrontRight() {
        return _initial_pos_tyre + Constants::DIM * Constants::FRONT_RIGHT;
    }
    double* getTyreInitialPositionRearLeft() {
        return _initial_pos_tyre + Constants::DIM * Constants::REAR_LEFT;
    }
    double* getTyreInitialPositionRearRight() {
        return _initial_pos_tyre + Constants::DIM * Constants::REAR_RIGHT;
    }

    double* getBodyInitialOrientation() { return _initial_angle; }
    bool getFlagInitialLeg() { return _initial_leg_flag; }

    bool getUseInterpolation() { return _interpolation; }

    MBDSolver getUsedSolverForMBD() { return _MBD_solver; }

    int getMaxNumberOfBroydenIterationForMBD() { return _max_num_iter; }

    int getToleranceBroydenIterationForMBD() { return _tolerance; }

    double getTimeStepSize() { return _timestep; }

    double getNumberOfTimeIterations() { return _num_time_iter; }

    double getSolutionVectorSize() { return _solution_dim; }

    double* getBodyExternalForce() { return _external_force_body; }

    double* getWheelExternalForceFrontLeft() {
        return _external_force_wheel + Constants::DIM * Constants::FRONT_LEFT;
    }
    double* getWheelExternalForceFrontRight() {
        return _external_force_wheel + Constants::DIM * Constants::FRONT_RIGHT;
    }
    double* getWheelExternalForceRearLeft() {
        return _external_force_wheel + Constants::DIM * Constants::REAR_LEFT;
    }
    double* getWheelExternalForceRearRight() {
        return _external_force_wheel + Constants::DIM * Constants::REAR_RIGHT;
    }

    double* getTyreExternalForceFrontLeft() {
        return _external_force_tyre + Constants::DIM * Constants::FRONT_LEFT;
    }
    double* getTyreExternalForceFrontRight() {
        return _external_force_tyre + Constants::DIM * Constants::FRONT_RIGHT;
    }
    double* getTyreExternalForceRearLeft() {
        return _external_force_tyre + Constants::DIM * Constants::REAR_LEFT;
    }
    double* getTyreExternalForceRearRight() {
        return _external_force_tyre + Constants::DIM * Constants::REAR_RIGHT;
    }
    arbitraryTrajectory<double>* getArbitraryTrajectory() { return _trajectory; }

    double getCircularRoadRadius() { return _profile_radius; }
    double* getCircularRoadCenter() { return _profile_center; }

    BoundaryConditionRoad getRoadConditions() { return _boundary_condition_road; }

    const EVAALookup<Constants::floatEVAA>& getLookupStiffness() const { return *_lookupStiffness; }

    const EVAALookup<Constants::floatEVAA>& getLookupDamping() const { return *_lookupDamping; }

private:
    /** Private constructor for the singleton instance. */
    MetaDataBase();

    /**
     * Reads the lookup table parameters from an XML file.
     * \param[in] filename The XML file.
     */
    void readLookupParameters(const std::string& filename);

    /**
     * \brief Reads 4 legs which contain each EVAA::Constants::DIM vectors.
     * Reads the vectors of the positions of the legs relative to the center of mass.
     * \param vec XML parser
     * \return storage all components in one vector with 12 elements [rr:XYZ,rl:XYZ,fl:XYZ,rl:XYZ]
     */
    template <typename T>
    void readVectorLegs(double* storage, T vec);

    /**
     * \brief Read 4 legs which contain each 1 double.
     * Reads paramaters which are given for all the legs.
     * \param vec XML parser
     * \return storage all components in one vector with 4 elements [rr,rl,fl,rl]
     */
    template <typename T>
    void readLegs(double* storage, T vec);

    /**
     * \brief Reads a vector with 3 doubles.
     * \param vec XML parser
     * \return storage all components in one vector with 3 elements [XYZ]
     */
    template <typename T>
    void readVector(double* storage, T vec);

    /**
     * \brief Reads a quaternion with 4 doubles.
     * \param vec XML parser
     * \return storage all components in one vector with 4 elements [XYZW]
     */
    template <typename T>
    void readAngles(double* storage, T vec);

    /*
    All general simulation parameters and car specific parameters (such as geometry and
    initial conditions)
    */
    double _k_tyre[Constants::NUM_LEGS];
    double _k_body[Constants::NUM_LEGS];
    double _c_tyre[Constants::NUM_LEGS];
    double _c_body[Constants::NUM_LEGS];
    double _l_long[Constants::NUM_LEGS];
    double _l_lat[Constants::NUM_LEGS];
    double _mass_body;
    double _I_body[9];
    double _mass_tyre[Constants::NUM_LEGS];
    double _mass_wheel[Constants::NUM_LEGS];
    double _lower_spring_length[Constants::NUM_LEGS];
    double _upper_spring_length[Constants::NUM_LEGS];
    double _initial_lower_spring_length[Constants::NUM_LEGS];
    double _initial_upper_spring_length[Constants::NUM_LEGS];
    double _initial_vel_body[Constants::DIM];
    double _initial_vel_wheel[Constants::DIM * Constants::NUM_LEGS];
    double _initial_vel_tyre[Constants::DIM * Constants::NUM_LEGS];
    double _initial_ang_vel_body[Constants::DIM];
    double _gravity[Constants::DIM];
    double _initial_pos_body[Constants::DIM];
    double _initial_angle[Constants::NUM_LEGS];
    double _initial_pos_wheel[Constants::DIM * Constants::NUM_LEGS];
    double _initial_pos_tyre[Constants::DIM * Constants::NUM_LEGS];  // this has to be removed or used
                                                                  // only if it is prescribed
    bool _initial_leg_flag = false;
    bool _interpolation = false;
    MBDSolver _MBD_solver;
    int _max_num_iter;
    double _tolerance;
    double _timestep;
    int _num_time_iter;
    int _solution_dim;

    arbitraryTrajectory<double>* _trajectory;

    /*
    Environment parameters (road conditions and external force fields)
    */
    double _external_force_body[Constants::DIM];
    double _external_force_wheel[Constants::DIM * Constants::NUM_LEGS];
    double _external_force_tyre[Constants::DIM * Constants::NUM_LEGS];

    /*
    Circular profile params
    */
    double _profile_radius;
    double _profile_center[Constants::DIM];
    BoundaryConditionRoad _boundary_condition_road;

    /** Lookup filename read from the car file. */
    std::string _lookup_filename;

    /** Stiffness lookup from the compute engine. */
    EVAALookup<Constants::floatEVAA>* _lookupStiffness;
    /** Damping lookup from the compute engine. */
    EVAALookup<Constants::floatEVAA>* _lookupDamping;
};

}  // namespace EVAA
