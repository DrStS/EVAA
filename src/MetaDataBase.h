// TODO: copyright header
#pragma once

// TODO: remove cstdio
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

#include "Constants.h"
#include "IO/Output.h"
#include "arbitraryTrajectory.h"

// TODO: U_Lookup feels artificial, find a way to not use it.
#ifndef U_Lookup
#define U_Lookup
#include "EVAALookup.h"
#endif

#include "IP_EVAA_XML.h"
#include "LOAD_EVAA_XML.h"
#include "LOOKUP_EVAA_XML.h"

namespace EVAA {

/** Road condition. */
enum class BoundaryConditionRoad { FIXED, NONFIXED, CIRCULAR, ARBITRARY };

/** Solver for the MBD system. */
enum class MBDSolver { EXPLICIT_EULER, RUNGE_KUTTA_4, BROYDEN_EULER, BROYDEN_CN, BROYDEN_BDF2 };

/** Singleton handling input data parsing from XML files. */
template <class T>
class MetaDataBase {
public:
    /**
     * \return The singleton instance.
     */
    static MetaDataBase& getDataBase() {
        static MetaDataBase database;
        return database;
    }

    /** Deleted copy constructor. */
    MetaDataBase(MetaDataBase const&) = delete;

    /** Deleted copy operator. */
    void operator=(MetaDataBase const&) = delete;

    /**
     * Reads the car, initial and simulation parameters from an XML file.
     * \param[in] filename The XML file.
     */
    void MetaDataBase::readParameters(const std::string& filename) {
        // Load car parameters

        const auto settings = EVAA_settings(filename, xml_schema::flags::dont_validate);
        const auto twoTrackModel = settings->VehicleXML().TwoTrackModelXML();

        _mass_body = twoTrackModel.MassXML().BodyXML();
        readLegs(_mass_wheel, twoTrackModel.MassXML().UnsprungMassXML());
        readLegs(_mass_tyre, twoTrackModel.MassXML().TyreXML());
        _I_body[0] = twoTrackModel.InertiaXML().XX();
        _I_body[1] = twoTrackModel.InertiaXML().XY();
        _I_body[2] = twoTrackModel.InertiaXML().XZ();
        _I_body[3] = twoTrackModel.InertiaXML().YX();
        _I_body[4] = twoTrackModel.InertiaXML().YY();
        _I_body[5] = twoTrackModel.InertiaXML().YZ();
        _I_body[6] = twoTrackModel.InertiaXML().ZX();
        _I_body[7] = twoTrackModel.InertiaXML().ZY();
        _I_body[8] = twoTrackModel.InertiaXML().ZZ();

        if (twoTrackModel.StiffnessXML().ConstantXML().present()) {
            std::cout << "Take constant stiffness without lookup table" << std::endl;
            readLegs(_k_tyre, twoTrackModel.StiffnessXML().ConstantXML().get().TyreXML());
            readLegs(_k_body, twoTrackModel.StiffnessXML().ConstantXML().get().BodyXML());
        }
        else {
            readLookupParameters(twoTrackModel.StiffnessXML().LookupTableXML().get().FilePathXML());
        }
        readLegs(_c_tyre, twoTrackModel.DampingCoefficientsXML().TyreXML());
        readLegs(_c_body, twoTrackModel.DampingCoefficientsXML().BodyXML());
        readLegs(_l_long, twoTrackModel.GeometryXML().LongitudinalReferenceToWheelXML());
        readLegs(_l_lat, twoTrackModel.GeometryXML().LateralReferenceToWheelXML());
        readVector(_vehicleCIR,
                   twoTrackModel.GeometryXML().RelativeCenterOfInstanteneousRotation());
        readLegs(_lower_spring_length, twoTrackModel.GeometryXML().SuspensionSpringsXML());
        readLegs(_upper_spring_length, twoTrackModel.GeometryXML().TyreSpringsXML());

        // Load initial parameters
        const auto initial = settings->InitialConditionsXML();

        readVector(_initial_vel_body, initial.VelocitiesXML().BodyXML());
        readVector(_initial_ang_vel_body, initial.VelocitiesXML().angularBodyXML());

        readVectorLegs(_initial_vel_wheel, initial.VelocitiesXML().UnsprungMassXML());
        readVectorLegs(_initial_vel_tyre, initial.VelocitiesXML().TyreXML());

        readLegs(_initial_lower_spring_length, initial.SpringElongationXML().TyreXML());
        readLegs(_initial_upper_spring_length, initial.SpringElongationXML().BodyXML());

        readVector(_initial_pos_body, initial.PositionXML().BodyXML());
        if (initial.PositionXML().UnsprungMassXML().present()) {
            _initial_leg_flag = 1;
            readVectorLegs(_initial_pos_wheel, initial.PositionXML().UnsprungMassXML().get());
            readVectorLegs(_initial_pos_tyre, initial.PositionXML().TyreXML().get());
        }

        readAngles(_initialAngleGlobal, initial.OrientationXML());

        // Load simulation parameters
        const auto simulation = settings->SimulationParametersXML();

        readVector(_gravity, simulation.GeneralSettingsXML().GravityXML());
        _num_time_iter = simulation.GeneralSettingsXML().NumberOfIterationsXML();
        _timestep = simulation.GeneralSettingsXML().TimestepSizeXML();

        _max_num_iter = simulation.MultyBodyDynamicsXML().MaximalIterationNumberXML();
        _tolerance = simulation.MultyBodyDynamicsXML().ToleranceXML();
        _solution_dim = simulation.MultyBodyDynamicsXML().SolutionDimensionXML();
        std::string solver = simulation.MultyBodyDynamicsXML().SolverXML();

        if (solver == "explicit_Euler") {
            _MBD_solver = MBDSolver::EXPLICIT_EULER;
        }
        else if (solver == "RK4") {
            _MBD_solver = MBDSolver::RUNGE_KUTTA_4;
        }
        else if (solver == "Broyden_Euler") {
            _MBD_solver = MBDSolver::BROYDEN_EULER;
        }
        else if (solver == "Broyden_CN") {
            _MBD_solver = MBDSolver::BROYDEN_CN;
        }
        else if (solver == "Broyden_BDF2") {
            _MBD_solver = MBDSolver::BROYDEN_BDF2;
        }
        else {
            throw std::logic_error("Wrong MBD-solver in XML: " + solver +
                                   ". Must be one of: "
                                   "explicit_Euler, RK4, Broyden_Euler, Broyden_CN, Broyden_BDF2");
        }
    }

    /**
     * Reads the load parameters from an XML file.
     * \param[in] loadFilename The XML file.
     */
    void readLoadParameters(const std::string& loadFilename) {
        //--------------------------------------------------
        // Load external parameters
        //--------------------------------------------------
        const auto load_data = EVAA_load_module(loadFilename, xml_schema::flags::dont_validate);
        const auto sinusoidalProfile = load_data->roadProfile()
                                           .arbitraryRoadProfile()
                                           .get()
                                           .verticalProfile()
                                           .sinusoidalProfile();
        const auto horizontalProfile =
            load_data->roadProfile().arbitraryRoadProfile().get().horizontalProfile();

        if (load_data->roadProfile().fixedTyre().present()) {
            _boundary_condition_road = BoundaryConditionRoad::FIXED;
            std::cout << "Run the simulation with fixed tyres" << std::endl;
        }
        else if (load_data->roadProfile().detachedTyre().present()) {
            _boundary_condition_road = BoundaryConditionRoad::NONFIXED;
            std::cout << "Run the simulation without any tyre constraints" << std::endl;
        }
        else if (load_data->roadProfile().circularRoadProfile().present()) {
            _boundary_condition_road = BoundaryConditionRoad::CIRCULAR;
            std::cout << "Run the simulation on a circular road" << std::endl;
            _profile_radius = load_data->roadProfile().circularRoadProfile()->radius();
            readVector(_profile_center, load_data->roadProfile().circularRoadProfile()->center());
        }
        else if (load_data->roadProfile().arbitraryRoadProfile().present()) {
            _boundary_condition_road = BoundaryConditionRoad::ARBITRARY;
            std::cout << "Run the simulation on an arbitrary road" << std::endl;
            _trajectory = new arbitraryTrajectory<T>(
                _num_time_iter, _timestep, sinusoidalProfile.rightTyre().amplitude(),
                sinusoidalProfile.leftTyre().amplitude(), sinusoidalProfile.rightTyre().period(),
                sinusoidalProfile.leftTyre().period(), sinusoidalProfile.rightTyre().shift(),
                sinusoidalProfile.leftTyre().shift());  // TODO: free memory without leaks

            T initialVelocity = sqrt(_initial_vel_body[0] * _initial_vel_body[0] +
                                     _initial_vel_body[1] * _initial_vel_body[1]);

            std::vector<T> wayPointsX;
            std::vector<T> wayPointsY;
            std::vector<T> wayPointsTimes;

            size_t numWayPoints = 0;

            for (auto i = horizontalProfile.wayPoint().begin();
                 i < horizontalProfile.wayPoint().end(); ++i) {
                numWayPoints++;
                wayPointsX.push_back(i->X());
                wayPointsY.push_back(i->Y());
                wayPointsTimes.push_back(i->time());
            }

            _trajectory->interpolateRoadPoints(numWayPoints, wayPointsX.data(), wayPointsY.data(),
                                               wayPointsTimes.data());
            _trajectory->calculateTyreShifts(
                _l_long[Constants::FRONT_LEFT], _l_long[Constants::FRONT_RIGHT],
                _l_long[Constants::REAR_LEFT], _l_long[Constants::REAR_RIGHT]);
            _trajectory->calculateTravelledDistance();
            _trajectory->calculateAngles();
            _trajectory->calculateAccelerationsCenterOfGravity();
            _trajectory->calculateAccelerationsLegs(
                _l_long[Constants::FRONT_LEFT], _l_long[Constants::FRONT_RIGHT],
                _l_long[Constants::REAR_LEFT], _l_long[Constants::REAR_RIGHT],
                _l_lat[Constants::FRONT_LEFT], _l_lat[Constants::FRONT_RIGHT],
                _l_lat[Constants::REAR_LEFT], _l_lat[Constants::REAR_RIGHT]);

            _trajectory->calculateVerticalAccelerations();
        }
        else {
            throw std::logic_error(
                "Wrong boundary conditions. Implemented so far: circle, fixed, nonfixed.");
        }

        readVectorLegs(_external_force_tyre, load_data->forces().forceTyre());
        readVectorLegs(_external_force_wheel, load_data->forces().forceWheel());
        readVector(_external_force_body, load_data->forces().forceBody());
    }

    /* Getters. TODO: make all "T *get<X>" return const ref, mark the getter as const. */
    T getTyreStiffnessFrontLeft() const { return _k_tyre[Constants::FRONT_LEFT]; }
    T getTyreStiffnessFrontRight() const { return _k_tyre[Constants::FRONT_RIGHT]; }
    T getTyreStiffnessRearLeft() const { return _k_tyre[Constants::REAR_LEFT]; }
    T getTyreStiffnessRearRight() const { return _k_tyre[Constants::REAR_RIGHT]; }

    T getBodyStiffnessFrontLeft() const { return _k_body[Constants::FRONT_LEFT]; }
    T getBodyStiffnessFrontRight() const { return _k_body[Constants::FRONT_RIGHT]; }
    T getBodyStiffnessRearLeft() const { return _k_body[Constants::REAR_LEFT]; }
    T getBodyStiffnessRearRight() const { return _k_body[Constants::REAR_RIGHT]; }

    T getTyreDampingFrontLeft() const { return _c_tyre[Constants::FRONT_LEFT]; }
    T getTyreDampingFrontRight() const { return _c_tyre[Constants::FRONT_RIGHT]; }
    T getTyreDampingRearLeft() const { return _c_tyre[Constants::REAR_LEFT]; }
    T getTyreDampingRearRight() const { return _c_tyre[Constants::REAR_RIGHT]; }

    T getBodyDampingFrontLeft() const { return _c_body[Constants::FRONT_LEFT]; }
    T getBodyDampingFrontRight() const { return _c_body[Constants::FRONT_RIGHT]; }
    T getBodyDampingRearLeft() const { return _c_body[Constants::REAR_LEFT]; }
    T getBodyDampingRearRight() const { return _c_body[Constants::REAR_RIGHT]; }

    T getLongitudalLegPositionFrontLeft() { return _l_long[Constants::FRONT_LEFT]; }
    T getLongitudalLegPositionFrontRight() { return _l_long[Constants::FRONT_RIGHT]; }
    T getLongitudalLegPositionRearLeft() { return _l_long[Constants::REAR_LEFT]; }
    T getLongitudalLegPositionRearRight() { return _l_long[Constants::REAR_RIGHT]; }
    // vector in the format fl fr rl rr
    T* getLongitudalLegPositionVector() { return _l_long; }

    T getLatidudalLegPositionFrontLeft() { return _l_lat[Constants::FRONT_LEFT]; }
    T getLatidudalLegPositionFrontRight() { return _l_lat[Constants::FRONT_RIGHT]; }
    T getLatidudalLegPositionRearLeft() { return _l_lat[Constants::REAR_LEFT]; }
    T getLatidudalLegPositionRearRight() { return _l_lat[Constants::REAR_RIGHT]; }
    // vector in the format fl fr rl rr
    T* getLatidudalLegPositionVector() { return _l_lat; }

    T* getPositionCenterOfInstantaneousRotation() { return _vehicleCIR; }

    T getTyreMassFrontLeft() { return _mass_tyre[Constants::FRONT_LEFT]; }
    T getTyreMassFrontRight() { return _mass_tyre[Constants::FRONT_RIGHT]; }
    T getTyreMassRearLeft() { return _mass_tyre[Constants::REAR_LEFT]; }
    T getTyreMassRearRight() { return _mass_tyre[Constants::REAR_RIGHT]; }
    // vector in the format fl fr rl rr
    T* getTyreMassVector() { return _mass_tyre; }

    T getWheelMassFrontLeft() { return _mass_wheel[Constants::FRONT_LEFT]; }
    T getWheelMassFrontRight() { return _mass_wheel[Constants::FRONT_RIGHT]; }
    T getWheelMassRearLeft() { return _mass_wheel[Constants::REAR_LEFT]; }
    T getWheelMassRearRight() { return _mass_wheel[Constants::REAR_RIGHT]; }
    // vector in the format fl fr rl rr
    T* getWheelMassVector() { return _mass_wheel; }

    T getBodyMass() { return _mass_body; }

    T getMomentOfInertiaXX() { return _I_body[0]; }
    T getMomentOfInertiaXY() { return _I_body[1]; }
    T getMomentOfInertiaXZ() { return _I_body[2]; }
    T getMomentOfInertiaYX() { return _I_body[3]; }
    T getMomentOfInertiaYY() { return _I_body[4]; }
    T getMomentOfInertiaYZ() { return _I_body[5]; }
    T getMomentOfInertiaZX() { return _I_body[6]; }
    T getMomentOfInertiaZY() { return _I_body[7]; }
    T getMomentOfInertiaZZ() { return _I_body[8]; }
    // contains all elements of the InertiaMatrix in the format [XX,XY,XZ,YX,YY,YZ,ZX,ZY,ZZ]
    T* getMomentOfInertiaVector() { return _I_body; }

    T getTyreSpringLengthFrontLeft() { return _lower_spring_length[Constants::FRONT_LEFT]; }
    T getTyreSpringLengthFrontRight() { return _lower_spring_length[Constants::FRONT_RIGHT]; }
    T getTyreSpringLengthRearLeft() { return _lower_spring_length[Constants::REAR_LEFT]; }
    T getTyreSpringLengthRearRight() { return _lower_spring_length[Constants::REAR_RIGHT]; }
    // vector in the format fl fr rl rr
    T* getTyreSpringLengthVector() { return _lower_spring_length; }

    T getBodySpringLengthFrontLeft() { return _upper_spring_length[Constants::FRONT_LEFT]; }
    T getBodySpringLengthFrontRight() { return _upper_spring_length[Constants::FRONT_RIGHT]; }
    T getBodySpringLengthRearLeft() { return _upper_spring_length[Constants::REAR_LEFT]; }
    T getBodySpringLengthRearRight() { return _upper_spring_length[Constants::REAR_RIGHT]; }
    // vector in the format fl fr rl rr
    T* getBodySpringLengthVector() { return _upper_spring_length; }

    T getTyreSpringInitialLengthFrontLeft() {
        return _initial_lower_spring_length[Constants::FRONT_LEFT];
    }
    T getTyreSpringInitialLengthFrontRight() {
        return _initial_lower_spring_length[Constants::FRONT_RIGHT];
    }
    T getTyreSpringInitialLengthRearLeft() {
        return _initial_lower_spring_length[Constants::REAR_LEFT];
    }
    T getTyreSpringInitialLengthRearRight() {
        return _initial_lower_spring_length[Constants::REAR_RIGHT];
    }
    // vector in the format fl fr rl rr
    T* getTyreSpringInitialLengthVector() { return _initial_lower_spring_length; }

    T getBodySpringInitialLengthFrontLeft() {
        return _initial_upper_spring_length[Constants::FRONT_LEFT];
    }
    T getBodySpringInitialLengthFrontRight() {
        return _initial_upper_spring_length[Constants::FRONT_RIGHT];
    }
    T getBodySpringInitialLengthRearLeft() {
        return _initial_upper_spring_length[Constants::REAR_LEFT];
    }
    T getBodySpringInitialLengthRearRight() {
        return _initial_upper_spring_length[Constants::REAR_RIGHT];
    }
    // vector in the format fl fr rl rr
    T* getBodySpringInitialLengthVector() { return _initial_upper_spring_length; }

    T* getBodyInitialVelocity() { return _initial_vel_body; }

    T* getWheelInitialVelocityFrontLeft() {
        return _initial_vel_wheel + Constants::DIM * Constants::FRONT_LEFT;
    }
    T* getWheelInitialVelocityFrontRight() {
        return _initial_vel_wheel + Constants::DIM * Constants::FRONT_RIGHT;
    }
    T* getWheelInitialVelocityRearLeft() {
        return _initial_vel_wheel + Constants::DIM * Constants::REAR_LEFT;
    }
    T* getWheelInitialVelocityRearRight() {
        return _initial_vel_wheel + Constants::DIM * Constants::REAR_RIGHT;
    }

    T* getTyreInitialVelocityFrontLeft() {
        return _initial_vel_tyre + Constants::DIM * Constants::FRONT_LEFT;
    }
    T* getTyreInitialVelocityFrontRight() {
        return _initial_vel_tyre + Constants::DIM * Constants::FRONT_RIGHT;
    }
    T* getTyreInitialVelocityRearLeft() {
        return _initial_vel_tyre + Constants::DIM * Constants::REAR_LEFT;
    }
    T* getTyreInitialVelocityRearRight() {
        return _initial_vel_tyre + Constants::DIM * Constants::REAR_RIGHT;
    }

    T* getBodyInitialAngularVelocity() { return _initial_ang_vel_body; }

    T* getGravityField() { return _gravity; }

    T* getBodyInitialPosition() { return _initial_pos_body; }

    T* getWheelInitialPositionFrontLeft() {
        return _initial_pos_wheel + Constants::DIM * Constants::FRONT_LEFT;
    }
    T* getWheelInitialPositionFrontRight() {
        return _initial_pos_wheel + Constants::DIM * Constants::FRONT_RIGHT;
    }
    T* getWheelInitialPositionRearLeft() {
        return _initial_pos_wheel + Constants::DIM * Constants::REAR_LEFT;
    }
    T* getWheelInitialPositionRearRight() {
        return _initial_pos_wheel + Constants::DIM * Constants::REAR_RIGHT;
    }

    T* getTyreInitialPositionFrontLeft() {
        return _initial_pos_tyre + Constants::DIM * Constants::FRONT_LEFT;
    }
    T* getTyreInitialPositionFrontRight() {
        return _initial_pos_tyre + Constants::DIM * Constants::FRONT_RIGHT;
    }
    T* getTyreInitialPositionRearLeft() {
        return _initial_pos_tyre + Constants::DIM * Constants::REAR_LEFT;
    }
    T* getTyreInitialPositionRearRight() {
        return _initial_pos_tyre + Constants::DIM * Constants::REAR_RIGHT;
    }

    T* getBodyInitialOrientation() { return _initialAngleGlobal; }
    bool getFlagInitialLeg() { return _initial_leg_flag; }

    bool getUseInterpolation() { return _interpolation; }

    MBDSolver getUsedSolverForMBD() { return _MBD_solver; }

    int getMaxNumberOfBroydenIterationForMBD() { return _max_num_iter; }

    int getToleranceBroydenIterationForMBD() { return _tolerance; }

    T getTimeStepSize() { return _timestep; }

    T getNumberOfTimeIterations() { return _num_time_iter; }

    T getSolutionVectorSize() { return _solution_dim; }

    T* getBodyExternalForce() { return _external_force_body; }

    T* getWheelExternalForceFrontLeft() {
        return _external_force_wheel + Constants::DIM * Constants::FRONT_LEFT;
    }
    T* getWheelExternalForceFrontRight() {
        return _external_force_wheel + Constants::DIM * Constants::FRONT_RIGHT;
    }
    T* getWheelExternalForceRearLeft() {
        return _external_force_wheel + Constants::DIM * Constants::REAR_LEFT;
    }
    T* getWheelExternalForceRearRight() {
        return _external_force_wheel + Constants::DIM * Constants::REAR_RIGHT;
    }

    T* getTyreExternalForceFrontLeft() {
        return _external_force_tyre + Constants::DIM * Constants::FRONT_LEFT;
    }
    T* getTyreExternalForceFrontRight() {
        return _external_force_tyre + Constants::DIM * Constants::FRONT_RIGHT;
    }
    T* getTyreExternalForceRearLeft() {
        return _external_force_tyre + Constants::DIM * Constants::REAR_LEFT;
    }
    T* getTyreExternalForceRearRight() {
        return _external_force_tyre + Constants::DIM * Constants::REAR_RIGHT;
    }
    arbitraryTrajectory<T>* getArbitraryTrajectory() { return _trajectory; }

    T getCircularRoadRadius() { return _profile_radius; }
    T* getCircularRoadCenter() { return _profile_center; }

    BoundaryConditionRoad getRoadConditions() { return _boundary_condition_road; }

    const EVAALookup<Constants::floatEVAA>& getLookupStiffness() const { return *_lookupStiffness; }

    const EVAALookup<Constants::floatEVAA>& getLookupDamping() const { return *_lookupDamping; }

private:
    /** Private constructor for the singleton instance. */
    MetaDataBase() : _lookup_filename("") {}

    /**
     * Reads the lookup table parameters from an XML file.
     * \param[in] filename The XML file.
     */
    void readLookupParameters(const std::string& filename) {
        IO::checkFileExists(filename);
        const auto lookupHandler = LookupHandler(filename, xml_schema::flags::dont_validate);
        if (lookupHandler->LookupTableGenerator().present()) {
            const auto lookupTable = lookupHandler->LookupTableGenerator().get();
            T *a, *_k_body, *_k_tyre;
            T b, c, l_min, l_max;
            int size, k, type, order;

            a = new (T[8]);
            _k_body = new (T[8]);
            _k_tyre = new (T[8]);

            std::cout << "Generate look up table from parameters." << std::endl;

            size = lookupTable.Size();
            b = lookupTable.TableParameters().b();
            c = lookupTable.TableParameters().c();
            l_min = lookupTable.Range().l_min();
            l_max = lookupTable.Range().l_max();
            k = lookupTable.InterpolationMethod().k();
            type = lookupTable.InterpolationMethod().type();
            order = lookupTable.InterpolationMethod().order();

            readLegs(_k_body, lookupTable.Magnitude().Body());
            readLegs(_k_tyre, lookupTable.Magnitude().Tyre());
            a[0] = _k_body[2];
            a[1] = _k_tyre[0];
            a[2] = _k_body[3];
            a[3] = _k_tyre[1];
            a[4] = _k_body[0];
            a[5] = _k_tyre[2];
            a[6] = _k_body[1];
            a[7] = _k_tyre[3];

            // TODO: release memory.
            _lookupStiffness =
                new EVAALookup<Constants::floatEVAA>(size, a, b, c, l_min, l_max, k, type, order);
            _interpolation = 1;  // to switch from constant to interpolation type

            // damping is /100 from the stiffness for the start
            for (auto j = 0; j < k; j++) {
                // TODO: extract magic number 100?
                a[j] /= 100;
            }

            // TODO: release memory.
            _lookupDamping =
                new EVAALookup<Constants::floatEVAA>(size, a, b, c, l_min, l_max, k, type, order);

            delete[] a;
            delete[] _k_body;
            delete[] _k_tyre;
        }
    }

    /**
     * \brief Reads 4 legs which contain each EVAA::Constants::DIM vectors.
     * Reads the vectors of the positions of the legs relative to the center of mass.
     * \param vec XML parser
     * \return storage all components in one vector with 12 elements [rr:XYZ,rl:XYZ,fl:XYZ,rl:XYZ]
     */
    template <typename S>
    void readVectorLegs(T* storage, S vec) {
        readVector(storage + Constants::DIM * Constants::FRONT_LEFT, vec.FrontLeft());
        readVector(storage + Constants::DIM * Constants::FRONT_RIGHT, vec.FrontRight());
        readVector(storage + Constants::DIM * Constants::REAR_LEFT, vec.ReerLeft());
        readVector(storage + Constants::DIM * Constants::REAR_RIGHT, vec.ReerRight());
    }

    /**
     * \brief Read 4 legs which contain each 1 T.
     * Reads paramaters which are given for all the legs.
     * \param vec XML parser
     * \return storage all components in one vector with 4 elements [rr,rl,fl,rl]
     */
    template <typename S>
    void readLegs(T* storage, S vec) {
        storage[Constants::FRONT_LEFT] = vec.FrontLeft();
        storage[Constants::FRONT_RIGHT] = vec.FrontRight();
        storage[Constants::REAR_LEFT] = vec.ReerLeft();
        storage[Constants::REAR_RIGHT] = vec.ReerRight();
    }

    /**
     * \brief Reads a vector with 3 Ts.
     * \param vec XML parser
     * \return storage all components in one vector with 3 elements [XYZ]
     */
    template <typename S>
    void readVector(T* storage, S vec) {
        storage[0] = vec.x();
        storage[1] = vec.y();
        storage[2] = vec.z();
    }

    /**
     * \brief Reads a quaternion with 4 Ts.
     * \param vec XML parser
     * \return storage all components in one vector with 4 elements [XYZW]
     */
    template <typename S>
    void readAngles(T* storage, S vec) {
        storage[0] = vec.x();
        storage[1] = vec.y();
        storage[2] = vec.z();
        storage[3] = vec.w();
    }

    /*
    All general simulation parameters and car specific parameters (such as geometry and
    initial conditions)
    */
    T _k_tyre[Constants::NUM_LEGS];
    T _k_body[Constants::NUM_LEGS];
    T _c_tyre[Constants::NUM_LEGS];
    T _c_body[Constants::NUM_LEGS];
    T _l_long[Constants::NUM_LEGS];
    T _l_lat[Constants::NUM_LEGS];
    T _vehicleCIR[Constants::DIM];
    T _mass_body;
    T _I_body[9];
    T _mass_tyre[Constants::NUM_LEGS];
    T _mass_wheel[Constants::NUM_LEGS];
    T _lower_spring_length[Constants::NUM_LEGS];
    T _upper_spring_length[Constants::NUM_LEGS];
    T _initial_lower_spring_length[Constants::NUM_LEGS];
    T _initial_upper_spring_length[Constants::NUM_LEGS];
    T _initial_vel_body[Constants::DIM];
    T _initial_vel_wheel[Constants::DIM * Constants::NUM_LEGS];
    T _initial_vel_tyre[Constants::DIM * Constants::NUM_LEGS];
    T _initial_ang_vel_body[Constants::DIM];
    T _gravity[Constants::DIM];
    T _initial_pos_body[Constants::DIM];
    T _initialAngleGlobal[Constants::DIM + 1];
    T _initial_pos_wheel[Constants::DIM * Constants::NUM_LEGS];
    T _initial_pos_tyre[Constants::DIM * Constants::NUM_LEGS];  // this has to be removed or used
                                                                // only if it is prescribed
    bool _initial_leg_flag = false;
    bool _interpolation = false;
    MBDSolver _MBD_solver;
    int _max_num_iter;
    T _tolerance;
    T _timestep;
    int _num_time_iter;
    int _solution_dim;

    arbitraryTrajectory<T>* _trajectory;

    /*
    Environment parameters (road conditions and external force fields)
    */
    T _external_force_body[Constants::DIM];
    T _external_force_wheel[Constants::DIM * Constants::NUM_LEGS];
    T _external_force_tyre[Constants::DIM * Constants::NUM_LEGS];

    /*
    Circular profile params
    */
    T _profile_radius;
    T _profile_center[Constants::DIM];
    BoundaryConditionRoad _boundary_condition_road;

    /** Lookup filename read from the car file. */
    std::string _lookup_filename;

    /** Stiffness lookup from the compute engine. */
    EVAALookup<Constants::floatEVAA>* _lookupStiffness;
    /** Damping lookup from the compute engine. */
    EVAALookup<Constants::floatEVAA>* _lookupDamping;
};

}  // namespace EVAA
