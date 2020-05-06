// TODO: copyright header
#pragma once

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

#include "ArbitraryTrajectory.h"
#include "Constants.h"
//#include "IO/Output.h"
#include "MathLibrary.h"

#ifdef INTERPOLATION
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

/** Solver for the ALE system */
enum class ALESolver {
    IMPLICIT_EULER, BDF2
};

/** Singleton handling input data parsing from XML files. */
template <class T>
class MetaDatabase {
public:
    /**
     * \return The singleton instance.
     */
    static MetaDatabase& getDatabase() {
        static MetaDatabase database;
        return database;
    }

    /** Deleted copy constructor. */
    MetaDatabase(MetaDatabase const&) = delete;

    /** Deleted copy operator. */
    void operator=(MetaDatabase const&) = delete;

    /*
    Checks if the xml has correct format with resect to compiler flag
    */
    void CheckStiffnessXML(const bool interpolation, const bool fixspring) {
#ifdef INTERPOLATION
        if (interpolation) {
            std::cout << "Using lookup tables for Stiffness and Damping" << std::endl;
        }
        else {
            throw "Lookup block are not available in the XML!";
        }
        if (fixspring) {
            throw "Cannot use constant block with interpolation. Check the "
                  "XML!";
        }
#else
        if (fixspring) {
            std::cout << "Using constant Stiffness and Damping" << std::endl;
        }
        else {
            throw "Constant Stiffness and damping blocks are not available in the XML!";
        }
        if (interpolation) {
            throw "Cannot use lookup block with constant compiler flag. Check the XML!";
        }
#endif
    }

    /**
     * Reads the car, initial and simulation parameters from an XML file.
     * \param[in] filename The XML file.
     */
    void MetaDatabase::readParameters(const std::string& filename) {
        // Load car parameters

        const auto settings = EVAA_settings(filename, xml_schema::flags::dont_validate);
        const auto twoTrackModel = settings->VehicleXML().TwoTrackModelXML();

        _mass_body = twoTrackModel.MassXML().BodyXML();
        readLegs(_mass, 2, twoTrackModel.MassXML().UnsprungMassXML());
        readLegs(_mass + 1, 2, twoTrackModel.MassXML().TyreXML());
        _I_body[0] = twoTrackModel.InertiaXML().XX();
        _I_body[1] = twoTrackModel.InertiaXML().XY();
        _I_body[2] = twoTrackModel.InertiaXML().XZ();
        _I_body[3] = twoTrackModel.InertiaXML().YX();
        _I_body[4] = twoTrackModel.InertiaXML().YY();
        _I_body[5] = twoTrackModel.InertiaXML().YZ();
        _I_body[6] = twoTrackModel.InertiaXML().ZX();
        _I_body[7] = twoTrackModel.InertiaXML().ZY();
        _I_body[8] = twoTrackModel.InertiaXML().ZZ();

        CheckStiffnessXML(twoTrackModel.StiffnessXML().LookupTableXML().present(), twoTrackModel.StiffnessXML().ConstantXML().present());

#ifdef INTERPOLATION
        readLookupParameters(twoTrackModel.StiffnessXML().LookupTableXML().get().FilePathXML());
#else
        readLegs(_k_tyre, twoTrackModel.StiffnessXML().ConstantXML().get().TyreXML());
        readLegs(_k_body, twoTrackModel.StiffnessXML().ConstantXML().get().BodyXML());
#endif
        readLegs(_c_tyre, twoTrackModel.DampingCoefficientsXML().TyreXML());
        readLegs(_c_body, twoTrackModel.DampingCoefficientsXML().BodyXML());
        readLegs(_l_long, twoTrackModel.GeometryXML().LongitudinalReferenceToWheelXML());
        readLegs(_l_lat, twoTrackModel.GeometryXML().LateralReferenceToWheelXML());
        readVector(_vehicleCIR, twoTrackModel.GeometryXML().RelativeCenterOfInstanteneousRotation());
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
        else if (solver == "BroydenEuler") {
            _MBD_solver = MBDSolver::BROYDEN_EULER;
        }
        else if (solver == "BroydenCN") {
            _MBD_solver = MBDSolver::BROYDEN_CN;
        }
        else if (solver == "BroydenBDF2") {
            _MBD_solver = MBDSolver::BROYDEN_BDF2;
        }
        else {
            throw std::logic_error("Wrong MBD-solver in XML: " + solver +
                                   ". Must be one of: "
                                   "explicit_Euler, RK4, BroydenEuler, BroydenCN, BroydenBDF2");
        }

        solver = simulation.LinearALEXML().Method();
        if (solver == "IMPLICIT_EULER") {
            _ALE_solver = ALESolver::IMPLICIT_EULER;
        }
        else if (solver == "BDF2") {
            _ALE_solver = ALESolver::BDF2;
        }
        else {
            throw std::logic_error("Wrong ALE-solver in XML: " + solver +
                                   ". Must be one of: "
                                   "IMPLICIT_EULER, BDF2");
        }
        _maxNewtonIterations = simulation.LinearALEXML().MaximumNumNewtonIterations();
        _newtonTolerance = simulation.LinearALEXML().Tolerance();
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
            const auto sinusoidalProfile = load_data->roadProfile().arbitraryRoadProfile().get().verticalProfile().sinusoidalProfile();
            const auto horizontalProfile = load_data->roadProfile().arbitraryRoadProfile().get().horizontalProfile();
            _boundary_condition_road = BoundaryConditionRoad::ARBITRARY;
            std::cout << "Run the simulation on an arbitrary road" << std::endl;
            _trajectory = new ArbitraryTrajectory<T>(
                _num_time_iter, _timestep, sinusoidalProfile.rightTyre().amplitude(),
                sinusoidalProfile.leftTyre().amplitude(), sinusoidalProfile.rightTyre().period(),
                sinusoidalProfile.leftTyre().period(), sinusoidalProfile.rightTyre().shift(),
                sinusoidalProfile.leftTyre().shift(), _initial_upper_spring_length,
                _initial_lower_spring_length, _l_lat, _l_long);  // TODO: free memory without leaks

            T initialVelocity = sqrt(_initial_vel_body[0] * _initial_vel_body[0] + _initial_vel_body[1] * _initial_vel_body[1]);

            std::vector<T> wayPointsX;
            std::vector<T> wayPointsY;
            std::vector<T> wayPointsTimes;

            size_t numWayPoints = 0;

            for (auto i = horizontalProfile.wayPoint().begin(); i < horizontalProfile.wayPoint().end(); ++i) {
                numWayPoints++;
                wayPointsX.push_back(i->X());
                wayPointsY.push_back(i->Y());
                wayPointsTimes.push_back(i->time());
            }

            _trajectory->calculateTyreShifts();
            _trajectory->calculateVerticalPositionsLegs();
            _trajectory->interpolateRoadPoints(numWayPoints, wayPointsX.data(), wayPointsY.data(),
                                               wayPointsTimes.data());
            _trajectory->calculateTravelledDistance();
            _trajectory->calculateAngles();
            _trajectory->calculateAccelerationsCenterOfGravity();
            _trajectory->calculateAccelerationsLegs();

            _trajectory->calculateVerticalAccelerations();
        }
        else {
            throw std::logic_error(
                "Wrong boundary conditions. Implemented so far: circle, fixed, "
                "nonfixed.");
        }

        readVectorLegs(_external_force_tyre, load_data->forces().forceTyre());
        readVectorLegs(_external_force_wheel, load_data->forces().forceWheel());
        readVector(_external_force_body, load_data->forces().forceBody());
    }

    /* Getters. TODO: make all "T *get<X>" return const ref, mark the getter as
     * const. */
    inline T getTyreStiffnessFrontLeft() const { return _k_tyre[Constants::FRONT_LEFT]; }
    inline T getTyreStiffnessFrontRight() const { return _k_tyre[Constants::FRONT_RIGHT]; }
    inline T getTyreStiffnessRearLeft() const { return _k_tyre[Constants::REAR_LEFT]; }
    inline T getTyreStiffnessRearRight() const { return _k_tyre[Constants::REAR_RIGHT]; }
    inline const T* getTyreStiffnessVector() const { return _k_tyre; }

    inline T getBodyStiffnessFrontLeft() const { return _k_body[Constants::FRONT_LEFT]; }
    inline T getBodyStiffnessFrontRight() const { return _k_body[Constants::FRONT_RIGHT]; }
    inline T getBodyStiffnessRearLeft() const { return _k_body[Constants::REAR_LEFT]; }
    inline T getBodyStiffnessRearRight() const { return _k_body[Constants::REAR_RIGHT]; }
    inline const T* getBodyStiffnessVector() const { return _k_body; }

    inline T getTyreDampingFrontLeft() const { return _c_tyre[Constants::FRONT_LEFT]; }
    inline T getTyreDampingFrontRight() const { return _c_tyre[Constants::FRONT_RIGHT]; }
    inline T getTyreDampingRearLeft() const { return _c_tyre[Constants::REAR_LEFT]; }
    inline T getTyreDampingRearRight() const { return _c_tyre[Constants::REAR_RIGHT]; }
    inline const T* getTyreDampingVector() const { return _c_tyre; }

    inline T getBodyDampingFrontLeft() const { return _c_body[Constants::FRONT_LEFT]; }
    inline T getBodyDampingFrontRight() const { return _c_body[Constants::FRONT_RIGHT]; }
    inline T getBodyDampingRearLeft() const { return _c_body[Constants::REAR_LEFT]; }
    inline T getBodyDampingRearRight() const { return _c_body[Constants::REAR_RIGHT]; }
    inline const T* getBodyDampingVector() const { return _c_body; }

    inline T getLongitudalLegPositionFrontLeft() const { return _l_long[Constants::FRONT_LEFT]; }
    inline T getLongitudalLegPositionFrontRight() const { return _l_long[Constants::FRONT_RIGHT]; }
    inline T getLongitudalLegPositionRearLeft() const { return _l_long[Constants::REAR_LEFT]; }
    inline T getLongitudalLegPositionRearRight() const { return _l_long[Constants::REAR_RIGHT]; }
    // vector in the format fl fr rl rr
    inline T* getLongitudalLegPositionVector() { return _l_long; }

    inline T getLatidudalLegPositionFrontLeft() const { return _l_lat[Constants::FRONT_LEFT]; }
    inline T getLatidudalLegPositionFrontRight() const { return _l_lat[Constants::FRONT_RIGHT]; }
    inline T getLatidudalLegPositionRearLeft() const { return _l_lat[Constants::REAR_LEFT]; }
    inline T getLatidudalLegPositionRearRight() const { return _l_lat[Constants::REAR_RIGHT]; }
    // vector in the format fl fr rl rr
    inline T* getLatidudalLegPositionVector() { return _l_lat; }

    inline T* getPositionCenterOfInstantaneousRotation() { return _vehicleCIR; }

    inline T getTyreMassFrontLeft() const { return _mass[2 * Constants::FRONT_LEFT + 1]; }
    inline T getTyreMassFrontRight() const { return _mass[2 * Constants::FRONT_RIGHT + 1]; }
    inline T getTyreMassRearLeft() const { return _mass[2 * Constants::REAR_LEFT + 1]; }
    inline T getTyreMassRearRight() const { return _mass[2 * Constants::REAR_RIGHT + 1]; }

    inline T getWheelMassFrontLeft() const { return _mass[2 * Constants::FRONT_LEFT]; }
    inline T getWheelMassFrontRight() const { return _mass[2 * Constants::FRONT_RIGHT]; }
    inline T getWheelMassRearLeft() const { return _mass[2 * Constants::REAR_LEFT]; }
    inline T getWheelMassRearRight() const { return _mass[2 * Constants::REAR_RIGHT]; }

    /**
     * Returns the leg mass vector.
     * \return A mass vector in format [W_fl, T_fl, W_fr, T_fr, W_rl, T_rl,
     * W_rr, T_rr]
     */
    inline const T* getWheelTyreMassVector() { return _mass; }

    inline T getBodyMass() const { return _mass_body; }

    inline T getMomentOfInertiaXX() const { return _I_body[0]; }
    inline T getMomentOfInertiaYY() const { return _I_body[4]; }
    // contains all elements of the InertiaMatrix in the format
    // [XX,XY,XZ,YX,YY,YZ,ZX,ZY,ZZ]
    inline T* getMomentOfInertiaVector() { return _I_body; }

    inline T getTyreSpringLengthFrontLeft() const { return _lower_spring_length[Constants::FRONT_LEFT]; }
    inline T getTyreSpringLengthFrontRight() const { return _lower_spring_length[Constants::FRONT_RIGHT]; }
    inline T getTyreSpringLengthRearLeft() const { return _lower_spring_length[Constants::REAR_LEFT]; }
    inline T getTyreSpringLengthRearRight() const { return _lower_spring_length[Constants::REAR_RIGHT]; }
    // vector in the format fl fr rl rr
    inline T* getTyreSpringLengthVector() { return _lower_spring_length; }

    inline T getBodySpringLengthFrontLeft() const { return _upper_spring_length[Constants::FRONT_LEFT]; }
    inline T getBodySpringLengthFrontRight() const { return _upper_spring_length[Constants::FRONT_RIGHT]; }
    inline T getBodySpringLengthRearLeft() const { return _upper_spring_length[Constants::REAR_LEFT]; }
    inline T getBodySpringLengthRearRight() const { return _upper_spring_length[Constants::REAR_RIGHT]; }
    // vector in the format fl fr rl rr
    inline T* getBodySpringLengthVector() { return _upper_spring_length; }

    inline T getTyreSpringInitialLengthFrontLeft() const { return _initial_lower_spring_length[Constants::FRONT_LEFT]; }
    inline T getTyreSpringInitialLengthFrontRight() const { return _initial_lower_spring_length[Constants::FRONT_RIGHT]; }
    inline T getTyreSpringInitialLengthRearLeft() const { return _initial_lower_spring_length[Constants::REAR_LEFT]; }
    inline T getTyreSpringInitialLengthRearRight() const { return _initial_lower_spring_length[Constants::REAR_RIGHT]; }
    // vector in the format fl fr rl rr
    inline T* getTyreSpringInitialLengthVector() { return _initial_lower_spring_length; }

    inline T getBodySpringInitialLengthFrontLeft() const { return _initial_upper_spring_length[Constants::FRONT_LEFT]; }
    inline T getBodySpringInitialLengthFrontRight() const { return _initial_upper_spring_length[Constants::FRONT_RIGHT]; }
    inline T getBodySpringInitialLengthRearLeft() const { return _initial_upper_spring_length[Constants::REAR_LEFT]; }
    inline T getBodySpringInitialLengthRearRight() const { return _initial_upper_spring_length[Constants::REAR_RIGHT]; }
    // vector in the format fl fr rl rr
    inline T* getBodySpringInitialLengthVector() { return _initial_upper_spring_length; }

    inline T* getBodyInitialVelocity() { return _initial_vel_body; }

    inline T* getWheelInitialVelocityFrontLeft() { return _initial_vel_wheel + Constants::DIM * Constants::FRONT_LEFT; }
    inline T* getWheelInitialVelocityFrontRight() { return _initial_vel_wheel + Constants::DIM * Constants::FRONT_RIGHT; }
    inline T* getWheelInitialVelocityRearLeft() { return _initial_vel_wheel + Constants::DIM * Constants::REAR_LEFT; }
    inline T* getWheelInitialVelocityRearRight() { return _initial_vel_wheel + Constants::DIM * Constants::REAR_RIGHT; }

    inline T* getTyreInitialVelocityFrontLeft() { return _initial_vel_tyre + Constants::DIM * Constants::FRONT_LEFT; }
    inline T* getTyreInitialVelocityFrontRight() { return _initial_vel_tyre + Constants::DIM * Constants::FRONT_RIGHT; }
    inline T* getTyreInitialVelocityRearLeft() { return _initial_vel_tyre + Constants::DIM * Constants::REAR_LEFT; }
    inline T* getTyreInitialVelocityRearRight() { return _initial_vel_tyre + Constants::DIM * Constants::REAR_RIGHT; }

    inline T* getBodyInitialAngularVelocity() { return _initial_ang_vel_body; }

    inline T* getGravityField() { return _gravity; }

    inline T* getBodyInitialPosition() { return _initial_pos_body; }

    inline T* getWheelInitialPositionFrontLeft() { return _initial_pos_wheel + Constants::DIM * Constants::FRONT_LEFT; }
    inline T* getWheelInitialPositionFrontRight() { return _initial_pos_wheel + Constants::DIM * Constants::FRONT_RIGHT; }
    inline T* getWheelInitialPositionRearLeft() { return _initial_pos_wheel + Constants::DIM * Constants::REAR_LEFT; }
    inline T* getWheelInitialPositionRearRight() { return _initial_pos_wheel + Constants::DIM * Constants::REAR_RIGHT; }

    inline T* getTyreInitialPositionFrontLeft() { return _initial_pos_tyre + Constants::DIM * Constants::FRONT_LEFT; }
    inline T* getTyreInitialPositionFrontRight() { return _initial_pos_tyre + Constants::DIM * Constants::FRONT_RIGHT; }
    inline T* getTyreInitialPositionRearLeft() { return _initial_pos_tyre + Constants::DIM * Constants::REAR_LEFT; }
    inline T* getTyreInitialPositionRearRight() { return _initial_pos_tyre + Constants::DIM * Constants::REAR_RIGHT; }

    inline T* getBodyInitialOrientation() { return _initialAngleGlobal; }
    bool getFlagInitialLeg() const { return _initial_leg_flag; }

    bool getUseInterpolation() const { return _interpolation; }

    ALESolver getALESolver() const { return _ALE_solver; }

    T getNewtonTolerance() const { return _newtonTolerance; }

    int getMaxNewtonIterations() const { return _maxNewtonIterations; }

    MBDSolver getUsedSolverForMBD() const { return _MBD_solver; }

    int getMaxNumberOfBroydenIterationForMBD() const { return _max_num_iter; }

    T getToleranceBroydenIterationForMBD() const { return _tolerance; }

    inline T getTimeStepSize() const { return _timestep; }

    inline T getNumberOfTimeIterations() const { return _num_time_iter; }

    inline T getSolutionVectorSize() const { return _solution_dim; }

    inline T* getBodyExternalForce() { return _external_force_body; }

    inline T* getWheelExternalForceFrontLeft() { return _external_force_wheel + Constants::DIM * Constants::FRONT_LEFT; }
    inline T* getWheelExternalForceFrontRight() { return _external_force_wheel + Constants::DIM * Constants::FRONT_RIGHT; }
    inline T* getWheelExternalForceRearLeft() { return _external_force_wheel + Constants::DIM * Constants::REAR_LEFT; }
    inline T* getWheelExternalForceRearRight() { return _external_force_wheel + Constants::DIM * Constants::REAR_RIGHT; }

    inline T* getTyreExternalForceFrontLeft() { return _external_force_tyre + Constants::DIM * Constants::FRONT_LEFT; }
    inline T* getTyreExternalForceFrontRight() { return _external_force_tyre + Constants::DIM * Constants::FRONT_RIGHT; }
    inline T* getTyreExternalForceRearLeft() { return _external_force_tyre + Constants::DIM * Constants::REAR_LEFT; }
    inline T* getTyreExternalForceRearRight() { return _external_force_tyre + Constants::DIM * Constants::REAR_RIGHT; }
    ArbitraryTrajectory<T>* getArbitraryTrajectory() { return _trajectory; }

    inline T getCircularRoadRadius() const { return _profile_radius; }
    inline T* getCircularRoadCenter() { return _profile_center; }

    BoundaryConditionRoad getRoadConditions() const { return _boundary_condition_road; }
#ifdef INTERPOLATION
    const EVAALookup<Constants::floatEVAA>& getLookupStiffness() const { return *_lookupStiffness; }

    const EVAALookup<Constants::floatEVAA>& getLookupDamping() const { return *_lookupDamping; }
#else
    const void* getLookupStiffness() const { return NULL; }

    const void* getLookupDamping() const { return NULL; }
#endif  // INTERPOLATION

    virtual ~MetaDatabase() {
#ifdef INTERPOLATION        
        delete _lookupDamping;
        delete _lookupStiffness;
#endif
        delete _trajectory;
    }

private:
    /** Private constructor for the singleton instance. */
    MetaDatabase() : _lookup_filename("") {}
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
            a[0] = _k_body[0];
            a[1] = _k_tyre[0];
            a[2] = _k_body[1];
            a[3] = _k_tyre[1];
            a[4] = _k_body[2];
            a[5] = _k_tyre[2];
            a[6] = _k_body[3];
            a[7] = _k_tyre[3];

            _lookupStiffness = new EVAALookup<Constants::floatEVAA>(size, a, b, c, l_min, l_max, k, type, order);
            _interpolation = 1;  // to switch from constant to interpolation type

            // damping is /100 from the stiffness for the start
            for (auto j = 0; j < k; j++) {
                // TODO: extract magic number 100?
                a[j] /= 100;
            }

            _lookupDamping = new EVAALookup<Constants::floatEVAA>(size, a, b, c, l_min, l_max, k, type, order);

            delete[] a;
            delete[] _k_body;
            delete[] _k_tyre;
        }
    }

    /**
     * \brief Reads 4 legs which contain each EVAA::Constants::DIM vectors.
     * Reads the vectors of the positions of the legs relative to the center of
     * mass. \param vec XML parser \return storage all components in one vector
     * with 12 elements [fl:XYZ,fr:XYZ,rl:XYZ,rr:XYZ]
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
     * \param[in] vec XML parser
     * \param[out] storage all components in one vector with 4 consecutive
     * elements [fl,fr,rl,rr]
     */
    template <typename S>
    void readLegs(T* storage, S vec) {
        readLegs(storage, 1, vec);
    }

    /**
     * \brief Read 4 legs which contain each 1 T with increment in output
     * storage. Reads paramaters which are given for all the legs. \param[in]
     * vec XML parser \param[in] inc increment in storage \param[out] storage
     * components in one vector at indices 0, inc, 2*inc, 3*inc [fl,fr,rl,rr]
     */
    template <typename S>
    void readLegs(T* storage, size_t inc, S vec) {
        storage[Constants::FRONT_LEFT * inc] = vec.FrontLeft();
        storage[Constants::FRONT_RIGHT * inc] = vec.FrontRight();
        storage[Constants::REAR_LEFT * inc] = vec.ReerLeft();
        storage[Constants::REAR_RIGHT * inc] = vec.ReerRight();
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
    All general simulation parameters and car specific parameters (such as
    geometry and initial conditions)
    */
    T _k_tyre[Constants::NUM_LEGS];
    T _k_body[Constants::NUM_LEGS];
    T _c_tyre[Constants::NUM_LEGS];
    T _c_body[Constants::NUM_LEGS];
    T _l_long[Constants::NUM_LEGS];
    T _l_lat[Constants::NUM_LEGS];
    T _vehicleCIR[Constants::DIM];
    T _mass_body;
    T _I_body[Constants::DIMDIM];
    T _mass[2 * Constants::NUM_LEGS];
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
    ALESolver _ALE_solver;
    T _newtonTolerance;
    int _maxNewtonIterations;

    ArbitraryTrajectory<T>* _trajectory = nullptr;

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

#ifdef INTERPOLATION
    /** Stiffness lookup from the compute engine. */
    EVAALookup<Constants::floatEVAA>* _lookupStiffness;
    /** Damping lookup from the compute engine. */
    EVAALookup<Constants::floatEVAA>* _lookupDamping;
#endif  // INTERPOLATION
};

}  // namespace EVAA