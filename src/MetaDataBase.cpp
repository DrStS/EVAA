/**
 * \file MetaDataBase.cpp
 * This file holds the function definitions of MetaDataBase.
 * \date 04/14/2020
 */

#include "MetaDataBase.h"

#include <stdio.h>
#include <stdlib.h>

#include <fstream>
#include <iostream>
#include <string>

#include "IP_EVAA_XML.h"
#include "LOAD_EVAA_XML.h"
#include "LOOKUP_EVAA_XML.h"

namespace EVAA {

MetaDataBase& MetaDataBase::getDataBase() {
    static MetaDataBase database;
    return database;
}

/**
 * \brief reads the vectors of the positions of the legs relative to the center of mass
 */
template <typename T>
void MetaDataBase::readVectorLegs(double* storage, T vec) {
    readVector(storage, vec.FrontLeft());
    readVector(storage + 3, vec.FrontRight());
    readVector(storage + 6, vec.ReerLeft());
    readVector(storage + 9, vec.ReerRight());
}

/**
 * \brief reads a 3 dim vector
 */
template <typename T>
void MetaDataBase::readVector(double* storage, T vec) {
    storage[0] = vec.x();
    storage[1] = vec.y();
    storage[2] = vec.z();
}

/**
 * \brief reads paramaters which are given for all the legs
 */
template <typename T>
void MetaDataBase::readLegs(double* storage, T vec) {
    storage[0] = vec.FrontLeft();
    storage[1] = vec.FrontRight();
    storage[2] = vec.ReerLeft();
    storage[3] = vec.ReerRight();
}

/**
 * \brief reads quaternion
 */
template <typename T>
void MetaDataBase::readAngles(double* storage, T vec) {
    storage[0] = vec.x();
    storage[1] = vec.y();
    storage[2] = vec.z();
    storage[3] = vec.w();
}

void MetaDataBase::readParameters(const std::string& filename) {
    // Load car parameters

    const auto settings = EVAA_settings(filename, xml_schema::flags::dont_validate);
    const auto twoTrackModel = settings->VehicleXML().TwoTrackModelXML();

    mass_body = twoTrackModel.MassXML().BodyXML();
    readLegs(mass_wheel, twoTrackModel.MassXML().UnsprungMassXML());
    readLegs(mass_tyre, twoTrackModel.MassXML().TyreXML());
    I_body[0] = twoTrackModel.InertiaXML().XX();
    I_body[1] = twoTrackModel.InertiaXML().XY();
    I_body[2] = twoTrackModel.InertiaXML().XZ();
    I_body[3] = twoTrackModel.InertiaXML().YX();
    I_body[4] = twoTrackModel.InertiaXML().YY();
    I_body[5] = twoTrackModel.InertiaXML().YZ();
    I_body[6] = twoTrackModel.InertiaXML().ZX();
    I_body[7] = twoTrackModel.InertiaXML().ZY();
    I_body[8] = twoTrackModel.InertiaXML().ZZ();

    if (twoTrackModel.StiffnessXML().ConstantXML().present()) {
        std::cout << "Take constant stiffness without lookup table" << std::endl;
        readLegs(k_tyre, twoTrackModel.StiffnessXML().ConstantXML().get().TyreXML());
        readLegs(k_body, twoTrackModel.StiffnessXML().ConstantXML().get().BodyXML());
    }
    else {
        readLookupParameters(twoTrackModel.StiffnessXML().LookupTableXML().get().FilePathXML());
    }
    readLegs(c_tyre, twoTrackModel.DampingCoefficientsXML().TyreXML());
    readLegs(c_body, twoTrackModel.DampingCoefficientsXML().BodyXML());
    readLegs(l_long, twoTrackModel.GeometryXML().LongitudinalReferenceToWheelXML());
    readLegs(l_lat, twoTrackModel.GeometryXML().LateralReferenceToWheelXML());
    readVector(vehicleCIR, twoTrackModel.GeometryXML().RelativeCenterOfInstanteneousRotation());
    readLegs(lower_spring_length, twoTrackModel.GeometryXML().SuspensionSpringsXML());
    readLegs(upper_spring_length, twoTrackModel.GeometryXML().TyreSpringsXML());

    // Load initial parameters
    const auto initial = settings->InitialConditionsXML();

    readVector(initial_vel_body, initial.VelocitiesXML().BodyXML());
    readVector(initial_ang_vel_body, initial.VelocitiesXML().angularBodyXML());

    readVectorLegs(initial_vel_wheel, initial.VelocitiesXML().UnsprungMassXML());
    readVectorLegs(initial_vel_tyre, initial.VelocitiesXML().TyreXML());

    readLegs(initial_lower_spring_length, initial.SpringElongationXML().TyreXML());
    readLegs(initial_upper_spring_length, initial.SpringElongationXML().BodyXML());

    readVector(initial_pos_body, initial.PositionXML().BodyXML());
    if (initial.PositionXML().UnsprungMassXML().present()) {
        initial_leg_flag = 1;
        readVectorLegs(initial_pos_wheel, initial.PositionXML().UnsprungMassXML().get());
        readVectorLegs(initial_pos_tyre, initial.PositionXML().TyreXML().get());
    }

    readAngles(initialAngleGlobal, initial.OrientationXML());

    // Load simulation parameters
    const auto simulation = settings->SimulationParametersXML();

    readVector(_gravity, simulation.GeneralSettingsXML().GravityXML());
    _num_time_iter = simulation.GeneralSettingsXML().NumberOfIterationsXML();
    _timestep = simulation.GeneralSettingsXML().TimestepSizeXML();

    max_num_iter = simulation.MultyBodyDynamicsXML().MaximalIterationNumberXML();
    tolerance = simulation.MultyBodyDynamicsXML().ToleranceXML();
    solution_dim = simulation.MultyBodyDynamicsXML().SolutionDimensionXML();
    std::string solver = simulation.MultyBodyDynamicsXML().SolverXML();

    if (solver == "explicit_Euler") {
        MBD_solver = EXPLICIT_EULER;
    }
    else if (solver == "RK4") {
        MBD_solver = RUNGE_KUTTA_4;
    }
    else if (solver == "Broyden_Euler") {
        MBD_solver = BROYDEN_EULER;
    }
    else if (solver == "Broyden_CN") {
        MBD_solver = BROYDEN_CN;
    }
    else if (solver == "Broyden_BDF2") {
        MBD_solver = BROYDEN_BDF2;
    }
    else {
        throw std::logic_error("Wrong MBD-solver in XML: " + solver +
                               ". Must be one of: "
                               "explicit_Euler, RK4, Broyden_Euler, Broyden_CN, Broyden_BDF2");
    }
}

void MetaDataBase::readLoadParameters(const std::string& filename) {
    // Load external parameters
    const auto load_data = EVAA_load_module(filename, xml_schema::flags::dont_validate);


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
       _trajectory = new arbitraryTrajectory<double>(_num_time_iter, _timestep,
           load_data->roadProfile().arbitraryRoadProfile().get().verticalProfile().sinusoidalProfile().rightTyre().amplitude(),
           load_data->roadProfile().arbitraryRoadProfile().get().verticalProfile().sinusoidalProfile().leftTyre().amplitude(),
           load_data->roadProfile().arbitraryRoadProfile().get().verticalProfile().sinusoidalProfile().rightTyre().period(),
           load_data->roadProfile().arbitraryRoadProfile().get().verticalProfile().sinusoidalProfile().leftTyre().period(),
           load_data->roadProfile().arbitraryRoadProfile().get().verticalProfile().sinusoidalProfile().rightTyre().shift(),
           load_data->roadProfile().arbitraryRoadProfile().get().verticalProfile().sinusoidalProfile().leftTyre().shift());  // TODO: free memory without leaks

       double initialVelocity = sqrt(_initial_vel_body[0] * _initial_vel_body[0] +
                                     _initial_vel_body[1] * _initial_vel_body[1]);

       std::vector<double> wayPointsX;
       std::vector<double> wayPointsY;
       std::vector<double> wayPointsTimes;

       size_t numWayPoints = 0;

       for (auto    i = load_data->roadProfile().arbitraryRoadProfile().get().horizontalProfile().wayPoint().begin();
                    i < load_data->roadProfile().arbitraryRoadProfile().get().horizontalProfile().wayPoint().end(); ++i) {
           numWayPoints++;
           wayPointsX.push_back(i->X());
           wayPointsY.push_back(i->Y());
           wayPointsTimes.push_back(i->time());
       }

       _trajectory->interpolateRoadPoints(numWayPoints, wayPointsX.data(), wayPointsY.data(), wayPointsTimes.data(), initialVelocity);

       _trajectory->calculateTyreShifts( _l_long[Constants::FRONT_LEFT], _l_long[Constants::FRONT_RIGHT],
                                         _l_long[Constants::REAR_LEFT], _l_long[Constants::REAR_RIGHT]);



   }

    else {
        throw std::logic_error(
            "Wrong boundary conditions. Implemented so far: circle, fixed, nonfixed.");
    }

    readVectorLegs(_external_force_tyre, load_data->forces().forceTyre());
    readVectorLegs(_external_force_wheel, load_data->forces().forceWheel());
    readVector(_external_force_body, load_data->forces().forceBody());
}

void MetaDataBase::readLookupParameters(const std::string& filename) {
    IO::checkFileExists(filename);
    const auto lookup = LookupHandler(filename, xml_schema::flags::dont_validate);

    const auto generator = lookup->LookupTableGenerator();
    if (generator.present()) {
        double *a, *k_body, *k_tyre;
        double b, c, l_min, l_max;
        int size, k, type, order;

        // TODO: statics or mkl_malloc;
        a = new (double[8]);
        k_body = new (double[4]);
        k_tyre = new (double[4]);

        std::cout << "Generate look up table from parameters." << std::endl;

        size = generator->Size();
        b = generator->TableParameters().b();
        c = generator->TableParameters().c();
        l_min = generator->Range().l_min();
        l_max = generator->Range().l_max();
        k = generator->InterpolationMethod().k();
        type = generator->InterpolationMethod().type();
        order = generator->InterpolationMethod().order();

        readLegs(k_body, generator->Magnitude().Body());
        readLegs(k_tyre, generator->Magnitude().Tyre());
        a[0] = k_body[0];
        a[1] = k_tyre[0];
        a[2] = k_body[1];
        a[3] = k_tyre[1];
        a[4] = k_body[2];
        a[5] = k_tyre[2];
        a[6] = k_body[3];
        a[7] = k_tyre[3];

        _lookupStiffness = std::make_unique<EVAALookup<Constants::floatEVAA>>(
            size, a, b, c, l_min, l_max, k, type, order);
        interpolation = 1;  // to switch from constant to interpolation type

        // damping is /100 from the stiffness for the start
        for (auto j = 0; j < k; j++) {
            a[j] /= 100;
        }
        _lookupDamping = std::make_unique<EVAALookup<Constants::floatEVAA>>(size, a, b, c, l_min,
                                                                            l_max, k, type, order);
        delete[] a;
        delete[] k_body;
        delete[] k_tyre;
    }
}

}  // namespace EVAA
