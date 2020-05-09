// TODO: Copyright header

#pragma once

#include <string>

#include "ArbitraryTrajectory.h"
#include "Car.h"
#include "Constants.h"
#include "MathLibrary.h"
#include "MetaDataBase.h"

namespace EVAA {

/**
 * Parent class for generic road profiles
 */
template <typename T>
class Profile {
protected:
    std::string Name;
    T* massOfComponents;  // [GC, W1, T1, ... ]						size: 9

public:
    Profile() { massOfComponents = Math::malloc<T>(Constants::VEC_DIM); }

    /**
     * In case the initial conditions have to be specific to a road profile
     * \param Car
     */
    virtual void ApplyProfileInitialCondition(Car<T>* carObj) = 0;

    std::string GetName() { return Name; }

    virtual ~Profile() { Math::free<T>(massOfComponents); };
};

/**Eulerian Profile*/
template <typename T>
class Euler : public Profile<T> {
protected:
    T* eulerianVelocity;  // [GC: Z, Angle: XY, W1: Z, T1: Z, ... ]   size: 11
    T gravity;

public:
    Euler(T& g) : Profile<T>() {
        eulerianVelocity = Math::malloc<T>(Constants::DOF);
        gravity = g;
    }

    virtual ~Euler() { Math::free<T>(eulerianVelocity); }

    /**
    Add Gravity force in the z direction
    */
    void AddGravity(Car<T>* carObj, T* profileInducedForce) {
        profileInducedForce[0] += carObj->getMassComponents()[0] * gravity;
#pragma loop(ivdep)
        for (auto i = Constants::DIM; i < Constants::DOF; ++i) {
            profileInducedForce[i] +=
                carObj->getMassComponents()[i - (Constants::DIM - 1)] * gravity;
        }
    }

    /**
     * Get external force acting on the car system in the Eulerian Frame
     * \param Car
     * \return profileInducedForce forces and torque acting on each component [GC: Z, GC(Torque):XY,
     * W1: Z, T1: Z, ...]
     */
    virtual void GetProfileForceEulerian(const size_t& _iterationCount, Car<T>* carObj,
                                         T* profileInducedForce) = 0;
};

/** Lagrangian profiles*/
template <typename T>
class Lagrange : public Profile<T> {
private:
    // bla bla

public:
    Lagrange() : Profile<T>() {}

    /**
     * Get external force acting on the car system in the Lagrangian Frame
     * \param Car
     * \return profileInducedForce forces acting on each component [GC: XY, W1: XY, T1: XY, ...]
     * \return reactionOnTyre reaction force on the tyre induced by profile forcce [T: XY, T2: XY,
     * T3: XY, T4: XY]
     */
    virtual void GetProfileForceLagrangian(const size_t& _iterationCount, Car<T>* carObj,
                                           T* profileInducedForce, T* reactionOnTyre) = 0;

    /**
     * Get external torque acting on the car system in Lagrangian Frame
     * \param Car
     * \return externalTorque torque acting on the car system [GC: Z]
     */
    virtual void GetProfileTorqueLagrangian(const size_t& _iterationCount, Car<T>* carObj,
                                            T* externalTorque) = 0;
    virtual ~Lagrange() {}
};

/** Follow a circular road profile with radius R and center C */
template <typename T>
class Circular : public Lagrange<T> {
private:
    /** Center of the Circle (if the motion is a Circle) */
    T* centerOfCircle = NULL;

    /** Radius of the Circle (if the motion is a Circle) */
    T inverseRadiusOfCircle = std::numeric_limits<T>::infinity();

    // vectors used in computations inside the methods
    T* tangentialVelocityDirection;  // arrow from center towards object

    // distance vector between center of circle and the car =
    // Positions of points from the car vs Center of Circle, seen as 0
    T* radiusVector;  // 3 x 9

public:
    Circular(T* circleCenter) : Lagrange<T>() {
        Name = "Circular";

        // Position
        centerOfCircle = Math::malloc<T>(Constants::DIM);
        Math::copy<T>(Constants::DIM, circleCenter, 1, centerOfCircle, 1);

        // auxiliary vectors
        tangentialVelocityDirection = Math::calloc<T>(Constants::DIM);
        radiusVector = Math::malloc<T>((Constants::DIM - 1) * Constants::VEC_DIM);
    }

    virtual ~Circular() {
        Math::free<T>(centerOfCircle);
        Math::free<T>(tangentialVelocityDirection);
        Math::free<T>(radiusVector);
    }

    /**
     * In case the initial conditions have to be specific to a road profile
     * \param Car
     */
    virtual void ApplyProfileInitialCondition(Car<T>* carObj) {
        // Compute Radius Vector based on the position of origin of euler frame
        radiusVector[0] = carObj->getCurrentPositionLagrangian()[0] - centerOfCircle[0];
        radiusVector[1] = carObj->getCurrentPositionLagrangian()[1] - centerOfCircle[1];
        inverseRadiusOfCircle =
            1.0 / Math::nrm2<T>(Constants::DIM - 1, radiusVector, Constants::INCX);
        Math::scal<T>(Constants::DIM - 1, inverseRadiusOfCircle, radiusVector, Constants::INCX);

        // Compute the tangential direction of the velocity to compute angular velocity
        tangentialVelocityDirection[0] = -radiusVector[1];
        tangentialVelocityDirection[1] = radiusVector[0];

        // Compute Angular velocity using w = v * sin(theta) / r = dot_product(v,
        // tangential_unit_vector) / r
        carObj->_currentAngularVelocityLagrangian =
            Math::dot<T>(Constants::DIM - 1, tangentialVelocityDirection, Constants::INCX,
                         carObj->getCurrentVelocityLagrangian(), Constants::INCX);
        carObj->_currentAngularVelocityLagrangian *= inverseRadiusOfCircle;
        // Copy this to initial vector too
        carObj->setInitialAngularVelocityGlobalZ(carObj->getCurrentAngularVelocityLagrangian());

        // apply angle condition on z angle
        carObj->setInitialAngleGlobalZ(carObj->_currentAngleLagrangian =
                                           std::atan2(carObj->_currentVelocityLagrangian[1],
                                                      carObj->_currentVelocityLagrangian[0]));
    }

    virtual void GetProfileTorqueLagrangian(const size_t& _iterationCount, Car<T>* carObj,
                                            T* externalTorque) {
        *externalTorque = 0;
    }

    /**
     * Computes centrifugal force component in the Lagrangian Frame
     * \param[in] distanceFromCenter: Distance of each component of the car from center of circle
     * \param[in] velocityLagrangian: Velocity vector of the components of the car in the Lagrangian
     * Frame \param[in] mass: Mass vector containing mass of each component of the car in the array
     * of form [CG, W1, T1, ... ] \param[out] centrifugalForce: Computed centrifugal force on each
     * component of the car
     */
    void ComputeCentrifugalForceLagrangian(const T* distanceFromCenter, const T* velocityLagrangian,
                                           const T* mass, T* centrifugalForce) {
        // Copy the radius vector to the force to extract direction
        Math::copy<T>((Constants::DIM - 1) * Constants::VEC_DIM, distanceFromCenter,
                      Constants::INCX, centrifugalForce, Constants::INCX);
        T velocityMagnitude;
        for (auto i = 0; i < Constants::VEC_DIM; ++i) {
            inverseRadiusOfCircle =
                1.0 / Math::nrm2<T>((Constants::DIM - 1),
                                    centrifugalForce + (Constants::DIM - 1) * i, Constants::INCX);

            Math::scal<T>((Constants::DIM - 1), inverseRadiusOfCircle,
                          centrifugalForce + (Constants::DIM - 1) * i, Constants::INCX);
            // Cross product to get tangential velocity component
            tangentialVelocityDirection[0] = *(centrifugalForce + (Constants::DIM - 1) * i + 1);
            tangentialVelocityDirection[1] = -*(centrifugalForce + (Constants::DIM - 1) * i);
            velocityMagnitude =
                Math::dot<T>((Constants::DIM - 1), velocityLagrangian + (Constants::DIM - 1) * i,
                             Constants::INCX, tangentialVelocityDirection, Constants::INCX);
            // Force computation
            Math::scal<T>(Constants::DIM - 1,
                          mass[i] * velocityMagnitude * velocityMagnitude * inverseRadiusOfCircle,
                          centrifugalForce + (Constants::DIM - 1) * i, Constants::INCX);
        }
    }

    /**
     * Get Lagrangian forces acting due to the road profile
     */
    virtual void GetProfileForceLagrangian(const size_t& _iterationCount, Car<T>* carObj,
                                           T* profileInducedForce, T* reactionOnTyre) {
        carObj->ComputeDisplacementToPointLagrangian(centerOfCircle, radiusVector);
        // compute centrifugal force on each component
        ComputeCentrifugalForceLagrangian(radiusVector, carObj->getCurrentVelocityLagrangian(),
                                          carObj->getMassComponents(), profileInducedForce);

        // Compute the reaction on tyre due to centrifugal force
        reactionOnTyre[0] = -profileInducedForce[0];
        reactionOnTyre[1] = -profileInducedForce[1];
        for (auto i = 1; i < Constants::VEC_DIM; ++i) {
            reactionOnTyre[0] -= profileInducedForce[(Constants::DIM - 1) * i];
            reactionOnTyre[1] -= profileInducedForce[(Constants::DIM - 1) * i + 1];
        }
    }
};  // Circular

template <typename T>
class Arbitrary : public Lagrange<T> {
private:
    ArbitraryTrajectory<T>* _trajectory;

public:
    Arbitrary(ArbitraryTrajectory<T>* trajectory) : Lagrange<T>(), _trajectory(trajectory) {
        Name = "Arbitrary";
    }

    /**
     * Get external force acting on the car system in the Lagrangian Frame
     * \param Car
     * \return profileInducedForce forces acting on each component [GC: XY, W1: XY, T1: XY, ...]
     * \return reactionOnTyre reaction force on the tyre induced by profile forcce [T: XY, T2: XY,
     * T3: XY, T4: XY]
     */
    virtual void GetProfileForceLagrangian(const size_t& _iterationCount, Car<T>* carObj,
                                           T* profileInducedForce, T* reactionOnTyre) {
        _trajectory->getLagrangianForcesCenterOfGravity(_iterationCount, profileInducedForce,
                                                        carObj->getMassComponents()[0]);

        _trajectory->getLagrangianForcesFrontLeft(
            _iterationCount, carObj->getMassComponents()[1 + 2 * Constants::FRONT_LEFT],
            profileInducedForce + (1 + 2 * Constants::FRONT_LEFT) * (Constants::DIM - 1));
        _trajectory->getLagrangianForcesFrontLeft(
            _iterationCount, carObj->getMassComponents()[2 + 2 * Constants::FRONT_LEFT],
            profileInducedForce + (2 + 2 * Constants::FRONT_LEFT) * (Constants::DIM - 1));

        _trajectory->getLagrangianForcesFrontRight(
            _iterationCount, carObj->getMassComponents()[1 + 2 * Constants::FRONT_RIGHT],
            profileInducedForce + (1 + 2 * Constants::FRONT_RIGHT) * (Constants::DIM - 1));
        _trajectory->getLagrangianForcesFrontRight(
            _iterationCount, carObj->getMassComponents()[2 + 2 * Constants::FRONT_RIGHT],
            profileInducedForce + (2 + 2 * Constants::FRONT_RIGHT) * (Constants::DIM - 1));

        _trajectory->getLagrangianForcesRearLeft(
            _iterationCount, carObj->getMassComponents()[1 + 2 * Constants::REAR_LEFT],
            profileInducedForce + (1 + 2 * Constants::REAR_LEFT) * (Constants::DIM - 1));
        _trajectory->getLagrangianForcesRearLeft(
            _iterationCount, carObj->getMassComponents()[2 + 2 * Constants::REAR_LEFT],
            profileInducedForce + (2 + 2 * Constants::REAR_LEFT) * (Constants::DIM - 1));

        _trajectory->getLagrangianForcesRearRight(
            _iterationCount, carObj->getMassComponents()[1 + 2 * Constants::REAR_RIGHT],
            profileInducedForce + (1 + 2 * Constants::REAR_RIGHT) * (Constants::DIM - 1));
        _trajectory->getLagrangianForcesRearRight(
            _iterationCount, carObj->getMassComponents()[2 + 2 * Constants::REAR_RIGHT],
            profileInducedForce + (2 + 2 * Constants::REAR_RIGHT) * (Constants::DIM - 1));

        // Compute the reaction on tyre due to centrifugal force
        reactionOnTyre[0] = profileInducedForce[0];
        reactionOnTyre[1] = profileInducedForce[1];
        for (auto i = 1; i < Constants::VEC_DIM; ++i) {
            reactionOnTyre[0] += profileInducedForce[(Constants::DIM - 1) * i];
            reactionOnTyre[1] += profileInducedForce[(Constants::DIM - 1) * i + 1];
        }

        Math::scal<T>((Constants::DIM - 1) * Constants::VEC_DIM, -1, profileInducedForce,
                      Constants::INCX);
    }

    /**
     * Get external torque acting on the car system in Lagrangian Frame
     * \param Car
     * \return externalTorque torque acting on the car system [GC: Z]
     */
    virtual void GetProfileTorqueLagrangian(const size_t& _iterationCount, Car<T>* carObj,
                                            T* externalTorque) {
        *externalTorque =
            _trajectory->getLagrangianTorque(_iterationCount, carObj->getMomentOfInertia()[8]);
    }
    virtual void ApplyProfileInitialCondition(Car<T>* carObj) {
        T* angle = &(carObj->_currentAngleLagrangian);
        T* wc = &(carObj->_currentAngularVelocityLagrangian);
        T* pcc = carObj->_currentPositionLagrangian;
        T* vcc = carObj->_currentVelocityLagrangian;
        _trajectory->updateInitialConditionsLagrange(angle, wc, pcc, vcc);
    }
    virtual ~Arbitrary() {}
};

/** Follow a circular road profile with radius R and center C */
template <typename T>
class Straight : public Lagrange<T> {
private:
public:
    Straight() : Lagrange<T>() { Name = "Straight"; }

    virtual ~Straight() {}

    /**
     * In case the initial conditions have to be specific to a road profile
     * \param Car
     */
    virtual void ApplyProfileInitialCondition(Car<T>* carObj) {}

    virtual void GetProfileTorqueLagrangian(const size_t& _iterationCount, Car<T>* carObj,
                                            T* externalTorque) {
        *externalTorque = 0;
    }

    /**
     * Get Lagrangian forces acting due to the road profile
     */
    virtual void GetProfileForceLagrangian(const size_t& _iterationCount, Car<T>* carObj,
                                           T* profileInducedForce, T* reactionOnTyre) {
        Math::scal<T>(Constants::VEC_DIM * (Constants::DIM - 1), 0, profileInducedForce,
                      Constants::INCX);
        Math::scal<T>((Constants::DIM - 1), 0, reactionOnTyre, Constants::INCX);
    }
};  // Straight Road

/* Eulerian profile type with fixed boundary condition */
template <typename T>
class Fixed : public Euler<T> {
private:
    // bla bla
public:
    Fixed(T& g) : Euler<T>(g) { Name = "Fixed"; }
    virtual ~Fixed() {}

    virtual void GetProfileForceEulerian(const size_t& _iterationCount, Car<T>* carObj,
                                         T* profileInducedForce) {
        // TODO optimize it
        Math::scal<T>(Constants::DOF, 0, profileInducedForce, Constants::INCX);
        AddGravity(carObj, profileInducedForce);
#pragma loop(ivdep)
        for (auto i = 0; i < Constants::NUM_LEGS; ++i) {
            profileInducedForce[Constants::TYRE_INDEX_EULER[i]] =  //
                carObj->getkVec()[2 * i + 1] *
                    (carObj->getCurrentDisplacementTwoTrackModel()[Constants::TYRE_INDEX_EULER[i]] -
                     carObj->getCurrentDisplacementTwoTrackModel()[Constants::TYRE_INDEX_EULER[i] -
                                                                   1])  //
                +
                carObj->getdVec()[2 * i + 1] *
                    (carObj->getCurrentVelocityTwoTrackModel()[Constants::TYRE_INDEX_EULER[i]] -
                     carObj->getCurrentVelocityTwoTrackModel()[Constants::TYRE_INDEX_EULER[i] - 1]);
        }
    }

    virtual void ApplyProfileInitialCondition(Car<T>* carObj) {


        T interpolationFactor[Constants::NUM_LEGS];
        for (auto i = 0; i < Constants::NUM_LEGS; ++i) {
            interpolationFactor[i] =
                0.5 * carObj->_lenLong[i] / (carObj->_lenLong[0] + carObj->_lenLong[2]);
        }

        // update displacements
        T& centerPosition =
            carObj->_currentDisplacementTwoTrackModel[0];
        T& unexcitedCenterPosition = carObj->_unexcitedPositionTwoTrackModel[0];
        for (auto i = 0; i < Constants::NUM_LEGS; ++i) {
            T& wheelPosition =
                carObj->_currentDisplacementTwoTrackModel[Constants::TYRE_INDEX_EULER[i] - 1];
            T& unexcitedWheelPosition =
                carObj->_unexcitedPositionTwoTrackModel[Constants::TYRE_INDEX_EULER[i] - 1];
            T& tyrePosition =
                carObj->_currentDisplacementTwoTrackModel[Constants::TYRE_INDEX_EULER[i]];
            T& unexcitedTyrePosition =
                carObj->_unexcitedPositionTwoTrackModel[Constants::TYRE_INDEX_EULER[i]];


            wheelPosition -= tyrePosition;
            centerPosition -= interpolationFactor[i] * tyrePosition;

            unexcitedTyrePosition += tyrePosition;
            unexcitedWheelPosition += tyrePosition;
            unexcitedCenterPosition += interpolationFactor[i] * tyrePosition;

            tyrePosition = 0;
        }

        // update velocities
        T& centerVelocity = carObj->_currentVelocityTwoTrackModel[0];
        for (auto i = 0; i < Constants::NUM_LEGS; ++i) {
            T& wheelVelocity =
                carObj->_currentVelocityTwoTrackModel[Constants::TYRE_INDEX_EULER[i] - 1];
            T& tyreVelocity = carObj->_currentVelocityTwoTrackModel[Constants::TYRE_INDEX_EULER[i]];
            wheelVelocity -= tyreVelocity;
            centerVelocity -= interpolationFactor[i] * tyreVelocity;
            tyreVelocity = 0;
        }

        // angle thingy from Shubham

        carObj->UpdateLengthsTwoTrackModel();
    }
};

/* Eulerian profile type with fixed boundary condition */
template <typename T>
class Nonfixed : public Euler<T> {
private:
    /// Bla Bla
public:
    Nonfixed(T& g) : Euler<T>(g) { Name = "Nonfixed"; }
    virtual ~Nonfixed() {}
    virtual void GetProfileForceEulerian(const size_t& _iterationCount, Car<T>* carObj,
                                         T* profileInducedForce) {
        Math::scal<T>(Constants::DOF, 0, profileInducedForce, Constants::INCX);
        AddGravity(carObj, profileInducedForce);
    }
    virtual void ApplyProfileInitialCondition(Car<T>* carObj) {}
};

template <typename T>
class Sinusoidal : public Euler<T> {
private:
    ArbitraryTrajectory<T>* _trajectory;

public:
    Sinusoidal(ArbitraryTrajectory<T>* trajectory, T& g) : Euler<T>(g), _trajectory(trajectory) {
        Name = "Sinusoidal";
    }
    virtual ~Sinusoidal() {}

    virtual void GetProfileForceEulerian(const size_t& _iterationCount, Car<T>* carObj,
                                         T* profileInducedForce) {
        // TODO optimize it
        Math::scal<T>(Constants::DOF, 0, profileInducedForce, Constants::INCX);
        AddGravity(carObj, profileInducedForce);
        profileInducedForce[Constants::TYRE_INDEX_EULER[Constants::FRONT_LEFT]] =
            _trajectory->getVerticalRoadForcesFrontLeft(
                _iterationCount, carObj->getMassComponents()[2 + 2 * Constants::FRONT_LEFT]);
        profileInducedForce[Constants::TYRE_INDEX_EULER[Constants::FRONT_RIGHT]] =
            _trajectory->getVerticalRoadForcesFrontRight(
                _iterationCount, carObj->getMassComponents()[2 + 2 * Constants::FRONT_RIGHT]);
        profileInducedForce[Constants::TYRE_INDEX_EULER[Constants::REAR_LEFT]] =
            _trajectory->getVerticalRoadForcesRearLeft(
                _iterationCount, carObj->getMassComponents()[2 + 2 * Constants::REAR_LEFT]);
        profileInducedForce[Constants::TYRE_INDEX_EULER[Constants::REAR_RIGHT]] =
            _trajectory->getVerticalRoadForcesRearRight(
                _iterationCount, carObj->getMassComponents()[2 + 2 * Constants::REAR_RIGHT]);

        for (auto i = 0; i < Constants::NUM_LEGS; ++i) {
            profileInducedForce[Constants::TYRE_INDEX_EULER[i]] +=  //
                carObj->getkVec()[2 * i + 1] *
                    (carObj->getCurrentDisplacementTwoTrackModel()[Constants::TYRE_INDEX_EULER[i]] -
                     carObj->getCurrentDisplacementTwoTrackModel()[Constants::TYRE_INDEX_EULER[i] -
                                                                   1])  //
                +
                carObj->getdVec()[2 * i + 1] *
                    (carObj->getCurrentVelocityTwoTrackModel()[Constants::TYRE_INDEX_EULER[i]] -
                     carObj->getCurrentVelocityTwoTrackModel()[Constants::TYRE_INDEX_EULER[i] - 1]);
        }
    }

    virtual void ApplyProfileInitialCondition(Car<T>* carObj) {
        T* angle = &(carObj->getCurrentDisplacementTwoTrackModel()[1]);
        T* wc = &(carObj->_currentVelocityTwoTrackModel[1]);

        T* pcc = &(carObj->getCurrentDisplacementTwoTrackModel()[0]);
        T* vcc = &(carObj->_currentVelocityTwoTrackModel[0]);

        T* pw_fl =
            &(carObj->getCurrentDisplacementTwoTrackModel()[Constants::TYRE_INDEX_EULER[0] - 1]);
        T* pw_fr =
            &(carObj->getCurrentDisplacementTwoTrackModel()[Constants::TYRE_INDEX_EULER[1] - 1]);
        T* pw_rl =
            &(carObj->getCurrentDisplacementTwoTrackModel()[Constants::TYRE_INDEX_EULER[2] - 1]);
        T* pw_rr =
            &(carObj->getCurrentDisplacementTwoTrackModel()[Constants::TYRE_INDEX_EULER[3] - 1]);

        T* pt_fl = &(carObj->getCurrentDisplacementTwoTrackModel()[Constants::TYRE_INDEX_EULER[0]]);
        T* pt_fr = &(carObj->getCurrentDisplacementTwoTrackModel()[Constants::TYRE_INDEX_EULER[1]]);
        T* pt_rl = &(carObj->getCurrentDisplacementTwoTrackModel()[Constants::TYRE_INDEX_EULER[2]]);
        T* pt_rr = &(carObj->getCurrentDisplacementTwoTrackModel()[Constants::TYRE_INDEX_EULER[3]]);

        T* vw_fl = &(carObj->_currentVelocityTwoTrackModel[Constants::TYRE_INDEX_EULER[0] - 1]);
        T* vw_fr = &(carObj->_currentVelocityTwoTrackModel[Constants::TYRE_INDEX_EULER[1] - 1]);
        T* vw_rl = &(carObj->_currentVelocityTwoTrackModel[Constants::TYRE_INDEX_EULER[2] - 1]);
        T* vw_rr = &(carObj->_currentVelocityTwoTrackModel[Constants::TYRE_INDEX_EULER[3] - 1]);

        T* vt_fl = &(carObj->_currentVelocityTwoTrackModel[Constants::TYRE_INDEX_EULER[0]]);
        T* vt_fr = &(carObj->_currentVelocityTwoTrackModel[Constants::TYRE_INDEX_EULER[1]]);
        T* vt_rl = &(carObj->_currentVelocityTwoTrackModel[Constants::TYRE_INDEX_EULER[2]]);
        T* vt_rr = &(carObj->_currentVelocityTwoTrackModel[Constants::TYRE_INDEX_EULER[3]]);

        _trajectory->updateInitialConditionsEuler(angle, wc, pcc, vcc, pt_fl, pt_fr, pt_rl, pt_rr,
                                                  pw_fl, pw_fr, pw_rl, pw_rr, vt_fl, vt_fr, vt_rl,
                                                  vt_rr, vw_fl, vw_fr, vw_rl, vw_rr);

        // calculate new unexcited positions
        //tyres
        carObj->_unexcitedPositionTwoTrackModel[Constants::TYRE_INDEX_EULER[Constants::FRONT_LEFT]] = *pt_fl;
        carObj->_unexcitedPositionTwoTrackModel[Constants::TYRE_INDEX_EULER[Constants::FRONT_RIGHT]] = *pt_fr;
        carObj->_unexcitedPositionTwoTrackModel[Constants::TYRE_INDEX_EULER[Constants::REAR_LEFT]] = *pt_rl;
        carObj->_unexcitedPositionTwoTrackModel[Constants::TYRE_INDEX_EULER[Constants::REAR_RIGHT]] = *pt_rr;

        //wheels
        carObj->_unexcitedPositionTwoTrackModel[Constants::TYRE_INDEX_EULER[Constants::FRONT_LEFT] - 1] =
            *pt_fl + carObj->_unexcitedSpringsLength[2*Constants::FRONT_LEFT+1];
        carObj->_unexcitedPositionTwoTrackModel[Constants::TYRE_INDEX_EULER[Constants::FRONT_RIGHT] - 1] =
            *pt_fr + carObj->_unexcitedSpringsLength[2*Constants::FRONT_RIGHT+1];
        carObj->_unexcitedPositionTwoTrackModel[Constants::TYRE_INDEX_EULER[Constants::REAR_LEFT] - 1] =
            *pt_rl + carObj->_unexcitedSpringsLength[2*Constants::REAR_LEFT+1];
        carObj->_unexcitedPositionTwoTrackModel[Constants::TYRE_INDEX_EULER[Constants::REAR_RIGHT] - 1] =
            *pt_rr + carObj->_unexcitedSpringsLength[2*Constants::REAR_RIGHT+1];

        // CoG
        carObj->_unexcitedPositionTwoTrackModel[0] =
            0.5 * carObj->_lenLong[Constants::FRONT_LEFT] /
                (carObj->_lenLong[Constants::FRONT_LEFT] + carObj->_lenLong[Constants::REAR_LEFT]) *
                (carObj->_unexcitedPositionTwoTrackModel
                     [Constants::TYRE_INDEX_EULER[Constants::FRONT_LEFT] - 1] +
                 carObj->_unexcitedSpringsLength[2 * Constants::FRONT_LEFT] +
                 carObj->_unexcitedPositionTwoTrackModel
                     [Constants::TYRE_INDEX_EULER[Constants::FRONT_RIGHT] - 1] +
                 carObj->_unexcitedSpringsLength[2 * Constants::FRONT_RIGHT]) +
            0.5 * carObj->_lenLong[Constants::REAR_LEFT] /
                (carObj->_lenLong[Constants::FRONT_LEFT] + carObj->_lenLong[Constants::REAR_LEFT]) *
                (carObj->_unexcitedPositionTwoTrackModel
                     [Constants::TYRE_INDEX_EULER[Constants::REAR_LEFT] - 1] +
                 carObj->_unexcitedSpringsLength[2 * Constants::REAR_LEFT] +
                 carObj->_unexcitedPositionTwoTrackModel
                     [Constants::TYRE_INDEX_EULER[Constants::REAR_RIGHT] - 1] +
                 carObj->_unexcitedSpringsLength[2 * Constants::REAR_RIGHT]);   


        // calculate the new displacements
        //tyres
        *pt_fl = 0;
        *pt_fr = 0;
        *pt_rl = 0;
        *pt_rr = 0;

        //wheels
        *pw_fl = carObj->_unexcitedPositionTwoTrackModel
                [Constants::TYRE_INDEX_EULER[Constants::FRONT_LEFT] - 1] - *pw_fl;
        *pw_fr = carObj->_unexcitedPositionTwoTrackModel
                [Constants::TYRE_INDEX_EULER[Constants::FRONT_RIGHT] - 1] - *pw_fr;
        *pw_rl = carObj->_unexcitedPositionTwoTrackModel
                [Constants::TYRE_INDEX_EULER[Constants::REAR_LEFT] - 1] - *pw_rl;
        *pw_rr = carObj->_unexcitedPositionTwoTrackModel
                [Constants::TYRE_INDEX_EULER[Constants::REAR_RIGHT] - 1] - *pw_rr;

        //CoG
        *pcc = carObj->_unexcitedPositionTwoTrackModel[0] - *pcc;

        
        // velocities
        // wheel
        *vw_fl -= *vt_fl;
        *vw_fr -= *vt_fr;
        *vw_rl -= *vt_rl;
        *vw_rr -= *vt_rr;

        // center
        *vcc -= 0.5 * carObj->_lenLong[Constants::FRONT_LEFT] /
                (carObj->_lenLong[Constants::FRONT_LEFT] + carObj->_lenLong[Constants::REAR_LEFT]) *
                (*vt_fl + *vt_fr) +
            0.5 * carObj->_lenLong[Constants::REAR_LEFT] /
                (carObj->_lenLong[Constants::FRONT_LEFT] + carObj->_lenLong[Constants::REAR_LEFT]) *
                (*vt_rl + *vt_rr); 

        // tyres
        *vt_fl = 0;
        *vt_fr = 0;
        *vt_rl = 0;
        *vt_rr = 0;

        // angle thingy from Shubham

        carObj->UpdateLengthsTwoTrackModel();

    }
};

}  // namespace EVAA
