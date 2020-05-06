// TODO: Copyright header

#pragma once

#include <string>

#include "Car.h"
#include "Constants.h"
#include "MathLibrary.h"

namespace EVAA {

/**
 * Parent class for generic road profiles
 */	
template <typename T>
class Profile {
protected:
    std::string Name;
	T* massOfComponents;            // [GC, W1, T1, ... ]						size: 9

public:

	Profile() {
		massOfComponents = Math::malloc<T>(Constants::VEC_DIM);
	}

	/**
	 * In case the initial conditions have to be specific to a road profile
	 * \param Car
	 */
	virtual void ApplyProfileInitialCondition(Car<T>* carObj) = 0;

	std::string GetName() { return Name; }

    virtual ~Profile(){		
		Math::free<T>(massOfComponents);
	};
};

/**Eulerian Profile*/
template <typename T>
class Euler : public Profile<T> {
protected:
	T* eulerianVelocity;			// [GC: Z, Angle: XY, W1: Z, T1: Z, ... ]   size: 11
	T gravity;
public:
	Euler(T& g) : Profile<T>() {
		eulerianVelocity = Math::malloc<T>(Constants::DOF);
		gravity = g;
	}

	virtual ~Euler() {
		Math::free<T>(eulerianVelocity);
	}

	/**
	Add Gravity force in the z direction
	*/
	void AddGravity(Car<T>* carObj, T* profileInducedForce) {
		profileInducedForce[0] += carObj->getMassComponents()[0] * gravity;
	#pragma loop(ivdep)
		for (auto i = Constants::DIM; i < Constants::DOF; ++i) {
			profileInducedForce[i] += carObj->getMassComponents()[i] * gravity;
		}
	}

	/**
	 * Get external force acting on the car system in the Eulerian Frame
	 * \param Car
	 * \return profileInducedForce forces and torque acting on each component [GC: Z, GC(Torque):XY, W1: Z, T1: Z, ...]
	 */
	virtual void GetProfileForceEulerian(const size_t& _iterationCount, Car<T>* carObj, T* profileInducedForce) = 0;

};

/** Lagrangian profiles*/
template <typename T>
class Lagrange : public Profile<T> {

public:
	Lagrange(): Profile<T>(){}

	/**
	 * Get external force acting on the car system in the Lagrangian Frame
	 * \param Car
	 * \return profileInducedForce forces acting on each component [GC: XY, W1: XY, T1: XY, ...]
	 * \return reactionOnTyre reaction force on the tyre induced by profile forcce [T: XY, T2: XY, T3: XY, T4: XY]
	 */
	virtual void GetProfileForceLagrangian(const size_t& _iterationCount, Car<T>* carObj, T* profileInducedForce, T* reactionOnTyre) = 0;

	/**
	 * Get external torque acting on the car system in Lagrangian Frame
	 * \param Car
	 * \return externalTorque torque acting on the car system [GC: Z]
	 */
	virtual void GetProfileTorqueLagrangian(const size_t& _iterationCount, Car<T>* carObj, T* externalTorque) = 0;
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
    Circular(T* circleCenter): Lagrange<T>(){
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
		inverseRadiusOfCircle = 1.0/Math::nrm2<T>(Constants::DIM - 1, radiusVector, Constants::INCX);
		Math::scal<T>(Constants::DIM - 1, inverseRadiusOfCircle, radiusVector, Constants::INCX);

		// Compute the tangential direction of the velocity to compute angular velocity
		tangentialVelocityDirection[0] = - radiusVector[1];
		tangentialVelocityDirection[1] =  radiusVector[0];

		// Compute Angular velocity using w = v * sin(theta) / r = dot_product(v, tangential_unit_vector) / r 
		carObj->_currentAngularVelocityLagrangian = Math::dot<T>(Constants::DIM - 1, tangentialVelocityDirection, Constants::INCX, carObj->getCurrentVelocityLagrangian(), Constants::INCX);
		carObj->_currentAngularVelocityLagrangian *= inverseRadiusOfCircle;
		// Copy this to initial vector too
        carObj->setInitialAngularVelocityGlobalZ(carObj->getCurrentAngularVelocityLagrangian());

		// apply angle condition on z angle


	}

	virtual void GetProfileTorqueLagrangian(const size_t& _iterationCount, Car<T>* carObj, T* externalTorque){
		*externalTorque = 0;
	}

	/**
	* Computes centrifugal force component in the Lagrangian Frame
	* \param[in] distanceFromCenter: Distance of each component of the car from center of circle
	* \param[in] velocityLagrangian: Velocity vector of the components of the car in the Lagrangian Frame
	* \param[in] mass: Mass vector containing mass of each component of the car in the array of form [CG, W1, T1, ... ]
	* \param[out] centrifugalForce: Computed centrifugal force on each component of the car
	*/
	void ComputeCentrifugalForceLagrangian(const T* distanceFromCenter, const T* velocityLagrangian, const T* mass, T* centrifugalForce) {
		// Copy the radius vector to the force to extract direction
		Math::copy<T>((Constants::DIM - 1) * Constants::VEC_DIM, distanceFromCenter, Constants::INCX, centrifugalForce, Constants::INCX);
		T velocityMagnitude;
		for (auto i = 0; i < Constants::VEC_DIM; ++i) {
			inverseRadiusOfCircle = 1.0 / Math::nrm2<T>((Constants::DIM - 1), centrifugalForce + (Constants::DIM - 1)*i, Constants::INCX);

			Math::scal<T>((Constants::DIM - 1), inverseRadiusOfCircle, centrifugalForce + (Constants::DIM - 1)*i, Constants::INCX);
			// Cross product to get tangential velocity component
			tangentialVelocityDirection[0] = *(centrifugalForce + (Constants::DIM - 1)*i + 1);
			tangentialVelocityDirection[1] = - *(centrifugalForce + (Constants::DIM - 1)*i);
			velocityMagnitude = Math::dot<T>((Constants::DIM - 1), velocityLagrangian + (Constants::DIM - 1)*i, Constants::INCX, tangentialVelocityDirection, Constants::INCX);
			// Force computation
			Math::scal<T>(Constants::DIM - 1, mass[i]* velocityMagnitude*velocityMagnitude*inverseRadiusOfCircle, centrifugalForce + (Constants::DIM - 1)*i, Constants::INCX);
		}
	}

	/**
	* Get Lagrangian forces acting due to the road profile
	*/
	virtual void GetProfileForceLagrangian(const size_t& _iterationCount, Car<T>* carObj, T* profileInducedForce, T* reactionOnTyre) {
		
		carObj->ComputeDisplacementToPointLagrangian(centerOfCircle, radiusVector);
		// compute centrifugal force on each component
        ComputeCentrifugalForceLagrangian(radiusVector, carObj->getCurrentVelocityLagrangian(), carObj->getMassComponents(), profileInducedForce);

		// Compute the reaction on tyre due to centrifugal force
		reactionOnTyre[0] = -profileInducedForce[0];
		reactionOnTyre[1] = -profileInducedForce[1];
		for (auto i = 1; i < Constants::VEC_DIM; ++i) {
			reactionOnTyre[0] -= profileInducedForce[(Constants::DIM - 1)*i];
			reactionOnTyre[1] -= profileInducedForce[(Constants::DIM - 1)*i + 1];
		}
	}
};  // Circular


/** Follow a circular road profile with radius R and center C */
template <typename T>
class Straight : public Lagrange<T> {
private:

public:
	Straight() : Lagrange<T>() {
		Name = "Straight";
	}

	virtual ~Straight() {}

	/**
	 * In case the initial conditions have to be specific to a road profile
	 * \param Car
	 */
	virtual void ApplyProfileInitialCondition(Car<T>* carObj) {}

	virtual void GetProfileTorqueLagrangian(const size_t& _iterationCount, Car<T>* carObj, T* externalTorque) {
		*externalTorque = 0;
	}

	/**
	* Get Lagrangian forces acting due to the road profile
	*/
	virtual void GetProfileForceLagrangian(const size_t& _iterationCount, Car<T>* carObj, T* profileInducedForce, T* reactionOnTyre) {
		Math::scal<T>(Constants::VEC_DIM * (Constants::DIM - 1), 0, profileInducedForce, Constants::INCX);
		Math::scal<T>((Constants::DIM - 1), 0, reactionOnTyre, Constants::INCX);
	}
};  // Straight Road



/* Eulerian profile type with fixed boundary condition */
template <typename T>
class Fixed : public Euler<T> {
private:
	// bla bla
public:
	Fixed(T &g) : Euler<T>(g) {
		Name = "Fixed";
	}
	virtual ~Fixed(){}

	virtual void GetProfileForceEulerian(const size_t& _iterationCount, Car<T>* carObj, T* profileInducedForce) {
		// TODO optimize it
		Math::scal<T>(Constants::DOF, 0, profileInducedForce, Constants::INCX);
		AddGravity(carObj, profileInducedForce);
		for (auto i = 0; i < Constants::NUM_LEGS; ++i) {
			profileInducedForce[Constants::TYRE_INDEX_EULER[i]] = //
                carObj->getkVec()[2 * i + 1] * (carObj->getCurrentDisplacementTwoTrackModel()[Constants::TYRE_INDEX_EULER[i]] - carObj->getCurrentDisplacementTwoTrackModel()[Constants::TYRE_INDEX_EULER[i] - 1])  //
              + carObj->getdVec()[2 * i + 1] * (carObj->getCurrentVelocityTwoTrackModel()[Constants::TYRE_INDEX_EULER[i]] - carObj->getCurrentVelocityTwoTrackModel()[Constants::TYRE_INDEX_EULER[i] - 1]);
		}
	}

	virtual void ApplyProfileInitialCondition(Car<T>* carObj) {		
		// The car has to start with fixed condition Velocity = 0; displacement = 0;
		for (auto i=0; i<Constants::NUM_LEGS; ++i){
			carObj->_currentDisplacementTwoTrackModel[Constants::TYRE_INDEX_EULER[i]] = 0;
		}
		for (auto i = 0; i < Constants::NUM_LEGS; ++i) {
			carObj->_currentVelocityTwoTrackModel[Constants::TYRE_INDEX_EULER[i]] = 0;
		}
		carObj->UpdateLengthsTwoTrackModel();
	}

};

/* Eulerian profile type with fixed boundary condition */
template <typename T>
class Nonfixed : public Euler<T> {
private:
	/// Bla Bla
public:
	Nonfixed(T& g) : Euler<T>(g) {
		Name = "Nonfixed";
	}
	virtual ~Nonfixed() {}
	virtual void GetProfileForceEulerian(const size_t& _iterationCount, Car<T>* carObj, T* profileInducedForce) {
		Math::scal<T>(Constants::DOF, 0, profileInducedForce, Constants::INCX);
		AddGravity(carObj, profileInducedForce);
	}
	virtual void ApplyProfileInitialCondition(Car<T>* carObj) {}
};
}  // namespace EVAA
