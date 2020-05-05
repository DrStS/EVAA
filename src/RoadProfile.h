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
	 * Get external force acting on the car system in the Eulerian Frame
	 * \param Car
	 * \return profileInducedForce forces and torque acting on each component [GC: Z, GC(Torque):XY, W1: Z, T1: Z, ...]
	 */
	virtual void GetProfileForceEulerian(Car<T>* carObj, T* profileInducedForce) = 0;

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
	virtual void GetProfileForceLagrangian(Car<T>* carObj, T* profileInducedForce, T* reactionOnTyre) = 0;

	/**
	 * Get external torque acting on the car system in Lagrangian Frame
	 * \param Car
	 * \return externalTorque torque acting on the car system [GC: Z]
	 */
	virtual void GetProfileTorqueLagrangian(Car<T>* carObj, T* externalTorque) = 0;
	virtual ~Lagrange() {}
};


/** Follow a circular road profile with radius R and center C */
template <typename T>
class Circular : public Lagrange<T> {
private:
    /** Center of the Circle (if the motion is a Circle) */
    T* centerOfCircle = NULL;

    /** Radius of the Circle (if the motion is a Circle) */
    T radiusOfCircle = 0;

    // vectors used in computations inside the methods
    T* tangentialVelocityDirection;  // arrow from center towards object

    // distance vector between center of circle and the car =
    // Positions of points from the car vs Center of Circle, seen as 0
    T* radiusVector;  // 3 x 9
    

public:
    Circular(T* circleCenter, T circleRadius): Lagrange<T>(){
        Name = "Circular";

        // Position
		centerOfCircle = Math::malloc<T>(Constants::DIM);
        Math::copy<T>(Constants::DIM, circleCenter, 1, centerOfCircle, 1);

		radiusOfCircle = circleRadius;

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

		T nrm;
		// Compute Radius Vector based on the position of origin of euler frame
		radiusVector[0] = carObj->currentPositionLagrangian[0] - centerOfCircle[0];
		radiusVector[1] = carObj->currentPositionLagrangian[1] - centerOfCircle[1];
		nrm = 1.0/Math::nrm2<T>(Constants::DIM - 1, radiusVector, Constants::INCX);
		Math::scal<T>(Constants::DIM - 1, nrm, radiusVector, Constants::INCX);

		// Compute the tangential direction of the velocity to compute angular velocity
		tangentialVelocityDirection[0] = - radiusVector[1];
		tangentialVelocityDirection[1] =  radiusVector[0];

		// Compute Angular velocity using w = v * sin(theta) / r = dot_product(v, tangential_unit_vector) / r 
		carObj->currentAngularVelocityLagrangian = Math::dot<T>(Constants::DIM - 1, tangentialVelocityDirection, Constants::INCX, carObj->currentVelocityLagrangian, Constants::INCX);
		carObj->currentAngularVelocityLagrangian = carObj->currentAngularVelocityLagrangian * nrm;
		// Copy this to initial vector too
		carObj->initialAngularVelocityGlobal[2] = carObj->currentAngularVelocityLagrangian;
		
	}

	virtual void GetProfileTorqueLagrangian(Car<T>* carObj, T* externalTorque){
		*externalTorque = 0;
	}

	/**
	* Computes centrifugal force component in the Lagrangian Frame
	* \param distanceFromCenter: Distance of each component of the car from center of circle
	* \param velocityLagrangian: Velocity vector of the components of the car in the Lagrangian Frame
	* \param mass: Mass vector containing mass of each component of the car in the array of form [CG, W1, T1, ... ]
	* \param centrifugalForce: Computed centrifugal force on each component of the car
	*/
	void ComputeCentrifugalForceLagrangian(T* distanceFromCenter, T* velocityLagrangian, T* mass, T* centrifugalForce) {
		// Copy the radius vector to the force to extract direction
		Math::copy<T>((Constants::DIM - 1) * Constants::VEC_DIM, distanceFromCenter, Constants::INCX, centrifugalForce, Constants::INCX);
		T inverseRadius, velocityMagnitude;
		for (auto i = 0; i < Constants::VEC_DIM; ++i) {
			inverseRadius = 1.0 / Math::nrm2<T>((Constants::DIM - 1), centrifugalForce + (Constants::DIM - 1)*i, Constants::INCX);

			Math::scal<T>((Constants::DIM - 1), inverseRadius, centrifugalForce + (Constants::DIM - 1)*i, Constants::INCX);
			// Cross product to get tangential velocity component
			tangentialVelocityDirection[0] = *(centrifugalForce + (Constants::DIM - 1)*i + 1);
			tangentialVelocityDirection[1] = - *(centrifugalForce + (Constants::DIM - 1)*i);
			velocityMagnitude = Math::dot<T>((Constants::DIM - 1), velocityLagrangian + (Constants::DIM - 1)*i, Constants::INCX, tangentialVelocityDirection, Constants::INCX);
			// Force computation
			Math::scal<T>(Constants::DIM - 1, mass[i]* velocityMagnitude*velocityMagnitude*inverseRadius, centrifugalForce + (Constants::DIM - 1)*i, Constants::INCX);
		}
	}

	/**
	* Get Lagrangian forces acting due to the road profile
	*/
	virtual void GetProfileForceLagrangian(Car<T>* carObj, T* profileInducedForce, T* reactionOnTyre) {
		
		carObj->ComputeDisplacementToPointLagrangian(centerOfCircle, radiusVector);
		// compute centrifugal force on each component
		ComputeCentrifugalForceLagrangian(radiusVector, carObj->currentVelocityLagrangian, carObj->massComponents, profileInducedForce);

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

	virtual void GetProfileTorqueLagrangian(Car<T>* carObj, T* externalTorque) {
		*externalTorque = 0;
	}

	/**
	* Get Lagrangian forces acting due to the road profile
	*/
	virtual void GetProfileForceLagrangian(Car<T>* carObj, T* profileInducedForce, T* reactionOnTyre) {	
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

	virtual void GetProfileForceEulerian(Car<T>* carObj, T* profileInducedForce) {
		Math::scal<T>(Constants::DOF, 0, profileInducedForce, Constants::INCX);
		for (auto i = 0; i < Constants::NUM_LEGS; ++i) {
			profileInducedForce[Constants::TYRE_INDEX_EULER[i]] = carObj->kVec[2 * i + 1] * (carObj->currentDisplacementTwoTrackModel[Constants::TYRE_INDEX_EULER[i]] - carObj->currentDisplacementTwoTrackModel[Constants::TYRE_INDEX_EULER[i] - 1]) 
				+ carObj->dVec[2 * i + 1] * (carObj->currentVelocityTwoTrackModel[Constants::TYRE_INDEX_EULER[i]] - carObj->currentVelocityTwoTrackModel[Constants::TYRE_INDEX_EULER[i] - 1]);
		}
	}

	virtual void ApplyProfileInitialCondition(Car<T>* carObj) {
		
		// The car has to start with fixed condition Velocity = 0; displacement = 0;
		for (auto i=0; i<Constants::NUM_LEGS; ++i){
			carObj->currentDisplacementTwoTrackModel[Constants::TYRE_INDEX_EULER[i]] = 0;
		}
		for (auto i = 0; i < Constants::NUM_LEGS; ++i) {
			carObj->currentVelocityTwoTrackModel[Constants::TYRE_INDEX_EULER[i]] = 0;
		}
		carObj->updateLengthsTwoTrackModel();
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
	virtual void GetProfileForceEulerian(Car<T>* carObj, T* profileInducedForce) { Math::scal<T>(Constants::DOF, 0, profileInducedForce, Constants::INCX); }
	virtual void ApplyProfileInitialCondition(Car<T>* carObj) {}
};
}  // namespace EVAA
