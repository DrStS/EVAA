// TODO: copyright header
#pragma once

// TODO: remove cstdio
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>

#include "Constants.h"


namespace EVAA {

/** Class that handles everything related to the arbitrary road trajetory */
template <class T>
class arbitraryTrajectory {
public:
    /**
     * \brief Constructor
     */
    arbitraryTrajectory() : _numIterations(0), _delta_t(0) {}
    /**
    * \brief Constructor
    */
    arbitraryTrajectory(size_t numIterations, T delta_t,
                        T amplitude_right, T amplitude_left,
                        T period_right, T period_left, 
                        T shift_right, T shift_left):
            _numIterations(numIterations), _delta_t(delta_t){

        _amplitudeRight = amplitude_right;
        _amplitudeLeft = amplitude_left;

        _frequencyRight = 2 * PI / period_right;
        _frequencyLeft = 2 * PI / period_left;

        _phaseShift_fl = 2 * PI / shift_left;
        _phaseShift_fr = 2 * PI / shift_right;
        _phaseShift_rl = 2 * PI / shift_left;
        _phaseShift_rr = 2 * PI / shift_right;


        // allocate memory
        _roadPointsX = Math::malloc<T>(numIterations);
        _roadPointsY = Math::malloc<T>(numIterations);
        _roadForcesX = Math::malloc<T>(numIterations);
        _roadForcesY = Math::malloc<T>(numIterations);
        _roadTorque = Math::malloc<T>(numIterations);
        _normedDistance = Math::malloc<T>(numIterations);
    }

    /**
    * Destructor
    */
    void ~arbitraryTrajectory() { 
        Math::free(_roadPointsX);
        Math::free(_roadPointsY);
        Math::free(_roadForcesX);
        Math::free(_roadForcesY);
        Math::free(_roadTorque);
        Math::free(_normedDistance);

    }

    
    /**
    * \brief calcuate the addition shift due to the geometry of the car
    */
    void calculateTyreShifts(long_fl, long_fr, long_rl, long_rr) {
        _phaseShift_fl += 2 * PI / long_fl;
        _phaseShift_fr += 2 * PI / long_fr;
        _phaseShift_rl -= 2 * PI / long_rl;
        _phaseShift_rr -= 2 * PI / long_rr;
    }


    /**
    * \brief interpolate all intermediate road points @Felix, help me out
    * \param numProvidedPoints number of road points from the XML
    * \param providedPointsX array of all X-coordinates of the points from the XML 
    * \param providedPointsY array of all Y-coordinates of the points from the XML
    * \param providedTimes times, at which the XML points should be reached 
    * \param initialVelocity such that in the first segment the acceleration is calculated correctly
    
    * This means, that between two XML points  [X(i), Y(i)] and [X(ì+1), Y(i+1)], there will be n=(times(í+1) - times(i)) / delta_t true road points
     */
    void interpolateRoadPoints(size_t numProvidedPoints, T* providedPointsX, T* providedPointsY, T* providedTimes, T initialVelocity) {


        // update _roadPointsX[j] and _roadPointsY[j]


    }


    /**
    * \brief calculate the travelled distance at each trajectory point 
    */
    void calculateTravelledDistance() {
        //TODO as private from the previous function
    }

    /**
    * \brief use a second order scheme to calculate all torques and forces between the trajectory points
    */
    void calculateForcesAndTorques() {
        //TODO
    }

    /**
    * \brief calculates the force acting on the tyre due to bumpy road
    * \param time at which this happens
    */
    T getVerticalRoadForcesFrontLeft(T time) {

    }

     /**
     * \brief calculates the force acting on the tyre due to bumpy road
     * \param time at which this happens
     */
    T getVerticalRoadForcesFrontRight(T time) {}


     /**
     * \brief calculates the force acting on the tyre due to bumpy road
     * \param time at which this happens
     */
    T getVerticalRoadForcesRearLeft(T time) {}


     /**
     * \brief calculates the force acting on the tyre due to bumpy road
     * \param time at which this happens
     */
    T getVerticalRoadForcesRearRight(T time) {}



     /**
     * \param iteration index of the current time iteration
     * \return vector the 2D vector to be written too [GC:XY]
     */
    void getLagrangianForces(size_t iteration, T* vector) {
        vector[0] = _roadForcesX[i];
        vector[1] = _roadForcesY[i];
    }
    /**
     * \param iteration index of the current time iteration
     * \return the Z component [GC:Z]
     */
    T getLagrangianTorque(iteration) { return _roadTorque[i]; }

    /**
     * \param iteration index of the current time iteration
     * \return vector the 2D vector to be written too [GC:XY]
     */
    void getLagrangianPosition(size_t iteration, T* vector) {
        vector[0] = _roadPointsX[i];
        vector[1] = _roadPointsY[i];
    }
 



private:
    const size_t _numIterations;    // total number of iterations

    const T _delta_t;               // time step size of the simulation

    // Eulerframe
    T _amplitudeRight;               // amplitude of the sinusoidal curve acting on the right tyres
    T _amplitudeLeft;                // amplitude of the sinusoidal curve acting on the left tyres  

    T _frequencyRight;               // frequency of the sinusoidal curve acting on the right tyres
    T _frequencyLeft;                // frequency of the sinusoidal curve acting on the left tyres

    T _phaseShift_fl;               // phase shift of the sinusoidal curve acting on the front_left tyres
    T _phaseShift_fr;               // phase shift of the sinusoidal curve acting on the front-right tyres
    T _phaseShift_rl;               // phase shift of the sinusoidal curve acting on the rear_left tyres
    T _phaseShift_rr;               // phase shift of the sinusoidal curve acting on the rear_right tyres
              

    // Lagrangian frame
    T* _roadPointsX;                // X-coordinates of the trajectory at all iterations
    T* _roadPointsY;                // Y-coordinates of the trajectory at all iterations

    T* _roadForcesX;                // X-components of the force acting on the center of gravity at all iterations
    T* _roadForcesY;                // Y-components of the force acting on the center of gravity at all iterations

    T* _roadTorque;                 // Z-components of the torque acting on the center of gravity at all iterations

    T* _normedDistance;             // absolute distance since the beginning of the road
    const double PI = 3.141592653589793238463   // TODO: add this into constants
};

}  // namespace EVAA
