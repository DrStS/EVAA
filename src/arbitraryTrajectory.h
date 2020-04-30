// TODO: copyright header
#pragma once

// TODO: remove cstdio
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <math.h>

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
    arbitraryTrajectory(size_t numIterations, T delta_t, T amplitude_right, T amplitude_left,
                        T period_right, T period_left, T shift_right, T shift_left) :
        _numIterations(numIterations), _delta_t(delta_t) {

        if (period_right == 0) {
            _frequencyRight = 0;
        }
        else {
            _frequencyRight = 2 * PI / period_right;
        }

        if (period_left == 0) {
            _frequencyLeft = 0;
        }
        else {
            _frequencyLeft = 2 * PI / period_left;
        }

        _amplitudeRight = amplitude_right;
        _amplitudeLeft = amplitude_left;

        _phaseShift_fl = 2 * PI / shift_left;
        _phaseShift_fr = 2 * PI / shift_right;
        _phaseShift_rl = 2 * PI / shift_left;
        _phaseShift_rr = 2 * PI / shift_right;

        // allocate memory
        _roadPointsX = Math::malloc<T>(numIterations);
        _roadPointsY = Math::malloc<T>(numIterations);
        _roadForcesX = Math::malloc<T>(numIterations);
        _roadForcesY = Math::malloc<T>(numIterations);
        _roadAngles = Math::malloc<T>(numIterations);
        _roadTorque = Math::malloc<T>(numIterations);
        _normedDistance = Math::malloc<T>(numIterations);
        _tyreForces_fl = Math::malloc<T>(numIterations);
        _tyreForces_fr = Math::malloc<T>(numIterations);
        _tyreForces_rl = Math::malloc<T>(numIterations);
        _tyreForces_rr = Math::malloc<T>(numIterations);
    }

    /**
     * Destructor
     */
    ~arbitraryTrajectory() {
        Math::free(_roadPointsX);
        Math::free(_roadPointsY);
        Math::free(_roadForcesX);
        Math::free(_roadForcesY);
        Math::free(_roadAngles);
        Math::free(_roadTorque);
        Math::free(_normedDistance);
        Math::free(_tyreForces_fl);
        Math::free(_tyreForces_fr);
        Math::free(_tyreForces_rl);
        Math::free(_tyreForces_rr);
    }

    /**
     * \brief calcuate the addition shift due to the geometry of the car
     */
    void calculateTyreShifts(T long_fl, T long_fr, T long_rl, T long_rr) {
        _phaseShift_fl += long_fl;
        _phaseShift_fr += long_fr;
        _phaseShift_rl -= long_rl;
        _phaseShift_rr -= long_rr;
    }

    /**
    * \brief interpolate all intermediate road points @Felix, help me out
    * \param numProvidedPoints number of road points from the XML
    * \param providedPointsX array of all X-coordinates of the points from the XML
    * \param providedPointsY array of all Y-coordinates of the points from the XML
    * \param providedTimes times, at which the XML points should be reached
    * \param initialVelocity such that in the first segment the acceleration is calculated correctly

    * This means, that between two XML points  [X(i), Y(i)] and [X(ì+1), Y(i+1)], there will be
    n=(times(í+1) - times(i)) / delta_t true road points
     */
    void interpolateRoadPoints(size_t numProvidedPoints, T* providedPointsX, T* providedPointsY,
                               T* providedTimes, T initialVelocity) {
        // update _roadPointsX[j] and _roadPointsY[j]
    }

    /**
     * \brief calculate the travelled distance at each trajectory point
     */
    void calculateTravelledDistance() {
        _normedDistance[0] = 0;
        for (int i = 1; i < _numIterations; ++i) {
            _normedDistance[i] = _normedDistance[i - 1] + 
                sqrt((_roadPointsX[i] - _roadPointsX[i - 1]) * (_roadPointsX[i] - _roadPointsX[i - 1]) +
                     (_roadPointsY[i] - _roadPointsY[i - 1]) * (_roadPointsY[i] - _roadPointsY[i - 1]));
        }
    }

    
     /**
     * \brief use a second order scheme to calculate all angles on the trajectory
     * points
     */
    void calculateAngles() {
        T invDeltaT = 1. / _delta_t;
        if (_numIterations > 3) {
            T directionX = -1.5 * _roadPointsX[0] + 2 * _roadPointsX[1] - 0.5 * _roadPointsX[2];
            T directionY = -1.5 * _roadPointsY[0] + 2 * _roadPointsY[1] - 0.5 * _roadPointsY[2];
            T directionNorm = directionX * directionX + directionY + directionY;
            if (directionNorm == 0) {
                _roadAngles[0] = 0;
            }
            else if (directionY >= 0) {
                _roadAngles[0] = std::acos(directionX / directionNorm);       
            }
            else {
                _roadAngles[0] = -std::acos(directionX / directionNorm);       
            }

            for (int i = 1; i < _numIterations - 1; ++i) {
                directionX = -0.5 * _roadPointsX[i - 1] + 0.5 * _roadPointsX[i + 1];
                directionY = -0.5 * _roadPointsY[i - 1] + 0.5 * _roadPointsY[i + 1];
                T directionNorm = directionX * directionX + directionY + directionY;
                if (directionNorm == 0) {
                    _roadAngles[i] = _roadAngles[i-1];
                }
                else if (directionY >= 0) {
                    _roadAngles[i] = std::acos(directionX / directionNorm);
                }
                else {
                    _roadAngles[i] = -std::acos(directionX / directionNorm);
                }

            }

            directionX = -1.5 * _roadPointsX[_numIterations - 1] +
                         2 * _roadPointsX[_numIterations - 2] -
                         0.5 * _roadPointsX[_numIterations - 3];
            directionY = -1.5 * _roadPointsY[_numIterations - 1] +
                         2 * _roadPointsY[_numIterations - 2] -
                         0.5 * _roadPointsY[_numIterations - 3];
            directionNorm = directionX * directionX + directionY + directionY;
            if (directionNorm == 0) {
                _roadAngles[_numIterations - 1] = _roadAngles[_numIterations - 2];
            }
            else if (directionY >= 0) {
                _roadAngles[_numIterations - 1] = std::acos(directionX / directionNorm);
            }
            else {
                _roadAngles[_numIterations - 1] = -std::acos(directionX / directionNorm);
            }
        }
        else {
            for (int i = 0; i < _numIterations; ++i) {
                _roadAngles[i] = 0;
            }
        }
    }


    /**
     * \brief use a second order scheme to calculate all torques and forces between the trajectory
     * points
     */
    void calculateForcesAndTorques() {
        T invDeltaT = 1. / (_delta_t * _delta_t);
        if (_numIterations > 4) {
            _roadForcesX[0] = (2 * _roadPointsX[0] + 5 * _roadPointsX[1] + 4 * _roadPointsX[2] + _roadPointsX[3]) * invDeltaT;
            _roadForcesY[0] = (2 * _roadPointsY[0] + 5 * _roadPointsY[1] + 4 * _roadPointsY[2] + _roadPointsY[3]) * invDeltaT;
            _roadTorque[0] =  (2 * _roadAngles[0] + 5 * _roadAngles[1] + 4 * _roadAngles[2] + _roadAngles[3]) * invDeltaT;

            for (int i = 1; i < _numIterations - 1; ++i) {
                _roadForcesX[i] =  (_roadPointsX[i - 1] - 2 * _roadPointsX[i] + _roadPointsX[i + 1]) * invDeltaT;
                _roadForcesY[i] =  (_roadPointsY[i - 1] - 2 * _roadPointsY[i] + _roadPointsY[i + 1]) * invDeltaT;
                _roadTorque[i] = (_roadAngles[i - 1] - 2 * _roadAngles[i] + _roadAngles[i + 1]) * invDeltaT;
            }

            _roadForcesX[_numIterations - 1] =
                (2 * _roadPointsX[_numIterations - 1] + 5 * _roadPointsX[_numIterations - 2] +
                 4 * _roadPointsX[_numIterations - 3] + _roadPointsX[_numIterations - 4]) * invDeltaT;
            _roadForcesY[_numIterations - 1] =
                (2 * _roadPointsY[_numIterations - 1] + 5 * _roadPointsY[_numIterations - 2] +
                 4 * _roadPointsY[_numIterations - 3] + _roadPointsY[_numIterations - 4]) * invDeltaT;
            _roadTorque[_numIterations - 1] =
                (2 * _roadAngles[_numIterations - 1] + 5 * _roadAngles[_numIterations - 2] +
                 4 * _roadAngles[_numIterations - 3] + _roadAngles[_numIterations - 3]) * invDeltaT;
        }
        else {
            for (int i = 0; i < _numIterations; ++i) {
                _roadForcesX[i] = 0;
                _roadForcesY[i] = 0;
                _roadTorque[i] = 0;
            }
        }
    }
 
  
       /**
     * \brief use a second order scheme to calculate the forces that act on the tyre on the trajectory
     * points
     */
    void calculateVerticalForces() {
        T invDeltaT = 1. / (_delta_t * _delta_t);
        if (_numIterations > 4) {
            T pos_fl_prev = _amplitudeLeft * std::sin(_frequencyLeft * (_normedDistance[0] + _phaseShift_fl));
            T pos_fr_prev = _amplitudeRight * std::sin(_frequencyRight * (_normedDistance[0] + _phaseShift_fr));
            T pos_rl_prev = _amplitudeLeft * std::sin(_frequencyLeft * (_normedDistance[0] + _phaseShift_rl));
            T pos_rr_prev = _amplitudeRight * std::sin(_frequencyRight * (_normedDistance[0] + _phaseShift_rr));

            T pos_fl = _amplitudeLeft * std::sin(_frequencyLeft * (_normedDistance[1] + _phaseShift_fl));
            T pos_fr = _amplitudeRight * std::sin(_frequencyRight * (_normedDistance[1] + _phaseShift_fr));
            T pos_rl = _amplitudeLeft * std::sin(_frequencyLeft * (_normedDistance[1] + _phaseShift_rl));
            T pos_rr = _amplitudeRight * std::sin(_frequencyRight * (_normedDistance[1] + _phaseShift_rr));

            T pos_fl_next = _amplitudeLeft * std::sin(_frequencyLeft * (_normedDistance[2] + _phaseShift_fl));
            T pos_fr_next = _amplitudeRight * std::sin(_frequencyRight * (_normedDistance[2] + _phaseShift_fr));
            T pos_rl_next = _amplitudeLeft * std::sin(_frequencyLeft * (_normedDistance[2] + _phaseShift_rl));
            T pos_rr_next = _amplitudeRight * std::sin(_frequencyRight * (_normedDistance[2] + _phaseShift_rr));

            T pos_fl_nnext = _amplitudeLeft * std::sin(_frequencyLeft * (_normedDistance[3] + _phaseShift_fl));
            T pos_fr_nnext = _amplitudeRight * std::sin(_frequencyRight * (_normedDistance[3] + _phaseShift_fr));
            T pos_rl_nnext = _amplitudeLeft * std::sin(_frequencyLeft * (_normedDistance[3] + _phaseShift_rl));
            T pos_rr_nnext = _amplitudeRight * std::sin(_frequencyRight * (_normedDistance[3] + _phaseShift_rr));

            _tyreForces_fl[0] = (2 * pos_fl_prev + 5 * pos_fl + 4 * pos_fl_next + pos_fl_nnext) * invDeltaT;
            _tyreForces_fr[0] = (2 * pos_fr_prev + 5 * pos_fr + 4 * pos_fr_next + pos_fr_nnext) * invDeltaT;
            _tyreForces_rl[0] = (2 * pos_rl_prev + 5 * pos_rl + 4 * pos_rl_next + pos_rl_nnext) * invDeltaT;
            _tyreForces_rr[0] = (2 * pos_rr_prev + 5 * pos_rr + 4 * pos_rr_next + pos_rr_nnext) * invDeltaT;

            for (int i = 1; i < _numIterations - 1; ++i) {
                pos_fl_prev = pos_fl;
                pos_fl = pos_fl_next;
                pos_fl_next = _amplitudeLeft * std::sin(_frequencyLeft * (_normedDistance[i + 1] + _phaseShift_fl));
                pos_fr_prev = pos_fr;
                pos_fr = pos_fr_next;
                pos_fr_next = _amplitudeRight * std::sin(_frequencyRight * (_normedDistance[i + 1] + _phaseShift_fr));
                pos_rl_prev = pos_rl;
                pos_rl = pos_rl_next;
                pos_rl_next = _amplitudeLeft * std::sin(_frequencyLeft * (_normedDistance[i + 1] + _phaseShift_rl));
                pos_rr_prev = pos_rr;
                pos_rr = pos_rr_next;
                pos_rr_next = _amplitudeRight * std::sin(_frequencyRight * (_normedDistance[i + 1] + _phaseShift_rr));

                _tyreForces_fl[i] = (pos_fl_prev - 2 * pos_fl + pos_fl_next) * invDeltaT;
                _tyreForces_fr[i] = (pos_fr_prev - 2 * pos_fr + pos_fr_next) * invDeltaT;
                _tyreForces_rl[i] = (pos_rl_prev - 2 * pos_rl + pos_rl_next) * invDeltaT;
                _tyreForces_rr[i] = (pos_rr_prev - 2 * pos_rr + pos_rr_next) * invDeltaT;
            }
            T pos_fl_nnext = _amplitudeLeft * std::sin(_frequencyLeft * (_normedDistance[_numIterations - 1] + _phaseShift_fl));
            T pos_fr_nnext = _amplitudeRight * std::sin(_frequencyRight * (_normedDistance[_numIterations - 1] + _phaseShift_fr));
            T pos_rl_nnext = _amplitudeLeft * std::sin(_frequencyLeft * (_normedDistance[_numIterations - 1] + _phaseShift_rl));
            T pos_rr_nnext = _amplitudeRight * std::sin(_frequencyRight * (_normedDistance[_numIterations - 1] + _phaseShift_rr));

            _tyreForces_fl[_numIterations - 1] = (2 * pos_fl_nnext + 5 * pos_fl_next + 4 * pos_fl + pos_fl_prev) * invDeltaT;
            _tyreForces_fr[_numIterations - 1] = (2 * pos_fr_nnext + 5 * pos_fr_next + 4 * pos_fr + pos_fr_prev) * invDeltaT;
            _tyreForces_rl[_numIterations - 1] = (2 * pos_rl_nnext + 5 * pos_rl_next + 4 * pos_rl + pos_rl_prev) * invDeltaT;
            _tyreForces_rr[_numIterations - 1] = (2 * pos_rr_nnext + 5 * pos_rr_next + 4 * pos_rr + pos_rr_prev) * invDeltaT;
        }
        else {
            for (int i = 0; i < _numIterations; ++i) {
                _tyreForces_fl[i] = 0;
                _tyreForces_fr[i] = 0;
                _tyreForces_rl[i] = 0;
                _tyreForces_rr[i] = 0;
            }
        }
    }
 

    /**
     * \brief calculates the force acting on the tyre due to bumpy road
     * \param time at which this happens
     */
    T getVerticalRoadForcesFrontLeft(size_t iteration) { return _tyreForces_fl[iteration]; }

    /**
     * \brief calculates the force acting on the tyre due to bumpy road
     * \param time at which this happens
     */
    T getVerticalRoadForcesFrontRight(size_t iteration) { return _tyreForces_rl[iteration]; }

    /**
     * \brief calculates the force acting on the tyre due to bumpy road
     * \param time at which this happens
     */
    T getVerticalRoadForcesRearLeft(size_t iteration) { return _tyreForces_fr[iteration]; }

    /**
     * \brief calculates the force acting on the tyre due to bumpy road
     * \param time at which this happens
     */
    T getVerticalRoadForcesRearRight(size_t iteration) { return _tyreForces_rr[iteration]; }

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
    T getLagrangianTorque(size_t iteration) { return _roadTorque[i]; }

    /**
     * \param iteration index of the current time iteration
     * \return vector the 2D vector to be written too [GC:XY]
     */
    void getLagrangianPosition(size_t iteration, T* vector) {
        vector[0] = _roadPointsX[i];
        vector[1] = _roadPointsY[i];
    }

private:
    const size_t _numIterations;  // total number of iterations

    const T _delta_t;  // time step size of the simulation

    // Eulerframe
    T _amplitudeRight;  // amplitude of the sinusoidal curve acting on the right tyres
    T _amplitudeLeft;   // amplitude of the sinusoidal curve acting on the left tyres

    T _frequencyRight;  // frequency of the sinusoidal curve acting on the right tyres
    T _frequencyLeft;   // frequency of the sinusoidal curve acting on the left tyres

    T _phaseShift_fl;  // phase shift of the sinusoidal curve acting on the front_left tyres
    T _phaseShift_fr;  // phase shift of the sinusoidal curve acting on the front-right tyres
    T _phaseShift_rl;  // phase shift of the sinusoidal curve acting on the rear_left tyres
    T _phaseShift_rr;  // phase shift of the sinusoidal curve acting on the rear_right tyres

    T* _tyreForces_fl; // force that acts on the tyre in Z direction
    T* _tyreForces_fr; // force that acts on the tyre in Z direction
    T* _tyreForces_rl; // force that acts on the tyre in Z direction
    T* _tyreForces_rr; // force that acts on the tyre in Z direction    

    // Lagrangian frame
    T* _roadPointsX;  // X-coordinates of the trajectory at all iterations
    T* _roadPointsY;  // Y-coordinates of the trajectory at all iterations

    T* _roadForcesX;  // X-components of the force acting on the center of gravity at all iterations
    T* _roadForcesY;  // Y-components of the force acting on the center of gravity at all iterations

    T* _roadAngles;   // contains all angles

    T* _roadTorque;  // Z-components of the torque acting on the center of gravity at all iterations

    T* _normedDistance;                        // absolute distance since the beginning of the road



    const double PI = 3.141592653589793238463; // TODO: add this into constants
};

}  // namespace EVAA
