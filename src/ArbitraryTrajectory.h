// TODO: copyright header
#pragma once

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

#include "Constants.h"

namespace EVAA {

/** Class that handles everything related to the arbitrary road trajetory */
template <class T>
class ArbitraryTrajectory {
public:
    /**
     * \brief Constructor
     */
    ArbitraryTrajectory() : _numIterations(0), _delta_t(0) {}
    /**
     * \brief Constructor
     */
    ArbitraryTrajectory(size_t numIterations, T delta_t, T amplitude_right, T amplitude_left,
                        T period_right, T period_left, T shift_right, T shift_left) :
        _numIterations(numIterations), _delta_t(delta_t) {
        if (period_right == 0) {
            _frequencyRight = 0;
        }
        else {
            _frequencyRight = 2 * Constants::PI / period_right;
        }

        if (period_left == 0) {
            _frequencyLeft = 0;
        }
        else {
            _frequencyLeft = 2 * Constants::PI / period_left;
        }

        _amplitudeRight = amplitude_right;
        _amplitudeLeft = amplitude_left;

        _phaseShift_fl = 2 * Constants::PI / shift_left;
        _phaseShift_fr = 2 * Constants::PI / shift_right;
        _phaseShift_rl = 2 * Constants::PI / shift_left;
        _phaseShift_rr = 2 * Constants::PI / shift_right;

        // allocate memory
        _roadPointsX = Math::malloc<T>(numIterations);
        _roadPointsY = Math::malloc<T>(numIterations);
        _roadAccelerationX = Math::malloc<T>(numIterations);
        _roadAccelerationY = Math::malloc<T>(numIterations);
        _roadAngles = Math::malloc<T>(numIterations);
        _roadAngularAcceleration = Math::malloc<T>(numIterations);
        _normedDistance = Math::malloc<T>(numIterations);
        _tyreAccelerationsX_fl = Math::malloc<T>(numIterations);
        _tyreAccelerationsX_fr = Math::malloc<T>(numIterations);
        _tyreAccelerationsX_rl = Math::malloc<T>(numIterations);
        _tyreAccelerationsX_rr = Math::malloc<T>(numIterations);
        _tyreAccelerationsY_fl = Math::malloc<T>(numIterations);
        _tyreAccelerationsY_fr = Math::malloc<T>(numIterations);
        _tyreAccelerationsY_rl = Math::malloc<T>(numIterations);
        _tyreAccelerationsY_rr = Math::malloc<T>(numIterations);
    }

    /**
     * Destructor
     */
    ~ArbitraryTrajectory() {
        Math::free<T>(_roadPointsX);
        Math::free<T>(_roadPointsY);
        Math::free<T>(_roadAccelerationX);
        Math::free<T>(_roadAccelerationY);
        Math::free<T>(_roadAngles);
        Math::free<T>(_roadAngularAcceleration);
        Math::free<T>(_normedDistance);
        Math::free<T>(_tyreAccelerationsX_fl);
        Math::free<T>(_tyreAccelerationsX_fr);
        Math::free<T>(_tyreAccelerationsX_rl);
        Math::free<T>(_tyreAccelerationsX_rr);
        Math::free<T>(_tyreAccelerationsY_fl);
        Math::free<T>(_tyreAccelerationsY_fr);
        Math::free<T>(_tyreAccelerationsY_rl);
        Math::free<T>(_tyreAccelerationsY_rr);
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

    * This means, that between two XML points  [X(i), Y(i)] and [X(�+1), Y(i+1)], there will be
    n=(times(�+1) - times(i)) / delta_t true road points
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
            _normedDistance[i] =
                _normedDistance[i - 1] + sqrt((_roadPointsX[i] - _roadPointsX[i - 1]) *
                                                  (_roadPointsX[i] - _roadPointsX[i - 1]) +
                                              (_roadPointsY[i] - _roadPointsY[i - 1]) *
                                                  (_roadPointsY[i] - _roadPointsY[i - 1]));
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
                    _roadAngles[i] = _roadAngles[i - 1];
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
     * \brief use a second order scheme to calculate all AngularAccelerations and Accelerations
     * between the trajectory points
     */
    void calculateAccelerationsCenterOfGravity() {
        calculateAccelerations(_roadAccelerationX, _roadPointsX);
        calculateAccelerations(_roadAccelerationX, _roadPointsY);
        calculateAccelerations(_roadAngularAcceleration, _roadAngles);
    }

    /**
     * \brief use a second order scheme to calculate all AngularAccelerations and Accelerations
     * between the trajectory points for the legs
     */
    void calculateAccelerationsLegs(T& l_long_fl, T& l_long_fr, T& l_long_rl, T& l_long_rr,
                                    T& l_lat_fl, T& l_lat_fr, T& l_lat_rl, T& l_lat_rr) {
        // get true local coordinates
        T localXcoordinates[Constants::NUM_LEGS] = {l_long_fl, l_long_fr, -l_long_rl, -l_long_rr};
        T localYcoordinates[Constants::NUM_LEGS] = {l_lat_fl, -l_lat_fr, l_lat_rl, -l_lat_rr};

        // get all leg positions
        T* legPointsX_fl = Math::malloc<T>(_numIterations);
        T* legPointsX_fr = Math::malloc<T>(_numIterations);
        T* legPointsX_rl = Math::malloc<T>(_numIterations);
        T* legPointsX_rr = Math::malloc<T>(_numIterations);

        T* legPointsY_fl = Math::malloc<T>(_numIterations);
        T* legPointsY_fr = Math::malloc<T>(_numIterations);
        T* legPointsY_rl = Math::malloc<T>(_numIterations);
        T* legPointsY_rr = Math::malloc<T>(_numIterations);

        for (int i = 0; i < _numIterations; ++i) {
            T c = std::cos(_roadAngles[i]);
            T s = std::sin(_roadAngles[i]);

            legPointsX_fl[i] =
                _roadPointsX[i] + localXcoordinates[0] * c - localYcoordinates[0] * s;
            legPointsY_fl[i] =
                _roadPointsY[i] + localXcoordinates[0] * s + localYcoordinates[0] * c;

            legPointsX_fr[i] =
                _roadPointsX[i] + localXcoordinates[1] * c - localYcoordinates[1] * s;
            legPointsY_fr[i] =
                _roadPointsY[i] + localXcoordinates[1] * s + localYcoordinates[1] * c;

            legPointsX_rl[i] =
                _roadPointsX[i] + localXcoordinates[2] * c - localYcoordinates[2] * s;
            legPointsY_rl[i] =
                _roadPointsY[i] + localXcoordinates[2] * s + localYcoordinates[2] * c;

            legPointsX_rr[i] =
                _roadPointsX[i] + localXcoordinates[3] * c - localYcoordinates[3] * s;
            legPointsY_rr[i] =
                _roadPointsY[i] + localXcoordinates[3] * s + localYcoordinates[3] * c;
        }

        // calculate tyre forces
        calculateAccelerations(_tyreAccelerationsX_fl, legPointsX_fl);
        calculateAccelerations(_tyreAccelerationsY_fl, legPointsY_fl);

        calculateAccelerations(_tyreAccelerationsX_fr, legPointsX_fr);
        calculateAccelerations(_tyreAccelerationsY_fr, legPointsY_fr);

        calculateAccelerations(_tyreAccelerationsX_rl, legPointsX_rl);
        calculateAccelerations(_tyreAccelerationsY_rl, legPointsY_rl);

        calculateAccelerations(_tyreAccelerationsX_rr, legPointsX_rr);
        calculateAccelerations(_tyreAccelerationsY_rr, legPointsY_rr);

        // delete leg positions
        Math::free<T>(legPointsX_fl);
        Math::free<T>(legPointsX_fr);
        Math::free<T>(legPointsX_rl);
        Math::free<T>(legPointsX_rr);

        Math::free<T>(legPointsY_fl);
        Math::free<T>(legPointsY_fr);
        Math::free<T>(legPointsY_rl);
        Math::free<T>(legPointsY_rr);
    }

    /**
     * \brief use a second order scheme to calculate the Accelerations that act on the tyre on the
     * trajectory points
     */
    void calculateVerticalAccelerations() {
        T invDeltaT = 1. / (_delta_t * _delta_t);
        if (_numIterations > 4) {
            T pos_fl_prev =
                _amplitudeLeft * std::sin(_frequencyLeft * (_normedDistance[0] + _phaseShift_fl));
            T pos_fr_prev =
                _amplitudeRight * std::sin(_frequencyRight * (_normedDistance[0] + _phaseShift_fr));
            T pos_rl_prev =
                _amplitudeLeft * std::sin(_frequencyLeft * (_normedDistance[0] + _phaseShift_rl));
            T pos_rr_prev =
                _amplitudeRight * std::sin(_frequencyRight * (_normedDistance[0] + _phaseShift_rr));

            T pos_fl =
                _amplitudeLeft * std::sin(_frequencyLeft * (_normedDistance[1] + _phaseShift_fl));
            T pos_fr =
                _amplitudeRight * std::sin(_frequencyRight * (_normedDistance[1] + _phaseShift_fr));
            T pos_rl =
                _amplitudeLeft * std::sin(_frequencyLeft * (_normedDistance[1] + _phaseShift_rl));
            T pos_rr =
                _amplitudeRight * std::sin(_frequencyRight * (_normedDistance[1] + _phaseShift_rr));

            T pos_fl_next =
                _amplitudeLeft * std::sin(_frequencyLeft * (_normedDistance[2] + _phaseShift_fl));
            T pos_fr_next =
                _amplitudeRight * std::sin(_frequencyRight * (_normedDistance[2] + _phaseShift_fr));
            T pos_rl_next =
                _amplitudeLeft * std::sin(_frequencyLeft * (_normedDistance[2] + _phaseShift_rl));
            T pos_rr_next =
                _amplitudeRight * std::sin(_frequencyRight * (_normedDistance[2] + _phaseShift_rr));

            T pos_fl_nnext =
                _amplitudeLeft * std::sin(_frequencyLeft * (_normedDistance[3] + _phaseShift_fl));
            T pos_fr_nnext =
                _amplitudeRight * std::sin(_frequencyRight * (_normedDistance[3] + _phaseShift_fr));
            T pos_rl_nnext =
                _amplitudeLeft * std::sin(_frequencyLeft * (_normedDistance[3] + _phaseShift_rl));
            T pos_rr_nnext =
                _amplitudeRight * std::sin(_frequencyRight * (_normedDistance[3] + _phaseShift_rr));

            _tyreAccelerations_fl[0] =
                (2 * pos_fl_prev + 5 * pos_fl + 4 * pos_fl_next + pos_fl_nnext) * invDeltaT;
            _tyreAccelerations_fr[0] =
                (2 * pos_fr_prev + 5 * pos_fr + 4 * pos_fr_next + pos_fr_nnext) * invDeltaT;
            _tyreAccelerations_rl[0] =
                (2 * pos_rl_prev + 5 * pos_rl + 4 * pos_rl_next + pos_rl_nnext) * invDeltaT;
            _tyreAccelerations_rr[0] =
                (2 * pos_rr_prev + 5 * pos_rr + 4 * pos_rr_next + pos_rr_nnext) * invDeltaT;

            for (int i = 1; i < _numIterations - 1; ++i) {
                pos_fl_prev = pos_fl;
                pos_fl = pos_fl_next;
                pos_fl_next = _amplitudeLeft *
                              std::sin(_frequencyLeft * (_normedDistance[i + 1] + _phaseShift_fl));
                pos_fr_prev = pos_fr;
                pos_fr = pos_fr_next;
                pos_fr_next = _amplitudeRight *
                              std::sin(_frequencyRight * (_normedDistance[i + 1] + _phaseShift_fr));
                pos_rl_prev = pos_rl;
                pos_rl = pos_rl_next;
                pos_rl_next = _amplitudeLeft *
                              std::sin(_frequencyLeft * (_normedDistance[i + 1] + _phaseShift_rl));
                pos_rr_prev = pos_rr;
                pos_rr = pos_rr_next;
                pos_rr_next = _amplitudeRight *
                              std::sin(_frequencyRight * (_normedDistance[i + 1] + _phaseShift_rr));

                _tyreAccelerations_fl[i] = (pos_fl_prev - 2 * pos_fl + pos_fl_next) * invDeltaT;
                _tyreAccelerations_fr[i] = (pos_fr_prev - 2 * pos_fr + pos_fr_next) * invDeltaT;
                _tyreAccelerations_rl[i] = (pos_rl_prev - 2 * pos_rl + pos_rl_next) * invDeltaT;
                _tyreAccelerations_rr[i] = (pos_rr_prev - 2 * pos_rr + pos_rr_next) * invDeltaT;
            }
            T pos_fl_nnext =
                _amplitudeLeft *
                std::sin(_frequencyLeft * (_normedDistance[_numIterations - 1] + _phaseShift_fl));
            T pos_fr_nnext =
                _amplitudeRight *
                std::sin(_frequencyRight * (_normedDistance[_numIterations - 1] + _phaseShift_fr));
            T pos_rl_nnext =
                _amplitudeLeft *
                std::sin(_frequencyLeft * (_normedDistance[_numIterations - 1] + _phaseShift_rl));
            T pos_rr_nnext =
                _amplitudeRight *
                std::sin(_frequencyRight * (_normedDistance[_numIterations - 1] + _phaseShift_rr));

            _tyreAccelerations_fl[_numIterations - 1] =
                (2 * pos_fl_nnext + 5 * pos_fl_next + 4 * pos_fl + pos_fl_prev) * invDeltaT;
            _tyreAccelerations_fr[_numIterations - 1] =
                (2 * pos_fr_nnext + 5 * pos_fr_next + 4 * pos_fr + pos_fr_prev) * invDeltaT;
            _tyreAccelerations_rl[_numIterations - 1] =
                (2 * pos_rl_nnext + 5 * pos_rl_next + 4 * pos_rl + pos_rl_prev) * invDeltaT;
            _tyreAccelerations_rr[_numIterations - 1] =
                (2 * pos_rr_nnext + 5 * pos_rr_next + 4 * pos_rr + pos_rr_prev) * invDeltaT;
        }
        else {
            for (int i = 0; i < _numIterations; ++i) {
                _tyreAccelerations_fl[i] = 0;
                _tyreAccelerations_fr[i] = 0;
                _tyreAccelerations_rl[i] = 0;
                _tyreAccelerations_rr[i] = 0;
            }
        }
    }

    /**
     * \brief calculates the Acceleration acting on the tyre due to bumpy road
     * \param time at which this happens
     */
    T getVerticalRoadAccelerationsFrontLeft(size_t iteration) {
        return _tyreAccelerations_fl[iteration];
    }

    /**
     * \brief calculates the Acceleration acting on the tyre due to bumpy road
     * \param time at which this happens
     */
    T getVerticalRoadAccelerationsFrontRight(size_t iteration) {
        return _tyreAccelerations_rl[iteration];
    }

    /**
     * \brief calculates the Acceleration acting on the tyre due to bumpy road
     * \param time at which this happens
     */
    T getVerticalRoadAccelerationsRearLeft(size_t iteration) {
        return _tyreAccelerations_fr[iteration];
    }

    /**
     * \brief calculates the Acceleration acting on the tyre due to bumpy road
     * \param time at which this happens
     */
    T getVerticalRoadAccelerationsRearRight(size_t iteration) {
        return _tyreAccelerations_rr[iteration];
    }

    /**
     * \param iteration index of the current time iteration
     * \return vector the 2D vector to be written too [GC:XY]
     */
    void getLagrangianAccelerations(size_t iteration, T* vector) {
        vector[0] = _roadAccelerationX[i];
        vector[1] = _roadAccelerationY[i];
    }
    /**
     * \param iteration index of the current time iteration
     * \return the Z component [GC:Z]
     */
    T getLagrangianAngularAcceleration(size_t iteration) { return _roadAngularAcceleration[i]; }

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

    T* _tyreAccelerationsX_fl;  // Acceleration that acts on the tyre in Z direction
    T* _tyreAccelerationsX_fr;  // Acceleration that acts on the tyre in Z direction
    T* _tyreAccelerationsX_rl;  // Acceleration that acts on the tyre in Z direction
    T* _tyreAccelerationsX_rr;  // Acceleration that acts on the tyre in Z direction

    T* _tyreAccelerationsY_fl;  // Acceleration that acts on the tyre in Z direction
    T* _tyreAccelerationsY_fr;  // Acceleration that acts on the tyre in Z direction
    T* _tyreAccelerationsY_rl;  // Acceleration that acts on the tyre in Z direction
    T* _tyreAccelerationsY_rr;  // Acceleration that acts on the tyre in Z direction

    // Lagrangian frame
    T* _roadPointsX;  // X-coordinates of the trajectory at all iterations
    T* _roadPointsY;  // Y-coordinates of the trajectory at all iterations

    T* _roadAccelerationX;  // X-components of the Acceleration acting on the center of gravity at
                            // all iterations
    T* _roadAccelerationY;  // Y-components of the Acceleration acting on the center of gravity at
                            // all iterations

    T* _roadAngles;  // contains all angles

    T* _roadAngularAcceleration;  // Z-components of the AngularAcceleration acting on the center of
                                  // gravity at all iterations

    T* _normedDistance;  // absolute distance since the beginning of the road

    /**
     * \brief use a second order scheme to calculate all Accelerations
     * between the trajectory points
     */
    void calculateAccelerations(T* acceleration, T* points) {
        T invDeltaT = 1. / (_delta_t * _delta_t);
        if (_numIterations > 4) {
            acceleration[0] =
                (2 * points[0] + 5 * points[1] + 4 * points[2] + points[3]) * invDeltaT;

            for (int i = 1; i < _numIterations - 1; ++i) {
                acceleration[i] = (points[i - 1] - 2 * points[i] + points[i + 1]) * invDeltaT;
            }

            acceleration[_numIterations - 1] =
                (2 * points[_numIterations - 1] + 5 * points[_numIterations - 2] +
                 4 * points[_numIterations - 3] + points[_numIterations - 4]) *
                invDeltaT;
        }
        else {
            for (int i = 0; i < _numIterations; ++i) {
                acceleration[i] = 0;
            }
        }
    }
};

}  // namespace EVAA
