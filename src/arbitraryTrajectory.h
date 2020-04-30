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

        _roadAngles = Math::malloc<T>(numIterations);
        _roadAngularAcceleration = Math::malloc<T>(numIterations);
        _normedDistance = Math::malloc<T>(numIterations);

        _roadAccelerationX = Math::malloc<T>(numIterations);
        _roadAccelerationY = Math::malloc<T>(numIterations);

        _tyreAccelerationsX_fl = Math::malloc<T>(numIterations);
        _tyreAccelerationsX_fr = Math::malloc<T>(numIterations);
        _tyreAccelerationsX_rl = Math::malloc<T>(numIterations);
        _tyreAccelerationsX_rr = Math::malloc<T>(numIterations);

        _tyreAccelerationsY_fl = Math::malloc<T>(numIterations);
        _tyreAccelerationsY_fr = Math::malloc<T>(numIterations);
        _tyreAccelerationsY_rl = Math::malloc<T>(numIterations);
        _tyreAccelerationsY_rr = Math::malloc<T>(numIterations);

        _tyreAccelerationsZ_fl = Math::malloc<T>(numIterations);
        _tyreAccelerationsZ_fr = Math::malloc<T>(numIterations);
        _tyreAccelerationsZ_rl = Math::malloc<T>(numIterations);
        _tyreAccelerationsZ_rr = Math::malloc<T>(numIterations);
    }

    /**
     * Destructor
     */
    ~arbitraryTrajectory() {
        Math::free(_roadPointsX);
        Math::free(_roadPointsY);
        Math::free(_roadAccelerationX);
        Math::free(_roadAccelerationY);
        Math::free(_roadAngles);
        Math::free(_roadAngularAcceleration);
        Math::free(_normedDistance);
        Math::free(_tyreAccelerationsX_fl);
        Math::free(_tyreAccelerationsX_fr);
        Math::free(_tyreAccelerationsX_rl);
        Math::free(_tyreAccelerationsX_rr);
        Math::free(_tyreAccelerationsY_fl);
        Math::free(_tyreAccelerationsY_fr);
        Math::free(_tyreAccelerationsY_rl);
        Math::free(_tyreAccelerationsY_rr);
        Math::free(_tyreAccelerationsZ_fl);
        Math::free(_tyreAccelerationsZ_fr);
        Math::free(_tyreAccelerationsZ_rl);
        Math::free(_tyreAccelerationsZ_rr);
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
        calculateAccelerations(_roadAccelerationY, _roadPointsY);
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
        Math::free(legPointsX_fl);
        Math::free(legPointsX_fr);
        Math::free(legPointsX_rl);
        Math::free(legPointsX_rr);

        Math::free(legPointsY_fl);
        Math::free(legPointsY_fr);
        Math::free(legPointsY_rl);
        Math::free(legPointsY_rr);
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

            _tyreAccelerationsZ_fl[0] =
                (2 * pos_fl_prev + 5 * pos_fl + 4 * pos_fl_next + pos_fl_nnext) * invDeltaT;
            _tyreAccelerationsZ_fr[0] =
                (2 * pos_fr_prev + 5 * pos_fr + 4 * pos_fr_next + pos_fr_nnext) * invDeltaT;
            _tyreAccelerationsZ_rl[0] =
                (2 * pos_rl_prev + 5 * pos_rl + 4 * pos_rl_next + pos_rl_nnext) * invDeltaT;
            _tyreAccelerationsZ_rr[0] =
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

                _tyreAccelerationsZ_fl[i] = (pos_fl_prev - 2 * pos_fl + pos_fl_next) * invDeltaT;
                _tyreAccelerationsZ_fr[i] = (pos_fr_prev - 2 * pos_fr + pos_fr_next) * invDeltaT;
                _tyreAccelerationsZ_rl[i] = (pos_rl_prev - 2 * pos_rl + pos_rl_next) * invDeltaT;
                _tyreAccelerationsZ_rr[i] = (pos_rr_prev - 2 * pos_rr + pos_rr_next) * invDeltaT;
            }
            pos_fl_nnext =
                _amplitudeLeft *
                std::sin(_frequencyLeft * (_normedDistance[_numIterations - 1] + _phaseShift_fl));
            pos_fr_nnext =
                _amplitudeRight *
                std::sin(_frequencyRight * (_normedDistance[_numIterations - 1] + _phaseShift_fr));
            pos_rl_nnext =
                _amplitudeLeft *
                std::sin(_frequencyLeft * (_normedDistance[_numIterations - 1] + _phaseShift_rl));
            pos_rr_nnext =
                _amplitudeRight *
                std::sin(_frequencyRight * (_normedDistance[_numIterations - 1] + _phaseShift_rr));

            _tyreAccelerationsZ_fl[_numIterations - 1] =
                (2 * pos_fl_nnext + 5 * pos_fl_next + 4 * pos_fl + pos_fl_prev) * invDeltaT;
            _tyreAccelerationsZ_fr[_numIterations - 1] =
                (2 * pos_fr_nnext + 5 * pos_fr_next + 4 * pos_fr + pos_fr_prev) * invDeltaT;
            _tyreAccelerationsZ_rl[_numIterations - 1] =
                (2 * pos_rl_nnext + 5 * pos_rl_next + 4 * pos_rl + pos_rl_prev) * invDeltaT;
            _tyreAccelerationsZ_rr[_numIterations - 1] =
                (2 * pos_rr_nnext + 5 * pos_rr_next + 4 * pos_rr + pos_rr_prev) * invDeltaT;
        }
        else {
            for (int i = 0; i < _numIterations; ++i) {
                _tyreAccelerationsZ_fl[i] = 0;
                _tyreAccelerationsZ_fr[i] = 0;
                _tyreAccelerationsZ_rl[i] = 0;
                _tyreAccelerationsZ_rr[i] = 0;
            }
        }
    }

    /**
     * \brief calculates the Acceleration acting on the tyre due to bumpy road
     * \param time at which this happens
     */
    T getVerticalRoadAccelerationsFrontLeft(size_t iteration) {
        return _tyreAccelerationsZ_fl[iteration];
    }

    /**
     * \brief calculates the Acceleration acting on the tyre due to bumpy road
     * \param time at which this happens
     */
    T getVerticalRoadAccelerationsFrontRight(size_t iteration) {
        return _tyreAccelerationsZ_rl[iteration];
    }

    /**
     * \brief calculates the Acceleration acting on the tyre due to bumpy road
     * \param time at which this happens
     */
    T getVerticalRoadAccelerationsRearLeft(size_t iteration) {
        return _tyreAccelerationsZ_fr[iteration];
    }

    /**
     * \brief calculates the Acceleration acting on the tyre due to bumpy road
     * \param time at which this happens
     */
    T getVerticalRoadAccelerationsRearRight(size_t iteration) {
        return _tyreAccelerationsZ_rr[iteration];
    }

    

    /**
     * \brief calculates the position of the t<ew
     * \param time at which this happens
     */
    T getVerticalPositionFrontLeft(size_t iteration) {
        return _amplitudeLeft * std::sin(_frequencyLeft * (_normedDistance[i] + _phaseShift_fl));
    }


     /**
     * \brief calculates the position of the t<ew
     * \param time at which this happens
     */
    T getVerticalPositionFrontRight(size_t iteration) {
        return _amplitudeRight * std::sin(_frequencyLeft * (_normedDistance[i] + _phaseShift_fr));
    }


     /**
     * \brief calculates the position of the t<ew
     * \param time at which this happens
     */
    T getVerticalPositionRearLeft(size_t iteration) {
        return _amplitudeLeft * std::sin(_frequencyLeft * (_normedDistance[i] + _phaseShift_rl));
    }

     /**
     * \brief calculates the position of the t<ew
     * \param time at which this happens
     */
    T getVerticalPositionRearRight(size_t iteration) {
        return _amplitudeRight * std::sin(_frequencyLeft * (_normedDistance[i] + _phaseShift_rr));
    }

    /**
     * \param iteration index of the current time iteration
     * \return vector the 2D vector to be written too [GC:XY]
     */
    void getLagrangianAccelerationsCenterOfGravity(size_t iteration, T* vector) {
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
     * \return vector the 2D vector to be written too [leg:XY]
     */
    void getLagrangianAccelerationsFrontLeft(size_t iteration, T* vector) {
        vector[0] = _tyreAccelerationsX_fl[i];
        vector[1] = _tyreAccelerationsY_fl[i];
    }

    /**
     * \param iteration index of the current time iteration
     * \return vector the 2D vector to be written too [leg:XY]
     */
    void getLagrangianAccelerationsFrontRight(size_t iteration, T* vector) {
        vector[0] = _tyreAccelerationsX_fr[i];
        vector[1] = _tyreAccelerationsY_fr[i];
    }

    /**
     * \param iteration index of the current time iteration
     * \return vector the 2D vector to be written too [leg:XY]
     */
    void getLagrangianAccelerationsRearLeft(size_t iteration, T* vector) {
        vector[0] = _tyreAccelerationsX_rl[i];
        vector[1] = _tyreAccelerationsY_rl[i];
    }

    /**
     * \param iteration index of the current time iteration
     * \return vector the 2D vector to be written too [leg:XY]
     */
    void getLagrangianAccelerationsRearRight(size_t iteration, T* vector) {
        vector[0] = _tyreAccelerationsX_rr[i];
        vector[1] = _tyreAccelerationsY_rr[i];
    }

    /**
     * \param iteration index of the current time iteration
     * \return vector the 2D vector to be written too [GC:XY]
     */
    void getLagrangianPositionCenterOfGravity(size_t iteration, T* vector) {
        vector[0] = _roadPointsX[i];
        vector[1] = _roadPointsY[i];
    }

    void testCircularTrajectory() {
        T center[2] = {2, -5};
        T radius = 30;
        T velocity = 3;

        T angularVelocity = velocity / radius;

        T x;

        // create circular road trajectory
        for (int i = 0; i < _numIterations; ++i) {
            x = _delta_t * i * velocity;
            _roadPointsX[i] = center[0] + radius * std::sin(x / radius);
            _roadPointsY[i] = center[1] + radius * std::cos(x / radius);
        }

        T longitude = 3;
        T latitude = 2;

        T tolerance = 1e-4;

         _amplitudeRight = 1;    // amplitude of the sinusoidal curve acting on the right tyres
         _amplitudeLeft = 1.5;   // amplitude of the sinusoidal curve acting on the left tyres

         _frequencyRight = 0.5;  // frequency of the sinusoidal curve acting on the right tyres
         _frequencyLeft = 0.8;   // frequency of the sinusoidal curve acting on the left tyres


        

        calculateTravelledDistance();
        calculateAngles();
        calculateAccelerationsCenterOfGravity();
        calculateAccelerationsLegs(longitude, longitude, longitude, longitude, latitude, latitude,
                                   latitude, latitude);

        for (int i = 0; i < _numIterations; ++i) {
            // check the correctness of the distance
            if (std::abs(_normedDistance[i] - velocity * _delta_t * i) > tolerance) {
                std::cout << "calculateTravelledDistance failed in arbitraryTrajectory"
                          << std::endl;
                std::cout << "difference: "
                          << std::abs(_normedDistance[i] - radius * velocity * _delta_t * i)
                          << _normedDistance[i] << std::endl;
            }

            // check correctness of the CoG acceleration
            T vectorToCenterX = (center[0] - _roadPointsX[i]) / radius;
            T vectorToCenterY = (center[1] - _roadPointsY[i]) / radius;

            T accelerationNorm = std::sqrt(_roadAccelerationX[i] * _roadAccelerationX[i] +
                                           _roadAccelerationY[i] * _roadAccelerationY[i]);

            T normedAccelerationX = _roadAccelerationX[i] / accelerationNorm;
            T normedAccelerationY = _roadAccelerationY[i] / accelerationNorm;

            // does it point to the center ?
            if ((i > 0) && (i < _numIterations - 1) &&
                ((std::abs(vectorToCenterX - normedAccelerationX) > tolerance) ||
                 (std::abs(vectorToCenterY - normedAccelerationY) > tolerance))) {
                std::cout << "calculateAccelerationsCenterOfGravity failed in arbitraryTrajectory "
                             "-- wrong direction"
                          << std::endl;
            }

            // check correctness of the fl acceleration
            accelerationNorm = std::sqrt(_tyreAccelerationsX_fl[i] * _tyreAccelerationsX_fl[i] +
                                         _tyreAccelerationsY_fl[i] * _tyreAccelerationsY_fl[i]);

            T legRadius = radius - latitude;
            T legVelocity = legRadius * angularVelocity;

            // is its norm correct ?
            if ((i > 0) && (i < _numIterations - 1) &&
                (std::abs(legVelocity * legVelocity / legRadius - accelerationNorm) > tolerance)) {
                std::cout << "calculateAccelerationsCenterOfGravity failed in arbitraryTrajectory "
                             "-- wrong norm of the front-left wheel"
                          << std::endl;
                std::cout << "acceleration: " << accelerationNorm
                          << " expected: " << velocity * velocity / radius << std::endl;
            }


            // is its norm correct ?
            if ((i > 0) && (i < _numIterations - 1) &&
                (std::abs(velocity * velocity / radius - accelerationNorm) > tolerance)) {
                std::cout << "calculateAccelerationsCenterOfGravity failed in arbitraryTrajectory "
                             "-- wrong norm"
                          << std::endl;
                std::cout << "acceleration: " << accelerationNorm
                          << " expected: " << velocity * velocity / radius << std::endl;
            }

            // check correctness of the fr acceleration
            accelerationNorm = std::sqrt(_tyreAccelerationsX_fr[i] * _tyreAccelerationsX_fr[i] +
                                           _tyreAccelerationsY_fr[i] * _tyreAccelerationsY_fr[i]);

            legRadius = radius + latitude;
            legVelocity = legRadius * angularVelocity;

            // is its norm correct ?
            if ((i > 0) && (i < _numIterations - 1) &&
                (std::abs(legVelocity * legVelocity / legRadius - accelerationNorm) > tolerance)) {
                std::cout << "calculateAccelerationsCenterOfGravity failed in arbitraryTrajectory "
                             "-- wrong norm of the front-left wheel"
                          << std::endl;
                std::cout << "acceleration: " << accelerationNorm
                          << " expected: " << velocity * velocity / radius << std::endl;
            }

            // check correctness of the rl acceleration
            accelerationNorm = std::sqrt(_tyreAccelerationsX_rl[i] * _tyreAccelerationsX_rl[i] +
                                         _tyreAccelerationsY_rl[i] * _tyreAccelerationsY_rl[i]);

            legRadius = radius - latitude;
            legVelocity = legRadius * angularVelocity;

            // is its norm correct ?
            if ((i > 0) && (i < _numIterations - 1) &&
                (std::abs(legVelocity * legVelocity / legRadius - accelerationNorm) > tolerance)) {
                std::cout << "calculateAccelerationsCenterOfGravity failed in arbitraryTrajectory "
                             "-- wrong norm of the front-left wheel"
                          << std::endl;
                std::cout << "acceleration: " << accelerationNorm
                          << " expected: " << velocity * velocity / radius << std::endl;
            }

            // check correctness of the rr acceleration
            accelerationNorm = std::sqrt(_tyreAccelerationsX_rr[i] * _tyreAccelerationsX_rr[i] +
                                         _tyreAccelerationsY_rr[i] * _tyreAccelerationsY_rr[i]);

            legRadius = radius + latitude;
            legVelocity = legRadius * angularVelocity;

            // is its norm correct ?
            if ((i > 0) && (i < _numIterations - 1) &&
                (std::abs(legVelocity * legVelocity / legRadius - accelerationNorm) > tolerance)) {
                std::cout << "calculateAccelerationsCenterOfGravity failed in arbitraryTrajectory "
                             "-- wrong norm of the front-left wheel"
                          << std::endl;
                std::cout << "acceleration: " << accelerationNorm
                          << " expected: " << velocity * velocity / radius << std::endl;
            }

        }
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

    T* _tyreAccelerationsZ_fl;  // Acceleration that acts on the tyre in Z direction
    T* _tyreAccelerationsZ_fr;  // Acceleration that acts on the tyre in Z direction
    T* _tyreAccelerationsZ_rl;  // Acceleration that acts on the tyre in Z direction
    T* _tyreAccelerationsZ_rr;  // Acceleration that acts on the tyre in Z direction

    // Lagrangian frame
    T* _roadPointsX;  // X-coordinates of the trajectory at all iterations
    T* _roadPointsY;  // Y-coordinates of the trajectory at all iterations

    T* _tyreAccelerationsX_fl;  // Acceleration that acts on the tyre in X direction
    T* _tyreAccelerationsX_fr;  // Acceleration that acts on the tyre in X direction
    T* _tyreAccelerationsX_rl;  // Acceleration that acts on the tyre in X direction
    T* _tyreAccelerationsX_rr;  // Acceleration that acts on the tyre in X direction

    T* _tyreAccelerationsY_fl;  // Acceleration that acts on the tyre in Y direction
    T* _tyreAccelerationsY_fr;  // Acceleration that acts on the tyre in Y direction
    T* _tyreAccelerationsY_rl;  // Acceleration that acts on the tyre in Y direction
    T* _tyreAccelerationsY_rr;  // Acceleration that acts on the tyre in Y direction



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
                -(2 * points[0] + 5 * points[1] + 4 * points[2] + points[3]) * invDeltaT;

            for (int i = 1; i < _numIterations - 1; ++i) {
                acceleration[i] = (points[i - 1] - 2 * points[i] + points[i + 1]) * invDeltaT;
            }

            acceleration[_numIterations - 1] =
                -(2 * points[_numIterations - 1] + 5 * points[_numIterations - 2] +
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
