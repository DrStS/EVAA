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
    ArbitraryTrajectory(size_t numIterations, T delta_t, T amplitude_right, T amplitude_left, T period_right, T period_left, T shift_right, T shift_left, T* initialUpperSpringLengths, T* initialLowerSpringLengths, T* latitudes, T* longitudes) : _numIterations(numIterations), _delta_t(delta_t) {
        // Read Eulerian force parameters
        _frequencyRight = (period_right == 0) ? 0 : 2 * Constants::PI / period_right;
        _frequencyLeft = (period_left == 0) ? 0 : 2 * Constants::PI / period_left;

        _amplitudeRight = amplitude_right;
        _amplitudeLeft = amplitude_left;

        _phaseShift_fl = (shift_left == 0) ? 0 : 2 * Constants::PI / shift_left;
        _phaseShift_fr = (shift_right == 0) ? 0 : 2 * Constants::PI / shift_right;
        _phaseShift_rl = (shift_left == 0) ? 0 : 2 * Constants::PI / shift_left;
        _phaseShift_rr = (shift_right == 0) ? 0 : 2 * Constants::PI / shift_right;

        // Read simulation parameters
        _initialUpperSpringLength_fl = initialUpperSpringLengths[Constants::FRONT_LEFT];
        _initialUpperSpringLength_fr = initialUpperSpringLengths[Constants::FRONT_RIGHT];
        _initialUpperSpringLength_rl = initialUpperSpringLengths[Constants::REAR_LEFT];
        _initialUpperSpringLength_rr = initialUpperSpringLengths[Constants::REAR_RIGHT];

        _initialLowerSpringLength_fl = initialLowerSpringLengths[Constants::FRONT_LEFT];
        _initialLowerSpringLength_fr = initialLowerSpringLengths[Constants::FRONT_RIGHT];
        _initialLowerSpringLength_rl = initialLowerSpringLengths[Constants::REAR_LEFT];
        _initialLowerSpringLength_rr = initialLowerSpringLengths[Constants::REAR_RIGHT];

        _l_lat_fl = latitudes[Constants::FRONT_LEFT];
        _l_lat_fr = latitudes[Constants::FRONT_RIGHT];
        _l_lat_rl = latitudes[Constants::REAR_LEFT];
        _l_lat_rr = latitudes[Constants::REAR_RIGHT];

        _l_long_fl = longitudes[Constants::FRONT_LEFT];
        _l_long_fr = longitudes[Constants::FRONT_RIGHT];
        _l_long_rl = longitudes[Constants::REAR_LEFT];
        _l_long_rr = longitudes[Constants::REAR_RIGHT];

        // allocate memory
        _roadPointsX = Math::malloc<T>(numIterations + 1);
        _roadPointsY = Math::malloc<T>(numIterations + 1);

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

        _legPointsX_fl = Math::malloc<T>(_numIterations + 1);
        _legPointsX_fr = Math::malloc<T>(_numIterations + 1);
        _legPointsX_rl = Math::malloc<T>(_numIterations + 1);
        _legPointsX_rr = Math::malloc<T>(_numIterations + 1);

        _legPointsY_fl = Math::malloc<T>(_numIterations + 1);
        _legPointsY_fr = Math::malloc<T>(_numIterations + 1);
        _legPointsY_rl = Math::malloc<T>(_numIterations + 1);
        _legPointsY_rr = Math::malloc<T>(_numIterations + 1);

        _legPointsZ_fl = Math::malloc<T>(_numIterations + 1);
        _legPointsZ_fr = Math::malloc<T>(_numIterations + 1);
        _legPointsZ_rl = Math::malloc<T>(_numIterations + 1);
        _legPointsZ_rr = Math::malloc<T>(_numIterations + 1);
    }

    /**
     * Destructor
     */
    ~ArbitraryTrajectory() {
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
        Math::free<T>(_legPointsX_fl);
        Math::free<T>(_legPointsX_fr);
        Math::free<T>(_legPointsX_rl);
        Math::free<T>(_legPointsX_rr);

        Math::free<T>(_legPointsY_fl);
        Math::free<T>(_legPointsY_fr);
        Math::free<T>(_legPointsY_rl);
        Math::free<T>(_legPointsY_rr);

        Math::free<T>(_legPointsZ_fl);
        Math::free<T>(_legPointsZ_fr);
        Math::free<T>(_legPointsZ_rl);
        Math::free<T>(_legPointsZ_rr);
    }

    /**
     * \brief calcuate the addition shift due to the geometry of the car
     */
    void calculateTyreShifts() {
        _phaseShift_fl += _l_long_fl;
        _phaseShift_fr += _l_long_fr;
        _phaseShift_rl -= _l_long_rl;
        _phaseShift_rr -= _l_long_rr;
    }

    /**
    * \brief interpolate all intermediate road points @Felix, help me out
    * \param numProvidedPoints number of road points from the XML
    * \param providedPointsX array of all X-coordinates of the points from the
    XML
    * \param providedPointsY array of all Y-coordinates of the points from the
    XML
    * \param providedTimes times, at which the XML points should be reached
    * \param initialVelocity such that in the first segment the acceleration is
    calculated correctly

    * This means, that between two XML points  [X(i), Y(i)] and [X(�+1),
    Y(i+1)], there will be n=(times(�+1) - times(i)) / delta_t true road points
     */
    void interpolateRoadPoints(size_t numProvidedPoints, T* providedPointsX, T* providedPointsY, T* providedTimes) {
        for (auto i = 0; i < numProvidedPoints - 1; i++) {
            if (providedTimes[i] >= providedTimes[i + 1]) {
                std::cout << providedTimes[i] << std::endl;
                std::cout << providedTimes[i + 1] << std::endl;
                throw std::domain_error("time is not increasing in road trajectory: " + std::to_string(providedTimes[i]) + " <= " + std::to_string(providedTimes[i + 1]));
            }
        }
        T maxTime = providedTimes[numProvidedPoints - 1];
        int numInterpolationPoints = maxTime / _delta_t + 1;
        T* timeInterpolationPoints = Math::malloc<T>(_numIterations + 1);
        T interpolationTime = _delta_t * _numIterations;
        // when interpolationTime is longer than road defined exit
        if (interpolationTime > maxTime) {
            throw std::domain_error("Simulation time is longer than road is defined: " + std::to_string(interpolationTime) + " > " + std::to_string(maxTime));
        }
        else {
            // fill in time axis
            for (auto i = 0; i < _numIterations + 1; i++) {
                timeInterpolationPoints[i] = i * _delta_t;
            }
            // do the interpolation
            interpolateAxis(numProvidedPoints, providedPointsX, providedTimes, timeInterpolationPoints, _roadPointsX);
            interpolateAxis(numProvidedPoints, providedPointsY, providedTimes, timeInterpolationPoints, _roadPointsY);
#ifdef WRITECSV
            T* noZComponent = nullptr;
            IO::writeRoadTrajectoryCSV(_roadPointsX, _roadPointsY, noZComponent, _numIterations + 1, "C:\\software\\repos\\EVAA\\output\\Trajectory.txt");
#endif  // WRITECSV
        }
        Math::free(timeInterpolationPoints);
    }

    /**
     * \brief interpolate the road at all time steps
     *
     * \param [in] numProvidedPoints
     * \param [in] providedPoints
     * \param [in] providedTimes
     * \param [in] interpolationAxis
     * \param [out] axis
     */
    void interpolateAxis(size_t numProvidedPoints, T* providedPoints, T* providedTimes, T* interpolationAxis, T* axis) {
        MKL_INT order = DF_PP_CUBIC;
        const MKL_INT dorder[1] = {1};
        DFTaskPtr Task;
        T* coeff = Math::malloc<T>((numProvidedPoints - 1) * order);
        int err = 0;
        err = Math::dfNewTask1D(&Task, numProvidedPoints, providedTimes, DF_NON_UNIFORM_PARTITION, 1, providedPoints, DF_NO_HINT);
        Math::dfCheckError(err);

        err = Math::dfEditPPSpline1D<T>(Task, order, DF_PP_NATURAL, DF_BC_FREE_END, nullptr, DF_NO_IC, nullptr, coeff, DF_NO_HINT);
        Math::dfCheckError(err);

        /* Construct spline using STD method */
        err = Math::dfConstruct1D<T>(Task, DF_PP_SPLINE, DF_METHOD_STD);
        Math::dfCheckError(err);

        for (auto i = 0; i < _numIterations + 1; i++) {
            err = Math::dfInterpolate1D<T>(Task, DF_INTERP, DF_METHOD_PP, 1, &interpolationAxis[i], DF_NO_HINT, 1, dorder, nullptr, &axis[i], DF_NO_HINT, 0);
            Math::dfCheckError(err);
        }
        Math::free(coeff);
        dfDeleteTask(&Task);
    }

    /**
     * \brief calculate the travelled distance at each trajectory point
     */
    void calculateTravelledDistance() {
        _normedDistance[0] = 0;
        for (int i = 1; i < _numIterations + 1; ++i) {
            _normedDistance[i] = _normedDistance[i - 1] + sqrt((_roadPointsX[i] - _roadPointsX[i - 1]) * (_roadPointsX[i] - _roadPointsX[i - 1]) + (_roadPointsY[i] - _roadPointsY[i - 1]) * (_roadPointsY[i] - _roadPointsY[i - 1]));
        }
    }

    /**
     * \brief use a second order scheme to calculate all angles on the
     * trajectory points
     */
    void calculateAngles() {
        T invDeltaT = 1. / _delta_t;
        if (_numIterations > 2) {
            T directionX = -1.5 * _roadPointsX[0] + 2 * _roadPointsX[1] - 0.5 * _roadPointsX[2];
            T directionY = -1.5 * _roadPointsY[0] + 2 * _roadPointsY[1] - 0.5 * _roadPointsY[2];
            T directionNorm = std::sqrt(directionX * directionX + directionY * directionY);
            if (directionNorm == 0) {
                _roadAngles[0] = 0;
            }
            else if (directionY >= 0) {
                _roadAngles[0] = std::acos(directionX / directionNorm);
            }
            else {
                _roadAngles[0] = -std::acos(directionX / directionNorm);
            }

            for (int i = 1; i < _numIterations; ++i) {
                directionX = -0.5 * _roadPointsX[i - 1] + 0.5 * _roadPointsX[i + 1];
                directionY = -0.5 * _roadPointsY[i - 1] + 0.5 * _roadPointsY[i + 1];
                directionNorm = std::sqrt(directionX * directionX + directionY * directionY);
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

            directionX = 1.5 * _roadPointsX[_numIterations] + -2 * _roadPointsX[_numIterations - 1] - -0.5 * _roadPointsX[_numIterations - 2];
            directionY = 1.5 * _roadPointsY[_numIterations] + -2 * _roadPointsY[_numIterations - 1] - -0.5 * _roadPointsY[_numIterations - 2];
            directionNorm = std::sqrt(directionX * directionX + directionY * directionY);
            if (directionNorm == 0) {
                _roadAngles[_numIterations] = _roadAngles[_numIterations - 1];
            }
            else if (directionY >= 0) {
                _roadAngles[_numIterations] = std::acos(directionX / directionNorm);
            }
            else {
                _roadAngles[_numIterations] = -std::acos(directionX / directionNorm);
            }
        }
        else {
            for (int i = 0; i < _numIterations + 1; ++i) {
                _roadAngles[i] = 0;
            }
        }
    }

    /**
     * \brief use a second order scheme to calculate all AngularAccelerations
     * and Accelerations between the trajectory points
     */
    void calculateAccelerationsCenterOfGravity() {
        calculateAccelerations(_roadAccelerationX, _roadPointsX);
        calculateAccelerations(_roadAccelerationY, _roadPointsY);
        calculateAccelerations(_roadAngularAcceleration, _roadAngles);
    }

    /**
     * \brief use a second order scheme to calculate all AngularAccelerations
     * and Accelerations between the trajectory points for the legs
     */
    void calculateAccelerationsLegs() {
        // get true local coordinates
        T localXcoordinates[Constants::NUM_LEGS] = {_l_long_fl, _l_long_fr, -_l_long_rl, -_l_long_rr};
        T localYcoordinates[Constants::NUM_LEGS] = {-_l_lat_fl, _l_lat_fr, -_l_lat_rl, _l_lat_rr};

        // get all leg positions
        for (int i = 0; i < _numIterations + 1; ++i) {
            T c = std::cos(_roadAngles[i]);
            T s = std::sin(_roadAngles[i]);

            _legPointsX_fl[i] = _roadPointsX[i] + localXcoordinates[0] * c - localYcoordinates[0] * s;
            _legPointsY_fl[i] = _roadPointsY[i] + localXcoordinates[0] * s + localYcoordinates[0] * c;

            _legPointsX_fr[i] = _roadPointsX[i] + localXcoordinates[1] * c - localYcoordinates[1] * s;
            _legPointsY_fr[i] = _roadPointsY[i] + localXcoordinates[1] * s + localYcoordinates[1] * c;

            _legPointsX_rl[i] = _roadPointsX[i] + localXcoordinates[2] * c - localYcoordinates[2] * s;
            _legPointsY_rl[i] = _roadPointsY[i] + localXcoordinates[2] * s + localYcoordinates[2] * c;

            _legPointsX_rr[i] = _roadPointsX[i] + localXcoordinates[3] * c - localYcoordinates[3] * s;
            _legPointsY_rr[i] = _roadPointsY[i] + localXcoordinates[3] * s + localYcoordinates[3] * c;
        }

        // calculate tyre forces
        calculateAccelerations(_tyreAccelerationsX_fl, _legPointsX_fl);
        calculateAccelerations(_tyreAccelerationsY_fl, _legPointsY_fl);

        calculateAccelerations(_tyreAccelerationsX_fr, _legPointsX_fr);
        calculateAccelerations(_tyreAccelerationsY_fr, _legPointsY_fr);

        calculateAccelerations(_tyreAccelerationsX_rl, _legPointsX_rl);
        calculateAccelerations(_tyreAccelerationsY_rl, _legPointsY_rl);

        calculateAccelerations(_tyreAccelerationsX_rr, _legPointsX_rr);
        calculateAccelerations(_tyreAccelerationsY_rr, _legPointsY_rr);

#ifdef WRITECSV
        IO::writeRoadTrajectoryCSV(_legPointsX_fl, _legPointsY_fl, _legPointsZ_fl, _numIterations + 1, "C:\\software\\repos\\EVAA\\output\\LegFl.txt");
        IO::writeRoadTrajectoryCSV(_legPointsX_fr, _legPointsY_fr, _legPointsZ_fr, _numIterations + 1, "C:\\software\\repos\\EVAA\\output\\LegFr.txt");
        IO::writeRoadTrajectoryCSV(_legPointsX_rl, _legPointsY_rl, _legPointsZ_rl, _numIterations + 1, "C:\\software\\repos\\EVAA\\output\\LegRl.txt");
        IO::writeRoadTrajectoryCSV(_legPointsX_rr, _legPointsY_rr, _legPointsZ_rr, _numIterations + 1, "C:\\software\\repos\\EVAA\\output\\LegRr.txt");
#endif  // WRITECSV
    }

    /**
     * \brief use a second order scheme to calculate the Accelerations that act
     * on the tyre on the trajectory points
     */
    void calculateVerticalAccelerations() {
        T invDeltaT = 1. / (_delta_t * _delta_t);
        if (_numIterations > 3) {
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

            _tyreAccelerationsZ_fl[0] = (2 * pos_fl_prev + 5 * pos_fl + 4 * pos_fl_next + pos_fl_nnext) * invDeltaT;
            _tyreAccelerationsZ_fr[0] = (2 * pos_fr_prev + 5 * pos_fr + 4 * pos_fr_next + pos_fr_nnext) * invDeltaT;
            _tyreAccelerationsZ_rl[0] = (2 * pos_rl_prev + 5 * pos_rl + 4 * pos_rl_next + pos_rl_nnext) * invDeltaT;
            _tyreAccelerationsZ_rr[0] = (2 * pos_rr_prev + 5 * pos_rr + 4 * pos_rr_next + pos_rr_nnext) * invDeltaT;

            for (int i = 1; i < _numIterations; ++i) {
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

                _tyreAccelerationsZ_fl[i] = (pos_fl_prev - 2 * pos_fl + pos_fl_next) * invDeltaT;
                _tyreAccelerationsZ_fr[i] = (pos_fr_prev - 2 * pos_fr + pos_fr_next) * invDeltaT;
                _tyreAccelerationsZ_rl[i] = (pos_rl_prev - 2 * pos_rl + pos_rl_next) * invDeltaT;
                _tyreAccelerationsZ_rr[i] = (pos_rr_prev - 2 * pos_rr + pos_rr_next) * invDeltaT;
            }
            pos_fl_nnext = _amplitudeLeft * std::sin(_frequencyLeft * (_normedDistance[_numIterations] + _phaseShift_fl));
            pos_fr_nnext = _amplitudeRight * std::sin(_frequencyRight * (_normedDistance[_numIterations] + _phaseShift_fr));
            pos_rl_nnext = _amplitudeLeft * std::sin(_frequencyLeft * (_normedDistance[_numIterations] + _phaseShift_rl));
            pos_rr_nnext = _amplitudeRight * std::sin(_frequencyRight * (_normedDistance[_numIterations] + _phaseShift_rr));

            _tyreAccelerationsZ_fl[_numIterations] = (2 * pos_fl_nnext + 5 * pos_fl_next + 4 * pos_fl + pos_fl_prev) * invDeltaT;
            _tyreAccelerationsZ_fr[_numIterations] = (2 * pos_fr_nnext + 5 * pos_fr_next + 4 * pos_fr + pos_fr_prev) * invDeltaT;
            _tyreAccelerationsZ_rl[_numIterations] = (2 * pos_rl_nnext + 5 * pos_rl_next + 4 * pos_rl + pos_rl_prev) * invDeltaT;
            _tyreAccelerationsZ_rr[_numIterations] = (2 * pos_rr_nnext + 5 * pos_rr_next + 4 * pos_rr + pos_rr_prev) * invDeltaT;
        }
        else {
            for (int i = 0; i < _numIterations + 1; ++i) {
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
     * \param mass of the object the force is applied
     */
    T getVerticalRoadForcesFrontLeft(size_t iteration, T mass) { return mass * _tyreAccelerationsZ_fl[iteration]; }

    /**
     * \brief calculates the force acting on the tyre due to bumpy road
     * \param time at which this happens
     * \param mass of the object the force is applied
     */
    T getVerticalRoadForcesFrontRight(size_t iteration, T mass) { return mass * _tyreAccelerationsZ_rl[iteration]; }

    /**
     * \brief calculates the Force acting on the tyre due to bumpy road
     * \param time at which this happens
     * \param mass of the object the force is applied
     */
    T getVerticalRoadForcesRearLeft(size_t iteration, T mass) { return mass * _tyreAccelerationsZ_fr[iteration]; }

    /**
     * \brief calculates the Force acting on the tyre due to bumpy road
     * \param time at which this happens
     * \param mass of the object the force is applied
     */
    T getVerticalRoadForcesRearRight(size_t iteration, T mass) { return mass * _tyreAccelerationsZ_rr[iteration]; }

    /**
     * \brief calculate vertical positions
     */
    void calculateVerticalPositionsLegs() {
        for (int i = 0; i < _numIterations + 1; ++i) {
            _legPointsZ_fl[i] = _amplitudeLeft * std::sin(_frequencyLeft * (_normedDistance[i] + _phaseShift_fl)) 
                - _initialUpperSpringLength_fl - _initialLowerSpringLength_fl;
            _legPointsZ_fr[i] = _amplitudeRight * std::sin(_frequencyLeft * (_normedDistance[i] + _phaseShift_fr))
            -_initialUpperSpringLength_fr - _initialLowerSpringLength_fr;
            _legPointsZ_rl[i] = _amplitudeLeft * std::sin(_frequencyLeft * (_normedDistance[i] + _phaseShift_rl))
            -_initialUpperSpringLength_rl - _initialLowerSpringLength_rl;
            _legPointsZ_rr[i] = _amplitudeRight * std::sin(_frequencyLeft * (_normedDistance[i] + _phaseShift_rr))
            -_initialUpperSpringLength_rr - _initialLowerSpringLength_rr;
        }
    }

    /**
     * \brief calculates the position of the t<ew
     * \param time at which this happens
     */
    T getVerticalPositionFrontLeft(size_t iteration) { return _legPointsZ_fl[iteration]; }

    /**
     * \brief calculates the position of the t<ew
     * \param time at which this happens
     */
    T getVerticalPositionFrontRight(size_t iteration) { return _legPointsZ_fr[iteration]; }

    /**
     * \brief calculates the position of the t<ew
     * \param time at which this happens
     */
    T getVerticalPositionRearLeft(size_t iteration) { return _legPointsZ_rl[iteration]; }

    /**
     * \brief calculates the position of the t<ew
     * \param time at which this happens
     */
    T getVerticalPositionRearRight(size_t iteration) { return _legPointsZ_rr[iteration]; }

    /**
     * \param iteration index of the current time iteration
     * \return vector the 2D vector to be written too [GC:XY]
     * \param mass of the object the force is applied
     */
    void getLagrangianForcesCenterOfGravity(size_t iteration, T* vector, T mass) {
        vector[0] = mass * _roadAccelerationX[i];
        vector[1] = mass * _roadAccelerationY[i];
    }
    /**
     * \param iteration index of the current time iteration
     * \param Z-component of the moment of inertia
     * \return the Z component [GC:Z]
     */
    T getLagrangianTorque(size_t iteration, T momentOfInertia) { return momentOfInertia * _roadAngularAcceleration[i]; }

    /**
     * \param iteration index of the current time iteration
     * \param mass of the object the force is applied
     * \return vector the 2D vector to be written too [leg:XY]
     */
    void getLagrangianForcesFrontLeft(size_t iteration, T mass, T* vector) {
        vector[0] = mass * _tyreAccelerationsX_fl[iteration];
        vector[1] = mass * _tyreAccelerationsY_fl[iteration];
    }

    /**
     * \param iteration index of the current time iteration
     * \param mass of the object the force is applied
     * \return vector the 2D vector to be written too [leg:XY]
     */
    void getLagrangianForcesFrontRight(size_t iteration, T mass, T* vector) {
        vector[0] = mass * _tyreAccelerationsX_fr[iteration];
        vector[1] = mass * _tyreAccelerationsY_fr[iteration];
    }

    /**
     * \param iteration index of the current time iteration
     * \param mass of the object the force is applied
     * \return vector the 2D vector to be written too [leg:XY]
     */
    void getLagrangianForcesRearLeft(size_t iteration, T mass, T* vector) {
        vector[0] = mass * _tyreAccelerationsX_rl[iteration];
        vector[1] = mass * _tyreAccelerationsY_rl[iteration];
    }

    /**
     * \param iteration index of the current time iteration
     * \param mass of the object the force is applied
     * \return vector the 2D vector to be written too [leg:XY]
     */
    void getLagrangianForcesRearRight(size_t iteration, T mass, T* vector) {
        vector[0] = mass * _tyreAccelerationsX_rr[iteration];
        vector[1] = mass * _tyreAccelerationsY_rr[iteration];
    }

    /**
     * \param iteration index of the current time iteration
     * \return vector the 2D vector to be written too [GC:XY]
     */
    void getLagrangianPositionCenterOfGravity(size_t iteration, T* vector) {
        vector[0] = _roadPointsX[iteration];
        vector[1] = _roadPointsY[iteration];
    }

    /**
     * \brief update the initial conditions of the car
     * \param[out] angles orientation of the car body [XYZ]
     * \param[out] wc angular velocity [XYZ]
     * \param[out] pcc position of the center of gravity [XYZ]
     * \param[out] vcc velocity of the center of gravity [XYZ]
     * \param[out] pt_fl tyre position front left [XYZ]
     * \param[out] pt_fr tyre position front right [XYZ]
     * \param[out] pt_rl tyre position rear left [XYZ]
     * \param[out] pt_rr tyre position rear right [XYZ]
     * \param[out] pw_fl wheel position front left [XYZ]
     * \param[out] pw_fr wheel position front right [XYZ]
     * \param[out] pw_rl wheel position rear left [XYZ]
     * \param[out] vw_rr wheel position rear right [XYZ]
     * \param[out] vt_fl tyre velocity front left [XYZ]
     * \param[out] vt_fr tyre velocity front right [XYZ]
     * \param[out] vt_rl tyre velocity rear left [XYZ]
     * \param[out] vt_rr tyre velocity rear right [XYZ]
     * \param[out] vw_fl wheel velocity front left [XYZ]
     * \param[out] vw_fr wheel velocity front right [XYZ]
     * \param[out] vw_rl wheel velocity rear left [XYZ]
     * \param[out] vw_rr wheel velocity rear right [XYZ]
     */
    void updateInitialConditions(T* angle, T* wc, T* pcc, T* vcc, T* pt_fl, T* pt_fr, T* pt_rl, T* pt_rr, T* pw_fl, T* pw_fr, T* pw_rl, T* pw_rr, T* vt_fl, T* vt_fr, T* vt_rl, T* vt_rr, T* vw_fl, T* vw_fr, T* vw_rl, T* vw_rr) {
        // angles
        angle[0] = 0;
        angle[1] = 0;
        angle[2] = _roadAngles[0];

        // angular velocity
        wc[0] = 0;
        wc[1] = 0;
        wc[2] = (1.5 * _roadAngles[0] - 2 * _roadAngles[1] + 0.5 * _roadAngles[2]) / _delta_t;

        // initial position
        pcc[0] = _roadPointsX[0];
        pcc[1] = _roadPointsY[0];
        pcc[2] = 0;

        // initial velocity
        vcc[0] = (-1.5 * _roadPointsX[0] + 2 * _roadPointsX[1] - 0.5 * _roadPointsX[2]) / _delta_t;
        vcc[1] = (-1.5 * _roadPointsY[0] + 2 * _roadPointsY[1] - 0.5 * _roadPointsY[2]) / _delta_t;
        vcc[2] = 0;

        // wheel positions
        pw_fl[0] = _legPointsX_fl[0];
        pw_fr[0] = _legPointsX_fr[0];
        pw_rl[0] = _legPointsX_rl[0];
        pw_rr[0] = _legPointsX_rr[0];

        pw_fl[1] = _legPointsY_fl[0];
        pw_fr[1] = _legPointsY_fr[0];
        pw_rl[1] = _legPointsY_rl[0];
        pw_rr[1] = _legPointsY_rr[0];

        pw_fl[2] = _initialUpperSpringLength_fl + _legPointsZ_fl[0];
        pw_fl[2] = _initialUpperSpringLength_fr + _legPointsZ_fr[0];
        pw_fl[2] = _initialUpperSpringLength_rl + _legPointsZ_rl[0];
        pw_fl[2] = _initialUpperSpringLength_rr + _legPointsZ_rr[0];

        // tyre positions
        pt_fl[0] = _legPointsX_fl[0];
        pt_fr[0] = _legPointsX_fr[0];
        pt_rl[0] = _legPointsX_rl[0];
        pt_rr[0] = _legPointsX_rr[0];

        pt_fl[1] = _legPointsY_fl[0];
        pt_fr[1] = _legPointsY_fr[0];
        pt_rl[1] = _legPointsY_rl[0];
        pt_rr[1] = _legPointsY_rr[0];

        pt_fl[2] = _legPointsZ_fl[0];
        pt_fr[2] = _legPointsZ_fr[0];
        pt_rl[2] = _legPointsZ_rl[0];
        pt_rr[2] = _legPointsZ_rr[0];

        // initial wheel velocity
        vw_fl[0] = (-1.5 * _legPointsX_fl[0] + 2 * _legPointsX_fl[1] - 0.5 * _legPointsX_fl[2]) / _delta_t;
        vw_fl[1] = (-1.5 * _legPointsY_fl[0] + 2 * _legPointsY_fl[1] - 0.5 * _legPointsY_fl[2]) / _delta_t;
        vw_fl[2] = (-1.5 * _legPointsZ_fl[0] + 2 * _legPointsZ_fl[1] - 0.5 * _legPointsZ_fl[2]) / _delta_t;

        vw_fr[0] = (-1.5 * _legPointsX_fr[0] + 2 * _legPointsX_fr[1] - 0.5 * _legPointsX_fr[2]) / _delta_t;
        vw_fr[1] = (-1.5 * _legPointsY_fr[0] + 2 * _legPointsY_fr[1] - 0.5 * _legPointsY_fr[2]) / _delta_t;
        vw_fr[2] = (-1.5 * _legPointsZ_fr[0] + 2 * _legPointsZ_fr[1] - 0.5 * _legPointsZ_fr[2]) / _delta_t;

        vw_rl[0] = (-1.5 * _legPointsX_rl[0] + 2 * _legPointsX_rl[1] - 0.5 * _legPointsX_rl[2]) / _delta_t;
        vw_rl[1] = (-1.5 * _legPointsY_rl[0] + 2 * _legPointsY_rl[1] - 0.5 * _legPointsY_rl[2]) / _delta_t;
        vw_rl[2] = (-1.5 * _legPointsZ_rl[0] + 2 * _legPointsZ_rl[1] - 0.5 * _legPointsZ_rl[2]) / _delta_t;

        vw_rr[0] = (-1.5 * _legPointsX_rr[0] + 2 * _legPointsX_rr[1] - 0.5 * _legPointsX_rr[2]) / _delta_t;
        vw_rr[1] = (-1.5 * _legPointsY_rr[0] + 2 * _legPointsY_rr[1] - 0.5 * _legPointsY_rr[2]) / _delta_t;
        vw_rr[2] = (-1.5 * _legPointsZ_rr[0] + 2 * _legPointsZ_rr[1] - 0.5 * _legPointsZ_rr[2]) / _delta_t;

        // initial wheel velocity
        vt_fl[0] = (-1.5 * _legPointsX_fl[0] + 2 * _legPointsX_fl[1] - 0.5 * _legPointsX_fl[2]) / _delta_t;
        vt_fl[1] = (-1.5 * _legPointsY_fl[0] + 2 * _legPointsY_fl[1] - 0.5 * _legPointsY_fl[2]) / _delta_t;
        vt_fl[2] = (-1.5 * _legPointsZ_fl[0] + 2 * _legPointsZ_fl[1] - 0.5 * _legPointsZ_fl[2]) / _delta_t;

        vt_fr[0] = (-1.5 * _legPointsX_fr[0] + 2 * _legPointsX_fr[1] - 0.5 * _legPointsX_fr[2]) / _delta_t;
        vt_fr[1] = (-1.5 * _legPointsY_fr[0] + 2 * _legPointsY_fr[1] - 0.5 * _legPointsY_fr[2]) / _delta_t;
        vt_fr[2] = (-1.5 * _legPointsZ_fr[0] + 2 * _legPointsZ_fr[1] - 0.5 * _legPointsZ_fr[2]) / _delta_t;

        vt_rl[0] = (-1.5 * _legPointsX_rl[0] + 2 * _legPointsX_rl[1] - 0.5 * _legPointsX_rl[2]) / _delta_t;
        vt_rl[1] = (-1.5 * _legPointsY_rl[0] + 2 * _legPointsY_rl[1] - 0.5 * _legPointsY_rl[2]) / _delta_t;
        vt_rl[2] = (-1.5 * _legPointsZ_rl[0] + 2 * _legPointsZ_rl[1] - 0.5 * _legPointsZ_rl[2]) / _delta_t;

        vt_rr[0] = (-1.5 * _legPointsX_rr[0] + 2 * _legPointsX_rr[1] - 0.5 * _legPointsX_rr[2]) / _delta_t;
        vt_rr[1] = (-1.5 * _legPointsY_rr[0] + 2 * _legPointsY_rr[1] - 0.5 * _legPointsY_rr[2]) / _delta_t;
        vt_rr[2] = (-1.5 * _legPointsZ_rr[0] + 2 * _legPointsZ_rr[1] - 0.5 * _legPointsZ_rr[2]) / _delta_t;
    }

    /**
     * \brief test functions for trajectory generation
     */
    void testCircularTrajectory() {
        T center[2] = {2, -5};
        T radius = 30;
        T velocity = 3;

        T angularVelocity = velocity / radius;

        T x;

        // create circular road trajectory
        for (int i = 0; i < _numIterations + 1; ++i) {
            x = _delta_t * i * velocity;
            _roadPointsX[i] = center[0] + radius * std::sin(x / radius);
            _roadPointsY[i] = center[1] + radius * std::cos(x / radius);
        }

        T longitude = 3;
        T latitude = 2;

        T tolerance = 1e-4;

        _amplitudeRight = 1;   // amplitude of the sinusoidal curve acting on the right tyres
        _amplitudeLeft = 1.5;  // amplitude of the sinusoidal curve acting on the left tyres

        _frequencyRight = 0.5;  // frequency of the sinusoidal curve acting on the right tyres
        _frequencyLeft = 0.8;   // frequency of the sinusoidal curve acting on the left tyres

        calculateTravelledDistance();
        calculateAngles();
        calculateAccelerationsCenterOfGravity();
        calculateAccelerationsLegs(longitude, longitude, longitude, longitude, latitude, latitude, latitude, latitude);

        for (int i = 0; i < _numIterations + 1; ++i) {
            // check the correctness of the distance
            if (std::abs(_normedDistance[i] - velocity * _delta_t * i) > tolerance) {
                std::cout << "calculateTravelledDistance failed in "
                             "arbitraryTrajectory"
                          << std::endl;
                std::cout << "difference: " << std::abs(_normedDistance[i] - radius * velocity * _delta_t * i) << _normedDistance[i] << std::endl;
            }

            // check correctness of the CoG acceleration
            T vectorToCenterX = (center[0] - _roadPointsX[i]) / radius;
            T vectorToCenterY = (center[1] - _roadPointsY[i]) / radius;

            T accelerationNorm = std::sqrt(_roadAccelerationX[i] * _roadAccelerationX[i] + _roadAccelerationY[i] * _roadAccelerationY[i]);

            T normedAccelerationX = _roadAccelerationX[i] / accelerationNorm;
            T normedAccelerationY = _roadAccelerationY[i] / accelerationNorm;

            // does it point to the center ?
            if ((i > 0) && (i < _numIterations) && ((std::abs(vectorToCenterX - normedAccelerationX) > tolerance) || (std::abs(vectorToCenterY - normedAccelerationY) > tolerance))) {
                std::cout << "calculateAccelerationsCenterOfGravity failed in "
                             "arbitraryTrajectory "
                             "-- wrong direction"
                          << std::endl;
            }

            // check correctness of the fl acceleration
            accelerationNorm = std::sqrt(_tyreAccelerationsX_fl[i] * _tyreAccelerationsX_fl[i] + _tyreAccelerationsY_fl[i] * _tyreAccelerationsY_fl[i]);

            T legRadius = radius - latitude;
            T legVelocity = legRadius * angularVelocity;

            // is its norm correct ?
            if ((i > 0) && (i < _numIterations) && (std::abs(legVelocity * legVelocity / legRadius - accelerationNorm) > tolerance)) {
                std::cout << "calculateAccelerationsCenterOfGravity failed in "
                             "arbitraryTrajectory "
                             "-- wrong norm of the front-left wheel"
                          << std::endl;
                std::cout << "acceleration: " << accelerationNorm << " expected: " << velocity * velocity / radius << std::endl;
            }

            // is its norm correct ?
            if ((i > 0) && (i < _numIterations) && (std::abs(velocity * velocity / radius - accelerationNorm) > tolerance)) {
                std::cout << "calculateAccelerationsCenterOfGravity failed in "
                             "arbitraryTrajectory "
                             "-- wrong norm"
                          << std::endl;
                std::cout << "acceleration: " << accelerationNorm << " expected: " << velocity * velocity / radius << std::endl;
            }

            // check correctness of the fr acceleration
            accelerationNorm = std::sqrt(_tyreAccelerationsX_fr[i] * _tyreAccelerationsX_fr[i] + _tyreAccelerationsY_fr[i] * _tyreAccelerationsY_fr[i]);

            legRadius = radius + latitude;
            legVelocity = legRadius * angularVelocity;

            // is its norm correct ?
            if ((i > 0) && (i < _numIterations) && (std::abs(legVelocity * legVelocity / legRadius - accelerationNorm) > tolerance)) {
                std::cout << "calculateAccelerationsCenterOfGravity failed in "
                             "arbitraryTrajectory "
                             "-- wrong norm of the front-left wheel"
                          << std::endl;
                std::cout << "acceleration: " << accelerationNorm << " expected: " << velocity * velocity / radius << std::endl;
            }

            // check correctness of the rl acceleration
            accelerationNorm = std::sqrt(_tyreAccelerationsX_rl[i] * _tyreAccelerationsX_rl[i] + _tyreAccelerationsY_rl[i] * _tyreAccelerationsY_rl[i]);

            legRadius = radius - latitude;
            legVelocity = legRadius * angularVelocity;

            // is its norm correct ?
            if ((i > 0) && (i < _numIterations) && (std::abs(legVelocity * legVelocity / legRadius - accelerationNorm) > tolerance)) {
                std::cout << "calculateAccelerationsCenterOfGravity failed in "
                             "arbitraryTrajectory "
                             "-- wrong norm of the front-left wheel"
                          << std::endl;
                std::cout << "acceleration: " << accelerationNorm << " expected: " << velocity * velocity / radius << std::endl;
            }

            // check correctness of the rr acceleration
            accelerationNorm = std::sqrt(_tyreAccelerationsX_rr[i] * _tyreAccelerationsX_rr[i] + _tyreAccelerationsY_rr[i] * _tyreAccelerationsY_rr[i]);

            legRadius = radius + latitude;
            legVelocity = legRadius * angularVelocity;

            // is its norm correct ?
            if ((i > 0) && (i < _numIterations) && (std::abs(legVelocity * legVelocity / legRadius - accelerationNorm) > tolerance)) {
                std::cout << "calculateAccelerationsCenterOfGravity failed in "
                             "arbitraryTrajectory "
                             "-- wrong norm of the front-left wheel"
                          << std::endl;
                std::cout << "acceleration: " << accelerationNorm << " expected: " << velocity * velocity / radius << std::endl;
            }
        }
    }

private:
    // Simulation parameters
    const size_t _numIterations;  // total number of iterations
    const T _delta_t;             // time step size of the simulation

    // Geometry parameters
    T _initialUpperSpringLength_fl;
    T _initialUpperSpringLength_fr;
    T _initialUpperSpringLength_rl;
    T _initialUpperSpringLength_rr;

    T _initialLowerSpringLength_fl;
    T _initialLowerSpringLength_fr;
    T _initialLowerSpringLength_rl;
    T _initialLowerSpringLength_rr;

    T _l_long_fl;
    T _l_long_fr;
    T _l_long_rl;
    T _l_long_rr;

    T _l_lat_fl;
    T _l_lat_fr;
    T _l_lat_rl;
    T _l_lat_rr;

    // Euler frame
    T* _legPointsZ_fl;
    T* _legPointsZ_fr;
    T* _legPointsZ_rl;
    T* _legPointsZ_rr;

    T _amplitudeRight;  // amplitude of the sinusoidal curve acting on the right
                        // tyres
    T _amplitudeLeft;   // amplitude of the sinusoidal curve acting on the left
                        // tyres

    T _frequencyRight;  // frequency of the sinusoidal curve acting on the right
                        // tyres
    T _frequencyLeft;   // frequency of the sinusoidal curve acting on the left
                        // tyres

    T _phaseShift_fl;  // phase shift of the sinusoidal curve acting on the
                       // front_left tyres
    T _phaseShift_fr;  // phase shift of the sinusoidal curve acting on the
                       // front-right tyres
    T _phaseShift_rl;  // phase shift of the sinusoidal curve acting on the
                       // rear_left tyres
    T _phaseShift_rr;  // phase shift of the sinusoidal curve acting on the
                       // rear_right tyres

    T* _tyreAccelerationsZ_fl;  // Acceleration that acts on the tyre in Z
                                // direction
    T* _tyreAccelerationsZ_fr;  // Acceleration that acts on the tyre in Z
                                // direction
    T* _tyreAccelerationsZ_rl;  // Acceleration that acts on the tyre in Z
                                // direction
    T* _tyreAccelerationsZ_rr;  // Acceleration that acts on the tyre in Z
                                // direction

    // Lagrangian frame
    T* _roadPointsX;  // X-coordinates of the trajectory at all iterations
    T* _roadPointsY;  // Y-coordinates of the trajectory at all iterations

    T* _legPointsX_fl;
    T* _legPointsX_fr;
    T* _legPointsX_rl;
    T* _legPointsX_rr;

    T* _legPointsY_fl;
    T* _legPointsY_fr;
    T* _legPointsY_rl;
    T* _legPointsY_rr;

    T* _tyreAccelerationsX_fl;  // Acceleration that acts on the tyre in X direction
    T* _tyreAccelerationsX_fr;  // Acceleration that acts on the tyre in X direction
    T* _tyreAccelerationsX_rl;  // Acceleration that acts on the tyre in X direction
    T* _tyreAccelerationsX_rr;  // Acceleration that acts on the tyre in X direction

    T* _tyreAccelerationsY_fl;  // Acceleration that acts on the tyre in Y direction
    T* _tyreAccelerationsY_fr;  // Acceleration that acts on the tyre in Y direction
    T* _tyreAccelerationsY_rl;  // Acceleration that acts on the tyre in Y direction
    T* _tyreAccelerationsY_rr;  // Acceleration that acts on the tyre in Y direction

    T* _roadAccelerationX;  // X-components of the Acceleration acting on the
                            // center of gravity at all iterations
    T* _roadAccelerationY;  // Y-components of the Acceleration acting on the
                            // center of gravity at all iterations

    T* _roadAngles;  // contains all angles

    T* _roadAngularAcceleration;  // Z-components of the AngularAcceleration
                                  // acting on the center of gravity at all
                                  // iterations

    T* _normedDistance;  // absolute distance since the beginning of the road

    /**
     * \brief use a second order scheme to calculate all Accelerations
     * between the trajectory points
     */
    void calculateAccelerations(T* acceleration, T* points) {
        T invDeltaT = 1. / (_delta_t * _delta_t);
        if (_numIterations > 3) {
            acceleration[0] = -(2 * points[0] + 5 * points[1] + 4 * points[2] + points[3]) * invDeltaT;

            for (int i = 1; i < _numIterations; ++i) {
                acceleration[i] = (points[i - 1] - 2 * points[i] + points[i + 1]) * invDeltaT;
            }

            acceleration[_numIterations] = -(2 * points[_numIterations] + 5 * points[_numIterations - 1] + 4 * points[_numIterations - 2] + points[_numIterations - 3]) * invDeltaT;
        }
        else {
            for (int i = 0; i < _numIterations + 1; ++i) {
                acceleration[i] = 0;
            }
        }
    }
};

}  // namespace EVAA
