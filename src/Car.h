// TODO: Copyright header

#pragma once

#include <mkl.h>

#include <cmath>

#include "BLAS.h"
#include "Constants.h"
#include "MathLibrary.h"
#include "MetaDataBase.h"

namespace EVAA {

template <typename T>
class Car {
private:
    /**
     * Copy all X and Y coordinates of the global vector to the local vector
     * \param Global_vector vector with coordinates [X,Y,Z,X,Y,Z,...]
     * \param local_vector vector with coordinates [X,Y,X,Y,...]
     */
    void ConstructALEVectors(T* globalVector, T* localVector)
    {
        for (size_t i = 0; i < Constants::VEC_DIM; ++i) {
            mkl<T>::copy(Constants::DIM - 1, globalVector + i * Constants::DIM, 1,
                         localVector + i * (Constants::DIM - 1), 1);
        }
    }

    /**
     * \brief compute the positions of the unexciting Euler system in the Euler frame
     */
    void ComputeUnexcitedTwoTrackPosition()
    {
        unexcitedPositionTwoTrackModel[0] = initialPositionGlobal[2];  // CG
        unexcitedPositionTwoTrackModel[1] = 0;                         // angles
        unexcitedPositionTwoTrackModel[2] = 0;
        unexcitedPositionTwoTrackModel[3] =
            initialPositionGlobal[2] - unexcitedSpringsLength[0];  // fl
        unexcitedPositionTwoTrackModel[4] =
            initialPositionGlobal[2] - unexcitedSpringsLength[0] - unexcitedSpringsLength[1];
        unexcitedPositionTwoTrackModel[5] =
            initialPositionGlobal[2] - unexcitedSpringsLength[2];  // fr
        unexcitedPositionTwoTrackModel[6] =
            initialPositionGlobal[2] - unexcitedSpringsLength[2] - unexcitedSpringsLength[3];
        unexcitedPositionTwoTrackModel[7] =
            initialPositionGlobal[2] - unexcitedSpringsLength[4];  // rl
        unexcitedPositionTwoTrackModel[8] =
            initialPositionGlobal[2] - unexcitedSpringsLength[4] - unexcitedSpringsLength[5];
        unexcitedPositionTwoTrackModel[9] =
            initialPositionGlobal[2] - unexcitedSpringsLength[6];  // rr
        unexcitedPositionTwoTrackModel[10] =
            initialPositionGlobal[2] - unexcitedSpringsLength[6] - unexcitedSpringsLength[7];
    }

    /**
     * \brief Initial displacements as required for the 11DOF
     */
    void ComputeInitialDisplacement()
    {
        currentDisplacementTwoTrackModel[0] = 0;
        currentDisplacementTwoTrackModel[1] = initialAngleGlobal[0];
        currentDisplacementTwoTrackModel[2] = initialAngleGlobal[1];
        currentDisplacementTwoTrackModel[3] = unexcitedSpringsLength[0] - currentSpringsLength[0] +
                                              currentCornerPositions[8] - initialPositionGlobal[2];
        currentDisplacementTwoTrackModel[4] = currentDisplacementTwoTrackModel[3] +
                                              unexcitedSpringsLength[1] - currentSpringsLength[1];
        currentDisplacementTwoTrackModel[5] = unexcitedSpringsLength[2] - currentSpringsLength[2] +
                                              currentCornerPositions[9] - initialPositionGlobal[2];
        currentDisplacementTwoTrackModel[6] = currentDisplacementTwoTrackModel[5] +
                                              unexcitedSpringsLength[3] - currentSpringsLength[3];
        currentDisplacementTwoTrackModel[7] = unexcitedSpringsLength[4] - currentSpringsLength[4] +
                                              currentCornerPositions[10] - initialPositionGlobal[2];
        currentDisplacementTwoTrackModel[8] = currentDisplacementTwoTrackModel[7] +
                                              unexcitedSpringsLength[5] - currentSpringsLength[5];
        currentDisplacementTwoTrackModel[9] = unexcitedSpringsLength[6] - currentSpringsLength[6] +
                                              currentCornerPositions[11] - initialPositionGlobal[2];
        currentDisplacementTwoTrackModel[10] = currentDisplacementTwoTrackModel[9] +
                                               unexcitedSpringsLength[7] - currentSpringsLength[7];
    }

    /**
     * Construct corner initilizer
     */
    void ConstructCornerRelativeToCG(T* corners)
    {
        corners[0] = +l_long[0];  // fl
        corners[4] = +l_lat[0];   // fl

        corners[1] = +l_long[1];  // fr
        corners[5] = -l_lat[1];   // fr

        corners[2] = -l_long[2];  // rl
        corners[6] = +l_lat[2];   // rl

        corners[3] = -l_long[3];  // rr
        corners[7] = -l_lat[3];   // rr
    }

    /**
     * Calculates the values of Corners for general angles.
     * TODO: Consider moving to MathLibrary or anonymous namespace.
     */
    void UpdateCorners11DOF(T* angles, T* rotation_mat_buffer, T* initial_corners,
                            T* updated_corners)
    {
        // zz, yy, xx
        MathLibrary::get_rotation_matrix(angles[2], angles[1], angles[0], rotation_mat_buffer);

        // do rotation: rotationMat * r
        // void cblas_dgemm(const CBLAS_LAYOUT Layout, const CBLAS_TRANSPOSE transa, const
        // CBLAS_TRANSPOSE transb, const MKL_INT m, const MKL_INT n, const MKL_INT k, const double
        // alpha, const double* a, const MKL_INT lda, const double* b, const MKL_INT ldb, const
        // double beta, double* c, const MKL_INT ldc);
        mkl<T>::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, Constants::DIM, Constants::NUM_LEGS,
                     Constants::DIM, 1, rotation_mat_buffer, Constants::DIM, initial_corners,
                     Constants::NUM_LEGS, 0, updated_corners, Constants::NUM_LEGS);
    }

    /**
     * \brief This function computes conter position about CG
     * \param corner: Array of corners about origin in column major format
     * \param center: Array of coordinate of center of mass
     */
    void CornerAboutCenter(T* corner, T* center)
    {
        for (size_t i = 0; i < Constants::NUM_LEGS; ++i) {
            mkl<T>::axpy(Constants::DIM, 1.0, center, 1, corner + i, Constants::NUM_LEGS);
        }
    }

    /**
     * Calculates the values of currentCornerPositions according to the current orientation
     */
    void UpdateCorners11DOF()
    {
        pos_buffer[0] = currentPositionLagrangian[0];
        pos_buffer[1] = currentPositionLagrangian[1];
        pos_buffer[2] = unexcitedPositionTwoTrackModel[0] + currentDisplacementTwoTrackModel[0];
        angle_buffer[0] = unexcitedPositionTwoTrackModel[1] + currentDisplacementTwoTrackModel[1];
        angle_buffer[1] = unexcitedPositionTwoTrackModel[2] + currentDisplacementTwoTrackModel[2];
        angle_buffer[2] = *currentAngleLagrangian;

        ConstructCornerRelativeToCG(relativeCornerPositions);
        UpdateCorners11DOF(angle_buffer, currentRotationMatrix, relativeCornerPositions,
                           currentCornerPositions);
        CornerAboutCenter(currentCornerPositions, pos_buffer);
    }

    // TODO: Was ist das? Maybe to MathLibrary
    void ConvertALEToGlobal(T* vect, T* global_vect)
    {
#pragma loop(ivdep)
        for (size_t i = 0; i < Constants::VEC_DIM; ++i) {
            global_vect[Constants::DIM * i] = vect[(Constants::DIM - 1) * i];
            global_vect[Constants::DIM * i + 1] = vect[(Constants::DIM - 1) * i + 1];
        }
    }

    // TODO: Was ist das? Maybe to MathLibrary
    void Convert11DOFToGlobal(T* vect, T* global_vect)
    {
        global_vect[2] += vect[0];
#pragma loop(ivdep)
        for (size_t i = 1; i < Constants::VEC_DIM; ++i) {
            global_vect[Constants::DIM * i + 2] += vect[(Constants::DIM - 1) + i];
        }
    }

public:
    /*
    Using Stefan's convention
    W1 = front left wheel
    W2 = front right wheel
    W3 = rear left wheel
    W4 = rear right wheel
    T1 = front left tyre
    T2 = front right tyre
    T3 = rear left tyre
    T4 = rear right tyre
    */

    // physical parameters
    T* massComponents;   // [CG,  W_fl, T_fl, W_fr, T_fr, W_rl, T_rl, W_rr, T_rr]    9 x 1
    T* massFullCar;      // 1 scalar
    T* momentOfInertia;  // [Ixx, Ixy, Ixz, Iyx, Iyy, Iyz, Izx, Izy, Izz]
    T* vehicleCIR;       // [XYZ]
    T* unexcitedPositionTwoTrackModel;  // [CG:Z, theta:XY, W_fl:Z, T_fl:Z, W_fr:Z, T_fr:Z, W_rl:Z,
                                        // T_rl:Z, W_rr:Z, T_rr:Z]
    T* unexcitedSpringsLength;  // [W_fl:Z, T_fl:Z, W_fr:Z, T_fr:Z, W_rl:Z, T_rl:Z, W_rr:Z, T_rr:Z]

    // only use to write to HDF5 DISCUSS
    T* Position_vec;  // [CG:XYZ, W_fl:XYZ, T_fl:XYZ, W_fr:XYZ, T_fr:XYZ, W_rl:XYZ, T_rl:XYZ,
                      // W_rr:XYZ, T_rr:XYZ] n x 9 x 3 !!! global Consider Constants::ALIGNMENT
                      // (3+1),(3+1),...
    T* Velocity_vec;  // [CG:XYZ, W_fl:XYZ, T_fl:XYZ, W_fr:XYZ, T_fr:XYZ, W_rl:XYZ, T_rl:XYZ,
                      // W_rr:XYZ, T_rr:XYZ] n x 9 x 3 !!! global
    T* angle_CG;      // [XYZ]
    T* w_CG;          // [XYZ]

    // Initial Conditions of the car
    T* initialPositionGlobal;         // [CG:XYZ, W_fl:XYZ, T_fl:XYZ, W_fr:XYZ, T_fr:XYZ, W_rl:XYZ,
                                      // T_rl:XYZ, W_rr:XYZ, T_rr:XYZ] !!! global
    T* initialVelocityGlobal;         // [CG:XYZ, W_fl:XYZ, T_fl:XYZ, W_fr:XYZ, T_fr:XYZ, W_rl:XYZ,
                                      // T_rl:XYZ, W_rr:XYZ, T_rr:XYZ] !!! global
    T* initialAngleGlobal;            // [XYZ]
    T* initialAngularVelocityGlobal;  // [XYZ]
    T* relativeCornerPositions;       // [X:fl,fr,rl,rr Y:fl,fr,rl,rr Z:fl,fr,rl,rr]

    // vectors to use in regular iteration
    // For 11DOF
    T* currentVelocityTwoTrackModel;  // [CG:Z, theta:XY, W_fl:Z, T_fl:Z, W_fr:Z, T_fr:Z, W_rl:Z,
                                      // T_rl:Z, W_rr:Z, T_rr:Z]
    T* currentPositionTwoTrackModel;  // [CG:Z, theta:XY, W_fl:Z, T_fl:Z, W_fr:Z, T_fr:Z, W_rl:Z,
                                      // T_rl:Z, W_rr:Z, T_rr:Z]
    T* currentDisplacementTwoTrackModel;  // [CG:Z, theta:XY, W_fl:Z, T_fl:Z, W_fr:Z, T_fr:Z,
                                          // W_rl:Z, T_rl:Z, W_rr:Z, T_rr:Z]
    T* currentCIRTwoTrackModel;  // [CG:Z, W_fl:Z, T_fl:Z, W_fr:Z, T_fr:Z, W_rl:Z, T_rl:Z, W_rr:Z,
                                 // T_rr:Z]
    T* currentSpringsLength;     // [W_fl:Z, T_fl:Z, W_fr:Z, T_fr:Z, W_rl:Z, T_rl:Z, W_rr:Z, T_rr:Z]

    // For ALE || suggestion: reduce it to only leg position (wheel == tyre)
    T* currentPositionLagrangian;  // [CG:XY, W_fl:XY, T_fl:XY, W_fr:XY, T_fr:XY, W_rl:XY, T_rl:XY,
                                   // W_rr:XY, T_rr:XY]
    T* currentAngleLagrangian;     // [CG:Z] -> change to T (pass by reference)
    T* currentVelocityLagrangian;  // [CG:XY, W_fl:XY, T_fl:XY, W_fr:XY, T_fr:XY, W_rl:XY, T_rl:XY,
                                   // W_rr:XY, T_rr:XY]
    T* currentAngularVelocityLagrangian;  // [CG:Z] -> change to T (pass by reference)

    // For global necessary updates
    T* currentRotationMatrix;   // [xx, xy, xz, yx, yy, yz, zx, zy, zz]
    T* currentCornerPositions;  // [X:fl,fr,rl,rr Y:fl,fr,rl,rr Z:fl,fr,rl,rr]         sorry for
                                // everyone

    EVAALookup<T>*_lookupStiffness, *_lookupDamping;

    // Members from 11 DOF system

    T *kVec, *dVec, *l_lat, *l_long;  // to be TEO
    size_t* tyre_index_set;           // to be something

    // Interpolator Members

    T *angle_buffer, *pos_buffer;  // to be removed

    /* Constructor */
    Car(EVAALookup<T>* lookupStiffness, EVAALookup<T>* lookupDamping)
    {
        _lookupStiffness = lookupStiffness;
        _lookupDamping = lookupDamping;

        // System Mono

        // Memory Allocation and matrix formulation

        const int positionAllocSize =
            Constants::DIM * Constants::VEC_DIM * sizeof(T);        // 27 dimensions
        const int pointAllocSize = Constants::VEC_DIM * sizeof(T);  // 9 dimensions
        const int dimAllocSize = Constants::DIM * sizeof(T);        // 3 dimensions

        // Params for Global coordinate

        Position_vec = (T*)mkl_malloc(positionAllocSize, Constants::ALIGNMENT);
        Velocity_vec = (T*)mkl_malloc(positionAllocSize, Constants::ALIGNMENT);
        massComponents = (T*)mkl_malloc(pointAllocSize, Constants::ALIGNMENT);
        angle_CG = (T*)mkl_malloc(dimAllocSize, Constants::ALIGNMENT);
        w_CG = (T*)mkl_malloc(dimAllocSize, Constants::ALIGNMENT);
        momentOfInertia = (T*)mkl_malloc(Constants::DIM * Constants::DIM * sizeof(T),
                                         Constants::ALIGNMENT);  // 9 Dimensions

        // Initial Params for Global coordinate

        initialPositionGlobal = (T*)mkl_malloc(positionAllocSize, Constants::ALIGNMENT);
        initialVelocityGlobal = (T*)mkl_malloc(positionAllocSize, Constants::ALIGNMENT);
        initialAngleGlobal =
            (T*)mkl_malloc(dimAllocSize, Constants::ALIGNMENT);  // 3 Constants::DIM
        initialAngularVelocityGlobal =
            (T*)mkl_malloc(dimAllocSize, Constants::ALIGNMENT);  // 3 Constants::DIM

        // Params for 11 DOF system

        currentDisplacementTwoTrackModel =
            (T*)mkl_malloc(Constants::DOF * sizeof(T), Constants::ALIGNMENT);  // 11 Constants::DIM
        currentVelocityTwoTrackModel =
            (T*)mkl_malloc(Constants::DOF * sizeof(T), Constants::ALIGNMENT);  // 3 Constants::DIM
        kVec = (T*)mkl_malloc(Constants::NUM_LEGS * 2 * sizeof(T),
                              Constants::ALIGNMENT);  // 4 Constants::DIM
        dVec = (T*)mkl_calloc(Constants::NUM_LEGS * 2, sizeof(T),
                              Constants::ALIGNMENT);  // 4 Constants::DIM
        l_lat = (T*)mkl_malloc(Constants::NUM_LEGS * sizeof(T),
                               Constants::ALIGNMENT);  // 4 Constants::DIM
        l_long = (T*)mkl_malloc(Constants::NUM_LEGS * sizeof(T),
                                Constants::ALIGNMENT);  // 4 Constants::DIM
        unexcitedSpringsLength = (T*)mkl_malloc(2 * Constants::NUM_LEGS * sizeof(T),
                                                Constants::ALIGNMENT);  // 8 Constants::DIM
        currentCIRTwoTrackModel =
            (T*)mkl_malloc((2 * Constants::NUM_LEGS + 1) * sizeof(T), Constants::ALIGNMENT);  // 9
        currentSpringsLength = (T*)mkl_malloc(2 * Constants::NUM_LEGS * sizeof(T),
                                              Constants::ALIGNMENT);  // 8 Constants::DIM
        tyre_index_set =
            (size_t*)mkl_malloc(Constants::NUM_LEGS * sizeof(size_t), Constants::ALIGNMENT);
        currentPositionTwoTrackModel =
            (T*)mkl_malloc(Constants::DOF * sizeof(T), Constants::ALIGNMENT);  // 11 Constants::DOF
        unexcitedPositionTwoTrackModel =
            (T*)mkl_malloc(Constants::DOF * sizeof(T), Constants::ALIGNMENT);  // 11 Constants::DOF

        // Memory allocation for interpolator

        currentCornerPositions = (T*)mkl_malloc(Constants::NUM_LEGS * dimAllocSize,
                                                Constants::ALIGNMENT);  // 12 Constants::DIM
        currentRotationMatrix = (T*)mkl_malloc(Constants::DIM * dimAllocSize,
                                               Constants::ALIGNMENT);  // 9 Constants::DIM
        relativeCornerPositions = (T*)mkl_malloc(Constants::NUM_LEGS * dimAllocSize,
                                                 Constants::ALIGNMENT);     // 12 Constants::DIM
        angle_buffer = (T*)mkl_malloc(dimAllocSize, Constants::ALIGNMENT);  // 3 Constants::DIM
        pos_buffer = (T*)mkl_malloc(dimAllocSize, Constants::ALIGNMENT);    // 3 Constants::DIM

        // ALE Buffer Allocation

        currentPositionLagrangian =
            (T*)mkl_malloc((Constants::DIM - 1) * pointAllocSize, Constants::ALIGNMENT);
        currentAngleLagrangian = new T;
        currentVelocityLagrangian =
            (T*)mkl_malloc((Constants::DIM - 1) * pointAllocSize, Constants::ALIGNMENT);
        currentAngularVelocityLagrangian = new T;
        massFullCar = new T;

        // Extract Data from parser

        mkl<T>::copy(Constants::NUM_LEGS,
                     MetaDataBase::DataBase()->getLongitudalLegPositionVector(), 1, l_long, 1);

        mkl<T>::copy(Constants::NUM_LEGS, MetaDataBase::DataBase()->getLatidudalLegPositionVector(),
                     1, l_lat, 1);

        mkl<T>::copy(Constants::DIM * Constants::DIM,
                     MetaDataBase::DataBase()->getMomentOfInertiaVector(), 1, momentOfInertia, 1);

        massComponents[0] = MetaDataBase::DataBase()->getBodyMass();
        massComponents[1] = MetaDataBase::DataBase()->getWheelMassFrontLeft();
        massComponents[2] = MetaDataBase::DataBase()->getTyreMassFrontLeft();
        massComponents[3] = MetaDataBase::DataBase()->getWheelMassFrontRight();
        massComponents[4] = MetaDataBase::DataBase()->getTyreMassFrontRight();
        massComponents[5] = MetaDataBase::DataBase()->getWheelMassRearLeft();
        massComponents[6] = MetaDataBase::DataBase()->getTyreMassRearLeft();
        massComponents[7] = MetaDataBase::DataBase()->getWheelMassRearRight();
        massComponents[8] = MetaDataBase::DataBase()->getTyreMassRearRight();
        *massFullCar = getMassFullCar();

        unexcitedSpringsLength[0] = MetaDataBase::DataBase()->getBodySpringLengthFrontLeft();
        unexcitedSpringsLength[1] = MetaDataBase::DataBase()->getTyreSpringLengthFrontLeft();
        unexcitedSpringsLength[2] = MetaDataBase::DataBase()->getBodySpringLengthFrontRight();
        unexcitedSpringsLength[3] = MetaDataBase::DataBase()->getTyreSpringLengthFrontRight();
        unexcitedSpringsLength[4] = MetaDataBase::DataBase()->getBodySpringLengthRearLeft();
        unexcitedSpringsLength[5] = MetaDataBase::DataBase()->getTyreSpringLengthRearLeft();
        unexcitedSpringsLength[6] = MetaDataBase::DataBase()->getBodySpringLengthRearRight();
        unexcitedSpringsLength[7] = MetaDataBase::DataBase()->getTyreSpringLengthRearRight();

        vehicleCIR = MetaDataBase::DataBase()->getPositionCenterOfInstantaneousRotation();

        // Initial Iteration vector

        // Initial Angles
        MathLibrary::ToEulerAngles<T>(MetaDataBase::DataBase()->getBodyInitialOrientation(),
                                      initialAngleGlobal);
        mkl<T>::copy(Constants::DIM, initialAngleGlobal, 1, angle_CG, 1);

        // Spring lengths
        currentSpringsLength[0] = MetaDataBase::DataBase()->getBodySpringInitialLengthFrontLeft();
        currentSpringsLength[1] = MetaDataBase::DataBase()->getTyreSpringInitialLengthFrontLeft();
        currentSpringsLength[2] = MetaDataBase::DataBase()->getBodySpringInitialLengthFrontRight();
        currentSpringsLength[3] = MetaDataBase::DataBase()->getTyreSpringInitialLengthFrontRight();
        currentSpringsLength[4] = MetaDataBase::DataBase()->getBodySpringInitialLengthRearLeft();
        currentSpringsLength[5] = MetaDataBase::DataBase()->getTyreSpringInitialLengthRearLeft();
        currentSpringsLength[6] = MetaDataBase::DataBase()->getBodySpringInitialLengthRearRight();
        currentSpringsLength[7] = MetaDataBase::DataBase()->getTyreSpringInitialLengthRearRight();

        // Filling the position vector with initial condition
        // CG
        mkl<T>::copy(Constants::DIM, MetaDataBase::DataBase()->getBodyInitialPosition(), 1,
                     initialPositionGlobal, 1);  // copy the center of mass position

        // Interpolator
        // Initialization: read init corners vectors into matrix

        ConstructCornerRelativeToCG(
            relativeCornerPositions);  // only CG position is used to construct corners
        mkl<T>::copy(Constants::DIM, initialAngleGlobal, 1, angle_buffer, 1);
        // angle_buffer[2] = 0;
        UpdateCorners11DOF(angle_buffer, currentRotationMatrix, relativeCornerPositions,
                           currentCornerPositions);
        CornerAboutCenter(currentCornerPositions, initialPositionGlobal);

        // Remaining position initialization

        const T* xml_start;
        T* position_start;
        if (MetaDataBase::DataBase()->getFlagInitialLeg()) {  // DEBUG to  be removed
            // if prescribed initial position (add a check for consistency with spring lengths)
            // W1 = W_fl
            xml_start = MetaDataBase::DataBase()->getWheelInitialPositionFrontLeft();
            position_start = initialPositionGlobal + 3;  //(end at 5)
            mkl<T>::copy(Constants::DIM, xml_start, 1, position_start, 1);
            // W2 = W_fr
            xml_start = MetaDataBase::DataBase()->getWheelInitialPositionFrontRight();
            position_start += 6;  // skip 3 for tyre (end at 11)
            mkl<T>::copy(Constants::DIM, xml_start, 1, position_start, 1);
            // W3 = W_rl
            xml_start = MetaDataBase::DataBase()->getWheelInitialPositionRearLeft();
            position_start += 6;  // skip 3 for tyre (end at 17)
            mkl<T>::copy(Constants::DIM, xml_start, 1, position_start, 1);
            // W2 = W_rr
            xml_start = MetaDataBase::DataBase()->getWheelInitialPositionRearRight();
            position_start += 6;  // skip 3 for tyre (end at 23)
            mkl<T>::copy(Constants::DIM, xml_start, 1, position_start, 1);

            // T1 = T_fl
            xml_start = MetaDataBase::DataBase()->getTyreInitialPositionFrontLeft();
            position_start =
                initialPositionGlobal + 6;  // skip 3 for center of mass and 3 for the wheel
            mkl<T>::copy(Constants::DIM, xml_start, 1, position_start, 1);  // (end at 8)
            // T2 = T_fr
            xml_start = MetaDataBase::DataBase()->getTyreInitialPositionFrontRight();
            position_start += 6;                                            // skip 3 for the wheel
            mkl<T>::copy(Constants::DIM, xml_start, 1, position_start, 1);  // (end at 14)
            // T3 = T_rl
            xml_start = MetaDataBase::DataBase()->getTyreInitialPositionRearLeft();
            position_start += 6;                                            // skip 3 for the wheel
            mkl<T>::copy(Constants::DIM, xml_start, 1, position_start, 1);  // (end at 20)
            // T4 = T_rr
            xml_start = MetaDataBase::DataBase()->getTyreInitialPositionRearRight();
            position_start += 6;                                            // skip 3 for the wheel
            mkl<T>::copy(Constants::DIM, xml_start, 1, position_start, 1);  // (end at 26)
            mkl<T>::copy(Constants::DIM * Constants::VEC_DIM, initialPositionGlobal, 1,
                         Position_vec, 1);
        }
        else {
            T* W_fl = initialPositionGlobal + 3;
            T* W_fr = initialPositionGlobal + 9;
            T* W_rl = initialPositionGlobal + 15;
            T* W_rr = initialPositionGlobal + 21;
            T* T_fl = initialPositionGlobal + 6;
            T* T_fr = initialPositionGlobal + 12;
            T* T_rl = initialPositionGlobal + 18;
            T* T_rr = initialPositionGlobal + 24;
            computeGlobalPositionWheelTyre(currentCornerPositions, currentSpringsLength, W_fl, T_fl,
                                           W_fr, T_fr, W_rl, T_rl, W_rr, T_rr);
            /* Update the mean position where changes are to be added */
            mkl<T>::copy(Constants::DIM * Constants::VEC_DIM, initialPositionGlobal, 1,
                         Position_vec, 1);
        }

        // Initial Velocity (Reuse the pointers)
        mkl<T>::copy(Constants::DIM, MetaDataBase::DataBase()->getBodyInitialVelocity(), 1,
                     initialVelocityGlobal, 1);
        // W1 = W_fl
        xml_start = MetaDataBase::DataBase()->getWheelInitialVelocityFrontLeft();
        position_start = initialVelocityGlobal + 3;
        mkl<T>::copy(Constants::DIM, xml_start, 1, position_start, 1);  // (end at 5)
        // W2 = W_fr
        xml_start = MetaDataBase::DataBase()->getWheelInitialVelocityFrontRight();
        position_start += 6;                                            // skip 3 for tyre
        mkl<T>::copy(Constants::DIM, xml_start, 1, position_start, 1);  // (end at 11)
        // W3 = W_rl
        xml_start = MetaDataBase::DataBase()->getWheelInitialVelocityRearLeft();
        position_start += 6;                                            // skip 3 for tyre
        mkl<T>::copy(Constants::DIM, xml_start, 1, position_start, 1);  // (end at 17)
        // W2 = W_rr
        xml_start = MetaDataBase::DataBase()->getWheelInitialVelocityRearRight();
        position_start += 6;                                            // skip 3 for tyre
        mkl<T>::copy(Constants::DIM, xml_start, 1, position_start, 1);  // (end at 23)

        // T1 = T_fl
        xml_start = MetaDataBase::DataBase()->getTyreInitialVelocityFrontLeft();
        position_start =
            initialVelocityGlobal + 6;  // skip 3 for center of mass and 3 for the wheel
        mkl<T>::copy(Constants::DIM, xml_start, 1, position_start, 1);  // (end at 8)

        // T2 = T_fr
        xml_start = MetaDataBase::DataBase()->getTyreInitialVelocityFrontRight();
        position_start += 6;                                            // skip 3 for the Tyre
        mkl<T>::copy(Constants::DIM, xml_start, 1, position_start, 1);  // (end at 14)
        // T3 = T_rl
        xml_start = MetaDataBase::DataBase()->getTyreInitialVelocityRearLeft();
        position_start += 6;                                            // skip 3 for the wheel
        mkl<T>::copy(Constants::DIM, xml_start, 1, position_start, 1);  // (end at 20)

        // T4 = T_rr
        xml_start = MetaDataBase::DataBase()->getTyreInitialVelocityRearRight();
        position_start += 6;                                            // skip 3 for the wheel
        mkl<T>::copy(Constants::DIM, xml_start, 1, position_start, 1);  // (end at 26)

        // copy the initial position to the position vector
        mkl<T>::copy(Constants::DIM * Constants::VEC_DIM, initialVelocityGlobal, 1, Velocity_vec,
                     1);

        // Initial Angular velocity
        mkl<T>::copy(Constants::DIM, MetaDataBase::DataBase()->getBodyInitialAngularVelocity(), 1,
                     initialAngularVelocityGlobal, 1);
        mkl<T>::copy(Constants::DIM, initialAngularVelocityGlobal, 1, w_CG, 1);

        /* Global assignments Done */

        // 11 DOF Buffer Initialization

        tyre_index_set[0] = 2;
        tyre_index_set[1] = 4;
        tyre_index_set[2] = 6;
        tyre_index_set[3] = 8;

        // unexcited state
        ComputeUnexcitedTwoTrackPosition();

        // solution vector of 11DOF
        ComputeInitialDisplacement();

        // initial velocity
        construct_11DOF_vector(initialVelocityGlobal, initialAngularVelocityGlobal,
                               currentVelocityTwoTrackModel);

        // initial position
        construct_11DOF_vector(initialPositionGlobal, initialAngleGlobal,
                               currentPositionTwoTrackModel);

        // we updated the initial global position, so we
        updateRadiusToCIR();

#ifdef INTERPOLATION
        _lookupStiffness->getInterpolation(currentSpringsLength, kVec);
//      _lookupDamping->getInterpolation(currentSpringsLength, dVec);
#else
        kVec[0] = MetaDataBase::DataBase()->getBodyStiffnessFrontLeft();
        kVec[1] = MetaDataBase::DataBase()->getTyreStiffnessFrontLeft();
        kVec[2] = MetaDataBase::DataBase()->getBodyStiffnessFrontRight();
        kVec[3] = MetaDataBase::DataBase()->getTyreStiffnessFrontRight();
        kVec[4] = MetaDataBase::DataBase()->getBodyStiffnessRearLeft();
        kVec[5] = MetaDataBase::DataBase()->getTyreStiffnessRearLeft();
        kVec[6] = MetaDataBase::DataBase()->getBodyStiffnessRearRight();
        kVec[7] = MetaDataBase::DataBase()->getTyreStiffnessRearRight();

        dVec[0] = MetaDataBase::DataBase()->getBodyDampingFrontLeft();
        dVec[1] = MetaDataBase::DataBase()->getTyreDampingFrontLeft();
        dVec[2] = MetaDataBase::DataBase()->getBodyDampingFrontRight();
        dVec[3] = MetaDataBase::DataBase()->getTyreDampingFrontRight();
        dVec[4] = MetaDataBase::DataBase()->getBodyDampingRearLeft();
        dVec[5] = MetaDataBase::DataBase()->getTyreDampingRearLeft();
        dVec[6] = MetaDataBase::DataBase()->getBodyDampingRearRight();
        dVec[7] = MetaDataBase::DataBase()->getTyreDampingRearRight();
#endif

        // ALE Buffer Initialization

        ConstructALEVectors(initialPositionGlobal, currentPositionLagrangian);
        ConstructALEVectors(initialVelocityGlobal, currentVelocityLagrangian);
        *currentAngleLagrangian = initialAngleGlobal[2];
        *currentAngularVelocityLagrangian = initialAngularVelocityGlobal[2];
    }

    /**
     * If forces and positiosn are negative, set them to zero, elsewise, keep them
     */
    void apply_normal_force(T* force, T* u, size_t* index, size_t n)
    {
#pragma loop(ivdep)
        for (int i = 0; i < n; ++i) {
            force[index[i]] = force[index[i]] > 0 ? force[index[i]] : 0;
        }
#pragma loop(ivdep)
        for (int i = 0; i < n; ++i) {
            u[index[i]] = u[index[i]] > 0 ? u[index[i]] : 0;
        }
    }

    /**
     * compute a reaction which is opposite to the internal force acting on the tyre
     */
    void compute_normal_force(T* K, T* u, T* force, size_t* index, size_t n)
    {
#pragma loop(ivdep)
        for (int i = 0; i < n; ++i) {
            force[index[i]] = -K[index[i] * Constants::DIM + index[i]] * u[index[i]];
        }
    }

    /**
     * Get the solution vector as required for the 11DOF system
     * \param Global_position in the format [GC:XYZ,W_fl:XYZ,T_fl:XYZ,W_fr:XYZ,T_fr:XYZ,...]
     * \param Global_angle with three angles [X,Y,Z]
     * \return Position_11dof in the format [GC:Y,angle:XY,W_fl:Y,T_fl:Y,W_fr:Y,T_fr:Y,...]
     */
    void construct_11DOF_vector(T* Global_position, T* Global_angle, T* Position_11dof)
    {
        Position_11dof[0] = Global_position[2];  // z coordinate of CG
        Position_11dof[1] = Global_angle[0];     // x angle of the CG
        Position_11dof[2] = Global_angle[1];     // y angle of the CG

        // copy y coordinate in order wheel, tyre, wheel, tyre, wheel, tyre, ...
        mkl<T>::copy(Constants::VEC_DIM - 1, Global_position + 5, 3, Position_11dof + 3, 1);
    }

    void computeGlobalPositionWheelTyre(T* Corners, T* curr_spring_len, T* W_fl, T* T_fl, T* W_fr,
                                        T* T_fr, T* W_rl, T* T_rl, T* W_rr, T* T_rr)
    {
        // W_fl & T_fl
        mkl<T>::copy(Constants::DIM, Corners, 4, W_fl, 1);
        W_fl[2] -= curr_spring_len[0];
        mkl<T>::copy(Constants::DIM, W_fl, 1, T_fl, 1);
        T_fl[2] -= curr_spring_len[1];

        // W_fr & T_fr
        mkl<T>::copy(Constants::DIM, Corners + 1, 4, W_fr, 1);
        W_fr[2] -= curr_spring_len[2];
        mkl<T>::copy(Constants::DIM, W_fr, 1, T_fr, 1);
        T_fr[2] -= curr_spring_len[3];

        // W_rl & T_rl
        mkl<T>::copy(Constants::DIM, Corners + 2, 4, W_rl, 1);
        W_rl[2] -= curr_spring_len[4];
        mkl<T>::copy(Constants::DIM, W_rl, 1, T_rl, 1);
        T_rl[2] -= curr_spring_len[5];

        // W_rr & T_rr
        mkl<T>::copy(Constants::DIM, Corners + 3, 4, W_rr, 1);
        W_rr[2] -= curr_spring_len[6];
        mkl<T>::copy(Constants::DIM, W_rr, 1, T_rr, 1);
        T_rr[2] -= curr_spring_len[7];
    }

    inline void compute_dx(const T* current_length, T* dx)
    {
        /*
        the dx follows the order
        [w1, t1, w2, t2, w3, t3, w4, t4]
        current_length has input of type
        [w1, t1, w2, t2, w3, t3, w4, t4]
        */
        mkl<T>::vSub(2 * (Constants::NUM_LEGS), unexcitedSpringsLength, current_length, dx);
    }

    /**
     * \brief From the current elongations, calculate the differences to the rest positions
     * \return length differences
     */
    inline void compute_dx(T* dx) { compute_dx(currentSpringsLength, dx); }

    inline void compute_dx_tyre(T* dx)
    {
        mkl<T>::copy(Constants::NUM_LEGS, currentDisplacementTwoTrackModel + 4, 2, dx, 1);
    }

    /**
     * \*brief compute current spring lengths
     */
    void computeCurrentSpringLength()
    {
        // upper spring_length = corner - wheel
        // lower spring length = wheel - tyre
        mkl<T>::copy(Constants::NUM_LEGS, currentCornerPositions + 8, 1, currentSpringsLength, 2);
        mkl<T>::axpy(Constants::NUM_LEGS, -1, currentPositionTwoTrackModel + 3, 2,
                     currentSpringsLength, 2);
        mkl<T>::copy(Constants::NUM_LEGS, currentPositionTwoTrackModel + 3, 2,
                     currentSpringsLength + 1, 2);
        mkl<T>::axpy(Constants::NUM_LEGS, -1, currentPositionTwoTrackModel + 4, 2,
                     currentSpringsLength + 1, 2);
    }

    /**
     * \brief updates the global position of the components based on euler frame displacement
     */
    void updateGlobalTwoTrackVectors()
    {
        // currentPositionTwoTrackModel = unexcitedPositionTwoTrackModel +
        // currentDisplacementTwoTrackModel;
        mkl<T>::copy(Constants::DOF, unexcitedPositionTwoTrackModel, 1,
                     currentPositionTwoTrackModel, 1);
        mkl<T>::axpy(Constants::DOF, +1, currentDisplacementTwoTrackModel, 1,
                     currentPositionTwoTrackModel, 1);
    }

    /**
     * \brief calculate the Z-distance from each car component to the CIR
     */
    void updateRadiusToCIR()
    {
        T globalNickpolPosition = currentPositionTwoTrackModel[0] +
                                  vehicleCIR[2];      // get global Z position of the Nickpol
        currentCIRTwoTrackModel[0] = -vehicleCIR[2];  // initialize in the constructor
        for (int i = 1; i < 2 * Constants::NUM_LEGS + 1; ++i) {
            currentCIRTwoTrackModel[i] =
                currentPositionTwoTrackModel[2 + i] - globalNickpolPosition;
        }
    }

    /**
     * \brief update global corner positions, compute global Euler vectors, compute current spring
     * lengths, compute CIR functions
     */
    void updateLengthsTwoTrackModel()
    {
        UpdateCorners11DOF();
        // global vector update
        updateGlobalTwoTrackVectors();

        // compute current spring lengths
        computeCurrentSpringLength();

        // Compute CIR distances
        updateRadiusToCIR();
    }

    /**
     * Fills the global vector with all entries (checkpointing)
     * \param ALE_vectors contains X and Y components [GC:XY,W_fl:XY,T_fl:XY,W_fr:XY,T_fr:XY,...]
     * \param vector 11DOF contains Z components [GC:Z,W_fl:Z,T_fl:Z,W_fr:Z,T_fr:Z,...]
     * \return global_vector [GC:XYZ,W_fl:XYZ,T_fl:XYZ,W_fr:XYZ,T_fr:XYZ,...]
     */
    void combineEulerianLagrangianVectors(T* ALE_vector, T* vector_11DOF, T* global_vector)
    {
        ConvertALEToGlobal(ALE_vector, global_vector);
        Convert11DOFToGlobal(vector_11DOF, global_vector);
    }

    /**
     * \brief flush the result to the output (checkpointing) || combine with function above
     */
    void combine_results()
    {
        mkl<T>::copy((Constants::VEC_DIM) * (Constants::DIM), initialPositionGlobal, 1,
                     Position_vec, 1);
        ConvertALEToGlobal(currentPositionLagrangian, Position_vec);
        Convert11DOFToGlobal(currentDisplacementTwoTrackModel, Position_vec);

        // Angles manually
        angle_CG[0] = initialAngleGlobal[0] + currentDisplacementTwoTrackModel[1];
        angle_CG[1] = initialAngleGlobal[1] + currentDisplacementTwoTrackModel[2];
        angle_CG[2] = *currentAngleLagrangian;
        w_CG[2] = *currentAngularVelocityLagrangian;
        ConvertALEToGlobal(currentVelocityLagrangian, Velocity_vec);
    }

    // Sums up all the 9 masses
    inline T getMassFullCar()
    {
        return (massComponents[0] +  // CG
                massComponents[1] + massComponents[2] + massComponents[3] + massComponents[4] +
                massComponents[5] + massComponents[6] + massComponents[7] + massComponents[8]);
    }

    void get_Velocity_vec_xy(T* Vel)
    {
        mkl<T>::copy((Constants::DIM - 1) * Constants::VEC_DIM, currentVelocityLagrangian, 1, Vel,
                     1);
    }

    void get_k_vec(T* k) { mkl<T>::copy(2 * Constants::NUM_LEGS, kVec, 1, k, 1); }
    void get_k_vec_tyre(T* k) { mkl<T>::copy(Constants::NUM_LEGS, kVec + 1, 2, k, 1); }
    void get_k_vec_wheel(T* k) { mkl<T>::copy(Constants::NUM_LEGS, kVec, 2, k, 1); }

    void get_Mass_vec(T* M) { mkl<T>::copy(Constants::VEC_DIM, massComponents, 1, M, 1); }

    /**
     * get distance vector from each important Point of the car (9: CG, 4*W_i, 4*T_i)
     * \param Point_P,
     * \return each entry from Position_vec
     */
    void get_dist_vector_xy(T* Point_P, T* dist_vector)
    {
        for (auto i = 0; i < Constants::VEC_DIM; ++i) {
            mkl<T>::copy(Constants::DIM - 1, Point_P, Constants::INCX,
                         &dist_vector[(Constants::DIM - 1) * i], Constants::INCX);
        }
        mkl<T>::vSub((Constants::DIM - 1) * Constants::VEC_DIM, currentPositionLagrangian,
                     dist_vector, dist_vector);
    }

    void get_ALE_change(T* current_ALE_vect, T* global_vect, T* change_vect)
    {
        change_vect[0] = current_ALE_vect[0] - global_vect[0];
        change_vect[1] = current_ALE_vect[1] - global_vect[1];
    }

    void apply_ALE_change()
    {
        /*
         * Now both vector are at current state. swap pointer and CG location in new previous will
         * be updated and following will be obselete which
         */
        // Suggestion: reuse the values from the corners (cheaper)

        T c = std::cos(*currentAngleLagrangian);
        T s = std::sin(*currentAngleLagrangian);

        // get the XY positions of all legs taking leg  = CG + R * r
        currentPositionLagrangian[2] =
            currentPositionLagrangian[0] + l_long[0] * c - l_lat[0] * s;  // fl
        currentPositionLagrangian[3] =
            currentPositionLagrangian[1] + l_lat[0] * c + l_long[0] * s;  // fl

        currentPositionLagrangian[6] =
            currentPositionLagrangian[0] + l_long[1] * c + l_lat[1] * s;  // fr
        currentPositionLagrangian[7] =
            currentPositionLagrangian[1] - l_lat[1] * c + l_long[1] * s;  // fr

        currentPositionLagrangian[10] =
            currentPositionLagrangian[0] - l_long[2] * c - l_lat[2] * s;  // rl
        currentPositionLagrangian[11] =
            currentPositionLagrangian[1] + l_lat[2] * c - l_long[2] * s;  // rl

        currentPositionLagrangian[14] =
            currentPositionLagrangian[0] - l_long[3] * c + l_lat[3] * s;  // rr
        currentPositionLagrangian[15] =
            currentPositionLagrangian[1] - l_lat[3] * c - l_long[3] * s;  // rr

        // copy all position value from the wheels to the tyres
        mkl<T>::copy(Constants::NUM_LEGS, currentPositionLagrangian + 2, 4,
                     currentPositionLagrangian + 4, 4);
        mkl<T>::copy(Constants::NUM_LEGS, currentPositionLagrangian + 3, 4,
                     currentPositionLagrangian + 5, 4);
    }

    ~Car()
    {
        mkl_free_buffers();

        /*
        mkl_free(Position_vec);
        mkl_free(Velocity_vec);
        Velocity_vec = nullptr;
        mkl_free(massComponents);
        massComponents = nullptr;
        mkl_free(angle_CG);
        mkl_free(w_CG);
        w_CG = nullptr;
        mkl_free(momentOfInertia);
        momentOfInertia = nullptr;


        // Initial Conditions of the car
        mkl_free(initialPositionGlobal);
        initialPositionGlobal = nullptr;
        mkl_free(initialVelocityGlobal);
        initialVelocityGlobal = nullptr;
        mkl_free(initialAngleGlobal);
        initialAngleGlobal = nullptr;
        mkl_free(initialAngularVelocityGlobal);
        initialAngularVelocityGlobal = nullptr;


        mkl_free(unexcitedSpringsLength);
        unexcitedSpringsLength = nullptr;
        mkl_free(currentCIRTwoTrackModel);
        currentCIRTwoTrackModel = nullptr;
        mkl_free(currentSpringsLength);
        currentSpringsLength = nullptr;
        mkl_free(currentDisplacementTwoTrackModel);
        currentDisplacementTwoTrackModel = nullptr;
        mkl_free(currentPositionTwoTrackModel);
        currentPositionTwoTrackModel = nullptr;
        mkl_free(unexcitedPositionTwoTrackModel);
        unexcitedPositionTwoTrackModel = nullptr;
        mkl_free(kVec);
        kVec = nullptr;
        mkl_free(dVec);
        dVec = nullptr;
        mkl_free(l_lat);
        mkl_free(l_long);
        l_long = nullptr;
        mkl_free(currentVelocityTwoTrackModel);
        currentVelocityTwoTrackModel = nullptr;
        mkl_free(tyre_index_set);

        // ALE Vectors

        mkl_free(currentPositionLagrangian);
        currentPositionLagrangian = nullptr;

        mkl_free(currentVelocityLagrangian);
        currentVelocityLagrangian = nullptr;

        delete currentAngleLagrangian;
        currentAngleLagrangian = nullptr;
        delete currentAngularVelocityLagrangian;
        currentAngularVelocityLagrangian = nullptr;
        delete massFullCar;
        massFullCar = nullptr;

        // Interpolator objects

        mkl_free(currentCornerPositions);
        currentCornerPositions = nullptr;
        mkl_free(currentRotationMatrix);
        currentRotationMatrix = nullptr;
        mkl_free(relativeCornerPositions);
        relativeCornerPositions = nullptr;
        mkl_free(angle_buffer);
        mkl_free(pos_buffer);
        */
    }
};

}  // namespace EVAA
