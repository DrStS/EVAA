// TODO: Copyright header

#pragma once

#include <cmath>

#include "Constants.h"
#include "MathLibrary.h"
#include "MetaDatabase.h"

//#define DAMPING 1  // TODO remove

namespace EVAA {

template <typename T>
class Car {
private:
    /**
     * Copy all X and Y coordinates of the global vector to the local vector
     * \param Global_vector vector with coordinates [X,Y,Z,X,Y,Z,...]
     * \param local_vector vector with coordinates [X,Y,X,Y,...]
     */
    void ConstructALEVectors(T* globalVector, T* localVector) {
        for (size_t i = 0; i < Constants::VEC_DIM; ++i) {
            Math::copy<T>(Constants::DIM - 1, globalVector + i * Constants::DIM, 1, localVector + i * (Constants::DIM - 1), 1);
        }
    }

	
	inline void RotationMatrixSmallAngle(const T arg1, const T arg2, const T arg3, T* arg4) { Math::GetRotationMatrixSmallAnglesXY<T>(arg1, arg2, arg3, arg4); }
	
	inline void RotationMatrixGeneral (const T arg1, const T arg2, const T arg3, T* arg4) { Math::GetRotationMatrix<T>(arg1, arg2, arg3, arg4); }

    /**
     * \brief compute the positions of the unexciting Euler system in the Euler
     * frame
     */
    void ComputeUnexcitedTwoTrackPosition() {
        _unexcitedPositionTwoTrackModel[0] = _initialPositionGlobal[2];  // CG
        _unexcitedPositionTwoTrackModel[1] = 0;                         // angles
        _unexcitedPositionTwoTrackModel[2] = 0;
        _unexcitedPositionTwoTrackModel[3] = _initialPositionGlobal[2] - _unexcitedSpringsLength[0];  // fl
        _unexcitedPositionTwoTrackModel[4] = _initialPositionGlobal[2] - _unexcitedSpringsLength[0] - _unexcitedSpringsLength[1];
        _unexcitedPositionTwoTrackModel[5] = _initialPositionGlobal[2] - _unexcitedSpringsLength[2];  // fr
        _unexcitedPositionTwoTrackModel[6] = _initialPositionGlobal[2] - _unexcitedSpringsLength[2] - _unexcitedSpringsLength[3];
        _unexcitedPositionTwoTrackModel[7] = _initialPositionGlobal[2] - _unexcitedSpringsLength[4];  // rl
        _unexcitedPositionTwoTrackModel[8] = _initialPositionGlobal[2] - _unexcitedSpringsLength[4] - _unexcitedSpringsLength[5];
        _unexcitedPositionTwoTrackModel[9] = _initialPositionGlobal[2] - _unexcitedSpringsLength[6];  // rr
        _unexcitedPositionTwoTrackModel[10] = _initialPositionGlobal[2] - _unexcitedSpringsLength[6] - _unexcitedSpringsLength[7];
    }

    /**
     * \brief Initial displacements as required for the 11DOF
     */
    void ComputeInitialDisplacement() {
        _currentDisplacementTwoTrackModel[0] = 0;
        _currentDisplacementTwoTrackModel[1] = _initialAngleGlobal[0];
        _currentDisplacementTwoTrackModel[2] = _initialAngleGlobal[1];
        _currentDisplacementTwoTrackModel[3] = _unexcitedSpringsLength[0] - _currentSpringsLengths[0] + _currentCornerPositions[8] - _initialPositionGlobal[2];
        _currentDisplacementTwoTrackModel[4] = _currentDisplacementTwoTrackModel[3] + _unexcitedSpringsLength[1] - _currentSpringsLengths[1];
        _currentDisplacementTwoTrackModel[5] = _unexcitedSpringsLength[2] - _currentSpringsLengths[2] + _currentCornerPositions[9] - _initialPositionGlobal[2];
        _currentDisplacementTwoTrackModel[6] = _currentDisplacementTwoTrackModel[5] + _unexcitedSpringsLength[3] - _currentSpringsLengths[3];
        _currentDisplacementTwoTrackModel[7] = _unexcitedSpringsLength[4] - _currentSpringsLengths[4] + _currentCornerPositions[10] - _initialPositionGlobal[2];
        _currentDisplacementTwoTrackModel[8] = _currentDisplacementTwoTrackModel[7] + _unexcitedSpringsLength[5] - _currentSpringsLengths[5];
        _currentDisplacementTwoTrackModel[9] = _unexcitedSpringsLength[6] - _currentSpringsLengths[6] + _currentCornerPositions[11] - _initialPositionGlobal[2];
        _currentDisplacementTwoTrackModel[10] = _currentDisplacementTwoTrackModel[9] + _unexcitedSpringsLength[7] - _currentSpringsLengths[7];
    }

    

     /**
     * Calculates the values of Corners for general angles.
     */
    void UpdateCorners11DOF(T* angles, T* rotation_mat_buffer, T* initial_corners, T* updated_corners) {
        // zz, yy, xx
		(this->*_RotationMatrix)(angles[2], angles[1], angles[0], rotation_mat_buffer);

        // do rotation: rotationMat * r
	    Math::gemm<T>(CblasRowMajor, CblasNoTrans, CblasNoTrans, Constants::DIM, Constants::NUM_LEGS, Constants::DIM, 1, rotation_mat_buffer, Constants::DIM, initial_corners, Constants::NUM_LEGS, 0, updated_corners, Constants::NUM_LEGS);
    }

    /**
     * Calculates the values of _currentCornerPositions according to the current
     * orientation
     */
    void UpdateCorners11DOF() {
        _positionBuffer[0] = _currentPositionLagrangian[0];
        _positionBuffer[1] = _currentPositionLagrangian[1];
        _positionBuffer[2] = _unexcitedPositionTwoTrackModel[0] + _currentDisplacementTwoTrackModel[0];
        _angleBuffer[0] = _unexcitedPositionTwoTrackModel[1] + _currentDisplacementTwoTrackModel[1];
        _angleBuffer[1] = _unexcitedPositionTwoTrackModel[2] + _currentDisplacementTwoTrackModel[2];
        _angleBuffer[2] = _currentAngleLagrangian;

        ConstructCornerRelativeToCG(_relativeCornerPositions);
        UpdateCorners11DOF(_angleBuffer, _currentRotationMatrix, _relativeCornerPositions, _currentCornerPositions);
        CornerAboutCenter(_currentCornerPositions, _positionBuffer);
    }

        /**
     * \brief This function computes conter position about CG
     * \param corner: Array of corners about origin in column major format
     * \param center: Array of coordinate of center of mass
     */
    void CornerAboutCenter(T* corner, T* center) {
        for (size_t i = 0; i < Constants::NUM_LEGS; ++i) {
            Math::axpy<T>(Constants::DIM, 1, center, 1, corner + i, Constants::NUM_LEGS);
        }
    }

    // TODO: Was ist das? Maybe to Math with name 2D to 3D
    void ConvertALEToGlobal(T* vector, T* globalVector) {
#pragma loop(ivdep)
        for (size_t i = 0; i < Constants::VEC_DIM; ++i) {
            globalVector[Constants::DIM * i] = vector[(Constants::DIM - 1) * i];
            globalVector[Constants::DIM * i + 1] = vector[(Constants::DIM - 1) * i + 1];
        }
    }

    // TODO: Was ist das? Maybe to Math Sure!!!!!!!
    void Convert11DOFToGlobal(T* vector, T* globalVector) {
        globalVector[2] = vector[0];
#pragma loop(ivdep)
        for (size_t i = 1; i < Constants::VEC_DIM; ++i) {
            globalVector[Constants::DIM * i + 2] = vector[(Constants::DIM - 1) + i];
        }
    }

    /**
     * Get the solution vector as required for the 11DOF system
     * \param Global_position in the format [GC:XYZ,W_fl:XYZ,T_fl:XYZ,W_fr:XYZ,T_fr:XYZ,...]
     * \param Global_angle with three angles [X,Y,Z]
     * \return Position_11dof in the format [GC:Y,angle:XY,W_fl:Y,T_fl:Y,W_fr:Y,T_fr:Y,...]
     */
    // TODO move to 11DOF ?
    void Construct11DOFVector(T* Global_position, T* Global_angle, T* Position_11dof) {
        Position_11dof[0] = Global_position[2];  // z coordinate of CG
        Position_11dof[1] = Global_angle[0];     // x angle of the CG
        Position_11dof[2] = Global_angle[1];     // y angle of the CG

        // copy y coordinate in order wheel, tyre, wheel, tyre, wheel, tyre, ...
        Math::copy<T>(Constants::VEC_DIM - 1, Global_position + 5, 3, Position_11dof + 3, 1);
    }

    void ComputeGlobalPositionWheelTyre(T* Corners, T* curr_spring_len, T* W_fl, T* T_fl, T* W_fr, T* T_fr, T* W_rl, T* T_rl, T* W_rr, T* T_rr) {
        // W_fl & T_fl
        Math::copy<T>(Constants::DIM, Corners, 4, W_fl, 1);
        W_fl[2] -= curr_spring_len[0];
        Math::copy<T>(Constants::DIM, W_fl, 1, T_fl, 1);
        T_fl[2] -= curr_spring_len[1];

        // W_fr & T_fr
        Math::copy<T>(Constants::DIM, Corners + 1, 4, W_fr, 1);
        W_fr[2] -= curr_spring_len[2];
        Math::copy<T>(Constants::DIM, W_fr, 1, T_fr, 1);
        T_fr[2] -= curr_spring_len[3];

        // W_rl & T_rl
        Math::copy<T>(Constants::DIM, Corners + 2, 4, W_rl, 1);
        W_rl[2] -= curr_spring_len[4];
        Math::copy<T>(Constants::DIM, W_rl, 1, T_rl, 1);
        T_rl[2] -= curr_spring_len[5];

        // W_rr & T_rr
        Math::copy<T>(Constants::DIM, Corners + 3, 4, W_rr, 1);
        W_rr[2] -= curr_spring_len[6];
        Math::copy<T>(Constants::DIM, W_rr, 1, T_rr, 1);
        T_rr[2] -= curr_spring_len[7];
    }
        
    /**
     * \*brief compute current spring lengths
     */
    void ComputeCurrentSpringLength() {
        // upper spring_length = corner - wheel
        // lower spring length = wheel - tyre
        Math::copy<T>(Constants::NUM_LEGS, _currentCornerPositions + 8, 1, _currentSpringsLengths, 2);
        Math::axpy<T>(Constants::NUM_LEGS, -1, _currentPositionTwoTrackModel + 3, 2, _currentSpringsLengths, 2);
        Math::copy<T>(Constants::NUM_LEGS, _currentPositionTwoTrackModel + 3, 2, _currentSpringsLengths + 1, 2);
        Math::axpy<T>(Constants::NUM_LEGS, -1, _currentPositionTwoTrackModel + 4, 2, _currentSpringsLengths + 1, 2);
    }

    /**
     * \brief updates the global position of the components based on euler frame
     * displacement
     */
    void UpdateGlobalTwoTrackVectors() {
        // _currentPositionTwoTrackModel = _unexcitedPositionTwoTrackModel +
        // _currentDisplacementTwoTrackModel;
        Math::copy<T>(Constants::DOF, _unexcitedPositionTwoTrackModel, 1, _currentPositionTwoTrackModel, 1);
        Math::axpy<T>(Constants::DOF, 1, _currentDisplacementTwoTrackModel, 1, _currentPositionTwoTrackModel, 1);
    }

    /**
     * \brief calculate the Z-distance from each car component to the CIR
     */
    void UpdateRadiusToCIR() {
        T globalNickpolPosition = _currentPositionTwoTrackModel[0] + _vehicleCIR[2];  // get global Z position of the Nickpol
        _currentCIRTwoTrackModel[0] = -_vehicleCIR[2];                                // initialize in the constructor
        for (int i = 1; i < 2 * Constants::NUM_LEGS + 1; ++i) {
            _currentCIRTwoTrackModel[i] = _currentPositionTwoTrackModel[2 + i] - globalNickpolPosition;
        }
    }

    // Sums up all the 9 masses (make it private)
    inline T getMassFullCar() {
        return (_massComponents[0] +  // CG
                _massComponents[1] + _massComponents[2] + _massComponents[3] + _massComponents[4] + _massComponents[5] + _massComponents[6] + _massComponents[7] + _massComponents[8]);
    }
    
    // physical parameters
    T* _massComponents;  // [CG,  W_fl, T_fl, W_fr, T_fr, W_rl, T_rl, W_rr, T_rr]
    T _massFullCar;      // 1 scalar
    T* _momentOfInertia;  // [Ixx, Ixy, Ixz, Iyx, Iyy, Iyz, Izx, Izy, Izz]
    T* _vehicleCIR;       // [XYZ]


    // TODO only use to write to HDF5 DISCUSS
    T* _PositionVector;  // [CG:XYZ, W_fl:XYZ, T_fl:XYZ, W_fr:XYZ, T_fr:XYZ,
                         // W_rl:XYZ, T_rl:XYZ, W_rr:XYZ, T_rr:XYZ] n x 9 x 3 !!!
                         // global Consider Constants::ALIGNMENT (3+1),(3+1),...
    T* _VelocityVector;  // [CG:XYZ, W_fl:XYZ, T_fl:XYZ, W_fr:XYZ, T_fr:XYZ,
                         // W_rl:XYZ, T_rl:XYZ, W_rr:XYZ, T_rr:XYZ] n x 9 x 3 !!!
                         // global
    T* _angleCG;         // [XYZ]
    T* _wCG;             // [XYZ]

    // Initial Conditions of the car
    T* _initialPositionGlobal;  // [CG:XYZ, W_fl:XYZ, T_fl:XYZ, W_fr:XYZ,
                                // T_fr:XYZ, W_rl:XYZ, T_rl:XYZ, W_rr:XYZ,
                                // T_rr:XYZ] !!! global
    T* _initialVelocityGlobal;  // [CG:XYZ, W_fl:XYZ, T_fl:XYZ, W_fr:XYZ,
                                // T_fr:XYZ, W_rl:XYZ, T_rl:XYZ, W_rr:XYZ,
                                // T_rr:XYZ] !!! global
    T* _initialAngleGlobal;     // [XYZ]
    T* _initialAngularVelocityGlobal;  // [XYZ]
    T* _relativeCornerPositions;       // [X:fl,fr,rl,rr Y:fl,fr,rl,rr Z:fl,fr,rl,rr]
    bool _smallAngleApprox = true;     // TODO Read this from MetaDatabase
    
    T* _currentPositionTwoTrackModel;  // [CG:Z, theta:XY, W_fl:Z, T_fl:Z,
                                       // W_fr:Z, T_fr:Z, W_rl:Z, T_rl:Z, W_rr:Z,
                                       // T_rr:Z]
    T* _currentCIRTwoTrackModel;       // [CG:Z, W_fl:Z, T_fl:Z, W_fr:Z, T_fr:Z,
                                       // W_rl:Z, T_rl:Z, W_rr:Z, T_rr:Z]
    T* _currentSpringsLengths;         // [W_fl:Z, T_fl:Z, W_fr:Z, T_fr:Z, W_rl:Z,
                                       // T_rl:Z, W_rr:Z, T_rr:Z]
    // For global necessary updates
    T* _currentRotationMatrix;   // [xx, xy, xz, yx, yy, yz, zx, zy, zz]
    

    void (Car<T>::*_RotationMatrix)(const T, const T, const T, T*);

    // Members from 11 DOF system
    T* _kVec;
    T* _dVec;

    // Interpolator Members
    T *_angleBuffer, *_positionBuffer;  // TODO to be removed

public:
    inline const T* getMassComponents() const { return _massComponents; }
    inline const T getMassCarFull() const { return _massFullCar; }
    inline const T* getMomentOfInertia() const { return _momentOfInertia; }
    inline const T* getPositionVector() const { return _PositionVector; }
    inline const T* getAngleCG() const { return _angleCG; }
    void setInitialAngularVelocityGlobalZ(T AngularVelocityLagrangian) { _initialAngularVelocityGlobal[2] = AngularVelocityLagrangian; }
	void setInitialAngleGlobalZ(T AngleLagrangian) { _initialAngleGlobal[2] = AngleLagrangian; }
    inline const T* getCurrentVelocityTwoTrackModel() const { return _currentVelocityTwoTrackModel; }
    inline T* getCurrentDisplacementTwoTrackModel() const { return _currentDisplacementTwoTrackModel; } // TODO consider friend
    inline const T* getCurrentCIRTwoTrackModel() const { return _currentCIRTwoTrackModel; }
    inline const T* getCurrentSpringsLengths() const { return _currentSpringsLengths; }
    inline const T* getCurrentPositionLagrangian() const { return _currentPositionLagrangian; }
	inline const T* getCurrentPositionTwoTrackModel() const { return _currentPositionTwoTrackModel; }
    inline const T getCurrentAngularVelocityLagrangian() const { return _currentAngularVelocityLagrangian; }
    inline const T* getCurrentVelocityLagrangian() const { return _currentVelocityLagrangian; }
    T getMomentOfInertiaLagrangian() const { return _lagrangianMomentOfInertia; }
    void setMomentOfInertiaLagrangian(T moment) { _lagrangianMomentOfInertia = moment; }
    inline const T* getCurrentCornerPositions() const { return _currentCornerPositions; }
	inline const T* getUnexcitedPositionTwoTrackModel() const { return _unexcitedPositionTwoTrackModel; }
    inline T* getkVec() const { return _kVec; } // TODO consider friend
    inline T* getdVec() const { return _dVec; } // TODO consider friend


public:
    /*
    * Convention:
    *   W1 = front left wheel
    *   W2 = front right wheel
    *   W3 = rear left wheel
    *   W4 = rear right wheel
    *   T1 = front left tyre
    *   T2 = front right tyre
    *   T3 = rear left tyre
    *   T4 = rear right tyre
    */

    // vectors to use in regular iteration
    // For 11DOF
    T* _currentVelocityTwoTrackModel;      // [CG:Z, theta:XY, W_fl:Z, T_fl:Z,
                                          // W_fr:Z, T_fr:Z, W_rl:Z, T_rl:Z, W_rr:Z,
                                          // T_rr:Z] (TODO: still public)
    
    T* _currentDisplacementTwoTrackModel;  // [CG:Z, theta:XY, W_fl:Z, T_fl:Z,
                                          // W_fr:Z, T_fr:Z, W_rl:Z, T_rl:Z,
                                          // W_rr:Z, T_rr:Z] (TODO: still public)

    T* _unexcitedPositionTwoTrackModel;    // [CG:Z, theta:XY, W_fl:Z, T_fl:Z,
                                           // W_fr:Z, T_fr:Z, W_rl:Z, T_rl:Z,
                                           // W_rr:Z, T_rr:Z] (TODO: still public, no getter)

    T* _unexcitedSpringsLength;				// [W_fl:Z, T_fl:Z, W_fr:Z, T_fr:Z, W_rl:Z,
											// T_rl:Z, W_rr:Z, T_rr:Z]

    
    // For ALE || suggestion: reduce it to only leg position (wheel == tyre)
    T* _currentPositionLagrangian;        // [CG:XY, W_fl:XY, T_fl:XY, W_fr:XY, T_fr:XY, W_rl:XY, T_rl:XY, W_rr:XY, T_rr:XY]
                                          // (TODO: still public) 

    T _currentAngularVelocityLagrangian;  // [CG:Z] (TODO: still public)  consider friend

    T _currentAngleLagrangian;            // [CG:Z] (TODO: still public)  consider friend

    T* _currentVelocityLagrangian;        // [CG:XY, W_fl:XY, T_fl:XY, W_fr:XY,
                                         // T_fr:XY, W_rl:XY, T_rl:XY, W_rr:XY,
                                         // T_rr:XY] (TODO: still public) consider friend

    T _lagrangianMomentOfInertia;         // corresponds to the whole car object (TODO:as before)
    
    
	T* _currentCornerPositions;  // [X:fl,fr,rl,rr Y:fl,fr,rl,rr Z:fl,fr,rl,rr]
								 // sorry for everyone

    // Members from 11 DOF system
    T* _lenLat;
    T *_lenLong;

    /* Constructor */
    Car() {
        // System Mono

        // Memory Allocation and matrix formulation

        const size_t positionAllocSize = Constants::DIM * Constants::VEC_DIM;  // 27 dimensions

		if (_smallAngleApprox) {
			_RotationMatrix = &Car<T>::RotationMatrixSmallAngle;
		}
		else {
			_RotationMatrix = &Car<T>::RotationMatrixGeneral;
		}

        // Params for Global coordinate

        _PositionVector = Math::malloc<T>(positionAllocSize);
        _VelocityVector = Math::malloc<T>(positionAllocSize);
        _massComponents = Math::malloc<T>(Constants::VEC_DIM);
        _angleCG = Math::malloc<T>(Constants::DIM);
        _wCG = Math::malloc<T>(Constants::DIM);
        _momentOfInertia = Math::malloc<T>(Constants::DIMDIM);

        // Initial Params for Global coordinate

        _initialPositionGlobal = Math::malloc<T>(positionAllocSize);
        _initialVelocityGlobal = Math::malloc<T>(positionAllocSize);
        _initialAngleGlobal = Math::malloc<T>(Constants::DIM);
        _initialAngularVelocityGlobal = Math::malloc<T>(Constants::DIM);

        // Params for 11 DOF system

        _currentDisplacementTwoTrackModel = Math::malloc<T>(Constants::DOF);
        _currentVelocityTwoTrackModel = Math::malloc<T>(Constants::DOF);
        _kVec = Math::calloc<T>(2 * Constants::NUM_LEGS); // TODO malloc
        _dVec = Math::calloc<T>(2 * Constants::NUM_LEGS);  // TODO malloc
        _lenLat = Math::malloc<T>(Constants::NUM_LEGS);
        _lenLong = Math::malloc<T>(Constants::NUM_LEGS);
        _unexcitedSpringsLength = Math::malloc<T>(2 * Constants::NUM_LEGS);
        _currentCIRTwoTrackModel = Math::malloc<T>(Constants::DIMDIM);
        _currentSpringsLengths = Math::malloc<T>(2 * Constants::NUM_LEGS);
        _currentPositionTwoTrackModel = Math::malloc<T>(Constants::DOF);
        _unexcitedPositionTwoTrackModel = Math::malloc<T>(Constants::DOF);

        // Memory allocation for interpolator

        _currentCornerPositions = Math::malloc<T>(Constants::NUM_LEGS * Constants::DIM);
        _currentRotationMatrix = Math::malloc<T>(Constants::DIMDIM);
        _relativeCornerPositions = Math::malloc<T>(Constants::NUM_LEGS * Constants::DIM);
        _angleBuffer = Math::malloc<T>(Constants::DIM);
        _positionBuffer = Math::malloc<T>(Constants::DIM);

        // ALE Buffer Allocation

        _currentPositionLagrangian = Math::malloc<T>((Constants::DIM - 1) * Constants::VEC_DIM);
        _currentVelocityLagrangian = Math::malloc<T>((Constants::DIM - 1) * Constants::VEC_DIM);

        // Extract Data from parser
        auto& db = MetaDatabase<T>::getDatabase();

        Math::copy<T>(Constants::NUM_LEGS, db.getLongitudalLegPositionVector(), 1, _lenLong, 1);

        Math::copy<T>(Constants::NUM_LEGS, db.getLatidudalLegPositionVector(), 1, _lenLat, 1);

        Math::copy<T>(Constants::DIMDIM, db.getMomentOfInertiaVector(), 1, _momentOfInertia, 1);
		
        // TODO: Current order [CG, W W W W T T T T] does not perform at its
        // best for the [W+T] part (alignment is not to 64). Rework the
        // formulation to put CG at the end.
        _massComponents[0] = db.getBodyMass();
        Math::copy<T>(2 * Constants::NUM_LEGS, db.getWheelTyreMassVector(), 1, _massComponents + 1, 1);

        _massFullCar = getMassFullCar();

        Math::copy<T>(Constants::NUM_LEGS, db.getBodySpringLengthVector(), 1, _unexcitedSpringsLength, 2);
        Math::copy<T>(Constants::NUM_LEGS, db.getTyreSpringLengthVector(), 1, _unexcitedSpringsLength + 1, 2);

        _vehicleCIR = db.getPositionCenterOfInstantaneousRotation();
		
        // Initial Iteration vector

        // Initial Angles
        Math::ToEulerAngles<T>(db.getBodyInitialOrientation(), _initialAngleGlobal);
        Math::copy<T>(Constants::DIM, _initialAngleGlobal, 1, _angleCG, 1);

		

        // Spring lengths
        Math::copy<T>(Constants::NUM_LEGS, db.getBodySpringInitialLengthVector(), 1, _currentSpringsLengths, 2);
        Math::copy<T>(Constants::NUM_LEGS, db.getTyreSpringInitialLengthVector(), 1, _currentSpringsLengths + 1, 2);

		
		

        // Filling the position vector with initial condition
        // CG
        Math::copy<T>(Constants::DIM, db.getBodyInitialPosition(), 1, _initialPositionGlobal, 1);  // copy the center of mass position
		


        // Interpolator
        // Initialization: read init corners vectors into matrix

        ConstructCornerRelativeToCG(_relativeCornerPositions);  // only CG position is used to construct corners
        Math::copy<T>(Constants::DIM, _initialAngleGlobal, 1, _angleBuffer, 1);
        // _angleBuffer[2] = 0;
        UpdateCorners11DOF(_angleBuffer, _currentRotationMatrix, _relativeCornerPositions, _currentCornerPositions);
        CornerAboutCenter(_currentCornerPositions, _initialPositionGlobal);
		
        // Remaining position initialization

        const T* _xmlStart;
        T* position_start;
        if (db.getFlagInitialLeg()) {  // TODO DEBUG to  be removed
            // if prescribed initial position (add a check for consistency with spring lengths)
            // W1 = W_fl
            _xmlStart = db.getWheelInitialPositionFrontLeft();
            position_start = _initialPositionGlobal + 3;  //(end at 5)
            Math::copy<T>(Constants::DIM, _xmlStart, 1, position_start, 1);
            // W2 = W_fr
            _xmlStart = db.getWheelInitialPositionFrontRight();
            position_start += 2 * Constants::DIM;  // skip 3 for tyre (end at 11)
            Math::copy<T>(Constants::DIM, _xmlStart, 1, position_start, 1);
            // W3 = W_rl
            _xmlStart = db.getWheelInitialPositionRearLeft();
            position_start += 2 * Constants::DIM;  // skip 3 for tyre (end at 17)
            Math::copy<T>(Constants::DIM, _xmlStart, 1, position_start, 1);
            // W2 = W_rr
            _xmlStart = db.getWheelInitialPositionRearRight();
            position_start += 2 * Constants::DIM;  // skip 3 for tyre (end at 23)
            Math::copy<T>(Constants::DIM, _xmlStart, 1, position_start, 1);

            // T1 = T_fl
            _xmlStart = db.getTyreInitialPositionFrontLeft();
            position_start = _initialPositionGlobal + 6;  // skip 3 for center of mass and 3 for the wheel
            Math::copy<T>(Constants::DIM, _xmlStart, 1, position_start, 1);  // (end at 8)
            // T2 = T_fr
            _xmlStart = db.getTyreInitialPositionFrontRight();
            position_start += 2 * Constants::DIM;  // skip 3 for the wheel
            Math::copy<T>(Constants::DIM, _xmlStart, 1, position_start,
                          1);  // (end at 14)
            // T3 = T_rl
            _xmlStart = db.getTyreInitialPositionRearLeft();
            position_start += 2 * Constants::DIM;  // skip 3 for the wheel
            Math::copy<T>(Constants::DIM, _xmlStart, 1, position_start,
                          1);  // (end at 20)
            // T4 = T_rr
            _xmlStart = db.getTyreInitialPositionRearRight();
            position_start += 2 * Constants::DIM;                                             // skip 3 for the wheel
            Math::copy<T>(Constants::DIM, _xmlStart, 1, position_start, 1);  // (end at 26)
            Math::copy<T>(Constants::DIM * Constants::VEC_DIM, _initialPositionGlobal, 1, _PositionVector, 1);
        }
        else {
            T* W_fl = _initialPositionGlobal + 3;
            T* W_fr = _initialPositionGlobal + 9;
            T* W_rl = _initialPositionGlobal + 15;
            T* W_rr = _initialPositionGlobal + 21;
            T* T_fl = _initialPositionGlobal + 6;
            T* T_fr = _initialPositionGlobal + 12;
            T* T_rl = _initialPositionGlobal + 18;
            T* T_rr = _initialPositionGlobal + 24;
            ComputeGlobalPositionWheelTyre(_currentCornerPositions, _currentSpringsLengths, W_fl, T_fl, W_fr, T_fr, W_rl, T_rl, W_rr, T_rr);
            /* Update the mean position where changes are to be added */
            Math::copy<T>(Constants::DIM * Constants::VEC_DIM, _initialPositionGlobal, 1, _PositionVector, 1);
        }
		
        // Initial Velocity (Reuse the pointers)
        Math::copy<T>(Constants::DIM, db.getBodyInitialVelocity(), 1, _initialVelocityGlobal, 1);
        // W1 = W_fl
        _xmlStart = db.getWheelInitialVelocityFrontLeft();
        position_start = _initialVelocityGlobal + 3;
        Math::copy<T>(Constants::DIM, _xmlStart, 1, position_start, 1);  // (end at 5)
        // W2 = W_fr
        _xmlStart = db.getWheelInitialVelocityFrontRight();
        position_start += 2 * Constants::DIM;  // skip 3 for tyre
        Math::copy<T>(Constants::DIM, _xmlStart, 1, position_start, 1);  // (end at 11)
        // W3 = W_rl
        _xmlStart = db.getWheelInitialVelocityRearLeft();
        position_start += 2 * Constants::DIM;  // skip 3 for tyre
        Math::copy<T>(Constants::DIM, _xmlStart, 1, position_start, 1);  // (end at 17)
        // W2 = W_rr
        _xmlStart = db.getWheelInitialVelocityRearRight();
        position_start += 2 * Constants::DIM;  // skip 3 for tyre
        Math::copy<T>(Constants::DIM, _xmlStart, 1, position_start, 1);  // (end at 23)

        // T1 = T_fl
        _xmlStart = db.getTyreInitialVelocityFrontLeft();
        position_start = _initialVelocityGlobal + 6;  // skip 3 for center of mass and 3 for the wheel
        Math::copy<T>(Constants::DIM, _xmlStart, 1, position_start, 1);  // (end at 8)

        // T2 = T_fr
        _xmlStart = db.getTyreInitialVelocityFrontRight();
        position_start += 2 * Constants::DIM;  // skip 3 for the Tyre
        Math::copy<T>(Constants::DIM, _xmlStart, 1, position_start, 1);  // (end at 14)
        // T3 = T_rl
        _xmlStart = db.getTyreInitialVelocityRearLeft();
        position_start += 2 * Constants::DIM;  // skip 3 for the wheel
        Math::copy<T>(Constants::DIM, _xmlStart, 1, position_start, 1);  // (end at 20)

        // T4 = T_rr
        _xmlStart = db.getTyreInitialVelocityRearRight();
        position_start += 2 * Constants::DIM;  // skip 3 for the wheel
        Math::copy<T>(Constants::DIM, _xmlStart, 1, position_start, 1);  // (end at 26)

        // copy the initial velocity to the velocity vector
        Math::copy<T>(Constants::DIM * Constants::VEC_DIM, _initialVelocityGlobal, 1, _VelocityVector, 1);
		Math::write_vector(_initialVelocityGlobal, 27);
		Math::write_vector(_VelocityVector, 27);
        // Initial Angular velocity
        Math::copy<T>(Constants::DIM, db.getBodyInitialAngularVelocity(), 1, _initialAngularVelocityGlobal, 1);
        Math::copy<T>(Constants::DIM, _initialAngularVelocityGlobal, 1, _wCG, 1);

        /* Global assignments Done */

        // 11 DOF Buffer Initialization

        // unexcited state
        ComputeUnexcitedTwoTrackPosition();

        // solution vector of 11DOF
        ComputeInitialDisplacement();

        // initial velocity
        Construct11DOFVector(_initialVelocityGlobal, _initialAngularVelocityGlobal, _currentVelocityTwoTrackModel);

        // initial position
        Construct11DOFVector(_initialPositionGlobal, _initialAngleGlobal, _currentPositionTwoTrackModel);

        // we updated the initial global position, so we TODO ??!
        UpdateRadiusToCIR();

#ifdef INTERPOLATION
        db.getLookupStiffness().getInterpolation(_currentSpringsLengths, _kVec);
#ifdef DAMPING
        db.getLookupDamping().getInterpolation(_currentSpringsLengths, _dVec);
#endif
#else
        // copy constant stiffnesses
        Math::copy<T>(Constants::NUM_LEGS, db.getBodyStiffnessVector(), 1, _kVec, 2);
        Math::copy<T>(Constants::NUM_LEGS, db.getTyreStiffnessVector(), 1, _kVec + 1, 2);

        // copy constant dampings
        Math::copy<T>(Constants::NUM_LEGS, db.getBodyDampingVector(), 1, _dVec, 2);
        Math::copy<T>(Constants::NUM_LEGS, db.getTyreDampingVector(), 1, _dVec + 1, 2);       
#endif

        // ALE Buffer Initialization

        ConstructALEVectors(_initialPositionGlobal, _currentPositionLagrangian);
        ConstructALEVectors(_initialVelocityGlobal, _currentVelocityLagrangian);
        _currentAngleLagrangian = _initialAngleGlobal[2];
        _currentAngularVelocityLagrangian = _initialAngularVelocityGlobal[2];
    }

    /**
     * \brief update global corner positions, compute global Euler vectors,
     * compute current spring lengths, compute CIR functions
     */
    void UpdateLengthsTwoTrackModel() {

        UpdateCorners11DOF();

        // global vector update
        UpdateGlobalTwoTrackVectors();
/*        std::cout << "new" << std::endl;
        IO::writeVector<T>(_unexcitedPositionTwoTrackModel, 11);
        IO::writeVector<T>(_currentDisplacementTwoTrackModel, 11);
        IO::writeVector<T>(_currentPositionTwoTrackModel, 11);
*/        // compute current spring lengths
        ComputeCurrentSpringLength();

        // Compute CIR distances
        UpdateRadiusToCIR();
    }

    /**
     * Fills the global vector with all entries (checkpointing)
     * \param[out] globalVector Contains the Lagrangian position vector on X and Y components and Eulerian position vector (from 11 DOF) for Z component
     */
    void CombineEulerianLagrangianVectors(T* globalVector) {
        ConvertALEToGlobal(_currentPositionLagrangian, globalVector);
        Convert11DOFToGlobal(_currentPositionTwoTrackModel, globalVector);
    }

    /**
     * \brief flush the result to the output (checkpointing) || combine with function above (TODO unclear documentation)
     */
    void CombineResults() {
        // Math::copy<T>(Constants::VEC_DIM * Constants::DIM, _initialPositionGlobal, 1, _PositionVector, 1);
        ConvertALEToGlobal(_currentPositionLagrangian, _PositionVector);
        Convert11DOFToGlobal(_currentPositionTwoTrackModel, _PositionVector);

        // Angles manually
        _angleCG[0] = _initialAngleGlobal[0] + _currentDisplacementTwoTrackModel[1];
        _angleCG[1] = _initialAngleGlobal[1] + _currentDisplacementTwoTrackModel[2];
        _angleCG[2] = _currentAngleLagrangian;
		_wCG[0] = _currentVelocityTwoTrackModel[1];
		_wCG[1] = _currentVelocityTwoTrackModel[2];
        _wCG[2] = _currentAngularVelocityLagrangian;
        ConvertALEToGlobal(_currentVelocityLagrangian, _VelocityVector);
		Convert11DOFToGlobal(_currentVelocityTwoTrackModel, _VelocityVector);
    }
	/**
	 * Construct corner initilizer
	 */
	void ConstructCornerRelativeToCG(T* corners) {
		corners[0] = +_lenLong[0];  // fl
		corners[4] = +_lenLat[0];   // fl

		corners[1] = +_lenLong[1];  // fr
		corners[5] = -_lenLat[1];   // fr

		corners[2] = -_lenLong[2];  // rl
		corners[6] = +_lenLat[2];   // rl

		corners[3] = -_lenLong[3];  // rr
		corners[7] = -_lenLat[3];   // rr

#pragma loop(ivdep)
		for (auto i = 8; i < 12; i++) {
			corners[i] = 0;
		}
	}
    /**
     * get distance vector from each important Point of the car (9: CG, 4*W_i, 4*T_i)
     * \param[in] Point_P, 
	 * \param[out] dist_vector compute distance of  each component from point P.
     */
    void ComputeDisplacementToPointLagrangian(T* Point_P, T* dist_vector) {
        for (auto i = 0; i < Constants::VEC_DIM; ++i) {
            Math::copy<T>(Constants::DIM - 1, Point_P, Constants::INCX, &dist_vector[(Constants::DIM - 1) * i], Constants::INCX);
        }
        Math::vSub<T>((Constants::DIM - 1) * Constants::VEC_DIM, _currentPositionLagrangian, dist_vector, dist_vector);
    }

	/*
	* Now both vector are at current state. swap pointer and CG location in
	* new previous will be updated and following will be obselete which
	*/
    void ApplyLagrangeChange() {
        // get the XY positions of all legs taking leg  = CG + R * r
        _currentPositionLagrangian[2] = _currentCornerPositions[0];  // fl
        _currentPositionLagrangian[3] = _currentCornerPositions[4];  // fl

        _currentPositionLagrangian[6] = _currentCornerPositions[1];  // fr
        _currentPositionLagrangian[7] = _currentCornerPositions[5];  // fr

        _currentPositionLagrangian[10] = _currentCornerPositions[2];  // rl
        _currentPositionLagrangian[11] = _currentCornerPositions[6];  // rl

        _currentPositionLagrangian[14] = _currentCornerPositions[3];  // rr
        _currentPositionLagrangian[15] = _currentCornerPositions[7];  // rr

        // copy all position value from the wheels to the tyres
        Math::copy<T>(Constants::NUM_LEGS, _currentPositionLagrangian + 2, 4, _currentPositionLagrangian + 4, 4);
        Math::copy<T>(Constants::NUM_LEGS, _currentPositionLagrangian + 3, 4, _currentPositionLagrangian + 5, 4);
    }

    ~Car() {
        // TODO: consider removing after checking performance.
        mkl_free_buffers();

        Math::free<T>(_PositionVector);
        Math::free<T>(_VelocityVector);
        Math::free<T>(_massComponents);
        Math::free<T>(_angleCG);
        Math::free<T>(_wCG);
        Math::free<T>(_momentOfInertia);

        // Initial Conditions of the car
        Math::free<T>(_initialPositionGlobal);
        Math::free<T>(_initialVelocityGlobal);
        Math::free<T>(_initialAngleGlobal);
        Math::free<T>(_initialAngularVelocityGlobal);

        Math::free<T>(_unexcitedSpringsLength);
        Math::free<T>(_currentCIRTwoTrackModel);
        Math::free<T>(_currentSpringsLengths);
        Math::free<T>(_currentDisplacementTwoTrackModel);
        Math::free<T>(_currentPositionTwoTrackModel);
        Math::free<T>(_unexcitedPositionTwoTrackModel);
        Math::free<T>(_kVec);
        Math::free<T>(_dVec);
        Math::free<T>(_lenLat);
        Math::free<T>(_lenLong);
        Math::free<T>(_currentVelocityTwoTrackModel);

        // ALE Vectors

        Math::free<T>(_currentPositionLagrangian);

        Math::free<T>(_currentVelocityLagrangian);

        // Interpolator objects

        Math::free<T>(_currentCornerPositions);
        Math::free<T>(_currentRotationMatrix);
        Math::free<T>(_relativeCornerPositions);
        Math::free<T>(_angleBuffer);
        Math::free<T>(_positionBuffer);
    }
};

}  // namespace EVAA
