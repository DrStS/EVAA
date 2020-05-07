// TODO: Copyright template

/**
 * \file 11DOF.h
 * This file holds the function declaration and definitions of the linear11dof
 * and the linear11dof_full class. \date 04/14/2020
 */

#pragma once

#include "Car.h"
#include "Constants.h"
#include "MathLibrary.h"
#include "LoadModule.h"
#include "MetaDatabase.h"

//#define DAMPING 1

namespace EVAA {

/**
 * \brief class to compute one timestep of the linear 11 dof system in small
 * angle approximation
 */
template <typename T>
class TwoTrackModelParent {
protected:
    /** pointer to car instance with all important car parameter */
    Car<T>* _car;
	LoadModule<T>* _loadModuleObj;

	/** Force Vector for the load module*/
	T* _twoTrackModelForce = nullptr;

    // time step related
    /** solution in next timestep */
    T _factor_h;
    /** solution in next timestep */
    T _h;
	size_t _currentIter;

    T _tolerance;
    int _maxNewtonIteration;
    int _newtonIteration;

    // solution in next timestep
    T* _u_n_p_1 = nullptr;

    // solution in previous timestep
    T* _u_n_m_1 = nullptr;

    // solution in current timestep
    /* pointers to matrices used for the backward euler */
    T* _u_n = nullptr;
    T* _A = nullptr;
    T* _B = nullptr;
    T* _C = nullptr; 
    
    T* _M_h2 = nullptr;
    T* _K = nullptr;
    T* _D = nullptr;
    
    /** used for the constructStiffnesMatrix and ConstructDampingMatrix. */
    T* _matrixTmp = nullptr;
    T* _kVec = nullptr;
    T* _dVec = nullptr;

    /** used for the constructStiffnesMatrix and ConstructDampingMatrix as well as for the BE. */
    T* _temp = nullptr;

#ifdef INTERPOLATION
    T* _J = nullptr;
    T* _dKdxx = nullptr;
#ifdef DAMPING
    T* _dDdxx = nullptr;
    T* _dddl = nullptr;
#endif // DAMPING
    T* _residual = nullptr;
    T _residualNorm;
    T* _dkdl = nullptr;
#endif // INTERPOLATION

public:
    /**
     * Constructor
     */
    TwoTrackModelParent(Car<T>* inputCar, LoadModule<T>* loadModel) : _car(inputCar), _loadModuleObj(loadModel) {}

    /**
     * Performs one timestep of the 11DOF solver
     * \param[in] load vector [angle:Z,GC:Y,W1:Y,T1:Y,W2:Y,T2:Y,...]
     * \oaram[out] solution of the following timestep
     * [angle:Z,GC:Y,W1:Y,T1:Y,W2:Y,T2:Y,...]
     */
    virtual void UpdateStep(size_t currentIter, T* solution) = 0;

    /**
     * Destructor
     */
    virtual ~TwoTrackModelParent() {}
};

/**
 * \brief class to compute one timestep of the linear 11 dof system in small
 * angle approximation
 */
template <typename T>
class TwoTrackModelBE : public TwoTrackModelParent<T> {
private:
	/** For pperformance benefit we keep it here */
	void (TwoTrackModelBE<T>::*_JacobianAdjustment)();

    /**
     * \brief construct _A Matrix
     *
     * this has to be called every time step for interpolation
     * A = 1/h^2 * M + 1/h * D + K
     */
    void ConstructAMatrix() {
        // A = M/h^2 (also acts as initialization of A)
        Math::copy<T>(Constants::DOFDOF, _M_h2, 1, _A, 1);
        // A += 1/h * _D => A = 1/h^2 * M + 1/h * D
        Math::axpy<T>(Constants::DOFDOF, _factor_h, _D, 1, _A, 1);
        // A += K => A = 1/h^2 * M + 1/h * D + K
        Math::axpy<T>(Constants::DOFDOF, 1, _K, 1, _A, 1);
    }
    /**
     * \brief construct _B
     *
     * this has to be called every time step for interpolation
     * _B = 2/(h*h) * M + 1/h * _D
     */
    void ConstructBMatrix() {
        Math::scal<T>(Constants::DOFDOF, 0, _B, 1);
        Math::axpy<T>(Constants::DOFDOF, 2, _M_h2, 1, _B, 1);
        Math::axpy<T>(Constants::DOFDOF, _factor_h, _D, 1, _B, 1);
    }
    /**
     * \brief construct _C
     *
     * this has to be called only once
     * _C = -1/h^2 M
     */
    void ConstructCMatrix() {
        Math::copy<T>(Constants::DOFDOF, _M_h2, 1, _C, 1);
        Math::scal<T>(Constants::DOFDOF, -1, _C, 1);
    }

protected:
    /**
     * \brief construct Mass matrix
     *
     * this is only done once. Therefore, there is no need to store Mass_vec or
     * I_CG
     */
    void ConstructMassMatrix() {
        // get values from MetaDatabase
        const T* Mass_vec = _car->getMassComponents();
        // construct M
        _M_h2[0] = Mass_vec[0];
        _M_h2[Constants::DOF + 1] = _car->getMomentOfInertia()[0];
        _M_h2[2 * Constants::DOF + 2] = _car->getMomentOfInertia()[4];
        // _M_h2 = M
        Math::copy<T>(Constants::DOF - 3, Mass_vec + 1, 1, _M_h2 + 3 * (Constants::DOF + 1),
                      Constants::DOF + 1);  // _M_h2 = diagonal matrix
        // _M_h2 = M / (h * h)
        Math::scal<T>(Constants::DOF, 1. / (_h * _h), _M_h2, Constants::DOF + 1);
    }

	/**
	* \breif Compute velocity of the car based on the displacement vectors
	* 
	* called inside UpdateStep function
	*/
	void UpdateVelocity() {
		Math::copy<T>(Constants::DOF, _u_n_p_1, Constants::INCX, _car->_currentVelocityTwoTrackModel, Constants::INCX);
		Math::axpy<T>(Constants::DOF, -1, _u_n, Constants::INCX, _car->_currentVelocityTwoTrackModel, Constants::INCX);
		Math::scal<T>(Constants::DOF, _factor_h, _car->_currentVelocityTwoTrackModel, Constants::INCX);
	}


    /**
     * \brief initialise the vectors _u_n, _u_n_m_1, _u_n_p_1
     *
     * called in the constructor
     * _u_n = ...
     * _u_n_m_1 = _u_n - _h * velocity
     * _u_n_p_1 = _u_n
     * tempVector is used for velocity
     */
    void InitPosVectors() {
        // _u_n
        Math::copy<T>(Constants::DOF, _u_n_p_1, 1, _u_n, 1);
        // _u_n_m_1 = _u_n - _h * velocity
        Math::copy<T>(Constants::DOF, _car->_currentVelocityTwoTrackModel, 1, _u_n_m_1, 1);
        Math::scal<T>(Constants::DOF, -_h, _u_n_m_1, 1);
        Math::axpy<T>(Constants::DOF, 1, _u_n, 1, _u_n_m_1, 1);
    }

public:
    /**
     * \brief Constructor
     */
    TwoTrackModelBE(Car<T>* inputCar, LoadModule<T>* loadModel) : TwoTrackModelParent(inputCar, loadModel) {
        _u_n_m_1 = Math::malloc<T>(Constants::DOF);					// velocity
        _u_n = Math::malloc<T>(Constants::DOF);						// position
        _u_n_p_1 = _car->getCurrentDisplacementTwoTrackModel();		// position
		_twoTrackModelForce = Math::malloc<T>(Constants::DOF);		// force

        _A = Math::malloc<T>(Constants::DOFDOF);
        _B = Math::calloc<T>(Constants::DOFDOF);
        _C = Math::calloc<T>(Constants::DOFDOF);
        _M_h2 = Math::calloc<T>(Constants::DOFDOF);
        _K = Math::calloc<T>(Constants::DOFDOF);
        _D = Math::calloc<T>(Constants::DOFDOF);

        _temp = Math::malloc<T>(Constants::DOF);
        _matrixTmp = Math::calloc<T>(Constants::DOFDOF);
#ifdef INTERPOLATION
		if (_loadModuleObj->GetEulerProfileName() == "Fixed") {
			_JacobianAdjustment = &TwoTrackModelBE<T>::ConstructFixedJacobian;
		}
		else if (_loadModuleObj->GetEulerProfileName() == "Nonfixed") {
			_JacobianAdjustment = &TwoTrackModelBE<T>::ConstructNonFixedJacobian;
		}
        _J = Math::calloc<T>(Constants::DOFDOF);
        _dKdxx = Math::calloc<T>(Constants::DOFDOF);
		_dkdl = Math::malloc<T>(2 * Constants::NUM_LEGS);
        
#ifdef DAMPING
		_dddl = Math::calloc<T>(2 * Constants::NUM_LEGS);
		_dDdxx = Math::calloc<T>(Constants::DOFDOF);
#endif
		_residual = Math::malloc<T>(Constants::DOF);
#endif
        _tolerance = MetaDatabase<T>::getDatabase().getNewtonTolerance();
        _maxNewtonIteration = MetaDatabase<T>::getDatabase().getMaxNewtonIterations();
        _h = MetaDatabase<T>::getDatabase().getTimeStepSize();
        _factor_h = 1 / _h;

        // if interpolation is used we need all the lengths to interpolate
        // damping and stiffness if not we have to define them populate kVec
        _kVec = _car->getkVec();
        // populate _dVec
        _dVec = _car->getdVec();

        InitPosVectors();
        ConstructMassMatrix();
        ConstructStiffnessMatrix();
#ifdef DAMPING
        ConstructDampingMatrix();
#endif
		ConstructAMatrix();
        ConstructBMatrix();
        ConstructCMatrix();
    }

    /**
     * Destructor
     */
    virtual ~TwoTrackModelBE() {
        Math::free<T>(_A);
        Math::free<T>(_B);
        Math::free<T>(_C);
        Math::free<T>(_u_n);
        Math::free<T>(_u_n_m_1);
		Math::free<T>(_twoTrackModelForce);
        Math::free<T>(_M_h2);
        Math::free<T>(_K);
        Math::free<T>(_D);
        Math::free<T>(_temp);
        Math::free<T>(_matrixTmp);
#ifdef INTERPOLATION
        Math::free<T>(_J);
		Math::free<T>(_dkdl);
		Math::free<T>(_dKdxx);
#ifdef DAMPING
		
        Math::free<T>(_dddl);
		Math::free<T>(_dDdxx);
#endif
        Math::free<T>(_residual);
#endif
    }

    /**
     * Performs one timestep of the 11DOF solver
     * \param load vector [angle:Z,GC:Y,W1:Y,T1:Y,W2:Y,T2:Y,...]
     */
    virtual void UpdateStep(size_t currentIter, T* solution) {
        // Math::scal<T>(Constants::DOF, -1, _u_n_m_1, Constants::INCX);
		// Get the force
		_currentIter = currentIter;
		_loadModuleObj->GetEulerianForce(currentIter, _twoTrackModelForce);
		Math::Solvers<T, TwoTrackModelBE<T> >::LinearBackwardEuler(_A, _B, _C, _u_n, _u_n_m_1, _twoTrackModelForce, _u_n_p_1, Constants::DOF);
        
#ifdef INTERPOLATION
		UpdateSystem();
		Math::Solvers<T, TwoTrackModelBE<T>>::Newton(this, _twoTrackModelForce, _J, _residual, &_residualNorm, _u_n_p_1, _temp, &_tolerance, &_maxNewtonIteration, &_newtonIteration);
#else
		ConstructAMatrix();
#endif
        Math::copy<T>(Constants::DOF, _u_n_p_1, 1, solution, 1);
        
		// update two track velocity in car
		UpdateVelocity();

        // _u_n_m_1 points to _u_n and _u_n points to _u_n_m_1
        Math::SwapAddress<T>(_u_n, _u_n_m_1);
        // _u_n points to _u_n_p_1 and
        // Math::SwapAddress<T>(_u_n_p_1, _u_n);
        // _u_n_p_1 point to _u_n_m_1 now
        // do not swap just copy
        Math::copy<T>(Constants::DOF, _u_n_p_1, 1, _u_n, 1);
    }

    /**
     * \brief calculate the _residual(Newton Function) + _residualNorm
     *
     * _residual = _A*x[n+1] - _B * x[n] + _M_h2 * x[n-1] - forces
     * _residualNorm = norm(_residual)
     * \param pointer to forces vector size of DOF
     */
    void CalculateResidual(T* force) {
        // _residual = _A*x[n+1]
        Math::gemm<T>(CblasRowMajor, CblasNoTrans, CblasNoTrans, Constants::DOF, 1, Constants::DOF, 1, _A, Constants::DOF, _u_n_p_1, 1, 0, _residual, 1);
        // _residual -= _B*x[n]
        Math::gemm<T>(CblasRowMajor, CblasNoTrans, CblasNoTrans, Constants::DOF, 1, Constants::DOF, -1, _B, Constants::DOF, _u_n, 1, 1, _residual, 1);
        // _residual += _M_h2 * x[n-1]
        Math::gemm<T>(CblasRowMajor, CblasNoTrans, CblasNoTrans, Constants::DOF, 1, Constants::DOF, 1, _M_h2, Constants::DOF, _u_n_m_1, 1, 1, _residual, 1);
        // _residual -= force
        Math::axpy<T>(Constants::DOF, -1, force, 1, _residual, 1);
        // res = norm(_residual)
        _residualNorm = Math::nrm2<T>(Constants::DOF, _residual, 1);
    }

    /**
     * \brief construct Jacobian for non fixed case
     *
     * this has to be called every newton iteraton
     * J = M_h2 + D / h + K + dKdx*x[n+1] + 1/h * dDdx ( x[n+1] - x[n] )
     */
    void ConstructJacobian() {
        // first update the derivative
        auto& db = MetaDatabase<T>::getDatabase();
        db.getLookupStiffness().getDerivative(_car->getCurrentSpringsLengths(), _dkdl);
        // construct the derivative (tensor) times a pos vector
        ConstructLookupDerivativeX(_dkdl, _u_n_p_1, _dKdxx);
        // J = A = M_h2 + D / h + K
        Math::copy<T>(Constants::DOFDOF, _A, 1, _J, 1);
        // _J += dKdx * x[n+1]
        Math::axpy<T>(Constants::DOFDOF, 1, _dKdxx, 1, _J, 1);
		(this->*_JacobianAdjustment)();
#ifdef DAMPING
                db.getLookupDamping().getDerivative(_car->getCurrentSpringsLengths(), _dddl);
		// temp = x[n+1]
        Math::copy<T>(Constants::DOF, _u_n_p_1, 1, _temp, 1);
        // temp += -x[n]
        Math::axpy<T>(Constants::DOF, -1, _u_n, 1, _temp, 1);
        // calc _dDdxx with (x[n+1] - x[n])
        ConstructLookupDerivativeX(_dddl, _temp, _dDdxx);
        // _J += 1/_h * _dDdxx
        Math::axpy<T>(Constants::DOFDOF, _factor_h, _dDdxx, 1, _J, 1);
#endif
    }
    /**
     * \brief construct Jacobian for fixed to road
     *
     * add the rows according to the tyres of [-(1/h dDdx + dKdx) x[n+1] - 1/h D - K + 1/h dDdx * x[n]] 
     * therefore J = M_h2 for the tyre positions
     */
    void ConstructFixedJacobian() {
        for (auto i = 0; i < Constants::NUM_LEGS; i++) {
            SetJacobianTyreLineToFixed(i);
        }
    }

    /**
     * \brief set a line in the jacobian according to a tyre to the fixed version
     *
     * \param[in] tyreNumber fl: 0, fr: 1, rl: 2, rr: 3
     */
    void SetJacobianTyreLineToFixed(int tyreNumber) {
        Math::copy<T>(Constants::DOF, _M_h2 + Constants::TYRE_INDEX_EULER[tyreNumber] * Constants::DOF, 1, _J + Constants::TYRE_INDEX_EULER[tyreNumber] * Constants::DOF, 1);
    }

    /**
     * \brief construct Jacobian for that one time step when the tyre hits the road
     *
     * _J = 1/h (_D + dD/dx (x[n+1]-x[n]) for the tyre positions
     */
    void ConstructTransitionJacobian() {
        // _matrixTmp = _D
        Math::copy(Constants::DOFDOF, _D, 1, _matrixTmp, 1);
        // Jacobian has already been called -> _dDdxx wit (x[n+1]-x[n]) can be directly accessed
        // _matrixTmp = _D + dD/dx (x[n+1]-x[n])
        Math::axpy(Constants::DOFDOF, 1, _dDdxx, 1, _matrixTmp, 1);
        // _matrixTmp /= h
        Math::scal(Constants::DOFDOF, fact_h, _matrixTmp, 1);
        for (auto i = 0; i < Constants::NUM_LEGS; i++) {
            SetJacobianTyreLineToTransition(i, M_temp);
        }
    }

    /**
     * \brief set a line in the jacobian according to a tyre to the transition versino
     *
     * \param[in] tyreNumber fl: 0, fr: 1, rl: 2, rr: 3
     */
    void SetJacobianTyreLineToTransition(int tyreNumber) {
        // _matrixTmp = _D
        Math::copy(Constants::DOFDOF, _D, 1, _matrixTmp, 1);
        // Jacobian has already been called -> _dDdxx wit (x[n+1]-x[n]) can be directly accessed
        // _matrixTmp = _D + dD/dx (x[n+1]-x[n])
        Math::axpy(Constants::DOFDOF, 1, _dDdxx, 1, _matrixTmp, 1);
        // _matrixTmp /= h
        Math::scal(Constants::DOFDOF, fact_h, _matrixTmp, 1);
        // write line into Jacobien
		SetJacobianTyreLineToTransition(tyreNumber, _matrixTmp);
    }
    /**
     * \brief set a line in the jacobian according to a tyre to the transition versino
     *
     * \param[in] tyreNumber fl: 0, fr: 1, rl: 2, rr: 3
     * \param[in] JNew precalculated Jacobian
     */
    void SetJacobianTyreLineToTransition(int tyreNumber, T* JNew) {
        Math::copy<T>(Constants::DOF, JNew + Constants::TYRE_INDEX_EULER[i] * Constants::DOF, 1, _J + Constants::TYRE_INDEX_EULER[i] * Constants::DOF, 1);
    }

	/**
	 * \brief construct Jacobian for non-fixed to road
	 *
	 * No Modification in the Jacobian*
	 */
	void ConstructNonFixedJacobian() {}

    /**
     * \brief update all dependent matrices on the position vector
     *
     * only needed in case of lookup Tables
     */
    void UpdateSystem() {
        _car->UpdateLengthsTwoTrackModel();
		_loadModuleObj->GetEulerianForce(_currentIter, _twoTrackModelForce);
        ConstructStiffnessMatrix();
#ifdef DAMPING
        ConstructDampingMatrix();
#endif
		ConstructAMatrix();
        ConstructBMatrix();
    }
    /**
     * \brief construct Stiffness Matrix
     *
     * this has to be called every newton iteraton and the spring lengtvector
     * has to be updated before
     */
    void ConstructStiffnessMatrix() {
#ifdef INTERPOLATION
        MetaDatabase<T>::getDatabase().getLookupStiffness().getInterpolation(_car->getCurrentSpringsLengths(), _kVec);
#endif
        _temp[0] = _kVec[0] + _kVec[2] + _kVec[4] + _kVec[6];
        _K[1] = _kVec[0] * _car->_lenLat[0] - _kVec[2] * _car->_lenLat[1] + _kVec[4] * _car->_lenLat[2] - _kVec[6] * _car->_lenLat[3];
        _K[2] = -_kVec[0] * _car->_lenLong[0] - _kVec[2] * _car->_lenLong[1] + _kVec[4] * _car->_lenLong[2] + _kVec[6] * _car->_lenLong[3];
        _K[3] = -_kVec[0];
        _K[4] = 0;
        _K[5] = -_kVec[2];
        _K[6] = 0;
        _K[7] = -_kVec[4];
        _K[8] = 0;
        _K[9] = -_kVec[6];
        _K[10] = 0;

        _temp[1] = _car->_lenLat[0] * _car->_lenLat[0] * _kVec[0] + _car->_lenLat[1] * _car->_lenLat[1] * _kVec[2] + _car->_lenLat[2] * _car->_lenLat[2] * _kVec[4] + _car->_lenLat[3] * _car->_lenLat[3] * _kVec[6];
        _K[Constants::DOF + 2] = -_car->_lenLong[0] * _car->_lenLat[0] * _kVec[0] + _car->_lenLat[1] * _car->_lenLong[1] * _kVec[2] + _car->_lenLong[2] * _car->_lenLat[2] * _kVec[4] - _car->_lenLong[3] * _car->_lenLat[3] * _kVec[6];
        _K[Constants::DOF + 3] = -_car->_lenLat[0] * _kVec[0];
        _K[Constants::DOF + 4] = 0;
        _K[Constants::DOF + 5] = _car->_lenLat[1] * _kVec[2];
        _K[Constants::DOF + 6] = 0;
        _K[Constants::DOF + 7] = -_car->_lenLat[2] * _kVec[4];
        _K[Constants::DOF + 8] = 0;
        _K[Constants::DOF + 9] = _car->_lenLat[3] * _kVec[6];
        _K[Constants::DOF + 10] = 0;

        _temp[2] = _car->_lenLong[0] * _car->_lenLong[0] * _kVec[0] + _car->_lenLong[1] * _car->_lenLong[1] * _kVec[2] + _car->_lenLong[2] * _car->_lenLong[2] * _kVec[4] + _car->_lenLong[3] * _car->_lenLong[3] * _kVec[6];
        _K[2 * Constants::DOF + 3] = _car->_lenLong[0] * _kVec[0];
        _K[2 * Constants::DOF + 4] = 0;
        _K[2 * Constants::DOF + 5] = _car->_lenLong[1] * _kVec[2];
        _K[2 * Constants::DOF + 6] = 0;
        _K[2 * Constants::DOF + 7] = -_car->_lenLong[2] * _kVec[4];
        _K[2 * Constants::DOF + 8] = 0;
        _K[2 * Constants::DOF + 9] = -_car->_lenLong[3] * _kVec[6];
        _K[2 * Constants::DOF + 10] = 0;

        _temp[3] = _kVec[0] + _kVec[1];
        _K[3 * Constants::DOF + 4] = -_kVec[1];
        _K[3 * Constants::DOF + 5] = 0;
        _K[3 * Constants::DOF + 6] = 0;
        _K[3 * Constants::DOF + 7] = 0;
        _K[3 * Constants::DOF + 8] = 0;
        _K[3 * Constants::DOF + 9] = 0;
        _K[3 * Constants::DOF + 10] = 0;
        // all others are zero

        _temp[4] = _kVec[1];
        _K[4 * Constants::DOF + 5] = 0;
        _K[4 * Constants::DOF + 6] = 0;
        _K[4 * Constants::DOF + 7] = 0;
        _K[4 * Constants::DOF + 8] = 0;
        _K[4 * Constants::DOF + 9] = 0;
        _K[4 * Constants::DOF + 10] = 0;

        _temp[5] = _kVec[2] + _kVec[3];
        _K[5 * Constants::DOF + 6] = -_kVec[3];
        _K[5 * Constants::DOF + 7] = 0;
        _K[5 * Constants::DOF + 8] = 0;
        _K[5 * Constants::DOF + 9] = 0;
        _K[5 * Constants::DOF + 10] = 0;

        _temp[6] = _kVec[3];
        _K[6 * Constants::DOF + 7] = 0;
        _K[6 * Constants::DOF + 8] = 0;
        _K[6 * Constants::DOF + 9] = 0;
        _K[6 * Constants::DOF + 10] = 0;

        _temp[7] = _kVec[4] + _kVec[5];
        _K[7 * Constants::DOF + 8] = -_kVec[5];
        _K[7 * Constants::DOF + 9] = 0;
        _K[7 * Constants::DOF + 10] = 0;

        _temp[8] = _kVec[5];
        _K[8 * Constants::DOF + 9] = 0;
        _K[8 * Constants::DOF + 10] = 0;

        _temp[9] = _kVec[6] + _kVec[7];
        _K[9 * Constants::DOF + 10] = -_kVec[7];

        _temp[10] = _kVec[7];

        // symmetrize _K
        // cblas_dcopy(DOF * DOF, _K, 1, K_trans, 1);
        Math::lacpy<T>(LAPACK_ROW_MAJOR, 'U', Constants::DOF, Constants::DOF, _K, Constants::DOF, _matrixTmp, Constants::DOF);

        Math::imatcopy<T>('R', 'T', Constants::DOF, Constants::DOF, 1, _matrixTmp, Constants::DOF,
                          Constants::DOF);  // get transpose of matrix

        Math::lacpy<T>(LAPACK_ROW_MAJOR, 'L', Constants::DOF, Constants::DOF, _matrixTmp, Constants::DOF, _K,
                       Constants::DOF);  // copy lower triangular in the orig matrix
        // cblas_daxpy(DOF * DOF,1, K_trans, 1, _K, 1); // _K = _K + _K'

        // add the diagonal to _K
        Math::CopyToDiagonal<T>(_K, _temp, Constants::DOF);  // _K = _K + _K'+ diag(_K)
    }
    /**
     * \brief construct Damping Matrix
     *
     * this has to be called every newton iteraton and the spring lengtvector
     * has to be updated before
     */
    void ConstructDampingMatrix() {
#ifdef INTERPOLATION
        MetaDatabase<T>::getDatabase().getLookupDamping().getInterpolation(_car->getCurrentSpringsLengths(), _dVec);
#endif
        _temp[0] = _dVec[0] + _dVec[2] + _dVec[4] + _dVec[6];
        _D[1] = _dVec[0] * _car->_lenLat[0] - _dVec[2] * _car->_lenLat[1] + _dVec[4] * _car->_lenLat[2] - _dVec[6] * _car->_lenLat[3];
        _D[2] = -_dVec[0] * _car->_lenLong[0] - _dVec[2] * _car->_lenLong[1] + _dVec[4] * _car->_lenLong[2] + _dVec[6] * _car->_lenLong[3];
        _D[3] = -_dVec[0];
        _D[4] = 0;
        _D[5] = -_dVec[2];
        _D[6] = 0;
        _D[7] = -_dVec[4];
        _D[8] = 0;
        _D[9] = -_dVec[6];
        _D[10] = 0;

        _temp[1] = _car->_lenLat[0] * _car->_lenLat[0] * _dVec[0] + _car->_lenLat[1] * _car->_lenLat[1] * _dVec[2] + _car->_lenLat[2] * _car->_lenLat[2] * _dVec[4] + _car->_lenLat[3] * _car->_lenLat[3] * _dVec[6];
        _D[Constants::DOF + 2] = -_car->_lenLong[0] * _car->_lenLat[0] * _dVec[0] + _car->_lenLat[1] * _car->_lenLong[1] * _dVec[2] + _car->_lenLong[2] * _car->_lenLat[2] * _dVec[4] - _car->_lenLong[3] * _car->_lenLat[3] * _dVec[6];
        _D[Constants::DOF + 3] = -_car->_lenLat[0] * _dVec[0];
        _D[Constants::DOF + 4] = 0;
        _D[Constants::DOF + 5] = _car->_lenLat[1] * _dVec[2];
        _D[Constants::DOF + 6] = 0;
        _D[Constants::DOF + 7] = -_car->_lenLat[2] * _dVec[4];
        _D[Constants::DOF + 8] = 0;
        _D[Constants::DOF + 9] = _car->_lenLat[3] * _dVec[6];
        _D[Constants::DOF + 10] = 0;

        _temp[2] = _car->_lenLong[0] * _car->_lenLong[0] * _dVec[0] + _car->_lenLong[1] * _car->_lenLong[1] * _dVec[2] + _car->_lenLong[2] * _car->_lenLong[2] * _dVec[4] + _car->_lenLong[3] * _car->_lenLong[3] * _dVec[6];
        _D[2 * Constants::DOF + 3] = _car->_lenLong[0] * _dVec[0];
        _D[2 * Constants::DOF + 4] = 0;
        _D[2 * Constants::DOF + 5] = _car->_lenLong[1] * _dVec[2];
        _D[2 * Constants::DOF + 6] = 0;
        _D[2 * Constants::DOF + 7] = -_car->_lenLong[2] * _dVec[4];
        _D[2 * Constants::DOF + 8] = 0;
        _D[2 * Constants::DOF + 9] = -_car->_lenLong[3] * _dVec[6];
        _D[2 * Constants::DOF + 10] = 0;

        _temp[3] = _dVec[0] + _dVec[1];
        _D[3 * Constants::DOF + 4] = -_dVec[1];
        _D[3 * Constants::DOF + 5] = 0;
        _D[3 * Constants::DOF + 6] = 0;
        _D[3 * Constants::DOF + 7] = 0;
        _D[3 * Constants::DOF + 8] = 0;
        _D[3 * Constants::DOF + 9] = 0;
        _D[3 * Constants::DOF + 10] = 0;
        // all others are zero

        _temp[4] = _dVec[1];
        _D[4 * Constants::DOF + 5] = 0;
        _D[4 * Constants::DOF + 6] = 0;
        _D[4 * Constants::DOF + 7] = 0;
        _D[4 * Constants::DOF + 8] = 0;
        _D[4 * Constants::DOF + 9] = 0;
        _D[4 * Constants::DOF + 10] = 0;

        _temp[5] = _dVec[2] + _dVec[3];
        _D[5 * Constants::DOF + 6] = -_dVec[3];
        _D[5 * Constants::DOF + 7] = 0;
        _D[5 * Constants::DOF + 8] = 0;
        _D[5 * Constants::DOF + 9] = 0;
        _D[5 * Constants::DOF + 10] = 0;

        _temp[6] = _dVec[3];
        _D[6 * Constants::DOF + 7] = 0;
        _D[6 * Constants::DOF + 8] = 0;
        _D[6 * Constants::DOF + 9] = 0;
        _D[6 * Constants::DOF + 10] = 0;

        _temp[7] = _dVec[4] + _dVec[5];
        _D[7 * Constants::DOF + 8] = -_dVec[5];
        _D[7 * Constants::DOF + 9] = 0;
        _D[7 * Constants::DOF + 10] = 0;

        _temp[8] = _dVec[5];
        _D[8 * Constants::DOF + 9] = 0;
        _D[8 * Constants::DOF + 10] = 0;

        _temp[9] = _dVec[6] + _dVec[7];
        _D[9 * Constants::DOF + 10] = -_dVec[7];

        _temp[10] = _dVec[7];

        // symmetrize _D
        // cblas_dcopy(DOF * DOF, _D, 1, D_trans, 1);
        Math::lacpy<T>(LAPACK_ROW_MAJOR, 'U', Constants::DOF, Constants::DOF, _D, Constants::DOF, _matrixTmp, Constants::DOF);

        Math::imatcopy<T>('R', 'T', Constants::DOF, Constants::DOF, 1, _matrixTmp, Constants::DOF,
                          Constants::DOF);  // get transpose of matrix

        Math::lacpy<T>(LAPACK_ROW_MAJOR, 'L', Constants::DOF, Constants::DOF, _matrixTmp, Constants::DOF, _D,
                       Constants::DOF);  // copy lower triangular in the orig matrix
        // cblas_daxpy(DOF * DOF,1, D_trans, 1, _D, 1); // _D = _D + _D'

        // add the diagonal to _K
        Math::CopyToDiagonal<T>(_D, _temp, Constants::DOF);  // _D = _D + _D'+ diag(_D)
    }

    /**
     * \brief construct The derivative of the stiffness lookupTable times a
     * position vector
     *
     * this has to be called every newton iteraton
     *
     * \param[in] der pointer to vector with the derivative of the lookup Table
     * of length 8 \param[in] x pointer to position vector of length DOF
     * \param[out] dMdxx pointer to matrix in which de derivative of the lookup
     * times a position vec is stored of size DOF * DOF.
     */
    void ConstructLookupDerivativeX(T* der, T* x, T* dMdxx) {
        _temp[0] = x[0] * (der[0] + der[2] + der[4] + der[6]) - der[2] * x[5] - der[4] * x[7] - der[6] * x[9] - der[0] * x[3] + x[1] * (der[0] * _car->_lenLat[0] - der[2] * _car->_lenLat[1] + der[4] * _car->_lenLat[2] - der[6] * _car->_lenLat[3]) - x[2] * (der[0] * _car->_lenLong[0] + der[2] * _car->_lenLong[1] - der[4] * _car->_lenLong[2] - der[6] * _car->_lenLong[3]);
        dMdxx[1] = der[2] * _car->_lenLat[1] * _car->_lenLat[1] * x[1] + der[4] * _car->_lenLat[2] * _car->_lenLat[2] * x[1] + der[6] * _car->_lenLat[3] * _car->_lenLat[3] * x[1] + der[0] * _car->_lenLat[0] * (x[0] - x[3] + _car->_lenLat[0] * x[1]) - der[2] * _car->_lenLat[1] * x[0] + der[2] * _car->_lenLat[1] * x[5] + der[4] * _car->_lenLat[2] * x[0] - der[4] * _car->_lenLat[2] * x[7] - der[6] * _car->_lenLat[3] * x[0] + der[6] * _car->_lenLat[3] * x[9] - der[0] * _car->_lenLat[0] * _car->_lenLong[0] * x[2] + der[2] * _car->_lenLat[1] * _car->_lenLong[1] * x[2] + der[4] * _car->_lenLat[2] * _car->_lenLong[2] * x[2] - der[6] * _car->_lenLat[3] * _car->_lenLong[3] * x[2];
        dMdxx[2] = der[0] * _car->_lenLong[0] * _car->_lenLong[0] * x[2] + der[4] * _car->_lenLong[2] * _car->_lenLong[2] * x[2] + der[6] * _car->_lenLong[3] * _car->_lenLong[3] * x[2] - der[0] * _car->_lenLong[0] * (x[0] - x[3] + _car->_lenLat[0] * x[1]) + der[2] * _car->_lenLong[1] * (x[5] - x[0] + _car->_lenLat[1] * x[1] + _car->_lenLong[1] * x[2]) + der[4] * _car->_lenLong[2] * x[0] - der[4] * _car->_lenLong[2] * x[7] + der[6] * _car->_lenLong[3] * x[0] - der[6] * _car->_lenLong[3] * x[9] + der[4] * _car->_lenLat[2] * _car->_lenLong[2] * x[1] - der[6] * _car->_lenLat[3] * _car->_lenLong[3] * x[1];
        dMdxx[3] = -der[0] * (x[0] - x[3] + _car->_lenLat[0] * x[1] - _car->_lenLong[0] * x[2]);
        dMdxx[4] = 0;
        dMdxx[5] = der[2] * (x[5] - x[0] + _car->_lenLat[1] * x[1] + _car->_lenLong[1] * x[2]);
        dMdxx[6] = 0;
        dMdxx[7] = -der[4] * (x[0] - x[7] + _car->_lenLat[2] * x[1] + _car->_lenLong[2] * x[2]);
        dMdxx[8] = 0;
        dMdxx[9] = -der[6] * (x[0] - x[9] - _car->_lenLat[3] * x[1] + _car->_lenLong[3] * x[2]);
        dMdxx[10] = 0;

        _temp[1] = der[2] * _car->_lenLat[1] * _car->_lenLat[1] * x[0] - der[2] * _car->_lenLat[1] * _car->_lenLat[1] * _car->_lenLat[1] * x[1] - der[2] * _car->_lenLat[1] * _car->_lenLat[1] * x[5] + der[4] * _car->_lenLat[2] * _car->_lenLat[2] * x[0] + der[4] * _car->_lenLat[2] * _car->_lenLat[2] * _car->_lenLat[2] * x[1] - der[4] * _car->_lenLat[2] * _car->_lenLat[2] * x[7] + der[6] * _car->_lenLat[3] * _car->_lenLat[3] * x[0] - der[6] * _car->_lenLat[3] * _car->_lenLat[3] * _car->_lenLat[3] * x[1] - der[6] * _car->_lenLat[3] * _car->_lenLat[3] * x[9] + der[0] * _car->_lenLat[0] * _car->_lenLat[0] * (x[0] - x[3] + _car->_lenLat[0] * x[1]) - der[0] * _car->_lenLat[0] * _car->_lenLat[0] * _car->_lenLong[0] * x[2] - der[2] * _car->_lenLat[1] * _car->_lenLat[1] * _car->_lenLong[1] * x[2] + der[4] * _car->_lenLat[2] * _car->_lenLat[2] * _car->_lenLong[2] * x[2] + der[6] * _car->_lenLat[3] * _car->_lenLat[3] * _car->_lenLong[3] * x[2];
        dMdxx[Constants::DOF + 2] = der[2] * _car->_lenLat[1] * _car->_lenLong[1] * x[0] - der[2] * _car->_lenLat[1] * _car->_lenLong[1] * x[5] + der[4] * _car->_lenLat[2] * _car->_lenLong[2] * x[0] - der[4] * _car->_lenLat[2] * _car->_lenLong[2] * x[7] - der[6] * _car->_lenLat[3] * _car->_lenLong[3] * x[0] + der[6] * _car->_lenLat[3] * _car->_lenLong[3] * x[9] + der[0] * _car->_lenLat[0] * _car->_lenLong[0] * _car->_lenLong[0] * x[2] - der[2] * _car->_lenLat[1] * _car->_lenLat[1] * _car->_lenLong[1] * x[1] - der[2] * _car->_lenLat[1] * _car->_lenLong[1] * _car->_lenLong[1] * x[2] + der[4] * _car->_lenLat[2] * _car->_lenLat[2] * _car->_lenLong[2] * x[1] + der[4] * _car->_lenLat[2] * _car->_lenLong[2] * _car->_lenLong[2] * x[2] + der[6] * _car->_lenLat[3] * _car->_lenLat[3] * _car->_lenLong[3] * x[1] - der[6] * _car->_lenLat[3] * _car->_lenLong[3] * _car->_lenLong[3] * x[2] - der[0] * _car->_lenLat[0] * _car->_lenLong[0] * (x[0] - x[3] + _car->_lenLat[0] * x[1]);
        dMdxx[Constants::DOF + 3] = -der[0] * _car->_lenLat[0] * (x[0] - x[3] + _car->_lenLat[0] * x[1] - _car->_lenLong[0] * x[2]);
        dMdxx[Constants::DOF + 4] = 0;
        dMdxx[Constants::DOF + 5] = -der[2] * _car->_lenLat[1] * (x[5] - x[0] + _car->_lenLat[1] * x[1] + _car->_lenLong[1] * x[2]);
        dMdxx[Constants::DOF + 6] = 0;
        dMdxx[Constants::DOF + 7] = -der[4] * _car->_lenLat[2] * (x[0] - x[7] + _car->_lenLat[2] * x[1] + _car->_lenLong[2] * x[2]);
        dMdxx[Constants::DOF + 8] = 0;
        dMdxx[Constants::DOF + 9] = der[6] * _car->_lenLat[3] * (x[0] - x[9] - _car->_lenLat[3] * x[1] + _car->_lenLong[3] * x[2]);
        dMdxx[Constants::DOF + 10] = 0;

        _temp[2] = der[2] * _car->_lenLong[1] * _car->_lenLong[1] * x[0] - der[0] * _car->_lenLong[0] * _car->_lenLong[0] * _car->_lenLong[0] * x[2] - der[2] * _car->_lenLong[1] * _car->_lenLong[1] * _car->_lenLong[1] * x[2] - der[2] * _car->_lenLong[1] * _car->_lenLong[1] * x[5] + der[4] * _car->_lenLong[2] * _car->_lenLong[2] * x[0] + der[4] * _car->_lenLong[2] * _car->_lenLong[2] * _car->_lenLong[2] * x[2] - der[4] * _car->_lenLong[2] * _car->_lenLong[2] * x[7] + der[6] * _car->_lenLong[3] * _car->_lenLong[3] * x[0] + der[6] * _car->_lenLong[3] * _car->_lenLong[3] * _car->_lenLong[3] * x[2] - der[6] * _car->_lenLong[3] * _car->_lenLong[3] * x[9] + der[0] * _car->_lenLong[0] * _car->_lenLong[0] * (x[0] - x[3] + _car->_lenLat[0] * x[1]) - der[2] * _car->_lenLat[1] * _car->_lenLong[1] * _car->_lenLong[1] * x[1] + der[4] * _car->_lenLat[2] * _car->_lenLong[2] * _car->_lenLong[2] * x[1] - der[6] * _car->_lenLat[3] * _car->_lenLong[3] * _car->_lenLong[3] * x[1];
        dMdxx[2 * Constants::DOF + 3] = der[0] * _car->_lenLong[0] * (x[0] - x[3] + _car->_lenLat[0] * x[1] - _car->_lenLong[0] * x[2]);
        dMdxx[2 * Constants::DOF + 4] = 0;
        dMdxx[2 * Constants::DOF + 5] = -der[2] * _car->_lenLong[1] * (x[5] - x[0] + _car->_lenLat[1] * x[1] + _car->_lenLong[1] * x[2]);
        dMdxx[2 * Constants::DOF + 6] = 0;
        dMdxx[2 * Constants::DOF + 7] = -der[2] * _car->_lenLong[1] * (x[5] - x[0] + _car->_lenLat[1] * x[1] + _car->_lenLong[1] * x[2]);
        dMdxx[2 * Constants::DOF + 8] = 0;
        dMdxx[2 * Constants::DOF + 9] = -der[2] * _car->_lenLong[1] * (x[5] - x[0] + _car->_lenLat[1] * x[1] + _car->_lenLong[1] * x[2]);
        dMdxx[2 * Constants::DOF + 10] = 0;

        _temp[3] = der[1] * (x[3] - x[4]) + der[0] * (x[0] - x[3] + _car->_lenLat[0] * x[1]) - der[0] * _car->_lenLong[0] * x[2];
        dMdxx[3 * Constants::DOF + 4] = -der[1] * (x[3] - x[4]);
        dMdxx[3 * Constants::DOF + 5] = 0;
        dMdxx[3 * Constants::DOF + 6] = 0;
        dMdxx[3 * Constants::DOF + 7] = 0;
        dMdxx[3 * Constants::DOF + 8] = 0;
        dMdxx[3 * Constants::DOF + 9] = 0;
        dMdxx[3 * Constants::DOF + 10] = 0;
        // all others are zero

        _temp[4] = der[1] * (x[3] - x[4]);
        dMdxx[4 * Constants::DOF + 5] = 0;
        dMdxx[4 * Constants::DOF + 6] = 0;
        dMdxx[4 * Constants::DOF + 7] = 0;
        dMdxx[4 * Constants::DOF + 8] = 0;
        dMdxx[4 * Constants::DOF + 9] = 0;
        dMdxx[4 * Constants::DOF + 10] = 0;

        _temp[5] = der[3] * (x[5] - x[6]) - der[2] * (x[5] - x[0] + _car->_lenLat[1] * x[1] + _car->_lenLong[1] * x[2]);
        dMdxx[5 * Constants::DOF + 6] = -der[3] * (x[5] - x[6]);
        dMdxx[5 * Constants::DOF + 7] = 0;
        dMdxx[5 * Constants::DOF + 8] = 0;
        dMdxx[5 * Constants::DOF + 9] = 0;
        dMdxx[5 * Constants::DOF + 10] = 0;

        _temp[6] = der[3] * (x[5] - x[6]);
        dMdxx[6 * Constants::DOF + 7] = 0;
        dMdxx[6 * Constants::DOF + 8] = 0;
        dMdxx[6 * Constants::DOF + 9] = 0;
        dMdxx[6 * Constants::DOF + 10] = 0;

        _temp[7] = der[5] * (x[7] - x[8]) + der[4] * (x[0] - x[7] + _car->_lenLat[2] * x[1] + _car->_lenLong[2] * x[2]);
        dMdxx[7 * Constants::DOF + 8] = -der[5] * (x[7] - x[8]);
        dMdxx[7 * Constants::DOF + 9] = 0;
        dMdxx[7 * Constants::DOF + 10] = 0;

        _temp[8] = der[5] * (x[7] - x[8]);
        dMdxx[8 * Constants::DOF + 9] = 0;
        dMdxx[8 * Constants::DOF + 10] = 0;

        _temp[9] = der[7] * (x[9] - x[10]) + der[6] * (x[0] - x[9] - _car->_lenLat[3] * x[1] + _car->_lenLong[3] * x[2]);
        dMdxx[9 * Constants::DOF + 10] = -der[7] * (x[9] - x[10]);

        _temp[10] = der[7] * (x[9] - x[10]);

        // symmetrize dMdxx
        // cblas_dcopy(DOF * DOF, dMdxx, 1, dMdxx_trans, 1);
        Math::lacpy<T>(LAPACK_ROW_MAJOR, 'U', Constants::DOF, Constants::DOF, dMdxx, Constants::DOF, _matrixTmp, Constants::DOF);

        Math::imatcopy<T>('R', 'T', Constants::DOF, Constants::DOF, 1, _matrixTmp, Constants::DOF,
                          Constants::DOF);  // get transpose of matrix

        Math::lacpy<T>(LAPACK_ROW_MAJOR, 'L', Constants::DOF, Constants::DOF, _matrixTmp, Constants::DOF, dMdxx,
                       Constants::DOF);  // copy lower triangular in the orig matrix
        // cblas_daxpy(DOF * DOF,1, dMdxx_trans, 1, dMdxx, 1); // dMdxx = dMdxx
        // + dMdxx'

        // add the diagonal to dM
        Math::CopyToDiagonal<T>(dMdxx, _temp, Constants::DOF);  // dMdxx = dMdxx + dMdxx'+ diag(dMdxx)
    }
};

template <typename T>
class TwoTrackModelBDF2 : public TwoTrackModelBE<T> {
private:
    T *_DMat, *_E;
    /** this vector is used to store the whole constants part on the rhs apart
     * from the force */
    T* _bVec;
    T *_u_n_m_2, *_u_n_m_3;
    size_t _timeStepCount = 0;
    void (TwoTrackModelBDF2<T>::*_activeExecutor)(T, T*);

    /**
     * \brief construct _A
     *
     * _A = 9/4h^2 * M + 3/2h * _D + _K
     */
    void ConstructAMatrix() {
        // _A = M/h^2
        Math::copy<T>(Constants::DOFDOF, _M_h2, 1, _A, 1);
        // _A *= 9/4
        Math::scal<T>(Constants::DOFDOF, 2.25, _A, 1);
        // _A += 3/2h * _D
        Math::axpy<T>(Constants::DOFDOF, 1.5 * _factor_h, _D, 1, _A, 1);
        // _A += _K
        Math::axpy<T>(Constants::DOFDOF, 1, _K, 1, _A, 1);
    }
    /**
     * \brief construct _bVec
     *
     * this has to be called every time step for interpolation and is only used
     * in case of interpolation _bVec = 1/h^2 * M * (6 * x[n] - 11/2 * x[n-1] + 2
     * * x[n-2] - 1/4 * x[n-3]) + 1/h
     * * _D * (2 * x[n] - 1/2 * x[n-1])
     */
    void ConstructbVec() {
        // temp = x[n]
        Math::copy<T>(Constants::DOF, _u_n, 1, _temp, 1);
        // temp *= 6
        Math::scal<T>(Constants::DOF, 6, _temp, 1);
        // temp += - 11/2 * x[n-1]
        Math::axpy<T>(Constants::DOF, -5.5, _u_n_m_1, 1, _temp, 1);
        // temp += 2 * x[n-2]
        Math::axpy<T>(Constants::DOF, 2, _u_n_m_2, 1, _temp, 1);
        // temp += -1/4 * x[n-3]
        Math::axpy<T>(Constants::DOF, -0.25, _u_n_m_3, 1, _temp, 1);
        // _bVec = _M_h2 * temp
        Math::gemv<T>(CblasRowMajor, CblasNoTrans, Constants::DOF, Constants::DOF, 1, _M_h2, Constants::DOF, _temp, 1, 0, _bVec, 1);
        // temp = x[n]
        Math::copy<T>(Constants::DOF, _u_n, 1, _temp, 1);
        // temp *= 2
        Math::scal<T>(Constants::DOF, 2, _temp, 1);
        // temp += - 1/2 * x[n-1]
        Math::axpy<T>(Constants::DOF, -0.5, _u_n_m_1, 1, _temp, 1);
        // _bVec += 1/h * _D * temp
        Math::gemv<T>(CblasRowMajor, CblasNoTrans, Constants::DOF, Constants::DOF, _factor_h, _D, Constants::DOF, _temp, 1, 1, _bVec, 1);
    }

    /**
     * \brief calculates a first guess for x[n+1]
     */
    void GetInitialGuess(T* forces) {
        ConstructbVec();
        lapack_int status;
        status = Math::potrf<T>(LAPACK_ROW_MAJOR, 'L', Constants::DOF, _A, Constants::DOF);
        Math::potrfCheckStatus(status);
        Math::vAdd<T>(Constants::DOF, _bVec, forces, _u_n_p_1);
        // _u_n_p_1=_A\(b+f)
        Math::potrs<T>(LAPACK_ROW_MAJOR, 'L', Constants::DOF, 1, _A, Constants::DOF, _u_n_p_1, 1);
        // reconstruct _A
        ConstructAMatrix();
    }

	/**
   * \brief update the velocity vector
   *
   * v = (3x[n+1]-4x[n]+x[n-1])/(2h)
   */
	void UpdateVelocity() {

		// v = x[n+1]
		Math::copy<T>(Constants::DOF, _u_n_p_1, Constants::INCX, _car->_currentVelocityTwoTrackModel, Constants::INCX);
		// v *= 3/2;
		Math::scal<T>(Constants::DOF, 1.5, _car->_currentVelocityTwoTrackModel, Constants::INCX);
		// v -= 2x[n]
		Math::axpy<T>(Constants::DOF, -2, _u_n, Constants::INCX, _car->_currentVelocityTwoTrackModel, Constants::INCX);
		// v += 0.5*x[n-1]
		Math::axpy<T>(Constants::DOF, 0.5, _u_n_m_1, Constants::INCX, _car->_currentVelocityTwoTrackModel, Constants::INCX);
		// v /= h
		Math::scal<T>(Constants::DOF, _factor_h, _car->_currentVelocityTwoTrackModel, Constants::INCX);
	}


public:
    /**
     * Constructor.
     */
    TwoTrackModelBDF2(Car<T>* inputCar, LoadModule<T>* loadModel) : TwoTrackModelBE<T>(inputCar, loadModel) {
        _DMat = Math::calloc<T>(Constants::DOFDOF);
        _E = Math::calloc<T>(Constants::DOFDOF);
        _u_n_m_2 = Math::malloc<T>(Constants::DOF);
        _u_n_m_3 = Math::malloc<T>(Constants::DOF);
        _bVec = Math::malloc<T>(Constants::DOF);
        _activeExecutor = &TwoTrackModelBDF2<T>::FirstTwoSteps;
    }

    void FirstTwoSteps(size_t currentIter, T* solution) {
        if (_timeStepCount == 0) {
            Math::copy<T>(Constants::DOF, _u_n_m_1, 1, _u_n_m_2, 1);
            TwoTrackModelBE<T>::UpdateStep(currentIter, solution);
        }
        else {
            Math::copy<T>(Constants::DOF, _u_n_m_2, 1, _u_n_m_3, 1);
            Math::copy<T>(Constants::DOF, _u_n_m_1, 1, _u_n_m_2, 1);
            TwoTrackModelBE<T>::UpdateStep(currentIter, solution);
            // construct _A
            ConstructAMatrix();
            ConstructAMatrix();
            _activeExecutor = &TwoTrackModelBDF2<T>::UpdateStepBDF2;
        }
    }
    virtual void UpdateStep(size_t currentIter, T* solution) {
		_currentIter = currentIter;
        (this->*_activeExecutor)(_currentIter, solution);
        _timeStepCount += 1;
    }

    /**
     * Performs one timestep of the 11DOF solver
     * \param load vector [angle:Z,GC:Y,W1:Y,T1:Y,W2:Y,T2:Y,...]
     * \return solution of the following timestep
     * [angle:Z,GC:Y,W1:Y,T1:Y,W2:Y,T2:Y,...]
     */
    void UpdateStepBDF2(size_t currentIter, T* solution) {
        // cblas_dscal(DOF,0, force, 1);
		// compute Force
		_loadModuleObj->GetEulerianForce(currentIter, _twoTrackModelForce);

        GetInitialGuess(_twoTrackModelForce);
#ifdef INTERPOLATION
        Math::Solvers<T, TwoTrackModelBDF2<T>>::Newton(this, _twoTrackModelForce, _J, _residual, &_residualNorm, _u_n_p_1, _temp, &_tolerance, &_maxNewtonIteration, &_newtonIteration);
#endif
        Math::copy<T>(Constants::DOF, _u_n_p_1, 1, solution, 1);
        // _u_n_m_2 points to _u_n_m_3 and _u_n_m_3 points to _u_n_m_2

		// update velocity inside car class
		UpdateVelocity();

        Math::SwapAddress<T>(_u_n_m_2, _u_n_m_3);
        // _u_n_m_2 points to _u_n_m_1 and _u_n_m_1 points to _u_n_m_3
        Math::SwapAddress<T>(_u_n_m_1, _u_n_m_2);
        // _u_n_m_1 points to _u_n and _u_n points to _u_n_m_3
        Math::SwapAddress<T>(_u_n, _u_n_m_1);
        Math::copy<T>(Constants::DOF, _u_n_p_1, 1, _u_n, 1);
    }
    /**
     * \brief calculate the _residual(Newton Function) + _residualNorm
     *
     * _residual = _A*x[n+1] - _bVec - forces
     * _residualNorm = norm(_residual)
     * \param force pointer to forces vector size of DOF
     */
    void CalculateResidual(T* force) {
        // _residual = _A*x[n+1]
        Math::gemm<T>(CblasRowMajor, CblasNoTrans, CblasNoTrans, Constants::DOF, 1, Constants::DOF, 1, _A, Constants::DOF, _u_n_p_1, 1, 0, _residual, 1);
        // _residual -= _bVec
        Math::axpy<T>(Constants::DOF, -1, _bVec, 1, _residual, 1);
        // _residual -= force
        Math::axpy<T>(Constants::DOF, -1, force, 1, _residual, 1);
        // res = norm(_residual)
        _residualNorm = Math::nrm2<T>(Constants::DOF, _residual, 1);
    }
    /**
     * \brief update all dependent matrices on the position vector
     *
     * only needed in case of lookup Tables
     */
    void UpdateSystem() {
        // constructSpringLengths();
        _car->UpdateLengthsTwoTrackModel();
		_loadModuleObj->GetEulerianForce(_currentIter, _twoTrackModelForce);
        ConstructStiffnessMatrix();
#ifdef DAMPING
        ConstructDampingMatrix();
#endif
        ConstructAMatrix();
        ConstructbVec();
    }
    /**
     * \brief construct Jacobian
     *
     * this has to be called every newton iteraton
     * _J = 9/4 * _M_h2 + 3/2h * _D + _K + dKdx*x[n+1] + 1/h * dDdx * (3/2 * x[n+1]
     * - 2 * x[n] + 1/2 * x[n-1])
     */
    void ConstructJacobian() {
        // first update the derivative
        auto& db = MetaDatabase<T>::getDatabase();
        db.getLookupStiffness().getDerivative(_car->getCurrentSpringsLengths(), _dkdl);
        ConstructLookupDerivativeX(_dkdl, _u_n_p_1, _dKdxx);
        // _J = _A
        Math::copy<T>(Constants::DOFDOF, _A, 1, _J, 1);
        // _J += dKdx * x[n+1]
        Math::axpy<T>(Constants::DOFDOF, 1, _dKdxx, 1, _J, 1);

#ifdef DAMPING
        db.getLookupDamping().getDerivative(_car->getCurrentSpringsLengths(), _dddl);
        // temp = x[n+1]
        Math::copy<T>(Constants::DOF, _u_n_p_1, 1, _temp, 1);
        // temp *= 3/2
        Math::scal<T>(Constants::DOF, 1.5, _temp, 1);
        // temp += -2 * x[n]
        Math::axpy<T>(Constants::DOF, -2, _u_n, 1, _temp, 1);
        // temp += 1/2 * x[n-1]
        Math::axpy<T>(Constants::DOF, 0.5, _u_n_m_1, 1, _temp, 1);
        // calc _dDdxx with (3/2 * x[n+1] - 2 * x[n] + 1/2 * x[n-1])
        ConstructLookupDerivativeX(_dddl, temp, _dDdxx);
        // _J += 1/_h * _dDdxx
        Math::axpy<T>(Constants::DOFDOF, _factor_h, _dDdxx, 1, _J, 1);
#endif
    }

    /**
     * \brief construct Jacobian for fixed to road
     *
     * add the rows according to the tyres of [-(1/h dDdx + dKdx) x[n+1] - 1/h _D
     * - _K + 1/h dDdx * x[n]] therefore _J = _M_h2 for the tyre positions
     */
    void ConstructFixedJacobian() {
        for (auto i = 0; i < Constants::NUM_LEGS; i++) {
            SetJacobianTyreLineToFixed(i);
        }
    }

    /**
     * \brief set a line in the jacobian according to a tyre to the fixed version
     *
     * \param[in] tyreNumber fl: 0, fr: 1, rl: 2, rr: 3
     */
    void SetJacobianTyreLineToFixed(int tyreNumber) {
        // _J = _M_h2
        Math::copy<T>(Constants::DOF, _M_h2 + Constants::TYRE_INDEX_EULER[i] * Constants::DOF, 1, _J + (4 + 2 * tyreNumber) * Constants::DOF, 1);
        // _J = 9/4 _M_h2
        Math::scal<T>(Constants::DOF, 2.25, _J + Constants::TYRE_INDEX_EULER[i] * Constants::DOF, 1);
    }

    /*
    Destructor
    */
    virtual ~TwoTrackModelBDF2() {
        Math::free<T>(_DMat);
        Math::free<T>(_E);
        Math::free<T>(_u_n_m_2);
        Math::free<T>(_u_n_m_3);
        Math::free<T>(_bVec);
    }
};

/**
 * For testing purposes.
 */
template <typename T>
class TwoTrackModelFull : public TwoTrackModelBE<T> {
public:
    TwoTrackModelFull(Car<T>* inputCar, LoadModule<T>* loadModel) : TwoTrackModelBE<T>(inputCar, loadModel) {
        auto& db = MetaDatabase<T>::getDatabase();
        _tend = db.getNumberOfTimeIterations() * _h;
        int sol_size = (floor(_tend / _h) + 1);
        _uSolution = Math::calloc<T>((sol_size + 1) * Constants::DOF);
    }

    void Solve(T* sol_vect) {
#ifdef WRITECSV
#ifdef INTERPOLATION
        IO::MyFile<T> newtonFile("C:\\software\\repos\\EVAA\\output\\newtonOutput.txt");
#endif // INTERPOLATION
#endif
        int iter = 1;
        T t = _h;
        double eps = _h / 100;
        T* solution_vect = _uSolution;
       
        while (std::abs(t - (_tend + _h)) > eps) {
            // solution_vect = _uSolution + iter * (DOF);
            UpdateStep(iter, sol_vect);
#ifdef WRITECSV
#ifdef INTERPOLATION
        newtonFile.writeSingleValue(_newtonIteration);
        newtonFile.writeSingleValue(_residualNorm);
        newtonFile.newLine();
#endif       
#endif // WRITECSV
            iter++;
            t += _h;
        }

        
    }

    void PrintFinalResults(T* sln) {
        std::cout.precision(15);
        std::cout << std::scientific;
        std::cout << "linear11DOF: orientation angles=\n\t[" << sln[1] << "\n\t " << sln[2] << "]" << std::endl;
        std::cout << "linear11DOF: car body position pc=\n\t[" << sln[0] << "]" << std::endl;
        std::cout << "linear11DOF: front-left wheel position pw3=\n\t[" << sln[3] << "]" << std::endl;
        std::cout << "linear11DOF: front-left tyre position pt3=\n\t[" << sln[4] << "]" << std::endl;
        std::cout << "linear11DOF: front-right wheel position pw4=\n\t[" << sln[5] << "]" << std::endl;
        std::cout << "linear11DOF: front-right tyre position pt4=\n\t[" << sln[6] << "]" << std::endl;
        std::cout << "linear11DOF: rear-left wheel position pw2=\n\t[" << sln[7] << "]" << std::endl;
        std::cout << "linear11DOF: rear-left tyre position pt2=\n\t[" << sln[8] << "]" << std::endl;
        std::cout << "linear11DOF: rear-right wheel position pw1=\n\t[" << sln[9] << "]" << std::endl;
        std::cout << "linear11DOF: rear-right tyre position pt1=\n\t[" << sln[10] << "]" << std::endl;
    }

    virtual ~TwoTrackModelFull() {
        Math::free<T>(_uSolution);
    }

private:
    T _tend;
    T *_uSolution;
};

}  // namespace EVAA
