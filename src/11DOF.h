/***********************************************************************************************//**
* \file 11DOF.h
* This file holds the function declaration and definitions of the linear11dof and the linear11dof_full class.
* \date 04/14/2020
**************************************************************************************************/
#pragma once
#include <mkl.h>
#include "car.h"
#include "MathLibrary.h"
#include "Constants.h"
#include "MetaDataBase.h"
#include "BLAS.h"

/**
* \brief class to compute one timestep of the linear 11 dof system in small angle approximation
*/
template <typename T>
class TwoTrackModelParent {
protected:
	// main car object
	Car<T>* _car;												/**< pointer to car instance with all important car parameter */

	// time step related
	T factor_h;													/**< solution in next timestep */
	T _h;														/**< solution in next timestep */

	// solution in next timestep
	T* u_n_p_1;

	// solution in previous timestep
	T* u_n_m_1;

	// solution in current timestep
	T* u_n;
	T* A, * B, *C;                                                  /**< pointer to matrices used for the backward euler */
	T* M_h2, * K, * D, * Mat_temp, * Mat_springLength;

	T* kVec, * dVec, * temp;
	T* springLengths, * springLengthsNormal;

	T* J, * dKdxx, * dDdxx;
	T* residual;
	T res_norm;
	T* dkdl, * dddl;

	// interpolation 
	EVAALookup<T>* lookupStiffness;
	EVAALookup<T>* lookupDamping;

public:
	TwoTrackModelParent(){}
	/*
	Constructor
	*/
	TwoTrackModelParent(Car<T>* input_car, EVAALookup<T>* stiffnessInterpolation, EVAALookup<T>* dampingInterpolation): lookupStiffness(stiffnessInterpolation), lookupDamping(dampingInterpolation), _car(input_car) {
	}

	/*
	Performs one timestep of the 11DOF solver
	\param load vector [angle:Z,GC:Y,W1:Y,T1:Y,W2:Y,T2:Y,...]
	\return solution of the following timestep [angle:Z,GC:Y,W1:Y,T1:Y,W2:Y,T2:Y,...]
	*/
	virtual void update_step(T* force, T* solution) = 0;
	/*
	Destructor
	*/
	virtual ~TwoTrackModelParent() {}
};

/**
* \brief class to compute one timestep of the linear 11 dof system in small angle approximation
*/
template <typename T>
class TwoTrackModelBE: public TwoTrackModelParent<T> {
private:
	/**
	* \brief construct A Matrix
	*
	* this has to be called every time step for interpolation
	* A = 1/h^2 * M + 1/h * D + K
	*/
	virtual void constructAMatrix() {
		// A = M/h^2 (also acts as initialization of A)
		mkl<T>::copy(Constants::DOFDOF, M_h2, 1, A, 1);
		// A += 1/h * D => A = 1/h^2 * M + 1/h * D
		mkl<T>::axpy(Constants::DOFDOF, factor_h, D, 1, A, 1);
		// A += K => A = 1/h^2 * M + 1/h * D + K
		mkl<T>::axpy(Constants::DOFDOF, 1, K, 1, A, 1);
	}
	/**
	* \brief construct B
	*
	* this has to be called every time step for interpolation
	* B = 2/(h*h) * M + 1/h * D
	*/
	void constructBMatrix() {
		mkl<T>::scal(Constants::DOFDOF, 0.0, B, 1);
		mkl<T>::axpy(Constants::DOFDOF, 2, M_h2, 1, B, 1);
		mkl<T>::axpy(Constants::DOFDOF, factor_h, D, 1, B, 1);
	}
	/**
	* \brief construct C
	*
	* this has to be called only once
	* C = -1/h^2 M
	*/
	void constructCMatrix() {
		mkl<T>::copy(Constants::DOFDOF,M_h2, 1, C,1);
		mkl<T>::scal(Constants::DOFDOF, -1, C, 1);
	}
protected:
	/**
	* \brief construct Mass matrix
	*
	* this is only done once. Therefore, there is no need to store Mass_vec or I_CG
	*/
	void constructMassMatrix() {
		// get values from MetaDataBase
		T* Mass_vec = _car->massComponents;
		T* I_CG = _car->momentOfInertia;
		// construct M
		M_h2[0] = Mass_vec[0];
		M_h2[Constants::DOF + 1] = I_CG[0];
		M_h2[2 * Constants::DOF + 2] = I_CG[4];
		// M_h2 = M
		mkl<T>::copy(Constants::DOF - 3, Mass_vec + 1, 1, M_h2 + 3 * (Constants::DOF + 1), Constants::DOF + 1); // M_h2 = diagonal matrix
		// M_h2 = M / (h * h)
		mkl<T>::scal(Constants::DOF, 1.0 / (_h * _h), M_h2, Constants::DOF + 1);
	}
	/**
	* \brief construct Mat_springLength
	*
	* this is only done once in the initalisation the matrix is used to calculated th springLenght with blas3
	*/
	/*
	void constructMatSpringLenghts() {
		Mat_springLength[0] = 1;
		Mat_springLength[1] = _car->l_lat[0];
		Mat_springLength[2] = -_car->l_long[0];
		Mat_springLength[3] = -1;
		Mat_springLength[Constants::DOF + 3] = 1;
		Mat_springLength[Constants::DOF + 4] = -1;
		Mat_springLength[2 * Constants::DOF] = 1;
		Mat_springLength[2 * Constants::DOF + 1] = -_car->l_lat[1];
		Mat_springLength[2 * Constants::DOF + 2] = -_car->l_long[1];
		Mat_springLength[2 * Constants::DOF + 5] = -1;
		Mat_springLength[3 * Constants::DOF + 5] = 1;
		Mat_springLength[3 * Constants::DOF + 6] = -1;
		Mat_springLength[4 * Constants::DOF] = 1;
		Mat_springLength[4 * Constants::DOF + 1] = _car->l_lat[2];
		Mat_springLength[4 * Constants::DOF + 2] = _car->l_long[2];
		Mat_springLength[4 * Constants::DOF + 7] = -1;
		Mat_springLength[5 * Constants::DOF + 7] = 1;
		Mat_springLength[5 * Constants::DOF + 8] = -1;
		Mat_springLength[6 * Constants::DOF] = 1;
		Mat_springLength[6 * Constants::DOF + 1] = -_car->l_lat[3];
		Mat_springLength[6 * Constants::DOF + 2] = _car->l_long[3];
		Mat_springLength[6 * Constants::DOF + 9] = -1;
		Mat_springLength[7 * Constants::DOF + 9] = 1;
		Mat_springLength[7 * Constants::DOF + 10] = -1;
	}
	*/
	/**
	* \brief stores the normal spring lengths
	*
	* called in the constructor
	*/
	/*
	void populateSpringLengthsNormal() {
		springLengthsNormal[0] = MetaDataBase::DataBase()->getBodySpringLengthFrontLeft();
		springLengthsNormal[1] = MetaDataBase::DataBase()->getBodySpringLengthFrontLeft();
		springLengthsNormal[2] = MetaDataBase::DataBase()->getBodySpringLengthFrontRight();
		springLengthsNormal[3] = MetaDataBase::DataBase()->getBodySpringLengthFrontRight();
		springLengthsNormal[4] = MetaDataBase::DataBase()->getBodySpringLengthRearLeft();
		springLengthsNormal[5] = MetaDataBase::DataBase()->getBodySpringLengthRearLeft();
		springLengthsNormal[6] = MetaDataBase::DataBase()->getBodySpringLengthRearRight();
		springLengthsNormal[7] = MetaDataBase::DataBase()->getBodySpringInitialLengthRearRight();
	}*/
	/**
	* \brief stores the initial spring lengths
	*
	* called in the constructor
	*/
	/*
	void initSpringLengths() {
		springLengths[0] = MetaDataBase::DataBase()->getBodySpringInitialLengthFrontLeft();
		springLengths[1] = MetaDataBase::DataBase()->getBodySpringInitialLengthFrontLeft();
		springLengths[2] = MetaDataBase::DataBase()->getBodySpringInitialLengthFrontRight();
		springLengths[3] = MetaDataBase::DataBase()->getBodySpringInitialLengthFrontRight();
		springLengths[4] = MetaDataBase::DataBase()->getBodySpringInitialLengthRearLeft();
		springLengths[5] = MetaDataBase::DataBase()->getBodySpringInitialLengthRearLeft();
		springLengths[6] = MetaDataBase::DataBase()->getBodySpringInitialLengthRearRight();
		springLengths[7] = MetaDataBase::DataBase()->getBodySpringInitialLengthRearRight();
	}
	
	*/
	/**
	* \brief initialise the vectors u_n, u_n_m_1, u_n_p_1
	*
	* called in the constructor
	* u_n = ...
	* u_n_m_1 = u_n - _h * velocity
	* u_n_p_1 = u_n
	* tempVector is used for velocity
	*/
	void initPosVectors() {
		// u_n
		/*
		u_n[0] = 0;
		u_n[1] = MetaDataBase::DataBase()->getBodyInitialOrientation()[0];
		u_n[2] = MetaDataBase::DataBase()->getBodyInitialOrientation()[1];
		u_n[3] = springLengthsNormal[0] + u_n[0] + _car->l_lat[0] * u_n[1] - _car->l_long[0] * u_n[2] - springLengths[0];
		u_n[4] = springLengthsNormal[1] + u_n[3] - springLengths[1];
		u_n[5] = springLengthsNormal[2] + u_n[0] - _car->l_lat[0] * u_n[1] - _car->l_long[0] * u_n[2] - springLengths[2];
		u_n[6] = springLengthsNormal[3] + u_n[5] - springLengths[3];
		u_n[7] = springLengthsNormal[4] + u_n[0] + _car->l_lat[0] * u_n[1] + _car->l_long[0] * u_n[2] - springLengths[4];
		u_n[8] = springLengthsNormal[5] + u_n[7] - springLengths[5];
		u_n[9] = springLengthsNormal[6] + u_n[0] - _car->l_lat[0] * u_n[1] + _car->l_long[0] * u_n[2] - springLengths[6];
		u_n[10] = springLengthsNormal[7] + u_n[9] - springLengths[7];
		
		// store velocity in temp
		temp[0] = MetaDataBase::DataBase()->getBodyInitialVelocity()[2];
		temp[1] = MetaDataBase::DataBase()->getBodyInitialAngularVelocity()[0];
		temp[2] = MetaDataBase::DataBase()->getBodyInitialAngularVelocity()[1];
		temp[3] = MetaDataBase::DataBase()->getWheelInitialVelocityFrontLeft()[2];
		temp[4] = MetaDataBase::DataBase()->getTyreInitialVelocityFrontLeft()[2];
		temp[5] = MetaDataBase::DataBase()->getWheelInitialVelocityFrontRight()[2];
		temp[6] = MetaDataBase::DataBase()->getTyreInitialVelocityFrontRight()[2];
		temp[7] = MetaDataBase::DataBase()->getWheelInitialVelocityRearLeft()[2];
		temp[8] = MetaDataBase::DataBase()->getTyreInitialVelocityRearLeft()[2];
		temp[9] = MetaDataBase::DataBase()->getWheelInitialVelocityRearRight()[2];
		temp[10] = MetaDataBase::DataBase()->getTyreInitialVelocityRearRight()[2];
		*/
		// u_n
		mkl<T>::copy(Constants::DOF, u_n_p_1, 1, u_n,1);
		// u_n_m_1 = u_n - _h * velocity
		mkl<T>::copy(Constants::DOF, _car->currentVelocityTwoTrackModel, 1, u_n_m_1, 1);
		mkl<T>::scal(Constants::DOF, -_h, u_n_m_1, 1);
		mkl<T>::axpy(Constants::DOF, 1, u_n, 1, u_n_m_1, 1);
	}

public:
	
	/**
	* \brief Constructor
	*/
	TwoTrackModelBE(Car<T>* input_car, EVAALookup<T>* stiffnessInterpolation, EVAALookup<T>* dampingInterpolation): TwoTrackModelParent<T>(input_car, stiffnessInterpolation, dampingInterpolation) {
		u_n_m_1 = (T*)mkl_malloc(Constants::DOF * sizeof(T), Constants::ALIGNMENT); // velocity
		u_n = (T*)mkl_malloc(Constants::DOF * sizeof(T), Constants::ALIGNMENT); // position
		u_n_p_1 = _car->currentDisplacementTwoTrackModel; // position

		A = (T*)mkl_malloc(Constants::DOFDOF * sizeof(T), Constants::ALIGNMENT);
		B = (T*)mkl_calloc(Constants::DOFDOF, sizeof(T), Constants::ALIGNMENT);
		C = (T*)mkl_calloc(Constants::DOFDOF, sizeof(T), Constants::ALIGNMENT);
		M_h2 = (T*)mkl_calloc(Constants::DOFDOF, sizeof(T), Constants::ALIGNMENT);
		K = (T*)mkl_calloc(Constants::DOFDOF, sizeof(T), Constants::ALIGNMENT);
		D = (T*)mkl_calloc(Constants::DOFDOF, sizeof(T), Constants::ALIGNMENT);

		temp = (T*)mkl_malloc(Constants::DOF * sizeof(T), Constants::ALIGNMENT);		/**< used for the constructStiffnesMatrix and constructDampingMatrix as well as for the BE*/
		Mat_temp = (T*)mkl_calloc(Constants::DOFDOF, sizeof(T), Constants::ALIGNMENT); /**< used for the constructStiffnesMatrix and constructDampingMatrix*/
		Mat_springLength = (T*)mkl_calloc( 2 * Constants::NUM_LEGS * Constants::DOF, sizeof(T), Constants::ALIGNMENT);
		springLengths = (T*)mkl_malloc(Constants::NUM_LEGS * 2 * sizeof(T), Constants::ALIGNMENT);
		springLengthsNormal = (T*)mkl_malloc(Constants::NUM_LEGS * 2 * sizeof(T), Constants::ALIGNMENT);
#ifdef INTERPOLATION
//		ABuffer = (T*)mkl_calloc(Constants::DOFDOF, sizeof(T), Constants::ALIGNMENT);
		J = (T*)mkl_calloc(Constants::DOFDOF, sizeof(T), Constants::ALIGNMENT);
		dKdxx = (T*)mkl_calloc(Constants::DOFDOF, sizeof(T), Constants::ALIGNMENT);
		dDdxx = (T*)mkl_calloc(Constants::DOFDOF, sizeof(T), Constants::ALIGNMENT);
		dkdl = (T*)mkl_malloc(Constants::NUM_LEGS * 2 * sizeof(T), Constants::ALIGNMENT);
		dddl = (T*)mkl_calloc(Constants::NUM_LEGS * 2, sizeof(T), Constants::ALIGNMENT);
		residual = (T*)mkl_malloc(Constants::DOF * sizeof(T), Constants::ALIGNMENT);
#endif
		_h = MetaDataBase::DataBase()->getTimeStepSize();
		factor_h = 1 / _h;
		// populateSpringLengthsNormal();
		// initSpringLengths();
		// constructMatSpringLenghts();

		// if interpolation is used we need alll the lengths to interpolate damping and stiffness if not we have to define them
		// populate kVec
		kVec = _car->kVec;
		// populate dVec
		dVec = _car->dVec;

		initPosVectors();
		constructMassMatrix();
		constructStiffnessMatrix();
		constructDampingMatrix();
		constructAMatrix();
		constructBMatrix();
		constructCMatrix();
	}
	/*
	Destructor
	*/
	virtual ~TwoTrackModelBE() {
		mkl_free(A);
		mkl_free(B);
		mkl_free(C);
		mkl_free(u_n);
		mkl_free(u_n_m_1);
//		mkl_free(u_n_p_1);
		mkl_free(M_h2);
		mkl_free(K);
		mkl_free(D);
		mkl_free(temp);
		mkl_free(Mat_temp);
		mkl_free(springLengths);
		mkl_free(springLengthsNormal);
		mkl_free(Mat_springLength);
#ifdef INTERPOLATION
//		mkl_free(ABuffer);
		mkl_free(J);
		mkl_free(dkdl);
		mkl_free(dddl);
		mkl_free(dKdxx);
		mkl_free(dDdxx);
		mkl_free(residual);
#endif
	}

	/*
	Performs one timestep of the 11DOF solver
	\param load vector [angle:Z,GC:Y,W1:Y,T1:Y,W2:Y,T2:Y,...]
	*/
	virtual void update_step(T* force, T* solution) {
		//mkl<T>::scal(Constants::DOF, -1.0, u_n_m_1, Constants::INCX);
		MathLibrary::Solvers<T, TwoTrackModelBE<T> >::Linear_Backward_Euler(A, B, C, u_n, u_n_m_1, force, u_n_p_1, Constants::DOF);
		constructAMatrix();
#ifdef INTERPOLATION
			MathLibrary::Solvers<T, TwoTrackModelBE<T> >::Newton(this, force, J, residual, &res_norm, u_n_p_1, temp);
#endif
		// solutio = solution[t_n] = u_n_p_1
		mkl<T>::copy(Constants::DOF, u_n_p_1, 1, solution, 1);
		// copy update to car

		MathLibrary::swap_address<T>(u_n, u_n_m_1); // u_n_m_1 points to u_n and u_n points to u_n_m_1
//		MathLibrary::swap_address<T>(u_n_p_1, u_n); // u_n points to u_n_p_1 and u_n_p_1 point to u_n_m_1 now
		// do not swap just copy
		mkl<T>::copy(Constants::DOF, u_n_p_1, 1, u_n, 1);
	}
	/**
	* \brief calculate the residual(Newton Function) + res_norm
	*
	* residual = A*x[n+1] - B * x[n] + M_h2 * x[n-1] - forces
	* res_norm = norm(residual)
	*/
	virtual void calcResidual(
		T* force /**< pointer to forces vector size of DOF*/
	) {
		
		// residual = A*x[n+1]
		mkl<T>::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, Constants::DOF, 1, Constants::DOF, 1, A, Constants::DOF, u_n_p_1, 1, 0, residual, 1);
		// residual -= B*x[n]
		mkl<T>::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, Constants::DOF, 1, Constants::DOF, -1, B, Constants::DOF, u_n, 1, 1, residual, 1);
		// residual += M_h2 * x[n-1]
		mkl<T>::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, Constants::DOF, 1, Constants::DOF, 1, M_h2, Constants::DOF, u_n_m_1, 1, 1, residual, 1);
		// residual -= force
		mkl<T>::axpy(Constants::DOF, -1, force, 1, residual, 1);
		// res = norm(residual)
		res_norm = mkl<T>::nrm2(Constants::DOF, residual, 1);
	}

	/**
	* \brief construct Jacobien for non fixed case
	*
	* this has to be called every newton iteraton
	* J = M_h2 + D / _h + K + dKdx*x[n+1] + 1/h_ * dDdx ( x[n+1] - x[n] )
	*/
	virtual void constructJacobien() {
		// first update the derivative 
		lookupStiffness->getDerivative(_car->currentSpringsLength, dkdl);
	//	lookupDamping->getDerivative(_car->currentSpringsLength, dddl);
		// construct the derivative (tensor) times a pos vector
		constructLookupDerivativeX(dkdl, u_n_p_1, dKdxx);
		// J = A =  M_h2 + D / _h + K
		mkl<T>::copy(Constants::DOFDOF, A, 1, J, 1);
		// J += dKdx * x[n+1]
		mkl<T>::axpy(Constants::DOFDOF, 1.0, dKdxx, 1, J, 1);
		// temp = x[n+1]
		mkl<T>::copy(Constants::DOF, u_n_p_1, 1, temp, 1);
		// temp += -x[n]
		mkl<T>::axpy(Constants::DOF, -1, u_n, 1, temp, 1);
		// calc dDdxx with (x[n+1] - x[n])
//		constructLookupDerivativeX(dddl, temp, dDdxx);
		// J += 1/_h * dDdxx
//		mkl<T>::axpy(Constants::DOFDOF, factor_h, dDdxx, 1, J, 1);
	}
	/**
	* \brief construct Jacobien for fixed to road
	*
	* add the rows according to the tyres of [-(1/h dDdx + dKdx) x[n+1] - 1/h D - K + 1/h dDdx * x[n]]
	* therefore J = M_h2 for the tyre positions
	*/
	virtual void constructFixedJacobien() {
		for (auto i = 0; i < Constants::NUM_LEGS; i++) {
			mkl<T>::copy(Constants::DOF, M_h2 + ( 4 + 2 * i ) * Constants::DOF,1, J + (4 + 2 * i) * Constants::DOF, 1);
		}
	}

	

	/**
	* \brief update all dependent matrices on the position vector
	*
	* only needed in case of lookup Tables
	*/
	virtual void updateSystem() {
		// constructSpringLengths();
		_car->updateLengthsTwoTrackModel();
		constructStiffnessMatrix();
		constructDampingMatrix();
		constructAMatrix();
		constructBMatrix();
	}
	/**
	* \brief construct Stiffness Matrix
	*
	* this has to be called every newton iteraton and the spring lengtvector has to be updated before
	*/
	void constructStiffnessMatrix() {
#ifdef INTERPOLATION
		lookupStiffness->getInterpolation(_car->currentSpringsLength, kVec);
#endif
		temp[0] = kVec[0] + kVec[2] + kVec[4] + kVec[6];
		K[1] = kVec[0] * _car->l_lat[0] - kVec[2] * _car->l_lat[1] + kVec[4] * _car->l_lat[2] - kVec[6] * _car->l_lat[3];
		K[2] = -kVec[0] * _car->l_long[0] - kVec[2] * _car->l_long[1] + kVec[4] * _car->l_long[2] + kVec[6] * _car->l_long[3];
		K[3] = -kVec[0];
		K[4] = 0.0;
		K[5] = -kVec[2];
		K[6] = 0.0;
		K[7] = -kVec[4];
		K[8] = 0.0;
		K[9] = -kVec[6];
		K[10] = 0.0;

		temp[1] = _car->l_lat[0] * _car->l_lat[0] * kVec[0] + _car->l_lat[1] * _car->l_lat[1] * kVec[2] + _car->l_lat[2] * _car->l_lat[2] * kVec[4] + _car->l_lat[3] * _car->l_lat[3] * kVec[6];
		K[Constants::DOF + 2] = -_car->l_long[0] * _car->l_lat[0] * kVec[0] + _car->l_lat[1] * _car->l_long[1] * kVec[2] + _car->l_long[2] * _car->l_lat[2] * kVec[4] - _car->l_long[3] * _car->l_lat[3] * kVec[6];
		K[Constants::DOF + 3] = -_car->l_lat[0] * kVec[0];
		K[Constants::DOF + 4] = 0;
		K[Constants::DOF + 5] = _car->l_lat[1] * kVec[2];
		K[Constants::DOF + 6] = 0;
		K[Constants::DOF + 7] = -_car->l_lat[2] * kVec[4];
		K[Constants::DOF + 8] = 0;
		K[Constants::DOF + 9] = _car->l_lat[3] * kVec[6];
		K[Constants::DOF + 10] = 0;

		temp[2] = _car->l_long[0] * _car->l_long[0] * kVec[0] + _car->l_long[1] * _car->l_long[1] * kVec[2] + _car->l_long[2] * _car->l_long[2] * kVec[4] + _car->l_long[3] * _car->l_long[3] * kVec[6];
		K[2 * Constants::DOF + 3] = _car->l_long[0] * kVec[0];
		K[2 * Constants::DOF + 4] = 0;
		K[2 * Constants::DOF + 5] = _car->l_long[1] * kVec[2];
		K[2 * Constants::DOF + 6] = 0;
		K[2 * Constants::DOF + 7] = -_car->l_long[2] * kVec[4];
		K[2 * Constants::DOF + 8] = 0;
		K[2 * Constants::DOF + 9] = -_car->l_long[3] * kVec[6];
		K[2 * Constants::DOF + 10] = 0;

		temp[3] = kVec[0] + kVec[1];
		K[3 * Constants::DOF + 4] = -kVec[1];
		K[3 * Constants::DOF + 5] = 0;
		K[3 * Constants::DOF + 6] = 0;
		K[3 * Constants::DOF + 7] = 0;
		K[3 * Constants::DOF + 8] = 0;
		K[3 * Constants::DOF + 9] = 0;
		K[3 * Constants::DOF + 10] = 0;
		// all others are zero

		temp[4] = kVec[1];
		K[4 * Constants::DOF + 5] = 0;
		K[4 * Constants::DOF + 6] = 0;
		K[4 * Constants::DOF + 7] = 0;
		K[4 * Constants::DOF + 8] = 0;
		K[4 * Constants::DOF + 9] = 0;
		K[4 * Constants::DOF + 10] = 0;

		temp[5] = kVec[2] + kVec[3];
		K[5 * Constants::DOF + 6] = -kVec[3];
		K[5 * Constants::DOF + 7] = 0;
		K[5 * Constants::DOF + 8] = 0;
		K[5 * Constants::DOF + 9] = 0;
		K[5 * Constants::DOF + 10] = 0;

		temp[6] = kVec[3];
		K[6 * Constants::DOF + 7] = 0;
		K[6 * Constants::DOF + 8] = 0;
		K[6 * Constants::DOF + 9] = 0;
		K[6 * Constants::DOF + 10] = 0;

		temp[7] = kVec[4] + kVec[5];
		K[7 * Constants::DOF + 8] = -kVec[5];
		K[7 * Constants::DOF + 9] = 0;
		K[7 * Constants::DOF + 10] = 0;


		temp[8] = kVec[5];
		K[8 * Constants::DOF + 9] = 0;
		K[8 * Constants::DOF + 10] = 0;


		temp[9] = kVec[6] + kVec[7];
		K[9 * Constants::DOF + 10] = -kVec[7];

		temp[10] = kVec[7];


		// symmetrize K
		//cblas_dcopy(DOF * DOF, K, 1, K_trans, 1);
		mkl<T>::lacpy(LAPACK_ROW_MAJOR, 'U', Constants::DOF, Constants::DOF, K, Constants::DOF, Mat_temp, Constants::DOF);

		mkl<T>::imatcopy('R', 'T', Constants::DOF, Constants::DOF, 1.0, Mat_temp, Constants::DOF, Constants::DOF); // get transpose of matrix

		mkl<T>::lacpy(LAPACK_ROW_MAJOR, 'L', Constants::DOF, Constants::DOF, Mat_temp, Constants::DOF, K, Constants::DOF); // copy lower triangular in the orig matrix
		//cblas_daxpy(DOF * DOF, 1.0, K_trans, 1, K, 1); // K = K + K'

		// add the diagonal to K
		MathLibrary::allocate_to_diagonal(K, temp, Constants::DOF); // K = K + K'+ diag(K)
	}
	/**
	* \brief construct Damping Matrix
	*
	* this has to be called every newton iteraton and the spring lengtvector has to be updated before
	*/
	void constructDampingMatrix() {
#ifdef INTERPOLATION
	//		lookupDamping->getInterpolation(_car->currentSpringsLength, dVec);
#endif
		temp[0] = dVec[0] + dVec[2] + dVec[4] + dVec[6];
		D[1] = dVec[0] * _car->l_lat[0] - dVec[2] * _car->l_lat[1] + dVec[4] * _car->l_lat[2] - dVec[6] * _car->l_lat[3];
		D[2] = -dVec[0] * _car->l_long[0] - dVec[2] * _car->l_long[1] + dVec[4] * _car->l_long[2] + dVec[6] * _car->l_long[3];
		D[3] = -dVec[0];
		D[4] = 0.0;
		D[5] = -dVec[2];
		D[6] = 0.0;
		D[7] = -dVec[4];
		D[8] = 0.0;
		D[9] = -dVec[6];
		D[10] = 0.0;

		temp[1] = _car->l_lat[0] * _car->l_lat[0] * dVec[0] + _car->l_lat[1] * _car->l_lat[1] * dVec[2] + _car->l_lat[2] * _car->l_lat[2] * dVec[4] + _car->l_lat[3] * _car->l_lat[3] * dVec[6];
		D[Constants::DOF + 2] = -_car->l_long[0] * _car->l_lat[0] * dVec[0] + _car->l_lat[1] * _car->l_long[1] * dVec[2] + _car->l_long[2] * _car->l_lat[2] * dVec[4] - _car->l_long[3] * _car->l_lat[3] * dVec[6];
		D[Constants::DOF + 3] = -_car->l_lat[0] * dVec[0];
		D[Constants::DOF + 4] = 0;
		D[Constants::DOF + 5] = _car->l_lat[1] * dVec[2];
		D[Constants::DOF + 6] = 0;
		D[Constants::DOF + 7] = -_car->l_lat[2] * dVec[4];
		D[Constants::DOF + 8] = 0;
		D[Constants::DOF + 9] = _car->l_lat[3] * dVec[6];
		D[Constants::DOF + 10] = 0;

		temp[2] = _car->l_long[0] * _car->l_long[0] * dVec[0] + _car->l_long[1] * _car->l_long[1] * dVec[2] + _car->l_long[2] * _car->l_long[2] * dVec[4] + _car->l_long[3] * _car->l_long[3] * dVec[6];
		D[2 * Constants::DOF + 3] = _car->l_long[0] * dVec[0];
		D[2 * Constants::DOF + 4] = 0;
		D[2 * Constants::DOF + 5] = _car->l_long[1] * dVec[2];
		D[2 * Constants::DOF + 6] = 0;
		D[2 * Constants::DOF + 7] = -_car->l_long[2] * dVec[4];
		D[2 * Constants::DOF + 8] = 0;
		D[2 * Constants::DOF + 9] = -_car->l_long[3] * dVec[6];
		D[2 * Constants::DOF + 10] = 0;

		temp[3] = dVec[0] + dVec[1];
		D[3 * Constants::DOF + 4] = -dVec[1];
		D[3 * Constants::DOF + 5] = 0;
		D[3 * Constants::DOF + 6] = 0;
		D[3 * Constants::DOF + 7] = 0;
		D[3 * Constants::DOF + 8] = 0;
		D[3 * Constants::DOF + 9] = 0;
		D[3 * Constants::DOF + 10] = 0;
		// all others are zero

		temp[4] = dVec[1];
		D[4 * Constants::DOF + 5] = 0;
		D[4 * Constants::DOF + 6] = 0;
		D[4 * Constants::DOF + 7] = 0;
		D[4 * Constants::DOF + 8] = 0;
		D[4 * Constants::DOF + 9] = 0;
		D[4 * Constants::DOF + 10] = 0;

		temp[5] = dVec[2] + dVec[3];
		D[5 * Constants::DOF + 6] = -dVec[3];
		D[5 * Constants::DOF + 7] = 0;
		D[5 * Constants::DOF + 8] = 0;
		D[5 * Constants::DOF + 9] = 0;
		D[5 * Constants::DOF + 10] = 0;

		temp[6] = dVec[3];
		D[6 * Constants::DOF + 7] = 0;
		D[6 * Constants::DOF + 8] = 0;
		D[6 * Constants::DOF + 9] = 0;
		D[6 * Constants::DOF + 10] = 0;

		temp[7] = dVec[4] + dVec[5];
		D[7 * Constants::DOF + 8] = -dVec[5];
		D[7 * Constants::DOF + 9] = 0;
		D[7 * Constants::DOF + 10] = 0;


		temp[8] = dVec[5];
		D[8 * Constants::DOF + 9] = 0;
		D[8 * Constants::DOF + 10] = 0;


		temp[9] = dVec[6] + dVec[7];
		D[9 * Constants::DOF + 10] = -dVec[7];

		temp[10] = dVec[7];


		// symmetrize D
		//cblas_dcopy(DOF * DOF, D, 1, D_trans, 1);
		mkl<T>::lacpy(LAPACK_ROW_MAJOR, 'U', Constants::DOF, Constants::DOF, D, Constants::DOF, Mat_temp, Constants::DOF);

		mkl<T>::imatcopy('R', 'T', Constants::DOF, Constants::DOF, 1.0, Mat_temp, Constants::DOF, Constants::DOF); // get transpose of matrix

		mkl<T>::lacpy(LAPACK_ROW_MAJOR, 'L', Constants::DOF, Constants::DOF, Mat_temp, Constants::DOF, D, Constants::DOF); // copy lower triangular in the orig matrix
		//cblas_daxpy(DOF * DOF, 1.0, D_trans, 1, D, 1); // D = D + D'

		// add the diagonal to K
		MathLibrary::allocate_to_diagonal(D, temp, Constants::DOF); // D = D + D'+ diag(D)
	}

	/**
	* \brief construct springLengths
	*
	* this has to be called every newton iteraton
	* springLengths = Mat_SpringLenghts * x + springLengthsNormal
	*/
	/*
	void constructSpringLengths() {
		// springLengths = springLengthsNormal
		mkl<T>::copy(Constants::NUM_LEGS * 2, springLengthsNormal, 1, springLengths, 1);
		// springLengths += Mat_SpringLenghts * x
		mkl<T>::gemv(CblasRowMajor, CblasNoTrans, Constants::NUM_LEGS * 2, Constants::DOF, 1.0, Mat_springLength, Constants::DOF, u_n_p_1, 1.0, 1.0, springLengths, 1.0);
	}
	*/
	/**
	* \brief construct The derivative of the stiffness lookupTable times a position vector
	*
	* this has to be called every newton iteraton
	*/
	void constructLookupDerivativeX(
		T* der /**< [in] pointer to vector with the derivative of the lookup Table of length 8 */,
		T* x /**< [in] pointer to position vector of length DOF */,
		T* dMdxx /**< [out] pointer to matrix in which de derivative of the lookup times a position vec is stored of size DOF * DOF */
	) {
		temp[0] = x[0] * (der[0] + der[2] + der[4] + der[6]) - der[2] * x[5] - der[4] * x[7] - der[6] * x[9] - der[0] * x[3] + x[1] * (der[0] * _car->l_lat[0] - der[2] * _car->l_lat[1] + der[4] * _car->l_lat[2] - der[6] * _car->l_lat[3]) - x[2] * (der[0] * _car->l_long[0] + der[2] * _car->l_long[1] - der[4] * _car->l_long[2] - der[6] * _car->l_long[3]);
		dMdxx[1] = der[2] * _car->l_lat[1] * _car->l_lat[1] * x[1] + der[4] * _car->l_lat[2] * _car->l_lat[2] * x[1] + der[6] * _car->l_lat[3] * _car->l_lat[3] * x[1] + der[0] * _car->l_lat[0] * (x[0] - x[3] + _car->l_lat[0] * x[1]) - der[2] * _car->l_lat[1] * x[0] + der[2] * _car->l_lat[1] * x[5] + der[4] * _car->l_lat[2] * x[0] - der[4] * _car->l_lat[2] * x[7] - der[6] * _car->l_lat[3] * x[0] + der[6] * _car->l_lat[3] * x[9] - der[0] * _car->l_lat[0] * _car->l_long[0] * x[2] + der[2] * _car->l_lat[1] * _car->l_long[1] * x[2] + der[4] * _car->l_lat[2] * _car->l_long[2] * x[2] - der[6] * _car->l_lat[3] * _car->l_long[3] * x[2];
		dMdxx[2] = der[0] * _car->l_long[0] * _car->l_long[0] * x[2] + der[4] * _car->l_long[2] * _car->l_long[2] * x[2] + der[6] * _car->l_long[3] * _car->l_long[3] * x[2] - der[0] * _car->l_long[0] * (x[0] - x[3] + _car->l_lat[0] * x[1]) + der[2] * _car->l_long[1] * (x[5] - x[0] + _car->l_lat[1] * x[1] + _car->l_long[1] * x[2]) + der[4] * _car->l_long[2] * x[0] - der[4] * _car->l_long[2] * x[7] + der[6] * _car->l_long[3] * x[0] - der[6] * _car->l_long[3] * x[9] + der[4] * _car->l_lat[2] * _car->l_long[2] * x[1] - der[6] * _car->l_lat[3] * _car->l_long[3] * x[1];
		dMdxx[3] = -der[0] * (x[0] - x[3] + _car->l_lat[0] * x[1] - _car->l_long[0] * x[2]);
		dMdxx[4] = 0.0;
		dMdxx[5] = der[2] * (x[5] - x[0] + _car->l_lat[1] * x[1] + _car->l_long[1] * x[2]);
		dMdxx[6] = 0.0;
		dMdxx[7] = -der[4] * (x[0] - x[7] + _car->l_lat[2] * x[1] + _car->l_long[2] * x[2]);
		dMdxx[8] = 0.0;
		dMdxx[9] = -der[6] * (x[0] - x[9] - _car->l_lat[3] * x[1] + _car->l_long[3] * x[2]);
		dMdxx[10] = 0.0;

		temp[1] = der[2] * _car->l_lat[1] * _car->l_lat[1] * x[0] - der[2] * _car->l_lat[1] * _car->l_lat[1] * _car->l_lat[1] * x[1] - der[2] * _car->l_lat[1] * _car->l_lat[1] * x[5] + der[4] * _car->l_lat[2] * _car->l_lat[2] * x[0] + der[4] * _car->l_lat[2] * _car->l_lat[2] * _car->l_lat[2] * x[1] - der[4] * _car->l_lat[2] * _car->l_lat[2] * x[7] + der[6] * _car->l_lat[3] * _car->l_lat[3] * x[0] - der[6] * _car->l_lat[3] * _car->l_lat[3] * _car->l_lat[3] * x[1] - der[6] * _car->l_lat[3] * _car->l_lat[3] * x[9] + der[0] * _car->l_lat[0] * _car->l_lat[0] * (x[0] - x[3] + _car->l_lat[0] * x[1]) - der[0] * _car->l_lat[0] * _car->l_lat[0] * _car->l_long[0] * x[2] - der[2] * _car->l_lat[1] * _car->l_lat[1] * _car->l_long[1] * x[2] + der[4] * _car->l_lat[2] * _car->l_lat[2] * _car->l_long[2] * x[2] + der[6] * _car->l_lat[3] * _car->l_lat[3] * _car->l_long[3] * x[2];
		dMdxx[Constants::DOF + 2] = der[2] * _car->l_lat[1] * _car->l_long[1] * x[0] - der[2] * _car->l_lat[1] * _car->l_long[1] * x[5] + der[4] * _car->l_lat[2] * _car->l_long[2] * x[0] - der[4] * _car->l_lat[2] * _car->l_long[2] * x[7] - der[6] * _car->l_lat[3] * _car->l_long[3] * x[0] + der[6] * _car->l_lat[3] * _car->l_long[3] * x[9] + der[0] * _car->l_lat[0] * _car->l_long[0] * _car->l_long[0] * x[2] - der[2] * _car->l_lat[1] * _car->l_lat[1] * _car->l_long[1] * x[1] - der[2] * _car->l_lat[1] * _car->l_long[1] * _car->l_long[1] * x[2] + der[4] * _car->l_lat[2] * _car->l_lat[2] * _car->l_long[2] * x[1] + der[4] * _car->l_lat[2] * _car->l_long[2] * _car->l_long[2] * x[2] + der[6] * _car->l_lat[3] * _car->l_lat[3] * _car->l_long[3] * x[1] - der[6] * _car->l_lat[3] * _car->l_long[3] * _car->l_long[3] * x[2] - der[0] * _car->l_lat[0] * _car->l_long[0] * (x[0] - x[3] + _car->l_lat[0] * x[1]);
		dMdxx[Constants::DOF + 3] = -der[0] * _car->l_lat[0] * (x[0] - x[3] + _car->l_lat[0] * x[1] - _car->l_long[0] * x[2]);
		dMdxx[Constants::DOF + 4] = 0;
		dMdxx[Constants::DOF + 5] = -der[2] * _car->l_lat[1] * (x[5] - x[0] + _car->l_lat[1] * x[1] + _car->l_long[1] * x[2]);
		dMdxx[Constants::DOF + 6] = 0;
		dMdxx[Constants::DOF + 7] = -der[4] * _car->l_lat[2] * (x[0] - x[7] + _car->l_lat[2] * x[1] + _car->l_long[2] * x[2]);
		dMdxx[Constants::DOF + 8] = 0;
		dMdxx[Constants::DOF + 9] = der[6] * _car->l_lat[3] * (x[0] - x[9] - _car->l_lat[3] * x[1] + _car->l_long[3] * x[2]);
		dMdxx[Constants::DOF + 10] = 0;

		temp[2] = der[2] * _car->l_long[1] * _car->l_long[1] * x[0] - der[0] * _car->l_long[0] * _car->l_long[0] * _car->l_long[0] * x[2] - der[2] * _car->l_long[1] * _car->l_long[1] * _car->l_long[1] * x[2] - der[2] * _car->l_long[1] * _car->l_long[1] * x[5] + der[4] * _car->l_long[2] * _car->l_long[2] * x[0] + der[4] * _car->l_long[2] * _car->l_long[2] * _car->l_long[2] * x[2] - der[4] * _car->l_long[2] * _car->l_long[2] * x[7] + der[6] * _car->l_long[3] * _car->l_long[3] * x[0] + der[6] * _car->l_long[3] * _car->l_long[3] * _car->l_long[3] * x[2] - der[6] * _car->l_long[3] * _car->l_long[3] * x[9] + der[0] * _car->l_long[0] * _car->l_long[0] * (x[0] - x[3] + _car->l_lat[0] * x[1]) - der[2] * _car->l_lat[1] * _car->l_long[1] * _car->l_long[1] * x[1] + der[4] * _car->l_lat[2] * _car->l_long[2] * _car->l_long[2] * x[1] - der[6] * _car->l_lat[3] * _car->l_long[3] * _car->l_long[3] * x[1];
		dMdxx[2 * Constants::DOF + 3] = der[0] * _car->l_long[0] * (x[0] - x[3] + _car->l_lat[0] * x[1] - _car->l_long[0] * x[2]);
		dMdxx[2 * Constants::DOF + 4] = 0;
		dMdxx[2 * Constants::DOF + 5] = -der[2] * _car->l_long[1] * (x[5] - x[0] + _car->l_lat[1] * x[1] + _car->l_long[1] * x[2]);
		dMdxx[2 * Constants::DOF + 6] = 0;
		dMdxx[2 * Constants::DOF + 7] = -der[2] * _car->l_long[1] * (x[5] - x[0] + _car->l_lat[1] * x[1] + _car->l_long[1] * x[2]);
		dMdxx[2 * Constants::DOF + 8] = 0;
		dMdxx[2 * Constants::DOF + 9] = -der[2] * _car->l_long[1] * (x[5] - x[0] + _car->l_lat[1] * x[1] + _car->l_long[1] * x[2]);
		dMdxx[2 * Constants::DOF + 10] = 0;

		temp[3] = der[1] * (x[3] - x[4]) + der[0] * (x[0] - x[3] + _car->l_lat[0] * x[1]) - der[0] * _car->l_long[0] * x[2];
		dMdxx[3 * Constants::DOF + 4] = -der[1] * (x[3] - x[4]);
		dMdxx[3 * Constants::DOF + 5] = 0;
		dMdxx[3 * Constants::DOF + 6] = 0;
		dMdxx[3 * Constants::DOF + 7] = 0;
		dMdxx[3 * Constants::DOF + 8] = 0;
		dMdxx[3 * Constants::DOF + 9] = 0;
		dMdxx[3 * Constants::DOF + 10] = 0;
		// all others are zero

		temp[4] = der[1] * (x[3] - x[4]);
		dMdxx[4 * Constants::DOF + 5] = 0;
		dMdxx[4 * Constants::DOF + 6] = 0;
		dMdxx[4 * Constants::DOF + 7] = 0;
		dMdxx[4 * Constants::DOF + 8] = 0;
		dMdxx[4 * Constants::DOF + 9] = 0;
		dMdxx[4 * Constants::DOF + 10] = 0;

		temp[5] = der[3] * (x[5] - x[6]) - der[2] * (x[5] - x[0] + _car->l_lat[1] * x[1] + _car->l_long[1] * x[2]);
		dMdxx[5 * Constants::DOF + 6] = -der[3] * (x[5] - x[6]);
		dMdxx[5 * Constants::DOF + 7] = 0;
		dMdxx[5 * Constants::DOF + 8] = 0;
		dMdxx[5 * Constants::DOF + 9] = 0;
		dMdxx[5 * Constants::DOF + 10] = 0;

		temp[6] = der[3] * (x[5] - x[6]);
		dMdxx[6 * Constants::DOF + 7] = 0;
		dMdxx[6 * Constants::DOF + 8] = 0;
		dMdxx[6 * Constants::DOF + 9] = 0;
		dMdxx[6 * Constants::DOF + 10] = 0;

		temp[7] = der[5] * (x[7] - x[8]) + der[4] * (x[0] - x[7] + _car->l_lat[2] * x[1] + _car->l_long[2] * x[2]);
		dMdxx[7 * Constants::DOF + 8] = -der[5] * (x[7] - x[8]);
		dMdxx[7 * Constants::DOF + 9] = 0;
		dMdxx[7 * Constants::DOF + 10] = 0;


		temp[8] = der[5] * (x[7] - x[8]);
		dMdxx[8 * Constants::DOF + 9] = 0;
		dMdxx[8 * Constants::DOF + 10] = 0;


		temp[9] = der[7] * (x[9] - x[10]) + der[6] * (x[0] - x[9] - _car->l_lat[3] * x[1] + _car->l_long[3] * x[2]);
		dMdxx[9 * Constants::DOF + 10] = -der[7] * (x[9] - x[10]);

		temp[10] = der[7] * (x[9] - x[10]);


		// symmetrize dMdxx
		//cblas_dcopy(DOF * DOF, dMdxx, 1, dMdxx_trans, 1);
		mkl<T>::lacpy(LAPACK_ROW_MAJOR, 'U', Constants::DOF, Constants::DOF, dMdxx, Constants::DOF, Mat_temp, Constants::DOF);

		mkl<T>::imatcopy('R', 'T', Constants::DOF, Constants::DOF, 1.0, Mat_temp, Constants::DOF, Constants::DOF); // get transpose of matrix

		mkl<T>::lacpy(LAPACK_ROW_MAJOR, 'L', Constants::DOF, Constants::DOF, Mat_temp, Constants::DOF, dMdxx, Constants::DOF); // copy lower triangular in the orig matrix
		//cblas_daxpy(DOF * DOF, 1.0, dMdxx_trans, 1, dMdxx, 1); // dMdxx = dMdxx + dMdxx'

		// add the diagonal to dM
		MathLibrary::allocate_to_diagonal(dMdxx, temp, Constants::DOF); // dMdxx = dMdxx + dMdxx'+ diag(dMdxx)
	}
};

template <typename T>
class TwoTrackModelBDF2 : public TwoTrackModelBE<T> {
private:
	T*Dmat, *E;
	T* bVec; /**< this vector is used to store the whole constants part on the rhs apart from the force */
	T* u_n_m_2, *u_n_m_3;
	size_t time_step_count = 0;
	void(TwoTrackModelBDF2<T>::*_active_executor)(T*, T*);
	
	/**
	* \brief construct A
	*
	* A = 9/4h^2 * M + 3/2h * D + K
	*/
	virtual void constructAMatrix() {
		// A = M/h^2
		mkl<T>::copy(Constants::DOFDOF, M_h2, 1, A, 1);
		// A *= 9/4
		mkl<T>::scal(Constants::DOFDOF, 2.25, A, 1);
		// A += 3/2h * D
		mkl<T>::axpy(Constants::DOFDOF, 1.5 * factor_h, D, 1, A, 1);
		// A += K
		mkl<T>::axpy(Constants::DOFDOF, 1, K, 1, A, 1);
	}
	/**
	* \brief construct bVec
	*
	* this has to be called every time step for interpolation and is only used in case of interpolation
	* bVec = 1/h^2 * M * (6 * x[n] - 11/2 * x[n-1] + 2 * x[n-2] - 1/4 * x[n-3]) + 1/h * D * (2 * x[n] - 1/2 * x[n-1])
	*/
	void constructbVec() {
		// temp = x[n]
		mkl<T>::copy(Constants::DOF, u_n, 1, temp, 1);
		// temp *= 6
		mkl<T>::scal(Constants::DOF, 6, temp, 1);
		// temp += - 11/2 * x[n-1]
		mkl<T>::axpy(Constants::DOF, -5.5, u_n_m_1, 1, temp, 1);
		// temp += 2 * x[n-2]
		mkl<T>::axpy(Constants::DOF, 2, u_n_m_2, 1, temp, 1);
		// temp += -1/4 * x[n-3]
		mkl<T>::axpy(Constants::DOF, -0.25, u_n_m_3, 1, temp, 1);
		// bVec = M_h2 * temp
		mkl<T>::gemv(CblasRowMajor, CblasNoTrans, Constants::DOF, Constants::DOF, 1, M_h2, Constants::DOF, temp, 1, 0, bVec, 1);
		// temp = x[n]
		mkl<T>::copy(Constants::DOF, u_n, 1, temp, 1);
		// temp *= 2
		mkl<T>::scal(Constants::DOF, 2, temp, 1);
		// temp += - 1/2 * x[n-1]
		mkl<T>::axpy(Constants::DOF, -0.5, u_n_m_1, 1, temp, 1);
		// bVec += 1/h * D * temp
		mkl<T>::gemv(CblasRowMajor, CblasNoTrans, Constants::DOF, Constants::DOF, factor_h, D, Constants::DOF, temp, 1, 1, bVec, 1);
	}

public:
	
	/*
	Constructor
	*/
	TwoTrackModelBDF2(Car<T>* input_car, EVAALookup<T>* stiffnessInterpolation, EVAALookup<T>* dampingInterpolation) : TwoTrackModelBE<T>(input_car, stiffnessInterpolation, dampingInterpolation) {
		Dmat = (T*)mkl_calloc(Constants::DOFDOF, sizeof(T), Constants::ALIGNMENT);
		E = (T*)mkl_calloc(Constants::DOFDOF, sizeof(T), Constants::ALIGNMENT);
		u_n_m_2 = (T*)mkl_malloc(Constants::DOF * sizeof(T), Constants::ALIGNMENT);
		u_n_m_3 = (T*)mkl_malloc(Constants::DOF * sizeof(T), Constants::ALIGNMENT);
		bVec = (T*)mkl_malloc(Constants::DOF * sizeof(T), Constants::ALIGNMENT);
		_active_executor = &Linear11dofBDF2<T>::first_two_steps;

#ifndef INTERPOLATION
		// C = (-11/2)*(1/(h*h))*M + (-1/(2*h))*D
		mkl<T>::axpy(Constants::DOFDOF, (-11.0 / 2.0), M_h2, 1, C, 1);
		mkl<T>::axpy(Constants::DOFDOF, (-1.0 / 2.0) * factor_h, D, 1, C, 1);

		// D = (2/(h*h))*M
		mkl<T>::axpy(Constants::DOFDOF, 2.0, M_h2, 1, Dmat, 1);

		// E = (-1/4)*(1/(h*h))*M
		mkl<T>::axpy(Constants::DOFDOF, (-1.0 / 4.0), M_h2, 1, E, 1);

		mkl<T>::copy(Constants::DOF, _car->u_prev_linear, 1, u_n, 1);
		mkl<T>::copy(Constants::DOF, _car->velocity_current_linear, 1, u_n_m_1, 1);

		mkl<T>::scal(Constants::DOF, -_h, u_n_m_1, 1);
		mkl<T>::axpy(Constants::DOF, 1, u_n, 1, u_n_m_1, 1);
#endif // INTERPOLATION
	}

	void first_two_steps(T* force, T* solution) {
		
		if (time_step_count == 0) {
			mkl<T>::copy(Constants::DOF, u_n_m_1, 1, u_n_m_2, 1);
			TwoTrackModelBE<T>::update_step(force, solution);
		}
		else if (time_step_count == 1) {
			mkl<T>::copy(Constants::DOF, u_n_m_2, 1, u_n_m_3, 1);
			mkl<T>::copy(Constants::DOF, u_n_m_1, 1, u_n_m_2, 1);
			TwoTrackModelBE<T>::update_step(force, solution);
		}
		else {
			mkl<T>::copy(Constants::DOF, u_n_m_2, 1, u_n_m_3, 1);
			mkl<T>::copy(Constants::DOF, u_n_m_1, 1, u_n_m_2, 1);
			// construct A
			constructAMatrix();
#ifdef INTERPOLATION
			constructbVec();
#else
			// B = (6/(h*h))*M + (2/h)*D
			mkl<T>::scal(Constants::DOFDOF, 0.0, B, 1);
			mkl<T>::axpy(Constants::DOFDOF, 6, M_h2, 1, B, 1);
			mkl<T>::axpy(Constants::DOFDOF, 2 * factor_h, D, 1, B, 1);
#endif // INTERPOLATION
			update_step_bdf2(force, solution);
			_active_executor = &Linear11dofBDF2<T>::update_step_bdf2;
		}
		
	}
	virtual void update_step(T* force, T* solution) {
		(this->*_active_executor)(force, solution);
		time_step_count += 1;
	}

	/*
	Performs one timestep of the 11DOF solver
	\param load vector [angle:Z,GC:Y,W1:Y,T1:Y,W2:Y,T2:Y,...]
	\return solution of the following timestep [angle:Z,GC:Y,W1:Y,T1:Y,W2:Y,T2:Y,...]
	*/
	void update_step_bdf2(T* force, T* solution) {
		//cblas_dscal(DOF, 0.0, force, 1);
		MathLibrary::Solvers<T, TwoTrackModelBDF2>::Linear_BDF2(A, B, C, Dmat, E, u_n, u_n_m_1, u_n_m_2, u_n_m_3, force, u_n_p_1, Constants::DOF);
#ifdef INTERPOLATION
			MathLibrary::Solvers<T, TwoTrackModelBDF2<T>>::Newton(this, force, J, residual, &res_norm, u_n_p_1, temp);
#endif
		/*compute_normal_force(K, u_n_p_1, f_n_p_1, tyre_index_set, DOF, num_tyre);
		apply_normal_force(f_n_p_1, u_n_p_1, tyre_index_set, num_tyre);*/
		mkl<T>::copy(Constants::DOF, u_n_p_1, 1, solution, 1);
		MathLibrary::swap_address<T>(u_n_m_2, u_n_m_3); // u_n_m_2 points to u_n_m_3 and u_n_m_3 points to u_n_m_2
		MathLibrary::swap_address<T>(u_n_m_1, u_n_m_2); // u_n_m_2 points to u_n_m_1 and u_n_m_1 points to u_n_m_3
		MathLibrary::swap_address<T>(u_n, u_n_m_1); // u_n_m_1 points to u_n and u_n points to u_n_m_3
		mkl<T>::copy(Constants::DOF, u_n_p_1, 1, u_n, 1);

	}
	/**
	* \brief calculate the residual(Newton Function) + res_norm
	*
	* residual = A*x[n+1] - bVec - forces
	* res_norm = norm(residual)
	*/
	virtual void calcResidual(
		T* force /**< pointer to forces vector size of DOF*/
	) {
		// residual = A*x[n+1]
		mkl<T>::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, Constants::DOF, 1, Constants::DOF, 1, A, Constants::DOF, u_n_p_1, 1, 0, residual, 1);
		// residual -= bVec
		mkl<T>::axpy(Constants::DOF, -1, bVec, 1, residual, 1);
		// residual -= force
		mkl<T>::axpy(Constants::DOF, -1, force, 1, residual, 1);
		// res = norm(residual)
		res_norm = mkl<T>::nrm2(Constants::DOF, residual, 1);
	}
	/**
	* \brief update all dependent matrices on the position vector
	*
	* only needed in case of lookup Tables
	*/
	virtual void updateSystem() {
		// constructSpringLengths();
		_car->updateLengthsTwoTrackModel();
		constructStiffnessMatrix();
		constructDampingMatrix();
		constructAMatrix();
		constructbVec();
	}
	/**
	* \brief construct Jacobien
	*
	* this has to be called every newton iteraton
	* J = 9/4 * M_h2 + 3/2h * D + K + dKdx*x[n+1] + 1/h * dDdx * (3/2 * x[n+1] - 2 * x[n] + 1/2 * x[n-1])
	*/
	virtual void constructJacobien() {
		// first update the derivative 
		lookupStiffness->getDerivative(_car->currentSpringsLength, dkdl);
//		lookupDamping->getDerivative(_car->currentSpringsLength, dddl);
		// construct the derivative (tensor) times a pos vector
		constructLookupDerivativeX(dddl, u_n_p_1, dDdxx);
		// J = A
		mkl<T>::copy(Constants::DOFDOF, A, 1, J, 1);
		// J += dKdx * x[n+1]
		mkl<T>::axpy(Constants::DOFDOF, 1.0, dKdxx, 1, J, 1);
		// temp = x[n+1]
		mkl<T>::copy(Constants::DOF, u_n_p_1, 1, temp, 1);
		// temp *= 3/2
		mkl<T>::scal(Constants::DOF, 1.5, temp, 1);
		// temp += -2 * x[n]
		mkl<T>::axpy(Constants::DOF, -2, u_n, 1, temp, 1);
		// temp += 1/2 * x[n-1]
		mkl<T>::axpy(Constants::DOF, 0.5, u_n_m_1, 1, temp, 1);
		// calc dDdxx with (3/2 * x[n+1] - 2 * x[n] + 1/2 * x[n-1])
		constructLookupDerivativeX(dddl, temp, dDdxx);
		// J += 1/_h * dDdxx
		mkl<T>::axpy(Constants::DOFDOF, factor_h, dDdxx, 1, J, 1);
	}
	/*
	Destructor
	*/
	virtual ~TwoTrackModelBDF2() {
		mkl_free(Dmat);
		mkl_free(E);
		mkl_free(u_n_m_2);
		mkl_free(u_n_m_3);
		mkl_free(bVec);
	}
};

/*
For testing purposes
*/
template <typename T>
class TwoTrackModelFull : public TwoTrackModelBE<T> {
public:
	TwoTrackModelFull(Car<T>* input_car, EVAALookup<T>* stiffnessInterpolation, EVAALookup<T>* dampingInterpolation): TwoTrackModelBE<T>(input_car, stiffnessInterpolation, dampingInterpolation) {

		tend_ = MetaDataBase::DataBase()->getNumberOfTimeIterations() * _h;
		int sol_size = (floor(tend_ / _h) + 1);
		f_n_p_1 = (T*)mkl_malloc(Constants::DOF * sizeof(T), Constants::ALIGNMENT);
		u_sol = (T*)mkl_calloc((sol_size + 1) * (Constants::DOF), sizeof(T), Constants::ALIGNMENT);


		f_n_p_1[0] = MetaDataBase::DataBase()->getBodyExternalForce()[2];
		f_n_p_1[3] = MetaDataBase::DataBase()->getWheelExternalForceFrontLeft()[2];
		f_n_p_1[4] = MetaDataBase::DataBase()->getTyreExternalForceFrontLeft()[2];
		f_n_p_1[5] = MetaDataBase::DataBase()->getWheelExternalForceFrontRight()[2];
		f_n_p_1[6] = MetaDataBase::DataBase()->getTyreExternalForceFrontRight()[2];
		f_n_p_1[7] = MetaDataBase::DataBase()->getWheelExternalForceFrontLeft()[2];
		f_n_p_1[8] = MetaDataBase::DataBase()->getTyreExternalForceFrontLeft()[2];
		f_n_p_1[9] = MetaDataBase::DataBase()->getWheelExternalForceFrontRight()[2];
		f_n_p_1[10] = MetaDataBase::DataBase()->getTyreExternalForceFrontRight()[2];


	}
	void apply_boundary_condition(int s) {
		condition_type = s;
		if (s == NONFIXED) {
			// don't do anything, proceed as usual (solve full 11DOF system)
		}
		else {
			throw "Incorrect boundary condition";
		}
	}
	void solve(T* sol_vect) {
		int iter = 1;
		T t = _h;
		double eps = _h / 100;
		T* solution_vect = u_sol;
		//cblas_dcopy(DOF, _car->u_prev_linear, 1, solution_vect, 1);
		while (std::abs(t - (tend_ + _h)) > eps) {
			//solution_vect = u_sol + iter * (DOF);
//			std::cout << "Running time step = " << iter << std::endl;
			update_step(f_n_p_1, sol_vect);
/*			if (iter > 3900) {
				MathLibrary::write_vector(sol_vect, 11);
			}
			*/
			iter++;
			t += _h;
		}
		MathLibrary::write_vector(_car->currentSpringsLength, 8);
		//cblas_dcopy(DOF, u_sol + (iter - 1)*(DOF), 1, sol_vect, 1);
	}

	void print_final_results(T* sln) {
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
		mkl_free(u_sol);
		mkl_free(f_n_p_1);
	}

private:
	T tend_;
	T* u_sol, * f_n_p_1;
	int condition_type;
};