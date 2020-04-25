/***********************************************************************************************//**
* \file 11DOF.h
* This file holds the function declaration and definitions of the linear11dof and the linear11dof_full class.
* \date 04/14/2020
**************************************************************************************************/
#pragma once
#define _CRTDBG_MAP_ALLOC
#include <cstdlib>
#include <crtdbg.h>

#ifdef _DEBUG
#define DBG_NEW new ( _NORMAL_BLOCK , __FILE__ , __LINE__ )
// Replace _NORMAL_BLOCK with _CLIENT_BLOCK if you want the
// allocations to be of _CLIENT_BLOCK type
#else
#define DBG_NEW new
#endif
#include <mkl.h>
#include "Car.h"
#include "MathLibrary.h"
#include "Constants.h"
#include "MetaDataBase.h"
#include "BLAS.h"

/**
* \brief class to compute one timestep of the linear 11 dof system in small angle approximation
*/
template <typename T>
class Linear11dofBE {
private:
	/**
	* \brief construct A Matrix
	*
	* this has to be called every time step for interpolation
	* A = 1/h^2 * M + 1/h * D + K
	*/
	void constructAMatrix() {
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
		mkl<T>::axpy(Constants::DOF, 2, M_h2, Constants::DOF + 1, B, Constants::DOF + 1);
		mkl<T>::axpy(Constants::DOFDOF, factor_h, D, 1, B, 1);
	}
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
	T* A, * B;                                                  /**< pointer to matrices used for the backward euler */
	T* M_h2, * K, * D, * Mat_temp, * Mat_springLength;

	T* kVec, * dVec, * temp;
	T* springLengths, * springLengthsNormal;
	size_t* tyre_index_set;

	T* J, * dKdxx, * dDdxx;
	T* residual;
	T res_norm;
	T* dkdl, * dddl;

	// interpolation 
	EVAALookup<T>* lookupStiffness;
	EVAALookup<T>* lookupDamping;


	void compute_normal_force(T* K, T* u, T* force, size_t* index, size_t dim, size_t n) {
#pragma loop( ivdep )
		for (int i = 0; i < n; ++i) {
			force[index[i]] = -K[index[i] * dim + index[i]] * u[index[i]];
		}
	}
	/**
	* \brief construct Mass matrix
	*
	* this is only done once. Therefore, there is no need to store Mass_vec or I_CG
	*/
	void constructMassMatrix() {
		// get values from MetaDataBase
		T Mass_vec[9];
		Mass_vec[0] = MetaDataBase::DataBase()->getBodyMass();
		Mass_vec[1] = MetaDataBase::DataBase()->getWheelMassFrontLeft();
		Mass_vec[2] = MetaDataBase::DataBase()->getTyreMassFrontLeft();
		Mass_vec[3] = MetaDataBase::DataBase()->getWheelMassFrontRight();
		Mass_vec[4] = MetaDataBase::DataBase()->getTyreMassFrontRight();
		Mass_vec[5] = MetaDataBase::DataBase()->getWheelMassRearLeft();
		Mass_vec[6] = MetaDataBase::DataBase()->getTyreMassRearLeft();
		Mass_vec[7] = MetaDataBase::DataBase()->getWheelMassRearRight();
		Mass_vec[8] = MetaDataBase::DataBase()->getTyreMassRearRight();
		T* I_CG = MetaDataBase::DataBase()->getMomentOfInertiaVector();
		// construct M
		M_h2[0] = Mass_vec[0];
		M_h2[Constants::DOF + 1] = I_CG[0];
		M_h2[2 * Constants::DOF + 2] = I_CG[4];
		mkl<T>::copy(Constants::DOF - 3, Mass_vec + 1, 1, M_h2 + 3 * (Constants::DOF + 1), Constants::DOF + 1); // M_h2 = diagonal matrix
		mkl<T>::scal(Constants::DOF, 1 / (_h * _h), M_h2, Constants::DOF + 1);
	}
	/**
	* \brief construct Mat_springLength
	*
	* this is only done once in the initalisation the matrix is used to calculated th springLenght with blas3
	*/
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
	/**
	* \brief stores the normal spring lengths
	*
	* called in the constructor
	*/
	void populateSpringLengthsNormal() {
		springLengthsNormal[0] = MetaDataBase::DataBase()->getBodySpringLengthFrontLeft();
		springLengthsNormal[1] = MetaDataBase::DataBase()->getBodySpringLengthFrontLeft();
		springLengthsNormal[2] = MetaDataBase::DataBase()->getBodySpringLengthFrontRight();
		springLengthsNormal[3] = MetaDataBase::DataBase()->getBodySpringLengthFrontRight();
		springLengthsNormal[4] = MetaDataBase::DataBase()->getBodySpringLengthRearLeft();
		springLengthsNormal[5] = MetaDataBase::DataBase()->getBodySpringLengthRearLeft();
		springLengthsNormal[6] = MetaDataBase::DataBase()->getBodySpringLengthRearRight();
		springLengthsNormal[7] = MetaDataBase::DataBase()->getBodySpringInitialLengthRearRight();
	}
	/**
	* \brief stores the initial spring lengths
	*
	* called in the constructor
	*/
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

	/**
	* \brief initialise the vectors u_n, u_n_m_1, u_n_p_1
	*
	* called in the constructor
	* u_n = lenght_init -length_normal
	* u_n_m_1 = u_n - _h * velocity
	* u_n_p_1 = u_n
	* tempVector is used for velocity
	*/
	void initPosVectors() {
		// u_n = lenght_init -length_normal
		mkl<T>::vSub(Constants::NUM_LEGS, springLengths, springLengthsNormal, u_n);
		
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
		// u_n_m_1 = u_n - _h * velocity
		mkl<T>::copy(Constants::DOF, temp, 1, u_n_m_1, 1);
		mkl<T>::scal(Constants::DOF, -_h, u_n_m_1, 1);
		mkl<T>::axpy(Constants::DOF, 1, u_n, 1, u_n_m_1, 1);
		// u_n_p_1 = u_n
		mkl<T>::copy(Constants::DOF, u_n, 1, u_n_p_1, 1);
	}

public:
	/**
	* \brief Constructor
	*/
	Linear11dofBE(Car<T>* input_car, EVAALookup<T>* stiffnessInterpolation, EVAALookup<T>* dampingInterpolation): lookupStiffness(stiffnessInterpolation), lookupDamping(dampingInterpolation), _car(input_car) {
		u_n_m_1 = (T*)mkl_malloc(Constants::DOF * sizeof(T), Constants::ALIGNMENT); // velocity
		u_n = (T*)mkl_malloc(Constants::DOF * sizeof(T), Constants::ALIGNMENT); // position
		u_n_p_1 = (T*)mkl_calloc(Constants::DOF, sizeof(T), Constants::ALIGNMENT);

		A = (T*)mkl_malloc(Constants::DOFDOF * sizeof(T), Constants::ALIGNMENT);
		B = (T*)mkl_calloc(Constants::DOFDOF, sizeof(T), Constants::ALIGNMENT);
		M_h2 = (T*)mkl_calloc(Constants::DOFDOF, sizeof(T), Constants::ALIGNMENT);
		K = (T*)mkl_calloc(Constants::DOFDOF, sizeof(T), Constants::ALIGNMENT);
		D = (T*)mkl_calloc(Constants::DOFDOF, sizeof(T), Constants::ALIGNMENT);

		kVec = (T*)mkl_malloc(Constants::NUM_LEGS * 2 * sizeof(T), Constants::ALIGNMENT);
		dVec = (T*)mkl_malloc(Constants::NUM_LEGS * 2 * sizeof(T), Constants::ALIGNMENT);
		temp = (T*)mkl_malloc(Constants::DOF * sizeof(T), Constants::ALIGNMENT);		/**< used for the constructStiffnesMatrix and constructDampingMatrix as well as for the BE*/
		Mat_temp = (T*)mkl_calloc(Constants::DOFDOF, sizeof(T), Constants::ALIGNMENT); /**< used for the constructStiffnesMatrix and constructDampingMatrix*/
		Mat_springLength = (T*)mkl_calloc( 2 * Constants::NUM_LEGS * Constants::DOF, sizeof(T), Constants::ALIGNMENT);
		springLengths = (T*)mkl_malloc(Constants::NUM_LEGS * 2 * sizeof(T), Constants::ALIGNMENT);
		springLengthsNormal = (T*)mkl_malloc(Constants::NUM_LEGS * 2 * sizeof(T), Constants::ALIGNMENT);
		if (Constants::USEINTERPOLATION) {
			J = (T*)mkl_calloc(Constants::DOFDOF, sizeof(T), Constants::ALIGNMENT);
			dKdxx = (T*)mkl_calloc(Constants::DOFDOF, sizeof(T), Constants::ALIGNMENT);
			dDdxx = (T*)mkl_calloc(Constants::DOFDOF, sizeof(T), Constants::ALIGNMENT);
			dkdl = (T*)mkl_malloc(Constants::NUM_LEGS * 2 * sizeof(T), Constants::ALIGNMENT);
			dddl = (T*)mkl_malloc(Constants::NUM_LEGS * 2 * sizeof(T), Constants::ALIGNMENT);
			residual = (T*)mkl_malloc(Constants::DOF * sizeof(T), Constants::ALIGNMENT);
		}

		_h = MetaDataBase::DataBase()->getTimeStepSize();
		factor_h = 1 / _h;
		populateSpringLengthsNormal();
		initSpringLengths();
		constructMatSpringLenghts();

		// if interpolation is used we need alll the lengths to interpolate damping and stiffness if not we have to define them
		if (!Constants::USEINTERPOLATION) {
			// populate kVec
			kVec[0] = MetaDataBase::DataBase()->getBodyStiffnessFrontLeft();
			kVec[1] = MetaDataBase::DataBase()->getTyreStiffnessFrontLeft();
			kVec[2] = MetaDataBase::DataBase()->getBodyStiffnessFrontRight();
			kVec[3] = MetaDataBase::DataBase()->getTyreStiffnessFrontRight();
			kVec[4] = MetaDataBase::DataBase()->getBodyStiffnessRearLeft();
			kVec[5] = MetaDataBase::DataBase()->getTyreStiffnessRearLeft();
			kVec[6] = MetaDataBase::DataBase()->getBodyStiffnessRearRight();
			kVec[7] = MetaDataBase::DataBase()->getTyreStiffnessRearRight();
			// populate dVec
			dVec[0] = MetaDataBase::DataBase()->getBodyDampingFrontLeft();
			dVec[1] = MetaDataBase::DataBase()->getBodyDampingFrontLeft();
			dVec[2] = MetaDataBase::DataBase()->getBodyDampingFrontRight();
			dVec[3] = MetaDataBase::DataBase()->getBodyDampingFrontRight();
			dVec[4] = MetaDataBase::DataBase()->getBodyDampingRearLeft();
			dVec[5] = MetaDataBase::DataBase()->getBodyDampingRearLeft();
			dVec[6] = MetaDataBase::DataBase()->getBodyDampingRearRight();
			dVec[7] = MetaDataBase::DataBase()->getBodyDampingRearRight();
		}

		initPosVectors();
		constructMassMatrix();
		constructStiffnessMatrix();
		constructDampingMatrix();
		constructAMatrix();
		constructBMatrix();
	}
	/*
	Destructor
	*/
	virtual ~Linear11dofBE() {
		mkl_free(A);
		mkl_free(B);
		mkl_free(u_n);
		mkl_free(u_n_m_1);
		mkl_free(u_n_p_1);
		mkl_free(M_h2);
		mkl_free(K);
		mkl_free(D);
		mkl_free(kVec);
		mkl_free(dVec);
		mkl_free(temp);
		mkl_free(Mat_temp);
		mkl_free(springLengths);
		mkl_free(springLengthsNormal);
		mkl_free(Mat_springLength);
		if (Constants::USEINTERPOLATION) {
			mkl_free(J);
			mkl_free(dkdl);
			mkl_free(dddl);
			mkl_free(dKdxx);
			mkl_free(dDdxx);
			mkl_free(residual);
		}
	}

	/*
	Performs one timestep of the 11DOF solver
	\param load vector [angle:Z,GC:Y,W1:Y,T1:Y,W2:Y,T2:Y,...]
	*/
	virtual void update_step(T* force, T* solution) {
		if (Constants::USEINTERPOLATION) {
			MathLibrary::Solvers<T, Linear11dofBE>::Newton(this, force, J, residual, &res_norm, u_n_p_1, temp);
		}
		else {
			MathLibrary::Solvers<T, Linear11dofBE>::Linear_Backward_Euler(A, B, M_h2, u_n, u_n_m_1, force, u_n_p_1, Constants::DOF);
		}
		/*compute_normal_force(K, u_n_p_1, f_n_p_1, tyre_index_set, DOF, num_tyre);
		apply_normal_force(f_n_p_1, u_n_p_1, tyre_index_set, num_tyre);*/
		// solutio = solution[t_n] = u_n_p_1
		mkl<T>::copy(Constants::DOF, u_n_p_1, 1, solution, 1);

		//std::cout << "sol vector:" << std::endl;
		//MathLibrary::write_vector(u_n_p_1,11);


		MathLibrary::swap_address<T>(u_n, u_n_m_1); // u_n_m_1 points to u_n and u_n points to u_n_m_1
		MathLibrary::swap_address<T>(u_n_p_1, u_n); // u_n points to u_n_p_1 and u_n_p_1 point to u_n_m_1 now
	}
	/**
	* \brief calculate the residual(Newton Function) + res_norm
	*
	* residual = A*x[n+1] - B * x[n] + M_h2 * x[n-1] - forces
	* res_norm = norm(residual)
	*/
	void calcResidual(
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
	* \brief construct Jacobien
	*
	* this has to be called every newton iteraton
	* J = M_h2 + D / _h + K + (dDdx + dKdx)*x[n+1] - dDdx / _h * x[n]
	*/
	void constructJacobien() {
		// first update the derivative 
		lookupStiffness->getDerivative(springLengths, dkdl);
		lookupDamping->getDerivative(springLengths, dddl);
		// construct the derivative (tensor) times a pos vector
		constructLookupDerivativeX(dkdl, u_n_p_1, dKdxx);
		constructLookupDerivativeX(dddl, u_n_p_1, dDdxx);
		// J = A
		mkl<T>::copy(Constants::DOFDOF, A, 1, J, 1);
		// J += dDdx * x[n+1]
		mkl<T>::axpy(Constants::DOFDOF, 1.0, dDdxx, 1, J, 1);
		// J += dKdx * x[n+1]
		mkl<T>::axpy(Constants::DOFDOF, 1.0, dKdxx, 1, J, 1);
		// calc dDdxx with x[n] 
		constructLookupDerivativeX(dddl, u_n, dDdxx);
		// J += -1/_h * dDdx * x[n]
		mkl<T>::axpy(Constants::DOFDOF, -factor_h, dDdxx, 1, J, 1);
	}
	/**
	* \brief update all dependent matrices on the position vector
	*
	* only needed in case of lookup Tables
	*/
	void updateSystem() {
		constructSpringLengths();
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
		if (Constants::USEINTERPOLATION) {
			lookupStiffness->getInterpolation(springLengths, kVec);
		}
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
		if (Constants::USEINTERPOLATION) {
			lookupDamping->getInterpolation(springLengths, dVec);
		}
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
	void constructSpringLengths() {
		// springLengths = springLengthsNormal
		mkl<T>::copy(Constants::NUM_LEGS * 2, springLengthsNormal, 1, springLengths, 1);
		// springLengths += Mat_SpringLenghts * x
		mkl<T>::gemv(CblasRowMajor, CblasNoTrans, Constants::NUM_LEGS * 2, Constants::DOF, 1.0, Mat_springLength, Constants::DOF,
			u_n_p_1, 1.0, 1.0, springLengths, 1.0);
	}
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
class Linear11dofBDF2 : public Linear11dofBE<T> {
private:
	T* C, *D, *E;
	T* bVec; /**< this vector is used to store the whole constants part on the rhs apart from the force */
	T* u_n_m_2, *u_n_m_3;
	size_t time_step_count = 0;
	void(Linear11dofBDF2<T>::*_active_executor)(T*, T*);
	
	/**
	* \brief construct A
	*
	* A = 9/4h^2 * M + 3/2h * D + K
	*/
	void constructAMatrix() {
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
	* \brief construct B
	*
	* this has to be called every time step for interpolation and is only used in case of interpolation
	* B = 1/h^2 * M * (6 * x[n] - 11/2 * x[n-1] + 2 * x[n-2] - 1/4 * x[n-3]) + 1/h * D * (2 * x[n] - 1/2 * x[n-1])
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
	Linear11dofBDF2(Car<T>* input_car, EVAALookup<T>* stiffnessInterpolation, EVAALookup<T>* dampingInterpolation) : Linear11dofBE<T>(input_car, stiffnessInterpolation, dampingInterpolation) {
		C = (T*)mkl_calloc(Constants::DOFDOF,  sizeof(T), Constants::ALIGNMENT);
		D = (T*)mkl_calloc(Constants::DOFDOF, sizeof(T), Constants::ALIGNMENT);
		E = (T*)mkl_calloc(Constants::DOFDOF, sizeof(T), Constants::ALIGNMENT);
		u_n_m_2 = (T*)mkl_malloc(Constants::DOF * sizeof(T), Constants::ALIGNMENT);
		u_n_m_3 = (T*)mkl_malloc(Constants::DOF * sizeof(T), Constants::ALIGNMENT);
		bVec = (T*)mkl_malloc(Constants::DOF * sizeof(T), Constants::ALIGNMENT);
		_active_executor = &Linear11dofBDF2<T>::first_two_steps;

		if (!Constants::USEINTERPOLATION) {
			// C = (-11/2)*(1/(h*h))*M + (-1/(2*h))*D
			mkl<T>::axpy(Constants::DOFDOF, (-11.0 / 2.0), M_h2, 1, C, 1);
			mkl<T>::axpy(Constants::DOFDOF, (-1.0 / 2.0) * factor_h, D, 1, C, 1);

			// D = (2/(h*h))*M
			mkl<T>::axpy(Constants::DOFDOF, 2.0, M_h2, 1, D, 1);

			// E = (-1/4)*(1/(h*h))*M
			mkl<T>::axpy(Constants::DOFDOF, (-1.0 / 4.0), M_h2, 1, E, 1);

			mkl<T>::copy(Constants::DOF, _car->u_prev_linear, 1, u_n, 1);
			mkl<T>::copy(Constants::DOF, _car->velocity_current_linear, 1, u_n_m_1, 1);

			mkl<T>::scal(Constants::DOF, -_h, u_n_m_1, 1);
			mkl<T>::axpy(Constants::DOF, 1, u_n, 1, u_n_m_1, 1);
		}
	}

	void first_two_steps(T* force, T* solution) {
		
		if (time_step_count == 0) {
			mkl<T>::copy(Constants::DOF, u_n_m_1, 1, u_n_m_2, 1);
			Linear11dofBE<T>::update_step(force, solution);
		}
		else if (time_step_count == 1) {
			mkl<T>::copy(Constants::DOF, u_n_m_2, 1, u_n_m_3, 1);
			mkl<T>::copy(Constants::DOF, u_n_m_1, 1, u_n_m_2, 1);
			Linear11dofBE<T>::update_step(force, solution);
		}
		else {
			mkl<T>::copy(Constants::DOF, u_n_m_2, 1, u_n_m_3, 1);
			mkl<T>::copy(Constants::DOF, u_n_m_1, 1, u_n_m_2, 1);
			// construct A
			constructAMatrix();
			if (Constants::USEINTERPOLATION) {
				constructbVec();
			}
			else {
				// B = (6/(h*h))*M + (2/h)*D
				mkl<T>::scal(Constants::DOFDOF, 0.0, B, 1);
				mkl<T>::axpy(Constants::DOFDOF, 6, M_h2, 1, B, 1);
				mkl<T>::axpy(Constants::DOFDOF, 2 * factor_h, D, 1, B, 1);
			}

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
		if (Constants::USEINTERPOLATION) {
			MathLibrary::Solvers<T, Linear11dofBDF2>::Newton(this, force, J, residual, &res_norm, u_n_p_1, temp);
		}
		else {
			MathLibrary::Solvers<T, Linear11dofBDF2>::Linear_BDF2(A, B, C, D, E, u_n, u_n_m_1, u_n_m_2, u_n_m_3, force, u_n_p_1, Constants::DOF);
		}
		/*compute_normal_force(K, u_n_p_1, f_n_p_1, tyre_index_set, DOF, num_tyre);
		apply_normal_force(f_n_p_1, u_n_p_1, tyre_index_set, num_tyre);*/
		mkl<T>::copy(Constants::DOF, u_n_p_1, 1, solution, 1);
		MathLibrary::swap_address<T>(u_n_m_2, u_n_m_3); // u_n_m_2 points to u_n_m_3 and u_n_m_3 points to u_n_m_2
		MathLibrary::swap_address<T>(u_n_m_1, u_n_m_2); // u_n_m_2 points to u_n_m_1 and u_n_m_1 points to u_n_m_3
		MathLibrary::swap_address<T>(u_n, u_n_m_1); // u_n_m_1 points to u_n and u_n points to u_n_m_3
		MathLibrary::swap_address<T>(u_n_p_1, u_n); // u_n points to u_n_p_1 and u_n_p_1 point to u_n_m_3 now

	}
	/**
	* \brief calculate the residual(Newton Function) + res_norm
	*
	* residual = A*x[n+1] - bVec - forces
	* res_norm = norm(residual)
	*/
	void calcResidual(
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
	void updateSystem() {
		constructSpringLengths();
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
	void constructJacobien() {
		// first update the derivative 
		lookupStiffness->getDerivative(springLengths, dkdl);
		lookupDamping->getDerivative(springLengths, dddl);
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
	virtual ~Linear11dofBDF2() {
		mkl_free(C);
		mkl_free(D);
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
class Linear11dofFull : public Linear11dofBDF2<T> {
public:
	Linear11dofFull(Car<T>* input_car, EVAALookup<T>* stiffnessInterpolation, EVAALookup<T>* dampingInterpolation) :Linear11dofBDF2<T>(input_car, stiffnessInterpolation, dampingInterpolation) {

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
			update_step(f_n_p_1, sol_vect);

			iter++;
			t += _h;
		}
		//cblas_dcopy(DOF, u_sol + (iter - 1)*(DOF), 1, sol_vect, 1);
	}

	void print_final_results(T* sln) {
		std::cout.precision(15);
		std::cout << std::scientific;
		std::cout << "linear11DOF: orientation angles=\n\t[" << sln[1] << "\n\t " << sln[2] << "]" << std::endl;
		std::cout << "linear11DOF: car body position pc=\n\t[" << sln[0] << "]" << std::endl;
		std::cout << "linear11DOF: rear-right wheel position pw1=\n\t[" << sln[9] << "]" << std::endl;
		std::cout << "linear11DOF: rear-left wheel position pw2=\n\t[" << sln[7] << "]" << std::endl;
		std::cout << "linear11DOF: front-left wheel position pw3=\n\t[" << sln[3] << "]" << std::endl;
		std::cout << "linear11DOF: front-right wheel position pw4=\n\t[" << sln[5] << "]" << std::endl;
		std::cout << "linear11DOF: rear-right tyre position pt1=\n\t[" << sln[10] << "]" << std::endl;
		std::cout << "linear11DOF: rear-left tyre position pt2=\n\t[" << sln[8] << "]" << std::endl;
		std::cout << "linear11DOF: front-left tyre position pt3=\n\t[" << sln[4] << "]" << std::endl;
		std::cout << "linear11DOF: front-right tyre position pt4=\n\t[" << sln[6] << "]" << std::endl;
	}


	virtual ~Linear11dofFull() {
		mkl_free(u_sol);
		mkl_free(f_n_p_1);
	}

private:
	T tend_;
	T* u_sol, * f_n_p_1;
	int condition_type;
};

/*
For testing purposes
*/
template <typename T>
class Linear11dofFullBE : public Linear11dofBE<T> {
public:
	Linear11dofFullBE(Car<T>* input_car, EVAALookup<T>* stiffnessInterpolation, EVAALookup<T>* dampingInterpolation): Linear11dofBE<T>(input_car, stiffnessInterpolation, dampingInterpolation) {

		tend_ = MetaDataBase::DataBase()->getNumberOfTimeIterations() * _h;
		f_n_p_1 = (T*)mkl_malloc(Constants::DOF * sizeof(T), Constants::ALIGNMENT);


		f_n_p_1[0] = MetaDataBase::DataBase()->getBodyExternalForce()[2];
		f_n_p_1[3] = MetaDataBase::DataBase()->getWheelExternalForceFrontLeft()[2];
		f_n_p_1[4] = MetaDataBase::DataBase()->getTyreExternalForceFrontLeft()[2];
		f_n_p_1[5] = MetaDataBase::DataBase()->getWheelExternalForceFrontRight()[2];
		f_n_p_1[6] = MetaDataBase::DataBase()->getTyreExternalForceFrontRight()[2];
		f_n_p_1[7] = MetaDataBase::DataBase()->getWheelExternalForceFrontLeft()[2];
		f_n_p_1[8] = MetaDataBase::DataBase()->getTyreExternalForceFrontLeft()[2];
		f_n_p_1[9] = MetaDataBase::DataBase()->getWheelExternalForceFrontRight()[2];
		f_n_p_1[10] = MetaDataBase::DataBase()->getTyreExternalForceFrontRight()[2];

		//std::cout << "force vector:" << std::endl;
		//MathLibrary::write_vector(f_n_p_1,11);

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

		while (std::abs(t - (tend_ + _h)) > eps) {
			//solution_vect = u_sol + iter * (DOF);
			update_step(f_n_p_1, sol_vect);

			iter++;
			t += _h;
		}
		//cblas_dcopy(DOF, u_sol + (iter - 1)*(DOF), 1, sol_vect, 1);
	}

	void print_final_results(T* sln) {
		std::cout.precision(15);
		std::cout << std::scientific;
		std::cout << "linear11DOF: orientation angles=\n\t[" << sln[1] << "\n\t " << sln[2] << "]" << std::endl;
		std::cout << "linear11DOF: car body position pc=\n\t[" << sln[0] << "]" << std::endl;
		std::cout << "linear11DOF: rear-right wheel position pw1=\n\t[" << sln[9] << "]" << std::endl;
		std::cout << "linear11DOF: rear-left wheel position pw2=\n\t[" << sln[7] << "]" << std::endl;
		std::cout << "linear11DOF: front-left wheel position pw3=\n\t[" << sln[3] << "]" << std::endl;
		std::cout << "linear11DOF: front-right wheel position pw4=\n\t[" << sln[5] << "]" << std::endl;
		std::cout << "linear11DOF: rear-right tyre position pt1=\n\t[" << sln[10] << "]" << std::endl;
		std::cout << "linear11DOF: rear-left tyre position pt2=\n\t[" << sln[8] << "]" << std::endl;
		std::cout << "linear11DOF: front-left tyre position pt3=\n\t[" << sln[4] << "]" << std::endl;
		std::cout << "linear11DOF: front-right tyre position pt4=\n\t[" << sln[6] << "]" << std::endl;
	}


	virtual ~Linear11dofFullBE() {
		mkl_free(f_n_p_1);
	}

private:
	T tend_;
	T* u_sol, * f_n_p_1;
	int condition_type;
};
