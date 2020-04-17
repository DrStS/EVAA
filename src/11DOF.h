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
class Linear11dofParent {
protected:
	// main car object
	Car<T>* car_;												/**< pointer to car instance with all important car parameter */

	// define constants
	const int dim = 3;
	const size_t num_tyre = 4;
	const int dim_x_dim = dim * dim;
	const int num_wheels_x_dim = num_tyre * dim;
	const int DOF_diag = 9; // the diagonal elements from A

	// time step related
	T factor_h2;												/**< solution in next timestep */
	T factor_h;													/**< solution in next timestep */
	T h_;														/**< solution in next timestep */
	size_t DOF;													/**< degrees of freedom; DOF=11 */
	size_t mat_len;												/**< number of elements in an 11 dof matrix; mat_len=121 */														/**< solution in current timestep */

	// solution in next timestep
	T* u_n_p_1;

	// solution in previous timestep
	T* u_n_m_1;

	// solution in current timestep
	T* u_n;
	T*A, *B;                                                  /**< pointer to matrices used for the backward euler */
	size_t* tyre_index_set;

	/**
	* \brief Debug only: Checks whether the matrix is SPD
	*/
	void check_status(lapack_int status) {
		if (status == 1) {
			std::cout << "Matrix non Positive Definite" << std::endl;
			exit(5);
		}
		else if (status == -1) {
			std::cout << "Matrix contain illegal value" << std::endl;
			exit(6);
		}
	}
	/**
	* \brief Debug only: Outputs any vector
	*
	* \param vect to be printed
	* \param count its length
	*/
	void write_vector(T* vect, int count) {
		std::cout << "Debug mode print" << std::endl;
		for (size_t i = 0; i < count; ++i) {
			//std::cout << vect[i] << std::endl;
			std::cout.precision(15);
			std::cout << std::scientific << vect[i] << std::endl;
			//printf("%1.15f\n", vect[i]);
		}
	}

	/**
	* \brief Debug only: Outputs any square matrix
	*
	* \param vect to be printed
	* \param count its size is count x count
	*/
	void write_matrix(T* vect, int count) {
		std::cout << "Debug mode print" << std::endl;
		for (size_t i = 0; i < count; ++i) {
			//std::cout << vect[i] << std::endl;
			std::cout.precision(5);
			for (size_t j = 0; j < count; ++j) {
				std::cout << std::scientific << vect[i * count + j] << "  ";
			}
			std::cout << "\n" << std::endl;
			//printf("%1.15f\n", vect[i]);
		}
	}

	void compute_normal_force(T* K, T* u, T* force, size_t* index, size_t dim, size_t n) {
#pragma loop( ivdep )
		for (int i = 0; i < n; ++i) {
			force[index[i]] = -K[index[i] * dim + index[i]] * u[index[i]];
		}
	}

public:
	/*
	Constructor
	*/

	Linear11dofParent(Car<T>* input_car) {
		car_ = input_car;
		DOF = car_->DOF;
		mat_len = (DOF) * (DOF);
		u_n_m_1 = (T*)mkl_malloc(DOF * sizeof(T), Constants::ALIGNMENT); // velocity
		u_n = (T*)mkl_malloc(DOF * sizeof(T), Constants::ALIGNMENT); // position
		u_n_p_1 = (T*)mkl_malloc(DOF * sizeof(T), Constants::ALIGNMENT);
		A = (T*)mkl_malloc(mat_len * sizeof(T), Constants::ALIGNMENT);
		B = (T*)mkl_calloc(mat_len, sizeof(T), Constants::ALIGNMENT);
	}

	/*
	Intializes the solution vector in the timestep -1 (before the simulation starts)
	*/
	virtual void initialize_solver(T h)=0; 
	/*
	Performs one timestep of the 11DOF solver
	\param load vector [angle:Z,GC:Y,W1:Y,T1:Y,W2:Y,T2:Y,...]
	\return solution of the following timestep [angle:Z,GC:Y,W1:Y,T1:Y,W2:Y,T2:Y,...]
	*/
	virtual void update_step(T* force, T* solution)=0; 
	/*
	Destructor
	*/
	virtual ~Linear11dofParent() {
		mkl_free(A);
		mkl_free(B);
		mkl_free(u_n);
		mkl_free(u_n_m_1);
		mkl_free(u_n_p_1);
	}
};

template <typename T>
class Linear11dof : public Linear11dofParent<T>{
public:
	/*
	Constructor
	*/
	Linear11dof(Car<T>* input_car): Linear11dofParent<T>(input_car){}

	/*
	Intializes the solution vector in the timestep -1 (before the simulation starts)
	*/
	virtual void initialize_solver(T h) {
		h_ = h;
		factor_h2 = (1 / (h_ * h_));
		factor_h = (1 / (h_));

		// B=((2/(h*h))*M+(1/h)*D);
		mkl<T>::axpy(this->DOF, 2 * factor_h2, car_->M_linear, this->DOF + 1, B, this->DOF + 1);
		mkl<T>::axpy(mat_len, factor_h, car_->D, 1, B, 1);
		
		mkl<T>::copy(DOF, car_->u_prev_linear, 1, u_n, 1);

		// u_n_m_1 = u_n - h_ * velocity_current_linear
		mkl<T>::copy(DOF, car_->velocity_current_linear, 1, u_n_m_1, 1);
		mkl<T>::scal(DOF, -h_, u_n_m_1, 1);
		mkl<T>::axpy(DOF, 1, u_n, 1, u_n_m_1, 1);
	}
	/*
	Performs one timestep of the 11DOF solver
	\param load vector [angle:Z,GC:Y,W1:Y,T1:Y,W2:Y,T2:Y,...]
	\return solution of the following timestep [angle:Z,GC:Y,W1:Y,T1:Y,W2:Y,T2:Y,...]
	*/
	virtual void update_step(T* force, T* solution) {
		// construct A
		// A = M_linear (also acts as initialization of A)
		mkl<T>::copy(mat_len, car_->M_linear, 1, A, 1);
		// A = 1/h^2 * A => A = 1/h^2 * M
		mkl<T>::scal(this->DOF, factor_h2, A, this->DOF + 1); // use the fact that M is diagonal
		// A += 1/h * D => A = 1/h^2 * M + 1/h * D
		mkl<T>::axpy(mat_len, factor_h, car_->D, 1, A, 1);
		// A += K => A = 1/h^2 * M + 1/h * D + K
		mkl<T>::axpy(mat_len, 1, car_->K, 1, A, 1);

		// u_n_m_1 = -1/h^2 * u_n_m_1
		mkl<T>::scal(DOF, -factor_h2, u_n_m_1, 1);
		//MathLibrary::write_vector(force, 11);
		//cblas_dscal(DOF, 0.0, force, 1);
		MathLibrary::Solvers<T, Linear11dof>::Linear_Backward_Euler(A, B, car_->M_linear, u_n, u_n_m_1, force, u_n_p_1, DOF);
		/*compute_normal_force(K, u_n_p_1, f_n_p_1, tyre_index_set, DOF, num_tyre);
		apply_normal_force(f_n_p_1, u_n_p_1, tyre_index_set, num_tyre);*/
		// solutio = solution[t_n] = u_n_p_1
		mkl<T>::copy(DOF, u_n_p_1, 1, solution, 1);
		MathLibrary::swap_address<T>(u_n, u_n_m_1); // u_n_m_1 points to u_n and u_n points to u_n_m_1
		MathLibrary::swap_address<T>(u_n_p_1, u_n); // u_n points to u_n_p_1 and u_n_p_1 point to u_n_m_1 now
	}
	/*
	Destructor
	*/
	virtual ~Linear11dof() {}
};

template <typename T>
class Linear11dofBDF2 : public Linear11dof<T> {
protected:
	T* C, *D, *E;
	T* u_n_m_2, *u_n_m_3;
	size_t time_step_count = 0;
	void(Linear11dofBDF2<T>::*_active_executor)(T*, T*);

public:
	/*
	Constructor
	*/
	Linear11dofBDF2(Car<T>* input_car) : Linear11dof<T>(input_car) {
		C = (T*)mkl_calloc(mat_len,  sizeof(T), Constants::ALIGNMENT);
		D = (T*)mkl_calloc(mat_len, sizeof(T), Constants::ALIGNMENT);
		E = (T*)mkl_calloc(mat_len, sizeof(T), Constants::ALIGNMENT);
		u_n_m_2 = (T*)mkl_malloc(DOF * sizeof(T), Constants::ALIGNMENT);
		u_n_m_3 = (T*)mkl_malloc(DOF * sizeof(T), Constants::ALIGNMENT);
		_active_executor = &Linear11dofBDF2<T>::first_two_steps;
	}

	/*
	Intializes the solution vector in the timestep -1 (before the simulation starts)
	*/
	virtual void initialize_solver(T h) {
		h_ = h;
		Linear11dof<T>::initialize_solver(h_);
	}

	void initialize_solver_bdf2(T h) {
		h_ = h;
		factor_h2 = (1 / (h_ * h_));
		factor_h = (1 / (h_));

		// B = (6/(h*h))*M + (2/h)*D
		mkl<T>::scal(mat_len, 0.0, B, 1);
		mkl<T>::axpy(mat_len, 6 * factor_h2, car_->M_linear, 1, B, 1);
		mkl<T>::axpy(mat_len, 2 * factor_h, car_->D, 1, B, 1);

		// C = (-11/2)*(1/(h*h))*M + (-1/(2*h))*D
		mkl<T>::axpy(mat_len, (-11.0/2.0) * factor_h2, car_->M_linear, 1, C, 1);
		mkl<T>::axpy(mat_len, (-1.0/2.0) * factor_h, car_->D, 1, C, 1);

		// D = (2/(h*h))*M
		mkl<T>::axpy(mat_len, 2.0 * factor_h2, car_->M_linear, 1, D, 1);

		// E = (-1/4)*(1/(h*h))*M
		mkl<T>::axpy(mat_len, (-1.0/4.0) * factor_h2, car_->M_linear, 1, E, 1);

		mkl<T>::copy(DOF, car_->u_prev_linear, 1, u_n, 1);
		mkl<T>::copy(DOF, car_->velocity_current_linear, 1, u_n_m_1, 1);

		mkl<T>::scal(DOF, -h_, u_n_m_1, 1);
		mkl<T>::axpy(DOF, 1, u_n, 1, u_n_m_1, 1);

	}

	void first_two_steps(T* force, T* solution) {
		
		if (time_step_count == 0) {
			mkl<T>::copy(DOF, u_n_m_1, 1, u_n_m_2, 1);
			Linear11dof<T>::update_step(force, solution);
			time_step_count += 1;
		}
		else if (time_step_count == 1) {
			mkl<T>::copy(DOF, u_n_m_2, 1, u_n_m_3, 1);
			mkl<T>::copy(DOF, u_n_m_1, 1, u_n_m_2, 1);
			Linear11dof<T>::update_step(force, solution);
			time_step_count += 1;
		}
		else {
			time_step_count += 1;
			mkl<T>::copy(DOF, u_n_m_2, 1, u_n_m_3, 1);
			mkl<T>::copy(DOF, u_n_m_1, 1, u_n_m_2, 1);
			initialize_solver_bdf2(h_);
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
		// construct A
		mkl<T>::copy(mat_len, car_->M_linear, 1, A, 1);
		mkl<T>::scal(mat_len, (9.0/4.0)*factor_h2, A, 1);
		mkl<T>::axpy(mat_len, (3.0/2.0)*factor_h, car_->D, 1, A, 1);
		mkl<T>::axpy(mat_len, 1, car_->K, 1, A, 1);
		//cblas_dscal(DOF, 0.0, force, 1);
		MathLibrary::Solvers<T, Linear11dofBDF2>::Linear_BDF2(A, B, C, D, E, u_n, u_n_m_1, u_n_m_2, u_n_m_3, force, u_n_p_1, DOF);
		/*compute_normal_force(K, u_n_p_1, f_n_p_1, tyre_index_set, DOF, num_tyre);
		apply_normal_force(f_n_p_1, u_n_p_1, tyre_index_set, num_tyre);*/
		mkl<T>::copy(DOF, u_n_p_1, 1, solution, 1);
		MathLibrary::swap_address<T>(u_n_m_2, u_n_m_3); // u_n_m_2 points to u_n_m_3 and u_n_m_3 points to u_n_m_2
		MathLibrary::swap_address<T>(u_n_m_1, u_n_m_2); // u_n_m_2 points to u_n_m_1 and u_n_m_1 points to u_n_m_3
		MathLibrary::swap_address<T>(u_n, u_n_m_1); // u_n_m_1 points to u_n and u_n points to u_n_m_3
		MathLibrary::swap_address<T>(u_n_p_1, u_n); // u_n points to u_n_p_1 and u_n_p_1 point to u_n_m_3 now

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
	}
};



/*
For testing purposes
*/
template <typename T>
class Linear11dofFull : public Linear11dofBDF2<T> {
public:
	Linear11dofFull(Car<T>* input_car) :Linear11dofBDF2<T>(input_car) {

		h_ = MetaDataBase::DataBase()->getTimeStepSize();
		tend_ = MetaDataBase::DataBase()->getNumberOfTimeIterations() * h_;
		int sol_size = (floor(tend_ / h_) + 1);
		f_n_p_1 = (T*)mkl_malloc(DOF * sizeof(T), Constants::ALIGNMENT);
		u_sol = (T*)mkl_calloc((sol_size + 1) * (DOF), sizeof(T), Constants::ALIGNMENT);


		interpolation_enabled = MetaDataBase::DataBase()->getUseInterpolation();

		f_n_p_1[0] = MetaDataBase::DataBase()->getBodyExternalForce()[1];
		f_n_p_1[3] = MetaDataBase::DataBase()->getWheelExternalForceFrontLeft()[1];
		f_n_p_1[4] = MetaDataBase::DataBase()->getTyreExternalForceFrontLeft()[1];
		f_n_p_1[5] = MetaDataBase::DataBase()->getWheelExternalForceFrontRight()[1];
		f_n_p_1[6] = MetaDataBase::DataBase()->getTyreExternalForceFrontRight()[1];
		f_n_p_1[7] = MetaDataBase::DataBase()->getWheelExternalForceFrontLeft()[1];
		f_n_p_1[8] = MetaDataBase::DataBase()->getTyreExternalForceFrontLeft()[1];
		f_n_p_1[9] = MetaDataBase::DataBase()->getWheelExternalForceFrontRight()[1];
		f_n_p_1[10] = MetaDataBase::DataBase()->getTyreExternalForceFrontRight()[1];


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
		initialize_solver(h_);
		int iter = 1;
		T t = h_;
		double eps = h_ / 100;
		T* solution_vect = u_sol;

		//cblas_dcopy(DOF, car_->u_prev_linear, 1, solution_vect, 1);
		while (std::abs(t - (tend_ + h_)) > eps) {
			if (interpolation_enabled) {
				(car_)->update_lengths_11DOF();
				(car_)->lookupStiffness->getInterpolation((car_)->current_spring_length, (car_)->k_vec);
				(car_)->update_K((car_)->k_vec);
			}
			//solution_vect = u_sol + iter * (DOF);
			update_step(f_n_p_1, sol_vect);

			iter++;
			t += h_;
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
	T h_;
	int condition_type;
	bool interpolation_enabled;
};
