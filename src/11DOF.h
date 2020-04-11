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
#include "ReadXML.h"

template <typename T>
class Linear11dof {
protected:
	// main car object
	Car<T>* car_;

	// define constants
	const int dim = 3;
	const size_t num_tyre = 4;
	const int dim_x_dim = dim * dim;
	const int num_wheels_x_dim = num_tyre * dim;
	const int DOF_diag = 9; // the diagonal elements from A

	// time step related
	T factor_h2;
	T factor_h;
	T h_;

	// DOF=11
	size_t DOF;

	// mat_len=121
	size_t mat_len;

	// solution in next timestep
	T* u_n_p_1;

	// solution in previous timestep
	T* u_n_m_1;

	// solution in current timestep
	T* u_n;

	T* A, * B;

	size_t* tyre_index_set;

	/*
	Debug only: Checks whether the matrix is SPD
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
	/*
	Debug only: Outputs any vector
	\param vect to be printed
	\param count its length
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

	/*
	Debug only: Outputs any square matrix
	\param vect to be printed
	\param count its size is count x count
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
	Linear11dof(Car<T>* input_car) {
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
	void initialize_solver(T h) {
		h_ = h;
		factor_h2 = (1 / (h_ * h_));
		factor_h = (1 / (h_));

		// B=((2/(h*h))*M+(1/h)*D);
		cblas_daxpy(this->DOF, 2 * factor_h2, car_->M_linear, this->DOF + 1, B, this->DOF + 1);
		cblas_daxpy(mat_len, factor_h, car_->D, 1, B, 1);
		
		cblas_dcopy(DOF, car_->u_prev_linear, 1, u_n, 1);

		// u_n_m_1 = u_n - h_ * velocity_current_linear
		cblas_dcopy(DOF, car_->velocity_current_linear, 1, u_n_m_1, 1);
		cblas_dscal(DOF, -h_, u_n_m_1, 1);
		cblas_daxpy(DOF, 1, u_n, 1, u_n_m_1, 1);
	}
	/*
	Performs one timestep of the 11DOF solver
	\param load vector [angle:Z,GC:Y,W1:Y,T1:Y,W2:Y,T2:Y,...]
	\return solution of the following timestep [angle:Z,GC:Y,W1:Y,T1:Y,W2:Y,T2:Y,...]
	*/
	void update_step(T* force, T* solution) {
		// construct A
		// A = M_linear (also acts as initialization of A)
		cblas_dcopy(mat_len, car_->M_linear, 1, A, 1);
		// A = 1/h^2 * A => A = 1/h^2 * M
		cblas_dscal(this->DOF, factor_h2, A, this->DOF + 1); // use the fact that M is diagonal
		// A += 1/h * D => A = 1/h^2 * M + 1/h * D
		cblas_daxpy(mat_len, factor_h, car_->D, 1, A, 1);
		// A += K => A = 1/h^2 * M + 1/h * D + K
		cblas_daxpy(mat_len, 1, car_->K, 1, A, 1);

		// u_n_m_1 = -1/h^2 * u_n_m_1
		cblas_dscal(DOF, -factor_h2, u_n_m_1, 1);
		//cblas_dscal(DOF, 0.0, force, 1);
		MathLibrary::Solvers<T, Linear11dof>::Linear_Backward_Euler(A, B, car_->M_linear, u_n, u_n_m_1, force, u_n_p_1, DOF);
		/*compute_normal_force(K, u_n_p_1, f_n_p_1, tyre_index_set, DOF, num_tyre);
		apply_normal_force(f_n_p_1, u_n_p_1, tyre_index_set, num_tyre);*/
		// solutio = solution[t_n] = u_n_p_1
		cblas_dcopy(DOF, u_n_p_1, 1, solution, 1);
		MathLibrary::swap_address<T>(u_n, u_n_m_1); // u_n_m_1 points to u_n and u_n points to u_n_m_1
		MathLibrary::swap_address<T>(u_n_p_1, u_n); // u_n points to u_n_p_1 and u_n_p_1 point to u_n_m_1 now
	}
	/*
	Destructor
	*/
	virtual ~Linear11dof() {
		mkl_free(A);
		mkl_free(B);
		mkl_free(u_n);
		mkl_free(u_n_m_1);
		mkl_free(u_n_p_1);
	}
};

/*
For testing purposes
*/
template <typename T>
class Linear11dofFull : public Linear11dof<T> {
public:
	Linear11dofFull(const Simulation_Parameters& params, const Load_Params& load_param, Car<T>* input_car) :Linear11dof<T>(input_car) {

		h_ = params.timestep;
		tend_ = params.num_time_iter * h_;
		int sol_size = (floor(tend_ / h_) + 1);
		f_n_p_1 = (T*)mkl_malloc(DOF * sizeof(T), Constants::ALIGNMENT);
		u_sol = (T*)mkl_calloc((sol_size + 1) * (DOF), sizeof(T), Constants::ALIGNMENT);


		interpolation_enabled = params.interpolation;

		f_n_p_1[0] = load_param.external_force_body[2];
		f_n_p_1[3] = load_param.external_force_wheel[2 * 3 + 2];
		f_n_p_1[4] = load_param.external_force_tyre[2 * 3 + 2];
		f_n_p_1[5] = load_param.external_force_wheel[3 * 3 + 2];
		f_n_p_1[6] = load_param.external_force_tyre[3 * 3 + 2];
		f_n_p_1[7] = load_param.external_force_wheel[1 * 3 + 2];
		f_n_p_1[8] = load_param.external_force_tyre[1 * 3 + 2];
		f_n_p_1[9] = load_param.external_force_wheel[0 * 3 + 2];
		f_n_p_1[10] = load_param.external_force_tyre[0 * 3 + 2];


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
				(car_)->lookupStiffness->getStiffness((car_)->current_spring_length, (car_)->k_vec);
				(car_)->update_K((car_)->k_vec);
			}
			//solution_vect = u_sol + iter * (DOF);
			update_step(f_n_p_1, sol_vect);

			iter++;
			t += h_;
		}
		//cblas_dcopy(DOF, u_sol + (iter - 1)*(DOF), 1, sol_vect, 1);


	}

	/*
	Prints all positions and angles in the car object
	*/
	void print_final_results(T* sln) {

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
