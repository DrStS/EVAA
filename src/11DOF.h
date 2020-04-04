#pragma once
#include <mkl.h>
#include "ReadXML.h"
#include "MathLibrary.h"
#include "car.h"

template <typename T>
class linear11dof {
private:
	Car* Car_
	// define constants
	const int alignment = 64;
	const int dim = 3;
	const int num_wheels = 4;
	const int dim_x_dim = dim * dim;
	const int num_wheels_x_dim = num_wheels * dim;
	const int DOF_diag = 9; // the diagonal elements from A

	//// Solver type selection based on type of boundary condition
	std::string condition_type;
	std::string cond1 = "fixed_to_road";
	std::string cond2 = "road_force";

	T* u_sol, *u_sol_red, *u_n_p_1, *u_n_p_1_red, *u_n_m_1, *u_n, *u_n_red, *u_n_m_1_red, *A, *Ared, *B, *Bred, *f_n_p_1, *f_n_p_1_red;
	size_t* tyre_index_set;
	size_t num_tyre = 4;
	void check_status(lapack_int status) {
		if (status == 1) {
			std::cout << "Matrix non Positive Definite" << std::endl;
			exit(5);
		}
		else if (status == -1) {
			std::cout << "Matrix contain illegal value" << std::endl;
			exit(5);
		}
	}
	void write_vector(T* vect, int count) {
		std::cout << "Debug mode print" << std::endl;
		for (size_t i = 0; i < count; ++i) {
			//std::cout << vect[i] << std::endl;
			std::cout.precision(15);
			std::cout << std::scientific << vect[i] << std::endl;
			//printf("%1.15f\n", vect[i]);
		}
	}

	void write_matrix(T* vect, int count) {
		std::cout << "Debug mode print" << std::endl;
		for (size_t i = 0; i < count; ++i) {
			//std::cout << vect[i] << std::endl;
			std::cout.precision(5);
			for (size_t j = 0; j < count; ++j) {
				std::cout << std::scientific << vect[i*count + j] << "  ";
			}
			std::cout << "\n" << std::endl;
			//printf("%1.15f\n", vect[i]);
		}
	}

public:
	linear11dof(Car* input_car) {

	}
};