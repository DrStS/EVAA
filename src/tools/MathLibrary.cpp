/*  Copyright &copy; 2019, Dr. Stefan Sicklinger, Munich \n
*
*  All rights reserved.
*
*  This file is part of EVAA.
*
*  EVAA is free software: you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation, either version 3 of the License, or
*  (at your option) any later version.
*
*  EVAA is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with EVAA.  If not, see http://www.gnu.org/licenses/.
*/
#include <iostream>
#include "MathLibrary.h"
#include  <mkl_cblas.h>
#include <mkl.h>

namespace MathLibrary {

	template <typename T>
	void transmutate_elements(T* arr, size_t start, size_t end, T value, T tol) {
		for (size_t i = start; i < end; ++i) {
			if (arr[i] * arr[i] < tol*tol) {
				arr[i] = value;
			}
		}
	}

	template <typename T>
	void diagonal_matrix(T* mat, size_t dim, T val) {
		for (size_t i = 0; i < dim; ++i) {
			mat[i*dim + i] = val;
		}
	}

	template <typename T>
	void elementwise_inversion(T* vect, size_t dim) {
		for (size_t i = 0; i < dim; ++i) {
			vect[i] = 1.0 / vect[i];
		}
	}

	template <typename T>
	void swap_address(T*& a, T*& b)
	{
		T* c = a;
		a = b;
		b = c;
	}

	void printMKLInfo() {
#ifdef USE_INTEL_MKL 
		int len = 198;
		char buf[198];
		mkl_get_version_string(buf, len);
		std::cout << buf << std::endl;
#endif
	}

	void computeDenseSymLUFactorisation(const int nElements, std::vector<double> &A, std::vector<int> &pivots) {
#ifdef USE_INTEL_MKL 
		LAPACKE_dgetrf(LAPACK_COL_MAJOR, nElements, nElements, A.data(), nElements, pivots.data());
#endif
	}

	void computeDenseSymSolution(const int nElements, std::vector<double> &A, std::vector<int> &pivots, std::vector<double> &rhs) {
#ifdef USE_INTEL_MKL 
		mkl_set_num_threads(EVAA::AuxiliaryParameters::denseVectorMatrixThreads);
		LAPACKE_dgetrs(LAPACK_COL_MAJOR, 'N', nElements, 1, A.data(), nElements, pivots.data(), rhs.data(), nElements);
#endif
	}

	void computeDenseVectorAddition(double *vec1, double *vec2, const double a, const int elements) {
#ifdef USE_INTEL_MKL
		mkl_set_num_threads(EVAA::AuxiliaryParameters::denseVectorMatrixThreads);
		cblas_daxpy(elements, a, vec1, 1, vec2, 1);
#endif
	}
	template <typename T>
	void Broyden_Euler(void(*f)(T, T*, T*), T* t, T* x_previous, T tol, T* x_vector_new, size_t max_iter) {
		int alignment = 64;
		size_t t_len = sizeof(t) / sizeof(t[0]);
		size_t x_len = sizeof(x_previous) / sizeof(x_previous[0]);
		T* f_old = (T*)mkl_malloc(sizeof(T) * x_len, alignment);
		T* f_new = (T*)mkl_malloc(sizeof(T) * x_len, alignment);
		T* dx = (T*)mkl_malloc(sizeof(T) * x_len, alignment);
		T* dx_inv = (T*)mkl_malloc(sizeof(T) * x_len, alignment);
		T* df = (T*)mkl_malloc(sizeof(T) * x_len, alignment);
		T* x = (T*)mkl_malloc(sizeof(T) * x_len, alignment);
		T* F = (T*)mkl_malloc(sizeof(T) * x_len, alignment);
		T* F_new = (T*)mkl_malloc(sizeof(T) * x_len, alignment);
		T* dF = (T*)mkl_malloc(sizeof(T) * x_len, alignment);
		T* J = (T*)mkl_calloc(x_len*x_len, sizeof(T), alignment);
		T* J_tmp = (T*)mkl_calloc(x_len*x_len, sizeof(T), alignment);
		lapack_int* piv = (lapack_int*)mkl_malloc(sizeof(lapack_int*) * x_len, alignment);

		T dt;
		T *it_start, *curr_time_start_pos;
		T* x_new = (T*)mkl_malloc(sizeof(T) * x_len, alignment);
		T* switcher;			// used for interchanging adresses between F and F_new
		T eps = 0.001;
		double nrm = 0.001;

		cblas_dcopy(x_len, x_previous, 1, x_vector_new, 1);
		for (size_t i = 1; i < t_len; ++i) {
			dt = t[i] - t[i - 1];
			// set the pointer to the start of previous timestep
			it_start = x_vector_new + (i - 1)*x_len;

			// 1. Initialize guess from previous time step
			// f_old = f(t(n-1), x_previous');
			f(t[i - 1], it_start, f_old);

			// in case the velocity is 0 add nuggets to avoid singular matrices (slow check for improvement)
			// f_old(abs(f_old) < 0.001) = 0.001;
			transmutate_elements(f_old, 0, x_len, eps, eps);

			// 2. Initial guess from Explicit Euler
			// x = x_previous + dt * f_old';
			// f_new = f(t(n), x');
			cblas_dcopy(x_len, it_start, 1, x, 1);
			cblas_daxpy(x_len, dt, f_old, 1, x, 1);
			f(t[i], x, f_new);

			// Initial approximation of the Jacobian
			// dx = x - x_previous;
			cblas_dcopy(x_len, x, 1, dx, 1);
			cblas_daxpy(x_len, -1.0, it_start, 1, dx, 1);
			// df = f_new' - f_old';
			cblas_dcopy(x_len, f_new, 1, df, 1);
			cblas_daxpy(x_len, -1.0, f_old, 1, df, 1);

			// approximate J(x_0)
			// J = eye(length(df)) - delta_t*((1./dx)'*df)';
			diagonal_matrix(J, x_len, 1.0);
			cblas_dcopy(x_len, dx, 1, dx_inv, 1);
			elementwise_inversion(dx_inv, x_len);
			cblas_dger(CblasRowMajor, x_len, x_len, -dt, dx_inv, 1, df, 1, J, x_len);

			// calculate initial F for stopping condition
			// x_dot = f_new';
			// F = dx - delta_t * x_dot;
			cblas_dcopy(x_len, dx, 1, F, 1);
			cblas_daxpy(x_len, -dt, f_new, 1, F, 1);

			// Broyden's Method
			for (size_t j = 0; j < max_iter; ++j) {
				if (cblas_dnrm2(x_len, F, 1) < tol) {
					break;
				}

				// x(i+1) = x(i) - J^(-1)*g(x(i))
				cblas_dcopy(x_len*x_len, J, 1, J_tmp, 1);
				LAPACKE_dgetrf(LAPACK_ROW_MAJOR, x_len, x_len, J_tmp, x_len, piv);
				cblas_dcopy(x_len, F, 1, x_new, 1);
				LAPACKE_dgetrs(LAPACK_ROW_MAJOR, "N", x_len, 1, J_tmp, x_len, piv, x_new, 1);
				cblas_daxpy(x_len, -1.0, x, 1, x_new, 1); // result here is -x_new
				cblas_dscal(x_len, -1.0, x_new, 1);

				// Calculate new derivative
				f(t[i - 1], x_new, f_new);

				// F_new = x_new - x_previous - delta_t * x_dot
				cblas_dcopy(x_len, x_new, 1, F_new, 1);
				cblas_daxpy(x_len, -1.0, it_start, 1, F_new, 1);
				cblas_daxpy(x_len, -dt, f_new, 1, F_new, 1);

				// dF = (F_new - F)'
				cblas_dcopy(x_len, F_new, 1, dF, 1);
				cblas_daxpy(x_len, -1.0, F, 1, dF, 1);

				// dx = (x_new - x)';
				cblas_dcopy(x_len, x_new, 1, dx, 1);
				cblas_daxpy(x_len, -1.0, x, 1, dx, 1);

				// J(n+1) = J(n) + (dF - J * dx) * dx' / norm(dx)^2
				// scaling dx and dF by norm(dx)
				nrm = cblas_dnrm2(x_len, dx, 1);
				cblas_dscal(x_len, 1.0 / nrm, dx, 1);
				cblas_dscal(x_len, 1.0 / nrm, dF, 1);
				// y := alpha*A*x + beta*y, dgemv operation
				cblas_dgemv(CblasRowMajor, CblasNoTrans, x_len, x_len, -1.0, J, x_len, dx, 1, 1.0, dF, 1); // using dF to store data
				cblas_dger(CblasRowMajor, x_len, x_len, 1.0, dF, 1, dx, 1, J, x_len);

				// F = F_new; interchanging pointers to avoid copy
				swap_address(F, F_new);
				// x = x_new; interchanging pointers to avoid copy
				swap_address(x, x_new);
			}
			// x_vector_new(n,:) = x;
			// start position of new time step
			curr_time_start_pos = x_vector_new + i * x_len;
			cblas_dcopy(x_len, x, 1, curr_time_start_pos, 1);
		}
	}
} /* namespace Math */


