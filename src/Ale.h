#pragma once
#include "Car.h"
#include "Profile.h"
#include "Load_module.h"
#include "EVAAComputeEngine.h"

template <class T>
class ALE {
private:
	Car* Car_obj = NULL; // suppose Interpolation in the Car
	Load_module* Load_module_obj = NULL; // needs Profile and Car
	linear11dof* linear11dof_obj = NULL;
	// interpolator member of EVAA

	time_vec;
		u_sol;

public:
	ALE(const Car* Car_obj_val,
		const Profile* Profile_obj_val,
		const Simulation_Parameters& params_val,
		const Load_Params& load_param_val,
		EVAAComputeStiffness* interpolator_val) {
		
		Car_obj = Car_obj_val;
		Profile_obj = Profile_obj_val;
		Load_module_obj(Profile_obj_val, Car_obj_val);

		params = params_val;
		load_param = load_param_val;
		interpolator = interpolator_val;
		linear11dof_obj = linear11dof(params_val, load_param_val, interpolator_val);
	}

	// class for boundary conditions (road forces maybe)

	// solver functions

	// get structure functions
};


void solve_full(T* sol_vect) {
	///////////////////////////////////// For 7 DOF System //////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////// Memory Allocation for Intermediate Solution Vectors //////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////
	int sol_size = (floor(tend_ / h_) + 1);
	T factor_h2 = (1 / (h_ * h_));
	T factor_h = (1 / (h_));
	const int mat_len = (DOF) * (DOF);
	time = (T*)mkl_calloc(sol_size, sizeof(T), alignment);
	u_sol = (T*)mkl_calloc(sol_size * (DOF), sizeof(T), alignment);

	// A=((1/(h*h))*M+(1/h)*D+K);

	// B=((2/(h*h))*M+(1/h)*D);
	cblas_daxpy(mat_len, 2 * factor_h2, M, 1, B, 1);
	cblas_daxpy(mat_len, factor_h, D, 1, B, 1);
	int iter = 1;
	T t = h_;
	double eps = h_ / 100;
	/*auto start = std::chrono::steady_clock::now();*/
	while (std::abs(t - (tend_ + h_)) > eps) {

		// K update here
		cblas_dcopy(mat_len, M, 1, A, 1);
		cblas_dscal(mat_len, factor_h2, A, 1);
		cblas_daxpy(mat_len, factor_h, D, 1, A, 1);
		cblas_daxpy(mat_len, 1, K, 1, A, 1);
		// Cholesky factorization of A
		lapack_int status;
		//lapack_int* piv = (lapack_int*)mkl_calloc((DOF + 4), sizeof(lapack_int), alignment);
		//write_matrix(A, DOF);
		status = LAPACKE_dpotrf(LAPACK_ROW_MAJOR, 'L', DOF, A, DOF);
		//status = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, DOF + 4, DOF + 4, A, DOF + 4, piv);
		check_status(status);

		// u_n_p_1=A\(B*u_n-((1/(h*h))*M)*u_n_m_1+f_n_p_1);
		// u_n_p_1 = B*u_n

		cblas_dgemv(CblasRowMajor, CblasNoTrans, DOF, DOF, 1, B, DOF, u_n, 1, 0, u_n_p_1, 1);
		// u_n_p_1 = -((1/(h*h))*M)*u_n_m_1 + u_n_p_1
		cblas_dgemv(CblasRowMajor, CblasNoTrans, DOF, DOF, -factor_h2, M, DOF, u_n_m_1, 1, 1, u_n_p_1, 1);
		// u_n_p_1 = f_n_p_1 + u_n_p_1
		cblas_daxpy(DOF, 1, f_n_p_1, 1, u_n_p_1, 1);
		// u_n_p_1=A\u_n_p_1    ====> u_n_p_1 computed
		LAPACKE_dpotrs(LAPACK_ROW_MAJOR, 'L', DOF, 1, A, DOF, u_n_p_1, 1);
		//LAPACKE_dgetrs(LAPACK_ROW_MAJOR, 'N', DOF + 4, 1, A, DOF + 4, piv, u_n_p_1, 1);

		///////////////// Debug print statements ///////////////////////////////
		/*if (iter==2){
			write_vector(u_n_p_1, DOF + 4);
		}*/
		time[iter] = t;
		/*
		f_n_p_1(idx) = -K(idx,idx)*u_n_p_1(idx);
		f_n_p_1(idx) = f_update(f_n_p_1, idx);
		u_n_p_1(idx) = 0;
		*/
		////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////// Normal force computation here /////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// compute_normal_force(K, u_n_p_1, f_n_p_1, tyre_index_set, DOF, num_tyre);


		///////////////// Debug print statements ///////////////////////////////
		/*if (iter==10){
			write_vector(u_n_p_1, DOF);
		}*/



		// apply_normal_force(f_n_p_1, u_n_p_1, tyre_index_set, num_tyre);

		////////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////////////////////////////

		// u_sol(j,:)=u_n_p_1;
		cblas_dcopy(DOF, u_n_p_1, 1, u_sol + iter * (DOF), 1);
		// u_n_m_1=u_n;
		cblas_dcopy(DOF, u_n, 1, u_n_m_1, 1);
		// u_n    =u_n_p_1;
		cblas_dcopy(DOF, u_n_p_1, 1, u_n, 1);
		iter++;
		t += h_;
	}
	/*auto end = std::chrono::steady_clock::now();
	std::cout << "Elapsed time in nanoseconds : "
		<< std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
		<< " ms" << std::endl;*/
		//std::cout << "iter = " << iter << " sol_size = "<< sol_size <<"\n\n" << std::endl;
	cblas_dcopy(DOF, u_sol + (iter - 1) * (DOF), 1, sol_vect, 1);
	clean("full");

}

void solve_reduced(T* sol_vect) {
	///////////////////////////////////// For 7 DOF System //////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////// Memory Allocation for Intermediate Solution Vectors //////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////
	int sol_size = (floor(tend_ / h_) + 1);
	T factor_h2 = (1 / (h_ * h_));
	T factor_h = (1 / (h_));
	int mat_len = (DOF) * (DOF);
	//u_sol = (T*)mkl_calloc(sol_size * (DOF + 4), sizeof(T), alignment);
	u_sol_red = (T*)mkl_calloc(sol_size * (DOF), sizeof(T), alignment);

	//u_n_m_1 = (T*)mkl_calloc((DOF + 4), sizeof(T), alignment);
	//u_n = (T*)mkl_calloc((DOF + 4), sizeof(T), alignment);


	// A=((1/(h*h))*M+(1/h)*D+K);
	cblas_daxpy(mat_len, factor_h2, M_red, 1, Ared, 1);
	cblas_daxpy(mat_len, factor_h, D_red, 1, Ared, 1);
	cblas_daxpy(mat_len, 1, K_red, 1, Ared, 1);
	// Cholesky factorization of A
	lapack_int status;
	//lapack_int* piv = (lapack_int*)mkl_calloc((DOF + 4), sizeof(lapack_int), alignment);
	status = LAPACKE_dpotrf(LAPACK_ROW_MAJOR, 'L', DOF, Ared, DOF);
	check_status(status);
	// B=((2/(h*h))*M+(1/h)*D);
	cblas_daxpy(mat_len, 2 * factor_h2, M_red, 1, Bred, 1);
	cblas_daxpy(mat_len, factor_h, D_red, 1, Bred, 1);

	int iter = 1;
	T t = h_;
	/*auto start = std::chrono::steady_clock::now();*/
	while (t < tend_) {
		// u_n_p_1=A\(B*u_n-((1/(h*h))*M)*u_n_m_1+f_n_p_1);
		// u_n_p_1 = B*u_n
		cblas_dgemv(CblasRowMajor, CblasNoTrans, DOF, DOF, 1, Bred, DOF, u_n_red, 1, 1, u_n_p_1_red, 1);
		// u_n_p_1 = -((1/(h*h))*M)*u_n_m_1 + u_n_p_1
		cblas_dgemv(CblasRowMajor, CblasNoTrans, DOF, DOF, -factor_h2, M_red, DOF, u_n_m_1_red, 1, 1, u_n_p_1_red, 1);
		// u_n_p_1 = f_n_p_1 + u_n_p_1
		cblas_daxpy(DOF, 1, f_n_p_1_red, 1, u_n_p_1_red, 1);
		// u_n_p_1=A\u_n_p_1
		LAPACKE_dpotrs(LAPACK_ROW_MAJOR, 'L', DOF, 1, Ared, DOF, u_n_p_1_red, 1);
		//LAPACKE_dgetrs(LAPACK_ROW_MAJOR, 'N', DOF + 4, 1, A, DOF + 4, piv, u_n_p_1, 1);
		///////////////// Debug print statements ///////////////////////////////
		/*if (iter==2){
			write_vector(u_n_p_1_red, DOF);
		}*/
		time[iter] = t;
		// u_sol(j,:)=u_n_p_1;
		cblas_dcopy(DOF, u_n_p_1_red, 1, u_sol_red + iter * (DOF), 1);
		// u_n_m_1=u_n;
		cblas_dcopy(DOF, u_n_red, 1, u_n_m_1_red, 1);
		// u_n    =u_n_p_1;
		cblas_dcopy(DOF, u_n_p_1_red, 1, u_n_red, 1);
		iter++;
		t += h_;
	}
	/*auto end = std::chrono::steady_clock::now();
	std::cout << "Elapsed time in nanoseconds : "
		<< std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
		<< " ms" << std::endl;*/
		//std::cout << "iter = " << iter << " sol_size = "<< sol_size <<"\n\n" << std::endl;
	cblas_dcopy(DOF, u_sol_red + (iter - 1) * (DOF), 1, sol_vect, 1);
	clean("reduced");

}
