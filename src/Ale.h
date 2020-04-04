#pragma once
//#include "Car.h"
#include "Profile.h"
#include "Load_module.h"
#include "EVAAComputeEngine.h"

template <class T>
class ALE {
private:
	// define constants
	const int alignment = 64;
	
	// necessary class objects
	Car* Car_obj; // suppose Interpolation in the Car
	Load_module* Load_module_obj; // needs Profile and Car
	linear11dof<T>* linear11dof_obj;
	EVAAComputeStiffness* interpolator;// interpolator member of EVAA
	Simulation_Parameters params;
	Load_Params load_param;
	Profile* profile_obj;

	// simulation parameters
	int DOF;
	T tend_;
	T h_;

	// time and solution vectors
	T* time_vec;
	T* u_sol, u_init;


public:
	ALE(Car* Car_obj_val,
		Load_module* Load_module_val,
		linear11dof<T>* linear11dof_val,
		Simulation_Parameters &params_val) {

		Car_obj = Car_obj_val;
		Load_module_obj = Load_module_val;
		linear11dof_obj = linear11dof_val;


		// general parameters of the simulation
		params = params_val;
		DOF = params.DOF;

		h_ = params.timestep;
		tend_ = params.num_time_iter * h_;


	}

	void apply_boundary_condition(int s) {
		linear11dof_obj->apply_boundary_condition(s);
	}

	void solve(T* sol_vect) {
		///////////////////////////////////// For 11 DOF System //////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////// Memory Allocation for Intermediate Solution Vectors //////////////////
		////////////////////////////////////////////////////////////////////////////////////////////////////////////

		//initialize solution vector
		int sol_size = (floor(tend_ / h_) + 1);
		time_vec = (T*)mkl_calloc(sol_size, sizeof(T), alignment);
		u_sol = (T*)mkl_calloc(sol_size * (DOF), sizeof(T), alignment);

		// start time iteration
		T t = h_;
		double eps = h_ / 100;

		// time iteration
		while (std::abs(t - (tend_ + h_)) > eps) {
			// TODO: ALE stuff
			
			// execute one time step of the linear solver
			linear11dof_obj->solve_full_one_step(u_sol);
		}
	}
};
