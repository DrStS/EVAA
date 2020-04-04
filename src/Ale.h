#pragma once
//#include "Car.h"
#include "Profile.h"
#include "Load_module.h"
#include "EVAAComputeEngine.h"
#include "11DOF.h"

template <class T>
class ALE {
private:
	// define constants
	const int alignment = 64;
	
	// necessary class objects
	Car<T>* Car_obj; // suppose Interpolation in the Car
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
	T* force_vector;
	T* weighted_forceXY;
	T* torque;
	T* posXY_vec;
	T* angleZ;
	T* velXY_vec;
	T* ang_velZ;

public:
	ALE(Car<T>* Car_obj_val,
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

	void global_frame_solver() {
		/* ??????????????????????????????????????? */

		// Compute weighted force sum TODO: Shubham
		function_from_shubham(weighted_forceXY);

		// 1. Update global X,Y velocities
		Car_obj->Velocity_vec_xy;

		// 2. Update global X,Y positions of the car
		Car_obj->Position_vect_xy;

		// 3. Update Z-angular velocities
		Car_obj->w_z;

		// 4. Update Z-rotation
		Car_obj->Angle_z;

		// Implement ALE solver!!!!!!

	}

	void solve(T* sol_vect) {

		//initialize solution vector
		int sol_size = (floor(tend_ / h_) + 1);
		time_vec = (T*)mkl_calloc(sol_size, sizeof(T), alignment);
		u_sol = (T*)mkl_calloc(sol_size * (DOF), sizeof(T), alignment);
		int force_dimensions = Car_obj->mkl_DIM * Car_obj->vec_DIM;
		int weighted_force_dimensions = 3;
		force_vector = (T*)mkl_calloc(force_dimensions, sizeof(T), alignment);
		weighted_forceXY = (T*)mkl_calloc(weighted_force_dimensions, sizeof(T), alignment);
		torque = new(T);

		// start time iteration
		T t = h_;
		double eps = h_ / 100;

		linear11dof_obj->solution_initialize(u_sol); // TODO: Shubham

		// time iteration
		while (std::abs(t - (tend_ + h_)) > eps) {

			Load_module_obj->update_force(t, force_vector, Delta_x_vec, External_force); // TODO: ask Teo
			Load_module_obj->update_torque(t, torque, Delta_x_vec, External_force); // TODO: ask Teo

			global_frame_solver();

			// execute one time step of the linear solver
			linear11dof_obj->solve_full_one_step(u_sol);
		}

		MKL_free(time_vec);
		MKL_free(u_sol);
		MKL_free(force_vector);
		MKL_free(weighted_force_vectorXY);
	}
};
