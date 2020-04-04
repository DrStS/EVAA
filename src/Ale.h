#pragma once
//#include "Car.h"
#include "Profile_class.h"
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
	T* new_force_vector;
	T* weighted_forceXY;
	T* new_weighted_forceXY;
	T* torque;
	T* new_torque;
	T* posXY_vec;
	T* angleZ;
	T* velXY_vec;
	T* ang_velZ;
	T global_inertial_Z;
	T global_mass;


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

		// 2. Update global X,Y positions of the car
		MathLibrary::Solvers<T, ALE>::Stoermer_Verlet_Position(Car_obj->Position_vect_xy[0], Car_obj->Velocity_vec_xy[0], weighted_forceXY[0], h_, global_mass);			;
		MathLibrary::Solvers<T, ALE>::Stoermer_Verlet_Position(Car_obj->Position_vect_xy[1], Car_obj->Velocity_vec_xy[1], weighted_forceXY[1], h_, global_mass); 

		// 4. Update Z-rotation
		MathLibrary::Solvers<T, ALE>::Stoermer_Verlet_Position(Car_obj->Angle_z, Car_obj->w_z, torque, h_, global_inertia_Z);

		// get forces 
		Load_module_obj->update_force(t, new_force_vector, Delta_x_vec); // TODO: ask Teo
		Load_module_obj->update_torque(t, new_torque, Delta_x_vec); // TODO: ask Teo

		// Compute weighted force sum TODO: Shubham
		function_from_shubham(new_weighted_forceXY);

		// 1. Update global X,Y velocities
		MathLibrary::Solvers<T, ALE>::Stoermer_Verlet_Velocity(Car_obj->Velocity_vec_xy[0], weighted_forceXY[0], new_weighted_forceXY[0], h_, global_mass);
		MathLibrary::Solvers<T, ALE>::Stoermer_Verlet_Velocity(Car_obj->Velocity_vec_xy[1], weighted_forceXY[1], new_weighted_forceXY[1], h_, global_mass);

		// 3. Update Z-angular velocities
		MathLibrary::Solvers<T, ALE>::Stoermer_Verlet_Velocity(Car_obj->w_z, torque, new_torque, h_, global_inertial_Z);


		// Implement ALE solver!!!!!!


		// update forces and torque
		weighted_forceXY[0] = new_weighted_forceXY[0];
		weighted_forceXY[1] = new_weighted_forceXY[1];

		torque[0] = new_torque[0];

	}

	void solve(T* sol_vect) {

		//initialize solution vector
		int sol_size = (floor(tend_ / h_) + 1);
		time_vec = (T*)mkl_calloc(sol_size, sizeof(T), alignment);
		u_sol = (T*)mkl_calloc(sol_size * (DOF), sizeof(T), alignment);
		int force_dimensions = Car_obj->DIM * Car_obj->vec_DIM;
		int weighted_force_dimensions = 3;
		force_vector = (T*)mkl_calloc(force_dimensions, sizeof(T), alignment);
		new_force_vector = (T*)mkl_calloc(force_dimensions, sizeof(T), alignment);
		weighted_forceXY = (T*)mkl_calloc(weighted_force_dimensions, sizeof(T), alignment);
		new_weighted_forceXY = (T*)mkl_calloc(weighted_force_dimensions, sizeof(T), alignment);
		torque = new(T);
		new_torque = new(T);

		calculate_global_inertia_Z();
		calculate_global_mass();

		Load_module_obj->update_force(t, force_vector, Delta_x_vec); // TODO: ask Teo
		Load_module_obj->update_torque(t, torque, Delta_x_vec); // TODO: ask Teo
		// Compute weighted force sum TODO: Shubham
		function_from_shubham(weighted_forceXY);


		// start time iteration
		T t = h_;
		double eps = h_ / 100;

		linear11dof_obj->solution_initialize(u_sol); // TODO: Shubham

		// time iteration
		while (std::abs(t - (tend_ + h_)) > eps) {

			global_frame_solver();

			// execute one time step of the linear solver
			linear11dof_obj->solve_full_one_step(u_sol);
		}

		MKL_free(time_vec);
		MKL_free(u_sol);
		MKL_free(force_vector);
		MKL_free(new_force_vector);
		MKL_free(weighted_forceXY);
		MKL_free(new_weighted_forceXY);

		delete torque;
		delete new_torque;
	}

	void calculate_global_inertia_Z() {
		// get the global inertia actiing in Z direction
		global_inertial_Z = Car_obj->I_CG[8];
		global_inertial_Z += (Car_obj->Mass_vec[1] + Car_obj->Mass_vec[2]) * 
			(Car_obj->l_lat[0] * Car_obj->l_lat[0] + Car_obj->l_long[0] * Car_obj->l_long[0]);
		global_inertial_Z += (Car_obj->Mass_vec[3] + Car_obj->Mass_vec[4]) *
			(Car_obj->l_lat[1] * Car_obj->l_lat[1] + Car_obj->l_long[1] * Car_obj->l_long[1]);
		global_inertial_Z += (Car_obj->Mass_vec[5] + Car_obj->Mass_vec[6]) *
			(Car_obj->l_lat[2] * Car_obj->l_lat[2] + Car_obj->l_long[2] * Car_obj->l_long[2]);
		global_inertial_Z += (Car_obj->Mass_vec[7] + Car_obj->Mass_vec[8]) *
			(Car_obj->l_lat[3] * Car_obj->l_lat[3] + Car_obj->l_long[3] * Car_obj->l_long[3]);
	}

	void calculate_global_mass() {
		global_mass = Car_obj->Mass_vec[0] + 
			Car_obj->Mass_vec[1] + Car_obj->Mass_vec[2] + Car_obj->Mass_vec[3] + Car_obj->Mass_vec[4] + 
			Car_obj->Mass_vec[5] + Car_obj->Mass_vec[6] + Car_obj->Mass_vec[7] + Car_obj->Mass_vec[8];
	}

};
