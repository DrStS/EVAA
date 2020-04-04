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

	// needed to solve the 11DOF system
	T* force_vector_11dof;
	T* k_vect;

	// needed to apply the load module
	T* force_vector;
	T* full_torque;
	T* Delta_x_vec;

	// needed for position integrator
	T* centripetal_force;
	T* new_centripetal_force;
	T* torque;
	T* new_torque;
	T* posXY_vec;
	T* angleZ;
	T* velXY_vec;
	T* ang_velZ;
	T global_inertia_Z;
	T global_mass;


public:
	ALE(Car<T>* Car_obj_val,
		Load_module* Load_module_val,
		linear11dof<T>* linear11dof_val,
		EVAAComputeStiffness* lookup_table,
		Simulation_Parameters &params_val) {

		Car_obj = Car_obj_val;
		Load_module_obj = Load_module_val;
		linear11dof_obj = linear11dof_val;
		interpolator = lookup_table;

		// general parameters of the simulation
		params = params_val;
		DOF = params.DOF;

		h_ = params.timestep;
		tend_ = params.num_time_iter * h_;
	}

	void global_frame_solver(T& t) {
		/* ??????????????????????????????????????? */

		// 2. Update global X,Y positions of the car
		MathLibrary::Solvers<T, ALE>::Stoermer_Verlet_Position(Car_obj->Position_vec_xy[0], Car_obj->Velocity_vec_xy[0], centripetal_force[0], h_, global_mass);			;
		MathLibrary::Solvers<T, ALE>::Stoermer_Verlet_Position(Car_obj->Position_vec_xy[1], Car_obj->Velocity_vec_xy[1], centripetal_force[1], h_, global_mass); 
		#pragma loop(ivdep)
		for (size_t i = 1; i < Car_obj->vec_DIM; ++i) {
			Car_obj->Position_vec_xy[2 * i] = Car_obj->Position_vec_xy[0];
		}
		#pragma loop(ivdep)
		for (size_t i = 1; i < Car_obj->vec_DIM; ++i) {
			Car_obj->Position_vec_xy[2 * i + 1] = Car_obj->Position_vec_xy[1];
		}

		// 4. Update Z-rotation
		MathLibrary::Solvers<T, ALE>::Stoermer_Verlet_Position(*Car_obj->Angle_z, *Car_obj->w_z, torque[2], h_, global_inertia_Z);

		// get forces 
		Car_obj->compute_dx(Delta_x_vec);
		Load_module_obj->update_force(t, force_vector, Delta_x_vec, new_centripetal_force);
		Load_module_obj->update_torque(t, new_torque, Delta_x_vec);


		// 1. Update global X,Y velocities
		MathLibrary::Solvers<T, ALE>::Stoermer_Verlet_Velocity(Car_obj->Velocity_vec_xy[0], centripetal_force[0], new_centripetal_force[0], h_, global_mass);
		MathLibrary::Solvers<T, ALE>::Stoermer_Verlet_Velocity(Car_obj->Velocity_vec_xy[1], centripetal_force[1], new_centripetal_force[1], h_, global_mass);
		#pragma loop(ivdep)
		for (size_t i = 1; i < Car_obj->vec_DIM; ++i) {
			Car_obj->Velocity_vec_xy[2 * i] = Car_obj->Velocity_vec_xy[0];
		}
		#pragma loop(ivdep)
		for (size_t i = 1; i < Car_obj->vec_DIM; ++i) {
			Car_obj->Velocity_vec_xy[2 * i + 1] = Car_obj->Velocity_vec_xy[1];
		}

		// 3. Update Z-angular velocities
		MathLibrary::Solvers<T, ALE>::Stoermer_Verlet_Velocity(*Car_obj->w_z, torque[2], new_torque[2], h_, global_inertia_Z);

		// update forces and torque
		centripetal_force[0] = new_centripetal_force[0];
		centripetal_force[1] = new_centripetal_force[1];

		torque[2] = new_torque[2]; // z - component

	}

	void solve(T* sol_vect) {

		//initialize solution vector
		int sol_size = (floor(tend_ / h_) + 1);
		int force_dimensions = Car_obj->DIM * Car_obj->vec_DIM;
		int centripetal_force_dimensions = 2;
		int full_torque_dimensions = 3;
		int num_springs = 2 * Car_obj->num_wheels;

		// allocate memory
		time_vec = (T*)mkl_calloc(sol_size, sizeof(T), alignment);
		u_sol = (T*)mkl_calloc(sol_size * (Car_obj->vec_DIM* Car_obj->DIM), sizeof(T), alignment);
		force_vector = (T*)mkl_calloc(force_dimensions, sizeof(T), alignment);
		force_vector_11dof = (T*)mkl_calloc(DOF, sizeof(T), alignment);
		full_torque = (T*)mkl_calloc(full_torque_dimensions, sizeof(T), alignment);
		centripetal_force = (T*)mkl_calloc(centripetal_force_dimensions, sizeof(T), alignment);
		new_centripetal_force = (T*)mkl_calloc(centripetal_force_dimensions, sizeof(T), alignment);
		k_vect = (T*)mkl_calloc(num_springs, sizeof(T), alignment);
		Delta_x_vec = (T*)mkl_calloc(num_springs, sizeof(T), alignment);

		torque = new T[3];
		new_torque = new T[3];

		// calculate characteristics of the whole car
		calculate_global_inertia_Z();
		calculate_global_mass();

		// start time iteration
		T t = h_;
		Car_obj->compute_dx(Delta_x_vec);
		Load_module_obj->update_force(t, force_vector, Delta_x_vec, centripetal_force);
		Load_module_obj->update_torque(t, torque, Delta_x_vec);

		// initialize the linear solver
		linear11dof_obj->initialize_solver(h_);
		T* solution_vect;
		int iter = 1;
		// time iteration
		double eps = h_ / 100;
		while (std::abs(t - (tend_ + h_)) > eps) {

			global_frame_solver(t);

			// translate 27 force vector + 3 torques into 11DOF
			Car_obj->construct_11DOF_vector(force_vector, new_torque, force_vector_11dof);
			
			linear11dof_obj->update_step(force_vector_11dof, Car_obj->u_current_linear);

			if (params.interpolation) {
				Car_obj->update_lengths_11DOF();
				interpolator->getStiffness(Car_obj->current_spring_length, k_vect);
				Car_obj->update_K(k_vect);
			}
			solution_vect = u_sol + iter * (Car_obj->vec_DIM*Car_obj->DIM);
			Car_obj->populate_results(Car_obj->Position_vec_xy, Car_obj->u_current_linear, solution_vect);
			MathLibrary::write_vector(solution_vect, (Car_obj->vec_DIM*Car_obj->DIM));
			exit(5);
			t += h_;
			iter++;
			
		}
		cblas_dcopy(DOF, u_sol + (iter - 1)*(DOF), 1, sol_vect, 1);



		MKL_free(time_vec);
		MKL_free(u_sol);
		MKL_free(force_vector);
		MKL_free(centripetal_force);
		MKL_free(new_centripetal_force);
		MKL_free(Delta_x_vec);

		delete[] torque;
		delete[] new_torque;
	}

	void calculate_global_inertia_Z() {
		// get the global inertia actiing in Z direction
		global_inertia_Z = Car_obj->I_CG[8];
		global_inertia_Z += (Car_obj->Mass_vec[1] + Car_obj->Mass_vec[2]) * 
			(Car_obj->l_lat[0] * Car_obj->l_lat[0] + Car_obj->l_long[0] * Car_obj->l_long[0]);
		global_inertia_Z += (Car_obj->Mass_vec[3] + Car_obj->Mass_vec[4]) *
			(Car_obj->l_lat[1] * Car_obj->l_lat[1] + Car_obj->l_long[1] * Car_obj->l_long[1]);
		global_inertia_Z += (Car_obj->Mass_vec[5] + Car_obj->Mass_vec[6]) *
			(Car_obj->l_lat[2] * Car_obj->l_lat[2] + Car_obj->l_long[2] * Car_obj->l_long[2]);
		global_inertia_Z += (Car_obj->Mass_vec[7] + Car_obj->Mass_vec[8]) *
			(Car_obj->l_lat[3] * Car_obj->l_lat[3] + Car_obj->l_long[3] * Car_obj->l_long[3]);
	}

	void calculate_global_mass() {
		global_mass = Car_obj->Mass_vec[0] + 
			Car_obj->Mass_vec[1] + Car_obj->Mass_vec[2] + Car_obj->Mass_vec[3] + Car_obj->Mass_vec[4] + 
			Car_obj->Mass_vec[5] + Car_obj->Mass_vec[6] + Car_obj->Mass_vec[7] + Car_obj->Mass_vec[8];
	}

};
