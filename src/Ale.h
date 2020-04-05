#pragma once
//#include "Car.h"
#include "Profile_class.h"
#include "Load_module.h"
#include "EVAAComputeEngine.h"
#include "11DOF.h"


/*
Implements the ALE method to extend the linear 11DOF system
*/
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
	Profile* profile_obj;

	// force and road parameters
	Load_Params load_param;

	// simulation parameters
	int DOF;
	T tend_;
	T h_;

	// time and solution vectors
	T* time_vec;
	T* u_sol;

	// needed to solve the 11DOF system
	T* force_vector_11dof;

	// Spring stiffnesses in Stefan's ordering
	T* k_vect;

	// Contains the forces in the ordering [GC:XYZ,W1:XYZ,T1:XYZ, ...]
	T* force_vector;

	// Contains he torque in the form [XYZ]
	T* full_torque;

	// Contains the dx in the spring elongation in Stefan's ordering
	T* Delta_x_vec;

	// global centripetal forces on the whole car [XYZ]
	T* centripetal_force;
	T* new_centripetal_force;

	// global torque on the car [XYZ]
	T* torque;
	T* new_torque;

	// global positions, velocities [XY] and angles [Z] of the center of mass of the car 
	T* posXY_vec;
	T* angleZ;
	T* velXY_vec;
	T* ang_velZ;

	// quantities for the whole car
	T global_inertia_Z;
	T global_mass;


public:
	/*
	Constructor
	*/
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

	/*
	Applies the Verlet_Stoermer algorithm to update the global XY position of the car and its Z orientation
	Store the global coordinates in the VelocityXY and PositionXY from the car object
	\param t current simulation time
	*/
	void global_frame_solver(T& t) {
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

	/*
	Executes the time iteration of the ALE solvers, switches from global position update to solving of the linear 11DOF system
	*/
	void solve(T* sol_vect) {

		//initialize solution vector
		int sol_size = (floor(tend_ / h_) + 1);
		int force_dimensions = Car_obj->DIM * Car_obj->vec_DIM;
		int centripetal_force_dimensions = 2;
		int full_torque_dimensions = 3;
		int num_springs = 2 * Car_obj->num_wheels;

		// allocate memory
		time_vec = (T*)mkl_calloc(sol_size, sizeof(T), alignment);
		u_sol = (T*)mkl_calloc(sol_size * (Car_obj->vec_DIM * Car_obj->DIM), sizeof(T), alignment);
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
			solution_vect = u_sol + iter * (Car_obj->vec_DIM * Car_obj->DIM);
			Car_obj->populate_results(Car_obj->Position_vec_xy, Car_obj->u_current_linear, solution_vect);
			
			t += h_;
			iter++;
			
		}
		cblas_dcopy(DOF, u_sol + (iter - 1) * (DOF), 1, sol_vect, 1);
		Car_obj->combine_results();


		MKL_free(time_vec);
		MKL_free(u_sol);
		MKL_free(force_vector);
		MKL_free(centripetal_force);
		MKL_free(new_centripetal_force);
		MKL_free(Delta_x_vec);

		delete[] torque;
		delete[] new_torque;
	}

	/*
	adds the contribution of the wheels and tyres to the inertia moment of the car
	*/
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

	/*
	Sums up all the 9 masses
	*/
	void calculate_global_mass() {
		global_mass = Car_obj->Mass_vec[0] + 
			Car_obj->Mass_vec[1] + Car_obj->Mass_vec[2] + Car_obj->Mass_vec[3] + Car_obj->Mass_vec[4] + 
			Car_obj->Mass_vec[5] + Car_obj->Mass_vec[6] + Car_obj->Mass_vec[7] + Car_obj->Mass_vec[8];
	}

	/*
	Prints all positions and angles in the car object
	*/
	void print_final_results() {
		T* sln = Car_obj->Position_vec;
		std::cout << "ALE: orientation angles=\n\t[" << Car_obj->angle_CG[0] << "\n\t " << Car_obj->angle_CG[1] << "\n\t " << Car_obj->angle_CG[2] << "]" << std::endl;
		std::cout << "ALE: car body position pc=\n\t[" << sln[0] << "\n\t " << sln[1] << "\n\t " << sln[2] << "]" << std::endl;
		std::cout << "ALE: rear-right wheel position pw1=\n\t[" << sln[21] << "\n\t " << sln[22] << "\n\t " << sln[23] << "]" << std::endl;
		std::cout << "ALE: rear-left wheel position pw2=\n\t[" << sln[15] << "\n\t " << sln[16] << "\n\t " << sln[17] << "]" << std::endl;
		std::cout << "ALE: front-left wheel position pw3=\n\t[" << sln[3] << "\n\t " << sln[4] << "\n\t " << sln[5] << "]" << std::endl;
		std::cout << "ALE: front-right wheel position pw4=\n\t[" << sln[9] << "\n\t " << sln[10] << "\n\t " << sln[11] << "]" << std::endl;
		std::cout << "ALE: rear-right tyre position pt1=\n\t[" << sln[24] << "\n\t " << sln[25] << "\n\t " << sln[26] << "]" << std::endl;
		std::cout << "ALE: rear-left tyre position pt2=\n\t[" << sln[18] << "\n\t " << sln[19] << "\n\t " << sln[20] << "]" << std::endl;
		std::cout << "ALE: front-left tyre position pt3=\n\t[" << sln[6] << "\n\t " << sln[7] << "\n\t " << sln[8] << "]" << std::endl;
		std::cout << "ALE: front-right tyre position pt4=\n\t[" << sln[12] << "\n\t " << sln[13] << "\n\t " << sln[14] << "]" << std::endl;

	}
};
