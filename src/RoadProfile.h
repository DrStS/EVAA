#pragma once

#include "BLAS.h"
#include "Car.h"
#include "Constants.h"
#include "MathLibrary.h"
#include <string>

/*
Parent class for generic road profiles
*/
template <typename T>
class Profile {
protected:
	std::string Name;

public:
	/*
	Get external force acting on the car system
	\param Car
	\return F_vec forces acting on each component [GC:XYZ,W1:XYZ,T1:XYZ, ...]
	\return Normal_ext on the car as a whole [XYZ]
	*/
	virtual void get_Profile_force(Car<T>* Car1, T* F_vec, T* Normal_ext) {};
	virtual void get_Profile_force_ALE(Car<T>* Car1, T* F_vec, T* Normal_ext) {};

	/*
	Get external torque acting on the car system
	\param Car
	\return F_vec torque acting on teh car system [XYZ]
	*/
	virtual void get_Profile_torque(Car<T>* Car1, T* Torque_vec) {};

	/*
	In case the initial conditions have to be specific to a road profile
	\param Car
	*/
	virtual void update_initial_condition(Car<T>* Car1) { std::cout << "No update of initial conditions" << std::endl; };
	virtual void set_fixed_index(size_t* index) {};

	virtual ~Profile() {};
};

// ===============================   Circular class =================================//
/*
Follow a circular road profile with radius R and center C
*/
template <typename T>
class Circular : public Profile<T> {
private:
	// Center of the Circle (if the motion is a Circle)
	T* Position = NULL; 
	
	// Radius of the Circle (if the motion is a Circle)
	T Radius = 0; 

	// vectors used in computations inside the methods
	T* unit_y_vector;
	T* velocity_direction; // arrow from center towards object

	// distance vector between center of circle and the car = 
	//		Positions of points from the car vs Center of Circle, seen as 0
	T* dist_car_center; // 3 x 9
	T* Velocity_vec;
	T* Mass_vec;

public:
	Circular(T* Pos, T Rad) {
		Name = "Circular";

		// Position
		Position = (T*)mkl_malloc(Constants::DIM * sizeof(T), Constants::ALIGNMENT);
		mkl<T>::copy(Constants::DIM, Pos, 1, Position, 1);

		Radius = Rad;

		// auxiliary vectors
		velocity_direction = (T*)mkl_calloc(Constants::DIM, sizeof(T), Constants::ALIGNMENT);
		Velocity_vec = (T*)mkl_malloc(sizeof(T) * Constants::DIM * Constants::VEC_DIM, Constants::ALIGNMENT);
		Mass_vec = (T*)mkl_malloc(sizeof(T) * Constants::VEC_DIM, Constants::ALIGNMENT);
		dist_car_center = (T*)mkl_malloc(sizeof(T) * (Constants::DIM - 1) * Constants::VEC_DIM, Constants::ALIGNMENT);
	}

	virtual ~Circular() {
		mkl_free(Position);
		mkl_free(velocity_direction);
		mkl_free(Velocity_vec);
		mkl_free(Mass_vec);
		mkl_free(dist_car_center);
	}

	void get_Position(T* Pos) {
		mkl<T>::copy(Constants::DIM, Position, 1, pos, 1);
	}

	inline T get_Radius() const {
		return Radius;
	}

	inline void set_Radius(const T& Rad) {
		Radius = rad;
	}

	void set_Position(const T* Pos) {
		mkl<T>::copy(Constants::DIM, pos, 1, Position, 1);
	}

	void get_centrifugal_force(T* Fr, T* v, T& m, T* p) {
		// calculates the force in a body only with respect to its velocity, mass and Position
		// v is the velocity of the body
		// m is the mass of the body
		// p is the global Position of the body
		// fr - the centripetal force
		// the rotation is always around the origin!

		mkl<T>::copy(Constants::DIM, p, 1, Fr, 1);

		Fr[2] = 0;		// path only in xy-plane

		T inv_radius = 1.0 / mkl<T>::nrm2(Constants::DIM, p, 1);		// corresponds to the (inverse) Radius of the trajectory at the considered body

		// Raffi: mkl<T>::scal(Constants::DIM, -inv_radius, fr, 1); - centripetal force 
		mkl<T>::scal(Constants::DIM, inv_radius, Fr, 1); // centrifugal force

		MathLibrary::crossProduct_unitvecZ(Fr, velocity_direction);

		//T velocity_magnitude = mkl<T>::dot(Constants::DIM, v, 1, velocity_direction, 1);
		T velocity_magnitude = mkl<T>::nrm2(Constants::DIM, v, 1);

		T force_magnitude = m * velocity_magnitude * velocity_magnitude * inv_radius;

		mkl<T>::scal(Constants::DIM, force_magnitude, Fr, 1);
	}

	void get_centrifugal_force_ALE(T* Fr, T* v, T& m, T* p) {
		// REFACTOR TO GENERAL DIRECTIONS *

		// Adapted for the 2D now!!!

		// calculates the force in a body only with respect to its velocity, mass and Position
		// v is the velocity of the body
		// m is the mass of the body
		// p is the global Position of the body (arrow to the body)
		// fr - the centripetal force
		// the rotation is always around the origin!

		mkl<T>::copy(Constants::DIM - 1, p, 1, Fr, 1);

		T inv_radius = 1.0 / mkl<T>::nrm2(Constants::DIM - 1, p, 1);		// corresponds to the (inverse) Radius of the trajectory at the considered body

		mkl<T>::scal(Constants::DIM - 1, inv_radius, Fr, 1); // centrifugal force

		// MathLibrary::crossProduct(Fr, unit_y_vector, velocity_direction);
		// * REFACTOR TO GENERAL DIRECTIONS * this is only for z direction
		velocity_direction[0] = Fr[1];
		velocity_direction[1] = -Fr[0];

		//const MKL_INT int_ddot = DIM - 1;
		//T velocity_magnitude = mkl<T>::dot(int_ddot, v, 1, velocity_direction, 1);
		T velocity_magnitude = mkl<T>::nrm2(Constants::DIM - 1, v, 1);

		T force_magnitude = m * velocity_magnitude * velocity_magnitude * inv_radius;

		mkl<T>::scal(Constants::DIM - 1, force_magnitude, Fr, 1);
	}


	/*
	Get external force acting on the car system
	\param Car
	\return F_vec forces acting on each component [GC:XYZ,W1:XYZ,T1:XYZ, ...] centripetals on each component
	\return Normal_ext on the car as a whole [XYZ]
	*/
	virtual void get_Profile_force(Car<T>* Car1, T* F_vec, T* Normal_ext) {
		// WRONG =========== - suitable for 2D (BROKEN)
		// REFACTOR (MAYBE another Profile class smth)

		// out: F_vec = [f_cg, f_w1, f_t1, f_w2, f_t2, f_w3, f_t3, f_w4, f_t4] (centripetal components of the force)
		// out: normal_ext - normal over the full body; updated in this function from the centripetal forces

		// get distance vector between center of circle and the car = positions of points from the car vs center of circle, which is seen as 0
		Car1->get_dist_vector_xy(Position, dist_car_center);

		//// the distance to [cg<->center of circle] must be equal with the Radius of the circle - place an assert here!!!!
		//if (abs(Car1->get_dist_vector_abs_val(Position) - Radius * Radius) > 1e-12) {
		//	std::cout << "\n\n error (in Circular.update_normal_force): the Radius of the circle is different than the distance  \\
		//		between the car and center of the circle!!!!please take care!!!";
		//}

		// vectors with velocities and masses
		Car1->get_Velocity_vec_xy(Velocity_vec);
		Velocity_vec[2] = 0;
		Car1->get_Mass_vec(Mass_vec);

		// compute each of the 9 centripetal forces
		for (int i = 0; i < Constants::VEC_DIM; ++i) {
			get_centrifugal_force(&F_vec[Constants::DIM * i], &Velocity_vec[Constants::DIM * i], Mass_vec[i], &dist_car_center[Constants::DIM * i]);
		}

		// compute centripetal part of the global normal force 
		// n = f_cg + f_w1 + f_t1 + f_w2 + f_t2 + f_w3 + f_t3 + f_w4 + f_t4
		mkl<T>::copy(Constants::DIM, F_vec, 1, Normal_ext, 1);
		for (auto i = 1; i < Constants::VEC_DIM; ++i) {
			mkl<T>::vAdd(Constants::DIM, Normal_ext, &F_vec[Constants::DIM * i], Normal_ext);
		}

	}

	virtual void get_Profile_force_ALE(Car<T>* Car1, T* F_vec, T* Normal_ext) {
		// out: F_vec = [f_cg, f_w1, f_t1, f_w2, f_t2, f_w3, f_t3, f_w4, f_t4] (centripetal components of the force)
		// out: normal_ext - normal over the full body; updated in this function from the centripetal forces

		// get distance vector between center of circle and the car = positions of points from the car vs center of circle, which is seen as 0
		Car1->get_dist_vector_xy(Position, dist_car_center);

		//// the distance to [cg<->center of circle] must be equal with the Radius of the circle - place an assert here!!!!
		//if (abs(Car1->get_dist_vector_abs_val(Position) - Radius * Radius) > 1e-12) {
		//	std::cout << "\n\n error (in Circular.update_normal_force): the Radius of the circle is different than the distance  \\
		//		between the car and center of the circle!!!!please take care!!!";
		//}

		// vectors with velocities and masses
		Car1->get_Velocity_vec_xy(Velocity_vec);

		//std::cout << "Velocity vec: \n\n";
		//MathLibrary::write_vector(Velocity_vec, 18);
		Car1->get_Mass_vec(Mass_vec);

		// compute each of the 9 centripetal forces
		for (int i = 0; i < Constants::VEC_DIM; ++i) {
			get_centrifugal_force_ALE(&F_vec[Constants::DIM * i], &Velocity_vec[(Constants::DIM - 1) * i], Mass_vec[i], &dist_car_center[(Constants::DIM - 1) * i]);
			F_vec[Constants::DIM * i + 2] = 0; // z direction is 0 !!! Have to be generalized to general directions
		}

		// compute centripetal part of the global normal force 
		get_centrifugal_force_ALE(Normal_ext, Velocity_vec, *(Car1->massFullCar), dist_car_center);
	}


	/*
	Get external torque acting on the car system
	\param Car
	\return F_vec torque acting on teh car system [XYZ]
	*/
	virtual void get_Profile_torque(Car<T>* Car1, T* Torque) {
		Torque[0] = 0;
		Torque[1] = 0;
		Torque[2] = 0; // Torque on z direction
	}

	/*
	overwrites all initial velocities but the one from the main car body
	Calculates the initial angular velocity such that the car perfectly rotates around its own axis as it follows the circle
	\param Car
	*/
	virtual void update_initial_condition(Car<T>* Car1) {

		std::cout << "Update initial conditions to circular motion" << std::endl;

		T* tangential_dir = (T*)mkl_malloc(Constants::DIM * sizeof(T), Constants::ALIGNMENT);
		T* radial_vector = (T*)mkl_malloc(Constants::DIM * sizeof(T), Constants::ALIGNMENT);
		radial_vector[0] = Car1->Position_vec[0] - this->Position[0];
		radial_vector[1] = Car1->Position_vec[1] - this->Position[1];
		radial_vector[2] = 0;

		T radius = mkl<T>::nrm2(Constants::DIM, radial_vector, 1);
		if (abs(radius - this->Radius) > 0.01)
			std::cout << "Warning! the initial position of the car is not on the trajectory provided in \
					the circular path. \n The expected radius is "
			<< this->Radius << ", but the car is at an initial distance of " << radius
			<< " from the center of the circle.\n The execution procedes with the current spatial "
			<< "configuration and with the current distance to the center of the circle." << std::endl;

		T inv_radius = 1. / radius;

		mkl<T>::scal(Constants::DIM, inv_radius, radial_vector, 1);
		MathLibrary::crossProduct_unitvecZ(radial_vector, tangential_dir);
		T magnitude = mkl<T>::dot(Constants::DIM, Car1->Velocity_vec, 1, tangential_dir, 1);
		mkl<T>::copy(Constants::DIM, tangential_dir, 1, Car1->Velocity_vec, 1);
		mkl<T>::scal(Constants::DIM, magnitude, Car1->Velocity_vec, 1);
		mkl<T>::scal(Constants::DIM, radius, radial_vector, 1);
		MathLibrary::crossProduct(radial_vector, Car1->Velocity_vec, Car1->w_CG);
		mkl<T>::scal(Constants::DIM, inv_radius * inv_radius, Car1->w_CG, 1);

		MKL_free(tangential_dir);
		MKL_free(radial_vector);
	}

};
// ============================ end of Circular class ============================= //


// ===============================   Fixed class ================================== //
template <typename T>
class Fixed : public Profile<T> {
public:
	Fixed(const T& g) {
		Name = "fixed";
		linear_idx = (size_t*)mkl_malloc(num_tyre * sizeof(size_t), Constants::ALIGNMENT);
		dx = (T*)mkl_malloc(sizeof(T) * num_tyre, Constants::ALIGNMENT);
		k_vec = (T*)mkl_malloc(sizeof(T) * num_tyre, Constants::ALIGNMENT);
		gravity = g;
	};
	virtual void get_Profile_force_ALE(Car<T>* Car1, T* F_vec, T* Normal_ext) {
		if (index_set) {
			Car1->compute_dx_tyre(dx);
			for (size_t i = 0; i < num_tyre; ++i) {
				F_vec[linear_idx[i] * Constants::DIM + 2] = -k_vec[i] * dx[i];
			}
		}
		else {
			std::cout << "Please Provide Tyre Index" << std::endl;
			exit(2);
		}
	}

	virtual void get_Profile_torque(Car<T>* Car1, T* Torque_vec) {
		Torque_vec[0] = 0;
		Torque_vec[1] = 0;
		Torque_vec[2] = 0; // Torque on z direction
	}

	virtual void set_fixed_index(size_t* index) {
		for (size_t i = 0; i < num_tyre; ++i) {
			linear_idx[i] = index[i];
		}
		index_set = 1;
	}

	~Fixed() {
		mkl_free(linear_idx);
		mkl_free(dx);
		mkl_free(k_vec);
	}
private:
	size_t* linear_idx;
	T *k_vec, *dx;
	T* external_force;
	size_t num_tyre = 4;
	bool index_set = 0;
	T gravity;
};


// =============================== Nonfixed class ===================================//
/*
Follow a circular road profile with radius R and center C
*/
template <typename T>
class Nonfixed : public Profile<T> {
private:

public:
	Nonfixed(T* Pos, T Rad) {
		Name = "Nonfixed";
	}

	virtual ~Nonfixed() {};

	void get_centrifugal_force(T* Fr, T* v, T& m, T* p) {
		// No external forces from a path are acting on the body
		Fr[0] = 0;
		Fr[1] = 0;
		Fr[2] = 0;
	}

	/*
	Get external force acting on the car system
	\param Car
	\return F_vec forces acting on each component [GC:XYZ,W1:XYZ,T1:XYZ, ...]
	\return Normal_ext on the car as a whole [XYZ]
	*/
	virtual void get_Profile_force(Car<T>* Car1, T* F_vec, T* Normal_ext) {
		for (int i = 0; i < Constants::VEC_DIM; ++i) {
			F_vec[Constants::DIM * i + 0] = 0;
			F_vec[Constants::DIM * i + 1] = 0;
			F_vec[Constants::DIM * i + 2] = 0;
		}

	}


	/*
	Get external torque acting on the car system
	\param Car
	\return F_vec torque acting on teh car system [XYZ]
	*/
	virtual void get_Profile_torque(Car<T>* Car1, T* Torque) {
		Torque[0] = 0; // Torque on x direction
		Torque[1] = 0; // Torque on y direction
		Torque[2] = 0; // Torque on z direction
	}

	/*
	overwrites all initial velocities but the one from the main car body
	Calculates the initial angular velocity such that the car perfectly rotates around its own axis as it follows the circle
	\param Car
	*/
	virtual void update_initial_condition(Car<T>* Car1) {
		// don't change them
	};
};

// =============================== end of Nonfixed class ==============================