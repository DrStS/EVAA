#include "profile_class.h"
#include "MathLibrary.h"


// ===============================   Circular class implementation ======================//
Circular::Circular(double* Pos, double Rad) {
	Name = "Circular";

	// Position
	Position = (double*)mkl_malloc(DIM * sizeof(double), Constants::ALIGNMENT);
	cblas_dcopy(DIM, Pos, 1, Position, 1);

	Radius = Rad;

	// auxiliary vectors
	std::cout << "Radius = " << Radius << "\n";
	velocity_direction = (double*)mkl_calloc(DIM, sizeof(double), Constants::ALIGNMENT);
	Velocity_vec = (double*)mkl_malloc(sizeof(double) * DIM * vec_DIM, Constants::ALIGNMENT);
	Mass_vec = (double*)mkl_malloc(sizeof(double) * vec_DIM, Constants::ALIGNMENT);
	dist_car_center = (double*)mkl_malloc(sizeof(double) * (DIM - 1) * vec_DIM, Constants::ALIGNMENT);
};

Circular::~Circular() {
	mkl_free(Position);
	mkl_free(velocity_direction);
	mkl_free(Velocity_vec);
	mkl_free(Mass_vec);
	mkl_free(dist_car_center);
}

void Circular::get_Position(double* pos) {
	cblas_dcopy(DIM, Position, 1, pos, 1);
}

double Circular::get_Radius() const {
	return Radius;
}

void Circular::set_Radius(const double& rad) {
	Radius = rad;
}

void Circular::set_Position(const double* pos) {
	cblas_dcopy(DIM, pos, 1, Position, 1);
}

void Circular::get_centrifugal_force_ALE(double* fr, double* v, double& m, double* p) {
	// REFACTOR TO GENERAL DIRECTIONS *

	// Adapted for the 2D now!!!

	// calculates the force in a body only with respect to its velocity, mass and Position
	// v is the velocity of the body
	// m is the mass of the body
	// p is the global Position of the body (arrow to the body)
	// fr - the centripetal force
	// the rotation is always around the origin!

	cblas_dcopy(DIM - 1, p, 1, fr, 1);

	double inv_radius = 1.0 / cblas_dnrm2(DIM - 1, p, 1);		// corresponds to the (inverse) Radius of the trajectory at the considered body

	cblas_dscal(DIM - 1, inv_radius, fr, 1); // centrifugal force

	// MathLibrary::crossProduct(fr, unit_y_vector, velocity_direction);
	// * REFACTOR TO GENERAL DIRECTIONS * this is only for z direction
	velocity_direction[0] = fr[1];
	velocity_direction[1] = -fr[0];	

	//const MKL_INT int_ddot = DIM - 1;
	//double velocity_magnitude = cblas_ddot(int_ddot, v, incx, velocity_direction, incx);
	double velocity_magnitude = cblas_dnrm2(DIM - 1, v, 1);

	double force_magnitude = m * velocity_magnitude * velocity_magnitude * inv_radius;

	cblas_dscal(DIM - 1, force_magnitude, fr, incx);
}

void Circular::get_centrifugal_force(double* fr, double* v, double& m, double* p) {
	// calculates the force in a body only with respect to its velocity, mass and Position
	// v is the velocity of the body
	// m is the mass of the body
	// p is the global Position of the body
	// fr - the centripetal force
	// the rotation is always around the origin!

	cblas_dcopy(DIM, p, 1, fr, 1);

	fr[2] = 0;		// path only in xy-plane

	double inv_radius = 1.0 / cblas_dnrm2(DIM, p, 1);		// corresponds to the (inverse) Radius of the trajectory at the considered body

	// Raffi: cblas_dscal(DIM, -inv_radius, fr, 1); - centripetal force 
	cblas_dscal(DIM, inv_radius, fr, 1); // centrifugal force

	MathLibrary::crossProduct_unitvecZ(fr, velocity_direction);

	 //double velocity_magnitude = cblas_ddot(DIM, v, incx, velocity_direction, incx);
	double velocity_magnitude = cblas_dnrm2(DIM, v, 1);

	double force_magnitude = m * velocity_magnitude * velocity_magnitude * inv_radius;

	cblas_dscal(DIM, force_magnitude, fr, incx);
}

void Circular::get_Profile_force(Car<double>* car1, double* f_vec, double* normal_ext) {
	// WRONG Implementation - suitable for 2D (BROKEN)
	// REFACTOR (MAYBE another Profile class smth)

	// out: f_vec = [f_cg, f_w1, f_t1, f_w2, f_t2, f_w3, f_t3, f_w4, f_t4] (centripetal components of the force)
	// out: normal_ext - normal over the full body; updated in this function from the centripetal forces

	// get distance vector between center of circle and the car = positions of points from the car vs center of circle, which is seen as 0
	car1->get_dist_vector_xy(Position, dist_car_center);

	//// the distance to [cg<->center of circle] must be equal with the Radius of the circle - place an assert here!!!!
	//if (abs(car1->get_dist_vector_abs_val(Position) - Radius * Radius) > 1e-12) {
	//	std::cout << "\n\n error (in Circular.update_normal_force): the Radius of the circle is different than the distance  \\
	//		between the car and center of the circle!!!!please take care!!!";
	//}

	// vectors with velocities and masses
	car1->get_Velocity_vec_xy(Velocity_vec);
	Velocity_vec[2] = 0;
	car1->get_Mass_vec(Mass_vec);

	// compute each of the 9 centripetal forces
	for (int i = 0; i < vec_DIM; ++i) {
		get_centrifugal_force(&f_vec[DIM * i], &Velocity_vec[DIM * i], Mass_vec[i], &dist_car_center[DIM * i]);
	}

	// compute centripetal part of the global normal force 
	// n = f_cg + f_w1 + f_t1 + f_w2 + f_t2 + f_w3 + f_t3 + f_w4 + f_t4
	cblas_dcopy(DIM, f_vec, incx, normal_ext, incx);
	for (auto i = 1; i < vec_DIM; ++i) {
		vdAdd(DIM, normal_ext, &f_vec[DIM * i], normal_ext);
	}

}

void Circular::get_Profile_force_ALE(Car<double>* car1, double* f_vec, double* normal_ext) {
	// out: f_vec = [f_cg, f_w1, f_t1, f_w2, f_t2, f_w3, f_t3, f_w4, f_t4] (centripetal components of the force)
	// out: normal_ext - normal over the full body; updated in this function from the centripetal forces

	// get distance vector between center of circle and the car = positions of points from the car vs center of circle, which is seen as 0
	car1->get_dist_vector_xy(Position, dist_car_center);

	//// the distance to [cg<->center of circle] must be equal with the Radius of the circle - place an assert here!!!!
	//if (abs(car1->get_dist_vector_abs_val(Position) - Radius * Radius) > 1e-12) {
	//	std::cout << "\n\n error (in Circular.update_normal_force): the Radius of the circle is different than the distance  \\
	//		between the car and center of the circle!!!!please take care!!!";
	//}

	// vectors with velocities and masses
	car1->get_Velocity_vec_xy(Velocity_vec);

	//std::cout << "Velocity vec: \n\n";
	//MathLibrary::write_vector(Velocity_vec, 18);
	car1->get_Mass_vec(Mass_vec);

	// compute each of the 9 centripetal forces
	for (int i = 0; i < vec_DIM; ++i) {
		get_centrifugal_force_ALE(&f_vec[DIM * i], &Velocity_vec[(DIM - 1) * i], Mass_vec[i], &dist_car_center[(DIM - 1) * i]);
		f_vec[DIM * i + 2] = 0; // z direction is 0 !!! Have to be generalized to general directions
	}

	// compute centripetal part of the global normal force 
	get_centrifugal_force_ALE(normal_ext, Velocity_vec, *(car1->global_mass), dist_car_center);
}

void Circular::get_Profile_torque(Car<double>* Car1, double* Torque) {
	Torque[2] = 0; // Torque on z direction
}

void Circular::update_initial_condition(Car<double>* Car1){

	std::cout << "Update initial conditions to circular motion" << std::endl;

	double* tangential_dir = (double*)mkl_malloc(Constants::DIM * sizeof(double), Constants::ALIGNMENT);
	double* radial_vector = (double*)mkl_malloc(Constants::DIM * sizeof(double), Constants::ALIGNMENT);
	radial_vector[0] = Car1->Position_vec[0] - this->Position[0];
	radial_vector[1] = Car1->Position_vec[1] - this->Position[1];
	radial_vector[2] = 0;

	double radius = cblas_dnrm2(Constants::DIM, radial_vector, 1);
	if (abs(radius - this->Radius) > 0.01)
		std::cout << "Warning! the initial position of the car is not on the trajectory provided in \
					the circular path. \n The expected radius is " 
				  << this->Radius << ", but the car is at an initial distance of " << radius 
				  << " from the center of the circle.\n The execution procedes with the current spatial "
				  << "configuration and with the current distance to the center of the circle." << std::endl;

	double inv_radius = 1. / radius;

	cblas_dscal(Constants::DIM, inv_radius, radial_vector, incx);
	MathLibrary::crossProduct_unitvecZ(radial_vector, tangential_dir);
	double magnitude = cblas_ddot(Constants::DIM, Car1->Velocity_vec, incx, tangential_dir, incx);
	cblas_dcopy(Constants::DIM, tangential_dir, 1, Car1->Velocity_vec, 1);
	cblas_dscal(Constants::DIM, magnitude, Car1->Velocity_vec, incx);
	cblas_dscal(Constants::DIM, radius, radial_vector, incx);
	MathLibrary::crossProduct(radial_vector, Car1->Velocity_vec, Car1->w_CG);
	cblas_dscal(Constants::DIM, inv_radius * inv_radius, Car1->w_CG, 1);

	MKL_free(tangential_dir);
	MKL_free(radial_vector);
}

// ============================ end of Circular class implementation ================== //



// ===============================   Fixed class implementation ======================= //
Fixed::Fixed(const double& g) {
	Name = "fixed";
	linear_idx = (size_t*)mkl_malloc(num_tyre * sizeof(size_t), Constants::ALIGNMENT);
	dx = (double*)mkl_malloc(sizeof(double) * num_tyre, Constants::ALIGNMENT);
	k_vec = (double*)mkl_malloc(sizeof(double) * num_tyre, Constants::ALIGNMENT);
	gravity = g;
};

void Fixed::get_Profile_force_ALE(Car<double>* Car1, double* F_vec, double* Normal_ext) {
	if (index_set) {
		Car1->get_k_vec_tyre(k_vec);
		Car1->compute_dx_tyre(dx);
		for (size_t i = 0; i < num_tyre; ++i) {
			F_vec[linear_idx[i] * DIM + 2] = -k_vec[i] * dx[i];
		}
	}
	else {
		std::cout << "Please Provide Tyre Index" << std::endl;
		exit(2);
	}
}

void Fixed::get_Profile_torque(Car<double>* Car1, double* Torque_vec) {
	Torque_vec[2] = 0; // Torque on z direction
}

void Fixed::set_fixed_index(size_t* index) {
	for (size_t i = 0; i < num_tyre; ++i) {
		linear_idx[i] = index[i];
	}
	index_set = 1;
}

Fixed::~Fixed() {
	mkl_free(linear_idx);
	mkl_free(dx);
	mkl_free(k_vec);
}




// =============================== Nonfixed class implementation ========================//
Nonfixed::Nonfixed(double* Pos, double Rad) {
	Name = "Nonfixed";
};

Nonfixed::~Nonfixed() {
}

void Nonfixed::get_centrifugal_force(double* fr, double* v, double& m, double* p) {
	// No external forces from a path are acting on the body
	fr[0] = 0;
	fr[1] = 0;
	fr[2] = 0;
}

void Nonfixed::get_Profile_force(Car<double>* car1, double* f_vec, double* normal_ext) {
	for (int i = 0; i < vec_DIM; ++i) {
		f_vec[DIM * i + 0] = 0;
		f_vec[DIM * i + 1] = 0;
		f_vec[DIM * i + 2] = 0;
	}

}

void Nonfixed::get_Profile_torque(Car<double>* Car1, double* Torque) {
	Torque[0] = 0; // Torque on x direction
	Torque[1] = 0; // Torque on y direction
	Torque[2] = 0; // Torque on z direction
}

void Nonfixed::update_initial_condition(Car<double>* Car1) {
	// don't change them
}

// =============================== end of Nonfixed class implementation ===================

