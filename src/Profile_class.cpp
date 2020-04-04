#include "profile_class.h"
#include "MathLibrary.h"

// ===============================   Circular class implementation ===================
Circular::Circular() {
	Name = "Circular";
	Position = (double*)mkl_malloc(sizeof(double) * 3, alignment);
	Position[2] = Position[1] = Position[0] = 0;
	Radius = 100;
	
	// auxiliary vectors
	unit_y_vector = (double*)mkl_calloc(DIM, sizeof(double), alignment);
	unit_y_vector[0] = 0.; unit_y_vector[1] = 0.; unit_y_vector[2] = 1.; // in z direction is 0

	velocity_direction = (double*)mkl_calloc(DIM, sizeof(double), alignment);
	dist_car_center = (double*)mkl_malloc(sizeof(double) * DIM * vec_DIM, alignment);
	Velocity_vec = (double*)mkl_malloc(sizeof(double) * DIM * vec_DIM, alignment);
	Mass_vec = (double*)mkl_malloc(sizeof(double) * vec_DIM, alignment);
}

Circular::Circular(double* pos) : Circular::Circular() {
	// Position
	if (pos != NULL) {
		cblas_dcopy(3, pos, 1, Position, 1);
	}
}

Circular::Circular(double* pos, double rad) : Circular::Circular(pos) {
	Radius = rad;
};

Circular::~Circular() {
	std::cout << "I am getting fucked in destruction of profile" << std::endl;
	mkl_free(Position);
	mkl_free(unit_y_vector);
	mkl_free(velocity_direction);
	mkl_free(Velocity_vec);
	mkl_free(Mass_vec);
	mkl_free(dist_car_center);
}

void Circular::get_Position(double* pos) {
	if (pos != NULL) {
		cblas_dcopy(3, Position, 1, pos, 1);
	}
}

double Circular::get_Radius() const {
	return Radius;
}

void Circular::set_Radius(const double& rad) {
	Radius = rad;
}

void Circular::set_Position(const double* pos) {
	if (pos != NULL) {
		cblas_dcopy(DIM, pos, 1, Position, 1);
	}
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

	double inv_radius = 1.0 / cblas_dnrm2(DIM, p, 1);		// corresponds to the (inverse) Radius of the trajectory at the considered tyre

	// Raffi: cblas_dscal(DIM, -inv_radius, fr, 1); - centripetal force 
	cblas_dscal(DIM, inv_radius, fr, 1); // centrifugal force

	MathLibrary::crossProduct(fr, unit_y_vector, velocity_direction);

	double velocity_magnitude = cblas_ddot(DIM, v, incx, velocity_direction, incx);

	double force_magnitude = m * velocity_magnitude * velocity_magnitude * inv_radius;

	cblas_dscal(DIM, force_magnitude, fr, incx);
}

void Circular::get_Profile_force(Car<double>* car1, double* f_vec, double* normal_ext) {
	// out: f_vec = [f_cg, f_w1, f_t1, f_w2, f_t2, f_w3, f_t3, f_w4, f_t4] (centripetal components of the force)
	// out: normal_ext - normal over the full body; updated in this function from the centripetal forces

	// get distance vector between center of circle and the car = positions of points from the car vs center of circle, which is seen as 0
	car1->get_dist_vector(Position, dist_car_center);
	
	//// the distance to [cg<->center of circle] must be equal with the Radius of the circle - place an assert here!!!!
	//if (abs(car1->get_dist_vector_abs_val(Position) - Radius * Radius) > 1e-12) {
	//	std::cout << "\n\n error (in Circular.update_normal_force): the Radius of the circle is different than the distance  \\
	//		between the car and center of the circle!!!!please take care!!!";
	//}

	// vectors with velocities and masses
	car1->get_Velocity_vec(Velocity_vec);
	car1->get_Mass_vec(Mass_vec);

	// compute each of the 9 centripetal forces
	for (int i = 0; i < vec_DIM; ++i) {
		get_centrifugal_force(&f_vec[DIM * i], &Velocity_vec[DIM * i], *(Mass_vec + i), &dist_car_center[DIM * i]);
	}

	// compute centripetal part of the global normal force 
	// n = f_cg + f_w1 + f_t1 + f_w2 + f_t2 + f_w3 + f_t3 + f_w4 + f_t4
	cblas_dcopy(DIM, f_vec, incx, normal_ext, incx);
	for (auto i = 1; i < vec_DIM; ++i) {
		vdAdd(DIM, normal_ext, &f_vec[DIM * i], normal_ext);
	}
	
}

void Circular::get_Profile_torque(Car<double>* Car1, double* Torque) {
	Torque[2] = 0; // Torque on z direction
}

void Circular::update_initial_condition(Car<double>* Car1){

	double* perpendicular_dir = (double*)mkl_calloc(Car1->DIM, sizeof(double), Car1->alignment);
	double* tangential_dir = (double*)mkl_calloc(Car1->alignment, sizeof(double), Car1->alignment);
	double* radial_vector = (double*)mkl_calloc(Car1->alignment, sizeof(double), Car1->alignment);
	radial_vector[0] = Car1->Position_vec[0] - this->Position[0];
	radial_vector[1] = Car1->Position_vec[1] - this->Position[1];
	radial_vector[2] = 0;

	double radius = cblas_dnrm2(Car1->DIM, radial_vector, 1);
	perpendicular_dir[2] = 1;
	if (abs(radius - this->Radius) > 0.01)
		std::cout << "Warning! the initial position of the car is not on the trajectory provided in the circular path. \n The expected radius is " << this->Radius << ", but the car is at an initial distance of " << radius << " from the center of the circle.\n The execution procedes with the current spatial configuration and with the current distance to the center of the circle." << std::endl;

	double inv_radius_squared = 1. / (radius * radius);

	const MKL_INT incx = 1;
	MathLibrary::crossProduct(radial_vector, perpendicular_dir, tangential_dir);
	double magnitude = cblas_ddot(Car1->DIM, Car1->Velocity_vec, incx, tangential_dir, incx);
	cblas_dcopy(Car1->DIM, tangential_dir, 1, Car1->Velocity_vec, 1);
	cblas_dscal(Car1->DIM, magnitude, Car1->Velocity_vec, incx);
	MathLibrary::crossProduct(radial_vector, Car1->Velocity_vec, Car1->w_CG);
	cblas_dscal(Car1->DIM, inv_radius_squared, Car1->w_CG, 1);

	MKL_free(perpendicular_dir);
	MKL_free(tangential_dir);
	MKL_free(radial_vector);
}

void Circular::this_is_a_test_fn(std::string s) {
	std::cout << "This is a test function I got " <<s<< std::endl;
}
// =============================== end of Circular class implementation ===================
