#include "profile_class.h"


// ===============================   Circular class implementation ===================
Circular::Circular() {
	Name = new char[strlen("Circular") + 1];
	Name = strcpy(Name, "Circular");
	Position = (double*)mkl_malloc(sizeof(double) * 3, alignment);
	Position[2] = Position[1] = Position[0] = 0;
	Radius = 100;

	// auxiliary vectors
	unit_y_vector = (double*)mkl_calloc(mkl_DIM, sizeof(double), alignment);
	unit_y_vector[0] = 0.; unit_y_vector[1] = 1.; unit_y_vector[2] = 0.;

	velocity_direction = (double*)mkl_calloc(mkl_DIM, sizeof(double), alignment);
	dist_car_center = (double*)mkl_malloc(sizeof(double) * mkl_DIM * vec_DIM, alignment);
	Velocity_vec = (double*)mkl_malloc(sizeof(double) * mkl_DIM * vec_DIM, alignment);
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

Circular::Circular(const Circular& circ1) {
	*this = circ1;
}

Circular& Circular::operator=(const Circular& circ1) {
	if (this != &circ1) {
		if (Name != NULL)
			delete[] Name;
		if (Position != NULL)
			delete[] Position;

		// Name
		if (circ1.Name != NULL) {
			Name = new char[strlen(circ1.Name) + 1];
			strcpy(Name, circ1.Name);
		}
		else
			Name = NULL;

		// Position
		if (circ1.Position != NULL) {
			Position = (double*)mkl_malloc(sizeof(double) * 3, alignment);
			cblas_dcopy(3, circ1.Position, 1, Position, 1);
		}
		else
			Position = NULL;

		// Radius
		Radius = circ1.Radius;
	}
	return *this;
}

Circular::~Circular() {
	if (Position != NULL) {
		mkl_free(Position);
		Position = NULL;
	}

	if (Name != NULL) {
		delete[] Name;
		Name = NULL;
	}

	if (unit_y_vector != NULL) {
		mkl_free(unit_y_vector);
		unit_y_vector = NULL;
	}

	if (velocity_direction != NULL) {
		mkl_free(velocity_direction);
		velocity_direction = NULL;
	}

	if (Velocity_vec != NULL) {
		mkl_free(Velocity_vec);
		Velocity_vec = NULL;
	}

	if (Mass_vec != NULL) {
		mkl_free(Mass_vec);
		Mass_vec = NULL;
	}

	if (dist_car_center != NULL) {
		mkl_free(dist_car_center);
		dist_car_center = NULL;
	}
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
		cblas_dcopy(mkl_DIM, pos, 1, Position, 1);
	}
}

void Circular::get_centrifugal_force(double* fr, double* v, double& m, double* p) {
	// calculates the force in a body only with respect to its velocity, mass and Position
	// v is the velocity of the body
	// m is the mass of the body
	// p is the global Position of the body
	// fr - the centripetal force
	// the rotation is always around the origin!

	cblas_dcopy(mkl_DIM, p, 1, fr, 1);

	fr[1] = 0;		// path only in xz-plane

	double inv_radius = 1.0 / cblas_dnrm2(mkl_DIM, p, 1);		// corresponds to the (inverse) Radius of the trajectory at the considered tyre

	// Raffi: cblas_dscal(mkl_DIM, -inv_radius, fr, 1);
	cblas_dscal(mkl_DIM, inv_radius, fr, 1);

	MathLibrary::crossProduct(fr, unit_y_vector, velocity_direction);

	double velocity_magnitude = cblas_ddot(mkl_DIM, v, incx, velocity_direction, incx);

	double force_magnitude = m * velocity_magnitude * velocity_magnitude * inv_radius;

	cblas_dscal(mkl_DIM, force_magnitude, fr, 1);
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
		get_centrifugal_force(&f_vec[mkl_DIM * i], &Velocity_vec[mkl_DIM * i], *(Mass_vec + i), &dist_car_center[mkl_DIM * i]);
	}

	// compute centripetal part of the global normal force 
	// n = f_cg + f_w1 + f_t1 + f_w2 + f_t2 + f_w3 + f_t3 + f_w4 + f_t4
	cblas_dcopy(mkl_DIM, f_vec, incx, normal_ext, incx);
	for (auto i = 1; i < vec_DIM; ++i) {
		vdAdd(mkl_DIM, normal_ext, &f_vec[mkl_DIM * i], normal_ext);
	}
}

void Circular::get_Profile_torque(Car<double>* Car1, double* Torque) {
	Torque[2] = 0; // Torque on z direction
}
// =============================== end of Circular class implementation ===================
