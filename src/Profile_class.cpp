#include "profile_class.h"


// ===============================   circular class implementation ===================
circular::circular() {
	name = new char[strlen("circular") + 1];
	name = strcpy(name, "circular");
	position = (double*)mkl_malloc(sizeof(double) * 3, alignment);
	position[2] = position[1] = position[0] = 0;
	radius = 100;

	// auxiliary vectors
	unit_y_vector = (double*)mkl_calloc(mkl_dim, sizeof(double), alignment);
	unit_y_vector[0] = 0.; unit_y_vector[1] = 1.; unit_y_vector[2] = 0.;

	velocity_direction = (double*)mkl_calloc(mkl_dim, sizeof(double), alignment);
	dist_car_center = (double*)mkl_malloc(sizeof(double) * mkl_dim * vec_dim, alignment);
	velocity_vec = (double*)mkl_malloc(sizeof(double) * mkl_dim * vec_dim, alignment);
	mass_vec = (double*)mkl_malloc(sizeof(double) * vec_dim, alignment);
}

circular::circular(double* pos) : circular::circular() {
	// position
	if (pos != null) {
		cblas_dcopy(3, pos, 1, position, 1);
	}
}

circular::circular(double* pos, double rad) : circular::circular(pos) {
	radius = rad;
};

circular::circular(const circular& circ1) {
	*this = circ1;
}

circular& circular::operator=(const circular& circ1) {
	if (this != &circ1) {
		if (name != null)
			delete[] name;
		if (position != null)
			delete[] position;

		// name
		if (circ1.name != null) {
			name = new char[strlen(circ1.name) + 1];
			strcpy(name, circ1.name);
		}
		else
			name = null;

		// position
		if (circ1.position != null) {
			position = (double*)mkl_malloc(sizeof(double) * 3, alignment);
			cblas_dcopy(3, circ1.position, 1, position, 1);
		}
		else
			position = null;

		// radius
		radius = circ1.radius;
	}
	return *this;
}

circular::~circular() {
	if (position != null) {
		mkl_free(position);
		position = null;
	}

	if (name != null) {
		delete[] name;
		name = null;
	}

	if (unit_y_vector != null) {
		mkl_free(unit_y_vector);
		unit_y_vector = null;
	}

	if (velocity_direction != null) {
		mkl_free(velocity_direction);
		velocity_direction = null;
	}

	if (velocity_vec != null) {
		mkl_free(velocity_vec);
		velocity_vec = null;
	}

	if (mass_vec != null) {
		mkl_free(mass_vec);
		mass_vec = null;
	}

	if (dist_car_center != null) {
		mkl_free(dist_car_center);
		dist_car_center = null;
	}
}

void circular::get_position(double* pos) {
	if (pos != null) {
		cblas_dcopy(3, position, 1, pos, 1);
	}
}

double circular::get_radius() const {
	return radius;
}

void circular::set_radius(const double& rad) {
	radius = rad;
}

void circular::set_position(const double* pos) {
	if (pos != null) {
		cblas_dcopy(mkl_dim, pos, 1, position, 1);
	}
}

void circular::get_centripet_force(double* fr, double* v, double& m, double* p) {
	// calculates the force in a body only with respect to its velocity, mass and position
	// v is the velocity of the body
	// m is the mass of the body
	// p is the global position of the body
	// fr - the centripetal force
	// the rotation is always around the origin!

	cblas_dcopy(mkl_dim, p, 1, fr, 1);

	fr[1] = 0;		// path only in xz-plane

	double inv_radius = 1.0 / cblas_dnrm2(mkl_dim, p, 1);		// corresponds to the (inverse) radius of the trajectory at the considered tyre

	cblas_dscal(mkl_dim, -inv_radius, fr, 1);

	mathlibrary::crossproduct(fr, unit_y_vector, velocity_direction);

	double velocity_magnitude = cblas_ddot(mkl_dim, v, incx, velocity_direction, incx);

	double force_magnitude = m * velocity_magnitude * velocity_magnitude * inv_radius;

	cblas_dscal(mkl_dim, force_magnitude, fr, 1);
}

void circular::update_profile_force(car<double>* car1, double* f_vec, double* normal_ext) {
	// out: f_vec = [f_cg, f_w1, f_t1, f_w2, f_t2, f_w3, f_t3, f_w4, f_t4] (centripetal components of the force)
	// out: normal_ext - normal over the full body; updated in this function from the centripetal forces

	// get distance vector between center of circle and the car = positions of points from the car vs center of circle, which is seen as 0
	car1->get_dist_vector(position, dist_car_center);

	//// the distance to [cg<->center of circle] must be equal with the radius of the circle - place an assert here!!!!
	//if (abs(car1->get_dist_vector_abs_val(position) - radius * radius) > 1e-12) {
	//	std::cout << "\n\n error (in circular.update_normal_force): the radius of the circle is different than the distance  \\
	//		between the car and center of the circle!!!!please take care!!!";
	//}

	// vectors with velocities and masses
	car1->get_velocity_vec(velocity_vec);
	car1->get_mass_vec(mass_vec);

	// compute each of the 9 centripetal forces
	for (int i = 0; i < vec_dim; ++i) {
		get_centripet_force(&f_vec[mkl_dim * i], &velocity_vec[mkl_dim * i], *(mass_vec + i), &dist_car_center[mkl_dim * i]);
	}

	// compute centripetal part of the global normal force 
	// n = f_cg + f_w1 + f_t1 + f_w2 + f_t2 + f_w3 + f_t3 + f_w4 + f_t4
	cblas_dcopy(mkl_dim, f_vec, incx, normal_ext, incx);
	for (auto i = 1; i < vec_dim; ++i) {
		vdadd(mkl_dim, normal_ext, &f_vec[mkl_dim * i], normal_ext);
	}
}

// =============================== end of circular class implementation ===================
