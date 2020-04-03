#import "Profile_class.h"

// ===============================   Circular class implementation ===================
Circular::Circular() {
	Name = new char[strlen("Circular") + 1];
	Name = strcpy(Name, "Circular");
	Position = (double*)mkl_malloc(sizeof(double) * 3, alignment);
	Position[2] = Position[1] = Position[0] = 0;
	Radius = 100;

	// Auxiliary vectors
	unit_y_vector = (double*)mkl_calloc(mkl_DIM, sizeof(double), alignment);
	unit_y_vector[0] = 0.; unit_y_vector[1] = 1.; unit_y_vector[2] = 0.;

	velocity_direction = (double*)mkl_calloc(mkl_DIM, sizeof(double), alignment);
	dist_car_center = (double*)mkl_malloc(sizeof(double) * mkl_DIM * vec_DIM, alignment);
	Velocity_vec = (double*)mkl_malloc(sizeof(double) * mkl_DIM * vec_DIM, alignment);
	Mass_vec = (double*)mkl_malloc(sizeof(double) * vec_DIM, alignment);
}

Circular::Circular(double* Pos) : Circular::Circular() {
	// Position
	if (Pos != NULL) {
		cblas_dcopy(3, Pos, 1, Position, 1);
	}
}

Circular::Circular(double* Pos, double Rad) : Circular::Circular(Pos) {
	Radius = Rad;
};

Circular::Circular(const Circular& Circ1) {
	*this = Circ1;
}

Circular& Circular::operator=(const Circular& Circ1) {
	if (this != &Circ1) {
		if (Name != NULL)
			delete[] Name;
		if (Position != NULL)
			delete[] Position;

		// Name
		if (Circ1.Name != NULL) {
			Name = new char[strlen(Circ1.Name) + 1];
			strcpy(Name, Circ1.Name);
		}
		else
			Name = NULL;

		// Position
		if (Circ1.Position != NULL) {
			Position = (double*)mkl_malloc(sizeof(double) * 3, alignment);
			cblas_dcopy(3, Circ1.Position, 1, Position, 1);
		}
		else
			Position = NULL;

		// Radius
		Radius = Circ1.Radius;
	}
	return *this;
}

Circular::~Circular() {
	if (Position != NULL) {
		MKL_free(Position);
		Position = NULL;
	}

	if (Name != NULL) {
		delete[] Name;
		Name = NULL;
	}

	if (unit_y_vector != NULL) {
		MKL_free(unit_y_vector);
		unit_y_vector = NULL;
	}

	if (velocity_direction != NULL) {
		MKL_free(velocity_direction);
		velocity_direction = NULL;
	}

	if (Velocity_vec != NULL) {
		MKL_free(Velocity_vec);
		Velocity_vec = NULL;
	}

	if (Mass_vec != NULL) {
		MKL_free(Mass_vec);
		Mass_vec = NULL;
	}

	if (dist_car_center != NULL) {
		MKL_free(dist_car_center);
		dist_car_center = NULL;
	}
}

void Circular::get_Position(double* Pos) {
	if (Pos != NULL) {
		cblas_dcopy(3, Position, 1, Pos, 1);
	}
}

double Circular::get_Radius() const {
	return Radius;
}

void Circular::set_Radius(const double& Rad) {
	Radius = Rad;
}

void Circular::set_Position(const double* Pos) {
	if (Pos != NULL) {
		cblas_dcopy(mkl_DIM, Pos, 1, Position, 1);
	}
}

void Circular::get_centripet_force(double* Fr, double* v, double& m, double* p) {
	// calculates the force in a body only with respect to its velocity, mass and position
	// v is the velocity of the body
	// m is the mass of the body
	// p is the global position of the body
	// Fr - the centripetal force
	// The rotation is always around the origin!

	cblas_dcopy(mkl_DIM, p, 1, Fr, 1);

	Fr[1] = 0;		// path only in XZ-plane

	double inv_radius = 1.0 / cblas_dnrm2(mkl_DIM, p, 1);		// corresponds to the (inverse) radius of the trajectory at the considered tyre

	cblas_dscal(mkl_DIM, -inv_radius, Fr, 1);

	MathLibrary::crossProduct(Fr, unit_y_vector, velocity_direction);

	double velocity_magnitude = cblas_ddot(mkl_DIM, v, incx, velocity_direction, incx);

	double force_magnitude = m * velocity_magnitude * velocity_magnitude * inv_radius;

	cblas_dscal(mkl_DIM, force_magnitude, Fr, 1);
}

void Circular::update_Profile_force(Car* Car1, double* F_vec, double* Normal_ext) {
	// out: F_vec = [F_CG, F_W1, F_T1, F_W2, F_T2, F_W3, F_T3, F_W4, F_T4] (Centripetal components of the force)
	// out: Normal_ext - normal over the full body; updated in this function from the centripetal forces

	// get distance vector between center of circle and the car = Positions of points from the car vs Center of Circle, which is seen as 0
	Car1->get_dist_vector(Position, dist_car_center);

	//// the distance to [CG<->Center of Circle] MUST BE EQUAL with the RADIUS of the Circle - place an ASSERT here!!!!
	//if (abs(Car1->get_dist_vector_abs_val(Position) - Radius * Radius) > 1e-12) {
	//	std::cout << "\n\n ERROR (in Circular.update_normal_force): The radius of the circle is different than the distance  \\
	//		between the car and center of the circle!!!!Please take care!!!";
	//}

	// Vectors with Velocities and Masses
	Car1->get_Velocity_vec(Velocity_vec);
	Car1->get_Mass_vec(Mass_vec);

	// Compute each of the 9 centripetal forces
	for (int i = 0; i < vec_DIM; ++i) {
		get_centripet_force(&F_vec[mkl_DIM * i], &Velocity_vec[mkl_DIM * i], *(Mass_vec + i), &dist_car_center[mkl_DIM * i]);
	}

	// Compute centripetal part of the global normal force 
	// N = F_CG + F_W1 + F_T1 + F_W2 + F_T2 + F_W3 + F_T3 + F_W4 + F_T4
	cblas_dcopy(mkl_DIM, F_vec, incx, Normal_ext, incx);
	for (auto i = 1; i < vec_DIM; ++i) {
		vdAdd(mkl_DIM, Normal_ext, &F_vec[mkl_DIM * i], Normal_ext);
	}
}

// =============================== END OF Circular class implementation ===================
