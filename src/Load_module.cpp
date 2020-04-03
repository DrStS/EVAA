#pragma once

#include <iostream>
#include "Load_module.h"
#include "MathLibrary.h"

#ifdef USE_INTEL_MKL
#include <mkl.h>
#define USE_GEMM
#endif

//// Circular, Car, Load_module etc.

// ===============================  Car class implementation =======================
Car::Car() {
	// Dimensions
	Length = 4.;
	Width = 2.;

	// Position
	Position_vec = (double*)mkl_malloc(sizeof(double) * mkl_DIM * vec_DIM, alignment);
	// Default Positions - CG = 100m (= On the circle(0., 100m))
	for (auto i = 0; i < mkl_DIM * vec_DIM; ++i) {
		Position_vec[3 * i] = 100. / 3;
		if (i == 1 || i == 2 || i == 3 || i == 4) // x-positions of {fl, fr} x {wheel, tyre}
			Position_vec[3 * i] = 100. / 3 + Length / 2;
		else if (i == 5 || i == 6 || i == 7 || i == 8) // x-positions of {rl, rr} x {wheel, tyre}
			Position_vec[3 * i] = 100. / 3 - Length / 2;
		else
			Position_vec[3 * i] = 100. / 3;

		Position_vec[3 * i + 1] = 200. / 3;

		Position_vec[3 * i + 2] = 200. / 3;

		if (i == 3 || i == 4 || i == 7 || i == 8) // z-positions of {fr, rr} x {wheel, tyre}
			Position_vec[3 * i + 2] = 200. / 3 + Width / 2;
		else if (i == 1 || i == 2 || i == 5 || i == 6) // z-positions of {fl, rl} x {wheel, tyre}
			Position_vec[3 * i + 2] = 200. / 3 + Width / 2;
		else
			Position_vec[3 * i + 2] = 200. / 3;
	}

	// Velocity - 10 m/s (random)
	Velocity_vec = (double*)mkl_malloc(sizeof(double) * mkl_DIM * vec_DIM, alignment);
	for (auto i = 0; i < 9; ++i) {
		Velocity_vec[3 * i] = 8;
		Velocity_vec[3 * i + 1] = 0;
		Velocity_vec[3 * i + 2] = 6;
	}

	// Mass
	Mass_vec = (double*)mkl_malloc(sizeof(double) * vec_DIM, alignment);
	Mass_vec[0] = 1936.; // Body
	Mass_vec[1] = 145. / 2; // wheel_fl
	Mass_vec[2] = 0.; // tire_fl
	Mass_vec[3] = 145. / 2; // wheel_fr
	Mass_vec[4] = 0.; // tire_fr
	Mass_vec[5] = 135. / 2; // wheel_rl
	Mass_vec[6] = 0.; // tire_rl
	Mass_vec[7] = 135. / 2; // wheel_rr
	Mass_vec[8] = 0.; // tire_rr

	// Springs
	// Stiffnesses
	k_vec = (double*)mkl_malloc(sizeof(double) * (vec_DIM - 1), alignment);
	for (auto i = 0; i < 4; ++i) {
		k_vec[2 * i] = 260e3; // Tyre
		k_vec[2 * i + 1] = 28e3 * 0.69; // Wheel
	}
	// Length of springs (?)
	L_body_fl = 1;
	L_tire_fl = 1;
	L_body_fr = 1;
	L_tire_fr = 1;
	L_body_rl = 1;
	L_tire_rl = 1;
	L_body_rr = 1;
	L_tire_rr = 1;

	// Distances from CG to wheels
	l_long_fl = 1.395;
	l_long_fr = 1.395;
	l_long_rl = 1.596;
	l_long_rr = 1.596;
	l_lat_fl = 2 * 0.8458;
	l_lat_fr = 2 * 0.8458;
	l_lat_rl = 2 * 0.84;
	l_lat_rr = 2 * 0.84;

	// 	// Moments of inertia
	I_body_xx = 640;
	I_body_yy = 4800;
};

Car::Car(double* Pos, double* Vel) : Car::Car() {
	// Position
	if (Pos != NULL) {
		cblas_dcopy(mkl_DIM * vec_DIM, Pos, 1, Position_vec, 1);
	}

	// Velocity
	if (Vel != NULL) {
		cblas_dcopy(mkl_DIM * vec_DIM, Vel, 1, Velocity_vec, 1);
	}
}

Car::Car(double* Pos, double* Vel, double* M) :Car::Car(Pos, Vel) {
	// Mass
	if (M != NULL) {
		cblas_dcopy(vec_DIM, M, 1, Mass_vec, 1);
	}
}

Car::Car(double* Pos, double* Vel, double* M, double* k) : Car::Car(Pos, Vel, M) {
	// Springs stiffnesses
	if (k != NULL) {
		cblas_dcopy(vec_DIM, k, 1, k_vec, 1);
	}
};

Car::Car(double* Pos, double* Vel, double* M, double* k,
	double Len, double Wid,
	double L_body_fl_val, double L_tire_fl_val,
	double L_body_fr_val, double L_tire_fr_val,
	double L_body_rl_val, double L_tire_rl_val,
	double L_body_rr_val, double L_tire_rr_val) : Car::Car(Pos, Vel, M, k) {
	Length = Len;
	Width = Wid;
	L_body_fl = L_body_fl_val;
	L_tire_fl = L_tire_fl_val;
	L_body_fr = L_body_fr_val;
	L_tire_fr = L_tire_fr_val;
	L_body_rl = L_body_rl_val;
	L_tire_rl = L_tire_rl_val;
	L_body_rr = L_body_rr_val;
	L_tire_rr = L_tire_rr_val;
};

Car::Car(double* Pos, double* Vel, double* M, double* k,
	double Len, double Wid,
	double L_body_fl_val, double L_tire_fl_val,
	double L_body_fr_val, double L_tire_fr_val,
	double L_body_rl_val, double L_tire_rl_val,
	double L_body_rr_val, double L_tire_rr_val,
	double l_long_fl_val, double l_long_fr_val,
	double l_long_rl_val, double l_long_rr_val,
	double l_lat_fl_val, double l_lat_fr_val,
	double l_lat_rl_val, double l_lat_rr_val) : Car::Car(Pos, Vel, M, k, Len, Wid,
		L_body_fl_val, L_tire_fl_val,
		L_body_fr_val, L_tire_fr_val,
		L_body_rl_val, L_tire_rl_val,
		L_body_rr_val, L_tire_rr_val) {

	l_long_fl = l_long_fl_val;
	l_long_fr = l_long_fr_val;
	l_long_rl = l_long_rl_val;
	l_long_rr = l_long_rr_val;
	l_lat_fl = l_lat_fl_val;
	l_lat_fr = l_lat_fr_val;
	l_lat_rl = l_lat_rl_val;
	l_lat_rr = l_lat_rr_val;
};

Car::Car(double* Pos, double* Vel, double* M, double* k,
	double Len, double Wid,
	double L_body_fl_val, double L_tire_fl_val,
	double L_body_fr_val, double L_tire_fr_val,
	double L_body_rl_val, double L_tire_rl_val,
	double L_body_rr_val, double L_tire_rr_val,
	double l_long_fl_val, double l_long_fr_val,
	double l_long_rl_val, double l_long_rr_val,
	double l_lat_fl_val, double l_lat_fr_val,
	double l_lat_rl_val, double l_lat_rr_val,
	double 	I_body_xx_val, double I_body_yy_val) : Car::Car(Pos, Vel, M, k, Len, Wid,
		L_body_fl_val, L_tire_fl_val,
		L_body_fr_val, L_tire_fr_val,
		L_body_rl_val, L_tire_rl_val,
		L_body_rr_val, L_tire_rr_val,
		l_long_fl_val, l_long_fr_val,
		l_long_rl_val, l_long_rr_val,
		l_lat_fl_val, l_lat_fr_val,
		l_lat_rl_val, l_lat_rr_val) {
	I_body_xx = I_body_xx_val;
	I_body_yy = I_body_yy_val;
};


Car::Car(const Car& Car1) {
	*this = Car1;
};

Car& Car::operator= (const Car& Car1) {
	if (this != &Car1) {
		// Position
		if (Car1.Position_vec != NULL) {
			Position_vec = (double*)mkl_malloc(sizeof(double) * mkl_DIM * vec_DIM, alignment);
			cblas_dcopy(mkl_DIM * vec_DIM, Car1.Position_vec, 1, Position_vec, 1);
		}

		// Velocity
		if (Car1.Velocity_vec != NULL) {
			Velocity_vec = (double*)mkl_malloc(sizeof(double) * mkl_DIM * vec_DIM, alignment);
			cblas_dcopy(mkl_DIM * vec_DIM, Car1.Velocity_vec, 1, Velocity_vec, 1);
		}

		// Mass vector
		if (Car1.Mass_vec != NULL) {
			Mass_vec = (double*)mkl_malloc(sizeof(double) * vec_DIM, alignment);
			cblas_dcopy(mkl_DIM * vec_DIM, Car1.Mass_vec, 1, Mass_vec, 1);
		}

		// Stifnesses vector
		if (Car1.k_vec != NULL) {
			k_vec = (double*)mkl_malloc(sizeof(double) * vec_DIM, alignment);
			cblas_dcopy(mkl_DIM * vec_DIM, Car1.k_vec, 1, k_vec, 1);
		}

		// Others
		Length = Car1.Length;
		Width = Car1.Width;

		L_body_fl = Car1.L_body_fl;
		L_tire_fl = Car1.L_tire_fl;
		L_body_fr = Car1.L_body_fr;
		L_tire_fr = Car1.L_tire_fr;
		L_body_rl = Car1.L_body_rl;
		L_tire_rl = Car1.L_tire_rl;
		L_body_rr = Car1.L_body_rr;
		L_tire_rr = Car1.L_tire_rr;

		l_long_fl = Car1.l_long_fl;
		l_long_fr = Car1.l_long_fr;
		l_long_rl = Car1.l_long_rl;
		l_long_rr = Car1.l_long_rr;
		l_lat_fl = Car1.l_lat_fl;
		l_lat_fr = Car1.l_lat_fr;
		l_lat_rl = Car1.l_lat_rl;
		l_lat_rr = Car1.l_lat_rr;

		I_body_xx = Car1.I_body_xx;
		I_body_yy = Car1.I_body_yy;
	}

	return *this;
}

Car::~Car() {
	if (Position_vec != NULL) {
		MKL_free(Position_vec);
		Position_vec = NULL;
	}

	if (Velocity_vec != NULL) {
		MKL_free(Velocity_vec);
		Velocity_vec = NULL;
	}

	if (Mass_vec != NULL) {
		MKL_free(Mass_vec);
		Mass_vec = NULL;
	}

	if (k_vec != NULL) {
		MKL_free(k_vec);
		k_vec = NULL;
	}
}

// Methods for Positions
void Car::get_Position_vec(double* Pos) {
	if (Pos != NULL) {
		cblas_dcopy(mkl_DIM * vec_DIM, Position_vec, 1, Pos, 1);
	}
}

void Car::set_Position_vec(const double* Pos) {
	if (Pos != NULL) {
		cblas_dcopy(mkl_DIM * vec_DIM, Pos, 1, Position_vec, 1);
	}
}

void Car::get_Position_vec_CG(double* Pos_CG) {
	if (Pos_CG != NULL) {
		cblas_dcopy(mkl_DIM, Position_vec, 1, Pos_CG, 1);
	}
}

void Car::set_Position_vec_CG(const double* Pos_CG) {
	if (Pos_CG != NULL) {
		cblas_dcopy(mkl_DIM, Pos_CG, 1, Position_vec, 1);
	}
}

// Methods for Velocities
void Car::get_Velocity_vec(double* Vel) {
	if (Vel != NULL) {
		// b=a, cblas_dcopy(n,a,inc,b,inc)
		cblas_dcopy(mkl_DIM * vec_DIM, Velocity_vec, 1, Vel, 1);
	}
}

void Car::set_Velocity_vec(const double* Vel) {
	if (Vel != NULL) {
		cblas_dcopy(mkl_DIM * vec_DIM, Vel, 1, Velocity_vec, 1);
	}
}

void Car::get_Velocity_vec_CG(double* Vel_CG) {
	if (Vel_CG != NULL) {
		// b=a, cblas_dcopy(n,a,inc,b,inc)
		cblas_dcopy(mkl_DIM, Velocity_vec, 1, Vel_CG, 1);
	}
}

void Car::set_Velocity_vec_CG(const double* Vel_CG) {
	if (Vel_CG != NULL) {
		cblas_dcopy(mkl_DIM, Vel_CG, 1, Velocity_vec, 1);
	}
}

// Methods for Stifnesses k_vec
void Car::get_k_vec(double* k) {
	if (k != NULL) {
		// b=a, cblas_dcopy(n,a,inc,b,inc)
		cblas_dcopy(mkl_DIM * vec_DIM, k_vec, 1, k, 1);
	}
}

void Car::set_k_vec(const double* k) {
	if (k != NULL) {
		cblas_dcopy(mkl_DIM * vec_DIM, k, 1, k_vec, 1);
	}
}


// Methods for Masses
void Car::get_Mass_vec(double* M) {
	if (M != NULL) {
		// b=a, cblas_dcopy(n,a,inc,b,inc)
		cblas_dcopy(mkl_DIM * vec_DIM, Mass_vec, 1, M, 1);
	}
}

void Car::set_Mass_vec(const double* M) {
	if (M != NULL) {
		cblas_dcopy(mkl_DIM * vec_DIM, M, 1, Mass_vec, 1);
	}
}

double Car::get_Mass_vec_CG() const {
	return Mass_vec[0];
}

// Methods to get distances
void Car::get_dist_vector(double* Point_P, double* dist_vector) {
	// get distance vector from each important Point of the car (9: CG, 4*W_i, 4*T_i)
	// source: Point_P, dest: each entry from Position_vec
	if (Point_P != NULL && dist_vector != NULL) {
		for (auto i = 0; i < vec_DIM; ++i) {
			cblas_dcopy(mkl_DIM, Point_P, incx, &dist_vector[mkl_DIM * i], incx);
		}
		// y=a-b, vdSub(n,a,b,y)
		vdSub(mkl_DIM * vec_DIM, Position_vec, dist_vector, dist_vector);
	}
}

void Car::get_dist_vector_CG(double* Point_P, double* dist_vector) {
	// get distance vector from Center of Gravity of the car to a Point P 
	// source: Point_P, dest: CG
	if (Point_P != NULL && dist_vector != NULL) {
		// y=a-b, vdSub(n,a,b,y)
		vdSub(mkl_DIM, Position_vec, Point_P, dist_vector);
	}
}

// Methods for moments of inertia
double Car::get_I_body_xx() const {
	return I_body_xx;
}

double Car::get_I_body_yy() const {
	return I_body_yy;
}

void Car::set_I_body_xx(const double& I_body_xx_val) {
	I_body_xx = I_body_xx_val;
}

void Car::set_I_body_yy(const double& I_body_yy_val) {
	I_body_yy = I_body_yy_val;
}


// =============================== END of Car class implementation =======================

// ===============================   Circular class implementation ===================
Circular::Circular() {
	Name = new char[strlen("Circular") + 1];
	Name = strcpy(Name, "Circular");
	Position = (double*)mkl_malloc(sizeof(double) * 3, alignment);
	Position[2] = Position[1] = Position[0] = 0;
	Radius = 100;
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
		free(Name);
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

	const MKL_INT mkl_DIM = 3;
	const MKL_INT mkl_incx = 1;
	const MKL_INT mkl_incy = 1;

	double* unit_y_vector = (double*)mkl_calloc(mkl_DIM, sizeof(double), this->alignment);
	double* velocity_direction = (double*)mkl_calloc(mkl_DIM, sizeof(double), this->alignment);

	double velocity_magnitude, inv_radius, force_magnitude;

	cblas_dcopy(mkl_DIM, p, 1, Fr, 1);

	Fr[1] = 0;		// path only in XZ-plane

	inv_radius = 1.0 / cblas_dnrm2(mkl_DIM, p, 1);		// corresponds to the (inverse) radius of the trajectory at the considered tyre

	cblas_dscal(mkl_DIM, -inv_radius, Fr, 1);

	unit_y_vector[0] = 0.;
	unit_y_vector[1] = 1.;
	unit_y_vector[2] = 0.;

	MathLibrary::crossProduct(Fr, unit_y_vector, velocity_direction);

	velocity_magnitude = cblas_ddot(mkl_DIM, v, mkl_incx, velocity_direction, mkl_incy);

	force_magnitude = m * velocity_magnitude * velocity_magnitude * inv_radius;

	cblas_dscal(mkl_DIM, force_magnitude, Fr, 1);

	MKL_free(unit_y_vector);
	MKL_free(velocity_direction);
}

void Circular::update_Profile_force(Car* Car1, double* F_vec, double* Normal_ext) {
	// out: F_vec = [F_CG, F_W1, F_T1, F_W2, F_T2, F_W3, F_T3, F_W4, F_T4] (Centripetal components of the force)
	// out: Normal_ext - normal over the full body; updated in this function from the centripetal forces

	// get distance vector between center of circle and the car = Positions of points from the car vs Center of Circle, which is seen as 0
	double* dist_car_center = (double*)mkl_malloc(sizeof(double) * mkl_DIM * vec_DIM, alignment);
	Car1->get_dist_vector(Position, dist_car_center);

	//// the distance to [CG<->Center of Circle] MUST BE EQUAL with the RADIUS of the Circle - place an ASSERT here!!!!
	//if (abs(Car1->get_dist_vector_abs_val(Position) - Radius * Radius) > 1e-12) {
	//	std::cout << "\n\n ERROR (in Circular.update_normal_force): The radius of the circle is different than the distance  \\
	//		between the car and center of the circle!!!!Please take care!!!";
	//}

	// Vectors with Velocities and Masses
	double* Velocity_vec = (double*)mkl_malloc(sizeof(double) * mkl_DIM * vec_DIM, alignment);
	double* Mass_vec = (double*)mkl_malloc(sizeof(double) * vec_DIM, alignment);

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

	MKL_free(dist_car_center);
	MKL_free(Velocity_vec);
	MKL_free(Mass_vec);
}

// =============================== END OF Circular class implementation ===================

// ===============================  Load_module class implementation =======================

// IN:  [F, T, W1, T1, W2, T2, W3, T3, W4, T4]
// OUT: [F, T, W1, T1, W2, T2, W3, T3, W4, T4]

Load_module::Load_module() {
	Active_Profile = new Profile();
	Car_obj = new Car();
}

Load_module::Load_module(Profile* Profile_type) : Load_module() {
	Active_Profile = Profile_type;
}

Load_module::Load_module(Profile* Profile_type, Car* Car1) : Load_module(Profile_type) {
	Car_obj = Car1;
}

Load_module::Load_module(const Load_module& Load_module_1) {
	*this = Load_module_1;
};

Load_module& Load_module::operator= (const Load_module& Load_module_1) {
	Active_Profile = Load_module_1.Active_Profile;
	Car_obj = Load_module_1.Car_obj;

	return *this;
}

void Load_module::set_Profile(Profile* Profile_type) {
	Active_Profile = Profile_type;
}

void Load_module::get_Profile(Profile* Profile_type) {
	Profile_type = Active_Profile;
}


void Load_module::update_force(double time_t, double* F_vec, double* Delta_x_vec, double* External_force) {
	double* Normal_ext = (double*)mkl_malloc(sizeof(double) * mkl_DIM, alignment); // Normal_force
	Active_Profile->update_Profile_force(Car_obj, F_vec, Normal_ext);

	// N += External_force
	if (External_force != NULL) {
		vdAdd(mkl_DIM, External_force, Normal_ext, Normal_ext);
	}
	// F_Ti += -0.25 * N; [F[2], F[4], F[6], F[8]]
	for (auto i = 2; i < vec_DIM; i += 2) {
		cblas_daxpy(mkl_DIM, -0.25, Normal_ext, incx, &F_vec[i * mkl_DIM], incx);
	}
	
	// get Stiffnesses vector k_vec
	double* k_vec = (double*)mkl_malloc(sizeof(double) * (vec_DIM - 1), alignment); // k_vec; vec_DIM-1 = 8
	Car_obj->get_k_vec(k_vec);

	for (auto i = 0; i < (vec_DIM - 1); i += 2) {
		// use the elastic forces at wheels
		// F_CG += k_wi * Delta_x_i
		cblas_daxpy(mkl_DIM, k_vec[i], &Delta_x_vec[mkl_DIM * i], incx, F_vec, incx);
		// F_w_i -= k_wi * Delta_x_i
		cblas_daxpy(mkl_DIM, -k_vec[i], &Delta_x_vec[mkl_DIM * i], incx, &F_vec[mkl_DIM * (i + 1)], incx);

		// use the elastic forces at tyres
		// F_w_i -= k_t_i * Delta_x_{i+1}
		cblas_daxpy(mkl_DIM, -k_vec[i + 1], &Delta_x_vec[mkl_DIM * (i + 1)], incx, &F_vec[mkl_DIM * (i + 1)], incx);
		// F_t_i += k_t_i * Delta_x_{i+1}
		cblas_daxpy(mkl_DIM, k_vec[i + 1], &Delta_x_vec[mkl_DIM * (i + 1)], incx, &F_vec[mkl_DIM * (i + 2)], incx);
	}


	MKL_free(Normal_ext);
	MKL_free(k_vec);

}

// =============================== END of Load_module class implementation =======================