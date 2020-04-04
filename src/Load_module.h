#pragma once

#include "car.h"
#include "Profile_class.h"

class Load_module {
private:
	// alignment and g are defined globally === REPLACE THEM!!!!
	const int alignment = 64;
	// MKL / vector constants:
	const int mkl_DIM = 3, vec_DIM = 9, incx = 1; // consider mkl_DIM = 4 for efficiency!!!

	Profile* Active_Profile = NULL;
	Car<double>* Car_obj;

	double* External_force = NULL;

	// Auxiliary vectors
	double* Normal_ext = NULL; // Normal_force computed inside update_Profile_force
	double* k_vec = NULL; // vector with springs' stiffnesses (copied from Car!!!)


public:
	Load_module(Profile* Profile_type, Car<double>* Car1, Load_Params load_param);
	~Load_module();
	Load_module(const Load_module& Load_module_1);
	Load_module& operator= (const Load_module& Load_module_1);
	void set_Profile(Profile* Profile_type);
	void get_Profile(Profile* Profile_type);
	void update_force(double time_t, double* F_vec, double* Delta_x_vec);
	void update_torque(double time_t, double* F_vec, double* Delta_x_vec);
	void set_External_force(Load_Params load_param);
};

