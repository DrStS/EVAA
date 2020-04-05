#pragma once

#include "car.h"
#include "Profile_class.h"

class Load_module {
private:
	// alignment and g are defined globally === REPLACE THEM!!!!
	const int alignment = 64;
	// MKL / vector constants:
	const int DIM = 3, vec_DIM = 9, incx = 1; // consider mkl_DIM = 4 for efficiency!!!

	Profile* Active_Profile;
	Car<double>* Car_obj;

	

	// Auxiliary vectors
	double* Normal_ext = NULL; // Normal_force updated with external forces
	double* k_vec = NULL; // vector with springs' stiffnesses (copied from Car!!!)
	double* Normal_from_profile = NULL; // Normal_force computed inside get_Profile_force
	double* External_force = NULL;

public:
	Load_module(Profile* Profile_type, Car<double>* Car1, Load_Params load_param);
	~Load_module();
	void set_Profile(Profile* Profile_type);
	void get_Profile(Profile* Profile_type);
	void update_force(double time_t, double* F_vec, double* Delta_x_vec, double* Normal_ext);
	void update_torque(double time_t, double* Torque, double* Delta_x_vec);
	void set_External_force(Load_Params load_param);
	void test() {
		std::cout << "External force:\n";
		MathLibrary::write_vector(External_force, 3);

		std::cout << "\n\nActive Profile test:\n";
		Active_Profile->test();
	}
};

