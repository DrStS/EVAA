#pragma once
#include "car.h"
// classes: Profile, Circular, Car, Load_module etc.


class Profile {
protected:
	// consider replacing double with floatEVAA
	// alignment and g are defined globally === REPLACE THEM!!!!
	const int alignment = 64, mkl_DIM = 3, vec_DIM = 9, incx = 1;

	char* Name = NULL;

public:
	Profile& operator=(const Profile& Prof1) {};
	void update_Profile_force(Car* Car1, double* F_vec, double* Normal_ext) {};
};

class Circular : public Profile {
private:
	double* Position = NULL; // Center of the Circle (if the motion is a Circle)
	double Radius = 0; // Radius of the Circle (if the motion is a Circle)
public:
	Circular();
	Circular(double* Pos);
	Circular(double* Pos, double Rad);
	Circular(const Circular& Circ1);
	virtual Circular& operator=(const Circular& Circ1);
	~Circular();
	void get_Position(double* Pos);
	double get_Radius() const;
	void set_Radius(const double& Rad);
	void set_Position(const double* Pos);
	void get_centripet_force(double* Fr, double* v, double& m, double* p);
	void update_Profile_force(Car* Car1, double* F_vec, double* Normal_ext);
};

class Load_module {
private:
	// alignment and g are defined globally === REPLACE THEM!!!!
	const int alignment = 64;
	// MKL / vector constants:
	const int mkl_DIM = 3, vec_DIM = 9, incx = 1; // consider mkl_DIM = 4 for efficiency!!!

	Profile* Active_Profile;
	Car* Car_obj;

public:
	Load_module();
	Load_module(Profile* Profile_type);
	Load_module(Profile* Profile_type, Car* Car1);
	Load_module(const Load_module& Load_module_1);
	Load_module& operator= (const Load_module& Load_module_1);
	void set_Profile(Profile* Profile_type);
	void get_Profile(Profile* Profile_type);
	void update_force(double time_t, double* F_vec, double* Delta_x_vec, double* External_force);
};

