#pragma once
#include <mkl.h>
#include <string>
#include "MathLibrary.h"
#include "car.h"

class Profile {
protected:
	// consider replacing double with floatEVAA
	// alignment and g are defined globally === REPLACE THEM!!!!
	const int alignment = 64, DIM = 3, vec_DIM = 9, incx = 1;

	std::string Name;

public:
	virtual void get_Profile_force(Car<double>* Car1, double* F_vec, double* Normal_ext) {};
	virtual void get_Profile_torque(Car<double>* Car1, double* Torque_vec) {};
	virtual void update_initial_condition(Car<double>* Car1) { std::cout << "I am fucked" << std::endl; };
	virtual void this_is_a_test_fn(std::string s) {};
	virtual ~Profile() {};
};

class Circular : public Profile {
private:
	double* Position = NULL; // Center of the Circle (if the motion is a Circle)
	double Radius = 0; // Radius of the Circle (if the motion is a Circle)

	// vectors used in computations inside the methods
	double* unit_y_vector;
	double* velocity_direction; // arrow from center towards object
	// distance vector between center of circle and the car = Positions of points from the car vs Center of Circle, which is seen as 0
	double* dist_car_center; // 3 x 9
	double* Velocity_vec;
	double* Mass_vec;

public:
	Circular(double* Pos, double Rad);
	virtual ~Circular();
	void get_Position(double* Pos);
	double get_Radius() const;
	void set_Radius(const double& Rad);
	void set_Position(const double* Pos);
	void get_centrifugal_force(double* Fr, double* v, double& m, double* p);
	virtual void get_Profile_force(Car<double>* Car1, double* F_vec, double* Normal_ext);
	virtual void get_Profile_torque(Car<double>* Car1, double* Torque);
	virtual void update_initial_condition(Car<double>* Car1);
	virtual void this_is_a_test_fn(std::string s);
};
