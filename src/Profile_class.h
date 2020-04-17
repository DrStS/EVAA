#pragma once
#include <mkl.h>
#include <string>
#include "MathLibrary.h"
#include "car.h"
/*
Parent class for generic road profiles
*/
class Profile {
protected:
	// consider replacing double with floatEVAA
	// alignment and g are defined globally === REPLACE THEM!!!!
	const int alignment = 64, DIM = 3, vec_DIM = 9, incx = 1;

	std::string Name;

public:
	/*
	Get external force acting on the car system
	\param Car
	\return F_vec forces acting on each component [GC:XYZ,W1:XYZ,T1:XYZ, ...]
	\return Normal_ext on the car as a whole [XYZ]
	*/
	virtual void get_Profile_force(Car<double>* Car1, double* F_vec, double* Normal_ext) {};
	virtual void get_Profile_force_ALE(Car<double>* Car1, double* F_vec, double* Normal_ext) {};

	/*
	Get external torque acting on the car system
	\param Car
	\return F_vec torque acting on teh car system [XYZ]
	*/
	virtual void get_Profile_torque(Car<double>* Car1, double* Torque_vec) {};

	/*
	In case the initial conditions have to be specific to a road profile
	\param Car
	*/
	virtual void update_initial_condition(Car<double>* Car1) { std::cout << "No update of initial conditions" << std::endl; };
	virtual void set_fixed_index(size_t* index) {};

	virtual ~Profile() {};
};


/*
Follow a circular road profile with radius R and center C
*/
class Circular : public Profile {
private:
	// Center of the Circle (if the motion is a Circle)
	double* Position = NULL; 
	
	// Radius of the Circle (if the motion is a Circle)
	double Radius = 0; 

	// vectors used in computations inside the methods
	double* unit_y_vector;
	double* velocity_direction; // arrow from center towards object

	// distance vector between center of circle and the car = 
	//		Positions of points from the car vs Center of Circle, seen as 0
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
	void get_centrifugal_force_ALE(double* Fr, double* v, double& m, double* p);

	/*
	Get external force acting on the car system
	\param Car
	\return F_vec forces acting on each component [GC:XYZ,W1:XYZ,T1:XYZ, ...] centripetals on each component
	\return Normal_ext on the car as a whole [XYZ]
	*/
	virtual void get_Profile_force(Car<double>* Car1, double* F_vec, double* Normal_ext);
	virtual void get_Profile_force_ALE(Car<double>* Car1, double* F_vec, double* Normal_ext);

	/*
	Get external torque acting on the car system
	\param Car
	\return F_vec torque acting on teh car system [XYZ]
	*/
	virtual void get_Profile_torque(Car<double>* Car1, double* Torque);

	/*
	overwrites all initial velocities but the one from the main car body
	Calculates the initial angular velocity such that the car perfectly rotates around its own axis as it follows the circle
	\param Car
	*/
	virtual void update_initial_condition(Car<double>* Car1);
};

class Fixed : public Profile {
public:
	Fixed(const double& g);
	virtual void get_Profile_force_ALE(Car<double>* Car1, double* F_vec, double* Normal_ext);
	virtual void get_Profile_torque(Car<double>* Car1, double* Torque_vec);
	virtual void set_fixed_index(size_t* index);
	~Fixed();
private:
	size_t* linear_idx;
	double *k_vec, *dx;
	double* external_force;
	size_t num_tyre = 4;
	bool index_set = 0;
	double gravity;
};

/*
Follow a circular road profile with radius R and center C
*/
class Nonfixed : public Profile {
private:

public:
	Nonfixed(double* Pos, double Rad);
	virtual ~Nonfixed();
	void get_centrifugal_force(double* Fr, double* v, double& m, double* p);

	/*
	Get external force acting on the car system
	\param Car
	\return F_vec forces acting on each component [GC:XYZ,W1:XYZ,T1:XYZ, ...]
	\return Normal_ext on the car as a whole [XYZ]
	*/
	virtual void get_Profile_force(Car<double>* Car1, double* F_vec, double* Normal_ext);

	/*
	Get external torque acting on the car system
	\param Car
	\return F_vec torque acting on teh car system [XYZ]
	*/
	virtual void get_Profile_torque(Car<double>* Car1, double* Torque);

	/*
	overwrites all initial velocities but the one from the main car body
	Calculates the initial angular velocity such that the car perfectly rotates around its own axis as it follows the circle
	\param Car
	*/
	virtual void update_initial_condition(Car<double>* Car1);
};

