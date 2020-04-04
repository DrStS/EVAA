#pragma once

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

	// vectors used in computations inside the methods
	double* unit_y_vector = NULL;
	double* velocity_direction = NULL; // arrow from center towards object
	// distance vector between center of circle and the car = Positions of points from the car vs Center of Circle, which is seen as 0
	double* dist_car_center = NULL; // 3 x 9
	double* Velocity_vec = NULL;
	double* Mass_vec = NULL;

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