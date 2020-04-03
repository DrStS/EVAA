#pragma once

// classes: Profile, Circular, Car, Load_module etc.

class Car {
private:
	// consider replacing double with floatEVAA
	// alignment and g are defined globally === REPLACE THEM!!!!
	const int alignment = 64;
	// MKL / vector constants:
	const int mkl_DIM = 3, vec_DIM = 10, incx = 1; // consider mkl_DIM = 4 for efficiency!!! // 10 dimension because of torque of the body

	double* Position_vec = NULL; // [CG, W1, T1, W2, T2, W3, T3, W4, T4] 9 x 3 !!! Consider alignment (3+1),(3+1),... 
	double* Velocity_vec = NULL; // [CG, W1, T1, W2, T2, W3, T3, W4, T4] 9 x 3
	double* Mass_vec = NULL; // [CG, W1, T1, W2, T2, W3, T3, W4, T4]     9

	// Length and Width of the car
	double Length;
	double Width;

	// Moments of inertia
	double I_body_xx = 640;
	double I_body_yy = 4800;

	// Springs
	// Stiffnesses
	double* k_vec = NULL;
	// Length of springs (?)
	double L_body_fl;
	double L_tire_fl;
	double L_body_fr;
	double L_tire_fr;
	double L_body_rl;
	double L_tire_rl;
	double L_body_rr;
	double L_tire_rr;

	// Distances from CG to wheels
	double l_long_fl;
	double l_long_fr;
	double l_long_rl;
	double l_long_rr;
	double l_lat_fl;
	double l_lat_fr;
	double l_lat_rl;
	double l_lat_rr;

public:
	Car();
	Car(double* Pos, double* Vel);
	Car(double* Pos, double* Vel, double* M);
	Car(double* Pos, double* Vel, double* M, double* k);
	Car(double* Pos, double* Vel, double* M, double* k,
		double Len, double Wid,
		double L_body_fl_val, double L_tire_fl_val,
		double L_body_fr_val, double L_tire_fr_val,
		double L_body_rl_val, double L_tire_rl_val,
		double L_body_rr_val, double L_tire_rr_val);
	Car(double* Pos, double* Vel, double* M, double* k,
		double Len, double Wid,
		double L_body_fl_val, double L_tire_fl_val,
		double L_body_fr_val, double L_tire_fr_val,
		double L_body_rl_val, double L_tire_rl_val,
		double L_body_rr_val, double L_tire_rr_val,
		double l_long_fl_val, double l_long_fr_val,
		double l_long_rl_val, double l_long_rr_val,
		double l_lat_fl_val, double l_lat_fr_val,
		double l_lat_rl_val, double l_lat_rr_val);
	Car(double* Pos, double* Vel, double* M, double* k,
		double Len, double Wid,
		double L_body_fl_val, double L_tire_fl_val,
		double L_body_fr_val, double L_tire_fr_val,
		double L_body_rl_val, double L_tire_rl_val,
		double L_body_rr_val, double L_tire_rr_val,
		double l_long_fl_val, double l_long_fr_val,
		double l_long_rl_val, double l_long_rr_val,
		double l_lat_fl_val, double l_lat_fr_val,
		double l_lat_rl_val, double l_lat_rr_val,
		double 	I_body_xx_val, double I_body_yy_val);
	Car(const Car& Car1);
	Car& operator= (const Car& Car1);
	~Car();
	void get_Position_vec(double* Pos);
	void set_Position_vec(const double* Pos);
	void get_Position_vec_CG(double* Pos);
	void set_Position_vec_CG(const double* Pos);

	void get_Velocity_vec(double* Vel);
	void set_Velocity_vec(const double* Vel);
	void get_Velocity_vec_CG(double* Vel);
	void set_Velocity_vec_CG(const double* Vel);

	void get_k_vec(double* k_vec);
	void set_k_vec(const double* k_vec);

	void get_Mass_vec(double* M);
	double get_Mass_vec_CG() const;
	void set_Mass_vec(const double* M);

	void get_dist_vector(double* Point_P, double* dist_vector); // 9 * 3 - from each important point to a fixed Point_P
	void get_dist_vector_CG(double* Point_P, double* dist_vector); // 3 - from Center of Gravity of a fixed Point_P

	double get_I_body_xx() const;
	void set_I_body_xx(const double& I_body_xx_val);
	double get_I_body_yy() const;
	void set_I_body_yy(const double& I_body_yy_val);
};


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

