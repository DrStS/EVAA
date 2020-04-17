#pragma once

#include <iostream>
#include "Load_module.h"
#include "MathLibrary.h"
#include "Constants.h"

#ifdef use_intel_mkl
#include <mkl.h>
#define use_gemm
#endif


// ===============================  Load_module class implementation =======================

// in:  [f, t, w1, t1, w2, t2, w3, t3, w4, t4]
// out: [f, t, w1, t1, w2, t2, w3, t3, w4, t4]

Load_module::Load_module(Profile* Profile_type, Car<double>* Car1) {
	Active_Profile = Profile_type;
	Car_obj = Car1;
	
	// auxiliary vectors
	Normal_ext = (double*)mkl_malloc(sizeof(double) * Constants::DIM, Constants::ALIGNMENT); // normal_force, with external forces
	k_vec = (double*)mkl_malloc(sizeof(double) * (Constants::VEC_DIM - 1), Constants::ALIGNMENT); // k_vec; Constants::VEC_DIM-1 = 8

	// read External_force
	External_force = (double*)mkl_calloc((Constants::VEC_DIM * Constants::DIM), sizeof(double), Constants::ALIGNMENT); // 3 * 9 
	//set_External_force();
	cblas_dcopy(Constants::DIM, MetaDataBase::DataBase()->getBodyExternalForce(), 1, External_force, 1); // copy the center of mass position
	double* xml_start, * position_start;
	xml_start = MetaDataBase::DataBase()->getWheelExternalForceFrontLeft();
	position_start = External_force + 3;
	cblas_dcopy(Constants::DIM, xml_start, 1, position_start, 1);
	// W2 = W_fr
	xml_start = MetaDataBase::DataBase()->getWheelExternalForceFrontRight();
	position_start += 6; // skip 3 for tyre
	cblas_dcopy(Constants::DIM, xml_start, 1, position_start, 1);
	// W3 = W_rl
	xml_start = MetaDataBase::DataBase()->getWheelExternalForceRearLeft();
	position_start += 6; // skip 3 for tyre
	cblas_dcopy(Constants::DIM, xml_start, 1, position_start, 1);
	// W2 = W_rr
	xml_start = MetaDataBase::DataBase()->getWheelExternalForceRearRight();
	position_start += 6; // skip 3 for tyre
	cblas_dcopy(Constants::DIM, xml_start, 1, position_start, 1);

	// T1 = T_fl
	xml_start = MetaDataBase::DataBase()->getTyreExternalForceFrontLeft();
	position_start = External_force + 6; // skip 3 for center of mass and 3 for the wheel
	cblas_dcopy(Constants::DIM, xml_start, 1, position_start, 1);
	// T2 = T_fr
	xml_start = MetaDataBase::DataBase()->getTyreExternalForceFrontRight();
	position_start += 6; // skip 3 for the wheel
	cblas_dcopy(Constants::DIM, xml_start, 1, position_start, 1);
	// T3 = T_rl
	xml_start = MetaDataBase::DataBase()->getTyreExternalForceRearLeft();
	position_start += 6; // skip 3 for the wheel
	cblas_dcopy(Constants::DIM, xml_start, 1, position_start, 1);
	// T4 = T_rr
	xml_start = MetaDataBase::DataBase()->getTyreExternalForceRearRight();
	position_start += 6; // skip 3 for the wheel
	cblas_dcopy(Constants::DIM, xml_start, 1, position_start, 1);
}

Load_module::~Load_module() {
	mkl_free(Normal_ext);
	mkl_free(k_vec);
	mkl_free(External_force);
}

void Load_module::set_External_force() {
	cblas_dcopy(DIM, MetaDataBase::DataBase()->getBodyExternalForce(), 1, External_force, 1); // copy the center of mass position
	double* xml_start, * position_start;
	xml_start = MetaDataBase::DataBase()->getWheelExternalForceFrontLeft();
	position_start = External_force + 3;
	cblas_dcopy(Constants::DIM, xml_start, 1, position_start, 1);
	// W2 = W_fr
	xml_start = MetaDataBase::DataBase()->getWheelExternalForceFrontRight();
	position_start += 6; // skip 3 for tyre
	cblas_dcopy(Constants::DIM, xml_start, 1, position_start, 1);
	// W3 = W_rl
	xml_start = MetaDataBase::DataBase()->getWheelExternalForceRearLeft();
	position_start += 6; // skip 3 for tyre
	cblas_dcopy(Constants::DIM, xml_start, 1, position_start, 1);
	// W2 = W_rr
	xml_start = MetaDataBase::DataBase()->getWheelExternalForceRearRight();
	position_start += 6; // skip 3 for tyre
	cblas_dcopy(Constants::DIM, xml_start, 1, position_start, 1);

	// T1 = T_fl
	xml_start = MetaDataBase::DataBase()->getTyreExternalForceFrontLeft();
	position_start = External_force + 6; // skip 3 for center of mass and 3 for the wheel
	cblas_dcopy(Constants::DIM, xml_start, 1, position_start, 1);
	// T2 = T_fr
	xml_start = MetaDataBase::DataBase()->getTyreExternalForceFrontRight();
	position_start += 6; // skip 3 for the wheel
	cblas_dcopy(Constants::DIM, xml_start, 1, position_start, 1);
	// T3 = T_rl
	xml_start = MetaDataBase::DataBase()->getTyreExternalForceRearLeft();
	position_start += 6; // skip 3 for the wheel
	cblas_dcopy(Constants::DIM, xml_start, 1, position_start, 1);
	// T4 = T_rr
	xml_start = MetaDataBase::DataBase()->getTyreExternalForceRearRight();
	position_start += 6; // skip 3 for the wheel
	cblas_dcopy(Constants::DIM, xml_start, 1, position_start, 1);
}

void Load_module::set_Profile(Profile* Profile_type) {
	Active_Profile = Profile_type;
}

void Load_module::get_Profile(Profile* Profile_type) {
	Profile_type = Active_Profile;
}

void Load_module::update_force(double time_t, double* F_vec, double* Delta_x_vec, double* Normal_ext) {
	cblas_dscal(Constants::DIM*Constants::VEC_DIM, 0.0, F_vec, 1);
	cblas_dscal(Constants::DIM, 0.0, Normal_ext, 1);
	cblas_dscal(2 * Constants::NUM_LEGS, 0.0, Delta_x_vec, 1);
	Active_Profile->get_Profile_force_ALE(Car_obj, F_vec, Normal_ext);
	/*
	Modify the profile to know where the ground is and apply normal force accordingly
	*/
	// F_Ti += -0.25 * N; [F[2], F[4], F[6], F[8]]
	for (auto i = 2; i < vec_DIM; i += 2) {
		cblas_daxpy(DIM, -0.25, Normal_ext, incx, &F_vec[i * DIM], incx); 
	}
	// =============== PAY ATTENTION to THIS ============================
	// N += external_force  /// this formulation is WRONG (!!!) if done before computing following steps, should be done at the end. external force doesn't necessarily have to create a normal force it can create acceleration, ex: when car flies
	vdAdd(DIM, External_force, Normal_ext, Normal_ext);
	// =============== PAY ATTENTION to THIS ============================

	// get stiffnesses vector k_vec
	Car_obj->get_k_vec(k_vec);

	for (auto i = 0; i < (vec_DIM - 1); i += 2) {
		// use the elastic forces at wheels
		// F_CG += k_wi * delta_x_i
		F_vec[2] += 0.0*k_vec[i] * Delta_x_vec[i];
		//cblas_daxpy(DIM, k_vec[i], &Delta_x_vec[DIM * i], incx, F_vec, incx);
		// F_W_i += -k_wi * delta_x_i
		//cblas_daxpy(DIM, -k_vec[i], &Delta_x_vec[DIM * i], incx, &F_vec[DIM * (i + 1)], incx);
		F_vec[DIM * (i + 1) + 2] -= 0.0*k_vec[i] * Delta_x_vec[i];

		// use the elastic forces at tyres
		// F_W_i += k_t_i * delta_x_{i+1}
		//cblas_daxpy(DIM, k_vec[i + 1], &Delta_x_vec[DIM * (i + 1)], incx, &F_vec[DIM * (i + 1)], incx);
		F_vec[DIM * (i + 1) + 2] += 0.0*k_vec[i + 1] * Delta_x_vec[(i + 1)];
		// F_T_i += -k_t_i * delta_x_{i+1}
		//cblas_daxpy(DIM, -k_vec[i + 1], &Delta_x_vec[DIM * (i + 1)], incx, &F_vec[DIM * (i + 2)], incx);
		F_vec[DIM * (i + 2) + 2] -= 0.0*k_vec[i + 1] * Delta_x_vec[(i + 1)];
	}
	// test add
	cblas_daxpy(DIM, 1.0, External_force, incx, F_vec, incx);

}

void Load_module::update_torque(double time_t, double* Torque, double* Delta_x_vec) {
	Active_Profile->get_Profile_torque(Car_obj, Torque);
}

// =============================== end of Load_module class implementation =======================