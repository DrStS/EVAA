#pragma once

#include <iostream>
#include "Load_module.h"
#include "MathLibrary.h"
#include "ReadXML.h"

#ifdef use_intel_mkl
#include <mkl.h>
#define use_gemm
#endif


// ===============================  Load_module class implementation =======================

// in:  [f, t, w1, t1, w2, t2, w3, t3, w4, t4]
// out: [f, t, w1, t1, w2, t2, w3, t3, w4, t4]

Load_module::Load_module(Profile* Profile_type, Car<double>* Car1, Load_Params load_param) {
	Active_Profile = Profile_type;
	Car_obj = Car1;
	
	// auxiliary vectors
	Normal_ext = (double*)mkl_malloc(sizeof(double) * Car1->DIM, alignment); // normal_force, with external forces
	k_vec = (double*)mkl_malloc(sizeof(double) * (Car1->vec_DIM - 1), alignment); // k_vec; Car1->vec_DIM-1 = 8

	// read External_force
	External_force = (double*)mkl_calloc((Car1->vec_DIM * Car1->DIM), sizeof(double), alignment); // 3 * 9 
	//set_External_force(load_param);
	cblas_dcopy(Car1->DIM, load_param.external_force_body, 1, External_force, 1); // copy the center of mass position
	double* xml_start, * position_start;
	xml_start = load_param.external_force_wheel + 2 * 3;
	position_start = External_force + 3;
	cblas_dcopy(Car1->DIM, xml_start, 1, position_start, 1);
	// W2 = W_fr
	xml_start = load_param.external_force_wheel + 3 * 3;
	position_start += 6; // skip 3 for tyre
	cblas_dcopy(Car1->DIM, xml_start, 1, position_start, 1);
	// W3 = W_rl
	xml_start = load_param.external_force_wheel + 1 * 3;
	position_start += 6; // skip 3 for tyre
	cblas_dcopy(Car1->DIM, xml_start, 1, position_start, 1);
	// W2 = W_rr
	xml_start = load_param.external_force_wheel + 0 * 3;
	position_start += 6; // skip 3 for tyre
	cblas_dcopy(Car1->DIM, xml_start, 1, position_start, 1);

	// T1 = T_fl
	xml_start = load_param.external_force_tyre + 2 * 3;
	position_start = External_force + 6; // skip 3 for center of mass and 3 for the wheel
	cblas_dcopy(Car1->DIM, xml_start, 1, position_start, 1);
	// T2 = T_fr
	xml_start = load_param.external_force_tyre + 3 * 3;
	position_start += 6; // skip 3 for the wheel
	cblas_dcopy(Car1->DIM, xml_start, 1, position_start, 1);
	// T3 = T_rl
	xml_start = load_param.external_force_tyre + 1 * 3;
	position_start += 6; // skip 3 for the wheel
	cblas_dcopy(Car1->DIM, xml_start, 1, position_start, 1);
	// T4 = T_rr
	xml_start = load_param.external_force_tyre + 0 * 3;
	position_start += 6; // skip 3 for the wheel
	cblas_dcopy(Car1->DIM, xml_start, 1, position_start, 1);
}

Load_module::~Load_module() {
	std::cout << "Load module destructor";
	mkl_free(Normal_ext);
	mkl_free(k_vec);
	mkl_free(External_force);
}

void Load_module::set_External_force(Load_Params load_param) {
	cblas_dcopy(DIM, load_param.external_force_body, 1, External_force, 1); // copy the center of mass position
	double* xml_start, * position_start;
	xml_start = load_param.external_force_wheel + 2 * 3;
	position_start = External_force + 3;
	cblas_dcopy(DIM, xml_start, 1, position_start, 1);
	// W2 = W_fr
	xml_start = load_param.external_force_wheel + 3 * 3;
	position_start += 6; // skip 3 for tyre
	cblas_dcopy(DIM, xml_start, 1, position_start, 1);
	// W3 = W_rl
	xml_start = load_param.external_force_wheel + 1 * 3;
	position_start += 6; // skip 3 for tyre
	cblas_dcopy(DIM, xml_start, 1, position_start, 1);
	// W2 = W_rr
	xml_start = load_param.external_force_wheel + 0 * 3;
	position_start += 6; // skip 3 for tyre
	cblas_dcopy(DIM, xml_start, 1, position_start, 1);

	// T1 = T_fl
	xml_start = load_param.external_force_tyre + 2 * 3;
	position_start = External_force + 6; // skip 3 for center of mass and 3 for the wheel
	cblas_dcopy(DIM, xml_start, 1, position_start, 1);
	// T2 = T_fr
	xml_start = load_param.external_force_tyre + 3 * 3;
	position_start += 6; // skip 3 for the wheel
	cblas_dcopy(DIM, xml_start, 1, position_start, 1);
	// T3 = T_rl
	xml_start = load_param.external_force_tyre + 1 * 3;
	position_start += 6; // skip 3 for the wheel
	cblas_dcopy(DIM, xml_start, 1, position_start, 1);
	// T4 = T_rr
	xml_start = load_param.external_force_tyre + 0 * 3;
	position_start += 6; // skip 3 for the wheel
	cblas_dcopy(DIM, xml_start, 1, position_start, 1);
}

void Load_module::set_Profile(Profile* Profile_type) {
	Active_Profile = Profile_type;
}

void Load_module::get_Profile(Profile* Profile_type) {
	Profile_type = Active_Profile;
}

void Load_module::update_force(double time_t, double* F_vec, double* Delta_x_vec, double* Normal_ext) {
	Active_Profile->get_Profile_force(Car_obj, F_vec, Normal_ext);
	

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
		cblas_daxpy(DIM, k_vec[i], &Delta_x_vec[DIM * i], incx, F_vec, incx);
		// F_W_i += -k_wi * delta_x_i
		cblas_daxpy(DIM, -k_vec[i], &Delta_x_vec[DIM * i], incx, &F_vec[DIM * (i + 1)], incx);

		// use the elastic forces at tyres
		// F_W_i += -k_t_i * delta_x_{i+1}
		cblas_daxpy(DIM, -k_vec[i + 1], &Delta_x_vec[DIM * (i + 1)], incx, &F_vec[DIM * (i + 1)], incx);
		// F_T_i += k_t_i * delta_x_{i+1}
		cblas_daxpy(DIM, k_vec[i + 1], &Delta_x_vec[DIM * (i + 1)], incx, &F_vec[DIM * (i + 2)], incx);
	}
	// test add
	
}

void Load_module::update_torque(double time_t, double* Torque, double* Delta_x_vec) {
	Active_Profile->get_Profile_torque(Car_obj, Torque);
}

// =============================== end of Load_module class implementation =======================