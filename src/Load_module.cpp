#pragma once

#include <iostream>
#include "load_module.h"
#include "mathlibrary.h"

#ifdef use_intel_mkl
#include <mkl.h>
#define use_gemm
#endif


// ===============================  load_module class implementation =======================

// in:  [f, t, w1, t1, w2, t2, w3, t3, w4, t4]
// out: [f, t, w1, t1, w2, t2, w3, t3, w4, t4]

load_module::load_module() {
	active_profile = new profile();
	car_obj = new car();

	// auxiliary vectors
	normal_ext = (double*)mkl_malloc(sizeof(double) * mkl_dim, alignment); // normal_force
	k_vec = (double*)mkl_malloc(sizeof(double) * (vec_dim - 1), alignment); // k_vec; vec_dim-1 = 8
}

load_module::load_module(profile* profile_type) {
	active_profile = profile_type;
	car_obj = new car();

	// auxiliary vectors
	normal_ext = (double*)mkl_malloc(sizeof(double) * mkl_dim, alignment); // normal_force
	k_vec = (double*)mkl_malloc(sizeof(double) * (vec_dim - 1), alignment); // k_vec; vec_dim-1 = 8
}

load_module::load_module(profile* profile_type, car* car1) {
	active_profile = profile_type;
	car_obj = car1;

	// auxiliary vectors
	normal_ext = (double*)mkl_malloc(sizeof(double) * mkl_dim, alignment); // normal_force
	k_vec = (double*)mkl_malloc(sizeof(double) * (vec_dim - 1), alignment); // k_vec; vec_dim-1 = 8
}

load_module::~load_module() {
	if (normal_ext != null) {
		mkl_free(normal_ext);
		normal_ext = null;
	}

	if (k_vec != null) {
		mkl_free(k_vec);
		k_vec = null;
	}
}

load_module::load_module(const load_module& load_module_1) {
	*this = load_module_1;
};

load_module& load_module::operator= (const load_module& load_module_1) {
	active_profile = load_module_1.active_profile;
	car_obj = load_module_1.car_obj;

	return *this;
}

void load_module::set_profile(profile* profile_type) {
	active_profile = profile_type;
}

void load_module::get_profile(profile* profile_type) {
	profile_type = active_profile;
}

void load_module::update_force(double time_t, double* f_vec, double* delta_x_vec, double* external_force) {
	active_profile->update_profile_force(car_obj, f_vec, normal_ext);

	// n += external_force
	if (external_force != null) {
		vdadd(mkl_dim, external_force, normal_ext, normal_ext);
	}
	// f_ti += -0.25 * n; [f[2], f[4], f[6], f[8]]
	for (auto i = 2; i < vec_dim; i += 2) {
		cblas_daxpy(mkl_dim, -0.25, normal_ext, incx, &f_vec[i * mkl_dim], incx);
	}
	
	// get stiffnesses vector k_vec
	car_obj->get_k_vec(k_vec);

	for (auto i = 0; i < (vec_dim - 1); i += 2) {
		// use the elastic forces at wheels
		// f_cg += k_wi * delta_x_i
		cblas_daxpy(mkl_dim, k_vec[i], &delta_x_vec[mkl_dim * i], incx, f_vec, incx);
		// f_w_i -= k_wi * delta_x_i
		cblas_daxpy(mkl_dim, -k_vec[i], &delta_x_vec[mkl_dim * i], incx, &f_vec[mkl_dim * (i + 1)], incx);

		// use the elastic forces at tyres
		// f_w_i -= k_t_i * delta_x_{i+1}
		cblas_daxpy(mkl_dim, -k_vec[i + 1], &delta_x_vec[mkl_dim * (i + 1)], incx, &f_vec[mkl_dim * (i + 1)], incx);
		// f_t_i += k_t_i * delta_x_{i+1}
		cblas_daxpy(mkl_dim, k_vec[i + 1], &delta_x_vec[mkl_dim * (i + 1)], incx, &f_vec[mkl_dim * (i + 2)], incx);
	}
}

// read external force
/*{
cblas_dcopy(mkl_dim, load_param.external_force_body, 1, force_vector, 1); // copy the center of mass position
t* xml_start, * position_start;
xml_start = load_param.external_force_wheel + 2 * 3;
position_start = force_vector + 3;
cblas_dcopy(mkl_dim, xml_start, 1, position_start, 1);
// w2 = w_fr
xml_start = load_param.external_force_wheel + 3 * 3;
position_start += 6; // skip 3 for tyre
cblas_dcopy(mkl_dim, xml_start, 1, position_start, 1);
// w3 = w_rl
xml_start = load_param.external_force_wheel + 1 * 3;
position_start += 6; // skip 3 for tyre
cblas_dcopy(mkl_dim, xml_start, 1, position_start, 1);
// w2 = w_rr
xml_start = load_param.external_force_wheel + 0 * 3;
position_start += 6; // skip 3 for tyre
cblas_dcopy(mkl_dim, xml_start, 1, position_start, 1);

// t1 = t_fl
xml_start = load_param.external_force_tyre + 2 * 3;
position_start = force_vector + 6; // skip 3 for center of mass and 3 for the wheel
cblas_dcopy(mkl_dim, xml_start, 1, position_start, 1);
// t2 = t_fr
xml_start = load_param.external_force_tyre + 3 * 3;
position_start += 6; // skip 3 for the wheel
cblas_dcopy(mkl_dim, xml_start, 1, position_start, 1);
// t3 = t_rl
xml_start = load_param.external_force_tyre + 1 * 3;
position_start += 6; // skip 3 for the wheel
cblas_dcopy(mkl_dim, xml_start, 1, position_start, 1);
// t4 = t_rr
xml_start = load_param.external_force_tyre + 0 * 3;
position_start += 6; // skip 3 for the wheel
cblas_dcopy(mkl_dim, xml_start, 1, position_start, 1);
}*/

// =============================== end of load_module class implementation =======================