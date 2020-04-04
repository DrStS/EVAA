#pragma once

#include <iostream>
#include "Load_module.h"
#include "mathlibrary.h"

#ifdef use_intel_mkl
#include <mkl.h>
#define use_gemm
#endif


// ===============================  Load_module class implementation =======================

// in:  [f, t, w1, t1, w2, t2, w3, t3, w4, t4]
// out: [f, t, w1, t1, w2, t2, w3, t3, w4, t4]

Load_module::Load_module(Profile* profile_type, Car<double>* car1) {
	Active_Profile = profile_type;
	Car_obj = car1;

	// auxiliary vectors
	Normal_ext = (double*)mkl_malloc(sizeof(double) * mkl_DIM, alignment); // normal_force
	k_vec = (double*)mkl_malloc(sizeof(double) * (vec_DIM - 1), alignment); // k_vec; vec_DIM-1 = 8
}

Load_module::~Load_module() {
	if (Normal_ext != NULL) {
		mkl_free(Normal_ext);
		Normal_ext = NULL;
	}

	if (k_vec != NULL) {
		mkl_free(k_vec);
		k_vec = NULL;
	}
}

Load_module::Load_module(const Load_module& load_module_1) {
	*this = load_module_1;
};

Load_module& Load_module::operator= (const Load_module& load_module_1) {
	Active_Profile = load_module_1.Active_Profile;
	Car_obj = load_module_1.Car_obj;

	return *this;
}

void Load_module::set_Profile(Profile* profile_type) {
	Active_Profile = profile_type;
}

void Load_module::get_Profile(Profile* profile_type) {
	profile_type = Active_Profile;
}

void Load_module::update_force(double time_t, double* f_vec, double* delta_x_vec, double* external_force) {
	Active_Profile->update_Profile_force(Car_obj, f_vec, Normal_ext);

	// n += external_force
	if (external_force != NULL) {
		vdAdd(mkl_DIM, external_force, Normal_ext, Normal_ext);
	}
	// f_ti += -0.25 * n; [f[2], f[4], f[6], f[8]]
	for (auto i = 2; i < vec_DIM; i += 2) {
		cblas_daxpy(mkl_DIM, -0.25, Normal_ext, incx, &f_vec[i * mkl_DIM], incx);
	}
	
	// get stiffnesses vector k_vec
	Car_obj->get_k_vec(k_vec);

	for (auto i = 0; i < (vec_DIM - 1); i += 2) {
		// use the elastic forces at wheels
		// f_cg += k_wi * delta_x_i
		cblas_daxpy(mkl_DIM, k_vec[i], &delta_x_vec[mkl_DIM * i], incx, f_vec, incx);
		// f_w_i -= k_wi * delta_x_i
		cblas_daxpy(mkl_DIM, -k_vec[i], &delta_x_vec[mkl_DIM * i], incx, &f_vec[mkl_DIM * (i + 1)], incx);

		// use the elastic forces at tyres
		// f_w_i -= k_t_i * delta_x_{i+1}
		cblas_daxpy(mkl_DIM, -k_vec[i + 1], &delta_x_vec[mkl_DIM * (i + 1)], incx, &f_vec[mkl_DIM * (i + 1)], incx);
		// f_t_i += k_t_i * delta_x_{i+1}
		cblas_daxpy(mkl_DIM, k_vec[i + 1], &delta_x_vec[mkl_DIM * (i + 1)], incx, &f_vec[mkl_DIM * (i + 2)], incx);
	}
}

// read external force
/*{
cblas_dcopy(mkl_DIM, load_param.external_force_body, 1, force_vector, 1); // copy the center of mass position
t* xml_start, * position_start;
xml_start = load_param.external_force_wheel + 2 * 3;
position_start = force_vector + 3;
cblas_dcopy(mkl_DIM, xml_start, 1, position_start, 1);
// w2 = w_fr
xml_start = load_param.external_force_wheel + 3 * 3;
position_start += 6; // skip 3 for tyre
cblas_dcopy(mkl_DIM, xml_start, 1, position_start, 1);
// w3 = w_rl
xml_start = load_param.external_force_wheel + 1 * 3;
position_start += 6; // skip 3 for tyre
cblas_dcopy(mkl_DIM, xml_start, 1, position_start, 1);
// w2 = w_rr
xml_start = load_param.external_force_wheel + 0 * 3;
position_start += 6; // skip 3 for tyre
cblas_dcopy(mkl_DIM, xml_start, 1, position_start, 1);

// t1 = t_fl
xml_start = load_param.external_force_tyre + 2 * 3;
position_start = force_vector + 6; // skip 3 for center of mass and 3 for the wheel
cblas_dcopy(mkl_DIM, xml_start, 1, position_start, 1);
// t2 = t_fr
xml_start = load_param.external_force_tyre + 3 * 3;
position_start += 6; // skip 3 for the wheel
cblas_dcopy(mkl_DIM, xml_start, 1, position_start, 1);
// t3 = t_rl
xml_start = load_param.external_force_tyre + 1 * 3;
position_start += 6; // skip 3 for the wheel
cblas_dcopy(mkl_DIM, xml_start, 1, position_start, 1);
// t4 = t_rr
xml_start = load_param.external_force_tyre + 0 * 3;
position_start += 6; // skip 3 for the wheel
cblas_dcopy(mkl_DIM, xml_start, 1, position_start, 1);
}*/

// =============================== end of Load_module class implementation =======================