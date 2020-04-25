#pragma once

#include "BLAS.h"
#include "car.h"
#include "Constants.h"
#include "MetaDataBase.h"
#include "RoadProfile.h"

template <typename T>
class LoadModule {
private:

	Profile<T>* Active_Profile;
	Car<T>* Car_obj;

	

	// Auxiliary vectors
	T* Normal_ext = NULL; // Normal_force updated with external forces
	T* k_vec = NULL; // vector with springs' stiffnesses (copied from Car!!!)
	T* Normal_from_profile = NULL; // Normal_force computed inside get_Profile_force
	T* External_force = NULL;

	/**
	* \brief adds the effect of the internal forces on the torque
	* \param Torque current torque of the body [XYZ]
	* \param F_vec vector of all internal XY-forces in the Eulerian frame at each position [GC:XYZ,W1:XYZ,T1:XYZ, ...]
	*/
	void ComputeInternalTorque(T* Torque, T* F_vec) {

		// compute torque around X-axis
//		#pragma loop(ivdep)
		for (int i = 0; i < 2 * Constants::NUM_LEGS + 1; ++i) {
			Torque[0] += F_vec[i * Constants::DIM + 1] * Car_obj->distance_nickpol[i];
		}

		// compute torque around Y-axis
	//	#pragma loop(ivdep)
		for (int i = 0; i < 2 * Constants::NUM_LEGS + 1; ++i) {
			Torque[1] += F_vec[i * Constants::DIM + 0] * Car_obj->distance_nickpol[i];
		}
	}


public:
	LoadModule(Profile<T>* Profile_type, Car<T>* Car1) {
		Active_Profile = Profile_type;
		Car_obj = Car1;

		// auxiliary vectors
		Normal_ext = (T*)mkl_malloc(sizeof(T) * Constants::DIM, Constants::ALIGNMENT); // normal_force, with external forces
		k_vec = (T*)mkl_malloc(sizeof(T) * (Constants::VEC_DIM - 1), Constants::ALIGNMENT); // k_vec; Constants::VEC_DIM-1 = 8

		// read External_force
		External_force = (T*)mkl_calloc((Constants::VEC_DIM * Constants::DIM), sizeof(T), Constants::ALIGNMENT); // 3 * 9 
		//set_External_force();
		mkl<T>::copy(Constants::DIM, MetaDataBase::DataBase()->getBodyExternalForce(), 1, External_force, 1); // copy the center of mass position
		T* xml_start, * position_start;
		xml_start = MetaDataBase::DataBase()->getWheelExternalForceFrontLeft();
		position_start = External_force + 3;
		mkl<T>::copy(Constants::DIM, xml_start, 1, position_start, 1);
		// W2 = W_fr
		xml_start = MetaDataBase::DataBase()->getWheelExternalForceFrontRight();
		position_start += 6; // skip 3 for tyre
		mkl<T>::copy(Constants::DIM, xml_start, 1, position_start, 1);
		// W3 = W_rl
		xml_start = MetaDataBase::DataBase()->getWheelExternalForceRearLeft();
		position_start += 6; // skip 3 for tyre
		mkl<T>::copy(Constants::DIM, xml_start, 1, position_start, 1);
		// W2 = W_rr
		xml_start = MetaDataBase::DataBase()->getWheelExternalForceRearRight();
		position_start += 6; // skip 3 for tyre
		mkl<T>::copy(Constants::DIM, xml_start, 1, position_start, 1);

		// T1 = T_fl
		xml_start = MetaDataBase::DataBase()->getTyreExternalForceFrontLeft();
		position_start = External_force + 6; // skip 3 for center of mass and 3 for the wheel
		mkl<T>::copy(Constants::DIM, xml_start, 1, position_start, 1);
		// T2 = T_fr
		xml_start = MetaDataBase::DataBase()->getTyreExternalForceFrontRight();
		position_start += 6; // skip 3 for the wheel
		mkl<T>::copy(Constants::DIM, xml_start, 1, position_start, 1);
		// T3 = T_rl
		xml_start = MetaDataBase::DataBase()->getTyreExternalForceRearLeft();
		position_start += 6; // skip 3 for the wheel
		mkl<T>::copy(Constants::DIM, xml_start, 1, position_start, 1);
		// T4 = T_rr
		xml_start = MetaDataBase::DataBase()->getTyreExternalForceRearRight();
		position_start += 6; // skip 3 for the wheel
		mkl<T>::copy(Constants::DIM, xml_start, 1, position_start, 1);
	}
	~LoadModule() {
		mkl_free(Normal_ext);
		mkl_free(k_vec);
		mkl_free(External_force);
	}
	void set_Profile(Profile<T>* Profile_type) {
		// not called
		Active_Profile = Profile_type;
	}
	void get_Profile(Profile<T>* Profile_type) {
		// not called
		Profile_type = Active_Profile;
	}

	/*
	Calculates the internal forces under certain conditions
	called in:  ALE.solve, global_frame_solver
	\param time_t current simulation time
	\param Delta_x_vec current dx of the spring lengths (in Stefan's ordering)
	\return F_vec vector of forces in the car [GC:XYZ,W1:XYZ,T1:XYZ, ...]
	\return Normal_ext external forces acting on the whole car system [XYZ]
	*/
	void update_force(T time_t, T* F_vec, T* Delta_x_vec, T* Normal_ext) {
		mkl<T>::scal(Constants::DIM * Constants::VEC_DIM, 0.0, F_vec, 1);
		mkl<T>::scal(Constants::DIM, 0.0, Normal_ext, 1);
		mkl<T>::scal(2 * Constants::NUM_LEGS, 0.0, Delta_x_vec, 1);
		Active_Profile->get_Profile_force_ALE(Car_obj, F_vec, Normal_ext);
		/*
		Modify the profile to know where the ground is and apply normal force accordingly
		*/
		// F_Ti += -0.25 * N; [F[2], F[4], F[6], F[8]]
		for (auto i = 2; i < Constants::VEC_DIM; i += 2) {
			mkl<T>::axpy(Constants::DIM, -0.25, Normal_ext, 1, &F_vec[i * Constants::DIM], 1);
		}
		// =============== PAY ATTENTION to THIS ============================
		// N += external_force  /// this formulation is WRONG (!!!) if done before computing following steps, should be done at the end. external force doesn't necessarily have to create a normal force it can create acceleration, ex: when car flies
		vdAdd(Constants::DIM, External_force, Normal_ext, Normal_ext);
		// =============== PAY ATTENTION to THIS ============================

		// get stiffnesses vector k_vec
		Car_obj->get_k_vec(k_vec);

		for (auto i = 0; i < (Constants::VEC_DIM - 1); i += 2) {
			// use the elastic forces at wheels
			// F_CG += k_wi * delta_x_i
			F_vec[2] += 0.0 * k_vec[i] * Delta_x_vec[i];
			//mkl<T>::axpy(DIM, k_vec[i], &Delta_x_vec[DIM * i], 1, F_vec, 1);
			// F_W_i += -k_wi * delta_x_i
			//mkl<T>::axpy(DIM, -k_vec[i], &Delta_x_vec[DIM * i], 1, &F_vec[DIM * (i + 1)], 1);
			F_vec[Constants::DIM * (i + 1) + 2] -= 0.0 * k_vec[i] * Delta_x_vec[i];

			// use the elastic forces at tyres
			// F_W_i += k_t_i * delta_x_{i+1}
			//mkl<T>::axpy(DIM, k_vec[i + 1], &Delta_x_vec[DIM * (i + 1)], 1, &F_vec[DIM * (i + 1)], 1);
			F_vec[Constants::DIM * (i + 1) + 2] += 0.0 * k_vec[i + 1] * Delta_x_vec[(i + 1)];
			// F_T_i += -k_t_i * delta_x_{i+1}
			//mkl<T>::axpy(DIM, -k_vec[i + 1], &Delta_x_vec[DIM * (i + 1)], 1, &F_vec[DIM * (i + 2)], 1);
			F_vec[Constants::DIM * (i + 2) + 2] -= 0.0 * k_vec[i + 1] * Delta_x_vec[(i + 1)];
		}

		// test add
		mkl<T>::axpy(Constants::DIM, 1.0, External_force, 1, F_vec, 1);

	}


	/*
	Calculates the torque acting on the whole car
	called in:  ALE.solve, global_frame_solver
	\param time_t current simulation time
	\param Delta_x_vec current dx of the spring lengths (in Stefan's ordering)
	\return Torque acting on the total car system [XYZ]
	*/
	void update_torque(T time_t, T* Torque, T* F_vec) {
		Active_Profile->get_Profile_torque(Car_obj, Torque);
		ComputeInternalTorque(Torque, F_vec);
	}


	void set_External_force() {
		// in LoadModule constructor commented

		mkl<T>::copy(Constants::DIM, MetaDataBase::DataBase()->getBodyExternalForce(), 1, External_force, 1); // copy the center of mass position
		T* xml_start, * position_start;
		xml_start = MetaDataBase::DataBase()->getWheelExternalForceFrontLeft();
		position_start = External_force + 3;
		mkl<T>::copy(Constants::DIM, xml_start, 1, position_start, 1);
		// W2 = W_fr
		xml_start = MetaDataBase::DataBase()->getWheelExternalForceFrontRight();
		position_start += 6; // skip 3 for tyre
		mkl<T>::copy(Constants::DIM, xml_start, 1, position_start, 1);
		// W3 = W_rl
		xml_start = MetaDataBase::DataBase()->getWheelExternalForceRearLeft();
		position_start += 6; // skip 3 for tyre
		mkl<T>::copy(Constants::DIM, xml_start, 1, position_start, 1);
		// W2 = W_rr
		xml_start = MetaDataBase::DataBase()->getWheelExternalForceRearRight();
		position_start += 6; // skip 3 for tyre
		mkl<T>::copy(Constants::DIM, xml_start, 1, position_start, 1);

		// T1 = T_fl
		xml_start = MetaDataBase::DataBase()->getTyreExternalForceFrontLeft();
		position_start = External_force + 6; // skip 3 for center of mass and 3 for the wheel
		mkl<T>::copy(Constants::DIM, xml_start, 1, position_start, 1);
		// T2 = T_fr
		xml_start = MetaDataBase::DataBase()->getTyreExternalForceFrontRight();
		position_start += 6; // skip 3 for the wheel
		mkl<T>::copy(Constants::DIM, xml_start, 1, position_start, 1);
		// T3 = T_rl
		xml_start = MetaDataBase::DataBase()->getTyreExternalForceRearLeft();
		position_start += 6; // skip 3 for the wheel
		mkl<T>::copy(Constants::DIM, xml_start, 1, position_start, 1);
		// T4 = T_rr
		xml_start = MetaDataBase::DataBase()->getTyreExternalForceRearRight();
		position_start += 6; // skip 3 for the wheel
		mkl<T>::copy(Constants::DIM, xml_start, 1, position_start, 1);
	}
	
};

