
// TODO: Copyright header

#pragma once

#include "Car.h"
#include "Constants.h"
#include "Math.h"
#include "MetaDataBase.h"
#include "RoadProfile.h"

namespace EVAA {

template <typename T>
class LoadModule {
private:
    Profile<T>* activeProfile;
	Lagrange<T>* lagrangianProfile;
	Euler<T>* eulerianProfile;
    Car<T>* carObj;
	bool aleCounter = false;

	void (LoadModule<T>::*_forceComputationModule)(T*, T*);
	void (LoadModule<T>::*_torqueComputationModule)(T*, T*);

    // Auxiliary vectors
	T* externalForce = NULL;
	T* reactionForce = NULL;
	T* _lagrangianForce = NULL;
	T* _eulerianForce = NULL;
	

    /**
     * \brief adds the effect of the internal forces on the torque
     * \param Torque current torque of the body [XYZ]
     * \param F_vec vector of all internal XY-forces in the Eulerian frame at each position
     * [GC:XYZ,W1:XYZ,T1:XYZ, ...]
     */
    void ComputeInternalTorqueFromLagrangianForce(T* Torque, T* F_vec) {
        // compute torque around X-axis
        // #pragma loop(ivdep)
		for (int i = 0; i < 2 * Constants::NUM_LEGS + 1; ++i) {
            Torque[0] += F_vec[i * (Constants::DIM - 1) + 1] * carObj->currentCIRTwoTrackModel[i];
        }

        // compute torque around Y-axis
        // #pragma loop(ivdep)
		for (int i = 0; i < 2 * Constants::NUM_LEGS + 1; ++i) {
            Torque[1] += F_vec[i * (Constants::DIM - 1) + 0] * carObj->currentCIRTwoTrackModel[i];
        }
    }
	void ReadExternalForce() {
		auto& db = MetaDataBase<T>::getDataBase();
		// copy the center of mass position
		Math::copy<T>(Constants::DIM, db.getBodyExternalForce(), 1, externalForce, 1);

		T *xml_start, *position_start;
		xml_start = db.getWheelExternalForceFrontLeft();
		position_start = externalForce + 3;
		Math::copy<T>(Constants::DIM, xml_start, 1, position_start, 1);
		// W2 = W_fr
		xml_start = db.getWheelExternalForceFrontRight();
		position_start += 6;  // skip 3 for tyre
		Math::copy<T>(Constants::DIM, xml_start, 1, position_start, 1);
		// W3 = W_rl
		xml_start = db.getWheelExternalForceRearLeft();
		position_start += 6;  // skip 3 for tyre
		Math::copy<T>(Constants::DIM, xml_start, 1, position_start, 1);
		// W2 = W_rr
		xml_start = db.getWheelExternalForceRearRight();
		position_start += 6;  // skip 3 for tyre
		Math::copy<T>(Constants::DIM, xml_start, 1, position_start, 1);

		// T1 = T_fl
		xml_start = db.getTyreExternalForceFrontLeft();
		position_start = externalForce + 6;  // skip 3 for center of mass and 3 for the wheel
		Math::copy<T>(Constants::DIM, xml_start, 1, position_start, 1);
		// T2 = T_fr
		xml_start = db.getTyreExternalForceFrontRight();
		position_start += 6;  // skip 3 for the wheel
		Math::copy<T>(Constants::DIM, xml_start, 1, position_start, 1);
		// T3 = T_rl
		xml_start = db.getTyreExternalForceRearLeft();
		position_start += 6;  // skip 3 for the wheel
		Math::copy<T>(Constants::DIM, xml_start, 1, position_start, 1);
		// T4 = T_rr
		xml_start = db.getTyreExternalForceRearRight();
		position_start += 6;  // skip 3 for the wheel
		Math::copy<T>(Constants::DIM, xml_start, 1, position_start, 1);
	}

public:
    LoadModule(Profile<T>* Profile_type, Car<T>* Car1) {
		throw "Not Implemented";
        activeProfile = Profile_type;
        carObj = Car1;

        auto& db = MetaDataBase<T>::getDataBase();

        // read External_force
		externalForce = Math::calloc<T>(Constants::VEC_DIM * Constants::DIM);
		reactionForce = Math::malloc<T>(Constants::DIM);
        // set_External_force();
		ReadExternalForce();
    }

	LoadModule(Lagrange<T>* lagrangeProfile, Euler<T>* eulerProfile, Car<T>* Car1) {
		lagrangianProfile = lagrangeProfile;
		eulerianProfile = eulerProfile;
		carObj = Car1;
		aleCounter = true;

		// read External_force
		externalForce = Math::calloc<T>(Constants::VEC_DIM * Constants::DIM);
		reactionForce = Math::malloc<T>(Constants::DIM);
		_lagrangianForce = Math::calloc<T>(Constants::VEC_DIM * (Constants::DIM - 1));
		_eulerianForce = Math::calloc<T>(Constants::DOF);
		// set_External_force();
		ReadExternalForce();

	}
    ~LoadModule() {
        Math::free<T>(externalForce);
		Math::free<T>(reactionForce);
		if (aleCounter) {
			Math::free<T>(_lagrangianForce);
			Math::free<T>(_eulerianForce);
		}
    }

	void GetLagrangianForce(T _time, T* lagrangeForce) {
		lagrangianProfile->GetProfileForceLagrangian(carObj, _lagrangianForce, lagrangeForce);

		// add external force
		for (auto i = 0; i < Constants::VEC_DIM; ++i) {
			lagrangeForce[0] += externalForce[i*Constants::DIM];
			lagrangeForce[1] += externalForce[i*Constants::DIM + 1];
		}

	}

	void GetEulerianForce(T _time, T* eulerForce) {
		// compute reaction force based on current computed _lagrangianForce always called after GetLagrangian Force otherwise gives wrong result
		reactionForce[0] = _lagrangianForce[0];
		reactionForce[1] = _lagrangianForce[1];
		for (auto i = 1; i < Constants::VEC_DIM; ++i) {
			reactionForce[0] += _lagrangianForce[(Constants::DIM - 1)*i];
			reactionForce[1] += _lagrangianForce[(Constants::DIM - 1)*i + 1];
		}

		eulerianProfile->GetProfileForceEulerian(carObj, eulerForce);
		// Add reaction component on the tyre
		for (auto i = 0; i < Constants::NUM_LEGS; ++i) {
			Math::axpy<T>(Constants::DIM, -0.25, reactionForce, Constants::INCX, _lagrangianForce + 4 + 2 * (Constants::DIM - 1)*i, Constants::INCX);
		}
		ComputeInternalTorqueFromLagrangianForce(eulerForce + 1, _lagrangianForce);
		// Add external forces
		eulerForce[0] += externalForce[2];
		for (auto i = 1; i < Constants::VEC_DIM; ++i) {
			eulerForce[2 + i] += externalForce[i*Constants::DIM + 2];
		}
	}

	void ComputeForceALE(T _time, T* lagrangeForce, T* eulerForce) {
		GetLagrangianForce(_time, lagrangeForce);
		GetEulerianForce(_time, eulerForce);
	}

	void GetTorqueLagrange(T _time, T* torque) {
		lagrangianProfile->GetProfileTorqueLagrangian(carObj, torque);
	}

    /**
     * Calculates the internal forces under certain conditions
     * called in:  ALE.solve, global_frame_solver
     * \param time_t current simulation time
     * \param Delta_x_vec current dx of the spring lengths (in Stefan's ordering)
     * \return F_vec vector of forces in the car [GC:XYZ,W1:XYZ,T1:XYZ, ...]
     * \return Normal_ext external forces acting on the whole car system [XYZ]
     */
//    void computeForce3D(T time_t, T* F_vec, T* Normal_ext) {
//        Math::scal<T>(Constants::DIM * Constants::VEC_DIM, 0, F_vec, 1);
//        Math::scal<T>(Constants::DIM, 0, Normal_ext, 1);
//        // fix this so TEOOOOOOOOOO
//        // why did Shubham tell me to comment everything out ?
//        activeProfile->get_Profile_force_ALE(carObj, F_vec, Normal_ext);
//        /* Modify the profile to know where the ground is and apply normal force accordingly */
//        // F_Ti += -0.25 * N; [F[2], F[4], F[6], F[8]]
//        for (auto i = 2; i < Constants::VEC_DIM; i += 2) {
//            Math::axpy<T>(Constants::DIM, -0.25, Normal_ext, 1, &F_vec[i * Constants::DIM], 1);
//        }
//        // PAY ATTENTION to THIS
//        // N += external_force
//        // this formulation is WRONG (!!!) if done before computing following steps, should be done
//        // at the end. external force doesn't necessarily have to create a normal force it can
//        // create acceleration, ex: when car flies
//        Math::vAdd<T>(Constants::DIM, External_force, Normal_ext, Normal_ext);
//        // PAY ATTENTION to THIS
//
//        // get stiffnesses vector k_vec
//#if MIGHT_BE_USEFUL
//        carObj->get_k_vec(k_vec);
//
//        for (auto i = 0; i < (Constants::VEC_DIM - 1); i += 2) {
//            // use the elastic forces at wheels F_CG += k_wi * delta_x_i
//            F_vec[2] += 0. * k_vec[i] * Delta_x_vec[i];
//            Math::axpy<T>(DIM, k_vec[i], &Delta_x_vec[DIM * i], 1, F_vec, 1);
//            F_W_i += -k_wi * delta_x_i Math::axpy<T>(DIM, -k_vec[i], &Delta_x_vec[DIM * i], 1,
//                                                     &F_vec[DIM(i + 1)], 1);
//            F_vec[Constants::DIM * (i + 1) + 2] -= 0. * k_vec[i] * Delta_x_vec[i];
//
//            // use the elastic forces at tyres F_W_i += k_t_i * delta_x_{i+1}
//            Math::axpy<T>(DIM, k_vec[i + 1], &Delta_x_vec[DIM * (i + 1)], 1, &F_vec[DIM * (i + 1)],
//                          1);
//            F_vec[Constants::DIM * (i + 1) + 2] += 0. * k_vec[i + 1] * Delta_x_vec[(i + 1)];
//            // F_T_i += -k_t_i * delta_x_{i+1}
//            Math::axpy<T>(DIM, -k_vec[i + 1], &Delta_x_vec[DIM * (i + 1)], 1, &F_vec[DIM * (i + 2)],
//                          1);
//            F_vec[Constants::DIM * (i + 2) + 2] -= 0. * k_vec[i + 1] * Delta_x_vec[(i + 1)];
//        }
//#endif
//
//        // test add
//        Math::axpy<T>(Constants::DIM, 1, External_force, 1, F_vec, 1);
//    }
//
//    /**
//     * Calculates the torque acting on the whole car
//     * called in:  ALE.solve, global_frame_solver
//     * \param time_t current simulation time
//     * \param Delta_x_vec current dx of the spring lengths (in Stefan's ordering)
//     * \return Torque acting on the total car system [XYZ]
//     */
//    void update_torque(T time_t, T* Torque, T* F_vec) {
//        activeProfile->get_Profile_torque(carObj, Torque);
//        ComputeInternalTorque(Torque, F_vec);
//    }
//
//    void set_External_force() {
//        // in LoadModule constructor commented
//        // TODO: remove method or code in constructor (call in constructor commented out).
//        // TODO: using xml_start/positon_start not be needed with all the getters in MetaDataBase?
//        // TODO: magic '3' and '6' below look like Constants::DIM.
//        auto& db = MetaDataBase<T>::getDataBase();
//
//        // copy the center of mass position
//        Math::copy<T>(Constants::DIM, db.getBodyExternalForce(), 1, External_force, 1);
//        T *xml_start, *position_start;
//        xml_start = db.getWheelExternalForceFrontLeft();
//        position_start = External_force + 3;
//        Math::copy<T>(Constants::DIM, xml_start, 1, position_start, 1);
//        // W2 = W_fr
//        xml_start = db.getWheelExternalForceFrontRight();
//        position_start += 6;  // skip 3 for tyre
//        Math::copy<T>(Constants::DIM, xml_start, 1, position_start, 1);
//        // W3 = W_rl
//        xml_start = db.getWheelExternalForceRearLeft();
//        position_start += 6;  // skip 3 for tyre
//        Math::copy<T>(Constants::DIM, xml_start, 1, position_start, 1);
//        // W2 = W_rr
//        xml_start = db.getWheelExternalForceRearRight();
//        position_start += 6;  // skip 3 for tyre
//        Math::copy<T>(Constants::DIM, xml_start, 1, position_start, 1);
//
//        // T1 = T_fl
//        xml_start = db.getTyreExternalForceFrontLeft();
//        position_start = External_force + 6;  // skip 3 for center of mass and 3 for the wheel
//        Math::copy<T>(Constants::DIM, xml_start, 1, position_start, 1);
//        // T2 = T_fr
//        xml_start = db.getTyreExternalForceFrontRight();
//        position_start += 6;  // skip 3 for the wheel
//        Math::copy<T>(Constants::DIM, xml_start, 1, position_start, 1);
//        // T3 = T_rl
//        xml_start = db.getTyreExternalForceRearLeft();
//        position_start += 6;  // skip 3 for the wheel
//        Math::copy<T>(Constants::DIM, xml_start, 1, position_start, 1);
//        // T4 = T_rr
//        xml_start = db.getTyreExternalForceRearRight();
//        position_start += 6;  // skip 3 for the wheel
//        Math::copy<T>(Constants::DIM, xml_start, 1, position_start, 1);
//    }

};

}  // namespace EVAA
