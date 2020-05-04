
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
	void (LoadModule<T>::*ComputeInternalTorqueFromLagrangianForce)(T*, T*);
	

    /**
     * \brief adds the effect of the internal forces on the torque
     * \param Torque current torque of the body [XYZ]
     * \param F_vec vector of all internal XY-forces in the Eulerian frame at each position
     * [GC:XYZ,W1:XYZ,T1:XYZ, ...]
     */
    void ComputeInternalTorqueFromLagrangianForceFixed(T* Torque, T* F_vec) {
        // compute torque around X-axis
		for (int i = 0; i < 2 * Constants::NUM_LEGS + 1; ++i) {
            Torque[0] += F_vec[i * (Constants::DIM - 1) + 1] * carObj->currentCIRTwoTrackModel[i];
        }

        // compute torque around Y-axis
		for (int i = 0; i < 2 * Constants::NUM_LEGS + 1; ++i) {
            Torque[1] += F_vec[i * (Constants::DIM - 1) + 0] * carObj->currentCIRTwoTrackModel[i];
        }
    }

	/**
	 * \brief adds the effect of the internal forces on the torque for non fixed
	 * \param Torque current torque of the body [XYZ]
	 * \param F_vec vector of all internal XY-forces in the Eulerian frame at each position
	 * [GC:XYZ,W1:XYZ,T1:XYZ, ...]
	 */
	void ComputeInternalTorqueFromLagrangianForceNonFixed(T* Torque, T* F_vec) {}

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
		if (eulerProfile->GetName() == "Fixed") {
			ComputeInternalTorqueFromLagrangianForce = &LoadModule<T>::ComputeInternalTorqueFromLagrangianForceFixed;
		}
		else if (eulerProfile->GetName() == "Nonfixed") {
			ComputeInternalTorqueFromLagrangianForce = &LoadModule<T>::ComputeInternalTorqueFromLagrangianForceNonFixed;
		}

	}
    ~LoadModule() {
        Math::free<T>(externalForce);
		Math::free<T>(reactionForce);
		if (aleCounter) {
			Math::free<T>(_lagrangianForce);
			Math::free<T>(_eulerianForce);
		}
    }

	std::string GetEulerProfileName() {
		return eulerianProfile->GetName();
	}

	std::string GetLagrangianProfileName() {
		return lagrangianProfile->GetName();
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
			Math::axpy<T>(Constants::DIM, -0.25, reactionForce, Constants::INCX, _lagrangianForce + Constants::TYRE_INDEX_LAGRANGE[i], Constants::INCX);
		}
		(this->*ComputeInternalTorqueFromLagrangianForce)(eulerForce + 1, _lagrangianForce);
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
};

}  // namespace EVAA
