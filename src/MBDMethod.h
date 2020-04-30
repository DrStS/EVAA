// TODO: Copyright header

#pragma once

#include <chrono>
#include <string>

#include "Constants.h"
#include "MathLibrary.h"
#include "MetaDataBase.h"

namespace EVAA {

/**
 * Handles the whole MBD simulation (from the Matlab code)
 */
template <class T>
class MBDMethod {
private:
    // Simulation Parameters

    T h;
    size_t num_iter;
    int max_iter;
    T tol;
    size_t solution_dim;  /// this is by the formulation
    std::string solver_name;
    MBDSolver used_solver;

    // Environment conditions

    BoundaryConditionRoad boundary_conditions;
    T radius_circular_path;
    T* center_of_circle;

    // Car Definition

    // bool use_interpolation;
    T k_body_fl;
    T k_tyre_fl;
    T k_body_fr;
    T k_tyre_fr;
    T k_body_rl;
    T k_tyre_rl;
    T k_body_rr;
    T k_tyre_rr;
    T k_body_rot_fl;
    T k_body_rot_fr;
    T k_body_rot_rl;
    T k_body_rot_rr;
    T k_tyre_rot_fl;
    T k_tyre_rot_fr;
    T k_tyre_rot_rl;
    T k_tyre_rot_rr;
    T c_body_fl;
    T c_tyre_fl;
    T c_body_fr;
    T c_tyre_fr;
    T c_body_rl;
    T c_tyre_rl;
    T c_body_rr;
    T c_tyre_rr;
    T l_long_fl;
    T l_long_fr;
    T l_long_rl;
    T l_long_rr;
    T l_lat_fl;
    T l_lat_fr;
    T l_lat_rl;
    T l_lat_rr;
    T mass;
    T I_body_xx;
    T I_body_yy;
    T I_body_zz;
    T mass_wheel_fl;
    T mass_tyre_fl;
    T mass_wheel_fr;
    T mass_tyre_fr;
    T mass_wheel_rl;
    T mass_tyre_rl;
    T mass_wheel_rr;
    T mass_tyre_rr;
    T upper_spring_length_rr;
    T upper_spring_length_rl;
    T upper_spring_length_fl;
    T upper_spring_length_fr;
    T lower_spring_length_rr;
    T lower_spring_length_rl;
    T lower_spring_length_fl;
    T lower_spring_length_fr;
    T* Ic;
    T *mass_wheel, *mass_tyre;
    T *upper_spring_length, *lower_spring_length;
    T *upper_spring_stiffness, *lower_spring_stiffness, *upper_spring_damping,
        *lower_spring_damping;
    T *upper_rotational_stiffness, *lower_rotational_stiffness;
    T *vc, *vw_fl, *vw_fr, *vw_rl, *vw_rr, *vt_fl, *vt_fr, *vt_rl,
        *vt_rr;  // velocity of the center of mass, wheel and tyre.

    // External Forces terms

    T g;
    T* FC;
    T *FT_fl, *FT_fr, *FT_rl, *FT_rr;
    T *FW_fl, *FW_fr, *FW_rl, *FW_rr;

    // Initial condition params

    T *initial_upper_spring_length, *initial_lower_spring_length, *initial_orientation,
        *initial_angular_velocity;

    // Arrays for intermediate steps

    T *r_fl, *r_fr, *r_rl, *r_rr;
    T* pcc;  // position of center of mass
    T* x_vector;

    // Auxillary Parameters for Compute_f function

    T *r_fl_tilda, *r_fr_tilda, *r_rl_tilda, *r_rr_tilda, *A_Ic, *A_rem;

    // Variables needed in compute_f function

    T *cf_C_cN, *cf_r_up_fl, *cf_r_up_fr, *cf_r_up_rl, *cf_r_up_rr, *cf_r_low_fl, *cf_r_low_fr,
        *cf_r_low_rl, *cf_r_low_rr;
    T *cf_upper_normal_fl, *cf_upper_normal_fr, *cf_upper_normal_rl, *cf_upper_normal_rr,
        *cf_lower_normal_fl, *cf_lower_normal_fr, *cf_lower_normal_rl, *cf_lower_normal_rr,
        *cf_col_dat;
    T *cf_upper_force_fl, *cf_upper_force_fr, *cf_upper_force_rl, *cf_upper_force_rr,
        *cf_lower_force_fl, *cf_lower_force_fr, *cf_lower_force_rl, *cf_lower_force_rr;
    T *cf_upper_dampf_fl, *cf_upper_dampf_fr, *cf_upper_dampf_rl, *cf_upper_dampf_rr,
        *cf_lower_dampf_fl, *cf_lower_dampf_fr, *cf_lower_dampf_rl, *cf_lower_dampf_rr, *cf_temp;
    T *cf_upper_angle_fl, *cf_upper_angle_fr, *cf_upper_angle_rl, *cf_upper_angle_rr,
        *cf_lower_angle_fl, *cf_lower_angle_fr, *cf_lower_angle_rl, *cf_lower_angle_rr;
    T *cf_upper_S_fl, *cf_upper_S_fr, *cf_upper_S_rl, *cf_upper_S_rr, *cf_lower_S_fl,
        *cf_lower_S_fr, *cf_lower_S_rl, *cf_lower_S_rr;
    T *cf_lower_rot_force_fl, *cf_lower_rot_force_fr, *cf_lower_rot_force_rl,
        *cf_lower_rot_force_rr, *cf_upper_rot_force_fl, *cf_upper_rot_force_fr,
        *cf_upper_rot_force_rl, *cf_upper_rot_force_rr;
    T *cf_car_rot_force_fl, *cf_car_rot_force_fr, *cf_car_rot_force_rl, *cf_car_rot_force_rr,
        *cf_sum_car_force_fl, *cf_sum_car_force_fr, *cf_sum_car_force_rl, *cf_sum_car_force_rr;
    T *cf_local_FR_fl, *cf_local_FR_fr, *cf_local_FR_rl, *cf_local_FR_rr;
    T *cf_Hc, *cf_sum_torque_spring_car, *cf_Tc, *cf_wc_tilda;
    T *cf_b_rem, *cf_Qc, *cf_qc_dot, *accelerations;
    T *current_spring_lengths, *stiffness_vector;
    T *wc_, *vc_, *vw_fl_, *vw_fr_, *vw_rl_, *vw_rr_, *vt_fl_, *vt_fr_, *vt_rl_, *vt_rr_, *qc_;
    T *pcc_, *pw_fl_, *pw_fr_, *pw_rl_, *pw_rr_, *pt_fl_, *pt_fr_, *pt_rl_, *pt_rr_;
    T inv_norm_r_up_fl, inv_norm_r_up_fr, inv_norm_r_up_rl, inv_norm_r_up_rr;
    T inv_norm_r_low_fl, inv_norm_r_low_fr, inv_norm_r_low_rl, inv_norm_r_low_rr;
    T norm_r_up_fl, norm_r_up_fr, norm_r_up_rl, norm_r_up_rr;
    T norm_r_low_fl, norm_r_low_fr, norm_r_low_rl, norm_r_low_rr;

    /**
     * Calculate the positions of the tyres and wheels according to the initial orientation of the
     * car The legs always form a 90� angle to the car body, such that the rotational springs are
     * at rest
     */
    void get_initial_length(T* initial_orientation_, const T* r_fl_, const T* r_fr_, const T* r_rl_,
                            const T* r_rr_, const T* pcc_, const T* initial_upper_spring_length_,
                            const T* initial_lower_spring_length_, T* wheel_coordinate_fl_,
                            T* wheel_coordinate_fr_, T* wheel_coordinate_rl_,
                            T* wheel_coordinate_rr_, T* tyre_coordinate_fl_, T* tyre_coordinate_fr_,
                            T* tyre_coordinate_rl_, T* tyre_coordinate_rr_)

    {
        /*
         * To reduce memory trace and better use cache this function is implemented in following
         * fashion:
         * Original steps for computation of one component:
         * 1. qc = qc/norm(qc);
         * 2. C_Nc = get_basis(qc);
         * 3. global_z = C_Nc(:,2);
         * 4. global_z = -global_z / norm(global_z);
         * 5. global_r_fl = pcc + C_Nc*r_fl;
         * 6. upper_global_spring__fl = upper_length(_fl) * global_z;
         * 7. lower_global_spring__fl = lower_length(_fl)*global_z;
         * 8. pw_fl = global_r_fl + upper_global_spring__fl;
         * 9. pt_fl = pw_fl + lower_global_spring__fl;
         * Modified steps for computation of one component:
         * 1. qc = qc/norm(qc);
         * 2. C_Nc = get_basis(qc);
         * 3. global_z = C_Nc(:,2);
         * 4. global_z = global_z / norm(global_z);
         * 5. pw_fl = pcc;
         * 6. pw_fl = pw_fl + C_Nc*r_fl;
         * 7. pw_fl = pw_fl + upper_length(1)*global_z;
         * 8. pt_fl = pw_fl
         * 9. pt_fl	= pt_fl + lower_length(_fl)*global_z;
         */

        T* global_z = Math::calloc<T>(Constants::DIM);
        T* C_Nc = Math::calloc<T>(Constants::DIM * Constants::DIM);

        // 1. qc = qc/norm(qc); This is in quaternions
        T nrm = Math::nrm2<T>(Constants::NUM_LEGS, initial_orientation_, 1);
        Math::scal<T>(Constants::NUM_LEGS, 1. / nrm, initial_orientation_, 1);

        // 2. C_Nc = get_basis(qc);
        Math::get_basis<T>(initial_orientation_, C_Nc);

        // 3. global_z = C_Nc(:,2);
        Math::copy<T>(Constants::DIM, C_Nc + 2, Constants::DIM, global_z, 1);

        // 4. global_z = -global_z / norm(global_z);
        nrm = Math::nrm2<T>(Constants::DIM, global_z, 1);
        Math::scal<T>(Constants::DIM, -1. / nrm, global_z, 1);

        // Leg 1

        // 5.	pw1 = pcc;
        Math::copy<T>(Constants::DIM, pcc_, 1, wheel_coordinate_fl_, 1);

        // 6.	pw1 = pw1 + C_Nc*r1;
        Math::gemv<T>(CblasRowMajor, CblasNoTrans, Constants::DIM, Constants::DIM, 1, C_Nc,
                      Constants::DIM, r_fl_, 1, 1, wheel_coordinate_fl_, 1);

        // 7.	pw1 = pw1 + upper_length(1)*global_z;
        Math::axpy<T>(Constants::DIM, initial_upper_spring_length_[0], global_z, 1,
                      wheel_coordinate_fl_, 1);

        // 8.	pt1 = pw1
        Math::copy<T>(Constants::DIM, wheel_coordinate_fl_, 1, tyre_coordinate_fl_, 1);

        // 9.	pt1 = pw1 + lower_length(1)*global_z;
        Math::axpy<T>(Constants::DIM, initial_lower_spring_length_[0], global_z, 1,
                      tyre_coordinate_fl_, 1);

        // Leg 2

        // 5.	pw2 = pcc;
        Math::copy<T>(Constants::DIM, pcc_, 1, wheel_coordinate_fr_, 1);

        // 6.	pw2 = pw2 + C_Nc*r2;
        Math::gemv<T>(CblasRowMajor, CblasNoTrans, Constants::DIM, Constants::DIM, 1, C_Nc,
                      Constants::DIM, r_fr_, 1, 1, wheel_coordinate_fr_, 1);

        // 7.	pw2 = pw2 + upper_length(2)*global_z;
        Math::axpy<T>(Constants::DIM, initial_upper_spring_length_[1], global_z, 1,
                      wheel_coordinate_fr_, 1);

        // 8.	pt2 = pw2
        Math::copy<T>(Constants::DIM, wheel_coordinate_fr_, 1, tyre_coordinate_fr_, 1);

        // 9.	pt2 = pw2 + lower_length(2)*global_z;
        Math::axpy<T>(Constants::DIM, initial_lower_spring_length_[1], global_z, 1,
                      tyre_coordinate_fr_, 1);

        // Leg 3

        // 5.	pw3 = pcc;
        Math::copy<T>(Constants::DIM, pcc_, 1, wheel_coordinate_rl_, 1);

        // 6.	pw3 = pw3 + C_Nc*r3;
        Math::gemv<T>(CblasRowMajor, CblasNoTrans, Constants::DIM, Constants::DIM, 1, C_Nc,
                      Constants::DIM, r_rl_, 1, 1, wheel_coordinate_rl_, 1);

        // 7.	pw3 = pw3 + upper_length(3)*global_z;
        Math::axpy<T>(Constants::DIM, initial_upper_spring_length_[2], global_z, 1,
                      wheel_coordinate_rl_, 1);

        // 8.	pt3 = pw3
        Math::copy<T>(Constants::DIM, wheel_coordinate_rl_, 1, tyre_coordinate_rl_, 1);

        // 9.	pt3 = pw3 + lower_length(3)*global_z;
        Math::axpy<T>(Constants::DIM, initial_lower_spring_length_[2], global_z, 1,
                      tyre_coordinate_rl_, 1);

        // Leg 4

        // 5.	pw4 = pcc;
        Math::copy<T>(Constants::DIM, pcc_, 1, wheel_coordinate_rr_, 1);

        // 6.	pw4 = pw4 + C_Nc*r4;
        Math::gemv<T>(CblasRowMajor, CblasNoTrans, Constants::DIM, Constants::DIM, 1, C_Nc,
                      Constants::DIM, r_rr_, 1, 1, wheel_coordinate_rr_, 1);

        // 7.	pw4 = pw4 + upper_length(4)*global_z;
        Math::axpy<T>(Constants::DIM, initial_upper_spring_length_[3], global_z, 1,
                      wheel_coordinate_rr_, 1);

        // 8.	pt4 = pw4
        Math::copy<T>(Constants::DIM, wheel_coordinate_rr_, 1, tyre_coordinate_rr_, 1);

        // 9.	pt4 = pw4 + lower_length(4)*global_z;
        Math::axpy<T>(Constants::DIM, initial_lower_spring_length_[3], global_z, 1,
                      tyre_coordinate_rr_, 1);

        Math::free<T>(global_z);
        Math::free<T>(C_Nc);
    }

    /**
     * For Debug purposes
     */
    void write_matrix(T* vect, int count) {
        std::cout << "Debug mode print" << std::endl;
        for (size_t i = 0; i < count; ++i) {
            // std::cout << vect[i] << std::endl;
            std::cout.precision(5);
            for (size_t j = 0; j < count; ++j) {
                std::cout << std::scientific << vect[i * count + j] << "  ";
            }
            std::cout << "\n" << std::endl;
        }
    }

    /**
     * For Debug purposes
     */
    void write_vector(T* vect, int count) {
        std::cout << "Debug mode print" << std::endl;
        for (size_t i = 0; i < count; ++i) {
            std::cout.precision(15);
            std::cout << std::scientific << vect[i] << std::endl;
        }
    }

    void MemoryAllocation() {
        // Memory Allocation and matrix formulation

        r_fl = Math::calloc<T>(Constants::DIM);
        r_fr = Math::calloc<T>(Constants::DIM);
        r_rl = Math::calloc<T>(Constants::DIM);
        r_rr = Math::calloc<T>(Constants::DIM);
        Ic = Math::calloc<T>(Constants::DIM * Constants::DIM);
        mass_wheel = Math::calloc<T>(Constants::NUM_LEGS);
        mass_tyre = Math::calloc<T>(Constants::NUM_LEGS);
        upper_spring_length = Math::calloc<T>(Constants::NUM_LEGS);
        lower_spring_length = Math::calloc<T>(Constants::NUM_LEGS);
        upper_spring_stiffness = Math::calloc<T>(Constants::NUM_LEGS);
        lower_spring_stiffness = Math::calloc<T>(Constants::NUM_LEGS);
        upper_rotational_stiffness = Math::calloc<T>(Constants::NUM_LEGS);
        lower_rotational_stiffness = Math::calloc<T>(Constants::NUM_LEGS);
        initial_upper_spring_length = Math::calloc<T>(Constants::NUM_LEGS);
        initial_lower_spring_length = Math::calloc<T>(Constants::NUM_LEGS);
        initial_orientation = Math::calloc<T>(Constants::NUM_LEGS);
        initial_angular_velocity = Math::calloc<T>(Constants::DIM);
        upper_spring_damping = Math::calloc<T>(Constants::NUM_LEGS);
        lower_spring_damping = Math::calloc<T>(Constants::NUM_LEGS);
        vc = Math::calloc<T>(Constants::DIM);
        vw_fl = Math::calloc<T>(Constants::DIM);
        vw_fr = Math::calloc<T>(Constants::DIM);
        vw_rl = Math::calloc<T>(Constants::DIM);
        vw_rr = Math::calloc<T>(Constants::DIM);
        vt_fl = Math::calloc<T>(Constants::DIM);
        vt_fr = Math::calloc<T>(Constants::DIM);
        vt_rl = Math::calloc<T>(Constants::DIM);
        vt_rr = Math::calloc<T>(Constants::DIM);
        pcc = Math::calloc<T>(Constants::DIM);
        FC = Math::calloc<T>(Constants::DIM);
        r_fl_tilda = Math::calloc<T>(Constants::DIM * Constants::DIM);
        r_fr_tilda = Math::calloc<T>(Constants::DIM * Constants::DIM);
        r_rl_tilda = Math::calloc<T>(Constants::DIM * Constants::DIM);
        r_rr_tilda = Math::calloc<T>(Constants::DIM * Constants::DIM);
        FT_fl = Math::calloc<T>(Constants::DIM);
        FT_rr = Math::calloc<T>(Constants::DIM);
        FT_rl = Math::calloc<T>(Constants::DIM);
        FT_fr = Math::calloc<T>(Constants::DIM);
        FW_fl = Math::calloc<T>(Constants::DIM);
        FW_rr = Math::calloc<T>(Constants::DIM);
        FW_rl = Math::calloc<T>(Constants::DIM);
        FW_fr = Math::calloc<T>(Constants::DIM);
        A_Ic = Math::calloc<T>(Constants::DIM * Constants::DIM);
        A_rem = Math::calloc<T>(Constants::DIM * Constants::DIM * Constants::DIM);
        center_of_circle = Math::calloc<T>(Constants::DIM);
        accelerations = Math::calloc<T>(Constants::DIM * Constants::DIM * Constants::DIM);
    }

    void ReadFromXML() {
        // Simulation Parameters
        auto& db = MetaDataBase<T>::getDataBase();

        h = db.getTimeStepSize();
        num_iter = db.getNumberOfTimeIterations();
        max_iter = db.getMaxNumberOfBroydenIterationForMBD();
        tol = db.getToleranceBroydenIterationForMBD();
        solution_dim = db.getSolutionVectorSize();  /// this is by the formulation
        used_solver = db.getUsedSolverForMBD();
        boundary_conditions = db.getRoadConditions();
        radius_circular_path = db.getCircularRoadRadius();
        // use_interpolation = db.getUseInterpolation();

        // Car Definition

        k_body_fl = db.getBodyStiffnessFrontLeft();
        k_tyre_fl = db.getTyreStiffnessFrontLeft();
        k_body_fr = db.getBodyStiffnessFrontRight();
        k_tyre_fr = db.getTyreStiffnessFrontRight();
        k_body_rl = db.getBodyStiffnessRearLeft();
        k_tyre_rl = db.getTyreStiffnessRearLeft();
        k_body_rr = db.getBodyStiffnessRearRight();
        k_tyre_rr = db.getTyreStiffnessRearRight();

        k_body_rot_fl = 1e5;
        k_body_rot_fr = 1e5;
        k_body_rot_rl = 1e5;
        k_body_rot_rr = 1e5;
        k_tyre_rot_fl = 1e5;
        k_tyre_rot_fr = 1e5;
        k_tyre_rot_rl = 1e5;
        k_tyre_rot_rr = 1e5;

        c_body_fl = db.getBodyDampingFrontLeft();
        c_tyre_fl = db.getTyreDampingFrontLeft();
        c_body_fr = db.getBodyDampingFrontRight();
        c_tyre_fr = db.getTyreDampingFrontRight();
        c_body_rl = db.getBodyDampingRearLeft();
        c_tyre_rl = db.getTyreDampingRearLeft();
        c_body_rr = db.getBodyDampingRearRight();
        c_tyre_rr = db.getTyreDampingRearRight();

        l_long_fl = db.getLongitudalLegPositionFrontLeft();
        l_long_fr = db.getLongitudalLegPositionFrontRight();
        l_long_rl = db.getLongitudalLegPositionRearLeft();
        l_long_rr = db.getLongitudalLegPositionRearRight();

        l_lat_fl = db.getLatidudalLegPositionFrontLeft();
        l_lat_fr = db.getLatidudalLegPositionFrontRight();
        l_lat_rl = db.getLatidudalLegPositionRearLeft();
        l_lat_rr = db.getLatidudalLegPositionRearRight();

        mass = db.getBodyMass();

        I_body_xx = db.getMomentOfInertiaXX();
        I_body_yy = db.getMomentOfInertiaYY();
        I_body_zz = db.getMomentOfInertiaZZ();

        mass_wheel_fl = db.getWheelMassFrontLeft();
        mass_tyre_fl = db.getTyreMassFrontLeft();
        mass_wheel_fr = db.getWheelMassFrontRight();
        mass_tyre_fr = db.getTyreMassFrontRight();
        mass_wheel_rl = db.getWheelMassRearLeft();
        mass_tyre_rl = db.getTyreMassRearLeft();
        mass_wheel_rr = db.getWheelMassRearRight();
        mass_tyre_rr = db.getTyreMassRearRight();

        upper_spring_length_fl = db.getBodySpringLengthFrontLeft();
        upper_spring_length_fr = db.getBodySpringLengthFrontRight();
        upper_spring_length_rl = db.getBodySpringLengthRearLeft();
        upper_spring_length_rr = db.getBodySpringLengthRearRight();

        lower_spring_length_fl = db.getTyreSpringLengthFrontLeft();
        lower_spring_length_fr = db.getTyreSpringLengthFrontRight();
        lower_spring_length_rl = db.getTyreSpringLengthRearLeft();
        lower_spring_length_rr = db.getTyreSpringLengthRearRight();

        g = db.getGravityField()[2];

        // Fill up vectors
        int i;

        i = 0;
        r_fl[i] = l_long_fl;
        r_fr[i] = l_long_fr;
        r_rl[i] = -l_long_rl;
        r_rr[i] = -l_long_rr;
        Ic[i * Constants::DIM + i] = I_body_xx;
        mass_wheel[i] = mass_wheel_fl;
        mass_tyre[i] = mass_tyre_fl;
        upper_spring_length[i] = upper_spring_length_fl;
        lower_spring_length[i] = lower_spring_length_fl;
        upper_spring_stiffness[i] = k_body_fl;
        lower_spring_stiffness[i] = k_tyre_fl;
        upper_rotational_stiffness[i] = k_body_rot_fl;
        lower_rotational_stiffness[i] = k_tyre_rot_fl;
        upper_spring_damping[i] = c_body_fl;
        lower_spring_damping[i] = c_tyre_fl;

        initial_upper_spring_length[i] = db.getBodySpringInitialLengthFrontLeft();
        initial_lower_spring_length[i] = db.getTyreSpringInitialLengthFrontLeft();
        initial_orientation[i] = db.getBodyInitialOrientation()[i];
        vc[i] = db.getBodyInitialVelocity()[i];
        vw_fl[i] = db.getWheelInitialVelocityFrontLeft()[i];
        vw_fr[i] = db.getWheelInitialVelocityFrontRight()[i];
        vw_rl[i] = db.getWheelInitialVelocityRearLeft()[i];
        vw_rr[i] = db.getWheelInitialVelocityRearRight()[i];
        vt_fl[i] = db.getTyreInitialVelocityFrontLeft()[i];
        vt_fr[i] = db.getTyreInitialVelocityFrontRight()[i];
        vt_rl[i] = db.getTyreInitialVelocityRearLeft()[i];
        vt_rr[i] = db.getTyreInitialVelocityRearRight()[i];
        pcc[i] = db.getBodyInitialPosition()[i];
        FC[i] = db.getBodyExternalForce()[i];
        FT_fl[i] = db.getTyreExternalForceFrontLeft()[i];
        FT_fr[i] = db.getTyreExternalForceFrontRight()[i];
        FT_rl[i] = db.getTyreExternalForceRearLeft()[i];
        FT_rr[i] = db.getTyreExternalForceRearRight()[i];
        FW_fl[i] = db.getWheelExternalForceFrontLeft()[i];
        FW_fr[i] = db.getWheelExternalForceFrontRight()[i];
        FW_rl[i] = db.getWheelExternalForceRearLeft()[i];
        FW_rr[i] = db.getWheelExternalForceRearRight()[i];
        center_of_circle[i] = db.getCircularRoadCenter()[i];

        i = 1;
        r_fl[i] = -l_lat_fl;
        r_fr[i] = l_lat_fr;
        r_rl[i] = -l_lat_rl;
        r_rr[i] = l_lat_rr;
        Ic[i * Constants::DIM + i] = I_body_yy;
        mass_wheel[i] = mass_wheel_fr;
        mass_tyre[i] = mass_tyre_fr;
        upper_spring_length[i] = upper_spring_length_fr;
        lower_spring_length[i] = lower_spring_length_fr;
        upper_spring_stiffness[i] = k_body_fr;
        lower_spring_stiffness[i] = k_tyre_fr;
        upper_rotational_stiffness[i] = k_body_rot_fr;
        lower_rotational_stiffness[i] = k_tyre_rot_fr;
        upper_spring_damping[i] = c_body_fr;
        lower_spring_damping[i] = c_tyre_fr;
        initial_upper_spring_length[i] = db.getBodySpringInitialLengthFrontRight();
        initial_lower_spring_length[i] = db.getTyreSpringInitialLengthFrontRight();
        initial_orientation[i] = db.getBodyInitialOrientation()[i];
        vc[i] = db.getBodyInitialVelocity()[i];
        vw_fl[i] = db.getWheelInitialVelocityFrontLeft()[i];
        vw_fr[i] = db.getWheelInitialVelocityFrontRight()[i];
        vw_rl[i] = db.getWheelInitialVelocityRearLeft()[i];
        vw_rr[i] = db.getWheelInitialVelocityRearRight()[i];
        vt_fl[i] = db.getTyreInitialVelocityFrontLeft()[i];
        vt_fr[i] = db.getTyreInitialVelocityFrontRight()[i];
        vt_rl[i] = db.getTyreInitialVelocityRearLeft()[i];
        vt_rr[i] = db.getTyreInitialVelocityRearRight()[i];
        pcc[i] = db.getBodyInitialPosition()[i];
        FC[i] = db.getBodyExternalForce()[i];
        FT_fl[i] = db.getTyreExternalForceFrontLeft()[i];
        FT_fr[i] = db.getTyreExternalForceFrontRight()[i];
        FT_rl[i] = db.getTyreExternalForceRearLeft()[i];
        FT_rr[i] = db.getTyreExternalForceRearRight()[i];
        FW_fl[i] = db.getWheelExternalForceFrontLeft()[i];
        FW_fr[i] = db.getWheelExternalForceFrontRight()[i];
        FW_rl[i] = db.getWheelExternalForceRearLeft()[i];
        FW_rr[i] = db.getWheelExternalForceRearRight()[i];
        center_of_circle[i] = db.getCircularRoadCenter()[i];

        i = 2;
        r_fl[i] = 0;
        r_fr[i] = 0;
        r_rl[i] = 0;
        r_rr[i] = 0;
        Ic[i * Constants::DIM + i] = I_body_zz;
        mass_wheel[i] = mass_wheel_rl;
        mass_tyre[i] = mass_tyre_rl;
        upper_spring_length[i] = upper_spring_length_rl;
        lower_spring_length[i] = lower_spring_length_rl;
        upper_spring_stiffness[i] = k_body_rl;
        lower_spring_stiffness[i] = k_tyre_rl;
        upper_rotational_stiffness[i] = k_body_rot_rl;
        lower_rotational_stiffness[i] = k_tyre_rot_rl;
        upper_spring_damping[i] = c_body_rl;
        lower_spring_damping[i] = c_tyre_rl;
        initial_upper_spring_length[i] = db.getBodySpringInitialLengthRearLeft();
        initial_lower_spring_length[i] = db.getTyreSpringInitialLengthRearLeft();
        initial_orientation[i] = db.getBodyInitialOrientation()[i];
        vc[i] = db.getBodyInitialVelocity()[i];
        vw_fl[i] = db.getWheelInitialVelocityFrontLeft()[i];
        vw_fr[i] = db.getWheelInitialVelocityFrontRight()[i];
        vw_rl[i] = db.getWheelInitialVelocityRearLeft()[i];
        vw_rr[i] = db.getWheelInitialVelocityRearRight()[i];
        vt_fl[i] = db.getTyreInitialVelocityFrontLeft()[i];
        vt_fr[i] = db.getTyreInitialVelocityFrontRight()[i];
        vt_rl[i] = db.getTyreInitialVelocityRearLeft()[i];
        vt_rr[i] = db.getTyreInitialVelocityRearRight()[i];
        pcc[i] = db.getBodyInitialPosition()[i];
        FC[i] = db.getBodyExternalForce()[i] - mass * g;
        FT_fl[i] = db.getTyreExternalForceFrontLeft()[i] - mass_tyre_fl * g;
        FT_fr[i] = db.getTyreExternalForceFrontRight()[i] - mass_tyre_fr * g;
        FT_rl[i] = db.getTyreExternalForceRearLeft()[i] - mass_tyre_rl * g;
        FT_rr[i] = db.getTyreExternalForceRearRight()[i] - mass_tyre_rr * g;
        FW_fl[i] = db.getWheelExternalForceFrontLeft()[i] - mass_wheel_fl * g;
        FW_fr[i] = db.getWheelExternalForceFrontRight()[i] - mass_wheel_fr * g;
        FW_rl[i] = db.getWheelExternalForceRearLeft()[i] - mass_wheel_rl * g;
        FW_rr[i] = db.getWheelExternalForceRearRight()[i] - mass_wheel_rr * g;
        center_of_circle[i] = db.getCircularRoadCenter()[i];

        i = 3;
        mass_wheel[i] = mass_wheel_rr;
        mass_tyre[i] = mass_tyre_rr;
        upper_spring_length[i] = upper_spring_length_rr;
        lower_spring_length[i] = lower_spring_length_rr;
        upper_spring_stiffness[i] = k_body_rr;
        lower_spring_stiffness[i] = k_tyre_rr;
        upper_rotational_stiffness[i] = k_body_rot_rr;
        lower_rotational_stiffness[i] = k_tyre_rot_rr;
        upper_spring_damping[i] = c_body_rr;
        lower_spring_damping[i] = c_tyre_rr;
        initial_upper_spring_length[i] = db.getBodySpringInitialLengthRearRight();
        initial_lower_spring_length[i] = db.getTyreSpringInitialLengthRearRight();
        initial_orientation[i] = db.getBodyInitialOrientation()[i];
    }

    void getCholeskyDecomposition() {
        // A_Ic has cholesky factorization of Ic
        Math::copy<T>(Constants::DIM * Constants::DIM, Ic, 1, A_Ic, 1);
        Math::potrf<T>(LAPACK_ROW_MAJOR, 'L', Constants::DIM, A_Ic, Constants::DIM);
    }

    /* Road Profile and Load module */

    /** Updates the velocity of the 4 wheels and tyres as well as the angular velocity, such that
     * the car is already in the trajectory of the circle overwrites the velocity in the tyres and
     * wheels and the angular velocities only keeps the tangential component of the velocity of the
     * car body
     */
    void circular_path_initialization(T* vc, T* vw_fl, T* vw_fr, T* vw_rl, T* vw_rr, T* vt_fl,
                                      T* vt_fr, T* vt_rl, T* vt_rr, T* omega, T* pcc, T* pt_fl,
                                      T* pt_fr, T* pt_rl, T* pt_rr, T& radius_param) {
        const MKL_INT dim = Constants::DIM;
        const MKL_INT incx = 1;

        // only consider circular motion in the XY-plane
        vc[2] = 0;

        // memory allocation
        T* perpendicular_dir = Math::calloc<T>(Constants::DIM);
        T* tangential_dir = Math::calloc<T>(Constants::DIM);
        T* radial_vector = Math::calloc<T>(Constants::DIM);

        // vector from the car to the center of the circle
        radial_vector[0] = pcc[0] - center_of_circle[0];
        radial_vector[1] = pcc[1] - center_of_circle[1];
        radial_vector[2] = 0;

        // absolute distance from the car to the center
        T radius = Math::nrm2<T>(Constants::DIM, radial_vector, 1);

        // vector out of the XY-plane (unit Z direction)
        perpendicular_dir[2] = 1;

        if (abs(radius - radius_param) > 0.1)
            std::cout << "Warning! the initial position of the car is not on the trajectory "
                         "provided in the circular path. \n The expected radius is "
                      << radius_circular_path << ", but the car is at an initial distance of "
                      << radius
                      << " from the center of the circle.\n The execution procedes with the "
                         "current spatial configuration and with the current distance to the "
                         "center of the circle."
                      << std::endl;

        T inv_radius = 1. / radius;

        // normalize the radial vector
        Math::scal<T>(dim, inv_radius, radial_vector, incx);

        // calculate the direction of the motion
        Math::crossProduct<T>(radial_vector, perpendicular_dir, tangential_dir);

        // calculate the velocity magnitude
        T magnitude = Math::dot<T>(dim, vc, incx, tangential_dir, incx);

        // get the velocity vector in tangential direction with computed magnitude
        Math::copy<T>(Constants::DIM, tangential_dir, 1, vc, 1);
        Math::scal<T>(dim, magnitude, vc, incx);

        // get the physical radial vector again
        Math::scal<T>(dim, radius, radial_vector, incx);

        // calculate the angular velocity of the car
        Math::crossProduct<T>(radial_vector, vc, omega);
        Math::scal<T>(Constants::DIM, inv_radius * inv_radius, omega, 1);

        // calculate the velocity in all legs
        Math::crossProduct<T>(omega, pt_fl, vt_fl);
        Math::crossProduct<T>(omega, pt_fr, vt_fr);
        Math::crossProduct<T>(omega, pt_rl, vt_rl);
        Math::crossProduct<T>(omega, pt_rr, vt_rr);

        Math::copy<T>(Constants::DIM, vt_fl, 1, vw_fl, 1);
        Math::copy<T>(Constants::DIM, vt_fr, 1, vw_fr, 1);
        Math::copy<T>(Constants::DIM, vt_rl, 1, vw_rl, 1);
        Math::copy<T>(Constants::DIM, vt_rr, 1, vw_rr, 1);

        Math::free<T>(perpendicular_dir);
        Math::free<T>(tangential_dir);
        Math::free<T>(radial_vector);
    }

    /**
     * Fixes the tyres to their initial position
     * The forces acting on the tyres are now always zero
     * \param v is the velocity of the tyre
     * \param m is the mass of the tyre
     * \param p is the global position of the tyre
     * \result Fr The only force acting on the tyre
     */
    void get_fixed_road_force(T* Fr_fl, T* Fr_fr, T* Fr_rl, T* Fr_rr) {
        Fr_fl[0] = 0;
        Fr_fl[1] = 0;
        Fr_fl[2] = 0;

        Fr_fr[0] = 0;
        Fr_fr[1] = 0;
        Fr_fr[2] = 0;

        Fr_rl[0] = 0;
        Fr_rl[1] = 0;
        Fr_rl[2] = 0;

        Fr_rr[0] = 0;
        Fr_rr[1] = 0;
        Fr_rr[2] = 0;
    }

    /**
     * No interaction with the road, no additional forces on the tyres
     * \param v is the velocity of the tyre
     * \param m is the mass of the tyre
     * \param p is the global position of the tyre
     * \result Fr The only force acting on the tyre
     */
    void get_nonfixed_road_force(T* Fr_fl, T* Fr_fr, T* Fr_rl, T* Fr_rr) {}

    /**
     * calculates the force in the tyre only with respect to its velocity, mass and position
     * \param v is the velocity of the tyre
     * \param m is the mass of the tyre
     * \param p is the global position of the tyre
     * \result Fr The only force acting on the tyre
     */
    void get_circular_road_force(T* Fr, T* v, T& m, T* p) {
        T *unit_z_vector, *velocity_direction_tyre;
        T velocity_magnitude_tyre, inv_radius, force_magnitude_tyre;

        // allocate memory
        unit_z_vector = Math::calloc<T>(Constants::DIM);
        velocity_direction_tyre = Math::calloc<T>(Constants::DIM);

        // some MKL constants (DIM=3)
        const MKL_INT mkl_DIM = Constants::DIM;
        const MKL_INT mkl_incx = 1;
        const MKL_INT mkl_incy = 1;

        // the force is in the same direction as the position vector
        // TODO: take care of situation when the center of the circle is not at the origin!! /////
        // RAFFI
        Math::copy<T>(Constants::DIM, p, 1, Fr, 1);

        Fr[2] = 0;  // path only in XZ-plane

        // inverse radius of the trajectory at the considered tyre (see TODO above)
        inv_radius = 1. / Math::nrm2<T>(Constants::DIM, p, 1);

        // normalize the force vector -- minus sign because the force points to the center
        Math::scal<T>(Constants::DIM, -inv_radius, Fr, 1);

        // perpendicular to the motion and radius
        unit_z_vector[0] = 0;
        unit_z_vector[1] = 0;
        unit_z_vector[2] = 1;

        // get normalized direction of motion
        Math::crossProduct<T>(Fr, unit_z_vector, velocity_direction_tyre);

        // get the physical velocity
        velocity_magnitude_tyre =
            Math::dot<T>(mkl_DIM, v, mkl_incx, velocity_direction_tyre, mkl_incy);

        // get the physical force
        force_magnitude_tyre = m * velocity_magnitude_tyre * velocity_magnitude_tyre * inv_radius;
        Math::scal<T>(Constants::DIM, force_magnitude_tyre, Fr, 1);

        // free memory
        Math::free<T>(unit_z_vector);
        Math::free<T>(velocity_direction_tyre);
    }

    /** Functions needed for compute_f */

    /**
     * Memory allocation of all the variables required in the solve function
     * To increase performance by removing repeting memory allocations
     * The same locations are overwritten at each timestep
     */
    void compute_f_mem_alloc() {
        // add to members
        ///// All the members belong to compute f function and are kind of nasty
        //// cf_* means parameter used in compute f function
        cf_C_cN = Math::calloc<T>(Constants::DIM * Constants::DIM);
        cf_r_up_fl = Math::calloc<T>(Constants::DIM);
        cf_r_up_fr = Math::calloc<T>(Constants::DIM);
        cf_r_up_rl = Math::calloc<T>(Constants::DIM);
        cf_r_up_rr = Math::calloc<T>(Constants::DIM);
        cf_r_low_fl = Math::calloc<T>(Constants::DIM);
        cf_r_low_fr = Math::calloc<T>(Constants::DIM);
        cf_r_low_rl = Math::calloc<T>(Constants::DIM);
        cf_r_low_rr = Math::calloc<T>(Constants::DIM);

        cf_upper_normal_fl = Math::calloc<T>(Constants::DIM);
        cf_upper_normal_fr = Math::calloc<T>(Constants::DIM);
        cf_upper_normal_rl = Math::calloc<T>(Constants::DIM);
        cf_upper_normal_rr = Math::calloc<T>(Constants::DIM);
        cf_lower_normal_fl = Math::calloc<T>(Constants::DIM);
        cf_lower_normal_fr = Math::calloc<T>(Constants::DIM);
        cf_lower_normal_rl = Math::calloc<T>(Constants::DIM);
        cf_lower_normal_rr = Math::calloc<T>(Constants::DIM);
        cf_col_dat = Math::calloc<T>(Constants::DIM);

        cf_upper_force_fl = Math::calloc<T>(Constants::DIM);
        cf_upper_force_fr = Math::calloc<T>(Constants::DIM);
        cf_upper_force_rl = Math::calloc<T>(Constants::DIM);
        cf_upper_force_rr = Math::calloc<T>(Constants::DIM);
        cf_lower_force_fl = Math::calloc<T>(Constants::DIM);
        cf_lower_force_fr = Math::calloc<T>(Constants::DIM);
        cf_lower_force_rl = Math::calloc<T>(Constants::DIM);
        cf_lower_force_rr = Math::calloc<T>(Constants::DIM);

        cf_upper_dampf_fl = Math::calloc<T>(Constants::DIM);
        cf_upper_dampf_fr = Math::calloc<T>(Constants::DIM);
        cf_upper_dampf_rl = Math::calloc<T>(Constants::DIM);
        cf_upper_dampf_rr = Math::calloc<T>(Constants::DIM);
        cf_lower_dampf_fl = Math::calloc<T>(Constants::DIM);
        cf_lower_dampf_fr = Math::calloc<T>(Constants::DIM);
        cf_lower_dampf_rl = Math::calloc<T>(Constants::DIM);
        cf_lower_dampf_rr = Math::calloc<T>(Constants::DIM);
        cf_temp = Math::calloc<T>(Constants::DIM);

        // TODO: :))) static vars, maybe an array/vector with all.
        cf_upper_angle_fl = new T;
        cf_upper_angle_fr = new T;
        cf_upper_angle_rl = new T;
        cf_upper_angle_rr = new T;
        cf_lower_angle_fl = new T;
        cf_lower_angle_fr = new T;
        cf_lower_angle_rl = new T;
        cf_lower_angle_rr = new T;

        cf_upper_S_fl = Math::calloc<T>(Constants::DIM);
        cf_upper_S_fr = Math::calloc<T>(Constants::DIM);
        cf_upper_S_rl = Math::calloc<T>(Constants::DIM);
        cf_upper_S_rr = Math::calloc<T>(Constants::DIM);
        cf_lower_S_fl = Math::calloc<T>(Constants::DIM);
        cf_lower_S_fr = Math::calloc<T>(Constants::DIM);
        cf_lower_S_rl = Math::calloc<T>(Constants::DIM);
        cf_lower_S_rr = Math::calloc<T>(Constants::DIM);

        cf_lower_rot_force_fl = Math::calloc<T>(Constants::DIM);
        cf_lower_rot_force_fr = Math::calloc<T>(Constants::DIM);
        cf_lower_rot_force_rl = Math::calloc<T>(Constants::DIM);
        cf_lower_rot_force_rr = Math::calloc<T>(Constants::DIM);
        cf_upper_rot_force_fl = Math::calloc<T>(Constants::DIM);
        cf_upper_rot_force_fr = Math::calloc<T>(Constants::DIM);
        cf_upper_rot_force_rl = Math::calloc<T>(Constants::DIM);
        cf_upper_rot_force_rr = Math::calloc<T>(Constants::DIM);
        cf_car_rot_force_fl = Math::calloc<T>(Constants::DIM);
        cf_car_rot_force_fr = Math::calloc<T>(Constants::DIM);
        cf_car_rot_force_rl = Math::calloc<T>(Constants::DIM);
        cf_car_rot_force_rr = Math::calloc<T>(Constants::DIM);
        cf_sum_car_force_fl = Math::calloc<T>(Constants::DIM);
        cf_sum_car_force_fr = Math::calloc<T>(Constants::DIM);
        cf_sum_car_force_rl = Math::calloc<T>(Constants::DIM);
        cf_sum_car_force_rr = Math::calloc<T>(Constants::DIM);

        cf_local_FR_fl = Math::calloc<T>(Constants::DIM);
        cf_local_FR_fr = Math::calloc<T>(Constants::DIM);
        cf_local_FR_rl = Math::calloc<T>(Constants::DIM);
        cf_local_FR_rr = Math::calloc<T>(Constants::DIM);

        cf_Hc = Math::calloc<T>(Constants::DIM);
        cf_sum_torque_spring_car = Math::calloc<T>(Constants::DIM);
        cf_Tc = Math::calloc<T>(Constants::DIM);
        cf_wc_tilda = Math::calloc<T>(Constants::DIM * Constants::DIM);

        cf_Tc[0] = 0;
        cf_Tc[1] = 0;
        cf_Tc[2] = 0;

        cf_b_rem = Math::calloc<T>(Constants::DIM * Constants::DIM * Constants::DIM);
        cf_Qc = Math::calloc<T>(Constants::NUM_LEGS * Constants::DIM);
        cf_qc_dot = Math::calloc<T>(Constants::NUM_LEGS);

        // required for the look up table
        current_spring_lengths = Math::calloc<T>(2 * Constants::NUM_LEGS * Constants::DIM);
        stiffness_vector = Math::calloc<T>(2 * Constants::NUM_LEGS * Constants::DIM);
    }

    /**
     * Clean the memory allocated for the main solver
     */
    void compute_f_clean() {
        // TODO: consider removing after checking performance.
        mkl_free_buffers();

        Math::free<T>(cf_C_cN);
        Math::free<T>(cf_r_up_fl);
        Math::free<T>(cf_r_up_fr);
        Math::free<T>(cf_r_up_rl);
        Math::free<T>(cf_r_up_rr);
        Math::free<T>(cf_r_low_fl);
        Math::free<T>(cf_r_low_fr);
        Math::free<T>(cf_r_low_rl);
        Math::free<T>(cf_r_low_rr);
        Math::free<T>(cf_upper_normal_fl);
        Math::free<T>(cf_upper_normal_fr);
        Math::free<T>(cf_upper_normal_rl);
        Math::free<T>(cf_upper_normal_rr);
        Math::free<T>(cf_lower_normal_fl);
        Math::free<T>(cf_lower_normal_fr);
        Math::free<T>(cf_lower_normal_rl);
        Math::free<T>(cf_lower_normal_rr);
        Math::free<T>(cf_col_dat);
        Math::free<T>(cf_upper_force_fl);
        Math::free<T>(cf_upper_force_fr);
        Math::free<T>(cf_upper_force_rl);
        Math::free<T>(cf_upper_force_rr);
        Math::free<T>(cf_lower_force_fl);
        Math::free<T>(cf_lower_force_fr);
        Math::free<T>(cf_lower_force_rl);
        Math::free<T>(cf_lower_force_rr);
        Math::free<T>(cf_upper_dampf_fl);
        Math::free<T>(cf_upper_dampf_fr);
        Math::free<T>(cf_upper_dampf_rl);
        Math::free<T>(cf_upper_dampf_rr);
        Math::free<T>(cf_lower_dampf_fl);
        Math::free<T>(cf_lower_dampf_fr);
        Math::free<T>(cf_lower_dampf_rl);
        Math::free<T>(cf_lower_dampf_rr);
        Math::free<T>(cf_temp);
        delete cf_upper_angle_fl;
        delete cf_upper_angle_fr;
        delete cf_upper_angle_rl;
        delete cf_upper_angle_rr;
        delete cf_lower_angle_fl;
        delete cf_lower_angle_fr;
        delete cf_lower_angle_rl;
        delete cf_lower_angle_rr;
        Math::free<T>(cf_upper_S_fl);
        Math::free<T>(cf_upper_S_fr);
        Math::free<T>(cf_upper_S_rl);
        Math::free<T>(cf_upper_S_rr);
        Math::free<T>(cf_lower_S_fl);
        Math::free<T>(cf_lower_S_fr);
        Math::free<T>(cf_lower_S_rl);
        Math::free<T>(cf_lower_S_rr);
        Math::free<T>(cf_lower_rot_force_fl);
        Math::free<T>(cf_lower_rot_force_fr);
        Math::free<T>(cf_lower_rot_force_rl);
        Math::free<T>(cf_lower_rot_force_rr);
        Math::free<T>(cf_upper_rot_force_fl);
        Math::free<T>(cf_upper_rot_force_fr);
        Math::free<T>(cf_upper_rot_force_rl);
        Math::free<T>(cf_upper_rot_force_rr);
        Math::free<T>(cf_car_rot_force_fl);
        Math::free<T>(cf_car_rot_force_fr);
        Math::free<T>(cf_car_rot_force_rl);
        Math::free<T>(cf_car_rot_force_rr);
        Math::free<T>(cf_sum_car_force_fl);
        Math::free<T>(cf_sum_car_force_fr);
        Math::free<T>(cf_sum_car_force_rl);
        Math::free<T>(cf_sum_car_force_rr);
        Math::free<T>(cf_local_FR_fl);
        Math::free<T>(cf_local_FR_fr);
        Math::free<T>(cf_local_FR_rl);
        Math::free<T>(cf_local_FR_rr);
        Math::free<T>(cf_Hc);
        Math::free<T>(cf_sum_torque_spring_car);
        Math::free<T>(cf_Tc);
        Math::free<T>(cf_wc_tilda);
        Math::free<T>(cf_b_rem);
        Math::free<T>(cf_Qc);
        Math::free<T>(cf_qc_dot);
        Math::free<T>(current_spring_lengths);
        Math::free<T>(stiffness_vector);
        Math::free<T>(accelerations);
    }

    void get_current_variables(T* x) {
        wc_ = x;
        vc_ = x + 3;
        vw_fl_ = x + 6;
        vw_fr_ = x + 9;
        vw_rl_ = x + 12;
        vw_rr_ = x + 15;
        vt_fl_ = x + 18;
        vt_fr_ = x + 21;
        vt_rl_ = x + 24;
        vt_rr_ = x + 27;
        qc_ = x + 30;
        pcc_ = x + 34;
        pw_fl_ = x + 37;
        pw_fr_ = x + 40;
        pw_rl_ = x + 43;
        pw_rr_ = x + 46;
        pt_fl_ = x + 49;
        pt_fr_ = x + 52;
        pt_rl_ = x + 55;
        pt_rr_ = x + 58;
    }

    void compute_spring_lengths(T* pcc_, T* pw_, T* pt_, T* cf_r_up_, T* cf_r_low_, T* r_,
                                T& norm_r_up, T& inv_norm_r_up, T& norm_r_low, T& inv_norm_r_low) {
        /*
        % global positions of the tyre connections
                pc1 = C_cN' * r1 + pcc;
                pc2 = C_cN' * r2 + pcc;
                pc3 = C_cN' * r3 + pcc;
                pc4 = C_cN' * r4 + pcc;

                % currrent length of the legs
                r_up1 = pc1 - pw1;
                r_up2 = pc2 - pw2;
                r_up3 = pc3 - pw3;
                r_up4 = pc4 - pw4;

                r_low1 = pw1 - pt1;
                r_low2 = pw2 - pt2;
                r_low3 = pw3 - pt3;
                r_low4 = pw4 - pt4;
        */

        Math::copy<T>(Constants::DIM, pcc_, 1, cf_r_up_, 1);
        Math::gemv<T>(CblasRowMajor, CblasNoTrans, Constants::DIM, Constants::DIM, 1, cf_C_cN,
                      Constants::DIM, r_, 1, 1, cf_r_up_, 1);
        Math::axpy<T>(Constants::DIM, -1, pw_, 1, cf_r_up_, 1);

        // r_low
        Math::copy<T>(Constants::DIM, pw_, 1, cf_r_low_, 1);
        Math::axpy<T>(Constants::DIM, -1, pt_, 1, cf_r_low_, 1);

        /* Compute spring lengths and their inverses */
        norm_r_up = Math::nrm2<T>(Constants::DIM, cf_r_up_, 1);

        inv_norm_r_up = 1. / norm_r_up;

        norm_r_low = Math::nrm2<T>(Constants::DIM, cf_r_low_, 1);

        inv_norm_r_low = 1. / norm_r_low;
    }

    void get_car_orientation() {
        /*
        // Basis //
        get cosine transforms (C_Nc means r_N = C_Nc * r_c)
        compute local base vectors
        basis_c = get_basis(qc);
        */
        Math::get_basis<T>(qc_, cf_C_cN);
    }

    void apply_stiffness_interpolation() {
        /* Compute stiffness from lookup table*/
#ifdef INTERPOLATION
        // populate the lenght_vector
        current_spring_lengths[0] = norm_r_up_fl;
        current_spring_lengths[1] = norm_r_low_fl;
        current_spring_lengths[2] = norm_r_up_fr;
        current_spring_lengths[3] = norm_r_low_fr;
        current_spring_lengths[4] = norm_r_up_rl;
        current_spring_lengths[5] = norm_r_low_rl;
        current_spring_lengths[6] = norm_r_up_rr;
        current_spring_lengths[7] = norm_r_low_rr;
        // calculate the new stiffnesses
        MetaDataBase<T>::getDataBase().getLookupStiffness().getInterpolation(current_spring_lengths,
                                                                             stiffness_vector);

        // overwrite stiffness values
        upper_spring_stiffness[0] = stiffness_vector[0];
        lower_spring_stiffness[0] = stiffness_vector[1];
        upper_spring_stiffness[1] = stiffness_vector[2];
        lower_spring_stiffness[1] = stiffness_vector[3];
        upper_spring_stiffness[2] = stiffness_vector[4];
        lower_spring_stiffness[2] = stiffness_vector[5];
        upper_spring_stiffness[3] = stiffness_vector[6];
        lower_spring_stiffness[3] = stiffness_vector[7];
#endif  // Interpolation compile flag, PLEASE USE THIS ONE ONLY !!!!!!!!!!!!!!!!!!!!!!!!!
    }

    void compute_angles() {
        /*
        get angle and normal vectors at the legs

        [~, upper_angle1, upper_normal1] = get_quaternion(r_up1, C_cN(2,:)');
        [~, upper_angle2, upper_normal2] = get_quaternion(r_up2, C_cN(2,:)');
        [~, upper_angle3, upper_normal3] = get_quaternion(r_up3, C_cN(2,:)');
        [~, upper_angle4, upper_normal4] = get_quaternion(r_up4, C_cN(2,:)');

        [~, lower_angle1, lower_normal1] = get_quaternion(r_low1, r_up1);
        [~, lower_angle2, lower_normal2] = get_quaternion(r_low2, r_up2);
        [~, lower_angle3, lower_normal3] = get_quaternion(r_low3, r_up3);
        [~, lower_angle4, lower_normal4] = get_quaternion(r_low4, r_up4);

        */

        Math::copy<T>(Constants::DIM, cf_C_cN + 2, Constants::DIM, cf_col_dat, 1);

        Math::get_quaternion<T>(cf_r_up_fl, cf_col_dat, cf_upper_angle_fl, cf_upper_normal_fl,
                                Constants::DIM);
        Math::get_quaternion<T>(cf_r_up_fr, cf_col_dat, cf_upper_angle_fr, cf_upper_normal_fr,
                                Constants::DIM);
        Math::get_quaternion<T>(cf_r_up_rl, cf_col_dat, cf_upper_angle_rl, cf_upper_normal_rl,
                                Constants::DIM);
        Math::get_quaternion<T>(cf_r_up_rr, cf_col_dat, cf_upper_angle_rr, cf_upper_normal_rr,
                                Constants::DIM);

        Math::get_quaternion<T>(cf_r_low_fl, cf_r_up_fl, cf_lower_angle_fl, cf_lower_normal_fl,
                                Constants::DIM);
        Math::get_quaternion<T>(cf_r_low_fr, cf_r_up_fr, cf_lower_angle_fr, cf_lower_normal_fr,
                                Constants::DIM);
        Math::get_quaternion<T>(cf_r_low_rl, cf_r_up_rl, cf_lower_angle_rl, cf_lower_normal_rl,
                                Constants::DIM);
        Math::get_quaternion<T>(cf_r_low_rr, cf_r_up_rr, cf_lower_angle_rr, cf_lower_normal_rr,
                                Constants::DIM);
    }

    void compute_elongational_forces() {
        // Forces and Torques
        // calculate the elongational spring forces (in global basis)

        /*
        upper_force1 =
            upper_spring_stiffness(1) * (r_up1) * (1 - upper_spring_length(1) * inv_norm_r_up1);
        upper_force2 =
            upper_spring_stiffness(2) * (r_up2) * (1 - upper_spring_length(2) * inv_norm_r_up2);
        upper_force3 =
            upper_spring_stiffness(3) * (r_up3) * (1 - upper_spring_length(3) * inv_norm_r_up3);
        upper_force4 =
            upper_spring_stiffness(4) * (r_up4) * (1 - upper_spring_length(4) * inv_norm_r_up4);

        lower_force1 =
            lower_spring_stiffness(1) * (r_low1) * (1 - lower_spring_length(1) * inv_norm_r_low1);
        lower_force2 =
            lower_spring_stiffness(2) * (r_low2) * (1 - lower_spring_length(2) * inv_norm_r_low2);
        lower_force3 =
            lower_spring_stiffness(3) * (r_low3) * (1 - lower_spring_length(3) * inv_norm_r_low3);
        lower_force4 =
            lower_spring_stiffness(4) * (r_low4) * (1 - lower_spring_length(4) * inv_norm_r_low4);
        */

        T scale;
        Math::copy<T>(Constants::DIM, cf_r_up_fl, 1, cf_upper_force_fl, 1);
        scale = this->upper_spring_stiffness[0] *
                (1. - this->upper_spring_length[0] * inv_norm_r_up_fl);
        Math::scal<T>(Constants::DIM, scale, cf_upper_force_fl, 1);
        Math::copy<T>(Constants::DIM, cf_r_up_fr, 1, cf_upper_force_fr, 1);
        scale = this->upper_spring_stiffness[1] *
                (1. - this->upper_spring_length[1] * inv_norm_r_up_fr);
        Math::scal<T>(Constants::DIM, scale, cf_upper_force_fr, 1);
        Math::copy<T>(Constants::DIM, cf_r_up_rl, 1, cf_upper_force_rl, 1);
        scale = this->upper_spring_stiffness[2] *
                (1. - this->upper_spring_length[2] * inv_norm_r_up_rl);
        Math::scal<T>(Constants::DIM, scale, cf_upper_force_rl, 1);
        Math::copy<T>(Constants::DIM, cf_r_up_rr, 1, cf_upper_force_rr, 1);
        scale = this->upper_spring_stiffness[3] *
                (1. - this->upper_spring_length[3] * inv_norm_r_up_rr);
        Math::scal<T>(Constants::DIM, scale, cf_upper_force_rr, 1);

        Math::copy<T>(Constants::DIM, cf_r_low_fl, 1, cf_lower_force_fl, 1);
        scale = this->lower_spring_stiffness[0] *
                (1. - this->lower_spring_length[0] * inv_norm_r_low_fl);
        Math::scal<T>(Constants::DIM, scale, cf_lower_force_fl, 1);
        Math::copy<T>(Constants::DIM, cf_r_low_fr, 1, cf_lower_force_fr, 1);
        scale = this->lower_spring_stiffness[1] *
                (1. - this->lower_spring_length[1] * inv_norm_r_low_fr);
        Math::scal<T>(Constants::DIM, scale, cf_lower_force_fr, 1);
        Math::copy<T>(Constants::DIM, cf_r_low_rl, 1, cf_lower_force_rl, 1);
        scale = this->lower_spring_stiffness[2] *
                (1. - this->lower_spring_length[2] * inv_norm_r_low_rl);
        Math::scal<T>(Constants::DIM, scale, cf_lower_force_rl, 1);
        Math::copy<T>(Constants::DIM, cf_r_low_rr, 1, cf_lower_force_rr, 1);
        scale = this->lower_spring_stiffness[3] *
                (1. - this->lower_spring_length[3] * inv_norm_r_low_rr);
        Math::scal<T>(Constants::DIM, scale, cf_lower_force_rr, 1);
    }

    void compute_damping_forces() {
        /*
        calculate forces from damping effects
        upper_vdiff1 = (dot((vc - C_cN' * (r1_tilda * wc)), r_up1) - dot(vw1, r_up1)) * r_up1 *
        inv_norm_r_up1 * inv_norm_r_up1; upper_vdiff2 = (dot((vc - C_cN' * (r2_tilda * wc)), r_up2)
        - dot(vw2, r_up2)) * r_up2 * inv_norm_r_up2 * inv_norm_r_up2; upper_vdiff3 = (dot((vc -
        C_cN' * (r3_tilda * wc)), r_up3) - dot(vw3, r_up3)) * r_up3 * inv_norm_r_up3 *
        inv_norm_r_up3; upper_vdiff4 = (dot((vc - C_cN' * (r4_tilda * wc)), r_up4) - dot(vw4,
        r_up4)) * r_up4 * inv_norm_r_up4 * inv_norm_r_up4;

        lower_vdiff1 = (dot(vw1, r_low1) - dot(vt1, r_low1)) * r_low1 * inv_norm_r_low1 *
        inv_norm_r_low1; lower_vdiff2 = (dot(vw2, r_low2) - dot(vt2, r_low2)) * r_low2 *
        inv_norm_r_low2 * inv_norm_r_low2; lower_vdiff3 = (dot(vw3, r_low3) - dot(vt3, r_low3)) *
        r_low3 * inv_norm_r_low3 * inv_norm_r_low3; lower_vdiff4 = (dot(vw4, r_low4) - dot(vt4,
        r_low4)) * r_low4 * inv_norm_r_low4 * inv_norm_r_low4;

        upper_dampf1 = upper_spring_damping1(upper_vdiff1);
        upper_dampf2 = upper_spring_damping2(upper_vdiff2);
        upper_dampf3 = upper_spring_damping3(upper_vdiff3);
        upper_dampf4 = upper_spring_damping4(upper_vdiff4);

        lower_dampf1 = lower_spring_damping1(lower_vdiff1);
        lower_dampf2 = lower_spring_damping2(lower_vdiff2);
        lower_dampf3 = lower_spring_damping3(lower_vdiff3);
        lower_dampf4 = lower_spring_damping4(lower_vdiff4);

*/

        const MKL_INT mkl_DIM = Constants::DIM;
        const MKL_INT mkl_incx = 1;
        const MKL_INT mkl_incy = 1;

        T scale;

        // upper_dampf1
        // compute: vc - C_cN' * (r1_tilda * wc))
        Math::gemv<T>(CblasRowMajor, CblasNoTrans, Constants::DIM, Constants::DIM, 1,
                      this->r_fl_tilda, Constants::DIM, wc_, 1, 0, cf_temp, 1);
        Math::copy<T>(Constants::DIM, vc_, 1, cf_upper_dampf_fl, 1);
        Math::gemv<T>(CblasRowMajor, CblasNoTrans, Constants::DIM, Constants::DIM, -1, cf_C_cN,
                      Constants::DIM, cf_temp, 1, 1, cf_upper_dampf_fl, 1);
        // dot((vc - C_cN' * (r1_tilda * wc)), r_up1)
        // std::cout << Math::dot<T>(Constants::DIM, upper_dampf1, 1, r_up1, 1) << std::endl;
        // std::cout << Math::dot_product<T>(upper_dampf1, r_up1, Constants::DIM) <<
        // std::endl;
        scale = Math::dot<T>(mkl_DIM, cf_upper_dampf_fl, mkl_incx, cf_r_up_fl, mkl_incy);
        // scale = Math::dot_product<T>(upper_dampf1, r_up1, Constants::DIM);

        // dot((vc - C_cN' * (r1_tilda * wc)), r_up1) - dot(vw1, r_up1)
        scale -= Math::dot<T>(mkl_DIM, vw_fl_, mkl_incx, cf_r_up_fl, mkl_incx);
        // scale -= Math::dot_product<T>(vw1_, r_up1, Constants::DIM);
        // (dot((vc - C_cN' * (r1_tilda * wc)), r_up1) - dot(vw1, r_up1))* inv_norm_r_up1 *
        // inv_norm_r_up1
        scale = scale * inv_norm_r_up_fl * inv_norm_r_up_fl * this->upper_spring_damping[0];
        Math::copy<T>(Constants::DIM, cf_r_up_fl, 1, cf_upper_dampf_fl, 1);
        Math::scal<T>(Constants::DIM, scale, cf_upper_dampf_fl, 1);

        //// upper_dampf_fr
        // compute: vc - C_cN' * (r_fr_tilda * wc))
        Math::gemv<T>(CblasRowMajor, CblasNoTrans, Constants::DIM, Constants::DIM, 1,
                      this->r_fr_tilda, Constants::DIM, wc_, 1, 0, cf_temp, 1);
        Math::copy<T>(Constants::DIM, vc_, 1, cf_upper_dampf_fr, 1);
        Math::gemv<T>(CblasRowMajor, CblasNoTrans, Constants::DIM, Constants::DIM, -1, cf_C_cN,
                      Constants::DIM, cf_temp, 1, 1, cf_upper_dampf_fr, 1);
        // dot((vc - C_cN' * (r_fr_tilda * wc)), r_up_fr)
        scale = Math::dot<T>(mkl_DIM, cf_upper_dampf_fr, mkl_incx, cf_r_up_fr, mkl_incy);
        // scale = Math::dot_product<T>(upper_dampf_fr, r_up_fr, Constants::DIM);
        // dot((vc - C_cN' * (r_fr_tilda * wc)), r_up_fr) - dot(vw_fr, r_up_fr)
        scale -= Math::dot<T>(mkl_DIM, vw_fr_, mkl_incx, cf_r_up_fr, mkl_incy);
        // scale -= Math::dot_product<T>(vw_fr_, r_up_fr, Constants::DIM);
        // (dot((vc - C_cN' * (r_fr_tilda * wc)), r_up_fr) - dot(vw_fr, r_up_fr))* inv_norm_r_up_fr
        // * inv_norm_r_up_fr
        scale = scale * inv_norm_r_up_fr * inv_norm_r_up_fr * this->upper_spring_damping[1];
        Math::copy<T>(Constants::DIM, cf_r_up_fr, 1, cf_upper_dampf_fr, 1);
        Math::scal<T>(Constants::DIM, scale, cf_upper_dampf_fr, 1);

        //// upper_dampf3
        // compute: vc - C_cN' * (r3_tilda * wc))
        Math::gemv<T>(CblasRowMajor, CblasNoTrans, Constants::DIM, Constants::DIM, 1,
                      this->r_rl_tilda, Constants::DIM, wc_, 1, 0, cf_temp, 1);
        Math::copy<T>(Constants::DIM, vc_, 1, cf_upper_dampf_rl, 1);
        Math::gemv<T>(CblasRowMajor, CblasNoTrans, Constants::DIM, Constants::DIM, -1, cf_C_cN,
                      Constants::DIM, cf_temp, 1, 1, cf_upper_dampf_rl, 1);
        // dot((vc - C_cN' * (r_rl_tilda * wc)), r_up_rl)
        scale = Math::dot<T>(mkl_DIM, cf_upper_dampf_rl, mkl_incx, cf_r_up_rl, mkl_incy);
        // scale = Math::dot_product<T>(upper_dampf_rl, r_up_rl, Constants::DIM);
        // dot((vc - C_cN' * (r_rl_tilda * wc)), r_up_rl) - dot(vw_rl, r_up_rl)
        scale -= Math::dot<T>(mkl_DIM, vw_rl_, mkl_incx, cf_r_up_rl, mkl_incy);
        // scale -= Math::dot_product<T>(vw_rl_, r_up_rl, Constants::DIM);
        // (dot((vc - C_cN' * (r_rl_tilda * wc)), r_up_rl) - dot(vw_rl, r_up_rl))* inv_norm_r_up_rl
        // * inv_norm_r_up_rl
        scale = scale * inv_norm_r_up_rl * inv_norm_r_up_rl * this->upper_spring_damping[2];
        Math::copy<T>(Constants::DIM, cf_r_up_rl, 1, cf_upper_dampf_rl, 1);
        Math::scal<T>(Constants::DIM, scale, cf_upper_dampf_rl, 1);

        //// upper_dampf4
        // compute: vc - C_cN' * (r4_tilda * wc))
        Math::gemv<T>(CblasRowMajor, CblasNoTrans, Constants::DIM, Constants::DIM, 1,
                      this->r_rr_tilda, Constants::DIM, wc_, 1, 0, cf_temp, 1);
        Math::copy<T>(Constants::DIM, vc_, 1, cf_upper_dampf_rr, 1);
        Math::gemv<T>(CblasRowMajor, CblasNoTrans, Constants::DIM, Constants::DIM, -1, cf_C_cN,
                      Constants::DIM, cf_temp, 1, 1, cf_upper_dampf_rr, 1);
        // dot((vc - C_cN' * (r4_tilda * wc)), r_up4)
        scale = Math::dot<T>(mkl_DIM, cf_upper_dampf_rr, mkl_incx, cf_r_up_rr, mkl_incy);
        // scale = Math::dot_product<T>(upper_dampf4, r_up4, Constants::DIM);
        // dot((vc - C_cN' * (r4_tilda * wc)), r_up4) - dot(vw4, r_up4)
        scale -= Math::dot<T>(mkl_DIM, vw_rr_, mkl_incx, cf_r_up_rr, mkl_incy);
        // scale -= Math::dot_product<T>(r_up4, vw4_, Constants::DIM);
        // (dot((vc - C_cN' * (r4_tilda * wc)), r_up4) - dot(vw4, r_up4))* inv_norm_r_up4 *
        // inv_norm_r_up4
        scale = scale * inv_norm_r_up_rr * inv_norm_r_up_rr * this->upper_spring_damping[3];
        Math::copy<T>(Constants::DIM, cf_r_up_rr, 1, cf_upper_dampf_rr, 1);
        Math::scal<T>(Constants::DIM, scale, cf_upper_dampf_rr, 1);

        //// lower_dampf1
        // dot(vw1, r_low1)
        scale = Math::dot<T>(mkl_DIM, vw_fl_, mkl_incx, cf_r_low_fl, mkl_incy);
        /*scale = Math::dot_product<T>(vw1_, r_low1, Constants::DIM);*/
        // (dot(vw1, r_low1) - dot(vt1, r_low1))
        scale -= Math::dot<T>(mkl_DIM, vt_fl_, mkl_incx, cf_r_low_fl, mkl_incx);
        // scale -= Math::dot_product<T>(vt1_, r_low1, Constants::DIM);
        //(dot(vw1, r_low1) - dot(vt1, r_low1)) * inv_norm_r_low1 * inv_norm_r_low1
        scale = scale * inv_norm_r_low_fl * inv_norm_r_low_fl * this->lower_spring_damping[0];
        Math::copy<T>(Constants::DIM, cf_r_low_fl, 1, cf_lower_dampf_fl, 1);
        Math::scal<T>(Constants::DIM, scale, cf_lower_dampf_fl, 1);

        //// lower_dampf2
        // dot(vw2, r_low2)
        scale = Math::dot<T>(mkl_DIM, vw_fr_, mkl_incx, cf_r_low_fr, mkl_incy);
        // scale = Math::dot_product<T>(vw2_, r_low2, Constants::DIM);
        // (dot(vw2, r_low2) - dot(vt2, r_low2))
        scale -= Math::dot<T>(mkl_DIM, vt_fr_, mkl_incx, cf_r_low_fr, mkl_incy);
        // scale -= Math::dot_product<T>(vt2_, r_low2, Constants::DIM);
        //(dot(vw2, r_low2) - dot(vt2, r_low2)) * inv_norm_r_low2 * inv_norm_r_low2
        scale = scale * inv_norm_r_low_fr * inv_norm_r_low_fr * this->lower_spring_damping[1];
        Math::copy<T>(Constants::DIM, cf_r_low_fr, 1, cf_lower_dampf_fr, 1);
        Math::scal<T>(Constants::DIM, scale, cf_lower_dampf_fr, 1);

        //// lower_dampf3
        // dot(vw3, r_low3)
        scale = Math::dot<T>(mkl_DIM, vw_rl_, mkl_incx, cf_r_low_rl, mkl_incy);
        // scale = Math::dot_product<T>(vw_rl_, r_low_rl, Constants::DIM);
        // (dot(vw_rl, r_low_rl) - dot(vt_rl, r_low_rl))
        scale -= Math::dot<T>(mkl_DIM, vt_rl_, mkl_incx, cf_r_low_rl, mkl_incy);
        // scale -= Math::dot_product<T>(vt_rl_, r_low_rl, Constants::DIM);
        //(dot(vw_rl, r_low_rl) - dot(vt_rl, r_low_rl)) * inv_norm_r_low_rl * inv_norm_r_low_rl
        scale = scale * inv_norm_r_low_rl * inv_norm_r_low_rl * this->lower_spring_damping[2];
        Math::copy<T>(Constants::DIM, cf_r_low_rl, 1, cf_lower_dampf_rl, 1);
        Math::scal<T>(Constants::DIM, scale, cf_lower_dampf_rl, 1);

        //// lower_dampf4
        // dot(vw4, r_low4)
        scale = Math::dot<T>(mkl_DIM, vw_rr_, mkl_incx, cf_r_low_rr, mkl_incy);
        // scale = Math::dot_product<T>(vw4_, r_low4, Constants::DIM);
        // (dot(vw4, r_low4) - dot(vt4, r_low4))
        scale -= Math::dot<T>(mkl_DIM, vt_rr_, mkl_incx, cf_r_low_rr, mkl_incy);
        // scale -= Math::dot_product<T>(vt4_, r_low4, Constants::DIM);
        //(dot(vw4, r_low4) - dot(vt4, r_low4)) * inv_norm_r_low4 * inv_norm_r_low4
        scale = scale * inv_norm_r_low_rr * inv_norm_r_low_rr * this->lower_spring_damping[3];
        Math::copy<T>(Constants::DIM, cf_r_low_rr, 1, cf_lower_dampf_rr, 1);
        Math::scal<T>(Constants::DIM, scale, cf_lower_dampf_rr, 1);
    }

    void compute_torques() {
        /*
torque from the rotational spring
upper_S1 = upper_rotational_stiffness(1) * upper_angle1 * upper_normal1;        % in global basis
upper_S2 = upper_rotational_stiffness(2) * upper_angle2 * upper_normal2;
upper_S3 = upper_rotational_stiffness(3) * upper_angle3 * upper_normal3;
upper_S4 = upper_rotational_stiffness(4) * upper_angle4 * upper_normal4;

lower_S1 = lower_rotational_stiffness(1) * lower_angle1 * lower_normal1;        % in global basis
lower_S2 = lower_rotational_stiffness(2) * lower_angle2 * lower_normal2;
lower_S3 = lower_rotational_stiffness(3) * lower_angle3 * lower_normal3;
lower_S4 = lower_rotational_stiffness(4) * lower_angle4 * lower_normal4;
*/

        T scale;

        // upper_S1 = upper_rotational_stiffness(1) * upper_angle1 * upper_normal1;
        scale = (this->upper_rotational_stiffness[0]) * (*cf_upper_angle_fl);
        Math::copy<T>(Constants::DIM, cf_upper_normal_fl, 1, cf_upper_S_fl, 1);
        Math::scal<T>(Constants::DIM, scale, cf_upper_S_fl, 1);

        // upper_S2 = upper_rotational_stiffness(2) * upper_angle2 * upper_normal2;
        scale = (this->upper_rotational_stiffness[1]) * (*cf_upper_angle_fr);
        Math::copy<T>(Constants::DIM, cf_upper_normal_fr, 1, cf_upper_S_fr, 1);
        Math::scal<T>(Constants::DIM, scale, cf_upper_S_fr, 1);

        // upper_S3 = upper_rotational_stiffness(3) * upper_angle3 * upper_normal3;
        scale = (this->upper_rotational_stiffness[2]) * (*cf_upper_angle_rl);
        Math::copy<T>(Constants::DIM, cf_upper_normal_rl, 1, cf_upper_S_rl, 1);
        Math::scal<T>(Constants::DIM, scale, cf_upper_S_rl, 1);

        // upper_S4 = upper_rotational_stiffness(4) * upper_angle4 * upper_normal4;
        scale = (this->upper_rotational_stiffness[3]) * (*cf_upper_angle_rr);
        Math::copy<T>(Constants::DIM, cf_upper_normal_rr, 1, cf_upper_S_rr, 1);
        Math::scal<T>(Constants::DIM, scale, cf_upper_S_rr, 1);

        // lower_S1 = lower_rotational_stiffness(1) * lower_angle1 * lower_normal1
        scale = (this->lower_rotational_stiffness[0]) * (*cf_lower_angle_fl);
        Math::copy<T>(Constants::DIM, cf_lower_normal_fl, 1, cf_lower_S_fl, 1);
        Math::scal<T>(Constants::DIM, scale, cf_lower_S_fl, 1);

        // lower_S2 = lower_rotational_stiffness(2) * lower_angle2 * lower_normal2
        scale = (this->lower_rotational_stiffness[1]) * (*cf_lower_angle_fr);
        Math::copy<T>(Constants::DIM, cf_lower_normal_fr, 1, cf_lower_S_fr, 1);
        Math::scal<T>(Constants::DIM, scale, cf_lower_S_fr, 1);

        // lower_S3 = lower_rotational_stiffness(3) * lower_angle3 * lower_normal3
        scale = (this->lower_rotational_stiffness[2]) * (*cf_lower_angle_rl);
        Math::copy<T>(Constants::DIM, cf_lower_normal_rl, 1, cf_lower_S_rl, 1);
        Math::scal<T>(Constants::DIM, scale, cf_lower_S_rl, 1);

        // lower_S_rr = lower_rotational_stiffness(4) * lower_angle4 * lower_normal4
        scale = (this->lower_rotational_stiffness[3]) * (*cf_lower_angle_rr);
        Math::copy<T>(Constants::DIM, cf_lower_normal_rr, 1, cf_lower_S_rr, 1);
        Math::scal<T>(Constants::DIM, scale, cf_lower_S_rr, 1);
    }

    void compute_torques_resultant_forces() {
        /*
         * calculate the effect of the rotational springs (in global basis)
         *
         * lower_rot_force1 = -cross(lower_S1, r_low1) / (r_low1'*r_low1);
         * lower_rot_force2 = -cross(lower_S2, r_low2) / (r_low2'*r_low2);
         * lower_rot_force3 = -cross(lower_S3, r_low3) / (r_low3'*r_low3);
         * lower_rot_force4 = -cross(lower_S4, r_low4) / (r_low4'*r_low4);
         *
         * upper_rot_force1 = -cross(upper_S1, r_up1) / (r_up1'*r_up1);
         * upper_rot_force2 = -cross(upper_S2, r_up2) / (r_up2'*r_up2);
         * upper_rot_force3 = -cross(upper_S3, r_up3) / (r_up3'*r_up3);
         * upper_rot_force4 = -cross(upper_S4, r_up4) / (r_up4'*r_up4);
         *
         * car_rot_force1 = -cross(lower_S1, r_up1) / (r_up1'*r_up1);
         * car_rot_force2 = -cross(lower_S2, r_up2) / (r_up2'*r_up2);
         * car_rot_force3 = -cross(lower_S3, r_up3) / (r_up3'*r_up3);
         * car_rot_force4 = -cross(lower_S4, r_up4) / (r_up4'*r_up4);
         *
         * sum_car_force1 = car_rot_force1 - upper_force1 - upper_dampf1 - upper_rot_force1;
         * sum_car_force2 = car_rot_force2 - upper_force2 - upper_dampf2 - upper_rot_force2;
         * sum_car_force3 = car_rot_force3 - upper_force3 - upper_dampf3 - upper_rot_force3;
         * sum_car_force4 = car_rot_force4 - upper_force4 - upper_dampf4 - upper_rot_force4;
         */
        // TODO: constants
        const MKL_INT mkl_DIM = Constants::DIM;
        const MKL_INT mkl_incx = 1;
        const MKL_INT mkl_incy = 1;

        T scale;
        T scale_u_fl, scale_u_fr, scale_u_rl, scale_u_rr;
        // lower_rot_force1 = -cross( lower_S1, r_low1) / (r_low1'*r_low1);
        scale = -1. / Math::dot<T>(mkl_DIM, cf_r_low_fl, mkl_incx, cf_r_low_fl, mkl_incy);
        // scale = -1. / Math::dot_product<T>(r_low1, r_low1, Constants::DIM);
        Math::crossProduct<T>(cf_lower_S_fl, cf_r_low_fl, cf_lower_rot_force_fl);
        Math::scal<T>(Constants::DIM, scale, cf_lower_rot_force_fl, 1);

        // lower_rot_force2 = -cross( lower_S2, r_low2) / (r_low2'*r_low2);
        scale = -1. / Math::dot<T>(mkl_DIM, cf_r_low_fr, mkl_incx, cf_r_low_fr, mkl_incy);
        // scale = -1. / Math::dot_product<T>(r_low_fr, r_low_fr, Constants::DIM);
        Math::crossProduct<T>(cf_lower_S_fr, cf_r_low_fr, cf_lower_rot_force_fr);
        Math::scal<T>(Constants::DIM, scale, cf_lower_rot_force_fr, 1);

        // lower_rot_force3 = -cross( lower_S3, r_low3) / (r_low3'*r_low3);
        scale = -1. / Math::dot<T>(mkl_DIM, cf_r_low_rl, mkl_incx, cf_r_low_rl, mkl_incy);
        // scale = -1. / Math::dot_product<T>(r_low3, r_low3, Constants::DIM);
        Math::crossProduct<T>(cf_lower_S_rl, cf_r_low_rl, cf_lower_rot_force_rl);
        Math::scal<T>(Constants::DIM, scale, cf_lower_rot_force_rl, 1);

        // lower_rot_force4 = -cross( lower_S4, r_low4) / (r_low4'*r_low4);
        scale = -1. / Math::dot<T>(mkl_DIM, cf_r_low_rr, mkl_incx, cf_r_low_rr, mkl_incx);
        // scale = -1. / Math::dot_product<T>(r_low4, r_low4, Constants::DIM);
        Math::crossProduct<T>(cf_lower_S_rr, cf_r_low_rr, cf_lower_rot_force_rr);
        Math::scal<T>(Constants::DIM, scale, cf_lower_rot_force_rr, 1);

        // upper_rot_force1 = -cross( upper_S1, r_up1) / (r_up1'*r_up1);
        scale_u_fl = -1. / Math::dot<T>(mkl_DIM, cf_r_up_fl, mkl_incx, cf_r_up_fl, mkl_incy);
        // scale_u1 = -1. / Math::dot_product<T>(r_up1, r_up1, Constants::DIM);
        Math::crossProduct<T>(cf_upper_S_fl, cf_r_up_fl, cf_upper_rot_force_fl);
        Math::scal<T>(Constants::DIM, scale_u_fl, cf_upper_rot_force_fl, 1);

        // upper_rot_force_fr = -cross( upper_S_fr, r_up_fr) / (r_up_fr'*r_up_fr);
        scale_u_fr = -1. / Math::dot<T>(mkl_DIM, cf_r_up_fr, mkl_incx, cf_r_up_fr, mkl_incy);
        // scale_u_fr = -1. / Math::dot_product<T>(r_up_fr, r_up_fr, Constants::DIM);
        Math::crossProduct<T>(cf_upper_S_fr, cf_r_up_fr, cf_upper_rot_force_fr);
        Math::scal<T>(Constants::DIM, scale_u_fr, cf_upper_rot_force_fr, 1);

        // upper_rot_force3 = -cross( upper_S3, r_up3) / (r_up3'*r_up3);
        scale_u_rl = -1. / Math::dot<T>(mkl_DIM, cf_r_up_rl, mkl_incx, cf_r_up_rl, mkl_incy);
        // scale_u3 = -1. / Math::dot_product<T>(r_up3, r_up3, Constants::DIM);
        Math::crossProduct<T>(cf_upper_S_rl, cf_r_up_rl, cf_upper_rot_force_rl);
        Math::scal<T>(Constants::DIM, scale_u_rl, cf_upper_rot_force_rl, 1);

        // upper_rot_force4 = -cross( upper_S4, r_up4) / (r_up4'*r_up4);
        scale_u_rr = -1. / Math::dot<T>(mkl_DIM, cf_r_up_rr, mkl_incx, cf_r_up_rr, mkl_incy);
        // scale_u4 = -1. / Math::dot_product<T>(r_up4, r_up4, Constants::DIM);
        Math::crossProduct<T>(cf_upper_S_rr, cf_r_up_rr, cf_upper_rot_force_rr);
        Math::scal<T>(Constants::DIM, scale_u_rr, cf_upper_rot_force_rr, 1);

        // car_rot_force1 = -cross( lower_S1, r_up1) / (r_up1'*r_up1);
        Math::crossProduct<T>(cf_lower_S_fl, cf_r_up_fl, cf_car_rot_force_fl);
        Math::scal<T>(Constants::DIM, scale_u_fl, cf_car_rot_force_fl, 1);

        // car_rot_force_fr = -cross( lower_S2, r_up2) / (r_up2'*r_up2);
        Math::crossProduct<T>(cf_lower_S_fr, cf_r_up_fr, cf_car_rot_force_fr);
        Math::scal<T>(Constants::DIM, scale_u_fr, cf_car_rot_force_fr, 1);

        // car_rot_force3 = -cross( lower_S3, r_up3) / (r_up3'*r_up3);
        Math::crossProduct<T>(cf_lower_S_rl, cf_r_up_rl, cf_car_rot_force_rl);
        Math::scal<T>(Constants::DIM, scale_u_rl, cf_car_rot_force_rl, 1);

        // car_rot_force4 = -cross( lower_S4, r_up4) / (r_up4'*r_up4);
        Math::crossProduct<T>(cf_lower_S_rr, cf_r_up_rr, cf_car_rot_force_rr);
        Math::scal<T>(Constants::DIM, scale_u_rr, cf_car_rot_force_rr, 1);
    }

    void compute_car_body_forces() {
        // sum_car_force1 = car_rot_force1 - upper_force1 - upper_dampf1 - upper_rot_force1;
        Math::copy<T>(Constants::DIM, cf_car_rot_force_fl, 1, cf_sum_car_force_fl, 1);
        Math::axpy<T>(Constants::DIM, -1, cf_upper_force_fl, 1, cf_sum_car_force_fl, 1);
        Math::axpy<T>(Constants::DIM, -1, cf_upper_dampf_fl, 1, cf_sum_car_force_fl, 1);
        Math::axpy<T>(Constants::DIM, -1, cf_upper_rot_force_fl, 1, cf_sum_car_force_fl, 1);

        // sum_car_force_fr = car_rot_force_fr - upper_force_fr - upper_dampf_fr -
        // upper_rot_force_fr;
        Math::copy<T>(Constants::DIM, cf_car_rot_force_fr, 1, cf_sum_car_force_fr, 1);
        Math::axpy<T>(Constants::DIM, -1, cf_upper_force_fr, 1, cf_sum_car_force_fr, 1);
        Math::axpy<T>(Constants::DIM, -1, cf_upper_dampf_fr, 1, cf_sum_car_force_fr, 1);
        Math::axpy<T>(Constants::DIM, -1, cf_upper_rot_force_fr, 1, cf_sum_car_force_fr, 1);

        // sum_car_force3 = car_rot_force3 - upper_force3 - upper_dampf3 - upper_rot_force3;
        Math::copy<T>(Constants::DIM, cf_car_rot_force_rl, 1, cf_sum_car_force_rl, 1);
        Math::axpy<T>(Constants::DIM, -1, cf_upper_force_rl, 1, cf_sum_car_force_rl, 1);
        Math::axpy<T>(Constants::DIM, -1, cf_upper_dampf_rl, 1, cf_sum_car_force_rl, 1);
        Math::axpy<T>(Constants::DIM, -1, cf_upper_rot_force_rl, 1, cf_sum_car_force_rl, 1);

        // sum_car_force4 = car_rot_force4 - upper_force4 - upper_dampf4 - upper_rot_force4;
        Math::copy<T>(Constants::DIM, cf_car_rot_force_rr, 1, cf_sum_car_force_rr, 1);
        Math::axpy<T>(Constants::DIM, -1, cf_upper_force_rr, 1, cf_sum_car_force_rr, 1);
        Math::axpy<T>(Constants::DIM, -1, cf_upper_dampf_rr, 1, cf_sum_car_force_rr, 1);
        Math::axpy<T>(Constants::DIM, -1, cf_upper_rot_force_rr, 1, cf_sum_car_force_rr, 1);
    }

    void compute_external_forces() {
        /*
         * get external forces for the legs
         *
         * local_FW1 = [0; FW(1); 0];      %in global basis
         * local_FW_fr = [0; FW(_fr); 0];
         * local_FW3 = [0; FW(3); 0];
         * local_FW4 = [0; FW(4); 0];
         *
         * local_FT1 = [0; FT(1); 0];      %in global basis
         * local_FT_fr = [0; FT(2); 0];
         * local_FT3 = [0; FT(3); 0];
         * local_FT4 = [0; FT(4); 0];
         */

        // Calculate the sum of all forces in the tyre (and keep it in local_FR)

        // local_FR1 = lower_force1 + lower_dampf1 + local_FT1 + local_FR1 + lower_rot_force1; ...
        // %vt1_dot
        Math::copy<T>(Constants::DIM, cf_lower_force_fl, 1, cf_local_FR_fl, 1);
        Math::axpy<T>(Constants::DIM, 1, cf_lower_dampf_fl, 1, cf_local_FR_fl, 1);
        Math::axpy<T>(Constants::DIM, 1, FT_fl, 1, cf_local_FR_fl, 1);
        Math::axpy<T>(Constants::DIM, 1, cf_local_FR_fl, 1, cf_local_FR_fl, 1);
        Math::axpy<T>(Constants::DIM, 1, cf_lower_rot_force_fl, 1, cf_local_FR_fl, 1);

        // local_FR2 = lower_force2 + lower_dampf2 + local_FT2 + local_FR2 + lower_rot_force2; ...
        // %vt2_dot
        Math::copy<T>(Constants::DIM, cf_lower_force_fr, 1, cf_local_FR_fr, 1);
        Math::axpy<T>(Constants::DIM, 1, cf_lower_dampf_fr, 1, cf_local_FR_fr, 1);
        Math::axpy<T>(Constants::DIM, 1, FT_fr, 1, cf_local_FR_fr, 1);
        Math::axpy<T>(Constants::DIM, 1, cf_local_FR_fr, 1, cf_local_FR_fr, 1);
        Math::axpy<T>(Constants::DIM, 1, cf_lower_rot_force_fr, 1, cf_local_FR_fr, 1);

        // local_FR3 = lower_force3 + lower_dampf3 + local_FT3 + local_FR3 + lower_rot_force3; ...
        // %vt3_dot
        Math::copy<T>(Constants::DIM, cf_lower_force_rl, 1, cf_local_FR_rl, 1);
        Math::axpy<T>(Constants::DIM, 1, cf_lower_dampf_rl, 1, cf_local_FR_rl, 1);
        Math::axpy<T>(Constants::DIM, 1, FT_rl, 1, cf_local_FR_rl, 1);
        Math::axpy<T>(Constants::DIM, 1, cf_local_FR_rl, 1, cf_local_FR_rl, 1);
        Math::axpy<T>(Constants::DIM, 1, cf_lower_rot_force_rl, 1, cf_local_FR_rl, 1);

        // local_FR4 = lower_force4 + lower_dampf4 + local_FT4 + local_FR4 + lower_rot_force4];
        // %vt4_dot
        Math::copy<T>(Constants::DIM, cf_lower_force_rr, 1, cf_local_FR_rr, 1);
        Math::axpy<T>(Constants::DIM, 1, cf_lower_dampf_rr, 1, cf_local_FR_rr, 1);
        Math::axpy<T>(Constants::DIM, 1, FT_rr, 1, cf_local_FR_rr, 1);
        Math::axpy<T>(Constants::DIM, 1, cf_local_FR_rr, 1, cf_local_FR_rr, 1);
        Math::axpy<T>(Constants::DIM, 1, cf_lower_rot_force_rr, 1, cf_local_FR_rr, 1);
    }

    void compute_car_body_total_torque() {
        /*
         * get H=I*w
         * Hc = A(1:3, 1:3) * wc;           %in local car basis
         *
         * external torque on the car body (use later for rotational damping)
         * Tc = zeros(3,1);     %in local car basis
         *
         * sum of all torques induced by forces
         * for the car body
         * sum_torque_spring_car = r1_tilda * (C_cN * sum_car_force1) + ... % from the elongational
         * springs r2_tilda * (C_cN * sum_car_force2) + ...
         *         r3_tilda * (C_cN * sum_car_force3) + ...
         *         r_rr_tilda * (C_cN * sum_car_force4) + ...
         *         -C_cN * upper_S1 - C_cN * upper_S2 - C_cN * upper_S3 - C_cN * upper_S4 + ...
         *           % ??from the rotational spring
         * -get_tilda(wc) * Hc + Tc;
         *    % from angular momentum and external torques
         */

        // Hc = A(1:3, 1:3) * wc;
        Math::gemv<T>(CblasRowMajor, CblasNoTrans, Constants::DIM, Constants::DIM, 1, this->Ic,
                      Constants::DIM, wc_, 1, 0, cf_Hc, 1);
        Math::get_tilda<T>(wc_, cf_wc_tilda);
        Math::copy<T>(Constants::DIM, cf_Hc, 1, cf_temp, 1);
        Math::gemv<T>(CblasRowMajor, CblasNoTrans, Constants::DIM, Constants::DIM, -1, cf_wc_tilda,
                      Constants::DIM, cf_temp, 1, 0, cf_Hc, 1);

        Math::copy<T>(Constants::DIM, cf_Hc, 1, cf_sum_torque_spring_car, 1);

        Math::axpy<T>(Constants::DIM, 1, cf_Tc, 1, cf_sum_torque_spring_car, 1);

        Math::gemv<T>(CblasRowMajor, CblasTrans, Constants::DIM, Constants::DIM, -1, cf_C_cN,
                      Constants::DIM, cf_upper_S_rr, 1, 1, cf_sum_torque_spring_car, 1);
        Math::gemv<T>(CblasRowMajor, CblasTrans, Constants::DIM, Constants::DIM, -1, cf_C_cN,
                      Constants::DIM, cf_upper_S_rl, 1, 1, cf_sum_torque_spring_car, 1);
        Math::gemv<T>(CblasRowMajor, CblasTrans, Constants::DIM, Constants::DIM, -1, cf_C_cN,
                      Constants::DIM, cf_upper_S_fr, 1, 1, cf_sum_torque_spring_car, 1);
        Math::gemv<T>(CblasRowMajor, CblasTrans, Constants::DIM, Constants::DIM, -1, cf_C_cN,
                      Constants::DIM, cf_upper_S_fl, 1, 1, cf_sum_torque_spring_car, 1);

        Math::gemv<T>(CblasRowMajor, CblasTrans, Constants::DIM, Constants::DIM, 1, cf_C_cN,
                      Constants::DIM, cf_sum_car_force_fl, 1, 0, cf_temp, 1);
        Math::gemv<T>(CblasRowMajor, CblasNoTrans, Constants::DIM, Constants::DIM, 1,
                      this->r_fl_tilda, Constants::DIM, cf_temp, 1, 1, cf_sum_torque_spring_car, 1);
        Math::gemv<T>(CblasRowMajor, CblasTrans, Constants::DIM, Constants::DIM, 1, cf_C_cN,
                      Constants::DIM, cf_sum_car_force_fr, 1, 0, cf_temp, 1);
        Math::gemv<T>(CblasRowMajor, CblasNoTrans, Constants::DIM, Constants::DIM, 1,
                      this->r_fr_tilda, Constants::DIM, cf_temp, 1, 1, cf_sum_torque_spring_car, 1);
        Math::gemv<T>(CblasRowMajor, CblasTrans, Constants::DIM, Constants::DIM, 1, cf_C_cN,
                      Constants::DIM, cf_sum_car_force_rl, 1, 0, cf_temp, 1);
        Math::gemv<T>(CblasRowMajor, CblasNoTrans, Constants::DIM, Constants::DIM, 1,
                      this->r_rl_tilda, Constants::DIM, cf_temp, 1, 1, cf_sum_torque_spring_car, 1);
        Math::gemv<T>(CblasRowMajor, CblasTrans, Constants::DIM, Constants::DIM, 1, cf_C_cN,
                      Constants::DIM, cf_sum_car_force_rr, 1, 0, cf_temp, 1);
        Math::gemv<T>(CblasRowMajor, CblasNoTrans, Constants::DIM, Constants::DIM, 1,
                      this->r_rr_tilda, Constants::DIM, cf_temp, 1, 1, cf_sum_torque_spring_car, 1);
    }

    void construct_right_hand_side() {
        /*
         * b = [sum_torque_spring_car; ...% w_dot_c
         * FC + sum_car_force1 + sum_car_force2 + sum_car_force3 + sum_car_force4; ...% vc_dot
         * upper_force1 - lower_force1 + upper_dampf1 - lower_dampf1 + local_FW1 + ...
         * upper_rot_force1 - car_rot_force1 - lower_rot_force1; ...% vw1_dot
         * upper_force2 - lower_force2 + upper_dampf2 - lower_dampf2 + local_FW2 + ...
         * upper_rot_force2 - car_rot_force2 - lower_rot_force2; ...% vw2_dot
         * upper_force3 - lower_force3 + upper_dampf3 - lower_dampf3 + local_FW3 + ...
         * upper_rot_force3 - car_rot_force3 - lower_rot_force3; ...% vw3_dot
         * upper_force4 - lower_force4 + upper_dampf4 - lower_dampf4 + local_FW4 + ...
         * upper_rot_force4 - car_rot_force4 - lower_rot_force4; ...% vw4_dot
         * lower_force1 + lower_dampf1 + local_FT1 + local_FR1 + lower_rot_force1; ...% vt1_dot
         * lower_force2 + lower_dampf2 + local_FT2 + local_FR2 + lower_rot_force2; ...% vt2_dot
         * lower_force3 + lower_dampf3 + local_FT3 + local_FR3 + lower_rot_force3; ...% vt3_dot
         * lower_force4 + lower_dampf4 + local_FT4 + local_FR4 + lower_rot_force4];% vt4_dot
         */
        // FC + sum_car_force1 + sum_car_force2 + sum_car_force3 + sum_car_force4; ... %vc_dot
        T* brem_start = cf_b_rem;
        Math::copy<T>(Constants::DIM, this->FC, 1, brem_start, 1);
        Math::axpy<T>(Constants::DIM, 1, cf_sum_car_force_fl, 1, brem_start, 1);
        Math::axpy<T>(Constants::DIM, 1, cf_sum_car_force_fr, 1, brem_start, 1);
        Math::axpy<T>(Constants::DIM, 1, cf_sum_car_force_rl, 1, brem_start, 1);
        Math::axpy<T>(Constants::DIM, 1, cf_sum_car_force_rr, 1, brem_start, 1);

        //  upper_force1 - lower_force1 + upper_dampf1 - lower_dampf1 + local_FW1 + upper_rot_force1
        //  - car_rot_force1 - lower_rot_force1;
        brem_start += Constants::DIM;
        Math::copy<T>(Constants::DIM, cf_upper_force_fl, 1, brem_start, 1);
        Math::axpy<T>(Constants::DIM, -1, cf_lower_force_fl, 1, brem_start, 1);
        Math::axpy<T>(Constants::DIM, 1, cf_upper_dampf_fl, 1, brem_start, 1);
        Math::axpy<T>(Constants::DIM, -1, cf_lower_dampf_fl, 1, brem_start, 1);
        Math::axpy<T>(Constants::DIM, 1, FW_fl, 1, brem_start, 1);
        Math::axpy<T>(Constants::DIM, 1, cf_upper_rot_force_fl, 1, brem_start, 1);
        Math::axpy<T>(Constants::DIM, -1, cf_car_rot_force_fl, 1, brem_start, 1);
        Math::axpy<T>(Constants::DIM, -1, cf_lower_rot_force_fl, 1, brem_start, 1);

        // upper_force2 - lower_force2 + upper_dampf2 - lower_dampf2 + local_FW2 + upper_rot_force2
        // - car_rot_force2 - lower_rot_force2;
        brem_start += Constants::DIM;
        Math::copy<T>(Constants::DIM, cf_upper_force_fr, 1, brem_start, 1);
        Math::axpy<T>(Constants::DIM, -1, cf_lower_force_fr, 1, brem_start, 1);
        Math::axpy<T>(Constants::DIM, 1, cf_upper_dampf_fr, 1, brem_start, 1);
        Math::axpy<T>(Constants::DIM, -1, cf_lower_dampf_fr, 1, brem_start, 1);
        Math::axpy<T>(Constants::DIM, 1, FW_fr, 1, brem_start, 1);
        Math::axpy<T>(Constants::DIM, 1, cf_upper_rot_force_fr, 1, brem_start, 1);
        Math::axpy<T>(Constants::DIM, -1, cf_car_rot_force_fr, 1, brem_start, 1);
        Math::axpy<T>(Constants::DIM, -1, cf_lower_rot_force_fr, 1, brem_start, 1);

        // upper_force3 - lower_force3 + upper_dampf3 - lower_dampf3 + local_FW3 + upper_rot_force3
        // - car_rot_force3 - lower_rot_force3;
        brem_start += Constants::DIM;
        Math::copy<T>(Constants::DIM, cf_upper_force_rl, 1, brem_start, 1);
        Math::axpy<T>(Constants::DIM, -1, cf_lower_force_rl, 1, brem_start, 1);
        Math::axpy<T>(Constants::DIM, 1, cf_upper_dampf_rl, 1, brem_start, 1);
        Math::axpy<T>(Constants::DIM, -1, cf_lower_dampf_rl, 1, brem_start, 1);
        Math::axpy<T>(Constants::DIM, 1, FW_rl, 1, brem_start, 1);
        Math::axpy<T>(Constants::DIM, 1, cf_upper_rot_force_rl, 1, brem_start, 1);
        Math::axpy<T>(Constants::DIM, -1, cf_car_rot_force_rl, 1, brem_start, 1);
        Math::axpy<T>(Constants::DIM, -1, cf_lower_rot_force_rl, 1, brem_start, 1);

        //  upper_force4 - lower_force4 + upper_dampf4 - lower_dampf4 + local_FW4 + upper_rot_force4
        //  - car_rot_force4 - lower_rot_force4;
        brem_start += Constants::DIM;
        Math::copy<T>(Constants::DIM, cf_upper_force_rr, 1, brem_start, 1);
        Math::axpy<T>(Constants::DIM, -1, cf_lower_force_rr, 1, brem_start, 1);
        Math::axpy<T>(Constants::DIM, 1, cf_upper_dampf_rr, 1, brem_start, 1);
        Math::axpy<T>(Constants::DIM, -1, cf_lower_dampf_rr, 1, brem_start, 1);
        Math::axpy<T>(Constants::DIM, 1, FW_rr, 1, brem_start, 1);
        Math::axpy<T>(Constants::DIM, 1, cf_upper_rot_force_rr, 1, brem_start, 1);
        Math::axpy<T>(Constants::DIM, -1, cf_car_rot_force_rr, 1, brem_start, 1);
        Math::axpy<T>(Constants::DIM, -1, cf_lower_rot_force_rr, 1, brem_start, 1);

        // local_FR1; ...					  %vt1_dot
        brem_start += Constants::DIM;
        Math::copy<T>(Constants::DIM, cf_local_FR_fl, 1, brem_start, 1);

        // local_FR2; ...			 	      %vt2_dot
        brem_start += Constants::DIM;
        Math::copy<T>(Constants::DIM, cf_local_FR_fr, 1, brem_start, 1);

        // local_FR3; ...			          %vt3_dot
        brem_start += Constants::DIM;
        Math::copy<T>(Constants::DIM, cf_local_FR_rl, 1, brem_start, 1);

        // local_FR4];			              %vt4_dot
        brem_start += Constants::DIM;
        Math::copy<T>(Constants::DIM, cf_local_FR_rr, 1, brem_start, 1);
    }

    void solve_accelerations() {
        Math::potrs<T>(LAPACK_ROW_MAJOR, 'L', Constants::DIM, 1, A_Ic, Constants::DIM,
                       cf_sum_torque_spring_car, 1);

        Math::vector_elem_wise_product<T>(A_rem, cf_b_rem, accelerations, 9 * Constants::DIM);
    }

    void compute_quaternion_change_rate() {
        /*
         * get the derivative of the altitude (expressed in quaternions) from the angular
         * velocities
         *
         * Qc = 0.5 * [qc(4) -qc(3) qc(2);
         *             qc(3) qc(4) -qc(1);
         *             -qc(2) qc(1) qc(4);
         *             -qc(1) -qc(2) -qc(3)];
         * qc_dot = Qc * wc;
         */

        cf_Qc[0] = 0.5 * qc_[3];
        cf_Qc[1] = -0.5 * qc_[2];
        cf_Qc[2] = 0.5 * qc_[1];
        cf_Qc[3] = 0.5 * qc_[2];
        cf_Qc[4] = 0.5 * qc_[3];
        cf_Qc[5] = -0.5 * qc_[0];
        cf_Qc[6] = -0.5 * qc_[1];
        cf_Qc[7] = 0.5 * qc_[0];
        cf_Qc[8] = 0.5 * qc_[3];
        cf_Qc[9] = -0.5 * qc_[0];
        cf_Qc[10] = -0.5 * qc_[1];
        cf_Qc[11] = -0.5 * qc_[2];
        Math::gemv<T>(CblasRowMajor, CblasNoTrans, Constants::NUM_LEGS, Constants::DIM, 1, cf_Qc,
                      Constants::DIM, wc_, 1, 0, cf_qc_dot, 1);
    }

    void construct_f_vector(T* f_) {
        Math::copy<T>(Constants::DIM, cf_sum_torque_spring_car, 1, f_, 1);
        T* start_next = f_ + Constants::DIM;

        Math::copy<T>(9 * Constants::DIM, accelerations, 1, start_next, 1);
        start_next += 9 * Constants::DIM;

        Math::copy<T>(Constants::NUM_LEGS, cf_qc_dot, 1, start_next, 1);
        start_next += Constants::NUM_LEGS;
        // add vc to f vector
        Math::copy<T>(Constants::DIM, vc_, 1, start_next, 1);
        start_next += Constants::DIM;

        // add vw1 to f vector

        Math::copy<T>(Constants::DIM, vw_fl_, 1, start_next, 1);
        start_next += Constants::DIM;

        // add vw2 to f vector
        Math::copy<T>(Constants::DIM, vw_fr_, 1, start_next, 1);
        start_next += Constants::DIM;

        // add vw3 to f vector
        Math::copy<T>(Constants::DIM, vw_rl_, 1, start_next, 1);
        start_next += Constants::DIM;

        // add vw_rr to f vector
        Math::copy<T>(Constants::DIM, vw_rr_, 1, start_next, 1);
        start_next += Constants::DIM;

        // add vt1 to f vector
        Math::copy<T>(Constants::DIM, vt_fl_, 1, start_next, 1);
        start_next += Constants::DIM;

        // add vt2 to f vector
        Math::copy<T>(Constants::DIM, vt_fr_, 1, start_next, 1);
        start_next += Constants::DIM;

        // add vt3 to f vector
        Math::copy<T>(Constants::DIM, vt_rl_, 1, start_next, 1);
        start_next += Constants::DIM;

        // add vt_rr to f vector
        Math::copy<T>(Constants::DIM, vt_rr_, 1, start_next, 1);
        start_next += Constants::DIM;
    }

public:
    /**
     * Constructor
     */
    MBDMethod() {
        MemoryAllocation();

        ReadFromXML();

        getCholeskyDecomposition();
    }

    /**
     * Initializes the time iteration and handles the numerical scheme
     */
    void solve(T* solution_vector) {
        // From the formulation we have 61 dimensions in the solution vector
        size_t solution_size = (this->num_iter + 1) * this->solution_dim;
        T* complete_vector = Math::calloc<T>(solution_size);
        x_vector = Math::calloc<T>(solution_dim);
        Math::get_tilda<T>(r_fl, r_fl_tilda);
        Math::get_tilda<T>(r_fr, r_fr_tilda);
        Math::get_tilda<T>(r_rl, r_rl_tilda);
        Math::get_tilda<T>(r_rr, r_rr_tilda);

        /*
         * Preparing x_vector in the form of
         * x_vector = [wc; ...     % 3 1:3
         *             vc; ...     % 3 4:6
         *             vw1; ...    % 3 7:9
         *             vw2; ...    % 3 10:12
         *             vw3; ...    % 3 13:15
         *             vw4; ...    % 3 16:18
         *             vt1; ...    % 3 19:21
         *             vt2; ...    % 3 22:24
         *             vt3; ...    % 3 25:27
         *             vt4; ...    % 3 28:30
         *             qc; ...     % 4 31:34
         *             pcc; ...    % 3 35:37
         *             pw1; ...    % 3 38:40
         *             pw2; ...    % 3 41:43
         *             pw3; ...    % 3 44:46
         *             pw4; ...    % 3 47:49
         *             pt1; ...    % 3 50:52
         *             pt2; ...    % 3 53:55
         *             pt3; ...    % 3 56:58
         *             pt4];       % 3 59:61
         */
        size_t i, j;
        i = 10;
        j = 0;

        // qc
        qc_ = x_vector + i * (Constants::DIM) + j * (Constants::NUM_LEGS);
        Math::copy<T>(Constants::NUM_LEGS, initial_orientation, 1,
                      x_vector + i * (Constants::DIM) + j * (Constants::NUM_LEGS), 1);
        j++;
        // pcc
        pcc_ = x_vector + i * (Constants::DIM) + j * (Constants::NUM_LEGS);
        Math::copy<T>(Constants::DIM, pcc, 1, pcc_, 1);
        i++;
        std::cout << "zoodle" << std::endl;
        // pw1
        pw_fl_ = x_vector + i * (Constants::DIM) + j * (Constants::NUM_LEGS);
        i++;
        // pw2
        pw_fr_ = x_vector + i * (Constants::DIM) + j * (Constants::NUM_LEGS);
        i++;
        // pw3
        pw_rl_ = x_vector + i * (Constants::DIM) + j * (Constants::NUM_LEGS);
        i++;
        // pw_rr
        pw_rr_ = x_vector + i * (Constants::DIM) + j * (Constants::NUM_LEGS);
        i++;
        // pt1
        pt_fl_ = x_vector + i * (Constants::DIM) + j * (Constants::NUM_LEGS);
        i++;
        // pt2
        pt_fr_ = x_vector + i * (Constants::DIM) + j * (Constants::NUM_LEGS);
        i++;
        // pt3
        pt_rl_ = x_vector + i * (Constants::DIM) + j * (Constants::NUM_LEGS);
        i++;
        // pt_rr
        pt_rr_ = x_vector + i * (Constants::DIM) + j * (Constants::NUM_LEGS);

        get_initial_length(qc_, r_fl, r_fr, r_rl, r_rr, pcc, initial_upper_spring_length,
                           initial_lower_spring_length, pw_fl_, pw_fr_, pw_rl_, pw_rr_, pt_fl_,
                           pt_fr_, pt_rl_, pt_rr_);

        // overwrites the initial velocity values
        if (boundary_conditions == BoundaryConditionRoad::CIRCULAR)
            circular_path_initialization(vc, vw_fl, vw_fr, vw_rl, vw_rr, vt_fl, vt_fr, vt_rl, vt_rr,
                                         initial_angular_velocity, pcc_, pt_fl_, pt_fr_, pt_rl_,
                                         pt_rr_, radius_circular_path);
        i = 0;
        j = 0;
        // wc
        wc_ = x_vector + i * (Constants::DIM) + j * (Constants::NUM_LEGS);
        Math::copy<T>(Constants::DIM, initial_angular_velocity, 1, wc_, 1);
        i++;
        // vc
        vc_ = x_vector + i * (Constants::DIM) + j * (Constants::NUM_LEGS);
        Math::copy<T>(Constants::DIM, vc, 1, vc_, 1);
        i++;
        // vw1
        vw_fl_ = x_vector + i * (Constants::DIM) + j * (Constants::NUM_LEGS);
        Math::copy<T>(Constants::DIM, vw_fl, 1, vw_fl_, 1);
        i++;
        // vw2
        vw_fr_ = x_vector + i * (Constants::DIM) + j * (Constants::NUM_LEGS);
        Math::copy<T>(Constants::DIM, vw_fr, 1, vw_fr_, 1);
        i++;
        // vw_rl
        vw_rl_ = x_vector + i * (Constants::DIM) + j * (Constants::NUM_LEGS);
        Math::copy<T>(Constants::DIM, vw_rl, 1, vw_rl_, 1);
        i++;
        // vw4
        vw_rr_ = x_vector + i * (Constants::DIM) + j * (Constants::NUM_LEGS);
        Math::copy<T>(Constants::DIM, vw_rr, 1, vw_rr_, 1);
        i++;
        // vt1
        vt_fl_ = x_vector + i * (Constants::DIM) + j * (Constants::NUM_LEGS);
        Math::copy<T>(Constants::DIM, vt_fl, 1, vt_fl_, 1);
        i++;
        // vt2
        vt_fr_ = x_vector + i * (Constants::DIM) + j * (Constants::NUM_LEGS);
        Math::copy<T>(Constants::DIM, vt_fr, 1, vt_fr_, 1);
        i++;
        // vt_rl
        vt_rl_ = x_vector + i * (Constants::DIM) + j * (Constants::NUM_LEGS);
        Math::copy<T>(Constants::DIM, vt_rl, 1, vt_rl_, 1);
        i++;
        // vt4
        vt_rr_ = x_vector + i * (Constants::DIM) + j * (Constants::NUM_LEGS);
        Math::copy<T>(Constants::DIM, vt_rr, 1, vt_rr_, 1);
        i++;

        j = 0;
        i = 0;
        while (i < Constants::DIM) {
            A_rem[j] = 1. / this->mass;
            i++;
            j++;
        }
        i = 0;
        while (i < Constants::DIM) {
            A_rem[j] = 1. / mass_wheel[0];
            i++;
            j++;
        }
        i = 0;
        while (i < Constants::DIM) {
            A_rem[j] = 1. / mass_wheel[1];
            i++;
            j++;
        }
        i = 0;
        while (i < Constants::DIM) {
            A_rem[j] = 1. / mass_wheel[2];
            i++;
            j++;
        }
        i = 0;
        while (i < Constants::DIM) {
            A_rem[j] = 1. / mass_wheel[3];
            i++;
            j++;
        }
        i = 0;
        while (i < Constants::DIM) {
            A_rem[j] = 1. / mass_tyre[0];
            i++;
            j++;
        }
        i = 0;
        while (i < Constants::DIM) {
            A_rem[j] = 1. / mass_tyre[1];
            i++;
            j++;
        }
        i = 0;
        while (i < Constants::DIM) {
            A_rem[j] = 1. / mass_tyre[2];
            i++;
            j++;
        }
        i = 0;
        while (i < Constants::DIM) {
            A_rem[j] = 1. / mass_tyre[3];
            i++;
            j++;
        }

        compute_f_mem_alloc();

        if (used_solver == MBDSolver::BROYDEN_CN) {
            Math::Solvers<T, MBDMethod>::Broyden_CN(this, x_vector, complete_vector, this->h,
                                                    this->num_iter, this->tol, this->max_iter);
        }
        else if (used_solver == MBDSolver::RUNGE_KUTTA_4) {
            Math::Solvers<T, MBDMethod>::RK4(this, x_vector, complete_vector, this->h,
                                             this->num_iter, this->tol, this->max_iter);
        }
        else if (used_solver == MBDSolver::BROYDEN_BDF2) {
            Math::Solvers<T, MBDMethod>::Broyden_PDF2(this, x_vector, complete_vector, this->h,
                                                      this->num_iter, this->tol, this->max_iter);
        }
        else if (used_solver == MBDSolver::BROYDEN_EULER) {
            Math::Solvers<T, MBDMethod>::Broyden_Euler(this, x_vector, complete_vector, this->h,
                                                       this->num_iter, this->tol, this->max_iter);
        }
        else if (used_solver == MBDSolver::EXPLICIT_EULER) {
            std::cout << "Explicit solver hasn't been implemented, you don't want to use it"
                      << std::endl;
        }
        else {
            std::cout << "sorry man, the solver you picked for MBD is weird and hasn't been "
                         "implemented yet"
                      << std::endl;
        }

        compute_f_clean();

        T* start = complete_vector + (this->num_iter) * this->solution_dim;
        Math::copy<T>(this->solution_dim, start, 1, solution_vector, 1);
        //	std::cout << "Solution copied!\n" << std::endl;
        Math::free<T>(complete_vector);
        Math::free<T>(x_vector);
    }

    /**
     * Solver which is called at each time step
     * Computes the forces and torques on each point mass and computes the right hand side of the
     * ODE
     * \param[in] x_ current solution of the system
     * \param[in] t_ current simulation time
     * \param[out] f_ the load vector
     */
    void compute_f3D_reduced(T* x_, T t_, T* f_) {
        /*
         * Small performance gain might be possible by transforming C_cN to column major
         * Note: corresponding MKL function call have to be changed too
         */
        const MKL_INT mkl_DIM = Constants::DIM;
        const MKL_INT mkl_incx = 1;
        const MKL_INT mkl_incy = 1;

        get_current_variables(x_);

        get_car_orientation();

        compute_spring_lengths(pcc_, pw_fl_, pt_fl_, cf_r_up_fl, cf_r_low_fl, this->r_fl,
                               norm_r_up_fl, inv_norm_r_up_fl, norm_r_low_fl, inv_norm_r_low_fl);
        compute_spring_lengths(pcc_, pw_fr_, pt_fr_, cf_r_up_fr, cf_r_low_fr, this->r_fr,
                               norm_r_up_fr, inv_norm_r_up_fr, norm_r_low_fr, inv_norm_r_low_fr);
        compute_spring_lengths(pcc_, pw_rl_, pt_rl_, cf_r_up_rl, cf_r_low_rl, this->r_rl,
                               norm_r_up_rl, inv_norm_r_up_rl, norm_r_low_rl, inv_norm_r_low_rl);
        compute_spring_lengths(pcc_, pw_rr_, pt_rr_, cf_r_up_rr, cf_r_low_rr, this->r_rr,
                               norm_r_up_rr, inv_norm_r_up_rr, norm_r_low_rr, inv_norm_r_low_rr);

        apply_stiffness_interpolation();

        compute_angles();

        compute_elongational_forces();

        compute_damping_forces();

        compute_torques();

        compute_torques_resultant_forces();

        compute_car_body_forces();

        compute_external_forces();

        if (boundary_conditions == BoundaryConditionRoad::FIXED) {
            get_fixed_road_force(cf_local_FR_fl, cf_local_FR_fr, cf_local_FR_rl, cf_local_FR_rr);
        }
        else if (boundary_conditions == BoundaryConditionRoad::NONFIXED) {
            get_nonfixed_road_force(cf_local_FR_fl, cf_local_FR_fr, cf_local_FR_rl, cf_local_FR_rr);
        }
        else if (boundary_conditions == BoundaryConditionRoad::CIRCULAR) {
            get_circular_road_force(cf_local_FR_fl, vt_fl_, mass_tyre_rr, pt_fl_);
            get_circular_road_force(cf_local_FR_fr, vt_fr_, mass_tyre_rl, pt_fr_);
            get_circular_road_force(cf_local_FR_rl, vt_rl_, mass_tyre_fl, pt_rl_);
            get_circular_road_force(cf_local_FR_rr, vt_rr_, mass_tyre_fr, pt_rr_);
        }

        compute_car_body_total_torque();

        construct_right_hand_side();

        solve_accelerations();

        compute_quaternion_change_rate();

        construct_f_vector(f_);
    }

    /**
     * Returns the dimension of the solution vector (=61)
     */
    size_t get_solution_dimension() { return this->solution_dim; }

    /**
     * Beatiful output of the result
     * \param sln solution vector
     */
    void print_final_result(T* sln) {
        std::cout << std::scientific;
        std::cout << std::setprecision(15);

        std::cout << "MBD: angular velocity w=\n\t[" << sln[0] << "\n\t " << sln[1] << "\n\t "
                  << sln[2] << "]" << std::endl;
        std::cout << "MBD: car body velocity vc=\n\t[" << sln[3] << "\n\t " << sln[4] << "\n\t "
                  << sln[5] << "]" << std::endl;
        std::cout << "MBD: front-left wheel velocity vw_fl=\n\t[" << sln[6] << "\n\t " << sln[7]
                  << "\n\t " << sln[8] << "]" << std::endl;
        std::cout << "MBD: front-right wheel velocity vw_fr=\n\t[" << sln[9] << "\n\t " << sln[10]
                  << "\n\t " << sln[11] << "]" << std::endl;
        std::cout << "MBD: rear-left wheel velocity vw_rl=\n\t[" << sln[12] << "\n\t " << sln[13]
                  << "\n\t " << sln[14] << "]" << std::endl;
        std::cout << "MBD: rear-right wheel velocity vw_rr=\n\t[" << sln[15] << "\n\t " << sln[16]
                  << "\n\t " << sln[17] << "]" << std::endl;
        std::cout << "MBD: front-left tyre velocity vt_fl=\n\t[" << sln[18] << "\n\t " << sln[19]
                  << "\n\t " << sln[20] << "]" << std::endl;
        std::cout << "MBD: front-right tyre velocity vt_fr=\n\t[" << sln[21] << "\n\t " << sln[22]
                  << "\n\t " << sln[23] << "]" << std::endl;
        std::cout << "MBD: rear-left tyre velocity vt_rl=\n\t[" << sln[24] << "\n\t " << sln[25]
                  << "\n\t " << sln[26] << "]" << std::endl;
        std::cout << "MBD: rear-right tyre velocity vt_rr=\n\t[" << sln[27] << "\n\t " << sln[28]
                  << "\n\t " << sln[29] << "]" << std::endl;
        std::cout << "MBD: orientation q=\n\t[" << sln[30] << "\n\t " << sln[31] << "\n\t "
                  << sln[32] << "\n\t " << sln[33] << "]" << std::endl;
        std::cout << "MBD: car body position pc=\n\t[" << sln[34] << "\n\t " << sln[35] << "\n\t "
                  << sln[36] << "]" << std::endl;
        std::cout << "MBD: front-left wheel position pw_fl=\n\t[" << sln[37] << "\n\t " << sln[38]
                  << "\n\t " << sln[39] << "]" << std::endl;
        std::cout << "MBD: front-right wheel position pw_fr=\n\t[" << sln[40] << "\n\t " << sln[41]
                  << "\n\t " << sln[42] << "]" << std::endl;
        std::cout << "MBD: rear-left wheel position pw_rl=\n\t[" << sln[43] << "\n\t " << sln[44]
                  << "\n\t " << sln[45] << "]" << std::endl;
        std::cout << "MBD: rear-right wheel position pw_rr=\n\t[" << sln[46] << "\n\t " << sln[47]
                  << "\n\t " << sln[48] << "]" << std::endl;
        std::cout << "MBD: front-left tyre position pt_fl=\n\t[" << sln[49] << "\n\t " << sln[50]
                  << "\n\t " << sln[51] << "]" << std::endl;
        std::cout << "MBD: front-right tyre position pt_fr=\n\t[" << sln[52] << "\n\t " << sln[53]
                  << "\n\t " << sln[54] << "]" << std::endl;
        std::cout << "MBD: rear-left tyre position pt_rl=\n\t[" << sln[55] << "\n\t " << sln[56]
                  << "\n\t " << sln[57] << "]" << std::endl;
        std::cout << "MBD: rear-right tyre position pt_rr=\n\t[" << sln[58] << "\n\t " << sln[59]
                  << "\n\t " << sln[60] << "]" << std::endl;
    }

    /**
     * Destructor
     */
    ~MBDMethod() {
        // TODO: consider removing after checking performance.
        mkl_free_buffers();

        Math::free<T>(r_fl);
        Math::free<T>(r_fr);
        Math::free<T>(r_rl);
        Math::free<T>(r_rr);
        Math::free<T>(Ic);
        Math::free<T>(mass_wheel);
        Math::free<T>(mass_tyre);
        Math::free<T>(upper_spring_length);
        Math::free<T>(lower_spring_length);
        Math::free<T>(upper_spring_stiffness);
        Math::free<T>(lower_spring_stiffness);
        Math::free<T>(upper_rotational_stiffness);
        Math::free<T>(lower_rotational_stiffness);
        Math::free<T>(initial_upper_spring_length);
        Math::free<T>(initial_lower_spring_length);
        Math::free<T>(initial_orientation);
        Math::free<T>(initial_angular_velocity);
        Math::free<T>(upper_spring_damping);
        Math::free<T>(lower_spring_damping);
        Math::free<T>(vc);
        Math::free<T>(vw_fl);
        Math::free<T>(vw_fr);
        Math::free<T>(vw_rl);
        Math::free<T>(vw_rr);
        Math::free<T>(vt_fl);
        Math::free<T>(vt_fr);
        Math::free<T>(vt_rl);
        Math::free<T>(vt_rr);
        Math::free<T>(pcc);
        Math::free<T>(FC);
        Math::free<T>(r_fl_tilda);
        Math::free<T>(r_fr_tilda);
        Math::free<T>(r_rl_tilda);
        Math::free<T>(r_rr_tilda);
        Math::free<T>(FW_fl);
        Math::free<T>(FW_fr);
        Math::free<T>(FW_rl);
        Math::free<T>(FW_rr);
        Math::free<T>(FT_fl);
        Math::free<T>(FT_fr);
        Math::free<T>(FT_rl);
        Math::free<T>(FT_rr);
        Math::free<T>(A_Ic);
        Math::free<T>(A_rem);
        Math::free<T>(center_of_circle);
    }
};

}  // namespace EVAA
