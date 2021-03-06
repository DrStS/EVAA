// TODO: Copyright header

#pragma once

#include <chrono>
#include <string>

#include "Constants.h"
#include "IO/Output.h"
#include "MathLibrary.h"
#include "MetaDatabase.h"

#ifdef USE_HDF5
#include "IO/OutputHDF5.h"
#endif

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
    std::string solver_name;
    MBDSolver used_solver;

    // Environment conditions
    LagrangianBoundaryConditionRoad lagrangian_boundary_conditions;
    EulerianBoundaryConditionRoad eulerian_boundary_conditions;

    // Car Definition
    T* Ic;
    T *upper_spring_length, *lower_spring_length;
    T *upper_spring_stiffness, *lower_spring_stiffness, *upper_spring_damping, *lower_spring_damping;
    T *upper_rotational_stiffness, *lower_rotational_stiffness;
    T *vc, *vw_fl, *vw_fr, *vw_rl, *vw_rr, *vt_fl, *vt_fr, *vt_rl,
        *vt_rr;  // velocity of the center of mass, wheel and tyre.

    // External Forces terms

    T* FC;
    T *FT_fl, *FT_fr, *FT_rl, *FT_rr;
    T *FW_fl, *FW_fr, *FW_rl, *FW_rr;

    // Initial condition params

    T* initial_angular_velocity;

    // Arrays for intermediate steps

    T *r_fl, *r_fr, *r_rl, *r_rr;
    T* pcc;  // position of center of mass
    T* x_vector;

    // Auxillary Parameters for Compute_f function

    T *r_fl_tilda, *r_fr_tilda, *r_rl_tilda, *r_rr_tilda, *A_Ic, *A_rem;

    // Variables needed in compute_f function

    T *cf_C_cN, *cf_r_up_fl, *cf_r_up_fr, *cf_r_up_rl, *cf_r_up_rr, *cf_r_low_fl, *cf_r_low_fr, *cf_r_low_rl, *cf_r_low_rr;
    T *cf_upper_normal_fl, *cf_upper_normal_fr, *cf_upper_normal_rl, *cf_upper_normal_rr, *cf_lower_normal_fl, *cf_lower_normal_fr, *cf_lower_normal_rl, *cf_lower_normal_rr, *cf_col_dat;
    T *cf_upper_force_fl, *cf_upper_force_fr, *cf_upper_force_rl, *cf_upper_force_rr, *cf_lower_force_fl, *cf_lower_force_fr, *cf_lower_force_rl, *cf_lower_force_rr;
    T *cf_upper_dampf_fl, *cf_upper_dampf_fr, *cf_upper_dampf_rl, *cf_upper_dampf_rr, *cf_lower_dampf_fl, *cf_lower_dampf_fr, *cf_lower_dampf_rl, *cf_lower_dampf_rr, *cf_temp;
    T *cf_upper_angle_fl, *cf_upper_angle_fr, *cf_upper_angle_rl, *cf_upper_angle_rr, *cf_lower_angle_fl, *cf_lower_angle_fr, *cf_lower_angle_rl, *cf_lower_angle_rr;
    T *cf_upper_S_fl, *cf_upper_S_fr, *cf_upper_S_rl, *cf_upper_S_rr, *cf_lower_S_fl, *cf_lower_S_fr, *cf_lower_S_rl, *cf_lower_S_rr;
    T *cf_lower_rot_force_fl, *cf_lower_rot_force_fr, *cf_lower_rot_force_rl, *cf_lower_rot_force_rr, *cf_upper_rot_force_fl, *cf_upper_rot_force_fr, *cf_upper_rot_force_rl, *cf_upper_rot_force_rr;
    T *cf_car_rot_force_fl, *cf_car_rot_force_fr, *cf_car_rot_force_rl, *cf_car_rot_force_rr, *cf_sum_car_force_fl, *cf_sum_car_force_fr, *cf_sum_car_force_rl, *cf_sum_car_force_rr;
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

    // solution vector from the previous iteration (for the flying car)
    T* _previousSolution;
    T* _previousTyreForces;
    bool _iterationBegin;
    size_t _currentIteration;

#ifdef USE_HDF5
#ifdef USE_CHECKPOINTS
    HDF5::OutputHDF5<T>* _checkpointsMBD;
    HDF5::OutputHDF5<T>* _checkpointsMBDFormatted;
    std::string _groupNameCheckpoints;  // basic name for a checkpoint group
#endif
#endif  //  USE_HDF5

    /**
     * Calculate the positions of the tyres and wheels according to the initial
     * orientation of the car The legs always form a 90 degrees angle to the car
     * body, such that the rotational springs are at rest
     */
    void get_initial_length(T* initial_orientation_, const T* r_fl_, const T* r_fr_, const T* r_rl_, const T* r_rr_, const T* pcc_, T* wheel_coordinate_fl_, T* wheel_coordinate_fr_, T* wheel_coordinate_rl_, T* wheel_coordinate_rr_, T* tyre_coordinate_fl_, T* tyre_coordinate_fr_, T* tyre_coordinate_rl_, T* tyre_coordinate_rr_)

    {
        /*
         * To reduce memory trace and better use cache this function is
         * implemented in following fashion: Original steps for computation of
         * one component:
         * 1. qc = qc/norm(qc);
         * 2. C_Nc = GetBasis(qc);
         * 3. global_z = C_Nc(:,2);
         * 4. global_z = -global_z / norm(global_z);
         * 5. global_r_fl = pcc + C_Nc*r_fl;
         * 6. upper_global_spring__fl = upper_length(_fl) * global_z;
         * 7. lower_global_spring__fl = lower_length(_fl)*global_z;
         * 8. pw_fl = global_r_fl + upper_global_spring__fl;
         * 9. pt_fl = pw_fl + lower_global_spring__fl;
         * Modified steps for computation of one component:
         * 1. qc = qc/norm(qc);
         * 2. C_Nc = GetBasis(qc);
         * 3. global_z = C_Nc(:,2);
         * 4. global_z = global_z / norm(global_z);
         * 5. pw_fl = pcc;
         * 6. pw_fl = pw_fl + C_Nc*r_fl;
         * 7. pw_fl = pw_fl + upper_length(1)*global_z;
         * 8. pt_fl = pw_fl
         * 9. pt_fl	= pt_fl + lower_length(_fl)*global_z;
         */
        auto& db = MetaDatabase<T>::getDatabase();
        const T* initial_upper_spring_length = db.getBodySpringInitialLengthVector();
        const T* initial_lower_spring_length = db.getTyreSpringInitialLengthVector();

        T* global_z = Math::calloc<T>(Constants::DIM);
        T* C_Nc = Math::calloc<T>(Constants::DIMDIM);

        // 1. qc = qc/norm(qc); This is in quaternions
        T nrm = Math::nrm2<T>(Constants::NUM_LEGS, initial_orientation_, 1);
        Math::scal<T>(Constants::NUM_LEGS, 1. / nrm, initial_orientation_, 1);

        // 2. C_Nc = GetBasis(qc);
        Math::GetBasis<T>(initial_orientation_, C_Nc);

        // 3. global_z = C_Nc(:,2);
        Math::copy<T>(Constants::DIM, C_Nc + 2, Constants::DIM, global_z, 1);

        // 4. global_z = -global_z / norm(global_z);
        nrm = Math::nrm2<T>(Constants::DIM, global_z, 1);
        Math::scal<T>(Constants::DIM, -1. / nrm, global_z, 1);

        // Leg 1

        // 5.	pw1 = pcc;
        Math::copy<T>(Constants::DIM, pcc_, 1, wheel_coordinate_fl_, 1);

        // 6.	pw1 = pw1 + C_Nc*r1;
        Math::gemv<T>(CblasRowMajor, CblasNoTrans, Constants::DIM, Constants::DIM, 1, C_Nc, Constants::DIM, r_fl_, 1, 1, wheel_coordinate_fl_, 1);

        // 7.	pw1 = pw1 + upper_length(1)*global_z;
        Math::axpy<T>(Constants::DIM, initial_upper_spring_length[0], global_z, 1, wheel_coordinate_fl_, 1);

        // 8.	pt1 = pw1
        Math::copy<T>(Constants::DIM, wheel_coordinate_fl_, 1, tyre_coordinate_fl_, 1);

        // 9.	pt1 = pw1 + lower_length(1)*global_z;
        Math::axpy<T>(Constants::DIM, initial_lower_spring_length[0], global_z, 1, tyre_coordinate_fl_, 1);

        // Leg 2

        // 5.	pw2 = pcc;
        Math::copy<T>(Constants::DIM, pcc_, 1, wheel_coordinate_fr_, 1);

        // 6.	pw2 = pw2 + C_Nc*r2;
        Math::gemv<T>(CblasRowMajor, CblasNoTrans, Constants::DIM, Constants::DIM, 1, C_Nc, Constants::DIM, r_fr_, 1, 1, wheel_coordinate_fr_, 1);

        // 7.	pw2 = pw2 + upper_length(2)*global_z;
        Math::axpy<T>(Constants::DIM, initial_upper_spring_length[1], global_z, 1, wheel_coordinate_fr_, 1);

        // 8.	pt2 = pw2
        Math::copy<T>(Constants::DIM, wheel_coordinate_fr_, 1, tyre_coordinate_fr_, 1);

        // 9.	pt2 = pw2 + lower_length(2)*global_z;
        Math::axpy<T>(Constants::DIM, initial_lower_spring_length[1], global_z, 1, tyre_coordinate_fr_, 1);

        // Leg 3

        // 5.	pw3 = pcc;
        Math::copy<T>(Constants::DIM, pcc_, 1, wheel_coordinate_rl_, 1);

        // 6.	pw3 = pw3 + C_Nc*r3;
        Math::gemv<T>(CblasRowMajor, CblasNoTrans, Constants::DIM, Constants::DIM, 1, C_Nc, Constants::DIM, r_rl_, 1, 1, wheel_coordinate_rl_, 1);

        // 7.	pw3 = pw3 + upper_length(3)*global_z;
        Math::axpy<T>(Constants::DIM, initial_upper_spring_length[2], global_z, 1, wheel_coordinate_rl_, 1);

        // 8.	pt3 = pw3
        Math::copy<T>(Constants::DIM, wheel_coordinate_rl_, 1, tyre_coordinate_rl_, 1);

        // 9.	pt3 = pw3 + lower_length(3)*global_z;
        Math::axpy<T>(Constants::DIM, initial_lower_spring_length[2], global_z, 1, tyre_coordinate_rl_, 1);

        // Leg 4

        // 5.	pw4 = pcc;
        Math::copy<T>(Constants::DIM, pcc_, 1, wheel_coordinate_rr_, 1);

        // 6.	pw4 = pw4 + C_Nc*r4;
        Math::gemv<T>(CblasRowMajor, CblasNoTrans, Constants::DIM, Constants::DIM, 1, C_Nc, Constants::DIM, r_rr_, 1, 1, wheel_coordinate_rr_, 1);

        // 7.	pw4 = pw4 + upper_length(4)*global_z;
        Math::axpy<T>(Constants::DIM, initial_upper_spring_length[3], global_z, 1, wheel_coordinate_rr_, 1);

        // 8.	pt4 = pw4
        Math::copy<T>(Constants::DIM, wheel_coordinate_rr_, 1, tyre_coordinate_rr_, 1);

        // 9.	pt4 = pw4 + lower_length(4)*global_z;
        Math::axpy<T>(Constants::DIM, initial_lower_spring_length[3], global_z, 1, tyre_coordinate_rr_, 1);

        Math::free<T>(global_z);
        Math::free<T>(C_Nc);
    }

    void MemoryAllocation() {
        // Memory Allocation and matrix formulation

        r_fl = Math::calloc<T>(Constants::DIM);
        r_fr = Math::calloc<T>(Constants::DIM);
        r_rl = Math::calloc<T>(Constants::DIM);
        r_rr = Math::calloc<T>(Constants::DIM);
        Ic = Math::calloc<T>(Constants::DIMDIM);
        upper_spring_length = Math::calloc<T>(Constants::NUM_LEGS);
        lower_spring_length = Math::calloc<T>(Constants::NUM_LEGS);
        upper_spring_stiffness = Math::calloc<T>(Constants::NUM_LEGS);
        lower_spring_stiffness = Math::calloc<T>(Constants::NUM_LEGS);
        upper_rotational_stiffness = Math::calloc<T>(Constants::NUM_LEGS);
        lower_rotational_stiffness = Math::calloc<T>(Constants::NUM_LEGS);
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
        r_fl_tilda = Math::calloc<T>(Constants::DIMDIM);
        r_fr_tilda = Math::calloc<T>(Constants::DIMDIM);
        r_rl_tilda = Math::calloc<T>(Constants::DIMDIM);
        r_rr_tilda = Math::calloc<T>(Constants::DIMDIM);
        FT_fl = Math::calloc<T>(Constants::DIM);
        FT_rr = Math::calloc<T>(Constants::DIM);
        FT_rl = Math::calloc<T>(Constants::DIM);
        FT_fr = Math::calloc<T>(Constants::DIM);
        FW_fl = Math::calloc<T>(Constants::DIM);
        FW_rr = Math::calloc<T>(Constants::DIM);
        FW_rl = Math::calloc<T>(Constants::DIM);
        FW_fr = Math::calloc<T>(Constants::DIM);
        A_Ic = Math::calloc<T>(Constants::DIMDIM);
        A_rem = Math::calloc<T>(Constants::DIMDIM * Constants::DIM);
        accelerations = Math::calloc<T>(Constants::DIMDIM * Constants::DIM);
        _previousSolution = Math::calloc<T>(Constants::MBD_SOLUTION_SIZE);
        _previousTyreForces = Math::calloc<T>(Constants::NUM_LEGS);
    }

    void ReadFromXML() {
        // Simulation Parameters
        auto& db = MetaDatabase<T>::getDatabase();

        h = db.getTimeStepSize();
        num_iter = db.getNumberOfTimeIterations();
        max_iter = db.getMaxNumberOfBroydenIterationForMBD();
        tol = db.getToleranceBroydenIterationForMBD();
        used_solver = db.getUsedSolverForMBD();
        lagrangian_boundary_conditions = db.getLagrangianRoadConditions();
        eulerian_boundary_conditions = db.getEulerianRoadConditions();

        // Car Definition
        T gx = db.getGravityField()[0];
        T gy = db.getGravityField()[1];
        T gz = db.getGravityField()[2];

        // Fill up vectors
        Math::copy(Constants::NUM_LEGS, db.getBodyStiffnessVector(), 1, upper_spring_stiffness, 1);
        Math::copy(Constants::NUM_LEGS, db.getTyreStiffnessVector(), 1, lower_spring_stiffness, 1);
        Math::copy(Constants::NUM_LEGS, db.getBodyDampingVector(), 1, upper_spring_damping, 1);
        Math::copy(Constants::NUM_LEGS, db.getTyreDampingVector(), 1, lower_spring_damping, 1);
        Math::copy(Constants::NUM_LEGS, db.getBodySpringLengthVector(), 1, upper_spring_length, 1);
        Math::copy(Constants::NUM_LEGS, db.getTyreSpringLengthVector(), 1, lower_spring_length, 1);

        // TODO: Extract constant or use MetaDatabase and vectorize copying.
        for (auto i = 0; i < Constants::NUM_LEGS; i++) {
            upper_rotational_stiffness[i] = 1e5;
            lower_rotational_stiffness[i] = 1e5;
        }

        Math::copy(Constants::DIM, db.getWheelInitialVelocityFrontLeft(), 1, vw_fl, 1);
        Math::copy(Constants::DIM, db.getWheelInitialVelocityFrontRight(), 1, vw_fr, 1);
        Math::copy(Constants::DIM, db.getWheelInitialVelocityRearLeft(), 1, vw_rl, 1);
        Math::copy(Constants::DIM, db.getWheelInitialVelocityRearRight(), 1, vw_rr, 1);

        Math::copy(Constants::DIM, db.getTyreInitialVelocityFrontLeft(), 1, vt_fl, 1);
        Math::copy(Constants::DIM, db.getTyreInitialVelocityFrontRight(), 1, vt_fr, 1);
        Math::copy(Constants::DIM, db.getTyreInitialVelocityRearLeft(), 1, vt_rl, 1);
        Math::copy(Constants::DIM, db.getTyreInitialVelocityRearRight(), 1, vt_rr, 1);

        // diag(Ic) = diag(MoI)
        Math::copy(Constants::DIM, db.getMomentOfInertiaVector(), Constants::DIM + 1, Ic, Constants::DIM + 1);

        Math::copy(Constants::DIM, db.getBodyInitialVelocity(), 1, vc, 1);
        Math::copy(Constants::DIM, db.getBodyInitialPosition(), 1, pcc, 1);
        Math::copy(Constants::DIM, db.getBodyExternalForce(), 1, FC, 1);
        Math::copy(Constants::DIM, db.getTyreExternalForceFrontLeft(), 1, FT_fl, 1);
        Math::copy(Constants::DIM, db.getTyreExternalForceFrontRight(), 1, FT_fr, 1);
        Math::copy(Constants::DIM, db.getTyreExternalForceRearLeft(), 1, FT_rl, 1);
        Math::copy(Constants::DIM, db.getTyreExternalForceRearRight(), 1, FT_rr, 1);
        Math::copy(Constants::DIM, db.getWheelExternalForceFrontLeft(), 1, FW_fl, 1);
        Math::copy(Constants::DIM, db.getWheelExternalForceFrontRight(), 1, FW_fr, 1);
        Math::copy(Constants::DIM, db.getWheelExternalForceRearLeft(), 1, FW_rl, 1);
        Math::copy(Constants::DIM, db.getWheelExternalForceRearRight(), 1, FW_rr, 1);

        int i;

        i = 0;
        r_fl[i] = db.getLongitudalLegPositionFrontLeft();
        r_fr[i] = db.getLongitudalLegPositionFrontRight();
        r_rl[i] = -db.getLongitudalLegPositionRearLeft();
        r_rr[i] = -db.getLongitudalLegPositionRearRight();
        FC[i] += db.getBodyMass() * gx;
        FT_fl[i] += db.getTyreMassFrontLeft() * gx;
        FT_fr[i] += db.getTyreMassFrontRight() * gx;
        FT_rl[i] += db.getTyreMassRearLeft() * gx;
        FT_rr[i] += db.getTyreMassRearRight() * gx;
        FW_fl[i] += db.getWheelMassFrontLeft() * gx;
        FW_fr[i] += db.getWheelMassFrontRight() * gx;
        FW_rl[i] += db.getWheelMassRearLeft() * gx;
        FW_rr[i] += db.getWheelMassRearRight() * gx;

        i = 1;
        r_fl[i] = -db.getLatidudalLegPositionFrontLeft();
        r_fr[i] = db.getLatidudalLegPositionFrontRight();
        r_rl[i] = -db.getLatidudalLegPositionRearLeft();
        r_rr[i] = db.getLatidudalLegPositionRearRight();
        FC[i] += db.getBodyMass() * gy;
        FT_fl[i] += db.getTyreMassFrontLeft() * gy;
        FT_fr[i] += db.getTyreMassFrontRight() * gy;
        FT_rl[i] += db.getTyreMassRearLeft() * gy;
        FT_rr[i] += db.getTyreMassRearRight() * gy;
        FW_fl[i] += db.getWheelMassFrontLeft() * gy;
        FW_fr[i] += db.getWheelMassFrontRight() * gy;
        FW_rl[i] += db.getWheelMassRearLeft() * gy;
        FW_rr[i] += db.getWheelMassRearRight() * gy;

        i = 2;
        r_fl[i] = 0;
        r_fr[i] = 0;
        r_rl[i] = 0;
        r_rr[i] = 0;
        FC[i] += db.getBodyMass() * gz;
        FT_fl[i] += db.getTyreMassFrontLeft() * gz;
        FT_fr[i] += db.getTyreMassFrontRight() * gz;
        FT_rl[i] += db.getTyreMassRearLeft() * gz;
        FT_rr[i] += db.getTyreMassRearRight() * gz;
        FW_fl[i] += db.getWheelMassFrontLeft() * gz;
        FW_fr[i] += db.getWheelMassFrontRight() * gz;
        FW_rl[i] += db.getWheelMassRearLeft() * gz;
        FW_rr[i] += db.getWheelMassRearRight() * gz;
    }

    void getCholeskyDecomposition() {
        // A_Ic has cholesky factorization of Ic
        Math::copy<T>(Constants::DIMDIM, Ic, 1, A_Ic, 1);
        Math::potrf<T>(LAPACK_ROW_MAJOR, 'L', Constants::DIM, A_Ic, Constants::DIM);
    }

    /* Road Profile and Load module */

    void circular_path_initialization_quaternion(T* vcc, T* q) {
        T* unit_x_vector = Math::calloc<T>(Constants::DIM);
        T* normal_quaternion_vector = Math::calloc<T>(Constants::DIM);
        T quaternion_angle = 0;

        // unit x vector
        unit_x_vector[0] = 1;
        unit_x_vector[1] = 0;
        unit_x_vector[2] = 0;

        // calculate quaternion
        Math::GetQuaternion<T>(unit_x_vector, vcc, &quaternion_angle, normal_quaternion_vector);
        q[0] = normal_quaternion_vector[0] * std::sin(0.5 * quaternion_angle);
        q[1] = normal_quaternion_vector[1] * std::sin(0.5 * quaternion_angle);
        q[2] = normal_quaternion_vector[2] * std::sin(0.5 * quaternion_angle);
        q[3] = std::cos(0.5 * quaternion_angle);

        Math::free<T>(normal_quaternion_vector);
        Math::free<T>(unit_x_vector);
    }

    /** Updates the velocity of the 4 wheels and tyres as well as the angular
     * velocity, such that the car is already in the trajectory of the circle
     * overwrites the velocity in the tyres and wheels and the angular
     * velocities only keeps the tangential component of the velocity of the car
     * body
     */
    void circular_path_initialization(T* vc, T* vw_fl, T* vw_fr, T* vw_rl, T* vw_rr, T* vt_fl, T* vt_fr, T* vt_rl, T* vt_rr, T* omega, T* pcc, T* pt_fl, T* pt_fr, T* pt_rl, T* pt_rr) {
        auto& db = MetaDatabase<T>::getDatabase();

        // only consider circular motion in the XY-plane
        vc[2] = 0;

        // memory allocation
        T* perpendicular_dir = Math::calloc<T>(Constants::DIM);
        T* tangential_dir = Math::calloc<T>(Constants::DIM);
        T* radial_vector = Math::calloc<T>(Constants::DIM);

        // vector from the car to the center of the circle
        const auto center_of_circle = db.getCircularRoadCenter();
        radial_vector[0] = pcc[0] - center_of_circle[0];
        radial_vector[1] = pcc[1] - center_of_circle[1];
        radial_vector[2] = 0;

        // absolute distance from the car to the center
        T radius = Math::nrm2<T>(Constants::DIM, radial_vector, 1);

        // vector out of the XY-plane (unit Z direction)
        perpendicular_dir[2] = 1;

        T inv_radius = 1. / radius;

        // normalize the radial vector
        Math::scal<T>(Constants::DIM, inv_radius, radial_vector, Constants::INCX);

        // calculate the direction of the motion
        Math::CrossProduct<T>(radial_vector, perpendicular_dir, tangential_dir);

        // calculate the velocity magnitude
        T magnitude = Math::dot<T>(Constants::DIM, vc, Constants::INCX, tangential_dir, Constants::INCX);

        // get the velocity vector in tangential direction with computed
        // magnitude
        Math::copy<T>(Constants::DIM, tangential_dir, 1, vc, 1);
        Math::scal<T>(Constants::DIM, magnitude, vc, Constants::INCX);

        // get the physical radial vector again
        Math::scal<T>(Constants::DIM, radius, radial_vector, Constants::INCX);

        // calculate the angular velocity of the car
        Math::CrossProduct<T>(radial_vector, vc, omega);
        Math::scal<T>(Constants::DIM, inv_radius * inv_radius, omega, 1);

        // calculate the velocity in all legs
        Math::CrossProduct<T>(omega, pt_fl, vt_fl);
        Math::CrossProduct<T>(omega, pt_fr, vt_fr);
        Math::CrossProduct<T>(omega, pt_rl, vt_rl);
        Math::CrossProduct<T>(omega, pt_rr, vt_rr);

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
     * \param[out] Fr_fl The force acting on the fl tyre
     * \param[out] Fr_fr The force acting on the fr tyre
     * \param[out] Fr_rl The force acting on the rl tyre
     * \param[out] Fr_rr The force acting on the rr tyre
     */
    void get_fixed_road_force(T* Fr_fl, T* Fr_fr, T* Fr_rl, T* Fr_rr) {
        Fr_fl[2] = 0;
        Fr_fr[2] = 0;
        Fr_rl[2] = 0;
        Fr_rr[2] = 0;
    }

    /**
     * No interaction with the road, no additional forces on the tyres
     * \param[out] Fr_fl The force acting on the fl tyre
     * \param[out] Fr_fr The force acting on the fr tyre
     * \param[out] Fr_rl The force acting on the rl tyre
     * \param[out] Fr_rr The force acting on the rr tyre
     * TODO: implement method
     */
    void get_nonfixed_road_force(T* Fr_fl, T* Fr_fr, T* Fr_rl, T* Fr_rr) {}

    /**
     * calculates the force in the tyre only with respect to its velocity, mass
     * and position \param[out] Fr The only force acting on the tyre \param[in]
     * v is the velocity of the tyre \param[in] m is the mass of the tyre
     * \param[in] p is the global position of the tyre
     */
    void get_circular_road_force(T* Fr, const T* v, const T& m, const T* p) {
        T *unit_z_vector, *velocity_direction_tyre;
        T velocity_magnitude_tyre, inv_radius, force_magnitude_tyre;

        // allocate memory
        unit_z_vector = Math::calloc<T>(Constants::DIM);
        velocity_direction_tyre = Math::calloc<T>(Constants::DIM);

        // the direction vector from the center to the car in Fr
        Math::copy<T>(Constants::DIM, p, 1, Fr, 1);
        Math::axpy<T>(Constants::DIM, -1, MetaDatabase<T>::getDatabase().getCircularRoadCenter(), 1, Fr, 1);

        Fr[2] = 0;  // path only in XY-plane

        // inverse radius of the trajectory at the considered tyre (see TODO
        // above)
        inv_radius = 1. / Math::nrm2<T>(Constants::DIM, Fr, 1);

        // normalize the force vector -- minus sign because the force points to
        // the center
        Math::scal<T>(Constants::DIM, -inv_radius, Fr, 1);

        // perpendicular to the motion and radius
        unit_z_vector[0] = 0;
        unit_z_vector[1] = 0;
        unit_z_vector[2] = 1;

        // get normalized direction of motion
        Math::CrossProduct<T>(Fr, unit_z_vector, velocity_direction_tyre);

        // get the physical velocity
        velocity_magnitude_tyre = Math::dot<T>(Constants::DIM, v, Constants::INCX, velocity_direction_tyre, Constants::INCX);

        // get the physical force
        force_magnitude_tyre = m * velocity_magnitude_tyre * velocity_magnitude_tyre * inv_radius;
        Math::scal<T>(Constants::DIM, force_magnitude_tyre, Fr, 1);

        // free memory
        Math::free<T>(unit_z_vector);
        Math::free<T>(velocity_direction_tyre);
    }

    /**
     * \brief set the road forces to the one from the trajectory
     */
    void getArbitraryRoadForces(T* Fr_fl, T* Fr_fr, T* Fr_rl, T* Fr_rr, size_t i) {
        if (_iterationBegin) {
            _previousTyreForces[Constants::FRONT_LEFT] = Fr_fl[2];
            _previousTyreForces[Constants::FRONT_RIGHT] = Fr_fr[2];
            _previousTyreForces[Constants::REAR_LEFT] = Fr_rl[2];
            _previousTyreForces[Constants::REAR_RIGHT] = Fr_rr[2];
            _iterationBegin = false;
        }

        auto& db = MetaDatabase<T>::getDatabase();

        db.getArbitraryTrajectory()->getLagrangianForcesFrontLeft(i, db.getTyreMassFrontLeft(), Fr_fl);
        db.getArbitraryTrajectory()->getLagrangianForcesFrontRight(i, db.getTyreMassFrontRight(), Fr_fr);
        db.getArbitraryTrajectory()->getLagrangianForcesRearLeft(i, db.getTyreMassRearLeft(), Fr_rl);
        db.getArbitraryTrajectory()->getLagrangianForcesRearRight(i, db.getTyreMassRearRight(), Fr_rr);

        T delta_t = db.getTimeStepSize();

        flyingCarRoadForces(db.getTyreMassFrontLeft(), Fr_fl[2], Fr_fl[2], Constants::FRONT_LEFT, db.getArbitraryTrajectory()->getVerticalPositionFrontLeft(i - 1), (db.getArbitraryTrajectory()->getVerticalPositionFrontLeft(i) - db.getArbitraryTrajectory()->getVerticalPositionFrontLeft(i - 1)) / delta_t, db.getArbitraryTrajectory()->getVerticalRoadForcesFrontLeft(i, db.getTyreMassFrontLeft()), delta_t);

        flyingCarRoadForces(db.getTyreMassFrontRight(), Fr_fr[2], Fr_fr[2], Constants::FRONT_RIGHT, db.getArbitraryTrajectory()->getVerticalPositionFrontRight(i - 1), (db.getArbitraryTrajectory()->getVerticalPositionFrontRight(i) - db.getArbitraryTrajectory()->getVerticalPositionFrontRight(i - 1)) / delta_t, db.getArbitraryTrajectory()->getVerticalRoadForcesFrontRight(i, db.getTyreMassFrontRight()), delta_t);

        flyingCarRoadForces(db.getTyreMassRearLeft(), Fr_rl[2], Fr_rl[2], Constants::REAR_LEFT, db.getArbitraryTrajectory()->getVerticalPositionRearLeft(i - 1), (db.getArbitraryTrajectory()->getVerticalPositionRearLeft(i) - db.getArbitraryTrajectory()->getVerticalPositionRearLeft(i - 1)) / delta_t, db.getArbitraryTrajectory()->getVerticalRoadForcesRearLeft(i, db.getTyreMassRearLeft()), delta_t);

        flyingCarRoadForces(db.getTyreMassRearRight(), Fr_rr[2], Fr_rr[2], Constants::REAR_RIGHT, db.getArbitraryTrajectory()->getVerticalPositionRearRight(i - 1), (db.getArbitraryTrajectory()->getVerticalPositionRearRight(i) - db.getArbitraryTrajectory()->getVerticalPositionRearRight(i - 1)) / delta_t, db.getArbitraryTrajectory()->getVerticalRoadForcesRearRight(i, db.getTyreMassRearRight()), delta_t);
    }

    void flyingCarRoadForces(T mass, T FT, T& FR, size_t leg_index, T traj_pos, T traj_vel, T traj_force, T& delta_t) {
        T& previous_velocity = _previousSolution[20 + Constants::DIM * leg_index];
        T& previous_position = _previousSolution[51 + Constants::DIM * leg_index];

        T nothing = 0.0;

        // tyre is below or on the road
        if (previous_position <= traj_pos) {
            FR = traj_force;

           //push the car up (NOT PHYSICAL, just for nice output) TODO: make the collision physical
            FR += mass * (traj_pos - previous_position) / delta_t;


           // add a bump if the car has a velocity inside the road
            if (previous_velocity < std::min(nothing, traj_vel)) {
                FR += mass * (std::min(nothing, traj_vel) - previous_velocity) / delta_t;
            }

            // do not pull the car downwards
            if (traj_force < _previousTyreForces[leg_index]) {
                FR = FT;
            }
        }
    }
    /** Functions needed for compute_f */

    /**
     * Memory allocation of all the variables required in the Solve function
     * To increase performance by removing repeting memory allocations
     * The same locations are overwritten at each timestep
     */
    void compute_f_mem_alloc() {
        // add to members
        ///// All the members belong to compute f function and are kind of nasty
        //// cf_* means parameter used in compute f function
        cf_C_cN = Math::calloc<T>(Constants::DIMDIM);
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
        cf_wc_tilda = Math::calloc<T>(Constants::DIMDIM);

        cf_Tc[0] = 0;
        cf_Tc[1] = 0;
        cf_Tc[2] = 0;

        cf_b_rem = Math::calloc<T>(Constants::DIMDIM * Constants::DIM);
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
        //        std::cout << "wc: " << wc_[0] << ", " << wc_[1] << ", " << wc_[2] << std::endl;
    }

    void compute_spring_lengths(T* pcc_, T* pw_, T* pt_, T* cf_r_up_, T* cf_r_low_, T* r_, T& norm_r_up, T& inv_norm_r_up, T& norm_r_low, T& inv_norm_r_low) {
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
        Math::gemv<T>(CblasRowMajor, CblasNoTrans, Constants::DIM, Constants::DIM, 1, cf_C_cN, Constants::DIM, r_, 1, 1, cf_r_up_, 1);
       // std::cout << cf_r_up_[2] << ", ";
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
        basis_c = GetBasis(qc);
        */
        Math::GetBasis<T>(qc_, cf_C_cN);
    }

    /**
     * Computes the stiffness from lookup table.
     */
    void apply_stiffness_interpolation() {
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
        MetaDatabase<T>::getDatabase().getLookupStiffness().getInterpolation(current_spring_lengths, stiffness_vector);

        // overwrite stiffness values
        Math::copy(Constants::NUM_LEGS, stiffness_vector, 2, upper_spring_stiffness, 1);
        Math::copy(Constants::NUM_LEGS, stiffness_vector + 1, 2, lower_spring_stiffness, 1);
    }

    void compute_angles() {
        /*
        get angle and normal vectors at the legs

        [~, upper_angle1, upper_normal1] = GetQuaternion(r_up1, C_cN(2,:)');
        [~, upper_angle2, upper_normal2] = GetQuaternion(r_up2, C_cN(2,:)');
        [~, upper_angle3, upper_normal3] = GetQuaternion(r_up3, C_cN(2,:)');
        [~, upper_angle4, upper_normal4] = GetQuaternion(r_up4, C_cN(2,:)');

        [~, lower_angle1, lower_normal1] = GetQuaternion(r_low1, r_up1);
        [~, lower_angle2, lower_normal2] = GetQuaternion(r_low2, r_up2);
        [~, lower_angle3, lower_normal3] = GetQuaternion(r_low3, r_up3);
        [~, lower_angle4, lower_normal4] = GetQuaternion(r_low4, r_up4);

        */

        Math::copy<T>(Constants::DIM, cf_C_cN + 2, Constants::DIM, cf_col_dat, 1);

        Math::GetQuaternion<T>(cf_r_up_fl, cf_col_dat, cf_upper_angle_fl, cf_upper_normal_fl);
        Math::GetQuaternion<T>(cf_r_up_fr, cf_col_dat, cf_upper_angle_fr, cf_upper_normal_fr);
        Math::GetQuaternion<T>(cf_r_up_rl, cf_col_dat, cf_upper_angle_rl, cf_upper_normal_rl);
        Math::GetQuaternion<T>(cf_r_up_rr, cf_col_dat, cf_upper_angle_rr, cf_upper_normal_rr);

        Math::GetQuaternion<T>(cf_r_low_fl, cf_r_up_fl, cf_lower_angle_fl, cf_lower_normal_fl);
        Math::GetQuaternion<T>(cf_r_low_fr, cf_r_up_fr, cf_lower_angle_fr, cf_lower_normal_fr);
        Math::GetQuaternion<T>(cf_r_low_rl, cf_r_up_rl, cf_lower_angle_rl, cf_lower_normal_rl);
        Math::GetQuaternion<T>(cf_r_low_rr, cf_r_up_rr, cf_lower_angle_rr, cf_lower_normal_rr);
    }

    void compute_elongational_forces() {
        // Forces and Torques
        // calculate the elongational spring forces (in global basis)

        /*
        upper_force1 =
            upper_spring_stiffness(1) * (r_up1) * (1 - upper_spring_length(1) *
        inv_norm_r_up1); upper_force2 = upper_spring_stiffness(2) * (r_up2) * (1
        - upper_spring_length(2) * inv_norm_r_up2); upper_force3 =
            upper_spring_stiffness(3) * (r_up3) * (1 - upper_spring_length(3) *
        inv_norm_r_up3); upper_force4 = upper_spring_stiffness(4) * (r_up4) * (1
        - upper_spring_length(4) * inv_norm_r_up4);

        lower_force1 =
            lower_spring_stiffness(1) * (r_low1) * (1 - lower_spring_length(1) *
        inv_norm_r_low1); lower_force2 = lower_spring_stiffness(2) * (r_low2) *
        (1 - lower_spring_length(2) * inv_norm_r_low2); lower_force3 =
            lower_spring_stiffness(3) * (r_low3) * (1 - lower_spring_length(3) *
        inv_norm_r_low3); lower_force4 = lower_spring_stiffness(4) * (r_low4) *
        (1 - lower_spring_length(4) * inv_norm_r_low4);
        */

        T scale;
        Math::copy<T>(Constants::DIM, cf_r_up_fl, 1, cf_upper_force_fl, 1);
        scale = this->upper_spring_stiffness[0] * (1. - this->upper_spring_length[0] * inv_norm_r_up_fl);
        Math::scal<T>(Constants::DIM, scale, cf_upper_force_fl, 1);
        Math::copy<T>(Constants::DIM, cf_r_up_fr, 1, cf_upper_force_fr, 1);
        scale = this->upper_spring_stiffness[1] * (1. - this->upper_spring_length[1] * inv_norm_r_up_fr);
        Math::scal<T>(Constants::DIM, scale, cf_upper_force_fr, 1);
        Math::copy<T>(Constants::DIM, cf_r_up_rl, 1, cf_upper_force_rl, 1);
        scale = this->upper_spring_stiffness[2] * (1. - this->upper_spring_length[2] * inv_norm_r_up_rl);
        Math::scal<T>(Constants::DIM, scale, cf_upper_force_rl, 1);
        Math::copy<T>(Constants::DIM, cf_r_up_rr, 1, cf_upper_force_rr, 1);
        scale = this->upper_spring_stiffness[3] * (1. - this->upper_spring_length[3] * inv_norm_r_up_rr);
        Math::scal<T>(Constants::DIM, scale, cf_upper_force_rr, 1);

        Math::copy<T>(Constants::DIM, cf_r_low_fl, 1, cf_lower_force_fl, 1);
        scale = this->lower_spring_stiffness[0] * (1. - this->lower_spring_length[0] * inv_norm_r_low_fl);
        Math::scal<T>(Constants::DIM, scale, cf_lower_force_fl, 1);
        Math::copy<T>(Constants::DIM, cf_r_low_fr, 1, cf_lower_force_fr, 1);
        scale = this->lower_spring_stiffness[1] * (1. - this->lower_spring_length[1] * inv_norm_r_low_fr);
        Math::scal<T>(Constants::DIM, scale, cf_lower_force_fr, 1);
        Math::copy<T>(Constants::DIM, cf_r_low_rl, 1, cf_lower_force_rl, 1);
        scale = this->lower_spring_stiffness[2] * (1. - this->lower_spring_length[2] * inv_norm_r_low_rl);
        Math::scal<T>(Constants::DIM, scale, cf_lower_force_rl, 1);
        Math::copy<T>(Constants::DIM, cf_r_low_rr, 1, cf_lower_force_rr, 1);
        scale = this->lower_spring_stiffness[3] * (1. - this->lower_spring_length[3] * inv_norm_r_low_rr);
        Math::scal<T>(Constants::DIM, scale, cf_lower_force_rr, 1);
    }

    void compute_damping_forces() {
        /*
        calculate forces from damping effects
        upper_vdiff1 = (dot((vc - C_cN' * (r1_tilda * wc)), r_up1) - dot(vw1,
        r_up1)) * r_up1 * inv_norm_r_up1 * inv_norm_r_up1; upper_vdiff2 =
        (dot((vc - C_cN' * (r2_tilda * wc)), r_up2)
        - dot(vw2, r_up2)) * r_up2 * inv_norm_r_up2 * inv_norm_r_up2;
        upper_vdiff3 = (dot((vc - C_cN' * (r3_tilda * wc)), r_up3) - dot(vw3,
        r_up3)) * r_up3 * inv_norm_r_up3 * inv_norm_r_up3; upper_vdiff4 =
        (dot((vc - C_cN' * (r4_tilda * wc)), r_up4) - dot(vw4, r_up4)) * r_up4 *
        inv_norm_r_up4 * inv_norm_r_up4;

        lower_vdiff1 = (dot(vw1, r_low1) - dot(vt1, r_low1)) * r_low1 *
        inv_norm_r_low1 * inv_norm_r_low1; lower_vdiff2 = (dot(vw2, r_low2) -
        dot(vt2, r_low2)) * r_low2 * inv_norm_r_low2 * inv_norm_r_low2;
        lower_vdiff3 = (dot(vw3, r_low3) - dot(vt3, r_low3)) * r_low3 *
        inv_norm_r_low3 * inv_norm_r_low3; lower_vdiff4 = (dot(vw4, r_low4) -
        dot(vt4, r_low4)) * r_low4 * inv_norm_r_low4 * inv_norm_r_low4;

        upper_dampf1 = upper_spring_damping1(upper_vdiff1);
        upper_dampf2 = upper_spring_damping2(upper_vdiff2);
        upper_dampf3 = upper_spring_damping3(upper_vdiff3);
        upper_dampf4 = upper_spring_damping4(upper_vdiff4);

        lower_dampf1 = lower_spring_damping1(lower_vdiff1);
        lower_dampf2 = lower_spring_damping2(lower_vdiff2);
        lower_dampf3 = lower_spring_damping3(lower_vdiff3);
        lower_dampf4 = lower_spring_damping4(lower_vdiff4);

*/

        T scale;

        // upper_dampf1
        // compute: vc - C_cN' * (r1_tilda * wc))
        Math::gemv<T>(CblasRowMajor, CblasNoTrans, Constants::DIM, Constants::DIM, 1, this->r_fl_tilda, Constants::DIM, wc_, 1, 0, cf_temp, 1);
        Math::copy<T>(Constants::DIM, vc_, 1, cf_upper_dampf_fl, 1);
        Math::gemv<T>(CblasRowMajor, CblasNoTrans, Constants::DIM, Constants::DIM, -1, cf_C_cN, Constants::DIM, cf_temp, 1, 1, cf_upper_dampf_fl, 1);

        // dot((vc - C_cN' * (r1_tilda * wc)), r_up1)
        scale = Math::dot<T>(Constants::DIM, cf_upper_dampf_fl, Constants::INCX, cf_r_up_fl, Constants::INCX);

        // scale = Math::dot_product<T>(upper_dampf1, r_up1, Constants::DIM);

        // dot((vc - C_cN' * (r1_tilda * wc)), r_up1) - dot(vw1, r_up1)
        scale -= Math::dot<T>(Constants::DIM, vw_fl_, Constants::INCX, cf_r_up_fl, Constants::INCX);
        // scale -= Math::dot_product<T>(vw1_, r_up1, Constants::DIM);
        // (dot((vc - C_cN' * (r1_tilda * wc)), r_up1) - dot(vw1, r_up1))*
        // inv_norm_r_up1 * inv_norm_r_up1
        scale = scale * inv_norm_r_up_fl * inv_norm_r_up_fl * this->upper_spring_damping[0];
        Math::copy<T>(Constants::DIM, cf_r_up_fl, 1, cf_upper_dampf_fl, 1);

        Math::scal<T>(Constants::DIM, scale, cf_upper_dampf_fl, 1);

        //// upper_dampf_fr
        // compute: vc - C_cN' * (r_fr_tilda * wc))
        Math::gemv<T>(CblasRowMajor, CblasNoTrans, Constants::DIM, Constants::DIM, 1, this->r_fr_tilda, Constants::DIM, wc_, 1, 0, cf_temp, 1);
        Math::copy<T>(Constants::DIM, vc_, 1, cf_upper_dampf_fr, 1);
        Math::gemv<T>(CblasRowMajor, CblasNoTrans, Constants::DIM, Constants::DIM, -1, cf_C_cN, Constants::DIM, cf_temp, 1, 1, cf_upper_dampf_fr, 1);
        // dot((vc - C_cN' * (r_fr_tilda * wc)), r_up_fr)
        scale = Math::dot<T>(Constants::DIM, cf_upper_dampf_fr, Constants::INCX, cf_r_up_fr, Constants::INCX);
        // scale = Math::dot_product<T>(upper_dampf_fr, r_up_fr,
        // Constants::DIM); dot((vc - C_cN' * (r_fr_tilda * wc)), r_up_fr) -
        // dot(vw_fr, r_up_fr)
        scale -= Math::dot<T>(Constants::DIM, vw_fr_, Constants::INCX, cf_r_up_fr, Constants::INCX);
        // scale -= Math::dot_product<T>(vw_fr_, r_up_fr, Constants::DIM);
        // (dot((vc - C_cN' * (r_fr_tilda * wc)), r_up_fr) - dot(vw_fr,
        // r_up_fr))* inv_norm_r_up_fr
        // * inv_norm_r_up_fr
        scale = scale * inv_norm_r_up_fr * inv_norm_r_up_fr * this->upper_spring_damping[1];
        Math::copy<T>(Constants::DIM, cf_r_up_fr, 1, cf_upper_dampf_fr, 1);
        Math::scal<T>(Constants::DIM, scale, cf_upper_dampf_fr, 1);

        //// upper_dampf3
        // compute: vc - C_cN' * (r3_tilda * wc))
        Math::gemv<T>(CblasRowMajor, CblasNoTrans, Constants::DIM, Constants::DIM, 1, this->r_rl_tilda, Constants::DIM, wc_, 1, 0, cf_temp, 1);
        Math::copy<T>(Constants::DIM, vc_, 1, cf_upper_dampf_rl, 1);
        Math::gemv<T>(CblasRowMajor, CblasNoTrans, Constants::DIM, Constants::DIM, -1, cf_C_cN, Constants::DIM, cf_temp, 1, 1, cf_upper_dampf_rl, 1);
        // dot((vc - C_cN' * (r_rl_tilda * wc)), r_up_rl)
        scale = Math::dot<T>(Constants::DIM, cf_upper_dampf_rl, Constants::INCX, cf_r_up_rl, Constants::INCX);
        // scale = Math::dot_product<T>(upper_dampf_rl, r_up_rl,
        // Constants::DIM); dot((vc - C_cN' * (r_rl_tilda * wc)), r_up_rl) -
        // dot(vw_rl, r_up_rl)
        scale -= Math::dot<T>(Constants::DIM, vw_rl_, Constants::INCX, cf_r_up_rl, Constants::INCX);
        // scale -= Math::dot_product<T>(vw_rl_, r_up_rl, Constants::DIM);
        // (dot((vc - C_cN' * (r_rl_tilda * wc)), r_up_rl) - dot(vw_rl,
        // r_up_rl))* inv_norm_r_up_rl
        // * inv_norm_r_up_rl
        scale = scale * inv_norm_r_up_rl * inv_norm_r_up_rl * this->upper_spring_damping[2];
        Math::copy<T>(Constants::DIM, cf_r_up_rl, 1, cf_upper_dampf_rl, 1);
        Math::scal<T>(Constants::DIM, scale, cf_upper_dampf_rl, 1);

        //// upper_dampf4
        // compute: vc - C_cN' * (r4_tilda * wc))
        Math::gemv<T>(CblasRowMajor, CblasNoTrans, Constants::DIM, Constants::DIM, 1, this->r_rr_tilda, Constants::DIM, wc_, 1, 0, cf_temp, 1);
        Math::copy<T>(Constants::DIM, vc_, 1, cf_upper_dampf_rr, 1);
        Math::gemv<T>(CblasRowMajor, CblasNoTrans, Constants::DIM, Constants::DIM, -1, cf_C_cN, Constants::DIM, cf_temp, 1, 1, cf_upper_dampf_rr, 1);
        // dot((vc - C_cN' * (r4_tilda * wc)), r_up4)
        scale = Math::dot<T>(Constants::DIM, cf_upper_dampf_rr, Constants::INCX, cf_r_up_rr, Constants::INCX);
        // scale = Math::dot_product<T>(upper_dampf4, r_up4, Constants::DIM);
        // dot((vc - C_cN' * (r4_tilda * wc)), r_up4) - dot(vw4, r_up4)
        scale -= Math::dot<T>(Constants::DIM, vw_rr_, Constants::INCX, cf_r_up_rr, Constants::INCX);
        // scale -= Math::dot_product<T>(r_up4, vw4_, Constants::DIM);
        // (dot((vc - C_cN' * (r4_tilda * wc)), r_up4) - dot(vw4, r_up4))*
        // inv_norm_r_up4 * inv_norm_r_up4
        scale = scale * inv_norm_r_up_rr * inv_norm_r_up_rr * this->upper_spring_damping[3];
        Math::copy<T>(Constants::DIM, cf_r_up_rr, 1, cf_upper_dampf_rr, 1);
        Math::scal<T>(Constants::DIM, scale, cf_upper_dampf_rr, 1);

        //// lower_dampf1
        // dot(vw1, r_low1)
        scale = Math::dot<T>(Constants::DIM, vw_fl_, Constants::INCX, cf_r_low_fl, Constants::INCX);
        /*scale = Math::dot_product<T>(vw1_, r_low1, Constants::DIM);*/
        // (dot(vw1, r_low1) - dot(vt1, r_low1))
        scale -= Math::dot<T>(Constants::DIM, vt_fl_, Constants::INCX, cf_r_low_fl, Constants::INCX);
        // scale -= Math::dot_product<T>(vt1_, r_low1, Constants::DIM);
        //(dot(vw1, r_low1) - dot(vt1, r_low1)) * inv_norm_r_low1 *
        // inv_norm_r_low1
        scale = scale * inv_norm_r_low_fl * inv_norm_r_low_fl * this->lower_spring_damping[0];
        Math::copy<T>(Constants::DIM, cf_r_low_fl, 1, cf_lower_dampf_fl, 1);
        Math::scal<T>(Constants::DIM, scale, cf_lower_dampf_fl, 1);

        //// lower_dampf2
        // dot(vw2, r_low2)
        scale = Math::dot<T>(Constants::DIM, vw_fr_, Constants::INCX, cf_r_low_fr, Constants::INCX);
        // scale = Math::dot_product<T>(vw2_, r_low2, Constants::DIM);
        // (dot(vw2, r_low2) - dot(vt2, r_low2))
        scale -= Math::dot<T>(Constants::DIM, vt_fr_, Constants::INCX, cf_r_low_fr, Constants::INCX);
        // scale -= Math::dot_product<T>(vt2_, r_low2, Constants::DIM);
        //(dot(vw2, r_low2) - dot(vt2, r_low2)) * inv_norm_r_low2 *
        // inv_norm_r_low2
        scale = scale * inv_norm_r_low_fr * inv_norm_r_low_fr * this->lower_spring_damping[1];
        Math::copy<T>(Constants::DIM, cf_r_low_fr, 1, cf_lower_dampf_fr, 1);
        Math::scal<T>(Constants::DIM, scale, cf_lower_dampf_fr, 1);

        //// lower_dampf3
        // dot(vw3, r_low3)
        scale = Math::dot<T>(Constants::DIM, vw_rl_, Constants::INCX, cf_r_low_rl, Constants::INCX);
        // scale = Math::dot_product<T>(vw_rl_, r_low_rl, Constants::DIM);
        // (dot(vw_rl, r_low_rl) - dot(vt_rl, r_low_rl))
        scale -= Math::dot<T>(Constants::DIM, vt_rl_, Constants::INCX, cf_r_low_rl, Constants::INCX);
        // scale -= Math::dot_product<T>(vt_rl_, r_low_rl, Constants::DIM);
        //(dot(vw_rl, r_low_rl) - dot(vt_rl, r_low_rl)) * inv_norm_r_low_rl *
        // inv_norm_r_low_rl
        scale = scale * inv_norm_r_low_rl * inv_norm_r_low_rl * this->lower_spring_damping[2];
        Math::copy<T>(Constants::DIM, cf_r_low_rl, 1, cf_lower_dampf_rl, 1);
        Math::scal<T>(Constants::DIM, scale, cf_lower_dampf_rl, 1);

        //// lower_dampf4
        // dot(vw4, r_low4)
        scale = Math::dot<T>(Constants::DIM, vw_rr_, Constants::INCX, cf_r_low_rr, Constants::INCX);
        // scale = Math::dot_product<T>(vw4_, r_low4, Constants::DIM);
        // (dot(vw4, r_low4) - dot(vt4, r_low4))
        scale -= Math::dot<T>(Constants::DIM, vt_rr_, Constants::INCX, cf_r_low_rr, Constants::INCX);
        // scale -= Math::dot_product<T>(vt4_, r_low4, Constants::DIM);
        //(dot(vw4, r_low4) - dot(vt4, r_low4)) * inv_norm_r_low4 *
        // inv_norm_r_low4
        scale = scale * inv_norm_r_low_rr * inv_norm_r_low_rr * this->lower_spring_damping[3];
        Math::copy<T>(Constants::DIM, cf_r_low_rr, 1, cf_lower_dampf_rr, 1);
        Math::scal<T>(Constants::DIM, scale, cf_lower_dampf_rr, 1);
    }

    void compute_torques() {
        /*
torque from the rotational spring
upper_S1 = upper_rotational_stiffness(1) * upper_angle1 * upper_normal1; % in
global basis upper_S2 = upper_rotational_stiffness(2) * upper_angle2 *
upper_normal2; upper_S3 = upper_rotational_stiffness(3) * upper_angle3 *
upper_normal3; upper_S4 = upper_rotational_stiffness(4) * upper_angle4 *
upper_normal4;

lower_S1 = lower_rotational_stiffness(1) * lower_angle1 * lower_normal1; % in
global basis lower_S2 = lower_rotational_stiffness(2) * lower_angle2 *
lower_normal2; lower_S3 = lower_rotational_stiffness(3) * lower_angle3 *
lower_normal3; lower_S4 = lower_rotational_stiffness(4) * lower_angle4 *
lower_normal4;
*/

        T scale;

        // upper_S1 = upper_rotational_stiffness(1) * upper_angle1 *
        // upper_normal1;
        scale = (this->upper_rotational_stiffness[0]) * (*cf_upper_angle_fl);
        Math::copy<T>(Constants::DIM, cf_upper_normal_fl, 1, cf_upper_S_fl, 1);
        Math::scal<T>(Constants::DIM, scale, cf_upper_S_fl, 1);

        // upper_S2 = upper_rotational_stiffness(2) * upper_angle2 *
        // upper_normal2;
        scale = (this->upper_rotational_stiffness[1]) * (*cf_upper_angle_fr);
        Math::copy<T>(Constants::DIM, cf_upper_normal_fr, 1, cf_upper_S_fr, 1);
        Math::scal<T>(Constants::DIM, scale, cf_upper_S_fr, 1);

        // upper_S3 = upper_rotational_stiffness(3) * upper_angle3 *
        // upper_normal3;
        scale = (this->upper_rotational_stiffness[2]) * (*cf_upper_angle_rl);
        Math::copy<T>(Constants::DIM, cf_upper_normal_rl, 1, cf_upper_S_rl, 1);
        Math::scal<T>(Constants::DIM, scale, cf_upper_S_rl, 1);

        // upper_S4 = upper_rotational_stiffness(4) * upper_angle4 *
        // upper_normal4;
        scale = (this->upper_rotational_stiffness[3]) * (*cf_upper_angle_rr);
        Math::copy<T>(Constants::DIM, cf_upper_normal_rr, 1, cf_upper_S_rr, 1);
        Math::scal<T>(Constants::DIM, scale, cf_upper_S_rr, 1);

        // lower_S1 = lower_rotational_stiffness(1) * lower_angle1 *
        // lower_normal1
        scale = (this->lower_rotational_stiffness[0]) * (*cf_lower_angle_fl);
        Math::copy<T>(Constants::DIM, cf_lower_normal_fl, 1, cf_lower_S_fl, 1);
        Math::scal<T>(Constants::DIM, scale, cf_lower_S_fl, 1);

        // lower_S2 = lower_rotational_stiffness(2) * lower_angle2 *
        // lower_normal2
        scale = (this->lower_rotational_stiffness[1]) * (*cf_lower_angle_fr);
        Math::copy<T>(Constants::DIM, cf_lower_normal_fr, 1, cf_lower_S_fr, 1);
        Math::scal<T>(Constants::DIM, scale, cf_lower_S_fr, 1);

        // lower_S3 = lower_rotational_stiffness(3) * lower_angle3 *
        // lower_normal3
        scale = (this->lower_rotational_stiffness[2]) * (*cf_lower_angle_rl);
        Math::copy<T>(Constants::DIM, cf_lower_normal_rl, 1, cf_lower_S_rl, 1);
        Math::scal<T>(Constants::DIM, scale, cf_lower_S_rl, 1);

        // lower_S_rr = lower_rotational_stiffness(4) * lower_angle4 *
        // lower_normal4
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
         * sum_car_force1 = car_rot_force1 - upper_force1 - upper_dampf1 -
         * upper_rot_force1; sum_car_force2 = car_rot_force2 - upper_force2 -
         * upper_dampf2 - upper_rot_force2; sum_car_force3 = car_rot_force3 -
         * upper_force3 - upper_dampf3 - upper_rot_force3; sum_car_force4 =
         * car_rot_force4 - upper_force4 - upper_dampf4 - upper_rot_force4;
         */
        // TODO: constants

        T scale;
        T scale_u_fl, scale_u_fr, scale_u_rl, scale_u_rr;
        // lower_rot_force1 = -cross( lower_S1, r_low1) / (r_low1'*r_low1);
        scale = -1. / Math::dot<T>(Constants::DIM, cf_r_low_fl, Constants::INCX, cf_r_low_fl, Constants::INCX);
        // scale = -1. / Math::dot_product<T>(r_low1, r_low1, Constants::DIM);
        Math::CrossProduct<T>(cf_lower_S_fl, cf_r_low_fl, cf_lower_rot_force_fl);
        Math::scal<T>(Constants::DIM, scale, cf_lower_rot_force_fl, 1);

        // lower_rot_force2 = -cross( lower_S2, r_low2) / (r_low2'*r_low2);
        scale = -1. / Math::dot<T>(Constants::DIM, cf_r_low_fr, Constants::INCX, cf_r_low_fr, Constants::INCX);
        // scale = -1. / Math::dot_product<T>(r_low_fr, r_low_fr,
        // Constants::DIM);
        Math::CrossProduct<T>(cf_lower_S_fr, cf_r_low_fr, cf_lower_rot_force_fr);
        Math::scal<T>(Constants::DIM, scale, cf_lower_rot_force_fr, 1);

        // lower_rot_force3 = -cross( lower_S3, r_low3) / (r_low3'*r_low3);
        scale = -1. / Math::dot<T>(Constants::DIM, cf_r_low_rl, Constants::INCX, cf_r_low_rl, Constants::INCX);
        // scale = -1. / Math::dot_product<T>(r_low3, r_low3, Constants::DIM);
        Math::CrossProduct<T>(cf_lower_S_rl, cf_r_low_rl, cf_lower_rot_force_rl);
        Math::scal<T>(Constants::DIM, scale, cf_lower_rot_force_rl, 1);

        // lower_rot_force4 = -cross( lower_S4, r_low4) / (r_low4'*r_low4);
        scale = -1. / Math::dot<T>(Constants::DIM, cf_r_low_rr, Constants::INCX, cf_r_low_rr, Constants::INCX);
        // scale = -1. / Math::dot_product<T>(r_low4, r_low4, Constants::DIM);
        Math::CrossProduct<T>(cf_lower_S_rr, cf_r_low_rr, cf_lower_rot_force_rr);
        Math::scal<T>(Constants::DIM, scale, cf_lower_rot_force_rr, 1);

        // upper_rot_force1 = -cross( upper_S1, r_up1) / (r_up1'*r_up1);
        scale_u_fl = -1. / Math::dot<T>(Constants::DIM, cf_r_up_fl, Constants::INCX, cf_r_up_fl, Constants::INCX);
        // scale_u1 = -1. / Math::dot_product<T>(r_up1, r_up1, Constants::DIM);
        Math::CrossProduct<T>(cf_upper_S_fl, cf_r_up_fl, cf_upper_rot_force_fl);
        Math::scal<T>(Constants::DIM, scale_u_fl, cf_upper_rot_force_fl, 1);

        // upper_rot_force_fr = -cross( upper_S_fr, r_up_fr) /
        // (r_up_fr'*r_up_fr);
        scale_u_fr = -1. / Math::dot<T>(Constants::DIM, cf_r_up_fr, Constants::INCX, cf_r_up_fr, Constants::INCX);
        // scale_u_fr = -1. / Math::dot_product<T>(r_up_fr, r_up_fr,
        // Constants::DIM);
        Math::CrossProduct<T>(cf_upper_S_fr, cf_r_up_fr, cf_upper_rot_force_fr);
        Math::scal<T>(Constants::DIM, scale_u_fr, cf_upper_rot_force_fr, 1);

        // upper_rot_force3 = -cross( upper_S3, r_up3) / (r_up3'*r_up3);
        scale_u_rl = -1. / Math::dot<T>(Constants::DIM, cf_r_up_rl, Constants::INCX, cf_r_up_rl, Constants::INCX);
        // scale_u3 = -1. / Math::dot_product<T>(r_up3, r_up3, Constants::DIM);
        Math::CrossProduct<T>(cf_upper_S_rl, cf_r_up_rl, cf_upper_rot_force_rl);
        Math::scal<T>(Constants::DIM, scale_u_rl, cf_upper_rot_force_rl, 1);

        // upper_rot_force4 = -cross( upper_S4, r_up4) / (r_up4'*r_up4);
        scale_u_rr = -1. / Math::dot<T>(Constants::DIM, cf_r_up_rr, Constants::INCX, cf_r_up_rr, Constants::INCX);
        // scale_u4 = -1. / Math::dot_product<T>(r_up4, r_up4, Constants::DIM);
        Math::CrossProduct<T>(cf_upper_S_rr, cf_r_up_rr, cf_upper_rot_force_rr);
        Math::scal<T>(Constants::DIM, scale_u_rr, cf_upper_rot_force_rr, 1);

        // car_rot_force1 = -cross( lower_S1, r_up1) / (r_up1'*r_up1);
        Math::CrossProduct<T>(cf_lower_S_fl, cf_r_up_fl, cf_car_rot_force_fl);
        Math::scal<T>(Constants::DIM, scale_u_fl, cf_car_rot_force_fl, 1);

        // car_rot_force_fr = -cross( lower_S2, r_up2) / (r_up2'*r_up2);
        Math::CrossProduct<T>(cf_lower_S_fr, cf_r_up_fr, cf_car_rot_force_fr);
        Math::scal<T>(Constants::DIM, scale_u_fr, cf_car_rot_force_fr, 1);

        // car_rot_force3 = -cross( lower_S3, r_up3) / (r_up3'*r_up3);
        Math::CrossProduct<T>(cf_lower_S_rl, cf_r_up_rl, cf_car_rot_force_rl);
        Math::scal<T>(Constants::DIM, scale_u_rl, cf_car_rot_force_rl, 1);

        // car_rot_force4 = -cross( lower_S4, r_up4) / (r_up4'*r_up4);
        Math::CrossProduct<T>(cf_lower_S_rr, cf_r_up_rr, cf_car_rot_force_rr);
        Math::scal<T>(Constants::DIM, scale_u_rr, cf_car_rot_force_rr, 1);
    }

    void compute_car_body_forces() {
        // sum_car_force1 = car_rot_force1 - upper_force1 - upper_dampf1 -
        // upper_rot_force1;
        Math::copy<T>(Constants::DIM, cf_car_rot_force_fl, 1, cf_sum_car_force_fl, 1);
        Math::axpy<T>(Constants::DIM, -1, cf_upper_force_fl, 1, cf_sum_car_force_fl, 1);
        Math::axpy<T>(Constants::DIM, -1, cf_upper_dampf_fl, 1, cf_sum_car_force_fl, 1);
        Math::axpy<T>(Constants::DIM, -1, cf_upper_rot_force_fl, 1, cf_sum_car_force_fl, 1);

        // sum_car_force_fr = car_rot_force_fr - upper_force_fr - upper_dampf_fr
        // - upper_rot_force_fr;
        Math::copy<T>(Constants::DIM, cf_car_rot_force_fr, 1, cf_sum_car_force_fr, 1);
        Math::axpy<T>(Constants::DIM, -1, cf_upper_force_fr, 1, cf_sum_car_force_fr, 1);
        Math::axpy<T>(Constants::DIM, -1, cf_upper_dampf_fr, 1, cf_sum_car_force_fr, 1);
        Math::axpy<T>(Constants::DIM, -1, cf_upper_rot_force_fr, 1, cf_sum_car_force_fr, 1);

        // sum_car_force3 = car_rot_force3 - upper_force3 - upper_dampf3 -
        // upper_rot_force3;
        Math::copy<T>(Constants::DIM, cf_car_rot_force_rl, 1, cf_sum_car_force_rl, 1);
        Math::axpy<T>(Constants::DIM, -1, cf_upper_force_rl, 1, cf_sum_car_force_rl, 1);
        Math::axpy<T>(Constants::DIM, -1, cf_upper_dampf_rl, 1, cf_sum_car_force_rl, 1);
        Math::axpy<T>(Constants::DIM, -1, cf_upper_rot_force_rl, 1, cf_sum_car_force_rl, 1);

        // sum_car_force4 = car_rot_force4 - upper_force4 - upper_dampf4 -
        // upper_rot_force4;
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

        // local_FR1 = lower_force1 + lower_dampf1 + local_FT1 + local_FR1 +
        // lower_rot_force1; ... %vt1_dot
        Math::copy<T>(Constants::DIM, cf_lower_force_fl, 1, cf_local_FR_fl, 1);
        Math::axpy<T>(Constants::DIM, 1, cf_lower_dampf_fl, 1, cf_local_FR_fl, 1);
        Math::axpy<T>(Constants::DIM, 1, FT_fl, 1, cf_local_FR_fl, 1);
        //        Math::axpy<T>(Constants::DIM, 1, cf_local_FR_fl, 1, cf_local_FR_fl, 1);
        Math::axpy<T>(Constants::DIM, 1, cf_lower_rot_force_fl, 1, cf_local_FR_fl, 1);

        // local_FR2 = lower_force2 + lower_dampf2 + local_FT2 + local_FR2 +
        // lower_rot_force2; ... %vt2_dot
        Math::copy<T>(Constants::DIM, cf_lower_force_fr, 1, cf_local_FR_fr, 1);
        Math::axpy<T>(Constants::DIM, 1, cf_lower_dampf_fr, 1, cf_local_FR_fr, 1);
        Math::axpy<T>(Constants::DIM, 1, FT_fr, 1, cf_local_FR_fr, 1);
        //        Math::axpy<T>(Constants::DIM, 1, cf_local_FR_fr, 1, cf_local_FR_fr, 1);
        Math::axpy<T>(Constants::DIM, 1, cf_lower_rot_force_fr, 1, cf_local_FR_fr, 1);

        // local_FR3 = lower_force3 + lower_dampf3 + local_FT3 + local_FR3 +
        // lower_rot_force3; ... %vt3_dot
        Math::copy<T>(Constants::DIM, cf_lower_force_rl, 1, cf_local_FR_rl, 1);
        Math::axpy<T>(Constants::DIM, 1, cf_lower_dampf_rl, 1, cf_local_FR_rl, 1);
        Math::axpy<T>(Constants::DIM, 1, FT_rl, 1, cf_local_FR_rl, 1);
        //        Math::axpy<T>(Constants::DIM, 1, cf_local_FR_rl, 1, cf_local_FR_rl, 1);
        Math::axpy<T>(Constants::DIM, 1, cf_lower_rot_force_rl, 1, cf_local_FR_rl, 1);

        // local_FR4 = lower_force4 + lower_dampf4 + local_FT4 + local_FR4 +
        // lower_rot_force4]; %vt4_dot
        Math::copy<T>(Constants::DIM, cf_lower_force_rr, 1, cf_local_FR_rr, 1);
        Math::axpy<T>(Constants::DIM, 1, cf_lower_dampf_rr, 1, cf_local_FR_rr, 1);
        Math::axpy<T>(Constants::DIM, 1, FT_rr, 1, cf_local_FR_rr, 1);
        //        Math::axpy<T>(Constants::DIM, 1, cf_local_FR_rr, 1, cf_local_FR_rr, 1);
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
         * sum_torque_spring_car = r1_tilda * (C_cN * sum_car_force1) + ... %
         * from the elongational springs r2_tilda * (C_cN * sum_car_force2) +
         * ... r3_tilda * (C_cN * sum_car_force3) + ... r_rr_tilda * (C_cN *
         * sum_car_force4) + ... -C_cN * upper_S1 - C_cN * upper_S2 - C_cN *
         * upper_S3 - C_cN * upper_S4 + ... % ??from the rotational spring
         * -GetTilda(wc) * Hc + Tc;
         *    % from angular momentum and external torques
         */

        // Hc = A(1:3, 1:3) * wc;
        Math::gemv<T>(CblasRowMajor, CblasNoTrans, Constants::DIM, Constants::DIM, 1, this->Ic, Constants::DIM, wc_, 1, 0, cf_Hc, 1);
        Math::GetTilda<T>(wc_, cf_wc_tilda);
        Math::copy<T>(Constants::DIM, cf_Hc, 1, cf_temp, 1);
        Math::gemv<T>(CblasRowMajor, CblasNoTrans, Constants::DIM, Constants::DIM, -1, cf_wc_tilda, Constants::DIM, cf_temp, 1, 0, cf_Hc, 1);

        Math::copy<T>(Constants::DIM, cf_Hc, 1, cf_sum_torque_spring_car, 1);

        Math::axpy<T>(Constants::DIM, 1, cf_Tc, 1, cf_sum_torque_spring_car, 1);

        Math::gemv<T>(CblasRowMajor, CblasTrans, Constants::DIM, Constants::DIM, -1, cf_C_cN, Constants::DIM, cf_upper_S_rr, 1, 1, cf_sum_torque_spring_car, 1);
        Math::gemv<T>(CblasRowMajor, CblasTrans, Constants::DIM, Constants::DIM, -1, cf_C_cN, Constants::DIM, cf_upper_S_rl, 1, 1, cf_sum_torque_spring_car, 1);
        Math::gemv<T>(CblasRowMajor, CblasTrans, Constants::DIM, Constants::DIM, -1, cf_C_cN, Constants::DIM, cf_upper_S_fr, 1, 1, cf_sum_torque_spring_car, 1);
        Math::gemv<T>(CblasRowMajor, CblasTrans, Constants::DIM, Constants::DIM, -1, cf_C_cN, Constants::DIM, cf_upper_S_fl, 1, 1, cf_sum_torque_spring_car, 1);

        Math::gemv<T>(CblasRowMajor, CblasTrans, Constants::DIM, Constants::DIM, 1, cf_C_cN, Constants::DIM, cf_sum_car_force_fl, 1, 0, cf_temp, 1);
        Math::gemv<T>(CblasRowMajor, CblasNoTrans, Constants::DIM, Constants::DIM, 1, this->r_fl_tilda, Constants::DIM, cf_temp, 1, 1, cf_sum_torque_spring_car, 1);
        Math::gemv<T>(CblasRowMajor, CblasTrans, Constants::DIM, Constants::DIM, 1, cf_C_cN, Constants::DIM, cf_sum_car_force_fr, 1, 0, cf_temp, 1);
        Math::gemv<T>(CblasRowMajor, CblasNoTrans, Constants::DIM, Constants::DIM, 1, this->r_fr_tilda, Constants::DIM, cf_temp, 1, 1, cf_sum_torque_spring_car, 1);
        Math::gemv<T>(CblasRowMajor, CblasTrans, Constants::DIM, Constants::DIM, 1, cf_C_cN, Constants::DIM, cf_sum_car_force_rl, 1, 0, cf_temp, 1);
        Math::gemv<T>(CblasRowMajor, CblasNoTrans, Constants::DIM, Constants::DIM, 1, this->r_rl_tilda, Constants::DIM, cf_temp, 1, 1, cf_sum_torque_spring_car, 1);
        Math::gemv<T>(CblasRowMajor, CblasTrans, Constants::DIM, Constants::DIM, 1, cf_C_cN, Constants::DIM, cf_sum_car_force_rr, 1, 0, cf_temp, 1);
        Math::gemv<T>(CblasRowMajor, CblasNoTrans, Constants::DIM, Constants::DIM, 1, this->r_rr_tilda, Constants::DIM, cf_temp, 1, 1, cf_sum_torque_spring_car, 1);
    }

    void construct_right_hand_side() {
        /*
         * b = [sum_torque_spring_car; ...% w_dot_c
         * FC + sum_car_force1 + sum_car_force2 + sum_car_force3 +
         * sum_car_force4; ...% vc_dot upper_force1 - lower_force1 +
         * upper_dampf1 - lower_dampf1 + local_FW1 + ... upper_rot_force1 -
         * car_rot_force1 - lower_rot_force1; ...% vw1_dot upper_force2 -
         * lower_force2 + upper_dampf2 - lower_dampf2 + local_FW2 + ...
         * upper_rot_force2 - car_rot_force2 - lower_rot_force2; ...% vw2_dot
         * upper_force3 - lower_force3 + upper_dampf3 - lower_dampf3 + local_FW3
         * + ... upper_rot_force3 - car_rot_force3 - lower_rot_force3; ...%
         * vw3_dot upper_force4 - lower_force4 + upper_dampf4 - lower_dampf4 +
         * local_FW4 + ... upper_rot_force4 - car_rot_force4 - lower_rot_force4;
         * ...% vw4_dot lower_force1 + lower_dampf1 + local_FT1 + local_FR1 +
         * lower_rot_force1; ...% vt1_dot lower_force2 + lower_dampf2 +
         * local_FT2 + local_FR2 + lower_rot_force2; ...% vt2_dot lower_force3 +
         * lower_dampf3 + local_FT3 + local_FR3 + lower_rot_force3; ...% vt3_dot
         * lower_force4 + lower_dampf4 + local_FT4 + local_FR4 +
         * lower_rot_force4];% vt4_dot
         */
        // FC + sum_car_force1 + sum_car_force2 + sum_car_force3 +
        // sum_car_force4; ... %vc_dot
        T* brem_start = cf_b_rem;
        Math::copy<T>(Constants::DIM, this->FC, 1, brem_start, 1);
        Math::axpy<T>(Constants::DIM, 1, cf_sum_car_force_fl, 1, brem_start, 1);
        Math::axpy<T>(Constants::DIM, 1, cf_sum_car_force_fr, 1, brem_start, 1);
        Math::axpy<T>(Constants::DIM, 1, cf_sum_car_force_rl, 1, brem_start, 1);
        Math::axpy<T>(Constants::DIM, 1, cf_sum_car_force_rr, 1, brem_start, 1);

        //  upper_force1 - lower_force1 + upper_dampf1 - lower_dampf1 +
        //  local_FW1 + upper_rot_force1
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

        // upper_force2 - lower_force2 + upper_dampf2 - lower_dampf2 + local_FW2
        // + upper_rot_force2
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

        // upper_force3 - lower_force3 + upper_dampf3 - lower_dampf3 + local_FW3
        // + upper_rot_force3
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

        //  upper_force4 - lower_force4 + upper_dampf4 - lower_dampf4 +
        //  local_FW4 + upper_rot_force4
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
        Math::potrs<T>(LAPACK_ROW_MAJOR, 'L', Constants::DIM, 1, A_Ic, Constants::DIM, cf_sum_torque_spring_car, 1);

        Math::vMul<T>(Constants::VEC_DIM * Constants::DIM, A_rem, cf_b_rem, accelerations);
    }

    void compute_quaternion_change_rate() {
        /*
         * get the derivative of the altitude (expressed in quaternions) from
         * the angular velocities
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
        Math::gemv<T>(CblasRowMajor, CblasNoTrans, Constants::NUM_LEGS, Constants::DIM, 1, cf_Qc, Constants::DIM, wc_, 1, 0, cf_qc_dot, 1);
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
     * Constructor (with default parameters for writing in Output HDF5)
     */
    MBDMethod(const std::string& filePath = "", const std::string& fileName = "MBD_Checkpoints.hdf5") {
        MemoryAllocation();

        ReadFromXML();

        getCholeskyDecomposition();

#ifdef USE_HDF5
#ifdef USE_CHECKPOINTS
        _groupNameCheckpoints = "MBD Checkpoint t = ";

        _checkpointsMBD = new HDF5::OutputHDF5<T>(filePath, fileName);
        _checkpointsMBD->CreateContainer(true);
        _checkpointsMBD->CloseContainer();

        _checkpointsMBDFormatted = new HDF5::OutputHDF5<T>(filePath, "MBD_Checkpoints_formatted.hdf5");
        _checkpointsMBDFormatted->CreateContainer(true);
        _checkpointsMBDFormatted->CloseContainer();
#endif
#endif
    }

    /**
     * Initializes the time iteration and handles the numerical scheme
     */
    void Solve(T* solution_vector) {
        // From the formulation we have 61 dimensions in the solution vector
        auto& db = MetaDatabase<T>::getDatabase();
        size_t solution_size = (this->num_iter + 1) * Constants::MBD_SOLUTION_SIZE;
        T* complete_vector = Math::calloc<T>(solution_size);
        x_vector = Math::calloc<T>(Constants::MBD_SOLUTION_SIZE);
        Math::GetTilda<T>(r_fl, r_fl_tilda);
        Math::GetTilda<T>(r_fr, r_fr_tilda);
        Math::GetTilda<T>(r_rl, r_rl_tilda);
        Math::GetTilda<T>(r_rr, r_rr_tilda);

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
        Math::copy<T>(Constants::NUM_LEGS, MetaDatabase<T>::getDatabase().getBodyInitialOrientation(), 1, x_vector + i * (Constants::DIM) + j * (Constants::NUM_LEGS), 1);
        j++;
        // pcc
        pcc_ = x_vector + i * (Constants::DIM) + j * (Constants::NUM_LEGS);
        Math::copy<T>(Constants::DIM, pcc, 1, pcc_, 1);
        i++;
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

        // compute correct quaternion
        if (lagrangian_boundary_conditions == LagrangianBoundaryConditionRoad::CIRCULAR) circular_path_initialization_quaternion(vc, qc_);

        get_initial_length(qc_, r_fl, r_fr, r_rl, r_rr, pcc, pw_fl_, pw_fr_, pw_rl_, pw_rr_, pt_fl_, pt_fr_, pt_rl_, pt_rr_);

        // overwrites the initial velocity values
        if (lagrangian_boundary_conditions == LagrangianBoundaryConditionRoad::CIRCULAR) circular_path_initialization(vc, vw_fl, vw_fr, vw_rl, vw_rr, vt_fl, vt_fr, vt_rl, vt_rr, initial_angular_velocity, pcc_, pt_fl_, pt_fr_, pt_rl_, pt_rr_);
        if (lagrangian_boundary_conditions == LagrangianBoundaryConditionRoad::ARBITRARY) {
            T angle[3];
            db.getArbitraryTrajectory()->updateInitialConditions(angle, initial_angular_velocity, pcc_, vc, pt_fl_, pt_fr_, pt_rl_, pt_rr_, pw_fl_, pw_fr_, pw_rl_, pw_rr_, vt_fl, vt_fr, vt_rl, vt_rr, vw_fl, vw_fr, vw_rl, vw_rr);
            // update quaternion
            qc_[0] = 0;
            qc_[1] = 0;
            qc_[2] = -std::sin(angle[2] * 0.5);
            qc_[3] = -std::cos(angle[2] * 0.5);
        }
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

        // A_rem = 1 / x, for x in [CG W1 T1 W2 T2 W3 T3 W4 T4].
        for (auto i = 0; i < Constants::DIM; i++) {
            A_rem[i] = db.getBodyMass();
            Math::copy<T>(Constants::NUM_LEGS, db.getWheelTyreMassVector(), 2, A_rem + Constants::DIM + i, 3);
            Math::copy<T>(Constants::NUM_LEGS, db.getWheelTyreMassVector() + 1, 2, A_rem + (Constants::NUM_LEGS + 1) * Constants::DIM + i, 3);
        }
        Math::vInv<T>(Constants::VEC_DIM * Constants::DIM, A_rem, A_rem);

        compute_f_mem_alloc();

#ifdef WRITECSV
        IO::MyFile<T> solutionCSV("C:\\software\\repos\\EVAA\\output\\mbdSolution.txt");
        IO::MyFile<T> parametersCSV("C:\\software\\repos\\EVAA\\output\\simulationParameters.txt");
        parametersCSV.writeParameters();
#endif  // WRITECSV

        Math::copy<T>(Constants::MBD_SOLUTION_SIZE, x_vector, 1, _previousSolution, 1);
        _currentIteration = 1;
        _iterationBegin = true;

        if (used_solver == MBDSolver::BROYDEN_CN) {
            Math::Solvers<T, MBDMethod>::BroydenCN(this, x_vector, complete_vector, this->h, this->num_iter, this->tol, this->max_iter);
        }
        else if (used_solver == MBDSolver::RUNGE_KUTTA_4) {
            Math::Solvers<T, MBDMethod>::RK4(this, x_vector, complete_vector, this->h, this->num_iter, this->tol, this->max_iter);
        }
        else if (used_solver == MBDSolver::BROYDEN_BDF2) {
            Math::Solvers<T, MBDMethod>::BroydenBDF2(this, x_vector, complete_vector, this->h, this->num_iter, this->tol, this->max_iter);
        }
        else if (used_solver == MBDSolver::BROYDEN_EULER) {
            Math::Solvers<T, MBDMethod>::BroydenEuler(this, x_vector, complete_vector, this->h, this->num_iter, this->tol, this->max_iter);
        }
        else if (used_solver == MBDSolver::EXPLICIT_EULER) {
            std::cout << "Explicit solver hasn't been implemented, you don't "
                         "want to use it"
                      << std::endl;
        }
        else {
            std::cout << "sorry man, the solver you picked for MBD is weird "
                         "and hasn't been implemented yet"
                      << std::endl;
        }

        compute_f_clean();

#ifdef WRITECSV
        solutionCSV.writeSolutionMatrixMBD(complete_vector, this->num_iter + 1);
#endif  // WRITECSV

#ifdef USE_HDF5
        WriteBulkResults(complete_vector, num_iter + 1, Constants::MBD_SOLUTION_SIZE);
#endif  // USE_HDF5

#ifdef USE_HDF5  // TODO should be USE_CHECKPOINTS
#ifdef USE_CHECKPOINTS
        // write at checkpoints
        // TODO checkpoints vector read from somewhere (contains the iteration numbers for checkpoints)
        const size_t numCheckpoints = 3;
        size_t checkpoints[numCheckpoints];  // to be read from somewhere
        checkpoints[0] = 0.25 * num_iter;
        checkpoints[1] = 0.5 * num_iter;
        checkpoints[2] = 0.75 * num_iter;

        for (auto i = 0; i < numCheckpoints; ++i) {
            _checkpointsMBD->CreateContainer(false, _groupNameCheckpoints + std::to_string((T)checkpoints[i] * h));
            _checkpointsMBD->WriteVector("MBD Checkpoint t = " + std::to_string((T)checkpoints[i] * h), &complete_vector[checkpoints[i]], Constants::MBD_SOLUTION_SIZE, HDF5FileHandle::GROUP);
            _checkpointsMBD->CloseContainer();
            // TODO : decide if we want with format or without
            _checkpointsMBDFormatted->CreateContainer(false, _groupNameCheckpoints + std::to_string((T)checkpoints[i] * h));
            WriteFormattedSolution(_checkpointsMBDFormatted, &complete_vector[checkpoints[i]], HDF5FileHandle::GROUP);
            _checkpointsMBDFormatted->CloseContainer();
        }
#endif
#endif
        
        T* start = complete_vector + (this->num_iter) * Constants::MBD_SOLUTION_SIZE;
        Math::copy<T>(Constants::MBD_SOLUTION_SIZE, start, 1, solution_vector, 1);
        Math::free<T>(complete_vector);
        Math::free<T>(x_vector);
    }

    /**
     * \brief here comes everything which has to be done at the end of one time iteration
     */
    void postprocessingTimeIteration(size_t iteration, T* solutionPreviousTimestep) {
        _currentIteration = iteration + 1;
        Math::copy<T>(Constants::MBD_SOLUTION_SIZE, solutionPreviousTimestep, 1, _previousSolution, 1);
        _iterationBegin = true;
    }

    /**
     * Solver which is called at each time step
     * Computes the forces and torques on each point mass and computes the right
     * hand side of the ODE \param[in] x_ current solution of the system
     * \param[in] t_ current simulation time
     * \param [in] iteration count of the current iteration
     * \param[out] f_ the load vector
     */
    void compute_f3D_reduced(T* x_, T t_, T* f_, size_t iteration) {
        /*
         * Small performance gain might be possible by transforming C_cN to
         * column major Note: corresponding MKL function call have to be changed
         * too
         */

        get_current_variables(x_);

        get_car_orientation();

  //      std::cout << "positions: ";
        compute_spring_lengths(pcc_, pw_fl_, pt_fl_, cf_r_up_fl, cf_r_low_fl, this->r_fl, norm_r_up_fl, inv_norm_r_up_fl, norm_r_low_fl, inv_norm_r_low_fl);
        compute_spring_lengths(pcc_, pw_fr_, pt_fr_, cf_r_up_fr, cf_r_low_fr, this->r_fr, norm_r_up_fr, inv_norm_r_up_fr, norm_r_low_fr, inv_norm_r_low_fr);
        compute_spring_lengths(pcc_, pw_rl_, pt_rl_, cf_r_up_rl, cf_r_low_rl, this->r_rl, norm_r_up_rl, inv_norm_r_up_rl, norm_r_low_rl, inv_norm_r_low_rl);
        compute_spring_lengths(pcc_, pw_rr_, pt_rr_, cf_r_up_rr, cf_r_low_rr, this->r_rr, norm_r_up_rr, inv_norm_r_up_rr, norm_r_low_rr, inv_norm_r_low_rr);
    //    std::cout << std::endl;

#ifdef INTERPOLATION
        apply_stiffness_interpolation();
#endif

        compute_angles();

        compute_elongational_forces();

        compute_damping_forces();

        compute_torques();

        compute_torques_resultant_forces();

        compute_car_body_forces();

        compute_external_forces();
        // TODO change this if else shit to something nice (function pointers maybe)
        if ((eulerian_boundary_conditions == EulerianBoundaryConditionRoad::FIXED) && (lagrangian_boundary_conditions == LagrangianBoundaryConditionRoad::STRAIGHT)) {
            get_fixed_road_force(cf_local_FR_fl, cf_local_FR_fr, cf_local_FR_rl, cf_local_FR_rr);
        }
        else if ((eulerian_boundary_conditions == EulerianBoundaryConditionRoad::NONFIXED) && (lagrangian_boundary_conditions == LagrangianBoundaryConditionRoad::STRAIGHT)) {
            get_nonfixed_road_force(cf_local_FR_fl, cf_local_FR_fr, cf_local_FR_rl, cf_local_FR_rr);
        }
        else if ((eulerian_boundary_conditions == EulerianBoundaryConditionRoad::FIXED) && (lagrangian_boundary_conditions == LagrangianBoundaryConditionRoad::CIRCULAR)) {
            auto& db = MetaDatabase<T>::getDatabase();
            get_circular_road_force(cf_local_FR_fl, vt_fl_, db.getTyreMassFrontLeft(), pt_fl_);
            get_circular_road_force(cf_local_FR_fr, vt_fr_, db.getTyreMassFrontRight(), pt_fr_);
            get_circular_road_force(cf_local_FR_rl, vt_rl_, db.getTyreMassRearLeft(), pt_rl_);
            get_circular_road_force(cf_local_FR_rr, vt_rr_, db.getTyreMassRearRight(), pt_rr_);
        }
        else if ((eulerian_boundary_conditions == EulerianBoundaryConditionRoad::SINUSOIDAL) && (lagrangian_boundary_conditions == LagrangianBoundaryConditionRoad::ARBITRARY)) {
            getArbitraryRoadForces(cf_local_FR_fl, cf_local_FR_fr, cf_local_FR_rl, cf_local_FR_rr, iteration);
        }
        else {
            throw std::logic_error(
                "The MBD simulation currently only supports following combinations of road conditions: \n"
                "   - fixed / circular \n"
                "   - fixed / straight\n"
                "   - detached / straight\n"
                "   - sinusoidal / arbitrary");
        }

        compute_car_body_total_torque();

        construct_right_hand_side();

        solve_accelerations();

        compute_quaternion_change_rate();

        construct_f_vector(f_);
    }

    size_t get_solution_dimension() { return Constants::MBD_SOLUTION_SIZE; }

#ifdef USE_HDF5
    /** Write a solution vector in formatted fashion.
     * \param[in] handle if HDF5FileHandle::FILE, then writes in file; if HDF5FileHandle::GROUP, then writes in group
     */
    void WriteFormattedSolution(HDF5::OutputHDF5<T>* writeMBD, T* sln, const HDF5FileHandle& handle = HDF5FileHandle::FILE) {
        std::string datasetName = "MBD: angular velocity w";
        writeMBD->WriteVector(datasetName, &sln[0], Constants::DIM, handle);
        datasetName = "MBD: car body velocity vc";
        writeMBD->WriteVector(datasetName, &sln[3], Constants::DIM, handle);
        // Wheels velocities
        datasetName = "MBD: front-left wheel velocity vw_fl";
        writeMBD->WriteVector(datasetName, &sln[6], Constants::DIM, handle);
        datasetName = "MBD: front-right wheel velocity vw_fr";
        writeMBD->WriteVector(datasetName, &sln[9], Constants::DIM, handle);
        datasetName = "MBD: rear-left wheel velocity vw_rl";
        writeMBD->WriteVector(datasetName, &sln[12], Constants::DIM, handle);
        datasetName = "MBD: rear-right wheel velocity vw_rr";
        writeMBD->WriteVector(datasetName, &sln[15], Constants::DIM, handle);
        // Tyres  velocities
        datasetName = "MBD: front-left tyre velocity vt_fl";
        writeMBD->WriteVector(datasetName, &sln[18], Constants::DIM, handle);
        datasetName = "MBD: front-right tyre velocity vt_fr";
        writeMBD->WriteVector(datasetName, &sln[21], Constants::DIM, handle);
        datasetName = "MBD: rear-left tyre velocity vt_rl";
        writeMBD->WriteVector(datasetName, &sln[24], Constants::DIM, handle);
        datasetName = "MBD: rear-right tyre velocity vt_rr";
        writeMBD->WriteVector(datasetName, &sln[27], Constants::DIM, handle);
        // Quaternion orientation
        datasetName = "MBD: orientation q";
        writeMBD->WriteVector(datasetName, &sln[30], Constants::DIM + 1, handle);
        datasetName = "MBD: car body position pc";
        writeMBD->WriteVector(datasetName, &sln[34], Constants::DIM, handle);
        // Wheels positions
        datasetName = "MBD: front-left wheel position pw_fl";
        writeMBD->WriteVector(datasetName, &sln[37], Constants::DIM, handle);
        datasetName = "MBD: front-right wheel position pw_f";
        writeMBD->WriteVector(datasetName, &sln[40], Constants::DIM, handle);
        datasetName = "MBD: rear-left wheel position pw_rl";
        writeMBD->WriteVector(datasetName, &sln[43], Constants::DIM, handle);
        datasetName = "MBD: rear-right wheel position pw_rr";
        writeMBD->WriteVector(datasetName, &sln[46], Constants::DIM, handle);
        // Tyres positions
        datasetName = "MBD: front-left tyre position pt_fl";
        writeMBD->WriteVector(datasetName, &sln[49], Constants::DIM, handle);
        datasetName = "MBD: front-right tyre position pt_fr";
        writeMBD->WriteVector(datasetName, &sln[52], Constants::DIM, handle);
        datasetName = "MBD: rear-left tyre position pt_rl";
        writeMBD->WriteVector(datasetName, &sln[55], Constants::DIM, handle);
        datasetName = "MBD: rear-right tyre position pt_rr";
        writeMBD->WriteVector(datasetName, &sln[58], Constants::DIM, handle);
    }
#endif

/** Write the Formatted Final Solution for MBD*/
#ifdef USE_HDF5
    void WriteFinalResultFormatted(T* sln, const std::string filePath = "", const std::string fileName = "MBD_final_solution_formatted.hdf5") {
        HDF5::OutputHDF5<T> finalMBD(filePath, fileName);
        finalMBD.CreateContainer(true);
        WriteFormattedSolution(&finalMBD, sln, HDF5FileHandle::FILE);
        finalMBD.CloseContainer();
    }
#endif  // USE_HDF5

/** Write the Bulk Final Solution for MBD*/
#ifdef USE_HDF5
    void WriteFinalResult(T* sln, const std::string filePath = "", const std::string fileName = "MBD_final_solution_bulk.hdf5") {
        HDF5::OutputHDF5<Constants::floatEVAA> finalMBD(filePath, fileName);
        std::string datasetName = "MBD final solution";
        finalMBD.CreateContainer(true);
        finalMBD.WriteVector(datasetName, &sln[0], Constants::MBD_SOLUTION_SIZE);
        finalMBD.CloseContainer();
    }
#endif  // USE_HDF5

#ifdef USE_HDF5
    void WriteBulkResults(T* fullSolution, const size_t& numIterations, const size_t& solutionDimension, const std::string& filePath = "", const std::string& fileName = "MBD_full_solution.hdf5") {
        HDF5::OutputHDF5<Constants::floatEVAA> finalMBD(filePath, fileName);
        std::string datasetName = "MBD bulk solution";
        finalMBD.CreateContainer(true);
        finalMBD.WriteMatrix(datasetName, fullSolution, numIterations, solutionDimension);
        finalMBD.CloseContainer();
    }
#endif  // USE_HDF5

    /**
     * Beautiful output of the result
     * \param sln solution vector
     */
    void PrintFinalResult(T* sln) {
        /*std::cout << std::scientific;
        std::cout << std::setprecision(15);*/

        std::cout << "MBD: angular velocity w =\n\t[" << sln[0] << "\n\t " << sln[1] << "\n\t " << sln[2] << "]" << std::endl;
        std::cout << "MBD: car body velocity vc =\n\t[" << sln[3] << "\n\t " << sln[4] << "\n\t " << sln[5] << "]" << std::endl;
        std::cout << "MBD: front-left wheel velocity vw_fl =\n\t[" << sln[6] << "\n\t " << sln[7] << "\n\t " << sln[8] << "]" << std::endl;
        std::cout << "MBD: front-right wheel velocity vw_fr =\n\t[" << sln[9] << "\n\t " << sln[10] << "\n\t " << sln[11] << "]" << std::endl;
        std::cout << "MBD: rear-left wheel velocity vw_rl =\n\t[" << sln[12] << "\n\t " << sln[13] << "\n\t " << sln[14] << "]" << std::endl;
        std::cout << "MBD: rear-right wheel velocity vw_rr =\n\t[" << sln[15] << "\n\t " << sln[16] << "\n\t " << sln[17] << "]" << std::endl;
        std::cout << "MBD: front-left tyre velocity vt_fl =\n\t[" << sln[18] << "\n\t " << sln[19] << "\n\t " << sln[20] << "]" << std::endl;
        std::cout << "MBD: front-right tyre velocity vt_fr =\n\t[" << sln[21] << "\n\t " << sln[22] << "\n\t " << sln[23] << "]" << std::endl;
        std::cout << "MBD: rear-left tyre velocity vt_rl =\n\t[" << sln[24] << "\n\t " << sln[25] << "\n\t " << sln[26] << "]" << std::endl;
        std::cout << "MBD: rear-right tyre velocity vt_rr =\n\t[" << sln[27] << "\n\t " << sln[28] << "\n\t " << sln[29] << "]" << std::endl;
        std::cout << "MBD: orientation q =\n\t[" << sln[30] << "\n\t " << sln[31] << "\n\t " << sln[32] << "\n\t " << sln[33] << "]" << std::endl;
        std::cout << "MBD: car body position pc =\n\t[" << sln[34] << "\n\t " << sln[35] << "\n\t " << sln[36] << "]" << std::endl;
        std::cout << "MBD: front-left wheel position pw_fl =\n\t[" << sln[37] << "\n\t " << sln[38] << "\n\t " << sln[39] << "]" << std::endl;
        std::cout << "MBD: front-right wheel position pw_fr =\n\t[" << sln[40] << "\n\t " << sln[41] << "\n\t " << sln[42] << "]" << std::endl;
        std::cout << "MBD: rear-left wheel position pw_rl =\n\t[" << sln[43] << "\n\t " << sln[44] << "\n\t " << sln[45] << "]" << std::endl;
        std::cout << "MBD: rear-right wheel position pw_rr =\n\t[" << sln[46] << "\n\t " << sln[47] << "\n\t " << sln[48] << "]" << std::endl;
        std::cout << "MBD: front-left tyre position pt_fl =\n\t[" << sln[49] << "\n\t " << sln[50] << "\n\t " << sln[51] << "]" << std::endl;
        std::cout << "MBD: front-right tyre position pt_fr =\n\t[" << sln[52] << "\n\t " << sln[53] << "\n\t " << sln[54] << "]" << std::endl;
        std::cout << "MBD: rear-left tyre position pt_rl =\n\t[" << sln[55] << "\n\t " << sln[56] << "\n\t " << sln[57] << "]" << std::endl;
        std::cout << "MBD: rear-right tyre position pt_rr =\n\t[" << sln[58] << "\n\t " << sln[59] << "\n\t " << sln[60] << "]" << std::endl;
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
        Math::free<T>(upper_spring_length);
        Math::free<T>(lower_spring_length);
        Math::free<T>(upper_spring_stiffness);
        Math::free<T>(lower_spring_stiffness);
        Math::free<T>(upper_rotational_stiffness);
        Math::free<T>(lower_rotational_stiffness);
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
        Math::free<T>(_previousTyreForces);
        Math::free<T>(_previousSolution);

#ifdef USE_HDF5  // TODO: use USE_CHECKPOINT
#ifdef USE_CHECKPOINTS
        delete _checkpointsMBD;
        delete _checkpointsMBDFormatted;
#endif
#endif
    }
};

}  // namespace EVAA
