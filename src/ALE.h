// TODO: Copyright header

#pragma once
#include "11DOF.h"
#include "BLAS.h"
#include "Constants.h"
#include "EVAAComputeEngine.h"
#include "LoadModule.h"
#include "RoadProfile.h"

namespace EVAA {

/** Implements the ALE method to extend the linear 11DOF system */
template <class T>
class ALE {
private:
    // necessary class objects
    Car<T>* Car_obj;                 // suppose Interpolation in the Car
    LoadModule<T>* Load_module_obj;  // needs Profile and Car
    TwoTrackModelParent<T>* TwoTrackModel_obj;
    EVAALookup<Constants::floatEVAA>* interpolator;  // interpolator member of EVAA
    Profile<T>* profile_obj;

    // simulation parameters
    int DOF;
    T tend_;
    T h_;

    // time and solution vectors
    T* time_vec;
    T* u_sol;

    // needed to solve the 11DOF system
    T* force_vector_11dof;

    // Spring stiffnesses in Stefan's ordering
    T* k_vect;

    // Contains the forces in the ordering [GC:XYZ,W1:XYZ,T1:XYZ, ...]
    T* force_vector;

    // Contains he torque in the form [XYZ]
    T* full_torque;

    // Contains the dx in the spring elongation in Stefan's ordering
    T* Delta_x_vec;

    // global centripetal forces on the whole car [XYZ]
    T* centripetal_force;
    T* new_centripetal_force;

    // global torque on the car [XYZ]
    T* torque;
    T* new_torque;

    // global positions, velocities [XY] and angles [Z] of the center of mass of the car
    T* posXY_vec;
    T* angleZ;
    T* velXY_vec;
    T* ang_velZ;

    // quantities for the whole car
    T global_inertia_Z;

public:
    /**
     * Constructor
     */
    ALE(Car<T>* Car_obj_val, LoadModule<T>* Load_module_val,
        TwoTrackModelParent<T>* TwoTrackModel_val, EVAALookup<T>* lookup_table)
    {
        Car_obj = Car_obj_val;
        Load_module_obj = Load_module_val;
        TwoTrackModel_obj = TwoTrackModel_val;
        interpolator = lookup_table;

        // general parameters of the simulation
        DOF = MetaDataBase::DataBase()->getDOF();

        h_ = MetaDataBase::DataBase()->getTimeStepSize();
        tend_ = MetaDataBase::DataBase()->getNumberOfTimeIterations() * h_;
    }

    /**
     * Applies the Verlet_Stoermer algorithm to update the global XY position of the car and its Z
     * orientation
     * Store the global coordinates in the VelocityXY and PositionXY from the car object
     * \param t current simulation time
     */
    void global_frame_solver(T& t)
    {
        // 2. Update global X,Y positions of the car
        MathLibrary::Solvers<T, ALE>::Stoermer_Verlet_Position(
            Car_obj->currentPositionLagrangian[0], Car_obj->currentVelocityLagrangian[0],
            centripetal_force[0], h_, *Car_obj->massFullCar);
        MathLibrary::Solvers<T, ALE>::Stoermer_Verlet_Position(
            Car_obj->currentPositionLagrangian[1], Car_obj->currentVelocityLagrangian[1],
            centripetal_force[1], h_, *Car_obj->massFullCar);

        // 4. Update Z-rotation
        MathLibrary::Solvers<T, ALE>::Stoermer_Verlet_Position(
            *Car_obj->currentAngleLagrangian, *Car_obj->currentAngularVelocityLagrangian, torque[2],
            h_, global_inertia_Z);

        // get forces
        Load_module_obj->update_force(t, force_vector, new_centripetal_force);
        Load_module_obj->update_torque(t, new_torque, force_vector);

        mkl<T>::scal(2, -1, new_centripetal_force, 1);

        // 1. Update global X,Y velocities
        MathLibrary::Solvers<T, ALE>::Stoermer_Verlet_Velocity(
            Car_obj->currentVelocityLagrangian[0], centripetal_force[0], new_centripetal_force[0],
            h_, *Car_obj->massFullCar);
        MathLibrary::Solvers<T, ALE>::Stoermer_Verlet_Velocity(
            Car_obj->currentVelocityLagrangian[1], centripetal_force[1], new_centripetal_force[1],
            h_, *Car_obj->massFullCar);

        // 3. Update Z-angular velocities
        MathLibrary::Solvers<T, ALE>::Stoermer_Verlet_Velocity(
            *Car_obj->currentAngularVelocityLagrangian, torque[2], new_torque[2], h_,
            global_inertia_Z);

        /*
         * Idea!! What if we do it at the end, since the displacement is a vector and by triangle
         * rule sum of all should add up force is computed using on
         */
        Car_obj->apply_ALE_change();

        // update forces and torque
        centripetal_force[0] = new_centripetal_force[0];
        centripetal_force[1] = new_centripetal_force[1];

        torque[2] = new_torque[2];  // z - component
    }
    void solve(T* sol_vect, T* u_sol)
    {
        // initialize solution vector
        int sol_size = (floor(tend_ / h_) + 1);
        const int force_dimensions = Constants::DIM * Constants::VEC_DIM;
        int centripetal_force_dimensions = 3;  // because update force needs it
        int full_torque_dimensions = 3;
        const int num_springs = 2 * Constants::NUM_LEGS;

        // allocate memory

        time_vec = (T*)mkl_calloc(sol_size, sizeof(T), Constants::ALIGNMENT);
        u_sol = (T*)mkl_calloc(sol_size * (Constants::VEC_DIM * Constants::DIM), sizeof(T),
                               Constants::ALIGNMENT);
        force_vector = (T*)mkl_calloc(force_dimensions, sizeof(T), Constants::ALIGNMENT);
        force_vector_11dof = (T*)mkl_calloc(DOF, sizeof(T), Constants::ALIGNMENT);
        full_torque = (T*)mkl_calloc(full_torque_dimensions, sizeof(T), Constants::ALIGNMENT);
        centripetal_force =
            (T*)mkl_calloc(centripetal_force_dimensions, sizeof(T), Constants::ALIGNMENT);
        new_centripetal_force =
            (T*)mkl_calloc(centripetal_force_dimensions, sizeof(T),
                           Constants::ALIGNMENT);  // this was 2 dimensional allocation and update
                                                   // force updates 3 dimension on this
        k_vect = (T*)mkl_calloc(num_springs, sizeof(T), Constants::ALIGNMENT);
        Delta_x_vec = (T*)mkl_calloc(num_springs, sizeof(T), Constants::ALIGNMENT);

        torque = new T[3];
        new_torque = new T[3];

        // calculate characteristics of the whole car
        calculate_global_inertia_Z();

        // start time iteration
        T t = h_;

        T* solution_vect;
        int iter = 1;
        // time iteration
        double eps = h_ / 100;
        while (std::abs(t - (tend_ + h_)) > eps) {
            // This has to be done at each time step
            //
            // update force vector
            // Car_obj->compute_dx(Delta_x_vec);
            Load_module_obj->update_force(t, force_vector, centripetal_force);
            Load_module_obj->update_torque(t, torque, force_vector);
            // convert centrifugal force to centripetal (only for x, y direction)
            mkl<T>::scal(centripetal_force_dimensions - 1, -1.0, centripetal_force, 1);
            if (iter == 100)
                MathLibrary::write_vector(centripetal_force, centripetal_force_dimensions);
            global_frame_solver(t);

            // translate 27 force vector + 3 torques into 11DOF
            Car_obj->construct_11DOF_vector(force_vector, new_torque, force_vector_11dof);

            TwoTrackModel_obj->update_step(force_vector_11dof,
                                           Car_obj->currentDisplacementTwoTrackModel);
            Car_obj->updateLengthsTwoTrackModel();
            solution_vect = u_sol + iter * (Constants::VEC_DIM * Constants::DIM);

            // only call this function at every checkpoint
            Car_obj->combineEulerianLagrangianVectors(Car_obj->currentPositionLagrangian,
                                                      Car_obj->currentDisplacementTwoTrackModel,
                                                      solution_vect);

            t += h_;
            iter++;
        }

        mkl<T>::copy((Constants::VEC_DIM * Constants::DIM),
                     u_sol + (iter - 1) * (Constants::VEC_DIM * Constants::DIM), 1, sol_vect, 1);
        Car_obj->combine_results();

        MKL_free(time_vec);

        MKL_free(force_vector);
        MKL_free(centripetal_force);
        MKL_free(new_centripetal_force);
        MKL_free(Delta_x_vec);

        delete[] torque;
        delete[] new_torque;
    }

    /**
     * Executes the time iteration of the ALE solvers, switches from global position update to
     * solving of the linear 11DOF system
     */
    void solve(T* sol_vect)
    {
        int sol_size = (floor(tend_ / h_) + 1);
        u_sol = (T*)mkl_calloc(sol_size * (Constants::VEC_DIM * Constants::DIM), sizeof(T),
                               Constants::ALIGNMENT);
        solve(sol_vect, u_sol);
        MKL_free(u_sol);
    }

    /**
     * Adds the contribution of the wheels and tyres to the inertia moment of the car
     */
    void calculate_global_inertia_Z()
    {
        // get the global inertia actiing in Z direction
        global_inertia_Z = Car_obj->momentOfInertia[8];
        global_inertia_Z +=
            (Car_obj->massComponents[1] + Car_obj->massComponents[2]) *
            (Car_obj->l_lat[0] * Car_obj->l_lat[0] + Car_obj->l_long[0] * Car_obj->l_long[0]);
        global_inertia_Z +=
            (Car_obj->massComponents[3] + Car_obj->massComponents[4]) *
            (Car_obj->l_lat[1] * Car_obj->l_lat[1] + Car_obj->l_long[1] * Car_obj->l_long[1]);
        global_inertia_Z +=
            (Car_obj->massComponents[5] + Car_obj->massComponents[6]) *
            (Car_obj->l_lat[2] * Car_obj->l_lat[2] + Car_obj->l_long[2] * Car_obj->l_long[2]);
        global_inertia_Z +=
            (Car_obj->massComponents[7] + Car_obj->massComponents[8]) *
            (Car_obj->l_lat[3] * Car_obj->l_lat[3] + Car_obj->l_long[3] * Car_obj->l_long[3]);
    }

    /**
     * Prints all positions and angles in the car object
     */
    void print_final_results()
    {
        T* sln = Car_obj->Position_vec;
        std::cout << "ALE: orientation angles=\n\t[" << Car_obj->angle_CG[0] << "\n\t "
                  << Car_obj->angle_CG[1] << "\n\t " << Car_obj->angle_CG[2] << "]" << std::endl;
        std::cout << "ALE: car body position pc=\n\t[" << sln[0] << "\n\t " << sln[1] << "\n\t "
                  << sln[2] << "]" << std::endl;
        std::cout << "ALE: rear-right wheel position pw1=\n\t[" << sln[21] << "\n\t " << sln[22]
                  << "\n\t " << sln[23] << "]" << std::endl;
        std::cout << "ALE: rear-left wheel position pw2=\n\t[" << sln[15] << "\n\t " << sln[16]
                  << "\n\t " << sln[17] << "]" << std::endl;
        std::cout << "ALE: front-left wheel position pw3=\n\t[" << sln[3] << "\n\t " << sln[4]
                  << "\n\t " << sln[5] << "]" << std::endl;
        std::cout << "ALE: front-right wheel position pw4=\n\t[" << sln[9] << "\n\t " << sln[10]
                  << "\n\t " << sln[11] << "]" << std::endl;
        std::cout << "ALE: rear-right tyre position pt1=\n\t[" << sln[24] << "\n\t " << sln[25]
                  << "\n\t " << sln[26] << "]" << std::endl;
        std::cout << "ALE: rear-left tyre position pt2=\n\t[" << sln[18] << "\n\t " << sln[19]
                  << "\n\t " << sln[20] << "]" << std::endl;
        std::cout << "ALE: front-left tyre position pt3=\n\t[" << sln[6] << "\n\t " << sln[7]
                  << "\n\t " << sln[8] << "]" << std::endl;
        std::cout << "ALE: front-right tyre position pt4=\n\t[" << sln[12] << "\n\t " << sln[13]
                  << "\n\t " << sln[14] << "]" << std::endl;
    }
};

}  // namespace EVAA
