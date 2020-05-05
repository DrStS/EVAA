// TODO: Copyright header

#pragma once
#include "11DOF.h"
#include "Constants.h"
#include "EVAAComputeEngine.h"
#include "LoadModule.h"
#include "MathLibrary.h"
#include "RoadProfile.h"

namespace EVAA {

/** Implements the ALE method to extend the linear 11DOF system */
template <class T>
class ALE {
private:
    // necessary class objects
    Car<T>* carObj;                 // suppose Interpolation in the Car
    LoadModule<T>* loadModuleObj;  // needs Profile and Car
    TwoTrackModelParent<T>* twoTrackModelObj;

    // simulation parameters
    T tend_;
    T h_;

    // time and solution vectors
    T* timeArray;
    T* solutionVector;

    // needed to solve the 11DOF system
    T* forceTwoTrackModel;

    // Contains the forces in the ordering [GC:XYZ,W1:XYZ,T1:XYZ, ...]
    T* combinedForceVector;

    // global centripetal forces on the whole car [XYZ]
    T* lagrangianForceVector;
    T* newLagrangianForceVector;

    // global torque on the car [XYZ]
    T* lagrangianTorque;
    T* newLagrangianTorque;

    // quantities for the whole car
    T momentOfInertiaZ;

public:
    /**
     * Constructor
     */
    ALE(Car<T>* Car_obj_val, LoadModule<T>* Load_module_val, TwoTrackModelParent<T>* TwoTrackModel_val) {
        carObj = Car_obj_val;
        loadModuleObj = Load_module_val;
        twoTrackModelObj = TwoTrackModel_val;

        h_ = MetaDataBase<T>::getDataBase().getTimeStepSize();
        tend_ = MetaDataBase<T>::getDataBase().getNumberOfTimeIterations() * h_;

        size_t solutionVectorSize = (floor(tend_ / h_) + 1);
        solutionVector = Math::malloc<T>(solutionVectorSize * (Constants::VEC_DIM * Constants::DIM));
    }

    ~ALE() { Math::free<T>(solutionVector); }

    /**
     * Applies the Verlet_Stoermer algorithm to update the global XY position of
     * the car and its Z orientation Store the global coordinates in the
     * VelocityXY and PositionXY from the car object \param t current simulation
     * time
     */
    void LagrangianUpdate(T& t) {
        // 2. Update global X,Y positions of the car
        Math::Solvers<T, ALE>::Stoermer_Verlet_Position(carObj->currentPositionLagrangian[0], carObj->currentVelocityLagrangian[0], lagrangianForceVector[0], h_, carObj->massFullCar);
        Math::Solvers<T, ALE>::Stoermer_Verlet_Position(carObj->currentPositionLagrangian[1], carObj->currentVelocityLagrangian[1], lagrangianForceVector[1], h_, carObj->massFullCar);

        // 4. Update Z-rotation
        Math::Solvers<T, ALE>::Stoermer_Verlet_Position(carObj->currentAngleLagrangian, carObj->currentAngularVelocityLagrangian, *lagrangianTorque, h_, momentOfInertiaZ);

        // get forces
        loadModuleObj->GetLagrangianForce(t, newLagrangianForceVector);
        loadModuleObj->GetTorqueLagrange(t, newLagrangianTorque);

        // 1. Update global X,Y velocities
        Math::Solvers<T, ALE>::Stoermer_Verlet_Velocity(carObj->currentVelocityLagrangian[0], lagrangianForceVector[0], newLagrangianForceVector[0], h_, carObj->massFullCar);
        Math::Solvers<T, ALE>::Stoermer_Verlet_Velocity(carObj->currentVelocityLagrangian[1], lagrangianForceVector[1], newLagrangianForceVector[1], h_, carObj->massFullCar);

        // 3. Update Z-angular velocities
        Math::Solvers<T, ALE>::Stoermer_Verlet_Velocity(carObj->currentAngularVelocityLagrangian, *lagrangianTorque, *newLagrangianTorque, h_, momentOfInertiaZ);
		/*
         * Idea!! What if we do it at the end, since the displacement is a vector and by triangle
         * rule sum of all should add up force is computed using on
         */
        carObj->ApplyLagrangeChange(); // Rethink how it can be combined with 11 dof computation

        // update forces and lagrangianTorque
        lagrangianForceVector[0] = newLagrangianForceVector[0];
        lagrangianForceVector[1] = newLagrangianForceVector[1];

        *lagrangianTorque = *newLagrangianTorque;  // z - component
    }

    void solve(T* sol_vect, T* u_sol_param) {
        // initialize solution vector
        int solutionVectorSize = (floor(tend_ / h_) + 1);
        int lagrangianForceDimension = Constants::DIM - 1;

        // allocate memory

        timeArray = Math::calloc<T>(solutionVectorSize);
        lagrangianForceVector = Math::calloc<T>(lagrangianForceDimension);

        // this was 2 dimensional allocation and update force updates 3 dimension on this
        newLagrangianForceVector = Math::calloc<T>(lagrangianForceDimension);

        //springElongation = Math::calloc<T>(2 * Constants::NUM_LEGS);

        lagrangianTorque = new T;
        newLagrangianTorque = new T;

        // calculate characteristics of the whole car
        CalculateGlobalMomentofInertiaZ();

        // start time iteration
        T t = h_;

#ifdef WRITECSV
        IO::MyFile<T> solutionCSV("C:\\software\\repos\\EVAA\\output\\aleSolution.txt");
        IO::MyFile<T> parametersCSV("C:\\software\\repos\\EVAA\\output\\simulationParameters.txt");
        parametersCSV.writeParameters();
        T* angleVecCSV = Math::malloc<T>(solutionVectorSize * Constants::DIM);
        T* posVecCSV = Math::malloc<T>(solutionVectorSize * Constants::DIM * Constants::VEC_DIM);
        T* velVecCSV = Math::malloc<T>(solutionVectorSize * (Constants::DIM - 1) * Constants::VEC_DIM);
#endif // WRITECSV

        T* solution_vect;
        int iter = 1;
        // time iteration
        double eps = h_ / 100;
        while (std::abs(t - (tend_ + h_)) > eps) {
            // This has to be done at each time step
            //
            // update force vector
			loadModuleObj->GetLagrangianForce(t, lagrangianForceVector);
			loadModuleObj->GetTorqueLagrange(t, lagrangianTorque);
			if (iter == 1000) IO::writeVector(lagrangianForceVector, lagrangianForceDimension); 
		    LagrangianUpdate(t);

            twoTrackModelObj->update_step(t, carObj->currentDisplacementTwoTrackModel);
            carObj->updateLengthsTwoTrackModel();
            solution_vect = u_sol_param + iter * (Constants::VEC_DIM * Constants::DIM);

            // only call this function at every checkpoint
            carObj->combineEulerianLagrangianVectors(solution_vect);

#ifdef WRITECSV
            carObj->combine_results();            
            Math::copy(Constants::VEC_DIM * Constants::DIM, solution_vect, 1, posVecCSV + iter * Constants::DIM * Constants::VEC_DIM, 1);
            Math::copy(Constants::DIM, carObj->angle_CG, 1, angleVecCSV + iter * Constants::DIM, 1);
            Math::copy(Constants::DIM, carObj->currentVelocityLagrangian, 1, velVecCSV + iter * (Constants::DIM-1) * Constants::VEC_DIM, 1);
#endif // WRITECSV

            t += h_;
            iter++;
        }

#ifdef WRITECSV
        solutionCSV.writeSolutionMatrix(posVecCSV, velVecCSV, angleVecCSV, solutionVectorSize);
        Math::free(angleVecCSV);
        Math::free(velVecCSV);
        Math::free(posVecCSV);
#endif  // WRITECSV


        Math::copy<T>(Constants::VEC_DIM * Constants::DIM, u_sol_param + (iter - 1) * (Constants::VEC_DIM * Constants::DIM), 1, sol_vect, 1);
        carObj->combine_results();

        Math::free<T>(timeArray);

        Math::free<T>(lagrangianForceVector);
        Math::free<T>(newLagrangianForceVector);
        delete lagrangianTorque;
		delete newLagrangianTorque;
    }

    /**
     * Executes the time iteration of the ALE solvers, switches from global
     * position update to solving of the linear 11DOF system
     */
    void solve(T* sol_vect) { solve(sol_vect, solutionVector); }

    /**
     * Adds the contribution of the wheels and tyres to the inertia moment of
     * the car
     */
    void CalculateGlobalMomentofInertiaZ() {
        // get the global inertia actiing in Z direction
        momentOfInertiaZ = carObj->momentOfInertia[8];
        momentOfInertiaZ +=
            (carObj->massComponents[1] + carObj->massComponents[2]) *
            (carObj->l_lat[0] * carObj->l_lat[0] + carObj->l_long[0] * carObj->l_long[0]);
        momentOfInertiaZ +=
            (carObj->massComponents[3] + carObj->massComponents[4]) *
            (carObj->l_lat[1] * carObj->l_lat[1] + carObj->l_long[1] * carObj->l_long[1]);
        momentOfInertiaZ +=
            (carObj->massComponents[5] + carObj->massComponents[6]) *
            (carObj->l_lat[2] * carObj->l_lat[2] + carObj->l_long[2] * carObj->l_long[2]);
        momentOfInertiaZ +=
            (carObj->massComponents[7] + carObj->massComponents[8]) *
            (carObj->l_lat[3] * carObj->l_lat[3] + carObj->l_long[3] * carObj->l_long[3]);
    }

    /**
     * Prints all positions and angles in the car object
     */
    void print_final_results() {
        T* sln = carObj->Position_vec;
        std::cout << "ALE: orientation angles=\n\t[" << carObj->angle_CG[0] << "\n\t " << carObj->angle_CG[1] << "\n\t " << carObj->angle_CG[2] << "]" << std::endl;
        std::cout << "ALE: car body position pc=\n\t[" << sln[0] << "\n\t " << sln[1] << "\n\t " << sln[2] << "]" << std::endl;
        std::cout << "ALE: rear-right wheel position pw1=\n\t[" << sln[21] << "\n\t " << sln[22] << "\n\t " << sln[23] << "]" << std::endl;
        std::cout << "ALE: rear-left wheel position pw2=\n\t[" << sln[15] << "\n\t " << sln[16] << "\n\t " << sln[17] << "]" << std::endl;
        std::cout << "ALE: front-left wheel position pw3=\n\t[" << sln[3] << "\n\t " << sln[4] << "\n\t " << sln[5] << "]" << std::endl;
        std::cout << "ALE: front-right wheel position pw4=\n\t[" << sln[9] << "\n\t " << sln[10] << "\n\t " << sln[11] << "]" << std::endl;
        std::cout << "ALE: rear-right tyre position pt1=\n\t[" << sln[24] << "\n\t " << sln[25] << "\n\t " << sln[26] << "]" << std::endl;
        std::cout << "ALE: rear-left tyre position pt2=\n\t[" << sln[18] << "\n\t " << sln[19] << "\n\t " << sln[20] << "]" << std::endl;
        std::cout << "ALE: front-left tyre position pt3=\n\t[" << sln[6] << "\n\t " << sln[7] << "\n\t " << sln[8] << "]" << std::endl;
        std::cout << "ALE: front-right tyre position pt4=\n\t[" << sln[12] << "\n\t " << sln[13] << "\n\t " << sln[14] << "]" << std::endl;
    }
};

}  // namespace EVAA
