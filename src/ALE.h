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
    Car<T>* _carObj;                 // suppose Interpolation in the Car
    LoadModule<T>* _loadModuleObj;  // needs Profile and Car
    TwoTrackModelParent<T>* _twoTrackModelObj;

    // simulation parameters
    T _tend;
    T _h;
    size_t _solutionVectorSize;

    // time and solution vectors
    T* _timeArray;
    T* _solutionVector;

    // global centripetal forces on the whole car [XYZ]
    T* _lagrangianForceVector;
    T* _newLagrangianForceVector;

    // global torque on the car [XYZ]
    T* _lagrangianTorque;
    T* _newLagrangianTorque;

    // quantities for the whole car
    T _momentOfInertiaZ;

public:
    /**
     * Constructor
     */
    ALE(Car<T>* carObjVal, LoadModule<T>* loadModuleVal, TwoTrackModelParent<T>* twoTrackModelVal) : //
        _carObj(carObjVal), _loadModuleObj(loadModuleVal), _twoTrackModelObj(twoTrackModelVal) {        
        _h = MetaDatabase<T>::getDatabase().getTimeStepSize();
        _tend = MetaDatabase<T>::getDatabase().getNumberOfTimeIterations() * _h;

        _solutionVectorSize = (floor(_tend / _h) + 1);
        _solutionVector = Math::malloc<T>(_solutionVectorSize * (Constants::VEC_DIM * Constants::DIM));
    }

    ~ALE() { Math::free<T>(_solutionVector); }

    /**
     * Applies the Verlet_Stoermer algorithm to update the global XY position of
     * the car and its Z orientation Store the global coordinates in the
     * VelocityXY and PositionXY from the car object \param t current simulation
     * time
     */
    void LagrangianUpdate(T& t) {
        // 2. Update global X,Y positions of the car
        Math::Solvers<T, ALE>::StoermerVerletPosition(_carObj->_currentPositionLagrangian[0], _carObj->getCurrentVelocityLagrangian()[0], _lagrangianForceVector[0], _h, _carObj->getMassCarFull());
        Math::Solvers<T, ALE>::StoermerVerletPosition(_carObj->_currentPositionLagrangian[1], _carObj->getCurrentVelocityLagrangian()[1], _lagrangianForceVector[1], _h, _carObj->getMassCarFull());

        // 4. Update Z-rotation
        Math::Solvers<T, ALE>::StoermerVerletPosition(_carObj->_currentAngleLagrangian, _carObj->getCurrentAngularVelocityLagrangian(), *_lagrangianTorque, _h, _momentOfInertiaZ);

        // get forces
        _loadModuleObj->GetLagrangianForce(t, _newLagrangianForceVector);
        _loadModuleObj->GetTorqueLagrange(t, _newLagrangianTorque);

        // 1. Update global X,Y velocities
        Math::Solvers<T, ALE>::StoermerVerletVelocity(_carObj->_currentVelocityLagrangian[0], _lagrangianForceVector[0], _newLagrangianForceVector[0], _h, _carObj->getMassCarFull());
        Math::Solvers<T, ALE>::StoermerVerletVelocity(_carObj->_currentVelocityLagrangian[1], _lagrangianForceVector[1], _newLagrangianForceVector[1], _h, _carObj->getMassCarFull());

        // 3. Update Z-angular velocities
        Math::Solvers<T, ALE>::StoermerVerletVelocity(_carObj->_currentAngularVelocityLagrangian, *_lagrangianTorque, *_newLagrangianTorque, _h, _momentOfInertiaZ);
		/*
         * Idea!! What if we do it at the end, since the displacement is a vector and by triangle
         * rule sum of all should add up force is computed using on
         */
        _carObj->ApplyLagrangeChange(); // TODO Rethink how it can be combined with 11 dof computation

        // update forces and _lagrangianTorque
        _lagrangianForceVector[0] = _newLagrangianForceVector[0];
        _lagrangianForceVector[1] = _newLagrangianForceVector[1];

        *_lagrangianTorque = *_newLagrangianTorque;  // z - component
    }

    void Solve(T* sol_vect, T* u_sol_param) {
        // initialize solution vector
        const int lagrangianForceDimension = Constants::DIM - 1;

        // allocate memory

        _timeArray = Math::calloc<T>(_solutionVectorSize);
        _lagrangianForceVector = Math::calloc<T>(lagrangianForceDimension);

        // this was 2 dimensional allocation and update force updates 3 dimension on this
        _newLagrangianForceVector = Math::calloc<T>(lagrangianForceDimension);

        //springElongation = Math::calloc<T>(2 * Constants::NUM_LEGS);

        _lagrangianTorque = new T;
        _newLagrangianTorque = new T;

        // calculate characteristics of the whole car
        CalculateGlobalMomentofInertiaZ();

        // start time iteration
        T t = _h;

#ifdef WRITECSV
        IO::MyFile<T> solutionCSV("C:\\software\\repos\\EVAA\\output\\aleSolution.txt");
        IO::MyFile<T> parametersCSV("C:\\software\\repos\\EVAA\\output\\simulationParameters.txt");
        parametersCSV.writeParameters();
        T* angleVecCSV = Math::malloc<T>(_solutionVectorSize * Constants::DIM);
        T* posVecCSV = Math::malloc<T>(_solutionVectorSize * Constants::DIM * Constants::VEC_DIM);
        T* velVecCSV = Math::malloc<T>(_solutionVectorSize * (Constants::DIM - 1) * Constants::VEC_DIM);
#endif // WRITECSV

        T* solution_vect;
        int iter = 1;
        // time iteration
        double eps = _h / 100;
        while (std::abs(t - (_tend + _h)) > eps) {
            // This has to be done at each time step
            //
            // update force vector
			_loadModuleObj->GetLagrangianForce(t, _lagrangianForceVector);
			_loadModuleObj->GetTorqueLagrange(t, _lagrangianTorque);
			//if (iter == 1000) IO::writeVector(_lagrangianForceVector, lagrangianForceDimension); 
		    LagrangianUpdate(t);

            _twoTrackModelObj->UpdateStep(t, _carObj->_currentDisplacementTwoTrackModel);
            _carObj->UpdateLengthsTwoTrackModel();
            _carObj->ApplyLagrangeChange();
            solution_vect = u_sol_param + iter * (Constants::VEC_DIM * Constants::DIM);

            // only call this function at every checkpoint
            _carObj->CombineEulerianLagrangianVectors(solution_vect);

#ifdef WRITECSV
            _carObj->CombineResults();            
            Math::copy(Constants::VEC_DIM * Constants::DIM, _carObj->getPositionVector(), 1, posVecCSV + iter * Constants::DIM * Constants::VEC_DIM, 1);
            Math::copy(Constants::DIM, _carObj->getAngleCG(), 1, angleVecCSV + iter * Constants::DIM, 1);
            Math::copy(Constants::DIM, _carObj->_currentVelocityLagrangian, 1, velVecCSV + iter * (Constants::DIM-1) * Constants::VEC_DIM, 1);
#endif // WRITECSV
            t += _h;
            iter++;
        }

#ifdef WRITECSV
        solutionCSV.writeSolutionMatrix(posVecCSV, velVecCSV, angleVecCSV, _solutionVectorSize);
        Math::free(angleVecCSV);
        Math::free(velVecCSV);
        Math::free(posVecCSV);
#endif  // WRITECSV


        Math::copy<T>(Constants::VEC_DIM * Constants::DIM, u_sol_param + (iter - 1) * (Constants::VEC_DIM * Constants::DIM), 1, sol_vect, 1);
        _carObj->CombineResults();

        Math::free<T>(_timeArray);

        Math::free<T>(_lagrangianForceVector);
        Math::free<T>(_newLagrangianForceVector);
        delete _lagrangianTorque;
		delete _newLagrangianTorque;
    }

    /**
     * Executes the time iteration of the ALE solvers, switches from global
     * position update to solving of the linear 11DOF system
     */
    void Solve(T* sol_vect) { Solve(sol_vect, _solutionVector); }

    /**
     * Adds the contribution of the wheels and tyres to the inertia moment of
     * the car
     */
    void CalculateGlobalMomentofInertiaZ() {
        // get the global inertia actiing in Z direction
        _momentOfInertiaZ = _carObj->getMomentOfInertia()[8];
        _momentOfInertiaZ +=
            (_carObj->getMassComponents()[1] + _carObj->getMassComponents()[2]) *
            (_carObj->_lenLat[0] * _carObj->_lenLat[0] + _carObj->_lenLong[0] * _carObj->_lenLong[0]);
        _momentOfInertiaZ +=
            (_carObj->getMassComponents()[3] + _carObj->getMassComponents()[4]) *
            (_carObj->_lenLat[1] * _carObj->_lenLat[1] + _carObj->_lenLong[1] * _carObj->_lenLong[1]);
        _momentOfInertiaZ +=
            (_carObj->getMassComponents()[5] + _carObj->getMassComponents()[6]) *
            (_carObj->_lenLat[2] * _carObj->_lenLat[2] + _carObj->_lenLong[2] * _carObj->_lenLong[2]);
        _momentOfInertiaZ +=
            (_carObj->getMassComponents()[7] + _carObj->getMassComponents()[8]) *
            (_carObj->_lenLat[3] * _carObj->_lenLat[3] + _carObj->_lenLong[3] * _carObj->_lenLong[3]);
    }

    /**
     * Prints all positions and angles in the car object
     */
    void PrintFinalResults() {
        const T* sln = _carObj->getPositionVector();
        std::cout << "ALE: orientation angles=\n\t[" << _carObj->getAngleCG()[0] << "\n\t " << _carObj->getAngleCG()[1] << "\n\t " << _carObj->getAngleCG()[2] << "]" << std::endl;
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
