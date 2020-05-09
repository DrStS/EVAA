// TODO: Copyright header

#pragma once
#include "11DOF.h"
#include "Constants.h"
#include "EVAAComputeEngine.h"
#include "LoadModule.h"
#include "MathLibrary.h"
#include "RoadProfile.h"

#ifdef USE_HDF5
#include "OutputHDF5.h"
#endif

namespace EVAA {

/** Implements the ALE method to extend the linear 11DOF system */
template <class T>
class ALE {
private:
    // necessary class objects
    Car<T>* _carObj;                // suppose Interpolation in the Car
    LoadModule<T>* _loadModuleObj;  // needs Profile and Car
    TwoTrackModelParent<T>* _twoTrackModelObj;

    // simulation parameters
    T _tend;
    T _h;
    size_t _solutionVectorSize;

    // time and solution vectors
    T* _timeArray = nullptr;
    T* _fullSolution = nullptr;

    // global centripetal forces on the whole car [XYZ]
    T* _lagrangianForceVector = nullptr;
    T* _newLagrangianForceVector = nullptr;

    // global torque on the car [XYZ]
    T _lagrangianTorque;
    T _newLagrangianTorque;

    // quantities for the whole car
    T _momentOfInertiaZ;

#ifdef WRITECSV
    IO::MyFile<T>* solutionCSV = nullptr;
    T* angleVecCSV;
    T* posVecCSV;
    T* velVecCSV;
#endif  // WRITECSV

#ifdef USE_HDF5
    HDF5::OutputHDF5<T>* _checkpointsALE;
    std::string _groupNameCheckpoints;  // basic name for a checkpoint group
    std::string _filePath;
    std::string _fileName;
#endif

public:
    /**
     * Constructor (without considering writing in Output HDF5)
     */
    ALE(Car<T>* carObjVal, LoadModule<T>* loadModuleVal, TwoTrackModelParent<T>* twoTrackModelVal) :
        //
        _carObj(carObjVal),
        _loadModuleObj(loadModuleVal),
        _twoTrackModelObj(twoTrackModelVal) {
        _h = MetaDatabase<T>::getDatabase().getTimeStepSize();
        _tend = MetaDatabase<T>::getDatabase().getNumberOfTimeIterations() * _h;

        _solutionVectorSize = (floor(_tend / _h) + 1);
        _fullSolution = Math::malloc<T>(_solutionVectorSize * (Constants::VEC_DIM * Constants::DIM));

        // initialize solution vector
        _timeArray = Math::calloc<T>(_solutionVectorSize);                              // TODO Allocation in constructor
        _lagrangianForceVector = Math::calloc<T>(Constants::lagrangianForceDimension);  // TODO Allocation in constructor

        // this was 2 dimensional allocation and update force updates 3 dimension on this
        _newLagrangianForceVector = Math::calloc<T>(Constants::lagrangianForceDimension);

#ifdef WRITECSV
        solutionCSV = new IO::MyFile<T>("C:\\software\\repos\\EVAA\\output\\aleSolution.txt");
        IO::MyFile<T> parametersCSV("C:\\software\\repos\\EVAA\\output\\simulationParameters.txt");
        parametersCSV.writeParameters();
        angleVecCSV = Math::malloc<T>(_solutionVectorSize * Constants::DIM);
        posVecCSV = Math::malloc<T>(_solutionVectorSize * Constants::DIM * Constants::VEC_DIM);
        velVecCSV = Math::malloc<T>(_solutionVectorSize * (Constants::DIM - 1) * Constants::VEC_DIM);
#endif  // WRITECSV
    }

    /**
     * Constructor (considering writing in Output HDF5)
     */
#ifdef USE_HDF5
    ALE(Car<T>* carObjVal, LoadModule<T>* loadModuleVal, TwoTrackModelParent<T>* twoTrackModelVal, std::string filePath, std::string fileName) : ALE(carObjVal, loadModuleVal, twoTrackModelVal) {
        _checkpointsALE = new HDF5::OutputHDF5<T>(filePath, fileName);
        _groupNameCheckpoints = "ALE Checkpoint t = ";
        _checkpointsALE->CreateContainer(true, _groupNameCheckpoints);
        _checkpointsALE->CloseContainer();
    }
#endif

    ~ALE() {
        Math::free<T>(_fullSolution);
        Math::free<T>(_timeArray);
        Math::free<T>(_lagrangianForceVector);
        Math::free<T>(_newLagrangianForceVector);
        // TODO : Move this to destructor, it is not part of the solver
#ifdef WRITECSV
        Math::free(angleVecCSV);
        Math::free(velVecCSV);
        Math::free(posVecCSV);
#endif  // WRITECSV
#ifdef USE_HDF5
        delete _checkpointsALE;
#endif
    }

    /**
     * Applies the Verlet_Stoermer algorithm to update the global XY position of
     * the car and its Z orientation Store the global coordinates in the
     * VelocityXY and PositionXY from the car object \param t current simulation
     * time
     */
    void LagrangianUpdate(const size_t iter) {
        // 1. Update global X,Y positions of the car
        Math::Solvers<T, ALE>::StoermerVerletPosition(_carObj->_currentPositionLagrangian[0], _carObj->getCurrentVelocityLagrangian()[0], _lagrangianForceVector[0], _h, _carObj->getMassCarFull());
        Math::Solvers<T, ALE>::StoermerVerletPosition(_carObj->_currentPositionLagrangian[1], _carObj->getCurrentVelocityLagrangian()[1], _lagrangianForceVector[1], _h, _carObj->getMassCarFull());

        // 2. Update Z-rotation
        Math::Solvers<T, ALE>::StoermerVerletPosition(_carObj->_currentAngleLagrangian, _carObj->getCurrentAngularVelocityLagrangian(), _lagrangianTorque, _h, _momentOfInertiaZ);

        // 3. get forces
        _loadModuleObj->GetLagrangianForce(iter, _newLagrangianForceVector);
        _loadModuleObj->GetTorqueLagrange(iter, &_newLagrangianTorque);

        // 4. Update global X,Y velocities
        Math::Solvers<T, ALE>::StoermerVerletVelocity(_carObj->_currentVelocityLagrangian[0], _lagrangianForceVector[0], _newLagrangianForceVector[0], _h, _carObj->getMassCarFull());
        Math::Solvers<T, ALE>::StoermerVerletVelocity(_carObj->_currentVelocityLagrangian[1], _lagrangianForceVector[1], _newLagrangianForceVector[1], _h, _carObj->getMassCarFull());

        // 5. Update Z-angular velocities
        Math::Solvers<T, ALE>::StoermerVerletVelocity(_carObj->_currentAngularVelocityLagrangian, _lagrangianTorque, _newLagrangianTorque, _h, _momentOfInertiaZ);

        // 6. update forces and _lagrangianTorque
        _lagrangianForceVector[0] = _newLagrangianForceVector[0];
        _lagrangianForceVector[1] = _newLagrangianForceVector[1];

        _lagrangianTorque = _newLagrangianTorque;  // z - component
    }

    /**
     * Executes the time iteration of the ALE solvers, switches from global
     * position update to solving of the linear 11DOF system
     */
    void Solve(T* finalSolution) {
        // calculate characteristics of the whole car
        CalculateGlobalMomentofInertiaZ();

        // time iteration
        size_t iter;
        T t;
        T eps = _h / 100;
        for (iter = 1, t = _h; std::abs(t - (_tend + _h)) > eps; iter++, t += _h) {
            // TODO : What is the stoping criterion? We have the number of iterations apriori!

            // This has to be done at each time step
            //
            // update force vector
            _loadModuleObj->GetLagrangianForce(iter, _lagrangianForceVector);
            _loadModuleObj->GetTorqueLagrange(iter, &_lagrangianTorque);
            // if (iter == 1000) IO::writeVector(_lagrangianForceVector,
            // Constants::lagrangianForceDimension);
            LagrangianUpdate(t);

            _twoTrackModelObj->UpdateStep(iter, _carObj->_currentDisplacementTwoTrackModel);
            _carObj->UpdateLengthsTwoTrackModel();
            _carObj->ApplyLagrangeChange();

#ifdef WRITECSV
            _carObj->CombineResults();
            Math::copy(Constants::VEC_DIM * Constants::DIM, _carObj->getPositionVector(), 1, posVecCSV + iter * Constants::DIM * Constants::VEC_DIM, 1);
            Math::copy(Constants::DIM, _carObj->getAngleCG(), 1, angleVecCSV + iter * Constants::DIM, 1);
            Math::copy(Constants::DIM, _carObj->_currentVelocityLagrangian, 1, velVecCSV + iter * (Constants::DIM - 1) * Constants::VEC_DIM, 1);
#endif  // WRITECSV

#ifdef USE_HDF5
            T* solution_vect = _fullSolution + iter * (Constants::VEC_DIM * Constants::DIM);
            _carObj->CombineEulerianLagrangianVectors(solution_vect);
            // Call this only at checkpoints
            _checkpointsALE->CreateContainer(false, _groupNameCheckpoints + std::to_string(t));
            // Write whatever vectors / matrices
            _checkpointsALE->WriteVector("Solution Vector; iter " + std::to_string(iter),     //
                                         solution_vect, Constants::VEC_DIM * Constants::DIM,  //
                                         EVAA::HDF5FileHandle::GROUP);
            // TODO Write specialised vectors
            // TODO Write vectors for points of interest
            _checkpointsALE->CloseContainer();
#endif  // USE_HDF5
        }

#ifdef WRITECSV
        solutionCSV->writeSolutionMatrix(posVecCSV, velVecCSV, angleVecCSV, _solutionVectorSize);
#endif  // WRITECSV

        // TODO How do we put the initial solution in _fullSolution?
        Math::copy<T>(Constants::VEC_DIM * Constants::DIM, _fullSolution + (iter - 1) * (Constants::VEC_DIM * Constants::DIM), 1, finalSolution, 1);
        _carObj->CombineResults();
    }

    /**
     * Adds the contribution of the wheels and tyres to the inertia moment of
     * the car
     */
    void CalculateGlobalMomentofInertiaZ() {
        // get the global inertia actiing in Z direction
        _momentOfInertiaZ = _carObj->getMomentOfInertia()[8];
        _momentOfInertiaZ += (_carObj->getMassComponents()[1] + _carObj->getMassComponents()[2]) * (_carObj->_lenLat[0] * _carObj->_lenLat[0] + _carObj->_lenLong[0] * _carObj->_lenLong[0]);
        _momentOfInertiaZ += (_carObj->getMassComponents()[3] + _carObj->getMassComponents()[4]) * (_carObj->_lenLat[1] * _carObj->_lenLat[1] + _carObj->_lenLong[1] * _carObj->_lenLong[1]);
        _momentOfInertiaZ += (_carObj->getMassComponents()[5] + _carObj->getMassComponents()[6]) * (_carObj->_lenLat[2] * _carObj->_lenLat[2] + _carObj->_lenLong[2] * _carObj->_lenLong[2]);
        _momentOfInertiaZ += (_carObj->getMassComponents()[7] + _carObj->getMassComponents()[8]) * (_carObj->_lenLat[3] * _carObj->_lenLat[3] + _carObj->_lenLong[3] * _carObj->_lenLong[3]);
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

    void WriteBulkResults(std::string filePath = "", std::string fileName = "ALE_full_solution.hdf5", std::string datasetName = "ALE Final Solution") {
#ifdef USE_HDF5
        HDF5::OutputHDF5<Constants::floatEVAA> fullALE(filePath, fileName);
        fullALE.CreateContainer(true);
        fullALE.WriteMatrix(datasetName, _fullSolution, _solutionVectorSize, Constants::VEC_DIM * Constants::DIM);
        fullALE.CloseContainer();
#endif  // USE_HDF5
    }
};

}  // namespace EVAA
