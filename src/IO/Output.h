// TODO: Copyright header

#pragma once

#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

#include "MetaDatabase.h"

namespace EVAA {

namespace IO {

/**
 * Puts a vector representation into standard output.
 * \param out The output stream.
 * \param vect The data
 * \param count The numver of elements in the vector.
 */
template <typename T>
void writeVector(const T* vector, const size_t count) {
    writeVector(std::cout, vector, count);
}

/**
 * Puts a vector representation into an output stream.
 * \param out The output stream.
 * \param vect The data
 * \param count The numver of elements in the vector.
 */
template <typename T>
void writeVector(std::ostream& out, const T* vect, const size_t count) {
    out << "[ ";
    for (auto i = 0; i < count; ++i) {
        out.precision(15);
        out << std::scientific << vect[i] << std::endl;
    }
    out << "]" << std::endl;
}

/**
 * Puts a matrix representation into an output stream.
 * \param out The output stream.
 * \param matrix The matrix liniarized into an array.
 * \param rows The number of rows in the matrix.
 * \param cols The number of columns in the matrix.
 */
template <typename T>
void writeMatrix(std::ostream& out, const T* matrix, const size_t rows, const size_t cols) {
    out << "[" << std::endl;
    for (auto i = 0; i < count; ++i) {
        out << "  ";
        writeVector(out, matrix + cols * i, cols);
    }
    out << std::endl;
}

/**
 * Writes a matrix representation into a file.
 * \param fname The file name.
 * \param mat The matrix liniarized into an array.
 * \param rows The number of rows in the matrix.
 * \param cols The number of columns in the matrix.
 */
template <typename T>
void writeMatrix(std::string& fname, const T* mat, const size_t rows, const size_t cols) {
    std::ofstream myfile(fname);
    writeMatrix(myfile, T, rows, cols);
    myfile.close();
};

/**
 * \brief function to print the values for the lookup table and the interpolated
 * values for debugging
 */

template <typename T>
void writeLookUpGridPlusInterpolateValues(T* axis /**< [in] pointer to axis */, T* grid /**< [in] pointer to grid */, int n /**< [in] number of elements in arrays */, T* interpolationPoints /**< [in] pointer to interpolation points */, T* interpolation /**< [in] pointer to interpolation points at points of axis
                                                                                                                                                                                                                                                              */
                                          ,
                                          int ni /**< [in] number of elements in arrays */, std::string fname /**< [in] name of output file */
) {
    std::ofstream myfile(fname);
    myfile << std::setprecision(15);
    for (auto i = 0; i < ni; i++) {
        myfile << interpolationPoints[i] << "," << interpolation[i];
        if (i < n) {
            myfile << "," << axis[i] << "," << grid[i];
        }
        myfile << "\n";
    }
    myfile.close();
}

/**
 * \brief function to print the values for the trajectory debugging
 */
template <typename T>
void writeRoadTrajectoryCSV(T* interpolationx /**< [in] pointer to interpolation
                                                 points at points of axis */
                            ,
                            T* interpolationy /**< [in] pointer to interpolation
                                                 points at points of axis */
                            ,
                            T* interpolationz /**< [in] pointer to interpolation
                                             points at points of axis */
                            ,
                            int ni /**< [in] number of elements in arrays */, std::string fname /**< [in] name of output file */
) {
    std::ofstream myfile(fname);
    myfile << std::setprecision(15);
    for (auto i = 0; i < ni; i++) {
        myfile << interpolationx[i] << "," << interpolationy[i] << "," << ((interpolationz == nullptr) ? 0 : interpolationz[i]) << "\n";
    }
    myfile.close();
}

/**
 * Checks if a file exists.
 * \param[in] filename The name of the file to be checked.
 * \throw std::system_error if the file does not exist.
 */
void checkFileExists(const std::string& filename);

template <typename T>
class MyFile {
public:
    MyFile(std::string filename){
        myfile.open(filename);
        myfile << std::setprecision(15);
    };
    ~MyFile(){
        myfile.close();
    };
    /**
    * \brief write one time step to the file
    */
    void writeSoltionVector(
        T* solVec, /**< of size Constants::VEC_DIM * Constants::DIM = 27 */
        T* lagVelVec, /**< Constants::DIM - 1) * Constants::VEC_DIM [CG:XY, W_fl:XY, T_fl:XY, W_fr:XY, T_fr:XY, W_rl:XY, T_rl:XY, W_rr:XY, T_rr:XY] */
        T* angVec /**< of size 3 */
    ) {
        myfile << 0 << ",";  // 1 wc[x]
        myfile << 0 << ",";  // 2 wc[z]
        myfile << 0 << ",";  // 3 wc[y]
        myfile << lagVelVec[0] << ",";  // 4 vc[x] currentVelocityLagrangian[0]
        myfile << 0 << ",";  // 5 vc[z]
        myfile << lagVelVec[1] << ",";  // 6 vc[y] currentVelocityLagrangian[1]
        myfile << lagVelVec[7 * (Constants::DIM - 1)] << ",";  // 7 vw_rr[x]
        myfile << 0 << ",";  // 8 vw1[z]
        myfile << lagVelVec[7 * (Constants::DIM - 1) + 1] << ",";  // 9 vw_rr[y]
        myfile << lagVelVec[5 * (Constants::DIM - 1)] << ",";  // 10 vw_rl[x]
        myfile << 0 << ",";  // 11 vw2[z]
        myfile << lagVelVec[5 * (Constants::DIM - 1) + 1] << ",";  // 12 vw_rl[y]
        myfile << lagVelVec[1 * (Constants::DIM - 1)] << ",";  // 13 vw_fl[x]
        myfile << 0 << ",";  // 14 vw3[z]
        myfile << lagVelVec[1 * (Constants::DIM - 1) + 1] << ",";  // 15 vw_fl[y]
        myfile << lagVelVec[3 * (Constants::DIM - 1)] << ",";  // 16 vw_fr[x]
        myfile << 0 << ",";  // 17 vw4[z]
        myfile << lagVelVec[3 * (Constants::DIM - 1) + 1] << ",";  // 18 vw_Fr[y]
        myfile << lagVelVec[8 * (Constants::DIM - 1)] << ",";  // 19 vt1[x]
        myfile << 0 << ",";  // 20 vt1[z]
        myfile << lagVelVec[8 * (Constants::DIM - 1) + 1] << ",";  // 21 vt1[y]
        myfile << lagVelVec[6 * (Constants::DIM - 1)] << ",";  // 22 vt2[x]
        myfile << 0 << ",";  // 23 vt2[z]
        myfile << lagVelVec[6 * (Constants::DIM - 1) + 1] << ",";  // 24 vt2[y]
        myfile << lagVelVec[2 * (Constants::DIM - 1)] << ",";  // 25 vt3[x]
        myfile << 0 << ",";  // 26 vt3[z]
        myfile << lagVelVec[2 * (Constants::DIM - 1) + 1] << ",";  // 27 vt3[y]
        myfile << lagVelVec[4 * (Constants::DIM - 1)] << ",";  // 28 vt4[x]
        myfile << 0 << ",";  // 29 vt4[z]
        myfile << lagVelVec[4 * (Constants::DIM - 1) + 1] << ",";  // 30 vt4[y]
        myfile << angVec[0] << ",";  // 31 qc[0] // here x
        myfile << angVec[2] << ",";  // 32 qc[1] // here z
        myfile << angVec[1] << ",";  // 33 qc[2] // here y
        myfile << 0 << ",";  // 34 qc[3]
        myfile << solVec[0] << ",";  // 35 pcc[x]
        myfile << solVec[2] << ",";  // 36 pcc[z]
        myfile << solVec[1] << ",";  // 37 pcc[y]
        myfile << solVec[5 * Constants::DIM + 0] << ",";  // 38 pw_rr[x]
        myfile << solVec[5 * Constants::DIM + 2] << ",";  // 39 pw_rr[z]
        myfile << solVec[5 * Constants::DIM + 1] << ",";  // 40 pw_rr[y]
        myfile << solVec[7 * Constants::DIM + 0] << ",";  // 41 pw_rl[x]
        myfile << solVec[7 * Constants::DIM + 2] << ",";  // 42 pw_rl[z]
        myfile << solVec[7 * Constants::DIM + 1] << ",";  // 43 pw_rl[y]
        myfile << solVec[3 * Constants::DIM + 0] << ",";  // 44 pw_fl[x]
        myfile << solVec[3 * Constants::DIM + 2] << ",";  // 45 pw_fl[z]
        myfile << solVec[3 * Constants::DIM + 1] << ",";  // 46 pw_fl[y]
        myfile << solVec[1 * Constants::DIM + 0] << ",";  // 47 pw_fr[x]
        myfile << solVec[1 * Constants::DIM + 2] << ",";  // 48 pw_fr[z]
        myfile << solVec[1 * Constants::DIM + 1] << ",";  // 49 pw_fr[y]
        myfile << solVec[6 * Constants::DIM + 0] << ",";  // 50 pt_rr[x]
        myfile << solVec[6 * Constants::DIM + 2] << ",";  // 51 pt_rr[z]
        myfile << solVec[6 * Constants::DIM + 1] << ",";  // 52 pt_rr[y]
        myfile << solVec[8 * Constants::DIM + 0] << ",";  // 53 pt_rl[x]
        myfile << solVec[8 * Constants::DIM + 2] << ",";  // 54 pt_rl[z]
        myfile << solVec[8 * Constants::DIM + 1] << ",";  // 55 pt_rl[y]
        myfile << solVec[4 * Constants::DIM + 0] << ",";  // 56 pt_fl[x]
        myfile << solVec[4 * Constants::DIM + 2] << ",";  // 57 pt_fl[z]
        myfile << solVec[4 * Constants::DIM + 1] << ",";  // 58 pt_fl[y]
        myfile << solVec[2 * Constants::DIM + 0] << ",";  // 59 pt_fr[x]
        myfile << solVec[2 * Constants::DIM + 2] << ",";  // 60 pt_fr[z]
        myfile << solVec[2 * Constants::DIM + 1] << "\n";  // 61 pt_fr[y]
    };

    /**
     * \brief write one time step to the file
     */
    void writeSoltionVectorMBD(T* solVec    /**< of size 61 */    ) {

        myfile << solVec[0] << ",";                                // 1 wc[x]
        myfile << solVec[2] << ",";                                // 2 wc[z]
        myfile << solVec[1] << ",";                                // 3 wc[y]
        myfile << solVec[3] << ",";                                // 4 vc[x]
        myfile << solVec[5] << ",";                                // 5 vc[z]
        myfile << solVec[4] << ",";                                // 6 vc[y] 
        myfile << solVec[6 + 3 * Constants::DIM + 0] << ",";       // 7 vw_rr[x]
        myfile << solVec[6 + 3 * Constants::DIM + 2] << ",";       // 8 vw_rr[z]
        myfile << solVec[6 + 3 * Constants::DIM + 1] << ",";       // 9 vw_rr[y]
        myfile << solVec[6 + 2 * Constants::DIM + 0] << ",";       // 10 vw_rl[x]
        myfile << solVec[6 + 2 * Constants::DIM + 2] << ",";       // 11 vw_rl[z]
        myfile << solVec[6 + 2 * Constants::DIM + 1] << ",";       // 12 vw_rl[y]
        myfile << solVec[6 + 0 * Constants::DIM + 0] << ",";       // 13 vw_fl[x]
        myfile << solVec[6 + 0 * Constants::DIM + 2] << ",";       // 14 vw_fl[z]
        myfile << solVec[6 + 0 * Constants::DIM + 1] << ",";       // 15 vw_fl[y]
        myfile << solVec[6 + 1 * Constants::DIM + 0] << ",";       // 16 vw_fr[x]
        myfile << solVec[6 + 1 * Constants::DIM + 2] << ",";       // 17 vw_fr[z]
        myfile << solVec[6 + 1 * Constants::DIM + 1] << ",";       // 18 vw_fr[y]
        myfile << solVec[18 + 3 * Constants::DIM + 0] << ",";      // 19 vt_rr[x]
        myfile << solVec[18 + 3 * Constants::DIM + 2] << ",";      // 20 vt_rr[z]
        myfile << solVec[18 + 3 * Constants::DIM + 1] << ",";      // 21 vt_rr[y]
        myfile << solVec[18 + 2 * Constants::DIM + 0] << ",";      // 22 vt_rl[x]
        myfile << solVec[18 + 2 * Constants::DIM + 2] << ",";      // 23 vt_rl[z]
        myfile << solVec[18 + 2 * Constants::DIM + 1] << ",";      // 24 vt_rl[y]
        myfile << solVec[18 + 0 * Constants::DIM + 0] << ",";      // 25 vt_fl[x]
        myfile << solVec[18 + 0 * Constants::DIM + 2] << ",";      // 26 vt_fl[z]
        myfile << solVec[18 + 0 * Constants::DIM + 1] << ",";      // 27 vt_fl[y]
        myfile << solVec[18 + 1 * Constants::DIM + 0] << ",";      // 28 vt_fr[x]
        myfile << solVec[18 + 1 * Constants::DIM + 2] << ",";      // 29 vt_fr[z]
        myfile << solVec[18 + 1 * Constants::DIM + 1] << ",";      // 30 vt_fr[y]
        myfile << solVec[30] << ",";                               // 31 qc[0] 
        myfile << -solVec[32] << ",";                               // 32 qc[1] 
        myfile << -solVec[31] << ",";                               // 33 qc[2] 
        myfile << solVec[33] << ",";                               // 34 qc[3]
        myfile << solVec[34] << ",";                               // 35 pcc[x]
        myfile << solVec[36] << ",";                               // 36 pcc[z]
        myfile << solVec[35] << ",";                               // 37 pcc[y]
        myfile << solVec[37 + 3 * Constants::DIM + 0] << ",";       // 38 pw_rr[x]
        myfile << solVec[37 + 3 * Constants::DIM + 2] << ",";       // 39 pw_rr[z]
        myfile << solVec[37 + 3 * Constants::DIM + 1] << ",";       // 40 pw_rr[y]
        myfile << solVec[37 + 2 * Constants::DIM + 0] << ",";       // 41 pw_rl[x]
        myfile << solVec[37 + 2 * Constants::DIM + 2] << ",";       // 42 pw_rl[z]
        myfile << solVec[37 + 2 * Constants::DIM + 1] << ",";       // 43 pw_rl[y]
        myfile << solVec[37 + 0 * Constants::DIM + 0] << ",";       // 44 pw_fl[x]
        myfile << solVec[37 + 0 * Constants::DIM + 2] << ",";       // 45 pw_fl[z]
        myfile << solVec[37 + 0 * Constants::DIM + 1] << ",";       // 46 pw_fl[y]
        myfile << solVec[37 + 1 * Constants::DIM + 0] << ",";       // 47 pw_fr[x]
        myfile << solVec[37 + 1 * Constants::DIM + 2] << ",";       // 48 pw_fr[z]
        myfile << solVec[37 + 1 * Constants::DIM + 1] << ",";       // 49 pw_fr[y]
        myfile << solVec[49 + 3 * Constants::DIM + 0] << ",";      // 50 pt_rr[x]
        myfile << solVec[49 + 3 * Constants::DIM + 2] << ",";      // 51 pt_rr[z]
        myfile << solVec[49 + 3 * Constants::DIM + 1] << ",";      // 52 pt_rr[y]
        myfile << solVec[49 + 2 * Constants::DIM + 0] << ",";      // 53 pt_rl[x]
        myfile << solVec[49 + 2 * Constants::DIM + 2] << ",";      // 54 pt_rl[z]
        myfile << solVec[49 + 2 * Constants::DIM + 1] << ",";      // 55 pt_rl[y]
        myfile << solVec[49 + 0 * Constants::DIM + 0] << ",";      // 56 pt_fl[x]
        myfile << solVec[49 + 0 * Constants::DIM + 2] << ",";      // 57 pt_fl[z]
        myfile << solVec[49 + 0 * Constants::DIM + 1] << ",";      // 58 pt_fl[y]
        myfile << solVec[49 + 1 * Constants::DIM + 0] << ",";      // 59 pt_fr[x]
        myfile << solVec[49 + 1 * Constants::DIM + 2] << ",";      // 60 pt_fr[z]
        myfile << solVec[49 + 1 * Constants::DIM + 1] << "\n";      // 61 pt_fr[y]
    }

    void writeSingleValue(T value) { myfile << value << ","; }
    void writeSingleValue(int value) { myfile << value << ","; }

    void newLine() { myfile << "\n"; }
    

    void writeParameters() { 
        writeSingleValue(MetaDatabase<T>::getDatabase().getTimeStepSize());

        writeSingleValue(-MetaDatabase<T>::getDatabase().getLongitudalLegPositionRearRight());
        writeSingleValue(-MetaDatabase<T>::getDatabase().getLongitudalLegPositionRearLeft());
        writeSingleValue(MetaDatabase<T>::getDatabase().getLongitudalLegPositionFrontLeft());
        writeSingleValue(MetaDatabase<T>::getDatabase().getLongitudalLegPositionFrontRight());

        writeSingleValue(MetaDatabase<T>::getDatabase().getLatidudalLegPositionRearRight());
        writeSingleValue(-MetaDatabase<T>::getDatabase().getLatidudalLegPositionRearLeft());
        writeSingleValue(-MetaDatabase<T>::getDatabase().getLatidudalLegPositionFrontLeft());
        writeSingleValue(MetaDatabase<T>::getDatabase().getLatidudalLegPositionFrontRight());

        writeSingleValue(MetaDatabase<T>::getDatabase().getBodySpringLengthRearRight());
        writeSingleValue(MetaDatabase<T>::getDatabase().getBodySpringLengthRearLeft());
        writeSingleValue(MetaDatabase<T>::getDatabase().getBodySpringLengthFrontLeft());
        writeSingleValue(MetaDatabase<T>::getDatabase().getBodySpringLengthFrontRight());

        writeSingleValue(MetaDatabase<T>::getDatabase().getTyreSpringLengthRearRight());
        writeSingleValue(MetaDatabase<T>::getDatabase().getTyreSpringLengthRearLeft());
        writeSingleValue(MetaDatabase<T>::getDatabase().getTyreSpringLengthFrontLeft());
        writeSingleValue(MetaDatabase<T>::getDatabase().getTyreSpringLengthFrontRight());

    }


    void writeSolutionMatrix(
         T* solVec, /**< of size sol_size * Constants::VEC_DIM * Constants::DIM = 27 */
        T* velVec,
        T* angVec,
        size_t sol_size
    ) {
        for (auto i = 0; i < sol_size; i++) {
            writeSoltionVector(
                solVec + i * Constants::VEC_DIM * Constants::DIM,
                velVec + i * (Constants::DIM - 1) * Constants::VEC_DIM,
                angVec + i * Constants::DIM
                );
        }
    };

    void writeSolutionMatrixMBD(T* solVec, size_t sol_size) {
        for (auto i = 0; i < sol_size; i++) {
            writeSoltionVectorMBD(solVec + i * 61);
        }
    };


private:
    std::ofstream myfile;
};

}  // namespace IO
}  // namespace EVAA

