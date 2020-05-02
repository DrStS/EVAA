// TODO: Copyright header

#pragma once

#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

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
                            int ni /**< [in] number of elements in arrays */, std::string fname /**< [in] name of output file */
) {
    std::ofstream myfile(fname);
    myfile << std::setprecision(15);
    for (auto i = 0; i < ni; i++) {
        myfile << interpolationx[i] << "," << interpolationy[i] << "\n";
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
        T* solVec /**< of size Constants::VEC_DIM * Constants::DIM = 27 */
    ) {
        myfile << solVec[0] << ",";  // 1 wc[x]
        myfile << solVec[0] << ",";  // 2 wc[z]
        myfile << solVec[0] << ",";  // 3 wc[y]
        myfile << solVec[0] << ",";  // 4 vc[x]
        myfile << solVec[0] << ",";  // 5 vc[z]
        myfile << solVec[0] << ",";  // 6 vc[y]
        myfile << solVec[0] << ",";  // 7 vw1[x]
        myfile << solVec[0] << ",";  // 8 vw1[z]
        myfile << solVec[0] << ",";  // 9 vw1[y]
        myfile << solVec[0] << ",";  // 10 vw2[x]
        myfile << solVec[0] << ",";  // 11 vw2[z]
        myfile << solVec[0] << ",";  // 12 vw2[y]
        myfile << solVec[0] << ",";  // 13 vw3[x]
        myfile << solVec[0] << ",";  // 14 vw3[z]
        myfile << solVec[0] << ",";  // 15 vw3[y]
        myfile << solVec[0] << ",";  // 16 vw4[x]
        myfile << solVec[0] << ",";  // 17 vw4[z]
        myfile << solVec[0] << ",";  // 18 vw4[y]
        myfile << solVec[0] << ",";  // 19 vt1[x]
        myfile << solVec[0] << ",";  // 20 vt1[z]
        myfile << solVec[0] << ",";  // 21 vt1[y]
        myfile << solVec[0] << ",";  // 22 vt2[x]
        myfile << solVec[0] << ",";  // 23 vt2[z]
        myfile << solVec[0] << ",";  // 24 vt2[y]
        myfile << solVec[0] << ",";  // 25 vt3[x]
        myfile << solVec[0] << ",";  // 26 vt3[z]
        myfile << solVec[0] << ",";  // 27 vt3[y]
        myfile << solVec[0] << ",";  // 28 vt4[x]
        myfile << solVec[0] << ",";  // 29 vt4[z]
        myfile << solVec[0] << ",";  // 30 vt4[y]
        myfile << solVec[0] << ",";  // 31 qc[0]
        myfile << solVec[0] << ",";  // 32 qc[1]
        myfile << solVec[0] << ",";  // 33 qc[2]
        myfile << solVec[0] << ",";  // 34 qc[3]
        myfile << solVec[0] << ",";  // 35 pcc[x]
        myfile << solVec[0] << ",";  // 36 pcc[z]
        myfile << solVec[0] << ",";  // 37 pcc[y]
        myfile << solVec[0] << ",";  // 38 pw1[x]
        myfile << solVec[0] << ",";  // 39 pw1[z]
        myfile << solVec[0] << ",";  // 40 pw1[y]
        myfile << solVec[0] << ",";  // 41 pw2[x]
        myfile << solVec[0] << ",";  // 42 pw2[z]
        myfile << solVec[0] << ",";  // 43 pw2[y]
        myfile << solVec[0] << ",";  // 44 pw3[x]
        myfile << solVec[0] << ",";  // 45 pw3[z]
        myfile << solVec[0] << ",";  // 46 pw3[y]
        myfile << solVec[0] << ",";  // 47 pw4[x]
        myfile << solVec[0] << ",";  // 48 pw4[z]
        myfile << solVec[0] << ",";  // 49 pw4[y]
        myfile << solVec[0] << ",";  // 50 pt1[x]
        myfile << solVec[0] << ",";  // 51 pt1[z]
        myfile << solVec[0] << ",";  // 52 pt1[y]
        myfile << solVec[0] << ",";  // 53 pt2[x]
        myfile << solVec[0] << ",";  // 54 pt2[z]
        myfile << solVec[0] << ",";  // 55 pt2[y]
        myfile << solVec[0] << ",";  // 56 pt3[x]
        myfile << solVec[0] << ",";  // 57 pt3[z]
        myfile << solVec[0] << ",";  // 58 pt3[y]
        myfile << solVec[0] << ",";  // 59 pt4[x]
        myfile << solVec[0] << ",";  // 60 pt4[z]
        myfile << solVec[0] << "\n";  // 61 pt4[y]
    };

    void writeSolutionMatrix(
         T* solVec, /**< of size sol_size * Constants::VEC_DIM * Constants::DIM = 27 */
        size_t sol_size
    ) {
        for (auto i = 0; i < sol_size; i++) {
            writeSoltionVector(solVec + i * Constants::VEC_DIM * Constants::DIM);
        }
    };

private:
    std::ofstream myfile;
};

}  // namespace IO
}  // namespace EVAA

