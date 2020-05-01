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
 * \brief function to print the values for the lookup table and the interpolated values for
 * debugging
 */
template <typename T>
void writeLookUpGridPlusInterpolateValues(
    T* axis /**< [in] pointer to axis */, T* grid /**< [in] pointer to grid */,
    int n /**< [in] number of elements in arrays */,
    T* interpolationPoints /**< [in] pointer to interpolation points */,
    T* interpolation /**< [in] pointer to interpolation points at points of axis */,
    int ni /**< [in] number of elements in arrays */,
    std::string fname /**< [in] name of output file */
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
void writeRoadTrajectoryCSV(
    T* interpolationx /**< [in] pointer to interpolation points at points of axis */,
    T* interpolationy /**< [in] pointer to interpolation points at points of axis */,
    int ni /**< [in] number of elements in arrays */,
    std::string fname /**< [in] name of output file */
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

}  // namespace IO
}  // namespace EVAA
