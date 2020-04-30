// TODO: Copyright header

#pragma once

#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

namespace EVAA {

namespace IO {
template <typename T>
void write_matrix(T* mat, std::string fname, size_t rows, size_t cols) {
    std::ofstream myfile(fname);
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            myfile << mat[i * cols + j] << std::setprecision(15) << " ";
        }
        myfile << "\n";
    }
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

};  // namespace IO
}  // namespace EVAA
