#include <array>

#include "OutputHDF5.h"
#include "gtest/gtest.h"

class OutputTest : public ::testing::Test {};

#ifdef USE_HDF5
TEST_F(OutputTest, writeVectorHDF5file) {
    const int vecDim = 6;
    std::array<double, vecDim> vec{5, 9, 91, 26, 5, 94};
    std::array<double, vecDim> vec_another;

    writeVectorToFile("outputVector.h5", "vecTOwrite", vec.data(), vecDim);

    readVectorFromFile("outputVector.h5", "vecTOwrite", vec_another.data(), vecDim);

    for (int i = 0; i < vec.size(); ++i) {
        EXPECT_EQ(vec[i], vec_another[i]) << "Vectors differ at index " << i;
    }
}

TEST_F(OutputTest, writeMatrixHDF5file) {
    const int matrixRows = 3;
    const int matrixColumns = 2;
    std::array<double, matrixRows * matrixColumns> matrix{5, 9, 91, 26, 5, 94};
    std::array<double, matrixRows * matrixColumns> matrix_another;

    writeMatrixToFile("outputMatrix.h5", "matrixTOwrite", matrix.data(), matrixRows, matrixColumns);
    
    readMatrixFromFile("outputMatrix.h5", "matrixTOwrite", matrix_another.data(), matrixRows,
                       matrixColumns);

    for (auto i = 0; i < matrixRows; ++i) {
        for (auto j = 0; j < matrixColumns; ++j) {
            EXPECT_EQ(matrix[i], matrix_another[i])
                << "Matrix differ at index (" << i << ", " << j << ")\n";
        }
    }
}

#endif