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

TEST_F(OutputTest, ConstructOutputHDF5Class) {
    std::string _fileName = "Output1.hdf5";
    std::string _filePath = "";
    EVAA::HDF5::OutputHDF5<double> testObj(_fileName, _filePath);
}

TEST_F(OutputTest, ContainersInOutputHDF5Class) {
    std::string _fileName = "Output1.hdf5";
    std::string _filePath = "";
    EVAA::HDF5::OutputHDF5<double> testObj(_fileName, _filePath);
    testObj.createContainer(1);
    testObj.openContainer(1);
    testObj.closeContainer();
}

TEST_F(OutputTest, WriteVectorOutputHDF5) {
    std::string _fileName = "OutputWriteVector.hdf5";
    std::string _filePath = "";
    EVAA::HDF5::OutputHDF5<double> testObj(_fileName, _filePath);
    testObj.createContainer(true);

    const size_t vecDim = 6;
    double vecWrite[vecDim]{5, 9, 91, 26, 5, 94};
    std::string _datasetName = "SolutionVector";
    testObj.WriteVector(_datasetName, vecWrite, vecDim);

    testObj.closeContainer();
}

TEST_F(OutputTest, ReadVectorFileHDF5) {
    std::string _fileName = "OutputWriteVector.hdf5";
    std::string _filePath = "";

    H5::H5File* myHDF5FileHandle = new H5::H5File(_filePath + _fileName, H5F_ACC_RDONLY);
    std::string _datasetName = "SolutionVector";
    H5::Exception::dontPrint();

    H5::DataSet dataset = myHDF5FileHandle->openDataSet(_datasetName);

    const size_t vecDim = 6;
    double vecWrite[vecDim]{5, 9, 91, 26, 5, 94};

    size_t vecDimRead = 0;
    double vecRead[vecDim];
    for (auto i = 0; i < vecDim; ++i) vecRead[i] = 0;

    dataset.read(vecRead, H5::PredType::NATIVE_DOUBLE);
    H5::Attribute attribute = dataset.openAttribute("Vector # elements");
    attribute.read(H5::PredType::NATIVE_INT, &vecDimRead);

    delete myHDF5FileHandle;

    // check
    EXPECT_EQ(vecDim, vecDimRead) << "Vectors length differ! ";
    for (int i = 0; i < vecDim; ++i) {
        EXPECT_EQ(vecWrite[i], vecRead[i]) << "Vectors differ at index " << i;
    }
}

TEST_F(OutputTest, ReadWritenVectorOutputHDF5) {
    H5::Exception::dontPrint();
    std::string _fileName = "OutputWriteVector.hdf5";
    std::string _filePath = "";

    EVAA::HDF5::OutputHDF5<double> testObj(_fileName, _filePath);
    testObj.openContainer(false);

    const size_t vecDim = 6;
    double vecWrite[vecDim]{5, 9, 91, 26, 5, 94};
    std::string _datasetName = "SolutionVector";

    size_t vecDimRead = 0;
    double vecRead[vecDim];
    for (auto i = 0; i < vecDim; ++i) vecRead[i] = 0;

    testObj.ReadVector(_datasetName, vecRead, vecDimRead);

    EXPECT_EQ(vecDim, vecDimRead) << "Vectors length differ! ";
    for (int i = 0; i < vecDim; ++i) {
        EXPECT_EQ(vecWrite[i], vecRead[i]) << "Vectors differ at index " << i;
    }

    testObj.closeContainer();
}

TEST_F(OutputTest, WriteMatrixOutputHDF5) {
    const int matrixRows = 3;
    const int matrixColumns = 2;
    double matrixWrite[matrixRows*matrixColumns]{5, 9, 91, 26, 5, 94};    

    // Create OutputHDF5 object
    std::string _fileName = "OutputWriteMatrix.hdf5";
    std::string _filePath = "";
    EVAA::HDF5::OutputHDF5<double> testObj(_fileName, _filePath);
    
    // Write part
    testObj.createContainer(true);
    std::string _datasetName = "SolutionMatrix";
    testObj.WriteMatrix(_datasetName, matrixWrite, matrixRows, matrixColumns);
    testObj.closeContainer();

    // Read part
    double matrixRead[matrixRows * matrixColumns];
    size_t matrixRowsRead = 0;
    size_t matrixColumnsRead = 0;
    testObj.openContainer(false);
    testObj.ReadMatrix(_datasetName, matrixRead, matrixRowsRead, matrixColumnsRead);
    
    EXPECT_EQ(matrixRows, matrixRowsRead) << "Matrix rows length differ! ";
    EXPECT_EQ(matrixColumns, matrixColumnsRead) << "Matrix rows length differ! ";
    for (auto i = 0; i < matrixRows; ++i) {
        for (auto j = 0; j < matrixColumns; ++j) {
            EXPECT_EQ(matrixWrite[i], matrixRead[i])
                << "Matrix differ at index (" << i << ", " << j << ")\n";
        }
    }
}

#endif
