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
    testObj.createContainer(true);
    testObj.closeContainer();

    testObj.openContainer(false);
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

TEST_F(OutputTest, TestHandles) {
    H5::Exception::dontPrint();
    std::string _fileName = "OutputWriteVector.hdf5";
    std::string _filePath = "";

    H5::H5File* myHDF5FileHandle = new H5::H5File(_fileName + _filePath, H5F_ACC_RDONLY);
    H5::Group* myHDF5GroupOperators = new H5::Group(myHDF5FileHandle->openGroup("Default"));

    myHDF5GroupOperators->close();
    myHDF5FileHandle->close();
    /*EVAA::HDF5::OutputHDF5<double> testObj(_fileName, _filePath);
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

    testObj.closeContainer();*/
}

TEST_F(OutputTest, WriteReadMatrixOutputHDF5) {
    const int matrixRows = 3;
    const int matrixColumns = 2;
    double matrixWrite[matrixRows * matrixColumns]{5, 9, 91, 26, 5, 94};

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

TEST_F(OutputTest, Write3Matrices1File1Group) {
    const int matrixRows1 = 3, matrixColumns1 = 2;
    const int matrixRows2 = 3, matrixColumns2 = 2;
    const int matrixRows3 = 3, matrixColumns3 = 2;
    double matrixWrite1[matrixRows1 * matrixColumns1]{0, 1, 2, 3, 4, 5};
    double matrixWrite2[matrixRows2 * matrixColumns2]{6, 7, 8, 9, 10, 11};
    double matrixWrite3[matrixRows3 * matrixColumns3]{12, 13, 14, 15, 16, 17};

    // Create OutputHDF5 object
    std::string _fileName = "Write3Matrices1File1Group.hdf5";
    std::string _filePath = "";
    EVAA::HDF5::OutputHDF5<double> testObj(_fileName, _filePath);

    // Write part

    // Matrix 1
    testObj.createContainer(true);
    std::string _datasetName1 = "Matrix 1";
    testObj.WriteMatrix(_datasetName1, matrixWrite1, matrixRows1, matrixColumns1);
    testObj.closeContainer();

    // Matrix 2
    testObj.openContainer(true);
    std::string _datasetName2 = "Matrix 2";
    testObj.WriteMatrix(_datasetName2, matrixWrite2, matrixRows2, matrixColumns2);
    testObj.closeContainer();

    // Matrix 3
    std::string groupName3 = "Group of Matrix 3";
    testObj.createContainer(false, groupName3);
    std::string _datasetName3 = "Matrix 3";
    testObj.WriteMatrix(_datasetName3, matrixWrite3, matrixRows3, matrixColumns3,
                        EVAA::HDF5FileHandle::GROUP);
    testObj.closeContainer();
}

TEST_F(OutputTest, WriteRead3Matrices1File1Group) {
    const int matrixRows1 = 3, matrixColumns1 = 2;
    const int matrixRows2 = 3, matrixColumns2 = 2;
    const int matrixRows3 = 3, matrixColumns3 = 2;
    double matrixWrite1[matrixRows1 * matrixColumns1]{0, 1, 2, 3, 4, 5};
    double matrixWrite2[matrixRows2 * matrixColumns2]{6, 7, 8, 9, 10, 11};
    double matrixWrite3[matrixRows3 * matrixColumns3]{12, 13, 14, 15, 16, 17};

    // Create OutputHDF5 object
    std::string _fileName = "WriteRead3Matrices1File1Group.hdf5";
    std::string _filePath = "";
    EVAA::HDF5::OutputHDF5<double> testObj(_fileName, _filePath);

    // Write part

    // Matrix 1
    testObj.createContainer(true);
    std::string _datasetName1 = "Matrix 1";
    testObj.WriteMatrix(_datasetName1, matrixWrite1, matrixRows1, matrixColumns1);
    testObj.closeContainer();

    // Matrix 2
    testObj.openContainer(true);
    std::string _datasetName2 = "Matrix 2";
    testObj.WriteMatrix(_datasetName2, matrixWrite2, matrixRows2, matrixColumns2);
    testObj.closeContainer();

    // Matrix 3
    std::string groupName3 = "Group of Matrix 3";
    testObj.createContainer(false, groupName3);
    std::string _datasetName3 = "Matrix 3";
    testObj.WriteMatrix(_datasetName3, matrixWrite3, matrixRows3, matrixColumns3,
                        EVAA::HDF5FileHandle::GROUP);
    testObj.closeContainer();

    // Read part
    double matrixRead1[matrixRows1 * matrixColumns1];
    double matrixRead2[matrixRows2 * matrixColumns2];
    double matrixRead3[matrixRows3 * matrixColumns3];
    size_t matrixRowsRead1 = 0, matrixColumnsRead1 = 0;
    size_t matrixRowsRead2 = 0, matrixColumnsRead2 = 0;
    size_t matrixRowsRead3 = 0, matrixColumnsRead3 = 0;
    testObj.openContainer(false);
    testObj.ReadMatrix(_datasetName1, matrixRead1, matrixRowsRead1, matrixColumnsRead1);
    testObj.ReadMatrix(_datasetName2, matrixRead2, matrixRowsRead2, matrixColumnsRead2);
    testObj.closeContainer();

    testObj.openContainer(false, groupName3);
    testObj.ReadMatrix(_datasetName3, matrixRead3, matrixRowsRead3, matrixColumnsRead3,
                       EVAA::HDF5FileHandle::GROUP);
    testObj.closeContainer();

    // Compare part
    EXPECT_EQ(matrixRows1, matrixRowsRead1) << "Matrix 1 rows length differ! ";
    EXPECT_EQ(matrixColumns1, matrixColumnsRead1) << "Matrix 1 columns length differ! ";
    for (auto i = 0; i < matrixRows1; ++i) {
        for (auto j = 0; j < matrixColumns1; ++j) {
            EXPECT_EQ(matrixWrite1[i], matrixRead1[i])
                << "Matrix 1 differ at index (" << i << ", " << j << ")\n";
        }
    }

    EXPECT_EQ(matrixRows2, matrixRowsRead2) << "Matrix 2 rows length differ! ";
    EXPECT_EQ(matrixColumns2, matrixColumnsRead2) << "Matrix 2 columns length differ! ";
    for (auto i = 0; i < matrixRows2; ++i) {
        for (auto j = 0; j < matrixColumns2; ++j) {
            EXPECT_EQ(matrixWrite2[i], matrixRead2[i])
                << "Matrix 2 differ at index (" << i << ", " << j << ")\n";
        }
    }

    EXPECT_EQ(matrixRows3, matrixRowsRead3) << "Matrix 3 rows length differ! ";
    EXPECT_EQ(matrixColumns3, matrixColumnsRead3) << "Matrix 3 columns length differ! ";
    for (auto i = 0; i < matrixRows3; ++i) {
        for (auto j = 0; j < matrixColumns3; ++j) {
            EXPECT_EQ(matrixWrite3[i], matrixRead3[i])
                << "Matrix 3 differ at index (" << i << ", " << j << ")\n";
        }
    }
}

TEST_F(OutputTest, WriteRead3Matrices2Groups) {
    const int matrixRows1 = 3, matrixColumns1 = 2;
    const int matrixRows2 = 3, matrixColumns2 = 2;
    const int matrixRows3 = 3, matrixColumns3 = 2;
    double matrixWrite1[matrixRows1 * matrixColumns1]{0, 1, 2, 3, 4, 5};
    double matrixWrite2[matrixRows2 * matrixColumns2]{6, 7, 8, 9, 10, 11};
    double matrixWrite3[matrixRows3 * matrixColumns3]{12, 13, 14, 15, 16, 17};

    // Create OutputHDF5 object
    std::string _fileName = "WriteRead3Matrices2Groups.hdf5";
    std::string _filePath = "";
    EVAA::HDF5::OutputHDF5<double> testObj(_fileName, _filePath);

    // Write part

    // Matrix 1
    testObj.createContainer(true);
    std::string _datasetName1 = "Matrix 1";
    testObj.WriteMatrix(_datasetName1, matrixWrite1, matrixRows1, matrixColumns1);
    testObj.closeContainer();

    // Matrix 2 in File
    testObj.openContainer(true);
    std::string _datasetName2 = "Matrix 2";
    testObj.WriteMatrix(_datasetName2, matrixWrite2, matrixRows2, matrixColumns2);
    testObj.closeContainer();

    // Matrix 2 in Group
    std::string groupName2 = "Group of Matrix 2";
    testObj.createContainer(false, groupName2);
    std::string _datasetName2G = "Matrix 2";
    testObj.WriteMatrix(_datasetName2G, matrixWrite2, matrixRows2, matrixColumns2,
                        EVAA::HDF5FileHandle::GROUP);
    testObj.closeContainer();

    // Matrix 3
    std::string groupName3 = "Group of Matrix 3";
    testObj.createContainer(false, groupName3);
    std::string _datasetName3 = "Matrix 3";
    testObj.WriteMatrix(_datasetName3, matrixWrite3, matrixRows3, matrixColumns3,
                        EVAA::HDF5FileHandle::GROUP);
    testObj.closeContainer();

    // Read part
    double matrixRead1[matrixRows1 * matrixColumns1];
    double matrixRead2[matrixRows2 * matrixColumns2];
    double matrixRead2G[matrixRows2 * matrixColumns2];
    double matrixRead3[matrixRows3 * matrixColumns3];
    size_t matrixRowsRead1 = 0, matrixColumnsRead1 = 0;
    size_t matrixRowsRead2 = 0, matrixColumnsRead2 = 0;
    size_t matrixRowsRead2G = 0, matrixColumnsRead2G = 0;
    size_t matrixRowsRead3 = 0, matrixColumnsRead3 = 0;
    testObj.openContainer(false);
    testObj.ReadMatrix(_datasetName1, matrixRead1, matrixRowsRead1, matrixColumnsRead1);
    testObj.ReadMatrix(_datasetName2, matrixRead2, matrixRowsRead2, matrixColumnsRead2);
    testObj.closeContainer();

    testObj.openContainer(false, groupName2);
    testObj.ReadMatrix(_datasetName2G, matrixRead2G, matrixRowsRead2G, matrixColumnsRead2G,
                       EVAA::HDF5FileHandle::GROUP);
    testObj.closeContainer();
    
    testObj.openContainer(false, groupName3);
    testObj.ReadMatrix(_datasetName3, matrixRead3, matrixRowsRead3, matrixColumnsRead3,
                       EVAA::HDF5FileHandle::GROUP);
    testObj.closeContainer();

    // Compare part

    // Matrix 1
    EXPECT_EQ(matrixRows1, matrixRowsRead1) << "Matrix 1 rows length differ! ";
    EXPECT_EQ(matrixColumns1, matrixColumnsRead1) << "Matrix 1 columns length differ! ";
    for (auto i = 0; i < matrixRows1; ++i) {
        for (auto j = 0; j < matrixColumns1; ++j) {
            EXPECT_EQ(matrixWrite1[i], matrixRead1[i])
                << "Matrix 1 differ at index (" << i << ", " << j << ")\n";
        }
    }

    // Matrix 2 File
    EXPECT_EQ(matrixRows2, matrixRowsRead2) << "Matrix 2 rows length differ! ";
    EXPECT_EQ(matrixColumns2, matrixColumnsRead2) << "Matrix 2 columns length differ! ";
    for (auto i = 0; i < matrixRows2; ++i) {
        for (auto j = 0; j < matrixColumns2; ++j) {
            EXPECT_EQ(matrixWrite2[i], matrixRead2[i])
                << "Matrix 2 differ at index (" << i << ", " << j << ")\n";
        }
    }

    // Matrix 2 Group
    EXPECT_EQ(matrixRows2, matrixRowsRead2G) << "Matrix 2 rows length differ! ";
    EXPECT_EQ(matrixColumns2, matrixColumnsRead2G) << "Matrix 2 columns length differ! ";
    for (auto i = 0; i < matrixRows2; ++i) {
        for (auto j = 0; j < matrixColumns2; ++j) {
            EXPECT_EQ(matrixWrite2[i], matrixRead2G[i])
                << "Matrix 2 differ at index (" << i << ", " << j << ")\n";
        }
    }

    // Matrix 3
    EXPECT_EQ(matrixRows3, matrixRowsRead3) << "Matrix 3 rows length differ! ";
    EXPECT_EQ(matrixColumns3, matrixColumnsRead3) << "Matrix 3 columns length differ! ";
    for (auto i = 0; i < matrixRows3; ++i) {
        for (auto j = 0; j < matrixColumns3; ++j) {
            EXPECT_EQ(matrixWrite3[i], matrixRead3[i])
                << "Matrix 3 differ at index (" << i << ", " << j << ")\n";
        }
    }
}

TEST_F(OutputTest, WriteReadMatrixFileGroup) {
    const int matrixRows = 3;
    const int matrixColumns = 2;
    double matrixWrite1[matrixRows * matrixColumns]{5, 9, 91, 26, 5, 94};
    double matrixWrite2[matrixRows * matrixColumns]{1, 2, 1, 2, 1, 2};

    // Create OutputHDF5 object
    std::string _fileName = "WriteReadMatrixFileGroup.hdf5";
    std::string _filePath = "";
    EVAA::HDF5::OutputHDF5<double> testObj(_fileName, _filePath);

    // Write part
    std::string groupName = "Group Test";
    testObj.createContainer(true, groupName);
    std::string _datasetName = "SolutionMatrix";
    testObj.WriteMatrix(_datasetName, matrixWrite1, matrixRows, matrixColumns);
    testObj.closeContainer();

    testObj.openContainer(true, groupName);
    testObj.WriteMatrix(_datasetName, matrixWrite2, matrixRows, matrixColumns,
                        EVAA::HDF5FileHandle::GROUP);
    testObj.closeContainer();

    // Read part
    double matrixRead1[matrixRows * matrixColumns];
    double matrixRead2[matrixRows * matrixColumns];
    size_t matrixRowsRead = 0;
    size_t matrixColumnsRead = 0;
    testObj.openContainer(false, groupName);
    testObj.ReadMatrix(_datasetName, matrixRead1, matrixRowsRead, matrixColumnsRead);
    testObj.ReadMatrix(_datasetName, matrixRead2, matrixRowsRead, matrixColumnsRead,
                       EVAA::HDF5FileHandle::GROUP);
    testObj.closeContainer();

    EXPECT_EQ(matrixRows, matrixRowsRead) << "Matrix 1 rows length differ! ";
    EXPECT_EQ(matrixColumns, matrixColumnsRead) << "Matrix 1 columns length differ! ";
    for (auto i = 0; i < matrixRows; ++i) {
        for (auto j = 0; j < matrixColumns; ++j) {
            EXPECT_EQ(matrixWrite1[i], matrixRead1[i])
                << "Matrix 1 differ at index (" << i << ", " << j << ")\n";
        }
    }

    EXPECT_EQ(matrixRows, matrixRowsRead) << "Matrix 2 rows length differ! ";
    EXPECT_EQ(matrixColumns, matrixColumnsRead) << "Matrix 2 columns length differ! ";
    for (auto i = 0; i < matrixRows; ++i) {
        for (auto j = 0; j < matrixColumns; ++j) {
            EXPECT_EQ(matrixWrite2[i], matrixRead2[i])
                << "Matrix 2 differ at index (" << i << ", " << j << ")\n";
        }
    }
}

#endif
