#include <gtest/gtest.h>
#include <sys/stat.h>

#include <array>
#include <cstdio>

#include "OutputHDF5.h"

#ifdef USE_HDF5

namespace Test {

class TestOutputHDF5 : public ::testing::Test {
protected:
    void SetUp() override {
        _filePath = "";
        _fileName = "TestOutput.hdf5";
    }

    void TearDown() override { std::remove(_fileName.c_str()); }

    std::string _filePath;
    std::string _fileName;
};

TEST_F(TestOutputHDF5, ConstructOutputHDF5Class) {
    EVAA::HDF5::OutputHDF5<double> testObj(_filePath, _fileName);
}

TEST_F(TestOutputHDF5, ContainersInOutputHDF5Class) {
    EVAA::HDF5::OutputHDF5<double> testObj(_filePath, _fileName);
    testObj.CreateContainer(true);
    testObj.CloseContainer();

    testObj.OpenContainer(false);
    testObj.CloseContainer();
}

TEST_F(TestOutputHDF5, WriteVectorOutputHDF5) {
    EVAA::HDF5::OutputHDF5<double> testObj(_filePath, _fileName);
    testObj.CreateContainer(true);

    const size_t vecDim = 6;
    double vecWrite[vecDim]{5, 9, 91, 26, 5, 94};
    std::string _datasetName = "SolutionVector";
    testObj.WriteVector(_datasetName, vecWrite, vecDim);

    testObj.CloseContainer();
}

TEST_F(TestOutputHDF5, WriteReadMatrixOutputHDF5) {
    const int matrixRows = 3;
    const int matrixColumns = 2;
    double matrixWrite[matrixRows * matrixColumns]{5, 9, 91, 26, 5, 94};

    // Create OutputHDF5 object
    EVAA::HDF5::OutputHDF5<double> testObj(_filePath, _fileName);

    // Write part
    testObj.CreateContainer(true);
    std::string _datasetName = "SolutionMatrix";
    testObj.WriteMatrix(_datasetName, matrixWrite, matrixRows, matrixColumns);
    testObj.CloseContainer();

    // Read part
    double matrixRead[matrixRows * matrixColumns];
    size_t matrixRowsRead = 0;
    size_t matrixColumnsRead = 0;
    testObj.OpenContainer(false);
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

TEST_F(TestOutputHDF5, Write3Matrices1File1Group) {
    const int matrixRows1 = 3, matrixColumns1 = 2;
    const int matrixRows2 = 3, matrixColumns2 = 2;
    const int matrixRows3 = 3, matrixColumns3 = 2;
    double matrixWrite1[matrixRows1 * matrixColumns1]{0, 1, 2, 3, 4, 5};
    double matrixWrite2[matrixRows2 * matrixColumns2]{6, 7, 8, 9, 10, 11};
    double matrixWrite3[matrixRows3 * matrixColumns3]{12, 13, 14, 15, 16, 17};

    // Create OutputHDF5 object
    EVAA::HDF5::OutputHDF5<double> testObj(_filePath, _fileName);

    // Write part

    // Matrix 1
    testObj.CreateContainer(true);
    std::string _datasetName1 = "Matrix 1";
    testObj.WriteMatrix(_datasetName1, matrixWrite1, matrixRows1, matrixColumns1);
    testObj.CloseContainer();

    // Matrix 2
    testObj.OpenContainer(true);
    std::string _datasetName2 = "Matrix 2";
    testObj.WriteMatrix(_datasetName2, matrixWrite2, matrixRows2, matrixColumns2);
    testObj.CloseContainer();

    // Matrix 3
    std::string groupName3 = "Group of Matrix 3";
    testObj.CreateContainer(false, groupName3);
    std::string _datasetName3 = "Matrix 3";
    testObj.WriteMatrix(_datasetName3, matrixWrite3, matrixRows3, matrixColumns3,
                        EVAA::HDF5FileHandle::GROUP);
    testObj.CloseContainer();
}

TEST_F(TestOutputHDF5, WriteRead3Matrices1File1Group) {
    const int matrixRows1 = 3, matrixColumns1 = 2;
    const int matrixRows2 = 3, matrixColumns2 = 2;
    const int matrixRows3 = 3, matrixColumns3 = 2;
    double matrixWrite1[matrixRows1 * matrixColumns1]{0, 1, 2, 3, 4, 5};
    double matrixWrite2[matrixRows2 * matrixColumns2]{6, 7, 8, 9, 10, 11};
    double matrixWrite3[matrixRows3 * matrixColumns3]{12, 13, 14, 15, 16, 17};

    // Create OutputHDF5 object
    EVAA::HDF5::OutputHDF5<double> testObj(_filePath, _fileName);

    // Write part

    // Matrix 1
    testObj.CreateContainer(true);
    std::string _datasetName1 = "Matrix 1";
    testObj.WriteMatrix(_datasetName1, matrixWrite1, matrixRows1, matrixColumns1);
    testObj.CloseContainer();

    // Matrix 2
    testObj.OpenContainer(true);
    std::string _datasetName2 = "Matrix 2";
    testObj.WriteMatrix(_datasetName2, matrixWrite2, matrixRows2, matrixColumns2);
    testObj.CloseContainer();

    // Matrix 3
    std::string groupName3 = "Group of Matrix 3";
    testObj.CreateContainer(false, groupName3);
    std::string _datasetName3 = "Matrix 3";
    testObj.WriteMatrix(_datasetName3, matrixWrite3, matrixRows3, matrixColumns3,
                        EVAA::HDF5FileHandle::GROUP);
    testObj.CloseContainer();

    // Read part
    double matrixRead1[matrixRows1 * matrixColumns1];
    double matrixRead2[matrixRows2 * matrixColumns2];
    double matrixRead3[matrixRows3 * matrixColumns3];
    size_t matrixRowsRead1 = 0, matrixColumnsRead1 = 0;
    size_t matrixRowsRead2 = 0, matrixColumnsRead2 = 0;
    size_t matrixRowsRead3 = 0, matrixColumnsRead3 = 0;
    testObj.OpenContainer(false);
    testObj.ReadMatrix(_datasetName1, matrixRead1, matrixRowsRead1, matrixColumnsRead1);
    testObj.ReadMatrix(_datasetName2, matrixRead2, matrixRowsRead2, matrixColumnsRead2);
    testObj.CloseContainer();

    testObj.OpenContainer(false, groupName3);
    testObj.ReadMatrix(_datasetName3, matrixRead3, matrixRowsRead3, matrixColumnsRead3,
                       EVAA::HDF5FileHandle::GROUP);
    testObj.CloseContainer();

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

TEST_F(TestOutputHDF5, WriteRead3Matrices2Groups) {
    const int matrixRows1 = 3, matrixColumns1 = 2;
    const int matrixRows2 = 3, matrixColumns2 = 2;
    const int matrixRows3 = 3, matrixColumns3 = 2;
    double matrixWrite1[matrixRows1 * matrixColumns1]{0, 1, 2, 3, 4, 5};
    double matrixWrite2[matrixRows2 * matrixColumns2]{6, 7, 8, 9, 10, 11};
    double matrixWrite3[matrixRows3 * matrixColumns3]{12, 13, 14, 15, 16, 17};

    // Create OutputHDF5 object
    EVAA::HDF5::OutputHDF5<double> testObj(_filePath, _fileName);

    // Write part

    // Matrix 1
    testObj.CreateContainer(true);
    std::string _datasetName1 = "Matrix 1";
    testObj.WriteMatrix(_datasetName1, matrixWrite1, matrixRows1, matrixColumns1);
    testObj.CloseContainer();

    // Matrix 2 in File
    testObj.OpenContainer(true);
    std::string _datasetName2 = "Matrix 2";
    testObj.WriteMatrix(_datasetName2, matrixWrite2, matrixRows2, matrixColumns2);
    testObj.CloseContainer();

    // Matrix 2 in Group
    std::string groupName2 = "Group of Matrix 2";
    testObj.CreateContainer(false, groupName2);
    std::string _datasetName2G = "Matrix 2";
    testObj.WriteMatrix(_datasetName2G, matrixWrite2, matrixRows2, matrixColumns2,
                        EVAA::HDF5FileHandle::GROUP);
    testObj.CloseContainer();

    // Matrix 3
    std::string groupName3 = "Group of Matrix 3";
    testObj.CreateContainer(false, groupName3);
    std::string _datasetName3 = "Matrix 3";
    testObj.WriteMatrix(_datasetName3, matrixWrite3, matrixRows3, matrixColumns3,
                        EVAA::HDF5FileHandle::GROUP);
    testObj.CloseContainer();

    // Read part
    double matrixRead1[matrixRows1 * matrixColumns1];
    double matrixRead2[matrixRows2 * matrixColumns2];
    double matrixRead2G[matrixRows2 * matrixColumns2];
    double matrixRead3[matrixRows3 * matrixColumns3];
    size_t matrixRowsRead1 = 0, matrixColumnsRead1 = 0;
    size_t matrixRowsRead2 = 0, matrixColumnsRead2 = 0;
    size_t matrixRowsRead2G = 0, matrixColumnsRead2G = 0;
    size_t matrixRowsRead3 = 0, matrixColumnsRead3 = 0;
    testObj.OpenContainer(false);
    testObj.ReadMatrix(_datasetName1, matrixRead1, matrixRowsRead1, matrixColumnsRead1);
    testObj.ReadMatrix(_datasetName2, matrixRead2, matrixRowsRead2, matrixColumnsRead2);
    testObj.CloseContainer();

    testObj.OpenContainer(false, groupName2);
    testObj.ReadMatrix(_datasetName2G, matrixRead2G, matrixRowsRead2G, matrixColumnsRead2G,
                       EVAA::HDF5FileHandle::GROUP);
    testObj.CloseContainer();

    testObj.OpenContainer(false, groupName3);
    testObj.ReadMatrix(_datasetName3, matrixRead3, matrixRowsRead3, matrixColumnsRead3,
                       EVAA::HDF5FileHandle::GROUP);
    testObj.CloseContainer();

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

TEST_F(TestOutputHDF5, WriteReadMatrixFileGroup) {
    const int matrixRows = 3;
    const int matrixColumns = 2;
    double matrixWrite1[matrixRows * matrixColumns]{5, 9, 91, 26, 5, 94};
    double matrixWrite2[matrixRows * matrixColumns]{1, 2, 1, 2, 1, 2};

    // Create OutputHDF5 object
    EVAA::HDF5::OutputHDF5<double> testObj(_filePath, _fileName);

    // Write part
    std::string groupName = "Group Test";
    testObj.CreateContainer(true, groupName);
    std::string _datasetName = "SolutionMatrix";
    testObj.WriteMatrix(_datasetName, matrixWrite1, matrixRows, matrixColumns);
    testObj.CloseContainer();

    testObj.OpenContainer(true, groupName);
    testObj.WriteMatrix(_datasetName, matrixWrite2, matrixRows, matrixColumns,
                        EVAA::HDF5FileHandle::GROUP);
    testObj.CloseContainer();

    // Read part
    double matrixRead1[matrixRows * matrixColumns];
    double matrixRead2[matrixRows * matrixColumns];
    size_t matrixRowsRead = 0;
    size_t matrixColumnsRead = 0;
    testObj.OpenContainer(false, groupName);
    testObj.ReadMatrix(_datasetName, matrixRead1, matrixRowsRead, matrixColumnsRead);
    testObj.ReadMatrix(_datasetName, matrixRead2, matrixRowsRead, matrixColumnsRead,
                       EVAA::HDF5FileHandle::GROUP);
    testObj.CloseContainer();

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

TEST_F(TestOutputHDF5, WriteString) {
    // Create OutputHDF5 object
    EVAA::HDF5::OutputHDF5<double> testObj(_filePath, _fileName);

    std::string dataWrite = "Data to write";
    const std::string datasetName = "Write String";
    const std::string groupName = "Group string";

    testObj.CreateContainer(true, groupName);
    testObj.WriteString(datasetName, dataWrite, EVAA::HDF5FileHandle::FILE);
    testObj.CloseContainer();
}

TEST_F(TestOutputHDF5, ReadWrittenString) {
    // Create OutputHDF5 object
    EVAA::HDF5::OutputHDF5<double> testObj(_filePath, _fileName);

    std::string dataRead;
    std::string dataWrite = "Data to write";
    const std::string datasetName = "Write String";
    const std::string groupName = "Group string";

    testObj.CreateContainer(true, groupName);
    testObj.WriteString(datasetName, dataWrite, EVAA::HDF5FileHandle::FILE);
    testObj.CloseContainer();

    testObj.OpenContainer(false, groupName);
    testObj.ReadString(datasetName, dataRead, EVAA::HDF5FileHandle::FILE);
    testObj.CloseContainer();

    EXPECT_EQ(dataRead, dataWrite) << "Read data differs from written!";
}

TEST_F(TestOutputHDF5, WriteMatrixVectorString) {
    EVAA::HDF5::OutputHDF5<double> testObj(_filePath, "write 3" + _fileName);

    // write vector
    const size_t vecDim = 6;
    double vecWrite[vecDim]{5, 9, 91, 26, 5, 94};
    std::string datasetNameVector = "SolutionVector";
    std::string groupName = "Group Test";
    testObj.CreateContainer(true, groupName);
    testObj.WriteVector(datasetNameVector, vecWrite, vecDim);
    testObj.WriteVector(datasetNameVector, vecWrite, vecDim, EVAA::HDF5FileHandle::GROUP);
    testObj.CloseContainer();

    // write matrix
    const int matrixRows = 3;
    const int matrixColumns = 2;
    double matrixWrite[matrixRows * matrixColumns]{5, 9, 91, 26, 5, 94};
    std::string datasetNameMatrix = "SolutionMatrix";
    testObj.OpenContainer(true, groupName);
    testObj.WriteMatrix(datasetNameMatrix, matrixWrite, matrixRows, matrixColumns);
    testObj.WriteMatrix(datasetNameMatrix, matrixWrite, matrixRows, matrixColumns,
                        EVAA::HDF5FileHandle::GROUP);
    testObj.CloseContainer();

    // write string
    std::string data = "Data to write";
    const std::string datasetNameString = "Write String";

    testObj.OpenContainer(true, groupName);
    testObj.WriteString(datasetNameString, data);
    testObj.WriteString(datasetNameString, data, EVAA::HDF5FileHandle::GROUP);
    testObj.CloseContainer();
}

}  // namespace Test
#endif
