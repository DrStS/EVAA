// TODO: Copyright header
#pragma once

#ifdef USE_HDF5
#include <cstdint>
#include <string>

#include "H5Cpp.h"
#include "H5FloatType.h"
#include "H5PredType.h"

namespace EVAA {
namespace HDF5 {
/**
 * Writes a single vector into a HDF5 file.
 * \param[in] fileName The HDF5 file to write to.
 * \param[in] vectorName The vector name in the file.
 * \param[in] vector The vector to write.
 * \param[in] size The vector size.
 */
template <typename T>
void writeVectorToFile(const std::string& fileName, const std::string& vectorName, T* vector,
                       size_t size);

/**
 * Reads a single vector from a HDF5 file.
 * \param[in] fileName The HDF5 file to read from.
 * \param[in] vectorName The vector name in the file.
 * \param[out] vector The vector to read.
 * \param[in] size The vector size.
 */
template <typename T>
void readVectorFromFile(const std::string& fileName, const std::string& vectorName, T* vector,
                        size_t size);

/**
 * Writes a single matrix to a HDF5 file.
 * \param[in] fileName The HDF5 file to write to.
 * \param[in] matrixName The matrix name in the file.
 * \param[in] matrix The matrix to write.
 * \param[in] numRows The matrix number of rows.
 * \param[in] numColumns The matrix number of columns.
 */
template <typename T>
void writeMatrixToFile(const std::string& fileName, const std::string& matrixName, T* matrix,
                       size_t numRows, size_t numColumns);

/**
 * Reads a single matrix from a HDF5 file.
 * \param[in] fileName The HDF5 file to read from.
 * \param[in] matrixName The matrix name in the file.
 * \param[out] matrix The vector to read.
 * \param[in] numRows The matrix number of rows.
 * \param[in] numColumns The matrix number of columns.
 */
template <typename T>
void readMatrixFromFile(const std::string& fileName, const std::string& matrixName, T* matrix,
                        size_t numRows, size_t numColumns);

// Templates for HDF5DataType
// see https://stackoverflow.com/a/15221077
template <typename T>
struct getHDF5DataType {
    static H5::PredType type() { return H5::PredType::NATIVE_DOUBLE; }
};
template <>
struct getHDF5DataType<float> {
    H5::FloatType type{H5::PredType::NATIVE_FLOAT};
};
template <>
struct getHDF5DataType<double> {
    H5::FloatType type{H5::PredType::NATIVE_DOUBLE};
};

template <class T>
class OutputHDF5 {
public:
    OutputHDF5(std::string _fileName, std::string _filePath) :
        myFileName(_fileName), myFilePath(_filePath) {
#ifdef USE_HDF5
        H5::Exception::dontPrint();
        myHDF5FileHandle = nullptr;
        myHDF5GroupOperators = nullptr;
#endif
    };
    /* Comment for Read part
    STACCATOComplexDouble* kdynRead = new STACCATOComplexDouble[LENGTH];
    dataset->read(kdynRead, mtype);
    std::cout << std::endl << "Field real : " << std::endl;
    for (i = 0; i < LENGTH; i++) std::cout << kdyn[i].real - kdynRead[i].real << " ";
    std::cout << std::endl;
    std::cout << std::endl << "Field imag : " << std::endl;
    for (i = 0; i < LENGTH; i++) std::cout << kdyn[i].imag - kdynRead[i].imag << " ";
    std::cout << std::endl;
    delete dataset;
    */

    ~OutputHDF5() {
        delete myHDF5FileHandle;
        delete myHDF5GroupOperators;
    }

    void createContainer(bool _forceWrite) {
#ifdef USE_HDF5

        try {
            if (_forceWrite) {
                myHDF5FileHandle = new H5::H5File(myFilePath + myFileName, H5F_ACC_TRUNC);
                myHDF5GroupOperators = new H5::Group(myHDF5FileHandle->createGroup("/EVAAOutput"));
            }
            else {
                myHDF5FileHandle = new H5::H5File(myFilePath + myFileName, H5F_ACC_EXCL);
                myHDF5GroupOperators = new H5::Group(myHDF5FileHandle->createGroup("/EVAAOutput"));
            }
        }
        catch (H5::FileIException error) {
            std::cout << "Error: Cannot create file" << std::endl;
            std::cout << "File already exists!" << std::endl;
        }
#endif  // USE_HDF5
    }

    void closeContainer(void) {
#ifdef USE_HDF5
        //myHDF5GroupOperators->close();
        myHDF5FileHandle->close();
#endif  // USE_HDF5
    }

    void openContainer(bool _writePermission) {
#ifdef USE_HDF5
        try {
            if (_writePermission) {
                myHDF5FileHandle = new H5::H5File(myFilePath + myFileName, H5F_ACC_RDWR);
            }
            else {
                myHDF5FileHandle = new H5::H5File(myFilePath + myFileName, H5F_ACC_RDONLY);
            }
        }
        catch (H5::FileIException error) {
            std::cout << "Error: Cannot open file" << std::endl;
            std::cout << "File already exists!" << std::endl;
        }
#endif  // USE_HDF5
    }

    /* \brief Templated over datatype, writes a Matrix in the container, together with its
     * dimensions
     */
    template <typename HDF5DataType>
    void WriteMatrix(const std::string& _datasetName, const T* _data, const size_t _numRows,
                     const size_t _numColumns, const HDF5DataType& datatype) {
        /*
         * Define the size of the array and create the data space for fixed
         * size dataset.
         */
        const size_t size[] = {_numRows, _numColumns};
        H5::DataSpace dataspace(2, size);

        /*
         * Create a new dataset within the file using defined dataspace and
         * datatype and default dataset creation properties.
         */
        H5::DataSet* dataset =
            new H5::DataSet(myHDF5FileHandle->createDataSet(_datasetName, datatype, dataspace));

        dataset->write(_data, datatype);

        // Add matrix dimension information to container
        H5::DataSpace attrDataspaceScalar(H5S_SCALAR);

        H5::Attribute attribute = dataset->createAttribute(
            "Matrix # columns", H5::PredType::STD_I32BE, attrDataspaceScalar);
        int attrDataScalar[1] = {_numColumns};
        attribute.write(H5::PredType::NATIVE_INT, attrDataScalar);

        attribute =
            dataset->createAttribute("Matrix # rows", H5::PredType::STD_I32BE, attrDataspaceScalar);
        attrDataScalar[0] = _numRows;
        attribute.write(H5::PredType::NATIVE_INT, attrDataScalar);

        delete dataset;
    }

    /* \brief Write a Matrix in the container, together with its dimensions
     *
     */
    void WriteMatrix(const std::string& _datasetName, const T* _data, const size_t _numRows,
                     const size_t _numColumns) {
        const getHDF5DataType<T> HDF5DataType;
        WriteMatrix(_datasetName, _data, _numRows, _numColumns, HDF5DataType.type);
    }

    /* \brief Templated over datatype, writes a Matrix in the container, together with its
     * dimensions
     */
    template <typename HDF5DataType>
    void ReadMatrix(const std::string& _datasetName, T* _data, size_t& _numRows,
                    size_t& _numColumns,
                    const HDF5DataType& datatype) {
        
        // Open dataset
        H5::DataSet dataset = myHDF5FileHandle->openDataSet(_datasetName);

        // Read Matrix
        dataset.read(_data, datatype);

        H5::Attribute attribute = dataset.openAttribute("Matrix # columns");
        attribute.read(H5::PredType::NATIVE_INT, &_numColumns);

        attribute = dataset.openAttribute("Matrix # rows");
        attribute.read(H5::PredType::NATIVE_INT, &_numRows);
    }

    /* \brief Write a Matrix in the container, together with its dimensions
     *
     */
    void ReadMatrix(const std::string& _datasetName, T* _data, size_t& _numRows,
                    size_t& _numColumns) {
        const getHDF5DataType<T> HDF5DataType;
        ReadMatrix(_datasetName, _data, _numRows, _numColumns, HDF5DataType.type);
    }

    /* \brief Templated over datatype, writes a Vector in the container, together with its
     * dimensions
     */
    template <typename HDF5DataType>
    void WriteVector(const std::string& _datasetName, const T* _data, const size_t _numElements,
                     const HDF5DataType& datatype) {
        /*
         * Define the size of the array and create the data space for fixed
         * size dataset.
         */
        const size_t size[] = {_numElements};
        H5::DataSpace dataspace(1, size);

        /*
         * Create a new dataset within the file using defined dataspace and
         * datatype and default dataset creation properties.
         */
        H5::DataSet* dataset =
            new H5::DataSet(myHDF5FileHandle->createDataSet(_datasetName, datatype, dataspace));

        dataset->write(_data, datatype);

        // Add vector dimension information to container
        H5::DataSpace attrDataspaceScalar(H5S_SCALAR);

        H5::Attribute attribute = dataset->createAttribute(
            "Vector # elements", H5::PredType::STD_I32BE, attrDataspaceScalar);
        int attrDataScalar[1] = {_numElements};
        attribute.write(H5::PredType::NATIVE_INT, attrDataScalar);

        delete dataset;
    }

    /* \brief Write a Vector in the container, together with its dimensions
     *
     */
    void WriteVector(const std::string& _datasetName, const T* _data, const size_t _numElements) {
        const getHDF5DataType<T> HDF5DataType;
        WriteVector(_datasetName, _data, _numElements, HDF5DataType.type);
    }

    /* \brief Templated over datatype, Reads a Vector in the container, together with its
     * dimensions
     */
    template <typename HDF5DataType>
    void ReadVector(const std::string& _datasetName, T* _data, size_t& _numElements,
                    const HDF5DataType& datatype) {
        /*
         * Create a new dataset within the file using defined dataspace and
         * datatype and default dataset creation properties.
         */
        H5::DataSet dataset = myHDF5FileHandle->openDataSet(_datasetName);

        dataset.read(_data, datatype);

        H5::Attribute attribute = dataset.openAttribute("Vector # elements");
        attribute.read(H5::PredType::NATIVE_INT, &_numElements);
    }

    /* \brief Reads a Vector in the container, together with its dimensions
     *
     */
    void ReadVector(const std::string& _datasetName, T* _data, size_t& _numElements) {
        const getHDF5DataType<T> HDF5DataType;
        ReadVector(_datasetName, _data, _numElements, HDF5DataType.type);
    }

public:
    /// my file name;
    std::string myFileName;
    /// my file path;
    std::string myFilePath;
    /// my file handle
    H5::H5File* myHDF5FileHandle;
    /// my group handle
    H5::Group* myHDF5GroupOperators;
};
}  // namespace HDF5
};  // namespace EVAA

#ifdef USE_HDF5
#include "OutputHDF5.cpp"
#endif  // USE_HDF5

#endif  // USE_HDF5
