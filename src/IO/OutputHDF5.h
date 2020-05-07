// TODO: Copyright header
#pragma once

#ifdef USE_HDF5
#include <cstdint>
#include <string>

#include "H5Cpp.h"
#include "H5FloatType.h"
#include "H5PredType.h"

namespace EVAA {

/** Handle for HDF5 write in files or groups fashion. */
enum class HDF5FileHandle { FILE, GROUP };

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

/** \brief Templates for HDF5DataType (see https://stackoverflow.com/a/15221077)
 * Allows writing float and doubles. It can be also extended for other data types.
 */
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

/** \brief OutputHDF5 is a class that handles the writing to and reading from HDF5 files.
 * The algorithms instantiate an OutputHDF5 object.
 */
template <class T>
class OutputHDF5 {
public:
    OutputHDF5(std::string _fileName, std::string _filePath) :
        myFileName(_fileName), myFilePath(_filePath) {
        H5::Exception::dontPrint();
        myHDF5FileHandle = nullptr;
        myHDF5GroupOperators = nullptr;
    };

    ~OutputHDF5() {
        delete myHDF5FileHandle;
        delete myHDF5GroupOperators;
    }

    /** \brief Function for creating a container. It allocates the handlers for the HDF5 File and
     * HDF5 Group
     * \param[in] _writePermission (true or false).
     *           - True is used for overwriting / creating a new file with a new container.
     *           - False is used for appending a different in a already defined file.
     * \param[groupName] defines the group name of the container. Default name is "/Default"
     */
    void createContainer(bool _forceWrite, std::string groupName = "/Default") {
        try {
            if (_forceWrite) {
                myHDF5FileHandle = new H5::H5File(myFilePath + myFileName, H5F_ACC_TRUNC);
                myHDF5GroupOperators = new H5::Group(myHDF5FileHandle->createGroup(groupName));
            }
            else {
                myHDF5FileHandle = new H5::H5File(myFilePath + myFileName, H5F_ACC_RDWR);
                myHDF5GroupOperators = new H5::Group(myHDF5FileHandle->createGroup(groupName));
            }
        }
        catch (H5::FileIException error) {
            std::cout << "Error: Cannot create file" << std::endl;
            std::cout << "File already exists!" << std::endl;
        }
    }

    /** \brief Function for opening a container. It allocates the handlers for the HDF5 File and
     * HDF5 Group
     * \param[in] _writePermission (true or false).
     *           - True is used for writing in a already created container.
     *           - False is used during reading.
     * \param[groupName] defines the group name of the container. Default name is "/Default"
     */
    void openContainer(bool _writePermission, std::string groupName = "/Default") {
        try {
            if (_writePermission) {
                myHDF5FileHandle = new H5::H5File(myFilePath + myFileName, H5F_ACC_RDWR);
                myHDF5GroupOperators = new H5::Group(myHDF5FileHandle->openGroup(groupName));
            }
            else {
                myHDF5FileHandle = new H5::H5File(myFilePath + myFileName, H5F_ACC_RDONLY);
                myHDF5GroupOperators = new H5::Group(myHDF5FileHandle->openGroup(groupName));
            }
        }
        catch (H5::FileIException error) {
            std::cout << "Error: Cannot open file" << std::endl;
            std::cout << "File already exists!" << std::endl;
        }
    }

    void closeContainer(void) {
        myHDF5GroupOperators->close();
        myHDF5FileHandle->close();
    }

    /** Functions for Matrices */

    /* \brief Write a Matrix in the container, in a group or outside of a group.
     * \param[in] _datasetName The name of the dataset containing the data.
     * \param[in] _data The matrix values.
     * \param[in] _numRows The matrix number of rows.
     * \param[in] _numColumns The matrix number of columns.
     * \param[in] handle Specifier for using a FILE (Default) or GROUP handle. If GROUP, then the
     * object use a (pre)defined group from createContainer or openContainer
     */
    void WriteMatrix(const std::string& _datasetName, const T* _data, const size_t _numRows,
                     const size_t _numColumns, HDF5FileHandle handle = HDF5FileHandle::FILE) {
        const getHDF5DataType<T> HDF5DataType;
        if (handle == HDF5FileHandle::FILE)
            WriteMatrix(_datasetName, _data, _numRows, _numColumns, HDF5DataType.type,
                        myHDF5FileHandle);
        else if (handle == HDF5FileHandle::GROUP)
            WriteMatrix(_datasetName, _data, _numRows, _numColumns, HDF5DataType.type,
                        myHDF5GroupOperators);
        else {
            throw std::logic_error("Wrong handle provided. Must be one of: FILE, GROUP");
        }
    }

    /* \brief Reads a Matrix from the container, in a group or outside of a group.
     * \param[in] _datasetName The name of the dataset containing the data.
     * \param[out] _data The matrix values.
     * \param[out] _numRows The matrix number of rows.
     * \param[out] _numColumns The matrix number of columns.
     * \param[in] handle Specifier for using a FILE (Default) or GROUP handle. If GROUP, then the
     * object use a (pre)defined group from createContainer or openContainer
     */
    void ReadMatrix(const std::string& _datasetName, T* _data, size_t& _numRows,
                    size_t& _numColumns, HDF5FileHandle handle = HDF5FileHandle::FILE) {
        const getHDF5DataType<T> HDF5DataType;
        if (handle == HDF5FileHandle::FILE)
            ReadMatrix(_datasetName, _data, _numRows, _numColumns, HDF5DataType.type,
                       myHDF5FileHandle);
        else if (handle == HDF5FileHandle::GROUP)
            ReadMatrix(_datasetName, _data, _numRows, _numColumns, HDF5DataType.type,
                       myHDF5GroupOperators);
        else {
            throw std::logic_error("Wrong handle provided. Must be one of: FILE, GROUP");
        }
    }

    /** Functions for Vectors */

    /* \brief Write a Vector in the container, in a group or outside of a group.
     * \param[in] _datasetName The name of the dataset containing the data.
     * \param[in] _data The vector values.
     * \param[in] _numElements The vector length.
     * \param[in] handle Specifier for using a FILE (Default) or GROUP handle. If GROUP, then the
     * object use a (pre)defined group from createContainer or openContainer
     */
    void WriteVector(const std::string& _datasetName, const T* _data, const size_t _numElements,
                     HDF5FileHandle handle = HDF5FileHandle::FILE) {
        const getHDF5DataType<T> HDF5DataType;
        if (handle == HDF5FileHandle::FILE)
            WriteVector(_datasetName, _data, _numElements, HDF5DataType.type, myHDF5FileHandle);
        else if (handle == HDF5FileHandle::GROUP)
            WriteVector(_datasetName, _data, _numElements, HDF5DataType.type, myHDF5GroupOperators);
        else {
            throw std::logic_error("Wrong handle provided. Must be one of: FILE, GROUP");
        }
    }

    /* \brief Reads a Vector from the container, in a group or outside of a group.
     * \param[in] _datasetName The name of the dataset containing the data.
     * \param[out] _data The vector values.
     * \param[out] _numElements The vector length.
     * \param[in] handle Specifier for using a FILE (Default) or GROUP handle. If GROUP, then the
     * object use a (pre)defined group from createContainer or openContainer
     */
    void ReadVector(const std::string& _datasetName, T* _data, size_t& _numElements,
                    HDF5FileHandle handle = HDF5FileHandle::FILE) {
        const getHDF5DataType<T> HDF5DataType;
        if (handle == HDF5FileHandle::FILE)
            ReadVector(_datasetName, _data, _numElements, HDF5DataType.type, myHDF5FileHandle);
        else if (handle == HDF5FileHandle::GROUP)
            ReadVector(_datasetName, _data, _numElements, HDF5DataType.type, myHDF5GroupOperators);
        else {
            throw std::logic_error("Wrong handle provided. Must be one of: FILE, GROUP");
        }
    }

private:
    /** my file name; */
    std::string myFileName;
    /** my file path; */
    std::string myFilePath;
    /** my file handle; */
    H5::H5File* myHDF5FileHandle;
    /** my group handle; */
    H5::Group* myHDF5GroupOperators;

private:
    /** Functions for Matrices */

    /* \brief Creates a dataset in a container and write a Matrix in a group or outside of a group.
     * Templated over datatype.
     * \param[in] _datasetName The name of the dataset containing the data.
     * \param[in] _data The matrix values.
     * \param[in] _numRows The matrix number of rows.
     * \param[in] _numColumns The matrix number of columns.
     * \param[in] datatype The type of data to be written, used for double / float.
     * \param[in] handle Can be of type H5File* or H5Group*.
     */
    template <typename HDF5DataType>
    void WriteMatrix(const std::string& _datasetName, const T* _data, const size_t _numRows,
                     const size_t _numColumns, const HDF5DataType& datatype, H5::Group* handle) {
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
            new H5::DataSet(handle->createDataSet(_datasetName, datatype, dataspace));

        dataset->write(_data, datatype);

        delete dataset;

        /*
        // Add matrix dimension information to container as attributes
        H5::DataSpace attrDataspaceScalar(H5S_SCALAR);

        H5::Attribute attribute = dataset->createAttribute(
            "Matrix # columns", H5::PredType::STD_I32BE, attrDataspaceScalar);
        int attrDataScalar[1] = {_numColumns};
        attribute.write(H5::PredType::NATIVE_INT, attrDataScalar);

        attribute =
            dataset->createAttribute("Matrix # rows", H5::PredType::STD_I32BE, attrDataspaceScalar);
        attrDataScalar[0] = _numRows;
        attribute.write(H5::PredType::NATIVE_INT, attrDataScalar);
        */
    }

    /* \brief Reads a Matrix from a dataset from a container, from a group or outside of a group.
     * Templated over datatype.
     * \param[in] _datasetName The name of the dataset containing the data.
     * \param[out] _data The matrix values.
     * \param[out] _numRows The matrix number of rows.
     * \param[out] _numColumns The matrix number of columns.
     * \param[in] datatype The type of data to be written, used for double / float.
     * \param[in] handle Can be of type H5File* or H5Group*.
     */
    template <typename HDF5DataType>
    void ReadMatrix(const std::string& _datasetName, T* _data, size_t& _numRows,
                    size_t& _numColumns, const HDF5DataType& datatype, H5::Group* handle) {
        // Open dataset
        H5::DataSet dataset = handle->openDataSet(_datasetName);

        // Read Matrix
        dataset.read(_data, datatype);

        size_t dims[2];
        // use this for simple data space
        H5Sget_simple_extent_dims(H5Dget_space(dataset.getId()), dims, NULL);
        _numRows = dims[0];
        _numColumns = dims[1];

        /*
        H5::Attribute attribute = dataset.openAttribute("Matrix # columns");
        attribute.read(H5::PredType::NATIVE_INT, &_numColumns);

        attribute = dataset.openAttribute("Matrix # rows");
        attribute.read(H5::PredType::NATIVE_INT, &_numRows);
        */
    }

    /** Functions for Vectors */

    /* \brief Creates a dataset in a container and write a Vector in a group or outside of a group.
     * Templated over datatype.
     * \param[in] _datasetName The name of the dataset containing the data.
     * \param[in] _data The vector values.
     * \param[in] _numElements The vector number of elements.
     * \param[in] datatype The type of data to be written, used for double / float.
     * \param[in] handle Can be of type H5File* or H5Group*.
     */
    template <typename HDF5DataType>
    void WriteVector(const std::string& _datasetName, const T* _data, const size_t _numElements,
                     const HDF5DataType& datatype, H5::Group* handle) {
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
            new H5::DataSet(handle->createDataSet(_datasetName, datatype, dataspace));

        dataset->write(_data, datatype);

        /*
        // Add vector dimension information to container (as attributes)
        H5::DataSpace attrDataspaceScalar(H5S_SCALAR);

        H5::Attribute attribute = dataset->createAttribute(
            "Vector # elements", H5::PredType::STD_I32BE, attrDataspaceScalar);
        int attrDataScalar[1] = {_numElements};
        attribute.write(H5::PredType::NATIVE_INT, attrDataScalar);
        */

        delete dataset;
    }

    /* \brief Reads a Vector from a dataset from a container, from a group or outside of a group.
     * Templated over datatype.
     * \param[in] _datasetName The name of the dataset containing the data.
     * \param[out] _data The vector values.
     * \param[out] _numElements The vector number of elements.
     * \param[in] datatype The type of data to be written, used for double / float.
     * \param[in] handle Can be of type H5File* or H5Group*.
     */
    template <typename HDF5DataType>
    void ReadVector(const std::string& _datasetName, T* _data, size_t& _numElements,
                    const HDF5DataType& datatype, H5::Group* handle) {
        /*
         * Create a new dataset within the file using defined dataspace and
         * datatype and default dataset creation properties.
         */
        H5::DataSet dataset = handle->openDataSet(_datasetName);

        dataset.read(_data, datatype);

        // use this for simple data space
        H5Sget_simple_extent_dims(H5Dget_space(dataset.getId()), &_numElements, NULL);
    }
};
}  // namespace HDF5
};  // namespace EVAA

#ifdef USE_HDF5
#include "OutputHDF5.cpp"
#endif  // USE_HDF5

#endif  // USE_HDF5
