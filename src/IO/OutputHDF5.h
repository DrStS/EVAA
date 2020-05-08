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
void W11riteVectorToFile(const std::string& fileName, const std::string& vectorName, T* vector,
                         size_t size);

/**
 * Reads a single vector from a HDF5 file.
 * \param[in] fileName The HDF5 file to read from.
 * \param[in] vectorName The vector name in the file.
 * \param[out] vector The vector to read.
 * \param[in] size The vector size.
 */
template <typename T>
void ReadVectorFromFile(const std::string& fileName, const std::string& vectorName, T* vector,
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
void WriteMatrixToFile(const std::string& fileName, const std::string& matrixName, T* matrix,
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
void ReadMatrixFromFile(const std::string& fileName, const std::string& matrixName, T* matrix,
                        size_t numRows, size_t numColumns);

/** \brief Templates for HDF5DataType (see https://stackoverflow.com/a/15221077)
 * Allows writing float and doubles. It can be also extended for other data types.
 */
template <typename T>
struct GetHDF5DataType {
    static H5::PredType type() { return H5::PredType::NATIVE_DOUBLE; }
};
template <>
struct GetHDF5DataType<float> {
    H5::FloatType type{H5::PredType::NATIVE_FLOAT};
};
template <>
struct GetHDF5DataType<double> {
    H5::FloatType type{H5::PredType::NATIVE_DOUBLE};
};

/** \brief OutputHDF5 is a class that handles the writing to and reading from HDF5 files.
 * The algorithms instantiate an OutputHDF5 object.
 */
template <class T>
class OutputHDF5 {
public:
    /** Constructor with parameters for OutputHDF5. Initialize _filePath and _fileName and set
     * nullptr to handles
     * \param[in] filePath Directory path where the files will be written. TODO: implement it
     * properly (create new folder with check)
     * \param[in] fileName Name of the file where the data is written
     */
    OutputHDF5(std::string filePath, std::string fileName) : _filePath(""), _fileName(fileName) {
        // TODO Use _filePath argument

        H5::Exception::dontPrint();
        _fileHandle = nullptr;
        _groupHandle = nullptr;
    };

    /** Destructor of OutputHDF5, deletes the handles. */
    ~OutputHDF5() {
        delete _fileHandle;
        delete _groupHandle;
    }

    /** \brief Function for creating a container. It allocates the handlers for the HDF5 File and
     * HDF5 Group
     * \param[in] _writePermission (true or false).
     *           - True is used for overwriting / creating a new file with a new container.
     *           - False is used for appending a different in a already defined file.
     * \param[groupName] defines the group name of the container. Default name is "/Default"
     */
    void CreateContainer(bool _forceWrite, std::string groupName = "/Default") {
        try {
            if (_forceWrite) {
                _fileHandle = new H5::H5File(_filePath + _fileName, H5F_ACC_TRUNC);
                _groupHandle = new H5::Group(_fileHandle->createGroup(groupName));
            }
            else {
                _fileHandle = new H5::H5File(_filePath + _fileName, H5F_ACC_RDWR);
                _groupHandle = new H5::Group(_fileHandle->createGroup(groupName));
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
    void OpenContainer(bool _writePermission, std::string groupName = "/Default") {
        try {
            if (_writePermission) {
                _fileHandle = new H5::H5File(_filePath + _fileName, H5F_ACC_RDWR);
                _groupHandle = new H5::Group(_fileHandle->openGroup(groupName));
            }
            else {
                _fileHandle = new H5::H5File(_filePath + _fileName, H5F_ACC_RDONLY);
                _groupHandle = new H5::Group(_fileHandle->openGroup(groupName));
            }
        }
        catch (H5::FileIException error) {
            std::cout << "Error: Cannot open file" << std::endl;
            std::cout << "File already exists!" << std::endl;
        }
    }

    /** \brief Close the container.
     */
    void CloseContainer(void) {
        _groupHandle->close();
        _fileHandle->close();
    }

    /** Functions for Matrices */

    /* \brief Write a Matrix in the container, in a group or outside of a group.
     * \param[in] datasetName The name of the dataset containing the data.
     * \param[in] data The matrix values.
     * \param[in] numRows The matrix number of rows.
     * \param[in] numColumns The matrix number of columns.
     * \param[in] handle Specifier for using a FILE (Default) or GROUP handle. If GROUP, then the
     * object use a (pre)defined group from createContainer or openContainer
     */
    void WriteMatrix(const std::string& datasetName, const T* data, const size_t numRows,
                     const size_t numColumns, const HDF5FileHandle handle = HDF5FileHandle::FILE) {
        const GetHDF5DataType<T> HDF5DataType;
        if (handle == HDF5FileHandle::FILE)
            WriteMatrix(datasetName, data, numRows, numColumns, HDF5DataType.type, _fileHandle);
        else if (handle == HDF5FileHandle::GROUP)
            WriteMatrix(datasetName, data, numRows, numColumns, HDF5DataType.type, _groupHandle);
        else {
            throw std::logic_error("Wrong handle provided. Must be one of: FILE, GROUP");
        }
    }

    /* \brief Reads a Matrix from the container, in a group or outside of a group.
     * \param[in] datasetName The name of the dataset containing the data.
     * \param[out] data The matrix values.
     * \param[out] numRows The matrix number of rows.
     * \param[out] numColumns The matrix number of columns.
     * \param[in] handle Specifier for using a FILE (Default) or GROUP handle. If GROUP, then the
     * object use a (pre)defined group from createContainer or openContainer
     */
    void ReadMatrix(const std::string& datasetName, T* data, size_t& numRows, size_t& numColumns,
                    const HDF5FileHandle handle = HDF5FileHandle::FILE) {
        const GetHDF5DataType<T> HDF5DataType;
        if (handle == HDF5FileHandle::FILE)
            ReadMatrix(datasetName, data, numRows, numColumns, HDF5DataType.type, _fileHandle);
        else if (handle == HDF5FileHandle::GROUP)
            ReadMatrix(datasetName, data, numRows, numColumns, HDF5DataType.type, _groupHandle);
        else {
            throw std::logic_error("Wrong handle provided. Must be one of: FILE, GROUP");
        }
    }

    /** Functions for Vectors */

    /* \brief Write a Vector in the container, in a group or outside of a group.
     * \param[in] datasetName The name of the dataset containing the data.
     * \param[in] data The vector values.
     * \param[in] numElements The vector length.
     * \param[in] handle Specifier for using a FILE (Default) or GROUP handle. If GROUP, then the
     * object use a (pre)defined group from createContainer or openContainer
     */
    void WriteVector(const std::string& datasetName, const T* data, const size_t numElements,
                     const HDF5FileHandle handle = HDF5FileHandle::FILE) {
        const GetHDF5DataType<T> HDF5DataType;
        if (handle == HDF5FileHandle::FILE)
            WriteVector(datasetName, data, numElements, HDF5DataType.type, _fileHandle);
        else if (handle == HDF5FileHandle::GROUP)
            WriteVector(datasetName, data, numElements, HDF5DataType.type, _groupHandle);
        else {
            throw std::logic_error(
                "Wrong handle provided. Must be one of: HDF5FileHandle::FILE, "
                "HDF5FileHandle::GROUP");
        }
    }

    /* \brief Reads a Vector from the container, in a group or outside of a group.
     * \param[in] datasetName The name of the dataset containing the data.
     * \param[out] data The vector values.
     * \param[out] numElements The vector length.
     * \param[in] handle Specifier for using a FILE (Default) or GROUP handle. If GROUP, then the
     * object use a (pre)defined group from createContainer or openContainer
     */
    void ReadVector(const std::string& datasetName, T* data, size_t& numElements,
                    const HDF5FileHandle handle = HDF5FileHandle::FILE) {
        const GetHDF5DataType<T> HDF5DataType;
        if (handle == HDF5FileHandle::FILE)
            ReadVector(datasetName, data, numElements, HDF5DataType.type, _fileHandle);
        else if (handle == HDF5FileHandle::GROUP)
            ReadVector(datasetName, data, numElements, HDF5DataType.type, _groupHandle);
        else {
            throw std::logic_error("Wrong handle provided. Must be one of: FILE, GROUP");
        }
    }

    /** Functions for Strings */
    void WriteString(const std::string& datasetName, const std::string data,
                     const HDF5FileHandle handle = HDF5FileHandle::FILE) {
        // allocate handle
        H5::Group* handleWrite;
        if (handle == HDF5FileHandle::FILE)
            handleWrite = _fileHandle;
        else if (handle == HDF5FileHandle::GROUP)
            handleWrite = _groupHandle;
        else {
            throw std::logic_error(
                "Wrong handle provided. Must be one of: HDF5FileHandle::FILE, "
                "HDF5FileHandle::GROUP");
        }

        // write the string in dataset

        // Create new dataspace for attribute
        H5::DataSpace dataspace(H5S_SCALAR);

        // Create new string datatype for attribute
        H5::StrType datatype(H5::PredType::C_S1, data.length());

        /*
         * Create a new dataset within the file using defined dataspace and
         * datatype and default dataset creation properties.
         */
        H5::DataSet* dataset =
            new H5::DataSet(handleWrite->createDataSet(datasetName, datatype, dataspace));

        dataset->write(data, datatype);

        delete dataset;
    }

    void ReadString(const std::string& datasetName, std::string& data,
                    const HDF5FileHandle handle = HDF5FileHandle::FILE) {
        // allocate handle
        H5::Group* handleRead;
        if (handle == HDF5FileHandle::FILE)
            handleRead = _fileHandle;
        else if (handle == HDF5FileHandle::GROUP)
            handleRead = _groupHandle;
        else {
            throw std::logic_error(
                "Wrong handle provided. Must be one of: HDF5FileHandle::FILE, "
                "HDF5FileHandle::GROUP");
        }

        /* Open a dataset */
        H5::DataSet dataset = handleRead->openDataSet(datasetName);

        /* Get datatype for dataset */
        H5::DataType dtype = dataset.getDataType();

        /* Get length of string */
        size_t dataLength = dataset.getDataType().getSize();

        /* Read dataset from disk */
        char* readData = new char[dataLength];
        dataset.read((void*)readData, dtype);

        /* Assign data values to data string */
        data.assign(readData, dataLength);

        delete[] readData;
    }

private:
    /** my file name; */
    std::string _fileName;
    /** my file path; */
    std::string _filePath;
    /** my file handle; */
    H5::H5File* _fileHandle;
    /** my group handle; */
    H5::Group* _groupHandle;

private:
    /** Functions for Matrices */

    /* \brief Creates a dataset in a container and write a Matrix in a group or outside of a group.
     * Templated over datatype.
     * \param[in] datasetName The name of the dataset containing the data.
     * \param[in] data The matrix values.
     * \param[in] numRows The matrix number of rows.
     * \param[in] numColumns The matrix number of columns.
     * \param[in] datatype The type of data to be written, used for double / float.
     * \param[in] handle Can be of type H5File* or H5Group*.
     */
    template <typename HDF5DataType>
    void WriteMatrix(const std::string& datasetName, const T* data, const size_t numRows,
                     const size_t numColumns, const HDF5DataType& datatype,
                     const H5::Group* handle) {
        /*
         * Define the size of the array and create the data space for fixed
         * size dataset.
         */
        const size_t size[] = {numRows, numColumns};
        H5::DataSpace dataspace(2, size);

        /*
         * Create a new dataset within the file using defined dataspace and
         * datatype and default dataset creation properties.
         */
        H5::DataSet* dataset =
            new H5::DataSet(handle->createDataSet(datasetName, datatype, dataspace));

        dataset->write(data, datatype);

        delete dataset;

        /*
        // if need to use attributes
        // Add matrix dimension information to container as attributes
        H5::DataSpace attrDataspaceScalar(H5S_SCALAR);

        H5::Attribute attribute = dataset->createAttribute(
            "Matrix # columns", H5::PredType::STD_I32BE, attrDataspaceScalar);
        int attrDataScalar[1] = {numColumns};
        attribute.write(H5::PredType::NATIVE_INT, attrDataScalar);

        attribute =
            dataset->createAttribute("Matrix # rows", H5::PredType::STD_I32BE, attrDataspaceScalar);
        attrDataScalar[0] = numRows;
        attribute.write(H5::PredType::NATIVE_INT, attrDataScalar);
        */
    }

    /* \brief Reads a Matrix from a dataset from a container, from a group or outside of a group.
     * Templated over datatype.
     * \param[in] datasetName The name of the dataset containing the data.
     * \param[out] data The matrix values.
     * \param[out] numRows The matrix number of rows.
     * \param[out] numColumns The matrix number of columns.
     * \param[in] datatype The type of data to be written, used for double / float.
     * \param[in] handle Can be of type H5File* or H5Group*.
     */
    template <typename HDF5DataType>
    void ReadMatrix(const std::string& datasetName, T* data, size_t& numRows, size_t& numColumns,
                    const HDF5DataType& datatype, const H5::Group* handle) {
        // Open dataset
        H5::DataSet dataset = handle->openDataSet(datasetName);

        // Read Matrix
        dataset.read(data, datatype);

        size_t dims[2];
        // use this for simple data space
        H5Sget_simple_extent_dims(H5Dget_space(dataset.getId()), dims, NULL);

        numRows = dims[0];
        numColumns = dims[1];

        /* // if need to use attributes
        H5::Attribute attribute = dataset.openAttribute("Matrix # columns");
        attribute.read(H5::PredType::NATIVE_INT, &numColumns);
        */
    }

    /** Functions for Vectors */

    /* \brief Creates a dataset in a container and write a Vector in a group or outside of a group.
     * Templated over datatype.
     * \param[in] datasetName The name of the dataset containing the data.
     * \param[in] data The vector values.
     * \param[in] numElements The vector number of elements.
     * \param[in] datatype The type of data to be written, used for double / float.
     * \param[in] handle Can be of type H5File* or H5Group*.
     */
    template <typename HDF5DataType>
    void WriteVector(const std::string& datasetName, const T* data, const size_t numElements,
                     const HDF5DataType& datatype, const H5::Group* handle) {
        /*
         * Define the size of the array and create the data space for fixed
         * size dataset.
         */
        const size_t size[] = {numElements};
        H5::DataSpace dataspace(1, size);

        /*
         * Create a new dataset within the file using defined dataspace and
         * datatype and default dataset creation properties.
         */
        H5::DataSet* dataset =
            new H5::DataSet(handle->createDataSet(datasetName, datatype, dataspace));

        dataset->write(data, datatype);

        /*
        // Add vector dimension information to container (as attributes)
        H5::DataSpace attrDataspaceScalar(H5S_SCALAR);

        H5::Attribute attribute = dataset->createAttribute(
            "Vector # elements", H5::PredType::STD_I32BE, attrDataspaceScalar);
        int attrDataScalar[1] = {numElements};
        attribute.write(H5::PredType::NATIVE_INT, attrDataScalar);
        */

        delete dataset;
    }

    /* \brief Reads a Vector from a dataset from a container, from a group or outside of a group.
     * Templated over datatype.
     * \param[in] datasetName The name of the dataset containing the data.
     * \param[out] data The vector values.
     * \param[out] numElements The vector number of elements.
     * \param[in] datatype The type of data to be written, used for double / float.
     * \param[in] handle Can be of type H5File* or H5Group*.
     */
    template <typename HDF5DataType>
    void ReadVector(const std::string& datasetName, T* data, size_t& numElements,
                    const HDF5DataType& datatype, const H5::Group* handle) {
        /*
         * Create a new dataset within the file using defined dataspace and
         * datatype and default dataset creation properties.
         */
        H5::DataSet dataset = handle->openDataSet(datasetName);

        dataset.read(data, datatype);

        // use this for simple 1D data space
        numElements = dataset.getDataType().getSize();
    }
};
}  // namespace HDF5
};  // namespace EVAA

#ifdef USE_HDF5
#include "OutputHDF5.cpp"
#endif  // USE_HDF5

#endif  // USE_HDF5
