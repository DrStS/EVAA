#include <string>

#ifdef USE_HDF5
#include "H5Cpp.h"
#include "H5FloatType.h"
#include "H5PredType.h"

template <typename T>
void writeVectorToFile(const std::string& fileName, const std::string& datasetName, T* vector,
                       size_t size) {
    /*
     * Create a new file using H5F_ACC_TRUNC access,
     * default file creation properties, and default file
     * access properties.
     */
    H5::H5File file(fileName, H5F_ACC_TRUNC);

    /*
     * Define the size of the array and create the data space for fixed
     * size dataset.
     */
    H5::DataSpace dataspace(1, &size);

    /*
     * Define datatype for the data in the file.
     * We will store little endian INT numbers.
     */
    H5::FloatType datatype(H5::PredType::NATIVE_DOUBLE);
    datatype.setOrder(H5T_ORDER_LE);

    /*
     * Create a new dataset within the file using defined dataspace and
     * datatype and default dataset creation properties.
     */
    H5::DataSet dataset = file.createDataSet(datasetName, datatype, dataspace);
    /*
     * Write the data to the dataset using default memory space, file
     * space, and transfer properties.
     */
    dataset.write(vector, H5::PredType::NATIVE_DOUBLE);
}

template <typename T>
void readVectorFromFile(const std::string& fileName, const std::string& datasetName, T* vector,
                        size_t size) {
    /*
     * Open the specified file and the specified dataset in the file.
     */
    H5::H5File file(fileName, H5F_ACC_RDONLY);
    H5::DataSet dataset = file.openDataSet(datasetName);
    /*
     * Read data from hyperslab in the file into the hyperslab in
     * memory and display the data.
     */
    dataset.read(vector, H5::PredType::NATIVE_DOUBLE);
}

template <typename T>
void writeMatrixToFile(const std::string& fileName, const std::string& matrixName, T* matrix,
                       size_t numRows, size_t numColumns) {
    /*
     * Create a new file using H5F_ACC_TRUNC access,
     * default file creation properties, and default file
     * access properties.
     */
    H5::H5File file(fileName, H5F_ACC_TRUNC);

    /*
     * Define the size of the array and create the data space for fixed
     * size dataset.
     */
    const int rankMatrix = 2;
    const size_t size[rankMatrix]{numRows, numColumns};
    H5::DataSpace dataspace(rankMatrix, size);

    /*
     * Define datatype for the data in the file.
     * We will store little endian INT numbers.
     */
    H5::FloatType datatype(H5::PredType::NATIVE_DOUBLE);
    datatype.setOrder(H5T_ORDER_LE);

    /*
     * Create a new dataset within the file using defined dataspace and
     * datatype and default dataset creation properties.
     */
    H5::DataSet dataset = file.createDataSet(matrixName, datatype, dataspace);
    /*
     * Write the data to the dataset using default memory space, file
     * space, and transfer properties.
     */
    dataset.write(matrix, H5::PredType::NATIVE_DOUBLE);
}

template <typename T>
void readMatrixFromFile(const std::string& fileName, const std::string& matrixName, T* matrix,
                        size_t numRows, size_t numColumns) {
    /*
     * Open the specified file and the specified dataset in the file.
     */
    H5::H5File file(fileName, H5F_ACC_RDONLY);
    H5::DataSet dataset = file.openDataSet(matrixName);
    /*
     * Read data from hyperslab in the file into the hyperslab in
     * memory and display the data.
     */
    dataset.read(matrix, H5::PredType::NATIVE_DOUBLE);
}

#endif  // USE_HDF5
