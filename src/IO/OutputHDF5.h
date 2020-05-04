// TODO: Copyright header

#pragma once

#include <string>

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

#ifdef USE_HDF5
#include "OutputHDF5.cpp"
#endif  // USE_HDF5
