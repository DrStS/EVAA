// TODO: Copyright header

#pragma once

#include <string>

/**
 * Writes a vector into a HDF5 file.
 * \param[in] fileName The HDF5 file to write to.
 * \param[in] vectorName The vector name in the file.
 * \param[in] vector The vector to write.
 * \param[in] size The vector size.
 */
template<typename T>
void writeVectorToFile(const std::string& fileName, const std::string& vectorName, T* vector, size_t size);

/**
 * Reads a vector from a HDF5 file.
 * \param[in] fileName The HDF5 file to read from.
 * \param[in] vectorName The vector name in the file.
 * \param[out] vector The vector to read.
 * \param[in] size The vector size.
 */
template<typename T>
void readVectorFromFile(const std::string& fileName, const std::string& vectorName, T* vector, size_t size);

#include "outputhdf5.cpp"
