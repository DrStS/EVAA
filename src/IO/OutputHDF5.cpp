#include <string>

#ifdef USE_HDF5
#include "H5Cpp.h"
#include "H5FloatType.h"
#include "H5PredType.h"

template<typename T>
void writeVectorToFile(const std::string& fileName, const std::string& vectorName, T* vector, size_t size) {
	H5::H5File file(fileName, H5F_ACC_TRUNC);
	H5::DataSpace dataspace(1, &size);

	H5::FloatType datatype(H5::PredType::NATIVE_DOUBLE);
	datatype.setOrder(H5T_ORDER_LE);
}

template<typename T>
void readVectorFromFile(const std::string& fileName, const std::string& vectorName, T* vector, size_t size) {

}

template <typename T> void writeVectorToFile(const std::string&, const std::string&, float, size_t);
template <typename T> void writeVectorToFile(const std::string&, const std::string&, double, size_t);
#endif  // USE_HDF5