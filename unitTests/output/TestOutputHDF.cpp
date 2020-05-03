#include "OutputHDF5.h"
#include "gtest/gtest.h"

#include <array>

class OutputTest : public ::testing::Test {};

#ifdef USE_HDF5
TEST_F(OutputTest, writeHDF5file) {
	const int vec_dim = 6;
	std::array<double, vec_dim> vec{ 5, 9, 91, 26, 5, 94 };
	std::array<double, vec_dim> vec_another;

	writeVectorToFile("output1.h5", "vec2write", vec.data(), vec_dim);
	readVectorFromFile("output1.h5", "vec2write", vec_another.data(), vec_dim);

	for (int i = 0; i < vec.size(); ++i) {
		EXPECT_EQ(vec[i], vec_another[i]) << "Vectors differ at index " << i;
	}
}
#endif