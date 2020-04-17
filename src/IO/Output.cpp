#include "Output.h"
#include <fstream>
#include <iostream>
#include <limits>
#include <iomanip>
#include <string>

void IO::write_matrix(double* T, std::string fname, size_t rows, size_t cols) {
	std::ofstream myfile(fname);
	for (size_t i = 0; i < rows; i++)
	{
		for (size_t j = 0; j < cols; j++)
		{
			myfile << T[i*cols + j] << std::setprecision(15) << " ";
		}
		myfile << "\n";
	}

}
