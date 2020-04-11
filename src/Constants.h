#pragma once

namespace Constants {
	constexpr int ALIGNMENT = 64;
	constexpr int DIM = 3; // consider Constants::DIM = 4 for efficiency!!! // 10 Constants::DIMensions because of torque of the body
	constexpr int VEC_DIM = 9; // [CG, Wi, Ti] (9 positions / vectors for CG, wheels and tyres)
	constexpr int INCX = 1; 
	constexpr int NUM_LEGS = 4;
};
