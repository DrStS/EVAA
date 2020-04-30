// TODO: Copyright header

#pragma once

namespace EVAA {
namespace Constants {
/** alignment for mkl malloc */
constexpr int ALIGNMENT = 64;

// consider Constants::DIM = 4 for efficiency!!!
// 10 Constants::DIMensions because of torque of the body
/** dimension of solution space */
constexpr int DIM = 3;

/** [CG, Wi, Ti] (9 positions / vectors for CG, wheels and tyres) */
constexpr int VEC_DIM = 9;

/** MKL constant incx */
constexpr int INCX = 1;

/** number of tyres, wheels, legs */
constexpr int NUM_LEGS = 4;

/** Degrees of Freedom in the Eulerian Frame */
constexpr int DOF = 11;

// TODO: Get rid of it.
using floatEVAA = double;

constexpr bool USEINTERPOLATION = true;

constexpr floatEVAA TOLERANCE = 1.e-8;

constexpr int DOFDOF = 121;

};  // namespace Constants
}  // namespace EVAA
