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

// TODO: FIX the single precision problem
#ifdef DOUBLE_PRECISION
using floatEVAA = double;
#elif SINGLE_PRECISION
using floatEVAA = float;
#endif

constexpr floatEVAA TOLERANCE = 1.e-8;

constexpr int DOFDOF = 121;

/** Convention order for front left leg/wheel/tire. */
constexpr int FRONT_LEFT = 0;
/** Convention order for front right leg/wheel/tire. */
constexpr int FRONT_RIGHT = 1;
/** Convention order for rear left leg/wheel/tire. */
constexpr int REAR_LEFT = 2;
/** Convention order for rear right leg/wheel/tire. */
constexpr int REAR_RIGHT = 3;
/** It's PI come on */
const double PI = 3.141592653589793238463;
};  // namespace Constants
}  // namespace EVAA
