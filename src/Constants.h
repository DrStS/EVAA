// TODO: Copyright header

#pragma once

namespace EVAA {
namespace Constants {
/** alignment for mkl malloc */
constexpr size_t ALIGNMENT = 64;

/** MBD solution vector size */
constexpr size_t MBD_SOLUTION_SIZE = 61;

// consider Constants::DIM = 4 for efficiency!!!
// 10 Constants::DIMensions because of torque of the body
/** dimension of solution space */
constexpr size_t DIM = 3;
constexpr size_t lagrangianForceDimension = 2;

/** [CG, Wi, Ti] (9 positions / vectors for CG, wheels and tyres) */
constexpr size_t VEC_DIM = 9;

/** MKL constant incx */
constexpr size_t INCX = 1;

/** number of tyres, wheels, legs */
constexpr size_t NUM_LEGS = 4;

/** Degrees of Freedom in the Eulerian Frame */
constexpr size_t DOF = 11;

/** Deciding over the floating point precision*/
#ifdef DOUBLE_PRECISION
using floatEVAA = double;
#elif SINGLE_PRECISION
using floatEVAA = float;
#endif

constexpr floatEVAA TOLERANCE = 1.e-8;

constexpr size_t DOFDOF = 121;
constexpr size_t DIMDIM = 9;

/** Convention order for front left leg/wheel/tire. */
constexpr size_t FRONT_LEFT = 0;
/** Convention order for front right leg/wheel/tire. */
constexpr size_t FRONT_RIGHT = 1;
/** Convention order for rear left leg/wheel/tire. */
constexpr size_t REAR_LEFT = 2;
/** Convention order for rear right leg/wheel/tire. */
constexpr int REAR_RIGHT = 3;

/** Tyre index in Euler Frame 11DOF */
constexpr int TYRE_INDEX_EULER[4]{4, 6, 8, 10};

/** Tyre index in the Lagrangian Frame Start position*/
constexpr int TYRE_INDEX_LAGRANGE[4]{4, 8, 12, 16};

/** It's PI come on */
constexpr floatEVAA PI = 3.141592653589793238463;
};  // namespace Constants
}  // namespace EVAA
