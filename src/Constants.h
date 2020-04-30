#pragma once

namespace Constants {
constexpr int ALIGNMENT = 64; /**< alignment for mkl malloc */
constexpr int DIM = 3;        /**< dimension of solution space  consider Constants::DIM = 4 for
                                 efficiency!!! // 10 Constants::DIMensions because of torque of the body */
constexpr int VEC_DIM = 9;    /**< [CG, Wi, Ti] (9 positions / vectors for CG, wheels and tyres) */
constexpr int INCX = 1;       /**< MKL constant incx */
constexpr int NUM_LEGS = 4;   /**< number of tyres, wheels, legs */
constexpr int DOF = 11;       /**< Degrees of Freedom in the Eulerian Frame */
using floatEVAA = double;
constexpr bool USEINTERPOLATION = true;
constexpr floatEVAA TOLERANCE = 1.0e-8;
constexpr int DOFDOF = 121;
// size_t alignment;											/**< alignment for mkl
// malloc */
// const int dim = 3;											/**< dimension of
// solution space */
// const size_t num_tyre = 4;									/**< number of tyres
// */
// const int dim_x_dim = dim * dim;							/**< number of elements of a
// matrix in solution space */
// const int num_wheels_x_dim = num_tyre * dim;				/**< number of elements of a
// matrix regarding the wheels */ const int DOF_diag = 9;
// /**< number of diagonal elements of A */
};  // namespace Constants
