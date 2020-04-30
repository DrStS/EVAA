// TODO: Copyright header

#include "MKLRoutines.h"

#include <stdexcept>

namespace EVAA {
namespace Math {
// namespace MKL {

void potrfCheckStatus(lapack_int status) {
    if (status == 1) {
        throw std::domain_error("Matrix non Positive Definite");
    }
    if (status == -1) {
        throw std::domain_error("Matrix contain illegal value");
    }
}

void dfCheckError(const int num) {
    switch (num) {

    case DF_ERROR_NULL_TASK_DESCRIPTOR:
        throw std::domain_error("Null task descriptor. Code: " + std::to_string(num));
        break;

    case DF_ERROR_MEM_FAILURE:
        throw std::domain_error("Memory allocation failure in DF functionality. Code: " +
                                std::to_string(num));
        break;

    case DF_ERROR_BAD_NX:
        throw std::domain_error("The number of breakpoints is invalid. Code: " +
                                std::to_string(num));
        break;

    case DF_ERROR_BAD_X:
        throw std::domain_error("The array which contains the breakpoints is not defined. Code: " +
                                std::to_string(num));
        break;

    case DF_ERROR_BAD_X_HINT:
        throw std::domain_error("Invalid flag describing structure of partition. Code: " +
                                std::to_string(num));
        break;

    case DF_ERROR_BAD_NY:
        throw std::domain_error("Invalid dimension of vector-valued function y. Code: " +
                                std::to_string(num));
        break;

    case DF_ERROR_BAD_Y:
        throw std::domain_error("The array which contains function values is invalid. Code: " +
                                std::to_string(num));
        break;

    case DF_ERROR_BAD_Y_HINT:
        throw std::domain_error("Invalid flag describing structure of function y. Code: " +
                                std::to_string(num));
        break;

    case DF_ERROR_BAD_SPLINE_ORDER:
        throw std::domain_error("Invalid spline order. Code: " + std::to_string(num));
        break;

    case DF_ERROR_BAD_SPLINE_TYPE:
        throw std::domain_error("Invalid type of the spline. Code: " + std::to_string(num));
        break;

    case DF_ERROR_BAD_IC_TYPE:
        throw std::domain_error(
            "Invalid type of internal conditions used in the spline construction. Code: " +
            std::to_string(num));
        break;

    case DF_ERROR_BAD_IC:
        throw std::domain_error(
            "Array of internal conditions for spline construction is not defined. Code: " +
            std::to_string(num));
        break;

    case DF_ERROR_BAD_BC_TYPE:
        throw std::domain_error(
            "Invalid type of boundary conditions used in the spline construction. Code: " +
            std::to_string(num));
        break;

    case DF_ERROR_BAD_BC:
        throw std::domain_error(
            "Array which presents boundary conditions for spline construction is not defined. "
            "Code: " +
            std::to_string(num));
        break;

    case DF_ERROR_BAD_PP_COEFF:
        throw std::domain_error(
            "Array of piece-wise polynomial spline coefficients is not defined. Code: " +
            std::to_string(num));
        break;

    case DF_ERROR_BAD_PP_COEFF_HINT:
        throw std::domain_error(
            "Invalid flag describing structure of the piece-wise polynomial spline coefficients. "
            "Code: " +
            std::to_string(num));
        break;

    case DF_ERROR_BAD_PERIODIC_VAL:
        throw std::domain_error(
            "Function values at the end points of the interpolation interval are not equal as "
            "required in periodic boundary conditions. Code: " +
            std::to_string(num));
        break;

    case DF_ERROR_BAD_DATA_ATTR:
        throw std::domain_error(
            "Invalid attribute of the pointer to be set or modified in Data Fitting task "
            "descriptor with EditIdxPtr editor. Code: " +
            std::to_string(num));
        break;

    case DF_ERROR_BAD_DATA_IDX:
        throw std::domain_error(
            "Index of pointer to be set or modified in Data Fitting task descriptor with "
            "EditIdxPtr editor is out of range. Code: " +
            std::to_string(num));
        break;

    case DF_ERROR_BAD_NSITE:
        throw std::domain_error("Invalid number of interpolation sites. Code: " +
                                std::to_string(num));
        break;

    case DF_ERROR_BAD_SITE:
        throw std::domain_error("Array of interpolation sites is not defined. Code: " +
                                std::to_string(num));
        break;

    case DF_ERROR_BAD_SITE_HINT:
        throw std::domain_error("Invalid flag describing structure of interpolation sites. Code: " +
                                std::to_string(num));
        break;

    case DF_ERROR_BAD_NDORDER:
        throw std::domain_error(
            "Invalid size of array that defines order of the derivatives to be computed at the "
            "interpolation sites. Code: " +
            std::to_string(num));
        break;

    case DF_ERROR_BAD_DORDER:
        throw std::domain_error(
            "Array defining derivative orders to be computed at interpolation sites is not "
            "defined. Code: " +
            std::to_string(num));
        break;

    case DF_ERROR_BAD_DATA_HINT:
        throw std::domain_error(
            "Invalid flag providing a-priori information about partition and/or interpolation "
            "sites. Code: " +
            std::to_string(num));
        break;

    case DF_ERROR_BAD_INTERP:
        throw std::domain_error(
            "Array of spline based interpolation results is not defined. Code: " +
            std::to_string(num));
        break;

    case DF_ERROR_BAD_INTERP_HINT:
        throw std::domain_error(
            "Invalid flag defining structure of spline based interpolation results. Code: " +
            std::to_string(num));
        break;

    case DF_ERROR_BAD_CELL_IDX:
        throw std::domain_error(
            "Array of indices of partition cells containing interpolation sites is not defined. "
            "Code: " +
            std::to_string(num));
        break;

    case DF_ERROR_BAD_NLIM:
        throw std::domain_error("Invalid size of arrays containing integration limits. Code: " +
                                std::to_string(num));
        break;

    case DF_ERROR_BAD_LLIM:
        throw std::domain_error("Array of left integration limits is not defined. Code: " +
                                std::to_string(num));
        break;

    case DF_ERROR_BAD_RLIM:
        throw std::domain_error("Array of right integration limits is not defined. Code: " +
                                std::to_string(num));
        break;

    case DF_ERROR_BAD_INTEGR:
        throw std::domain_error("Array of spline based integration results is not defined. Code: " +
                                std::to_string(num));
        break;

    case DF_ERROR_BAD_INTEGR_HINT:
        throw std::domain_error(
            "Invalid flag defining structure of spline based integration results. Code: " +
            std::to_string(num));
        break;

    case DF_ERROR_BAD_LOOKUP_INTERP_SITE:
        throw std::domain_error(
            "Bad site provided for interpolation with look-up interpolator. Code: " +
            std::to_string(num));
        break;

    case DF_ERROR_NULL_PTR:
        throw std::domain_error("Bad pointer provided in DF function. Code: " +
                                std::to_string(num));
        break;

    default:
        if (num < 0) {
            throw std::domain_error("Unknown df error: Code: " + std::to_string(num));
        }
        break;
    }
}

// }  // namespace MKL
}  // namespace Math
}  // namespace EVAA
