/**
 * \file EVAALookup.h
 * This file holds the class of EVAALookup and EVAAComputeGrid.
 * \date 04/14/2020
 */

#pragma once

#include <Output.h>

#include <iostream>
#include <vector>

#include "Constants.h"
#include "MathLibrary.h"

namespace EVAA {

/**
 * function returns a readable output for mkl errors
 */
// void CheckDfError(int num);

/**
 * \brief class to generate a grid for interpolation
 *
 * this class is only used for testing purposes
 */
template <typename T>
class EVAAComputeGrid {
private:
    /**
     * \brief Function to feed values to the grid, used for the lookup table
     * \param[in] a coefficient for const part
     * \param[in] b coefficient for linear part
     * \param[in] c coefficient for quadratic part
     * \param[in] length docs for input parameter v
     * \return y
     */
    static T responseFunction(const T a, const T b, const T c, const T length) {
        return c * length * length + b * length + a;
    }

public:
    /**
     * \brief build equidistant grid and the corresponding axis
     * \param[out] grid pointer to grid of size size*k to store y values
     * \param[out] axis pointer to axis of size size to write x value
     * \param[in] size size of one grid
     * \param[in] l_min min length of spring in lookup
     * \param[in] l_max max length of spring in lookup
     * \param[in] a pointer to array of size k with constant coefficient for each spring
     * \param[in] b coefficient for linear part
     * \param[in] c coefficient for quadratic part
     * \param[in] k number of springs - TODO: use Constants.
     */
    static void buildLinearGrid(T* grid, T* axis, int size, T l_min, T l_max, T* a, T b, T c,
                                int k) {
        T density = (l_max - l_min) / (size - 1);
        for (auto i = 0; i < size; i++) {
            axis[i] = l_min + i * density;
        }
        for (auto j = 0; j < k; j++) {
            for (auto i = 0; i < size; i++) {
                grid[i + j * size] = responseFunction(a[j], b, c, axis[i]);
            }
        }
    }

    /**
     * \brief build Chebyshev grid and the corresponding axis
     * \param[out] grid pointer to grid of size size*k to store y values
     * \param[out] axis pointer to axis of size size to write x value
     * \param[in] size size of one grid
     * \param[in] l_min min length of spring in lookup
     * \param[in] l_max max length of spring in lookup
     * \param[in] a pointer to array of size k with constant coefficient for each spring
     * \param[in] b coefficient for linear part
     * \param[in] c coefficient for quadratic part
     * \param[in] k number of springs - TODO: use Constants.
     */
    static void buildChebyshevGrid(T* grid, T* axis, int size, T l_min, T l_max, T* a, T b, T c,
                                   int k) {
        for (auto i = 0; i < size; i++) {
            axis[i] = (1 + cos((2 * i + 1) / (2 * size) * M_PI)) / 2 * (l_max - l_min) + l_min;
        }
        for (auto j = 0; j < k; j++) {
            for (auto i = 0; i < size; i++) {
                grid[i + j * size] = responseFunction(a[j], b, c, axis[i]);
            }
        }
    }
};

/**
 * \brief Class in case there are values which are not available in analytical form
 *
 * generate an instance of this class before the simulation and call the interpolate function during
 * the simulation
 */
template <typename T>
class EVAALookup {
private:
    DFTaskPtr* task = nullptr; /**< Data Fitting task descriptor */
    MKL_INT nx;                /**< number of break points (data points in the grid)*/
    MKL_INT xhint = DF_NON_UNIFORM_PARTITION; /**< additional info about break points*/
    MKL_INT ny;                 /**< number of functions ( in our case 1 per task but 8 in total)*/
    MKL_INT yhint = DF_NO_HINT; /**< additional info about function*/
    MKL_INT scoeffhint = DF_NO_HINT; /**< additional info about spline coefficients*/
    MKL_INT bc_type = DF_NO_BC;      /**< boundary conditions type*/
    MKL_INT ic_type = DF_NO_IC;      /**< internal conditions type*/
    MKL_INT rhint = DF_NO_HINT;      /**< interpolation results storage format*/
    T* axis;               /**< array of break points (axis with length values of the grid points)*/
    T* grid;               /**< function values */
    T* ic = nullptr;       /**< internal conditions*/
    T* bc = nullptr;       /**< boundary conditions*/
    T* scoeff = nullptr;   /**< array of spline coefficients*/
    T* datahint = nullptr; /**< additional info about the structure*/
    MKL_INT stype;         /**< spline type: linear = DF_PP_DEFAULT, spline = DF_PP_NATURAL*/
    MKL_INT sorder;        /**< spline order: linear = DF_PP_LINEAR, spline = DF_PP_CUBIC*/
    T l_min, l_max;        /**< to test wheather its inbound */

public:
    /**
     * \brief Constructor of lookup class
     *
     * for linear interpolation:
     * type = DF_PP_DEFAULT, order = DF_PP_LINEAR
     * for spline interpolation:
     * type = DF_PP_NATURAL, order = DF_PP_CUBIC
     */
    EVAALookup(
        int size /**< [in] size of one grid */,
        T* a /**< [in] pointer to array of size k with constant coefficient for each spring */,
        T b /**< [in] coefficient for linear part of grid function */,
        T c /**< [in] coefficient for quadratic part of grid function */,
        T l_min /**< [in] min length of spring in lookup */,
        T l_max /**< [in] max length of spring in lookup */, int k /**< [in] number of springs */,
        int type /**< [in] int which corersponds to a certain type of interpolation */,
        int order /**< [in] order of the interpolation: depends on type */
        ) :
        nx(size), ny(k), sorder(order), stype(type), l_min(l_min), l_max(l_max) {
        if (sorder == DF_PP_CUBIC) {
            bc_type = DF_BC_FREE_END;
        }
        // we will get the grid aftwards from a file. That why I do not directly write size, l_min,
        // l_max into the variables
        task = (DFTaskPtr*)mkl_malloc(ny * sizeof(DFTaskPtr), Constants::ALIGNMENT);
        grid = (T*)mkl_calloc(nx * ny, sizeof(T), Constants::ALIGNMENT);
        axis = (T*)mkl_calloc(nx, sizeof(T), Constants::ALIGNMENT);
        scoeff = (T*)mkl_calloc(ny * (nx - 1) * sorder, sizeof(T), Constants::ALIGNMENT);

        /* create grid */
        EVAAComputeGrid<T>::buildLinearGrid(grid, axis, nx, l_min, l_max, a, b, c, ny);

        int err = 0;
        for (auto i = 0; i < ny; i++) {
            /* Create Data Fitting task */
            err = dfdNewTask1D(&task[i], nx, axis, xhint, 1, &grid[i * nx], yhint);

            // CheckDfError(err);

            /* Edit task parameters for look up interpolant */
            err = dfdEditPPSpline1D(task[i], sorder, stype, bc_type, bc, ic_type, ic,
                                    &scoeff[i * (nx - 1) * sorder], scoeffhint);

            // CheckDfError(err);

            /* Construct linear spline using STD method */
            err = dfdConstruct1D(task[i], DF_PP_SPLINE, DF_METHOD_STD);
            // CheckDfError(err);
        }

        // for debugging purposes
        generateLookupOutputFile(l_min, l_max, a[0]);

        mkl_free(grid);
        grid = nullptr;
        mkl_free(axis);
        axis = nullptr;
        delete ic;
        ic = nullptr;
        delete bc;
        bc = nullptr;
    }
    /**
     * \brief Destructor of lookup class
     *
     * free all the allocated space
     */
    ~EVAALookup() {
        /* Delete Data Fitting task */
        for (auto i = 0; i < ny; i++) {
            dfDeleteTask(&task[i]);
        }
        mkl_free(scoeff);
        scoeff = nullptr;
        delete datahint;
        datahint = nullptr;
    }
    /**
     * \brief interpolates the ny = k grids ob the lookuptable
     *
     * The lookup table has been generated in the initialisation of the object
     * this function uses the calculated coefficients to interpolate certain values
     */
    void getInterpolation(
        T* length /**< [in] pointer to array of size k with lenght values of springs*/,
        T* inter /**< [out] pointer to array of size k to store interpolation values*/
    ) {
        const MKL_INT ndorder =
            1;  // size of array describing derivative (dorder), which is definde two lines below
        const MKL_INT dorder[1] = {1};  // only the values are computed

        for (auto i = 0; i < ny; i++) {
            if (length[i] > l_max) {
                std::cout << "spring length to big for lookup: " << length[i] << std::endl;
                exit(100);
            }
            if (length[i] < l_min) {
                std::cout << "spring length to small for lookup: " << length[i] << std::endl;
                exit(100);
            }
            dfdInterpolate1D(task[i], DF_INTERP, DF_METHOD_PP, 1, &length[i], DF_NO_HINT, ndorder,
                             dorder, datahint, &inter[i], rhint, 0);
        }
    }
    /*
     * \brief derivative from stifftness k after length
     *
     * The lookup table has been generated in the initialisation of the object
     * this function uses the calculated coefficients to interpolate certain values
     */
    void getDerivative(
        T* length /**< [in] pointer to array of size k with lenght values of springs*/,
        T* deriv /**< [out] pointer to array of size k to store values of the derivative*/
    ) {
        // size of array describing derivative (dorder), which is defined two lines below
        const MKL_INT ndorder = 2;
        const MKL_INT dorder[2] = {0, 1};  // only the derivative values are computed
        for (auto i = 0; i < ny; i++) {
            dfdInterpolate1D(task[i], DF_INTERP, DF_METHOD_PP, 1, &length[i], DF_NO_HINT, ndorder,
                             dorder, datahint, &deriv[i], rhint, 0);
        }
    }

    /**
     * \brief Calculate Matrix for the 11 Dof system with fl,fr,rl,rr notation
     */
    void getMatrixInterpolation(
        T* length /**< [in] pointer to array of size k with lenght values of springs*/,
        T* inter /**< [out] pointer to array of size k to store interpolation values*/,
        T* mat /**< [out] pointer to array of size k * k to store interpolation values*/
    ) {
        getInterpolation(length, inter);
    }

    /**
     * \brief interpolate the first task on every point of axis to check for correctnes
     */
    void generateLookupOutputFile(T l_min, T l_max, T add) {
        // size of array describing derivative (dorder), which is defined two lines below
        const MKL_INT ndorder = 1;
        const MKL_INT dorder[1] = {1};  // only the values are computed
        T* interpolation = (T*)mkl_malloc((2 * nx - 1) * sizeof(T), Constants::ALIGNMENT);
        T* interpolationPoints = (T*)mkl_malloc((2 * nx - 1) * sizeof(T), Constants::ALIGNMENT);
        for (auto i = 0; i < (2 * nx - 1); i++) {
            interpolationPoints[i] = l_min + i * (l_max - l_min) / (2 * nx - 2);
            dfdInterpolate1D(task[0], DF_INTERP, DF_METHOD_PP, 1, &interpolationPoints[i],
                             DF_NO_HINT, ndorder, dorder, datahint, &interpolation[i], rhint, 0);
        }

        IO::writeLookUpGridPlusInterpolateValues<T>(
            axis, grid, nx, interpolationPoints, interpolation, 2 * nx - 1,
            "LookupTablePlusInterpolation" + std::to_string(add) + ".txt");
        mkl_free(interpolation);
    }
};

#if MIGHT_BE_USEFUL
void CheckDfError(int num) {
    switch (num) {
    case DF_ERROR_NULL_TASK_DESCRIPTOR: {
        printf("Error: null task descriptor (code %d).\n", num);
        break;
    }
    case DF_ERROR_MEM_FAILURE: {
        printf("Error: memory allocation failure in DF functionality (code %d).\n", num);
        break;
    }
    case DF_ERROR_BAD_NX: {
        printf("Error: the number of breakpoints is invalid (code %d).\n", num);
        break;
    }
    case DF_ERROR_BAD_X: {
        printf("Error: the array which contains the breakpoints is not defined (code
               % d)
            .\n ", num); break;
    }
    case DF_ERROR_BAD_X_HINT: {
        printf("Error: invalid flag describing structure of partition (code %d).\n", num);
        break;
    }
    case DF_ERROR_BAD_NY: {
        printf("Error: invalid dimension of vector-valued function y (code %d).\n", num);
        break;
    }
    case DF_ERROR_BAD_Y: {
        printf("Error: the array which contains function values is invalid (code %d).\n", num);
        break;
    }
    case DF_ERROR_BAD_Y_HINT: {
        printf("Error: invalid flag describing structure of function y (code %d).\n", num);
        break;
    }
    case DF_ERROR_BAD_SPLINE_ORDER: {
        printf("Error: invalid spline order (code %d).\n", num);
        break;
    }
    case DF_ERROR_BAD_SPLINE_TYPE: {
        printf("Error: invalid type of the spline (code %d).\n", num);
        break;
    }
    case DF_ERROR_BAD_IC_TYPE: {
                printf("Error: invalid type of internal conditions used in the spline construction
(code %d).\n", num); break;
    }
    case DF_ERROR_BAD_IC: {
                printf("Error: array of internal conditions for spline construction is not defined
(code %d).\n", num); break;
    }
    case DF_ERROR_BAD_BC_TYPE: {
                printf("Error: invalid type of boundary conditions used in the spline construction
(code %d).\n", num); break;
    }
    case DF_ERROR_BAD_BC: {
                printf("Error: array which presents boundary conditions for spline construction is
not defined (code %d).\n", num); break;
    }
    case DF_ERROR_BAD_PP_COEFF: {
                printf("Error: array of piece-wise polynomial spline coefficients is not defined
(code %d).\n", num); break;
    }
    case DF_ERROR_BAD_PP_COEFF_HINT: {
                printf("Error: invalid flag describing structure of the piece-wise polynomial spline
coefficients (code %d).\n", num); break;
    }
    case DF_ERROR_BAD_PERIODIC_VAL: {
                printf("Error: function values at the end points of the interpolation interval are
not equal as required in periodic boundary conditions (code %d).\n", num); break;
    }
    case DF_ERROR_BAD_DATA_ATTR: {
                printf("Error: invalid attribute of the pointer to be set or modified in Data
Fitting task descriptor with EditIdxPtr editor (code %d).\n", num); break;
    }
    case DF_ERROR_BAD_DATA_IDX: {
                printf("Error: index of pointer to be set or modified in Data Fitting task
descriptor with EditIdxPtr editor is out of range (code %d).\n", num); break;
    }
    case DF_ERROR_BAD_NSITE: {
        printf("Error: invalid number of interpolation sites (code %d).\n", num);
        break;
    }
    case DF_ERROR_BAD_SITE: {
        printf("Error: array of interpolation sites is not defined (code %d).\n", num);
        break;
    }
    case DF_ERROR_BAD_SITE_HINT: {
        printf("Error: invalid flag describing structure of interpolation sites (code
               % d)
            .\n ", num); break;
    }
    case DF_ERROR_BAD_NDORDER: {
                printf("Error: invalid size of array that defines order of the derivatives to be
computed at the interpolation sites (code %d).\n", num); break;
    }
    case DF_ERROR_BAD_DORDER: {
                printf("Error: array defining derivative orders to be computed at interpolation
sites is not defined (code %d).\n", num); break;
    }
    case DF_ERROR_BAD_DATA_HINT: {
                printf("Error: invalid flag providing a-priori information about partition and/or
interpolation sites (code %d).\n", num); break;
    }
    case DF_ERROR_BAD_INTERP: {
        printf("Error: array of spline based interpolation results is not defined (code
               % d)
            .\n ", num); break;
    }
    case DF_ERROR_BAD_INTERP_HINT: {
                printf("Error: invalid flag defining structure of spline based interpolation results
(code %d).\n", num); break;
    }
    case DF_ERROR_BAD_CELL_IDX: {
                printf("Error: array of indices of partition cells containing interpolation sites is
not defined (code %d).\n", num); break;
    }
    case DF_ERROR_BAD_NLIM: {
        printf("Error: invalid size of arrays containing integration limits (code %d).\n", num);
        break;
    }
    case DF_ERROR_BAD_LLIM: {
        printf("Error: array of left integration limits is not defined (code %d).\n", num);
        break;
    }
    case DF_ERROR_BAD_RLIM: {
        printf("Error: array of right integration limits is not defined (code %d).\n", num);
        break;
    }
    case DF_ERROR_BAD_INTEGR: {
        printf("Error: array of spline based integration results is not defined (code
               % d)
            .\n ", num); break;
    }
    case DF_ERROR_BAD_INTEGR_HINT: {
                printf("Error: invalid flag defining structure of spline based integration results
(code %d).\n", num); break;
    }
    case DF_ERROR_BAD_LOOKUP_INTERP_SITE: {
        printf("Error: bad site provided for interpolation with look-up interpolator (code
               % d)
            .\n ", num); break;
    }
    case DF_ERROR_NULL_PTR: {
        printf("Error: bad pointer provided in DF function (code %d).\n", num);
        break;
    }
    default:
        break;
    }

    if (num < 0) {
        exit(1);
    }
}
#endif

}  // namespace EVAA
