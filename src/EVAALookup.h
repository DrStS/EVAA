// TODO: Copyright check

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
    static T responseFunction(const T& a, const T& b, const T& c, const T& length) { return c * length * length + b * length + a; }

public:
    /**
     * \brief build equidistant grid and the corresponding axis
     * \param[out] grid pointer to grid of size size*k to store y values
     * \param[out] axis pointer to axis of size size to write x value
     * \param[in] size size of one grid
     * \param[in] l_min min length of spring in lookup
     * \param[in] l_max max length of spring in lookup
     * \param[in] a pointer to array of size k with constant coefficient for
     * each spring \param[in] b coefficient for linear part \param[in] c
     * coefficient for quadratic part \param[in] k number of springs
     */
    static void buildLinearGrid(T* grid, T* axis, const size_t& size, const T& l_min, const T& l_max, const T a[8], const T& b, const T& c, const size_t& k) {
        T density = (l_max - l_min) / (size - 1);
#pragma loop(ivdep)
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
     * \param[in] a pointer to array of size k with constant coefficient for
     * each spring \param[in] b coefficient for linear part \param[in] c
     * coefficient for quadratic part \param[in] k number of springs
     */
    static void buildChebyshevGrid(T* grid, T* axis, const size_t& size, const T& l_min, const T& l_max, const T a[8], const T& b, const T& c, const size_t& k) {
#pragma loop(ivdep)
        for (auto i = 0; i < size; i++) {
            axis[i] = (1 + cos((2 * i + 1) / (2 * size) * Constants::PI)) / 2 * (l_max - l_min) + l_min;
        }
        for (auto j = 0; j < k; j++) {
            for (auto i = 0; i < size; i++) {
                grid[i + j * size] = responseFunction(a[j], b, c, axis[i]);
            }
        }
    }
};

/** additional info about break points*/
constexpr MKL_INT xhint = DF_NON_UNIFORM_PARTITION;
/** additional info about function*/
constexpr MKL_INT yhint = DF_NO_HINT;
/** additional info about spline coefficients*/
constexpr MKL_INT scoeffhint = DF_NO_HINT;
/** interpolation results storage format*/
constexpr MKL_INT rhint = DF_NO_HINT;

/**
 * \brief Class in case there are values which are not available in analytical
 * form
 *
 * generate an instance of this class before the simulation and call the
 * interpolate function during the simulation
 */
template <typename T>
class EVAALookup {
private:
    DFTaskPtr* task = nullptr;  /**< Data Fitting task descriptor */
    size_t nx;                  /**< number of break points (data points in the grid)*/
    size_t ny;                  /**< number of functions ( in our case 1 per task but 8 in total)*/
    MKL_INT bc_type = DF_NO_BC; /**< boundary conditions type*/
    MKL_INT ic_type = DF_NO_IC; /**< internal conditions type*/

    T* axis = nullptr;   /**< array of break points (axis with length values of the grid points)*/
    T* grid = nullptr;   /**< function values */
    T* scoeff = nullptr; /**< array of spline coefficients*/
    MKL_INT stype;       /**< spline type: linear = DF_PP_DEFAULT, spline = DF_PP_NATURAL*/
    MKL_INT sorder;      /**< spline order: linear = DF_PP_LINEAR, spline = DF_PP_CUBIC*/
    T l_min, l_max;      /**< to test wheather its inbound */

public:
    /**
     * \brief Constructor of lookup class
     *
     * for linear interpolation:
     * type = DF_PP_DEFAULT, order = DF_PP_LINEAR
     * for spline interpolation:
     * type = DF_PP_NATURAL, order = DF_PP_CUBIC
     *
     * \param size size of one grid
     * \param a pointer to array of size k with constant coefficient for each
     * spring \param b coefficient for linear part of grid function \param c
     * coefficient for quadratic part of grid function \param l_min min length
     * of spring in lookup \param l_max max length of spring in lookup \param k
     * number of springs \param type int which corersponds to a certain type of
     * interpolation \param order order of the interpolation: depends on type
     */
    EVAALookup(const size_t size, const T a[8], const T& b, const T& c, const T& l_min, const T& l_max, const size_t k, const size_t type, const size_t order) :
        //
        nx(size),
        l_min(l_min),
        l_max(l_max),
        ny(k),
        stype(type),
        sorder(order) {
        if (sorder == DF_PP_CUBIC) {
            bc_type = DF_BC_FREE_END;
        }
        // we will get the grid aftwards from a file. That why I do not directly
        // write size, l_min, l_max into the variables
        task = Math::malloc<DFTaskPtr>(ny);
        scoeff = Math::malloc<T>(ny * (nx - 1) * sorder);

        grid = Math::malloc<T>(nx * ny);
        axis = Math::malloc<T>(nx);

        /* create grid */
        EVAAComputeGrid<T>::buildLinearGrid(grid, axis, nx, l_min, l_max, a, b, c, ny);

        int err = 0;
        for (auto i = 0; i < ny; i++) {
            /* Create Data Fitting task */
            err = Math::dfNewTask1D<T>(&task[i], nx, axis, xhint, 1, &grid[i * nx], yhint);
            Math::dfCheckError(err);

            /* Edit task parameters for look up interpolant */
            err = Math::dfEditPPSpline1D<T>(task[i], sorder, stype, bc_type, nullptr, ic_type, nullptr, &scoeff[i * (nx - 1) * sorder], scoeffhint);
            Math::dfCheckError(err);

            /* Construct linear spline using STD method */
            err = Math::dfConstruct1D<T>(task[i], DF_PP_SPLINE, DF_METHOD_STD);
            Math::dfCheckError(err);
        }

        // for debugging purposes
#ifdef WRITECSV
        generateLookupOutputFile(l_min, l_max, a[0]);
#endif  // WRITECSV

        Math::free<T>(grid);
    }

    /**
     * \brief Destructor of lookup class
     *
     * free all the allocated space
     */
    ~EVAALookup() {
        Math::free<T>(axis);
        Math::free<T>(scoeff);

        /* Delete Data Fitting task */
        for (auto i = 0; i < ny; i++) {
            dfDeleteTask(&task[i]);
        }
        Math::free<DFTaskPtr>(task);
    }
    /**
     * \brief interpolates the ny = k grids ob the lookuptable
     *
     * The lookup table has been generated in the initialisation of the object
     * this function uses the calculated coefficients to interpolate certain
     * values
     * \param[in] length  pointer to array of size k with length values of springs
     * \param[out] inter pointer to array of size k to store interpolation values
     */
    void getInterpolation(const T* length, T* inter) const {
        count_interp_debug++;

        // size of array describing derivative (dorder), which is definde two lines below
        const MKL_INT ndorder = 1;
        const MKL_INT dorder[1] = {1};  // only the values are computed

        for (auto i = 0; i < ny; i++) {
            if (length[i] > l_max) {
                std::cout << "Number of getInterpolation calls: " << count_interp_debug << "\n";
                throw std::domain_error("spring length to big for lookup: " + std::to_string(length[i]) + " > " + std::to_string(l_max));
            }
            if (length[i] < l_min) {
                std::cout << "Number of getInterpolation calls: " << count_interp_debug << "\n";
                throw std::domain_error("spring length to small for lookup: " + std::to_string(length[i]) + " < " + std::to_string(l_min));
            }
            Math::dfInterpolate1D<T>(task[i], DF_INTERP, DF_METHOD_PP, 1, &length[i], DF_NO_HINT, ndorder, dorder, nullptr, &inter[i], rhint, 0);
        }
    }
    /*
     * \brief derivative from stifftness k after length
     *
     * The lookup table has been generated in the initialisation of the object
     * this function uses the calculated coefficients to interpolate certain
     * values
     * \param[in] length pointer to array of size k with length values of springs
     * \param[out] deriv pointer to array of size k to store values of the derivative
     */
    void getDerivative(const T* length, T* deriv) const {
        // size of array describing derivative (dorder), which is defined two
        // lines below
        const MKL_INT ndorder = 2;
        const MKL_INT dorder[2] = {0, 1};  // only the derivative values are computed
        for (auto i = 0; i < ny; i++) {
            Math::dfInterpolate1D<T>(task[i], DF_INTERP, DF_METHOD_PP, 1, &length[i], DF_NO_HINT, ndorder, dorder, nullptr, &deriv[i], rhint, 0);
        }
    }

    /**
     * \brief interpolate the first task on every point of axis to check for
     * correctness
     */
    void generateLookupOutputFile(T l_min, T l_max, T add) {
        // size of array describing derivative (order), which is defined two
        // lines below
        const MKL_INT ndorder = 1;
        const MKL_INT dorder[1] = {1};  // only the values are computed
        T* interpolation = Math::malloc<T>(2 * nx - 1);
        T* interpolationPoints = Math::malloc<T>(2 * nx - 1);
        for (auto i = 0; i < (2 * nx - 1); i++) {
            interpolationPoints[i] = l_min + i * (l_max - l_min) / (2 * nx - 2);
            // For Debugging
            Math::dfInterpolate1D<T>(task[0], DF_INTERP, DF_METHOD_PP, 1, &interpolationPoints[i], DF_NO_HINT, ndorder, dorder, nullptr, &interpolation[i], rhint, 0);
        }

        IO::writeLookUpGridPlusInterpolateValues<T>(axis, grid, nx, interpolationPoints, interpolation, 2 * nx - 1, "C:\\software\\repos\\EVAA\\output\\LookupTablePlusInterpolation" + std::to_string(add) + ".txt");
        Math::free<T>(interpolation);
        Math::free<T>(interpolationPoints);
    }
};

}  // namespace EVAA
