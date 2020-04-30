/*
 * Copyright &copy; 2019, Stefan Sicklinger, Munich
 *
 *  All rights reserved.
 *
 *  This file is part of EVAA.
 *
 *  EVAA is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  EVAA is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with EVAA.  If not, see http://www.gnu.org/licenses/.
 */

/**
 * \file MathLibrary.h
 * This file holds the math funcktion of EVAA.
 * \date 6/13/2019
 */

#pragma once

#include <cmath>
#include <iostream>
#include <vector>

#include "AuxiliaryParameters.h"
#ifdef USE_INTEL_MKL
#include <mkl.h>

#include "BLAS.h"
#endif

namespace EVAA {

/**
 * Implements several useful math functions and solvers
 */
namespace MathLibrary {

/**
 * \brief Compute dense symmetrix matrix LU factorisation
 * \param[in] _nElements number of rows = number of columns
 * \param[in] _A matrix
 * \param[in] _pivot elements
 * \author Stefan Sicklinger
 */
void computeDenseSymLUFactorisation(const int _nElements, std::vector<double>& _A,
                                    std::vector<int>& _pivots);

/**
 * \brief Compute backward/forward substitution
 * \param[in] _nElements number of rows = number of columns
 * \param[in] _A matrix
 * \param[in] _pivot elements
 * \author Stefan Sicklinger
 */
void computeDenseSymSolution(const int _nElements, std::vector<double>& _A,
                             std::vector<int>& _pivots, std::vector<double>& _rhs);

/**
 * \brief Computes a vector-scalar product and adds the result to a vector. vec2 <- a*vec1 + vec2
 * \param[in] _vec1 the 1st vector
 * \param[in] _vec2 the 2nd vector
 * \param[in] _alpha   scalar
 * \param[in] _nElements number of elements in vec1
 * \author Stefan Sicklinger
 */
void computeDenseVectorAddition(double* vec1, double* vec2, const double _alpha,
                                const int _nElements);

/**
 * \brief Print info of the Intel MKL
 * \author Stefan Sicklinger
 */
void printMKLInfo(void);

/**
 * Replaces all values in an array which are below a threshold to a value
 * \param arr vector to be adapted
 * \param start index to start the update
 * \param end index to end the update
 * \param value to be written to the elements
 * \param tol threshold
 */
template <typename T>
void transmutate_elements(T* arr, size_t start, size_t end, T value, T tol)
{
#pragma loop(ivdep)
    for (size_t i = start; i < end; ++i) {
        if (arr[i] * arr[i] < tol * tol) {
            arr[i] = value;
        }
    }
}

/**
 * Performs the elementwise multiplication of two vectors
 * \param v1 first vector
 * \param v2 second vector
 * \param dim length of the vector
 * \return result vector
 */
template <typename T>
void vector_elem_wise_product(T* v1, T* v2, T* result, size_t dim)
{
#pragma loop(ivdep)
    for (size_t i = 0; i < dim; i++) {
        result[i] = v1[i] * v2[i];
    }
}

/**
 * Creates a diagonal matrix with one same value on all diagonal entries (or overwrites the
 * diagonal elements)
 * \param dim size of the matrix (dim x dim)
 * \param val value to be written to all diagonal entries
 * \return mat the diagonal matrix
 */
template <typename T>
void diagonal_matrix(T* mat, size_t dim, T val)
{
#pragma loop(ivdep)
    for (size_t i = 0; i < dim; ++i) {
        mat[i * dim + i] = val;
    }
}

/**
 * outputs the determinant
 */
template <typename T>
void check_status(T status)
{
    if (status == 1) {
        std::cout << "Matrix non Positive Definite" << std::endl;
        exit(5);
    }
    else if (status == -1) {
        std::cout << "Matrix contain illegal value" << std::endl;
        exit(5);
    }
}

/**
 * Broken function (not safe to use)
 */
template <typename T>
void elementwise_inversion(T* vect, size_t dim)
{
    int alignment = 64;
    T* temp = (T*)mkl_calloc(dim, sizeof(T), alignment);
    mkl<T>::vInv(dim, vect, temp);
    mkl<T>::copy(dim, temp, 1, vect, 1);
    // for (size_t i = 0; i < dim; ++i) {
    //	vect[i] = 1.0 / vect[i];
    //}
    mkl_free(temp);
}

/**
 * Swaps the adresses of the elements a and b
 */
template <typename T>
void swap_address(T*& a, T*& b)
{
    T* c = a;
    a = b;
    b = c;
}

/**
 * Creates a diagonal matrix with diagonal entries from a vector(or overwrites the diagonal
 * elements)
 * \param dim size of the matrix (dim x dim)
 * \param vector to be written to all diagonal entries
 * \return mat the diagonal matrix
 */
template <typename T>
void allocate_to_diagonal(T* matrix, T* vector, size_t dim)
{
    mkl<T>::copy(dim, vector, 1, matrix, dim + 1);
    /*for (int i = 0; i < dim; ++i) {
            matrix[i*dim + i] = vector[i];
    }*/
}

/**
 * Performs the cross product of two vectors of size 3
 * \param vect_A first vector
 * \param vect_B second vector
 * \return cross_P vect_A x vect_B
 */
template <typename T>
void crossProduct(const T* vect_A, const T* vect_B, T* cross_P)

{
    cross_P[0] = vect_A[1] * vect_B[2] - vect_A[2] * vect_B[1];
    cross_P[1] = vect_A[2] * vect_B[0] - vect_A[0] * vect_B[2];
    cross_P[2] = vect_A[0] * vect_B[1] - vect_A[1] * vect_B[0];
}

/**
 * Performs the cross product of a vector of size 3 with [0,0,1]
 * \param vect_A first vector; A = [x, y, z]
 * \return cross_P = vect_A x [0, 0, 1] = [y, -x, 0]
 */
template <typename T>
void crossProduct_unitvecZ(const T* vect_A, T* cross_P)

{
    cross_P[0] = vect_A[1];
    cross_P[1] = -vect_A[0];
    cross_P[2] = 0;
}

/**
 * Performs the cross product of a vector of size 3 with [0,0,1]
 * ignoring the third component
 * \param vect_A first vector; A = [x, y, z]
 * \return cross_P = vect_A x [0, 0, 1] = [y, -x, 0] i.e. [y,-x]
 */
template <typename T>
void crossProduct_unitvecZ2D(const T* vect_A, T* cross_P)

{
    cross_P[0] = vect_A[1];
    cross_P[1] = -vect_A[0];
}

/**
 * Performs the cross product of 2 vector A and B with no component on Z direction
 * \param vect_A first vector; A = [xA, yA, 0]
 * \param vect_B second vector; B = [xB, yB, 0]
 * \return cross_P = vect_A x vect_B = [0, 0, xA * yB - xB * yA]
 */
template <typename T>
void crossProduct_unitvecXY(const T* vect_A, const T* vect_B, T* cross_P)

{
    cross_P[0] = 0;
    cross_P[1] = 0;
    cross_P[2] = vect_A[0] * vect_B[1] - vect_A[1] * vect_B[0];
}

/**
 * Converts quaternions to Euler angles from
 * https://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles
 * \param q quaterion in format [XYZW]
 * \return E Euler angles in format [XYZ]
 */
template <typename T>
void ToEulerAngles(const T* q, T* E)
{
    // Converts Quaternion of form w,x,y,z in euler angle

    // roll (x-axis rotation)
    T sinr_cosp = 2 * (q[3] * q[0] + q[1] * q[2]);
    T cosr_cosp = 1 - 2 * (q[0] * q[0] + q[1] * q[1]);
    E[0] = std::atan2(sinr_cosp, cosr_cosp);
    double M_PI = 3.14159265358979323846;
    // pitch (y-axis rotation)
    T sinp = 2 * (q[3] * q[1] - q[2] * q[0]);
    if (std::abs(sinp) >= 1)
        E[1] = std::copysign(M_PI / 2, sinp);  // use 90 degrees if out of range
    else
        E[1] = std::asin(sinp);

    // yaw (z-axis rotation)
    T siny_cosp = 2 * (q[3] * q[2] + q[0] * q[1]);
    T cosy_cosp = 1 - 2 * (q[1] * q[1] + q[2] * q[2]);
    E[2] = std::atan2(siny_cosp, cosy_cosp);
}

/**
 * Computes the rotation matrix out of 3 angle representation
 * \result R rotation matrix
 */
template <typename T>
void get_rotation_matrix(const T yaw, const T pitch, const T roll, T* R)
{
    R[0] = (std::cos(yaw)) * (std::cos(pitch));
    R[1] = ((std::cos(yaw)) * (std::sin(pitch)) * (std::sin(roll)) -
            (std::sin(yaw)) * (std::cos(roll)));
    R[2] = ((std::cos(yaw)) * (std::sin(pitch)) * (std::cos(roll)) +
            (std::sin(yaw)) * (std::sin(roll)));
    R[3] = (std::sin(yaw)) * (std::cos(pitch));
    R[4] = ((std::sin(yaw)) * (std::sin(pitch)) * (std::sin(roll)) +
            (std::cos(yaw)) * (std::cos(roll)));
    R[5] = ((std::sin(yaw)) * (std::sin(pitch)) * (std::cos(roll)) -
            (std::cos(yaw)) * (std::sin(roll)));
    R[6] = -(std::sin(pitch));
    R[7] = (std::cos(pitch)) * (std::sin(roll));
    R[8] = (std::cos(pitch)) * (std::cos(roll));
}

/**
 * Computes the angle between two vectors
 * \param v1 first vector
 * \param v2 second vector
 * \param dim size of the vectors
 * \return angle in radians
 */
template <typename T>
T compute_angle_using_dot_product(const T* v1, const T* v2, size_t dim)
{
    T angle;
    T nrm_v1, nrm_v2;
    nrm_v1 = mkl<T>::nrm2(dim, v1, 1);
    nrm_v2 = mkl<T>::nrm2(dim, v2, 1);
    const MKL_INT DIM = dim;
    const MKL_INT incx = 1;
    angle = mkl<T>::dot(DIM, v1, incx, v2, incx);
    // angle = dot_product<T>(v1, v2, dim);
    angle = angle / (nrm_v1 * nrm_v2);
    angle = angle >= 1 ? 1.0 : angle;
    angle = std::acosl(angle);
    return angle;
}

/**
 * print a vector
 * \param vect the data
 * \param count its size
 */
template <typename T>
void write_vector(T* vect, int count)
{
    std::cout << "Debug mode print" << std::endl;
    for (size_t i = 0; i < count; ++i) {
        // std::cout << vect[i] << std::endl;
        std::cout.precision(15);
        std::cout << std::scientific << vect[i] << std::endl;
        // printf("%1.15f\n", vect[i]);
    }
    // exit(5);
}

/**
 * print a matrix
 * \param vect the data
 * \param count its size=count x count
 */
template <typename T>
void write_matrix(T* vect, int count)
{
    std::cout << "Debug mode print" << std::endl;
    for (size_t i = 0; i < count; ++i) {
        // std::cout << vect[i] << std::endl;
        std::cout.precision(5);
        for (size_t j = 0; j < count; ++j) {
            std::cout << std::scientific << vect[i * count + j] << "  ";
        }
        std::cout << "\n" << std::endl;
        // printf("%1.15f\n", vect[i]);
    }
    // exit(5);
}

/**
 * Calculates the rotation quaternion between two angles
 * \param v1 first vector
 * \param v2 second vector
 * \param dim size of the vector
 * \return angle the angle between the two vectors (in radians)
 * \return rotation_axis unit axis
 */
template <typename T>
void get_quaternion(const T* v1, const T* v2, T* angle, T* rotation_axis, size_t dim_)
{
    T nrm_v1, nrm_v2, ra_nrm, q_nrm;
    nrm_v1 = mkl<T>::nrm2(dim_, v1, 1);
    nrm_v2 = mkl<T>::nrm2(dim_, v2, 1);
    crossProduct(v1, v2, rotation_axis);
    mkl<T>::scal(dim_, 1.0 / (nrm_v1 * nrm_v2), rotation_axis, 1);
    *angle = compute_angle_using_dot_product<T>(v1, v2, dim_);
    ra_nrm = mkl<T>::nrm2(dim_, rotation_axis, 1);
    T eps = 1e-8;
    if (ra_nrm < eps) {
        *angle = 0.0;
        if ((v1[2] >= 0 ? v1[2] : -v1[2]) > eps) {
            rotation_axis[0] = -1;
            rotation_axis[1] = -1;
            rotation_axis[2] = (v1[0] + v1[1]) / v1[2];
        }
        else if ((v1[1] >= 0 ? v1[1] : -v1[1]) > eps) {
            rotation_axis[0] = -1;
            rotation_axis[1] = (v1[0] + v1[2]) / v1[1];
            rotation_axis[2] = -1;
        }
        else if ((v1[0] >= 0 ? v1[0] : -v1[0]) > eps) {
            rotation_axis[0] = (v1[2] + v1[1]) / v1[0];
            rotation_axis[1] = -1;
            rotation_axis[2] = -1;
        }
        else {
            rotation_axis[0] = 1;
            rotation_axis[1] = 0;
            rotation_axis[2] = 0;
        }
        ra_nrm = mkl<T>::nrm2(dim_, rotation_axis, 1);
    }
    mkl<T>::scal(dim_, 1.0 / (ra_nrm), rotation_axis, 1);
}

/**
 * Calculates the rotation quaternion between two angles
 * \param v1 first vector
 * \param v2 second vector
 * \param dim size of the vector
 * \return angle the angle between the two vectors (in radians)
 * \return rotation_axis unit axis
 * \return q the quaternion
 */
template <typename T>
void get_quaternion(const T* v1, const T* v2, T* q, T* angle, T* rotation_axis, size_t dim_)
{
    get_quaternion(v1, v2, angle, rotation_axis, dim_);
    q[0] = rotation_axis[0] * std::sin(*angle / 2.0);
    q[1] = rotation_axis[1] * std::sin(*angle / 2.0);
    q[2] = rotation_axis[2] * std::sin(*angle / 2.0);
    q[3] = std::cos(*angle / 2.0);
    q_nrm = mkl<T>::nrm2(4, q, 1);
    mkl<T>::scal(4, 1.0 / (q_nrm), q, 1);
}

/**
 * get the basis vectors connected to an orientation
 * \param initial_orientation the quaternion
 * \return transformed_basis the basis vectors in matrix form
 */
template <typename T>
void get_basis(const T* initial_orientation, T* transfrormed_basis)
{
    // quaternion initial_orientation yields basis(calculates basically the matrix rotation
    // the euclidian basis to the local basis)

    // get local basis vectors using unit quaternion rotation
    // s = 1 / norm(q) ^ 2;      %normalizer, only to remove numerical stuffy stuff
    size_t quad_dim = 4;
    size_t dim = 3;
    T nrm = mkl<T>::nrm2(quad_dim, initial_orientation, 1);
    T s = 1.0 / (nrm * nrm);
    mkl<T>::scal(dim * dim, 0.0, transfrormed_basis, 1);
    size_t i, j;

    T quad_sum = 0.0;
    for (i = 0; i < dim; ++i) {
        quad_sum += initial_orientation[i] * initial_orientation[i];
    }
    i = 0;
    j = 0;
    transfrormed_basis[i * dim + j] =
        1 - 2 * s * (quad_sum - initial_orientation[i] * initial_orientation[i]);
    j = 1;
    transfrormed_basis[i * dim + j] = 2 * s *
                                      (initial_orientation[i] * initial_orientation[j] -
                                       initial_orientation[i + 2] * initial_orientation[j + 2]);
    j = 2;
    transfrormed_basis[i * dim + j] = 2 * s *
                                      (initial_orientation[i] * initial_orientation[j] +
                                       initial_orientation[i + 1] * initial_orientation[j + 1]);
    i = 1;
    j = 0;
    transfrormed_basis[i * dim + j] = 2 * s *
                                      (initial_orientation[i] * initial_orientation[j] +
                                       initial_orientation[i + 2] * initial_orientation[j + 2]);
    j = 1;
    transfrormed_basis[i * dim + j] =
        1 - 2 * s * (quad_sum - initial_orientation[i] * initial_orientation[i]);
    j = 2;
    transfrormed_basis[i * dim + j] = 2 * s *
                                      (initial_orientation[i] * initial_orientation[j] -
                                       initial_orientation[i - 1] * initial_orientation[j + 1]);
    i = 2;
    j = 0;
    transfrormed_basis[i * dim + j] = 2 * s *
                                      (initial_orientation[i] * initial_orientation[j] -
                                       initial_orientation[i + 1] * initial_orientation[j + 1]);
    j = 1;
    transfrormed_basis[i * dim + j] = 2 * s *
                                      (initial_orientation[i] * initial_orientation[j] +
                                       initial_orientation[i + 1] * initial_orientation[j - 1]);
    j = 2;
    transfrormed_basis[i * dim + j] =
        1 - 2 * s * (quad_sum - initial_orientation[i] * initial_orientation[i]);
}

/**
 * Performs the cross product in matrix form
 * a x b = a_tilda * b
 * \param input_vector a
 * \return tilda the matrix a_tilda
 */
template <typename T>
void get_tilda(const T* input_vector, T* tilda_output)
{
    /*
    This function is only suitable for a 3 dimensional system and renders unusable, might throw
    exceptions when used with other dimensions.
    given y: x_tilda*y = cross(x,y) [Stoneking, page 3 bottom]
    x_tilda = [	0, -x(3), x(2) x(3), 0, -x(1) -x(2), x(1), 0	]
    */
    tilda_output[0] = 0;
    tilda_output[1] = -input_vector[2];
    tilda_output[2] = input_vector[1];
    tilda_output[3] = input_vector[2];
    tilda_output[4] = 0;
    tilda_output[5] = -input_vector[0];
    tilda_output[6] = -input_vector[1];
    tilda_output[7] = input_vector[0];
    tilda_output[8] = 0;
}

/**
Implements different numerical schemes
*/
template <typename T, class C>
class Solvers {
public:
    /**
     *  Implicit Euler with Broyden scheme
     *  \param tol tolerance of the Broyden iteration
     *  \param max_iter of the Broyden iteration
     *  \param num_time_iter number of time steps to perform
     *  \param dt timestep
     *  \param x_previous previous solution
     *  \return new solution
     */
    static void Broyden_Euler(C* obj, T* x_previous, T* x_vector_new, T dt, size_t num_time_iter,
                              T tol, size_t max_iter)
    {
        size_t alignment = obj->get_alignment();
        size_t x_len = obj->get_solution_dimension();
        T* f_old = (T*)mkl_malloc(sizeof(T) * x_len, alignment);
        T* f_new = (T*)mkl_malloc(sizeof(T) * x_len, alignment);
        T* dx = (T*)mkl_malloc(sizeof(T) * x_len, alignment);
        T* dx_inv = (T*)mkl_malloc(sizeof(T) * x_len, alignment);
        T* df = (T*)mkl_malloc(sizeof(T) * x_len, alignment);
        T* x = (T*)mkl_malloc(sizeof(T) * x_len, alignment);
        T* F = (T*)mkl_malloc(sizeof(T) * x_len, alignment);
        T* F_new = (T*)mkl_malloc(sizeof(T) * x_len, alignment);
        T* dF = (T*)mkl_malloc(sizeof(T) * x_len, alignment);
        T* J = (T*)mkl_calloc(x_len * x_len, sizeof(T), alignment);
        T* J_tmp = (T*)mkl_calloc(x_len * x_len, sizeof(T), alignment);
        lapack_int* piv = (lapack_int*)mkl_malloc(sizeof(lapack_int*) * x_len, alignment);

        T *it_start, *curr_time_start_pos;
        T* x_new = (T*)mkl_malloc(sizeof(T) * x_len, alignment);
        T* switcher;  // used for interchanging adresses between F and F_new
        T eps = 1e-4;
        double nrm = 0.001;
        T t = 0.0;
        T val = 0.0;
        mkl<T>::copy(x_len, x_previous, 1, x_vector_new, 1);
        for (size_t i = 1; i <= num_time_iter; ++i) {
            // set the pointer to the start of previous timestep
            it_start = x_vector_new + (i - 1) * x_len;

            // 1. Initialize guess from previous time step
            // f_old = f(t(n-1), x_previous');
            obj->compute_f3D_reduced(it_start, t, f_old);

            // in case the velocity is 0 add nuggets to avoid singular matrices (slow check for
            // improvement) f_old(abs(f_old) < 0.001) = 0.001;
            val = eps * (2 * (1 + rand() % 2) - 3);
            transmutate_elements(f_old, 0, x_len, val, eps);

            // 2. Initial guess from Explicit Euler
            // x = x_previous + dt * f_old';
            // f_new = f(t(n), x');
            t += dt;
            mkl<T>::copy(x_len, it_start, 1, x, 1);
            mkl<T>::axpy(x_len, dt, f_old, 1, x, 1);
            obj->compute_f3D_reduced(x, t, f_new);

            // Initial approximation of the Jacobian
            // dx = x - x_previous;
            mkl<T>::copy(x_len, x, 1, dx, 1);
            mkl<T>::axpy(x_len, -1.0, it_start, 1, dx, 1);
            // df = f_new' - f_old';
            mkl<T>::copy(x_len, f_new, 1, df, 1);
            mkl<T>::axpy(x_len, -1.0, f_old, 1, df, 1);

            // approximate J(x_0)
            // J = eye(length(df)) - delta_t*((1./dx)'*df)';
            mkl<T>::scal(x_len * x_len, 0.0, J, 1);
            diagonal_matrix(J, x_len, 1.0);
            mkl<T>::copy(x_len, dx, 1, dx_inv, 1);
            elementwise_inversion(dx_inv, x_len);
            mkl<T>::ger(CblasRowMajor, x_len, x_len, -dt, df, 1, dx_inv, 1, J, x_len);

            // calculate initial F for stopping condition
            // x_dot = f_new';
            // F = dx - delta_t * x_dot;
            mkl<T>::copy(x_len, dx, 1, F, 1);
            mkl<T>::axpy(x_len, -dt, f_new, 1, F, 1);

            // Broyden's Method
            for (size_t j = 0; j < max_iter; ++j) {
                if (mkl<T>::nrm2(x_len, F, 1) < tol) {
                    break;
                }

                // x(i+1) = x(i) - J^(-1)*g(x(i))
                mkl<T>::copy(x_len * x_len, J, 1, J_tmp, 1);
                mkl<T>::getrf(LAPACK_ROW_MAJOR, x_len, x_len, J_tmp, x_len, piv);
                mkl<T>::copy(x_len, F, 1, x_new, 1);
                mkl<T>::getrs(LAPACK_ROW_MAJOR, 'N', x_len, 1, J_tmp, x_len, piv, x_new, 1);
                mkl<T>::axpy(x_len, -1.0, x, 1, x_new, 1);  // result here is -x_new
                mkl<T>::scal(x_len, -1.0, x_new, 1);

                // Calculate new derivative
                obj->compute_f3D_reduced(x_new, t, f_new);

                // F_new = x_new - x_previous - delta_t * x_dot
                mkl<T>::copy(x_len, x_new, 1, F_new, 1);
                mkl<T>::axpy(x_len, -1.0, it_start, 1, F_new, 1);
                mkl<T>::axpy(x_len, -dt, f_new, 1, F_new, 1);

                // dF = (F_new - F)'
                mkl<T>::copy(x_len, F_new, 1, dF, 1);
                mkl<T>::axpy(x_len, -1.0, F, 1, dF, 1);

                // dx = (x_new - x)';
                mkl<T>::copy(x_len, x_new, 1, dx, 1);
                mkl<T>::axpy(x_len, -1.0, x, 1, dx, 1);

                // J(n+1) = J(n) + (dF - J * dx) * dx' / norm(dx)^2
                // scaling dx and dF by norm(dx)
                nrm = mkl<T>::nrm2(x_len, dx, 1);
                mkl<T>::scal(x_len, 1.0 / nrm, dx, 1);
                mkl<T>::scal(x_len, 1.0 / nrm, dF, 1);
                // y := alpha*A*x + beta*y, dgemv operation
                mkl<T>::gemv(CblasRowMajor, CblasNoTrans, x_len, x_len, -1.0, J, x_len, dx, 1, 1.0,
                             dF, 1);  // using dF to store data
                mkl<T>::ger(CblasRowMajor, x_len, x_len, 1.0, dF, 1, dx, 1, J, x_len);

                // F = F_new; interchanging pointers to avoid copy
                swap_address(F, F_new);
                // x = x_new; interchanging pointers to avoid copy
                swap_address(x, x_new);
            }
            // x_vector_new(n,:) = x;
            // start position of new time step
            curr_time_start_pos = x_vector_new + i * x_len;
            mkl<T>::copy(x_len, x, 1, curr_time_start_pos, 1);
        }

        MKL_free(f_old);
        MKL_free(f_new);
        MKL_free(dx);
        MKL_free(dx_inv);
        MKL_free(df);
        MKL_free(x);
        MKL_free(F);
        MKL_free(F_new);
        MKL_free(dF);
        MKL_free(J);
        MKL_free(J_tmp);
        MKL_free(piv);
        MKL_free(x_new);
    }

    /**
     * Backward-Differences2 (BDF2) with Broyden scheme
     * \param tol tolerance of the Broyden iteration
     * \param max_iter of the Broyden iteration
     * \param num_time_iter number of time steps to perform
     * \param dt timestep
     * \param x_previous previous solution
     * \return new solution
     */
    static void Broyden_PDF2(C* obj, T* x_previous, T* x_vector_new, T dt, size_t num_time_iter,
                             T tol, size_t max_iter)
    {
        size_t alignment = obj->get_alignment();
        size_t x_len = obj->get_solution_dimension();
        T* f_old = (T*)mkl_malloc(sizeof(T) * x_len, alignment);
        T* f_new = (T*)mkl_malloc(sizeof(T) * x_len, alignment);
        T* dx = (T*)mkl_malloc(sizeof(T) * x_len, alignment);
        T* dx_inv = (T*)mkl_malloc(sizeof(T) * x_len, alignment);
        T* df = (T*)mkl_malloc(sizeof(T) * x_len, alignment);
        T* x = (T*)mkl_malloc(sizeof(T) * x_len, alignment);
        T* F = (T*)mkl_malloc(sizeof(T) * x_len, alignment);
        T* F_new = (T*)mkl_malloc(sizeof(T) * x_len, alignment);
        T* dF = (T*)mkl_malloc(sizeof(T) * x_len, alignment);
        T* J = (T*)mkl_calloc(x_len * x_len, sizeof(T), alignment);
        T* J_tmp = (T*)mkl_calloc(x_len * x_len, sizeof(T), alignment);
        lapack_int* piv = (lapack_int*)mkl_malloc(sizeof(lapack_int*) * x_len, alignment);

        T *it_start, *curr_time_start_pos, *prev_prev_pos;
        T* x_new = (T*)mkl_malloc(sizeof(T) * x_len, alignment);
        T* switcher;  // used for interchanging adresses between F and F_new
        T eps = 1e-4;
        double nrm = 0.001;
        T t = 0.0;
        T val = 0.0;
        mkl<T>::copy(x_len, x_previous, 1, x_vector_new, 1);

        size_t i = 1;

        // set the pointer to the start of previous timestep
        it_start = x_vector_new + (i - 1) * x_len;

        // 1. Initialize guess from previous time step
        // f_old = f(t(n-1), x_previous');
        obj->compute_f3D_reduced(it_start, t, f_old);
        // in case the velocity is 0 add nuggets to avoid singular matrices (slow check for
        // improvement) f_old(abs(f_old) < 0.001) = 0.001;
        val = eps * (2 * (1 + rand() % 2) - 3);
        transmutate_elements(f_old, 0, x_len, val, eps);

        // 2. Initial guess from Explicit Euler
        // x = x_previous + dt * f_old';
        // f_new = f(t(n), x');
        t += dt;
        mkl<T>::copy(x_len, it_start, 1, x, 1);
        mkl<T>::axpy(x_len, dt, f_old, 1, x, 1);
        obj->compute_f3D_reduced(x, t, f_new);

        // Initial approximation of the Jacobian
        // dx = x - x_previous;
        mkl<T>::copy(x_len, x, 1, dx, 1);
        mkl<T>::axpy(x_len, -1.0, it_start, 1, dx, 1);
        // df = f_new' - f_old';
        mkl<T>::copy(x_len, f_new, 1, df, 1);
        mkl<T>::axpy(x_len, -1.0, f_old, 1, df, 1);

        // approximate J(x_0)
        // J = eye(length(df)) - delta_t*((1./dx)'*df)';
        mkl<T>::scal(x_len * x_len, 0.0, J, 1);
        diagonal_matrix(J, x_len, 1.0);
        mkl<T>::copy(x_len, dx, 1, dx_inv, 1);
        elementwise_inversion(dx_inv, x_len);
        mkl<T>::ger(CblasRowMajor, x_len, x_len, -dt, df, 1, dx_inv, 1, J, x_len);

        // calculate initial F for stopping condition
        // x_dot = f_new';
        // F = dx - delta_t * x_dot;
        mkl<T>::copy(x_len, dx, 1, F, 1);
        mkl<T>::axpy(x_len, -dt, f_new, 1, F, 1);

        // Broyden's Method
        for (size_t j = 0; j < max_iter; ++j) {
            if (mkl<T>::nrm2(x_len, F, 1) < tol) {
                break;
            }

            // x(i+1) = x(i) - J^(-1)*g(x(i))
            mkl<T>::copy(x_len * x_len, J, 1, J_tmp, 1);
            mkl<T>::getrf(LAPACK_ROW_MAJOR, x_len, x_len, J_tmp, x_len, piv);
            mkl<T>::copy(x_len, F, 1, x_new, 1);
            mkl<T>::getrs(LAPACK_ROW_MAJOR, 'N', x_len, 1, J_tmp, x_len, piv, x_new, 1);
            mkl<T>::axpy(x_len, -1.0, x, 1, x_new, 1);  // result here is -x_new
            mkl<T>::scal(x_len, -1.0, x_new, 1);

            // Calculate new derivative
            obj->compute_f3D_reduced(x_new, t, f_new);

            // F_new = x_new - x_previous - delta_t * x_dot
            mkl<T>::copy(x_len, x_new, 1, F_new, 1);
            mkl<T>::axpy(x_len, -1.0, it_start, 1, F_new, 1);
            mkl<T>::axpy(x_len, -dt, f_new, 1, F_new, 1);

            // dF = (F_new - F)'
            mkl<T>::copy(x_len, F_new, 1, dF, 1);
            mkl<T>::axpy(x_len, -1.0, F, 1, dF, 1);

            // dx = (x_new - x)';
            mkl<T>::copy(x_len, x_new, 1, dx, 1);
            mkl<T>::axpy(x_len, -1.0, x, 1, dx, 1);

            // J(n+1) = J(n) + (dF - J * dx) * dx' / norm(dx)^2
            // scaling dx and dF by norm(dx)
            nrm = mkl<T>::nrm2(x_len, dx, 1);
            mkl<T>::scal(x_len, 1.0 / nrm, dx, 1);
            mkl<T>::scal(x_len, 1.0 / nrm, dF, 1);
            // y := alpha*A*x + beta*y, dgemv operation
            mkl<T>::gemv(CblasRowMajor, CblasNoTrans, x_len, x_len, -1.0, J, x_len, dx, 1, 1.0, dF,
                         1);  // using dF to store data
            mkl<T>::ger(CblasRowMajor, x_len, x_len, 1.0, dF, 1, dx, 1, J, x_len);

            // F = F_new; interchanging pointers to avoid copy
            swap_address(F, F_new);
            // x = x_new; interchanging pointers to avoid copy
            swap_address(x, x_new);
        }
        // x_vector_new(n,:) = x;
        // start position of new time step
        curr_time_start_pos = x_vector_new + i * x_len;
        mkl<T>::copy(x_len, x, 1, curr_time_start_pos, 1);

        if (num_time_iter >= 2) {
            // PDF 2 Method
            for (size_t i = 2; i <= num_time_iter; ++i) {
                it_start = x_vector_new + (i - 1) * x_len;
                prev_prev_pos = x_vector_new + (i - 2) * x_len;

                // 1. Initialize guess from previous time step
                // f_old = f(t(n-1), x_previous');
                obj->compute_f3D_reduced(it_start, t, f_old);
                // in case the velocity is 0 add nuggets to avoid singular matrices (slow check for
                // improvement) f_old(abs(f_old) < 0.001) = 0.001;
                val = eps * (2 * (1 + rand() % 2) - 3);
                transmutate_elements(f_old, 0, x_len, val, eps);

                // 2. Initial guess from Explicit Euler
                // x = x_previous + dt * f_old';
                // f_new = f(t(n), x');
                t += dt;
                mkl<T>::copy(x_len, it_start, 1, x, 1);
                mkl<T>::axpy(x_len, dt, f_old, 1, x, 1);
                obj->compute_f3D_reduced(x, t, f_new);

                // Initial approximation of the Jacobian
                // dx = x - x_previous;
                mkl<T>::copy(x_len, x, 1, dx, 1);
                mkl<T>::axpy(x_len, -1.0, it_start, 1, dx, 1);
                // df = f_new' - f_old';
                mkl<T>::copy(x_len, f_new, 1, df, 1);
                mkl<T>::axpy(x_len, -1.0, f_old, 1, df, 1);

                // approximate J(x_0)
                // J = eye(length(df)) - delta_t*((1./dx)'*df)';
                mkl<T>::scal(x_len * x_len, 0.0, J, 1);
                diagonal_matrix(J, x_len, 1.0);
                mkl<T>::copy(x_len, dx, 1, dx_inv, 1);
                elementwise_inversion(dx_inv, x_len);
                mkl<T>::ger(CblasRowMajor, x_len, x_len, -dt, df, 1, dx_inv, 1, J, x_len);

                // calculate initial F for stopping condition
                // x_dot = f_new';
                // F = x - 4/3 * x_previous + 1/3 * x_previous_previous - 2/3 * delta_t * x_dot;
                mkl<T>::copy(x_len, dx, 1, F, 1);
                mkl<T>::axpy(x_len, -1.0 / 3.0, it_start, 1, F, 1);
                mkl<T>::axpy(x_len, 1.0 / 3.0, prev_prev_pos, 1, F, 1);
                mkl<T>::axpy(x_len, (-2.0 / 3.0) * dt, f_new, 1, F, 1);

                // Broyden's Method
                for (size_t j = 0; j < max_iter; ++j) {
                    if (mkl<T>::nrm2(x_len, F, 1) < tol) {
                        break;
                    }

                    // x(i+1) = x(i) - J^(-1)*g(x(i))
                    mkl<T>::copy(x_len * x_len, J, 1, J_tmp, 1);
                    mkl<T>::getrf(LAPACK_ROW_MAJOR, x_len, x_len, J_tmp, x_len, piv);
                    mkl<T>::copy(x_len, F, 1, x_new, 1);
                    mkl<T>::getrs(LAPACK_ROW_MAJOR, 'N', x_len, 1, J_tmp, x_len, piv, x_new, 1);
                    mkl<T>::axpy(x_len, -1.0, x, 1, x_new, 1);  // result here is -x_new
                    mkl<T>::scal(x_len, -1.0, x_new, 1);

                    // Calculate new derivative
                    obj->compute_f3D_reduced(x_new, t - dt, f_new);

                    // F_new = x_new - 4/3 * x_previous + 1/3 * x_previous_previous - 2/3 * delta_t
                    // * x_dot;
                    mkl<T>::copy(x_len, x_new, 1, F_new, 1);
                    mkl<T>::axpy(x_len, -4.0 / 3.0, it_start, 1, F_new, 1);
                    mkl<T>::axpy(x_len, 1.0 / 3.0, prev_prev_pos, 1, F_new, 1);
                    mkl<T>::axpy(x_len, (-2.0 / 3.0) * dt, f_new, 1, F_new, 1);

                    // dF = (F_new - F)'
                    mkl<T>::copy(x_len, F_new, 1, dF, 1);
                    mkl<T>::axpy(x_len, -1.0, F, 1, dF, 1);

                    // dx = (x_new - x)';
                    mkl<T>::copy(x_len, x_new, 1, dx, 1);
                    mkl<T>::axpy(x_len, -1.0, x, 1, dx, 1);

                    // J(n+1) = J(n) + (dF - J * dx) * dx' / norm(dx)^2
                    // scaling dx and dF by norm(dx)
                    nrm = mkl<T>::nrm2(x_len, dx, 1);
                    mkl<T>::scal(x_len, 1.0 / nrm, dx, 1);
                    mkl<T>::scal(x_len, 1.0 / nrm, dF, 1);
                    // y := alpha*A*x + beta*y, dgemv operation
                    mkl<T>::gemv(CblasRowMajor, CblasNoTrans, x_len, x_len, -1.0, J, x_len, dx, 1,
                                 1.0, dF, 1);  // using dF to store data
                    mkl<T>::ger(CblasRowMajor, x_len, x_len, 1.0, dF, 1, dx, 1, J, x_len);

                    // F = F_new; interchanging pointers to avoid copy
                    swap_address(F, F_new);
                    // x = x_new; interchanging pointers to avoid copy
                    swap_address(x, x_new);
                }
                // x_vector_new(n,:) = x;
                // start position of new time step
                curr_time_start_pos = x_vector_new + i * x_len;
                mkl<T>::copy(x_len, x, 1, curr_time_start_pos, 1);
            }
        }
        MKL_free(f_old);
        MKL_free(f_new);
        MKL_free(dx);
        MKL_free(dx_inv);
        MKL_free(df);
        MKL_free(x);
        MKL_free(F);
        MKL_free(F_new);
        MKL_free(dF);
        MKL_free(J);
        MKL_free(J_tmp);
        MKL_free(piv);
        MKL_free(x_new);
    }

    /**
     * Crank-Nicholson with Broyden scheme
     * \param tol tolerance of the Broyden iteration
     * \param max_iter of the Broyden iteration
     * \param num_time_iter number of time steps to perform
     * \param dt timestep
     * \param x_previous previous solution
     * \return new solution
     */
    static void Broyden_CN(C* obj, T* x_previous, T* x_vector_new, T dt, size_t num_time_iter,
                           T tol, size_t max_iter)
    {
        /*
        method requires that the object using this solver has following public member functions
        1. get_alignment()
        2. get_solution_dimension()
        */
        // std::cout << "Broyden started!\n" << std::endl;

        size_t alignment = obj->get_alignment();
        size_t x_len = obj->get_solution_dimension();
        T* f_old = (T*)mkl_malloc(sizeof(T) * x_len, alignment);
        T* f_new = (T*)mkl_malloc(sizeof(T) * x_len, alignment);
        T* dx = (T*)mkl_malloc(sizeof(T) * x_len, alignment);
        T* dx_inv = (T*)mkl_malloc(sizeof(T) * x_len, alignment);
        T* df = (T*)mkl_malloc(sizeof(T) * x_len, alignment);
        T* x = (T*)mkl_malloc(sizeof(T) * x_len, alignment);
        T* F = (T*)mkl_malloc(sizeof(T) * x_len, alignment);
        T* F_new = (T*)mkl_malloc(sizeof(T) * x_len, alignment);
        T* dF = (T*)mkl_malloc(sizeof(T) * x_len, alignment);
        T* J = (T*)mkl_calloc(x_len * x_len, sizeof(T), alignment);
        T* J_tmp = (T*)mkl_calloc(x_len * x_len, sizeof(T), alignment);
        lapack_int* piv = (lapack_int*)mkl_malloc(sizeof(lapack_int*) * x_len, alignment);

        T *it_start, *curr_time_start_pos;
        T* x_new = (T*)mkl_malloc(sizeof(T) * x_len, alignment);
        T* switcher;  // used for interchanging adresses between F and F_new
        T eps = 1e-4;
        double nrm = 0.001;
        T t = 0.0;
        T val = 0.0;
        mkl<T>::copy(x_len, x_previous, 1, x_vector_new, 1);
        for (size_t i = 1; i <= num_time_iter; ++i) {
            // set the pointer to the start of previous timestep
            // std::cout << "executing iter = " << i << std::endl;
            it_start = x_vector_new + (i - 1) * x_len;

            // 1. Initialize guess from previous time step
            // f_old = f(t(n-1), x_previous');
            obj->compute_f3D_reduced(it_start, t, f_old);

            // in case the velocity is 0 add nuggets to avoid singular matrices (slow check for
            // improvement) f_old(abs(f_old) < 0.0001) = 0.0001;
            val = eps * (2 * (1 + rand() % 2) - 3);
            transmutate_elements(f_old, 0, x_len, val, eps);

            // 2. Initial guess from Explicit Euler
            // x = x_previous + dt * f_old';
            // f_new = f(t(n), x');
            t += dt;
            mkl<T>::copy(x_len, it_start, 1, x, 1);
            mkl<T>::axpy(x_len, dt, f_old, 1, x, 1);
            obj->compute_f3D_reduced(x, t, f_new);

            // Initial approximation of the Jacobian
            // dx = x - x_previous;
            mkl<T>::copy(x_len, x, 1, dx, 1);
            mkl<T>::axpy(x_len, -1.0, it_start, 1, dx, 1);

            // df = f_new' - f_old';
            mkl<T>::copy(x_len, f_new, 1, df, 1);
            mkl<T>::axpy(x_len, -1.0, f_old, 1, df, 1);

            // approximate J(x_0)
            // J = eye(length(df)) - delta_t*0.5 *((1./dx)'*df)';
            mkl<T>::scal(x_len * x_len, 0.0, J, 1);
            diagonal_matrix(J, x_len, 1.0);
            mkl<T>::copy(x_len, dx, 1, dx_inv, 1);
            elementwise_inversion(dx_inv, x_len);
            // write_vector(dx_inv, x_len);
            mkl<T>::ger(CblasRowMajor, x_len, x_len, -dt / 2.0, df, 1, dx_inv, 1, J, x_len);

            // calculate initial F for stopping condition
            // x_dot = f_new';
            // x_dot_previous = f_old';
            // F = dx - delta_t * 0.5 * (x_dot + x_dot_previous);
            mkl<T>::copy(x_len, dx, 1, F, 1);
            mkl<T>::axpy(x_len, -dt / 2.0, f_new, 1, F, 1);
            mkl<T>::axpy(x_len, -dt / 2.0, f_old, 1, F, 1);

            // Broyden's Method
            for (size_t j = 0; j < max_iter; ++j) {
                if (mkl<T>::nrm2(x_len, F, 1) < tol) {
                    break;
                }

                // x(i+1) = x(i) - J^(-1)*g(x(i))
                mkl<T>::copy(x_len * x_len, J, 1, J_tmp, 1);
                mkl<T>::getrf(LAPACK_ROW_MAJOR, x_len, x_len, J_tmp, x_len, piv);
                mkl<T>::copy(x_len, F, 1, x_new, 1);
                mkl<T>::getrs(LAPACK_ROW_MAJOR, 'N', x_len, 1, J_tmp, x_len, piv, x_new, 1);
                mkl<T>::axpy(x_len, -1.0, x, 1, x_new, 1);  // result here is -x_new
                mkl<T>::scal(x_len, -1.0, x_new, 1);

                // Calculate new derivative
                obj->compute_f3D_reduced(x_new, t, f_new);

                // F_new = x_new - x_previous - delta_t * 0.5 * (x_dot + x_dot_previous);
                mkl<T>::copy(x_len, x_new, 1, F_new, 1);
                mkl<T>::axpy(x_len, -1.0, it_start, 1, F_new, 1);
                mkl<T>::axpy(x_len, -dt / 2.0, f_new, 1, F_new, 1);
                mkl<T>::axpy(x_len, -dt / 2.0, f_old, 1, F_new, 1);

                // dF = (F_new - F)'
                mkl<T>::copy(x_len, F_new, 1, dF, 1);
                mkl<T>::axpy(x_len, -1.0, F, 1, dF, 1);

                // dx = (x_new - x)';
                mkl<T>::copy(x_len, x_new, 1, dx, 1);
                mkl<T>::axpy(x_len, -1.0, x, 1, dx, 1);

                // J(n+1) = J(n) + (dF - J * dx) * dx' / norm(dx)^2
                // scaling dx and dF by norm(dx)
                nrm = mkl<T>::nrm2(x_len, dx, 1);
                mkl<T>::scal(x_len, 1.0 / nrm, dx, 1);
                mkl<T>::scal(x_len, 1.0 / nrm, dF, 1);
                // y := alpha*A*x + beta*y, dgemv operation
                mkl<T>::gemv(CblasRowMajor, CblasNoTrans, x_len, x_len, -1.0, J, x_len, dx, 1, 1.0,
                             dF, 1);  // using dF to store data
                mkl<T>::ger(CblasRowMajor, x_len, x_len, 1.0, dF, 1, dx, 1, J, x_len);

                // F = F_new; interchanging pointers to avoid copy
                swap_address(F, F_new);

                // x = x_new; interchanging pointers to avoid copy
                swap_address(x, x_new);
            }
            // x_vector_new(n,:) = x;
            // start position of new time step
            curr_time_start_pos = x_vector_new + i * x_len;
            mkl<T>::copy(x_len, x, 1, curr_time_start_pos, 1);
        }
        // std::cout << "Broyden cleaning!\n" << std::endl;
        MKL_free(f_old);
        MKL_free(f_new);
        MKL_free(dx);
        MKL_free(dx_inv);
        MKL_free(df);
        MKL_free(x);
        MKL_free(F);
        MKL_free(F_new);
        MKL_free(dF);
        MKL_free(J);
        MKL_free(J_tmp);
        MKL_free(piv);
        MKL_free(x_new);

        // std::cout << "Exiting Broyden!\n" << std::endl;
    }

    /**
     * \brief newton loop for the 11 Dof system
     */
    static void Newton(C* obj /**< instance of 11dof class*/,
                       T* force /**< pointer to force vector (not from 11dof class)*/,
                       T* J /**< pointer to Jacobien from 11dofClass*/,
                       T* res /**< residual from 11Dof class*/,
                       T* res_norm /**< norm of the residual from the 11 Dof class*/,
                       T* u_n_p_1 /**< current position vector*/,
                       T* temp /**< pointer to temp vector from 11 dof class for stopping criteria*/
    )
    {
        int count = 0;
        T delta_norm = 1, delta_norm2 = 0;
        lapack_int status;
        obj->calcResidual(force);

        do {
            count++;
            obj->constructJacobien();

            // get J=LL^T
            status = mkl<T>::potrf(LAPACK_ROW_MAJOR, 'L', Constants::DOF, J, Constants::DOF);
            check_status<lapack_int>(status);
            // residual=J\residual -> residual = delta
            mkl<T>::potrs(LAPACK_ROW_MAJOR, 'L', Constants::DOF, 1, J, Constants::DOF, res, 1);
            // u_n_p_1 -= delta
            mkl<T>::axpy(Constants::DOF, -1, res, 1, u_n_p_1, 1);
            delta_norm = mkl<T>::nrm2(Constants::DOF, res, 1);

            // update all the matrices
            obj->updateSystem();
            // calculate the newton function
            obj->calcResidual(force);
            // copy the new residual to a temp to check if newton has converged without losing
            // residual
            mkl<T>::copy(Constants::DOF, res, 1, temp, 1);
            mkl<T>::potrs(LAPACK_ROW_MAJOR, 'L', Constants::DOF, 1, J, Constants::DOF, temp, 1);
            delta_norm2 = mkl<T>::nrm2(Constants::DOF, temp, 1);
        } while (*res_norm > Constants::TOLERANCE && delta_norm2 < delta_norm && count < 10);
    }

    /**
     * explicit Runge-Kutta-4
     * \param num_time_iter number of time steps to perform
     * \param dt timestep
     * \param x_previous previous solution
     * \return new solution
     */
    static void RK4(C* obj, T* x_previous, T* x_vector_new, T dt, size_t num_time_iter, T tol,
                    size_t max_iter)
    {
        size_t alignment = obj->get_alignment();
        size_t x_len = obj->get_solution_dimension();
        T* f_old = (T*)mkl_malloc(sizeof(T) * x_len, alignment);
        T* k1 = (T*)mkl_malloc(sizeof(T) * x_len, alignment);
        T* k2 = (T*)mkl_malloc(sizeof(T) * x_len, alignment);
        T* k3 = (T*)mkl_malloc(sizeof(T) * x_len, alignment);
        T* k4 = (T*)mkl_malloc(sizeof(T) * x_len, alignment);
        T* x = (T*)mkl_malloc(sizeof(T) * x_len, alignment);

        T *it_start, *curr_time_start_pos;
        T* switcher;  // used for interchanging adresses between F and F_new
        T t = 0.0;
        T coeff1, coeff2, coeff3, coeff4;
        coeff1 = 1.0 / 6.0;
        coeff2 = 2.0 / 6.0;
        coeff3 = 2.0 / 6.0;
        coeff4 = 1.0 / 6.0;
        mkl<T>::copy(x_len, x_previous, 1, x_vector_new, 1);
        for (size_t i = 1; i <= num_time_iter; ++i) {
            it_start = x_vector_new + (i - 1) * x_len;
            curr_time_start_pos = x_vector_new + (i)*x_len;
            // k1 = delta_t * f(t_prev, x_previous')';
            obj->compute_f3D_reduced(it_start, t, f_old);
            mkl<T>::copy(x_len, f_old, 1, k1, 1);
            mkl<T>::scal(x_len, dt, k1, 1);

            // k2 = delta_t * f(t_prev + delta_t/2, x_previous' + k1'/2)';
            mkl<T>::copy(x_len, k1, 1, x, 1);
            mkl<T>::scal(x_len, 1.0 / 2.0, x, 1);
            mkl<T>::axpy(x_len, 1.0, it_start, 1, x, 1);
            obj->compute_f3D_reduced(x, t + dt / 2.0, f_old);
            mkl<T>::copy(x_len, f_old, 1, k2, 1);
            mkl<T>::scal(x_len, dt, k2, 1);

            // k3 = delta_t * f(t_prev + delta_t/2, x_previous' + k2'/2)';
            mkl<T>::copy(x_len, k2, 1, x, 1);
            mkl<T>::scal(x_len, 1.0 / 2.0, x, 1);
            mkl<T>::axpy(x_len, 1.0, it_start, 1, x, 1);
            obj->compute_f3D_reduced(x, t + dt / 2.0, f_old);
            mkl<T>::copy(x_len, f_old, 1, k3, 1);
            mkl<T>::scal(x_len, dt, k3, 1);

            // k4 = delta_t * f(t_prev + delta_t, x_previous' + k3')';
            mkl<T>::copy(x_len, k3, 1, x, 1);
            mkl<T>::axpy(x_len, 1.0, it_start, 1, x, 1);
            obj->compute_f3D_reduced(x, t + dt, f_old);
            mkl<T>::copy(x_len, f_old, 1, k4, 1);
            mkl<T>::scal(x_len, dt, k4, 1);

            // x_vector_new(n,:) = x_previous + 1/6 * (k1 + 2*k2 + 2*k3 + k4);
            mkl<T>::copy(x_len, it_start, 1, curr_time_start_pos, 1);
            mkl<T>::axpy(x_len, coeff1, k1, 1, curr_time_start_pos, 1);
            mkl<T>::axpy(x_len, coeff2, k2, 1, curr_time_start_pos, 1);
            mkl<T>::axpy(x_len, coeff3, k3, 1, curr_time_start_pos, 1);
            mkl<T>::axpy(x_len, coeff4, k4, 1, curr_time_start_pos, 1);
            t += dt;
        }
        mkl_free(f_old);
        mkl_free(k1);
        mkl_free(k2);
        mkl_free(k3);
        mkl_free(k4);
        mkl_free(x);
    }

    /**
     * linear backward Euler
     * \param num_time_iter number of time steps to perform
     * \param dt timestep
     * \param x_previous previous solution
     * \return new solution
     */
    static void Linear_Backward_Euler(T* A, T* B, T* C, T* x_prev, T* x_prev_prev, T* b, T* x,
                                      size_t dim)
    {
        /*
         * This works for only symmetric positive definite A the provided matrix A would be
         * overwritten and results stored in x computes the backward euler step of following form:
         * Ax = B*x_prev + C*x_prev_prev + b
         * A, B, C are coefficient of the euler formation matrix not implemented
         */
        lapack_int status;
        // get A=LL^T
        status = mkl<T>::potrf(LAPACK_ROW_MAJOR, 'L', dim, A, dim);

        check_status<lapack_int>(status);
        // u_n_p_1 = B * u_n
        mkl<T>::gemv(CblasRowMajor, CblasNoTrans, dim, dim, 1.0, B, dim, x_prev, 1, 0, x, 1);
        // u_n_p_1 += C * u_n_m_1 <=> u_n_p_1 += ((1/(h*h))*M)*(-u_n_m_1)
        mkl<T>::gemv(CblasRowMajor, CblasNoTrans, dim, dim, 1.0, C, dim, x_prev_prev, 1, 1, x, 1);
        // u_n_p_1 += x <=> u_n_p_1 += f_n_p_1
        mkl<T>::axpy(dim, 1, b, 1, x, 1);
        // u_n_p_1=A\u_n_p_1
        mkl<T>::potrs(LAPACK_ROW_MAJOR, 'L', dim, 1, A, dim, x, 1);
    }

    /**
     * linear backward Euler - C (Mass_matrix) diagonal matrix, considered as vector
     * \param num_time_iter number of time steps to perform
     * \param dt timestep
     * \param x_previous previous solution
     * \return new solution
     */
    static void Linear_Backward_Euler_diag(T* A, T* B, T* C, T* x_prev, T* x_prev_prev, T* b, T* x,
                                           size_t dim)
    {
        /*
         * This works for only symmetric positive definite A the provided matrix A would be
         * overwritten and results stored in x computes the backward euler step of following form:
         * Ax = B*x_prev + C*x_prev_prev + b
         * A, B, C are coefficient of the euler formation matrix not implemented
         */
        lapack_int status;
        // get A=LL^T
        status = mkl<T>::potrf(LAPACK_ROW_MAJOR, 'L', dim, A, dim);
        check_status<lapack_int>(status);
        // u_n_p_1 += B * u_n
        mkl<T>::gemv(CblasRowMajor, CblasNoTrans, dim, dim, 1.0, B, dim, x_prev, 1, 0, x, 1);

        // u_n_p_1 += C * u_n_m_1 <=> u_n_p_1 += ((1/(h*h))*M)*(-u_n_m_1)
        // cblas_dgemv(CblasRowMajor, CblasNoTrans, dim, dim, 1.0, C, dim, x_prev_prev, 1, 1, x, 1);
        T vec_tmp[11];
        mkl<T>::vMul(dim, C, x_prev_prev, vec_tmp);
        mkl<T>::axpy(dim, 1, vec_tmp, 1, x, 1);

        // u_n_p_1 += x <=> u_n_p_1 += f_n_p_1
        mkl<T>::axpy(dim, 1, b, 1, x, 1);
        // u_n_p_1=A\u_n_p_1
        mkl<T>::potrs(LAPACK_ROW_MAJOR, 'L', dim, 1, A, dim, x, 1);
    }

    /**
     * 2nd order Stoermer-Verlet algorithm to update the position of one scalar
     * \param x position
     * \param v velocity
     * \param F force
     * \param delta_t timestep
     * \param mass
     * \return x position
     */
    static void Stoermer_Verlet_Position(T& x, T& v, T& F, T& delta_t, T& mass)
    {
        x += delta_t * v + delta_t * delta_t / (2 * mass) * F;
    }

    /**
     * 2nd order Stoermer-Verlet algorithm to update the velocity of one scalar
     * \param v velocity
     * \param F force
     * \param F_new force at the following time step
     * \param delta_t timestep
     * \param mass
     * \return v velocity
     */
    static void Stoermer_Verlet_Velocity(T& v, T& F, T& F_new, T& delta_t, T& mass)
    {
        v += delta_t / (2 * mass) * (F + F_new);
    }

    static void Linear_BDF2(T* A, T* B, T* C, T* D, T* E, T* x_n, T* x_n_m_1, T* x_n_m_2,
                            T* x_n_m_3, T* b, T* x_n_p_1, size_t dim)
    {
        /*
         * This works for only symmetric positive definite A the provided matrix A would be
         * overwritten and results stored in x computes the backward euler step of following form:
         * A*x_n_p_1 = B*x_n + C*x_n_m_1 + D*x_n_m_2 + E*x_n_m_3 + b
         * A, B, C are coefficient of the euler formation matrix not implemented
         */
        lapack_int status;
        // get A=LL^T
        status = mkl<T>::potrf(LAPACK_ROW_MAJOR, 'L', dim, A, dim);
        check_status<lapack_int>(status);
        // u_n_p_1 += B * u_n
        mkl<T>::gemv(CblasRowMajor, CblasNoTrans, dim, dim, 1.0, B, dim, x_n, 1, 0, x_n_p_1, 1);
        // u_n_p_1 += C * u_n_m_1 <=> u_n_p_1 += ((1/(h*h))*M)*(-u_n_m_1)
        mkl<T>::gemv(CblasRowMajor, CblasNoTrans, dim, dim, 1.0, C, dim, x_n_m_1, 1, 1, x_n_p_1, 1);
        mkl<T>::gemv(CblasRowMajor, CblasNoTrans, dim, dim, 1.0, D, dim, x_n_m_2, 1, 1, x_n_p_1, 1);
        mkl<T>::gemv(CblasRowMajor, CblasNoTrans, dim, dim, 1.0, E, dim, x_n_m_3, 1, 1, x_n_p_1, 1);

        // u_n_p_1 += x <=> u_n_p_1 += f_n_p_1
        mkl<T>::axpy(dim, 1, b, 1, x_n_p_1, 1);
        // u_n_p_1=A\u_n_p_1
        mkl<T>::potrs(LAPACK_ROW_MAJOR, 'L', dim, 1, A, dim, x_n_p_1, 1);
    }
};

}  // namespace MathLibrary
}  // namespace EVAA
