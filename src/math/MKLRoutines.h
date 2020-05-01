// TODO: copyright header

#pragma once

// Force enable gemm for now.
#define USE_GEMM

#include <mkl.h>

#include <string>

#include "Constants.h"

namespace EVAA {
namespace Math {

// Until MKL does not get an alternative for this project's purpose, keep the
// methods in the EVAA::Math namespace to write cleaner code. namespace MKL {

/**
 * Template magic to invoke one of two calls, based on matching argument types.
 * Used for all float/double alternatives in MKL.
 *
 * TODO: Possible improvements:
 * - make use of std::forward to stop writing all arguments. The downside is
 * that you might not get autocomplete features.
 * - Copy EIGEN's approach: declare all [sb] versions of the functions using the
 * same base name since C++ allows overloads, e.g. define apxy(float) and
 * apxy(double) letting the type to be inferred when used.
 * */
namespace {

template <typename F1, typename F2, typename... Args>
inline auto invokeFD(F1 f1, F2 f2, Args&&... args) -> decltype(f1(args...)) {
    return f1(args...);
}

template <typename F1, typename F2, typename... Args>
inline auto invokeFD(F1 f1, F2 f2, Args&&... args) -> decltype(f2(args...)) {
    return f2(args...);
}

}  // namespace

/*
 * Reduce verbosity for template function definition below.
 * NOTE: std::is_floating_point does not filter out long double.
 */
#define TEMPLATE_FLOAT_TYPE template <typename T, typename = std::enable_if_t<std::is_floating_point<T>::value>>

/** Returns information about MKL. */
std::string getInfo();

/**
 * \brief Allocates a buffer memory for count elements.
 * Wraps a call to mkl_malloc.
 * \param count Number of elments.
 * \return The newly allocated buffer.
 */
template <typename T>
T* malloc(size_t count) {
    if (count == 0) {
        throw std::invalid_argument("count cannot be 0 for malloc");
    }
    // TODO: signal allocation error.
    // See https://software.intel.com/en-us/mkl-developer-reference-c-mkl-malloc
    T* ptr = (T*)mkl_malloc(count * sizeof(T), Constants::ALIGNMENT);
    if (ptr == nullptr) {
        throw "Memory allocation Error in malloc!";
    }
    return ptr;
}

/**
 * \brief Allocates a buffer memory for count elements and initializes them to
 * 0. Wraps a call to mkl_calloc. \param count Number of elments. \return The
 * newly allocated buffer.
 */
template <typename T>
T* calloc(size_t count) {
    if (count == 0) {
        throw std::invalid_argument("count cannot be 0 for calloc");
    }
    // TODO: signal allocation error.
    // See https://software.intel.com/en-us/mkl-developer-reference-c-mkl-calloc
    T* ptr = (T*)mkl_calloc(count, sizeof(T), Constants::ALIGNMENT);
    if (ptr == nullptr) {
        throw "Memory allocation Error in calloc!";
    }
    return ptr;
}

/**
 * \brief Frees buffer memory previously allocated with Math::malloc or
 * Math::calloc. Wraps a call to mkl_free. \param buffer The buffer to be freed.
 */
template <typename T>
void free(T*& buffer) {
    if (buffer == nullptr) {
        throw std::invalid_argument("Attempt to free a nullptr");
    }
    mkl_free(buffer);
    buffer = nullptr;
}

// BLAS 1

/** Wraps cblas_saxpy and cblas_daxpy. */
TEMPLATE_FLOAT_TYPE
inline void axpy(const MKL_INT n, const T a, const T* x, const MKL_INT incx, T* y, const MKL_INT incy) { invokeFD(cblas_saxpy, cblas_daxpy, n, a, x, incx, y, incy); }

/** Wraps cblas_scopy and cblas_dcopy. */
TEMPLATE_FLOAT_TYPE
inline void copy(const MKL_INT n, const T* x, const MKL_INT incx, T* y, const MKL_INT incy) { invokeFD(cblas_scopy, cblas_dcopy, n, x, incx, y, incy); }

/** Wraps cblas_sdot and cblas_ddot. */
TEMPLATE_FLOAT_TYPE
inline T dot(const MKL_INT n, const T* x, const MKL_INT incx, const T* y, const MKL_INT incy) { return invokeFD(cblas_sdot, cblas_ddot, n, x, incx, y, incy); }

/** Wraps cblas_snrm2 and cblas_dnrm2. */
TEMPLATE_FLOAT_TYPE
inline T nrm2(const MKL_INT n, const T* x, const MKL_INT incx) { return invokeFD(cblas_snrm2, cblas_dnrm2, n, x, incx); }

/** Wraps cblas_sscal and cblas_dscal. */
TEMPLATE_FLOAT_TYPE
inline void scal(const MKL_INT n, const T a, T* x, const MKL_INT incx) { invokeFD(cblas_sscal, cblas_dscal, n, a, x, incx); }

/** Wraps cblas_sswap and cblas_dswap. */
TEMPLATE_FLOAT_TYPE
inline void swap(const MKL_INT n, T* x, const MKL_INT incx, T* y, const MKL_INT incy) { invokeFD(cblas_sswap, cblas_dswap, n, x, incx, y, incy); }

/** Wraps cblas_srot and cblas_drot. */
TEMPLATE_FLOAT_TYPE
inline void rot(const MKL_INT n, T* x, const MKL_INT incx, T* y, const MKL_INT incy, const T c, const T s) { invokeFD(cblas_srot, cblas_drot, n, x, incx, y, incy, c, s); }

// BLAS 2

/** Wraps cblas_sgbmv and cblas_dgbmv. */
TEMPLATE_FLOAT_TYPE
inline void gbmv(const CBLAS_LAYOUT Layout, const CBLAS_TRANSPOSE trans, const MKL_INT m, const MKL_INT n, const MKL_INT kl, const MKL_INT ku, const T alpha, const T* a, const MKL_INT lda, const T* x, const MKL_INT incx, const T beta, T* y, const MKL_INT incy) { invokeFD(cblas_sgbmv, cblas_dgbmv, Layout, trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy); }

/** Wraps cblas_sgemv and cblas_dgemv. */
TEMPLATE_FLOAT_TYPE
inline void gemv(const CBLAS_LAYOUT Layout, const CBLAS_TRANSPOSE trans, const MKL_INT m, const MKL_INT n, const T alpha, const T* a, const MKL_INT lda, const T* x, const MKL_INT incx, const T beta, T* y, const MKL_INT incy) { invokeFD(cblas_sgemv, cblas_dgemv, Layout, trans, m, n, alpha, a, lda, x, incx, beta, y, incy); }

/** Wraps cblas_sger and cblas_dger. */
TEMPLATE_FLOAT_TYPE
inline void ger(const CBLAS_LAYOUT Layout, const MKL_INT m, const MKL_INT n, const T alpha, const T* x, const MKL_INT incx, const T* y, const MKL_INT incy, T* a, const MKL_INT lda) { invokeFD(cblas_sger, cblas_dger, Layout, m, n, alpha, x, incx, y, incy, a, lda); }

/** Wraps cblas_ssbmv and cblas_dsbmv. */
TEMPLATE_FLOAT_TYPE
inline void sbmv(const CBLAS_LAYOUT Layout, const CBLAS_UPLO uplo, const MKL_INT n, const MKL_INT k, const T alpha, const T* a, const MKL_INT lda, const T* x, const MKL_INT incx, const T beta, T* y, const MKL_INT incy) { invokeFD(cblas_ssbmv, cblas_dsbmv, Layout, uplo, n, k, alpha, a, lda, x, incx, beta, y, incy); }

/** Wraps cblas_simatcopy and cblas_dimatcopy. */
TEMPLATE_FLOAT_TYPE
inline void imatcopy(const char ordering, const char trans, size_t rows, size_t cols, const T alpha, T* AB, size_t lda, size_t ldb) { invokeFD(mkl_simatcopy, mkl_dimatcopy, ordering, trans, rows, cols, alpha, AB, lda, ldb); }

/** Wraps cblas_ssymv and cblas_dsymv. */
TEMPLATE_FLOAT_TYPE
inline void symv(const CBLAS_LAYOUT Layout, const CBLAS_UPLO uplo, const MKL_INT n, const T alpha, const T* a, const MKL_INT lda, const T* x, const MKL_INT incx, const T beta, T* y, const MKL_INT incy) { invokeFD(cblas_ssymv, cblas_dsymv, Layout, uplo, n, alpha, a, lda, x, incx, beta, y, incy); }

/** Wraps cblas_ssyr and cblas_dsyr. */
TEMPLATE_FLOAT_TYPE
inline void syr(const CBLAS_LAYOUT Layout, const CBLAS_UPLO uplo, const MKL_INT n, const T alpha, const T* x, const MKL_INT incx, T* a, const MKL_INT lda) { invokeFD(cblas_ssyr, cblas_dsyr, Layout, uplo, n, alpha, x, incx, a, lda); }

/** Wraps cblas_ssyr2 and cblas_dsyr2. */
TEMPLATE_FLOAT_TYPE
inline void syr2(const CBLAS_LAYOUT Layout, const CBLAS_UPLO uplo, const MKL_INT n, const T alpha, const T* x, const MKL_INT incx, const T* y, const MKL_INT incy, T* a, const MKL_INT lda) { invokeFD(cblas_ssyr2, cblas_dsyr2, Layout, uplo, n, alpha, x, incx, y, incy, a, lda); }

/** Wraps cblas_stbmv and cblas_dtbmv. */
TEMPLATE_FLOAT_TYPE
inline void tbmv(const CBLAS_LAYOUT Layout, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE trans, const CBLAS_DIAG diag, const MKL_INT n, const MKL_INT k, const T* a, const MKL_INT lda, T* x, const MKL_INT incx) { invokeFD(cblas_stbmv, cblas_dtbmv, Layout, uplo, trans, diag, n, k, a, lda, x, incx); }

// BLAS 3

/** Wraps cblas_sgemm and cblas_dgemm. */
TEMPLATE_FLOAT_TYPE
inline void gemm(const CBLAS_LAYOUT Layout, const CBLAS_TRANSPOSE transa, const CBLAS_TRANSPOSE transb, const MKL_INT m, const MKL_INT n, const MKL_INT k, const T alpha, const T* a, const MKL_INT lda, const T* b, const MKL_INT ldb, const T beta, T* c, const MKL_INT ldc) { invokeFD(cblas_sgemm, cblas_dgemm, Layout, transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc); }

/** Wraps cblas_ssymm and cblas_dsymm. */
TEMPLATE_FLOAT_TYPE
inline void symm(const CBLAS_LAYOUT Layout, const CBLAS_SIDE side, const CBLAS_UPLO uplo, const MKL_INT m, const MKL_INT n, const T alpha, const T* a, const MKL_INT lda, const T* b, const MKL_INT ldb, const T beta, T* c, const MKL_INT ldc) { invokeFD(cblas_ssymm, cblas_dsymm, Layout, side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc); }

/** Wraps cblas_ssyrk and cblas_dsyrk. */
TEMPLATE_FLOAT_TYPE
inline void syrk(const CBLAS_LAYOUT Layout, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE trans, const MKL_INT n, const MKL_INT k, const T alpha, const T* a, const MKL_INT lda, const T beta, T* c, const MKL_INT ldc) { invokeFD(cblas_ssyrk, cblas_dsyrk, Layout, uplo, trans, n, k, alpha, a, lda, beta, c, ldc); }

/** Wraps cblas_strmm and cblas_dtrmm. */
TEMPLATE_FLOAT_TYPE
inline void trmm(const CBLAS_LAYOUT Layout, const CBLAS_SIDE side, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE transa, const CBLAS_DIAG diag, const MKL_INT m, const MKL_INT n, const T alpha, const T* a, const MKL_INT lda, T* b, const MKL_INT ldb) { invokeFD(cblas_strmm, cblas_dtrmm, Layout, side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb); }

// LAPACK

/** Wraps LAPACKE_sgetrf and LAPACKE_dgetrf. */
TEMPLATE_FLOAT_TYPE
inline lapack_int getrf(int matrix_layout, lapack_int m, lapack_int n, T* a, lapack_int lda, lapack_int* ipiv) { return invokeFD(LAPACKE_sgetrf, LAPACKE_dgetrf, matrix_layout, m, n, a, lda, ipiv); }

/** Wraps LAPACKE_sgetrs and LAPACKE_dgetrs. */
TEMPLATE_FLOAT_TYPE
inline lapack_int getrs(int matrix_layout, char trans, lapack_int n, lapack_int nrhs, const T* a, lapack_int lda, const lapack_int* ipiv, T* b, lapack_int ldb) { return invokeFD(LAPACKE_sgetrs, LAPACKE_dgetrs, matrix_layout, trans, n, nrhs, a, lda, ipiv, b, ldb); }

/** Wraps LAPACKE_sgecon and LAPACKE_dgecon. */
TEMPLATE_FLOAT_TYPE
inline lapack_int gecon(int matrix_layout, char norm, lapack_int n, const T* a, lapack_int lda, T anorm, T* rcond) { return invokeFD(LAPACKE_sgecon, LAPACKE_dgecon, matrix_layout, norm, n, a, lda, anorm, rcond); }

/** Wraps LAPACKE_sgbtrf and LAPACKE_dgbtrf. */
TEMPLATE_FLOAT_TYPE
inline lapack_int gbtrf(int matrix_layout, lapack_int m, lapack_int n, lapack_int kl, lapack_int ku, T* ab, lapack_int ldab, lapack_int* ipiv) { return invokeFD(LAPACKE_sgbtrf, LAPACKE_dgbtrf, matrix_layout, m, n, kl, ku, ab, ldab, ipiv); }

/** Wraps LAPACKE_sgbtrs and LAPACKE_dgbtrs. */
TEMPLATE_FLOAT_TYPE
inline lapack_int gbtrs(int matrix_layout, char trans, lapack_int n, lapack_int kl, lapack_int ku, lapack_int nrhs, const T* ab, lapack_int ldab, const lapack_int* ipiv, T* b, lapack_int ldb) { return invokeFD(LAPACKE_sgbtrs, LAPACKE_dgbtrs, matrix_layout, trans, n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb); }

/** Wraps LAPACKE_sgbcon and LAPACKE_dgbcon. */
TEMPLATE_FLOAT_TYPE
inline lapack_int gbcon(int matrix_layout, char norm, lapack_int n, lapack_int kl, lapack_int ku, const T* ab, lapack_int ldab, const lapack_int* ipiv, T anorm, T* rcond) { return invokeFD(LAPACKE_sgbcon, LAPACKE_dgbcon, matrix_layout, norm, n, kl, ku, ab, ldab, ipiv, anorm, rcond); }

/** Wraps LAPACKE_spotrf and LAPACKE_dpotrf. */
TEMPLATE_FLOAT_TYPE
inline lapack_int potrf(int matrix_layout, char uplo, lapack_int n, T* a, lapack_int lda) { return invokeFD(LAPACKE_spotrf, LAPACKE_dpotrf, matrix_layout, uplo, n, a, lda); }

/**
 * Checks if a status is an error for Math::potrf.
 * \param status the status to be checked.
 * \throw std::domain_error If the status reports a potrf failure.
 * \note: consider putting the status check inside Math::potrf.
 */
void potrfCheckStatus(lapack_int status);

/** Wraps LAPACKE_spotrs and LAPACKE_dpotrs. */
TEMPLATE_FLOAT_TYPE
inline lapack_int potrs(int matrix_layout, char uplo, lapack_int n, lapack_int nrhs, const T* a, lapack_int lda, T* b, lapack_int ldb) { return invokeFD(LAPACKE_spotrs, LAPACKE_dpotrs, matrix_layout, uplo, n, nrhs, a, lda, b, ldb); }

/** Wraps LAPACKE_spocon and LAPACKE_dpocon. */
TEMPLATE_FLOAT_TYPE
inline lapack_int pocon(int matrix_layout, char uplo, lapack_int n, const T* a, lapack_int lda, T anorm, T* rcond) { return invokeFD(LAPACKE_spocon, LAPACKE_dpocon, matrix_layout, uplo, n, a, lda, anorm, rcond); }

/** Wraps LAPACKE_spttrf and LAPACKE_dpttrf. */
TEMPLATE_FLOAT_TYPE
inline lapack_int pttrf(lapack_int n, T* d, T* e) { return invokeFD(LAPACKE_spttrf, LAPACKE_dpttrf, n, d, e); }

/** Wraps LAPACKE_spttrs and LAPACKE_dpttrs. */
TEMPLATE_FLOAT_TYPE
inline lapack_int pttrs(int matrix_layout, lapack_int n, lapack_int nrhs, const T* d, const T* e, T* b, lapack_int ldb) { return invokeFD(LAPACKE_spttrs, LAPACKE_dpttrs, matrix_layout, n, nrhs, d, e, b, ldb); }

/** Wraps LAPACKE_sptcon and LAPACKE_dptcon. */
TEMPLATE_FLOAT_TYPE
inline lapack_int ptcon(lapack_int n, const T* d, const T* e, T anorm, T* rcond) { return invokeFD(LAPACKE_sptcon, LAPACKE_dptcon, n, d, e, anorm, rcond); }

// Orthogonal Factorizations

/** Wraps LAPACKE_sgeqrf and LAPACKE_dgeqrf. */
TEMPLATE_FLOAT_TYPE
inline lapack_int geqrf(int matrix_layout, lapack_int m, lapack_int n, T* a, lapack_int lda, T* tau) { return invokeFD(LAPACKE_sgeqrf, LAPACKE_dgeqrf, matrix_layout, m, n, a, lda, tau); }

/** Wraps LAPACKE_sgeqpf and LAPACKE_dgeqpf. */
TEMPLATE_FLOAT_TYPE
inline lapack_int geqpf(int matrix_layout, lapack_int m, lapack_int n, T* a, lapack_int lda, lapack_int* jpvt, T* tau) { return invokeFD(LAPACKE_sgeqpf, LAPACKE_dgeqpf, matrix_layout, m, n, a, lda, jpvt, tau); }

/** Wraps LAPACKE_sorgqr and LAPACKE_dorgqr. */
TEMPLATE_FLOAT_TYPE
inline lapack_int orgqr(int matrix_layout, lapack_int m, lapack_int n, lapack_int k, T* a, lapack_int lda, const T* tau) { return invokeFD(LAPACKE_sorgqr, LAPACKE_dorgqr, matrix_layout, m, n, k, a, lda, tau); }

/** Wraps LAPACKE_sormqr and LAPACKE_dormqr. */
TEMPLATE_FLOAT_TYPE
inline lapack_int ormqr(int matrix_layout, char side, char trans, lapack_int m, lapack_int n, lapack_int k, const T* a, lapack_int lda, const T* tau, T* c, lapack_int ldc) { return invokeFD(LAPACKE_sormqr, LAPACKE_dormqr, matrix_layout, side, trans, m, n, k, a, lda, tau, c, ldc); }

/** Wraps LAPACKE_sgebrd and LAPACKE_dgebrd. */
TEMPLATE_FLOAT_TYPE
inline lapack_int gebrd(int matrix_layout, lapack_int m, lapack_int n, T* a, lapack_int lda, T* d, T* e, T* tauq, T* taup) { return invokeFD(LAPACKE_sgebrd, LAPACKE_dgebrd, matrix_layout, m, n, a, lda, d, e, tauq, taup); }

/** Wraps LAPACKE_sbdsqr and LAPACKE_dbdsqr. */
TEMPLATE_FLOAT_TYPE
inline lapack_int bdsqr(int matrix_layout, char uplo, lapack_int n, lapack_int ncvt, lapack_int nru, lapack_int ncc, T* d, T* e, T* vt, lapack_int ldvt, T* u, lapack_int ldu, T* c, lapack_int ldc) { return invokeFD(LAPACKE_sbdsqr, LAPACKE_dbdsqr, matrix_layout, uplo, n, ncvt, nru, ncc, d, e, vt, ldvt, u, ldu, c, ldc); }

/** Wraps LAPACKE_sbdsdc and LAPACKE_dbdsdc. */
TEMPLATE_FLOAT_TYPE
inline lapack_int bdsdc(int matrix_layout, char uplo, char compq, lapack_int n, T* d, T* e, T* u, lapack_int ldu, T* vt, lapack_int ldvt, T* q, lapack_int* iq) { return invokeFD(LAPACKE_sbdsdc, LAPACKE_dbdsdc, matrix_layout, uplo, compq, n, d, e, u, ldu, vt, ldvt, q, iq); }

/** Wraps LAPACKE_ssytrd and LAPACKE_dsytrd. */
TEMPLATE_FLOAT_TYPE
inline lapack_int sytrd(int matrix_layout, char uplo, lapack_int n, T* a, lapack_int lda, T* d, T* e, T* tau) { invokeFD(LAPACKE_ssytrd, LAPACKE_dsytrd, matrix_layout, uplo, n, a, lda, d, e, tau); }

/** Wraps LAPACKE_ssbtrd and LAPACKE_dsbtrd. */
TEMPLATE_FLOAT_TYPE
inline lapack_int sbtrd(int matrix_layout, char vect, char uplo, lapack_int n, lapack_int kd, T* ab, lapack_int ldab, T* d, T* e, T* q, lapack_int ldq) { return invokeFD(LAPACKE_ssbtrd, LAPACKE_dsbtrd, matrix_layout, vect, uplo, n, kd, ab, ldab, d, e, q, ldq); }

/** Wraps LAPACKE_sorgtr and LAPACKE_dorgtr. */
TEMPLATE_FLOAT_TYPE
inline lapack_int orgtr(int matrix_layout, char uplo, lapack_int n, T* a, lapack_int lda, const T* tau) { return invokeFD(LAPACKE_sorgtr, LAPACKE_dorgtr, matrix_layout, uplo, n, a, lda, tau); }

/** Wraps LAPACKE_sormtr and LAPACKE_dormtr. */
TEMPLATE_FLOAT_TYPE
inline lapack_int ormtr(int matrix_layout, char side, char uplo, char trans, lapack_int m, lapack_int n, const T* a, lapack_int lda, const T* tau, T* c, lapack_int ldc) { return invokeFD(LAPACKE_sormtr, LAPACKE_dormtr, matrix_layout, side, uplo, trans, m, n, a, lda, tau, c, ldc); }

/** Wraps LAPACKE_ssteqr and LAPACKE_dsteqr. */
TEMPLATE_FLOAT_TYPE
inline lapack_int steqr(int matrix_layout, char compz, lapack_int n, T* d, T* e, T* z, lapack_int ldz) { return invokeFD(LAPACKE_ssteqr, LAPACKE_dsteqr, matrix_layout, compz, n, d, e, z, ldz); }

/** Wraps LAPACKE_sstedc and LAPACKE_dstedc. */
TEMPLATE_FLOAT_TYPE
inline lapack_int stedc(int matrix_layout, char compz, lapack_int n, T* d, T* e, T* z, lapack_int ldz) { return invokeFD(LAPACKE_sstedc, LAPACKE_dstedc, matrix_layout, compz, n, d, e, z, ldz); }

/** Wraps LAPACKE_spteqr and LAPACKE_dpteqr. */
TEMPLATE_FLOAT_TYPE
inline lapack_int pteqr(int matrix_layout, char compz, lapack_int n, T* d, T* e, T* z, lapack_int ldz) { return invokeFD(LAPACKE_spteqr, LAPACKE_dpteqr, matrix_layout, compz, n, d, e, z, ldz); }

/** Wraps LAPACKE_sggsvp and LAPACKE_dggsvp. */
TEMPLATE_FLOAT_TYPE
inline lapack_int ggsvp(int matrix_layout, char jobu, char jobv, char jobq, lapack_int m, lapack_int p, lapack_int n, T* a, lapack_int lda, T* b, lapack_int ldb, T tola, T tolb, lapack_int* k, lapack_int* l, T* u, lapack_int ldu, T* v, lapack_int ldv, T* q, lapack_int ldq) { return invokeFD(LAPACKE_sggsvp, LAPACKE_dggsvp, matrix_layout, jobu, jobv, jobq, m, p, n, a, lda, b, ldb, tola, tolb, k, l, u, ldu, v, ldv, q, ldq); }

/** Wraps LAPACKE_sggsvp3 and LAPACKE_dggsvp3. */
TEMPLATE_FLOAT_TYPE
inline lapack_int ggsvp3(int matrix_layout, char jobu, char jobv, char jobq, lapack_int m, lapack_int p, lapack_int n, T* a, lapack_int lda, T* b, lapack_int ldb, T tola, T tolb, lapack_int* k, lapack_int* l, T* u, lapack_int ldu, T* v, lapack_int ldv, T* q, lapack_int ldq) { return invokeFD(LAPACKE_sggsvp3, LAPACKE_dggsvp3, matrix_layout, jobu, jobv, jobq, m, p, n, a, lda, b, ldb, tola, tolb, k, l, u, ldu, v, ldv, q, ldq); }

/** Wraps LAPACKE_sggsvd3 and LAPACKE_dggsvd3. */
TEMPLATE_FLOAT_TYPE
inline lapack_int ggsvd3(int matrix_layout, char jobu, char jobv, char jobq, lapack_int m, lapack_int n, lapack_int p, lapack_int* k, lapack_int* l, T* a, lapack_int lda, T* b, lapack_int ldb, T* alpha, T* beta, T* u, lapack_int ldu, T* v, lapack_int ldv, T* q, lapack_int ldq, lapack_int* iwork) { return invokeFD(LAPACKE_sggsvd3, LAPACKE_dggsvd3, matrix_layout, jobu, jobv, jobq, m, n, p, k, l, a, lda, b, ldb, alpha, beta, u, ldu, v, ldv, q, ldq, iwork); }

/** Wraps LAPACKE_stgsja and LAPACKE_dtgsja. */
TEMPLATE_FLOAT_TYPE
inline lapack_int tgsja(int matrix_layout, char jobu, char jobv, char jobq, lapack_int m, lapack_int p, lapack_int n, lapack_int k, lapack_int l, T* a, lapack_int lda, T* b, lapack_int ldb, T tola, T tolb, T* alpha, T* beta, T* u, lapack_int ldu, T* v, lapack_int ldv, T* q, lapack_int ldq, lapack_int* ncycle) { return invokeFD(LAPACKE_stgsja, LAPACKE_dtgsja, matrix_layout, jobu, jobv, jobq, m, p, n, k, l, a, lda, b, ldb, tola, tolb, alpha, beta, u, ldu, v, ldv, q, ldq, ncycle); }

/** Wraps LAPACKE_slacpy and LAPACKE_dlacpy. */
TEMPLATE_FLOAT_TYPE
inline lapack_int lacpy(int matrix_layout, char uplo, lapack_int m, lapack_int n, const T* a, lapack_int lda, T* b, lapack_int ldb) { return invokeFD(LAPACKE_slacpy, LAPACKE_dlacpy, matrix_layout, uplo, m, n, a, lda, b, ldb); }

// Vectorized Ops

/** Wraps vsAdd and vdAdd. */
TEMPLATE_FLOAT_TYPE
inline void vAdd(const MKL_INT n, const T* a, const T* b, T* y) { invokeFD(vsAdd, vdAdd, n, a, b, y); }

/** Wraps vsSub and vdSub. */
TEMPLATE_FLOAT_TYPE
inline void vSub(const MKL_INT n, const T* a, const T* b, T* y) { invokeFD(vsSub, vdSub, n, a, b, y); }

/** Wraps vsMul and vdMul. */
TEMPLATE_FLOAT_TYPE
inline void vMul(const MKL_INT n, const T* a, const T* b, T* y) { invokeFD(vsMul, vdMul, n, a, b, y); }

/** Wraps vsSqr and vdSqr. */
TEMPLATE_FLOAT_TYPE
inline void vSqr(const MKL_INT n, const T* a, T* y) { invokeFD(vsSqr, vdSqr, n, a, y); }

/** Wraps vsAbs and vdAbs. */
TEMPLATE_FLOAT_TYPE
inline void vAbs(const MKL_INT n, const T* a, T* y) { invokeFD(vsAbs, vdAbs, n, a, y); }

/** Wraps vsLinearFrac and vdLinearFrac. */
TEMPLATE_FLOAT_TYPE
inline void vLinearFrac(const MKL_INT n, const T* a, const T* b, const T scalea, const T shifta, const T scaleb, const T shiftb, const T* y) { invokeFD(vsLinearFrac, vdLinearFrac, n, a, b, scalea, shifta, scaleb, shiftb, y); }

/** Wraps vsInv and vdInv. */
TEMPLATE_FLOAT_TYPE
inline void vInv(const MKL_INT n, const T* a, T* y) { invokeFD(vsInv, vdInv, n, a, y); }

/** Wraps vsDiv and vdDiv. */
TEMPLATE_FLOAT_TYPE
inline void vDiv(const MKL_INT n, const T* a, const T* b, T* y) { invokeFD(vsDiv, vdDiv, n, a, b, y); }

/** Wraps vsSqrt and vdSqrt. */
TEMPLATE_FLOAT_TYPE
inline void vSqrt(const MKL_INT n, const T* a, T* y) { invokeFD(vsSqrt, vdSqrt, n, a, y); }

/** Wraps vsInvSqrt and vdInvSqrt. */
TEMPLATE_FLOAT_TYPE
inline void vInvSqrt(const MKL_INT n, const T* a, T* y) { invokeFD(vsInvSqrt, vdInvSqrt, n, a, y); }

/** Wraps vsCbrt and vdCbrt. */
TEMPLATE_FLOAT_TYPE
inline void vCbrt(const MKL_INT n, const T* a, T* y) { invokeFD(vsCbrt, vdCbrt, n, a, y); }

/** Wraps vsInvCbrt and vdInvCbrt. */
TEMPLATE_FLOAT_TYPE
inline void vInvCbrt(const MKL_INT n, const T* a, T* y) { invokeFD(vsInvCbrt, vdInvCbrt, n, a, y); }

/** Wraps vsPow and vdPow. */
TEMPLATE_FLOAT_TYPE
inline void vPow(const MKL_INT n, const T* a, const T* b, T* y) { invokeFD(vsPow, vdPow, n, a, b, y); }

/** Wraps vsExp and vdExp. */
TEMPLATE_FLOAT_TYPE
inline void vExp(const MKL_INT n, const T* a, T* y) { invokeFD(vsExp, vdExp, n, a, y); }

/** Wraps vsExp2 and vdExp2. */
TEMPLATE_FLOAT_TYPE
inline void vExp2(const MKL_INT n, const T* a, T* y) { invokeFD(vsExp2, vdExp2, n, a, y); }

/** Wraps vsExp10 and vdExp10. */
TEMPLATE_FLOAT_TYPE
inline void vExp10(const MKL_INT n, const T* a, T* y) { invokeFD(vsExp10, vdExp10, n, a, y); }

/** Wraps vsLn and vdLn. */
TEMPLATE_FLOAT_TYPE
inline void vLn(const MKL_INT n, const T* a, T* y) { invokeFD(vsLn, vdLn, n, a, y); }

/** Wraps vsLog2 and vdLog2. */
TEMPLATE_FLOAT_TYPE
inline void vLog2(const MKL_INT n, const T* a, T* y) { invokeFD(vsLog2, vdLog2, n, a, y); }

/** Wraps vsLog10 and vdLog10. */
TEMPLATE_FLOAT_TYPE
inline void vLog10(const MKL_INT n, const T* a, T* y) { invokeFD(vsLog10, vdLog10, n, a, y); }

/** Wraps vsLogb and vdLogb. */
TEMPLATE_FLOAT_TYPE
inline void vLogb(const MKL_INT n, const T* a, T* y) { invokeFD(vsLogb, vdLogb, n, a, y); }

/** Wraps vsCos and vdCos. */
TEMPLATE_FLOAT_TYPE
inline void vCos(const MKL_INT n, const T* a, T* y) { invokeFD(vsCos, vdCos, n, a, y); }

/** Wraps vsSin and vdSin. */
TEMPLATE_FLOAT_TYPE
inline void vSin(const MKL_INT n, const T* a, T* y) { invokeFD(vsSin, vdSin, n, a, y); }

/** Wraps vsTan and vdTan. */
TEMPLATE_FLOAT_TYPE
inline void vTan(const MKL_INT n, const T* a, T* y) { invokeFD(vsTan, vdTan, n, a, y); }

/** Wraps vsAcos and vdAcos. */
TEMPLATE_FLOAT_TYPE
inline void vAcos(const MKL_INT n, const T* a, T* y) { invokeFD(vsAcos, vdAcos, n, a, y); }

/** Wraps vsAsin and vdAsin. */
TEMPLATE_FLOAT_TYPE
inline void vAsin(const MKL_INT n, const T* a, T* y) { invokeFD(vsAsin, vdAsin, n, a, y); }

/** Wraps vsAtan and vdAtan. */
TEMPLATE_FLOAT_TYPE
inline void vAtan(const MKL_INT n, const T* a, T* y) { invokeFD(vsAtan, vdAtan, n, a, y); }

/** Wraps vsAtan2 and vdAtan2. */
TEMPLATE_FLOAT_TYPE
inline void vAtan2(const MKL_INT n, const T* a, const T* b, T* y) { invokeFD(vsAtan2, vdAtan2, n, a, b, y); }

/** Wraps vsCosh and vdCosh. */
TEMPLATE_FLOAT_TYPE
inline void vCosh(const MKL_INT n, const T* a, T* y) { invokeFD(vsCosh, vdCosh, n, a, y); }

/** Wraps vsSinh and vdSinh. */
TEMPLATE_FLOAT_TYPE
inline void vSinh(const MKL_INT n, const T* a, T* y) { invokeFD(vsSinh, vdSinh, n, a, y); }

/** Wraps vsTanh and vdTanh. */
TEMPLATE_FLOAT_TYPE
inline void vTanh(const MKL_INT n, const T* a, T* y) { invokeFD(vsTanh, vdTanh, n, a, y); }

/** Wraps vsAcosh and vdAcosh. */
TEMPLATE_FLOAT_TYPE
inline void vAcosh(const MKL_INT n, const T* a, T* y) { invokeFD(vsAcosh, vdAcosh, n, a, y); }

/** Wraps vsAsinh and vdAsinh. */
TEMPLATE_FLOAT_TYPE
inline void vAsinh(const MKL_INT n, const T* a, T* y) { invokeFD(vsAsinh, vdAsinh, n, a, y); }

/** Wraps vsAtanh and vdAtanh. */
TEMPLATE_FLOAT_TYPE
inline void vAtanh(const MKL_INT n, const T* a, T* y) { invokeFD(vsAtanh, vdAtanh, n, a, y); }

/** Wraps vsFloor and vdFloor. */
TEMPLATE_FLOAT_TYPE
inline void vFloor(const MKL_INT n, const T* a, T* y) { invokeFD(vsFloor, vdFloor, n, a, y); }

/** Wraps vsCeil and vdCeil. */
TEMPLATE_FLOAT_TYPE
inline void vCeil(const MKL_INT n, const T* a, T* y) { invokeFD(vsCeil, vdCeil, n, a, y); }

/** Wraps vsTrunc and vdTrunc. */
TEMPLATE_FLOAT_TYPE
inline void vTrunc(const MKL_INT n, const T* a, T* y) { invokeFD(vsTrunc, vdTrunc, n, a, y); }

/** Wraps vsRound and vdRound. */
TEMPLATE_FLOAT_TYPE
inline void vRound(const MKL_INT n, const T* a, T* y) { invokeFD(vsRound, vdRound, n, a, y); }

/** Wraps vsNearbyInt and vdNearbyInt. */
TEMPLATE_FLOAT_TYPE
inline void vNearbyInt(const MKL_INT n, const T* a, T* y) { invokeFD(vsNearbyInt, vdNearbyInt, n, a, y); }

/** Wraps vsModf and vdModf. */
TEMPLATE_FLOAT_TYPE
inline void vModf(const MKL_INT n, const T* a, T* y, T* z) { invokeFD(vsModf, vdModf, n, a, y, z); }

/** Wraps vsFmax and vdFmax. */
TEMPLATE_FLOAT_TYPE
inline void vFmax(const MKL_INT n, const T* a, const T* b, T* y) { invokeFD(vsFmax, vdFmax, n, a, b, y); }

/** Wraps vsFmin and vdFmin. */
TEMPLATE_FLOAT_TYPE
inline void vFmin(const MKL_INT n, const T* a, const T* b, T* y) { invokeFD(vsFmin, vdFmin, n, a, b, y); }

/** Wraps vsMaxMag and vdMaxMag. */
TEMPLATE_FLOAT_TYPE
inline void vMaxMag(const MKL_INT n, const T* a, const T* b, T* y) { invokeFD(vsMaxMag, vdMaxMag, n, a, b, y); }

/** Wraps vsMinMag and vdMinMag. */
TEMPLATE_FLOAT_TYPE
inline void vMinMag(const MKL_INT n, const T* a, const T* b, T* y) { invokeFD(vsMinMag, vdMinMag, n, a, b, y); }

// Data Fitting

/**
 * Checks if an MKL data fitting return is signling an error.
 * \throw std::domain_error if num is a known MKL error code or less than 0.
 */
void dfCheckError(const int num);

/** Wraps dfsNewTask1D and dfdNewTask1D. */
TEMPLATE_FLOAT_TYPE
inline int dfNewTask1D(DFTaskPtr* task, const MKL_INT nx, const T* x, const MKL_INT xhint, const MKL_INT ny, const T* y, const MKL_INT yhint) { return invokeFD(dfsNewTask1D, dfdNewTask1D, task, nx, x, xhint, ny, y, yhint); }

/** Wraps dfsEditPPSpline1D and dfdEditPPSpline1D. */
TEMPLATE_FLOAT_TYPE
inline int dfEditPPSpline1D(DFTaskPtr task, const MKL_INT s_order, const MKL_INT s_type, const MKL_INT bc_type, const T* bc, const MKL_INT ic_type, const T* ic, const T* scoeff, const MKL_INT scoeffhint) { return invokeFD(dfsEditPPSpline1D, dfdEditPPSpline1D, task, s_order, s_type, bc_type, bc, ic_type, ic, scoeff, scoeffhint); }

/** Wraps dfsConstruct1D and dfdConstruct1D. */
template <typename T>
inline int dfConstruct1D(DFTaskPtr task, const MKL_INT s_format, const MKL_INT method);

template <>
inline int dfConstruct1D<float>(DFTaskPtr task, const MKL_INT s_format, const MKL_INT method) {
    return dfsConstruct1D(task, s_format, method);
}

template <>
inline int dfConstruct1D<double>(DFTaskPtr task, const MKL_INT s_format, const MKL_INT method) {
    return dfdConstruct1D(task, s_format, method);
}

/** Wraps dfsInterpolate1D and dfdInterpolate1D. */
TEMPLATE_FLOAT_TYPE
inline int dfInterpolate1D(DFTaskPtr task, const MKL_INT type, const MKL_INT method, const MKL_INT nsite, const T* site, const MKL_INT sitehint, const MKL_INT ndorder, const MKL_INT* dorder, const T* datahint, T* r, const MKL_INT rhint, MKL_INT* cell) { return invokeFD(dfsInterpolate1D, dfdInterpolate1D, task, type, method, nsite, site, sitehint, ndorder, dorder, datahint, r, rhint, cell); }

/** Wraps dfsInterpolateEx1D and dfdInterpolateEx1D. */
TEMPLATE_FLOAT_TYPE
inline int dfInterpolateEx1D(DFTaskPtr task, const MKL_INT type, const MKL_INT method, const MKL_INT nsite, const T* site, const MKL_INT sitehint, const MKL_INT ndorder, const MKL_INT* dorder, const T* datahint, T* r, const MKL_INT rhint, MKL_INT* cell, const void* le_params, const void* re_params, const void* i_params, const void* search_params) {
    /* TODO: This needs modification to incorporate callbacks. */
    return invokeFD(dfsInterpolateEx1D, dfdInterpolateEx1D, task, type, method, nsite, site, sitehint, ndorder, dorder, datahint, r, rhint, cell, NULL, le_params, NULL, re_params, NULL, i_params, NULL, search_params);
}

// }  // namespace MKL
}  // namespace Math
}  // namespace EVAA
