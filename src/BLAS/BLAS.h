// TODO: Copyright header

#pragma once
#include <mkl.h>

template <typename T>
class mkl {
public:
    static void mkl_error()
    {
        throw "Incorrect Datatype selected, supported types are Doubles and Floats";
    }

    // BLAS 1

    // DAXPY and SAXPY
    static void axpy(const MKL_INT n, const T a, const T* x, const MKL_INT incx, T* y,
                     const MKL_INT incy)
    {
#ifdef SINGLE_PRECISION
        cblas_saxpy(n, a, x, incx, y, incy);
#elif DOUBLE_PRECISION
        cblas_daxpy(n, a, x, incx, y, incy);
#else
        mkl::mkl_error();
#endif
    }

    // DCOPY and SCOPY
    static void copy(const MKL_INT n, const T* x, const MKL_INT incx, T* y, const MKL_INT incy)
    {
#ifdef SINGLE_PRECISION
        cblas_scopy(n, x, incx, y, incy);
#elif DOUBLE_PRECISION
        cblas_dcopy(n, x, incx, y, incy);
#else
        mkl::mkl_error();
#endif
    }

    // DDOT and SDOT
    static T dot(const MKL_INT n, const T* x, const MKL_INT incx, const T* y, const MKL_INT incy)
    {
#ifdef SINGLE_PRECISION
        return cblas_sdot(n, x, incx, y, incy);
#elif DOUBLE_PRECISION
        return cblas_ddot(n, x, incx, y, incy);
#else
        mkl::mkl_error();
#endif
    }

    // DNRM2 and SNRM2
    static T nrm2(const MKL_INT n, const T* x, const MKL_INT incx)
    {
#ifdef SINGLE_PRECISION
        return cblas_snrm2(n, x, incx);
#elif DOUBLE_PRECISION
        return cblas_dnrm2(n, x, incx);
#else
        mkl::mkl_error();
#endif
    }

    // DSCAL and SSCAL
    static void scal(const MKL_INT n, const T a, T* x, const MKL_INT incx)
    {
#ifdef SINGLE_PRECISION
        cblas_sscal(n, a, x, incx);
#elif DOUBLE_PRECISION
        cblas_dscal(n, a, x, incx);
#else
        mkl::mkl_error();
#endif
    }

    // DSWAP and SSWAP
    static void swap(const MKL_INT n, T* x, const MKL_INT incx, T* y, const MKL_INT incy)
    {
#ifdef SINGLE_PRECISION
        cblas_sswap(n, x, incx, y, incy);
#elif DOUBLE_PRECISION
        cblas_dswap(n, x, incx, y, incy);
#else
        mkl::mkl_error();
#endif
    }

    // DROT and SROT
    static void rot(const MKL_INT n, T* x, const MKL_INT incx, T* y, const MKL_INT incy, const T c,
                    const T s)
    {
#ifdef SINGLE_PRECISION
        cblas_srot(n, x, incx, y, incy, c, s);
#elif DOUBLE_PRECISION
        cblas_drot(n, x, incx, y, incy, c, s);
#else
        mkl::mkl_error();
#endif
    }

    // BLAS 2

    static void gbmv(const CBLAS_LAYOUT Layout, const CBLAS_TRANSPOSE trans, const MKL_INT m,
                     const MKL_INT n, const MKL_INT kl, const MKL_INT ku, const T alpha, const T* a,
                     const MKL_INT lda, const T* x, const MKL_INT incx, const T beta, T* y,
                     const MKL_INT incy)
    {
#ifdef SINGLE_PRECISION
        cblas_sgbmv(Layout, trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy);
#elif DOUBLE_PRECISION
        cblas_dgbmv(Layout, trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy);
#else
        mkl::mkl_error();
#endif
    }

    static void gemv(const CBLAS_LAYOUT Layout, const CBLAS_TRANSPOSE trans, const MKL_INT m,
                     const MKL_INT n, const T alpha, const T* a, const MKL_INT lda, const T* x,
                     const MKL_INT incx, const T beta, T* y, const MKL_INT incy)
    {
#ifdef SINGLE_PRECISION
        cblas_sgemv(Layout, trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
#elif DOUBLE_PRECISION
        cblas_dgemv(Layout, trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
#else
        mkl::mkl_error();
#endif
    }

    static void ger(const CBLAS_LAYOUT Layout, const MKL_INT m, const MKL_INT n, const T alpha,
                    const T* x, const MKL_INT incx, const T* y, const MKL_INT incy, T* a,
                    const MKL_INT lda)
    {
#ifdef SINGLE_PRECISION
        cblas_sger(Layout, m, n, alpha, x, incx, y, incy, a, lda);
#elif DOUBLE_PRECISION
        cblas_dger(Layout, m, n, alpha, x, incx, y, incy, a, lda);
#else
        mkl::mkl_error();
#endif
    }

    static void sbmv(const CBLAS_LAYOUT Layout, const CBLAS_UPLO uplo, const MKL_INT n,
                     const MKL_INT k, const T alpha, const T* a, const MKL_INT lda, const T* x,
                     const MKL_INT incx, const T beta, T* y, const MKL_INT incy)
    {
#ifdef SINGLE_PRECISION
        cblas_ssbmv(Layout, uplo, n, k, alpha, a, lda, x, incx, beta, y, incy);
#elif DOUBLE_PRECISION
        cblas_dsbmv(Layout, uplo, n, k, alpha, a, lda, x, incx, beta, y, incy);
#else
        mkl::mkl_error();
#endif
    }

    static void imatcopy(const char ordering, const char trans, size_t rows, size_t cols,
                         const T alpha, T* AB, size_t lda, size_t ldb)
    {
#ifdef SINGLE_PRECISION
        mkl_simatcopy(ordering, trans, rows, cols, alpha, AB, lda, ldb);
#elif DOUBLE_PRECISION
        mkl_dimatcopy(ordering, trans, rows, cols, alpha, AB, lda, ldb);
#else
        mkl::mkl_error();
#endif
    }

    static void symv(const CBLAS_LAYOUT Layout, const CBLAS_UPLO uplo, const MKL_INT n,
                     const T alpha, const T* a, const MKL_INT lda, const T* x, const MKL_INT incx,
                     const T beta, T* y, const MKL_INT incy)
    {
#ifdef SINGLE_PRECISION
        cblas_ssymv(Layout, uplo, n, alpha, a, lda, x, incx, beta, y, incy);
#elif DOUBLE_PRECISION
        cblas_dsymv(Layout, uplo, n, alpha, a, lda, x, incx, beta, y, incy);
#else
        mkl::mkl_error();
#endif
    }

    static void syr(const CBLAS_LAYOUT Layout, const CBLAS_UPLO uplo, const MKL_INT n,
                    const T alpha, const T* x, const MKL_INT incx, T* a, const MKL_INT lda)
    {
#ifdef SINGLE_PRECISION
        cblas_ssyr(Layout, uplo, n, alpha, x, incx, a, lda);
#elif DOUBLE_PRECISION
        cblas_dsyr(Layout, uplo, n, alpha, x, incx, a, lda);
#else
        mkl::mkl_error();
#endif
    }

    static void syr2(const CBLAS_LAYOUT Layout, const CBLAS_UPLO uplo, const MKL_INT n,
                     const T alpha, const T* x, const MKL_INT incx, const T* y, const MKL_INT incy,
                     T* a, const MKL_INT lda)
    {
#ifdef SINGLE_PRECISION
        cblas_ssyr2(Layout, uplo, n, alpha, x, incx, y, incy, a, lda);
#elif DOUBLE_PRECISION
        cblas_dsyr2(Layout, uplo, n, alpha, x, incx, y, incy, a, lda);
#else
        mkl::mkl_error();
#endif
    }

    static void tbmv(const CBLAS_LAYOUT Layout, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE trans,
                     const CBLAS_DIAG diag, const MKL_INT n, const MKL_INT k, const T* a,
                     const MKL_INT lda, T* x, const MKL_INT incx)
    {
#ifdef SINGLE_PRECISION
        cblas_stbmv(Layout, uplo, trans, diag, n, k, a, lda, x, incx);
#elif DOUBLE_PRECISION
        cblas_dtbmv(Layout, uplo, trans, diag, n, k, a, lda, x, incx);
#else
        mkl::mkl_error();
#endif
    }

    // BLAS 3

    static void gemm(const CBLAS_LAYOUT Layout, const CBLAS_TRANSPOSE transa,
                     const CBLAS_TRANSPOSE transb, const MKL_INT m, const MKL_INT n,
                     const MKL_INT k, const T alpha, const T* a, const MKL_INT lda, const T* b,
                     const MKL_INT ldb, const T beta, T* c, const MKL_INT ldc)
    {
#ifdef SINGLE_PRECISION
        cblas_sgemm(Layout, transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
#elif DOUBLE_PRECISION
        cblas_dgemm(Layout, transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
#else
        mkl::mkl_error();
#endif
    }

    static void symm(const CBLAS_LAYOUT Layout, const CBLAS_SIDE side, const CBLAS_UPLO uplo,
                     const MKL_INT m, const MKL_INT n, const T alpha, const T* a, const MKL_INT lda,
                     const T* b, const MKL_INT ldb, const T beta, T* c, const MKL_INT ldc)
    {
#ifdef SINGLE_PRECISION
        cblas_ssymm(Layout, side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc);
#elif DOUBLE_PRECISION
        cblas_dsymm(Layout, side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc);
#else
        mkl::mkl_error();
#endif
    }

    static void syrk(const CBLAS_LAYOUT Layout, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE trans,
                     const MKL_INT n, const MKL_INT k, const T alpha, const T* a, const MKL_INT lda,
                     const T beta, T* c, const MKL_INT ldc)
    {
#ifdef SINGLE_PRECISION
        cblas_ssyrk(Layout, uplo, trans, n, k, alpha, a, lda, beta, c, ldc);
#elif DOUBLE_PRECISION
        cblas_dsyrk(Layout, uplo, trans, n, k, alpha, a, lda, beta, c, ldc);
#else
        mkl::mkl_error();
#endif
    }

    static void trmm(const CBLAS_LAYOUT Layout, const CBLAS_SIDE side, const CBLAS_UPLO uplo,
                     const CBLAS_TRANSPOSE transa, const CBLAS_DIAG diag, const MKL_INT m,
                     const MKL_INT n, const T alpha, const T* a, const MKL_INT lda, T* b,
                     const MKL_INT ldb)
    {
#ifdef SINGLE_PRECISION
        cblas_strmm(Layout, side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
#elif DOUBLE_PRECISION
        cblas_dtrmm(Layout, side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
#else
        mkl::mkl_error();
#endif
    }

    // LAPACK

    static lapack_int getrf(int matrix_layout, lapack_int m, lapack_int n, T* a, lapack_int lda,
                            lapack_int* ipiv)
    {
#ifdef SINGLE_PRECISION
        return LAPACKE_sgetrf(matrix_layout, m, n, a, lda, ipiv);
#elif DOUBLE_PRECISION
        return LAPACKE_dgetrf(matrix_layout, m, n, a, lda, ipiv);
#else
        mkl::mkl_error();
#endif
    }

    static lapack_int getrs(int matrix_layout, char trans, lapack_int n, lapack_int nrhs,
                            const T* a, lapack_int lda, const lapack_int* ipiv, T* b,
                            lapack_int ldb)
    {
#ifdef SINGLE_PRECISION
        return LAPACKE_sgetrs(matrix_layout, trans, n, nrhs, a, lda, ipiv, b, ldb);
#elif DOUBLE_PRECISION
        return LAPACKE_dgetrs(matrix_layout, trans, n, nrhs, a, lda, ipiv, b, ldb);
#else
        mkl::mkl_error();
#endif
    }

    static lapack_int gecon(int matrix_layout, char norm, lapack_int n, const T* a, lapack_int lda,
                            T anorm, T* rcond)
    {
#ifdef SINGLE_PRECISION
        return LAPACKE_sgecon(matrix_layout, norm, n, a, lda, anorm, rcond);
#elif DOUBLE_PRECISION
        return LAPACKE_dgecon(matrix_layout, norm, n, a, lda, anorm, rcond);
#else
        mkl::mkl_error();
#endif
    }

    static lapack_int gbtrf(int matrix_layout, lapack_int m, lapack_int n, lapack_int kl,
                            lapack_int ku, T* ab, lapack_int ldab, lapack_int* ipiv)
    {
#ifdef SINGLE_PRECISION
        return LAPACKE_sgbtrf(matrix_layout, m, n, kl, ku, ab, ldab, ipiv);
#elif DOUBLE_PRECISION
        return LAPACKE_dgbtrf(matrix_layout, m, n, kl, ku, ab, ldab, ipiv);
#else
        mkl::mkl_error();
#endif
    }

    static lapack_int gbtrs(int matrix_layout, char trans, lapack_int n, lapack_int kl,
                            lapack_int ku, lapack_int nrhs, const T* ab, lapack_int ldab,
                            const lapack_int* ipiv, T* b, lapack_int ldb)
    {
#ifdef SINGLE_PRECISION
        return LAPACKE_sgbtrs(matrix_layout, trans, n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb);
#elif DOUBLE_PRECISION
        return LAPACKE_dgbtrs(matrix_layout, trans, n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb);
#else
        mkl::mkl_error();
#endif
    }

    static lapack_int gbcon(int matrix_layout, char norm, lapack_int n, lapack_int kl,
                            lapack_int ku, const T* ab, lapack_int ldab, const lapack_int* ipiv,
                            T anorm, T* rcond)
    {
#ifdef SINGLE_PRECISION
        return LAPACKE_sgbcon(matrix_layout, norm, n, kl, ku, ab, ldab, ipiv, anorm, rcond);
#elif DOUBLE_PRECISION
        return LAPACKE_dgbcon(matrix_layout, norm, n, kl, ku, ab, ldab, ipiv, anorm, rcond);
#else
        mkl::mkl_error();
#endif
    }

    static lapack_int potrf(int matrix_layout, char uplo, lapack_int n, T* a, lapack_int lda)
    {
#ifdef SINGLE_PRECISION
        return LAPACKE_spotrf(matrix_layout, uplo, n, a, lda);
#elif DOUBLE_PRECISION
        return LAPACKE_dpotrf(matrix_layout, uplo, n, a, lda);
#else
        mkl::mkl_error();
#endif
    }

    static lapack_int potrs(int matrix_layout, char uplo, lapack_int n, lapack_int nrhs, const T* a,
                            lapack_int lda, T* b, lapack_int ldb)
    {
#ifdef SINGLE_PRECISION
        return LAPACKE_spotrs(matrix_layout, uplo, n, nrhs, a, lda, b, ldb);
#elif DOUBLE_PRECISION
        return LAPACKE_dpotrs(matrix_layout, uplo, n, nrhs, a, lda, b, ldb);
#else
        mkl::mkl_error();
#endif
    }

    static lapack_int pocon(int matrix_layout, char uplo, lapack_int n, const T* a, lapack_int lda,
                            T anorm, T* rcond)
    {
#ifdef SINGLE_PRECISION
        return LAPACKE_spocon(matrix_layout, uplo, n, a, lda, anorm, rcond);
#elif DOUBLE_PRECISION
        return LAPACKE_dpocon(matrix_layout, uplo, n, a, lda, anorm, rcond);
#else
        mkl::mkl_error();
#endif
    }

    static lapack_int pttrf(lapack_int n, T* d, T* e)
    {
#ifdef SINGLE_PRECISION
        return LAPACKE_spttrf(n, d, e);
#elif DOUBLE_PRECISION
        return LAPACKE_dpttrf(n, d, e);
#else
        mkl::mkl_error();
#endif
    }

    static lapack_int pttrs(int matrix_layout, lapack_int n, lapack_int nrhs, const T* d,
                            const T* e, T* b, lapack_int ldb)
    {
#ifdef SINGLE_PRECISION
        return LAPACKE_spttrs(matrix_layout, n, nrhs, d, e, b, ldb);
#elif DOUBLE_PRECISION
        return LAPACKE_dpttrs(matrix_layout, n, nrhs, d, e, b, ldb);
#else
        mkl::mkl_error();
#endif
    }

    static lapack_int ptcon(lapack_int n, const T* d, const T* e, T anorm, T* rcond)
    {
#ifdef SINGLE_PRECISION
        return LAPACKE_sptcon(n, d, e, anorm, rcond);
#elif DOUBLE_PRECISION
        return LAPACKE_dptcon(n, d, e, anorm, rcond);
#else
        mkl::mkl_error();
#endif
    }

    // Orthogonal Factorizations

    static lapack_int geqrf(int matrix_layout, lapack_int m, lapack_int n, T* a, lapack_int lda,
                            T* tau)
    {
#ifdef SINGLE_PRECISION
        return LAPACKE_sgeqrf(matrix_layout, m, n, a, lda, tau);
#elif DOUBLE_PRECISION
        return LAPACKE_dgeqrf(matrix_layout, m, n, a, lda, tau);
#else
        mkl::mkl_error();
#endif
    }

    static lapack_int geqpf(int matrix_layout, lapack_int m, lapack_int n, T* a, lapack_int lda,
                            lapack_int* jpvt, T* tau)
    {
#ifdef SINGLE_PRECISION
        return LAPACKE_sgeqpf(matrix_layout, m, n, a, lda, jpvt, tau);
#elif DOUBLE_PRECISION
        return LAPACKE_dgeqpf(matrix_layout, m, n, a, lda, jpvt, tau);
#else
        mkl::mkl_error();
#endif
    }

    static lapack_int orgqr(int matrix_layout, lapack_int m, lapack_int n, lapack_int k, T* a,
                            lapack_int lda, const T* tau)
    {
#ifdef SINGLE_PRECISION
        return LAPACKE_sorgqr(matrix_layout, m, n, k, a, lda, tau);
#elif DOUBLE_PRECISION
        return LAPACKE_dorgqr(matrix_layout, m, n, k, a, lda, tau);
#else
        mkl::mkl_error();
#endif
    }

    static lapack_int ormqr(int matrix_layout, char side, char trans, lapack_int m, lapack_int n,
                            lapack_int k, const T* a, lapack_int lda, const T* tau, T* c,
                            lapack_int ldc)
    {
#ifdef SINGLE_PRECISION
        return LAPACKE_sormqr(matrix_layout, side, trans, m, n, k, a, lda, tau, c, ldc);
#elif DOUBLE_PRECISION
        return LAPACKE_dormqr(matrix_layout, side, trans, m, n, k, a, lda, tau, c, ldc);
#else
        mkl::mkl_error();
#endif
    }

    static lapack_int gebrd(int matrix_layout, lapack_int m, lapack_int n, T* a, lapack_int lda,
                            T* d, T* e, T* tauq, T* taup)
    {
#ifdef SINGLE_PRECISION
        return LAPACKE_sgebrd(matrix_layout, m, n, a, lda, d, e, tauq, taup);
#elif DOUBLE_PRECISION
        return LAPACKE_dgebrd(matrix_layout, m, n, a, lda, d, e, tauq, taup);
#else
        mkl::mkl_error();
#endif
    }

    static lapack_int bdsqr(int matrix_layout, char uplo, lapack_int n, lapack_int ncvt,
                            lapack_int nru, lapack_int ncc, T* d, T* e, T* vt, lapack_int ldvt,
                            T* u, lapack_int ldu, T* c, lapack_int ldc)
    {
#ifdef SINGLE_PRECISION
        return LAPACKE_sbdsqr(matrix_layout, uplo, n, ncvt, nru, ncc, d, e, vt, ldvt, u, ldu, c,
                              ldc);
#elif DOUBLE_PRECISION
        return LAPACKE_dbdsqr(matrix_layout, uplo, n, ncvt, nru, ncc, d, e, vt, ldvt, u, ldu, c,
                              ldc);
#else
        mkl::mkl_error();
#endif
    }

    static lapack_int bdsdc(int matrix_layout, char uplo, char compq, lapack_int n, T* d, T* e,
                            T* u, lapack_int ldu, T* vt, lapack_int ldvt, T* q, lapack_int* iq)
    {
#ifdef SINGLE_PRECISION
        return LAPACKE_sbdsdc(matrix_layout, uplo, compq, n, d, e, u, ldu, vt, ldvt, q, iq);
#elif DOUBLE_PRECISION
        return LAPACKE_dbdsdc(matrix_layout, uplo, compq, n, d, e, u, ldu, vt, ldvt, q, iq);
#else
        mkl::mkl_error();
#endif
    }

    static lapack_int sytrd(int matrix_layout, char uplo, lapack_int n, T* a, lapack_int lda, T* d,
                            T* e, T* tau)
    {
#ifdef SINGLE_PRECISION
        return LAPACKE_ssytrd(matrix_layout, uplo, n, a, lda, d, e, tau);
        ;
#elif DOUBLE_PRECISION
        return LAPACKE_dsytrd(matrix_layout, uplo, n, a, lda, d, e, tau);
#else
        mkl::mkl_error();
#endif
    }

    static lapack_int sbtrd(int matrix_layout, char vect, char uplo, lapack_int n, lapack_int kd,
                            T* ab, lapack_int ldab, T* d, T* e, T* q, lapack_int ldq)
    {
#ifdef SINGLE_PRECISION
        return LAPACKE_ssbtrd(matrix_layout, vect, uplo, n, kd, ab, ldab, d, e, q, ldq);
#elif DOUBLE_PRECISION
        return LAPACKE_dsbtrd(matrix_layout, vect, uplo, n, kd, ab, ldab, d, e, q, ldq);
#else
        mkl::mkl_error();
#endif
    }

    static lapack_int orgtr(int matrix_layout, char uplo, lapack_int n, T* a, lapack_int lda,
                            const T* tau)
    {
#ifdef SINGLE_PRECISION
        return LAPACKE_sorgtr(matrix_layout, uplo, n, a, lda, tau);
#elif DOUBLE_PRECISION
        return LAPACKE_dorgtr(matrix_layout, uplo, n, a, lda, tau);
#else
        mkl::mkl_error();
#endif
    }

    static lapack_int ormtr(int matrix_layout, char side, char uplo, char trans, lapack_int m,
                            lapack_int n, const T* a, lapack_int lda, const T* tau, T* c,
                            lapack_int ldc)
    {
#ifdef SINGLE_PRECISION
        return LAPACKE_sormtr(matrix_layout, side, uplo, trans, m, n, a, lda, tau, c, ldc);
#elif DOUBLE_PRECISION
        return LAPACKE_dormtr(matrix_layout, side, uplo, trans, m, n, a, lda, tau, c, ldc);
#else
        mkl::mkl_error();
#endif
    }

    static lapack_int steqr(int matrix_layout, char compz, lapack_int n, T* d, T* e, T* z,
                            lapack_int ldz)
    {
#ifdef SINGLE_PRECISION
        return LAPACKE_ssteqr(matrix_layout, compz, n, d, e, z, ldz);
#elif DOUBLE_PRECISION
        return LAPACKE_dsteqr(matrix_layout, compz, n, d, e, z, ldz);
#else
        mkl::mkl_error();
#endif
    }

    static lapack_int stedc(int matrix_layout, char compz, lapack_int n, T* d, T* e, T* z,
                            lapack_int ldz)
    {
#ifdef SINGLE_PRECISION
        return LAPACKE_sstedc(matrix_layout, compz, n, d, e, z, ldz);
#elif DOUBLE_PRECISION
        return LAPACKE_dstedc(matrix_layout, compz, n, d, e, z, ldz);
#else
        mkl::mkl_error();
#endif
    }

    static lapack_int pteqr(int matrix_layout, char compz, lapack_int n, T* d, T* e, T* z,
                            lapack_int ldz)
    {
#ifdef SINGLE_PRECISION
        return LAPACKE_spteqr(matrix_layout, compz, n, d, e, z, ldz);
#elif DOUBLE_PRECISION
        return LAPACKE_dpteqr(matrix_layout, compz, n, d, e, z, ldz);
#else
        mkl::mkl_error();
#endif
    }

    static lapack_int ggsvp(int matrix_layout, char jobu, char jobv, char jobq, lapack_int m,
                            lapack_int p, lapack_int n, T* a, lapack_int lda, T* b, lapack_int ldb,
                            T tola, T tolb, lapack_int* k, lapack_int* l, T* u, lapack_int ldu,
                            T* v, lapack_int ldv, T* q, lapack_int ldq)
    {
#ifdef SINGLE_PRECISION
        return LAPACKE_sggsvp(matrix_layout, jobu, jobv, jobq, m, p, n, a, lda, b, ldb, tola, tolb,
                              k, l, u, ldu, v, ldv, q, ldq);
#elif DOUBLE_PRECISION
        return LAPACKE_dggsvp(matrix_layout, jobu, jobv, jobq, m, p, n, a, lda, b, ldb, tola, tolb,
                              k, l, u, ldu, v, ldv, q, ldq);
#else
        mkl::mkl_error();
#endif
    }

    static lapack_int ggsvp3(int matrix_layout, char jobu, char jobv, char jobq, lapack_int m,
                             lapack_int p, lapack_int n, T* a, lapack_int lda, T* b, lapack_int ldb,
                             T tola, T tolb, lapack_int* k, lapack_int* l, T* u, lapack_int ldu,
                             T* v, lapack_int ldv, T* q, lapack_int ldq)
    {
#ifdef SINGLE_PRECISION
        return LAPACKE_sggsvp3(matrix_layout, jobu, jobv, jobq, m, p, n, a, lda, b, ldb, tola, tolb,
                               k, l, u, ldu, v, ldv, q, ldq);
#elif DOUBLE_PRECISION
        return LAPACKE_dggsvp3(matrix_layout, jobu, jobv, jobq, m, p, n, a, lda, b, ldb, tola, tolb,
                               k, l, u, ldu, v, ldv, q, ldq);
#else
        mkl::mkl_error();
#endif
    }

    static lapack_int ggsvd3(int matrix_layout, char jobu, char jobv, char jobq, lapack_int m,
                             lapack_int n, lapack_int p, lapack_int* k, lapack_int* l, T* a,
                             lapack_int lda, T* b, lapack_int ldb, T* alpha, T* beta, T* u,
                             lapack_int ldu, T* v, lapack_int ldv, T* q, lapack_int ldq,
                             lapack_int* iwork)
    {
#ifdef SINGLE_PRECISION
        return LAPACKE_sggsvd3(matrix_layout, jobu, jobv, jobq, m, n, p, k, l, a, lda, b, ldb,
                               alpha, beta, u, ldu, v, ldv, q, ldq, iwork);
#elif DOUBLE_PRECISION
        return LAPACKE_dggsvd3(matrix_layout, jobu, jobv, jobq, m, n, p, k, l, a, lda, b, ldb,
                               alpha, beta, u, ldu, v, ldv, q, ldq, iwork);
#else
        mkl::mkl_error();
#endif
    }

    static lapack_int tgsja(int matrix_layout, char jobu, char jobv, char jobq, lapack_int m,
                            lapack_int p, lapack_int n, lapack_int k, lapack_int l, T* a,
                            lapack_int lda, T* b, lapack_int ldb, T tola, T tolb, T* alpha, T* beta,
                            T* u, lapack_int ldu, T* v, lapack_int ldv, T* q, lapack_int ldq,
                            lapack_int* ncycle)
    {
#ifdef SINGLE_PRECISION
        return LAPACKE_stgsja(matrix_layout, jobu, jobv, jobq, m, p, n, k, l, a, lda, b, ldb, tola,
                              tolb, alpha, beta, u, ldu, v, ldv, q, ldq, ncycle);
#elif DOUBLE_PRECISION
        return LAPACKE_dtgsja(matrix_layout, jobu, jobv, jobq, m, p, n, k, l, a, lda, b, ldb, tola,
                              tolb, alpha, beta, u, ldu, v, ldv, q, ldq, ncycle);
#else
        mkl::mkl_error();
#endif
    }

    static lapack_int lacpy(int matrix_layout, char uplo, lapack_int m, lapack_int n, const T* a,
                            lapack_int lda, T* b, lapack_int ldb)
    {
#ifdef SINGLE_PRECISION
        return LAPACKE_slacpy(matrix_layout, uplo, m, n, a, lda, b, ldb);
#elif DOUBLE_PRECISION
        return LAPACKE_dlacpy(matrix_layout, uplo, m, n, a, lda, b, ldb);
#else
        mkl::mkl_error();
#endif
    }

    // Vectorized Ops

    static void vAdd(const MKL_INT n, const T* a, const T* b, T* y)
    {
#ifdef SINGLE_PRECISION
        vsAdd(n, a, b, y);
#elif DOUBLE_PRECISION
        vdAdd(n, a, b, y);
#else
        mkl::mkl_error();
#endif
    }

    static void vSub(const MKL_INT n, const T* a, const T* b, T* y)
    {
#ifdef SINGLE_PRECISION
        vsSub(n, a, b, y);
#elif DOUBLE_PRECISION
        vdSub(n, a, b, y);
#else
        mkl::mkl_error();
#endif
    }

    static void vMul(const MKL_INT n, const T* a, const T* b, T* y)
    {
#ifdef SINGLE_PRECISION
        vsMul(n, a, b, y);
#elif DOUBLE_PRECISION
        vdMul(n, a, b, y);
#else
        mkl::mkl_error();
#endif
    }

    static void vSqr(const MKL_INT n, const T* a, T* y)
    {
#ifdef SINGLE_PRECISION
        vsSqr(n, a, y);
#elif DOUBLE_PRECISION
        vdSqr(n, a, y);
#else
        mkl::mkl_error();
#endif
    }

    static void vAbs(const MKL_INT n, const T* a, T* y)
    {
#ifdef SINGLE_PRECISION
        vsAbs(n, a, y);
#elif DOUBLE_PRECISION
        vdAbs(n, a, y);
#else
        mkl::mkl_error();
#endif
    }

    static void vLinearFrac(const MKL_INT n, const T* a, const T* b, const T scalea, const T shifta,
                            const T scaleb, const T shiftb, const T* y)
    {
#ifdef SINGLE_PRECISION
        vsLinearFrac(n, a, b, scalea, shifta, scaleb, shiftb, y);
#elif DOUBLE_PRECISION
        vdLinearFrac(n, a, b, scalea, shifta, scaleb, shiftb, y);
#else
        mkl::mkl_error();
#endif
    }

    static void vInv(const MKL_INT n, const T* a, T* y)
    {
#ifdef SINGLE_PRECISION
        vsInv(n, a, y);
#elif DOUBLE_PRECISION
        vdInv(n, a, y);
#else
        mkl::mkl_error();
#endif
    }

    static void vDiv(const MKL_INT n, const T* a, const T* b, T* y)
    {
#ifdef SINGLE_PRECISION
        vsDiv(n, a, b, y);
#elif DOUBLE_PRECISION
        vdDiv(n, a, b, y);
#else
        mkl::mkl_error();
#endif
    }

    static void vSqrt(const MKL_INT n, const T* a, T* y)
    {
#ifdef SINGLE_PRECISION
        vsSqrt(n, a, y);
#elif DOUBLE_PRECISION
        vdSqrt(n, a, y);
#else
        mkl::mkl_error();
#endif
    }

    static void vInvSqrt(const MKL_INT n, const T* a, T* y)
    {
#ifdef SINGLE_PRECISION
        vsInvSqrt(n, a, y);
#elif DOUBLE_PRECISION
        vdInvSqrt(n, a, y);
#else
        mkl::mkl_error();
#endif
    }

    static void vCbrt(const MKL_INT n, const T* a, T* y)
    {
#ifdef SINGLE_PRECISION
        vsCbrt(n, a, y);
#elif DOUBLE_PRECISION
        vdCbrt(n, a, y);
#else
        mkl::mkl_error();
#endif
    }

    static void vInvCbrt(const MKL_INT n, const T* a, T* y)
    {
#ifdef SINGLE_PRECISION
        vsInvCbrt(n, a, y);
#elif DOUBLE_PRECISION
        vdInvCbrt(n, a, y);
#else
        mkl::mkl_error();
#endif
    }

    static void vPow(const MKL_INT n, const T* a, const T* b, T* y)
    {
#ifdef SINGLE_PRECISION
        vsPow(n, a, b, y);
#elif DOUBLE_PRECISION
        vdPow(n, a, b, y);
#else
        mkl::mkl_error();
#endif
    }

    static void vExp(const MKL_INT n, const T* a, T* y)
    {
#ifdef SINGLE_PRECISION
        vsExp(n, a, y);
#elif DOUBLE_PRECISION
        vdExp(n, a, y);
#else
        mkl::mkl_error();
#endif
    }

    static void vExp2(const MKL_INT n, const T* a, T* y)
    {
#ifdef SINGLE_PRECISION
        vsExp2(n, a, y);
#elif DOUBLE_PRECISION
        vdExp2(n, a, y);
#else
        mkl::mkl_error();
#endif
    }

    static void vExp10(const MKL_INT n, const T* a, T* y)
    {
#ifdef SINGLE_PRECISION
        vsExp10(n, a, y);
#elif DOUBLE_PRECISION
        vdExp10(n, a, y);
#else
        mkl::mkl_error();
#endif
    }

    static void vLn(const MKL_INT n, const T* a, T* y)
    {
#ifdef SINGLE_PRECISION
        vsLn(n, a, y);
#elif DOUBLE_PRECISION
        vdLn(n, a, y);
#else
        mkl::mkl_error();
#endif
    }

    static void vLog2(const MKL_INT n, const T* a, T* y)
    {
#ifdef SINGLE_PRECISION
        vsLog2(n, a, y);
#elif DOUBLE_PRECISION
        vdLog2(n, a, y);
#else
        mkl::mkl_error();
#endif
    }

    static void vLog10(const MKL_INT n, const T* a, T* y)
    {
#ifdef SINGLE_PRECISION
        vsLog10(n, a, y);
#elif DOUBLE_PRECISION
        vdLog10(n, a, y);
#else
        mkl::mkl_error();
#endif
    }

    static void vLogb(const MKL_INT n, const T* a, T* y)
    {
#ifdef SINGLE_PRECISION
        vsLogb(n, a, y);
#elif DOUBLE_PRECISION
        vdLogb(n, a, y);
#else
        mkl::mkl_error();
#endif
    }

    static void vCos(const MKL_INT n, const T* a, T* y)
    {
#ifdef SINGLE_PRECISION
        vsCos(n, a, y);
#elif DOUBLE_PRECISION
        vdCos(n, a, y);
#else
        mkl::mkl_error();
#endif
    }

    static void vSin(const MKL_INT n, const T* a, T* y)
    {
#ifdef SINGLE_PRECISION
        vsSin(n, a, y);
#elif DOUBLE_PRECISION
        vdSin(n, a, y);
#else
        mkl::mkl_error();
#endif
    }

    static void vTan(const MKL_INT n, const T* a, T* y)
    {
#ifdef SINGLE_PRECISION
        vsTan(n, a, y);
#elif DOUBLE_PRECISION
        vdTan(n, a, y);
#else
        mkl::mkl_error();
#endif
    }

    static void vAcos(const MKL_INT n, const T* a, T* y)
    {
#ifdef SINGLE_PRECISION
        vsAcos(n, a, y);
#elif DOUBLE_PRECISION
        vdAcos(n, a, y);
#else
        mkl::mkl_error();
#endif
    }

    static void vAsin(const MKL_INT n, const T* a, T* y)
    {
#ifdef SINGLE_PRECISION
        vsAsin(n, a, y);
#elif DOUBLE_PRECISION
        vdAsin(n, a, y);
#else
        mkl::mkl_error();
#endif
    }

    static void vAtan(const MKL_INT n, const T* a, T* y)
    {
#ifdef SINGLE_PRECISION
        vsAtan(n, a, y);
#elif DOUBLE_PRECISION
        vdAtan(n, a, y);
#else
        mkl::mkl_error();
#endif
    }

    static void vAtan2(const MKL_INT n, const T* a, const T* b, T* y)
    {
#ifdef SINGLE_PRECISION
        vsAtan2(n, a, b, y);
#elif DOUBLE_PRECISION
        vdAtan2(n, a, b, y);
#else
        mkl::mkl_error();
#endif
    }

    static void vCosh(const MKL_INT n, const T* a, T* y)
    {
#ifdef SINGLE_PRECISION
        vsCosh(n, a, y);
#elif DOUBLE_PRECISION
        vdCosh(n, a, y);
#else
        mkl::mkl_error();
#endif
    }

    static void vSinh(const MKL_INT n, const T* a, T* y)
    {
#ifdef SINGLE_PRECISION
        vsSinh(n, a, y);
#elif DOUBLE_PRECISION
        vdSinh(n, a, y);
#else
        mkl::mkl_error();
#endif
    }

    static void vTanh(const MKL_INT n, const T* a, T* y)
    {
#ifdef SINGLE_PRECISION
        vsTanh(n, a, y);
#elif DOUBLE_PRECISION
        vdTanh(n, a, y);
#else
        mkl::mkl_error();
#endif
    }

    static void vAcosh(const MKL_INT n, const T* a, T* y)
    {
#ifdef SINGLE_PRECISION
        vsAcosh(n, a, y);
#elif DOUBLE_PRECISION
        vdAcosh(n, a, y);
#else
        mkl::mkl_error();
#endif
    }

    static void vAsinh(const MKL_INT n, const T* a, T* y)
    {
#ifdef SINGLE_PRECISION
        vsAsinh(n, a, y);
#elif DOUBLE_PRECISION
        vdAsinh(n, a, y);
#else
        mkl::mkl_error();
#endif
    }

    static void vAtanh(const MKL_INT n, const T* a, T* y)
    {
#ifdef SINGLE_PRECISION
        vsAtanh(n, a, y);
#elif DOUBLE_PRECISION
        vdAtanh(n, a, y);
#else
        mkl::mkl_error();
#endif
    }

    static void vFloor(const MKL_INT n, const T* a, T* y)
    {
#ifdef SINGLE_PRECISION
        vsFloor(n, a, y);
#elif DOUBLE_PRECISION
        vdFloor(n, a, y);
#else
        mkl::mkl_error();
#endif
    }

    static void vCeil(const MKL_INT n, const T* a, T* y)
    {
#ifdef SINGLE_PRECISION
        vsCeil(n, a, y);
#elif DOUBLE_PRECISION
        vdCeil(n, a, y);
#else
        mkl::mkl_error();
#endif
    }

    static void vTrunc(const MKL_INT n, const T* a, T* y)
    {
#ifdef SINGLE_PRECISION
        vsTrunc(n, a, y);
#elif DOUBLE_PRECISION
        vdTrunc(n, a, y);
#else
        mkl::mkl_error();
#endif
    }

    static void vRound(const MKL_INT n, const T* a, T* y)
    {
#ifdef SINGLE_PRECISION
        vsRound(n, a, y);
#elif DOUBLE_PRECISION
        vdRound(n, a, y);
#else
        mkl::mkl_error();
#endif
    }

    static void vNearbyInt(const MKL_INT n, const T* a, T* y)
    {
#ifdef SINGLE_PRECISION
        vsNearbyInt(n, a, y);
#elif DOUBLE_PRECISION
        vdNearbyInt(n, a, y);
#else
        mkl::mkl_error();
#endif
    }

    static void vModf(const MKL_INT n, const T* a, T* y, T* z)
    {
#ifdef SINGLE_PRECISION
        vsModf(n, a, y, z);
#elif DOUBLE_PRECISION
        vdModf(n, a, y, z);
#else
        mkl::mkl_error();
#endif
    }

    static void vFmax(const MKL_INT n, const T* a, const T* b, T* y)
    {
#ifdef SINGLE_PRECISION
        vsFmax(n, a, b, y);
#elif DOUBLE_PRECISION
        vdFmax(n, a, b, y);
#else
        mkl::mkl_error();
#endif
    }

    static void vFmin(const MKL_INT n, const T* a, const T* b, T* y)
    {
#ifdef SINGLE_PRECISION
        vsFmin(n, a, b, y);
#elif DOUBLE_PRECISION
        vdFmin(n, a, b, y);
#else
        mkl::mkl_error();
#endif
    }

    static void vMaxMag(const MKL_INT n, const T* a, const T* b, T* y)
    {
#ifdef SINGLE_PRECISION
        vsMaxMag(n, a, b, y);
#elif DOUBLE_PRECISION
        vdMaxMag(n, a, b, y);
#else
        mkl::mkl_error();
#endif
    }

    static void vMinMag(const MKL_INT n, const T* a, const T* b, T* y)
    {
#ifdef SINGLE_PRECISION
        vsMinMag(n, a, b, y);
#elif DOUBLE_PRECISION
        vdMinMag(n, a, b, y);
#else
        mkl::mkl_error();
#endif
    }

    // Data Fitting

    static int dfNewTask1D(DFTaskPtr* task, const MKL_INT nx, const T* x, const MKL_INT xhint,
                           const MKL_INT ny, const T* y, const MKL_INT yhint)
    {
#ifdef SINGLE_PRECISION
        return dfsNewTask1D(task, nx, x, xhint, ny, y, yhint);
#elif DOUBLE_PRECISION
        return dfdNewTask1D(task, nx, x, xhint, ny, y, yhint);
#else
        mkl::mkl_error();
#endif
    }

    static int dfEditPPSpline1D(DFTaskPtr task, const MKL_INT s_order, const MKL_INT s_type,
                                const MKL_INT bc_type, const T* bc, const MKL_INT ic_type,
                                const T* ic, const T* scoeff, const MKL_INT scoeffhint)
    {
#ifdef SINGLE_PRECISION
        return dfsEditPPSpline1D(task, s_order, s_type, bc_type, bc, ic_type, ic, scoeff,
                                 scoeffhint);
#elif DOUBLE_PRECISION
        return dfdEditPPSpline1D(task, s_order, s_type, bc_type, bc, ic_type, ic, scoeff,
                                 scoeffhint);
#else
        mkl::mkl_error();
#endif
    }

    static int dfConstruct1D(DFTaskPtr task, const MKL_INT s_format, const MKL_INT method)
    {
#ifdef SINGLE_PRECISION
        return dfsConstruct1D(task, s_format, method);
#elif DOUBLE_PRECISION
        return dfsConstruct1D(task, s_format, method);
#else
        mkl::mkl_error();
#endif
    }

    static int dfInterpolate1D(DFTaskPtr task, const MKL_INT type, const MKL_INT method,
                               const MKL_INT nsite, const T* site, const MKL_INT sitehint,
                               const MKL_INT ndorder, const MKL_INT* dorder, const T* datahint,
                               T* r, const MKL_INT rhint, MKL_INT* cell)
    {
#ifdef SINGLE_PRECISION
        return dfsInterpolate1D(task, type, method, nsite, site, sitehint, ndorder, dorder,
                                datahint, r, rhint, cell);
#elif DOUBLE_PRECISION
        return dfdInterpolate1D(task, type, method, nsite, site, sitehint, ndorder, dorder,
                                datahint, r, rhint, cell);
#else
        mkl::mkl_error();
#endif
    }

    static int dfInterpolateEx1D(DFTaskPtr task, const MKL_INT type, const MKL_INT method,
                                 const MKL_INT nsite, const T* site, const MKL_INT sitehint,
                                 const MKL_INT ndorder, const MKL_INT* dorder, const T* datahint,
                                 T* r, const MKL_INT rhint, MKL_INT* cell, const void* le_params,
                                 const void* re_params, const void* i_params,
                                 const void* search_params)
    {
        /* This needs modification to incorporate callbacks */
#ifdef SINGLE_PRECISION
        return dfsInterpolateEx1D(task, type, method, nsite, site, sitehint, ndorder, dorder,
                                  datahint, r, rhint, cell, NULL, le_params, NULL, re_params, NULL,
                                  i_params, NULL, search_params);
#elif DOUBLE_PRECISION
        return dfdInterpolateEx1D(task, type, method, nsite, site, sitehint, ndorder, dorder,
                                  datahint, r, rhint, cell, NULL, le_params, NULL, re_params, NULL,
                                  i_params, NULL, search_params);
#else
        mkl::mkl_error();
#endif
    }
};
