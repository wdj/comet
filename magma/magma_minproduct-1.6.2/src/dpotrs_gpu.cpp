/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from zpotrs_gpu.cpp normal z -> d, Fri Jan 30 19:00:12 2015

*/
#include "common_magma_minproduct.h"

/**
    Purpose
    -------
    DPOTRS solves a system of linear equations A*X = B with a symmetric
    positive definite matrix A using the Cholesky factorization
    A = U**H*U or A = L*L**H computed by DPOTRF.

    Arguments
    ---------
    @param[in]
    uplo    magma_minproduct_uplo_t
      -     = Magma_minproductUpper:  Upper triangle of A is stored;
      -     = Magma_minproductLower:  Lower triangle of A is stored.

    @param[in]
    n       INTEGER
            The order of the matrix A.  N >= 0.

    @param[in]
    nrhs    INTEGER
            The number of right hand sides, i.e., the number of columns
            of the matrix B.  NRHS >= 0.

    @param[in]
    dA      DOUBLE_PRECISION array on the GPU, dimension (LDDA,N)
            The triangular factor U or L from the Cholesky factorization
            A = U**H*U or A = L*L**H, as computed by DPOTRF.

    @param[in]
    ldda    INTEGER
            The leading dimension of the array A.  LDDA >= max(1,N).

    @param[in,out]
    dB      DOUBLE_PRECISION array on the GPU, dimension (LDDB,NRHS)
            On entry, the right hand side matrix B.
            On exit, the solution matrix X.

    @param[in]
    lddb    INTEGER
            The leading dimension of the array B.  LDDB >= max(1,N).

    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument had an illegal value

    @ingroup magma_minproduct_dposv_comp
    ********************************************************************/
extern "C" magma_minproduct_int_t
magma_minproduct_dpotrs_gpu(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDouble_ptr dB, magma_minproduct_int_t lddb,
    magma_minproduct_int_t *info)
{
    double c_one = MAGMA_minproduct_D_ONE;

    *info = 0;
    if ( uplo != Magma_minproductUpper && uplo != Magma_minproductLower )
        *info = -1;
    if ( n < 0 )
        *info = -2;
    if ( nrhs < 0)
        *info = -3;
    if ( ldda < max(1, n) )
        *info = -5;
    if ( lddb < max(1, n) )
        *info = -7;
    if (*info != 0) {
        magma_minproduct_xerbla( __func__, -(*info) );
        return *info;
    }

    /* Quick return if possible */
    if ( (n == 0) || (nrhs == 0) ) {
        return *info;
    }

    if ( uplo == Magma_minproductUpper ) {
        if ( nrhs == 1) {
            magma_minproduct_dtrsv(Magma_minproductUpper, Magma_minproductConjTrans, Magma_minproductNonUnit, n, dA, ldda, dB, 1 );
            magma_minproduct_dtrsv(Magma_minproductUpper, Magma_minproductNoTrans,   Magma_minproductNonUnit, n, dA, ldda, dB, 1 );
        } else {
            magma_minproduct_dtrsm(Magma_minproductLeft, Magma_minproductUpper, Magma_minproductConjTrans, Magma_minproductNonUnit, n, nrhs, c_one, dA, ldda, dB, lddb);
            magma_minproduct_dtrsm(Magma_minproductLeft, Magma_minproductUpper, Magma_minproductNoTrans,   Magma_minproductNonUnit, n, nrhs, c_one, dA, ldda, dB, lddb);
        }
    }
    else {
        if ( nrhs == 1) {
            magma_minproduct_dtrsv(Magma_minproductLower, Magma_minproductNoTrans,   Magma_minproductNonUnit, n, dA, ldda, dB, 1 );
            magma_minproduct_dtrsv(Magma_minproductLower, Magma_minproductConjTrans, Magma_minproductNonUnit, n, dA, ldda, dB, 1 );
        } else {
            magma_minproduct_dtrsm(Magma_minproductLeft, Magma_minproductLower, Magma_minproductNoTrans,   Magma_minproductNonUnit, n, nrhs, c_one, dA, ldda, dB, lddb);
            magma_minproduct_dtrsm(Magma_minproductLeft, Magma_minproductLower, Magma_minproductConjTrans, Magma_minproductNonUnit, n, nrhs, c_one, dA, ldda, dB, lddb);
        }
    }

    return *info;
}
