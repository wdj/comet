/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions normal z -> s d c
*/
#include "common_magma_minproduct.h"

/**
    Purpose
    -------
    ZPOSV computes the solution to a complex system of linear equations
       A * X = B,
    where A is an N-by-N Hermitian positive definite matrix and X and B
    are N-by-NRHS matrices.
    The Cholesky decomposition is used to factor A as
       A = U**H * U,  if UPLO = Magma_minproductUpper, or
       A = L * L**H,  if UPLO = Magma_minproductLower,
    where U is an upper triangular matrix and  L is a lower triangular
    matrix.  The factored form of A is then used to solve the system of
    equations A * X = B.

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

    @param[in,out]
    dA      COMPLEX_16 array on the GPU, dimension (LDDA,N)
            On entry, the Hermitian matrix dA.  If UPLO = Magma_minproductUpper, the leading
            N-by-N upper triangular part of dA contains the upper
            triangular part of the matrix dA, and the strictly lower
            triangular part of dA is not referenced.  If UPLO = Magma_minproductLower, the
            leading N-by-N lower triangular part of dA contains the lower
            triangular part of the matrix dA, and the strictly upper
            triangular part of dA is not referenced.
    \n
            On exit, if INFO = 0, the factor U or L from the Cholesky
            factorization dA = U**H*U or dA = L*L**H.

    @param[in]
    ldda    INTEGER
            The leading dimension of the array A.  LDA >= max(1,N).

    @param[in,out]
    dB      COMPLEX_16 array on the GPU, dimension (LDB,NRHS)
            On entry, the right hand side matrix B.
            On exit, the solution matrix X.

    @param[in]
    lddb    INTEGER
            The leading dimension of the array B.  LDB >= max(1,N).

    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument had an illegal value

    @ingroup magma_minproduct_zposv_driver
    ********************************************************************/
extern "C" magma_minproduct_int_t
magma_minproduct_zposv_gpu(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex_ptr dB, magma_minproduct_int_t lddb,
    magma_minproduct_int_t *info )
{
    *info = 0;
    if ( uplo != Magma_minproductUpper && uplo != Magma_minproductLower )
        *info = -1;
    if ( n < 0 )
        *info = -2;
    if ( nrhs < 0 )
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

    magma_minproduct_zpotrf_gpu( uplo, n, dA, ldda, info );
    if ( *info == 0 ) {
        magma_minproduct_zpotrs_gpu( uplo, n, nrhs, dA, ldda, dB, lddb, info );
    }

    return *info;
}
