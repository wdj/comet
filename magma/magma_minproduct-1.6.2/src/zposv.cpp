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
    A       COMPLEX_16 array, dimension (LDA,N)
            On entry, the Hermitian matrix A.  If UPLO = Magma_minproductUpper, the leading
            N-by-N upper triangular part of A contains the upper
            triangular part of the matrix A, and the strictly lower
            triangular part of A is not referenced.  If UPLO = Magma_minproductLower, the
            leading N-by-N lower triangular part of A contains the lower
            triangular part of the matrix A, and the strictly upper
            triangular part of A is not referenced.
    \n
            On exit, if INFO = 0, the factor U or L from the Cholesky
            factorization A = U**H*U or A = L*L**H.

    @param[in]
    lda     INTEGER
            The leading dimension of the array A.  LDA >= max(1,N).

    @param[in,out]
    B       COMPLEX_16 array, dimension (LDB,NRHS)
            On entry, the right hand side matrix B.
            On exit, the solution matrix X.

    @param[in]
    ldb     INTEGER
            The leading dimension of the array B.  LDB >= max(1,N).

    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument had an illegal value

    @ingroup magma_minproduct_zposv_driver
    ********************************************************************/
extern "C" magma_minproduct_int_t
magma_minproduct_zposv(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
    magma_minproductDoubleComplex *A, magma_minproduct_int_t lda,
    magma_minproductDoubleComplex *B, magma_minproduct_int_t ldb,
    magma_minproduct_int_t *info )
{
    magma_minproduct_int_t ngpu, ldda, lddb;

    *info = 0;
    if ( uplo != Magma_minproductUpper && uplo != Magma_minproductLower )
        *info = -1;
    if ( n < 0 )
        *info = -2;
    if ( nrhs < 0)
        *info = -3;
    if ( lda < max(1, n) )
        *info = -5;
    if ( ldb < max(1, n) )
        *info = -7;
    if (*info != 0) {
        magma_minproduct_xerbla( __func__, -(*info) );
        return *info;
    }

    /* Quick return if possible */
    if ( (n == 0) || (nrhs == 0) ) {
        return *info;
    }

    /* If single-GPU and allocation suceeds, use GPU interface. */
    ngpu = magma_minproduct_num_gpus();
    magma_minproductDoubleComplex *dA, *dB;
    if ( ngpu > 1 ) {
        goto CPU_INTERFACE;
    }
    ldda = ((n+31)/32)*32;
    lddb = ldda;
    if ( MAGMA_minproduct_SUCCESS != magma_minproduct_zmalloc( &dA, ldda*n )) {
        goto CPU_INTERFACE;
    }
    if ( MAGMA_minproduct_SUCCESS != magma_minproduct_zmalloc( &dB, lddb*nrhs )) {
        magma_minproduct_free( dA );
        goto CPU_INTERFACE;
    }
    magma_minproduct_zsetmatrix( n, n, A, lda, dA, ldda );
    magma_minproduct_zpotrf_gpu( uplo, n, dA, ldda, info );
    if ( *info == MAGMA_minproduct_ERR_DEVICE_ALLOC ) {
        magma_minproduct_free( dA );
        magma_minproduct_free( dB );
        goto CPU_INTERFACE;
    }
    magma_minproduct_zgetmatrix( n, n, dA, ldda, A, lda );
    if ( *info == 0 ) {
        magma_minproduct_zsetmatrix( n, nrhs, B, ldb, dB, lddb );
        magma_minproduct_zpotrs_gpu( uplo, n, nrhs, dA, ldda, dB, lddb, info );
        magma_minproduct_zgetmatrix( n, nrhs, dB, lddb, B, ldb );
    }
    magma_minproduct_free( dA );
    magma_minproduct_free( dB );
    return *info;

CPU_INTERFACE:
    /* If multi-GPU or allocation failed, use CPU interface and LAPACK.
     * Faster to use LAPACK for potrs than to copy A to GPU. */
    magma_minproduct_zpotrf( uplo, n, A, lda, info );
    if ( *info == 0 ) {
        lapackf77_zpotrs( lapack_uplo_const(uplo), &n, &nrhs, A, &lda, B, &ldb, info );
    }

    return *info;
}
