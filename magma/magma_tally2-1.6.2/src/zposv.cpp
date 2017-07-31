/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions normal z -> s d c

*/
#include "common_magma_tally2.h"

/**
    Purpose
    -------
    ZPOSV computes the solution to a complex system of linear equations
       A * X = B,
    where A is an N-by-N Hermitian positive definite matrix and X and B
    are N-by-NRHS matrices.
    The Cholesky decomposition is used to factor A as
       A = U**H * U,  if UPLO = Magma_tally2Upper, or
       A = L * L**H,  if UPLO = Magma_tally2Lower,
    where U is an upper triangular matrix and  L is a lower triangular
    matrix.  The factored form of A is then used to solve the system of
    equations A * X = B.

    Arguments
    ---------
    @param[in]
    uplo    magma_tally2_uplo_t
      -     = Magma_tally2Upper:  Upper triangle of A is stored;
      -     = Magma_tally2Lower:  Lower triangle of A is stored.

    @param[in]
    n       INTEGER
            The order of the matrix A.  N >= 0.

    @param[in]
    nrhs    INTEGER
            The number of right hand sides, i.e., the number of columns
            of the matrix B.  NRHS >= 0.

    @param[in,out]
    A       COMPLEX_16 array, dimension (LDA,N)
            On entry, the Hermitian matrix A.  If UPLO = Magma_tally2Upper, the leading
            N-by-N upper triangular part of A contains the upper
            triangular part of the matrix A, and the strictly lower
            triangular part of A is not referenced.  If UPLO = Magma_tally2Lower, the
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

    @ingroup magma_tally2_zposv_driver
    ********************************************************************/
extern "C" magma_tally2_int_t
magma_tally2_zposv(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n, magma_tally2_int_t nrhs,
    magma_tally2DoubleComplex *A, magma_tally2_int_t lda,
    magma_tally2DoubleComplex *B, magma_tally2_int_t ldb,
    magma_tally2_int_t *info )
{
    magma_tally2_int_t ngpu, ldda, lddb;

    *info = 0;
    if ( uplo != Magma_tally2Upper && uplo != Magma_tally2Lower )
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
        magma_tally2_xerbla( __func__, -(*info) );
        return *info;
    }

    /* Quick return if possible */
    if ( (n == 0) || (nrhs == 0) ) {
        return *info;
    }

    /* If single-GPU and allocation suceeds, use GPU interface. */
    ngpu = magma_tally2_num_gpus();
    magma_tally2DoubleComplex *dA, *dB;
    if ( ngpu > 1 ) {
        goto CPU_INTERFACE;
    }
    ldda = ((n+31)/32)*32;
    lddb = ldda;
    if ( MAGMA_tally2_SUCCESS != magma_tally2_zmalloc( &dA, ldda*n )) {
        goto CPU_INTERFACE;
    }
    if ( MAGMA_tally2_SUCCESS != magma_tally2_zmalloc( &dB, lddb*nrhs )) {
        magma_tally2_free( dA );
        goto CPU_INTERFACE;
    }
    magma_tally2_zsetmatrix( n, n, A, lda, dA, ldda );
    magma_tally2_zpotrf_gpu( uplo, n, dA, ldda, info );
    if ( *info == MAGMA_tally2_ERR_DEVICE_ALLOC ) {
        magma_tally2_free( dA );
        magma_tally2_free( dB );
        goto CPU_INTERFACE;
    }
    magma_tally2_zgetmatrix( n, n, dA, ldda, A, lda );
    if ( *info == 0 ) {
        magma_tally2_zsetmatrix( n, nrhs, B, ldb, dB, lddb );
        magma_tally2_zpotrs_gpu( uplo, n, nrhs, dA, ldda, dB, lddb, info );
        magma_tally2_zgetmatrix( n, nrhs, dB, lddb, B, ldb );
    }
    magma_tally2_free( dA );
    magma_tally2_free( dB );
    return *info;

CPU_INTERFACE:
    /* If multi-GPU or allocation failed, use CPU interface and LAPACK.
     * Faster to use LAPACK for potrs than to copy A to GPU. */
    magma_tally2_zpotrf( uplo, n, A, lda, info );
    if ( *info == 0 ) {
        lapackf77_zpotrs( lapack_uplo_const_tally2(uplo), &n, &nrhs, A, &lda, B, &ldb, info );
    }

    return *info;
}
