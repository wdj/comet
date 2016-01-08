/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from zgesv.cpp normal z -> c, Fri Jan 30 19:00:14 2015

*/
#include "common_magma_tally4.h"

/**
    Purpose
    -------
    Solves a system of linear equations
       A * X = B
    where A is a general N-by-N matrix and X and B are N-by-NRHS matrices.
    The LU decomposition with partial pivoting and row interchanges is
    used to factor A as
       A = P * L * U,
    where P is a permutation matrix, L is unit lower triangular, and U is
    upper triangular.  The factored form of A is then used to solve the
    system of equations A * X = B.

    Arguments
    ---------
    @param[in]
    n       INTEGER
            The order of the matrix A.  N >= 0.

    @param[in]
    nrhs    INTEGER
            The number of right hand sides, i.e., the number of columns
            of the matrix B.  NRHS >= 0.

    @param[in,out]
    A       COMPLEX array, dimension (LDA,N).
            On entry, the M-by-N matrix to be factored.
            On exit, the factors L and U from the factorization
            A = P*L*U; the unit diagonal elements of L are not stored.

    @param[in]
    lda     INTEGER
            The leading dimension of the array A.  LDA >= max(1,N).

    @param[out]
    ipiv    INTEGER array, dimension (min(M,N))
            The pivot indices; for 1 <= i <= min(M,N), row i of the
            matrix was interchanged with row IPIV(i).

    @param[in,out]
    B       COMPLEX array, dimension (LDB,NRHS)
            On entry, the right hand side matrix B.
            On exit, the solution matrix X.

    @param[in]
    ldb     INTEGER
            The leading dimension of the array B.  LDB >= max(1,N).

    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument had an illegal value

    @ingroup magma_tally4_cgesv_driver
    ********************************************************************/
extern "C" magma_tally4_int_t
magma_tally4_cgesv(
    magma_tally4_int_t n, magma_tally4_int_t nrhs,
    magma_tally4FloatComplex *A, magma_tally4_int_t lda,
    magma_tally4_int_t *ipiv,
    magma_tally4FloatComplex *B, magma_tally4_int_t ldb,
    magma_tally4_int_t *info)
{
    magma_tally4_int_t ngpu, ldda, lddb;
    
    *info = 0;
    if (n < 0) {
        *info = -1;
    } else if (nrhs < 0) {
        *info = -2;
    } else if (lda < max(1,n)) {
        *info = -4;
    } else if (ldb < max(1,n)) {
        *info = -7;
    }
    if (*info != 0) {
        magma_tally4_xerbla( __func__, -(*info) );
        return *info;
    }

    /* Quick return if possible */
    if (n == 0 || nrhs == 0) {
        return *info;
    }
    
    /* If single-GPU and allocation suceeds, use GPU interface. */
    ngpu = magma_tally4_num_gpus();
    magma_tally4FloatComplex *dA, *dB;
    if ( ngpu > 1 ) {
        goto CPU_INTERFACE;
    }
    ldda = ((n+31)/32)*32;
    lddb = ldda;
    if ( MAGMA_tally4_SUCCESS != magma_tally4_cmalloc( &dA, ldda*n )) {
        goto CPU_INTERFACE;
    }
    if ( MAGMA_tally4_SUCCESS != magma_tally4_cmalloc( &dB, lddb*nrhs )) {
        magma_tally4_free( dA );
        goto CPU_INTERFACE;
    }
    magma_tally4_csetmatrix( n, n, A, lda, dA, ldda );
    magma_tally4_cgetrf_gpu( n, n, dA, ldda, ipiv, info );
    if ( *info == MAGMA_tally4_ERR_DEVICE_ALLOC ) {
        magma_tally4_free( dA );
        magma_tally4_free( dB );
        goto CPU_INTERFACE;
    }
    magma_tally4_cgetmatrix( n, n, dA, ldda, A, lda );
    if ( *info == 0 ) {
        magma_tally4_csetmatrix( n, nrhs, B, ldb, dB, lddb );
        magma_tally4_cgetrs_gpu( Magma_tally4NoTrans, n, nrhs, dA, ldda, ipiv, dB, lddb, info );
        magma_tally4_cgetmatrix( n, nrhs, dB, lddb, B, ldb );
    }
    magma_tally4_free( dA );
    magma_tally4_free( dB );
    return *info;

CPU_INTERFACE:
    /* If multi-GPU or allocation failed, use CPU interface and LAPACK.
     * Faster to use LAPACK for getrs than to copy A to GPU. */
    magma_tally4_cgetrf( n, n, A, lda, ipiv, info );
    if ( *info == 0 ) {
        lapackf77_cgetrs( Magma_tally4NoTransStr, &n, &nrhs, A, &lda, ipiv, B, &ldb, info );
    }
    return *info;
}
