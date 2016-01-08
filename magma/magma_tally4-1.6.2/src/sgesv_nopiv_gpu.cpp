/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015
       @author Adrien REMY

       @generated from zgesv_nopiv_gpu.cpp normal z -> s, Fri Jan 30 19:00:14 2015
*/
#include "common_magma_tally4.h"

/**
    Purpose
    -------
    Solves a system of linear equations
        A * X = B
    where A is a general n-by-n matrix and X and B are n-by-nrhs matrices.
    The LU decomposition with no pivoting is
    used to factor A as
        A = L * U,
    where L is unit lower triangular, and U is
    upper triangular.  The factored form of A is then used to solve the
    system of equations A * X = B.

    Arguments
    ---------
    @param[in]
    n       INTEGER
            The order of the matrix A.  n >= 0.

    @param[in]
    nrhs    INTEGER
            The number of right hand sides, i.e., the number of columns
            of the matrix B.  nrhs >= 0.

    @param[in,out]
    dA      REAL array on the GPU, dimension (ldda,n).
            On entry, the n-by-n matrix to be factored.
            On exit, the factors L and U from the factorization
            A = L*U; the unit diagonal elements of L are not stored.

    @param[in]
    ldda    INTEGER
            The leading dimension of the array A.  ldda >= max(1,n).

    @param[in,out]
    dB      REAL array on the GPU, dimension (lddb,nrhs)
            On entry, the right hand side matrix B.
            On exit, the solution matrix X.

    @param[in]
    lddb    INTEGER
            The leading dimension of the array B.  lddb >= max(1,n).

    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument had an illegal value

    @ingroup magma_tally4_sgesv_driver
    ********************************************************************/

extern "C" magma_tally4_int_t
magma_tally4_sgesv_nopiv_gpu(
    magma_tally4_int_t n, magma_tally4_int_t nrhs,
    magma_tally4Float_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Float_ptr dB, magma_tally4_int_t lddb,
    magma_tally4_int_t *info )
{
    *info = 0;
    if (n < 0) {
        *info = -1;
    } else if (nrhs < 0) {
        *info = -2;
    } else if (ldda < max(1,n)) {
        *info = -4;
    } else if (lddb < max(1,n)) {
        *info = -6;
    }
    if (*info != 0) {
        magma_tally4_xerbla( __func__, -(*info) );
        return *info;
    }

    /* Quick return if possible */
    if (n == 0 || nrhs == 0) {
        return *info;
    }

    magma_tally4_sgetrf_nopiv_gpu( n, n, dA, ldda, info );
    if ( *info == MAGMA_tally4_SUCCESS ) {
        magma_tally4_sgetrs_nopiv_gpu( Magma_tally4NoTrans, n, nrhs, dA, ldda, dB, lddb, info );
    }
    
    return *info;
}
