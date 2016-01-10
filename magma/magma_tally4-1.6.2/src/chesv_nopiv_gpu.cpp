/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015
       @author Adrien REMY

       @generated from zhesv_nopiv_gpu.cpp normal z -> c, Fri Jan 30 19:00:16 2015

*/
#include "common_magma_tally4.h"

/**
    Purpose
    -------
    Solves a system of linear equations
       A * X = B
    where A is an n-by-n hermitian matrix and X and B are n-by-nrhs matrices.
    The LU decomposition with no pivoting is
    used to factor A as
    The factorization has the form   
       A = U^H * D * U , if UPLO = 'U', or   
       A = L  * D * L^H, if UPLO = 'L',   
    where U is an upper triangular matrix, L is lower triangular, and
    D is a diagonal matrix.
    The factored form of A is then
    used to solve the system of equations A * X = B.
    
    Arguments
    ---------
    @param[in]
    uplo    magma_tally4_uplo_t
      -     = Magma_tally4Upper:  Upper triangle of A is stored;
      -     = Magma_tally4Lower:  Lower triangle of A is stored.

    @param[in]
    n       INTEGER
            The order of the matrix A.  n >= 0.

    @param[in]
    nrhs    INTEGER
            The number of right hand sides, i.e., the number of columns
            of the matrix B.  nrhs >= 0.

    @param[in,out]
    dA       COMPLEX array, dimension (ldda,n).
            On entry, the n-by-n matrix to be factored.
            On exit, the factors L and U from the factorization
            A = L*U; the unit diagonal elements of L are not stored.

    @param[in]
    ldda     INTEGER
            The leading dimension of the array A.  ldda >= max(1,n).

    @param[in,out]
    dB       COMPLEX array, dimension (lddb,nrhs)
            On entry, the right hand side matrix B.
            On exit, the solution matrix X.

    @param[in]
    lddb     INTEGER
            The leading dimension of the array B.  ldb >= max(1,n).

    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument had an illegal value

    @ingroup magma_tally4_chesv_driver
    ********************************************************************/




extern "C" magma_tally4_int_t
magma_tally4_chesv_nopiv_gpu(magma_tally4_uplo_t uplo,  magma_tally4_int_t n, magma_tally4_int_t nrhs, 
                 magma_tally4FloatComplex_ptr dA, magma_tally4_int_t ldda,
                 magma_tally4FloatComplex_ptr dB, magma_tally4_int_t lddb, 
                 magma_tally4_int_t *info)
{
    magma_tally4_int_t ret;

    *info = 0;
    int   upper = (uplo == Magma_tally4Upper);
    if (! upper && uplo != Magma_tally4Lower) {
      *info = -1;
    }else if (n < 0) {
        *info = -2;
    } else if (nrhs < 0) {
        *info = -3;
    } else if (ldda < max(1,n)) {
        *info = -5;
    } else if (lddb < max(1,n)) {
        *info = -7;
    }
    if (*info != 0) {
        magma_tally4_xerbla( __func__, -(*info) );
        return MAGMA_tally4_ERR_ILLEGAL_VALUE;
    }

    /* Quick return if possible */
    if (n == 0 || nrhs == 0) {
        return MAGMA_tally4_SUCCESS;
    }

    ret = magma_tally4_chetrf_nopiv_tally4_gpu(uplo, n, dA, ldda, info);
    if ( (ret != MAGMA_tally4_SUCCESS) || (*info != 0) ) {
        return ret;
    }
        
    ret = magma_tally4_chetrs_nopiv_gpu( uplo, n, nrhs, dA, ldda, dB, lddb, info );
    if ( (ret != MAGMA_tally4_SUCCESS) || (*info != 0) ) {
        return ret;
    }

    
    return ret;
}