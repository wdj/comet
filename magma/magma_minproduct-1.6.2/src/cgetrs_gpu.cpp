/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from zgetrs_gpu.cpp normal z -> c, Fri Jan 30 19:00:14 2015

*/
#include "common_magma_minproduct.h"

/**
    Purpose
    -------
    Solves a system of linear equations
      A * X = B,  A**T * X = B,  or  A**H * X = B
    with a general N-by-N matrix A using the LU factorization computed by CGETRF_GPU.

    Arguments
    ---------
    @param[in]
    trans   magma_minproduct_trans_t
            Specifies the form of the system of equations:
      -     = Magma_minproductNoTrans:    A    * X = B  (No transpose)
      -     = Magma_minproductTrans:      A**T * X = B  (Transpose)
      -     = Magma_minproductConjTrans:  A**H * X = B  (Conjugate transpose)

    @param[in]
    n       INTEGER
            The order of the matrix A.  N >= 0.

    @param[in]
    nrhs    INTEGER
            The number of right hand sides, i.e., the number of columns
            of the matrix B.  NRHS >= 0.

    @param[in]
    dA      COMPLEX array on the GPU, dimension (LDA,N)
            The factors L and U from the factorization A = P*L*U as computed
            by CGETRF_GPU.

    @param[in]
    ldda    INTEGER
            The leading dimension of the array A.  LDA >= max(1,N).

    @param[in]
    ipiv    INTEGER array, dimension (N)
            The pivot indices from CGETRF; for 1 <= i <= N, row i of the
            matrix was interchanged with row IPIV(i).

    @param[in,out]
    dB      COMPLEX array on the GPU, dimension (LDB,NRHS)
            On entry, the right hand side matrix B.
            On exit, the solution matrix X.

    @param[in]
    lddb    INTEGER
            The leading dimension of the array B.  LDB >= max(1,N).

    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument had an illegal value

    @ingroup magma_minproduct_cgesv_comp
    ********************************************************************/
extern "C" magma_minproduct_int_t
magma_minproduct_cgetrs_gpu(
    magma_minproduct_trans_t trans, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda, magma_minproduct_int_t *ipiv,
    magma_minproductFloatComplex_ptr dB, magma_minproduct_int_t lddb,
    magma_minproduct_int_t *info)
{
    magma_minproductFloatComplex c_one = MAGMA_minproduct_C_ONE;
    magma_minproductFloatComplex *work = NULL;
    int notran = (trans == Magma_minproductNoTrans);
    magma_minproduct_int_t i1, i2, inc;

    *info = 0;
    if ( (! notran) &&
         (trans != Magma_minproductTrans) &&
         (trans != Magma_minproductConjTrans) ) {
        *info = -1;
    } else if (n < 0) {
        *info = -2;
    } else if (nrhs < 0) {
        *info = -3;
    } else if (ldda < max(1,n)) {
        *info = -5;
    } else if (lddb < max(1,n)) {
        *info = -8;
    }
    if (*info != 0) {
        magma_minproduct_xerbla( __func__, -(*info) );
        return *info;
    }

    /* Quick return if possible */
    if (n == 0 || nrhs == 0) {
        return *info;
    }

    magma_minproduct_cmalloc_cpu( &work, n * nrhs );
    if ( work == NULL ) {
        *info = MAGMA_minproduct_ERR_HOST_ALLOC;
        return *info;
    }
      
    i1 = 1;
    i2 = n;
    if (notran) {
        inc = 1;

        /* Solve A * X = B. */
        magma_minproduct_cgetmatrix( n, nrhs, dB, lddb, work, n );
        lapackf77_claswp(&nrhs, work, &n, &i1, &i2, ipiv, &inc);
        magma_minproduct_csetmatrix( n, nrhs, work, n, dB, lddb );

        if ( nrhs == 1) {
            magma_minproduct_ctrsv(Magma_minproductLower, Magma_minproductNoTrans, Magma_minproductUnit,    n, dA, ldda, dB, 1 );
            magma_minproduct_ctrsv(Magma_minproductUpper, Magma_minproductNoTrans, Magma_minproductNonUnit, n, dA, ldda, dB, 1 );
        } else {
            magma_minproduct_ctrsm(Magma_minproductLeft, Magma_minproductLower, Magma_minproductNoTrans, Magma_minproductUnit,    n, nrhs, c_one, dA, ldda, dB, lddb );
            magma_minproduct_ctrsm(Magma_minproductLeft, Magma_minproductUpper, Magma_minproductNoTrans, Magma_minproductNonUnit, n, nrhs, c_one, dA, ldda, dB, lddb );
        }
    } else {
        inc = -1;

        /* Solve A**T * X = B  or  A**H * X = B. */
        if ( nrhs == 1) {
            magma_minproduct_ctrsv(Magma_minproductUpper, trans, Magma_minproductNonUnit, n, dA, ldda, dB, 1 );
            magma_minproduct_ctrsv(Magma_minproductLower, trans, Magma_minproductUnit,    n, dA, ldda, dB, 1 );
        } else {
            magma_minproduct_ctrsm(Magma_minproductLeft, Magma_minproductUpper, trans, Magma_minproductNonUnit, n, nrhs, c_one, dA, ldda, dB, lddb );
            magma_minproduct_ctrsm(Magma_minproductLeft, Magma_minproductLower, trans, Magma_minproductUnit,    n, nrhs, c_one, dA, ldda, dB, lddb );
        }

        magma_minproduct_cgetmatrix( n, nrhs, dB, lddb, work, n );
        lapackf77_claswp(&nrhs, work, &n, &i1, &i2, ipiv, &inc);
        magma_minproduct_csetmatrix( n, nrhs, work, n, dB, lddb );
    }
    magma_minproduct_free_cpu(work);

    return *info;
}
