/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from zgetrs_gpu.cpp normal z -> c, Fri Jan 30 19:00:14 2015

*/
#include "common_magma_tally4.h"

/**
    Purpose
    -------
    Solves a system of linear equations
      A * X = B,  A**T * X = B,  or  A**H * X = B
    with a general N-by-N matrix A using the LU factorization computed by CGETRF_GPU.

    Arguments
    ---------
    @param[in]
    trans   magma_tally4_trans_t
            Specifies the form of the system of equations:
      -     = Magma_tally4NoTrans:    A    * X = B  (No transpose)
      -     = Magma_tally4Trans:      A**T * X = B  (Transpose)
      -     = Magma_tally4ConjTrans:  A**H * X = B  (Conjugate transpose)

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

    @ingroup magma_tally4_cgesv_comp
    ********************************************************************/
extern "C" magma_tally4_int_t
magma_tally4_cgetrs_gpu(
    magma_tally4_trans_t trans, magma_tally4_int_t n, magma_tally4_int_t nrhs,
    magma_tally4FloatComplex_ptr dA, magma_tally4_int_t ldda, magma_tally4_int_t *ipiv,
    magma_tally4FloatComplex_ptr dB, magma_tally4_int_t lddb,
    magma_tally4_int_t *info)
{
    magma_tally4FloatComplex c_one = MAGMA_tally4_C_ONE;
    magma_tally4FloatComplex *work = NULL;
    int notran = (trans == Magma_tally4NoTrans);
    magma_tally4_int_t i1, i2, inc;

    *info = 0;
    if ( (! notran) &&
         (trans != Magma_tally4Trans) &&
         (trans != Magma_tally4ConjTrans) ) {
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
        magma_tally4_xerbla( __func__, -(*info) );
        return *info;
    }

    /* Quick return if possible */
    if (n == 0 || nrhs == 0) {
        return *info;
    }

    magma_tally4_cmalloc_cpu( &work, n * nrhs );
    if ( work == NULL ) {
        *info = MAGMA_tally4_ERR_HOST_ALLOC;
        return *info;
    }
      
    i1 = 1;
    i2 = n;
    if (notran) {
        inc = 1;

        /* Solve A * X = B. */
        magma_tally4_cgetmatrix( n, nrhs, dB, lddb, work, n );
        lapackf77_claswp(&nrhs, work, &n, &i1, &i2, ipiv, &inc);
        magma_tally4_csetmatrix( n, nrhs, work, n, dB, lddb );

        if ( nrhs == 1) {
            magma_tally4_ctrsv(Magma_tally4Lower, Magma_tally4NoTrans, Magma_tally4Unit,    n, dA, ldda, dB, 1 );
            magma_tally4_ctrsv(Magma_tally4Upper, Magma_tally4NoTrans, Magma_tally4NonUnit, n, dA, ldda, dB, 1 );
        } else {
            magma_tally4_ctrsm(Magma_tally4Left, Magma_tally4Lower, Magma_tally4NoTrans, Magma_tally4Unit,    n, nrhs, c_one, dA, ldda, dB, lddb );
            magma_tally4_ctrsm(Magma_tally4Left, Magma_tally4Upper, Magma_tally4NoTrans, Magma_tally4NonUnit, n, nrhs, c_one, dA, ldda, dB, lddb );
        }
    } else {
        inc = -1;

        /* Solve A**T * X = B  or  A**H * X = B. */
        if ( nrhs == 1) {
            magma_tally4_ctrsv(Magma_tally4Upper, trans, Magma_tally4NonUnit, n, dA, ldda, dB, 1 );
            magma_tally4_ctrsv(Magma_tally4Lower, trans, Magma_tally4Unit,    n, dA, ldda, dB, 1 );
        } else {
            magma_tally4_ctrsm(Magma_tally4Left, Magma_tally4Upper, trans, Magma_tally4NonUnit, n, nrhs, c_one, dA, ldda, dB, lddb );
            magma_tally4_ctrsm(Magma_tally4Left, Magma_tally4Lower, trans, Magma_tally4Unit,    n, nrhs, c_one, dA, ldda, dB, lddb );
        }

        magma_tally4_cgetmatrix( n, nrhs, dB, lddb, work, n );
        lapackf77_claswp(&nrhs, work, &n, &i1, &i2, ipiv, &inc);
        magma_tally4_csetmatrix( n, nrhs, work, n, dB, lddb );
    }
    magma_tally4_free_cpu(work);

    return *info;
}
