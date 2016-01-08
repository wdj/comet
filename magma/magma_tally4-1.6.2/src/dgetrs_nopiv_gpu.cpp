/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015
       @author Adrien REMY

       @generated from zgetrs_nopiv_gpu.cpp normal z -> d, Fri Jan 30 19:00:14 2015

*/
#include "common_magma_tally4.h"

/**
    Purpose
    -------
    Solves a system of linear equations
      A * X = B,  A**T * X = B,  or  A**H * X = B
    with a general N-by-N matrix A using the LU factorization computed by DGETRF_NOPIV_GPU.

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
    dA      DOUBLE_PRECISION array on the GPU, dimension (LDA,N)
            The factors L and U from the factorization A = P*L*U as computed
            by DGETRF_GPU.

    @param[in]
    ldda    INTEGER
            The leading dimension of the array A.  LDA >= max(1,N).

    param[in,out]
    dB      DOUBLE_PRECISION array on the GPU, dimension (LDB,NRHS)
            On entry, the right hand side matrix B.
            On exit, the solution matrix X.

    @param[in]
    lddb    INTEGER
            The leading dimension of the array B.  LDB >= max(1,N).

    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument had an illegal value

    @ingroup magma_tally4_dgesv_comp
    ********************************************************************/
extern "C" magma_tally4_int_t
magma_tally4_dgetrs_nopiv_gpu(
    magma_tally4_trans_t trans, magma_tally4_int_t n, magma_tally4_int_t nrhs,
    magma_tally4Double_ptr dA, magma_tally4_int_t ldda,
    magma_tally4Double_ptr dB, magma_tally4_int_t lddb,
    magma_tally4_int_t *info)
{
    double c_one = MAGMA_tally4_D_ONE;
    int notran = (trans == Magma_tally4NoTrans);

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

    if (notran) {
        /* Solve A * X = B. */
        if ( nrhs == 1) {
            magma_tally4_dtrsv(Magma_tally4Lower, Magma_tally4NoTrans, Magma_tally4Unit,    n, dA, ldda, dB, 1 );
            magma_tally4_dtrsv(Magma_tally4Upper, Magma_tally4NoTrans, Magma_tally4NonUnit, n, dA, ldda, dB, 1 );
        } else {
            magma_tally4_dtrsm(Magma_tally4Left, Magma_tally4Lower, Magma_tally4NoTrans, Magma_tally4Unit,    n, nrhs, c_one, dA, ldda, dB, lddb );
            magma_tally4_dtrsm(Magma_tally4Left, Magma_tally4Upper, Magma_tally4NoTrans, Magma_tally4NonUnit, n, nrhs, c_one, dA, ldda, dB, lddb );
        }
    } else {
        /* Solve A**T * X = B  or  A**H * X = B. */
        if ( nrhs == 1) {
            magma_tally4_dtrsv(Magma_tally4Upper, trans, Magma_tally4NonUnit, n, dA, ldda, dB, 1 );
            magma_tally4_dtrsv(Magma_tally4Lower, trans, Magma_tally4Unit,    n, dA, ldda, dB, 1 );
        } else {
            magma_tally4_dtrsm(Magma_tally4Left, Magma_tally4Upper, trans, Magma_tally4NonUnit, n, nrhs, c_one, dA, ldda, dB, lddb );
            magma_tally4_dtrsm(Magma_tally4Left, Magma_tally4Lower, trans, Magma_tally4Unit,    n, nrhs, c_one, dA, ldda, dB, lddb );
        }
    }

    return *info;
}
