/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015
       @author Adrien REMY

       @generated from zgetrs_nopiv_gpu.cpp normal z -> s, Fri Jan 30 19:00:14 2015

*/
#include "common_magma_tally3.h"

/**
    Purpose
    -------
    Solves a system of linear equations
      A * X = B,  A**T * X = B,  or  A**H * X = B
    with a general N-by-N matrix A using the LU factorization computed by SGETRF_NOPIV_GPU.

    Arguments
    ---------
    @param[in]
    trans   magma_tally3_trans_t
            Specifies the form of the system of equations:
      -     = Magma_tally3NoTrans:    A    * X = B  (No transpose)
      -     = Magma_tally3Trans:      A**T * X = B  (Transpose)
      -     = Magma_tally3ConjTrans:  A**H * X = B  (Conjugate transpose)

    @param[in]
    n       INTEGER
            The order of the matrix A.  N >= 0.

    @param[in]
    nrhs    INTEGER
            The number of right hand sides, i.e., the number of columns
            of the matrix B.  NRHS >= 0.

    @param[in]
    dA      REAL array on the GPU, dimension (LDA,N)
            The factors L and U from the factorization A = P*L*U as computed
            by SGETRF_GPU.

    @param[in]
    ldda    INTEGER
            The leading dimension of the array A.  LDA >= max(1,N).

    param[in,out]
    dB      REAL array on the GPU, dimension (LDB,NRHS)
            On entry, the right hand side matrix B.
            On exit, the solution matrix X.

    @param[in]
    lddb    INTEGER
            The leading dimension of the array B.  LDB >= max(1,N).

    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument had an illegal value

    @ingroup magma_tally3_sgesv_comp
    ********************************************************************/
extern "C" magma_tally3_int_t
magma_tally3_sgetrs_nopiv_gpu(
    magma_tally3_trans_t trans, magma_tally3_int_t n, magma_tally3_int_t nrhs,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Float_ptr dB, magma_tally3_int_t lddb,
    magma_tally3_int_t *info)
{
    float c_one = MAGMA_tally3_S_ONE;
    int notran = (trans == Magma_tally3NoTrans);

    *info = 0;
    if ( (! notran) &&
         (trans != Magma_tally3Trans) &&
         (trans != Magma_tally3ConjTrans) ) {
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
        magma_tally3_xerbla( __func__, -(*info) );
        return *info;
    }

    /* Quick return if possible */
    if (n == 0 || nrhs == 0) {
        return *info;
    }

    if (notran) {
        /* Solve A * X = B. */
        if ( nrhs == 1) {
            magma_tally3_strsv(Magma_tally3Lower, Magma_tally3NoTrans, Magma_tally3Unit,    n, dA, ldda, dB, 1 );
            magma_tally3_strsv(Magma_tally3Upper, Magma_tally3NoTrans, Magma_tally3NonUnit, n, dA, ldda, dB, 1 );
        } else {
            magma_tally3_strsm(Magma_tally3Left, Magma_tally3Lower, Magma_tally3NoTrans, Magma_tally3Unit,    n, nrhs, c_one, dA, ldda, dB, lddb );
            magma_tally3_strsm(Magma_tally3Left, Magma_tally3Upper, Magma_tally3NoTrans, Magma_tally3NonUnit, n, nrhs, c_one, dA, ldda, dB, lddb );
        }
    } else {
        /* Solve A**T * X = B  or  A**H * X = B. */
        if ( nrhs == 1) {
            magma_tally3_strsv(Magma_tally3Upper, trans, Magma_tally3NonUnit, n, dA, ldda, dB, 1 );
            magma_tally3_strsv(Magma_tally3Lower, trans, Magma_tally3Unit,    n, dA, ldda, dB, 1 );
        } else {
            magma_tally3_strsm(Magma_tally3Left, Magma_tally3Upper, trans, Magma_tally3NonUnit, n, nrhs, c_one, dA, ldda, dB, lddb );
            magma_tally3_strsm(Magma_tally3Left, Magma_tally3Lower, trans, Magma_tally3Unit,    n, nrhs, c_one, dA, ldda, dB, lddb );
        }
    }

    return *info;
}
