/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from zcgetrs_gpu.cpp mixed zc -> ds, Fri Jan 30 19:00:13 2015

*/
#include "common_magma_tally2.h"

/**
    Purpose
    -------
    DSGETRS solves a system of linear equations
       A * X = B,  A**T * X = B,  or  A**H * X = B
    with a general N-by-N matrix A using the LU factorization computed
    by MAGMA_tally2_SGETRF_GPU. B and X are in DOUBLE PRECISION, and A is in SINGLE PRECISION.
    This routine is used in the mixed precision iterative solver
    magma_tally2_dsgesv.

    Arguments
    ---------
    @param[in]
    trans   magma_tally2_trans_t
            Specifies the form of the system of equations:
      -     = Magma_tally2NoTrans:    A * X = B     (No transpose)
      -     = Magma_tally2Trans:      A**T * X = B  (Transpose)
      -     = Magma_tally2ConjTrans:  A**H * X = B  (Conjugate transpose)

    @param[in]
    n       INTEGER
            The order of the matrix A.  N >= 0.

    @param[in]
    nrhs    INTEGER
            The number of right hand sides, i.e., the number of columns
            of the matrix B.  NRHS >= 0.

    @param[in]
    dA      SINGLE PRECISION array on the GPU, dimension (LDDA,N)
            The factors L and U from the factorization A = P*L*U
            as computed by CGETRF_GPU.

    @param[in]
    ldda    INTEGER
            The leading dimension of the array dA.  LDDA >= max(1,N).

    @param[in]
    dipiv   INTEGER array on the GPU, dimension (N)
            The pivot indices; for 1 <= i <= N, after permuting, row i of the
            matrix was moved to row dIPIV(i).
            Note this is different than IPIV from DGETRF, where interchanges
            are applied one-after-another.

    @param[in]
    dB      DOUBLE PRECISION array on the GPU, dimension (LDDB,NRHS)
            On entry, the right hand side matrix B.

    @param[in]
    lddb    INTEGER
            The leading dimension of the arrays X and B.  LDDB >= max(1,N).

    @param[out]
    dX      DOUBLE PRECISION array on the GPU, dimension (LDDX, NRHS)
            On exit, the solution matrix dX.

    @param[in]
    lddx    INTEGER
            The leading dimension of the array dX, LDDX >= max(1,N).

    @param
    dSX     (workspace) SINGLE PRECISION array on the GPU used as workspace,
            dimension (N, NRHS)

    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument had an illegal value

    @ingroup magma_tally2_dgesv_comp
    ********************************************************************/
extern "C" magma_tally2_int_t
magma_tally2_dsgetrs_gpu(
    magma_tally2_trans_t trans, magma_tally2_int_t n, magma_tally2_int_t nrhs,
    magma_tally2Float_ptr  dA, magma_tally2_int_t ldda,
    magma_tally2Int_ptr        dipiv,
    magma_tally2Double_ptr dB, magma_tally2_int_t lddb,
    magma_tally2Double_ptr dX, magma_tally2_int_t lddx,
    magma_tally2Float_ptr dSX,
    magma_tally2_int_t *info)
{
    float c_one = MAGMA_tally2_S_ONE;
    int notran = (trans == Magma_tally2NoTrans);
    magma_tally2_int_t inc;

    *info = 0;
    if ( (! notran) &&
         (trans != Magma_tally2Trans) &&
         (trans != Magma_tally2ConjTrans) ) {
        *info = -1;
    } else if (n < 0) {
        *info = -2;
    } else if (nrhs < 0) {
        *info = -3;
    } else if (ldda < n) {
        *info = -5;
    } else if (lddb < n) {
        *info = -8;
    } else if (lddx < n) {
        *info = -10;
    } else if (lddx != lddb) { /* TODO: remove it when dslaswp will have the correct interface */
        *info = -10;
    }
    if (*info != 0) {
        magma_tally2_xerbla( __func__, -(*info) );
        return *info;
    }

    /* Quick return if possible */
    if (n == 0 || nrhs == 0) {
        return *info;
    }
    
    if (notran) {
        inc = 1;
        
        /* Get X by row applying interchanges to B and cast to single */
        /*
         * TODO: clean dslaswp interface to have interface closer to zlaswp
         */
        //magma_tally2blas_dslaswp(nrhs, dB, lddb, dSX, lddbx, 1, n, dipiv);
        magma_tally2blas_dslaswp(nrhs, dB, lddb, dSX, n, dipiv, inc);
        
        /* Solve L*X = B, overwriting B with SX. */
        magma_tally2_strsm( Magma_tally2Left, Magma_tally2Lower, Magma_tally2NoTrans, Magma_tally2Unit,
                     n, nrhs, c_one, dA, ldda, dSX, n);
        
        /* Solve U*X = B, overwriting B with X. */
        magma_tally2_strsm( Magma_tally2Left, Magma_tally2Upper, Magma_tally2NoTrans, Magma_tally2NonUnit,
                     n, nrhs, c_one, dA, ldda, dSX, n);
        
        magma_tally2blas_slag2d( n, nrhs, dSX, n, dX, lddx, info );
    }
    else {
        inc = -1;
        
        /* Cast the DOUBLE PRECISION RHS to SINGLE PRECISION */
        magma_tally2blas_dlag2s( n, nrhs, dB, lddb, dSX, n, info );
        
        /* Solve A**T * X = B, or A**H * X = B */
        magma_tally2_strsm( Magma_tally2Left, Magma_tally2Upper, trans, Magma_tally2NonUnit,
                     n, nrhs, c_one, dA, ldda, dSX, n );
        magma_tally2_strsm( Magma_tally2Left, Magma_tally2Lower, trans, Magma_tally2Unit,
                     n, nrhs, c_one, dA, ldda, dSX, n );
        
        magma_tally2blas_dslaswp( nrhs, dX, lddx, dSX, n, dipiv, inc );
    }

    return *info;
} /* magma_tally2_dsgetrs */