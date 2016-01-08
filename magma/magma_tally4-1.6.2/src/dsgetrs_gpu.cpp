/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from zcgetrs_gpu.cpp mixed zc -> ds, Fri Jan 30 19:00:13 2015

*/
#include "common_magma_tally4.h"

/**
    Purpose
    -------
    DSGETRS solves a system of linear equations
       A * X = B,  A**T * X = B,  or  A**H * X = B
    with a general N-by-N matrix A using the LU factorization computed
    by MAGMA_tally4_SGETRF_GPU. B and X are in DOUBLE PRECISION, and A is in SINGLE PRECISION.
    This routine is used in the mixed precision iterative solver
    magma_tally4_dsgesv.

    Arguments
    ---------
    @param[in]
    trans   magma_tally4_trans_t
            Specifies the form of the system of equations:
      -     = Magma_tally4NoTrans:    A * X = B     (No transpose)
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

    @ingroup magma_tally4_dgesv_comp
    ********************************************************************/
extern "C" magma_tally4_int_t
magma_tally4_dsgetrs_gpu(
    magma_tally4_trans_t trans, magma_tally4_int_t n, magma_tally4_int_t nrhs,
    magma_tally4Float_ptr  dA, magma_tally4_int_t ldda,
    magma_tally4Int_ptr        dipiv,
    magma_tally4Double_ptr dB, magma_tally4_int_t lddb,
    magma_tally4Double_ptr dX, magma_tally4_int_t lddx,
    magma_tally4Float_ptr dSX,
    magma_tally4_int_t *info)
{
    float c_one = MAGMA_tally4_S_ONE;
    int notran = (trans == Magma_tally4NoTrans);
    magma_tally4_int_t inc;

    *info = 0;
    if ( (! notran) &&
         (trans != Magma_tally4Trans) &&
         (trans != Magma_tally4ConjTrans) ) {
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
        magma_tally4_xerbla( __func__, -(*info) );
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
        //magma_tally4blas_dslaswp(nrhs, dB, lddb, dSX, lddbx, 1, n, dipiv);
        magma_tally4blas_dslaswp(nrhs, dB, lddb, dSX, n, dipiv, inc);
        
        /* Solve L*X = B, overwriting B with SX. */
        magma_tally4_strsm( Magma_tally4Left, Magma_tally4Lower, Magma_tally4NoTrans, Magma_tally4Unit,
                     n, nrhs, c_one, dA, ldda, dSX, n);
        
        /* Solve U*X = B, overwriting B with X. */
        magma_tally4_strsm( Magma_tally4Left, Magma_tally4Upper, Magma_tally4NoTrans, Magma_tally4NonUnit,
                     n, nrhs, c_one, dA, ldda, dSX, n);
        
        magma_tally4blas_slag2d( n, nrhs, dSX, n, dX, lddx, info );
    }
    else {
        inc = -1;
        
        /* Cast the DOUBLE PRECISION RHS to SINGLE PRECISION */
        magma_tally4blas_dlag2s( n, nrhs, dB, lddb, dSX, n, info );
        
        /* Solve A**T * X = B, or A**H * X = B */
        magma_tally4_strsm( Magma_tally4Left, Magma_tally4Upper, trans, Magma_tally4NonUnit,
                     n, nrhs, c_one, dA, ldda, dSX, n );
        magma_tally4_strsm( Magma_tally4Left, Magma_tally4Lower, trans, Magma_tally4Unit,
                     n, nrhs, c_one, dA, ldda, dSX, n );
        
        magma_tally4blas_dslaswp( nrhs, dX, lddx, dSX, n, dipiv, inc );
    }

    return *info;
} /* magma_tally4_dsgetrs */
