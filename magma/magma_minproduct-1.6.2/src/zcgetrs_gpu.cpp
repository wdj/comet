/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions mixed zc -> ds

*/
#include "common_magma_minproduct.h"

/**
    Purpose
    -------
    ZCGETRS solves a system of linear equations
       A * X = B,  A**T * X = B,  or  A**H * X = B
    with a general N-by-N matrix A using the LU factorization computed
    by MAGMA_minproduct_CGETRF_GPU. B and X are in COMPLEX_16, and A is in COMPLEX.
    This routine is used in the mixed precision iterative solver
    magma_minproduct_zcgesv.

    Arguments
    ---------
    @param[in]
    trans   magma_minproduct_trans_t
            Specifies the form of the system of equations:
      -     = Magma_minproductNoTrans:    A * X = B     (No transpose)
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
    dA      COMPLEX array on the GPU, dimension (LDDA,N)
            The factors L and U from the factorization A = P*L*U
            as computed by CGETRF_GPU.

    @param[in]
    ldda    INTEGER
            The leading dimension of the array dA.  LDDA >= max(1,N).

    @param[in]
    dipiv   INTEGER array on the GPU, dimension (N)
            The pivot indices; for 1 <= i <= N, after permuting, row i of the
            matrix was moved to row dIPIV(i).
            Note this is different than IPIV from ZGETRF, where interchanges
            are applied one-after-another.

    @param[in]
    dB      COMPLEX_16 array on the GPU, dimension (LDDB,NRHS)
            On entry, the right hand side matrix B.

    @param[in]
    lddb    INTEGER
            The leading dimension of the arrays X and B.  LDDB >= max(1,N).

    @param[out]
    dX      COMPLEX_16 array on the GPU, dimension (LDDX, NRHS)
            On exit, the solution matrix dX.

    @param[in]
    lddx    INTEGER
            The leading dimension of the array dX, LDDX >= max(1,N).

    @param
    dSX     (workspace) COMPLEX array on the GPU used as workspace,
            dimension (N, NRHS)

    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument had an illegal value

    @ingroup magma_minproduct_zgesv_comp
    ********************************************************************/
extern "C" magma_minproduct_int_t
magma_minproduct_zcgetrs_gpu(
    magma_minproduct_trans_t trans, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
    magma_minproductFloatComplex_ptr  dA, magma_minproduct_int_t ldda,
    magma_minproductInt_ptr        dipiv,
    magma_minproductDoubleComplex_ptr dB, magma_minproduct_int_t lddb,
    magma_minproductDoubleComplex_ptr dX, magma_minproduct_int_t lddx,
    magma_minproductFloatComplex_ptr dSX,
    magma_minproduct_int_t *info)
{
    magma_minproductFloatComplex c_one = MAGMA_minproduct_C_ONE;
    int notran = (trans == Magma_minproductNoTrans);
    magma_minproduct_int_t inc;

    *info = 0;
    if ( (! notran) &&
         (trans != Magma_minproductTrans) &&
         (trans != Magma_minproductConjTrans) ) {
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
    } else if (lddx != lddb) { /* TODO: remove it when zclaswp will have the correct interface */
        *info = -10;
    }
    if (*info != 0) {
        magma_minproduct_xerbla( __func__, -(*info) );
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
         * TODO: clean zclaswp interface to have interface closer to zlaswp
         */
        //magma_minproductblas_zclaswp(nrhs, dB, lddb, dSX, lddbx, 1, n, dipiv);
        magma_minproductblas_zclaswp(nrhs, dB, lddb, dSX, n, dipiv, inc);
        
        /* Solve L*X = B, overwriting B with SX. */
        magma_minproduct_ctrsm( Magma_minproductLeft, Magma_minproductLower, Magma_minproductNoTrans, Magma_minproductUnit,
                     n, nrhs, c_one, dA, ldda, dSX, n);
        
        /* Solve U*X = B, overwriting B with X. */
        magma_minproduct_ctrsm( Magma_minproductLeft, Magma_minproductUpper, Magma_minproductNoTrans, Magma_minproductNonUnit,
                     n, nrhs, c_one, dA, ldda, dSX, n);
        
        magma_minproductblas_clag2z( n, nrhs, dSX, n, dX, lddx, info );
    }
    else {
        inc = -1;
        
        /* Cast the COMPLEX_16 RHS to COMPLEX */
        magma_minproductblas_zlag2c( n, nrhs, dB, lddb, dSX, n, info );
        
        /* Solve A**T * X = B, or A**H * X = B */
        magma_minproduct_ctrsm( Magma_minproductLeft, Magma_minproductUpper, trans, Magma_minproductNonUnit,
                     n, nrhs, c_one, dA, ldda, dSX, n );
        magma_minproduct_ctrsm( Magma_minproductLeft, Magma_minproductLower, trans, Magma_minproductUnit,
                     n, nrhs, c_one, dA, ldda, dSX, n );
        
        magma_minproductblas_zclaswp( nrhs, dX, lddx, dSX, n, dipiv, inc );
    }

    return *info;
} /* magma_minproduct_zcgetrs */
