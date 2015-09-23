/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Hatem Ltaief
       @author Mathieu Faverge

       @generated from zgessm_gpu.cpp normal z -> s, Fri Jan 30 19:00:14 2015

*/
#include "common_magma_minproduct.h"

/**
    Purpose
    -------
    SGESSM applies the factors L computed by SGETRF_INCPIV to
    a real M-by-N tile A.
    
    Arguments
    ---------
    @param[in]
    m       INTEGER
            The number of rows of the matrix A.  M >= 0.

    @param[in]
    n       INTEGER
            The number of columns of the matrix A.  N >= 0.

    @param[in]
    k       INTEGER
            The number of columns of the matrix L.  K >= 0.

    @param[in]
    ib      INTEGER
            The inner-blocking size.  IB >= 0.

    @param[in]
    ipiv    INTEGER array on the cpu.
            The pivot indices array of size K as returned by
            SGETRF_INCPIV.

    @param[in]
    dL1     REAL array, dimension(LDDL1, N)
            The IB-by-K matrix in which is stored L^(-1) as returned by GETRF_INCPIV

    @param[in]
    lddl1   INTEGER
            The leading dimension of the array L1.  LDDL1 >= max(1,2*IB).

    @param[in]
    dL      REAL array, dimension(LDDL, N)
            The M-by-K lower triangular tile on the gpu.

    @param[in]
    lddl    INTEGER
            The leading dimension of the array L.  LDDL >= max(1,M).

    @param[in,out]
    dA      REAL array, dimension (LDDA, N)
            On entry, the M-by-N tile A on the gpu.
            On exit, updated by the application of L on the gpu.

    @param[in]
    ldda    INTEGER
            The leading dimension of the array A.  LDDA >= max(1,M).

    @ingroup magma_minproduct_sgesv_tile
    ********************************************************************/
extern "C" magma_minproduct_int_t
magma_minproduct_sgessm_gpu(
    magma_minproduct_order_t order, magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k, magma_minproduct_int_t ib,
    magma_minproduct_int_t *ipiv,
    magma_minproductFloat_ptr dL1, magma_minproduct_int_t lddl1,
    magma_minproductFloat_ptr dL,  magma_minproduct_int_t lddl,
    magma_minproductFloat_ptr dA,  magma_minproduct_int_t ldda,
    magma_minproduct_int_t *info)
{
#define AT(i,j) (dAT + (i)*ldda + (j)      )
#define L(i,j)  (dL  + (i)      + (j)*lddl )
#define dL1(j)  (dL1            + (j)*lddl1)

    float c_one     = MAGMA_minproduct_S_ONE;
    float c_neg_one = MAGMA_minproduct_S_NEG_ONE;

    int i, sb;
    magma_minproductFloat_ptr dAT;

    /* Check arguments */
    *info = 0;
    if (m < 0)
        *info = -1;
    else if (n < 0)
        *info = -2;
    else if (ldda < max(1,m))
        *info = -4;

    if (*info != 0) {
        magma_minproduct_xerbla( __func__, -(*info) );
        return *info;
    }

    /* Quick return if possible */
    if (m == 0 || n == 0)
        return *info;

    if ( order == Magma_minproductColMajor ) {
        magma_minproductblas_sgetmo_in( dA, dAT, ldda, m, n );
    } else {
        dAT = dA;
    }

    for (i = 0; i < k; i += ib) {
        sb = min(ib, k-i);

        magma_minproductblas_slaswp( n, dAT, ldda, i+1, i+sb, ipiv, 1 );

#ifndef WITHOUTTRTRI
        magma_minproduct_strmm( Magma_minproductRight, Magma_minproductLower, Magma_minproductTrans, Magma_minproductUnit,
                     n, sb,
                     c_one, dL1(i),   lddl1,
                            AT(i, 0), ldda);
#else
        magma_minproduct_strsm( Magma_minproductRight, Magma_minproductLower, Magma_minproductTrans, Magma_minproductUnit,
                     n, sb,
                     c_one, L( i, i), lddl,
                            AT(i, 0), ldda);
#endif

        if ( (i+sb) < m) {
            magma_minproduct_sgemm( Magma_minproductNoTrans, Magma_minproductTrans,
                         n, m-(i+sb), sb,
                         c_neg_one, AT(i,    0), ldda,
                                    L( i+sb, i), lddl,
                         c_one,     AT(i+sb, 0), ldda );
        }
    }

    if ( order == Magma_minproductColMajor ) {
        magma_minproductblas_sgetmo_in( dA, dAT, ldda, m, n );
    }

    return *info;
} /* magma_minproduct_sgessm_gpu */
