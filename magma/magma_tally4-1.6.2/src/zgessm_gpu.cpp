/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Hatem Ltaief
       @author Mathieu Faverge

       @precisions normal z -> c d s

*/
#include "common_magma_tally4.h"

/**
    Purpose
    -------
    ZGESSM applies the factors L computed by ZGETRF_INCPIV to
    a complex M-by-N tile A.
    
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
            ZGETRF_INCPIV.

    @param[in]
    dL1     COMPLEX_16 array, dimension(LDDL1, N)
            The IB-by-K matrix in which is stored L^(-1) as returned by GETRF_INCPIV

    @param[in]
    lddl1   INTEGER
            The leading dimension of the array L1.  LDDL1 >= max(1,2*IB).

    @param[in]
    dL      COMPLEX_16 array, dimension(LDDL, N)
            The M-by-K lower triangular tile on the gpu.

    @param[in]
    lddl    INTEGER
            The leading dimension of the array L.  LDDL >= max(1,M).

    @param[in,out]
    dA      COMPLEX_16 array, dimension (LDDA, N)
            On entry, the M-by-N tile A on the gpu.
            On exit, updated by the application of L on the gpu.

    @param[in]
    ldda    INTEGER
            The leading dimension of the array A.  LDDA >= max(1,M).

    @ingroup magma_tally4_zgesv_tile
    ********************************************************************/
extern "C" magma_tally4_int_t
magma_tally4_zgessm_gpu(
    magma_tally4_order_t order, magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k, magma_tally4_int_t ib,
    magma_tally4_int_t *ipiv,
    magma_tally4DoubleComplex_ptr dL1, magma_tally4_int_t lddl1,
    magma_tally4DoubleComplex_ptr dL,  magma_tally4_int_t lddl,
    magma_tally4DoubleComplex_ptr dA,  magma_tally4_int_t ldda,
    magma_tally4_int_t *info)
{
#define AT(i,j) (dAT + (i)*ldda + (j)      )
#define L(i,j)  (dL  + (i)      + (j)*lddl )
#define dL1(j)  (dL1            + (j)*lddl1)

    magma_tally4DoubleComplex c_one     = MAGMA_tally4_Z_ONE;
    magma_tally4DoubleComplex c_neg_one = MAGMA_tally4_Z_NEG_ONE;

    int i, sb;
    magma_tally4DoubleComplex_ptr dAT;

    /* Check arguments */
    *info = 0;
    if (m < 0)
        *info = -1;
    else if (n < 0)
        *info = -2;
    else if (ldda < max(1,m))
        *info = -4;

    if (*info != 0) {
        magma_tally4_xerbla( __func__, -(*info) );
        return *info;
    }

    /* Quick return if possible */
    if (m == 0 || n == 0)
        return *info;

    if ( order == Magma_tally4ColMajor ) {
        magma_tally4blas_zgetmo_in( dA, dAT, ldda, m, n );
    } else {
        dAT = dA;
    }

    for (i = 0; i < k; i += ib) {
        sb = min(ib, k-i);

        magma_tally4blas_zlaswp( n, dAT, ldda, i+1, i+sb, ipiv, 1 );

#ifndef WITHOUTTRTRI
        magma_tally4_ztrmm( Magma_tally4Right, Magma_tally4Lower, Magma_tally4Trans, Magma_tally4Unit,
                     n, sb,
                     c_one, dL1(i),   lddl1,
                            AT(i, 0), ldda);
#else
        magma_tally4_ztrsm( Magma_tally4Right, Magma_tally4Lower, Magma_tally4Trans, Magma_tally4Unit,
                     n, sb,
                     c_one, L( i, i), lddl,
                            AT(i, 0), ldda);
#endif

        if ( (i+sb) < m) {
            magma_tally4_zgemm( Magma_tally4NoTrans, Magma_tally4Trans,
                         n, m-(i+sb), sb,
                         c_neg_one, AT(i,    0), ldda,
                                    L( i+sb, i), lddl,
                         c_one,     AT(i+sb, 0), ldda );
        }
    }

    if ( order == Magma_tally4ColMajor ) {
        magma_tally4blas_zgetmo_in( dA, dAT, ldda, m, n );
    }

    return *info;
} /* magma_tally4_zgessm_gpu */