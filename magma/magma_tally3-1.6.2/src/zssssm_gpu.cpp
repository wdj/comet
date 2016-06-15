/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Hatem Ltaief
       @author Mathieu Faverge

       @precisions normal z -> c d s

*/
#include "common_magma_tally3.h"

/**
    Purpose
    -------
    ZSSSSM applies the LU factorization update from a complex
    matrix formed by a lower triangular IB-by-K tile L1 on top of a
    M2-by-K tile L2 to a second complex matrix formed by a M1-by-N1
    tile A1 on top of a M2-by-N2 tile A2 (N1 == N2).

    This is the right-looking Level 2.5 BLAS version of the algorithm.

    Arguments
    ---------
    @param[in]
    m1      INTEGER
            The number of rows of the matrix A1.  M1 >= 0.

    @param[in]
    n1      INTEGER
            The number of columns of the matrix A1.  N1 >= 0.

    @param[in]
    m2      INTEGER
            The number of rows of the matrix A2.  M2 >= 0.

    @param[in]
    n2      INTEGER
            The number of columns of the matrix A2.  N2 >= 0.

    @param[in]
    k       INTEGER
            The number of columns of the matrix L1 and L2.  K >= 0.

    @param[in]
    ib      INTEGER
            The inner-blocking size.  IB >= 0.

    @param[in,out]
    dA1     COMPLEX_16 array, dimension(LDDA1, N), on gpu.
            On entry, the M1-by-N1 tile dA1.
            On exit, dA1 is updated by the application of dL (dL1 dL2).

    @param[in]
    ldda1   INTEGER
            The leading dimension of the array dA1.  LDDA1 >= max(1,M1).

    @param[in,out]
    dA2     COMPLEX_16 array, dimension(LDDA2, N), on gpu.
            On entry, the M2-by-N2 tile dA2.
            On exit, dA2 is updated by the application of dL (dL1 dL2).

    @param[in]
    ldda2   INTEGER
            The leading dimension of the array dA2.  LDDA2 >= max(1,M2).

    @param[in]
    dL1     COMPLEX_16 array, dimension(LDDL1, K), on gpu.
            The inverse of the IB-by-K lower triangular tile as returned by
            ZTSTRF.

    @param[in]
    lddl1   INTEGER
            The leading dimension of the array L1.  LDDL1 >= max(1,2*IB).

    @param[in]
    dL2     COMPLEX_16 array, dimension(LDDL2, K)
            The M2-by-K tile as returned by ZTSTRF.

    @param[in]
    lddl2   INTEGER
            The leading dimension of the array L2.  LDDL2 >= max(1,M2).

    @param[in]
    ipiv    INTEGER array on the cpu.
            The pivot indices array of size K as returned by ZTSTRF

    @ingroup magma_tally3_zgesv_tile
    ********************************************************************/
extern "C" magma_tally3_int_t
magma_tally3_zssssm_gpu(
    magma_tally3_order_t order, magma_tally3_int_t m1, magma_tally3_int_t n1,
    magma_tally3_int_t m2, magma_tally3_int_t n2, magma_tally3_int_t k, magma_tally3_int_t ib,
    magma_tally3DoubleComplex_ptr dA1, magma_tally3_int_t ldda1,
    magma_tally3DoubleComplex_ptr dA2, magma_tally3_int_t ldda2,
    magma_tally3DoubleComplex_ptr dL1, magma_tally3_int_t lddl1,
    magma_tally3DoubleComplex_ptr dL2, magma_tally3_int_t lddl2,
    magma_tally3_int_t *ipiv,
    magma_tally3_int_t *info)
{
#define A1T(i,j) (dA1T + (i)*ldda1 + (j))
#define A2T(i,j) (dA2T + (i)*ldda2 + (j))
#define L1(i)    (dL1  + (i)*lddl1      )
#define L2(i,j)  (dL2  + (i)*lddl2i + (j)*lddl2j)

    magma_tally3DoubleComplex c_one     = MAGMA_tally3_Z_ONE;
    magma_tally3DoubleComplex c_neg_one = MAGMA_tally3_Z_NEG_ONE;

    int ip, ii, sb;
    magma_tally3DoubleComplex_ptr dA1T, dA2T;
    magma_tally3_trans_t transL;
    int lddl2i, lddl2j;
    
    MAGMA_tally3_UNUSED( ip );  // used only if NOSWAPBLK

    /* Check input arguments */
    *info = 0;
    if (m1 < 0) {
        *info = -1;
    }
    else if (n1 < 0) {
        *info = -2;
    }
    else if (m2 < 0) {
        *info = -3;
    }
    else if (n2 < 0) {
        *info = -4;
    }
    else if (k < 0) {
        *info = -5;
    }
    else if (ib < 0) {
        *info = -6;
    }
    else if (ldda1 < max(1,m1)) {
        *info = -8;
    }
    else if (ldda2 < max(1,m2)) {
        *info = -10;
    }
    else if (lddl1 < max(1,ib)) {
        *info = -12;
    }
    else if (lddl2 < max(1,m2)) {
        *info = -14;
    }

    if (*info != 0) {
        magma_tally3_xerbla( __func__, -(*info) );
        return *info;
    }

    /* Quick return */
    if ((m1 == 0) || (n1 == 0) || (m2 == 0) || (n2 == 0) || (k == 0) || (ib == 0))
        return *info;

    if ( order == Magma_tally3ColMajor ) {
        magma_tally3blas_zgetmo_in( dA1, dA1T, ldda1, m1, n1 );
        magma_tally3blas_zgetmo_in( dA2, dA2T, ldda2, m2, n2 );
        transL = Magma_tally3Trans;
        lddl2i = 1; lddl2j = lddl2;
    } else {
        dA1T = dA1;
        dA2T = dA2;
        transL = Magma_tally3NoTrans;
        lddl2i = lddl2; lddl2j = 1;
    }

    ip = 0;
    for( ii=0; ii < k; ii += ib ) {
        sb = min( k-ii, ib);

#ifndef NOSWAPBLK
        magma_tally3blas_zswapblk( Magma_tally3RowMajor, n1,
                            A1T(0, 0), ldda1,
                            A2T(0, 0), ldda2,
                            ii+1, ii+ib, ipiv, 1, m1 );
#else
        {
            int im;
            for (i=0; i < ib; i++) {
                im = ipiv[ip]-1;

                if (im != (ii+i)) {
                    im = im - m1;

                    assert( (im >= 0) && (im < m1) && (im < m2) );
                    magma_tally3blas_zswap( n1, A1T(ii+i, 0), 1, A2T(im, 0), 1 );
                }
                ip++;
            }
        }
#endif

#ifndef WITHOUTTRTRI
        /* Lower, Trans, because L1 is not transposed */
        magma_tally3_ztrmm( Magma_tally3Right, Magma_tally3Lower, Magma_tally3Trans, Magma_tally3Unit,
                     n1, sb,
                     c_one, L1( ii),    lddl1,
                            A1T(ii, 0), ldda1);
#else
        /* Lower, Trans, because L1 is not transposed */
        magma_tally3_ztrsm( Magma_tally3Right, Magma_tally3Lower, Magma_tally3Trans, Magma_tally3Unit,
                     n1, sb,
                     c_one, L1( ii),    lddl1,
                            A1T(ii, 0), ldda1);
#endif

        /* Second parameter is trans because L2 is not transposed */
        magma_tally3_zgemm( Magma_tally3NoTrans, transL,
                     n2, m2, sb,
                     c_neg_one, A1T(ii, 0), ldda1,
                                L2( 0, ii), lddl2,
                     c_one,     A2T(0, 0 ), ldda2 );
    }

    if ( order == Magma_tally3ColMajor ) {
        magma_tally3blas_zgetmo_out( dA1, dA1T, ldda1, m1, n1 );
        magma_tally3blas_zgetmo_out( dA2, dA2T, ldda2, m2, n2 );
    }
    return *info;
}