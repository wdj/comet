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
#ifdef MAGMA_tally4_WITH_PLASMA

#include <plasma.h>
#include <core_blas.h>
#include "common_magma_tally4.h"

/**
    Purpose
    -------
    ZGETRF_INCPIV computes an LU factorization of a general M-by-N tile A
    using partial pivoting with row interchanges.

    The factorization has the form

      A = P * L * U

    where P is a permutation matrix, L is lower triangular with unit
    diagonal elements (lower trapezoidal if m > n), and U is upper
    triangular (upper trapezoidal if m < n).

    This is the right-looking Level 2.5 BLAS version of the algorithm.

    Arguments
    ---------
    @param[in]
    m       INTEGER
            The number of rows of the matrix A.  M >= 0.

    @param[in]
    n       INTEGER
            The number of columns of the matrix A.  N >= 0.

    @param[in]
    ib      INTEGER
            The inner-blocking size.  IB >= 0.

    @param[in,out]
    hA      COMPLEX_16 array, dimension(LDHA, N), on cpu.
            On entry, only the M-by-IB first panel needs to be identical to dA(1..M, 1..IB).
            On exit, the content is incomplete. Shouldn't be used.

    @param[in]
    ldha    INTEGER
            The leading dimension of the array hA.  LDHA >= max(1,M).

    @param[in,out]
    dA      COMPLEX_16 array, dimension(LDDA, N), on gpu.
            On entry, the M-by-N tile to be factored.
            On exit, the factors L and U from the factorization
            A = P*L*U; the unit diagonal elements of L are not stored.

    @param[in]
    ldda    INTEGER
            The leading dimension of the array dA.  LDDA >= max(1,M).

    @param[out]
    hL      COMPLEX_16 array, dimension(LDHL, min(M,N)), on vpu.
            On exit, contains in the upper part the IB-by-K lower triangular tile,
            and in the lower part IB-by-min(M,N) the inverse of the top part.

    @param[in]
    ldhl    INTEGER
            The leading dimension of the array hL.  LDHL >= max(1,2*IB).

    @param[out]
    dL      COMPLEX_16 array, dimension(LDDL, K), on gpu.
            On exit, contains in the upper part the IB-by-min(M,N) lower triangular tile,
            and in the lower part IB-by-min(M,N) the inverse of the top part.

    @param[in]
    lddl    INTEGER
            The leading dimension of the array dL.  LDDL >= max(1,2*IB).

    @param[out]
    ipiv    INTEGER array, dimension min(M,N), on the cpu.
            The pivot indices array.

    @param[out]
    dWORK   COMPLEX_16 array, dimension(LDDWORK, 2*IB), on gpu.
            Workspace.

    @param[in]
    lddwork INTEGER
            The leading dimension of the array dWORK.  LDDWORK >= max(NB, 1).

    @param[out]
    info    INTEGER
            - PLASMA_SUCCESS successful exit
            - < 0 if INFO = -k, the k-th argument had an illegal value
            - > 0 if INFO = k, U(k,k) is exactly zero. The factorization
                has been completed, but the factor U is exactly
                singular, and division by zero will occur if it is used
                to solve a system of equations.

    @ingroup magma_tally4_zgesv_comp
    ********************************************************************/
extern "C" magma_tally4_int_t
magma_tally4_zgetrf_incpiv_gpu(
    magma_tally4_order_t order, magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t ib,
    magma_tally4DoubleComplex    *hA, magma_tally4_int_t ldha,
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex    *hL, magma_tally4_int_t ldhl,
    magma_tally4DoubleComplex_ptr dL, magma_tally4_int_t lddl,
    magma_tally4_int_t *ipiv,
    magma_tally4DoubleComplex_ptr dwork, magma_tally4_int_t lddwork,
    magma_tally4_int_t *info)
{
#define AT(i,j) (dAT + (i)*ib*ldda + (j)*ib)
#define hA(i,j) (hA  + (i)*ib + (j)*ib*ldha)
#define hL(j)   (hL  + (j)*ib*ldhl         )
#define hL2(j)  (hL2 + (j)*ib*ldhl         )
#define dL(j)   (dL  + (j)*ib*lddl         )
#define dL2(j)  (dL2 + (j)*ib*lddl         )

    magma_tally4DoubleComplex c_one     = MAGMA_tally4_Z_ONE;
    magma_tally4DoubleComplex c_neg_one = MAGMA_tally4_Z_NEG_ONE;

    magma_tally4_int_t iinfo;
    magma_tally4_int_t maxm, mindim;
    magma_tally4_int_t i, rows, cols, s, ii, sb;
    magma_tally4DoubleComplex_ptr dAT;
#ifndef WITHOUTTRTRI
    magma_tally4DoubleComplex_ptr dL2 = dL + ib;
    magma_tally4DoubleComplex *hL2 = hL + ib;
#endif

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

    /* Function Body */
    mindim = min(m, n);
    s      = mindim / ib;

    if ( ib >= mindim ) {
        /* Use CPU code. */
        lapackf77_zgetrf(&m, &n, hA, &ldha, ipiv, info);

#ifndef WITHOUTTRTRI
        CORE_zlacpy(PlasmaUpperLower, mindim, mindim,
                    (PLASMA_Complex64_t*)hA, ldha,
                    (PLASMA_Complex64_t*)hL2, ldhl );

        CORE_ztrtri( PlasmaLower, PlasmaUnit, mindim,
                     (PLASMA_Complex64_t*)hL2, ldhl, info );
        if (*info != 0 ) {
            fprintf(stderr, "ERROR, trtri returned with info = %d\n", *info);
        }

        magma_tally4_zsetmatrix( mindim, mindim, hL2, ldhl, dL2, lddl );
#endif

        if ( order == Magma_tally4RowMajor ) {
            magma_tally4_zsetmatrix( m, n, hA, ldha, dwork, lddwork );
            magma_tally4blas_ztranspose( m, n, dwork, lddwork, dA, ldda );
        } else {
            magma_tally4_zsetmatrix( m, n, hA, ldha, dA, ldda );
        }
    }
    else {
        /* Use hybrid blocked code. */
        maxm = ((m + 31)/32)*32;

        if ( order == Magma_tally4ColMajor ) {
            magma_tally4blas_zgetmo_in( dA, dAT, ldda, m, n );
        } else {
            dAT = dA;
        }

        for( i=0; i < s; i++ ) {
            ii = i * ib;
            sb = min(ib, mindim-ii);
            cols = maxm - ii;

            if ( i > 0 ) {
                // download i-th panel
                magma_tally4blas_ztranspose( sb, m, AT(0,i), ldda, dwork, maxm );
                magma_tally4_zgetmatrix( m, sb, dwork, maxm, hA(0, i), ldha );

                // make sure that gpu queue is empty
                //magma_tally4_device_sync();
#ifndef WITHOUTTRTRI
                magma_tally4_ztrmm( Magma_tally4Right, Magma_tally4Lower, Magma_tally4Trans, Magma_tally4Unit,
                             n - (ii+sb), ib,
                             c_one, dL2(i-1),    lddl,
                                    AT(i-1,i+1), ldda );
#else
                magma_tally4_ztrsm( Magma_tally4Right, Magma_tally4Upper, Magma_tally4NoTrans, Magma_tally4Unit,
                             n - (ii+sb), ib,
                             c_one, AT(i-1,i-1), ldda,
                                    AT(i-1,i+1), ldda );
#endif
                magma_tally4_zgemm( Magma_tally4NoTrans, Magma_tally4NoTrans,
                             n-(ii+sb), m-ii, ib,
                             c_neg_one, AT(i-1,i+1), ldda,
                                        AT(i,  i-1), ldda,
                             c_one,     AT(i,  i+1), ldda );
            }

            // do the cpu part
            rows = m - ii;
            lapackf77_zgetrf( &rows, &sb, hA(i, i), &ldha, ipiv+ii, &iinfo);
            if ( (*info == 0) && (iinfo > 0) )
                *info = iinfo + ii;

            {
                int j;
                int fin = ii + sb;
                for (j=ii; j < fin; j++) {
                    ipiv[j] = ii + ipiv[j];
                }
            }
            magma_tally4blas_zlaswp( n-ii, AT(0, i), ldda, ii+1, ii+sb, ipiv, 1 );

#ifndef WITHOUTTRTRI
            CORE_zlacpy(PlasmaLower, sb, sb,
                        (PLASMA_Complex64_t*)hA(i, i), ldha,
                        (PLASMA_Complex64_t*)hL2(i), ldhl );

            CORE_ztrtri( PlasmaLower, PlasmaUnit, sb,
                         (PLASMA_Complex64_t*)hL2(i), ldhl, info );
            if (*info != 0 ) {
                fprintf(stderr, "ERROR, trtri returned with info = %d\n", *info);
            }
            magma_tally4_zsetmatrix( sb, sb, hL2(i), ldhl, dL2(i), lddl );
#endif
            // upload i-th panel
            magma_tally4_zsetmatrix( rows, sb, hA(i, i), ldha, dwork, cols );
            magma_tally4blas_ztranspose( rows, sb, dwork, cols, AT(i,i), ldda );

            // do the small non-parallel computations
            if ( s > (i+1) ) {
#ifndef WITHOUTTRTRI
                magma_tally4_ztrmm( Magma_tally4Right, Magma_tally4Lower, Magma_tally4Trans, Magma_tally4Unit,
                             sb, sb,
                             c_one, dL2(i),     lddl,
                                    AT(i, i+1), ldda);
#else
                magma_tally4_ztrsm( Magma_tally4Right, Magma_tally4Upper, Magma_tally4NoTrans, Magma_tally4Unit,
                             sb, sb,
                             c_one, AT(i, i  ), ldda,
                                    AT(i, i+1), ldda);
#endif
                magma_tally4_zgemm( Magma_tally4NoTrans, Magma_tally4NoTrans,
                             sb, m-(ii+sb), sb,
                             c_neg_one, AT(i,   i+1), ldda,
                                        AT(i+1, i  ), ldda,
                             c_one,     AT(i+1, i+1), ldda );
            }
            else {
                /* Update of the last panel */
#ifndef WITHOUTTRTRI
                magma_tally4_ztrmm( Magma_tally4Right, Magma_tally4Lower, Magma_tally4Trans, Magma_tally4Unit,
                             n-mindim, sb,
                             c_one, dL2(i),     lddl,
                                    AT(i, i+1), ldda);
#else
                magma_tally4_ztrsm( Magma_tally4Right, Magma_tally4Upper, Magma_tally4NoTrans, Magma_tally4Unit,
                             n-mindim, sb,
                             c_one, AT(i, i  ), ldda,
                                    AT(i, i+1), ldda);
#endif
                /* m-(ii+sb) should be always 0 */
                magma_tally4_zgemm( Magma_tally4NoTrans, Magma_tally4NoTrans,
                             n-mindim, m-(ii+sb), sb,
                             c_neg_one, AT(i,   i+1), ldda,
                                        AT(i+1, i  ), ldda,
                             c_one,     AT(i+1, i+1), ldda );
            }
        }

        if ( order == Magma_tally4ColMajor ) {
            magma_tally4blas_zgetmo_out( dA, dAT, ldda, m, n );
        }
    }
    return *info;
}

#endif