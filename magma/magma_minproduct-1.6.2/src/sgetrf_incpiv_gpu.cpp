/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Hatem Ltaief
       @author Mathieu Faverge

       @generated from zgetrf_incpiv_gpu.cpp normal z -> s, Fri Jan 30 19:00:14 2015

*/
#ifdef MAGMA_minproduct_WITH_PLASMA

#include <plasma.h>
#include <core_blas.h>
#include "common_magma_minproduct.h"

/**
    Purpose
    -------
    SGETRF_INCPIV computes an LU factorization of a general M-by-N tile A
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
    hA      REAL array, dimension(LDHA, N), on cpu.
            On entry, only the M-by-IB first panel needs to be identical to dA(1..M, 1..IB).
            On exit, the content is incomplete. Shouldn't be used.

    @param[in]
    ldha    INTEGER
            The leading dimension of the array hA.  LDHA >= max(1,M).

    @param[in,out]
    dA      REAL array, dimension(LDDA, N), on gpu.
            On entry, the M-by-N tile to be factored.
            On exit, the factors L and U from the factorization
            A = P*L*U; the unit diagonal elements of L are not stored.

    @param[in]
    ldda    INTEGER
            The leading dimension of the array dA.  LDDA >= max(1,M).

    @param[out]
    hL      REAL array, dimension(LDHL, min(M,N)), on vpu.
            On exit, contains in the upper part the IB-by-K lower triangular tile,
            and in the lower part IB-by-min(M,N) the inverse of the top part.

    @param[in]
    ldhl    INTEGER
            The leading dimension of the array hL.  LDHL >= max(1,2*IB).

    @param[out]
    dL      REAL array, dimension(LDDL, K), on gpu.
            On exit, contains in the upper part the IB-by-min(M,N) lower triangular tile,
            and in the lower part IB-by-min(M,N) the inverse of the top part.

    @param[in]
    lddl    INTEGER
            The leading dimension of the array dL.  LDDL >= max(1,2*IB).

    @param[out]
    ipiv    INTEGER array, dimension min(M,N), on the cpu.
            The pivot indices array.

    @param[out]
    dWORK   REAL array, dimension(LDDWORK, 2*IB), on gpu.
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

    @ingroup magma_minproduct_sgesv_comp
    ********************************************************************/
extern "C" magma_minproduct_int_t
magma_minproduct_sgetrf_incpiv_gpu(
    magma_minproduct_order_t order, magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t ib,
    float    *hA, magma_minproduct_int_t ldha,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda,
    float    *hL, magma_minproduct_int_t ldhl,
    magma_minproductFloat_ptr dL, magma_minproduct_int_t lddl,
    magma_minproduct_int_t *ipiv,
    magma_minproductFloat_ptr dwork, magma_minproduct_int_t lddwork,
    magma_minproduct_int_t *info)
{
#define AT(i,j) (dAT + (i)*ib*ldda + (j)*ib)
#define hA(i,j) (hA  + (i)*ib + (j)*ib*ldha)
#define hL(j)   (hL  + (j)*ib*ldhl         )
#define hL2(j)  (hL2 + (j)*ib*ldhl         )
#define dL(j)   (dL  + (j)*ib*lddl         )
#define dL2(j)  (dL2 + (j)*ib*lddl         )

    float c_one     = MAGMA_minproduct_S_ONE;
    float c_neg_one = MAGMA_minproduct_S_NEG_ONE;

    magma_minproduct_int_t iinfo;
    magma_minproduct_int_t maxm, mindim;
    magma_minproduct_int_t i, rows, cols, s, ii, sb;
    magma_minproductFloat_ptr dAT;
#ifndef WITHOUTTRTRI
    magma_minproductFloat_ptr dL2 = dL + ib;
    float *hL2 = hL + ib;
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
        magma_minproduct_xerbla( __func__, -(*info) );
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
        lapackf77_sgetrf(&m, &n, hA, &ldha, ipiv, info);

#ifndef WITHOUTTRTRI
        CORE_slacpy(PlasmaUpperLower, mindim, mindim,
                    (float*)hA, ldha,
                    (float*)hL2, ldhl );

        CORE_strtri( PlasmaLower, PlasmaUnit, mindim,
                     (float*)hL2, ldhl, info );
        if (*info != 0 ) {
            fprintf(stderr, "ERROR, trtri returned with info = %d\n", *info);
        }

        magma_minproduct_ssetmatrix( mindim, mindim, hL2, ldhl, dL2, lddl );
#endif

        if ( order == Magma_minproductRowMajor ) {
            magma_minproduct_ssetmatrix( m, n, hA, ldha, dwork, lddwork );
            magma_minproductblas_stranspose( m, n, dwork, lddwork, dA, ldda );
        } else {
            magma_minproduct_ssetmatrix( m, n, hA, ldha, dA, ldda );
        }
    }
    else {
        /* Use hybrid blocked code. */
        maxm = ((m + 31)/32)*32;

        if ( order == Magma_minproductColMajor ) {
            magma_minproductblas_sgetmo_in( dA, dAT, ldda, m, n );
        } else {
            dAT = dA;
        }

        for( i=0; i < s; i++ ) {
            ii = i * ib;
            sb = min(ib, mindim-ii);
            cols = maxm - ii;

            if ( i > 0 ) {
                // download i-th panel
                magma_minproductblas_stranspose( sb, m, AT(0,i), ldda, dwork, maxm );
                magma_minproduct_sgetmatrix( m, sb, dwork, maxm, hA(0, i), ldha );

                // make sure that gpu queue is empty
                //magma_minproduct_device_sync();
#ifndef WITHOUTTRTRI
                magma_minproduct_strmm( Magma_minproductRight, Magma_minproductLower, Magma_minproductTrans, Magma_minproductUnit,
                             n - (ii+sb), ib,
                             c_one, dL2(i-1),    lddl,
                                    AT(i-1,i+1), ldda );
#else
                magma_minproduct_strsm( Magma_minproductRight, Magma_minproductUpper, Magma_minproductNoTrans, Magma_minproductUnit,
                             n - (ii+sb), ib,
                             c_one, AT(i-1,i-1), ldda,
                                    AT(i-1,i+1), ldda );
#endif
                magma_minproduct_sgemm( Magma_minproductNoTrans, Magma_minproductNoTrans,
                             n-(ii+sb), m-ii, ib,
                             c_neg_one, AT(i-1,i+1), ldda,
                                        AT(i,  i-1), ldda,
                             c_one,     AT(i,  i+1), ldda );
            }

            // do the cpu part
            rows = m - ii;
            lapackf77_sgetrf( &rows, &sb, hA(i, i), &ldha, ipiv+ii, &iinfo);
            if ( (*info == 0) && (iinfo > 0) )
                *info = iinfo + ii;

            {
                int j;
                int fin = ii + sb;
                for (j=ii; j < fin; j++) {
                    ipiv[j] = ii + ipiv[j];
                }
            }
            magma_minproductblas_slaswp( n-ii, AT(0, i), ldda, ii+1, ii+sb, ipiv, 1 );

#ifndef WITHOUTTRTRI
            CORE_slacpy(PlasmaLower, sb, sb,
                        (float*)hA(i, i), ldha,
                        (float*)hL2(i), ldhl );

            CORE_strtri( PlasmaLower, PlasmaUnit, sb,
                         (float*)hL2(i), ldhl, info );
            if (*info != 0 ) {
                fprintf(stderr, "ERROR, trtri returned with info = %d\n", *info);
            }
            magma_minproduct_ssetmatrix( sb, sb, hL2(i), ldhl, dL2(i), lddl );
#endif
            // upload i-th panel
            magma_minproduct_ssetmatrix( rows, sb, hA(i, i), ldha, dwork, cols );
            magma_minproductblas_stranspose( rows, sb, dwork, cols, AT(i,i), ldda );

            // do the small non-parallel computations
            if ( s > (i+1) ) {
#ifndef WITHOUTTRTRI
                magma_minproduct_strmm( Magma_minproductRight, Magma_minproductLower, Magma_minproductTrans, Magma_minproductUnit,
                             sb, sb,
                             c_one, dL2(i),     lddl,
                                    AT(i, i+1), ldda);
#else
                magma_minproduct_strsm( Magma_minproductRight, Magma_minproductUpper, Magma_minproductNoTrans, Magma_minproductUnit,
                             sb, sb,
                             c_one, AT(i, i  ), ldda,
                                    AT(i, i+1), ldda);
#endif
                magma_minproduct_sgemm( Magma_minproductNoTrans, Magma_minproductNoTrans,
                             sb, m-(ii+sb), sb,
                             c_neg_one, AT(i,   i+1), ldda,
                                        AT(i+1, i  ), ldda,
                             c_one,     AT(i+1, i+1), ldda );
            }
            else {
                /* Update of the last panel */
#ifndef WITHOUTTRTRI
                magma_minproduct_strmm( Magma_minproductRight, Magma_minproductLower, Magma_minproductTrans, Magma_minproductUnit,
                             n-mindim, sb,
                             c_one, dL2(i),     lddl,
                                    AT(i, i+1), ldda);
#else
                magma_minproduct_strsm( Magma_minproductRight, Magma_minproductUpper, Magma_minproductNoTrans, Magma_minproductUnit,
                             n-mindim, sb,
                             c_one, AT(i, i  ), ldda,
                                    AT(i, i+1), ldda);
#endif
                /* m-(ii+sb) should be always 0 */
                magma_minproduct_sgemm( Magma_minproductNoTrans, Magma_minproductNoTrans,
                             n-mindim, m-(ii+sb), sb,
                             c_neg_one, AT(i,   i+1), ldda,
                                        AT(i+1, i  ), ldda,
                             c_one,     AT(i+1, i+1), ldda );
            }
        }

        if ( order == Magma_minproductColMajor ) {
            magma_minproductblas_sgetmo_out( dA, dAT, ldda, m, n );
        }
    }
    return *info;
}

#endif
