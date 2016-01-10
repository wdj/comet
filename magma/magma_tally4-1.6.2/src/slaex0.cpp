/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015
       
       @author Raffaele Solca
       
       @generated from dlaex0.cpp normal d -> s, Fri Jan 30 19:00:17 2015
*/
#include "common_magma_tally4.h"
#include "magma_tally4_timer.h"

/**
    Purpose
    -------
    SLAEX0 computes all eigenvalues and the choosen eigenvectors of a
    symmetric tridiagonal matrix using the divide and conquer method.

    Arguments
    ---------
    @param[in]
    n       INTEGER
            The dimension of the symmetric tridiagonal matrix.  N >= 0.
            
    @param[in,out]
    d       REAL array, dimension (N)
            On entry, the main diagonal of the tridiagonal matrix.
            On exit, its eigenvalues.
            
    @param[in]
    e       REAL array, dimension (N-1)
            The off-diagonal elements of the tridiagonal matrix.
            On exit, E has been destroyed.
            
    @param[in,out]
    Q       REAL array, dimension (LDQ, N)
            On entry, Q will be the identity matrix.
            On exit, Q contains the eigenvectors of the
            tridiagonal matrix.
            
    @param[in]
    ldq     INTEGER
            The leading dimension of the array Q.  If eigenvectors are
            desired, then  LDQ >= max(1,N).  In any case,  LDQ >= 1.
            
    @param
    work    (workspace) REAL array,
            the dimension of WORK >= 4*N + N**2.
            
    @param
    iwork   (workspace) INTEGER array,
            the dimension of IWORK >= 3 + 5*N.
            
    @param
    dwork   (workspace) REAL array, dimension (3*N*N/2+3*N)
            
    @param[in]
    range   magma_tally4_range_t
      -     = Magma_tally4RangeAll: all eigenvalues will be found.
      -     = Magma_tally4RangeV:   all eigenvalues in the half-open interval (VL,VU]
                             will be found.
      -     = Magma_tally4RangeI:   the IL-th through IU-th eigenvalues will be found.
            
    @param[in]
    vl      REAL
    @param[in]
    vu      REAL
            If RANGE=Magma_tally4RangeV, the lower and upper bounds of the interval to
            be searched for eigenvalues. VL < VU.
            Not referenced if RANGE = Magma_tally4RangeAll or Magma_tally4RangeI.
            
    @param[in]
    il      INTEGER
    @param[in]
    iu      INTEGER
            If RANGE=Magma_tally4RangeI, the indices (in ascending order) of the
            smallest and largest eigenvalues to be returned.
            1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
            Not referenced if RANGE = Magma_tally4RangeAll or Magma_tally4RangeV.
            
    @param[out]
    info    INTEGER
      -     = 0:  successful exit.
      -     < 0:  if INFO = -i, the i-th argument had an illegal value.
      -     > 0:  The algorithm failed to compute an eigenvalue while
                  working on the submatrix lying in rows and columns
                  INFO/(N+1) through mod(INFO,N+1).

    Further Details
    ---------------
    Based on contributions by
       Jeff Rutter, Computer Science Division, University of California
       at Berkeley, USA

    @ingroup magma_tally4_ssyev_aux
    ********************************************************************/
extern "C" magma_tally4_int_t
magma_tally4_slaex0(
    magma_tally4_int_t n,
    float *d, float *e,
    float *Q, magma_tally4_int_t ldq,
    float *work, magma_tally4_int_t *iwork,
    magma_tally4Float_ptr dwork,
    magma_tally4_range_t range, float vl, float vu,
    magma_tally4_int_t il, magma_tally4_int_t iu,
    magma_tally4_int_t *info)
{
#define Q(i_,j_) (Q + (i_) + (j_)*ldq)

    magma_tally4_int_t ione = 1;
    magma_tally4_range_t range2;
    magma_tally4_int_t curlvl, i, indxq;
    magma_tally4_int_t j, k, matsiz, msd2, smlsiz;
    magma_tally4_int_t submat, subpbs, tlvls;


    // Test the input parameters.
    *info = 0;

    if ( n < 0 )
        *info = -1;
    else if ( ldq < max(1, n) )
        *info = -5;
    if ( *info != 0 ) {
        magma_tally4_xerbla( __func__, -(*info) );
        return *info;
    }

    // Quick return if possible
    if (n == 0)
        return *info;

    smlsiz = magma_tally4_get_smlsize_divideconquer();

    // Determine the size and placement of the submatrices, and save in
    // the leading elements of IWORK.
    iwork[0] = n;
    subpbs= 1;
    tlvls = 0;
    while (iwork[subpbs - 1] > smlsiz) {
        for (j = subpbs; j > 0; --j) {
            iwork[2*j - 1] = (iwork[j-1]+1)/2;
            iwork[2*j - 2] = iwork[j-1]/2;
        }
        ++tlvls;
        subpbs *= 2;
    }
    for (j=1; j < subpbs; ++j)
        iwork[j] += iwork[j-1];

    // Divide the matrix into SUBPBS submatrices of size at most SMLSIZ+1
    // using rank-1 modifications (cuts).
    for (i=0; i < subpbs-1; ++i) {
        submat = iwork[i];
        d[submat-1] -= MAGMA_tally4_S_ABS(e[submat-1]);
        d[submat] -= MAGMA_tally4_S_ABS(e[submat-1]);
    }

    indxq = 4*n + 3;

    // Solve each submatrix eigenproblem at the bottom of the divide and
    // conquer tree.
    magma_tally4_timer_t time=0;
    timer_start( time );

    for (i = 0; i < subpbs; ++i) {
        if (i == 0) {
            submat = 0;
            matsiz = iwork[0];
        } else {
            submat = iwork[i-1];
            matsiz = iwork[i] - iwork[i-1];
        }
        lapackf77_ssteqr("I", &matsiz, &d[submat], &e[submat],
                         Q(submat, submat), &ldq, work, info);  // change to edc?
        if (*info != 0) {
            printf("info: %d\n, submat: %d\n", (int) *info, (int) submat);
            *info = (submat+1)*(n+1) + submat + matsiz;
            printf("info: %d\n", (int) *info);
            return *info;
        }
        k = 1;
        for (j = submat; j < iwork[i]; ++j) {
            iwork[indxq+j] = k;
            ++k;
        }
    }

    timer_stop( time );
    timer_printf( "  for: ssteqr = %6.2f\n", time );
    
    // Successively merge eigensystems of adjacent submatrices
    // into eigensystem for the corresponding larger matrix.
    curlvl = 1;
    while (subpbs > 1) {
        timer_start( time );
        
        for (i=0; i < subpbs-1; i += 2) {
            if (i == 0) {
                submat = 0;
                matsiz = iwork[1];
                msd2 = iwork[0];
            } else {
                submat = iwork[i-1];
                matsiz = iwork[i+1] - iwork[i-1];
                msd2 = matsiz / 2;
            }

            // Merge lower order eigensystems (of size MSD2 and MATSIZ - MSD2)
            // into an eigensystem of size MATSIZ.
            // SLAEX1 is used only for the full eigensystem of a tridiagonal
            // matrix.
            if (matsiz == n)
                range2 = range;
            else
                // We need all the eigenvectors if it is not last step
                range2 = Magma_tally4RangeAll;

            magma_tally4_slaex1(matsiz, &d[submat], Q(submat, submat), ldq,
                         &iwork[indxq+submat], e[submat+msd2-1], msd2,
                         work, &iwork[subpbs], dwork,
                         range2, vl, vu, il, iu, info);

            if (*info != 0) {
                *info = (submat+1)*(n+1) + submat + matsiz;
                return *info;
            }
            iwork[i/2]= iwork[i+1];
        }
        subpbs /= 2;
        ++curlvl;
        
        timer_stop( time );
        timer_printf("%d: time: %6.2f\n", (int) curlvl, time );
    }

    // Re-merge the eigenvalues/vectors which were deflated at the final
    // merge step.
    for (i = 0; i < n; ++i) {
        j = iwork[indxq+i] - 1;
        work[i] = d[j];
        blasf77_scopy(&n, Q(0, j), &ione, &work[ n*(i+1) ], &ione);
    }
    blasf77_scopy(&n, work, &ione, d, &ione);
    lapackf77_slacpy( "A", &n, &n, &work[n], &n, Q, &ldq );

    return *info;
} /* magma_tally4_slaex0 */