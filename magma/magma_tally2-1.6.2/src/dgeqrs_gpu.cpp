/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from zgeqrs_gpu.cpp normal z -> d, Fri Jan 30 19:00:15 2015

*/
#include "common_magma_tally2.h"

/**
    Purpose
    -------
    Solves the least squares problem
           min || A*X - C ||
    using the QR factorization A = Q*R computed by DGEQRF_GPU.

    Arguments
    ---------
    @param[in]
    m       INTEGER
            The number of rows of the matrix A. M >= 0.

    @param[in]
    n       INTEGER
            The number of columns of the matrix A. M >= N >= 0.

    @param[in]
    nrhs    INTEGER
            The number of columns of the matrix C. NRHS >= 0.

    @param[in]
    dA      DOUBLE_PRECISION array on the GPU, dimension (LDDA,N)
            The i-th column must contain the vector which defines the
            elementary reflector H(i), for i = 1,2,...,n, as returned by
            DGEQRF_GPU in the first n columns of its array argument A.

    @param[in]
    ldda    INTEGER
            The leading dimension of the array A, LDDA >= M.

    @param[in]
    tau     DOUBLE_PRECISION array, dimension (N)
            TAU(i) must contain the scalar factor of the elementary
            reflector H(i), as returned by MAGMA_tally2_DGEQRF_GPU.

    @param[in,out]
    dB      DOUBLE_PRECISION array on the GPU, dimension (LDDB,NRHS)
            On entry, the M-by-NRHS matrix C.
            On exit, the N-by-NRHS solution matrix X.

    @param[in]
    dT      DOUBLE_PRECISION array that is the output (the 6th argument)
            of magma_tally2_dgeqrf_gpu of size
            2*MIN(M, N)*NB + ((N+31)/32*32 )* MAX(NB, NRHS).
            The array starts with a block of size MIN(M,N)*NB that stores
            the triangular T matrices used in the QR factorization,
            followed by MIN(M,N)*NB block storing the diagonal block
            inverses for the R matrix, followed by work space of size
            ((N+31)/32*32 )* MAX(NB, NRHS).

    @param[in]
    lddb    INTEGER
            The leading dimension of the array dB. LDDB >= M.

    @param[out]
    hwork   (workspace) DOUBLE_PRECISION array, dimension (LWORK)
            On exit, if INFO = 0, WORK[0] returns the optimal LWORK.

    @param[in]
    lwork   INTEGER
            The dimension of the array WORK,
            LWORK >= (M - N + NB)*(NRHS + NB) + NRHS*NB,
            where NB is the blocksize given by magma_tally2_get_dgeqrf_nb( M ).
    \n
            If LWORK = -1, then a workspace query is assumed; the routine
            only calculates the optimal size of the HWORK array, returns
            this value as the first entry of the WORK array.

    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument had an illegal value

    @ingroup magma_tally2_dgels_comp
    ********************************************************************/
extern "C" magma_tally2_int_t
magma_tally2_dgeqrs_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t nrhs,
    magma_tally2Double_ptr dA,    magma_tally2_int_t ldda,
    double *tau,
    magma_tally2Double_ptr dT,
    magma_tally2Double_ptr dB,    magma_tally2_int_t lddb,
    double *hwork, magma_tally2_int_t lwork,
    magma_tally2_int_t *info)
{
    #define dA(a_1,a_2) (dA + (a_2)*(ldda) + (a_1))
    #define dT(a_1)     (dT + (lddwork+(a_1))*nb)

    double c_zero    = MAGMA_tally2_D_ZERO;
    double c_one     = MAGMA_tally2_D_ONE;
    double c_neg_one = MAGMA_tally2_D_NEG_ONE;
    magma_tally2Double_ptr dwork;
    magma_tally2_int_t i, k, lddwork, rows, ib;
    magma_tally2_int_t ione = 1;

    magma_tally2_int_t nb     = magma_tally2_get_dgeqrf_nb(m);
    magma_tally2_int_t lwkopt = (m - n + nb)*(nrhs + nb) + nrhs*nb;
    int lquery = (lwork == -1);

    hwork[0] = MAGMA_tally2_D_MAKE( (double)lwkopt, 0. );

    *info = 0;
    if (m < 0)
        *info = -1;
    else if (n < 0 || m < n)
        *info = -2;
    else if (nrhs < 0)
        *info = -3;
    else if (ldda < max(1,m))
        *info = -5;
    else if (lddb < max(1,m))
        *info = -9;
    else if (lwork < lwkopt && ! lquery)
        *info = -11;

    if (*info != 0) {
        magma_tally2_xerbla( __func__, -(*info) );
        return *info;
    }
    else if (lquery)
        return *info;

    k = min(m,n);
    if (k == 0) {
        hwork[0] = c_one;
        return *info;
    }

    /* B := Q' * B */
    magma_tally2_dormqr_gpu( Magma_tally2Left, Magma_tally2Trans,
                      m, nrhs, n,
                      dA(0,0), ldda, tau,
                      dB, lddb, hwork, lwork, dT, nb, info );
    if ( *info != 0 ) {
        return *info;
    }

    /* Solve R*X = B(1:n,:) */
    lddwork= k;
    if (nb < k)
        dwork = dT+2*lddwork*nb;
    else
        dwork = dT;
    // To do: Why did we have this line originally; seems to be a bug (Stan)?
    // dwork = dT;

    i    = (k-1)/nb * nb;
    ib   = n-i;
    rows = m-i;

    // TODO: this assumes that, on exit from magma_tally2_dormqr_gpu, hwork contains
    // the last block of A and B (i.e., C in dormqr). This should be fixed.
    // Seems this data should already be on the GPU, so could switch to
    // magma_tally2_dtrsm and drop the dsetmatrix.
    if ( nrhs == 1 ) {
        blasf77_dtrsv( Magma_tally2UpperStr, Magma_tally2NoTransStr, Magma_tally2NonUnitStr,
                       &ib, hwork,         &rows,
                            hwork+rows*ib, &ione);
    } else {
        blasf77_dtrsm( Magma_tally2LeftStr, Magma_tally2UpperStr, Magma_tally2NoTransStr, Magma_tally2NonUnitStr,
                       &ib, &nrhs,
                       &c_one, hwork,         &rows,
                               hwork+rows*ib, &rows);
    }
    
    // update the solution vector
    magma_tally2_dsetmatrix( ib, nrhs, hwork+rows*ib, rows, dwork+i, lddwork );

    // update c
    if (nrhs == 1)
        magma_tally2_dgemv( Magma_tally2NoTrans, i, ib,
                     c_neg_one, dA(0, i), ldda,
                                dwork + i,   1,
                     c_one,     dB,           1);
    else
        magma_tally2_dgemm( Magma_tally2NoTrans, Magma_tally2NoTrans,
                     i, nrhs, ib,
                     c_neg_one, dA(0, i), ldda,
                                dwork + i,   lddwork,
                     c_one,     dB,           lddb);

    int start = i-nb;
    if (nb < k) {
        for (i = start; i >= 0; i -= nb) {
            ib = min(k-i, nb);
            rows = m -i;

            if (i + ib < n) {
                if (nrhs == 1) {
                    magma_tally2_dgemv( Magma_tally2NoTrans, ib, ib,
                                 c_one,  dT(i), ib,
                                         dB+i,      1,
                                 c_zero, dwork+i,  1);
                    magma_tally2_dgemv( Magma_tally2NoTrans, i, ib,
                                 c_neg_one, dA(0, i), ldda,
                                            dwork + i,   1,
                                 c_one,     dB,           1);
                }
                else {
                    magma_tally2_dgemm( Magma_tally2NoTrans, Magma_tally2NoTrans,
                                 ib, nrhs, ib,
                                 c_one,  dT(i), ib,
                                         dB+i,      lddb,
                                 c_zero, dwork+i,  lddwork);
                    magma_tally2_dgemm( Magma_tally2NoTrans, Magma_tally2NoTrans,
                                 i, nrhs, ib,
                                 c_neg_one, dA(0, i), ldda,
                                            dwork + i,   lddwork,
                                 c_one,     dB,          lddb);
                }
            }
        }
    }

    magma_tally2_dcopymatrix( (n), nrhs,
                       dwork, lddwork,
                       dB,    lddb );
    
    return *info;
}

#undef dA
#undef dT
