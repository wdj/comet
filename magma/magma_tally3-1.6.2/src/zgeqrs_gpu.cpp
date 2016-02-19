/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions normal z -> s d c

*/
#include "common_magma_tally3.h"

/**
    Purpose
    -------
    Solves the least squares problem
           min || A*X - C ||
    using the QR factorization A = Q*R computed by ZGEQRF_GPU.

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
    dA      COMPLEX_16 array on the GPU, dimension (LDDA,N)
            The i-th column must contain the vector which defines the
            elementary reflector H(i), for i = 1,2,...,n, as returned by
            ZGEQRF_GPU in the first n columns of its array argument A.

    @param[in]
    ldda    INTEGER
            The leading dimension of the array A, LDDA >= M.

    @param[in]
    tau     COMPLEX_16 array, dimension (N)
            TAU(i) must contain the scalar factor of the elementary
            reflector H(i), as returned by MAGMA_tally3_ZGEQRF_GPU.

    @param[in,out]
    dB      COMPLEX_16 array on the GPU, dimension (LDDB,NRHS)
            On entry, the M-by-NRHS matrix C.
            On exit, the N-by-NRHS solution matrix X.

    @param[in]
    dT      COMPLEX_16 array that is the output (the 6th argument)
            of magma_tally3_zgeqrf_gpu of size
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
    hwork   (workspace) COMPLEX_16 array, dimension (LWORK)
            On exit, if INFO = 0, WORK[0] returns the optimal LWORK.

    @param[in]
    lwork   INTEGER
            The dimension of the array WORK,
            LWORK >= (M - N + NB)*(NRHS + NB) + NRHS*NB,
            where NB is the blocksize given by magma_tally3_get_zgeqrf_nb( M ).
    \n
            If LWORK = -1, then a workspace query is assumed; the routine
            only calculates the optimal size of the HWORK array, returns
            this value as the first entry of the WORK array.

    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument had an illegal value

    @ingroup magma_tally3_zgels_comp
    ********************************************************************/
extern "C" magma_tally3_int_t
magma_tally3_zgeqrs_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t nrhs,
    magma_tally3DoubleComplex_ptr dA,    magma_tally3_int_t ldda,
    magma_tally3DoubleComplex *tau,
    magma_tally3DoubleComplex_ptr dT,
    magma_tally3DoubleComplex_ptr dB,    magma_tally3_int_t lddb,
    magma_tally3DoubleComplex *hwork, magma_tally3_int_t lwork,
    magma_tally3_int_t *info)
{
    #define dA(a_1,a_2) (dA + (a_2)*(ldda) + (a_1))
    #define dT(a_1)     (dT + (lddwork+(a_1))*nb)

    magma_tally3DoubleComplex c_zero    = MAGMA_tally3_Z_ZERO;
    magma_tally3DoubleComplex c_one     = MAGMA_tally3_Z_ONE;
    magma_tally3DoubleComplex c_neg_one = MAGMA_tally3_Z_NEG_ONE;
    magma_tally3DoubleComplex_ptr dwork;
    magma_tally3_int_t i, k, lddwork, rows, ib;
    magma_tally3_int_t ione = 1;

    magma_tally3_int_t nb     = magma_tally3_get_zgeqrf_nb(m);
    magma_tally3_int_t lwkopt = (m - n + nb)*(nrhs + nb) + nrhs*nb;
    int lquery = (lwork == -1);

    hwork[0] = MAGMA_tally3_Z_MAKE( (double)lwkopt, 0. );

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
        magma_tally3_xerbla( __func__, -(*info) );
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
    magma_tally3_zunmqr_gpu( Magma_tally3Left, Magma_tally3_ConjTrans,
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

    // TODO: this assumes that, on exit from magma_tally3_zunmqr_gpu, hwork contains
    // the last block of A and B (i.e., C in zunmqr). This should be fixed.
    // Seems this data should already be on the GPU, so could switch to
    // magma_tally3_ztrsm and drop the zsetmatrix.
    if ( nrhs == 1 ) {
        blasf77_ztrsv( Magma_tally3UpperStr, Magma_tally3NoTransStr, Magma_tally3NonUnitStr,
                       &ib, hwork,         &rows,
                            hwork+rows*ib, &ione);
    } else {
        blasf77_ztrsm( Magma_tally3LeftStr, Magma_tally3UpperStr, Magma_tally3NoTransStr, Magma_tally3NonUnitStr,
                       &ib, &nrhs,
                       &c_one, hwork,         &rows,
                               hwork+rows*ib, &rows);
    }
    
    // update the solution vector
    magma_tally3_zsetmatrix( ib, nrhs, hwork+rows*ib, rows, dwork+i, lddwork );

    // update c
    if (nrhs == 1)
        magma_tally3_zgemv( Magma_tally3NoTrans, i, ib,
                     c_neg_one, dA(0, i), ldda,
                                dwork + i,   1,
                     c_one,     dB,           1);
    else
        magma_tally3_zgemm( Magma_tally3NoTrans, Magma_tally3NoTrans,
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
                    magma_tally3_zgemv( Magma_tally3NoTrans, ib, ib,
                                 c_one,  dT(i), ib,
                                         dB+i,      1,
                                 c_zero, dwork+i,  1);
                    magma_tally3_zgemv( Magma_tally3NoTrans, i, ib,
                                 c_neg_one, dA(0, i), ldda,
                                            dwork + i,   1,
                                 c_one,     dB,           1);
                }
                else {
                    magma_tally3_zgemm( Magma_tally3NoTrans, Magma_tally3NoTrans,
                                 ib, nrhs, ib,
                                 c_one,  dT(i), ib,
                                         dB+i,      lddb,
                                 c_zero, dwork+i,  lddwork);
                    magma_tally3_zgemm( Magma_tally3NoTrans, Magma_tally3NoTrans,
                                 i, nrhs, ib,
                                 c_neg_one, dA(0, i), ldda,
                                            dwork + i,   lddwork,
                                 c_one,     dB,          lddb);
                }
            }
        }
    }

    magma_tally3_zcopymatrix( (n), nrhs,
                       dwork, lddwork,
                       dB,    lddb );
    
    return *info;
}

#undef dA
#undef dT
