/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from zgeqrs_gpu.cpp normal z -> c, Fri Jan 30 19:00:16 2015

*/
#include "common_magma_tally4.h"

/**
    Purpose
    -------
    Solves the least squares problem
           min || A*X - C ||
    using the QR factorization A = Q*R computed by CGEQRF_GPU.

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
    dA      COMPLEX array on the GPU, dimension (LDDA,N)
            The i-th column must contain the vector which defines the
            elementary reflector H(i), for i = 1,2,...,n, as returned by
            CGEQRF_GPU in the first n columns of its array argument A.

    @param[in]
    ldda    INTEGER
            The leading dimension of the array A, LDDA >= M.

    @param[in]
    tau     COMPLEX array, dimension (N)
            TAU(i) must contain the scalar factor of the elementary
            reflector H(i), as returned by MAGMA_tally4_CGEQRF_GPU.

    @param[in,out]
    dB      COMPLEX array on the GPU, dimension (LDDB,NRHS)
            On entry, the M-by-NRHS matrix C.
            On exit, the N-by-NRHS solution matrix X.

    @param[in]
    dT      COMPLEX array that is the output (the 6th argument)
            of magma_tally4_cgeqrf_gpu of size
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
    hwork   (workspace) COMPLEX array, dimension (LWORK)
            On exit, if INFO = 0, WORK[0] returns the optimal LWORK.

    @param[in]
    lwork   INTEGER
            The dimension of the array WORK,
            LWORK >= (M - N + NB)*(NRHS + NB) + NRHS*NB,
            where NB is the blocksize given by magma_tally4_get_cgeqrf_nb( M ).
    \n
            If LWORK = -1, then a workspace query is assumed; the routine
            only calculates the optimal size of the HWORK array, returns
            this value as the first entry of the WORK array.

    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument had an illegal value

    @ingroup magma_tally4_cgels_comp
    ********************************************************************/
extern "C" magma_tally4_int_t
magma_tally4_cgeqrs_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t nrhs,
    magma_tally4FloatComplex_ptr dA,    magma_tally4_int_t ldda,
    magma_tally4FloatComplex *tau,
    magma_tally4FloatComplex_ptr dT,
    magma_tally4FloatComplex_ptr dB,    magma_tally4_int_t lddb,
    magma_tally4FloatComplex *hwork, magma_tally4_int_t lwork,
    magma_tally4_int_t *info)
{
    #define dA(a_1,a_2) (dA + (a_2)*(ldda) + (a_1))
    #define dT(a_1)     (dT + (lddwork+(a_1))*nb)

    magma_tally4FloatComplex c_zero    = MAGMA_tally4_C_ZERO;
    magma_tally4FloatComplex c_one     = MAGMA_tally4_C_ONE;
    magma_tally4FloatComplex c_neg_one = MAGMA_tally4_C_NEG_ONE;
    magma_tally4FloatComplex_ptr dwork;
    magma_tally4_int_t i, k, lddwork, rows, ib;
    magma_tally4_int_t ione = 1;

    magma_tally4_int_t nb     = magma_tally4_get_cgeqrf_nb(m);
    magma_tally4_int_t lwkopt = (m - n + nb)*(nrhs + nb) + nrhs*nb;
    int lquery = (lwork == -1);

    hwork[0] = MAGMA_tally4_C_MAKE( (float)lwkopt, 0. );

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
        magma_tally4_xerbla( __func__, -(*info) );
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
    magma_tally4_cunmqr_gpu( Magma_tally4Left, Magma_tally4_ConjTrans,
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

    // TODO: this assumes that, on exit from magma_tally4_cunmqr_gpu, hwork contains
    // the last block of A and B (i.e., C in cunmqr). This should be fixed.
    // Seems this data should already be on the GPU, so could switch to
    // magma_tally4_ctrsm and drop the csetmatrix.
    if ( nrhs == 1 ) {
        blasf77_ctrsv( Magma_tally4UpperStr, Magma_tally4NoTransStr, Magma_tally4NonUnitStr,
                       &ib, hwork,         &rows,
                            hwork+rows*ib, &ione);
    } else {
        blasf77_ctrsm( Magma_tally4LeftStr, Magma_tally4UpperStr, Magma_tally4NoTransStr, Magma_tally4NonUnitStr,
                       &ib, &nrhs,
                       &c_one, hwork,         &rows,
                               hwork+rows*ib, &rows);
    }
    
    // update the solution vector
    magma_tally4_csetmatrix( ib, nrhs, hwork+rows*ib, rows, dwork+i, lddwork );

    // update c
    if (nrhs == 1)
        magma_tally4_cgemv( Magma_tally4NoTrans, i, ib,
                     c_neg_one, dA(0, i), ldda,
                                dwork + i,   1,
                     c_one,     dB,           1);
    else
        magma_tally4_cgemm( Magma_tally4NoTrans, Magma_tally4NoTrans,
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
                    magma_tally4_cgemv( Magma_tally4NoTrans, ib, ib,
                                 c_one,  dT(i), ib,
                                         dB+i,      1,
                                 c_zero, dwork+i,  1);
                    magma_tally4_cgemv( Magma_tally4NoTrans, i, ib,
                                 c_neg_one, dA(0, i), ldda,
                                            dwork + i,   1,
                                 c_one,     dB,           1);
                }
                else {
                    magma_tally4_cgemm( Magma_tally4NoTrans, Magma_tally4NoTrans,
                                 ib, nrhs, ib,
                                 c_one,  dT(i), ib,
                                         dB+i,      lddb,
                                 c_zero, dwork+i,  lddwork);
                    magma_tally4_cgemm( Magma_tally4NoTrans, Magma_tally4NoTrans,
                                 i, nrhs, ib,
                                 c_neg_one, dA(0, i), ldda,
                                            dwork + i,   lddwork,
                                 c_one,     dB,          lddb);
                }
            }
        }
    }

    magma_tally4_ccopymatrix( (n), nrhs,
                       dwork, lddwork,
                       dB,    lddb );
    
    return *info;
}

#undef dA
#undef dT