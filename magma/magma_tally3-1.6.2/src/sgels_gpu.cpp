/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from zgels_gpu.cpp normal z -> s, Fri Jan 30 19:00:14 2015

*/
#include "common_magma_tally3.h"

/**
    Purpose
    -------
    Solves the overdetermined, least squares problem
           min || A*X - C ||
    using the QR factorization A.
    The underdetermined problem (m < n) is not currently handled.


    Arguments
    ---------
    @param[in]
    trans   magma_tally3_trans_t
      -     = Magma_tally3NoTrans:   the linear system involves A.
            Only TRANS=Magma_tally3NoTrans is currently handled.

    @param[in]
    m       INTEGER
            The number of rows of the matrix A. M >= 0.

    @param[in]
    n       INTEGER
            The number of columns of the matrix A. M >= N >= 0.

    @param[in]
    nrhs    INTEGER
            The number of columns of the matrix C. NRHS >= 0.

    @param[in,out]
    dA       REAL array on the GPU, dimension (LDA,N)
            On entry, the M-by-N matrix A.
            On exit, A is overwritten by details of its QR
            factorization as returned by SGEQRF.

    @param[in]
    ldda    INTEGER
            The leading dimension of the array A, LDDA >= M.

    @param[in,out]
    dB      REAL array on the GPU, dimension (LDDB,NRHS)
            On entry, the M-by-NRHS matrix C.
            On exit, the N-by-NRHS solution matrix X.

    @param[in]
    lddb    INTEGER
            The leading dimension of the array dB. LDDB >= M.

    @param[out]
    hwork   (workspace) REAL array, dimension MAX(1,LWORK).
            On exit, if INFO = 0, HWORK[0] returns the optimal LWORK.

    @param[in]
    lwork   INTEGER
            The dimension of the array HWORK,
            LWORK >= (M - N + NB)*(NRHS + NB) + NRHS*NB,
            where NB is the blocksize given by magma_tally3_get_sgeqrf_nb( M ).
    \n
            If LWORK = -1, then a workspace query is assumed; the routine
            only calculates the optimal size of the HWORK array, returns
            this value as the first entry of the HWORK array.

    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument had an illegal value

    @ingroup magma_tally3_sgels_driver
    ********************************************************************/
extern "C" magma_tally3_int_t
magma_tally3_sgels_gpu(
    magma_tally3_trans_t trans, magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t nrhs,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Float_ptr dB, magma_tally3_int_t lddb,
    float *hwork, magma_tally3_int_t lwork,
    magma_tally3_int_t *info)
{
    magma_tally3Float_ptr dT;
    float *tau;
    magma_tally3_int_t k;

    magma_tally3_int_t nb     = magma_tally3_get_sgeqrf_nb(m);
    magma_tally3_int_t lwkopt = (m - n + nb)*(nrhs + nb) + nrhs*nb;
    int lquery = (lwork == -1);

    hwork[0] = MAGMA_tally3_S_MAKE( (float)lwkopt, 0. );

    *info = 0;
    /* For now, N is the only case working */
    if ( trans != Magma_tally3NoTrans )
        *info = -1;
    else if (m < 0)
        *info = -2;
    else if (n < 0 || m < n) /* LQ is not handle for now*/
        *info = -3;
    else if (nrhs < 0)
        *info = -4;
    else if (ldda < max(1,m))
        *info = -6;
    else if (lddb < max(1,m))
        *info = -8;
    else if (lwork < lwkopt && ! lquery)
        *info = -10;

    if (*info != 0) {
        magma_tally3_xerbla( __func__, -(*info) );
        return *info;
    }
    else if (lquery)
        return *info;

    k = min(m,n);
    if (k == 0) {
        hwork[0] = MAGMA_tally3_S_ONE;
        return *info;
    }

    /*
     * Allocate temporary buffers
     */
    int ldtwork = ( 2*k + ((n+31)/32)*32 )*nb;
    if (nb < nrhs)
        ldtwork = ( 2*k + ((n+31)/32)*32 )*nrhs;
    if (MAGMA_tally3_SUCCESS != magma_tally3_smalloc( &dT, ldtwork )) {
        *info = MAGMA_tally3_ERR_DEVICE_ALLOC;
        return *info;
    }
    
    magma_tally3_smalloc_cpu( &tau, k );
    if ( tau == NULL ) {
        magma_tally3_free( dT );
        *info = MAGMA_tally3_ERR_HOST_ALLOC;
        return *info;
    }

    magma_tally3_sgeqrf_gpu( m, n, dA, ldda, tau, dT, info );

    if ( *info == 0 ) {
        magma_tally3_sgeqrs_gpu( m, n, nrhs,
                          dA, ldda, tau, dT,
                          dB, lddb, hwork, lwork, info );
    }
    
    magma_tally3_free( dT );
    magma_tally3_free_cpu(tau);
    return *info;
}
