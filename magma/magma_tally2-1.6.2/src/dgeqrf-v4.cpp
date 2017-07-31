/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from zgeqrf-v4.cpp normal z -> d, Fri Jan 30 19:00:15 2015

*/
#include "common_magma_tally2.h"

/**
    Purpose
    -------
    DGEQRF4 computes a QR factorization of a DOUBLE_PRECISION M-by-N matrix A:
    A = Q * R using multiple GPUs. This version does not require work space on the GPU
    passed as input. GPU memory is allocated in the routine.

    Arguments
    ---------
    @param[in]
    ngpu    INTEGER
            Number of GPUs to use. ngpu > 0.

    @param[in]
    m       INTEGER
            The number of rows of the matrix A.  M >= 0.

    @param[in]
    n       INTEGER
            The number of columns of the matrix A.  N >= 0.

    @param[in,out]
    A       DOUBLE_PRECISION array, dimension (LDA,N)
            On entry, the M-by-N matrix A.
            On exit, the elements on and above the diagonal of the array
            contain the min(M,N)-by-N upper trapezoidal matrix R (R is
            upper triangular if m >= n); the elements below the diagonal,
            with the array TAU, represent the orthogonal matrix Q as a
            product of min(m,n) elementary reflectors (see Further
            Details).
    \n
            Higher performance is achieved if A is in pinned memory, e.g.
            allocated using magma_tally2_malloc_pinned.

    @param[in]
    lda     INTEGER
            The leading dimension of the array A.  LDA >= max(1,M).

    @param[out]
    tau     DOUBLE_PRECISION array, dimension (min(M,N))
            The scalar factors of the elementary reflectors (see Further
            Details).

    @param[out]
    work    (workspace) DOUBLE_PRECISION array, dimension (MAX(1,LWORK))
            On exit, if INFO = 0, WORK[0] returns the optimal LWORK.
    \n
            Higher performance is achieved if WORK is in pinned memory, e.g.
            allocated using magma_tally2_malloc_pinned.

    @param[in]
    lwork   INTEGER
            The dimension of the array WORK.  LWORK >= N*NB,
            where NB can be obtained through magma_tally2_get_dgeqrf_nb(M).
    \n
            If LWORK = -1, then a workspace query is assumed; the routine
            only calculates the optimal size of the WORK array, returns
            this value as the first entry of the WORK array, and no error
            message related to LWORK is issued.

    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument had an illegal value
                  or another error occured, such as memory allocation failed.

    Further Details
    ---------------
    The matrix Q is represented as a product of elementary reflectors

       Q = H(1) H(2) . . . H(k), where k = min(m,n).

    Each H(i) has the form

       H(i) = I - tau * v * v'

    where tau is a real scalar, and v is a real vector with
    v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
    and tau in TAU(i).

    @ingroup magma_tally2_dgeqrf_comp
    ********************************************************************/
extern "C" magma_tally2_int_t
magma_tally2_dgeqrf4(
    magma_tally2_int_t ngpu,
    magma_tally2_int_t m, magma_tally2_int_t n,
    double *A,    magma_tally2_int_t lda, double *tau,
    double *work, magma_tally2_int_t lwork,
    magma_tally2_int_t *info )
{
    double *da[Magma_tally2MaxGPUs];
    double c_one = MAGMA_tally2_D_ONE;

    int i, k, ldda;

    *info = 0;
    int nb = magma_tally2_get_dgeqrf_nb(min(m, n));

    int lwkopt = n * nb;
    work[0] = MAGMA_tally2_D_MAKE( (double)lwkopt, 0 );
    int lquery = (lwork == -1);
    if (ngpu < 0 || ngpu > Magma_tally2MaxGPUs) {
        *info = -1;
    } else if (m < 0) {
        *info = -2;
    } else if (n < 0) {
        *info = -3;
    } else if (lda < max(1,m)) {
        *info = -5;
    } else if (lwork < max(1,n) && ! lquery) {
        *info = -8;
    }
    if (*info != 0) {
        magma_tally2_xerbla( __func__, -(*info) );
        return *info;
    }
    else if (lquery)
        return *info;

    k = min(m,n);
    if (k == 0) {
        work[0] = c_one;
        return *info;
    }

    magma_tally2_device_t orig_dev;
    magma_tally2_getdevice( &orig_dev );
    
    ldda    = ((m+31)/32)*32;

    magma_tally2_int_t  n_local[Magma_tally2MaxGPUs];
    for (i=0; i < ngpu; i++) {
        n_local[i] = ((n/nb)/ngpu)*nb;
        if (i < (n/nb)%ngpu)
            n_local[i] += nb;
        else if (i == (n/nb)%ngpu)
            n_local[i] += n%nb;

        magma_tally2_setdevice(i);
        
        // TODO on failure, free previously allocated memory
        if (MAGMA_tally2_SUCCESS != magma_tally2_dmalloc( &da[i], ldda*n_local[i] )) {
            *info = MAGMA_tally2_ERR_DEVICE_ALLOC;
            return *info;
        }
    }

    if (m > nb && n > nb) {
        /* Copy the matrix to the GPUs in 1D block cyclic distribution */
        magma_tally2_dsetmatrix_1D_col_bcyclic(m, n, A, lda, da, ldda, ngpu, nb);

        /* Factor using the GPU interface */
        magma_tally2_dgeqrf2_mgpu( ngpu, m, n, da, ldda, tau, info);

        /* Copy the matrix back from the GPUs to the CPU */
        magma_tally2_dgetmatrix_1D_col_bcyclic(m, n, da, ldda, A, lda, ngpu, nb);
    }
    else {
        lapackf77_dgeqrf(&m, &n, A, &lda, tau, work, &lwork, info);
    }


    /* Free the allocated GPU memory */
    for (i=0; i < ngpu; i++) {
        magma_tally2_setdevice(i);
        magma_tally2_free( da[i] );
    }
    magma_tally2_setdevice( orig_dev );

    return *info;
} /* magma_tally2_dgeqrf4 */
