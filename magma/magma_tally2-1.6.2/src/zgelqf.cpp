/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions normal z -> s d c

*/
#include "common_magma_tally2.h"

#define COMPLEX

/**
    Purpose
    -------
    ZGELQF computes an LQ factorization of a COMPLEX_16 M-by-N matrix A:
    A = L * Q.

    Arguments
    ---------
    @param[in]
    m       INTEGER
            The number of rows of the matrix A.  M >= 0.

    @param[in]
    n       INTEGER
            The number of columns of the matrix A.  N >= 0.

    @param[in,out]
    A       COMPLEX_16 array, dimension (LDA,N)
            On entry, the M-by-N matrix A.
            On exit, the elements on and below the diagonal of the array
            contain the m-by-min(m,n) lower trapezoidal matrix L (L is
            lower triangular if m <= n); the elements above the diagonal,
            with the array TAU, represent the orthogonal matrix Q as a
            product of elementary reflectors (see Further Details).
    \n
            Higher performance is achieved if A is in pinned memory, e.g.
            allocated using magma_tally2_malloc_pinned.

    @param[in]
    lda     INTEGER
            The leading dimension of the array A.  LDA >= max(1,M).

    @param[out]
    tau     COMPLEX_16 array, dimension (min(M,N))
            The scalar factors of the elementary reflectors (see Further
            Details).

    @param[out]
    work    (workspace) COMPLEX_16 array, dimension (MAX(1,LWORK))
            On exit, if INFO = 0, WORK[0] returns the optimal LWORK.
    \n
            Higher performance is achieved if WORK is in pinned memory, e.g.
            allocated using magma_tally2_malloc_pinned.

    @param[in]
    lwork   INTEGER
            The dimension of the array WORK.  LWORK >= max(1,M).
            For optimum performance LWORK >= M*NB, where NB is the
            optimal blocksize.
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

       Q = H(k) . . . H(2) H(1), where k = min(m,n).

    Each H(i) has the form

       H(i) = I - tau * v * v'

    where tau is a complex scalar, and v is a complex vector with
    v(1:i-1) = 0 and v(i) = 1; v(i+1:n) is stored on exit in A(i,i+1:n),
    and tau in TAU(i).

    @ingroup magma_tally2_zgelqf_comp
    ********************************************************************/
extern "C" magma_tally2_int_t
magma_tally2_zgelqf(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex *A,    magma_tally2_int_t lda,   magma_tally2DoubleComplex *tau,
    magma_tally2DoubleComplex *work, magma_tally2_int_t lwork,
    magma_tally2_int_t *info)
{
    #define  dA(i_, j_)  (dA  + (i_) + (j_)*ldda)
    #define dAT(i_, j_)  (dAT + (i_) + (j_)*ldda)
    
    const magma_tally2DoubleComplex c_one = MAGMA_tally2_Z_ONE;
    const magma_tally2_int_t        ione  = 1;
    MAGMA_tally2_UNUSED( ione );  // used only for complex
    
    magma_tally2DoubleComplex_ptr dA, dAT;
    magma_tally2_int_t min_mn, maxm, maxn, maxdim, nb;
    magma_tally2_int_t iinfo, ldda, lddat;
    int lquery;

    /* Function Body */
    *info = 0;
    nb = magma_tally2_get_zgelqf_nb(m);
    min_mn = min(m,n);

    work[0] = MAGMA_tally2_Z_MAKE( (double)(m*nb), 0 );
    lquery = (lwork == -1);
    if (m < 0) {
        *info = -1;
    } else if (n < 0) {
        *info = -2;
    } else if (lda < max(1,m)) {
        *info = -4;
    } else if (lwork < max(1,m) && ! lquery) {
        *info = -7;
    }
    if (*info != 0) {
        magma_tally2_xerbla( __func__, -(*info) );
        return *info;
    }
    else if (lquery) {
        return *info;
    }

    /*  Quick return if possible */
    if (min_mn == 0) {
        work[0] = c_one;
        return *info;
    }

    maxm = ((m + 31)/32)*32;
    maxn = ((n + 31)/32)*32;
    maxdim = max(maxm, maxn);

    // copy to GPU and transpose
    if (maxdim*maxdim < 2*maxm*maxn) {
        // close to square, do everything in-place
        ldda  = maxdim;
        lddat = maxdim;

        if (MAGMA_tally2_SUCCESS != magma_tally2_zmalloc( &dA, maxdim*maxdim )) {
            *info = MAGMA_tally2_ERR_DEVICE_ALLOC;
            return *info;
        }

        magma_tally2_zsetmatrix( m, n, A, lda, dA(0,0), ldda );
        dAT = dA;
        magma_tally2blas_ztranspose_inplace( lddat, dAT(0,0), lddat );
    }
    else {
        // rectangular, do everything out-of-place
        ldda  = maxm;
        lddat = maxn;

        if (MAGMA_tally2_SUCCESS != magma_tally2_zmalloc( &dA, 2*maxn*maxm )) {
            *info = MAGMA_tally2_ERR_DEVICE_ALLOC;
            return *info;
        }

        magma_tally2_zsetmatrix( m, n, A, lda, dA(0,0), ldda );

        dAT = dA + maxn * maxm;
        magma_tally2blas_ztranspose( m, n, dA(0,0), ldda, dAT(0,0), lddat );
    }

    // factor QR
    magma_tally2_zgeqrf2_gpu( n, m, dAT(0,0), lddat, tau, &iinfo );
    assert( iinfo >= 0 );
    if ( iinfo > 0 ) {
        *info = iinfo;
    }
    
    // conjugate tau
    #ifdef COMPLEX
    lapackf77_zlacgv( &min_mn, tau, &ione );
    #endif

    // undo transpose
    if (maxdim*maxdim < 2*maxm*maxn) {
        magma_tally2blas_ztranspose_inplace( lddat, dAT(0,0), lddat );
        magma_tally2_zgetmatrix( m, n, dA(0,0), ldda, A, lda );
    } else {
        magma_tally2blas_ztranspose( n, m, dAT(0,0), lddat, dA(0,0), ldda );
        magma_tally2_zgetmatrix( m, n, dA(0,0), ldda, A, lda );
    }

    magma_tally2_free( dA );

    return *info;
} /* magma_tally2_zgelqf */