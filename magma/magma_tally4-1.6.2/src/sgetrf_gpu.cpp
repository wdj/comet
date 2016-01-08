/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Stan Tomov
       @generated from zgetrf_gpu.cpp normal z -> s, Fri Jan 30 19:00:14 2015
*/
#include "common_magma_tally4.h"


/**
    Purpose
    -------
    SGETRF computes an LU factorization of a general M-by-N matrix A
    using partial pivoting with row interchanges.

    The factorization has the form
        A = P * L * U
    where P is a permutation matrix, L is lower triangular with unit
    diagonal elements (lower trapezoidal if m > n), and U is upper
    triangular (upper trapezoidal if m < n).

    This is the right-looking Level 3 BLAS version of the algorithm.
    
    If the current stream is NULL, this version replaces it with a new
    stream to overlap computation with communication.

    Arguments
    ---------
    @param[in]
    m       INTEGER
            The number of rows of the matrix A.  M >= 0.

    @param[in]
    n       INTEGER
            The number of columns of the matrix A.  N >= 0.

    @param[in,out]
    dA      REAL array on the GPU, dimension (LDDA,N).
            On entry, the M-by-N matrix to be factored.
            On exit, the factors L and U from the factorization
            A = P*L*U; the unit diagonal elements of L are not stored.

    @param[in]
    ldda     INTEGER
            The leading dimension of the array A.  LDDA >= max(1,M).

    @param[out]
    ipiv    INTEGER array, dimension (min(M,N))
            The pivot indices; for 1 <= i <= min(M,N), row i of the
            matrix was interchanged with row IPIV(i).

    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument had an illegal value
                  or another error occured, such as memory allocation failed.
      -     > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
                  has been completed, but the factor U is exactly
                  singular, and division by zero will occur if it is used
                  to solve a system of equations.

    @ingroup magma_tally4_sgesv_comp
    ********************************************************************/
extern "C" magma_tally4_int_t
magma_tally4_sgetrf_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4Float_ptr dA, magma_tally4_int_t ldda, magma_tally4_int_t *ipiv,
    magma_tally4_int_t *info)
{
    #define dAT(i_, j_) (dAT + (i_)*nb*lddat + (j_)*nb)

    float c_one     = MAGMA_tally4_S_ONE;
    float c_neg_one = MAGMA_tally4_S_NEG_ONE;

    magma_tally4_int_t iinfo, nb;
    magma_tally4_int_t maxm, maxn, mindim;
    magma_tally4_int_t i, j, rows, cols, s, lddat, ldwork;
    float *dAT, *dAP, *work;

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
    nb     = magma_tally4_get_sgetrf_nb(m);
    s      = mindim / nb;

    if (nb <= 1 || nb >= min(m,n)) {
        /* Use CPU code. */
        magma_tally4_smalloc_cpu( &work, m * n );
        if ( work == NULL ) {
            *info = MAGMA_tally4_ERR_HOST_ALLOC;
            return *info;
        }
        magma_tally4_sgetmatrix( m, n, dA, ldda, work, m );
        lapackf77_sgetrf(&m, &n, work, &m, ipiv, info);
        magma_tally4_ssetmatrix( m, n, work, m, dA, ldda );
        magma_tally4_free_cpu(work);
    }
    else {
        /* Use hybrid blocked code. */
        maxm = ((m + 31)/32)*32;
        maxn = ((n + 31)/32)*32;

        if (MAGMA_tally4_SUCCESS != magma_tally4_smalloc( &dAP, nb*maxm )) {
            *info = MAGMA_tally4_ERR_DEVICE_ALLOC;
            return *info;
        }

        // square matrices can be done in place;
        // rectangular requires copy to transpose
        if ( m == n ) {
            dAT = dA;
            lddat = ldda;
            magma_tally4blas_stranspose_inplace( m, dAT, ldda );
        }
        else {
            lddat = maxn;  // N-by-M
            if (MAGMA_tally4_SUCCESS != magma_tally4_smalloc( &dAT, lddat*maxm )) {
                magma_tally4_free( dAP );
                *info = MAGMA_tally4_ERR_DEVICE_ALLOC;
                return *info;
            }
            magma_tally4blas_stranspose( m, n, dA, ldda, dAT, lddat );
        }

        ldwork = maxm;
        if (MAGMA_tally4_SUCCESS != magma_tally4_smalloc_pinned( &work, ldwork*nb )) {
            magma_tally4_free( dAP );
            if ( ! (m == n))
                magma_tally4_free( dAT );
            *info = MAGMA_tally4_ERR_HOST_ALLOC;
            return *info;
        }

        /* Define user stream if current stream is NULL */
        magma_tally4_queue_t stream[2];
        
        magma_tally4_queue_t orig_stream;
        magma_tally4blasGetKernelStream( &orig_stream );

        magma_tally4_queue_create( &stream[0] );
        if (orig_stream == NULL) {
            magma_tally4_queue_create( &stream[1] );
            magma_tally4blasSetKernelStream(stream[1]);
        }
        else {
            stream[1] = orig_stream;
        }
  
        for( j=0; j < s; j++ ) {
            // download j-th panel
            cols = maxm - j*nb;
            magma_tally4blas_stranspose( nb, m-j*nb, dAT(j,j), lddat, dAP, cols );

            // make sure that the transpose has completed
            magma_tally4_queue_sync( stream[1] );
            magma_tally4_sgetmatrix_async( m-j*nb, nb, dAP, cols, work, ldwork,
                                    stream[0]);

            if ( j > 0 ) {
                magma_tally4_strsm( Magma_tally4Right, Magma_tally4Upper, Magma_tally4NoTrans, Magma_tally4Unit,
                             n - (j+1)*nb, nb,
                             c_one, dAT(j-1,j-1), lddat,
                                    dAT(j-1,j+1), lddat );
                magma_tally4_sgemm( Magma_tally4NoTrans, Magma_tally4NoTrans,
                             n-(j+1)*nb, m-j*nb, nb,
                             c_neg_one, dAT(j-1,j+1), lddat,
                                        dAT(j,  j-1), lddat,
                             c_one,     dAT(j,  j+1), lddat );
            }

            // do the cpu part
            rows = m - j*nb;
            magma_tally4_queue_sync( stream[0] );
            lapackf77_sgetrf( &rows, &nb, work, &ldwork, ipiv+j*nb, &iinfo);
            if ( *info == 0 && iinfo > 0 )
                *info = iinfo + j*nb;

            // upload j-th panel
            magma_tally4_ssetmatrix_async( m-j*nb, nb, work, ldwork, dAP, maxm,
                                    stream[0]);

            for( i=j*nb; i < j*nb + nb; ++i ) {
                ipiv[i] += j*nb;
            }
            magma_tally4blas_slaswp( n, dAT, lddat, j*nb + 1, j*nb + nb, ipiv, 1 );

            magma_tally4_queue_sync( stream[0] );
            magma_tally4blas_stranspose( m-j*nb, nb, dAP, maxm, dAT(j,j), lddat );

            // do the small non-parallel computations (next panel update)
            if ( s > (j+1) ) {
                magma_tally4_strsm( Magma_tally4Right, Magma_tally4Upper, Magma_tally4NoTrans, Magma_tally4Unit,
                             nb, nb,
                             c_one, dAT(j, j  ), lddat,
                                    dAT(j, j+1), lddat);
                magma_tally4_sgemm( Magma_tally4NoTrans, Magma_tally4NoTrans,
                             nb, m-(j+1)*nb, nb,
                             c_neg_one, dAT(j,   j+1), lddat,
                                        dAT(j+1, j  ), lddat,
                             c_one,     dAT(j+1, j+1), lddat );
            }
            else {
                magma_tally4_strsm( Magma_tally4Right, Magma_tally4Upper, Magma_tally4NoTrans, Magma_tally4Unit,
                             n-s*nb, nb,
                             c_one, dAT(j, j  ), lddat,
                                    dAT(j, j+1), lddat);
                magma_tally4_sgemm( Magma_tally4NoTrans, Magma_tally4NoTrans,
                             n-(j+1)*nb, m-(j+1)*nb, nb,
                             c_neg_one, dAT(j,   j+1), lddat,
                                        dAT(j+1, j  ), lddat,
                             c_one,     dAT(j+1, j+1), lddat );
            }
        }

        magma_tally4_int_t nb0 = min(m - s*nb, n - s*nb);
        if ( nb0 > 0 ) {
            rows = m - s*nb;
            cols = maxm - s*nb;
    
            magma_tally4blas_stranspose( nb0, rows, dAT(s,s), lddat, dAP, maxm );
            magma_tally4_sgetmatrix( rows, nb0, dAP, maxm, work, ldwork );
    
            // do the cpu part
            lapackf77_sgetrf( &rows, &nb0, work, &ldwork, ipiv+s*nb, &iinfo);
            if ( *info == 0 && iinfo > 0 )
                *info = iinfo + s*nb;
                
            for( i=s*nb; i < s*nb + nb0; ++i ) {
                ipiv[i] += s*nb;
            }
            magma_tally4blas_slaswp( n, dAT, lddat, s*nb + 1, s*nb + nb0, ipiv, 1 );
    
            // upload j-th panel
            magma_tally4_ssetmatrix( rows, nb0, work, ldwork, dAP, maxm );
            magma_tally4blas_stranspose( rows, nb0, dAP, maxm, dAT(s,s), lddat );
    
            magma_tally4_strsm( Magma_tally4Right, Magma_tally4Upper, Magma_tally4NoTrans, Magma_tally4Unit,
                         n-s*nb-nb0, nb0,
                         c_one, dAT(s,s),     lddat,
                                dAT(s,s)+nb0, lddat);
        }
        
        // undo transpose
        if ( m == n ) {
            magma_tally4blas_stranspose_inplace( m, dAT, lddat );
        }
        else {
            magma_tally4blas_stranspose( n, m, dAT, lddat, dA, ldda );
            magma_tally4_free( dAT );
        }

        magma_tally4_free( dAP );
        magma_tally4_free_pinned( work );
        
        magma_tally4_queue_destroy( stream[0] );
        if (orig_stream == NULL) {
            magma_tally4_queue_destroy( stream[1] );
        }
        magma_tally4blasSetKernelStream( orig_stream );
    }
    
    return *info;
} /* magma_tally4_sgetrf_gpu */

#undef dAT
