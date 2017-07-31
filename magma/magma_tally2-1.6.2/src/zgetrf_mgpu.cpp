/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions normal z -> s d c

*/
#include "common_magma_tally2.h"


/**
    Purpose
    -------
    ZGETRF computes an LU factorization of a general M-by-N matrix A
    using partial pivoting with row interchanges.

    The factorization has the form
        A = P * L * U
    where P is a permutation matrix, L is lower triangular with unit
    diagonal elements (lower trapezoidal if m > n), and U is upper
    triangular (upper trapezoidal if m < n).

    This is the right-looking Level 3 BLAS version of the algorithm.

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
    d_lA    COMPLEX_16 array of pointers on the GPU, dimension (ngpu).
            On entry, the M-by-N matrix A distributed over GPUs
            (d_lA[d] points to the local matrix on d-th GPU).
            It uses 1D block column cyclic format with the block size of nb,
            and each local matrix is stored by column.
            On exit, the factors L and U from the factorization
            A = P*L*U; the unit diagonal elements of L are not stored.

    @param[in]
    ldda     INTEGER
            The leading dimension of the array d_lA.  LDDA >= max(1,M).

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

    @ingroup magma_tally2_zgesv_comp
    ********************************************************************/
extern "C" magma_tally2_int_t
magma_tally2_zgetrf_mgpu(
    magma_tally2_int_t ngpu,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex_ptr d_lA[], magma_tally2_int_t ldda, magma_tally2_int_t *ipiv,
    magma_tally2_int_t *info)
{
    magma_tally2_int_t nb, n_local[Magma_tally2MaxGPUs];
    magma_tally2_int_t maxm;
    magma_tally2_int_t i, j, d, lddat, lddwork;
    magma_tally2DoubleComplex *d_lAT[Magma_tally2MaxGPUs];
    magma_tally2DoubleComplex *d_panel[Magma_tally2MaxGPUs], *work;
    magma_tally2_queue_t streaml[Magma_tally2MaxGPUs][2];

    /* Check arguments */
    *info = 0;
    if (m < 0)
        *info = -2;
    else if (n < 0)
        *info = -3;
    else if (ldda < max(1,m))
        *info = -5;

    if (*info != 0) {
        magma_tally2_xerbla( __func__, -(*info) );
        return *info;
    }

    /* Quick return if possible */
    if (m == 0 || n == 0)
        return *info;

    /* Function Body */
    nb = magma_tally2_get_zgetrf_nb(m);

    if (nb <= 1 || nb >= n) {
        /* Use CPU code. */
        magma_tally2_zmalloc_cpu( &work, m * n );
        if ( work == NULL ) {
            *info = MAGMA_tally2_ERR_HOST_ALLOC;
            return *info;
        }
        magma_tally2_zgetmatrix( m, n, d_lA[0], ldda, work, m );
        lapackf77_zgetrf(&m, &n, work, &m, ipiv, info);
        magma_tally2_zsetmatrix( m, n, work, m, d_lA[0], ldda );
        magma_tally2_free_cpu(work);
    } else {
        /* Use hybrid blocked code. */
        magma_tally2_device_t orig_dev;
        magma_tally2_getdevice( &orig_dev );
        magma_tally2_queue_t orig_stream;
        magma_tally2blasGetKernelStream( &orig_stream );
        
        maxm = ((m + 31)/32)*32;
        if ( ngpu > ceil((double)n/nb) ) {
            printf( " * too many GPUs for the matrix size, using %d GPUs\n", (int) ngpu );
            *info = -1;
            return *info;
        }

        /* allocate workspace for each GPU */
        lddat = ((((((n+nb-1)/nb)/ngpu)*nb)+31)/32)*32;
        lddat = (n+nb-1)/nb;                 /* number of block columns         */
        lddat = (lddat+ngpu-1)/ngpu;         /* number of block columns per GPU */
        lddat = nb*lddat;                    /* number of columns per GPU       */
        lddat = ((lddat+31)/32)*32;          /* make it a multiple of 32        */
        for (i=0; i < ngpu; i++) {
            magma_tally2_setdevice(i);
            
            /* local-n and local-ld */
            n_local[i] = ((n/nb)/ngpu)*nb;
            if (i < (n/nb)%ngpu)
                n_local[i] += nb;
            else if (i == (n/nb)%ngpu)
                n_local[i] += n%nb;
            
            /* workspaces */
            if (MAGMA_tally2_SUCCESS != magma_tally2_zmalloc( &d_panel[i], (3+ngpu)*nb*maxm )) {
                for( j=0; j <= i; j++ ) {
                    magma_tally2_setdevice(j);
                }
                for( j=0; j < i; j++ ) {
                    magma_tally2_setdevice(j);
                    magma_tally2_free( d_panel[j] );
                    magma_tally2_free( d_lAT[j]   );
                }
                *info = MAGMA_tally2_ERR_DEVICE_ALLOC;
                return *info;
            }
            
            /* local-matrix storage */
            if (MAGMA_tally2_SUCCESS != magma_tally2_zmalloc( &d_lAT[i], lddat*maxm )) {
                for( j=0; j <= i; j++ ) {
                    magma_tally2_setdevice(j);
                    magma_tally2_free( d_panel[j] );
                }
                for( j=0; j < i; j++ ) {
                    magma_tally2_setdevice(j);
                    magma_tally2_free( d_lAT[j] );
                }
                *info = MAGMA_tally2_ERR_DEVICE_ALLOC;
                return *info;
            }
            
            /* create the streams */
            magma_tally2_queue_create( &streaml[i][0] );
            magma_tally2_queue_create( &streaml[i][1] );
            
            magma_tally2blasSetKernelStream(streaml[i][1]);
            magma_tally2blas_ztranspose( m, n_local[i], d_lA[i], ldda, d_lAT[i], lddat );
        }
        for (i=0; i < ngpu; i++) {
            magma_tally2_setdevice(i);
            magma_tally2_queue_sync(streaml[i][0]);
            magma_tally2blasSetKernelStream(NULL);
        }
        magma_tally2_setdevice(0);

        /* cpu workspace */
        lddwork = maxm;
        if (MAGMA_tally2_SUCCESS != magma_tally2_zmalloc_pinned( &work, lddwork*nb*ngpu )) {
            for (i=0; i < ngpu; i++ ) {
                magma_tally2_setdevice(i);
                magma_tally2_free( d_panel[i] );
                magma_tally2_free( d_lAT[i]   );
            }
            *info = MAGMA_tally2_ERR_HOST_ALLOC;
            return *info;
        }

        /* calling multi-gpu interface with allocated workspaces and streams */
        magma_tally2_zgetrf2_mgpu(ngpu, m, n, nb, 0, d_lAT, lddat, ipiv, d_panel, work, maxm,
                           streaml, info);

        /* clean up */
        for( d=0; d < ngpu; d++ ) {
            magma_tally2_setdevice(d);
            
            /* save on output */
            magma_tally2blas_ztranspose( n_local[d], m, d_lAT[d], lddat, d_lA[d], ldda );
            magma_tally2_device_sync();
            magma_tally2_free( d_lAT[d]   );
            magma_tally2_free( d_panel[d] );
            magma_tally2_queue_destroy( streaml[d][0] );
            magma_tally2_queue_destroy( streaml[d][1] );
        } /* end of for d=1,..,ngpu */
        magma_tally2_setdevice( orig_dev );
        magma_tally2blasSetKernelStream( orig_stream );
        magma_tally2_free_pinned( work );
    }
        
    return *info;
}
