/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from zpotrf_mgpu.cpp normal z -> c, Fri Jan 30 19:00:13 2015

*/
#include "common_magma_minproduct.h"


/**
    Purpose
    -------
    CPOTRF computes the Cholesky factorization of a complex Hermitian
    positive definite matrix dA.

    The factorization has the form
       dA = U**H * U,   if UPLO = Magma_minproductUpper, or
       dA = L  * L**H,  if UPLO = Magma_minproductLower,
    where U is an upper triangular matrix and L is lower triangular.

    This is the block version of the algorithm, calling Level 3 BLAS.

    Arguments
    ---------
    @param[in]
    ngpu    INTEGER
            Number of GPUs to use. ngpu > 0.

    @param[in]
    uplo    magma_minproduct_uplo_t
      -     = Magma_minproductUpper:  Upper triangle of dA is stored;
      -     = Magma_minproductLower:  Lower triangle of dA is stored.

    @param[in]
    n       INTEGER
            The order of the matrix dA.  N >= 0.

    @param[in,out]
    d_lA    COMPLEX array of pointers on the GPU, dimension (ngpu)
            On entry, the Hermitian matrix dA distributed over GPUs
            (d_lA[d] points to the local matrix on the d-th GPU).
            It is distributed in 1D block column or row cyclic (with the
            block size of nb) if UPLO = Magma_minproductUpper or Magma_minproductLower, respectively.
            If UPLO = Magma_minproductUpper, the leading N-by-N upper triangular
            part of dA contains the upper triangular part of the matrix dA,
            and the strictly lower triangular part of dA is not referenced.
            If UPLO = Magma_minproductLower, the leading N-by-N lower triangular part
            of dA contains the lower triangular part of the matrix dA, and
            the strictly upper triangular part of dA is not referenced.
    \n
            On exit, if INFO = 0, the factor U or L from the Cholesky
            factorization dA = U**H * U or dA = L * L**H.

    @param[in]
    ldda     INTEGER
            The leading dimension of the array d_lA. LDDA >= max(1,N).
            To benefit from coalescent memory accesses LDDA must be
            divisible by 16.

    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument had an illegal value
      -     > 0:  if INFO = i, the leading minor of order i is not
                  positive definite, and the factorization could not be
                  completed.

    @ingroup magma_minproduct_cposv_comp
    ********************************************************************/
extern "C" magma_minproduct_int_t
magma_minproduct_cpotrf_mgpu(
    magma_minproduct_int_t ngpu,
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductFloatComplex_ptr d_lA[], magma_minproduct_int_t ldda,
    magma_minproduct_int_t *info)
{
    magma_minproduct_int_t     j, nb, d, lddp, h;
    const char* uplo_ = lapack_uplo_const( uplo );
    magma_minproductFloatComplex *work;
    int upper = (uplo == Magma_minproductUpper);
    magma_minproductFloatComplex *dwork[Magma_minproductMaxGPUs];
    magma_minproduct_queue_t    stream[Magma_minproductMaxGPUs][3];
    magma_minproduct_event_t     event[Magma_minproductMaxGPUs][5];

    *info = 0;
    nb = magma_minproduct_get_cpotrf_nb(n);
    if (! upper && uplo != Magma_minproductLower) {
        *info = -1;
    } else if (n < 0) {
        *info = -2;
    } else if (!upper) {
        lddp = nb*(n/(nb*ngpu));
        if ( n%(nb*ngpu) != 0 ) lddp += min(nb, n-ngpu*lddp);
        if ( ldda < lddp ) *info = -4;
    } else if ( ldda < n ) {
        *info = -4;
    }
    if (*info != 0) {
        magma_minproduct_xerbla( __func__, -(*info) );
        return *info;
    }

    magma_minproduct_device_t orig_dev;
    magma_minproduct_getdevice( &orig_dev );
    
    if (ngpu == 1 && ((nb <= 1) || (nb >= n)) ) {
        /*  Use unblocked code. */
        magma_minproduct_setdevice(0);
        if (MAGMA_minproduct_SUCCESS != magma_minproduct_cmalloc_pinned( &work, n*nb )) {
            *info = MAGMA_minproduct_ERR_HOST_ALLOC;
            return *info;
        }
        magma_minproduct_cgetmatrix( n, n, d_lA[0], ldda, work, n );
        lapackf77_cpotrf(uplo_, &n, work, &n, info);
        magma_minproduct_csetmatrix( n, n, work, n, d_lA[0], ldda );
        magma_minproduct_free_pinned( work );
    }
    else {
        lddp = nb*((n+nb-1)/nb);
        for( d=0; d < ngpu; d++ ) {
            magma_minproduct_setdevice(d);
            if (MAGMA_minproduct_SUCCESS != magma_minproduct_cmalloc( &dwork[d], ngpu*nb*lddp )) {
                for( j=0; j < d; j++ ) {
                    magma_minproduct_setdevice(j);
                    magma_minproduct_free( dwork[j] );
                }
                *info = MAGMA_minproduct_ERR_DEVICE_ALLOC;
                return *info;
            }
            for( j=0; j < 3; j++ )
                magma_minproduct_queue_create( &stream[d][j] );
            for( j=0; j < 5; j++ )
                magma_minproduct_event_create( &event[d][j]  );
        }
        magma_minproduct_setdevice(0);
        h = 1; //ngpu; //(n+nb-1)/nb;
        if (MAGMA_minproduct_SUCCESS != magma_minproduct_cmalloc_pinned( &work, n*nb*h )) {
            *info = MAGMA_minproduct_ERR_HOST_ALLOC;
            return *info;
        }
        if (upper) {
            /* with three streams */
            magma_minproduct_cpotrf3_mgpu(ngpu, uplo, n, n, 0, 0, nb, d_lA, ldda, dwork, lddp, work, n,
                               h, stream, event, info);
        } else {
            /* with three streams */
            magma_minproduct_cpotrf3_mgpu(ngpu, uplo, n, n, 0, 0, nb, d_lA, ldda, dwork, lddp, work, nb*h,
                               h, stream, event, info);
        }

        /* clean up */
        for( d=0; d < ngpu; d++ ) {
            magma_minproduct_setdevice(d);
            for( j=0; j < 3; j++ ) {
                magma_minproduct_queue_sync( stream[d][j] );
                magma_minproduct_queue_destroy( stream[d][j] );
            }
            
            for( j=0; j < 5; j++ )
                magma_minproduct_event_destroy( event[d][j] );
            
            magma_minproduct_free( dwork[d] );
        }
        magma_minproduct_free_pinned( work );
    } /* end of not lapack */

    magma_minproduct_setdevice( orig_dev );
    
    return *info;
} /* magma_minproduct_cpotrf_mgpu */
