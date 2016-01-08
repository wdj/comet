/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions normal z -> s d c

*/
#include "common_magma_tally4.h"

/**
    Purpose
    -------
    ZUNGQR generates an M-by-N COMPLEX_16 matrix Q with orthonormal columns,
    which is defined as the first N columns of a product of K elementary
    reflectors of order M

          Q  =  H(1) H(2) . . . H(k)

    as returned by ZGEQRF_GPU.

    Arguments
    ---------
    @param[in]
    m       INTEGER
            The number of rows of the matrix Q. M >= 0.

    @param[in]
    n       INTEGER
            The number of columns of the matrix Q. M >= N >= 0.

    @param[in]
    k       INTEGER
            The number of elementary reflectors whose product defines the
            matrix Q. N >= K >= 0.

    @param[in,out]
    dA      COMPLEX_16 array A on the GPU device,
            dimension (LDDA,N). On entry, the i-th column must contain
            the vector which defines the elementary reflector H(i), for
            i = 1,2,...,k, as returned by ZGEQRF_GPU in the first k
            columns of its array argument A.
            On exit, the M-by-N matrix Q.

    @param[in]
    ldda    INTEGER
            The first dimension of the array A. LDDA >= max(1,M).

    @param[in]
    tau     COMPLEX_16 array, dimension (K)
            TAU(i) must contain the scalar factor of the elementary
            reflector H(i), as returned by ZGEQRF_GPU.

    @param[in]
    dT      COMPLEX_16 work space array on the GPU device,
            dimension (MIN(M, N) )*NB.
            This must be the 6th argument of magma_tally4_zgeqrf_gpu
            [ note that if N here is bigger than N in magma_tally4_zgeqrf_gpu,
              the workspace requirement DT in magma_tally4_zgeqrf_gpu must be
              as specified in this routine ].

    @param[in]
    nb      INTEGER
            This is the block size used in ZGEQRF_GPU, and correspondingly
            the size of the T matrices, used in the factorization, and
            stored in DT.

    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument has an illegal value

    @ingroup magma_tally4_zheev_2stage
    ********************************************************************/
extern "C" magma_tally4_int_t
magma_tally4_zungqr_2stage_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex *tau,
    magma_tally4DoubleComplex_ptr dT, magma_tally4_int_t nb,
    magma_tally4_int_t *info)
{
    #define dA(a_1,a_2) (dA + (a_2)*(ldda) + (a_1))
    #define dT(a_1)     (dT + (a_1)*nb)

    magma_tally4DoubleComplex c_zero = MAGMA_tally4_Z_ZERO;
    magma_tally4DoubleComplex c_one  = MAGMA_tally4_Z_ONE;
    
    magma_tally4_int_t  i__1, i__2, i__3;
    //magma_tally4_int_t lwork;
    magma_tally4_int_t i, ib, ki, kk;  //, iinfo;
    //magma_tally4_int_t lddwork = min(m, n);
    //magma_tally4DoubleComplex *work, *panel;
    magma_tally4DoubleComplex_ptr dwork;
    //magma_tally4_queue_t stream[2];
    magma_tally4_int_t ldt=nb; // need to be an input parameter

    *info = 0;
    if (m < 0) {
        *info = -1;
    } else if ((n < 0) || (n > m)) {
        *info = -2;
    } else if ((k < 0) || (k > n)) {
        *info = -3;
    } else if (ldda < max(1,m)) {
        *info = -5;
    }
    if (*info != 0) {
        magma_tally4_xerbla( __func__, -(*info) );
        return *info;
    }

    if (n <= 0)
        return *info;

    if (MAGMA_tally4_SUCCESS != magma_tally4_zmalloc( &dwork, n*nb )) {
        printf ("!!!! zungqr_2stage magma_tally4_alloc failed for: dwork\n" );
        exit(-1);
    }

    if ( (nb > 1) && (nb < k) ) {
        /*  Use blocked code after the last block.
            The first kk columns are handled by the block method.
            ki is start of 2nd-to-last block. */
        ki = (k - nb - 1) / nb * nb;
        kk = min(k, ki + nb);

        /* Set A(1:kk,kk+1:n) to zero. */
        /* and A(kk+1:m, kk+1:n) = I */
        magma_tally4blas_zlaset( Magma_tally4Full, kk,   n-kk, c_zero, c_zero, dA(0, kk), ldda );
        magma_tally4blas_zlaset( Magma_tally4Full, m-kk, n-kk, c_zero, c_one,  dA(kk,kk), ldda );
    }
    else {
        ki = 0;
        kk = 0;
    }
    
    /* Allocate work space on CPU in pinned memory */
    //lwork = (n+m) * nb;
    //if (kk < n)
    //  lwork = max(lwork, n * nb + (m-kk)*(n-kk));

    //if (MAGMA_tally4_SUCCESS != magma_tally4_zmalloc_pinned( &work, (lwork) )) {
    //    *info = MAGMA_tally4_ERR_HOST_ALLOC;
    //    return *info;
    //}
    //panel = work + n * nb;

    //magma_tally4_queue_create( &stream[0] );
    //magma_tally4_queue_create( &stream[1] );
    /* Use unblocked code for the last or only block. */
    if (kk < n) {
        i__1 = m - kk;
        i__2 = n - kk;
        i__3 = k - kk;
        //magma_tally4_zgetmatrix(i__1, i__2, dA(kk, kk), ldda, panel, i__1);
        //lapackf77_zungqr(&i__1, &i__2, &i__3, panel, &i__1, &tau[kk],
        //                 work, &lwork, &iinfo);
        //
        //magma_tally4_zsetmatrix(i__1, i__2, panel, i__1, dA(kk, kk), ldda);
        
        magma_tally4_zlarfb_gpu( Magma_tally4Left, Magma_tally4NoTrans, Magma_tally4Forward, Magma_tally4Columnwise,
                          i__1, i__2, i__3,
                          dA(kk, kk-nb), ldda, dT(kk-nb), ldt,
                          dA(kk, kk), ldda, dwork, i__2);
        
        //magma_tally4blas_zlaset(Magma_tally4Full, kk-nb,     nb, c_zero, c_zero, dA(0,kk-nb),     ldda);
        //magma_tally4blas_zlaset(Magma_tally4Full, m-(kk-nb), nb, c_zero, c_one,  dA(kk-nb,kk-nb), ldda);
    }

    if (kk > 0) {
        /* Use blocked code */
        for (i = ki; i >= nb; i -= nb) {
            ib = min(nb, k - i);
            /* Send current panel to the CPU for update */
            i__2 = m - i;
            //magma_tally4_zgetmatrix_async( i__2, ib, dA(i,i), ldda, panel, i__2, stream[0] );  // verify
            if (i + ib < n) {
                /* Apply H to A(i:m,i+ib:n) from the left */
                i__3 = n - i;

                magma_tally4blas_zlaset( Magma_tally4Full, i,   ib, c_zero, c_zero, dA(0,i), ldda );
                magma_tally4blas_zlaset( Magma_tally4Full, m-i, ib, c_zero, c_one,  dA(i,i), ldda );

                magma_tally4_zlarfb_gpu( Magma_tally4Left, Magma_tally4NoTrans, Magma_tally4Forward, Magma_tally4Columnwise,
                                  i__2, i__3, ib,
                                  dA(i, i-nb), ldda, dT(i-nb),             ldt,
                                  dA(i, i), ldda, dwork, i__3);
            }

            /* Apply H to rows i:m of current block on the CPU */
            //magma_tally4_queue_sync( stream[0] );
            //lapackf77_zungqr(&i__2, &ib, &ib, panel, &i__2, &tau[i],
            //                 work, &lwork, &iinfo);
            //magma_tally4_zsetmatrix_async( i__2, ib, panel, i__2, dA(i,i), ldda, stream[1] );  // verify

            /* Set rows 1:i-1 of current block to zero */
            i__2 = i + ib;
            //magma_tally4blas_zlaset(Magma_tally4Full, i-ib,     ib, c_zero, c_zero, dA(0,i-ib),    ldda);
            //magma_tally4blas_zlaset(Magma_tally4Full, m-(i-ib), ib, c_zero, c_one,  dA(i-ib,i-ib), ldda);
        }
    }

    magma_tally4blas_zlaset( Magma_tally4Full, m, nb, c_zero, c_one, dA(0,0), ldda );

    magma_tally4_free( dwork );
    //magma_tally4_free_pinned( work );
    //magma_tally4_queue_destroy( stream[0] );
    //magma_tally4_queue_destroy( stream[1] );

    return *info;
} /* magma_tally4_zungqr_gpu */

#undef dA
#undef dT
