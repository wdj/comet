/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Stan Tomov
       @precisions normal z -> s d c

*/
#include "common_magma_tally2.h"

/**
    Purpose
    -------
    ZGEQRF computes a QR factorization of a complex M-by-N matrix A:
    A = Q * R.
    
    This version has LAPACK-complaint arguments.
    
    If the current stream is NULL, this version replaces it with a new
    stream to overlap computation with communication.

    Other versions (magma_tally2_zgeqrf_gpu and magma_tally2_zgeqrf3_gpu) store the
    intermediate T matrices.

    Arguments
    ---------
    @param[in]
    m       INTEGER
            The number of rows of the matrix A.  M >= 0.

    @param[in]
    n       INTEGER
            The number of columns of the matrix A.  N >= 0.

    @param[in,out]
    dA      COMPLEX_16 array on the GPU, dimension (LDDA,N)
            On entry, the M-by-N matrix A.
            On exit, the elements on and above the diagonal of the array
            contain the min(M,N)-by-N upper trapezoidal matrix R (R is
            upper triangular if m >= n); the elements below the diagonal,
            with the array TAU, represent the orthogonal matrix Q as a
            product of min(m,n) elementary reflectors (see Further
            Details).

    @param[in]
    ldda    INTEGER
            The leading dimension of the array dA.  LDDA >= max(1,M).
            To benefit from coalescent memory accesses LDDA must be
            divisible by 16.

    @param[out]
    tau     COMPLEX_16 array, dimension (min(M,N))
            The scalar factors of the elementary reflectors (see Further
            Details).

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

    where tau is a complex scalar, and v is a complex vector with
    v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
    and tau in TAU(i).

    @ingroup magma_tally2_zgeqrf_comp
    ********************************************************************/
extern "C" magma_tally2_int_t
magma_tally2_zgeqrf2_gpu(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda,
    magma_tally2DoubleComplex *tau,
    magma_tally2_int_t *info )
{
    #define dA(a_1,a_2)    ( dA+(a_2)*(ldda) + (a_1))
    #define work_ref(a_1)  ( work + (a_1))
    #define hwork          ( work + (nb)*(m))

    magma_tally2DoubleComplex_ptr dwork;
    magma_tally2DoubleComplex *work;
    magma_tally2_int_t i, k, ldwork, lddwork, old_i, old_ib, rows;
    magma_tally2_int_t nbmin, nx, ib, nb;
    magma_tally2_int_t lhwork, lwork;

    /* Function Body */
    *info = 0;
    if (m < 0) {
        *info = -1;
    } else if (n < 0) {
        *info = -2;
    } else if (ldda < max(1,m)) {
        *info = -4;
    }
    if (*info != 0) {
        magma_tally2_xerbla( __func__, -(*info) );
        return *info;
    }

    k = min(m,n);
    if (k == 0)
        return *info;

    nb = magma_tally2_get_zgeqrf_nb(m);

    lwork  = (m+n) * nb;
    lhwork = lwork - (m)*nb;

    if (MAGMA_tally2_SUCCESS != magma_tally2_zmalloc( &dwork, n*nb )) {
        *info = MAGMA_tally2_ERR_DEVICE_ALLOC;
        return *info;
    }

    if (MAGMA_tally2_SUCCESS != magma_tally2_zmalloc_pinned( &work, lwork )) {
        magma_tally2_free( dwork );
        *info = MAGMA_tally2_ERR_HOST_ALLOC;
        return *info;
    }

    /* Define user stream if current stream is NULL */
    magma_tally2_queue_t stream[2];
    
    magma_tally2_queue_t orig_stream;
    magma_tally2blasGetKernelStream( &orig_stream );

    magma_tally2_queue_create( &stream[0] );
    if (orig_stream == NULL) {
        magma_tally2_queue_create( &stream[1] );
        magma_tally2blasSetKernelStream(stream[1]);
    }
    else {
        stream[1] = orig_stream;
    }

    nbmin = 2;
    nx    = nb;
    ldwork = m;
    lddwork= n;

    if (nb >= nbmin && nb < k && nx < k) {
        /* Use blocked code initially */
        old_i = 0; old_ib = nb;
        for (i = 0; i < k-nx; i += nb) {
            ib = min(k-i, nb);
            rows = m -i;

            /* download i-th panel */
            magma_tally2_queue_sync( stream[1] );
            magma_tally2_zgetmatrix_async( rows, ib,
                                    dA(i,i),       ldda,
                                    work_ref(i), ldwork, stream[0] );
            if (i > 0) {
                /* Apply H' to A(i:m,i+2*ib:n) from the left */
                magma_tally2_zlarfb_gpu( Magma_tally2Left, Magma_tally2ConjTrans, Magma_tally2Forward, Magma_tally2Columnwise,
                                  m-old_i, n-old_i-2*old_ib, old_ib,
                                  dA(old_i, old_i         ), ldda, dwork,        lddwork,
                                  dA(old_i, old_i+2*old_ib), ldda, dwork+old_ib, lddwork);

                magma_tally2_zsetmatrix_async( old_ib, old_ib,
                                        work_ref(old_i),  ldwork,
                                        dA(old_i, old_i), ldda, stream[1] );
            }

            magma_tally2_queue_sync( stream[0] );
            lapackf77_zgeqrf(&rows, &ib, work_ref(i), &ldwork, tau+i, hwork, &lhwork, info);
            /* Form the triangular factor of the block reflector
               H = H(i) H(i+1) . . . H(i+ib-1) */
            lapackf77_zlarft( Magma_tally2ForwardStr, Magma_tally2ColumnwiseStr,
                              &rows, &ib,
                              work_ref(i), &ldwork, tau+i, hwork, &ib);

            zpanel_to_q_tally2( Magma_tally2Upper, ib, work_ref(i), ldwork, hwork+ib*ib );

            /* download the i-th V matrix */
            magma_tally2_zsetmatrix_async( rows, ib, work_ref(i), ldwork, dA(i,i), ldda, stream[0] );

            /* download the T matrix */
            magma_tally2_queue_sync( stream[1] );
            magma_tally2_zsetmatrix_async( ib, ib, hwork, ib, dwork, lddwork, stream[0] );
            magma_tally2_queue_sync( stream[0] );

            if (i + ib < n) {
                if (i+nb < k-nx) {
                    /* Apply H' to A(i:m,i+ib:i+2*ib) from the left */
                    magma_tally2_zlarfb_gpu( Magma_tally2Left, Magma_tally2ConjTrans, Magma_tally2Forward, Magma_tally2Columnwise,
                                      rows, ib, ib,
                                      dA(i, i   ), ldda, dwork,    lddwork,
                                      dA(i, i+ib), ldda, dwork+ib, lddwork);
                    zq_to_panel_tally2( Magma_tally2Upper, ib, work_ref(i), ldwork, hwork+ib*ib );
                }
                else {
                    magma_tally2_zlarfb_gpu( Magma_tally2Left, Magma_tally2ConjTrans, Magma_tally2Forward, Magma_tally2Columnwise,
                                      rows, n-i-ib, ib,
                                      dA(i, i   ), ldda, dwork,    lddwork,
                                      dA(i, i+ib), ldda, dwork+ib, lddwork);
                    zq_to_panel_tally2( Magma_tally2Upper, ib, work_ref(i), ldwork, hwork+ib*ib );
                    magma_tally2_zsetmatrix_async( ib, ib,
                                            work_ref(i), ldwork,
                                            dA(i,i),     ldda, stream[1] );
                }
                old_i  = i;
                old_ib = ib;
            }
        }
    } else {
        i = 0;
    }
    magma_tally2_free( dwork );

    /* Use unblocked code to factor the last or only block. */
    if (i < k) {
        ib   = n-i;
        rows = m-i;
        magma_tally2_zgetmatrix_async( rows, ib, dA(i, i), ldda, work, rows, stream[1] );
        magma_tally2_queue_sync( stream[1] );
        lhwork = lwork - rows*ib;
        lapackf77_zgeqrf(&rows, &ib, work, &rows, tau+i, work+ib*rows, &lhwork, info);
        
        magma_tally2_zsetmatrix_async( rows, ib, work, rows, dA(i, i), ldda, stream[1] );
    }

    magma_tally2_free_pinned( work );

    magma_tally2_queue_destroy( stream[0] );
    if (orig_stream == NULL) {
        magma_tally2_queue_destroy( stream[1] );
    }
    magma_tally2blasSetKernelStream( orig_stream );

    return *info;
} /* magma_tally2_zgeqrf2_gpu */
