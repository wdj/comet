/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Stan Tomov
       @generated from zgeqrf2_gpu.cpp normal z -> d, Fri Jan 30 19:00:15 2015

*/
#include "common_magma_tally3.h"

/**
    Purpose
    -------
    DGEQRF computes a QR factorization of a real M-by-N matrix A:
    A = Q * R.
    
    This version has LAPACK-complaint arguments.
    
    If the current stream is NULL, this version replaces it with a new
    stream to overlap computation with communication.

    Other versions (magma_tally3_dgeqrf_gpu and magma_tally3_dgeqrf3_gpu) store the
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
    dA      DOUBLE_PRECISION array on the GPU, dimension (LDDA,N)
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
    tau     DOUBLE_PRECISION array, dimension (min(M,N))
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

    where tau is a real scalar, and v is a real vector with
    v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
    and tau in TAU(i).

    @ingroup magma_tally3_dgeqrf_comp
    ********************************************************************/
extern "C" magma_tally3_int_t
magma_tally3_dgeqrf2_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3Double_ptr dA, magma_tally3_int_t ldda,
    double *tau,
    magma_tally3_int_t *info )
{
    #define dA(a_1,a_2)    ( dA+(a_2)*(ldda) + (a_1))
    #define work_ref(a_1)  ( work + (a_1))
    #define hwork          ( work + (nb)*(m))

    magma_tally3Double_ptr dwork;
    double *work;
    magma_tally3_int_t i, k, ldwork, lddwork, old_i, old_ib, rows;
    magma_tally3_int_t nbmin, nx, ib, nb;
    magma_tally3_int_t lhwork, lwork;

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
        magma_tally3_xerbla( __func__, -(*info) );
        return *info;
    }

    k = min(m,n);
    if (k == 0)
        return *info;

    nb = magma_tally3_get_dgeqrf_nb(m);

    lwork  = (m+n) * nb;
    lhwork = lwork - (m)*nb;

    if (MAGMA_tally3_SUCCESS != magma_tally3_dmalloc( &dwork, n*nb )) {
        *info = MAGMA_tally3_ERR_DEVICE_ALLOC;
        return *info;
    }

    if (MAGMA_tally3_SUCCESS != magma_tally3_dmalloc_pinned( &work, lwork )) {
        magma_tally3_free( dwork );
        *info = MAGMA_tally3_ERR_HOST_ALLOC;
        return *info;
    }

    /* Define user stream if current stream is NULL */
    magma_tally3_queue_t stream[2];
    
    magma_tally3_queue_t orig_stream;
    magma_tally3blasGetKernelStream( &orig_stream );

    magma_tally3_queue_create( &stream[0] );
    if (orig_stream == NULL) {
        magma_tally3_queue_create( &stream[1] );
        magma_tally3blasSetKernelStream(stream[1]);
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
            magma_tally3_queue_sync( stream[1] );
            magma_tally3_dgetmatrix_async( rows, ib,
                                    dA(i,i),       ldda,
                                    work_ref(i), ldwork, stream[0] );
            if (i > 0) {
                /* Apply H' to A(i:m,i+2*ib:n) from the left */
                magma_tally3_dlarfb_gpu( Magma_tally3Left, Magma_tally3ConjTrans, Magma_tally3Forward, Magma_tally3Columnwise,
                                  m-old_i, n-old_i-2*old_ib, old_ib,
                                  dA(old_i, old_i         ), ldda, dwork,        lddwork,
                                  dA(old_i, old_i+2*old_ib), ldda, dwork+old_ib, lddwork);

                magma_tally3_dsetmatrix_async( old_ib, old_ib,
                                        work_ref(old_i),  ldwork,
                                        dA(old_i, old_i), ldda, stream[1] );
            }

            magma_tally3_queue_sync( stream[0] );
            lapackf77_dgeqrf(&rows, &ib, work_ref(i), &ldwork, tau+i, hwork, &lhwork, info);
            /* Form the triangular factor of the block reflector
               H = H(i) H(i+1) . . . H(i+ib-1) */
            lapackf77_dlarft( Magma_tally3ForwardStr, Magma_tally3ColumnwiseStr,
                              &rows, &ib,
                              work_ref(i), &ldwork, tau+i, hwork, &ib);

            dpanel_to_q_tally3( Magma_tally3Upper, ib, work_ref(i), ldwork, hwork+ib*ib );

            /* download the i-th V matrix */
            magma_tally3_dsetmatrix_async( rows, ib, work_ref(i), ldwork, dA(i,i), ldda, stream[0] );

            /* download the T matrix */
            magma_tally3_queue_sync( stream[1] );
            magma_tally3_dsetmatrix_async( ib, ib, hwork, ib, dwork, lddwork, stream[0] );
            magma_tally3_queue_sync( stream[0] );

            if (i + ib < n) {
                if (i+nb < k-nx) {
                    /* Apply H' to A(i:m,i+ib:i+2*ib) from the left */
                    magma_tally3_dlarfb_gpu( Magma_tally3Left, Magma_tally3ConjTrans, Magma_tally3Forward, Magma_tally3Columnwise,
                                      rows, ib, ib,
                                      dA(i, i   ), ldda, dwork,    lddwork,
                                      dA(i, i+ib), ldda, dwork+ib, lddwork);
                    dq_to_panel_tally3( Magma_tally3Upper, ib, work_ref(i), ldwork, hwork+ib*ib );
                }
                else {
                    magma_tally3_dlarfb_gpu( Magma_tally3Left, Magma_tally3ConjTrans, Magma_tally3Forward, Magma_tally3Columnwise,
                                      rows, n-i-ib, ib,
                                      dA(i, i   ), ldda, dwork,    lddwork,
                                      dA(i, i+ib), ldda, dwork+ib, lddwork);
                    dq_to_panel_tally3( Magma_tally3Upper, ib, work_ref(i), ldwork, hwork+ib*ib );
                    magma_tally3_dsetmatrix_async( ib, ib,
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
    magma_tally3_free( dwork );

    /* Use unblocked code to factor the last or only block. */
    if (i < k) {
        ib   = n-i;
        rows = m-i;
        magma_tally3_dgetmatrix_async( rows, ib, dA(i, i), ldda, work, rows, stream[1] );
        magma_tally3_queue_sync( stream[1] );
        lhwork = lwork - rows*ib;
        lapackf77_dgeqrf(&rows, &ib, work, &rows, tau+i, work+ib*rows, &lhwork, info);
        
        magma_tally3_dsetmatrix_async( rows, ib, work, rows, dA(i, i), ldda, stream[1] );
    }

    magma_tally3_free_pinned( work );

    magma_tally3_queue_destroy( stream[0] );
    if (orig_stream == NULL) {
        magma_tally3_queue_destroy( stream[1] );
    }
    magma_tally3blasSetKernelStream( orig_stream );

    return *info;
} /* magma_tally3_dgeqrf2_gpu */
