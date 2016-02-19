/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Stan Tomov
       @author Raffaele Solca
       @author Ichitaro Yamazaki
       @author Mark Gates

       @generated from zlatrd_mgpu.cpp normal z -> s, Fri Jan 30 19:00:17 2015

*/
#include "common_magma_tally3.h"
#include "trace.h"

#define PRECISION_s
#define REAL

/**
    Purpose
    -------
    SLATRD reduces NB rows and columns of a real symmetric matrix A to
    symmetric tridiagonal form by an orthogonal similarity
    transformation Q' * A * Q, and returns the matrices V and W which are
    needed to apply the transformation to the unreduced part of A.

    If UPLO = Magma_tally3Upper, SLATRD reduces the last NB rows and columns of a
    matrix, of which the upper triangle is supplied;
    if UPLO = Magma_tally3Lower, SLATRD reduces the first NB rows and columns of a
    matrix, of which the lower triangle is supplied.

    This is an auxiliary routine called by SSYTRD.

    Arguments
    ---------
    @param[in]
    uplo    magma_tally3_uplo_t
            Specifies whether the upper or lower triangular part of the
            symmetric matrix A is stored:
      -     = Magma_tally3Upper: Upper triangular
      -     = Magma_tally3Lower: Lower triangular

    @param[in]
    n       INTEGER
            The order of the matrix A.

    @param[in]
    nb      INTEGER
            The number of rows and columns to be reduced.

    @param[in,out]
    A       REAL array, dimension (LDA,N)
            On entry, the symmetric matrix A.  If UPLO = Magma_tally3Upper, the leading
            n-by-n upper triangular part of A contains the upper
            triangular part of the matrix A, and the strictly lower
            triangular part of A is not referenced.  If UPLO = Magma_tally3Lower, the
            leading n-by-n lower triangular part of A contains the lower
            triangular part of the matrix A, and the strictly upper
            triangular part of A is not referenced.
            On exit:
      -     if UPLO = Magma_tally3Upper, the last NB columns have been reduced to
              tridiagonal form, with the diagonal elements overwriting
              the diagonal elements of A; the elements above the diagonal
              with the array TAU, represent the orthogonal matrix Q as a
              product of elementary reflectors;
      -     if UPLO = Magma_tally3Lower, the first NB columns have been reduced to
              tridiagonal form, with the diagonal elements overwriting
              the diagonal elements of A; the elements below the diagonal
              with the array TAU, represent the  orthogonal matrix Q as a
              product of elementary reflectors.
            See Further Details.

    @param[in]
    lda     INTEGER
            The leading dimension of the array A.  LDA >= (1,N).

    @param[out]
    e       REAL array, dimension (N-1)
            If UPLO = Magma_tally3Upper, E(n-nb:n-1) contains the superdiagonal
            elements of the last NB columns of the reduced matrix;
            if UPLO = Magma_tally3Lower, E(1:nb) contains the subdiagonal elements of
            the first NB columns of the reduced matrix.

    @param[out]
    tau     REAL array, dimension (N-1)
            The scalar factors of the elementary reflectors, stored in
            TAU(n-nb:n-1) if UPLO = Magma_tally3Upper, and in TAU(1:nb) if UPLO = Magma_tally3Lower.
            See Further Details.

    @param[out]
    W       REAL array, dimension (LDW,NB)
            The n-by-nb matrix W required to update the unreduced part
            of A.

    @param[in]
    ldw     INTEGER
            The leading dimension of the array W. LDW >= max(1,N).

    @param
    dA

    @param[in]
    ldda

    @param[in]
    offset

    @param
    dW

    @param[in]
    lddw

    @param
    hwork

    @param[in]
    lhwork

    @param
    dwork

    @param[in]
    ldwork
             
    @param[in]
    queues  magma_tally3_queue_t array of dimension (ngpu).
            queues[dev] is an execution queue on GPU dev.
    
    Further Details
    ---------------
    If UPLO = Magma_tally3Upper, the matrix Q is represented as a product of elementary
    reflectors

       Q = H(n) H(n-1) . . . H(n-nb+1).

    Each H(i) has the form

       H(i) = I - tau * v * v'

    where tau is a real scalar, and v is a real vector with
    v(i:n) = 0 and v(i-1) = 1; v(1:i-1) is stored on exit in A(1:i-1,i),
    and tau in TAU(i-1).

    If UPLO = Magma_tally3Lower, the matrix Q is represented as a product of elementary
    reflectors

       Q = H(1) H(2) . . . H(nb).

    Each H(i) has the form

       H(i) = I - tau * v * v'

    where tau is a real scalar, and v is a real vector with
    v(1:i) = 0 and v(i+1) = 1; v(i+1:n) is stored on exit in A(i+1:n,i),
    and tau in TAU(i).

    The elements of the vectors v together form the n-by-nb matrix V
    which is needed, with W, to apply the transformation to the unreduced
    part of the matrix, using a symmetric rank-2k update of the form:
    A := A - V*W' - W*V'.

    The contents of A on exit are illustrated by the following examples
    with n = 5 and nb = 2:

    if UPLO = Magma_tally3Upper:                       if UPLO = Magma_tally3Lower:

      (  a   a   a   v4  v5 )              (  d                  )
      (      a   a   v4  v5 )              (  1   d              )
      (          a   1   v5 )              (  v1  1   a          )
      (              d   1  )              (  v1  v2  a   a      )
      (                  d  )              (  v1  v2  a   a   a  )

    where d denotes a diagonal element of the reduced matrix, a denotes
    an element of the original matrix that is unchanged, and vi denotes
    an element of the vector defining H(i).

    @ingroup magma_tally3_ssyev_aux
    ********************************************************************/
extern "C" magma_tally3_int_t
magma_tally3_slatrd_mgpu(
    magma_tally3_int_t ngpu,
    magma_tally3_uplo_t uplo,
    magma_tally3_int_t n, magma_tally3_int_t nb, magma_tally3_int_t nb0,
    float *A,  magma_tally3_int_t lda,
    float *e, float *tau,
    float *W,          magma_tally3_int_t ldw,
    magma_tally3Float_ptr dA[],    magma_tally3_int_t ldda, magma_tally3_int_t offset,
    magma_tally3Float_ptr dW[],    magma_tally3_int_t lddw,
    float    *hwork,   magma_tally3_int_t lhwork,
    magma_tally3Float_ptr dwork[], magma_tally3_int_t ldwork,
    magma_tally3_queue_t queues[] )
{
#define A(i, j) (A + (j)*lda + (i))
#define W(i, j) (W + (j)*ldw + (i))

#define dA(dev, i, j)  (dA[(dev)] + ((j)+loffset)*ldda + (i) + offset)
#define dW(dev, i, j)  (dW[(dev)] + (j)          *lddw + (i))
#define dW1(dev, i, j) (dW[(dev)] + ((j)+nb)     *lddw + (i))

    const float c_neg_one = MAGMA_tally3_S_NEG_ONE;
    const float c_one     = MAGMA_tally3_S_ONE;
    const float c_zero    = MAGMA_tally3_S_ZERO;
    const magma_tally3_int_t ione = 1;

    float alpha, value;
    magma_tally3_int_t dev;
    magma_tally3_int_t i, n_i, n_i_1, ip1, iw;

    // TODO check arguments
    magma_tally3_int_t info = 0;
    if (n <= 0) {
        return info;
    }
    
    // TODO allocate f in ssytrd and pass into slatrd. (e.g., expand hwork a bit)
    float *f;
    magma_tally3_smalloc_cpu( &f, n );
    if ( f == NULL ) {
        info = MAGMA_tally3_ERR_HOST_ALLOC;
        return info;
    }

    magma_tally3_device_t orig_dev;
    magma_tally3_getdevice( &orig_dev );
    magma_tally3_queue_t orig_stream;
    magma_tally3blasGetKernelStream( &orig_stream );
    
    if (uplo == Magma_tally3Upper) {
        /* Reduce last NB columns of upper triangle */
        for (i = n-1; i >= n - nb; --i) {
            ip1 = i + 1;
            n_i_1 = n - i - 1;
            iw = i - n + nb;
            if (i < n-1) {
                /* Update A(1:i,i) */
                float wii = -conj( *W(i, iw+1) );
                blasf77_saxpy( &ip1, &wii, A(0, i+1), &ione, A(0, i), &ione );

                wii = -conj( *A(i, i+1) );
                blasf77_saxpy( &ip1, &wii, W(0, iw+1), &ione, A(0, i), &ione );
            }
            if (i > 0) {
                /* Generate elementary reflector H(i) to annihilate A(1:i-2,i) */
                alpha = *A(i-1, i);
                lapackf77_slarfg( &i, &alpha, A(0, i), &ione, &tau[i - 1] );

                e[i-1] = MAGMA_tally3_S_REAL( alpha );
                *A(i-1,i) = MAGMA_tally3_S_ONE;
                
                // TODO Previously, this set dx2[dev] = dW1(dev, 0, iw); and used dx2 in ssymv.
                // TODO Now ssymv handles broadcasting x to the GPUs, but data in dW1 is
                // TODO apparently still used in ssytrd_mgpu / ssyr2k_mgpu.
                for( dev=0; dev < ngpu; dev++ ) {
                    magma_tally3_setdevice( dev );
                    magma_tally3_ssetvector_async( n, A(0,i), 1, dW1(dev, 0, iw), 1, queues[dev] );
                }
                magma_tally3blas_ssymv_mgpu(
                    Magma_tally3Upper, i, c_one, dA, ldda, 0,
                    A(0,i), 1, c_zero, W(0, iw), 1,
                    hwork, lhwork, dwork, ldwork, ngpu, nb0, queues );

                if (i < n-1) {
                    blasf77_sgemv( Magma_tally3ConjTransStr, &i, &n_i_1, &c_one,
                                   W(0,   iw+1), &ldw,
                                   A(0,   i),    &ione, &c_zero,
                                   W(i+1, iw),   &ione );
                }

                /* overlap update */
                if ( i < n-1 && i-1 >= n - nb ) {
                    /* Update A(1:i,i) */
                    #ifdef COMPLEX
                    lapackf77_slacgv( &n_i_1, W(i-1, iw+1), &ldw );
                    #endif
                    blasf77_sgemv( "No transpose", &i, &n_i_1, &c_neg_one,
                                   A(0,   i+1),  &lda,
                                   W(i-1, iw+1), &ldw, &c_one,
                                   A(0,   i-1),  &ione );
                    #ifdef COMPLEX
                    lapackf77_slacgv( &n_i_1, W(i-1, iw+1), &ldw );
                    lapackf77_slacgv( &n_i_1, A(i-1, i +1), &lda );
                    #endif
                    blasf77_sgemv( "No transpose", &i, &n_i_1, &c_neg_one,
                                   W(0,   iw+1), &ldw,
                                   A(i-1, i+1),  &lda, &c_one,
                                   A(0,   i-1),  &ione );
                    #ifdef COMPLEX
                    lapackf77_slacgv( &n_i_1, A(i-1, i+1), &lda );
                    #endif
                }

                // synchronize to get ssymv result W(0, iw)
                magma_tally3blas_ssymv_mgpu_sync(
                    Magma_tally3Upper, i, c_one, dA, ldda, 0,
                    A(0,i), 1, c_zero, W(0, iw), 1,
                    hwork, lhwork, dwork, ldwork, ngpu, nb0, queues );

                if (i < n-1) {
                    blasf77_sgemv( "No transpose", &i, &n_i_1, &c_neg_one,
                                   A(0,   i+1), &lda,
                                   W(i+1, iw),  &ione, &c_one,
                                   W(0,   iw),  &ione );

                    blasf77_sgemv( Magma_tally3ConjTransStr, &i, &n_i_1, &c_one,
                                   A(0,   i+1), &lda,
                                   A(0,   i),   &ione, &c_zero,
                                   W(i+1, iw),  &ione );

                    blasf77_sgemv( "No transpose", &i, &n_i_1, &c_neg_one,
                                   W(0,   iw+1), &ldw,
                                   W(i+1, iw),   &ione, &c_one,
                                   W(0,   iw),   &ione );
                }

                blasf77_sscal( &i, &tau[i - 1], W(0, iw), &ione );

                value = magma_tally3_cblas_sdot( i, W(0,iw), ione, A(0,i), ione );
                alpha = tau[i - 1] * -0.5f * value;
                blasf77_saxpy( &i, &alpha, A(0, i), &ione, W(0, iw), &ione );

                for( dev=0; dev < ngpu; dev++ ) {
                    magma_tally3_setdevice( dev );
                    magma_tally3_ssetvector_async( n, W(0,iw), 1, dW(dev, 0, iw), 1, queues[dev] );
                }
            }
        }
    } else {
        /*  Reduce first NB columns of lower triangle */
        for (i = 0; i < nb; ++i) {
            /* Update A(i:n,i) */
            n_i = n - i;
            //idw = ((offset+i)/nb)%ngpu;
            if ( i > 0 ) {
                trace_cpu_start( 0, "gemv", "gemv" );
                float wii = -conj( *W(i, i-1) );
                blasf77_saxpy( &n_i, &wii, A(i, i-1), &ione, A(i, i), &ione );

                wii = -conj( *A(i, i-1) );
                blasf77_saxpy( &n_i, &wii, W(i, i-1), &ione, A(i, i), &ione );
            }

            if (i < n-1) {
                /* Generate elementary reflector H(i) to annihilate A(i+2:n,i) */
                n_i_1 = n - i - 1;
                trace_cpu_start( 0, "larfg", "larfg" );
                alpha = *A(i+1, i);
                lapackf77_slarfg( &n_i_1, &alpha, A(min(i+2,n-1), i), &ione, &tau[i] );
                e[i] = MAGMA_tally3_S_REAL( alpha );
                *A(i+1,i) = MAGMA_tally3_S_ONE;
                trace_cpu_end( 0 );

                /* Compute W(i+1:n,i) */
                // TODO Previously, this set dx2[id] = dW1(id, 0, i)-offset; and used dx2 in ssymv.
                // TODO Now ssymv handles broadcasting x to the GPUs, but data in dW1 is
                // TODO apparently still used in ssytrd_mgpu / ssyr2k_mgpu.
                for( dev=0; dev < ngpu; dev++ ) {
                    magma_tally3_setdevice( dev );
                    magma_tally3_ssetvector_async( n, A(0,i), 1, dW1(dev, 0, i), 1, queues[dev] );
                }
                
                magma_tally3blas_ssymv_mgpu(
                    Magma_tally3Lower, n_i_1, c_one, dA, ldda, offset+i+1,
                    A(i+1, i), 1, c_zero, W(i+1, i), 1,
                    hwork, lhwork, dwork, ldwork, ngpu, nb0, queues );
                
                trace_cpu_start( 0, "gemv", "gemv" );
                blasf77_sgemv( Magma_tally3ConjTransStr, &n_i_1, &i, &c_one,
                               W(i+1, 0), &ldw,
                               A(i+1, i), &ione, &c_zero,
                               W(0,   i), &ione );
                
                blasf77_sgemv( "No transpose", &n_i_1, &i, &c_neg_one,
                               A(i+1, 0), &lda,
                               W(0,   i), &ione, &c_zero,
                               f,         &ione );
                
                blasf77_sgemv( Magma_tally3ConjTransStr, &n_i_1, &i, &c_one,
                               A(i+1, 0), &lda,
                               A(i+1, i), &ione, &c_zero,
                               W(0,   i), &ione );
                trace_cpu_end( 0 );

                /* overlap update */
                if ( i > 0 && i+1 < n ) {
                    trace_cpu_start( 0, "gemv", "gemv" );
                    #ifdef COMPLEX
                    lapackf77_slacgv( &i, W(i+1, 0), &ldw );
                    #endif
                    blasf77_sgemv( "No transpose", &n_i_1, &i, &c_neg_one,
                                   A(i+1, 0),   &lda,
                                   W(i+1, 0),   &ldw, &c_one,
                                   A(i+1, i+1), &ione );
                    #ifdef COMPLEX
                    lapackf77_slacgv( &i, W(i+1, 0), &ldw );
                    lapackf77_slacgv( &i, A(i+1, 0), &lda );
                    #endif
                    blasf77_sgemv( "No transpose", &n_i_1, &i, &c_neg_one,
                                   W(i+1, 0),   &ldw,
                                   A(i+1, 0),   &lda, &c_one,
                                   A(i+1, i+1), &ione );
                    #ifdef COMPLEX
                    lapackf77_slacgv( &i, A(i+1, 0), &lda );
                    #endif
                    trace_cpu_end( 0 );
                }

                // synchronize to get ssymv result W(i+1, i)
                magma_tally3blas_ssymv_mgpu_sync(
                    Magma_tally3Lower, n_i_1, c_one, dA, ldda, offset+i+1,
                    A(i+1, i), 1, c_zero, W(i+1, i), 1,
                    hwork, lhwork, dwork, ldwork, ngpu, nb0, queues );
                
                trace_cpu_start( 0, "axpy", "axpy" );
                if (i != 0) {
                    blasf77_saxpy( &n_i_1, &c_one, f, &ione, W(i+1, i), &ione );
                }

                blasf77_sgemv( "No transpose", &n_i_1, &i, &c_neg_one,
                               W(i+1, 0), &ldw,
                               W(0,   i), &ione, &c_one,
                               W(i+1, i), &ione );
                blasf77_sscal( &n_i_1, &tau[i], W(i+1,i), &ione );

                value = magma_tally3_cblas_sdot( n_i_1, W(i+1,i), ione, A(i+1,i), ione );
                alpha = tau[i] * -0.5f * value;
                blasf77_saxpy( &n_i_1, &alpha, A(i+1, i), &ione, W(i+1,i), &ione );
                trace_cpu_end( 0 );
                for( dev=0; dev < ngpu; dev++ ) {
                    magma_tally3_setdevice( dev );
                    magma_tally3_ssetvector_async( n, W(0,i), 1, dW(dev, 0, i), 1, queues[dev] );
                }
            }
        }
    }

    magma_tally3_free_cpu( f );

    magma_tally3_setdevice( orig_dev );
    magma_tally3blasSetKernelStream( orig_stream );
    
    return info;
} /* magma_tally3_slatrd_mgpu */

#undef A
#undef W
#undef dA
#undef dW
#undef dW1
