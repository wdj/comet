/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Stan Tomov
       @author Raffaele Solca
       @author Ichitaro Yamazaki
       @author Mark Gates

       @generated from zlatrd_mgpu.cpp normal z -> d, Fri Jan 30 19:00:17 2015

*/
#include "common_magma_minproduct.h"
#include "trace.h"

#define PRECISION_d
#define REAL

/**
    Purpose
    -------
    DLATRD reduces NB rows and columns of a real symmetric matrix A to
    symmetric tridiagonal form by an orthogonal similarity
    transformation Q' * A * Q, and returns the matrices V and W which are
    needed to apply the transformation to the unreduced part of A.

    If UPLO = Magma_minproductUpper, DLATRD reduces the last NB rows and columns of a
    matrix, of which the upper triangle is supplied;
    if UPLO = Magma_minproductLower, DLATRD reduces the first NB rows and columns of a
    matrix, of which the lower triangle is supplied.

    This is an auxiliary routine called by DSYTRD.

    Arguments
    ---------
    @param[in]
    uplo    magma_minproduct_uplo_t
            Specifies whether the upper or lower triangular part of the
            symmetric matrix A is stored:
      -     = Magma_minproductUpper: Upper triangular
      -     = Magma_minproductLower: Lower triangular

    @param[in]
    n       INTEGER
            The order of the matrix A.

    @param[in]
    nb      INTEGER
            The number of rows and columns to be reduced.

    @param[in,out]
    A       DOUBLE_PRECISION array, dimension (LDA,N)
            On entry, the symmetric matrix A.  If UPLO = Magma_minproductUpper, the leading
            n-by-n upper triangular part of A contains the upper
            triangular part of the matrix A, and the strictly lower
            triangular part of A is not referenced.  If UPLO = Magma_minproductLower, the
            leading n-by-n lower triangular part of A contains the lower
            triangular part of the matrix A, and the strictly upper
            triangular part of A is not referenced.
            On exit:
      -     if UPLO = Magma_minproductUpper, the last NB columns have been reduced to
              tridiagonal form, with the diagonal elements overwriting
              the diagonal elements of A; the elements above the diagonal
              with the array TAU, represent the orthogonal matrix Q as a
              product of elementary reflectors;
      -     if UPLO = Magma_minproductLower, the first NB columns have been reduced to
              tridiagonal form, with the diagonal elements overwriting
              the diagonal elements of A; the elements below the diagonal
              with the array TAU, represent the  orthogonal matrix Q as a
              product of elementary reflectors.
            See Further Details.

    @param[in]
    lda     INTEGER
            The leading dimension of the array A.  LDA >= (1,N).

    @param[out]
    e       DOUBLE_PRECISION array, dimension (N-1)
            If UPLO = Magma_minproductUpper, E(n-nb:n-1) contains the superdiagonal
            elements of the last NB columns of the reduced matrix;
            if UPLO = Magma_minproductLower, E(1:nb) contains the subdiagonal elements of
            the first NB columns of the reduced matrix.

    @param[out]
    tau     DOUBLE_PRECISION array, dimension (N-1)
            The scalar factors of the elementary reflectors, stored in
            TAU(n-nb:n-1) if UPLO = Magma_minproductUpper, and in TAU(1:nb) if UPLO = Magma_minproductLower.
            See Further Details.

    @param[out]
    W       DOUBLE_PRECISION array, dimension (LDW,NB)
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
    queues  magma_minproduct_queue_t array of dimension (ngpu).
            queues[dev] is an execution queue on GPU dev.
    
    Further Details
    ---------------
    If UPLO = Magma_minproductUpper, the matrix Q is represented as a product of elementary
    reflectors

       Q = H(n) H(n-1) . . . H(n-nb+1).

    Each H(i) has the form

       H(i) = I - tau * v * v'

    where tau is a real scalar, and v is a real vector with
    v(i:n) = 0 and v(i-1) = 1; v(1:i-1) is stored on exit in A(1:i-1,i),
    and tau in TAU(i-1).

    If UPLO = Magma_minproductLower, the matrix Q is represented as a product of elementary
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

    if UPLO = Magma_minproductUpper:                       if UPLO = Magma_minproductLower:

      (  a   a   a   v4  v5 )              (  d                  )
      (      a   a   v4  v5 )              (  1   d              )
      (          a   1   v5 )              (  v1  1   a          )
      (              d   1  )              (  v1  v2  a   a      )
      (                  d  )              (  v1  v2  a   a   a  )

    where d denotes a diagonal element of the reduced matrix, a denotes
    an element of the original matrix that is unchanged, and vi denotes
    an element of the vector defining H(i).

    @ingroup magma_minproduct_dsyev_aux
    ********************************************************************/
extern "C" magma_minproduct_int_t
magma_minproduct_dlatrd_mgpu(
    magma_minproduct_int_t ngpu,
    magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t n, magma_minproduct_int_t nb, magma_minproduct_int_t nb0,
    double *A,  magma_minproduct_int_t lda,
    double *e, double *tau,
    double *W,          magma_minproduct_int_t ldw,
    magma_minproductDouble_ptr dA[],    magma_minproduct_int_t ldda, magma_minproduct_int_t offset,
    magma_minproductDouble_ptr dW[],    magma_minproduct_int_t lddw,
    double    *hwork,   magma_minproduct_int_t lhwork,
    magma_minproductDouble_ptr dwork[], magma_minproduct_int_t ldwork,
    magma_minproduct_queue_t queues[] )
{
#define A(i, j) (A + (j)*lda + (i))
#define W(i, j) (W + (j)*ldw + (i))

#define dA(dev, i, j)  (dA[(dev)] + ((j)+loffset)*ldda + (i) + offset)
#define dW(dev, i, j)  (dW[(dev)] + (j)          *lddw + (i))
#define dW1(dev, i, j) (dW[(dev)] + ((j)+nb)     *lddw + (i))

    const double c_neg_one = MAGMA_minproduct_D_NEG_ONE;
    const double c_one     = MAGMA_minproduct_D_ONE;
    const double c_zero    = MAGMA_minproduct_D_ZERO;
    const magma_minproduct_int_t ione = 1;

    double alpha, value;
    magma_minproduct_int_t dev;
    magma_minproduct_int_t i, n_i, n_i_1, ip1, iw;

    // TODO check arguments
    magma_minproduct_int_t info = 0;
    if (n <= 0) {
        return info;
    }
    
    // TODO allocate f in dsytrd and pass into dlatrd. (e.g., expand hwork a bit)
    double *f;
    magma_minproduct_dmalloc_cpu( &f, n );
    if ( f == NULL ) {
        info = MAGMA_minproduct_ERR_HOST_ALLOC;
        return info;
    }

    magma_minproduct_device_t orig_dev;
    magma_minproduct_getdevice( &orig_dev );
    magma_minproduct_queue_t orig_stream;
    magma_minproductblasGetKernelStream( &orig_stream );
    
    if (uplo == Magma_minproductUpper) {
        /* Reduce last NB columns of upper triangle */
        for (i = n-1; i >= n - nb; --i) {
            ip1 = i + 1;
            n_i_1 = n - i - 1;
            iw = i - n + nb;
            if (i < n-1) {
                /* Update A(1:i,i) */
                double wii = -conj( *W(i, iw+1) );
                blasf77_daxpy( &ip1, &wii, A(0, i+1), &ione, A(0, i), &ione );

                wii = -conj( *A(i, i+1) );
                blasf77_daxpy( &ip1, &wii, W(0, iw+1), &ione, A(0, i), &ione );
            }
            if (i > 0) {
                /* Generate elementary reflector H(i) to annihilate A(1:i-2,i) */
                alpha = *A(i-1, i);
                lapackf77_dlarfg( &i, &alpha, A(0, i), &ione, &tau[i - 1] );

                e[i-1] = MAGMA_minproduct_D_REAL( alpha );
                *A(i-1,i) = MAGMA_minproduct_D_ONE;
                
                // TODO Previously, this set dx2[dev] = dW1(dev, 0, iw); and used dx2 in dsymv.
                // TODO Now dsymv handles broadcasting x to the GPUs, but data in dW1 is
                // TODO apparently still used in dsytrd_mgpu / dsyr2k_mgpu.
                for( dev=0; dev < ngpu; dev++ ) {
                    magma_minproduct_setdevice( dev );
                    magma_minproduct_dsetvector_async( n, A(0,i), 1, dW1(dev, 0, iw), 1, queues[dev] );
                }
                magma_minproductblas_dsymv_mgpu(
                    Magma_minproductUpper, i, c_one, dA, ldda, 0,
                    A(0,i), 1, c_zero, W(0, iw), 1,
                    hwork, lhwork, dwork, ldwork, ngpu, nb0, queues );

                if (i < n-1) {
                    blasf77_dgemv( Magma_minproductConjTransStr, &i, &n_i_1, &c_one,
                                   W(0,   iw+1), &ldw,
                                   A(0,   i),    &ione, &c_zero,
                                   W(i+1, iw),   &ione );
                }

                /* overlap update */
                if ( i < n-1 && i-1 >= n - nb ) {
                    /* Update A(1:i,i) */
                    #ifdef COMPLEX
                    lapackf77_dlacgv( &n_i_1, W(i-1, iw+1), &ldw );
                    #endif
                    blasf77_dgemv( "No transpose", &i, &n_i_1, &c_neg_one,
                                   A(0,   i+1),  &lda,
                                   W(i-1, iw+1), &ldw, &c_one,
                                   A(0,   i-1),  &ione );
                    #ifdef COMPLEX
                    lapackf77_dlacgv( &n_i_1, W(i-1, iw+1), &ldw );
                    lapackf77_dlacgv( &n_i_1, A(i-1, i +1), &lda );
                    #endif
                    blasf77_dgemv( "No transpose", &i, &n_i_1, &c_neg_one,
                                   W(0,   iw+1), &ldw,
                                   A(i-1, i+1),  &lda, &c_one,
                                   A(0,   i-1),  &ione );
                    #ifdef COMPLEX
                    lapackf77_dlacgv( &n_i_1, A(i-1, i+1), &lda );
                    #endif
                }

                // synchronize to get dsymv result W(0, iw)
                magma_minproductblas_dsymv_mgpu_sync(
                    Magma_minproductUpper, i, c_one, dA, ldda, 0,
                    A(0,i), 1, c_zero, W(0, iw), 1,
                    hwork, lhwork, dwork, ldwork, ngpu, nb0, queues );

                if (i < n-1) {
                    blasf77_dgemv( "No transpose", &i, &n_i_1, &c_neg_one,
                                   A(0,   i+1), &lda,
                                   W(i+1, iw),  &ione, &c_one,
                                   W(0,   iw),  &ione );

                    blasf77_dgemv( Magma_minproductConjTransStr, &i, &n_i_1, &c_one,
                                   A(0,   i+1), &lda,
                                   A(0,   i),   &ione, &c_zero,
                                   W(i+1, iw),  &ione );

                    blasf77_dgemv( "No transpose", &i, &n_i_1, &c_neg_one,
                                   W(0,   iw+1), &ldw,
                                   W(i+1, iw),   &ione, &c_one,
                                   W(0,   iw),   &ione );
                }

                blasf77_dscal( &i, &tau[i - 1], W(0, iw), &ione );

                value = magma_minproduct_cblas_ddot( i, W(0,iw), ione, A(0,i), ione );
                alpha = tau[i - 1] * -0.5f * value;
                blasf77_daxpy( &i, &alpha, A(0, i), &ione, W(0, iw), &ione );

                for( dev=0; dev < ngpu; dev++ ) {
                    magma_minproduct_setdevice( dev );
                    magma_minproduct_dsetvector_async( n, W(0,iw), 1, dW(dev, 0, iw), 1, queues[dev] );
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
                double wii = -conj( *W(i, i-1) );
                blasf77_daxpy( &n_i, &wii, A(i, i-1), &ione, A(i, i), &ione );

                wii = -conj( *A(i, i-1) );
                blasf77_daxpy( &n_i, &wii, W(i, i-1), &ione, A(i, i), &ione );
            }

            if (i < n-1) {
                /* Generate elementary reflector H(i) to annihilate A(i+2:n,i) */
                n_i_1 = n - i - 1;
                trace_cpu_start( 0, "larfg", "larfg" );
                alpha = *A(i+1, i);
                lapackf77_dlarfg( &n_i_1, &alpha, A(min(i+2,n-1), i), &ione, &tau[i] );
                e[i] = MAGMA_minproduct_D_REAL( alpha );
                *A(i+1,i) = MAGMA_minproduct_D_ONE;
                trace_cpu_end( 0 );

                /* Compute W(i+1:n,i) */
                // TODO Previously, this set dx2[id] = dW1(id, 0, i)-offset; and used dx2 in dsymv.
                // TODO Now dsymv handles broadcasting x to the GPUs, but data in dW1 is
                // TODO apparently still used in dsytrd_mgpu / dsyr2k_mgpu.
                for( dev=0; dev < ngpu; dev++ ) {
                    magma_minproduct_setdevice( dev );
                    magma_minproduct_dsetvector_async( n, A(0,i), 1, dW1(dev, 0, i), 1, queues[dev] );
                }
                
                magma_minproductblas_dsymv_mgpu(
                    Magma_minproductLower, n_i_1, c_one, dA, ldda, offset+i+1,
                    A(i+1, i), 1, c_zero, W(i+1, i), 1,
                    hwork, lhwork, dwork, ldwork, ngpu, nb0, queues );
                
                trace_cpu_start( 0, "gemv", "gemv" );
                blasf77_dgemv( Magma_minproductConjTransStr, &n_i_1, &i, &c_one,
                               W(i+1, 0), &ldw,
                               A(i+1, i), &ione, &c_zero,
                               W(0,   i), &ione );
                
                blasf77_dgemv( "No transpose", &n_i_1, &i, &c_neg_one,
                               A(i+1, 0), &lda,
                               W(0,   i), &ione, &c_zero,
                               f,         &ione );
                
                blasf77_dgemv( Magma_minproductConjTransStr, &n_i_1, &i, &c_one,
                               A(i+1, 0), &lda,
                               A(i+1, i), &ione, &c_zero,
                               W(0,   i), &ione );
                trace_cpu_end( 0 );

                /* overlap update */
                if ( i > 0 && i+1 < n ) {
                    trace_cpu_start( 0, "gemv", "gemv" );
                    #ifdef COMPLEX
                    lapackf77_dlacgv( &i, W(i+1, 0), &ldw );
                    #endif
                    blasf77_dgemv( "No transpose", &n_i_1, &i, &c_neg_one,
                                   A(i+1, 0),   &lda,
                                   W(i+1, 0),   &ldw, &c_one,
                                   A(i+1, i+1), &ione );
                    #ifdef COMPLEX
                    lapackf77_dlacgv( &i, W(i+1, 0), &ldw );
                    lapackf77_dlacgv( &i, A(i+1, 0), &lda );
                    #endif
                    blasf77_dgemv( "No transpose", &n_i_1, &i, &c_neg_one,
                                   W(i+1, 0),   &ldw,
                                   A(i+1, 0),   &lda, &c_one,
                                   A(i+1, i+1), &ione );
                    #ifdef COMPLEX
                    lapackf77_dlacgv( &i, A(i+1, 0), &lda );
                    #endif
                    trace_cpu_end( 0 );
                }

                // synchronize to get dsymv result W(i+1, i)
                magma_minproductblas_dsymv_mgpu_sync(
                    Magma_minproductLower, n_i_1, c_one, dA, ldda, offset+i+1,
                    A(i+1, i), 1, c_zero, W(i+1, i), 1,
                    hwork, lhwork, dwork, ldwork, ngpu, nb0, queues );
                
                trace_cpu_start( 0, "axpy", "axpy" );
                if (i != 0) {
                    blasf77_daxpy( &n_i_1, &c_one, f, &ione, W(i+1, i), &ione );
                }

                blasf77_dgemv( "No transpose", &n_i_1, &i, &c_neg_one,
                               W(i+1, 0), &ldw,
                               W(0,   i), &ione, &c_one,
                               W(i+1, i), &ione );
                blasf77_dscal( &n_i_1, &tau[i], W(i+1,i), &ione );

                value = magma_minproduct_cblas_ddot( n_i_1, W(i+1,i), ione, A(i+1,i), ione );
                alpha = tau[i] * -0.5f * value;
                blasf77_daxpy( &n_i_1, &alpha, A(i+1, i), &ione, W(i+1,i), &ione );
                trace_cpu_end( 0 );
                for( dev=0; dev < ngpu; dev++ ) {
                    magma_minproduct_setdevice( dev );
                    magma_minproduct_dsetvector_async( n, W(0,i), 1, dW(dev, 0, i), 1, queues[dev] );
                }
            }
        }
    }

    magma_minproduct_free_cpu( f );

    magma_minproduct_setdevice( orig_dev );
    magma_minproductblasSetKernelStream( orig_stream );
    
    return info;
} /* magma_minproduct_dlatrd_mgpu */

#undef A
#undef W
#undef dA
#undef dW
#undef dW1
