/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Raffaele Solca
       @author Stan Tomov
       @author Mark Gates

       @generated from zhetrd_mgpu.cpp normal z -> c, Fri Jan 30 19:00:18 2015

*/
#include "common_magma_tally3.h"
#include "trace.h"

/**
    Purpose
    -------
    CHETRD reduces a complex Hermitian matrix A to real symmetric
    tridiagonal form T by an orthogonal similarity transformation:
    Q**H * A * Q = T.

    Arguments
    ---------
    @param[in]
    ngpu    INTEGER
            Number of GPUs to use. ngpu > 0.

    @param[in]
    nqueue  INTEGER
            The number of GPU streams used for update.  10 >= nqueue > 0.

    @param[in]
    uplo    magma_tally3_uplo_t
      -     = Magma_tally3Upper:  Upper triangle of A is stored;
      -     = Magma_tally3Lower:  Lower triangle of A is stored.

    @param[in]
    n       INTEGER
            The order of the matrix A.  N >= 0.

    @param[in,out]
    A       COMPLEX array, dimension (LDA,N)
            On entry, the Hermitian matrix A.  If UPLO = Magma_tally3Upper, the leading
            N-by-N upper triangular part of A contains the upper
            triangular part of the matrix A, and the strictly lower
            triangular part of A is not referenced.  If UPLO = Magma_tally3Lower, the
            leading N-by-N lower triangular part of A contains the lower
            triangular part of the matrix A, and the strictly upper
            triangular part of A is not referenced.
            On exit, if UPLO = Magma_tally3Upper, the diagonal and first superdiagonal
            of A are overwritten by the corresponding elements of the
            tridiagonal matrix T, and the elements above the first
            superdiagonal, with the array TAU, represent the orthogonal
            matrix Q as a product of elementary reflectors; if UPLO
            = Magma_tally3Lower, the diagonal and first subdiagonal of A are over-
            written by the corresponding elements of the tridiagonal
            matrix T, and the elements below the first subdiagonal, with
            the array TAU, represent the orthogonal matrix Q as a product
            of elementary reflectors. See Further Details.

    @param[in]
    lda     INTEGER
            The leading dimension of the array A.  LDA >= max(1,N).

    @param[out]
    d       COMPLEX array, dimension (N)
            The diagonal elements of the tridiagonal matrix T:
            D(i) = A(i,i).
 
    @param[out]
    e       COMPLEX array, dimension (N-1)
            The off-diagonal elements of the tridiagonal matrix T:
            E(i) = A(i,i+1) if UPLO = Magma_tally3Upper, E(i) = A(i+1,i) if UPLO = Magma_tally3Lower.

    @param[out]
    tau     COMPLEX array, dimension (N-1)
            The scalar factors of the elementary reflectors (see Further
            Details).

    @param[out]
    work    (workspace) COMPLEX array, dimension (MAX(1,LWORK))
            On exit, if INFO = 0, WORK[0] returns the optimal LWORK.

    @param[in]
    lwork   INTEGER
            The dimension of the array WORK.  LWORK >= N*NB, where NB is the
            optimal blocksize given by magma_tally3_get_chetrd_nb().
    \n
            If LWORK = -1, then a workspace query is assumed; the routine
            only calculates the optimal size of the WORK array, returns
            this value as the first entry of the WORK array, and no error
            message related to LWORK is issued by XERBLA.

    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument had an illegal value

    Further Details
    ---------------
    If UPLO = Magma_tally3Upper, the matrix Q is represented as a product of elementary
    reflectors

        Q = H(n-1) . . . H(2) H(1).

    Each H(i) has the form

        H(i) = I - tau * v * v'

    where tau is a complex scalar, and v is a complex vector with
    v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in
    A(1:i-1,i+1), and tau in TAU(i).

    If UPLO = Magma_tally3Lower, the matrix Q is represented as a product of elementary
    reflectors

        Q = H(1) H(2) . . . H(n-1).

    Each H(i) has the form

        H(i) = I - tau * v * v'

    where tau is a complex scalar, and v is a complex vector with
    v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in A(i+2:n,i),
    and tau in TAU(i).

    The contents of A on exit are illustrated by the following examples
    with n = 5:

    if UPLO = Magma_tally3Upper:                if UPLO = Magma_tally3Lower:

        (  d   e   v2  v3  v4 )              (  d                  )
        (      d   e   v3  v4 )              (  e   d              )
        (          d   e   v4 )              (  v1  e   d          )
        (              d   e  )              (  v1  v2  e   d      )
        (                  d  )              (  v1  v2  v3  e   d  )

    where d and e denote diagonal and off-diagonal elements of T, and vi
    denotes an element of the vector defining H(i).

    @ingroup magma_tally3_cheev_comp
    ********************************************************************/
extern "C" magma_tally3_int_t
magma_tally3_chetrd_mgpu(
    magma_tally3_int_t ngpu,
    magma_tally3_int_t nqueue, magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3FloatComplex *A, magma_tally3_int_t lda,
    float *d, float *e, magma_tally3FloatComplex *tau,
    magma_tally3FloatComplex *work, magma_tally3_int_t lwork,
    magma_tally3_int_t *info)
{
#define  A(i, j)     (A           + (j)*lda  + (i))
#define dA(id, i, j) (dA[(id)]    + (j)*ldda + (i))
#define dW(id, i, j) (dW[(id)] + (j)*ldda + (i))

    const char* uplo_ = lapack_uplo_const_tally3( uplo );
    
    magma_tally3_int_t nlocal, ldda;
    magma_tally3_int_t nb = magma_tally3_get_chetrd_nb(n), ib, ib2;

    const magma_tally3FloatComplex c_neg_one = MAGMA_tally3_C_NEG_ONE;
    const magma_tally3FloatComplex c_one     = MAGMA_tally3_C_ONE;
    const float             d_one     = MAGMA_tally3_D_ONE;
    
    #ifdef PROFILE_SY2RK
    float mv_time = 0.0;
    float up_time = 0.0;
    #endif

    magma_tally3_int_t kk, nx;
    magma_tally3_int_t i, ii, iii, j, dev, i_n;
    magma_tally3_int_t iinfo;
    magma_tally3_int_t ldwork, lddw, lwkopt, ldwork2, lhwork;
    magma_tally3_int_t lquery;
    
    // set pointers to NULL so it is safe to goto CLEANUP if any malloc fails.
    magma_tally3_queue_t queues[Magma_tally3MaxGPUs][10] = { { NULL, NULL } };
    magma_tally3_queue_t queues0[Magma_tally3MaxGPUs]    = { NULL };
    magma_tally3FloatComplex *hwork = NULL;
    magma_tally3FloatComplex_ptr dwork2[Magma_tally3MaxGPUs] = { NULL };
    magma_tally3FloatComplex_ptr dA[Magma_tally3MaxGPUs]     = { NULL };
    magma_tally3FloatComplex_ptr dW[Magma_tally3MaxGPUs]     = { NULL };

    *info = 0;
    int upper = (uplo == Magma_tally3Upper);
    lquery = (lwork == -1);
    if (! upper && uplo != Magma_tally3Lower) {
        *info = -1;
    } else if (n < 0) {
        *info = -2;
    } else if (lda < max(1,n)) {
        *info = -4;
    } else if (lwork < nb*n && ! lquery) {
        *info = -9;
    } else if ( nqueue > 2 ) {
        *info = 2;  // TODO fix
    }

    /* Determine the block size. */
    ldwork = n;
    lwkopt = n * nb;
    if (*info == 0) {
        work[0] = MAGMA_tally3_C_MAKE( lwkopt, 0 );
    }

    if (*info != 0) {
        magma_tally3_xerbla( __func__, -(*info) );
        return *info;
    }
    else if (lquery) {
        return *info;
    }

    /* Quick return if possible */
    if (n == 0) {
        work[0] = c_one;
        return *info;
    }

    magma_tally3_device_t orig_dev;
    magma_tally3_getdevice( &orig_dev );
    magma_tally3_queue_t orig_stream;
    magma_tally3blasGetKernelStream( &orig_stream );

    //#define PROFILE_SY2RK
    #ifdef PROFILE_SY2RK
    float times[11] = { 0 };
    magma_tally3_event_t start, stop;
    float etime;
    magma_tally3_setdevice( 0 );
    magma_tally3_event_create( &start );
    magma_tally3_event_create( &stop  );
    #endif

    ldda = roundup( lda, 32 );
    lddw = ldda;
    nlocal = nb*(1 + n/(nb*ngpu));
    ldwork2 = ldda*( ((n - 1)/nb + 1) + 1);  // i.e., ldda*(blocks + 1)
    for( dev=0; dev < ngpu; dev++ ) {
        magma_tally3_setdevice( dev );
        // TODO fix memory leak
        if ( MAGMA_tally3_SUCCESS != magma_tally3_cmalloc( &dA[dev],     nlocal*ldda + 3*lddw*nb ) ||
             MAGMA_tally3_SUCCESS != magma_tally3_cmalloc( &dwork2[dev], ldwork2 ) ) {
            *info = MAGMA_tally3_ERR_DEVICE_ALLOC;
            goto CLEANUP;
        }
        dW[dev] = dA[dev] + nlocal*ldda;
        
        for( kk=0; kk < nqueue; kk++ ) {
            magma_tally3_queue_create( &queues[dev][kk] );
        }
        queues0[dev] = queues[dev][0];
    }
    
    lhwork = nqueue*ngpu*n;
    if ( MAGMA_tally3_SUCCESS != magma_tally3_cmalloc_pinned( &hwork, lhwork ) ) {
        *info = MAGMA_tally3_ERR_HOST_ALLOC;
        goto CLEANUP;
    }

    // crossover point: use CPU code for last nx columns
    //if (n < 2048)
    //    nx = n;
    //else
    //    nx = 512;
    nx = min( 128, n );  // nx <= n is required

    if (upper) {
        /* Copy the matrix to the GPU */
        if (1 <= n-nx) {
            magma_tally3_chtodhe( ngpu, uplo, n, nb, A, lda, dA, ldda, queues, &iinfo );
        }

        /*  Reduce the upper triangle of A.
            Columns 1:kk are handled by the unblocked method. */
        for (i = nb*((n-1)/nb); i >= nx; i -= nb) {
            ib = min(nb, n-i);

            ii  = nb*(i/(nb*ngpu));
            dev = (i/nb)%ngpu;

            /* wait for the next panel */
            if (i != nb*((n-1)/nb)) {
                magma_tally3_setdevice( dev );
                magma_tally3_queue_sync( queues[dev][0] );
            }

            magma_tally3_clatrd_mgpu( ngpu, uplo, i+ib, ib, nb,
                               A(0, 0), lda, e, tau,
                               work, ldwork,
                               dA, ldda, 0,
                               dW, i+ib,
                               hwork,  lhwork,
                               dwork2, ldwork2,
                               queues0 );

            magma_tally3_cher2k_mgpu( ngpu, Magma_tally3Upper, Magma_tally3NoTrans, nb, i, ib,
                               c_neg_one, dW, i+ib, 0,
                               d_one,     dA, ldda, 0,
                               nqueue, queues);

            /* get the next panel */
            if (i-nb >= nx ) {
                ib2 = min(nb, n-(i-nb));
                
                ii  = nb*((i-nb)/(nb*ngpu));
                dev = ((i-nb)/nb)%ngpu;
                magma_tally3_setdevice( dev );
                
                magma_tally3_cgetmatrix_async( (i-nb)+ib2, ib2,
                                        dA(dev, 0, ii), ldda,
                                        A(0, i-nb),     lda,
                                        queues[dev][0] );
            }

            /* Copy superdiagonal elements back into A, and diagonal
               elements into D */
            for (j = i; j < i+ib; ++j) {
                if ( j > 0 ) {
                    *A(j-1,j) = MAGMA_tally3_C_MAKE( e[j - 1], 0 );
                }
                d[j] = MAGMA_tally3_C_REAL( *A(j, j) );
            }
        } /* end of for i=... */
      
        if ( nx > 0 ) {
            if (1 <= n-nx) { /* else A is already on CPU */
                for (i=0; i < nx; i += nb) {
                    ib = min(nb, n-i);
                    ii  = nb*(i/(nb*ngpu));
                    dev = (i/nb)%ngpu;
                
                    magma_tally3_setdevice( dev );
                    magma_tally3_cgetmatrix_async( nx, ib,
                                            dA(dev, 0, ii), ldda,
                                            A(0, i),        lda,
                                            queues[dev][0] );
                }
            }
            
            for( dev=0; dev < ngpu; dev++ ) {
                magma_tally3_setdevice( dev );
                magma_tally3_queue_sync( queues[dev][0] );
            }
            /* Use CPU code to reduce the last or only block */
            lapackf77_chetrd( uplo_, &nx, A(0, 0), &lda, d, e, tau,
                              work, &lwork, &iinfo );
        }
    }
    else {
        trace_init( 1, ngpu, nqueue, (CUstream_st**)queues );
        /* Copy the matrix to the GPU */
        if (1 <= n-nx) {
            magma_tally3_chtodhe( ngpu, uplo, n, nb, A, lda, dA, ldda, queues, &iinfo );
        }

        /* Reduce the lower triangle of A */
        for (i = 0; i < n-nx; i += nb) {
            ib = min(nb, n-i);

            ii  = nb*(i/(nb*ngpu));
            dev = (i/nb)%ngpu;
            /* Reduce columns i:i+ib-1 to tridiagonal form and form the
               matrix W which is needed to update the unreduced part of
               the matrix */

            /*   Get the current panel (no need for the 1st iteration) */
            if (i != 0) {
                magma_tally3_setdevice( dev );
                trace_gpu_start( dev, 0, "comm", "get" );
                magma_tally3_cgetmatrix_async( n-i, ib,
                                        dA(dev, i, ii), ldda,
                                        A(i,i),         lda,
                                        queues[dev][0] );
                trace_gpu_end( dev, 0 );
                magma_tally3_queue_sync( queues[dev][0] );
                magma_tally3_setdevice( 0 );
            }
            
            magma_tally3_clatrd_mgpu( ngpu, uplo, n-i, ib, nb,
                               A(i, i), lda, &e[i], &tau[i],
                               work, ldwork,
                               dA, ldda, i,
                               dW, n-i,
                               hwork,  lhwork,
                               dwork2, ldwork2,
                               queues0 );

            #ifdef PROFILE_SY2RK
            magma_tally3_setdevice( 0 );
            if ( i > 0 ) {
                cudaEventElapsedTime( &etime, start, stop );
                up_time += (etime/1000.0);
            }
            magma_tally3_event_record( start, 0 );
            #endif
            
            magma_tally3_cher2k_mgpu( ngpu, Magma_tally3Lower, Magma_tally3NoTrans, nb, n-i-ib, ib,
                               c_neg_one, dW, n-i, ib,
                               d_one, dA, ldda, i+ib, nqueue, queues);
            
            #ifdef PROFILE_SY2RK
            magma_tally3_setdevice( 0 );
            magma_tally3_event_record( stop, 0 );
            #endif

            /* Copy subdiagonal elements back into A, and diagonal
               elements into D */
            for (j = i; j < i+ib; ++j) {
                if ( j+1 < n ) {
                    *A(j+1,j) = MAGMA_tally3_C_MAKE( e[j], 0 );
                }
                d[j] = MAGMA_tally3_C_REAL( *A(j, j) );
            }
        } /* for i=... */

        /* Use CPU code to reduce the last or only block */
        if ( i < n ) {
            iii = i;
            i_n = n-i;
            if ( i > 0 ) {
                for (; i < n; i += nb) {
                    ib = min(nb, n-i);
                    ii  = nb*(i/(nb*ngpu));
                    dev = (i/nb)%ngpu;
                
                    magma_tally3_setdevice( dev );
                    magma_tally3_cgetmatrix_async( i_n, ib,
                                            dA(dev, iii, ii), ldda,
                                            A(iii, i),        lda,
                                            queues[dev][0] );
                }
                for( dev=0; dev < ngpu; dev++ ) {
                    magma_tally3_setdevice( dev );
                    magma_tally3_queue_sync( queues[dev][0] );
                }
            }
            lapackf77_chetrd( uplo_, &i_n, A(iii, iii), &lda, &d[iii], &e[iii],
                              &tau[iii], work, &lwork, &iinfo );
        }
    }
    
    for( dev=0; dev < ngpu; dev++ ) {
        magma_tally3_setdevice( dev );
        for( kk=0; kk < nqueue; kk++ ) {
            magma_tally3_queue_sync( queues[dev][kk] );
        }
    }
    
    #ifdef PROFILE_SY2RK
    magma_tally3_setdevice( 0 );
    if ( n > nx ) {
        cudaEventElapsedTime( &etime, start, stop );
        up_time += (etime/1000.0);
    }
    magma_tally3_event_destroy( start );
    magma_tally3_event_destroy( stop  );
    #endif

    trace_finalize( "chetrd.svg", "trace.css" );
    
    #ifdef PROFILE_SY2RK
    printf( " n=%d nb=%d\n", n, nb );
    printf( " Time in CLARFG: %.2e seconds\n", times[0] );
    //printf( " Time in CHEMV : %.2e seconds\n", mv_time );
    printf( " Time in CHER2K: %.2e seconds\n", up_time );
    #endif
    
CLEANUP:
    for( dev=0; dev < ngpu; dev++ ) {
        magma_tally3_setdevice( dev );
        for( kk=0; kk < nqueue; kk++ ) {
            magma_tally3_queue_destroy( queues[dev][kk] );
        }
        magma_tally3_free( dA[dev] );
        magma_tally3_free( dwork2[dev] );
    }
    magma_tally3_free_pinned( hwork );
    
    magma_tally3_setdevice( orig_dev );
    magma_tally3blasSetKernelStream( orig_stream );
    
    work[0] = MAGMA_tally3_C_MAKE( lwkopt, 0 );
    
    return *info;
} /* magma_tally3_chetrd */


// ----------------------------------------------------------------------
// TODO info is unused
extern "C" magma_tally3_int_t
magma_tally3_chtodhe(
    magma_tally3_int_t ngpu,
    magma_tally3_uplo_t uplo, magma_tally3_int_t n, magma_tally3_int_t nb,
    magma_tally3FloatComplex     *A,   magma_tally3_int_t lda,
    magma_tally3FloatComplex_ptr dA[], magma_tally3_int_t ldda,
    magma_tally3_queue_t queues[][10],
    magma_tally3_int_t *info)
{
    magma_tally3_device_t orig_dev;
    magma_tally3_getdevice( &orig_dev );
    
    magma_tally3_int_t k;
    if (uplo == Magma_tally3Lower) {
        /* go through each block-column */
        magma_tally3_int_t j, jj, jb, mj;
        for (j=0; j < n; j += nb) {
            jj =  j/(nb*ngpu);
            k  = (j/nb)%ngpu;
            
            jb = min(nb, (n-j));
            mj = n-j;
            
            magma_tally3_setdevice( k );
            magma_tally3_csetmatrix_async( mj, jb,
                                     A(j,j),         lda,
                                    dA(k, j, jj*nb), ldda,
                                    queues[k][0] );
        }
    }
    else {
        /* go through each block-column */
        magma_tally3_int_t j, jj, jb, mj;
        for (j=0; j < n; j += nb) {
            jj =  j/(nb*ngpu);
            k  = (j/nb)%ngpu;
            
            jb = min(nb, (n-j));
            mj = j+jb;
            
            magma_tally3_setdevice( k );
            magma_tally3_csetmatrix_async( mj, jb,
                                     A(0, j),        lda,
                                    dA(k, 0, jj*nb), ldda,
                                    queues[k][0] );
        }
    }
    for( k=0; k < ngpu; k++ ) {
        magma_tally3_setdevice( k );
        magma_tally3_queue_sync( queues[k][0] );
    }
    magma_tally3_setdevice( orig_dev );
    
    return *info;
}


// ----------------------------------------------------------------------
extern "C" void
magma_tally3_cher2k_mgpu(
    magma_tally3_int_t ngpu,
    magma_tally3_uplo_t uplo, magma_tally3_trans_t trans, magma_tally3_int_t nb, magma_tally3_int_t n, magma_tally3_int_t k,
    magma_tally3FloatComplex alpha,
    magma_tally3FloatComplex_ptr dB[], magma_tally3_int_t lddb, magma_tally3_int_t b_offset,
    float beta,
    magma_tally3FloatComplex_ptr dC[], magma_tally3_int_t lddc, magma_tally3_int_t c_offset,
    magma_tally3_int_t nqueue, magma_tally3_queue_t queues[][10])
{
#define dB(id, i, j)  (dB[(id)] + (j)*lddb + (i)+b_offset)
#define dB1(id, i, j) (dB[(id)] + (j)*lddb + (i)+b_offset) + k*lddb
#define dC(id, i, j)  (dC[(id)] + (j)*lddc + (i))

    magma_tally3_int_t i, id, ib, ii, kk, n1;
    magma_tally3FloatComplex c_one = MAGMA_tally3_C_ONE;

    magma_tally3_device_t orig_dev;
    magma_tally3_getdevice( &orig_dev );
    magma_tally3_queue_t orig_stream;
    magma_tally3blasGetKernelStream( &orig_stream );
    
    /* diagonal update */
    for( i=0; i < n; i += nb ) {
        id = ((i+c_offset)/nb)%ngpu;
        kk = (i/(nb*ngpu))%nqueue;
        magma_tally3_setdevice( id );
        magma_tally3blasSetKernelStream( queues[id][kk] );

        ib = min(nb, n-i);
        ii = nb*((i+c_offset)/(nb*ngpu));

        /* cher2k on diagonal block */
        trace_gpu_start( id, kk, "syr2k", "syr2k" );
        magma_tally3_cher2k( uplo, trans, ib, k,
                      alpha, dB1(id, i,          0 ), lddb,
                             dB(id,  i,          0 ), lddb,
                      beta,  dC(id,  i+c_offset, ii), lddc );
        trace_gpu_end( id, kk );
    }

    /* off-diagonal update */
    if (uplo == Magma_tally3Upper) {
        for( i=nb; i < n; i += nb ) {
            id = ((i+c_offset)/nb)%ngpu;
            kk = (i/(nb*ngpu))%nqueue;
            magma_tally3_setdevice( id );
            magma_tally3blasSetKernelStream( queues[id][kk] );
            
            ib = min(nb, n-i);
            ii = nb*((i+c_offset)/(nb*ngpu));
            magma_tally3_cgemm( Magma_tally3NoTrans, Magma_tally3ConjTrans, i, ib, k,
                         alpha, dB1(id, 0, 0 ), lddb,
                                dB(id,  i, 0 ), lddb,
                         c_one, dC(id,  0, ii), lddc );
        }
    }
    else {
        for( i=0; i < n-nb; i += nb ) {
            id = ((i+c_offset)/nb)%ngpu;
            kk = (i/(nb*ngpu))%nqueue;
            magma_tally3_setdevice( id );
            magma_tally3blasSetKernelStream( queues[id][kk] );
            
            ib = min(nb, n-i);
            ii = nb*((i+c_offset)/(nb*ngpu));
            n1 = n-i-ib;
            
            // cgemm on off-diagonal blocks
            trace_gpu_start( id, kk, "gemm_up", "gemm_up" );
            magma_tally3_cgemm( Magma_tally3NoTrans, Magma_tally3ConjTrans, n1, ib, k,
                         alpha, dB1(id, i+ib,          0 ), lddb,
                                dB(id,  i,             0 ), lddb,
                         c_one, dC(id,  i+c_offset+ib, ii), lddc );
            trace_gpu_end( id, kk );
        }
    }

    if (uplo == Magma_tally3Upper) {
        for( i=nb; i < n; i += nb ) {
            id = ((i+c_offset)/nb)%ngpu;
            kk = (i/(nb*ngpu))%nqueue;
            magma_tally3_setdevice( id );
            magma_tally3blasSetKernelStream( queues[id][kk] );
            
            ib = min(nb, n-i);
            ii = nb*((i+c_offset)/(nb*ngpu));
            magma_tally3_cgemm( Magma_tally3NoTrans, Magma_tally3ConjTrans, i, ib, k,
                         alpha, dB( id, 0, 0 ), lddb,
                                dB1(id, i, 0 ), lddb,
                         c_one, dC(id,  0, ii), lddc );
        }
    } else {
        for( i=0; i < n-nb; i += nb ) {
            id = ((i+c_offset)/nb)%ngpu;
            kk = (i/(nb*ngpu))%nqueue;
            magma_tally3_setdevice( id );
            magma_tally3blasSetKernelStream( queues[id][kk] );
            
            ib = min(nb, n-i);
            ii = nb*((i+c_offset)/(nb*ngpu));
            n1 = n-i-ib;
            
            /* cgemm on off-diagonal blocks */
            trace_gpu_start( id, kk, "gemm_up", "gemm_up" );
            magma_tally3_cgemm( Magma_tally3NoTrans, Magma_tally3ConjTrans, n1, ib, k,
                         alpha, dB(id,  i+ib,          0 ), lddb,
                                dB1(id, i,             0 ), lddb,
                         c_one, dC(id,  i+c_offset+ib, ii), lddc );
            trace_gpu_end( id, kk );
        }
    }

    for( id=0; id < ngpu; id++ ) {
        magma_tally3_setdevice( id );
        for( kk=0; kk < nqueue; kk++ ) {
            magma_tally3_queue_sync( queues[id][kk] );
        }
    }
    magma_tally3_setdevice( orig_dev );
    magma_tally3blasSetKernelStream( orig_stream );
}
