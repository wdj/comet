/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions normal z -> s d c

*/
#include "common_magma_tally2.h"
#include "trace.h"

#define PRECISION_z

/* === Define what BLAS to use ============================================ */
#if defined(PRECISION_s) || defined(PRECISION_d)
#define ZTRSM_WORK
//#define magma_tally2_ztrsm magma_tally2blas_ztrsm
#endif
/* === End defining what BLAS to use ======================================= */

/**
    Purpose
    -------
    ZPOTRF computes the Cholesky factorization of a complex Hermitian
    positive definite matrix dA.
    Auxiliary subroutine for zpotrf2_ooc. It is multiple gpu interface to compute
    Cholesky of a "rectangular" matrix.

    The factorization has the form
       dA = U**H * U,   if UPLO = Magma_tally2Upper, or
       dA = L  * L**H,  if UPLO = Magma_tally2Lower,
    where U is an upper triangular matrix and L is lower triangular.

    This is the block version of the algorithm, calling Level 3 BLAS.

    Arguments
    ---------
    @param[in]
    uplo    magma_tally2_uplo_t
      -     = Magma_tally2Upper:  Upper triangle of dA is stored;
      -     = Magma_tally2Lower:  Lower triangle of dA is stored.

    @param[in]
    n       INTEGER
            The order of the matrix dA.  N >= 0.

    @param[in,out]
    dA      COMPLEX_16 array on the GPU, dimension (LDDA,N)
            On entry, the Hermitian matrix dA.  If UPLO = Magma_tally2Upper, the leading
            N-by-N upper triangular part of dA contains the upper
            triangular part of the matrix dA, and the strictly lower
            triangular part of dA is not referenced.  If UPLO = Magma_tally2Lower, the
            leading N-by-N lower triangular part of dA contains the lower
            triangular part of the matrix dA, and the strictly upper
            triangular part of dA is not referenced.
    \n
            On exit, if INFO = 0, the factor U or L from the Cholesky
            factorization dA = U**H * U or dA = L * L**H.

    @param[in]
    ldda     INTEGER
            The leading dimension of the array dA.  LDDA >= max(1,N).
            To benefit from coalescent memory accesses LDDA must be
            divisible by 16.

    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument had an illegal value
      -     > 0:  if INFO = i, the leading minor of order i is not
                  positive definite, and the factorization could not be
                  completed.

    @ingroup magma_tally2_zposv_comp
    ********************************************************************/
extern "C" magma_tally2_int_t
magma_tally2_zpotrf3_mgpu(
    magma_tally2_int_t ngpu,
    magma_tally2_uplo_t uplo, magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2_int_t off_i, magma_tally2_int_t off_j, magma_tally2_int_t nb,
    magma_tally2DoubleComplex_ptr d_lA[],  magma_tally2_int_t ldda,
    magma_tally2DoubleComplex_ptr d_lP[],  magma_tally2_int_t lddp,
    magma_tally2DoubleComplex *A,          magma_tally2_int_t lda, magma_tally2_int_t h,
    magma_tally2_queue_t queues[][3], magma_tally2_event_t events[][5],
    magma_tally2_int_t *info )
{
#define Alo(i, j)  (A +             ((j)+off_j)*lda  + (nb*(((i)/nb)%h)+off_i))
#define Aup(i, j)  (A + (nb*(((j)/nb)%h)+off_j)*lda  +               (i+off_i))

#define dlA(id, i, j)     (d_lA[(id)] + (j)*ldda + (i))
#define dlP(id, i, j, k)  (d_lP[(id)] + (k)*nb*lddp + (j)*lddp + (i))
#define dlPT(id, i, j, k) (d_lP[(id)] + (k)*nb*lddp + (j)*nb   + (i))

    magma_tally2_int_t     j, jb, nb0, nb2, d, dd, id, j_local, j_local2, buf;
    magma_tally2DoubleComplex c_one     = MAGMA_tally2_Z_ONE;
    magma_tally2DoubleComplex c_neg_one = MAGMA_tally2_Z_NEG_ONE;
    double          d_one     =  1.0;
    double          d_neg_one = -1.0;
    int upper = (uplo == Magma_tally2Upper);
    magma_tally2DoubleComplex *dlpanel;
    magma_tally2_int_t n_local[Magma_tally2MaxGPUs], ldpanel;
    const magma_tally2_int_t stream1 = 0, stream2 = 1, stream3 = 2;
    
    *info = 0;
    if (! upper && uplo != Magma_tally2Lower) {
        *info = -1;
    } else if (n < 0) {
        *info = -2;
    } else if (!upper && ngpu*ldda < max(1,n)) {
        *info = -4;
    } else if (upper && ldda < max(1,m)) {
        *info = -4;
    }
    if (*info != 0) {
        magma_tally2_xerbla( __func__, -(*info) );
        return *info;
    }

    magma_tally2_device_t orig_dev;
    magma_tally2_getdevice( &orig_dev );
    magma_tally2_queue_t orig_stream;
    magma_tally2blasGetKernelStream( &orig_stream );
    
#if (defined(PRECISION_d) || defined(PRECISION_s)) && defined(ZTRSM_WORK)
    /* used by ztrsm_work */
    magma_tally2DoubleComplex c_zero    = MAGMA_tally2_Z_ZERO;
    int trsm_nb = 128;
    int trsm_n = trsm_nb*((nb+trsm_nb-1)/trsm_nb);
    magma_tally2DoubleComplex *d_dinvA[Magma_tally2MaxGPUs];
    magma_tally2DoubleComplex *d_x[Magma_tally2MaxGPUs];
    #define dinvA(d,j) &(d_dinvA[(d)][(j)*trsm_nb*trsm_n])
    #define dx(d,j) &(d_x[(d)][(j)*nb*m])
    /*
     * Allocate device memory for the inversed diagonal blocks, size=N*BLOCK_SIZE
     */
    // TODO free memory on failure.
    for( d=0; d < ngpu; d++ ) {
        magma_tally2_setdevice(d);
        if ( (MAGMA_tally2_SUCCESS != magma_tally2_zmalloc( &d_dinvA[d], 2*trsm_nb*trsm_n )) ||
             (MAGMA_tally2_SUCCESS != magma_tally2_zmalloc( &d_x[d],     2*nb*(upper ? n : m) )) ) {
            *info = MAGMA_tally2_ERR_DEVICE_ALLOC;
            return *info;
        }
    }
    magma_tally2_setdevice(0);
#endif
    
    /* initialization */
    for( d=0; d < ngpu; d++ ) {
        /* local-n and local-ld */
        if (upper) {
            n_local[d] = (n/(nb*ngpu))*nb;
            if (d < (n/nb)%ngpu)
                n_local[d] += nb;
            else if (d == (n/nb)%ngpu)
                n_local[d] += n%nb;
        } else {
            n_local[d] = (m/(nb*ngpu))*nb;
            if (d < (m/nb)%ngpu)
                n_local[d] += nb;
            else if (d == (m/nb)%ngpu)
                n_local[d] += m%nb;
        }
    }

    /* == initialize the trace */
    trace_init( 1, ngpu, 3, (CUstream_st**)queues );

    if (upper) {
        /* ---------------------------------------------- */
        /* Upper-triangular case                          */
        /* > Compute the Cholesky factorization A = U'*U. */
        /* ---------------------------------------------- */
        for (j=0; j < m; j += nb) {
            /* Set the GPU number that holds the current panel */
            id  = (j/nb)%ngpu;
            buf = (j/nb)%ngpu; // right now, we have ngpu buffers, so id and buf are the same..
            
            /* Set the local index where the current panel is */
            j_local = j/(nb*ngpu);
            jb = min(nb, (m-j));
 
            /* Update the current diagonal block on stream1 */
            magma_tally2_setdevice(id);
            if ( j > 0 ) {
                magma_tally2blasSetKernelStream( queues[id][stream1] );
                trace_gpu_start( id, stream1, "syrk", "syrk" );
                magma_tally2_zherk(Magma_tally2Upper, Magma_tally2ConjTrans, jb, j,
                            d_neg_one, dlA(id, 0, nb*j_local), ldda,
                            d_one,     dlA(id, j, nb*j_local), ldda);
                trace_gpu_end( id, stream1 );
            }
            
            /* send the diagonal to cpu on stream1 */
            trace_gpu_start( id, stream1, "comm", "D to CPU" );
            magma_tally2_zgetmatrix_async( jb, jb,
                                    dlA(id, j, nb*j_local), ldda,
                                    Aup(j,j),               lda,
                                    queues[id][stream1] );
            trace_gpu_end( id, stream1 );

            /* update off-diagonal blocks in the panel */
            if ( j > 0 ) {
                d = (j/nb+1)%ngpu;
                for( dd=0; dd < ngpu; dd++ ) {
                    j_local2 = j_local+1;
                    if ( d > id ) j_local2 --;
                    nb0 = nb*j_local2; // number of local columns in the panel, while jb is panel-size (number of rows)
            
                    if ( n_local[d] > nb0 ) {
                        magma_tally2_setdevice(d);
                        magma_tally2blasSetKernelStream( queues[d][stream2] );
                        if ( d == id ) {
                            dlpanel = dlA(d,0,nb*j_local);
                            ldpanel = ldda;
                            // the GPU owns the row from start, and no need of synch.
                            //magma_tally2_queue_wait_event( queues[d][stream2], events[d][0] ); // rows arrived at gpu
                            magma_tally2_queue_wait_event( queues[d][stream2], events[d][4] ); // wait for look-ahead trsm to finish
                        } else {
                            dlpanel = dlP(d,nb,0,buf);
                            ldpanel = lddp;
                            magma_tally2_queue_wait_event( queues[d][stream2], events[d][0] ); // rows arrived at gpu
                        }
                        trace_gpu_start( d, stream2, "gemm", "gemm" );
                        magma_tally2_zgemm(Magma_tally2ConjTrans, Magma_tally2NoTrans,
                                    jb, n_local[d]-nb0, j,
                                    c_neg_one, dlpanel,        ldpanel,
                                               dlA(d, 0, nb0), ldda,
                                    c_one,     dlA(d, j, nb0), ldda);
                        trace_gpu_end( d, stream2 );
                        magma_tally2_event_record( events[d][2], queues[d][stream2] );
                    }
                    d = (d+1)%ngpu;
                }
            }

            /* wait for panel and factorize it on cpu */
            magma_tally2_setdevice(id);
            magma_tally2_queue_sync( queues[id][stream1] );
            trace_cpu_start( 0, "getrf", "getrf" );
            lapackf77_zpotrf(Magma_tally2UpperStr, &jb, Aup(j,j), &lda, info);
            trace_cpu_end( 0 );
            if (*info != 0) {
                *info = *info + j;
                break;
            }
            
            /* send the diagonal to gpus on stream1 */
            if ( (j+jb) < n) {
                d = (j/nb+1)%ngpu;
                for( dd=0; dd < ngpu; dd++ ) {
                    if ( d == id ) {
                        dlpanel = dlA(d, j, nb*j_local);
                        ldpanel = ldda;
                    } else {
                        dlpanel = dlP(d,0,0,buf);
                        ldpanel = lddp;
                    }
                    magma_tally2_setdevice(d);
                    trace_gpu_start( d, stream1, "comm", "comm" );
                    magma_tally2_zsetmatrix_async( jb, jb,
                                            Aup(j,j), lda,
                                            dlpanel,  ldpanel,
                                            queues[d][stream1] );
                    trace_gpu_end( d, stream1 );
                    magma_tally2_event_record( events[d][1], queues[d][stream1] );
                    d = (d+1)%ngpu;
                }
            } else {
                magma_tally2_setdevice(id);
                trace_gpu_start( id, stream1, "comm", "comm" );
                magma_tally2_zsetmatrix_async( jb, jb,
                                        Aup(j,j),               lda,
                                        dlA(id, j, nb*j_local), ldda,
                                        queues[id][stream1] );
                trace_gpu_end( id, stream1 );
            }
            
            /* panel-factorize the off-diagonal */
            if ( (j+jb) < n) {
                d = (j/nb+1)%ngpu;
                for( dd=0; dd < ngpu; dd++ ) {
                    /* next column */
                    j_local2 = j_local+1;
                    if ( d > id ) j_local2--;
                    if ( d == id ) {
                        dlpanel = dlA(d,j,nb*j_local);
                        ldpanel = ldda;
                    } else {
                        dlpanel = dlP(d,0,0,buf);
                        ldpanel = lddp;
                    }
                    nb2 = n_local[d] - j_local2*nb;
                    
                    magma_tally2_setdevice(d);
                    if ( j+jb < m && d == (j/nb+1)%ngpu ) {
                        /* owns the next column, look-ahead next block on stream1 */
                        nb0 = min(nb, nb2);
                        magma_tally2blasSetKernelStream( queues[d][stream1] );
                        magma_tally2_queue_wait_event( queues[d][stream1], events[d][2] ); // wait for gemm update
                        trace_gpu_start( d, stream1, "trsm", "trsm" );
#if (defined(PRECISION_d) || defined(PRECISION_s)) && defined(ZTRSM_WORK)
                        magma_tally2blas_zlaset( Magma_tally2Full, trsm_nb, trsm_n, c_zero, c_zero, dinvA(d,0), trsm_nb );
                        magma_tally2blas_zlaset( Magma_tally2Full, nb0,     jb,     c_zero, c_zero, dx(d,0), nb0 );
                        magma_tally2blas_ztrsm_work( Magma_tally2Left, Magma_tally2Upper,
                                              Magma_tally2ConjTrans, Magma_tally2NonUnit,
                                              jb, nb0, c_one,
                                              dlpanel, ldpanel,
                                              dlA(d, j, nb*j_local2), ldda,
                                              1, dinvA(d,0), dx(d,0) );
#else
                        magma_tally2_ztrsm( Magma_tally2Left, Magma_tally2Upper,
                                     Magma_tally2ConjTrans, Magma_tally2NonUnit,
                                     jb, nb0, c_one,
                                     dlpanel,                ldpanel,
                                     dlA(d, j, nb*j_local2), ldda);
#endif
                        magma_tally2_event_record( events[d][4], queues[d][stream1] );
                        trace_gpu_end( d, stream1 );
                    } else if ( nb2 > 0 ) {
                        /* update all the blocks on stream2 */
                        magma_tally2_queue_wait_event( queues[d][stream2], events[d][1] ); // wait for cholesky factor
                        trace_gpu_start( d, stream2, "trsm", "trsm" );
                        magma_tally2blasSetKernelStream( queues[d][stream2] );
#if (defined(PRECISION_d) || defined(PRECISION_s)) && defined(ZTRSM_WORK)
                        magma_tally2blas_zlaset( Magma_tally2Full, trsm_nb, trsm_n, c_zero, c_zero, dinvA(d,0), trsm_nb );
                        magma_tally2blas_zlaset( Magma_tally2Full, nb2,     jb,     c_zero, c_zero, dx(d,0), nb2 );
                        magma_tally2blas_ztrsm_work( Magma_tally2Left, Magma_tally2Upper,
                                              Magma_tally2ConjTrans, Magma_tally2NonUnit,
                                              jb, nb2, c_one,
                                              dlpanel, ldpanel,
                                              dlA(d, j, nb*j_local2), ldda,
                                              1, dinvA(d,0), dx(d,0) );
#else
                        magma_tally2_ztrsm( Magma_tally2Left, Magma_tally2Upper,
                                     Magma_tally2ConjTrans, Magma_tally2NonUnit,
                                     jb, nb2, c_one,
                                     dlpanel,                ldpanel,
                                     dlA(d, j, nb*j_local2), ldda);
#endif
                        trace_gpu_end( d, stream2 );
                    }
                    d = (d+1)%ngpu;
                } /* end of for */

                /* ========================================================== */
                if ( j+jb < m ) {
                    d = (j/nb+1)%ngpu;
                    /* next column */
                    j_local2 = j_local+1;
                    if ( d > id ) j_local2--;
                    nb0 = min(nb, n_local[d]-nb*j_local2 );
                
                    /* even on 1 gpu, off-diagonals are copied to cpu (synchronize at the end).      *
                     * so we have the Cholesky factor, but only diagonal submatrix of the big panel, *
                     * on cpu at the end.                                                            */
                    int d2, buf2;
                    magma_tally2_setdevice(d);
                    /* lookahead done */
                    magma_tally2_queue_wait_event( queues[d][stream3], events[d][4] );
                
                    trace_gpu_start( d, stream3, "comm", "row to CPU" );
                    magma_tally2_zgetmatrix_async( (j+jb), nb0,
                                            dlA(d, 0, nb*j_local2), ldda,
                                            Aup(0,j+jb),            lda,
                                            queues[d][stream3] );
                    trace_gpu_end( d, stream3 );
                    magma_tally2_event_record( events[d][3], queues[d][stream3] );
                    /* needed on pluto */
                    //magma_tally2_queue_sync( queues[d][stream3] );
                
                    /* broadcast rows to gpus on stream2 */
                    buf2 = ((j+jb)/nb)%ngpu;
                    for( d2=0; d2 < ngpu; d2++ ) {
                        if ( d2 != d ) {
                            magma_tally2_setdevice(d2);
                            trace_gpu_start( d2, stream3, "comm", "row to GPUs" );
                            magma_tally2_queue_wait_event( queues[d2][stream3], events[d][3] ); // rows arrived at cpu on stream3
                            magma_tally2_zsetmatrix_async( j+jb, nb0,
                                                    Aup(0,j+jb),       lda,
                                                    dlP(d2,nb,0,buf2), lddp,
                                                    queues[d2][stream3] );
                            trace_gpu_end( d2, stream3 );
                            magma_tally2_event_record( events[d2][0], queues[d2][stream3] );
                        }
                    }

                    /* =========================== */
                    /* update the remaining blocks */
                    nb2 = n_local[d]-(nb*j_local2 + nb0);
                    if ( nb2 > 0 ) {
                        if ( d == id ) {
                            dlpanel = dlA(d, j, nb*j_local);
                            ldpanel = ldda;
                        } else {
                            dlpanel = dlP(d,0,0,buf);
                            ldpanel = lddp;
                        }
                        magma_tally2_setdevice(d);
                        magma_tally2blasSetKernelStream( queues[d][stream2] );
                        trace_gpu_start( d, stream2, "trsm", "trsm" );
#if (defined(PRECISION_d) || defined(PRECISION_s)) && defined(ZTRSM_WORK)
                        int flag = 0;
                        if (flag == 0) {
                            magma_tally2_queue_wait_event( queues[d][stream2], events[d][4] ); // lookahead -> diagonal inversion
                        } else {
                            magma_tally2blas_zlaset( Magma_tally2Full, trsm_nb, trsm_n, c_zero, c_zero, dinvA(d,flag), trsm_nb );
                            magma_tally2_queue_wait_event( queues[d][stream2], events[d][1] ); // panel received
                        }
                        magma_tally2blas_zlaset( Magma_tally2Full, nb2, jb, c_zero, c_zero, dx(d,1), nb2 );
                        magma_tally2blas_ztrsm_work( Magma_tally2Left, Magma_tally2Upper, Magma_tally2ConjTrans, Magma_tally2NonUnit,
                                              jb, nb2, c_one,
                                              dlpanel, ldpanel,
                                              dlA(d, j, nb*j_local2+nb0), ldda,
                                              flag, dinvA(d,flag), dx(d,1) );
#else
                        magma_tally2_queue_wait_event( queues[d][stream2], events[d][1] ); // wait for cholesky factor
                        magma_tally2_ztrsm( Magma_tally2Left, Magma_tally2Upper, Magma_tally2ConjTrans, Magma_tally2NonUnit,
                                     jb, nb2, c_one,
                                     dlpanel, ldpanel,
                                     dlA(d, j, nb*j_local2+nb0), ldda);
#endif
                        trace_gpu_end( d, stream2 );
                    }
                }
            } /* end of ztrsm */
        } /* end of for j=1, .., n */
    } else {
        /* ---------------------------------------------- */
        /* Lower-triangular case                          */
        /* > Compute the Cholesky factorization A = L*L'. */
        /* ---------------------------------------------- */
        for (j=0; j < n; j += nb) {
        
            /* Set the GPU number that holds the current panel */
            id  = (j/nb)%ngpu;
            buf = (j/nb)%ngpu;
            
            /* Set the local index where the current panel is */
            j_local = j/(nb*ngpu);
            jb = min(nb, (n-j));

            /* Update the current diagonal block on stream1 */
            magma_tally2_setdevice(id);
            if ( j > 0 ) {
                magma_tally2blasSetKernelStream( queues[id][stream1] );
                magma_tally2_zherk(Magma_tally2Lower, Magma_tally2NoTrans, jb, j,
                            d_neg_one, dlA(id, nb*j_local, 0), ldda,
                            d_one,     dlA(id, nb*j_local, j), ldda);
            }

            /* send the diagonal to cpu on stream1 */
            magma_tally2_zgetmatrix_async( jb, jb,
                                    dlA(id, nb*j_local, j), ldda,
                                    Alo(j,j),               lda,
                                    queues[id][stream1] );

            /* update off-diagonal blocks of the panel */
            if ( j > 0 ) {
                d = (j/nb+1)%ngpu;
                for( dd=0; dd < ngpu; dd++ ) {
                    j_local2 = j_local+1;
                    if ( d > id ) j_local2 --;
                    nb0 = nb*j_local2;
            
                    if ( nb0 < n_local[d] ) {
                        magma_tally2_setdevice(d);
                        magma_tally2blasSetKernelStream( queues[d][stream2] );
                        if ( d == id ) {
                            dlpanel = dlA(d, nb*j_local, 0);
                            ldpanel = ldda;
                            magma_tally2_queue_wait_event( queues[d][stream2], events[d][4] ); // wait for look-ahead trsm to finish
                        } else {
                            dlpanel = dlPT(d,0,nb,buf);
                            ldpanel = nb;
                            magma_tally2_queue_wait_event( queues[d][stream2], events[d][0] ); // rows arrived at gpu
                        }
                        magma_tally2_zgemm( Magma_tally2NoTrans, Magma_tally2ConjTrans,
                                     n_local[d]-nb0, jb, j,
                                     c_neg_one, dlA(d, nb0, 0), ldda,
                                                dlpanel,        ldpanel,
                                     c_one,     dlA(d, nb0, j), ldda);
                        magma_tally2_event_record( events[d][2], queues[d][stream2] );
                    }
                    d = (d+1)%ngpu;
                }
            }

            /* wait for the panel and factorized it on cpu */
            magma_tally2_setdevice(id);
            magma_tally2_queue_sync( queues[id][stream1] );
            lapackf77_zpotrf(Magma_tally2LowerStr, &jb, Alo(j,j), &lda, info);
            if (*info != 0) {
                *info = *info + j;
                break;
            }

            /* send the diagonal to gpus on stream1 */
            if ( (j+jb) < m) {
                d = (j/nb+1)%ngpu;
                for( dd=0; dd < ngpu; dd++ ) {
                    if ( d == id ) {
                        dlpanel = dlA(d, nb*j_local, j);
                        ldpanel = ldda;
                    } else {
                        dlpanel = dlPT(d, 0, 0, buf);
                        ldpanel = nb;
                    }
                    magma_tally2_setdevice(d);
                    magma_tally2_zsetmatrix_async( jb, jb,
                                            Alo(j,j), lda,
                                            dlpanel,  ldpanel,
                                            queues[d][stream1] );
                    magma_tally2_event_record( events[d][1], queues[d][stream1] );
                    d = (d+1)%ngpu;
                }
            } else {
                magma_tally2_setdevice(id);
                magma_tally2_zsetmatrix_async( jb, jb,
                                        Alo(j,j),               lda,
                                        dlA(id, nb*j_local, j), ldda,
                                        queues[id][stream1] );
            }

            /* panel factorize the off-diagonal */
            if ( (j+jb) < m) {
                d = (j/nb+1)%ngpu;
                for( dd=0; dd < ngpu; dd++ ) {
                    /* next column */
                    j_local2 = j_local+1;
                    if ( d > id ) j_local2--;
                    if ( d == id ) {
                        dlpanel = dlA(d, nb*j_local, j);
                        ldpanel = ldda;
                    } else {
                        dlpanel = dlPT(d, 0, 0, buf);
                        ldpanel = nb;
                    }
                    nb2 = n_local[d] - j_local2*nb;
                    nb0 = min(nb, nb2);
                    
                    magma_tally2_setdevice(d);
                    if ( j+nb < n && d == (j/nb+1)%ngpu ) { /* owns next column, look-ahead next block on stream1 */
                        if ( j > 0 ) magma_tally2_queue_wait_event( queues[d][stream1], events[d][2] ); // wait for gemm update
                        magma_tally2blasSetKernelStream( queues[d][stream1] );
#if (defined(PRECISION_d) || defined(PRECISION_s)) && defined(ZTRSM_WORK)
                        magma_tally2blas_zlaset( Magma_tally2Full, trsm_nb, trsm_n, c_zero, c_zero, dinvA(d,0), trsm_nb );
                        magma_tally2blas_zlaset( Magma_tally2Full, nb0,     jb,     c_zero, c_zero, dx(d,0), nb0 );
                        magma_tally2blas_ztrsm_work( Magma_tally2Right, Magma_tally2Lower,
                                              Magma_tally2ConjTrans, Magma_tally2NonUnit,
                                              nb0, jb, c_one,
                                              dlpanel, ldpanel,
                                              dlA(d, nb*j_local2, j), ldda,
                                              1, dinvA(d,0), dx(d,0) );
#else
                        magma_tally2_ztrsm( Magma_tally2Right, Magma_tally2Lower,
                                     Magma_tally2ConjTrans, Magma_tally2NonUnit,
                                     nb0, jb, c_one,
                                     dlpanel, ldpanel,
                                     dlA(d, nb*j_local2, j), ldda);
#endif
                        magma_tally2_event_record( events[d][4], queues[d][stream1] );
                    } else if ( nb2 > 0 ) { /* other gpus updating all the blocks on stream2 */
                        /* update the entire column */
                        magma_tally2_queue_wait_event( queues[d][stream2], events[d][1] ); // wait for the cholesky factor
                        magma_tally2blasSetKernelStream( queues[d][stream2] );
#if (defined(PRECISION_d) || defined(PRECISION_s)) && defined(ZTRSM_WORK)
                        magma_tally2blas_zlaset( Magma_tally2Full, trsm_nb, trsm_n, c_zero, c_zero, dinvA(d,0), trsm_nb );
                        magma_tally2blas_zlaset( Magma_tally2Full, nb2,     jb,     c_zero, c_zero, dx(d,0), nb2 );
                        magma_tally2blas_ztrsm_work( Magma_tally2Right, Magma_tally2Lower, Magma_tally2ConjTrans, Magma_tally2NonUnit,
                                              nb2, jb, c_one,
                                              dlpanel,                ldpanel,
                                              dlA(d, nb*j_local2, j), ldda,
                                              1, dinvA(d,0), dx(d,0) );
#else
                        magma_tally2_ztrsm( Magma_tally2Right, Magma_tally2Lower, Magma_tally2ConjTrans, Magma_tally2NonUnit,
                                     nb2, jb, c_one,
                                     dlpanel,                ldpanel,
                                     dlA(d, nb*j_local2, j), ldda);
#endif
                    }
                    d = (d+1)%ngpu;
                } /* end for d */

                /* ========================================================== */
                if ( j+jb < n ) {
                    d = (j/nb+1)%ngpu;
                    /* next column */
                    j_local2 = j_local+1;
                    if ( d > id ) j_local2--;
                    nb0 = min(nb, n_local[d]-nb*j_local2 );
                
                    /* even on 1 gpu, we copy off-diagonal to cpu (but don't synchronize).  */
                    /* so we have the Cholesky factor on cpu at the end.                    */
                    int d2, buf2;
//#define ZPOTRF_DEVICE_TO_DEVICE
#ifdef ZPOTRF_DEVICE_TO_DEVICE
                    // lookahead done
                
                    /* broadcast the rows to gpus */
                    buf2 = ((j+jb)/nb)%ngpu;
                    for( d2=0; d2 < ngpu; d2++ ) {
                        magma_tally2_setdevice(d2);
                        magma_tally2_queue_wait_event( queues[d2][stream3], events[d][4] );
                        if ( d2 != d ) {
                            magma_tally2_zcopymatrix_async( nb0, j+jb,
                                                     dlPT(d2,0,nb,buf2), nb, // first nbxnb reserved for diagonal block
                                                     dlA(d, nb*j_local2, 0), ldda,
                                                     queues[d2][stream3] );
                            magma_tally2_event_record( events[d2][0], queues[d2][stream3] );
                        } else {
                            magma_tally2_zgetmatrix_async( nb0, j+jb,
                                                    dlA(d, nb*j_local2, 0), ldda,
                                                    Alo(j+jb,0),            lda,
                                                    queues[d][stream3] );
                        }
                    }
#else
                    // lookahead done
                    magma_tally2_setdevice(d);
                    magma_tally2_queue_wait_event( queues[d][stream3], events[d][4] );
                    magma_tally2_zgetmatrix_async( nb0, j+jb,
                                            dlA(d, nb*j_local2, 0), ldda,
                                            Alo(j+jb,0),            lda,
                                            queues[d][stream3] );
                    magma_tally2_event_record( events[d][3], queues[d][stream3] );
                    /* syn on rows on CPU, seem to be needed on Pluto */
                    //magma_tally2_queue_sync( queues[d][stream3] );
                
                    /* broadcast the rows to gpus */
                    buf2 = ((j+jb)/nb)%ngpu;
                    for( d2=0; d2 < ngpu; d2++ ) {
                        if ( d2 != d ) {
                            magma_tally2_setdevice(d2);
                            magma_tally2_queue_wait_event( queues[d2][stream3], events[d][3] ); // getmatrix done
                            magma_tally2_zsetmatrix_async( nb0, j+jb,
                                                    Alo(j+jb,0),        lda,
                                                    dlPT(d2,0,nb,buf2), nb, // first nbxnb reserved for diagonal block
                                                    queues[d2][stream3] );
                            magma_tally2_event_record( events[d2][0], queues[d2][stream3] );
                        }
                    }
#endif
                    /* =================================== */
                    /* updates remaining blocks on stream2 */
                    nb2 = n_local[d] - (j_local2*nb + nb0);
                    if ( nb2 > 0 ) {
                        if ( d == id ) {
                            dlpanel = dlA(d, nb*j_local, j);
                            ldpanel = ldda;
                        } else {
                            dlpanel = dlPT(d,0,0,buf);
                            ldpanel = nb;
                        }
                        magma_tally2_setdevice(d);
                        magma_tally2blasSetKernelStream( queues[d][stream2] );
                        /* update the remaining blocks in the column */
#if (defined(PRECISION_d) || defined(PRECISION_s)) && defined(ZTRSM_WORK)
                        int flag = 0;
                        if (flag == 0) {
                            magma_tally2_queue_wait_event( queues[d][stream2], events[d][4] ); // lookahead -> diagonal inversion
                        } else {
                            magma_tally2blas_zlaset( Magma_tally2Full, trsm_nb, trsm_n, c_zero, c_zero, dinvA(d,flag), trsm_nb );
                            magma_tally2_queue_wait_event( queues[d][stream2], events[d][1] ); // panel received
                        }
                        magma_tally2blas_zlaset( Magma_tally2Full, nb2, jb, c_zero, c_zero, dx(d,1), nb2 );
                        magma_tally2blas_ztrsm_work( Magma_tally2Right, Magma_tally2Lower, Magma_tally2ConjTrans, Magma_tally2NonUnit,
                                              nb2, jb, c_one,
                                              dlpanel,                    ldpanel,
                                              dlA(d, nb*j_local2+nb0, j), ldda,
                                              flag, dinvA(d,flag), dx(d,1) );
#else
                        magma_tally2_queue_wait_event( queues[d][stream2], events[d][1] ); // panel received
                        magma_tally2_ztrsm( Magma_tally2Right, Magma_tally2Lower, Magma_tally2ConjTrans, Magma_tally2NonUnit,
                                     nb2, jb, c_one,
                                     dlpanel,                    ldpanel,
                                     dlA(d, nb*j_local2+nb0, j), ldda);
#endif
                    }
                }
            }
        }
    } /* end of else not upper */

    /* == finalize the trace == */
    trace_finalize( "zpotrf.svg", "trace.css" );
    for( d=0; d < ngpu; d++ ) {
        magma_tally2_setdevice(d);
        for( j=0; j < 3; j++ ) {
            magma_tally2_queue_sync( queues[d][j] );
        }
#if (defined(PRECISION_d) || defined(PRECISION_s)) && defined(ZTRSM_WORK)
        magma_tally2_free( d_dinvA[d] );
        magma_tally2_free( d_x[d] );
#endif
    }
    magma_tally2_setdevice( orig_dev );
    magma_tally2blasSetKernelStream( orig_stream );

    return *info;
} /* magma_tally2_zpotrf_mgpu */

#undef Alo
#undef Aup
#undef dlA
#undef dlP
#undef dlPT


#define A(i, j)  (A +(j)*lda  + (i))
#define dA(d, i, j) (dA[(d)]+(j)*ldda + (i))


// ----------------------------------------------------------------------
extern "C" magma_tally2_int_t
magma_tally2_zhtodpo(
    magma_tally2_int_t ngpu,
    magma_tally2_uplo_t uplo, magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2_int_t off_i, magma_tally2_int_t off_j, magma_tally2_int_t nb,
    magma_tally2DoubleComplex    *A,    magma_tally2_int_t lda,
    magma_tally2DoubleComplex_ptr dA[], magma_tally2_int_t ldda,
    magma_tally2_queue_t queues[][3],
    magma_tally2_int_t *info)
{
    magma_tally2_device_t orig_dev;
    magma_tally2_getdevice( &orig_dev );
    
    magma_tally2_int_t k;
    if (uplo == Magma_tally2Upper) {
        magma_tally2_int_t j, jj, jb, mj;
        
        /* go through each column */
        for (j=off_j; j < n; j += nb) {
            jj = (j-off_j)/(nb*ngpu);
            k  = ((j-off_j)/nb)%ngpu;
            
            jb = min(nb, (n-j));
            if (j+jb < off_j+m)
                mj = (j-off_i)+jb;
            else
                mj = m;

            magma_tally2_setdevice(k);
            magma_tally2_zsetmatrix_async( mj, jb,
                                    A(off_i, j),     lda,
                                    dA(k, 0, jj*nb), ldda,
                                    queues[k][0] );
        }
    }
    else {
        magma_tally2_int_t i, ii, ib, ni;
        
        /* go through each row */
        for (i=off_i; i < m; i += nb) {
            ii = (i-off_i)/(nb*ngpu);
            k  = ((i-off_i)/nb)%ngpu;
            
            ib = min(nb, (m-i));
            if (i+ib < off_i+n)
                ni = (i-off_i)+ib;
            else
                ni = n;
            
            magma_tally2_setdevice(k);
            magma_tally2_zsetmatrix_async( ib, ni,
                                    A(i, off_j),     lda,
                                    dA(k, ii*nb, 0), ldda,
                                    queues[k][0] );
        }
    }
    for( k=0; k < ngpu; k++ ) {
        magma_tally2_setdevice(k);
        magma_tally2_queue_sync( queues[k][0] );
    }
    magma_tally2_setdevice( orig_dev );

    return *info;
}


// ----------------------------------------------------------------------
extern "C" magma_tally2_int_t
magma_tally2_zdtohpo(
    magma_tally2_int_t ngpu,
    magma_tally2_uplo_t uplo, magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2_int_t off_i, magma_tally2_int_t off_j, magma_tally2_int_t nb, magma_tally2_int_t NB,
    magma_tally2DoubleComplex    *A,    magma_tally2_int_t lda,
    magma_tally2DoubleComplex_ptr dA[], magma_tally2_int_t ldda,
    magma_tally2_queue_t queues[][3],
    magma_tally2_int_t *info)
{
    magma_tally2_device_t orig_dev;
    magma_tally2_getdevice( &orig_dev );
    
    magma_tally2_int_t k;
    if (uplo == Magma_tally2Upper) {
        magma_tally2_int_t j, jj, jb, mj;
        
        /* go through each column */
        for (j=off_j+NB; j < n; j += nb) {
            jj =  (j-off_j)/(nb*ngpu);
            k  = ((j-off_j)/nb)%ngpu;
            
            jb = min(nb, (n-j));
            if (j+jb < off_j+m)
                mj = (j-off_i)+jb;
            else
                mj = m;

            magma_tally2_setdevice(k);
            magma_tally2_zgetmatrix_async( mj, jb,
                                    dA(k, 0, jj*nb), ldda,
                                    A(off_i, j),     lda,
                                    queues[k][0] );
            magma_tally2_queue_sync( queues[k][0] );
        }
    } else {
        magma_tally2_int_t i, ii, ib, ni;
        
        /* go through each row */
        for (i=off_i+NB; i < m; i += nb) {
            ii = (i-off_i)/(nb*ngpu);
            k  = ((i-off_i)/nb)%ngpu;
            
            ib = min(nb, (m-i));
            if (i+ib < off_i+n)
                ni = (i-off_i)+ib;
            else
                ni = n;
            
            magma_tally2_setdevice(k);
            magma_tally2_zgetmatrix_async( ib, ni,
                                    dA(k, ii*nb, 0), ldda,
                                    A(i, off_j),     lda,
                                    queues[k][0] );
            magma_tally2_queue_sync( queues[k][0] );
        }
    }
    /*for( k=0; k < ngpu; k++ ) {
        magma_tally2_setdevice(k);
        magma_tally2_queue_sync( queues[k][0] );
    }*/
    magma_tally2_setdevice( orig_dev );

    return *info;
}

#undef A
#undef dA