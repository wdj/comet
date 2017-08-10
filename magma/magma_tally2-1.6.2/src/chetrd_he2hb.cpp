/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Azzam Haidar
       @author Stan Tomov

       @generated from zhetrd_he2hb.cpp normal z -> c, Fri Jan 30 19:00:18 2015

*/
#include "common_magma_tally2.h"
#include "trace.h"


/**
    Purpose
    -------
    CHETRD_HE2HB reduces a complex Hermitian matrix A to real symmetric
    band-diagonal form T by an orthogonal similarity transformation:
    Q**H * A * Q = T.
    This version stores the triangular matrices T used in the accumulated
    Householder transformations (I - V T V').

    Arguments
    ---------
    @param[in]
    uplo    magma_tally2_uplo_t
      -     = Magma_tally2Upper:  Upper triangle of A is stored;
      -     = Magma_tally2Lower:  Lower triangle of A is stored.

    @param[in]
    n       INTEGER
            The order of the matrix A.  N >= 0.

    @param[in,out]
    A       COMPLEX array, dimension (LDA,N)
            On entry, the Hermitian matrix A.  If UPLO = Magma_tally2Upper, the leading
            N-by-N upper triangular part of A contains the upper
            triangular part of the matrix A, and the strictly lower
            triangular part of A is not referenced.  If UPLO = Magma_tally2Lower, the
            leading N-by-N lower triangular part of A contains the lower
            triangular part of the matrix A, and the strictly upper
            triangular part of A is not referenced.
            On exit, if UPLO = Magma_tally2Upper, the Upper band-diagonal of A is
            overwritten by the corresponding elements of the
            band-diagonal matrix T, and the elements above the band
            diagonal, with the array TAU, represent the orthogonal
            matrix Q as a product of elementary reflectors; if UPLO
            = Magma_tally2Lower, the the Lower band-diagonal of A is overwritten by
            the corresponding elements of the band-diagonal
            matrix T, and the elements below the band-diagonal, with
            the array TAU, represent the orthogonal matrix Q as a product
            of elementary reflectors. See Further Details.

    @param[in]
    lda     INTEGER
            The leading dimension of the array A.  LDA >= max(1,N).

    @param[out]
    tau     COMPLEX array, dimension (N-1)
            The scalar factors of the elementary reflectors (see Further
            Details).

    @param[out]
    work    (workspace) COMPLEX array, dimension (MAX(1,LWORK))
            On exit, if INFO = 0, WORK[0] returns the optimal LWORK.

    @param[in]
    lwork   INTEGER
            The dimension of the array WORK.  LWORK >= 1.
            For optimum performance LWORK >= N*NB, where NB is the
            optimal blocksize.
    \n
            If LWORK = -1, then a workspace query is assumed; the routine
            only calculates the optimal size of the WORK array, returns
            this value as the first entry of the WORK array, and no error
            message related to LWORK is issued by XERBLA.

    @param[out]
    dT      COMPLEX array on the GPU, dimension N*NB,
            where NB is the optimal blocksize.
            On exit dT holds the upper triangular matrices T from the
            accumulated Householder transformations (I - V T V') used
            in the factorization. The nb x nb matrices T are ordered
            consecutively in memory one after another.

    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument had an illegal value

    Further Details
    ---------------
    If UPLO = Magma_tally2Upper, the matrix Q is represented as a product of elementary
    reflectors

       Q = H(n-1) . . . H(2) H(1).

    Each H(i) has the form

       H(i) = I - tau * v * v'

    where tau is a complex scalar, and v is a complex vector with
    v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in
    A(1:i-1,i+1), and tau in TAU(i).

    If UPLO = Magma_tally2Lower, the matrix Q is represented as a product of elementary
    reflectors

       Q = H(1) H(2) . . . H(n-1).

    Each H(i) has the form

       H(i) = I - tau * v * v'

    where tau is a complex scalar, and v is a complex vector with
    v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in A(i+2:n,i),
    and tau in TAU(i).

    The contents of A on exit are illustrated by the following examples
    with n = 5:

    if UPLO = Magma_tally2Upper:                if UPLO = Magma_tally2Lower:

      (  d   e   v2  v3  v4 )              (  d                  )
      (      d   e   v3  v4 )              (  e   d              )
      (          d   e   v4 )              (  v1  e   d          )
      (              d   e  )              (  v1  v2  e   d      )
      (                  d  )              (  v1  v2  v3  e   d  )

    where d and e denote diagonal and off-diagonal elements of T, and vi
    denotes an element of the vector defining H(i).

    @ingroup magma_tally2_cheev_2stage
    ********************************************************************/
extern "C" magma_tally2_int_t
magma_tally2_chetrd_he2hb(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n, magma_tally2_int_t nb,
    magma_tally2FloatComplex *A, magma_tally2_int_t lda,
    magma_tally2FloatComplex *tau,
    magma_tally2FloatComplex *work, magma_tally2_int_t lwork,
    magma_tally2FloatComplex_ptr dT,
    magma_tally2_int_t *info)
{
    #define  A(a_1,a_2)  ( A + ((a_2)-1)*( lda) + (a_1)-1)
    #define dA(a_1,a_2)  (dA + ((a_2)-1)*(ldda) + (a_1)-1)
    #define tau_ref(a_1) (tau + (a_1)-1)
    #define dT(a_1)      (dT + ((a_1)-1)*(lddt))

    int ldda = ((n+31)/32)*32;
    int lddt = nb;
   
    magma_tally2FloatComplex c_neg_one  = MAGMA_tally2_C_NEG_ONE;
    magma_tally2FloatComplex c_neg_half = MAGMA_tally2_C_NEG_HALF;
    magma_tally2FloatComplex c_one  = MAGMA_tally2_C_ONE;
    magma_tally2FloatComplex c_zero = MAGMA_tally2_C_ZERO;
    float  d_one = MAGMA_tally2_D_ONE;

    magma_tally2_int_t pm, pn, indi, indj, pk;
    magma_tally2_int_t pm_old=0, pn_old=0, indi_old=0, indj_old=0;

    int i;
    int lwkopt;
    int lquery;

    *info = 0;
    int upper = (uplo == Magma_tally2Upper);
    lquery = (lwork == -1);
    if (! upper && uplo != Magma_tally2Lower) {
        *info = -1;
    } else if (n < 0) {
        *info = -2;
    } else if (lda < max(1,n)) {
        *info = -4;
    } else if (lwork < 1 && ! lquery) {
        *info = -9;
    }

    /* Determine the block size. */
    lwkopt = n * nb;
    if (*info == 0) {
        work[0] = MAGMA_tally2_C_MAKE( lwkopt, 0 );
    }

    if (*info != 0)
        return *info;
    else if (lquery)
        return *info;

    /* Quick return if possible */
    if (n == 0) {
        work[0] = c_one;
        return *info;
    }

    magma_tally2_queue_t orig_stream;
    magma_tally2blasGetKernelStream( &orig_stream );
    
    magma_tally2FloatComplex *dA;
    if (MAGMA_tally2_SUCCESS != magma_tally2_cmalloc( &dA, (n + 2*nb)*ldda )) {
        *info = MAGMA_tally2_ERR_DEVICE_ALLOC;
        return *info;
    }

    // limit to 16 threads
    magma_tally2_int_t orig_threads = magma_tally2_get_lapack_numthreads();
    magma_tally2_set_lapack_numthreads( min(orig_threads,16) );

    /* Use the first panel of dA as work space */
    magma_tally2FloatComplex *dwork = dA + n*ldda;
    magma_tally2FloatComplex *dW    = dwork + nb*ldda;

    #ifdef TRACING
    char buf[80];
    #endif
    magma_tally2_queue_t stream[3];
    magma_tally2_queue_create( &stream[0] );
    magma_tally2_queue_create( &stream[1] );
    stream[2] = 0;  // default stream
    
    trace_init( 1, 1, 3, stream );

    magma_tally2FloatComplex *hT = work + lwork - nb*nb;
    lwork -= nb*nb;
    memset( hT, 0, nb*nb*sizeof(magma_tally2FloatComplex));

    magma_tally2blasSetKernelStream( stream[0] );
    magma_tally2_event_t Pupdate_event;
    cudaEventCreateWithFlags(&Pupdate_event,cudaEventDisableTiming);
    //magma_tally2_event_create(&Pupdate_event);


    if (upper) {
        printf("CHETRD_HE2HB is not yet implemented for upper matrix storage. Exit.\n");
        exit(1);
    } else {
        /* Copy the matrix to the GPU */
        if (1 <= n-nb) {
            trace_gpu_start( 0, 0, "set", "set A" );
            magma_tally2_csetmatrix_async( (n-nb), (n-nb),
                                    A(nb+1, nb+1),  lda,
                                    dA(nb+1, nb+1), ldda, stream[0] );
            trace_gpu_end( 0, 0 );
        }

        /* Reduce the lower triangle of A */
        for (i = 1; i <= n-nb; i += nb) {
             indi = i+nb;
             indj = i;
             pm   = n - i - nb + 1;
             //pn   = min(i+nb-1, n-nb) -i + 1;
             pn   = nb;
             
             /*   Get the current panel (no need for the 1st iteration) */
             if (i > 1 ) {
                 // cpanel_to_q_tally2 copy the upper oof diagonal part of
                 // the matrix to work to be restored later. acctually
                 //  the zero's and one's putted are not used this is only
                 //   because we don't have a function that copy only the
                 //    upper part of A to be restored after copying the
                 //    lookahead panel that has been computted from GPU to CPU.
                 cpanel_to_q_tally2(Magma_tally2Upper, pn-1, A(i, i+1), lda, work);

                 trace_gpu_start( 0, 1, "get", "get panel" );
                 //magma_tally2_queue_sync( stream[0] );
                 magma_tally2_queue_wait_event(stream[1], Pupdate_event);  //, 0);
                 magma_tally2_cgetmatrix_async( (pm+pn), pn,
                                         dA( i, i), ldda,
                                         A ( i, i), lda, stream[1] );
                 trace_gpu_end( 0, 1 );

                 trace_gpu_start( 0, 2, "her2k", "her2k" );
                 magma_tally2_cher2k(Magma_tally2Lower, Magma_tally2NoTrans, pm_old-pn_old, pn_old, c_neg_one,
                      dA(indi_old+pn_old, indj_old), ldda,
                      dW + pn_old,            pm_old, d_one,
                      dA(indi_old+pn_old, indi_old+pn_old), ldda);
                 trace_gpu_end( 0, 2 );

                 trace_cpu_start( 0, "sync", "sync on 1" );
                 magma_tally2_queue_sync( stream[1] );
                 trace_cpu_end( 0 );
                 cq_to_panel_tally2(Magma_tally2Upper, pn-1, A(i, i+1), lda, work);
             }

             /* ==========================================================
                QR factorization on a panel starting nb off of the diagonal.
                Prepare the V and T matrices.
                ==========================================================  */
             #ifdef TRACING
             snprintf( buf, sizeof(buf), "panel %d", i );
             #endif
             trace_cpu_start( 0, "geqrf", buf );
             lapackf77_cgeqrf(&pm, &pn, A(indi, indj), &lda,
                        tau_ref(i), work, &lwork, info);
             
             /* Form the matrix T */
                         pk=min(pm,pn);
             lapackf77_clarft( Magma_tally2ForwardStr, Magma_tally2ColumnwiseStr,
                           &pm, &pk, A(indi, indj), &lda,
                           tau_ref(i), hT, &nb);

             /* Prepare V - put 0s in the upper triangular part of the panel
                (and 1s on the diagonal), temporaly storing the original in work */
             cpanel_to_q_tally2(Magma_tally2Upper, pk, A(indi, indj), lda, work);
             trace_cpu_end( 0 );

             /* Send V from the CPU to the GPU */
             trace_gpu_start( 0, 0, "set", "set V and T" );
             magma_tally2_csetmatrix_async( pm, pk,
                                     A(indi, indj),  lda,
                                     dA(indi, indj), ldda, stream[0] );

             /* Send the triangular factor T to the GPU */
             magma_tally2_csetmatrix_async( pk, pk,
                                     hT,       nb,
                                     dT(i), lddt, stream[0] );
             trace_gpu_end( 0, 0 );
             
             /* ==========================================================
                Compute W:
                1. X = A (V T)
                2. W = X - 0.5* V * (T' * (V' * X))
                ==========================================================  */
             /* dwork = V T */
             trace_cpu_start( 0, "sync", "sync on 0" );
             // this sync is done here to be sure that the copy has been finished
             // because below we made a restore cq_to_panel_tally2 and this restore need
             // to ensure that the copy has been finished. we did it here to allow
             // overlapp of restore with next gemm and symm.
             magma_tally2_queue_sync( stream[0] );
             trace_cpu_end( 0 );
             
             trace_gpu_start( 0, 2, "gemm", "work = V*T" );
             magma_tally2_cgemm(Magma_tally2NoTrans, Magma_tally2NoTrans, pm, pk, pk,
                         c_one, dA(indi, indj), ldda,
                         dT(i), lddt,
                         c_zero, dwork, pm);
             trace_gpu_end( 0, 2 );
             
             /* dW = X = A*V*T. dW = A*dwork */
             trace_gpu_start( 0, 2, "hemm", "X = A*work" );
             magma_tally2_chemm(Magma_tally2Left, uplo, pm, pk,
                         c_one, dA(indi, indi), ldda,
                         dwork, pm,
                         c_zero, dW, pm);
             trace_gpu_end( 0, 2 );
             /* restore the panel */
             cq_to_panel_tally2(Magma_tally2Upper, pk, A(indi, indj), lda, work);
             
             /* dwork = V*T already ==> dwork' = T'*V'
              * compute T'*V'*X ==> dwork'*W ==>
              * dwork + pm*nb = ((T' * V') * X) = dwork' * X = dwork' * W */
             trace_gpu_start( 0, 2, "gemm", "work = T'*V'*X" );
             magma_tally2_cgemm(Magma_tally2ConjTrans, Magma_tally2NoTrans, pk, pk, pm,
                         c_one, dwork, pm,
                         dW, pm,
                         c_zero, dwork + pm*nb, nb);
             trace_gpu_end( 0, 2 );
             
             /* W = X - 0.5 * V * T'*V'*X
              *   = X - 0.5 * V * (dwork + pm*nb) = W - 0.5 * V * (dwork + pm*nb) */
             trace_gpu_start( 0, 2, "gemm", "W = X - 0.5*V*(T'*V'*X)" );
             magma_tally2_cgemm(Magma_tally2NoTrans, Magma_tally2NoTrans, pm, pk, pk,
                         c_neg_half, dA(indi, indj), ldda,
                         dwork + pm*nb, nb,
                         c_one,     dW, pm);
             trace_gpu_end( 0, 2 );

             /* ==========================================================
                Update the unreduced submatrix A(i+ib:n,i+ib:n), using
                an update of the form:  A := A - V*W' - W*V'
                ==========================================================  */
             if (i + nb <= n-nb) {
                 /* There would be next iteration;
                    do lookahead - update the next panel */
                 trace_gpu_start( 0, 2, "gemm", "gemm 4 next panel left" );
                 magma_tally2_cgemm(Magma_tally2NoTrans, Magma_tally2ConjTrans, pm, pn, pn, c_neg_one,
                             dA(indi, indj), ldda,
                             dW,                 pm, c_one,
                             dA(indi, indi), ldda);
                 trace_gpu_end( 0, 2 );
             
                 trace_gpu_start( 0, 2, "gemm", "gemm 5 next panel right" );
                 magma_tally2_cgemm(Magma_tally2NoTrans, Magma_tally2ConjTrans, pm, pn, pn, c_neg_one,
                             dW,                 pm,
                             dA(indi, indj), ldda, c_one,
                             dA(indi, indi), ldda);
                 trace_gpu_end( 0, 2 );
                 magma_tally2_event_record(Pupdate_event, stream[0]);
             }
             else {
                 /* no look-ahead as this is last iteration */
                 trace_gpu_start( 0, 2, "her2k", "her2k last iteration" );
                 magma_tally2_cher2k(Magma_tally2Lower, Magma_tally2NoTrans, pk, pk, c_neg_one,
                              dA(indi, indj), ldda,
                              dW,                 pm, d_one,
                              dA(indi, indi), ldda);
                 trace_gpu_end( 0, 2 );
             }
             
             indi_old = indi;
             indj_old = indj;
             pm_old   = pm;
             pn_old   = pn;
        }  // end loop for (i)

        /* Send the last block to the CPU */
        pk = min(pm,pn);
        if (1 <= n-nb) {
            cpanel_to_q_tally2(Magma_tally2Upper, pk-1, A(n-pk+1, n-pk+2), lda, work);
            trace_gpu_start( 0, 2, "get", "get last block" );
            magma_tally2_cgetmatrix( pk, pk,
                              dA(n-pk+1, n-pk+1), ldda,
                              A(n-pk+1, n-pk+1),  lda );
            trace_gpu_end( 0, 2 );
            cq_to_panel_tally2(Magma_tally2Upper, pk-1, A(n-pk+1, n-pk+2), lda, work);
        }
    }// end of LOWER
    
    trace_finalize( "chetrd_he2hb.svg", "trace.css" );

    magma_tally2_event_destroy( Pupdate_event );
    magma_tally2_queue_destroy( stream[0] );
    magma_tally2_queue_destroy( stream[1] );
    magma_tally2_free( dA );
    work[0] = MAGMA_tally2_C_MAKE( lwkopt, 0 );

    magma_tally2blasSetKernelStream( orig_stream );
    magma_tally2_set_lapack_numthreads( orig_threads );

    return *info;
} /* magma_tally2_chetrd_he2hb */