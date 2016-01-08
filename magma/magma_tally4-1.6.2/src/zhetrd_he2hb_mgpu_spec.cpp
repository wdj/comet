/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Azzam Haidar
       @author Stan Tomov

       @precisions normal z -> s d c

*/
#include "common_magma_tally4.h"
#include "magma_tally4_bulge.h"
#include "trace.h"

/**
    Purpose
    -------
    ZHETRD_HE2HB reduces a complex Hermitian matrix A to real symmetric
    band-diagonal form T by an orthogonal similarity transformation:
    Q**H * A * Q = T.
    This version stores the triangular matrices T used in the accumulated
    Householder transformations (I - V T V').

    Arguments
    ---------
    @param[in]
    uplo    magma_tally4_uplo_t
      -     = Magma_tally4Upper:  Upper triangle of A is stored;
      -     = Magma_tally4Lower:  Lower triangle of A is stored.

    @param[in]
    n       INTEGER
            The order of the matrix A.  N >= 0.

    @param[in,out]
    A       COMPLEX_16 array, dimension (LDA,N)
            On entry, the Hermitian matrix A.  If UPLO = Magma_tally4Upper, the leading
            N-by-N upper triangular part of A contains the upper
            triangular part of the matrix A, and the strictly lower
            triangular part of A is not referenced.  If UPLO = Magma_tally4Lower, the
            leading N-by-N lower triangular part of A contains the lower
            triangular part of the matrix A, and the strictly upper
            triangular part of A is not referenced.
            On exit, if UPLO = Magma_tally4Upper, the Upper band-diagonal of A is
            overwritten by the corresponding elements of the
            band-diagonal matrix T, and the elements above the band
            diagonal, with the array TAU, represent the orthogonal
            matrix Q as a product of elementary reflectors; if UPLO
            = Magma_tally4Lower, the the Lower band-diagonal of A is overwritten by
            the corresponding elements of the band-diagonal
            matrix T, and the elements below the band-diagonal, with
            the array TAU, represent the orthogonal matrix Q as a product
            of elementary reflectors. See Further Details.

    @param[in]
    lda     INTEGER
            The leading dimension of the array A.  LDA >= max(1,N).

    @param[out]
    tau     COMPLEX_16 array, dimension (N-1)
            The scalar factors of the elementary reflectors (see Further
            Details).

    @param[out]
    work    (workspace) COMPLEX_16 array, dimension (MAX(1,LWORK))
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
    dT      COMPLEX_16 array on the GPU, dimension N*NB,
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
    If UPLO = Magma_tally4Upper, the matrix Q is represented as a product of elementary
    reflectors

        Q = H(n-1) . . . H(2) H(1).

    Each H(i) has the form

        H(i) = I - tau * v * v'

    where tau is a complex scalar, and v is a complex vector with
    v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in
    A(1:i-1,i+1), and tau in TAU(i).

    If UPLO = Magma_tally4Lower, the matrix Q is represented as a product of elementary
    reflectors

        Q = H(1) H(2) . . . H(n-1).

    Each H(i) has the form

        H(i) = I - tau * v * v'

    where tau is a complex scalar, and v is a complex vector with
    v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in A(i+2:n,i),
    and tau in TAU(i).

    The contents of A on exit are illustrated by the following examples
    with n = 5:

    if UPLO = Magma_tally4Upper:                if UPLO = Magma_tally4Lower:

      (  d   e   v2  v3  v4 )              (  d                  )
      (      d   e   v3  v4 )              (  e   d              )
      (          d   e   v4 )              (  v1  e   d          )
      (              d   e  )              (  v1  v2  e   d      )
      (                  d  )              (  v1  v2  v3  e   d  )

    where d and e denote diagonal and off-diagonal elements of T, and vi
    denotes an element of the vector defining H(i).

    @ingroup magma_tally4_zheev_2stage
    ********************************************************************/
extern "C" magma_tally4_int_t
magma_tally4_zhetrd_he2hb_mgpu_spec(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n, magma_tally4_int_t nb,
    magma_tally4DoubleComplex *A, magma_tally4_int_t lda,
    magma_tally4DoubleComplex *tau,
    magma_tally4DoubleComplex *work, magma_tally4_int_t lwork,
    magma_tally4DoubleComplex_ptr dAmgpu[], magma_tally4_int_t ldda,
    magma_tally4DoubleComplex_ptr dTmgpu[], magma_tally4_int_t lddt,
    magma_tally4_int_t ngpu, magma_tally4_int_t distblk,
    magma_tally4_queue_t queues[][20], magma_tally4_int_t nqueue,
    magma_tally4_int_t *info)
{
    #define A(a_1,a_2)         ( A  + ((a_2)-1)*( lda) + (a_1)-1)
    #define tau_ref(a_1)       (tau + (a_1)-1)
    #define dT(a_0, a_1, a_2)  (dTmgpu[a_0]  + ((a_2)-1)*(lddt) + (a_1)-1)
    #define dA(a_0, a_1, a_2)  (dAmgpu[a_0]  + ((a_2)-1)*(ldda) + (a_1)-1)

    magma_tally4DoubleComplex c_neg_one  = MAGMA_tally4_Z_NEG_ONE;
    magma_tally4DoubleComplex c_neg_half = MAGMA_tally4_Z_NEG_HALF;
    magma_tally4DoubleComplex c_one  = MAGMA_tally4_Z_ONE;
    magma_tally4DoubleComplex c_zero = MAGMA_tally4_Z_ZERO;
    double  d_one = MAGMA_tally4_D_ONE;

    magma_tally4_int_t pm, pn, indi, indj, pk;
    magma_tally4_int_t pm_old=0, pn_old=0, indi_old=0, flipV=-1;
    magma_tally4_int_t iblock, idev, di;
    int i;
    int lwkopt;
    int lquery;


    assert (nqueue >= (ngpu+1));


    *info = 0;
    int upper = (uplo == Magma_tally4Upper);
    lquery = (lwork == -1);
    if (! upper && uplo != Magma_tally4Lower) {
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
        work[0] = MAGMA_tally4_Z_MAKE( lwkopt, 0 );
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

    magma_tally4_device_t orig_dev;
    magma_tally4_getdevice( &orig_dev );
    magma_tally4_queue_t orig_stream;
    magma_tally4blasGetKernelStream( &orig_stream );

    // limit to 16 threads
    magma_tally4_int_t orig_threads = magma_tally4_get_lapack_numthreads();
    magma_tally4_set_lapack_numthreads( min(orig_threads,16) );

    magma_tally4_int_t gnode[Magma_tally4MaxGPUs][Magma_tally4MaxGPUs+2];
    magma_tally4_int_t nbcmplx=0;
    magma_tally4_buildconnection_mgpu(gnode, &nbcmplx,  ngpu);
    #ifdef ENABLE_DEBUG
    printf(" Initializing communication pattern.... GPU-ncmplx %d\n\n", nbcmplx);
    #endif

    magma_tally4DoubleComplex *dspace[Magma_tally4MaxGPUs];
    magma_tally4DoubleComplex *dwork[Magma_tally4MaxGPUs], *dworkbis[Magma_tally4MaxGPUs];
    magma_tally4DoubleComplex *dvall[Magma_tally4MaxGPUs], *dv[Magma_tally4MaxGPUs], *dw[Magma_tally4MaxGPUs];
    magma_tally4DoubleComplex *workngpu[Magma_tally4MaxGPUs+1];
    magma_tally4_event_t     redevents[Magma_tally4MaxGPUs][Magma_tally4MaxGPUs*Magma_tally4MaxGPUs+10];
    magma_tally4_int_t nbevents = Magma_tally4MaxGPUs*Magma_tally4MaxGPUs;

    magma_tally4_int_t lddv        = ldda;
    magma_tally4_int_t lddw        = lddv;
    magma_tally4_int_t dwrk2siz    = ldda*nb*3;
    magma_tally4_int_t worksiz     = n*nb;
    magma_tally4_int_t devworksiz  = 2*nb*lddv + nb*lddw + nb*ldda + dwrk2siz; // 2*dv(dv0+dv1) + dw + dwork +dworkbis

    // local allocation and stream creation
    // TODO check malloc
    for( magma_tally4_int_t dev = 0; dev < ngpu; ++dev ) {
        magma_tally4_setdevice( dev );
        magma_tally4_zmalloc( &dspace[dev], devworksiz );
        magma_tally4_zmalloc_pinned ( &workngpu[dev], worksiz);
        dvall[dev]    = dspace[dev];
        dw[dev]       = dvall[dev]   + 2*nb*lddv;
        dwork[dev]    = dw[dev]      + nb*lddw;
        dworkbis[dev] = dwork[dev]   + nb*ldda;
        magma_tally4blasSetKernelStream( queues[ dev ][ 0 ] );
        for( magma_tally4_int_t i = 0; i < nbevents; ++i ) {
            cudaEventCreateWithFlags(&redevents[dev][i],cudaEventDisableTiming);
        }
    }
    magma_tally4_zmalloc_pinned ( &workngpu[ngpu], worksiz);
    magma_tally4DoubleComplex *worktest = NULL;
    //magma_tally4_zmalloc_cpu( &worktest, n*nb ); // not used
    // ======================
  

    magma_tally4DoubleComplex *hT = work + lwork - nb*nb;
    lwork -= nb*nb;
    memset( hT, 0, nb*nb*sizeof(magma_tally4DoubleComplex));

    if (upper) {
        printf("ZHETRD_HE2HB is not yet implemented for upper matrix storage. Exit.\n");
        exit(1);
    } else {
        /* Reduce the lower triangle of A */
        for (i = 1; i <= n-nb; i += nb) {
            indi = i+nb;
            indj = i;
            pm   = n - i - nb + 1;
            //pn   = min(i+nb-1, n-nb) -i + 1;
            pn   = nb;
            
            /*   Get the current panel (no need for the 1st iteration) */
            if (i > 1 ) {
                // zpanel_to_q copy the upper oof diagonal part of
                // the matrix to work to be restored later. acctually
                // the zero's and one's putted are not used this is only
                // because we don't have a function that copy only the
                // upper part of A to be restored after copying the
                // lookahead panel that has been computted from GPU to CPU.
                zpanel_to_q(Magma_tally4Upper, pn-1, A(i, i+1), lda, work);

                // find the device who own the panel then send it to the CPU.
                // below a -1 was added and then a -1 was done on di because of the fortran indexing
                iblock = ((i-1) / distblk) / ngpu;          // local block id
                di     = iblock*distblk + (i-1)%distblk;     // local index in parent matrix
                idev   = ((i-1) / distblk) % ngpu;          // device with this block


                //printf("Receiving panel ofsize %d %d from idev %d A(%d,%d) \n",(pm+pn), pn,idev,i-1,di);
                magma_tally4_setdevice( idev );

                //magma_tally4_device_sync();
                magma_tally4_zgetmatrix_async( (pm+pn), pn,
                                        dA(idev, i, di+1), ldda,
                                        A ( i, i), lda, queues[ idev ][ nqueue-1 ] );
              
                //magma_tally4_setdevice( 0 );
                //printf("updating zher2k on A(%d,%d) of size %d %d \n",indi_old+pn_old-1,indi_old+pn_old-1,pm_old-pn_old,pn_old);
                // compute ZHER2K_MGPU
                magma_tally4blas_zher2k_mgpu_spec(
                    Magma_tally4Lower, Magma_tally4NoTrans, pm_old-pn_old, pn_old,
                    c_neg_one, dv, pm_old, pn_old,
                               dw, pm_old, pn_old,
                    d_one,     dAmgpu, ldda, indi_old+pn_old-1,
                    ngpu, distblk, queues, 2 );
                //magma_tally4_setdevice( 0 );
             
                magma_tally4_setdevice( idev );
                magma_tally4_queue_sync( queues[idev][ nqueue-1 ] );
                //magma_tally4_setdevice( 0 );
                zq_to_panel(Magma_tally4Upper, pn-1, A(i, i+1), lda, work);
            }

            /* ==========================================================
               QR factorization on a panel starting nb off of the diagonal.
               Prepare the V and T matrices.
               ==========================================================  */
            lapackf77_zgeqrf(&pm, &pn, A(indi, indj), &lda,
                             tau_ref(i), work, &lwork, info);
            
            /* Form the matrix T */
            pk=min(pm,pn);
            lapackf77_zlarft( Magma_tally4ForwardStr, Magma_tally4ColumnwiseStr,
                              &pm, &pk, A(indi, indj), &lda,
                              tau_ref(i), hT, &nb);

            /* Prepare V - put 0s in the upper triangular part of the panel
               (and 1s on the diagonal), temporaly storing the original in work */
            zpanel_to_q(Magma_tally4Upper, pk, A(indi, indj), lda, work);



            /* Send V and T from the CPU to the GPU */
            // To be able to overlap the GET with the ZHER2K
            // it should be done on last stream.
            // TO Avoid a BUG that is overwriting the old_V
            // used atthis moment by zher2k with the new_V
            // send it now, we decide to have a flipflop
            // vector of Vs. if step%2=0 use V[0] else use V[nb*n]
            flipV = ((i-1)/nb)%2;
            for( magma_tally4_int_t dev = 0; dev < ngpu; ++dev ) {
                dv[dev] = dvall[dev] + flipV*nb*lddv;
            }

            for( magma_tally4_int_t dev = 0; dev < ngpu; ++dev ) {
                magma_tally4_setdevice( dev );
                // send V
                magma_tally4_zsetmatrix_async( pm, pk,
                                        A(indi, indj), lda,
                                        dv[dev], pm, queues[dev][nqueue-1] );

                // Send the triangular factor T to the GPU
                magma_tally4_zsetmatrix_async( pk, pk,
                                        hT, nb,
                                        dT(dev, 1, i), lddt, queues[dev][nqueue-1] );
            }

            /* ==========================================================
               Compute W:
               1. X = A (V T)
               2. W = X - 0.5* V * (T' * (V' * X))
               ==========================================================  */
            for( magma_tally4_int_t dev = 0; dev < ngpu; ++dev ) {
                // dwork = V T
                magma_tally4_setdevice( dev );
                magma_tally4blasSetKernelStream( queues[ dev ][ nqueue-1 ] );
                magma_tally4_queue_sync( queues[dev][nqueue-1] );
                magma_tally4_zgemm(Magma_tally4NoTrans, Magma_tally4NoTrans, pm, pk, pk,
                        c_one, dv[dev], pm,
                        dT(dev, 1, i), lddt,
                        c_zero, dwork[dev], pm);
            }

            // ===============================================
            //   SYNC TO BE SURE THAT BOTH V AND T WERE
            //   RECEIVED AND VT IS COMPUTED and SYR2K is done
            // ===============================================
            for( magma_tally4_int_t dev = 0; dev < ngpu; ++dev ) {
                magma_tally4_setdevice( dev );
                for( magma_tally4_int_t s = 0; s < nqueue; ++s )
                magma_tally4_queue_sync( queues[dev][s] );
            }


            // compute ZHEMM_MGPU
            // The broadcast of the result done inside this function
            // should be done in stream [0] because i am assuming this
            // for the GEMMs below otherwise I have to SYNC over the
            // Broadcasting stream.
            if (ngpu == 1) {
                magma_tally4blasSetKernelStream( queues[ 0 ][ 0 ] );
                magma_tally4_zhemm(
                    Magma_tally4Left, uplo, pm, pk,
                    c_one, dAmgpu[0]+(indi-1)*ldda+(indi-1), ldda,
                    dwork[0], pm,
                    c_zero, dw[0], pm);
            } else {
                magma_tally4blas_zhemm_mgpu_spec(
                    Magma_tally4Left, uplo, pm, pk,
                    c_one, dAmgpu, ldda, indi-1,
                                dwork, pm,
                    c_zero,     dw, pm, dworkbis, dwrk2siz, worktest, pm, workngpu, pm,
                    ngpu, distblk, queues, nqueue-1, redevents, nbevents, gnode, nbcmplx);
            }

            
            /* dwork = V*T already ==> dwork' = T'*V'
             * compute T'*V'*X ==> dwork'*W ==>
             * dwork + pm*nb = ((T' * V') * X) = dwork' * X = dwork' * W */
            for( magma_tally4_int_t dev = 0; dev < ngpu; ++dev ) {
                // Here we have to wait until the broadcast of ZHEMM has been done.
                // Note that the broadcast should be done on stream[0] so in a way
                // we can continue here on the same stream and avoid a sync
                magma_tally4_setdevice( dev );
                magma_tally4blasSetKernelStream( queues[ dev ][ dev ] );
                // magma_tally4_queue_sync( queues[dev][0] );
                magma_tally4_zgemm(Magma_tally4ConjTrans, Magma_tally4NoTrans, pk, pk, pm,
                            c_one, dwork[dev], pm,
                            dw[dev], pm,
                            c_zero, dworkbis[dev], nb);
                
                /* W = X - 0.5 * V * T'*V'*X
                 *   = X - 0.5 * V * (dwork + pm*nb) = W - 0.5 * V * (dwork + pm*nb) */
                magma_tally4_zgemm(Magma_tally4NoTrans, Magma_tally4NoTrans, pm, pk, pk,
                            c_neg_half, dv[dev], pm,
                            dworkbis[dev], nb,
                            c_one,     dw[dev], pm);
            }
            /* restore the panel it is put here to overlap with the previous GEMM*/
            zq_to_panel(Magma_tally4Upper, pk, A(indi, indj), lda, work);
            // ===============================================
            //   SYNC TO BE SURE THAT BOTH V AND W ARE DONE
            // ===============================================
            // Synchronise to be sure that W has been computed
            // because next ZHER2K use streaming and may happen
            // that lunch a gemm on stream 2 while stream 0
            // which compute those 2 GEMM above has not been
            // computed and also used for the same reason in
            // the panel update below and also for the last HER2K
            for( magma_tally4_int_t dev = 0; dev < ngpu; ++dev ) {
                magma_tally4_setdevice( dev );
                magma_tally4_queue_sync( queues[dev][dev] );
            }

            /* ==========================================================
               Update the unreduced submatrix A(i+ib:n,i+ib:n), using
               an update of the form:  A := A - V*W' - W*V'
               ==========================================================  */
            if (i + nb <= n-nb) {
                /* There would be next iteration;
                   do lookahead - update the next panel */
                // below a -1 was added and then a -1 was done on di because of the fortran indexing
                iblock = ((indi-1) / distblk) / ngpu;          // local block id
                di     = iblock*distblk + (indi-1)%distblk;     // local index in parent matrix
                idev   = ((indi-1) / distblk) % ngpu;          // device with this block
                magma_tally4_setdevice( idev );
                magma_tally4blasSetKernelStream( queues[ idev ][ nqueue-1 ] );
                //magma_tally4_queue_sync( queues[idev][0] ); removed because the sync has been done in the loop above
                magma_tally4_zgemm(Magma_tally4NoTrans, Magma_tally4ConjTrans, pm, pn, pn, c_neg_one,
                            dv[idev], pm,
                            dw[idev], pm, c_one,
                            dA(idev, indi, di+1), ldda);
            
                magma_tally4_zgemm(Magma_tally4NoTrans, Magma_tally4ConjTrans, pm, pn, pn, c_neg_one,
                            dw[idev], pm,
                            dv[idev], pm, c_one,
                            dA(idev, indi, di+1), ldda);
                //printf("updating next panel distblk %d  idev %d  on A(%d,%d) of size %d %d %d \n",distblk,idev,indi-1,di,pm,pn,pn);
            }
            else {
                /* no look-ahead as this is last iteration */
                // below a -1 was added and then a -1 was done on di because of the fortran indexing
                iblock = ((indi-1) / distblk) / ngpu;          // local block id
                di     = iblock*distblk + (indi-1)%distblk;     // local index in parent matrix
                idev   = ((indi-1) / distblk) % ngpu;          // device with this block
                magma_tally4_setdevice( idev );
                magma_tally4blasSetKernelStream( queues[ idev ][ 0 ] );
                //printf("LAST ZHER2K idev %d on A(%d,%d) of size %d \n",idev, indi-1,di,pk);
                magma_tally4_zher2k(Magma_tally4Lower, Magma_tally4NoTrans, pk, pk, c_neg_one,
                             dv[idev], pm,
                             dw[idev], pm, d_one,
                             dA(idev, indi, di+1), ldda);


                /* Send the last block to the CPU */
                zpanel_to_q(Magma_tally4Upper, pk-1, A(n-pk+1, n-pk+2), lda, work);
                magma_tally4_zgetmatrix( pk, pk,
                                  dA(idev, indi, di+1), ldda,
                                  A(n-pk+1, n-pk+1),  lda );
                zq_to_panel(Magma_tally4Upper, pk-1, A(n-pk+1, n-pk+2), lda, work);
            }

            indi_old = indi;
            //indj_old = indj;
            pm_old   = pm;
            pn_old   = pn;
        }  // end loop for (i)
    }// end of LOWER
    //magma_tally4_setdevice( 0 );

    for( magma_tally4_int_t dev = 0; dev < ngpu; ++dev ) {
        magma_tally4_setdevice( dev );
        magma_tally4_free( dspace[dev]);
        magma_tally4_free_pinned(workngpu[dev]);
        for( magma_tally4_int_t e = 0; e < nbevents; ++e ) {
            magma_tally4_event_destroy( redevents[dev][e] );
        }
    }
    magma_tally4_free_pinned(workngpu[ngpu]);
    magma_tally4_free_cpu(worktest);

    magma_tally4_setdevice( orig_dev );
    magma_tally4blasSetKernelStream( orig_stream );
    magma_tally4_set_lapack_numthreads( orig_threads );

    work[0] = MAGMA_tally4_Z_MAKE( lwkopt, 0 );
    return *info;
} /* magma_tally4_zhetrd_he2hb_mgpu_spec */
