/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Raffaele Solca
       @author Azzam Haidar
       @author Stan Tomov

       @precisions normal z -> s d c

*/
#include "common_magma_tally4.h"

/**
    Purpose
    -------
    ZUNMQR overwrites the general complex M-by-N matrix C with

    @verbatim
                                SIDE = Magma_tally4Left    SIDE = Magma_tally4Right
    TRANS = Magma_tally4NoTrans:       Q * C               C * Q
    TRANS = Magma_tally4_ConjTrans:    Q**H * C            C * Q**H
    @endverbatim

    where Q is a complex unitary matrix defined as the product of k
    elementary reflectors

          Q = H(1) H(2) . . . H(k)

    as returned by ZGEQRF. Q is of order M if SIDE = Magma_tally4Left and of order N
    if SIDE = Magma_tally4Right.

    Arguments
    ---------
    @param[in]
    ngpu    INTEGER
            Number of GPUs to use. ngpu > 0.

    @param[in]
    side    magma_tally4_side_t
      -     = Magma_tally4Left:      apply Q or Q**H from the Left;
      -     = Magma_tally4Right:     apply Q or Q**H from the Right.

    @param[in]
    trans   magma_tally4_trans_t
      -     = Magma_tally4NoTrans:    No transpose, apply Q;
      -     = Magma_tally4_ConjTrans: Conjugate transpose, apply Q**H.

    @param[in]
    m       INTEGER
            The number of rows of the matrix C. M >= 0.

    @param[in]
    n       INTEGER
            The number of columns of the matrix C. N >= 0.

    @param[in]
    k       INTEGER
            The number of elementary reflectors whose product defines
            the matrix Q.
            If SIDE = Magma_tally4Left,  M >= K >= 0;
            if SIDE = Magma_tally4Right, N >= K >= 0.

    @param[in]
    A       COMPLEX_16 array, dimension (LDA,K)
            The i-th column must contain the vector which defines the
            elementary reflector H(i), for i = 1,2,...,k, as returned by
            ZGEQRF in the first k columns of its array argument A.

    @param[in]
    lda     INTEGER
            The leading dimension of the array A.
            If SIDE = Magma_tally4Left,  LDA >= max(1,M);
            if SIDE = Magma_tally4Right, LDA >= max(1,N).

    @param[in]
    tau     COMPLEX_16 array, dimension (K)
            TAU(i) must contain the scalar factor of the elementary
            reflector H(i), as returned by ZGEQRF.

    @param[in,out]
    C       COMPLEX_16 array, dimension (LDC,N)
            On entry, the M-by-N matrix C.
            On exit, C is overwritten by Q*C or Q**H*C or C*Q**H or C*Q.

    @param[in]
    ldc     INTEGER
            The leading dimension of the array C. LDC >= max(1,M).

    @param[out]
    work    (workspace) COMPLEX_16 array, dimension (MAX(1,LWORK))
            On exit, if INFO = 0, WORK[0] returns the optimal LWORK.

    @param[in]
    lwork   INTEGER
            The dimension of the array WORK.
            If SIDE = Magma_tally4Left,  LWORK >= max(1,N);
            if SIDE = Magma_tally4Right, LWORK >= max(1,M).
            For optimum performance LWORK >= N*NB if SIDE = Magma_tally4Left, and
            LWORK >= M*NB if SIDE = Magma_tally4Right, where NB is the optimal
            blocksize.
    \n
            If LWORK = -1, then a workspace query is assumed; the routine
            only calculates the optimal size of the WORK array, returns
            this value as the first entry of the WORK array, and no error
            message related to LWORK is issued by XERBLA.

    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument had an illegal value

    @ingroup magma_tally4_zgeqrf_comp
    ********************************************************************/
extern "C" magma_tally4_int_t
magma_tally4_zunmqr_m(
    magma_tally4_int_t ngpu,
    magma_tally4_side_t side, magma_tally4_trans_t trans,
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4DoubleComplex *A,    magma_tally4_int_t lda,
    magma_tally4DoubleComplex *tau,
    magma_tally4DoubleComplex *C,    magma_tally4_int_t ldc,
    magma_tally4DoubleComplex *work, magma_tally4_int_t lwork,
    magma_tally4_int_t *info)
{
#define  A(i, j) (A + (j)*lda  + (i))
#define  C(i, j) (C + (j)*ldc  + (i))

#define    dC(gpui,      i, j) (dw[gpui] + (j)*lddc + (i))
#define  dA_c(gpui, ind, i, j) (dw[gpui] + maxnlocal*lddc + (ind)*lddar*lddac + (i) + (j)*lddac)
#define  dA_r(gpui, ind, i, j) (dw[gpui] + maxnlocal*lddc + (ind)*lddar*lddac + (i) + (j)*lddar)
#define    dT(gpui, ind)       (dw[gpui] + maxnlocal*lddc + 2*lddac*lddar + (ind)*((nb+1)*nb))
#define dwork(gpui, ind)       (dw[gpui] + maxnlocal*lddc + 2*lddac*lddar + 2*((nb+1)*nb) + (ind)*(lddwork*nb))

    magma_tally4DoubleComplex c_zero = MAGMA_tally4_Z_ZERO;
    magma_tally4DoubleComplex c_one  = MAGMA_tally4_Z_ONE;

    const char* side_  = lapack_side_const( side );
    const char* trans_ = lapack_trans_const( trans );

    // TODO fix memory leak (alloc after argument checks)
    magma_tally4_int_t nb = 128;
    magma_tally4DoubleComplex *T;
    magma_tally4_zmalloc_pinned(&T, nb*nb);
    //printf("calling zunmqr_m with nb=%d\n", (int) nb);

    magma_tally4DoubleComplex* dw[Magma_tally4MaxGPUs];
    magma_tally4_queue_t stream [Magma_tally4MaxGPUs][2];
    magma_tally4_event_t  event [Magma_tally4MaxGPUs][2];

    magma_tally4_int_t ind_c;
    magma_tally4_device_t igpu;
    
    magma_tally4_device_t orig_dev;
    magma_tally4_getdevice( &orig_dev );
    magma_tally4_queue_t orig_stream;
    magma_tally4blasGetKernelStream( &orig_stream );

    *info = 0;

    magma_tally4_int_t left   = (side == Magma_tally4Left);
    magma_tally4_int_t notran = (trans == Magma_tally4NoTrans);
    magma_tally4_int_t lquery = (lwork == -1);

    /* NQ is the order of Q and NW is the minimum dimension of WORK */
    magma_tally4_int_t nq, nw;
    if (left) {
        nq = m;
        nw = n;
    } else {
        nq = n;
        nw = m;
    }


    if (! left && side != Magma_tally4Right) {
        *info = -1;
    } else if (! notran && trans != Magma_tally4_ConjTrans) {
        *info = -2;
    } else if (m < 0) {
        *info = -3;
    } else if (n < 0) {
        *info = -4;
    } else if (k < 0 || k > nq) {
        *info = -5;
    } else if (lda < max(1,nq)) {
        *info = -7;
    } else if (ldc < max(1,m)) {
        *info = -10;
    } else if (lwork < max(1,nw) && ! lquery) {
        *info = -12;
    }

    magma_tally4_int_t lwkopt = max(1,nw) * nb;
    if (*info == 0) {
        work[0] = MAGMA_tally4_Z_MAKE( lwkopt, 0 );
    }

    if (*info != 0) {
        magma_tally4_xerbla( __func__, -(*info) );
        return *info;
    }
    else if (lquery) {
        return *info;
    }

    /* Quick return if possible */
    if (m == 0 || n == 0 || k == 0) {
        work[0] = c_one;
        return *info;
    }

    if (nb >= k) {
        /* Use CPU code */
        lapackf77_zunmqr(side_, trans_, &m, &n, &k, A, &lda, tau,
                         C, &ldc, work, &lwork, info);
        return *info;
    }

    magma_tally4_int_t lddc = (m+63)/64*64;
    magma_tally4_int_t lddac = nq;
    magma_tally4_int_t lddar = nb;
    magma_tally4_int_t lddwork = nw;

    magma_tally4_int_t nlocal[ Magma_tally4MaxGPUs ] = { 0 };

    magma_tally4_int_t nb_l=256;
    magma_tally4_int_t nbl = (n-1)/nb_l+1; // number of blocks
    magma_tally4_int_t maxnlocal = (nbl+ngpu-1)/ngpu*nb_l;

    ngpu = min(ngpu, (n+nb_l-1)/nb_l); // Don't use GPU that will not have data.

    magma_tally4_int_t ldw = maxnlocal*lddc // dC
                    + 2*lddac*lddar // 2*dA
                    + 2*(nb + 1 + lddwork)*nb; // 2*(dT and dwork)

    for (igpu = 0; igpu < ngpu; ++igpu) {
        magma_tally4_setdevice(igpu);
        if (MAGMA_tally4_SUCCESS != magma_tally4_zmalloc( &dw[igpu], ldw )) {
            *info = MAGMA_tally4_ERR_DEVICE_ALLOC;
            magma_tally4_xerbla( __func__, -(*info) );
            return *info;
        }
        magma_tally4_queue_create( &stream[igpu][0] );
        magma_tally4_queue_create( &stream[igpu][1] );
        magma_tally4_event_create( &event[igpu][0] );
        magma_tally4_event_create( &event[igpu][1] );
    }

    /* Use hybrid CPU-MGPU code */
    if (left) {
        //copy C to mgpus
        for (magma_tally4_int_t i = 0; i < nbl; ++i) {
            magma_tally4_int_t igpu = i%ngpu;
            magma_tally4_setdevice(igpu);
            magma_tally4_int_t kb = min(nb_l, n-i*nb_l);
            magma_tally4_zsetmatrix_async( m, kb,
                                   C(0, i*nb_l), ldc,
                                   dC(igpu, 0, i/ngpu*nb_l), lddc, stream[igpu][0] );
            nlocal[igpu] += kb;
        }

        magma_tally4_int_t i1, i2, i3;
        if ( !notran ) {
            i1 = 0;
            i2 = k;
            i3 = nb;
        } else {
            i1 = (k - 1) / nb * nb;
            i2 = 0;
            i3 = -nb;
        }

        ind_c = 0;

        for (magma_tally4_int_t i = i1; (i3 < 0 ? i >= i2 : i < i2); i += i3) {
            // start the copy of A panel
            magma_tally4_int_t kb = min(nb, k - i);
            for (igpu = 0; igpu < ngpu; ++igpu) {
                magma_tally4_setdevice(igpu);
                magma_tally4_event_sync(event[igpu][ind_c]); // check if the new data can be copied
                magma_tally4_zsetmatrix_async(nq-i, kb,
                                       A(i, i),                 lda,
                                       dA_c(igpu, ind_c, i, 0), lddac, stream[igpu][0] );
                // set upper triangular part of dA to identity
                magma_tally4blas_zlaset_band_q( Magma_tally4Upper, kb, kb, kb, c_zero, c_one, dA_c(igpu, ind_c, i, 0), lddac, stream[igpu][0] );
            }

            /* Form the triangular factor of the block reflector
             H = H(i) H(i+1) . . . H(i+ib-1) */
            magma_tally4_int_t nqi = nq - i;
            lapackf77_zlarft("F", "C", &nqi, &kb, A(i, i), &lda,
                             &tau[i], T, &kb);

            /* H or H' is applied to C(1:m,i:n) */

            /* Apply H or H'; First copy T to the GPU */
            for (igpu = 0; igpu < ngpu; ++igpu) {
                magma_tally4_setdevice(igpu);
                magma_tally4_zsetmatrix_async(kb, kb,
                                       T,               kb,
                                       dT(igpu, ind_c), kb, stream[igpu][0] );
            }

            for (igpu = 0; igpu < ngpu; ++igpu) {
                magma_tally4_setdevice(igpu);
                magma_tally4_queue_sync( stream[igpu][0] ); // check if the data was copied
                magma_tally4blasSetKernelStream(stream[igpu][1]);
                magma_tally4_zlarfb_gpu( side, trans, Magma_tally4Forward, Magma_tally4Columnwise,
                                 m-i, nlocal[igpu], kb,
                                 dA_c(igpu, ind_c, i, 0), lddac, dT(igpu, ind_c), kb,
                                 dC(igpu, i, 0), lddc,
                                 dwork(igpu, ind_c), lddwork);
                magma_tally4_event_record(event[igpu][ind_c], stream[igpu][1] );
            }

            ind_c = (ind_c+1)%2;
        }

        for (igpu = 0; igpu < ngpu; ++igpu) {
            magma_tally4_setdevice(igpu);
            magma_tally4_queue_sync( stream[igpu][1] );
        }

        //copy C from mgpus
        for (magma_tally4_int_t i = 0; i < nbl; ++i) {
            magma_tally4_int_t igpu = i%ngpu;
            magma_tally4_setdevice(igpu);
            magma_tally4_int_t kb = min(nb_l, n-i*nb_l);
            magma_tally4_zgetmatrix( m, kb,
                              dC(igpu, 0, i/ngpu*nb_l), lddc,
                              C(0, i*nb_l), ldc );
//            magma_tally4_zgetmatrix_async( m, kb,
//                                   dC(igpu, 0, i/ngpu*nb_l), lddc,
//                                   C(0, i*nb_l), ldc, stream[igpu][0] );
        }
    } else {
        // TODO fix memory leak T, dw, event, stream
        fprintf(stderr, "The case (side == right) is not implemented\n");
        *info = MAGMA_tally4_ERR_NOT_IMPLEMENTED;
        magma_tally4_xerbla( __func__, -(*info) );
        return *info;
        /*
        if ( notran ) {
            i1 = 0;
            i2 = k;
            i3 = nb;
        } else {
            i1 = (k - 1) / nb * nb;
            i2 = 0;
            i3 = -nb;
        }

        mi = m;
        ic = 0;

        for (i = i1; (i3 < 0 ? i >= i2 : i < i2); i += i3) {
            ib = min(nb, k - i);
            
            // Form the triangular factor of the block reflector
            // H = H(i) H(i+1) . . . H(i+ib-1)
            i__4 = nq - i;
            lapackf77_zlarft("F", "C", &i__4, &ib, A(i, i), &lda,
            &tau[i], T, &ib);
            
            // 1) copy the panel from A to the GPU, and
            // 2) set upper triangular part of dA to identity
            magma_tally4_zsetmatrix( i__4, ib, A(i, i), lda, dA(i, 0), ldda );
            magma_tally4blas_zlaset_band( Magma_tally4Upper, ib, ib, ib, c_zero, c_one, dA(i, 0), ldda );
            
            // H or H' is applied to C(1:m,i:n)
            ni = n - i;
            jc = i;
            
            // Apply H or H'; First copy T to the GPU
            magma_tally4_zsetmatrix( ib, ib, T, ib, dT, ib );
            magma_tally4_zlarfb_gpu( side, trans, Magma_tally4Forward, Magma_tally4Columnwise,
            mi, ni, ib,
            dA(i, 0), ldda, dT, ib,
            dC(ic, jc), lddc,
            dwork, lddwork);
        }
        */
    }

    work[0] = MAGMA_tally4_Z_MAKE( lwkopt, 0 );

    for (igpu = 0; igpu < ngpu; ++igpu) {
        magma_tally4_setdevice(igpu);
        magma_tally4_event_destroy( event[igpu][0] );
        magma_tally4_event_destroy( event[igpu][1] );
        magma_tally4_queue_destroy( stream[igpu][0] );
        magma_tally4_queue_destroy( stream[igpu][1] );
        magma_tally4_free( dw[igpu] );
    }
    magma_tally4_setdevice( orig_dev );
    magma_tally4blasSetKernelStream( orig_stream );

    return *info;
} /* magma_tally4_zunmqr */
