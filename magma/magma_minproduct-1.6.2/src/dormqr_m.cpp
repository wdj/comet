/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Raffaele Solca
       @author Azzam Haidar
       @author Stan Tomov

       @generated from zunmqr_m.cpp normal z -> d, Fri Jan 30 19:00:16 2015

*/
#include "common_magma_minproduct.h"

/**
    Purpose
    -------
    DORMQR overwrites the general real M-by-N matrix C with

    @verbatim
                                SIDE = Magma_minproductLeft    SIDE = Magma_minproductRight
    TRANS = Magma_minproductNoTrans:       Q * C               C * Q
    TRANS = Magma_minproductTrans:    Q**H * C            C * Q**H
    @endverbatim

    where Q is a real unitary matrix defined as the product of k
    elementary reflectors

          Q = H(1) H(2) . . . H(k)

    as returned by DGEQRF. Q is of order M if SIDE = Magma_minproductLeft and of order N
    if SIDE = Magma_minproductRight.

    Arguments
    ---------
    @param[in]
    ngpu    INTEGER
            Number of GPUs to use. ngpu > 0.

    @param[in]
    side    magma_minproduct_side_t
      -     = Magma_minproductLeft:      apply Q or Q**H from the Left;
      -     = Magma_minproductRight:     apply Q or Q**H from the Right.

    @param[in]
    trans   magma_minproduct_trans_t
      -     = Magma_minproductNoTrans:    No transpose, apply Q;
      -     = Magma_minproductTrans: Conjugate transpose, apply Q**H.

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
            If SIDE = Magma_minproductLeft,  M >= K >= 0;
            if SIDE = Magma_minproductRight, N >= K >= 0.

    @param[in]
    A       DOUBLE_PRECISION array, dimension (LDA,K)
            The i-th column must contain the vector which defines the
            elementary reflector H(i), for i = 1,2,...,k, as returned by
            DGEQRF in the first k columns of its array argument A.

    @param[in]
    lda     INTEGER
            The leading dimension of the array A.
            If SIDE = Magma_minproductLeft,  LDA >= max(1,M);
            if SIDE = Magma_minproductRight, LDA >= max(1,N).

    @param[in]
    tau     DOUBLE_PRECISION array, dimension (K)
            TAU(i) must contain the scalar factor of the elementary
            reflector H(i), as returned by DGEQRF.

    @param[in,out]
    C       DOUBLE_PRECISION array, dimension (LDC,N)
            On entry, the M-by-N matrix C.
            On exit, C is overwritten by Q*C or Q**H*C or C*Q**H or C*Q.

    @param[in]
    ldc     INTEGER
            The leading dimension of the array C. LDC >= max(1,M).

    @param[out]
    work    (workspace) DOUBLE_PRECISION array, dimension (MAX(1,LWORK))
            On exit, if INFO = 0, WORK[0] returns the optimal LWORK.

    @param[in]
    lwork   INTEGER
            The dimension of the array WORK.
            If SIDE = Magma_minproductLeft,  LWORK >= max(1,N);
            if SIDE = Magma_minproductRight, LWORK >= max(1,M).
            For optimum performance LWORK >= N*NB if SIDE = Magma_minproductLeft, and
            LWORK >= M*NB if SIDE = Magma_minproductRight, where NB is the optimal
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

    @ingroup magma_minproduct_dgeqrf_comp
    ********************************************************************/
extern "C" magma_minproduct_int_t
magma_minproduct_dormqr_m(
    magma_minproduct_int_t ngpu,
    magma_minproduct_side_t side, magma_minproduct_trans_t trans,
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    double *A,    magma_minproduct_int_t lda,
    double *tau,
    double *C,    magma_minproduct_int_t ldc,
    double *work, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info)
{
#define  A(i, j) (A + (j)*lda  + (i))
#define  C(i, j) (C + (j)*ldc  + (i))

#define    dC(gpui,      i, j) (dw[gpui] + (j)*lddc + (i))
#define  dA_c(gpui, ind, i, j) (dw[gpui] + maxnlocal*lddc + (ind)*lddar*lddac + (i) + (j)*lddac)
#define  dA_r(gpui, ind, i, j) (dw[gpui] + maxnlocal*lddc + (ind)*lddar*lddac + (i) + (j)*lddar)
#define    dT(gpui, ind)       (dw[gpui] + maxnlocal*lddc + 2*lddac*lddar + (ind)*((nb+1)*nb))
#define dwork(gpui, ind)       (dw[gpui] + maxnlocal*lddc + 2*lddac*lddar + 2*((nb+1)*nb) + (ind)*(lddwork*nb))

    double c_zero = MAGMA_minproduct_D_ZERO;
    double c_one  = MAGMA_minproduct_D_ONE;

    const char* side_  = lapack_side_const( side );
    const char* trans_ = lapack_trans_const( trans );

    // TODO fix memory leak (alloc after argument checks)
    magma_minproduct_int_t nb = 128;
    double *T;
    magma_minproduct_dmalloc_pinned(&T, nb*nb);
    //printf("calling dormqr_m with nb=%d\n", (int) nb);

    double* dw[Magma_minproductMaxGPUs];
    magma_minproduct_queue_t stream [Magma_minproductMaxGPUs][2];
    magma_minproduct_event_t  event [Magma_minproductMaxGPUs][2];

    magma_minproduct_int_t ind_c;
    magma_minproduct_device_t igpu;
    
    magma_minproduct_device_t orig_dev;
    magma_minproduct_getdevice( &orig_dev );
    magma_minproduct_queue_t orig_stream;
    magma_minproductblasGetKernelStream( &orig_stream );

    *info = 0;

    magma_minproduct_int_t left   = (side == Magma_minproductLeft);
    magma_minproduct_int_t notran = (trans == Magma_minproductNoTrans);
    magma_minproduct_int_t lquery = (lwork == -1);

    /* NQ is the order of Q and NW is the minimum dimension of WORK */
    magma_minproduct_int_t nq, nw;
    if (left) {
        nq = m;
        nw = n;
    } else {
        nq = n;
        nw = m;
    }


    if (! left && side != Magma_minproductRight) {
        *info = -1;
    } else if (! notran && trans != Magma_minproductTrans) {
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

    magma_minproduct_int_t lwkopt = max(1,nw) * nb;
    if (*info == 0) {
        work[0] = MAGMA_minproduct_D_MAKE( lwkopt, 0 );
    }

    if (*info != 0) {
        magma_minproduct_xerbla( __func__, -(*info) );
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
        lapackf77_dormqr(side_, trans_, &m, &n, &k, A, &lda, tau,
                         C, &ldc, work, &lwork, info);
        return *info;
    }

    magma_minproduct_int_t lddc = (m+63)/64*64;
    magma_minproduct_int_t lddac = nq;
    magma_minproduct_int_t lddar = nb;
    magma_minproduct_int_t lddwork = nw;

    magma_minproduct_int_t nlocal[ Magma_minproductMaxGPUs ] = { 0 };

    magma_minproduct_int_t nb_l=256;
    magma_minproduct_int_t nbl = (n-1)/nb_l+1; // number of blocks
    magma_minproduct_int_t maxnlocal = (nbl+ngpu-1)/ngpu*nb_l;

    ngpu = min(ngpu, (n+nb_l-1)/nb_l); // Don't use GPU that will not have data.

    magma_minproduct_int_t ldw = maxnlocal*lddc // dC
                    + 2*lddac*lddar // 2*dA
                    + 2*(nb + 1 + lddwork)*nb; // 2*(dT and dwork)

    for (igpu = 0; igpu < ngpu; ++igpu) {
        magma_minproduct_setdevice(igpu);
        if (MAGMA_minproduct_SUCCESS != magma_minproduct_dmalloc( &dw[igpu], ldw )) {
            *info = MAGMA_minproduct_ERR_DEVICE_ALLOC;
            magma_minproduct_xerbla( __func__, -(*info) );
            return *info;
        }
        magma_minproduct_queue_create( &stream[igpu][0] );
        magma_minproduct_queue_create( &stream[igpu][1] );
        magma_minproduct_event_create( &event[igpu][0] );
        magma_minproduct_event_create( &event[igpu][1] );
    }

    /* Use hybrid CPU-MGPU code */
    if (left) {
        //copy C to mgpus
        for (magma_minproduct_int_t i = 0; i < nbl; ++i) {
            magma_minproduct_int_t igpu = i%ngpu;
            magma_minproduct_setdevice(igpu);
            magma_minproduct_int_t kb = min(nb_l, n-i*nb_l);
            magma_minproduct_dsetmatrix_async( m, kb,
                                   C(0, i*nb_l), ldc,
                                   dC(igpu, 0, i/ngpu*nb_l), lddc, stream[igpu][0] );
            nlocal[igpu] += kb;
        }

        magma_minproduct_int_t i1, i2, i3;
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

        for (magma_minproduct_int_t i = i1; (i3 < 0 ? i >= i2 : i < i2); i += i3) {
            // start the copy of A panel
            magma_minproduct_int_t kb = min(nb, k - i);
            for (igpu = 0; igpu < ngpu; ++igpu) {
                magma_minproduct_setdevice(igpu);
                magma_minproduct_event_sync(event[igpu][ind_c]); // check if the new data can be copied
                magma_minproduct_dsetmatrix_async(nq-i, kb,
                                       A(i, i),                 lda,
                                       dA_c(igpu, ind_c, i, 0), lddac, stream[igpu][0] );
                // set upper triangular part of dA to identity
                magma_minproductblas_dlaset_band_q( Magma_minproductUpper, kb, kb, kb, c_zero, c_one, dA_c(igpu, ind_c, i, 0), lddac, stream[igpu][0] );
            }

            /* Form the triangular factor of the block reflector
             H = H(i) H(i+1) . . . H(i+ib-1) */
            magma_minproduct_int_t nqi = nq - i;
            lapackf77_dlarft("F", "C", &nqi, &kb, A(i, i), &lda,
                             &tau[i], T, &kb);

            /* H or H' is applied to C(1:m,i:n) */

            /* Apply H or H'; First copy T to the GPU */
            for (igpu = 0; igpu < ngpu; ++igpu) {
                magma_minproduct_setdevice(igpu);
                magma_minproduct_dsetmatrix_async(kb, kb,
                                       T,               kb,
                                       dT(igpu, ind_c), kb, stream[igpu][0] );
            }

            for (igpu = 0; igpu < ngpu; ++igpu) {
                magma_minproduct_setdevice(igpu);
                magma_minproduct_queue_sync( stream[igpu][0] ); // check if the data was copied
                magma_minproductblasSetKernelStream(stream[igpu][1]);
                magma_minproduct_dlarfb_gpu( side, trans, Magma_minproductForward, Magma_minproductColumnwise,
                                 m-i, nlocal[igpu], kb,
                                 dA_c(igpu, ind_c, i, 0), lddac, dT(igpu, ind_c), kb,
                                 dC(igpu, i, 0), lddc,
                                 dwork(igpu, ind_c), lddwork);
                magma_minproduct_event_record(event[igpu][ind_c], stream[igpu][1] );
            }

            ind_c = (ind_c+1)%2;
        }

        for (igpu = 0; igpu < ngpu; ++igpu) {
            magma_minproduct_setdevice(igpu);
            magma_minproduct_queue_sync( stream[igpu][1] );
        }

        //copy C from mgpus
        for (magma_minproduct_int_t i = 0; i < nbl; ++i) {
            magma_minproduct_int_t igpu = i%ngpu;
            magma_minproduct_setdevice(igpu);
            magma_minproduct_int_t kb = min(nb_l, n-i*nb_l);
            magma_minproduct_dgetmatrix( m, kb,
                              dC(igpu, 0, i/ngpu*nb_l), lddc,
                              C(0, i*nb_l), ldc );
//            magma_minproduct_dgetmatrix_async( m, kb,
//                                   dC(igpu, 0, i/ngpu*nb_l), lddc,
//                                   C(0, i*nb_l), ldc, stream[igpu][0] );
        }
    } else {
        // TODO fix memory leak T, dw, event, stream
        fprintf(stderr, "The case (side == right) is not implemented\n");
        *info = MAGMA_minproduct_ERR_NOT_IMPLEMENTED;
        magma_minproduct_xerbla( __func__, -(*info) );
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
            lapackf77_dlarft("F", "C", &i__4, &ib, A(i, i), &lda,
            &tau[i], T, &ib);
            
            // 1) copy the panel from A to the GPU, and
            // 2) set upper triangular part of dA to identity
            magma_minproduct_dsetmatrix( i__4, ib, A(i, i), lda, dA(i, 0), ldda );
            magma_minproductblas_dlaset_band( Magma_minproductUpper, ib, ib, ib, c_zero, c_one, dA(i, 0), ldda );
            
            // H or H' is applied to C(1:m,i:n)
            ni = n - i;
            jc = i;
            
            // Apply H or H'; First copy T to the GPU
            magma_minproduct_dsetmatrix( ib, ib, T, ib, dT, ib );
            magma_minproduct_dlarfb_gpu( side, trans, Magma_minproductForward, Magma_minproductColumnwise,
            mi, ni, ib,
            dA(i, 0), ldda, dT, ib,
            dC(ic, jc), lddc,
            dwork, lddwork);
        }
        */
    }

    work[0] = MAGMA_minproduct_D_MAKE( lwkopt, 0 );

    for (igpu = 0; igpu < ngpu; ++igpu) {
        magma_minproduct_setdevice(igpu);
        magma_minproduct_event_destroy( event[igpu][0] );
        magma_minproduct_event_destroy( event[igpu][1] );
        magma_minproduct_queue_destroy( stream[igpu][0] );
        magma_minproduct_queue_destroy( stream[igpu][1] );
        magma_minproduct_free( dw[igpu] );
    }
    magma_minproduct_setdevice( orig_dev );
    magma_minproductblasSetKernelStream( orig_stream );

    return *info;
} /* magma_minproduct_dormqr */
