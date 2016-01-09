/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Stan Tomov
       @author Mark Gates

       @precisions normal z -> s d c

*/
#include "common_magma_tally4.h"

/**
    Purpose
    -------
    ZUNMQR overwrites the general complex M-by-N matrix C with

    @verbatim
                              SIDE = Magma_tally4Left   SIDE = Magma_tally4Right
    TRANS = Magma_tally4NoTrans:     Q * C              C * Q
    TRANS = Magma_tally4_ConjTrans:  Q**H * C           C * Q**H
    @endverbatim

    where Q is a complex unitary matrix defined as the product of k
    elementary reflectors

        Q = H(1) H(2) . . . H(k)

    as returned by ZGEQRF. Q is of order M if SIDE = Magma_tally4Left and of order N
    if SIDE = Magma_tally4Right.

    Arguments
    ---------
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
            A is modified by the routine but restored on exit.

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
            On exit, C is overwritten by Q*C or Q**H * C or C * Q**H or C*Q.

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
            For optimum performance
            if SIDE = Magma_tally4Left,  LWORK >= N*NB;
            if SIDE = Magma_tally4Right, LWORK >= M*NB,
            where NB is the optimal blocksize.
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
magma_tally4_zunmqr(
    magma_tally4_side_t side, magma_tally4_trans_t trans,
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4DoubleComplex *A,    magma_tally4_int_t lda,
    magma_tally4DoubleComplex *tau,
    magma_tally4DoubleComplex *C,    magma_tally4_int_t ldc,
    magma_tally4DoubleComplex *work, magma_tally4_int_t lwork,
    magma_tally4_int_t *info)
{
    #define  A(i_,j_) ( A + (i_) + (j_)*lda)
    #define dC(i_,j_) (dC + (i_) + (j_)*lddc)
    #define dV(i_,j_) (dV + (i_) + (j_)*nq_i)
    #define dT(i_,j_) (dT + (i_) + (j_)*ib)
    #define dwork(i_) (dwork + (i_))
    
    magma_tally4DoubleComplex *T, *T2;
    magma_tally4_int_t i, i1, i2, ib, ic, jc, nb, mi, ni, nq, nq_i, nw, step;
    magma_tally4_int_t iinfo, ldwork, lwkopt;
    magma_tally4_int_t left, notran, lquery;

    *info = 0;
    left   = (side == Magma_tally4Left);
    notran = (trans == Magma_tally4NoTrans);
    lquery = (lwork == -1);

    /* NQ is the order of Q and NW is the minimum dimension of WORK */
    if (left) {
        nq = m;
        nw = n;
    } else {
        nq = n;
        nw = m;
    }
    
    /* Test the input arguments */
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

    if (*info == 0) {
        nb = magma_tally4_get_zgelqf_nb( min( m, n ));
        lwkopt = max(1,nw)*nb;
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
        work[0] = MAGMA_tally4_Z_ONE;
        return *info;
    }

    ldwork = nw;

    if (nb >= k) {
        /* Use CPU code */
        lapackf77_zunmqr( lapack_side_const_tally4(side), lapack_trans_const_tally4(trans),
            &m, &n, &k, A, &lda, tau, C, &ldc, work, &lwork, &iinfo);
    }
    else {
        /* Use hybrid CPU-GPU code */
        /* Allocate work space on the GPU.
         * nw*nb  for dwork (m or n) by nb
         * nq*nb  for dV    (n or m) by nb
         * nb*nb  for dT
         * lddc*n for dC.
         */
        magma_tally4_int_t lddc = ((m+31)/32)*32;
        magma_tally4DoubleComplex_ptr dwork, dV, dT, dC;
        magma_tally4_zmalloc( &dwork, (nw + nq + nb)*nb + lddc*n );
        if ( dwork == NULL ) {
            *info = MAGMA_tally4_ERR_DEVICE_ALLOC;
            return *info;
        }
        dV = dwork + nw*nb;
        dT = dV    + nq*nb;
        dC = dT    + nb*nb;
        
        /* work space on CPU.
         * nb*nb for T
         * nb*nb for T2, used to save and restore diagonal block of panel */
        magma_tally4_zmalloc_cpu( &T, 2*nb*nb );
        if ( T == NULL ) {
            magma_tally4_free( dwork );
            *info = MAGMA_tally4_ERR_HOST_ALLOC;
            return *info;
        }
        T2 = T + nb*nb;
        
        /* Copy matrix C from the CPU to the GPU */
        magma_tally4_zsetmatrix( m, n, C, ldc, dC(0,0), lddc );
        
        if ( (left && ! notran) ||  (! left && notran) ) {
            i1 = 0;
            i2 = k;
            step = nb;
        } else {
            i1 = ((k - 1) / nb) * nb;
            i2 = 0;
            step = -nb;
        }

        // silence "uninitialized" warnings
        mi = 0;
        ni = 0;
        
        if (left) {
            ni = n;
            jc = 0;
        } else {
            mi = m;
            ic = 0;
        }
        
        for (i = i1; (step < 0 ? i >= i2 : i < i2); i += step) {
            ib = min(nb, k - i);

            /* Form the triangular factor of the block reflector
               H = H(i) H(i+1) . . . H(i+ib-1) */
            nq_i = nq - i;
            lapackf77_zlarft("Forward", "Columnwise", &nq_i, &ib,
                             A(i,i), &lda, &tau[i], T, &ib);

            /* 1) set upper triangle of panel in A to identity,
               2) copy the panel from A to the GPU, and
               3) restore A                                      */
            zpanel_to_q_tally4( Magma_tally4Upper, ib, A(i,i), lda, T2 );
            magma_tally4_zsetmatrix( nq_i,  ib, A(i,i), lda, dV(0,0), nq_i );
            zq_to_panel_tally4( Magma_tally4Upper, ib, A(i,i), lda, T2 );

            if (left) {
                /* H or H**H is applied to C(i:m,1:n) */
                mi = m - i;
                ic = i;
            }
            else {
                /* H or H**H is applied to C(1:m,i:n) */
                ni = n - i;
                jc = i;
            }

            /* Apply H or H**H; First copy T to the GPU */
            magma_tally4_zsetmatrix( ib, ib, T, ib, dT(0,0), ib );
            magma_tally4_zlarfb_gpu( side, trans, Magma_tally4Forward, Magma_tally4Columnwise,
                              mi, ni, ib,
                              dV(0,0), nq_i,
                              dT(0,0), ib,
                              dC(ic,jc), lddc,
                              dwork(0), ldwork );
        }
        magma_tally4_zgetmatrix( m, n, dC(0,0), lddc, C, ldc );

        magma_tally4_free( dwork );
        magma_tally4_free_cpu( T );
    }
    work[0] = MAGMA_tally4_Z_MAKE( lwkopt, 0 );

    return *info;
} /* magma_tally4_zunmqr */
