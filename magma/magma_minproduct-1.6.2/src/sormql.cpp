/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Raffaele Solca
       @author Mark Gates

       @generated from zunmql.cpp normal z -> s, Fri Jan 30 19:00:16 2015

*/
#include "common_magma_minproduct.h"

/**
    Purpose
    -------
    SORMQL overwrites the general real M-by-N matrix C with

    @verbatim
                              SIDE = Magma_minproductLeft   SIDE = Magma_minproductRight
    TRANS = Magma_minproductNoTrans:     Q * C              C * Q
    TRANS = Magma_minproductTrans:  Q**H * C           C * Q**H
    @endverbatim

    where Q is a real unitary matrix defined as the product of k
    elementary reflectors

          Q = H(k) . . . H(2) H(1)

    as returned by SGEQLF. Q is of order M if SIDE = Magma_minproductLeft and of order N
    if SIDE = Magma_minproductRight.

    Arguments
    ---------
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
    A       REAL array, dimension (LDA,K)
            The i-th column must contain the vector which defines the
            elementary reflector H(i), for i = 1,2,...,k, as returned by
            SGEQLF in the last k columns of its array argument A.
            A is modified by the routine but restored on exit.

    @param[in]
    lda     INTEGER
            The leading dimension of the array A.
            If SIDE = Magma_minproductLeft,  LDA >= max(1,M);
            if SIDE = Magma_minproductRight, LDA >= max(1,N).

    @param[in]
    tau     REAL array, dimension (K)
            TAU(i) must contain the scalar factor of the elementary
            reflector H(i), as returned by SGEQLF.

    @param[in,out]
    C       REAL array, dimension (LDC,N)
            On entry, the M-by-N matrix C.
            On exit, C is overwritten by Q*C or Q**H*C or C*Q**H or C*Q.

    @param[in]
    ldc     INTEGER
            The leading dimension of the array C. LDC >= max(1,M).

    @param[out]
    work    (workspace) REAL array, dimension (MAX(1,LWORK))
            On exit, if INFO = 0, WORK[0] returns the optimal LWORK.

    @param[in]
    lwork   INTEGER
            The dimension of the array WORK.
            If SIDE = Magma_minproductLeft,  LWORK >= max(1,N);
            if SIDE = Magma_minproductRight, LWORK >= max(1,M).
            For optimum performance
            if SIDE = Magma_minproductLeft,  LWORK >= N*NB;
            if SIDE = Magma_minproductRight, LWORK >= M*NB,
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

    @ingroup magma_minproduct_sgeqlf_comp
    ********************************************************************/
extern "C" magma_minproduct_int_t
magma_minproduct_sormql(
    magma_minproduct_side_t side, magma_minproduct_trans_t trans,
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    float *A, magma_minproduct_int_t lda,
    float *tau,
    float *C, magma_minproduct_int_t ldc,
    float *work, magma_minproduct_int_t lwork,
    magma_minproduct_int_t *info)
{
    #define  A(i_,j_) ( A + (i_) + (j_)*lda)
    #define dC(i_,j_) (dC + (i_) + (j_)*lddc)
    
    float *T, *T2;
    magma_minproduct_int_t i, i1, i2, ib, nb, mi, ni, nq, nq_i, nw, step;
    magma_minproduct_int_t iinfo, ldwork, lwkopt;
    magma_minproduct_int_t left, notran, lquery;

    *info  = 0;
    left   = (side == Magma_minproductLeft);
    notran = (trans == Magma_minproductNoTrans);
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

    if (*info == 0) {
        nb = magma_minproduct_get_sgelqf_nb( min( m, n ));
        lwkopt = max(1,nw)*nb;
        work[0] = MAGMA_minproduct_S_MAKE( lwkopt, 0 );
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
        work[0] = MAGMA_minproduct_S_ONE;
        return *info;
    }

    ldwork = nw;

    if ( nb >= k ) {
        /* Use CPU code */
        lapackf77_sormql( lapack_side_const(side), lapack_trans_const(trans),
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
        magma_minproduct_int_t lddc = ((m+31)/32)*32;
        float *dwork, *dV, *dT, *dC;
        magma_minproduct_smalloc( &dwork, (nw + nq + nb)*nb + lddc*n );
        if ( dwork == NULL ) {
            *info = MAGMA_minproduct_ERR_DEVICE_ALLOC;
            return *info;
        }
        dV = dwork + nw*nb;
        dT = dV    + nq*nb;
        dC = dT    + nb*nb;
        
        /* work space on CPU.
         * nb*nb for T
         * nb*nb for T2, used to save and restore diagonal block of panel */
        magma_minproduct_smalloc_pinned( &T, 2*nb*nb );
        if ( T == NULL ) {
            magma_minproduct_free( dwork );
            *info = MAGMA_minproduct_ERR_HOST_ALLOC;
            return *info;
        }
        T2 = T + nb*nb;
    
        /* Copy matrix C from the CPU to the GPU */
        magma_minproduct_ssetmatrix( m, n, C, ldc, dC, lddc );
        
        if ( (left && notran) || (! left && ! notran) ) {
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
        } else {
            mi = m;
        }

        for (i = i1; (step < 0 ? i >= i2 : i < i2); i += step) {
            ib = min(nb, k - i);
            
            /* Form the triangular factor of the block reflector
               H = H(i+ib-1) . . . H(i+1) H(i) */
            nq_i = nq - k + i + ib;
            lapackf77_slarft("Backward", "Columnwise", &nq_i, &ib,
                             A(0,i), &lda, &tau[i], T, &ib);
            
            /* 1) set lower triangle of panel in A to identity,
               2) copy the panel from A to the GPU, and
               3) restore A                                      */
            spanel_to_q( Magma_minproductLower, ib, A(nq_i-ib,i), lda, T2 );
            magma_minproduct_ssetmatrix( nq_i,  ib, A(0,      i), lda, dV, nq_i );
            sq_to_panel( Magma_minproductLower, ib, A(nq_i-ib,i), lda, T2 );
            
            if (left) {
                /* H or H**H is applied to C(1:m-k+i+ib-1,1:n) */
                mi = m - k + i + ib;
            }
            else {
                /* H or H**H is applied to C(1:m,1:n-k+i+ib-1) */
                ni = n - k + i + ib;
            }
            
            /* Apply H or H**H; First copy T to the GPU */
            magma_minproduct_ssetmatrix( ib, ib, T, ib, dT, ib );
            magma_minproduct_slarfb_gpu( side, trans, Magma_minproductBackward, Magma_minproductColumnwise,
                              mi, ni, ib,
                              dV, nq_i,
                              dT, ib,
                              dC, lddc,
                              dwork, ldwork );
        }
        magma_minproduct_sgetmatrix( m, n, dC, lddc, C, ldc );

        magma_minproduct_free( dwork );
        magma_minproduct_free_pinned( T );
    }
    work[0] = MAGMA_minproduct_S_MAKE( lwkopt, 0 );

    return *info;
} /* magma_minproduct_sormql */
