/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Mark Gates

       @precisions normal z -> s d c

*/
#include "common_magma_tally3.h"

/*
 * Version 1 - LAPACK
 * Version 2 - MAGMA_tally3
 */
#define VERSION 2

/**
    Purpose
    -------
    If VECT = Magma_tally3Q, ZUNMBR overwrites the general complex M-by-N matrix C with
                                 SIDE = Magma_tally3Left     SIDE = Magma_tally3Right
    TRANS = Magma_tally3NoTrans:        Q*C                  C*Q
    TRANS = Magma_tally3_ConjTrans:     Q**H*C               C*Q**H
    
    If VECT = Magma_tally3P, ZUNMBR overwrites the general complex M-by-N matrix C with
                                 SIDE = Magma_tally3Left     SIDE = Magma_tally3Right
    TRANS = Magma_tally3NoTrans:        P*C                  C*P
    TRANS = Magma_tally3_ConjTrans:     P**H*C               C*P**H
    
    Here Q and P**H are the unitary matrices determined by ZGEBRD when
    reducing A complex matrix A to bidiagonal form: A = Q*B * P**H. Q
    and P**H are defined as products of elementary reflectors H(i) and
    G(i) respectively.
    
    Let nq = m if SIDE = Magma_tally3Left and nq = n if SIDE = Magma_tally3Right. Thus nq is the
    order of the unitary matrix Q or P**H that is applied.
    
    If VECT = Magma_tally3Q, A is assumed to have been an NQ-by-K matrix:
    if nq >= k, Q = H(1) H(2) . . . H(k);
    if nq <  k, Q = H(1) H(2) . . . H(nq-1).
    
    If VECT = Magma_tally3P, A is assumed to have been A K-by-NQ matrix:
    if k <  nq, P = G(1) G(2) . . . G(k);
    if k >= nq, P = G(1) G(2) . . . G(nq-1).
    
    Arguments
    ---------
    @param[in]
    vect    magma_tally3_vect_t
      -     = Magma_tally3Q: apply Q or Q**H;
      -     = Magma_tally3P: apply P or P**H.
    
    @param[in]
    side    magma_tally3_side_t
      -     = Magma_tally3Left:  apply Q, Q**H, P or P**H from the Left;
      -     = Magma_tally3Right: apply Q, Q**H, P or P**H from the Right.
    
    @param[in]
    trans   magma_tally3_trans_t
      -     = Magma_tally3NoTrans:    No transpose, apply Q or P;
      -     = Magma_tally3_ConjTrans: Conjugate transpose, apply Q**H or P**H.
    
    @param[in]
    m       INTEGER
            The number of rows of the matrix C. M >= 0.
    
    @param[in]
    n       INTEGER
            The number of columns of the matrix C. N >= 0.
    
    @param[in]
    k       INTEGER
            If VECT = Magma_tally3Q, the number of columns in the original
            matrix reduced by ZGEBRD.
            If VECT = Magma_tally3P, the number of rows in the original
            matrix reduced by ZGEBRD.
            K >= 0.
    
    @param[in]
    A       COMPLEX_16 array, dimension
                                  (LDA,min(nq,K)) if VECT = Magma_tally3Q
                                  (LDA,nq)        if VECT = Magma_tally3P
            The vectors which define the elementary reflectors H(i) and
            G(i), whose products determine the matrices Q and P, as
            returned by ZGEBRD.
    
    @param[in]
    lda     INTEGER
            The leading dimension of the array A.
            If VECT = Magma_tally3Q, LDA >= max(1,nq);
            if VECT = Magma_tally3P, LDA >= max(1,min(nq,K)).
    
    @param[in]
    tau     COMPLEX_16 array, dimension (min(nq,K))
            TAU(i) must contain the scalar factor of the elementary
            reflector H(i) or G(i) which determines Q or P, as returned
            by ZGEBRD in the array argument TAUQ or TAUP.
    
    @param[in,out]
    C       COMPLEX_16 array, dimension (LDC,N)
            On entry, the M-by-N matrix C.
            On exit, C is overwritten by Q*C or Q**H*C or C*Q**H or C*Q
            or P*C or P**H*C or C*P or C*P**H.
    
    @param[in]
    ldc     INTEGER
            The leading dimension of the array C. LDC >= max(1,M).
    
    @param[out]
    work    (workspace) COMPLEX_16 array, dimension (MAX(1,LWORK))
            On exit, if INFO = 0, WORK[0] returns the optimal LWORK.
    
    @param[in]
    lwork   INTEGER
            The dimension of the array WORK.
            If SIDE = Magma_tally3Left,  LWORK >= max(1,N);
            if SIDE = Magma_tally3Right, LWORK >= max(1,M);
            if N = 0 or M = 0, LWORK >= 1.
            For optimum performance
            if SIDE = Magma_tally3Left,  LWORK >= max(1,N*NB);
            if SIDE = Magma_tally3Right, LWORK >= max(1,M*NB),
            where NB is the optimal blocksize. (NB = 0 if M = 0 or N = 0.)
    \n
            If LWORK = -1, then a workspace query is assumed; the routine
            only calculates the optimal size of the WORK array, returns
            this value as the first entry of the WORK array, and no error
            message related to LWORK is issued by XERBLA.
    
    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument had an illegal value
    

    @ingroup magma_tally3_zgesvd_comp
    ********************************************************************/
extern "C" magma_tally3_int_t
magma_tally3_zunmbr(
    magma_tally3_vect_t vect, magma_tally3_side_t side, magma_tally3_trans_t trans,
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    magma_tally3DoubleComplex *A, magma_tally3_int_t lda,
    magma_tally3DoubleComplex *tau,
    magma_tally3DoubleComplex *C, magma_tally3_int_t ldc,
    magma_tally3DoubleComplex *work, magma_tally3_int_t lwork,
    magma_tally3_int_t *info)
{
    #define A(i,j)  (A + (i) + (j)*lda)
    #define C(i,j)  (C + (i) + (j)*ldc)
            
    magma_tally3_int_t i1, i2, nb, mi, ni, nq, nq_1, nw, iinfo, lwkopt;
    magma_tally3_int_t left, notran, applyq, lquery;
    magma_tally3_trans_t transt;
    
    MAGMA_tally3_UNUSED( nq_1 );  // used only in version 1

    *info = 0;
    applyq = (vect  == Magma_tally3Q);
    left   = (side  == Magma_tally3Left);
    notran = (trans == Magma_tally3NoTrans);
    lquery = (lwork == -1);

    /* NQ is the order of Q or P and NW is the minimum dimension of WORK */
    if (left) {
        nq = m;
        nw = n;
    }
    else {
        nq = n;
        nw = m;
    }
    if (m == 0 || n == 0) {
        nw = 0;
    }
    
    /* check arguments */
    if (! applyq && vect != Magma_tally3P) {
        *info = -1;
    }
    else if (! left && side != Magma_tally3Right) {
        *info = -2;
    }
    else if (! notran && trans != Magma_tally3_ConjTrans) {
        *info = -3;
    }
    else if (m < 0) {
        *info = -4;
    }
    else if (n < 0) {
        *info = -5;
    }
    else if (k < 0) {
        *info = -6;
    }
    else if ( (   applyq && lda < max(1,nq)        ) ||
              ( ! applyq && lda < max(1,min(nq,k)) ) ) {
        *info = -8;
    }
    else if (ldc < max(1,m)) {
        *info = -11;
    }
    else if (lwork < max(1,nw) && ! lquery) {
        *info = -13;
    }

    if (*info == 0) {
        if (nw > 0) {
            // TODO have get_zunmqr_nb and get_zunmlq_nb routines? see original LAPACK zunmbr.
            // TODO make them dependent on m, n, and k?
            nb = magma_tally3_get_zgebrd_nb( min( m, n ));
            lwkopt = max(1, nw*nb);
        }
        else {
            lwkopt = 1;
        }
        work[0] = MAGMA_tally3_Z_MAKE( lwkopt, 0 );
    }

    if (*info != 0) {
        magma_tally3_xerbla( __func__, -(*info) );
        return *info;
    }
    else if (lquery) {
        return *info;
    }

    /* Quick return if possible */
    if (m == 0 || n == 0) {
        return *info;
    }

    if (applyq) {
        /* Apply Q */
        if (nq >= k) {
            /* Q was determined by a call to ZGEBRD with nq >= k */
            #if VERSION == 1
            lapackf77_zunmqr( lapack_side_const_tally3(side), lapack_trans_const_tally3(trans),
                              &m, &n, &k, A, &lda, tau, C, &ldc, work, &lwork, &iinfo);
            #else
            magma_tally3_zunmqr( side, trans,
                          m, n, k, A, lda, tau, C, ldc, work, lwork, &iinfo);
            #endif
        }
        else if (nq > 1) {
            /* Q was determined by a call to ZGEBRD with nq < k */
            if (left) {
                mi = m - 1;
                ni = n;
                i1 = 1;
                i2 = 0;
            }
            else {
                mi = m;
                ni = n - 1;
                i1 = 0;
                i2 = 1;
            }
            #if VERSION == 1
            nq_1 = nq - 1;
            lapackf77_zunmqr( lapack_side_const_tally3(side), lapack_trans_const_tally3(trans),
                              &mi, &ni, &nq_1, A(1,0), &lda, tau, C(i1,i2), &ldc, work, &lwork, &iinfo);
            #else
            magma_tally3_zunmqr( side, trans,
                          mi, ni, nq-1, A(1,0), lda, tau, C(i1,i2), ldc, work, lwork, &iinfo);
            #endif
        }
    }
    else {
        /* Apply P */
        if (notran) {
            transt = Magma_tally3_ConjTrans;
        }
        else {
            transt = Magma_tally3NoTrans;
        }
        if (nq > k) {
            /* P was determined by a call to ZGEBRD with nq > k */
            #if VERSION == 1
            lapackf77_zunmlq( lapack_side_const_tally3(side), lapack_trans_const_tally3(transt),
                              &m, &n, &k, A, &lda, tau, C, &ldc, work, &lwork, &iinfo);
            #else
            magma_tally3_zunmlq( side, transt,
                          m, n, k, A, lda, tau, C, ldc, work, lwork, &iinfo);
            #endif
        }
        else if (nq > 1) {
            /* P was determined by a call to ZGEBRD with nq <= k */
            if (left) {
                mi = m - 1;
                ni = n;
                i1 = 1;
                i2 = 0;
            }
            else {
                mi = m;
                ni = n - 1;
                i1 = 0;
                i2 = 1;
            }
            #if VERSION == 1
            nq_1 = nq - 1;
            lapackf77_zunmlq( lapack_side_const_tally3(side), lapack_trans_const_tally3(transt),
                              &mi, &ni, &nq_1, A(0,1), &lda, tau, C(i1,i2), &ldc, work, &lwork, &iinfo);
            #else
            magma_tally3_zunmlq( side, transt,
                          mi, ni, nq-1, A(0,1), lda, tau, C(i1,i2), ldc, work, lwork, &iinfo);
            #endif
        }
    }
    work[0] = MAGMA_tally3_Z_MAKE( lwkopt, 0 );
    return *info;
} /* magma_tally3_zunmbr */
