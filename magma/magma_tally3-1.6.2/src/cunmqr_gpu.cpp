/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Stan Tomov
       @author Mark Gates

       @generated from zunmqr_gpu.cpp normal z -> c, Fri Jan 30 19:00:15 2015

*/
#include "common_magma_tally3.h"

/**
    Purpose
    -------
    CUNMQR_GPU overwrites the general complex M-by-N matrix C with

    @verbatim
                               SIDE = Magma_tally3Left    SIDE = Magma_tally3Right
    TRANS = Magma_tally3NoTrans:      Q * C               C * Q
    TRANS = Magma_tally3_ConjTrans:   Q**H * C            C * Q**H
    @endverbatim

    where Q is a complex unitary matrix defined as the product of k
    elementary reflectors

          Q = H(1) H(2) . . . H(k)

    as returned by CGEQRF. Q is of order M if SIDE = Magma_tally3Left and of order N
    if SIDE = Magma_tally3Right.

    Arguments
    ---------
    @param[in]
    side    magma_tally3_side_t
      -     = Magma_tally3Left:   apply Q or Q**H from the Left;
      -     = Magma_tally3Right:  apply Q or Q**H from the Right.

    @param[in]
    trans   magma_tally3_trans_t
      -     = Magma_tally3NoTrans:    No transpose, apply Q;
      -     = Magma_tally3_ConjTrans: Conjugate transpose, apply Q**H.

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
            If SIDE = Magma_tally3Left,  M >= K >= 0;
            if SIDE = Magma_tally3Right, N >= K >= 0.

    @param[in]
    dA      COMPLEX array on the GPU, dimension (LDDA,K)
            The i-th column must contain the vector which defines the
            elementary reflector H(i), for i = 1,2,...,k, as returned by
            CGEQRF in the first k columns of its array argument DA.
            DA is modified by the routine but restored on exit.

    @param[in]
    ldda    INTEGER
            The leading dimension of the array DA.
            If SIDE = Magma_tally3Left,  LDDA >= max(1,M);
            if SIDE = Magma_tally3Right, LDDA >= max(1,N).

    @param[in]
    tau     COMPLEX array, dimension (K)
            TAU(i) must contain the scalar factor of the elementary
            reflector H(i), as returned by CGEQRF.

    @param[in,out]
    dC      COMPLEX array on the GPU, dimension (LDDC,N)
            On entry, the M-by-N matrix C.
            On exit, C is overwritten by Q*C or Q**H * C or C * Q**H or C*Q.

    @param[in]
    lddc    INTEGER
            The leading dimension of the array DC. LDDC >= max(1,M).

    @param[out]
    hwork   (workspace) COMPLEX array, dimension (MAX(1,LWORK))
    \n
            Currently, cgetrs_gpu assumes that on exit, hwork contains the last
            block of A and C. This will change and *should not be relied on*!

    @param[in]
    lwork   INTEGER
            The dimension of the array HWORK.
            LWORK >= (M-K+NB)*(N+NB) + N*NB if SIDE = Magma_tally3Left, and
            LWORK >= (N-K+NB)*(M+NB) + M*NB if SIDE = Magma_tally3Right,
            where NB is the given blocksize.
    \n
            If LWORK = -1, then a workspace query is assumed; the routine
            only calculates the optimal size of the HWORK array, returns
            this value as the first entry of the HWORK array, and no error
            message related to LWORK is issued by XERBLA.

    @param[in]
    dT      COMPLEX array on the GPU that is the output
            (the 9th argument) of magma_tally3_cgeqrf_gpu.

    @param[in]
    nb      INTEGER
            This is the blocking size that was used in pre-computing DT, e.g.,
            the blocking size used in magma_tally3_cgeqrf_gpu.

    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument had an illegal value

    @ingroup magma_tally3_cgeqrf_comp
    ********************************************************************/
extern "C" magma_tally3_int_t
magma_tally3_cunmqr_gpu(
    magma_tally3_side_t side, magma_tally3_trans_t trans,
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    magma_tally3FloatComplex_ptr dA,    magma_tally3_int_t ldda,
    magma_tally3FloatComplex    *tau,
    magma_tally3FloatComplex_ptr dC,    magma_tally3_int_t lddc,
    magma_tally3FloatComplex    *hwork, magma_tally3_int_t lwork,
    magma_tally3FloatComplex_ptr dT,    magma_tally3_int_t nb,
    magma_tally3_int_t *info)
{
    #define dA(a_1,a_2) (dA + (a_1) + (a_2)*ldda)
    #define dC(a_1,a_2) (dC + (a_1) + (a_2)*lddc)
    #define dT(a_1)     (dT + (a_1)*nb)

    magma_tally3FloatComplex c_one = MAGMA_tally3_C_ONE;

    const char* side_  = lapack_side_const_tally3( side  );
    const char* trans_ = lapack_trans_const_tally3( trans );

    magma_tally3FloatComplex_ptr dwork;
    magma_tally3_int_t i, lddwork;
    magma_tally3_int_t i1, i2, step, ib, ic, jc, ma, mi, ni, nq, nw;
    int left, notran, lquery;
    magma_tally3_int_t lwkopt;

    *info = 0;
    left   = (side == Magma_tally3Left);
    notran = (trans == Magma_tally3NoTrans);
    lquery = (lwork == -1);

    /* NQ is the order of Q and NW is the minimum dimension of WORK */
    if (left) {
        nq = m;
        nw = n;
    } else {
        nq = n;
        nw = m;
    }
    lwkopt = (nq - k + nb)*(nw + nb) + nw*nb;
    hwork[0] = MAGMA_tally3_C_MAKE( lwkopt, 0 );
    
    if ( ! left && side != Magma_tally3Right ) {
        *info = -1;
    } else if ( ! notran && trans != Magma_tally3_ConjTrans ) {
        *info = -2;
    } else if (m < 0) {
        *info = -3;
    } else if (n < 0) {
        *info = -4;
    } else if (k < 0 || k > nq) {
        *info = -5;
    } else if (ldda < max(1,nq)) {
        *info = -7;
    } else if (lddc < max(1,m)) {
        *info = -10;
    } else if (lwork < lwkopt && ! lquery) {
        *info = -12;
    }

    if (*info != 0) {
        magma_tally3_xerbla( __func__, -(*info) );
        return *info;
    }
    else if (lquery) {
        return *info;
    }

    /* Quick return if possible */
    if (m == 0 || n == 0 || k == 0) {
        hwork[0] = c_one;
        return *info;
    }

    lddwork = k;
    dwork = dT(2*lddwork);

    if ( (left && (! notran)) || ((! left) && notran) ) {
        // left  trans:    Q^T C
        // right notrans:  C Q
        // multiply from first block, i = 0, to next-to-last block, i < k-nb
        i1 = 0;
        i2 = k-nb;
        step = nb;
    } else {
        // left  notrans:  Q C
        // right trans:    C Q^T
        // multiply from next-to-last block, i = floor((k-1-nb)/nb)*nb, to first block, i = 0
        i1 = ((k - 1 - nb) / nb) * nb;
        i2 = 0;
        step = -nb;
    }

    if (left) {
        ni = n;
        jc = 0;
    } else {
        mi = m;
        ic = 0;
    }
    
    /* Use unblocked code to multiply last or only block (cases Q*C or C*Q^T). */
    // workspace left:  A(mi*nb) + C(mi*ni) + work(ni*nb_la) = (m-k-nb)*nb + (m-k-nb)*n + n*nb
    // workspace right: A(ni*nb) + C(mi*ni) + work(mi*nb_la) = (n-k-nb)*nb + m*(n-k-nb) + m*nb
    if ( step < 0 ) {
        // i is beginning of last block
        i = i1 - step;
        if ( i >= k ) {
            i = i1;
        }
        ib = k - i;
        if (left) {
            // ni=n, jc=0, H or H^T is applied to C(i:m-1,0:n-1)
            mi = m - i;
            ma = mi;
            ic = i;
        }
        else {
            // mi=m, ic=0, H or H^T is applied to C(0:m-1,i:n-1)
            ni = n - i;
            ma = ni;
            jc = i;
        }
        
        magma_tally3FloatComplex* hA = hwork;
        magma_tally3FloatComplex* hC = hwork + ma*ib;
        magma_tally3FloatComplex* hW = hwork + ma*ib + mi*ni;
        magma_tally3_int_t lhwork = lwork - (ma*ib + mi*ni);
        
        magma_tally3_cgetmatrix( ma, ib, dA(i,  i ), ldda, hA, ma );
        magma_tally3_cgetmatrix( mi, ni, dC(ic, jc), lddc, hC, mi );

        lapackf77_cunmqr( side_, trans_,
                          &mi, &ni, &ib,
                          hA, &ma, tau+i,
                          hC, &mi,
                          hW, &lhwork, info );

        // send the updated part of C back to the GPU
        magma_tally3_csetmatrix( mi, ni, hC, mi, dC(ic, jc), lddc );
    }

    /* Use blocked code to multiply blocks */
    if (nb < k) {
        for( i=i1; (step < 0 ? i >= i2 : i < i2); i += step ) {
            ib = min(nb, k - i);
            if (left) {
                // ni=n, jc=0, H or H^T is applied to C(i:m-1,0:n-1)
                mi = m - i;
                ic = i;
            }
            else {
                // mi=m, ic=0, H or H^T is applied to C(0:m-1,i:n-1)
                ni = n - i;
                jc = i;
            }
            
            magma_tally3_clarfb_gpu( side, trans, Magma_tally3Forward, Magma_tally3Columnwise,
                              mi, ni, ib,
                              dA(i,  i ), ldda, dT(i), nb,
                              dC(ic, jc), lddc, dwork, nw );
        }
    }
    else {
        i = i1;
    }

    /* Use unblocked code to multiply the last or only block (cases Q^T*C or C*Q). */
    if ( step > 0 ) {
        ib = k-i;
        if (left) {
            // ni=n, jc=0, H or H^T is applied to C(i:m-1,0:n-1)
            mi = m - i;
            ma = mi;
            ic = i;
        }
        else {
            // mi=m, ic=0, H or H^T is applied to C(0:m-1,i:n-1)
            ni = n - i;
            ma = ni;
            jc = i;
        }
        
        magma_tally3FloatComplex* hA = hwork;
        magma_tally3FloatComplex* hC = hwork + ma*ib;
        magma_tally3FloatComplex* hW = hwork + ma*ib + mi*ni;
        magma_tally3_int_t lhwork = lwork - (ma*ib + mi*ni);
        
        magma_tally3_cgetmatrix( ma, ib, dA(i,  i ), ldda, hA, ma );
        magma_tally3_cgetmatrix( mi, ni, dC(ic, jc), lddc, hC, mi );

        lapackf77_cunmqr( side_, trans_,
                          &mi, &ni, &ib,
                          hA, &ma, tau+i,
                          hC, &mi,
                          hW, &lhwork, info );
        
        // send the updated part of C back to the GPU
        magma_tally3_csetmatrix( mi, ni, hC, mi, dC(ic, jc), lddc );
    }
    
    // TODO sync. For cases Q*C and C*Q^T, last call is magma_tally3_clarfb_gpu,
    // which is async magma_tally3_gemm calls, so cunmqr can be unfinished.

    // TODO: cgeqrs_gpu ASSUMES that hwork contains the last block of A and C.
    // That needs to be fixed, but until then, don't modify hwork[0] here.
    // In LAPACK: On exit, if INFO = 0, HWORK[0] returns the optimal LWORK.
    //hwork[0] = MAGMA_tally3_C_MAKE( lwkopt, 0 );
    return *info;
} /* magma_tally3_cunmqr_gpu */
