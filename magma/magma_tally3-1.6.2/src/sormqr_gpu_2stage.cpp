/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Azzam Haidar
       @author Stan Tomov
       @author Raffaele Solca

       @generated from zunmqr_gpu_2stage.cpp normal z -> s, Fri Jan 30 19:00:17 2015

*/
#include "common_magma_tally3.h"

/**
    Purpose
    -------
    SORMQR_GPU overwrites the general real M-by-N matrix C with

    @verbatim
                               SIDE = Magma_tally3Left    SIDE = Magma_tally3Right
    TRANS = Magma_tally3NoTrans:      Q * C               C * Q
    TRANS = Magma_tally3Trans:   Q**H * C            C * Q**H
    @endverbatim

    where Q is a real unitary matrix defined as the product of k
    elementary reflectors

        Q = H(1) H(2) . . . H(k)

    as returned by SGEQRF. Q is of order M if SIDE = Magma_tally3Left and of order N
    if SIDE = Magma_tally3Right.

    Arguments
    ---------
    @param[in]
    side    magma_tally3_side_t
      -      = Magma_tally3Left:      apply Q or Q**H from the Left;
      -      = Magma_tally3Right:     apply Q or Q**H from the Right.

    @param[in]
    trans   magma_tally3_trans_t
      -     = Magma_tally3NoTrans:    No transpose, apply Q;
      -     = Magma_tally3Trans: Conjugate transpose, apply Q**H.

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
    dA      REAL array on the GPU, dimension (LDDA,K)
            The i-th column must contain the vector which defines the
            elementary reflector H(i), for i = 1,2,...,k, as returned by
            SGEQRF in the first k columns of its array argument DA.
            DA is modified by the routine but restored on exit.

    @param[in]
    ldda    INTEGER
            The leading dimension of the array DA.
            If SIDE = Magma_tally3Left,  LDDA >= max(1,M);
            if SIDE = Magma_tally3Right, LDDA >= max(1,N).

    @param[in,out]
    dC      REAL array on the GPU, dimension (LDDC,N)
            On entry, the M-by-N matrix C.
            On exit, C is overwritten by Q*C or Q**H * C or C * Q**H or C*Q.

    @param[in]
    lddc     INTEGER
            The leading dimension of the array DC. LDDC >= max(1,M).

    @param[in]
    dT      REAL array on the GPU that is the output
            (the 9th argument) of magma_tally3_sgeqrf_gpu.

    @param[in]
    nb      INTEGER
            This is the blocking size that was used in pre-computing DT, e.g.,
            the blocking size used in magma_tally3_sgeqrf_gpu.

    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument had an illegal value

    @ingroup magma_tally3_ssyev_2stage
    ********************************************************************/
extern "C" magma_tally3_int_t
magma_tally3_sormqr_gpu_2stages(
    magma_tally3_side_t side, magma_tally3_trans_t trans,
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    magma_tally3Float_ptr dA,   magma_tally3_int_t ldda,
    magma_tally3Float_ptr dC,    magma_tally3_int_t lddc,
    magma_tally3Float_ptr dT,    magma_tally3_int_t nb,
    magma_tally3_int_t *info)
{
    #define dA(i_,j_) (dA + (i_) + (j_)*ldda)
    #define dC(i_,j_) (dC + (i_) + (j_)*lddc)
    
    magma_tally3Float_ptr dwork;

    magma_tally3_int_t i, i1, i2, step, ib, ic, jc, mi, ni, nq, nw;
    int left, notran;
    //magma_tally3_int_t lwkopt;

    *info = 0;
    left   = (side == Magma_tally3Left);
    notran = (trans == Magma_tally3NoTrans);

    /* NQ is the order of Q and NW is the minimum dimension of WORK */
    if (left) {
        nq = m;
        nw = n;
    } else {
        nq = n;
        nw = m;
    }
    if ( ! left && side != Magma_tally3Right ) {
        *info = -1;
    } else if ( ! notran && trans != Magma_tally3Trans ) {
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
    }

    // TODO alloc after xerbla & quick return, else memory leak
    if (MAGMA_tally3_SUCCESS != magma_tally3_smalloc( &dwork, n*nb )) {
        printf ("!!!! sorgqr_2stage magma_tally3_alloc failed for: dwork\n" );
        exit(-1);
    }

    if (*info != 0) {
        magma_tally3_xerbla( __func__, -(*info) );
        return *info;
    }

    /* Quick return if possible */
    if (m == 0 || n == 0 || k == 0) {
        return *info;
    }

    if ( (left && (! notran)) || ( (! left) && notran ) ) {
        i1 = 0;
        i2 = k;
        step = nb;
    } else {
        i1 = (k - 1) / nb * nb;
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

    for (i=i1; (step < 0 ? i >= i2 : i < i2); i += step) {
        ib = min(nb, k - i);
        if (left) {
            mi = m - i;
            ic = i;
        }
        else {
            ni = n - i;
            jc = i;
        }
        magma_tally3_slarfb_gpu( Magma_tally3Left, trans, Magma_tally3Forward, Magma_tally3Columnwise,
                          mi, ni, ib, dA(i,i), ldda, dT+i*nb, nb,
                          dC(ic,jc), lddc, dwork, nw );
    }
    
    magma_tally3_free( dwork );
    return *info;
} /* magma_tally3_sormqr_gpu_2stages */
