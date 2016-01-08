/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from zgeqr2x_gpu-v3.cpp normal z -> c, Fri Jan 30 19:00:15 2015

*/
#include "common_magma_tally4.h"

extern "C" magma_tally4_int_t
magma_tally4_clarfb2_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4FloatComplex_const_ptr dV,    magma_tally4_int_t lddv,
    magma_tally4FloatComplex_const_ptr dT,    magma_tally4_int_t lddt,
    magma_tally4FloatComplex_ptr       dC,    magma_tally4_int_t lddc,
    magma_tally4FloatComplex_ptr       dwork, magma_tally4_int_t ldwork )
{
    magma_tally4FloatComplex c_zero    = MAGMA_tally4_C_ZERO;
    magma_tally4FloatComplex c_one     = MAGMA_tally4_C_ONE;
    magma_tally4FloatComplex c_neg_one = MAGMA_tally4_C_NEG_ONE;

    if (m <= 0 || n <= 0)
        return MAGMA_tally4_SUCCESS;

    // W = C^H V
    magma_tally4_cgemm( Magma_tally4ConjTrans, Magma_tally4NoTrans,
    //magma_tally4blas_cgemm_reduce(
                           n, k, m,
                           c_one,  dC,    lddc,
                                   dV,    lddv,
                           c_zero, dwork, ldwork);

    // W = W T^H = C^H V T^H
    magma_tally4_ctrmm( Magma_tally4Right, Magma_tally4Upper, Magma_tally4NoTrans, Magma_tally4NonUnit,
                 n, k,
                 c_one, dT,    lddt,
                        dwork, ldwork);

    // C = C - V W^H = C - V T V^H C = (I - V T V^H) C = H C
    magma_tally4_cgemm( Magma_tally4NoTrans, Magma_tally4ConjTrans,
                 m, n, k,
                 c_neg_one, dV,    lddv,
                            dwork, ldwork,
                 c_one,     dC,    lddc);
    
    return MAGMA_tally4_SUCCESS;
}

//////////////////////////////////////////////////////////////////////////////

/**
    Purpose
    -------
    CGEQR2 computes a QR factorization of a complex m by n matrix A:
    A = Q * R.

    This expert routine requires two more arguments than the standard
    cgeqr2, namely, dT and ddA, explained below. The storage for A is
    also not as in the LAPACK's cgeqr2 routine (see below).

    The first is used to output the triangular
    n x n factor T of the block reflector used in the factorization.
    The second holds the diagonal nxn blocks of A, i.e., the diagonal
    submatrices of R. This routine implements the left looking QR.

    This version adds internal blocking.

    Arguments
    ---------
    @param[in]
    m       INTEGER
            The number of rows of the matrix A.  M >= 0.

    @param[in]
    n       INTEGER
            The number of columns of the matrix A.  N >= 0.

    @param[in,out]
    dA      COMPLEX array, dimension (LDA,N)
            On entry, the m by n matrix A.
            On exit, the unitary matrix Q as a
            product of elementary reflectors (see Further Details).
    \n
            the elements on and above the diagonal of the array
            contain the min(m,n) by n upper trapezoidal matrix R (R is
            upper triangular if m >= n); the elements below the diagonal,
            with the array TAU, represent the unitary matrix Q as a
            product of elementary reflectors (see Further Details).

    @param[in]
    ldda    INTEGER
            The leading dimension of the array A.  LDA >= max(1,M).

    @param[out]
    dtau    COMPLEX array, dimension (min(M,N))
            The scalar factors of the elementary reflectors (see Further
            Details).

    @param[out]
    dT      COMPLEX array, dimension N x N.
            Stores the triangular N x N factor T of the block reflector
            used in the factorization. The lower triangular part is 0.

    @param[out]
    ddA     COMPLEX array, dimension N x N.
            Stores the elements of the upper N x N diagonal block of A.
            LAPACK stores this array in A. There are 0s below the diagonal.

    @param
    dwork   (workspace) DOUBLE_PRECISION array, dimension (3 N)

    @param[out]
    info    INTEGER
      -     = 0: successful exit
      -     < 0: if INFO = -i, the i-th argument had an illegal value

    Further Details
    ---------------
    The matrix Q is represented as a product of elementary reflectors

       Q = H(1) H(2) . . . H(k), where k = min(m,n).

    Each H(i) has the form

       H(i) = I - tau * v * v'

    where tau is a complex scalar, and v is a complex vector with
    v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
    and tau in TAU(i).

    @ingroup magma_tally4_cgeqrf_comp
    ********************************************************************/
extern "C" magma_tally4_int_t
magma_tally4_cgeqr2x3_gpu(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4FloatComplex_ptr dA, magma_tally4_int_t ldda,
    magma_tally4FloatComplex_ptr dtau,
    magma_tally4FloatComplex_ptr dT,
    magma_tally4FloatComplex_ptr ddA,
    magma_tally4Float_ptr        dwork,
    magma_tally4_int_t *info)
{
    #define dA(i_,j_) (dA + (j_)*(ldda) + (i_))
    #define BLOCK_SIZE 32

    magma_tally4_int_t i, k;

    magma_tally4Float_ptr dnorm = dwork;
    magma_tally4FloatComplex *work = (magma_tally4FloatComplex *)(dwork+2*(n));

    *info = 0;
    if (m < 0) {
        *info = -1;
    } else if (n < 0) {
        *info = -2;
    } else if (ldda < max(1,m)) {
        *info = -4;
    }
    if (*info != 0) {
        magma_tally4_xerbla( __func__, -(*info) );
        return *info;
    }

    /* Compute the norms of the trailing columns */
    k = min(m,n);
    // magma_tally4blas_scnrm2_cols(m, k, dA(0,0), ldda, dnorm);

    for (int b=0; b < k; b += BLOCK_SIZE) {
        for (i = b; i < min(k, b+BLOCK_SIZE); ++i) {
            /*   Apply H' to A(:,i) from the left */
            if ( i-b > 0)
                magma_tally4_clarfbx_gpu(m-b, i-b, dA(b, b), ldda,
                                  dT+b+b*k, k, dA(b, i), work);

            /*   Adjust the dnorm[i] to hold the norm of A(i:m,i) */
            //if ( i > 0 )
            //    magma_tally4blas_scnrm2_adjust(i, dnorm+i, dA(0, i));
            magma_tally4blas_scnrm2_cols(m-i, 1, dA(i,i), ldda, dnorm+i);
            
            /*  Generate elementary reflector H(i) to annihilate A(i+1:m,i)
                1. 1 is not yet put on the diagonal of A
                2. Elements above the diagonal are copied in ddA and
                   the ones in A are set to zero
                3. update T */
            magma_tally4_clarfgtx_gpu(m-i, dA(i, i), dA(min(i+1,m), i), dtau+i,
                               dnorm+i, ddA + i + i*(n), i,
                               dA(i,0), ldda,  dT, k, work);
        }
        
        /* Apply the transformations to the trailing matrix. */
        //magma_tally4_clarfb2_gpu( Magma_tally4Left, Magma_tally4ConjTrans, Magma_tally4Forward, Magma_tally4Columnwise,
        magma_tally4_clarfb2_gpu(
                           m-b, k-i, BLOCK_SIZE,
                           dA(b, b), ldda, dT+b+b*k, k,
                           dA(b, i), ldda, work, k-i);
    }

    return *info;
} /* magma_tally4_cgeqr2 */
