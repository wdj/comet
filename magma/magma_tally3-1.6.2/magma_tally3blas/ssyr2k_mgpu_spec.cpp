/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from zher2k_mgpu_spec.cpp normal z -> s, Fri Jan 30 19:00:10 2015
       @author Mark Gates
       @author Azzam Haidar 
*/
#include "common_magma_tally3.h"

/**
    Purpose
    -------
    SSYR2K performs one of the symmetric rank 2k operations

       C := alpha*A*B**H + conjg( alpha )*B*A**H + beta*C,

    or

       C := alpha*A**H*B + conjg( alpha )*B**H*A + beta*C,

    where alpha and beta are scalars with beta real, C is an n by n
    symmetric matrix and A and B are n by k matrices in the first case
    and k by n matrices in the second case.
    
    This version assumes C has been symmetrized, so both upper and lower are
    stored, and it maintains the symmetry, doing twice the operations.

    Arguments
    ----------

    @param[in]
    uplo     magma_tally3_uplo_t.
             On entry, UPLO specifies whether the upper or lower
             triangular part of the array C is to be referenced as
             follows:
      -     = Magma_tally3Upper:  Only the upper triangular part of C is to be referenced.
      -     = Magma_tally3Lower:  Only the lower triangular part of C is to be referenced.

             **** current only Lower case is implemented.

    @param[in]
    trans    magma_tally3_trans_t.
             On entry, TRANS specifies the operation to be performed as
             follows:
      -     = Magma_tally3NoTrans:     C := alpha*A*B**H + conj( alpha )*B*A**H + beta*C.
      -     = Magma_tally3Trans:  C := alpha*A**H*B + conj( alpha )*B**H*A + beta*C.

             **** current only NoTrans case is implemented.

    @param[in]
    n        INTEGER.
             On entry, N specifies the order of the matrix C. N must be
             at least zero.

    @param[in]
    k        INTEGER.
             On entry with TRANS = Magma_tally3NoTrans, K specifies the number
             of columns of the matrices A and B, and on entry with
             TRANS = Magma_tally3Trans, K specifies the number of rows of the
             matrices A and B. K must be at least zero.

    @param[in]
    alpha    REAL.
             On entry, ALPHA specifies the scalar alpha.

    @param[in]
    dA       REAL array of DIMENSION ( LDA, ka ), where ka is
             k when TRANS = Magma_tally3NoTrans, and is n otherwise.
             Before entry with TRANS = Magma_tally3NoTrans, the leading n by k
             part of the array A must contain the matrix A, otherwise
             the leading k by n part of the array A must contain the
             matrix A.
             
             [TODO: describe distribution: duplicated on all GPUs.]

    @param[in]
    ldda     INTEGER.
             On entry, LDA specifies the first dimension of A as declared
             in the calling (sub) program. When TRANS = Magma_tally3NoTrans
             then LDA must be at least max( 1, n ), otherwise LDA must
             be at least max( 1, k ).

    @param[in]
    a_offset INTEGER
             Row offset to start sub-matrix of dA. Uses dA(a_offset:a_offset+n, :).
             0 <= a_offset < ldda.

    @param[in]
    dB       REAL array of DIMENSION ( LDB, kb ), where kb is
             k when TRANS = Magma_tally3NoTrans, and is n otherwise.
             Before entry with TRANS = Magma_tally3NoTrans, the leading n by k
             part of the array B must contain the matrix B, otherwise
             the leading k by n part of the array B must contain the
             matrix B.
             
             [TODO: describe distribution: duplicated on all GPUs.]

    @param[in]
    lddb     INTEGER.
             On entry, LDB specifies the first dimension of B as declared
             in the calling (sub) program. When TRANS = Magma_tally3NoTrans
             then LDB must be at least max( 1, n ), otherwise LDB must
             be at least max( 1, k ).

    @param[in]
    b_offset INTEGER
             Row offset to start sub-matrix of dB. Uses dB(b_offset:b_offset+n, :).
             0 <= b_offset < lddb.

    @param[in]
    beta     REAL.
             On entry, BETA specifies the scalar beta.

    @param[in,out]
    dC       REAL array of DIMENSION ( LDC, n ).
             Before entry with UPLO = Magma_tally3Upper, the leading n by n
             upper triangular part of the array C must contain the upper
             triangular part of the symmetric matrix and the strictly
             lower triangular part of C is not referenced. On exit, the
             upper triangular part of the array C is overwritten by the
             upper triangular part of the updated matrix.
    \n
             Before entry with UPLO = Magma_tally3Lower, the leading n by n
             lower triangular part of the array C must contain the lower
             triangular part of the symmetric matrix and the strictly
             upper triangular part of C is not referenced. On exit, the
             lower triangular part of the array C is overwritten by the
             lower triangular part of the updated matrix.
    \n
             Note that the imaginary parts of the diagonal elements need
             not be set, they are assumed to be zero, and on exit they
             are set to zero. [TODO: verify]
             
             [TODO: describe distribution: 1D column block-cyclic across GPUs.]

    @param[in]
    lddc     INTEGER.
             On entry, LDC specifies the first dimension of C as declared
             in the calling (sub) program. LDC must be at least max( 1, n ).

    @param[in]
    c_offset INTEGER.
             Row and column offset to start sub-matrix of dC.
             Uses dC(c_offset:c_offset+n, c_offset:c_offset+n).
             0 <= c_offset < lddc.

    @param[in]
    ngpu     INTEGER.
             Number of GPUs over which matrix C is distributed.

    @param[in]
    nb       INTEGER.
             Block size used for distribution of C.

    @param[in]
    queues   array of CUDA queues, of dimension NGPU by 20.
             Streams to use for running multiple GEMMs in parallel.
             Only up to NSTREAM queues are used on each GPU.

    @param[in]
    nqueue   INTEGER.
             Number of queues to use on each device

    @ingroup magma_tally3_sblas3
    ********************************************************************/
extern "C"
void magma_tally3blas_ssyr2k_mgpu_spec(
    magma_tally3_uplo_t uplo, magma_tally3_trans_t trans, magma_tally3_int_t n, magma_tally3_int_t k,
    float alpha,
    magma_tally3Float_ptr dA[], magma_tally3_int_t ldda, magma_tally3_int_t a_offset,
    magma_tally3Float_ptr dB[], magma_tally3_int_t lddb, magma_tally3_int_t b_offset,
    float beta,
    magma_tally3Float_ptr dC[], magma_tally3_int_t lddc, magma_tally3_int_t c_offset,
    magma_tally3_int_t ngpu, magma_tally3_int_t nb, magma_tally3_queue_t queues[][20], magma_tally3_int_t nqueue )
{
    #define dA(dev, i, j) (dA[dev] + (i) + (j)*ldda + (a_offset) )
    #define dB(dev, i, j) (dB[dev] + (i) + (j)*lddb + (b_offset) )
    #define dC(dev, i, j) (dC[dev] + (i) + (j)*lddc)
    
    /* Check arguments */
    magma_tally3_int_t info = 0;
    if ( uplo != Magma_tally3Lower ) {
        info = -1;  // 'u' not yet handled
    } else if ( trans != Magma_tally3NoTrans ) {
        info = -2;  // 'c' not yet handled
    } else if ( n < 0 ) {
        info = -3;
    } else if ( k < 0 ) {
        info = -4;
    } else if ( ((trans == Magma_tally3NoTrans)    && ldda < max(1,n)) ||
                ((trans == Magma_tally3Trans) && ldda < max(1,k)) ) {
        info = -7;
    } else if ( a_offset < 0 || a_offset > ldda ) {
        info = -8;
    } else if ( ((trans == Magma_tally3NoTrans)    && lddb < max(1,n)) ||
                ((trans == Magma_tally3Trans) && lddb < max(1,k)) ) {
        info = -10;
    } else if ( b_offset < 0 || b_offset > lddb ) {
        info = -11;
    } else if ( lddc < max(1,n) ) {
        info = -13;
    } else if ( c_offset < 0 || c_offset > lddc ) {
        info = -14;
    } else if ( ngpu <= 0 ) {
        info = -15;
    } else if ( nb <= 0 ) {
        info = -16;
    } else if ( nqueue <= 0 ) {
        info = -18;
    }
    if ( info != 0 ) {
        magma_tally3_xerbla( __func__, -(info) );
        return;
    }
    
    const float c_one = MAGMA_tally3_S_ONE;
    float cbeta = MAGMA_tally3_S_MAKE( beta, 0. );
    
    magma_tally3_int_t ib, ioff, iblock, idev, di, s;
    
    magma_tally3_device_t cdev;
    magma_tally3_queue_t cqueue;
    magma_tally3_getdevice( &cdev );
    magma_tally3blasGetKernelStream( &cqueue );
    
    // loop over all blocks
    // Faster to have two loops: first loop does C_hat = alpha*A*B**H + beta*C
    // blockoffset is offset within first block; for subsequent blocks it is 0
    magma_tally3_int_t blockoffset = c_offset % nb;
    for( magma_tally3_int_t i = 0; i < n; i += ib ) {
        ib     = min( nb-blockoffset, n-i );  // block size
        ioff   = i + c_offset;                 // global index in parent matrix
        iblock = (ioff / nb) / ngpu;          // local block id
        idev   = (ioff / nb) % ngpu;          // device with this block
        di     = iblock*nb + blockoffset;     // local index in parent matrix
        
        magma_tally3_setdevice( idev );
        s = iblock % nqueue;
        magma_tally3blasSetKernelStream( queues[ idev ][ s ] );
        
        // C[i:n,i] = alpha * A[i:n,0] * B[i,0]' + beta*C[i:n,i]
        //printf( "sgemm  n=%4d, ib=%4d, k=%4d, i=%4d\n", n-i, ib, k, i );
        magma_tally3_sgemm( Magma_tally3NoTrans, Magma_tally3Trans, n, ib, k,
                     alpha, dA(idev,0,0), ldda,
                            dB(idev,i,0), lddb,
                     cbeta, dC(idev,c_offset,di), lddc );
        blockoffset = 0;
    }
    
    // second loop does C = conj(alpha)*B*A**H + C_hat
    alpha = MAGMA_tally3_S_CNJG( alpha );
    blockoffset = c_offset % nb;
    for( magma_tally3_int_t i = 0; i < n; i += ib ) {
        ib     = min( nb-blockoffset, n-i );  // block size
        ioff   = i + c_offset;                 // global index in parent matrix
        iblock = (ioff / nb) / ngpu;          // local block id
        idev   = (ioff / nb) % ngpu;          // device with this block
        di     = iblock*nb + blockoffset;     // local index in parent matrix
        
        magma_tally3_setdevice( idev );
        s = iblock % nqueue;
        magma_tally3blasSetKernelStream( queues[ idev ][ s ] );
        
        // C[i:n,i] += conj(alpha) * B[i:n,0] * A[i,0]'
        //printf( "sgemm  n=%4d, ib=%4d, k=%4d, i=%4d\n", n-i, ib, k, i );
        magma_tally3_sgemm( Magma_tally3NoTrans, Magma_tally3Trans, n, ib, k,
                     alpha, dB(idev,0,0), lddb,
                            dA(idev,i,0), ldda,
                     c_one, dC(idev,c_offset,di), lddc );
        blockoffset = 0;
    }
    
    magma_tally3_setdevice( cdev );
    magma_tally3blasSetKernelStream( cqueue );
}
