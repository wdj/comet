/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions normal z -> s d c

*/
#include "common_magma_tally3.h"


/*********************************************************/
/*
 *  Swap diagonal blocks of two matrices.
 *  Each thread block swaps one diagonal block.
 *  Each thread iterates across one row of the block.
 */

__global__ void 
zswapdblk_kernel( int nb,
                  magma_tally3DoubleComplex *dA, int ldda, int inca,
                  magma_tally3DoubleComplex *dB, int lddb, int incb )
{
    const int tx = threadIdx.x;
    const int bx = blockIdx.x;

    dA += tx + bx * nb * (ldda + inca);
    dB += tx + bx * nb * (lddb + incb);

    magma_tally3DoubleComplex tmp;

    #pragma unroll
    for( int i = 0; i < nb; i++ ){
        tmp        = dA[i*ldda];
        dA[i*ldda] = dB[i*lddb];
        dB[i*lddb] = tmp;
    }
}


/**
    Purpose
    -------
    zswapdblk swaps diagonal blocks of size nb x nb between matrices
    dA and dB on the GPU. It swaps nblocks = n/nb blocks.
    For i = 1 .. nblocks, submatrices
    dA( i*nb*inca, i*nb ) and
    dB( i*nb*incb, i*nb ) are swapped.
    
    Arguments
    ---------
    @param[in]
    n       INTEGER
            The number of columns of the matrices dA and dB.  N >= 0.

    @param[in]
    nb      INTEGER
            The size of diagonal blocks.
            NB > 0 and NB <= maximum threads per CUDA block (512 or 1024).

    @param[in,out]
    dA      COMPLEX_16 array, dimension (LDDA,N)
            The matrix dA.

    @param[in]
    ldda    INTEGER
            The leading dimension of the array dA.
            LDDA >= (nblocks - 1)*nb*inca + nb.

    @param[in]
    inca    INTEGER
            The row increment between diagonal blocks of dA. inca >= 0. For example,
            inca = 1 means blocks are stored on the diagonal at dA(i*nb, i*nb),
            inca = 0 means blocks are stored side-by-side    at dA(0,    i*nb).

    @param[in,out]
    dB      COMPLEX_16 array, dimension (LDDB,N)
            The matrix dB.

    @param[in]
    lddb    INTEGER
            The leading dimension of the array db.
            LDDB >= (nblocks - 1)*nb*incb + nb.

    @param[in]
    incb    INTEGER
            The row increment between diagonal blocks of dB. incb >= 0. See inca.
    
    @param[in]
    queue   magma_tally3_queue_t
            Queue to execute in.

    @ingroup magma_tally3_zaux2
    ********************************************************************/
extern "C" void 
magma_tally3blas_zswapdblk_q(
    magma_tally3_int_t n, magma_tally3_int_t nb,
    magma_tally3DoubleComplex_ptr dA, magma_tally3_int_t ldda, magma_tally3_int_t inca,
    magma_tally3DoubleComplex_ptr dB, magma_tally3_int_t lddb, magma_tally3_int_t incb,
    magma_tally3_queue_t queue )
{
    magma_tally3_int_t nblocks = n / nb;
    
    magma_tally3_int_t info = 0;
    if (n < 0) {
        info = -1;
    } else if (nb < 1 || nb > 1024) {
        info = -2;
    } else if (ldda < (nblocks-1)*nb*inca + nb) {
        info = -4;
    } else if (inca < 0) {
        info = -5;
    } else if (lddb < (nblocks-1)*nb*incb + nb) {
        info = -7;
    } else if (incb < 0) {
        info = -8;
    }

    if (info != 0) {
        magma_tally3_xerbla( __func__, -(info) );
        return;  //info;
    }

    if ( nblocks > 0 ) {
        zswapdblk_kernel<<< nblocks, nb, 0, queue >>>
            ( nb, dA, ldda, inca,
                  dB, lddb, incb );
    }
}


/**
    @see magma_tally3blas_zswapdblk_q
    @ingroup magma_tally3_zaux2
    ********************************************************************/
extern "C" void 
magma_tally3blas_zswapdblk(
    magma_tally3_int_t n, magma_tally3_int_t nb,
    magma_tally3DoubleComplex_ptr dA, magma_tally3_int_t ldda, magma_tally3_int_t inca,
    magma_tally3DoubleComplex_ptr dB, magma_tally3_int_t lddb, magma_tally3_int_t incb )
{
    magma_tally3blas_zswapdblk_q( n, nb, dA, ldda, inca, dB, lddb, incb, magma_tally3_stream );
}
