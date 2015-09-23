/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from zswapdblk.cu normal z -> c, Fri Jan 30 19:00:09 2015

*/
#include "common_magma_minproduct.h"


/*********************************************************/
/*
 *  Swap diagonal blocks of two matrices.
 *  Each thread block swaps one diagonal block.
 *  Each thread iterates across one row of the block.
 */

__global__ void 
cswapdblk_kernel( int nb,
                  magma_minproductFloatComplex *dA, int ldda, int inca,
                  magma_minproductFloatComplex *dB, int lddb, int incb )
{
    const int tx = threadIdx.x;
    const int bx = blockIdx.x;

    dA += tx + bx * nb * (ldda + inca);
    dB += tx + bx * nb * (lddb + incb);

    magma_minproductFloatComplex tmp;

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
    cswapdblk swaps diagonal blocks of size nb x nb between matrices
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
    dA      COMPLEX array, dimension (LDDA,N)
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
    dB      COMPLEX array, dimension (LDDB,N)
            The matrix dB.

    @param[in]
    lddb    INTEGER
            The leading dimension of the array db.
            LDDB >= (nblocks - 1)*nb*incb + nb.

    @param[in]
    incb    INTEGER
            The row increment between diagonal blocks of dB. incb >= 0. See inca.
    
    @param[in]
    queue   magma_minproduct_queue_t
            Queue to execute in.

    @ingroup magma_minproduct_caux2
    ********************************************************************/
extern "C" void 
magma_minproductblas_cswapdblk_q(
    magma_minproduct_int_t n, magma_minproduct_int_t nb,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda, magma_minproduct_int_t inca,
    magma_minproductFloatComplex_ptr dB, magma_minproduct_int_t lddb, magma_minproduct_int_t incb,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t nblocks = n / nb;
    
    magma_minproduct_int_t info = 0;
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
        magma_minproduct_xerbla( __func__, -(info) );
        return;  //info;
    }

    if ( nblocks > 0 ) {
        cswapdblk_kernel<<< nblocks, nb, 0, queue >>>
            ( nb, dA, ldda, inca,
                  dB, lddb, incb );
    }
}


/**
    @see magma_minproductblas_cswapdblk_q
    @ingroup magma_minproduct_caux2
    ********************************************************************/
extern "C" void 
magma_minproductblas_cswapdblk(
    magma_minproduct_int_t n, magma_minproduct_int_t nb,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda, magma_minproduct_int_t inca,
    magma_minproductFloatComplex_ptr dB, magma_minproduct_int_t lddb, magma_minproduct_int_t incb )
{
    magma_minproductblas_cswapdblk_q( n, nb, dA, ldda, inca, dB, lddb, incb, magma_minproduct_stream );
}
