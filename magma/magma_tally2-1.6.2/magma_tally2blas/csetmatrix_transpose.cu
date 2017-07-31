/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from zsetmatrix_transpose.cu normal z -> c, Fri Jan 30 19:00:09 2015

*/
#include "common_magma_tally2.h"

#define PRECISION_c


//
//      m, n - dimensions in the source (input) matrix.
//             This routine copies the hA matrix from the CPU
//             to dAT on the GPU. In addition, the output matrix
//             is transposed. The routine uses a buffer of size
//             2*lddwork*nb pointed to by dwork (lddwork > m) on the GPU. 
//             Note that lda >= m and lddat >= n.
//
extern "C" void 
magma_tally2blas_csetmatrix_transpose_q(
    magma_tally2_int_t m, magma_tally2_int_t n,
    const magma_tally2FloatComplex     *hA, magma_tally2_int_t lda, 
    magma_tally2FloatComplex_ptr       dAT, magma_tally2_int_t ldda,
    magma_tally2FloatComplex_ptr     dwork, magma_tally2_int_t lddwork, magma_tally2_int_t nb,
    magma_tally2_queue_t queues[2] )
{
#define    hA(i_, j_)    (hA + (i_) + (j_)*lda)
#define   dAT(i_, j_)   (dAT + (i_) + (j_)*ldda)
#define dwork(i_, j_) (dwork + (i_) + (j_)*lddwork)

    magma_tally2_int_t i = 0, j = 0, ib;

    /* Quick return */
    if ( (m == 0) || (n == 0) )
        return;

    // TODO standard check arguments
    if (lda < m || ldda < n || lddwork < m){
        printf("Wrong arguments in %s.\n", __func__);
        return;
    }

    /* Move data from CPU to GPU in the first panel in the dwork buffer */
    ib = min(n-i, nb);
    magma_tally2_csetmatrix_async( m, ib,
                            hA(0,i), lda,
                            dwork(0,(j%2)*nb), lddwork, queues[j%2] );
    j++;

    for(i=nb; i < n; i += nb) {
        /* Move data from CPU to GPU in the second panel in the dwork buffer */
        ib = min(n-i, nb);
        magma_tally2_csetmatrix_async( m, ib,
                                hA(0,i), lda,
                                dwork(0,(j%2)*nb), lddwork, queues[j%2] );
        j++;
        
        /* Note that the previous panel (i.e., j%2) comes through the queue
           for the kernel so there is no need to synchronize.             */
        // TODO should this be ib not nb?
        magma_tally2blas_ctranspose_q( m, nb, dwork(0,(j%2)*nb), lddwork, dAT(i-nb,0), ldda, queues[j%2] );
    }

    /* Transpose the last part of the matrix.                            */
    j++;
    magma_tally2blas_ctranspose_q( m, ib, dwork(0,(j%2)*nb), lddwork, dAT(i-nb,0), ldda, queues[j%2] );
}


// @see magma_tally2blas_csetmatrix_transpose_q
extern "C" void 
magma_tally2blas_csetmatrix_transpose(
    magma_tally2_int_t m, magma_tally2_int_t n,
    const magma_tally2FloatComplex     *hA, magma_tally2_int_t lda, 
    magma_tally2FloatComplex_ptr       dAT, magma_tally2_int_t ldda,
    magma_tally2FloatComplex_ptr     dwork, magma_tally2_int_t lddwork, magma_tally2_int_t nb )
{
    magma_tally2_queue_t queues[2];
    magma_tally2_queue_create( &queues[0] );
    magma_tally2_queue_create( &queues[1] );

    magma_tally2blas_csetmatrix_transpose_q( m, n, hA, lda, dAT, ldda, dwork, lddwork, nb, queues );
    
    magma_tally2_queue_destroy( queues[0] );
    magma_tally2_queue_destroy( queues[1] );
}
