/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from zsetmatrix_transpose.cu normal z -> s, Fri Jan 30 19:00:09 2015

*/
#include "common_magma_tally4.h"

#define PRECISION_s


//
//      m, n - dimensions in the source (input) matrix.
//             This routine copies the hA matrix from the CPU
//             to dAT on the GPU. In addition, the output matrix
//             is transposed. The routine uses a buffer of size
//             2*lddwork*nb pointed to by dwork (lddwork > m) on the GPU. 
//             Note that lda >= m and lddat >= n.
//
extern "C" void 
magma_tally4blas_ssetmatrix_transpose_q(
    magma_tally4_int_t m, magma_tally4_int_t n,
    const float     *hA, magma_tally4_int_t lda, 
    magma_tally4Float_ptr       dAT, magma_tally4_int_t ldda,
    magma_tally4Float_ptr     dwork, magma_tally4_int_t lddwork, magma_tally4_int_t nb,
    magma_tally4_queue_t queues[2] )
{
#define    hA(i_, j_)    (hA + (i_) + (j_)*lda)
#define   dAT(i_, j_)   (dAT + (i_) + (j_)*ldda)
#define dwork(i_, j_) (dwork + (i_) + (j_)*lddwork)

    magma_tally4_int_t i = 0, j = 0, ib;

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
    magma_tally4_ssetmatrix_async( m, ib,
                            hA(0,i), lda,
                            dwork(0,(j%2)*nb), lddwork, queues[j%2] );
    j++;

    for(i=nb; i < n; i += nb) {
        /* Move data from CPU to GPU in the second panel in the dwork buffer */
        ib = min(n-i, nb);
        magma_tally4_ssetmatrix_async( m, ib,
                                hA(0,i), lda,
                                dwork(0,(j%2)*nb), lddwork, queues[j%2] );
        j++;
        
        /* Note that the previous panel (i.e., j%2) comes through the queue
           for the kernel so there is no need to synchronize.             */
        // TODO should this be ib not nb?
        magma_tally4blas_stranspose_q( m, nb, dwork(0,(j%2)*nb), lddwork, dAT(i-nb,0), ldda, queues[j%2] );
    }

    /* Transpose the last part of the matrix.                            */
    j++;
    magma_tally4blas_stranspose_q( m, ib, dwork(0,(j%2)*nb), lddwork, dAT(i-nb,0), ldda, queues[j%2] );
}


// @see magma_tally4blas_ssetmatrix_transpose_q
extern "C" void 
magma_tally4blas_ssetmatrix_transpose(
    magma_tally4_int_t m, magma_tally4_int_t n,
    const float     *hA, magma_tally4_int_t lda, 
    magma_tally4Float_ptr       dAT, magma_tally4_int_t ldda,
    magma_tally4Float_ptr     dwork, magma_tally4_int_t lddwork, magma_tally4_int_t nb )
{
    magma_tally4_queue_t queues[2];
    magma_tally4_queue_create( &queues[0] );
    magma_tally4_queue_create( &queues[1] );

    magma_tally4blas_ssetmatrix_transpose_q( m, n, hA, lda, dAT, ldda, dwork, lddwork, nb, queues );
    
    magma_tally4_queue_destroy( queues[0] );
    magma_tally4_queue_destroy( queues[1] );
}
