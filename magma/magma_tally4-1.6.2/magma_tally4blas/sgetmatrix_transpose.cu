/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from zgetmatrix_transpose.cu normal z -> s, Fri Jan 30 19:00:09 2015

*/
#include "common_magma_tally4.h"

#define PRECISION_s


//
//      m, n - dimensions in the output (hA) matrix.
//             This routine copies the dAT matrix from the GPU
//             to hA on the CPU. In addition, the output matrix
//             is transposed. The routine uses a buffer of size
//             2*lddwork*nb pointed to by dwork (lddwork > m) on the GPU. 
//             Note that lda >= m and lddat >= n.
//
extern "C" void 
magma_tally4blas_sgetmatrix_transpose_q(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4Float_const_ptr dAT, magma_tally4_int_t ldda,
    float          *hA,  magma_tally4_int_t lda,
    magma_tally4Float_ptr       dwork,  magma_tally4_int_t lddwork, magma_tally4_int_t nb,
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
        printf("Wrong arguments in sgetmatrix_transpose.\n");
        return;
    }

    for(i=0; i < n; i += nb) {
        /* Move data from GPU to CPU using 2 buffers; 1st transpose the data on the GPU */
        ib = min(n-i, nb);
        
        magma_tally4blas_stranspose_q( ib, m, dAT(i,0), ldda, dwork(0,(j%2)*nb), lddwork, queues[j%2] );
        magma_tally4_sgetmatrix_async( m, ib,
                                dwork(0,(j%2)*nb), lddwork,
                                hA(0,i), lda, queues[j%2] );
        j++;
    }
}


// @see magma_tally4blas_sgetmatrix_transpose_q
extern "C" void 
magma_tally4blas_sgetmatrix_transpose(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4Float_const_ptr dAT, magma_tally4_int_t ldda,
    float          *hA,  magma_tally4_int_t lda,
    magma_tally4Float_ptr       dwork,  magma_tally4_int_t lddwork, magma_tally4_int_t nb )
{
    magma_tally4_queue_t queues[2];
    magma_tally4_queue_create( &queues[0] );
    magma_tally4_queue_create( &queues[1] );

    magma_tally4blas_sgetmatrix_transpose_q( m, n, dAT, ldda, hA, lda, dwork, lddwork, nb, queues );

    magma_tally4_queue_destroy( queues[0] );
    magma_tally4_queue_destroy( queues[1] );
}
