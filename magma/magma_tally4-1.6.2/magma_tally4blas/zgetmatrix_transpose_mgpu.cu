/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions normal z -> s d c
       @author Ichitaro Yamazaki
*/
#include "common_magma_tally4.h"

#define PRECISION_z

//
//    m, n - dimensions in the output (hA) matrix.
//             This routine copies the dAT matrix from the GPU
//             to hA on the CPU. In addition, the output matrix
//             is transposed. The routine uses a buffer of size
//             2*lddw*nb pointed to by dwork (lddw > m) on the GPU. 
//             Note that lda >= m and lddat >= n.
//
extern "C" void 
magma_tally4blas_zgetmatrix_transpose_mgpu(
    magma_tally4_int_t ngpu, magma_tally4_queue_t queues[][2],
    magma_tally4DoubleComplex_const_ptr const dAT[],   magma_tally4_int_t ldda,
    magma_tally4DoubleComplex                *hA,      magma_tally4_int_t lda,
    magma_tally4DoubleComplex_ptr             dwork[], magma_tally4_int_t lddw,
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t nb)
{
#define hA(j)       (hA         + (j)*lda)
#define dwork(d, j) (dwork[(d)] + (j)*nb*lddw)
#define dAT(d, j)   (dAT[(d)]   + (j)*nb)

    magma_tally4_int_t nstreams = 2, d, j, j_local, id, ib;

    /* Quick return */
    if ( (m == 0) || (n == 0) )
        return;

    if (lda < m || ngpu*ldda < n || lddw < m){
        printf( "Wrong arguments in magma_tally4blas_zgetmatrix_transpose_mgpu (%d<%d), (%d*%d<%d), or (%d<%d).\n",
                (int) lda, (int) m, (int) ngpu, (int) ldda, (int) n, (int) lddw, (int) m );
        return;
    }
    
    /* Move data from GPU to CPU using two buffers; first transpose the data on the GPU */
    for(j=0; j<n; j+=nb){
       d       = (j/nb)%ngpu;
       j_local = (j/nb)/ngpu;
       id      = j_local%nstreams;
       magma_tally4_setdevice(d);

       ib = min(n-j, nb);
       magma_tally4blas_ztranspose_q( ib, m, dAT(d,j_local), ldda, dwork(d,id), lddw, queues[d][id] );
       magma_tally4_zgetmatrix_async( m, ib,
                               dwork(d, id), lddw,
                               hA(j),        lda, 
                               queues[d][id] );
    }
}
