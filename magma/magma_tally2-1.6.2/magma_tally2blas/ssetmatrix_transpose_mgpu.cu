/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from zsetmatrix_transpose_mgpu.cu normal z -> s, Fri Jan 30 19:00:10 2015
       @author Ichitaro Yamazaki
*/
#include "common_magma_tally2.h"

#define PRECISION_s

//
//    m, n - dimensions in the source (input) matrix.
//             This routine copies the hA matrix from the CPU
//             to dAT on the GPU. In addition, the output matrix
//             is transposed. The routine uses a buffer of size
//             2*lddw*nb pointed to by dwork (lddw > m) on the GPU. 
//             Note that lda >= m and lddat >= n.
//
extern "C" void 
magma_tally2blas_ssetmatrix_transpose_mgpu(
    magma_tally2_int_t ngpu, magma_tally2_queue_t queues[][2],
    const float *hA,       magma_tally2_int_t lda, 
    magma_tally2Float_ptr    dAT[],    magma_tally2_int_t ldda, 
    magma_tally2Float_ptr    dwork[],  magma_tally2_int_t lddw,
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t nb)
{
#define hA(j)       (hA         + (j)*lda)
#define dwork(d, j) (dwork[(d)] + (j)*nb*lddw)
#define dAT(d, j)   (dAT[(d)]   + (j)*nb)

    magma_tally2_int_t nqueues = 2, d, j, j_local, id, ib;

    /* Quick return */
    if ( (m == 0) || (n == 0) )
        return;

    if (lda < m || ngpu*ldda < n || lddw < m){
        printf( "Wrong arguments in magma_tally2blas_ssetmatrix_transpose_mgpu (%d<%d), (%d*%d<%d), or (%d<%d).\n",
                (int) lda, (int) m, (int) ngpu, (int) ldda, (int) n, (int) lddw, (int) m );
        return;
    }
    
    /* Move data from CPU to GPU by block columns and transpose it */
    for(j=0; j<n; j+=nb){
       d       = (j/nb)%ngpu;
       j_local = (j/nb)/ngpu;
       id      = j_local%nqueues;
       magma_tally2_setdevice(d);

       ib = min(n-j, nb);
       magma_tally2_ssetmatrix_async( m, ib,
                               hA(j),        lda,
                               dwork(d, id), lddw, 
                               queues[d][id] );

       magma_tally2blas_stranspose_q( m, ib, dwork(d,id), lddw, dAT(d,j_local), ldda, queues[d][id] );
    }
}
