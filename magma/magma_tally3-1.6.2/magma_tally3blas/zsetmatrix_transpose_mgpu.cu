/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions normal z -> s d c
       @author Ichitaro Yamazaki
*/
#include "common_magma_tally3.h"

#define PRECISION_z

//
//    m, n - dimensions in the source (input) matrix.
//             This routine copies the hA matrix from the CPU
//             to dAT on the GPU. In addition, the output matrix
//             is transposed. The routine uses a buffer of size
//             2*lddw*nb pointed to by dwork (lddw > m) on the GPU. 
//             Note that lda >= m and lddat >= n.
//
extern "C" void 
magma_tally3blas_zsetmatrix_transpose_mgpu(
    magma_tally3_int_t ngpu, magma_tally3_queue_t queues[][2],
    const magma_tally3DoubleComplex *hA,       magma_tally3_int_t lda, 
    magma_tally3DoubleComplex_ptr    dAT[],    magma_tally3_int_t ldda, 
    magma_tally3DoubleComplex_ptr    dwork[],  magma_tally3_int_t lddw,
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t nb)
{
#define hA(j)       (hA         + (j)*lda)
#define dwork(d, j) (dwork[(d)] + (j)*nb*lddw)
#define dAT(d, j)   (dAT[(d)]   + (j)*nb)

    magma_tally3_int_t nqueues = 2, d, j, j_local, id, ib;

    /* Quick return */
    if ( (m == 0) || (n == 0) )
        return;

    if (lda < m || ngpu*ldda < n || lddw < m){
        printf( "Wrong arguments in magma_tally3blas_zsetmatrix_transpose_mgpu (%d<%d), (%d*%d<%d), or (%d<%d).\n",
                (int) lda, (int) m, (int) ngpu, (int) ldda, (int) n, (int) lddw, (int) m );
        return;
    }
    
    /* Move data from CPU to GPU by block columns and transpose it */
    for(j=0; j<n; j+=nb){
       d       = (j/nb)%ngpu;
       j_local = (j/nb)/ngpu;
       id      = j_local%nqueues;
       magma_tally3_setdevice(d);

       ib = min(n-j, nb);
       magma_tally3_zsetmatrix_async( m, ib,
                               hA(j),        lda,
                               dwork(d, id), lddw, 
                               queues[d][id] );

       magma_tally3blas_ztranspose_q( m, ib, dwork(d,id), lddw, dAT(d,j_local), ldda, queues[d][id] );
    }
}