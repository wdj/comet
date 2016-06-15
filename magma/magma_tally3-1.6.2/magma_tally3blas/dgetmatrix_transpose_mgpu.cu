/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from zgetmatrix_transpose_mgpu.cu normal z -> d, Fri Jan 30 19:00:09 2015
       @author Ichitaro Yamazaki
*/
#include "common_magma_tally3.h"

#define PRECISION_d

//
//    m, n - dimensions in the output (hA) matrix.
//             This routine copies the dAT matrix from the GPU
//             to hA on the CPU. In addition, the output matrix
//             is transposed. The routine uses a buffer of size
//             2*lddw*nb pointed to by dwork (lddw > m) on the GPU. 
//             Note that lda >= m and lddat >= n.
//
extern "C" void 
magma_tally3blas_dgetmatrix_transpose_mgpu(
    magma_tally3_int_t ngpu, magma_tally3_queue_t queues[][2],
    magma_tally3Double_const_ptr const dAT[],   magma_tally3_int_t ldda,
    double                *hA,      magma_tally3_int_t lda,
    magma_tally3Double_ptr             dwork[], magma_tally3_int_t lddw,
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t nb)
{
#define hA(j)       (hA         + (j)*lda)
#define dwork(d, j) (dwork[(d)] + (j)*nb*lddw)
#define dAT(d, j)   (dAT[(d)]   + (j)*nb)

    magma_tally3_int_t nstreams = 2, d, j, j_local, id, ib;

    /* Quick return */
    if ( (m == 0) || (n == 0) )
        return;

    if (lda < m || ngpu*ldda < n || lddw < m){
        printf( "Wrong arguments in magma_tally3blas_dgetmatrix_transpose_mgpu (%d<%d), (%d*%d<%d), or (%d<%d).\n",
                (int) lda, (int) m, (int) ngpu, (int) ldda, (int) n, (int) lddw, (int) m );
        return;
    }
    
    /* Move data from GPU to CPU using two buffers; first transpose the data on the GPU */
    for(j=0; j<n; j+=nb){
       d       = (j/nb)%ngpu;
       j_local = (j/nb)/ngpu;
       id      = j_local%nstreams;
       magma_tally3_setdevice(d);

       ib = min(n-j, nb);
       magma_tally3blas_dtranspose_q( ib, m, dAT(d,j_local), ldda, dwork(d,id), lddw, queues[d][id] );
       magma_tally3_dgetmatrix_async( m, ib,
                               dwork(d, id), lddw,
                               hA(j),        lda, 
                               queues[d][id] );
    }
}