/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2011

       @author Azzam Haidar
       @author Tingxing Dong

       @precisions normal z -> s d c
*/



#include "common_magma_tally3.h"
#include "batched_kernel_param.h"



static    magma_tally3DoubleComplex neg_one = MAGMA_tally3_Z_NEG_ONE;
static    magma_tally3DoubleComplex one  = MAGMA_tally3_Z_ONE;
static    magma_tally3DoubleComplex zero  = MAGMA_tally3_Z_ZERO;

__global__ void
zgeqrf_copy_upper_kernel_batched(                
                  int n, int nb,
                  magma_tally3DoubleComplex **dV_array,    int ldv,
                  magma_tally3DoubleComplex **dR_array,    int ldr)
{

    magma_tally3DoubleComplex *dV = dV_array[blockIdx.x];
    magma_tally3DoubleComplex *dR = dR_array[blockIdx.x];

    int tid = threadIdx.x;

    int column = (tid / nb + 1) * nb; 
    
    if( tid < n && column < n) 
    {
       for(int i=column; i<n; i++)
       {
          dR[tid + i * ldr]  =  dV[tid + i * ldv];  
       }
    }
}

void zgeqrf_copy_upper_batched(                
                  magma_tally3_int_t n, magma_tally3_int_t nb,
                  magma_tally3DoubleComplex **dV_array,    magma_tally3_int_t ldv,
                  magma_tally3DoubleComplex **dR_array,    magma_tally3_int_t ldr,
          magma_tally3_int_t batchCount, magma_tally3_queue_t queue)
{
   /* 
        copy some data in dV to dR
   */

      if( nb >= n) return ;

      zgeqrf_copy_upper_kernel_batched<<<batchCount, n, 0, queue>>>(n, nb, dV_array, ldv, dR_array, ldr);

}



extern "C" magma_tally3_int_t
magma_tally3_zlarfb_zgemm_batched(
                  cublasHandle_t myhandle,
                  magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
                  magma_tally3DoubleComplex **dV_array,    magma_tally3_int_t ldv,
                  magma_tally3DoubleComplex **dT_array,    magma_tally3_int_t ldt,
                  magma_tally3DoubleComplex **dA_array,    magma_tally3_int_t lda,
                  magma_tally3DoubleComplex **W_array,     magma_tally3_int_t ldw,
                  magma_tally3DoubleComplex **W2_array,    magma_tally3_int_t ldw2,
                  magma_tally3_int_t batchCount, magma_tally3_queue_t queue)

{

    // W is workspace size of W is nb * n 
    // W = V^H * A. V is stored in A(i:m, i:ib)

    
    if( m <=0 || n <= 0 || k <=0 ) return 1;

#if 1  // CUBLAS is faster than MAGMA_tally3BLAS by 17GFLOP/S at size 512 batchCount = 2000
    cublasZgemmBatched(myhandle, CUBLAS_OP_C, CUBLAS_OP_N, k, n, m,
                             &one, (const magma_tally3DoubleComplex**) dV_array, ldv,
                                    (const magma_tally3DoubleComplex**) dA_array, lda,
                             &zero,  W_array, ldw, batchCount );



    // W2 = T^H * W        
    cublasZgemmBatched(myhandle, CUBLAS_OP_C, CUBLAS_OP_N, k, n, k,
                             &one, (const magma_tally3DoubleComplex**) dT_array, ldt,
                                    (const magma_tally3DoubleComplex**) W_array, ldw,
                             &zero,  W2_array, ldw2, batchCount );

        
    // A = A - V * W2 
    cublasZgemmBatched(myhandle, CUBLAS_OP_N, CUBLAS_OP_N, m, n, k,
                             &neg_one, (const magma_tally3DoubleComplex**) dV_array, ldv,
                                    (const magma_tally3DoubleComplex**) W2_array, ldw2,
                             &one,  dA_array, lda, batchCount );

#else 

    magma_tally3blas_zgemm_batched(Magma_tally3ConjTrans, Magma_tally3NoTrans, k, n, m,
                             one, (const magma_tally3DoubleComplex**) dV_array, ldv,
                                    (const magma_tally3DoubleComplex**) dA_array, lda,
                             zero,  W_array, ldw, batchCount );



    // W2 = T^H * W        
    magma_tally3blas_zgemm_batched(Magma_tally3ConjTrans, Magma_tally3NoTrans, k, n, k,
                             one, (const magma_tally3DoubleComplex**) dT_array, ldt,
                                    (const magma_tally3DoubleComplex**) W_array, ldw,
                             zero,  W2_array, ldw2, batchCount );

        
    // A = A - V * W2 
    magma_tally3blas_zgemm_batched(Magma_tally3NoTrans, Magma_tally3NoTrans, m, n, k,
                             neg_one, (const magma_tally3DoubleComplex**) dV_array, ldv,
                                    (const magma_tally3DoubleComplex**) W2_array, ldw2,
                             one,  dA_array, lda, batchCount );
          
#endif       
    return 0;

}



