/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2011

       @author Azzam Haidar
       @author Tingxing Dong

       @generated from zgeqr2_kernels.cu normal z -> s, Fri Jan 30 19:00:10 2015
*/



#include "common_magma_tally2.h"
#include "batched_kernel_param.h"



static    float neg_one = MAGMA_tally2_S_NEG_ONE;
static    float one  = MAGMA_tally2_S_ONE;
static    float zero  = MAGMA_tally2_S_ZERO;

__global__ void
sgeqrf_copy_upper_kernel_batched(                
                  int n, int nb,
                  float **dV_array,    int ldv,
                  float **dR_array,    int ldr)
{

    float *dV = dV_array[blockIdx.x];
    float *dR = dR_array[blockIdx.x];

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

void sgeqrf_copy_upper_batched(                
                  magma_tally2_int_t n, magma_tally2_int_t nb,
                  float **dV_array,    magma_tally2_int_t ldv,
                  float **dR_array,    magma_tally2_int_t ldr,
          magma_tally2_int_t batchCount, magma_tally2_queue_t queue)
{
   /* 
        copy some data in dV to dR
   */

      if( nb >= n) return ;

      sgeqrf_copy_upper_kernel_batched<<<batchCount, n, 0, queue>>>(n, nb, dV_array, ldv, dR_array, ldr);

}



extern "C" magma_tally2_int_t
magma_tally2_slarfb_sgemm_batched(
                  cublasHandle_t myhandle,
                  magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
                  float **dV_array,    magma_tally2_int_t ldv,
                  float **dT_array,    magma_tally2_int_t ldt,
                  float **dA_array,    magma_tally2_int_t lda,
                  float **W_array,     magma_tally2_int_t ldw,
                  float **W2_array,    magma_tally2_int_t ldw2,
                  magma_tally2_int_t batchCount, magma_tally2_queue_t queue)

{

    // W is workspace size of W is nb * n 
    // W = V^H * A. V is stored in A(i:m, i:ib)

    
    if( m <=0 || n <= 0 || k <=0 ) return 1;

#if 1  // CUBLAS is faster than MAGMA_tally2BLAS by 17GFLOP/S at size 512 batchCount = 2000
    cublasSgemmBatched(myhandle, CUBLAS_OP_C, CUBLAS_OP_N, k, n, m,
                             &one, (const float**) dV_array, ldv,
                                    (const float**) dA_array, lda,
                             &zero,  W_array, ldw, batchCount );



    // W2 = T^H * W        
    cublasSgemmBatched(myhandle, CUBLAS_OP_C, CUBLAS_OP_N, k, n, k,
                             &one, (const float**) dT_array, ldt,
                                    (const float**) W_array, ldw,
                             &zero,  W2_array, ldw2, batchCount );

        
    // A = A - V * W2 
    cublasSgemmBatched(myhandle, CUBLAS_OP_N, CUBLAS_OP_N, m, n, k,
                             &neg_one, (const float**) dV_array, ldv,
                                    (const float**) W2_array, ldw2,
                             &one,  dA_array, lda, batchCount );

#else 

    magma_tally2blas_sgemm_batched(Magma_tally2ConjTrans, Magma_tally2NoTrans, k, n, m,
                             one, (const float**) dV_array, ldv,
                                    (const float**) dA_array, lda,
                             zero,  W_array, ldw, batchCount );



    // W2 = T^H * W        
    magma_tally2blas_sgemm_batched(Magma_tally2ConjTrans, Magma_tally2NoTrans, k, n, k,
                             one, (const float**) dT_array, ldt,
                                    (const float**) W_array, ldw,
                             zero,  W2_array, ldw2, batchCount );

        
    // A = A - V * W2 
    magma_tally2blas_sgemm_batched(Magma_tally2NoTrans, Magma_tally2NoTrans, m, n, k,
                             neg_one, (const float**) dV_array, ldv,
                                    (const float**) W2_array, ldw2,
                             one,  dA_array, lda, batchCount );
          
#endif       
    return 0;

}



