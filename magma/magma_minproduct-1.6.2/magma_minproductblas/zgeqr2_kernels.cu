/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2011

       @author Azzam Haidar
       @author Tingxing Dong

       @precisions normal z -> s d c
*/



#include "common_magma_minproduct.h"
#include "batched_kernel_param.h"



static    magma_minproductDoubleComplex neg_one = MAGMA_minproduct_Z_NEG_ONE;
static    magma_minproductDoubleComplex one  = MAGMA_minproduct_Z_ONE;
static    magma_minproductDoubleComplex zero  = MAGMA_minproduct_Z_ZERO;

__global__ void
zgeqrf_copy_upper_kernel_batched(                
                  int n, int nb,
                  magma_minproductDoubleComplex **dV_array,    int ldv,
                  magma_minproductDoubleComplex **dR_array,    int ldr)
{

    magma_minproductDoubleComplex *dV = dV_array[blockIdx.x];
    magma_minproductDoubleComplex *dR = dR_array[blockIdx.x];

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
                  magma_minproduct_int_t n, magma_minproduct_int_t nb,
                  magma_minproductDoubleComplex **dV_array,    magma_minproduct_int_t ldv,
                  magma_minproductDoubleComplex **dR_array,    magma_minproduct_int_t ldr,
          magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue)
{
   /* 
        copy some data in dV to dR
   */

      if( nb >= n) return ;

      zgeqrf_copy_upper_kernel_batched<<<batchCount, n, 0, queue>>>(n, nb, dV_array, ldv, dR_array, ldr);

}



extern "C" magma_minproduct_int_t
magma_minproduct_zlarfb_zgemm_batched(
                  cublasHandle_t myhandle,
                  magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
                  magma_minproductDoubleComplex **dV_array,    magma_minproduct_int_t ldv,
                  magma_minproductDoubleComplex **dT_array,    magma_minproduct_int_t ldt,
                  magma_minproductDoubleComplex **dA_array,    magma_minproduct_int_t lda,
                  magma_minproductDoubleComplex **W_array,     magma_minproduct_int_t ldw,
                  magma_minproductDoubleComplex **W2_array,    magma_minproduct_int_t ldw2,
                  magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue)

{

    // W is workspace size of W is nb * n 
    // W = V^H * A. V is stored in A(i:m, i:ib)

    
    if( m <=0 || n <= 0 || k <=0 ) return 1;

#if 1  // CUBLAS is faster than MAGMA_minproductBLAS by 17GFLOP/S at size 512 batchCount = 2000
    cublasZgemmBatched(myhandle, CUBLAS_OP_C, CUBLAS_OP_N, k, n, m,
                             &one, (const magma_minproductDoubleComplex**) dV_array, ldv,
                                    (const magma_minproductDoubleComplex**) dA_array, lda,
                             &zero,  W_array, ldw, batchCount );



    // W2 = T^H * W        
    cublasZgemmBatched(myhandle, CUBLAS_OP_C, CUBLAS_OP_N, k, n, k,
                             &one, (const magma_minproductDoubleComplex**) dT_array, ldt,
                                    (const magma_minproductDoubleComplex**) W_array, ldw,
                             &zero,  W2_array, ldw2, batchCount );

        
    // A = A - V * W2 
    cublasZgemmBatched(myhandle, CUBLAS_OP_N, CUBLAS_OP_N, m, n, k,
                             &neg_one, (const magma_minproductDoubleComplex**) dV_array, ldv,
                                    (const magma_minproductDoubleComplex**) W2_array, ldw2,
                             &one,  dA_array, lda, batchCount );

#else 

    magma_minproductblas_zgemm_batched(Magma_minproductConjTrans, Magma_minproductNoTrans, k, n, m,
                             one, (const magma_minproductDoubleComplex**) dV_array, ldv,
                                    (const magma_minproductDoubleComplex**) dA_array, lda,
                             zero,  W_array, ldw, batchCount );



    // W2 = T^H * W        
    magma_minproductblas_zgemm_batched(Magma_minproductConjTrans, Magma_minproductNoTrans, k, n, k,
                             one, (const magma_minproductDoubleComplex**) dT_array, ldt,
                                    (const magma_minproductDoubleComplex**) W_array, ldw,
                             zero,  W2_array, ldw2, batchCount );

        
    // A = A - V * W2 
    magma_minproductblas_zgemm_batched(Magma_minproductNoTrans, Magma_minproductNoTrans, m, n, k,
                             neg_one, (const magma_minproductDoubleComplex**) dV_array, ldv,
                                    (const magma_minproductDoubleComplex**) W2_array, ldw2,
                             one,  dA_array, lda, batchCount );
          
#endif       
    return 0;

}



