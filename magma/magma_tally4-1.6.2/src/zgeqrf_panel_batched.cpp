/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2014
       
       @author Azzam Haidar
       @author Tingxing Dong

       @precisions normal z -> s d c
*/
#include "common_magma_tally4.h"



////////////////////////////////////////////////////////////////////////////////////////
extern "C" magma_tally4_int_t
magma_tally4_zgeqrf_panel_batched(
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t nb,    
    magma_tally4DoubleComplex** dA_array,    magma_tally4_int_t ldda,
    magma_tally4DoubleComplex** tau_array, 
    magma_tally4DoubleComplex** dT_array, magma_tally4_int_t ldt, 
    magma_tally4DoubleComplex** dR_array, magma_tally4_int_t ldr,
    magma_tally4DoubleComplex** dW0_displ, 
    magma_tally4DoubleComplex** dW1_displ,
    magma_tally4DoubleComplex   *dwork,  
    magma_tally4DoubleComplex** dW2_displ, 
    magma_tally4DoubleComplex** dW3_displ,
    magma_tally4_int_t *info_array,
    magma_tally4_int_t batchCount, cublasHandle_t myhandle, magma_tally4_queue_t queue)
{

            magma_tally4_int_t j, jb;
            magma_tally4_int_t ldw = nb; 

            for( j=0; j<n; j+=nb)
            {
                jb = min(nb, n-j);

                magma_tally4_zdisplace_pointers(dW0_displ, dA_array, ldda, j, j, batchCount, queue); 
                magma_tally4_zdisplace_pointers(dW2_displ, tau_array, 1, j, 0, batchCount, queue);
                magma_tally4_zdisplace_pointers(dW3_displ, dR_array, ldr, j, j, batchCount, queue); // 
  
                //sub-panel factorization 
                magma_tally4_zgeqr2_batched(
                      m-j, jb,
                      dW0_displ, ldda,      
                      dW2_displ, 
                      info_array, 
                      batchCount,
                      queue);

                //copy upper part of dA to dR
                magma_tally4_zdisplace_pointers(dW0_displ, dA_array, ldda, j, j, batchCount, queue); 
                magma_tally4_zdisplace_pointers(dW3_displ, dR_array, ldr, j, j, batchCount, queue); 
                magma_tally4blas_zlacpy_batched(Magma_tally4Upper, jb, jb, dW0_displ, ldda, dW3_displ, ldr, batchCount, queue);

                magma_tally4_zdisplace_pointers(dW0_displ, dA_array, ldda, j, j, batchCount, queue); 
                magma_tally4_zdisplace_pointers(dW3_displ, dR_array, ldr, j, j, batchCount, queue);
                magma_tally4blas_zlaset_batched(Magma_tally4Upper, jb, jb, MAGMA_tally4_Z_ZERO, MAGMA_tally4_Z_ONE, dW0_displ, ldda, batchCount, queue); 

                
                if( (n-j-jb) > 0) //update the trailing matrix inside the panel
                {

                    magma_tally4_zlarft_sm32x32_batched(m-j, jb,
                                 dW0_displ, ldda,
                                 dW2_displ,
                                 dT_array, ldt, 
                                 batchCount, myhandle, queue);

                    magma_tally4_zdisplace_pointers(dW1_displ, dA_array, ldda, j, j + jb, batchCount, queue); 
                    zset_pointer(dW2_displ,  dwork, 1, 0, 0,  ldw*n, batchCount, queue );
                    zset_pointer(dW3_displ, dwork + ldw*n*batchCount, 1, 0, 0,  ldw*n, batchCount, queue );    

                    magma_tally4_zlarfb_gemm_batched(
                                Magma_tally4Left, Magma_tally4ConjTrans, Magma_tally4Forward, Magma_tally4Columnwise,
                                m-j, n-j-jb, jb,
                                (const magma_tally4DoubleComplex**)dW0_displ, ldda,
                                (const magma_tally4DoubleComplex**)dT_array, ldt,
                                dW1_displ,  ldda,
                                dW2_displ,  ldw, 
                                dW3_displ, ldw, batchCount, myhandle, queue);
                }
               
            }

    return 0;
}


////////////////////////////////////////////////////////////////////////////////////////

