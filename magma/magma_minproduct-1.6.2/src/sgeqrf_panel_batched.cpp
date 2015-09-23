/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2014
       
       @author Azzam Haidar
       @author Tingxing Dong

       @generated from zgeqrf_panel_batched.cpp normal z -> s, Fri Jan 30 19:00:19 2015
*/
#include "common_magma_minproduct.h"



////////////////////////////////////////////////////////////////////////////////////////
extern "C" magma_minproduct_int_t
magma_minproduct_sgeqrf_panel_batched(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t nb,    
    float** dA_array,    magma_minproduct_int_t ldda,
    float** tau_array, 
    float** dT_array, magma_minproduct_int_t ldt, 
    float** dR_array, magma_minproduct_int_t ldr,
    float** dW0_displ, 
    float** dW1_displ,
    float   *dwork,  
    float** dW2_displ, 
    float** dW3_displ,
    magma_minproduct_int_t *info_array,
    magma_minproduct_int_t batchCount, cublasHandle_t myhandle, magma_minproduct_queue_t queue)
{

            magma_minproduct_int_t j, jb;
            magma_minproduct_int_t ldw = nb; 

            for( j=0; j<n; j+=nb)
            {
                jb = min(nb, n-j);

                magma_minproduct_sdisplace_pointers(dW0_displ, dA_array, ldda, j, j, batchCount, queue); 
                magma_minproduct_sdisplace_pointers(dW2_displ, tau_array, 1, j, 0, batchCount, queue);
                magma_minproduct_sdisplace_pointers(dW3_displ, dR_array, ldr, j, j, batchCount, queue); // 
  
                //sub-panel factorization 
                magma_minproduct_sgeqr2_batched(
                      m-j, jb,
                      dW0_displ, ldda,      
                      dW2_displ, 
                      info_array, 
                      batchCount,
                      queue);

                //copy upper part of dA to dR
                magma_minproduct_sdisplace_pointers(dW0_displ, dA_array, ldda, j, j, batchCount, queue); 
                magma_minproduct_sdisplace_pointers(dW3_displ, dR_array, ldr, j, j, batchCount, queue); 
                magma_minproductblas_slacpy_batched(Magma_minproductUpper, jb, jb, dW0_displ, ldda, dW3_displ, ldr, batchCount, queue);

                magma_minproduct_sdisplace_pointers(dW0_displ, dA_array, ldda, j, j, batchCount, queue); 
                magma_minproduct_sdisplace_pointers(dW3_displ, dR_array, ldr, j, j, batchCount, queue);
                magma_minproductblas_slaset_batched(Magma_minproductUpper, jb, jb, MAGMA_minproduct_S_ZERO, MAGMA_minproduct_S_ONE, dW0_displ, ldda, batchCount, queue); 

                
                if( (n-j-jb) > 0) //update the trailing matrix inside the panel
                {

                    magma_minproduct_slarft_sm32x32_batched(m-j, jb,
                                 dW0_displ, ldda,
                                 dW2_displ,
                                 dT_array, ldt, 
                                 batchCount, myhandle, queue);

                    magma_minproduct_sdisplace_pointers(dW1_displ, dA_array, ldda, j, j + jb, batchCount, queue); 
                    sset_pointer(dW2_displ,  dwork, 1, 0, 0,  ldw*n, batchCount, queue );
                    sset_pointer(dW3_displ, dwork + ldw*n*batchCount, 1, 0, 0,  ldw*n, batchCount, queue );    

                    magma_minproduct_slarfb_gemm_batched(
                                Magma_minproductLeft, Magma_minproductConjTrans, Magma_minproductForward, Magma_minproductColumnwise,
                                m-j, n-j-jb, jb,
                                (const float**)dW0_displ, ldda,
                                (const float**)dT_array, ldt,
                                dW1_displ,  ldda,
                                dW2_displ,  ldw, 
                                dW3_displ, ldw, batchCount, myhandle, queue);
                }
               
            }

    return 0;
}


////////////////////////////////////////////////////////////////////////////////////////

