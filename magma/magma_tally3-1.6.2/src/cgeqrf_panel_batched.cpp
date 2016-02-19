/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2014
       
       @author Azzam Haidar
       @author Tingxing Dong

       @generated from zgeqrf_panel_batched.cpp normal z -> c, Fri Jan 30 19:00:19 2015
*/
#include "common_magma_tally3.h"



////////////////////////////////////////////////////////////////////////////////////////
extern "C" magma_tally3_int_t
magma_tally3_cgeqrf_panel_batched(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t nb,    
    magma_tally3FloatComplex** dA_array,    magma_tally3_int_t ldda,
    magma_tally3FloatComplex** tau_array, 
    magma_tally3FloatComplex** dT_array, magma_tally3_int_t ldt, 
    magma_tally3FloatComplex** dR_array, magma_tally3_int_t ldr,
    magma_tally3FloatComplex** dW0_displ, 
    magma_tally3FloatComplex** dW1_displ,
    magma_tally3FloatComplex   *dwork,  
    magma_tally3FloatComplex** dW2_displ, 
    magma_tally3FloatComplex** dW3_displ,
    magma_tally3_int_t *info_array,
    magma_tally3_int_t batchCount, cublasHandle_t myhandle, magma_tally3_queue_t queue)
{

            magma_tally3_int_t j, jb;
            magma_tally3_int_t ldw = nb; 

            for( j=0; j<n; j+=nb)
            {
                jb = min(nb, n-j);

                magma_tally3_cdisplace_pointers(dW0_displ, dA_array, ldda, j, j, batchCount, queue); 
                magma_tally3_cdisplace_pointers(dW2_displ, tau_array, 1, j, 0, batchCount, queue);
                magma_tally3_cdisplace_pointers(dW3_displ, dR_array, ldr, j, j, batchCount, queue); // 
  
                //sub-panel factorization 
                magma_tally3_cgeqr2_batched(
                      m-j, jb,
                      dW0_displ, ldda,      
                      dW2_displ, 
                      info_array, 
                      batchCount,
                      queue);

                //copy upper part of dA to dR
                magma_tally3_cdisplace_pointers(dW0_displ, dA_array, ldda, j, j, batchCount, queue); 
                magma_tally3_cdisplace_pointers(dW3_displ, dR_array, ldr, j, j, batchCount, queue); 
                magma_tally3blas_clacpy_batched(Magma_tally3Upper, jb, jb, dW0_displ, ldda, dW3_displ, ldr, batchCount, queue);

                magma_tally3_cdisplace_pointers(dW0_displ, dA_array, ldda, j, j, batchCount, queue); 
                magma_tally3_cdisplace_pointers(dW3_displ, dR_array, ldr, j, j, batchCount, queue);
                magma_tally3blas_claset_batched(Magma_tally3Upper, jb, jb, MAGMA_tally3_C_ZERO, MAGMA_tally3_C_ONE, dW0_displ, ldda, batchCount, queue); 

                
                if( (n-j-jb) > 0) //update the trailing matrix inside the panel
                {

                    magma_tally3_clarft_sm32x32_batched(m-j, jb,
                                 dW0_displ, ldda,
                                 dW2_displ,
                                 dT_array, ldt, 
                                 batchCount, myhandle, queue);

                    magma_tally3_cdisplace_pointers(dW1_displ, dA_array, ldda, j, j + jb, batchCount, queue); 
                    cset_pointer(dW2_displ,  dwork, 1, 0, 0,  ldw*n, batchCount, queue );
                    cset_pointer(dW3_displ, dwork + ldw*n*batchCount, 1, 0, 0,  ldw*n, batchCount, queue );    

                    magma_tally3_clarfb_gemm_batched(
                                Magma_tally3Left, Magma_tally3ConjTrans, Magma_tally3Forward, Magma_tally3Columnwise,
                                m-j, n-j-jb, jb,
                                (const magma_tally3FloatComplex**)dW0_displ, ldda,
                                (const magma_tally3FloatComplex**)dT_array, ldt,
                                dW1_displ,  ldda,
                                dW2_displ,  ldw, 
                                dW3_displ, ldw, batchCount, myhandle, queue);
                }
               
            }

    return 0;
}


////////////////////////////////////////////////////////////////////////////////////////

