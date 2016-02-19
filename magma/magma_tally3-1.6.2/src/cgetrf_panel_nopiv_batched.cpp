/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2013
       
   @author Azzam Haidar
   @author Adrien Remy

   @generated from zgetrf_panel_nopiv_batched.cpp normal z -> c, Fri Jan 30 19:00:19 2015
*/
#include "common_magma_tally3.h"
////////////////////////////////////////////////////////////////////////////////////////
extern "C" magma_tally3_int_t
magma_tally3_cgetrf_panel_nopiv_batched(
    magma_tally3_int_t m, magma_tally3_int_t nb,    
    magma_tally3FloatComplex** dA_array,    magma_tally3_int_t ldda,
    magma_tally3FloatComplex** dX_array,    magma_tally3_int_t dX_length,
    magma_tally3FloatComplex** dinvA_array, magma_tally3_int_t dinvA_length,
    magma_tally3FloatComplex** dW0_displ, magma_tally3FloatComplex** dW1_displ,  
    magma_tally3FloatComplex** dW2_displ, magma_tally3FloatComplex** dW3_displ,
    magma_tally3FloatComplex** dW4_displ,     
    magma_tally3_int_t *info_array, magma_tally3_int_t gbstep,  
    magma_tally3_int_t batchCount, cublasHandle_t myhandle, magma_tally3_queue_t queue )
{
     magma_tally3_int_t arginfo = 0;
    //===============================================
    //  panel factorization
    //===============================================
    if(m < nb){
        printf("magma_tally3_cgetrf_panel_nopiv_batched_q m < nb %d < %d \n",(int) m, (int) nb);
        return -101;
    }

#if 0
    arginfo = magma_tally3_cgetf2_nopiv_batched(
                       m, nb,
                       dA_array, ldda,
                       dW1_displ, dW2_displ, dW3_displ,
                       info_array, gbstep, batchCount, myhandle);
    if (arginfo != 0) return arginfo;
#else
    arginfo = magma_tally3_cgetf2_nopiv_batched(
                       nb, nb,
                       dA_array, ldda,
                       dW1_displ, dW2_displ, dW3_displ,
                       info_array, gbstep, batchCount, myhandle, queue);
    if (arginfo != 0) return arginfo;
    if((m-nb) > 0){
        magma_tally3_cdisplace_pointers(dW0_displ, dA_array, ldda, nb, 0, batchCount, queue);
        magma_tally3blas_ctrsm_work_batched(Magma_tally3Right, Magma_tally3Upper, Magma_tally3NoTrans, Magma_tally3NonUnit,
                              1, m-nb, nb, 
                              MAGMA_tally3_C_ONE,
                              dA_array,    ldda, 
                              dW0_displ,   ldda, 
                              dX_array,    m-nb, 
                              dinvA_array, dinvA_length, 
                              dW1_displ,   dW2_displ, 
                              dW3_displ,   dW4_displ,
                              1, batchCount, queue);
    }
#endif
    return 0;
}


////////////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////////////
extern "C" magma_tally3_int_t
magma_tally3_cgetrf_recpanel_nopiv_batched(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t min_recpnb,    
    magma_tally3FloatComplex** dA_array,    magma_tally3_int_t ldda,
    magma_tally3FloatComplex** dX_array,    magma_tally3_int_t dX_length,
    magma_tally3FloatComplex** dinvA_array, magma_tally3_int_t dinvA_length,
    magma_tally3FloatComplex** dW1_displ, magma_tally3FloatComplex** dW2_displ,  
    magma_tally3FloatComplex** dW3_displ, magma_tally3FloatComplex** dW4_displ,
    magma_tally3FloatComplex** dW5_displ, 
    magma_tally3_int_t *info_array, magma_tally3_int_t gbstep,
    magma_tally3_int_t batchCount, cublasHandle_t myhandle, magma_tally3_queue_t queue )
{
    // Quick return if possible
    if (m == 0 || n == 0) {
        return 0;
    }
    magma_tally3_int_t arginfo = 0;


    magma_tally3FloatComplex **dA_displ  = NULL;
    magma_tally3_malloc((void**)&dA_displ,   batchCount * sizeof(*dA_displ));
    
    magma_tally3_int_t panel_nb = n;
    if(panel_nb <= min_recpnb){
        // if(DEBUG>0)printf("calling bottom panel recursive with m=%d nb=%d\n",m,n);
        //  panel factorization
        //magma_tally3_cdisplace_pointers(dA_displ, dA_array, ldda, 0, 0, batchCount);
        arginfo = magma_tally3_cgetrf_panel_nopiv_batched(
                           m, panel_nb, 
                           dA_array, ldda,
                           dX_array, dX_length,
                           dinvA_array, dinvA_length,
                           dW1_displ, dW2_displ,
                           dW3_displ, dW4_displ, dW5_displ,
                           info_array, gbstep, batchCount, myhandle, queue);
        if (arginfo != 0) return arginfo;
    }
    else{
        // split A over two [A A2]
        // panel on A1, update on A2 then panel on A1    
        magma_tally3_int_t n1 = n/2;
        magma_tally3_int_t n2 = n-n1;
        magma_tally3_int_t m1 = m;
        magma_tally3_int_t m2 = m-n1;
        magma_tally3_int_t p1 = 0;
        magma_tally3_int_t p2 = n1;
        // panel on A1
        //printf("calling recursive panel on A1 with m=%d nb=%d min_recpnb %d\n",m1,n1,min_recpnb);
        magma_tally3_cdisplace_pointers(dA_displ, dA_array, ldda, p1, p1, batchCount, queue); 
        arginfo = magma_tally3_cgetrf_recpanel_nopiv_batched(
                           m1, n1, min_recpnb,
                           dA_displ, ldda,
                           dX_array, dX_length,
                           dinvA_array, dinvA_length,
                           dW1_displ, dW2_displ,
                           dW3_displ, dW4_displ, dW5_displ,
                           info_array, gbstep, batchCount, myhandle, queue);
        if (arginfo != 0) return arginfo;

        // update A2
        //printf("calling update A2 with             m=%d n=%d k=%d\n",m2,n2,n1);
        
        magma_tally3_cdisplace_pointers(dW5_displ, dA_array, ldda, p1, p2, batchCount, queue); 
        magma_tally3blas_ctrsm_work_batched(Magma_tally3Left, Magma_tally3Lower, Magma_tally3NoTrans, Magma_tally3Unit, 1,
                              n1, n2,
                              MAGMA_tally3_C_ONE,
                              dA_displ,    ldda, // dA
                              dW5_displ,   ldda, // dB
                              dX_array,  n1, // dX
                              dinvA_array, dinvA_length,
                              dW1_displ,   dW2_displ, 
                              dW3_displ,   dW4_displ,
                              1, batchCount, queue);

        magma_tally3_cdisplace_pointers(dW1_displ, dA_array, ldda, p2, 0, batchCount, queue); 
        magma_tally3_cdisplace_pointers(dA_displ, dA_array, ldda, p2, p2, batchCount, queue); 
        magma_tally3blas_cgemm_batched( Magma_tally3NoTrans, Magma_tally3NoTrans, m2, n2, n1, 
                              MAGMA_tally3_C_NEG_ONE, dW1_displ, ldda, 
                              dW5_displ, ldda, 
                              MAGMA_tally3_C_ONE,  dA_displ, ldda, 
                              batchCount, queue);
        // panel on A2
        //printf("calling recursive panel on A2 with m=%d nb=%d min_recpnb %d\n",m2,n2,min_recpnb);
        arginfo = magma_tally3_cgetrf_recpanel_nopiv_batched(
                           m2, n2, min_recpnb,
                           dA_displ, ldda,
                           dX_array, dX_length,
                           dinvA_array, dinvA_length,
                           dW1_displ, dW2_displ,
                           dW3_displ, dW4_displ, dW5_displ,
                           info_array, gbstep+p2, batchCount, myhandle, queue);
        if (arginfo != 0) return arginfo;
    }

    magma_tally3_free(dA_displ);
    return 0;
}


////////////////////////////////////////////////////////////////////////////////////////

