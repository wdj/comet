/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2013
       
   @author Azzam Haidar
   @author Adrien Remy

   @generated from zgetrf_panel_nopiv_batched.cpp normal z -> d, Fri Jan 30 19:00:19 2015
*/
#include "common_magma_tally2.h"
////////////////////////////////////////////////////////////////////////////////////////
extern "C" magma_tally2_int_t
magma_tally2_dgetrf_panel_nopiv_batched(
    magma_tally2_int_t m, magma_tally2_int_t nb,    
    double** dA_array,    magma_tally2_int_t ldda,
    double** dX_array,    magma_tally2_int_t dX_length,
    double** dinvA_array, magma_tally2_int_t dinvA_length,
    double** dW0_displ, double** dW1_displ,  
    double** dW2_displ, double** dW3_displ,
    double** dW4_displ,     
    magma_tally2_int_t *info_array, magma_tally2_int_t gbstep,  
    magma_tally2_int_t batchCount, cublasHandle_t myhandle, magma_tally2_queue_t queue )
{
     magma_tally2_int_t arginfo = 0;
    //===============================================
    //  panel factorization
    //===============================================
    if(m < nb){
        printf("magma_tally2_dgetrf_panel_nopiv_batched_q m < nb %d < %d \n",(int) m, (int) nb);
        return -101;
    }

#if 0
    arginfo = magma_tally2_dgetf2_nopiv_batched(
                       m, nb,
                       dA_array, ldda,
                       dW1_displ, dW2_displ, dW3_displ,
                       info_array, gbstep, batchCount, myhandle);
    if (arginfo != 0) return arginfo;
#else
    arginfo = magma_tally2_dgetf2_nopiv_batched(
                       nb, nb,
                       dA_array, ldda,
                       dW1_displ, dW2_displ, dW3_displ,
                       info_array, gbstep, batchCount, myhandle, queue);
    if (arginfo != 0) return arginfo;
    if((m-nb) > 0){
        magma_tally2_ddisplace_pointers(dW0_displ, dA_array, ldda, nb, 0, batchCount, queue);
        magma_tally2blas_dtrsm_work_batched(Magma_tally2Right, Magma_tally2Upper, Magma_tally2NoTrans, Magma_tally2NonUnit,
                              1, m-nb, nb, 
                              MAGMA_tally2_D_ONE,
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
extern "C" magma_tally2_int_t
magma_tally2_dgetrf_recpanel_nopiv_batched(
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t min_recpnb,    
    double** dA_array,    magma_tally2_int_t ldda,
    double** dX_array,    magma_tally2_int_t dX_length,
    double** dinvA_array, magma_tally2_int_t dinvA_length,
    double** dW1_displ, double** dW2_displ,  
    double** dW3_displ, double** dW4_displ,
    double** dW5_displ, 
    magma_tally2_int_t *info_array, magma_tally2_int_t gbstep,
    magma_tally2_int_t batchCount, cublasHandle_t myhandle, magma_tally2_queue_t queue )
{
    // Quick return if possible
    if (m == 0 || n == 0) {
        return 0;
    }
    magma_tally2_int_t arginfo = 0;


    double **dA_displ  = NULL;
    magma_tally2_malloc((void**)&dA_displ,   batchCount * sizeof(*dA_displ));
    
    magma_tally2_int_t panel_nb = n;
    if(panel_nb <= min_recpnb){
        // if(DEBUG>0)printf("calling bottom panel recursive with m=%d nb=%d\n",m,n);
        //  panel factorization
        //magma_tally2_ddisplace_pointers(dA_displ, dA_array, ldda, 0, 0, batchCount);
        arginfo = magma_tally2_dgetrf_panel_nopiv_batched(
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
        magma_tally2_int_t n1 = n/2;
        magma_tally2_int_t n2 = n-n1;
        magma_tally2_int_t m1 = m;
        magma_tally2_int_t m2 = m-n1;
        magma_tally2_int_t p1 = 0;
        magma_tally2_int_t p2 = n1;
        // panel on A1
        //printf("calling recursive panel on A1 with m=%d nb=%d min_recpnb %d\n",m1,n1,min_recpnb);
        magma_tally2_ddisplace_pointers(dA_displ, dA_array, ldda, p1, p1, batchCount, queue); 
        arginfo = magma_tally2_dgetrf_recpanel_nopiv_batched(
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
        
        magma_tally2_ddisplace_pointers(dW5_displ, dA_array, ldda, p1, p2, batchCount, queue); 
        magma_tally2blas_dtrsm_work_batched(Magma_tally2Left, Magma_tally2Lower, Magma_tally2NoTrans, Magma_tally2Unit, 1,
                              n1, n2,
                              MAGMA_tally2_D_ONE,
                              dA_displ,    ldda, // dA
                              dW5_displ,   ldda, // dB
                              dX_array,  n1, // dX
                              dinvA_array, dinvA_length,
                              dW1_displ,   dW2_displ, 
                              dW3_displ,   dW4_displ,
                              1, batchCount, queue);

        magma_tally2_ddisplace_pointers(dW1_displ, dA_array, ldda, p2, 0, batchCount, queue); 
        magma_tally2_ddisplace_pointers(dA_displ, dA_array, ldda, p2, p2, batchCount, queue); 
        magma_tally2blas_dgemm_batched( Magma_tally2NoTrans, Magma_tally2NoTrans, m2, n2, n1, 
                              MAGMA_tally2_D_NEG_ONE, dW1_displ, ldda, 
                              dW5_displ, ldda, 
                              MAGMA_tally2_D_ONE,  dA_displ, ldda, 
                              batchCount, queue);
        // panel on A2
        //printf("calling recursive panel on A2 with m=%d nb=%d min_recpnb %d\n",m2,n2,min_recpnb);
        arginfo = magma_tally2_dgetrf_recpanel_nopiv_batched(
                           m2, n2, min_recpnb,
                           dA_displ, ldda,
                           dX_array, dX_length,
                           dinvA_array, dinvA_length,
                           dW1_displ, dW2_displ,
                           dW3_displ, dW4_displ, dW5_displ,
                           info_array, gbstep+p2, batchCount, myhandle, queue);
        if (arginfo != 0) return arginfo;
    }

    magma_tally2_free(dA_displ);
    return 0;
}


////////////////////////////////////////////////////////////////////////////////////////

