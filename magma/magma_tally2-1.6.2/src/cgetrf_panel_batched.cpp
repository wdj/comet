/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2013
       
       @author Azzam Haidar
       @author Tingxing Dong

       @generated from zgetrf_panel_batched.cpp normal z -> c, Fri Jan 30 19:00:19 2015
*/
#include "common_magma_tally2.h"



////////////////////////////////////////////////////////////////////////////////////////
extern "C" magma_tally2_int_t
magma_tally2_cgetrf_recpanel_batched(
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t min_recpnb,    
    magma_tally2FloatComplex** dA_array,    magma_tally2_int_t ldda,
    magma_tally2_int_t** dipiv_array, magma_tally2_int_t** dpivinfo_array,
    magma_tally2FloatComplex** dX_array,    magma_tally2_int_t dX_length,
    magma_tally2FloatComplex** dinvA_array, magma_tally2_int_t dinvA_length,
    magma_tally2FloatComplex** dW1_displ, magma_tally2FloatComplex** dW2_displ,  
    magma_tally2FloatComplex** dW3_displ, magma_tally2FloatComplex** dW4_displ,
    magma_tally2FloatComplex** dW5_displ,
    magma_tally2_int_t *info_array, magma_tally2_int_t gbstep,  
    magma_tally2_int_t batchCount,  cublasHandle_t myhandle, magma_tally2_queue_t queue)
{

    //magma_tally2_int_t DEBUG = 3;
    // Quick return if possible
    if (m ==0 || n == 0) {
        return 0;
    }


    magma_tally2FloatComplex **dA_displ  = NULL;
    magma_tally2_malloc((void**)&dA_displ,   batchCount * sizeof(*dA_displ));
    magma_tally2_int_t **dipiv_displ = NULL;
    magma_tally2_malloc((void**)&dipiv_displ, batchCount * sizeof(*dipiv_displ));
    
    magma_tally2_int_t panel_nb = n;
    if(panel_nb <= min_recpnb){
        //if(DEBUG>0)printf("calling bottom panel recursive with m=%d nb=%d\n",m,n);
        //  panel factorization
        //magma_tally2_cdisplace_pointers(dA_displ, dA_array, ldda, 0, 0, batchCount);
        magma_tally2_cgetf2_batched(m, panel_nb,
                           dA_array, ldda,
                           dW1_displ, dW2_displ, dW3_displ,
                           dipiv_array, info_array, gbstep, batchCount, myhandle, queue);
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
        //if(DEBUG>0)printf("calling recursive panel on A1 with m=%d nb=%d min_recpnb %d\n",m1,n1,min_recpnb);
        magma_tally2_cdisplace_pointers(dA_displ, dA_array, ldda, p1, p1, batchCount, queue); 
        magma_tally2_idisplace_pointers(dipiv_displ, dipiv_array, 1, p1, 0, batchCount, queue);
        magma_tally2_cgetrf_recpanel_batched(
                           m1, n1, min_recpnb,
                           dA_displ, ldda,
                           dipiv_displ, dpivinfo_array,
                           dX_array, dX_length,
                           dinvA_array, dinvA_length,
                           dW1_displ, dW2_displ,
                           dW3_displ, dW4_displ, dW5_displ,
                           info_array, gbstep, batchCount, myhandle, queue);

        // update A2
        //if(DEBUG>0)printf("calling TRSM  with             m=%d n=%d \n",m1,n2);
        
        // setup pivinfo 
        setup_pivinfo_batched(dpivinfo_array, dipiv_displ, m1, n1, batchCount, queue);
        magma_tally2_cdisplace_pointers(dW5_displ, dA_array, ldda, p1, p2, batchCount, queue); 
        magma_tally2_claswp_rowparallel_batched( n2, dW5_displ, ldda,
                           dX_array, n1,
                           0, n1,
                           dpivinfo_array, batchCount, queue);
        magma_tally2blas_ctrsm_outofplace_batched(Magma_tally2Left, Magma_tally2Lower, Magma_tally2NoTrans, Magma_tally2Unit, 1,
                              n1, n2,
                              MAGMA_tally2_C_ONE,
                              dA_displ,    ldda, // dA
                              dX_array,  n1, // dB
                              dW5_displ,   ldda, // dX
                              dinvA_array, dinvA_length,
                              dW1_displ,   dW2_displ, 
                              dW3_displ,   dW4_displ,
                              0, batchCount, queue);

        magma_tally2_cdisplace_pointers(dW1_displ, dA_array, ldda, p2, 0, batchCount, queue); 
        magma_tally2_cdisplace_pointers(dA_displ, dA_array, ldda, p2, p2, batchCount, queue); 

        //if(DEBUG>0)printf("calling update A2(%d,%d) -= A(%d,%d)*A(%d,%d)  with             m=%d n=%d k=%d ldda %d\n",p2,p2,p2,0,p1,p2,m2,n2,n1,ldda);

#if 0
        magma_tally2FloatComplex neg_one = MAGMA_tally2_C_NEG_ONE;
        magma_tally2FloatComplex one  = MAGMA_tally2_C_ONE;
        cublasCgemmBatched(myhandle, CUBLAS_OP_N, CUBLAS_OP_N, m2, n2, n1,
                                     &neg_one, (const magma_tally2FloatComplex**) dW1_displ, ldda,
                                               (const magma_tally2FloatComplex**) dW5_displ, ldda,
                                     &one,  dA_displ, ldda, batchCount );


#else

        magma_tally2blas_cgemm_batched( Magma_tally2NoTrans, Magma_tally2NoTrans, m2, n2, n1, 
                              MAGMA_tally2_C_NEG_ONE, dW1_displ, ldda, 
                              dW5_displ, ldda, 
                              MAGMA_tally2_C_ONE,  dA_displ, ldda, 
                              batchCount, queue);
#endif

        // panel on A2
        //if(DEBUG>0)printf("calling recursive panel on A2 with m=%d nb=%d min_recpnb %d\n",m2,n2,min_recpnb);
        magma_tally2_idisplace_pointers(dipiv_displ, dipiv_array, 1, p2, 0, batchCount, queue);
        magma_tally2_cgetrf_recpanel_batched(
                           m2, n2, min_recpnb,
                           dA_displ, ldda,
                           dipiv_displ, dpivinfo_array,
                           dX_array, dX_length,
                           dinvA_array, dinvA_length,
                           dW1_displ, dW2_displ,
                           dW3_displ, dW4_displ, dW5_displ,
                           info_array, gbstep+p2, batchCount, myhandle, queue);

        // setup pivinfo
        setup_pivinfo_batched(dpivinfo_array, dipiv_displ, m2, n2, batchCount, queue);
        adjust_ipiv_batched(dipiv_displ, n2, n1, batchCount, queue);
        
        magma_tally2_cdisplace_pointers(dW1_displ, dA_array, ldda, p2, 0, batchCount, queue); // no need since it is above
        magma_tally2_claswp_rowparallel_batched( n1, dW1_displ, ldda,
                           dW1_displ, ldda,
                           n1, n,
                           dpivinfo_array, batchCount, queue);

        
    }

    magma_tally2_free(dA_displ);
    magma_tally2_free(dipiv_displ);
    return 0;
}


////////////////////////////////////////////////////////////////////////////////////////

