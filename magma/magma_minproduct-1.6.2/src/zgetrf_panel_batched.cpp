/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2013
       
       @author Azzam Haidar
       @author Tingxing Dong

       @precisions normal z -> s d c
*/
#include "common_magma_minproduct.h"



////////////////////////////////////////////////////////////////////////////////////////
extern "C" magma_minproduct_int_t
magma_minproduct_zgetrf_recpanel_batched(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t min_recpnb,    
    magma_minproductDoubleComplex** dA_array,    magma_minproduct_int_t ldda,
    magma_minproduct_int_t** dipiv_array, magma_minproduct_int_t** dpivinfo_array,
    magma_minproductDoubleComplex** dX_array,    magma_minproduct_int_t dX_length,
    magma_minproductDoubleComplex** dinvA_array, magma_minproduct_int_t dinvA_length,
    magma_minproductDoubleComplex** dW1_displ, magma_minproductDoubleComplex** dW2_displ,  
    magma_minproductDoubleComplex** dW3_displ, magma_minproductDoubleComplex** dW4_displ,
    magma_minproductDoubleComplex** dW5_displ,
    magma_minproduct_int_t *info_array, magma_minproduct_int_t gbstep,  
    magma_minproduct_int_t batchCount,  cublasHandle_t myhandle, magma_minproduct_queue_t queue)
{

    //magma_minproduct_int_t DEBUG = 3;
    // Quick return if possible
    if (m ==0 || n == 0) {
        return 0;
    }


    magma_minproductDoubleComplex **dA_displ  = NULL;
    magma_minproduct_malloc((void**)&dA_displ,   batchCount * sizeof(*dA_displ));
    magma_minproduct_int_t **dipiv_displ = NULL;
    magma_minproduct_malloc((void**)&dipiv_displ, batchCount * sizeof(*dipiv_displ));
    
    magma_minproduct_int_t panel_nb = n;
    if(panel_nb <= min_recpnb){
        //if(DEBUG>0)printf("calling bottom panel recursive with m=%d nb=%d\n",m,n);
        //  panel factorization
        //magma_minproduct_zdisplace_pointers(dA_displ, dA_array, ldda, 0, 0, batchCount);
        magma_minproduct_zgetf2_batched(m, panel_nb,
                           dA_array, ldda,
                           dW1_displ, dW2_displ, dW3_displ,
                           dipiv_array, info_array, gbstep, batchCount, myhandle, queue);
    }
    else{
        // split A over two [A A2]
        // panel on A1, update on A2 then panel on A1    
        magma_minproduct_int_t n1 = n/2;
        magma_minproduct_int_t n2 = n-n1;
        magma_minproduct_int_t m1 = m;
        magma_minproduct_int_t m2 = m-n1;
        magma_minproduct_int_t p1 = 0;
        magma_minproduct_int_t p2 = n1;
        // panel on A1
        //if(DEBUG>0)printf("calling recursive panel on A1 with m=%d nb=%d min_recpnb %d\n",m1,n1,min_recpnb);
        magma_minproduct_zdisplace_pointers(dA_displ, dA_array, ldda, p1, p1, batchCount, queue); 
        magma_minproduct_idisplace_pointers(dipiv_displ, dipiv_array, 1, p1, 0, batchCount, queue);
        magma_minproduct_zgetrf_recpanel_batched(
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
        magma_minproduct_zdisplace_pointers(dW5_displ, dA_array, ldda, p1, p2, batchCount, queue); 
        magma_minproduct_zlaswp_rowparallel_batched( n2, dW5_displ, ldda,
                           dX_array, n1,
                           0, n1,
                           dpivinfo_array, batchCount, queue);
        magma_minproductblas_ztrsm_outofplace_batched(Magma_minproductLeft, Magma_minproductLower, Magma_minproductNoTrans, Magma_minproductUnit, 1,
                              n1, n2,
                              MAGMA_minproduct_Z_ONE,
                              dA_displ,    ldda, // dA
                              dX_array,  n1, // dB
                              dW5_displ,   ldda, // dX
                              dinvA_array, dinvA_length,
                              dW1_displ,   dW2_displ, 
                              dW3_displ,   dW4_displ,
                              0, batchCount, queue);

        magma_minproduct_zdisplace_pointers(dW1_displ, dA_array, ldda, p2, 0, batchCount, queue); 
        magma_minproduct_zdisplace_pointers(dA_displ, dA_array, ldda, p2, p2, batchCount, queue); 

        //if(DEBUG>0)printf("calling update A2(%d,%d) -= A(%d,%d)*A(%d,%d)  with             m=%d n=%d k=%d ldda %d\n",p2,p2,p2,0,p1,p2,m2,n2,n1,ldda);

#if 0
        magma_minproductDoubleComplex neg_one = MAGMA_minproduct_Z_NEG_ONE;
        magma_minproductDoubleComplex one  = MAGMA_minproduct_Z_ONE;
        cublasZgemmBatched(myhandle, CUBLAS_OP_N, CUBLAS_OP_N, m2, n2, n1,
                                     &neg_one, (const magma_minproductDoubleComplex**) dW1_displ, ldda,
                                               (const magma_minproductDoubleComplex**) dW5_displ, ldda,
                                     &one,  dA_displ, ldda, batchCount );


#else

        magma_minproductblas_zgemm_batched( Magma_minproductNoTrans, Magma_minproductNoTrans, m2, n2, n1, 
                              MAGMA_minproduct_Z_NEG_ONE, dW1_displ, ldda, 
                              dW5_displ, ldda, 
                              MAGMA_minproduct_Z_ONE,  dA_displ, ldda, 
                              batchCount, queue);
#endif

        // panel on A2
        //if(DEBUG>0)printf("calling recursive panel on A2 with m=%d nb=%d min_recpnb %d\n",m2,n2,min_recpnb);
        magma_minproduct_idisplace_pointers(dipiv_displ, dipiv_array, 1, p2, 0, batchCount, queue);
        magma_minproduct_zgetrf_recpanel_batched(
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
        
        magma_minproduct_zdisplace_pointers(dW1_displ, dA_array, ldda, p2, 0, batchCount, queue); // no need since it is above
        magma_minproduct_zlaswp_rowparallel_batched( n1, dW1_displ, ldda,
                           dW1_displ, ldda,
                           n1, n,
                           dpivinfo_array, batchCount, queue);

        
    }

    magma_minproduct_free(dA_displ);
    magma_minproduct_free(dipiv_displ);
    return 0;
}


////////////////////////////////////////////////////////////////////////////////////////

