/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Azzam Haidar
       @generated from zlarft_batched.cpp normal z -> c, Fri Jan 30 19:00:19 2015
*/
#include "common_magma_tally3.h"
#define  max_shared_bsiz 32


static    magma_tally3FloatComplex one   = MAGMA_tally3_C_ONE;
static    magma_tally3FloatComplex zero  = MAGMA_tally3_C_ZERO;

//===================================================================================================================
//===================================================================================================================
//===================================================================================================================
extern "C" void
magma_tally3_clarft_sm32x32_batched(magma_tally3_int_t n, magma_tally3_int_t k, magma_tally3FloatComplex **v_array, magma_tally3_int_t ldv,
                    magma_tally3FloatComplex **tau_array, magma_tally3FloatComplex **T_array, magma_tally3_int_t ldt, magma_tally3_int_t batchCount, cublasHandle_t myhandle, magma_tally3_queue_t queue)
{

    if( k <= 0) return;

     //==================================
     //          GEMV
     //==================================
#define USE_GEMV2
#define use_gemm_larft_sm32

#if defined(use_gemm_larft_sm32)
    //magma_tally3blas_cgemm_batched( Magma_tally3ConjTrans, Magma_tally3NoTrans, k, k, n, MAGMA_tally3_C_ONE, v_array, ldv, v_array, ldv, MAGMA_tally3_C_ZERO, T_array, ldt, batchCount, queue);
    cublasCgemmBatched(myhandle, CUBLAS_OP_C, CUBLAS_OP_N, k, k, n,
                             &one, (const magma_tally3FloatComplex**) v_array, ldv,
                                    (const magma_tally3FloatComplex**) v_array, ldv,
                             &zero,  T_array, ldt, batchCount);

    magma_tally3blas_claset_batched(Magma_tally3Lower, k, k, MAGMA_tally3_C_ZERO, MAGMA_tally3_C_ZERO, T_array, ldt, batchCount, queue);
#else
    #if 1
    for(int i=0; i<k; i++)
    {
        //W(1:i-1) := - tau(i) * V(i:n,1:i-1)' * V(i:n,i)
        //T( i, i ) = tau( i ) 
        //custom implementation.
        #ifdef USE_GEMV2
        magma_tally3blas_clarft_gemvrowwise_batched( n-i, i, 
                            tau_array,
                            v_array, ldv, 
                            T_array, ldt,
                            batchCount, queue);
                            
        #else       
        magma_tally3blas_clarft_gemvcolwise_batched( n-i, i, v_array, ldv, T_array, ldt, tau_array, batchCount, queue);
        #endif
    }
    #else
        //seems to be very slow when k=32 while the one by one loop above is faster
        clarft_gemv_loop_inside_kernel_batched(n, k, tau_array, v_array, ldv, T_array, ldt, batchCount, queue); 
    #endif
#endif
     //==================================
     //          TRMV
     //==================================
     //T(1:i-1,i) := T(1:i-1,1:i-1) * W(1:i-1) i=[1:k]
     magma_tally3blas_clarft_ctrmv_sm32x32_batched(k, k, tau_array, T_array, ldt, T_array, ldt, batchCount, queue);
}





//===================================================================================================================
//===================================================================================================================
//===================================================================================================================
#define RFT_MAG_GEM
extern "C" magma_tally3_int_t
magma_tally3_clarft_batched(magma_tally3_int_t n, magma_tally3_int_t k, magma_tally3_int_t stair_T, 
                magma_tally3FloatComplex **v_array, magma_tally3_int_t ldv,
                magma_tally3FloatComplex **tau_array, magma_tally3FloatComplex **T_array, magma_tally3_int_t ldt, 
                magma_tally3FloatComplex **work_array, magma_tally3_int_t lwork, magma_tally3_int_t batchCount, cublasHandle_t myhandle, 
                magma_tally3_queue_t queue)
{
    if( k <= 0) return 0;
    if( stair_T > 0 && k <= stair_T) return 0;

    magma_tally3_int_t maxnb = max_shared_bsiz;

    if( lwork < k*ldt) 
    {
        magma_tally3_xerbla( __func__, -(10) );
        return -10;
    }

    if( stair_T > 0 && stair_T > maxnb)
    { 
        magma_tally3_xerbla( __func__, -(3) );
        return -3;
    }
    magma_tally3_int_t DEBUG=0;
    magma_tally3_int_t nb = stair_T == 0 ? min(k,maxnb) : stair_T;

    magma_tally3_int_t i, j, prev_n, mycol, rows;

    magma_tally3FloatComplex **dW1_displ  = NULL;
    magma_tally3FloatComplex **dW2_displ  = NULL;
    magma_tally3FloatComplex **dW3_displ  = NULL;
    magma_tally3FloatComplex **dTstep_array  = NULL;

    magma_tally3_malloc((void**)&dW1_displ,  batchCount * sizeof(*dW1_displ));
    magma_tally3_malloc((void**)&dW2_displ,  batchCount * sizeof(*dW2_displ));
    magma_tally3_malloc((void**)&dW3_displ,  batchCount * sizeof(*dW3_displ));
    magma_tally3_malloc((void**)&dTstep_array,  batchCount * sizeof(*dTstep_array));

    //magma_tally3FloatComplex *Tstep =  k > nb ? work : T;
    if(k > nb)
    {
        magma_tally3_cdisplace_pointers(dTstep_array, work_array, lwork, 0, 0, batchCount, queue);
    }
    else
    {
        magma_tally3_cdisplace_pointers(dTstep_array, T_array, ldt, 0, 0, batchCount, queue);
    }

    //magma_tally3_int_t ldtstep = k > nb ? k : ldt;
    magma_tally3_int_t ldtstep = ldt; //a enlever
    // stair_T = 0 meaning all T
    // stair_T > 0 meaning the triangular portion of T has been computed. 
    //                    the value of stair_T is the nb of these triangulars
   

    //GEMV compute the whole triangular upper portion of T (phase 1)
    // TODO addcublas to check perf

#ifdef RFT_MAG_GEM
    magma_tally3blas_cgemm_batched( Magma_tally3ConjTrans, Magma_tally3NoTrans, 
            k, k, n, 
            one,  v_array, ldv, 
                  v_array, ldv, 
            zero, dTstep_array, ldtstep, 
            batchCount, queue);
#else
    cublasCgemmBatched(myhandle, CUBLAS_OP_C, CUBLAS_OP_N, k, k, n,
                             &one, (const magma_tally3FloatComplex**) v_array, ldv,
                                    (const magma_tally3FloatComplex**) v_array, ldv,
                             &zero, dTstep_array, ldtstep, batchCount);
#endif

    magma_tally3blas_claset_batched(Magma_tally3Lower, k, k, MAGMA_tally3_C_ZERO, MAGMA_tally3_C_ZERO, dTstep_array, ldtstep, batchCount, queue);
    // no need for it as T is expected to be lower zero
    //if(k > nb) magma_tally3blas_claset_batched(Magma_tally3Lower, k, k, MAGMA_tally3_C_ZERO, MAGMA_tally3_C_ZERO, dTstep_array, ldtstep, batchCount);
    

    //TRMV
    //T(1:i-1,i) := T(1:i-1,1:i-1) * W(1:i-1) i=[1:k]
    // TRMV is split over block of column of size nb 
    // the update should be done from top to bottom so:
    // 1- a gemm using the previous computed columns
    //    of T to update rectangular upper protion above 
    //    the triangle of my columns 
    // 2- the columns need to be updated by a serial 
    //    loop over of gemv over itself. since we limit the
    //    shared memory to nb, this nb column 
    //    are split vertically by chunk of nb rows

    dim3 grid(1, 1, batchCount);

    for(j=0; j<k; j+=nb)
    {
        prev_n =  j;
        mycol  =  min(nb, k-j);
        // note that myrow = prev_n + mycol;
        if(prev_n>0 && mycol>0){

            if(DEBUG==3) printf("doing gemm on the rectangular portion of size %d %d of T(%d,%d)\n",prev_n,mycol,0,j);

            magma_tally3_cdisplace_pointers(dW1_displ, dTstep_array, ldtstep, 0, j, batchCount, queue);
            magma_tally3_cdisplace_pointers(dW2_displ, T_array,     ldt, 0, j, batchCount, queue);
#ifdef RFT_MAG_GEM
            magma_tally3blas_cgemm_batched( Magma_tally3NoTrans, Magma_tally3NoTrans, 
                    prev_n, mycol, prev_n, 
                    one,  T_array, ldt, 
                          dW1_displ, ldtstep, 
                    zero, dW2_displ, ldt, 
                    batchCount, queue );
#else
            cublasCgemmBatched(myhandle, CUBLAS_OP_N, CUBLAS_OP_N, 
                    prev_n, mycol, prev_n,
                    &one, (const magma_tally3FloatComplex**) T_array, ldt,
                          (const magma_tally3FloatComplex**) dW1_displ, ldtstep,
                    &zero, dW2_displ, ldt, batchCount);
#endif

            // update my rectangular portion (prev_n,mycol) using sequence of gemv 
            magma_tally3_cdisplace_pointers(dW1_displ, dTstep_array, ldtstep, j, j, batchCount, queue);
            magma_tally3_cdisplace_pointers(dW3_displ, tau_array,  1, j, 0, batchCount, queue);

            for(i=0; i<prev_n; i+=nb)
            {
                rows = min(nb,prev_n-i);
                if(DEBUG==3) printf("        doing recctrmv on the rectangular portion of size %d %d of T(%d,%d)\n",rows,mycol,i,j);

                if(rows>0 && mycol>0)
                {
                    magma_tally3_cdisplace_pointers(dW2_displ, T_array,     ldt, i, j, batchCount, queue);
                    magma_tally3blas_clarft_recctrmv_sm32x32_batched(rows, mycol, dW3_displ, dW2_displ, ldt, dW1_displ, ldtstep, batchCount, queue);
                }
            }
        }

        // the upper rectangular protion is updated, now if needed update the triangular portion
        if(stair_T == 0){
            if(DEBUG==3) printf("doing ctrmv on the triangular portion of size %d %d of T(%d,%d)\n",mycol,mycol,j,j);

            if(mycol>0)
            {
                magma_tally3_cdisplace_pointers(dW1_displ, dTstep_array, ldtstep, j, j, batchCount, queue);
                magma_tally3_cdisplace_pointers(dW3_displ, tau_array,  1, j, 0, batchCount, queue);
                magma_tally3_cdisplace_pointers(dW2_displ, T_array,     ldt, j, j, batchCount, queue);
                magma_tally3blas_clarft_ctrmv_sm32x32_batched(mycol, mycol, dW3_displ, dW1_displ, ldtstep, dW2_displ, ldt, batchCount, queue);

            }
        }
    }// end of j

    magma_tally3_free(dW1_displ);
    magma_tally3_free(dW2_displ);
    magma_tally3_free(dW3_displ);
    magma_tally3_free(dTstep_array);

    return 0;
}


/*===============================================================================================================================*/

