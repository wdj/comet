/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2013
       
       @author Azzam Haidar
       @author Tingxing Dong

       @generated from zpotf2_batched.cpp normal z -> c, Fri Jan 30 19:00:19 2015
*/
#include "common_magma_tally4.h"
#include "batched_kernel_param.h"

#define PRECISION_c
/////////////////////////////////////////////////////////////////
extern "C" magma_tally4_int_t
magma_tally4_cpotf2_ctrsm_batched(
    magma_tally4_uplo_t uplo, magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4FloatComplex **dA_array, magma_tally4_int_t lda,
    magma_tally4FloatComplex **dA_displ, 
    magma_tally4FloatComplex **dB_displ, 
    magma_tally4FloatComplex **dC_displ,
    magma_tally4_int_t *info_array, magma_tally4_int_t gbstep,  
    magma_tally4_int_t batchCount, magma_tally4_queue_t queue)
{


    magma_tally4_int_t j;
    magma_tally4_int_t arginfo = 0;
    if( m > MAX_NTHREADS )
    {
        printf("magma_tally4_cpotf2_ctrsm_batched m=%d > %d not supported today \n", (int) m, (int) MAX_NTHREADS);
        arginfo =-13;
        return arginfo;
    }

    // Quick return if possible
    if (n == 0) {
        return arginfo;
    }


    magma_tally4FloatComplex alpha = MAGMA_tally4_C_NEG_ONE;
    magma_tally4FloatComplex beta  = MAGMA_tally4_C_ONE;

    if (uplo == Magma_tally4Upper) {
        printf("Upper side is unavailable \n");
    }
    else {
        for(j = 0; j < n; j++) {
            magma_tally4_cpotf2_cdotc_batched(j, dA_array, lda, j, info_array, gbstep, batchCount, queue); // including cdotc product and update a(j,j)
            if (j < n) {
                #if defined(PRECISION_z) || defined(PRECISION_c)
                magma_tally4_clacgv_batched(j, dA_array, lda, j, batchCount, queue);
                #endif

                magma_tally4_cdisplace_pointers(dA_displ, dA_array, lda, j+1, 0, batchCount, queue);
                magma_tally4_cdisplace_pointers(dB_displ, dA_array, lda, j, 0, batchCount, queue);
                magma_tally4_cdisplace_pointers(dC_displ, dA_array, lda, j+1, j, batchCount, queue);

                // Compute elements J+1:N of column J = A(j+1:n,1:j-1) * A(j,1:j-1) (row).
                magma_tally4blas_cgemv_batched(Magma_tally4NoTrans, m-j-1, j,
                                 alpha, dA_displ, lda,
                                        dB_displ,    lda,
                                 beta,  dC_displ, 1,
                                 batchCount, queue);// 


                #if defined(PRECISION_z) || defined(PRECISION_c)
                magma_tally4_clacgv_batched(j, dA_array, lda, j, batchCount, queue);
                #endif
                magma_tally4_cpotf2_csscal_batched(m-j, dA_array, 1, j+j*lda, info_array, batchCount, queue);
            }
        }
    }

    return arginfo;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
extern "C" magma_tally4_int_t
magma_tally4_cpotf2_batched(
    magma_tally4_uplo_t uplo, magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4FloatComplex **dA_array, magma_tally4_int_t lda,
    magma_tally4FloatComplex **dA_displ, 
    magma_tally4FloatComplex **dW_displ,
    magma_tally4FloatComplex **dB_displ, 
    magma_tally4FloatComplex **dC_displ, 
    magma_tally4_int_t *info_array, magma_tally4_int_t gbstep, 
    magma_tally4_int_t batchCount, cublasHandle_t myhandle, magma_tally4_queue_t queue)
{

    magma_tally4_int_t j;

    // Quick return if possible
    if (n == 0) {
        return 1;
    }

    magma_tally4FloatComplex alpha = MAGMA_tally4_C_NEG_ONE;
    magma_tally4FloatComplex beta  = MAGMA_tally4_C_ONE;


    int nb = POTF2_NB;
    int ib, rows;

    if (uplo == Magma_tally4Upper) {
       printf("Upper side is unavailable \n");
    }
    else {
        for(j = 0; j < n; j+= nb) {
            ib   = min(nb, n-j);
            rows = m-j;
            if( (rows <= POTF2_TILE_SIZE) && (ib <= POTF2_TILE_SIZE) ){
                magma_tally4_cdisplace_pointers(dA_displ, dA_array, lda, j, j, batchCount, queue);
                magma_tally4_cpotf2_tile_batched(
                               uplo, rows, ib,
                               dA_displ, lda,
                               info_array, gbstep, batchCount, queue);
            }
            else{
                 magma_tally4_cdisplace_pointers(dA_displ, dA_array, lda, j, j, batchCount, queue); 
                 magma_tally4_cpotf2_ctrsm_batched(
                           uplo, rows, ib,
                           dA_displ, lda,
                           dW_displ, dB_displ, dC_displ, 
                           info_array, gbstep, batchCount, queue);

            }
#if 1

//#define RIGHT_LOOKING
            if( (n-j-ib) > 0){
#ifdef RIGHT_LOOKING
                magma_tally4_cdisplace_pointers(dA_displ, dA_array, lda, j+ib, j, batchCount, queue);
                magma_tally4_cdisplace_pointers(dC_displ, dA_array, lda, j+ib, j+ib, batchCount, queue);
                #if 1
                magma_tally4blas_cgemm_batched(Magma_tally4NoTrans, Magma_tally4ConjTrans, m-j-ib, n-j-ib, ib,
                             alpha, dA_displ, lda,
                                    dA_displ, lda,
                             beta,  dC_displ, lda, batchCount, queue );
                #else
                cublasCgemmBatched(myhandle, CUBLAS_OP_N, CUBLAS_OP_C, m-j-ib, n-j-ib, ib,
                             &alpha, (const magma_tally4FloatComplex**) dA_displ, lda,
                                    (const magma_tally4FloatComplex**) dA_displ, lda,
                             &beta,  dC_displ, lda, batchCount, queue );
                #endif
#else
                // update next subpanel
                magma_tally4_cdisplace_pointers(dA_displ, dA_array, lda, j+ib, 0, batchCount, queue);
                magma_tally4_cdisplace_pointers(dC_displ, dA_array, lda, j+ib, j+ib, batchCount, queue);
                #if 1
                magma_tally4blas_cgemm_batched(Magma_tally4NoTrans, Magma_tally4ConjTrans, m-j-ib, min((n-j-ib),ib), j+ib,
                             alpha, dA_displ, lda,
                                    dA_displ, lda,
                             beta,  dC_displ, lda, batchCount, queue );
                #else
                cublasCgemmBatched(myhandle, CUBLAS_OP_N, CUBLAS_OP_C, m-j-ib, min((n-j-ib),ib), j+ib,
                             &alpha, (const magma_tally4FloatComplex**) dA_displ, lda,
                                    (const magma_tally4FloatComplex**) dA_displ, lda,
                             &beta,  dC_displ, lda, batchCount );
                #endif
#endif
             } // end of if( (n-j-ib) > 0)
#endif
        }
    }

    return 0;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



