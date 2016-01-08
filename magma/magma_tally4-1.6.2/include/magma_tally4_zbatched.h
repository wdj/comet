/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Azzam Haidar
       @author Tingxing Dong

       @precisions normal z -> s d c
*/

#ifndef MAGMA_tally4_ZBATCHED_H
#define MAGMA_tally4_ZBATCHED_H

#include "magma_tally4_types.h"

#define COMPLEX

#ifdef __cplusplus
extern "C" {
#endif
  /*
   *  local auxiliary routines
   */
void 
zset_pointer_int(magma_tally4_int_t **output_array,
        magma_tally4_int_t *input,
        magma_tally4_int_t lda,
        magma_tally4_int_t row, magma_tally4_int_t column, 
        magma_tally4_int_t batchSize,
        magma_tally4_int_t batchCount, magma_tally4_queue_t queue);

void 
zset_pointer(magma_tally4DoubleComplex **output_array,
                 magma_tally4DoubleComplex *input,
                 magma_tally4_int_t lda,
                 magma_tally4_int_t row, magma_tally4_int_t column,
                 magma_tally4_int_t batchSize,
                 magma_tally4_int_t batchCount, magma_tally4_queue_t queue);


void 
zset_array(magma_tally4DoubleComplex **output_array,
               magma_tally4DoubleComplex **input_array, magma_tally4_int_t lda,
               magma_tally4_int_t row, magma_tally4_int_t column,
               magma_tally4_int_t batchCount, magma_tally4_queue_t queue);

void 
magma_tally4_zdisplace_pointers(magma_tally4DoubleComplex **output_array,
               magma_tally4DoubleComplex **input_array, magma_tally4_int_t lda,
               magma_tally4_int_t row, magma_tally4_int_t column, 
               magma_tally4_int_t batchCount, magma_tally4_queue_t queue);

  /*
   *  LAPACK batched routines
   */

  /*
   *  BLAS batched routines
   */

void
magma_tally4blas_zgemm_batched(
    magma_tally4_trans_t transA, magma_tally4_trans_t transB, magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4DoubleComplex alpha,
    magma_tally4DoubleComplex const * const * dA_array, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex const * const * dB_array, magma_tally4_int_t lddb,
    magma_tally4DoubleComplex beta,
    magma_tally4DoubleComplex **dC_array, magma_tally4_int_t lddc, magma_tally4_int_t batchCount, magma_tally4_queue_t queue );

void
magma_tally4blas_zgemm_batched_lg(
    magma_tally4_trans_t transA, magma_tally4_trans_t transB, magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4DoubleComplex alpha,
    magma_tally4DoubleComplex const * const * dA_array, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex const * const * dB_array, magma_tally4_int_t lddb,
    magma_tally4DoubleComplex beta,
    magma_tally4DoubleComplex **dC_array, magma_tally4_int_t lddc, magma_tally4_int_t batchCount, magma_tally4_queue_t queue );
void
magma_tally4blas_zgemm_batched_k32(
    magma_tally4_trans_t transA, magma_tally4_trans_t transB, magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4DoubleComplex alpha,
    magma_tally4DoubleComplex const * const * dA_array, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex const * const * dB_array, magma_tally4_int_t lddb,
    magma_tally4DoubleComplex beta,
    magma_tally4DoubleComplex **dC_array, magma_tally4_int_t lddc, magma_tally4_int_t batchCount, magma_tally4_queue_t queue );



void 
magma_tally4blas_zherk_NC_batched( magma_tally4_trans_t TRANSA, magma_tally4_trans_t TRANSB, int m , int n , int k , 
                       magma_tally4DoubleComplex alpha, magma_tally4DoubleComplex **dA_array, int lda, 
                       magma_tally4DoubleComplex **B_array, int ldb, 
                       magma_tally4DoubleComplex beta,        magma_tally4DoubleComplex **C_array, int ldc, 
                       magma_tally4_int_t batchCount);

void
magma_tally4blas_zherk_batched(
    magma_tally4_uplo_t uplo, magma_tally4_trans_t trans, magma_tally4_int_t n, magma_tally4_int_t k,
    double alpha,
    magma_tally4DoubleComplex const * const * dA_array, magma_tally4_int_t ldda,
    double beta,
    magma_tally4DoubleComplex **dC_array, magma_tally4_int_t lddc, magma_tally4_int_t batchCount, magma_tally4_queue_t queue );

void
magma_tally4blas_zherk_batched_lg(
    magma_tally4_uplo_t uplo, magma_tally4_trans_t trans, magma_tally4_int_t n, magma_tally4_int_t k,
    double alpha,
    magma_tally4DoubleComplex const * const * dA_array, magma_tally4_int_t ldda,
    double beta,
    magma_tally4DoubleComplex **dC_array, magma_tally4_int_t lddc, magma_tally4_int_t batchCount, magma_tally4_queue_t queue );

void
magma_tally4blas_zherk_batched_k32(
    magma_tally4_uplo_t uplo, magma_tally4_trans_t trans, magma_tally4_int_t n, magma_tally4_int_t k,
    double alpha,
    magma_tally4DoubleComplex const * const * dA_array, magma_tally4_int_t ldda,
    double beta,
    magma_tally4DoubleComplex **dC_array, magma_tally4_int_t lddc, magma_tally4_int_t batchCount, magma_tally4_queue_t queue );






magma_tally4_int_t 
magma_tally4_zpotf2_tile_batched(
    magma_tally4_uplo_t uplo, magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex **dA_array, magma_tally4_int_t lda,
    magma_tally4_int_t *info_array, magma_tally4_int_t gbstep, magma_tally4_int_t batchCount, magma_tally4_queue_t queue);

magma_tally4_int_t 
magma_tally4_zpotrf_panel(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n, magma_tally4_int_t nb,     
    magma_tally4DoubleComplex *A, magma_tally4_int_t lda,
    magma_tally4DoubleComplex* dinvA, magma_tally4DoubleComplex** dinvA_array, magma_tally4_int_t invA_msize,
    magma_tally4DoubleComplex* x, magma_tally4DoubleComplex** x_array,  magma_tally4_int_t x_msize,
    magma_tally4_int_t *info_array,  magma_tally4_int_t batchCount, magma_tally4_int_t matrixSize, magma_tally4_queue_t queue);


void 
magma_tally4blas_ztrtri_diag_batched(
    magma_tally4_uplo_t uplo, magma_tally4_diag_t diag, magma_tally4_int_t n,
    magma_tally4DoubleComplex const * const *dA_array, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex **dinvA_array, 
    magma_tally4_int_t resetozero, magma_tally4_int_t batchCount, magma_tally4_queue_t queue);



void 
magma_tally4blas_ztrsm_batched(
    magma_tally4_side_t side, magma_tally4_uplo_t uplo, magma_tally4_trans_t transA, magma_tally4_diag_t diag,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex alpha,
    magma_tally4DoubleComplex** dA_array,    magma_tally4_int_t ldda,
    magma_tally4DoubleComplex** dB_array,    magma_tally4_int_t lddb,
    magma_tally4_int_t batchCount, magma_tally4_queue_t queue);


void 
magma_tally4blas_ztrsm_work_batched(
    magma_tally4_side_t side, magma_tally4_uplo_t uplo, magma_tally4_trans_t transA, magma_tally4_diag_t diag,
    magma_tally4_int_t flag, magma_tally4_int_t m, magma_tally4_int_t n, 
    magma_tally4DoubleComplex alpha, 
    magma_tally4DoubleComplex** dA_array,    magma_tally4_int_t ldda,
    magma_tally4DoubleComplex** dB_array,    magma_tally4_int_t lddb,
    magma_tally4DoubleComplex** dX_array,    magma_tally4_int_t lddx, 
    magma_tally4DoubleComplex** dinvA_array, magma_tally4_int_t dinvA_length,
    magma_tally4DoubleComplex** dA_displ, magma_tally4DoubleComplex** dB_displ, 
    magma_tally4DoubleComplex** dX_displ, magma_tally4DoubleComplex** dinvA_displ,
    magma_tally4_int_t resetozero, magma_tally4_int_t batchCount, magma_tally4_queue_t queue);

void 
magma_tally4blas_ztrsm_outofplace_batched(
    magma_tally4_side_t side, magma_tally4_uplo_t uplo, magma_tally4_trans_t transA, magma_tally4_diag_t diag,
    magma_tally4_int_t flag, magma_tally4_int_t m, magma_tally4_int_t n, 
    magma_tally4DoubleComplex alpha, 
    magma_tally4DoubleComplex** dA_array,    magma_tally4_int_t ldda,
    magma_tally4DoubleComplex** dB_array,    magma_tally4_int_t lddb,
    magma_tally4DoubleComplex** dX_array,    magma_tally4_int_t lddx, 
    magma_tally4DoubleComplex** dinvA_array, magma_tally4_int_t dinvA_length,
    magma_tally4DoubleComplex** dA_displ, magma_tally4DoubleComplex** dB_displ, 
    magma_tally4DoubleComplex** dX_displ, magma_tally4DoubleComplex** dinvA_displ,
    magma_tally4_int_t resetozero, magma_tally4_int_t batchCount, magma_tally4_queue_t queue);

magma_tally4_int_t 
magma_tally4_zpotrf_batched(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    magma_tally4DoubleComplex **dA_array, magma_tally4_int_t lda,
    magma_tally4_int_t *info_array,  magma_tally4_int_t batchCount, magma_tally4_queue_t queue);


magma_tally4_int_t 
magma_tally4_zpotf2_batched(
    magma_tally4_uplo_t uplo, magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex **dA_array, magma_tally4_int_t lda,
    magma_tally4DoubleComplex **dA_displ, 
    magma_tally4DoubleComplex **dW_displ,
    magma_tally4DoubleComplex **dB_displ, 
    magma_tally4DoubleComplex **dC_displ, 
    magma_tally4_int_t *info_array, magma_tally4_int_t gbstep, 
    magma_tally4_int_t batchCount, cublasHandle_t myhandle, magma_tally4_queue_t queue);

magma_tally4_int_t 
magma_tally4_zpotrf_panel_batched(
    magma_tally4_uplo_t uplo, magma_tally4_int_t n, magma_tally4_int_t nb,     
    magma_tally4DoubleComplex** dA_array,    magma_tally4_int_t ldda,
    magma_tally4DoubleComplex** dX_array,    magma_tally4_int_t dX_length,
    magma_tally4DoubleComplex** dinvA_array, magma_tally4_int_t dinvA_length,
    magma_tally4DoubleComplex** dW0_displ, magma_tally4DoubleComplex** dW1_displ, 
    magma_tally4DoubleComplex** dW2_displ, magma_tally4DoubleComplex** dW3_displ,
    magma_tally4DoubleComplex** dW4_displ, 
    magma_tally4_int_t *info_array, magma_tally4_int_t gbstep,
    magma_tally4_int_t batchCount, cublasHandle_t myhandle, magma_tally4_queue_t queue);

magma_tally4_int_t 
magma_tally4_zpotrf_recpanel_batched(
    magma_tally4_uplo_t uplo, magma_tally4_int_t m, magma_tally4_int_t n, 
    magma_tally4_int_t min_recpnb,    
    magma_tally4DoubleComplex** dA_array,    magma_tally4_int_t ldda,
    magma_tally4DoubleComplex** dX_array,    magma_tally4_int_t dX_length,
    magma_tally4DoubleComplex** dinvA_array, magma_tally4_int_t dinvA_length,
    magma_tally4DoubleComplex** dW0_displ, magma_tally4DoubleComplex** dW1_displ,  
    magma_tally4DoubleComplex** dW2_displ, magma_tally4DoubleComplex** dW3_displ,
    magma_tally4DoubleComplex** dW4_displ,
    magma_tally4_int_t *info_array, magma_tally4_int_t gbstep, 
    magma_tally4_int_t batchCount, cublasHandle_t myhandle, magma_tally4_queue_t queue);

magma_tally4_int_t 
magma_tally4_zpotrf_rectile_batched(
    magma_tally4_uplo_t uplo, magma_tally4_int_t m, magma_tally4_int_t n, 
    magma_tally4_int_t min_recpnb,    
    magma_tally4DoubleComplex** dA_array,    magma_tally4_int_t ldda,
    magma_tally4DoubleComplex** dX_array,    magma_tally4_int_t dX_length,
    magma_tally4DoubleComplex** dinvA_array, magma_tally4_int_t dinvA_length,
    magma_tally4DoubleComplex** dW0_displ, magma_tally4DoubleComplex** dW1_displ,  
    magma_tally4DoubleComplex** dW2_displ, magma_tally4DoubleComplex** dW3_displ,
    magma_tally4DoubleComplex** dW4_displ,
    magma_tally4_int_t *info_array, magma_tally4_int_t gbstep,
    magma_tally4_int_t batchCount, cublasHandle_t myhandle, magma_tally4_queue_t queue);

magma_tally4_int_t 
magma_tally4_zpotrs_batched(
                  magma_tally4_uplo_t uplo, magma_tally4_int_t n, magma_tally4_int_t nrhs,
                  magma_tally4DoubleComplex **dA_array, magma_tally4_int_t ldda,
                  magma_tally4DoubleComplex **dB_array, magma_tally4_int_t lddb,
                  magma_tally4_int_t batchCount, magma_tally4_queue_t queue);

magma_tally4_int_t 
magma_tally4_zposv_batched(
                  magma_tally4_uplo_t uplo, magma_tally4_int_t n, magma_tally4_int_t nrhs,
                  magma_tally4DoubleComplex **dA_array, magma_tally4_int_t ldda,
                  magma_tally4DoubleComplex **dB_array, magma_tally4_int_t lddb,
                  magma_tally4_int_t *dinfo_array,
                  magma_tally4_int_t batchCount, magma_tally4_queue_t queue);

magma_tally4_int_t 
magma_tally4_zgetrs_batched(
                  magma_tally4_trans_t trans, magma_tally4_int_t n, magma_tally4_int_t nrhs,
                  magma_tally4DoubleComplex **dA_array, magma_tally4_int_t ldda,
                  magma_tally4_int_t **dipiv_array, 
                  magma_tally4DoubleComplex **dB_array, magma_tally4_int_t lddb,
                  magma_tally4_int_t batchCount, magma_tally4_queue_t queue);



void 
magma_tally4_zlaswp_rowparallel_batched( magma_tally4_int_t n, magma_tally4DoubleComplex** input_array, magma_tally4_int_t ldi,
                   magma_tally4DoubleComplex** output_array, magma_tally4_int_t ldo,
                   magma_tally4_int_t k1, magma_tally4_int_t k2,
                   magma_tally4_int_t **pivinfo_array, 
                   magma_tally4_int_t batchCount, magma_tally4_queue_t queue );

void 
magma_tally4_zlaswp_rowserial_batched(magma_tally4_int_t n, magma_tally4DoubleComplex** dA_array, magma_tally4_int_t lda,
                   magma_tally4_int_t k1, magma_tally4_int_t k2,
                   magma_tally4_int_t **ipiv_array, 
                   magma_tally4_int_t batchCount, magma_tally4_queue_t queue);

void 
magma_tally4_zlaswp_columnserial_batched(magma_tally4_int_t n, magma_tally4DoubleComplex** dA_array, magma_tally4_int_t lda,
                   magma_tally4_int_t k1, magma_tally4_int_t k2,
                   magma_tally4_int_t **ipiv_array, 
                   magma_tally4_int_t batchCount, magma_tally4_queue_t queue);

void 
magma_tally4blas_ztranspose_batched(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex **dA_array,  magma_tally4_int_t ldda,
    magma_tally4DoubleComplex **dAT_array, magma_tally4_int_t lddat, magma_tally4_int_t batchCount );


void 
magma_tally4blas_zlaset_batched(
    magma_tally4_uplo_t uplo, magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex offdiag, magma_tally4DoubleComplex diag,
    magma_tally4DoubleComplex_ptr dAarray[], magma_tally4_int_t ldda,
    magma_tally4_int_t batchCount, magma_tally4_queue_t queue);

void 
magma_tally4blas_zmemset_batched(magma_tally4_int_t length, 
        magma_tally4DoubleComplex_ptr dAarray[], magma_tally4DoubleComplex val, 
        magma_tally4_int_t batchCount, magma_tally4_queue_t queue);

magma_tally4_int_t 
magma_tally4_zgetf2_batched(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex **dA_array, magma_tally4_int_t lda,
    magma_tally4DoubleComplex **GERA_array,
    magma_tally4DoubleComplex **GERB_array,
    magma_tally4DoubleComplex **GERC_array,
    magma_tally4_int_t **ipiv_array,
    magma_tally4_int_t *info_array, 
    magma_tally4_int_t gbstep,            
    magma_tally4_int_t batchCount,
    cublasHandle_t myhandle, magma_tally4_queue_t queue);

magma_tally4_int_t 
magma_tally4_zgetrf_recpanel_batched(
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t min_recpnb,    
    magma_tally4DoubleComplex** dA_array,    magma_tally4_int_t ldda,
    magma_tally4_int_t** dipiv_array, magma_tally4_int_t** dpivinfo_array,
    magma_tally4DoubleComplex** dX_array,    magma_tally4_int_t dX_length,
    magma_tally4DoubleComplex** dinvA_array, magma_tally4_int_t dinvA_length,
    magma_tally4DoubleComplex** dW1_displ, magma_tally4DoubleComplex** dW2_displ,  
    magma_tally4DoubleComplex** dW3_displ, magma_tally4DoubleComplex** dW4_displ,
    magma_tally4DoubleComplex** dW5_displ,
    magma_tally4_int_t *info_array, magma_tally4_int_t gbstep,
    magma_tally4_int_t batchCount, cublasHandle_t myhandle, magma_tally4_queue_t queue);

magma_tally4_int_t 
magma_tally4_zgetrf_batched(
        magma_tally4_int_t m, magma_tally4_int_t n,
        magma_tally4DoubleComplex **dA_array, 
        magma_tally4_int_t lda,
        magma_tally4_int_t **ipiv_array, 
        magma_tally4_int_t *info_array, 
        magma_tally4_int_t batchCount, magma_tally4_queue_t queue);

magma_tally4_int_t 
magma_tally4_zgetri_outofplace_batched( magma_tally4_int_t n, 
                  magma_tally4DoubleComplex **dA_array, magma_tally4_int_t ldda,
                  magma_tally4_int_t **dipiv_array, 
                  magma_tally4DoubleComplex **dinvA_array, magma_tally4_int_t lddia,
                  magma_tally4_int_t *info_array,
                  magma_tally4_int_t batchCount, magma_tally4_queue_t queue);

void 
magma_tally4_zdisplace_intpointers(magma_tally4_int_t **output_array,
               magma_tally4_int_t **input_array, magma_tally4_int_t lda,
               magma_tally4_int_t row, magma_tally4_int_t column, 
               magma_tally4_int_t batchCount, magma_tally4_queue_t queue);





void 
magma_tally4blas_izamax_atomic_batched(magma_tally4_int_t n, magma_tally4DoubleComplex** x_array, magma_tally4_int_t incx, magma_tally4_int_t **max_id_array, magma_tally4_int_t batchCount);

void 
magma_tally4blas_izamax_tree_batched(magma_tally4_int_t n, magma_tally4DoubleComplex** x_array, magma_tally4_int_t incx, magma_tally4_int_t **max_id_array, magma_tally4_int_t batchCount);



void 
magma_tally4blas_izamax_batched(magma_tally4_int_t n, magma_tally4DoubleComplex** x_array, magma_tally4_int_t incx, magma_tally4_int_t **max_id_array, magma_tally4_int_t batchCount);

void 
magma_tally4blas_izamax(magma_tally4_int_t n, magma_tally4DoubleComplex* x, magma_tally4_int_t incx, magma_tally4_int_t *max_id);


magma_tally4_int_t 
magma_tally4_izamax_batched(magma_tally4_int_t length, 
        magma_tally4DoubleComplex **x_array, magma_tally4_int_t incx, magma_tally4_int_t step,  magma_tally4_int_t lda,
        magma_tally4_int_t** ipiv_array, magma_tally4_int_t *info_array, magma_tally4_int_t gbstep, magma_tally4_int_t batchCount, magma_tally4_queue_t queue);

magma_tally4_int_t 
magma_tally4_zswap_batched(magma_tally4_int_t n, magma_tally4DoubleComplex **x_array, magma_tally4_int_t incx, magma_tally4_int_t j, 
                 magma_tally4_int_t** ipiv_array, magma_tally4_int_t batchCount, magma_tally4_queue_t queue);

magma_tally4_int_t 
magma_tally4_zscal_zgeru_batched(magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t step,
                                      magma_tally4DoubleComplex **dA_array, magma_tally4_int_t lda,
                                      magma_tally4_int_t *info_array, magma_tally4_int_t gbstep, 
                                      magma_tally4_int_t batchCount, magma_tally4_queue_t queue);

magma_tally4_int_t 
magma_tally4_zcomputecolumn_batched(magma_tally4_int_t m, magma_tally4_int_t paneloffset, magma_tally4_int_t step, 
                                        magma_tally4DoubleComplex **dA_array,  magma_tally4_int_t lda,
                                        magma_tally4_int_t **ipiv_array, 
                                        magma_tally4_int_t *info_array, magma_tally4_int_t gbstep, 
                                        magma_tally4_int_t batchCount, magma_tally4_queue_t queue);

void 
magma_tally4_zgetf2trsm_batched(magma_tally4_int_t ib, magma_tally4_int_t n, magma_tally4DoubleComplex **dA_array,  magma_tally4_int_t j, magma_tally4_int_t lda,
                       magma_tally4_int_t batchCount, magma_tally4_queue_t queue);


magma_tally4_int_t 
magma_tally4_zgetf2_nopiv_batched(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex **dA_array, magma_tally4_int_t lda,
    magma_tally4DoubleComplex **dW0_displ,
    magma_tally4DoubleComplex **dW1_displ,
    magma_tally4DoubleComplex **dW2_displ,
    magma_tally4_int_t *info_array,            
    magma_tally4_int_t gbstep, 
    magma_tally4_int_t batchCount,
    cublasHandle_t myhandle, magma_tally4_queue_t queue);

magma_tally4_int_t 
magma_tally4_zgetrf_panel_nopiv_batched(
    magma_tally4_int_t m, magma_tally4_int_t nb,    
    magma_tally4DoubleComplex** dA_array,    magma_tally4_int_t ldda,
    magma_tally4DoubleComplex** dX_array,    magma_tally4_int_t dX_length,
    magma_tally4DoubleComplex** dinvA_array, magma_tally4_int_t dinvA_length,
    magma_tally4DoubleComplex** dW0_displ, magma_tally4DoubleComplex** dW1_displ,  
    magma_tally4DoubleComplex** dW2_displ, magma_tally4DoubleComplex** dW3_displ,
    magma_tally4DoubleComplex** dW4_displ,     
    magma_tally4_int_t *info_array, magma_tally4_int_t gbstep, 
    magma_tally4_int_t batchCount, cublasHandle_t myhandle, magma_tally4_queue_t queue);

magma_tally4_int_t 
magma_tally4_zgetrf_recpanel_nopiv_batched(
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t min_recpnb,    
    magma_tally4DoubleComplex** dA_array,    magma_tally4_int_t ldda,
    magma_tally4DoubleComplex** dX_array,    magma_tally4_int_t dX_length,
    magma_tally4DoubleComplex** dinvA_array, magma_tally4_int_t dinvA_length,
    magma_tally4DoubleComplex** dW1_displ, magma_tally4DoubleComplex** dW2_displ,  
    magma_tally4DoubleComplex** dW3_displ, magma_tally4DoubleComplex** dW4_displ,
    magma_tally4DoubleComplex** dW5_displ, 
    magma_tally4_int_t *info_array, magma_tally4_int_t gbstep,
    magma_tally4_int_t batchCount, cublasHandle_t myhandle, magma_tally4_queue_t queue);


magma_tally4_int_t 
magma_tally4_zgetrf_nopiv_batched(
        magma_tally4_int_t m, magma_tally4_int_t n,
        magma_tally4DoubleComplex **dA_array, 
        magma_tally4_int_t lda,
        magma_tally4_int_t *info_array, 
        magma_tally4_int_t batchCount, magma_tally4_queue_t queue);

magma_tally4_int_t 
magma_tally4_zgetrs_nopiv_batched(
                  magma_tally4_trans_t trans, magma_tally4_int_t n, magma_tally4_int_t nrhs,
                  magma_tally4DoubleComplex **dA_array, magma_tally4_int_t ldda,
                  magma_tally4DoubleComplex **dB_array, magma_tally4_int_t lddb,
                  magma_tally4_int_t *info_array,
                  magma_tally4_int_t batchCount, magma_tally4_queue_t queue);

magma_tally4_int_t 
magma_tally4_zgesv_nopiv_batched(
                  magma_tally4_int_t n, magma_tally4_int_t nrhs,
                  magma_tally4DoubleComplex **dA_array, magma_tally4_int_t ldda,
                  magma_tally4DoubleComplex **dB_array, magma_tally4_int_t lddb,
                  magma_tally4_int_t *info_array,
                  magma_tally4_int_t batchCount, magma_tally4_queue_t queue);

magma_tally4_int_t 
magma_tally4_zgesv_rbt_batched(
                  magma_tally4_int_t n, magma_tally4_int_t nrhs,
                  magma_tally4DoubleComplex **dA_array, magma_tally4_int_t ldda,
                  magma_tally4DoubleComplex **dB_array, magma_tally4_int_t lddb,
                  magma_tally4_int_t *info_array,
                  magma_tally4_int_t batchCount, magma_tally4_queue_t queue);

magma_tally4_int_t 
magma_tally4_zgesv_batched(
                  magma_tally4_int_t n, magma_tally4_int_t nrhs,
                  magma_tally4DoubleComplex **dA_array, magma_tally4_int_t ldda,
                  magma_tally4_int_t **dipiv_array, 
                  magma_tally4DoubleComplex **dB_array, magma_tally4_int_t lddb,
                  magma_tally4_int_t *dinfo_array,
                  magma_tally4_int_t batchCount, magma_tally4_queue_t queue);

magma_tally4_int_t
magma_tally4_zgerbt_batched(
    magma_tally4_bool_t gen, magma_tally4_int_t n, magma_tally4_int_t nrhs,
    magma_tally4DoubleComplex **dA_array, magma_tally4_int_t ldda,
    magma_tally4DoubleComplex **dB_array, magma_tally4_int_t lddb,
    magma_tally4DoubleComplex *U, magma_tally4DoubleComplex *V,
    magma_tally4_int_t *info, magma_tally4_int_t batchCount, magma_tally4_queue_t queue);

void 
magma_tally4blas_zprbt_batched(
    magma_tally4_int_t n, 
    magma_tally4DoubleComplex **dA_array, magma_tally4_int_t ldda, 
    magma_tally4DoubleComplex *du, magma_tally4DoubleComplex *dv,
    magma_tally4_int_t batchCount, magma_tally4_queue_t queue);

void
magma_tally4blas_zprbt_mv_batched(
    magma_tally4_int_t n, 
    magma_tally4DoubleComplex *dv, magma_tally4DoubleComplex **db_array, magma_tally4_int_t batchCount, magma_tally4_queue_t queue);


void
magma_tally4blas_zprbt_mtv_batched(
    magma_tally4_int_t n, 
    magma_tally4DoubleComplex *du, magma_tally4DoubleComplex **db_array, magma_tally4_int_t batchCount, magma_tally4_queue_t queue);





void 
magma_tally4_zlacgv_batched(magma_tally4_int_t n, magma_tally4DoubleComplex **x_array, magma_tally4_int_t incx, int offset, int batchCount, magma_tally4_queue_t queue);

void 
magma_tally4_zpotf2_zdscal_batched(magma_tally4_int_t n, magma_tally4DoubleComplex **x_array, magma_tally4_int_t incx, magma_tally4_int_t offset, magma_tally4_int_t *info_array, magma_tally4_int_t batchCount, magma_tally4_queue_t queue);

void 
magma_tally4_zpotf2_zdotc_batched(magma_tally4_int_t n, magma_tally4DoubleComplex **x_array, magma_tally4_int_t incx, magma_tally4_int_t offset, magma_tally4_int_t *info_array, magma_tally4_int_t gbstep, magma_tally4_int_t batchCount, magma_tally4_queue_t queue);


void 
setup_pivinfo( magma_tally4_int_t *pivinfo, magma_tally4_int_t *ipiv, 
                      magma_tally4_int_t m, magma_tally4_int_t nb, 
                      magma_tally4_queue_t queue);


void
magma_tally4blas_zgeadd_batched_q(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex alpha,
    magma_tally4DoubleComplex_const_ptr  const dAarray[], magma_tally4_int_t ldda,
    magma_tally4DoubleComplex_ptr              dBarray[], magma_tally4_int_t lddb,
    magma_tally4_int_t batchCount, magma_tally4_queue_t queue );

void
magma_tally4blas_zlacpy_batched(
    magma_tally4_uplo_t uplo, magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex_const_ptr  const dAarray[], magma_tally4_int_t ldda,
    magma_tally4DoubleComplex_ptr              dBarray[], magma_tally4_int_t lddb,
    magma_tally4_int_t batchCount, magma_tally4_queue_t queue );

void
magma_tally4blas_zgeadd_batched(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex alpha,
    magma_tally4DoubleComplex_const_ptr  const dAarray[], magma_tally4_int_t ldda,
    magma_tally4DoubleComplex_ptr              dBarray[], magma_tally4_int_t lddb,
    magma_tally4_int_t batchCount );


void
magma_tally4blas_zgemv_batched(
    magma_tally4_trans_t trans, magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex alpha,
    magma_tally4DoubleComplex_ptr dA_array[], magma_tally4_int_t ldda,
    magma_tally4DoubleComplex_ptr dx_array[], magma_tally4_int_t incx,
    magma_tally4DoubleComplex beta,
    magma_tally4DoubleComplex_ptr dy_array[], magma_tally4_int_t incy,
    magma_tally4_int_t batchCount, magma_tally4_queue_t queue);


magma_tally4_int_t 
magma_tally4_zgeqrf_batched(
    magma_tally4_int_t m, magma_tally4_int_t n, 
    magma_tally4DoubleComplex **dA_array,
    magma_tally4_int_t lda, 
    magma_tally4DoubleComplex **tau_array,
    magma_tally4_int_t *info_array, magma_tally4_int_t batchCount, magma_tally4_queue_t queue);

magma_tally4_int_t 
magma_tally4_zgeqrf_batched_v4(
    magma_tally4_int_t m, magma_tally4_int_t n, 
    magma_tally4DoubleComplex **dA_array,
    magma_tally4_int_t lda, 
    magma_tally4DoubleComplex **tau_array,
    magma_tally4_int_t *info_array, magma_tally4_int_t batchCount);

magma_tally4_int_t
magma_tally4_zgeqrf_panel_batched(
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t nb,    
    magma_tally4DoubleComplex** dA_array,    magma_tally4_int_t ldda,
    magma_tally4DoubleComplex** tau_array, 
    magma_tally4DoubleComplex** dT_array, magma_tally4_int_t ldt, 
    magma_tally4DoubleComplex** dR_array, magma_tally4_int_t ldr,
    magma_tally4DoubleComplex** dW0_displ, 
    magma_tally4DoubleComplex** dW1_displ,
    magma_tally4DoubleComplex *dwork,  
    magma_tally4DoubleComplex** W_array, 
    magma_tally4DoubleComplex** W2_array,
    magma_tally4_int_t *info_array,
    magma_tally4_int_t batchCount, cublasHandle_t myhandle, magma_tally4_queue_t queue);

magma_tally4_int_t
magma_tally4_zgeqrf_panel_batched_v4(
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t nb,    
    magma_tally4DoubleComplex** dA_array,    magma_tally4_int_t ldda,
    magma_tally4DoubleComplex** tau_array, 
    magma_tally4DoubleComplex** dT_array, magma_tally4_int_t ldt, 
    magma_tally4DoubleComplex** dR_array, magma_tally4_int_t ldr,
    double** dnorm_array,  
    magma_tally4DoubleComplex** dW0_displ, 
    magma_tally4DoubleComplex** dW1_displ,
    magma_tally4DoubleComplex *dwork,  
    magma_tally4DoubleComplex** W_array, 
    magma_tally4DoubleComplex** W2_array,
    magma_tally4_int_t *info_array,
    magma_tally4_int_t batchCount, cublasHandle_t myhandle);

magma_tally4_int_t
magma_tally4_zgeqr2x_batched_v4(magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4DoubleComplex **dA_array,
                  magma_tally4_int_t lda, 
                  magma_tally4DoubleComplex **tau_array,
                  magma_tally4DoubleComplex **dR_array, magma_tally4_int_t ldr,
                  double **dwork_array,  
                  magma_tally4_int_t *info, magma_tally4_int_t batchCount);

magma_tally4_int_t
magma_tally4_zgeqr2_batched(magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4DoubleComplex **dA_array,
                  magma_tally4_int_t lda, 
                  magma_tally4DoubleComplex **tau_array,
                  magma_tally4_int_t *info, magma_tally4_int_t batchCount, magma_tally4_queue_t queue);

magma_tally4_int_t
magma_tally4_zlarfb_zgemm_batched(
                  cublasHandle_t myhandle,
                  magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
                  magma_tally4DoubleComplex **dV_array,    magma_tally4_int_t ldv,
                  magma_tally4DoubleComplex **dT_array,    magma_tally4_int_t ldt,
                  magma_tally4DoubleComplex **dA_array,    magma_tally4_int_t lda,
                  magma_tally4DoubleComplex **W_array,     magma_tally4_int_t ldw,
                  magma_tally4DoubleComplex **W2_array,    magma_tally4_int_t ldw2,
                  magma_tally4_int_t batchCount, magma_tally4_queue_t queue);

magma_tally4_int_t
magma_tally4_zlarfb_gemm_batched(
    magma_tally4_side_t side, magma_tally4_trans_t trans, magma_tally4_direct_t direct, magma_tally4_storev_t storev,
    magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k,
    magma_tally4DoubleComplex_const_ptr dV_array[],    magma_tally4_int_t lddv,
    magma_tally4DoubleComplex_const_ptr dT_array[],    magma_tally4_int_t lddt,
    magma_tally4DoubleComplex_ptr dC_array[],          magma_tally4_int_t lddc,
    magma_tally4DoubleComplex_ptr dwork_array[],       magma_tally4_int_t ldwork,
    magma_tally4DoubleComplex_ptr dworkvt_array[],     magma_tally4_int_t ldworkvt,
    magma_tally4_int_t batchCount, cublasHandle_t myhandle, magma_tally4_queue_t queue);




void
magma_tally4_zlarft_batched_vold(magma_tally4_int_t n, magma_tally4_int_t k, magma_tally4DoubleComplex **v_array, magma_tally4_int_t ldv,
                    magma_tally4DoubleComplex **tau_array,
                    magma_tally4DoubleComplex **T_array, magma_tally4_int_t ldt, 
                    magma_tally4_int_t batchCount);





magma_tally4_int_t
magma_tally4_zlarft_batched(magma_tally4_int_t n, magma_tally4_int_t k, magma_tally4_int_t stair_T, 
                magma_tally4DoubleComplex **v_array, magma_tally4_int_t ldv,
                magma_tally4DoubleComplex **tau_array, magma_tally4DoubleComplex **T_array, magma_tally4_int_t ldt, 
                magma_tally4DoubleComplex **work_array, magma_tally4_int_t lwork, magma_tally4_int_t batchCount, cublasHandle_t myhandle, magma_tally4_queue_t queue);

void
magma_tally4_zlarft_sm32x32_batched(magma_tally4_int_t n, magma_tally4_int_t k, magma_tally4DoubleComplex **v_array, magma_tally4_int_t ldv,
                    magma_tally4DoubleComplex **tau_array, magma_tally4DoubleComplex **T_array, magma_tally4_int_t ldt, 
                    magma_tally4_int_t batchCount, cublasHandle_t myhandle, magma_tally4_queue_t queue);



void magma_tally4blas_zlarft_recztrmv_sm32x32(
    magma_tally4_int_t m, magma_tally4_int_t n, 
    magma_tally4DoubleComplex *tau, 
    magma_tally4DoubleComplex *Trec, magma_tally4_int_t ldtrec, 
    magma_tally4DoubleComplex *Ttri, magma_tally4_int_t ldttri);


void magma_tally4blas_zlarft_recztrmv_sm32x32_batched(
    magma_tally4_int_t m, magma_tally4_int_t n, 
    magma_tally4DoubleComplex **tau_array, 
    magma_tally4DoubleComplex **Trec_array, magma_tally4_int_t ldtrec, 
    magma_tally4DoubleComplex **Ttri_array, magma_tally4_int_t ldttri,
    magma_tally4_int_t batchCount, magma_tally4_queue_t queue);

void magma_tally4blas_zlarft_ztrmv_sm32x32(
    magma_tally4_int_t m, magma_tally4_int_t n, 
    magma_tally4DoubleComplex *tau, 
    magma_tally4DoubleComplex *Tin, magma_tally4_int_t ldtin, 
    magma_tally4DoubleComplex *Tout, magma_tally4_int_t ldtout);

void magma_tally4blas_zlarft_ztrmv_sm32x32_batched(
    magma_tally4_int_t m, magma_tally4_int_t n, 
    magma_tally4DoubleComplex **tau_array, 
    magma_tally4DoubleComplex **Tin_array, magma_tally4_int_t ldtin, 
    magma_tally4DoubleComplex **Tout_array, magma_tally4_int_t ldtout,
    magma_tally4_int_t batchCount, magma_tally4_queue_t queue);

void magma_tally4blas_zlarft_gemv_loop_inside(
    int n, int k, 
    magma_tally4DoubleComplex *tau, 
    magma_tally4DoubleComplex *v, int ldv, 
    magma_tally4DoubleComplex *T, int ldt);

void magma_tally4blas_zlarft_gemv_loop_inside_batched(
    int n, int k, 
    magma_tally4DoubleComplex **tau_array, 
    magma_tally4DoubleComplex **v_array, int ldv, 
    magma_tally4DoubleComplex **T_array, int ldt, magma_tally4_int_t batchCount, magma_tally4_queue_t queue);

void magma_tally4blas_zlarft_gemvrowwise(
    magma_tally4_int_t m, magma_tally4_int_t i, 
    magma_tally4DoubleComplex *tau, 
    magma_tally4DoubleComplex *v, magma_tally4_int_t ldv, 
    magma_tally4DoubleComplex *T, magma_tally4_int_t ldt,
    magma_tally4DoubleComplex *W);

void magma_tally4blas_zlarft_gemvrowwise_batched(
    magma_tally4_int_t m, magma_tally4_int_t i, 
    magma_tally4DoubleComplex **tau_array, 
    magma_tally4DoubleComplex **v_array, magma_tally4_int_t ldv, 
    magma_tally4DoubleComplex **T_array, magma_tally4_int_t ldt,
    magma_tally4_int_t batchCount, magma_tally4_queue_t queue);

void magma_tally4blas_zlarft_gemvcolwise(
    magma_tally4_int_t m,  magma_tally4_int_t step,
    magma_tally4DoubleComplex *v, magma_tally4_int_t ldv, 
    magma_tally4DoubleComplex *T,  magma_tally4_int_t ldt,
    magma_tally4DoubleComplex *tau);

void magma_tally4blas_zlarft_gemvcolwise_batched(
    magma_tally4_int_t m,  magma_tally4_int_t step,
    magma_tally4DoubleComplex **v_array, magma_tally4_int_t ldv, 
    magma_tally4DoubleComplex **T_array,  magma_tally4_int_t ldt,
    magma_tally4DoubleComplex **tau_array, magma_tally4_int_t batchCount, magma_tally4_queue_t queue);




void zgeqrf_copy_upper_batched(                
                  magma_tally4_int_t n, magma_tally4_int_t nb,
                  magma_tally4DoubleComplex **dV_array,    magma_tally4_int_t ldv,
                  magma_tally4DoubleComplex **dR_array,    magma_tally4_int_t ldr,
          magma_tally4_int_t batchCount, magma_tally4_queue_t queue);



void 
magma_tally4blas_dznrm2_cols_batched(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex **dA_array, magma_tally4_int_t lda, 
    double **dxnorm_array, magma_tally4_int_t batchCount);
 
void 
magma_tally4_zlarfgx_batched(magma_tally4_int_t n, magma_tally4DoubleComplex **dx0_array, magma_tally4DoubleComplex **dx_array, 
                  magma_tally4DoubleComplex **dtau_array, double **dxnorm_array, 
                  magma_tally4DoubleComplex **dR_array, magma_tally4_int_t it, magma_tally4_int_t batchCount);


void 
magma_tally4_zlarfx_batched_v4(magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4DoubleComplex **v_array, magma_tally4DoubleComplex **tau_array,
                magma_tally4DoubleComplex **C_array, magma_tally4_int_t ldc, double **xnorm_array, 
                magma_tally4_int_t step, 
                magma_tally4_int_t batchCount);


void 
magma_tally4blas_zlarfg_batched(
    magma_tally4_int_t n,
    magma_tally4DoubleComplex** dalpha_array, magma_tally4DoubleComplex** dx_array, magma_tally4_int_t incx,
    magma_tally4DoubleComplex** dtau_array, magma_tally4_int_t batchCount );





// for debugging purpose
void 
zset_stepinit_ipiv(magma_tally4_int_t **ipiv_array,
                 magma_tally4_int_t pm,
                 magma_tally4_int_t batchCount);



#ifdef __cplusplus
}
#endif

#undef COMPLEX

#endif  /* MAGMA_tally4_ZBATCHED_H */
