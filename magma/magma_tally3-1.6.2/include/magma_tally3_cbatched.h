/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Azzam Haidar
       @author Tingxing Dong

       @generated from magma_tally3_zbatched.h normal z -> c, Fri Jan 30 19:00:06 2015
*/

#ifndef MAGMA_tally3_CBATCHED_H
#define MAGMA_tally3_CBATCHED_H

#include "magma_tally3_types.h"

#define COMPLEX

#ifdef __cplusplus
extern "C" {
#endif
  /*
   *  local auxiliary routines
   */
void 
cset_pointer_int(magma_tally3_int_t **output_array,
        magma_tally3_int_t *input,
        magma_tally3_int_t lda,
        magma_tally3_int_t row, magma_tally3_int_t column, 
        magma_tally3_int_t batchSize,
        magma_tally3_int_t batchCount, magma_tally3_queue_t queue);

void 
cset_pointer(magma_tally3FloatComplex **output_array,
                 magma_tally3FloatComplex *input,
                 magma_tally3_int_t lda,
                 magma_tally3_int_t row, magma_tally3_int_t column,
                 magma_tally3_int_t batchSize,
                 magma_tally3_int_t batchCount, magma_tally3_queue_t queue);


void 
cset_array(magma_tally3FloatComplex **output_array,
               magma_tally3FloatComplex **input_array, magma_tally3_int_t lda,
               magma_tally3_int_t row, magma_tally3_int_t column,
               magma_tally3_int_t batchCount, magma_tally3_queue_t queue);

void 
magma_tally3_cdisplace_pointers(magma_tally3FloatComplex **output_array,
               magma_tally3FloatComplex **input_array, magma_tally3_int_t lda,
               magma_tally3_int_t row, magma_tally3_int_t column, 
               magma_tally3_int_t batchCount, magma_tally3_queue_t queue);

  /*
   *  LAPACK batched routines
   */

  /*
   *  BLAS batched routines
   */

void
magma_tally3blas_cgemm_batched(
    magma_tally3_trans_t transA, magma_tally3_trans_t transB, magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    magma_tally3FloatComplex alpha,
    magma_tally3FloatComplex const * const * dA_array, magma_tally3_int_t ldda,
    magma_tally3FloatComplex const * const * dB_array, magma_tally3_int_t lddb,
    magma_tally3FloatComplex beta,
    magma_tally3FloatComplex **dC_array, magma_tally3_int_t lddc, magma_tally3_int_t batchCount, magma_tally3_queue_t queue );

void
magma_tally3blas_cgemm_batched_lg(
    magma_tally3_trans_t transA, magma_tally3_trans_t transB, magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    magma_tally3FloatComplex alpha,
    magma_tally3FloatComplex const * const * dA_array, magma_tally3_int_t ldda,
    magma_tally3FloatComplex const * const * dB_array, magma_tally3_int_t lddb,
    magma_tally3FloatComplex beta,
    magma_tally3FloatComplex **dC_array, magma_tally3_int_t lddc, magma_tally3_int_t batchCount, magma_tally3_queue_t queue );
void
magma_tally3blas_cgemm_batched_k32(
    magma_tally3_trans_t transA, magma_tally3_trans_t transB, magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    magma_tally3FloatComplex alpha,
    magma_tally3FloatComplex const * const * dA_array, magma_tally3_int_t ldda,
    magma_tally3FloatComplex const * const * dB_array, magma_tally3_int_t lddb,
    magma_tally3FloatComplex beta,
    magma_tally3FloatComplex **dC_array, magma_tally3_int_t lddc, magma_tally3_int_t batchCount, magma_tally3_queue_t queue );



void 
magma_tally3blas_cherk_NC_batched( magma_tally3_trans_t TRANSA, magma_tally3_trans_t TRANSB, int m , int n , int k , 
                       magma_tally3FloatComplex alpha, magma_tally3FloatComplex **dA_array, int lda, 
                       magma_tally3FloatComplex **B_array, int ldb, 
                       magma_tally3FloatComplex beta,        magma_tally3FloatComplex **C_array, int ldc, 
                       magma_tally3_int_t batchCount);

void
magma_tally3blas_cherk_batched(
    magma_tally3_uplo_t uplo, magma_tally3_trans_t trans, magma_tally3_int_t n, magma_tally3_int_t k,
    float alpha,
    magma_tally3FloatComplex const * const * dA_array, magma_tally3_int_t ldda,
    float beta,
    magma_tally3FloatComplex **dC_array, magma_tally3_int_t lddc, magma_tally3_int_t batchCount, magma_tally3_queue_t queue );

void
magma_tally3blas_cherk_batched_lg(
    magma_tally3_uplo_t uplo, magma_tally3_trans_t trans, magma_tally3_int_t n, magma_tally3_int_t k,
    float alpha,
    magma_tally3FloatComplex const * const * dA_array, magma_tally3_int_t ldda,
    float beta,
    magma_tally3FloatComplex **dC_array, magma_tally3_int_t lddc, magma_tally3_int_t batchCount, magma_tally3_queue_t queue );

void
magma_tally3blas_cherk_batched_k32(
    magma_tally3_uplo_t uplo, magma_tally3_trans_t trans, magma_tally3_int_t n, magma_tally3_int_t k,
    float alpha,
    magma_tally3FloatComplex const * const * dA_array, magma_tally3_int_t ldda,
    float beta,
    magma_tally3FloatComplex **dC_array, magma_tally3_int_t lddc, magma_tally3_int_t batchCount, magma_tally3_queue_t queue );






magma_tally3_int_t 
magma_tally3_cpotf2_tile_batched(
    magma_tally3_uplo_t uplo, magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex **dA_array, magma_tally3_int_t lda,
    magma_tally3_int_t *info_array, magma_tally3_int_t gbstep, magma_tally3_int_t batchCount, magma_tally3_queue_t queue);

magma_tally3_int_t 
magma_tally3_cpotrf_panel(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n, magma_tally3_int_t nb,     
    magma_tally3FloatComplex *A, magma_tally3_int_t lda,
    magma_tally3FloatComplex* dinvA, magma_tally3FloatComplex** dinvA_array, magma_tally3_int_t invA_msize,
    magma_tally3FloatComplex* x, magma_tally3FloatComplex** x_array,  magma_tally3_int_t x_msize,
    magma_tally3_int_t *info_array,  magma_tally3_int_t batchCount, magma_tally3_int_t matrixSize, magma_tally3_queue_t queue);


void 
magma_tally3blas_ctrtri_diag_batched(
    magma_tally3_uplo_t uplo, magma_tally3_diag_t diag, magma_tally3_int_t n,
    magma_tally3FloatComplex const * const *dA_array, magma_tally3_int_t ldda,
    magma_tally3FloatComplex **dinvA_array, 
    magma_tally3_int_t resetozero, magma_tally3_int_t batchCount, magma_tally3_queue_t queue);



void 
magma_tally3blas_ctrsm_batched(
    magma_tally3_side_t side, magma_tally3_uplo_t uplo, magma_tally3_trans_t transA, magma_tally3_diag_t diag,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex alpha,
    magma_tally3FloatComplex** dA_array,    magma_tally3_int_t ldda,
    magma_tally3FloatComplex** dB_array,    magma_tally3_int_t lddb,
    magma_tally3_int_t batchCount, magma_tally3_queue_t queue);


void 
magma_tally3blas_ctrsm_work_batched(
    magma_tally3_side_t side, magma_tally3_uplo_t uplo, magma_tally3_trans_t transA, magma_tally3_diag_t diag,
    magma_tally3_int_t flag, magma_tally3_int_t m, magma_tally3_int_t n, 
    magma_tally3FloatComplex alpha, 
    magma_tally3FloatComplex** dA_array,    magma_tally3_int_t ldda,
    magma_tally3FloatComplex** dB_array,    magma_tally3_int_t lddb,
    magma_tally3FloatComplex** dX_array,    magma_tally3_int_t lddx, 
    magma_tally3FloatComplex** dinvA_array, magma_tally3_int_t dinvA_length,
    magma_tally3FloatComplex** dA_displ, magma_tally3FloatComplex** dB_displ, 
    magma_tally3FloatComplex** dX_displ, magma_tally3FloatComplex** dinvA_displ,
    magma_tally3_int_t resetozero, magma_tally3_int_t batchCount, magma_tally3_queue_t queue);

void 
magma_tally3blas_ctrsm_outofplace_batched(
    magma_tally3_side_t side, magma_tally3_uplo_t uplo, magma_tally3_trans_t transA, magma_tally3_diag_t diag,
    magma_tally3_int_t flag, magma_tally3_int_t m, magma_tally3_int_t n, 
    magma_tally3FloatComplex alpha, 
    magma_tally3FloatComplex** dA_array,    magma_tally3_int_t ldda,
    magma_tally3FloatComplex** dB_array,    magma_tally3_int_t lddb,
    magma_tally3FloatComplex** dX_array,    magma_tally3_int_t lddx, 
    magma_tally3FloatComplex** dinvA_array, magma_tally3_int_t dinvA_length,
    magma_tally3FloatComplex** dA_displ, magma_tally3FloatComplex** dB_displ, 
    magma_tally3FloatComplex** dX_displ, magma_tally3FloatComplex** dinvA_displ,
    magma_tally3_int_t resetozero, magma_tally3_int_t batchCount, magma_tally3_queue_t queue);

magma_tally3_int_t 
magma_tally3_cpotrf_batched(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3FloatComplex **dA_array, magma_tally3_int_t lda,
    magma_tally3_int_t *info_array,  magma_tally3_int_t batchCount, magma_tally3_queue_t queue);


magma_tally3_int_t 
magma_tally3_cpotf2_batched(
    magma_tally3_uplo_t uplo, magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex **dA_array, magma_tally3_int_t lda,
    magma_tally3FloatComplex **dA_displ, 
    magma_tally3FloatComplex **dW_displ,
    magma_tally3FloatComplex **dB_displ, 
    magma_tally3FloatComplex **dC_displ, 
    magma_tally3_int_t *info_array, magma_tally3_int_t gbstep, 
    magma_tally3_int_t batchCount, cublasHandle_t myhandle, magma_tally3_queue_t queue);

magma_tally3_int_t 
magma_tally3_cpotrf_panel_batched(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n, magma_tally3_int_t nb,     
    magma_tally3FloatComplex** dA_array,    magma_tally3_int_t ldda,
    magma_tally3FloatComplex** dX_array,    magma_tally3_int_t dX_length,
    magma_tally3FloatComplex** dinvA_array, magma_tally3_int_t dinvA_length,
    magma_tally3FloatComplex** dW0_displ, magma_tally3FloatComplex** dW1_displ, 
    magma_tally3FloatComplex** dW2_displ, magma_tally3FloatComplex** dW3_displ,
    magma_tally3FloatComplex** dW4_displ, 
    magma_tally3_int_t *info_array, magma_tally3_int_t gbstep,
    magma_tally3_int_t batchCount, cublasHandle_t myhandle, magma_tally3_queue_t queue);

magma_tally3_int_t 
magma_tally3_cpotrf_recpanel_batched(
    magma_tally3_uplo_t uplo, magma_tally3_int_t m, magma_tally3_int_t n, 
    magma_tally3_int_t min_recpnb,    
    magma_tally3FloatComplex** dA_array,    magma_tally3_int_t ldda,
    magma_tally3FloatComplex** dX_array,    magma_tally3_int_t dX_length,
    magma_tally3FloatComplex** dinvA_array, magma_tally3_int_t dinvA_length,
    magma_tally3FloatComplex** dW0_displ, magma_tally3FloatComplex** dW1_displ,  
    magma_tally3FloatComplex** dW2_displ, magma_tally3FloatComplex** dW3_displ,
    magma_tally3FloatComplex** dW4_displ,
    magma_tally3_int_t *info_array, magma_tally3_int_t gbstep, 
    magma_tally3_int_t batchCount, cublasHandle_t myhandle, magma_tally3_queue_t queue);

magma_tally3_int_t 
magma_tally3_cpotrf_rectile_batched(
    magma_tally3_uplo_t uplo, magma_tally3_int_t m, magma_tally3_int_t n, 
    magma_tally3_int_t min_recpnb,    
    magma_tally3FloatComplex** dA_array,    magma_tally3_int_t ldda,
    magma_tally3FloatComplex** dX_array,    magma_tally3_int_t dX_length,
    magma_tally3FloatComplex** dinvA_array, magma_tally3_int_t dinvA_length,
    magma_tally3FloatComplex** dW0_displ, magma_tally3FloatComplex** dW1_displ,  
    magma_tally3FloatComplex** dW2_displ, magma_tally3FloatComplex** dW3_displ,
    magma_tally3FloatComplex** dW4_displ,
    magma_tally3_int_t *info_array, magma_tally3_int_t gbstep,
    magma_tally3_int_t batchCount, cublasHandle_t myhandle, magma_tally3_queue_t queue);

magma_tally3_int_t 
magma_tally3_cpotrs_batched(
                  magma_tally3_uplo_t uplo, magma_tally3_int_t n, magma_tally3_int_t nrhs,
                  magma_tally3FloatComplex **dA_array, magma_tally3_int_t ldda,
                  magma_tally3FloatComplex **dB_array, magma_tally3_int_t lddb,
                  magma_tally3_int_t batchCount, magma_tally3_queue_t queue);

magma_tally3_int_t 
magma_tally3_cposv_batched(
                  magma_tally3_uplo_t uplo, magma_tally3_int_t n, magma_tally3_int_t nrhs,
                  magma_tally3FloatComplex **dA_array, magma_tally3_int_t ldda,
                  magma_tally3FloatComplex **dB_array, magma_tally3_int_t lddb,
                  magma_tally3_int_t *dinfo_array,
                  magma_tally3_int_t batchCount, magma_tally3_queue_t queue);

magma_tally3_int_t 
magma_tally3_cgetrs_batched(
                  magma_tally3_trans_t trans, magma_tally3_int_t n, magma_tally3_int_t nrhs,
                  magma_tally3FloatComplex **dA_array, magma_tally3_int_t ldda,
                  magma_tally3_int_t **dipiv_array, 
                  magma_tally3FloatComplex **dB_array, magma_tally3_int_t lddb,
                  magma_tally3_int_t batchCount, magma_tally3_queue_t queue);



void 
magma_tally3_claswp_rowparallel_batched( magma_tally3_int_t n, magma_tally3FloatComplex** input_array, magma_tally3_int_t ldi,
                   magma_tally3FloatComplex** output_array, magma_tally3_int_t ldo,
                   magma_tally3_int_t k1, magma_tally3_int_t k2,
                   magma_tally3_int_t **pivinfo_array, 
                   magma_tally3_int_t batchCount, magma_tally3_queue_t queue );

void 
magma_tally3_claswp_rowserial_batched(magma_tally3_int_t n, magma_tally3FloatComplex** dA_array, magma_tally3_int_t lda,
                   magma_tally3_int_t k1, magma_tally3_int_t k2,
                   magma_tally3_int_t **ipiv_array, 
                   magma_tally3_int_t batchCount, magma_tally3_queue_t queue);

void 
magma_tally3_claswp_columnserial_batched(magma_tally3_int_t n, magma_tally3FloatComplex** dA_array, magma_tally3_int_t lda,
                   magma_tally3_int_t k1, magma_tally3_int_t k2,
                   magma_tally3_int_t **ipiv_array, 
                   magma_tally3_int_t batchCount, magma_tally3_queue_t queue);

void 
magma_tally3blas_ctranspose_batched(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex **dA_array,  magma_tally3_int_t ldda,
    magma_tally3FloatComplex **dAT_array, magma_tally3_int_t lddat, magma_tally3_int_t batchCount );


void 
magma_tally3blas_claset_batched(
    magma_tally3_uplo_t uplo, magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex offdiag, magma_tally3FloatComplex diag,
    magma_tally3FloatComplex_ptr dAarray[], magma_tally3_int_t ldda,
    magma_tally3_int_t batchCount, magma_tally3_queue_t queue);

void 
magma_tally3blas_cmemset_batched(magma_tally3_int_t length, 
        magma_tally3FloatComplex_ptr dAarray[], magma_tally3FloatComplex val, 
        magma_tally3_int_t batchCount, magma_tally3_queue_t queue);

magma_tally3_int_t 
magma_tally3_cgetf2_batched(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex **dA_array, magma_tally3_int_t lda,
    magma_tally3FloatComplex **GERA_array,
    magma_tally3FloatComplex **GERB_array,
    magma_tally3FloatComplex **GERC_array,
    magma_tally3_int_t **ipiv_array,
    magma_tally3_int_t *info_array, 
    magma_tally3_int_t gbstep,            
    magma_tally3_int_t batchCount,
    cublasHandle_t myhandle, magma_tally3_queue_t queue);

magma_tally3_int_t 
magma_tally3_cgetrf_recpanel_batched(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t min_recpnb,    
    magma_tally3FloatComplex** dA_array,    magma_tally3_int_t ldda,
    magma_tally3_int_t** dipiv_array, magma_tally3_int_t** dpivinfo_array,
    magma_tally3FloatComplex** dX_array,    magma_tally3_int_t dX_length,
    magma_tally3FloatComplex** dinvA_array, magma_tally3_int_t dinvA_length,
    magma_tally3FloatComplex** dW1_displ, magma_tally3FloatComplex** dW2_displ,  
    magma_tally3FloatComplex** dW3_displ, magma_tally3FloatComplex** dW4_displ,
    magma_tally3FloatComplex** dW5_displ,
    magma_tally3_int_t *info_array, magma_tally3_int_t gbstep,
    magma_tally3_int_t batchCount, cublasHandle_t myhandle, magma_tally3_queue_t queue);

magma_tally3_int_t 
magma_tally3_cgetrf_batched(
        magma_tally3_int_t m, magma_tally3_int_t n,
        magma_tally3FloatComplex **dA_array, 
        magma_tally3_int_t lda,
        magma_tally3_int_t **ipiv_array, 
        magma_tally3_int_t *info_array, 
        magma_tally3_int_t batchCount, magma_tally3_queue_t queue);

magma_tally3_int_t 
magma_tally3_cgetri_outofplace_batched( magma_tally3_int_t n, 
                  magma_tally3FloatComplex **dA_array, magma_tally3_int_t ldda,
                  magma_tally3_int_t **dipiv_array, 
                  magma_tally3FloatComplex **dinvA_array, magma_tally3_int_t lddia,
                  magma_tally3_int_t *info_array,
                  magma_tally3_int_t batchCount, magma_tally3_queue_t queue);

void 
magma_tally3_cdisplace_intpointers(magma_tally3_int_t **output_array,
               magma_tally3_int_t **input_array, magma_tally3_int_t lda,
               magma_tally3_int_t row, magma_tally3_int_t column, 
               magma_tally3_int_t batchCount, magma_tally3_queue_t queue);





void 
magma_tally3blas_icamax_atomic_batched(magma_tally3_int_t n, magma_tally3FloatComplex** x_array, magma_tally3_int_t incx, magma_tally3_int_t **max_id_array, magma_tally3_int_t batchCount);

void 
magma_tally3blas_icamax_tree_batched(magma_tally3_int_t n, magma_tally3FloatComplex** x_array, magma_tally3_int_t incx, magma_tally3_int_t **max_id_array, magma_tally3_int_t batchCount);



void 
magma_tally3blas_icamax_batched(magma_tally3_int_t n, magma_tally3FloatComplex** x_array, magma_tally3_int_t incx, magma_tally3_int_t **max_id_array, magma_tally3_int_t batchCount);

void 
magma_tally3blas_icamax(magma_tally3_int_t n, magma_tally3FloatComplex* x, magma_tally3_int_t incx, magma_tally3_int_t *max_id);


magma_tally3_int_t 
magma_tally3_icamax_batched(magma_tally3_int_t length, 
        magma_tally3FloatComplex **x_array, magma_tally3_int_t incx, magma_tally3_int_t step,  magma_tally3_int_t lda,
        magma_tally3_int_t** ipiv_array, magma_tally3_int_t *info_array, magma_tally3_int_t gbstep, magma_tally3_int_t batchCount, magma_tally3_queue_t queue);

magma_tally3_int_t 
magma_tally3_cswap_batched(magma_tally3_int_t n, magma_tally3FloatComplex **x_array, magma_tally3_int_t incx, magma_tally3_int_t j, 
                 magma_tally3_int_t** ipiv_array, magma_tally3_int_t batchCount, magma_tally3_queue_t queue);

magma_tally3_int_t 
magma_tally3_cscal_cgeru_batched(magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t step,
                                      magma_tally3FloatComplex **dA_array, magma_tally3_int_t lda,
                                      magma_tally3_int_t *info_array, magma_tally3_int_t gbstep, 
                                      magma_tally3_int_t batchCount, magma_tally3_queue_t queue);

magma_tally3_int_t 
magma_tally3_ccomputecolumn_batched(magma_tally3_int_t m, magma_tally3_int_t paneloffset, magma_tally3_int_t step, 
                                        magma_tally3FloatComplex **dA_array,  magma_tally3_int_t lda,
                                        magma_tally3_int_t **ipiv_array, 
                                        magma_tally3_int_t *info_array, magma_tally3_int_t gbstep, 
                                        magma_tally3_int_t batchCount, magma_tally3_queue_t queue);

void 
magma_tally3_cgetf2trsm_batched(magma_tally3_int_t ib, magma_tally3_int_t n, magma_tally3FloatComplex **dA_array,  magma_tally3_int_t j, magma_tally3_int_t lda,
                       magma_tally3_int_t batchCount, magma_tally3_queue_t queue);


magma_tally3_int_t 
magma_tally3_cgetf2_nopiv_batched(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex **dA_array, magma_tally3_int_t lda,
    magma_tally3FloatComplex **dW0_displ,
    magma_tally3FloatComplex **dW1_displ,
    magma_tally3FloatComplex **dW2_displ,
    magma_tally3_int_t *info_array,            
    magma_tally3_int_t gbstep, 
    magma_tally3_int_t batchCount,
    cublasHandle_t myhandle, magma_tally3_queue_t queue);

magma_tally3_int_t 
magma_tally3_cgetrf_panel_nopiv_batched(
    magma_tally3_int_t m, magma_tally3_int_t nb,    
    magma_tally3FloatComplex** dA_array,    magma_tally3_int_t ldda,
    magma_tally3FloatComplex** dX_array,    magma_tally3_int_t dX_length,
    magma_tally3FloatComplex** dinvA_array, magma_tally3_int_t dinvA_length,
    magma_tally3FloatComplex** dW0_displ, magma_tally3FloatComplex** dW1_displ,  
    magma_tally3FloatComplex** dW2_displ, magma_tally3FloatComplex** dW3_displ,
    magma_tally3FloatComplex** dW4_displ,     
    magma_tally3_int_t *info_array, magma_tally3_int_t gbstep, 
    magma_tally3_int_t batchCount, cublasHandle_t myhandle, magma_tally3_queue_t queue);

magma_tally3_int_t 
magma_tally3_cgetrf_recpanel_nopiv_batched(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t min_recpnb,    
    magma_tally3FloatComplex** dA_array,    magma_tally3_int_t ldda,
    magma_tally3FloatComplex** dX_array,    magma_tally3_int_t dX_length,
    magma_tally3FloatComplex** dinvA_array, magma_tally3_int_t dinvA_length,
    magma_tally3FloatComplex** dW1_displ, magma_tally3FloatComplex** dW2_displ,  
    magma_tally3FloatComplex** dW3_displ, magma_tally3FloatComplex** dW4_displ,
    magma_tally3FloatComplex** dW5_displ, 
    magma_tally3_int_t *info_array, magma_tally3_int_t gbstep,
    magma_tally3_int_t batchCount, cublasHandle_t myhandle, magma_tally3_queue_t queue);


magma_tally3_int_t 
magma_tally3_cgetrf_nopiv_batched(
        magma_tally3_int_t m, magma_tally3_int_t n,
        magma_tally3FloatComplex **dA_array, 
        magma_tally3_int_t lda,
        magma_tally3_int_t *info_array, 
        magma_tally3_int_t batchCount, magma_tally3_queue_t queue);

magma_tally3_int_t 
magma_tally3_cgetrs_nopiv_batched(
                  magma_tally3_trans_t trans, magma_tally3_int_t n, magma_tally3_int_t nrhs,
                  magma_tally3FloatComplex **dA_array, magma_tally3_int_t ldda,
                  magma_tally3FloatComplex **dB_array, magma_tally3_int_t lddb,
                  magma_tally3_int_t *info_array,
                  magma_tally3_int_t batchCount, magma_tally3_queue_t queue);

magma_tally3_int_t 
magma_tally3_cgesv_nopiv_batched(
                  magma_tally3_int_t n, magma_tally3_int_t nrhs,
                  magma_tally3FloatComplex **dA_array, magma_tally3_int_t ldda,
                  magma_tally3FloatComplex **dB_array, magma_tally3_int_t lddb,
                  magma_tally3_int_t *info_array,
                  magma_tally3_int_t batchCount, magma_tally3_queue_t queue);

magma_tally3_int_t 
magma_tally3_cgesv_rbt_batched(
                  magma_tally3_int_t n, magma_tally3_int_t nrhs,
                  magma_tally3FloatComplex **dA_array, magma_tally3_int_t ldda,
                  magma_tally3FloatComplex **dB_array, magma_tally3_int_t lddb,
                  magma_tally3_int_t *info_array,
                  magma_tally3_int_t batchCount, magma_tally3_queue_t queue);

magma_tally3_int_t 
magma_tally3_cgesv_batched(
                  magma_tally3_int_t n, magma_tally3_int_t nrhs,
                  magma_tally3FloatComplex **dA_array, magma_tally3_int_t ldda,
                  magma_tally3_int_t **dipiv_array, 
                  magma_tally3FloatComplex **dB_array, magma_tally3_int_t lddb,
                  magma_tally3_int_t *dinfo_array,
                  magma_tally3_int_t batchCount, magma_tally3_queue_t queue);

magma_tally3_int_t
magma_tally3_cgerbt_batched(
    magma_tally3_bool_t gen, magma_tally3_int_t n, magma_tally3_int_t nrhs,
    magma_tally3FloatComplex **dA_array, magma_tally3_int_t ldda,
    magma_tally3FloatComplex **dB_array, magma_tally3_int_t lddb,
    magma_tally3FloatComplex *U, magma_tally3FloatComplex *V,
    magma_tally3_int_t *info, magma_tally3_int_t batchCount, magma_tally3_queue_t queue);

void 
magma_tally3blas_cprbt_batched(
    magma_tally3_int_t n, 
    magma_tally3FloatComplex **dA_array, magma_tally3_int_t ldda, 
    magma_tally3FloatComplex *du, magma_tally3FloatComplex *dv,
    magma_tally3_int_t batchCount, magma_tally3_queue_t queue);

void
magma_tally3blas_cprbt_mv_batched(
    magma_tally3_int_t n, 
    magma_tally3FloatComplex *dv, magma_tally3FloatComplex **db_array, magma_tally3_int_t batchCount, magma_tally3_queue_t queue);


void
magma_tally3blas_cprbt_mtv_batched(
    magma_tally3_int_t n, 
    magma_tally3FloatComplex *du, magma_tally3FloatComplex **db_array, magma_tally3_int_t batchCount, magma_tally3_queue_t queue);





void 
magma_tally3_clacgv_batched(magma_tally3_int_t n, magma_tally3FloatComplex **x_array, magma_tally3_int_t incx, int offset, int batchCount, magma_tally3_queue_t queue);

void 
magma_tally3_cpotf2_csscal_batched(magma_tally3_int_t n, magma_tally3FloatComplex **x_array, magma_tally3_int_t incx, magma_tally3_int_t offset, magma_tally3_int_t *info_array, magma_tally3_int_t batchCount, magma_tally3_queue_t queue);

void 
magma_tally3_cpotf2_cdotc_batched(magma_tally3_int_t n, magma_tally3FloatComplex **x_array, magma_tally3_int_t incx, magma_tally3_int_t offset, magma_tally3_int_t *info_array, magma_tally3_int_t gbstep, magma_tally3_int_t batchCount, magma_tally3_queue_t queue);


void 
setup_pivinfo( magma_tally3_int_t *pivinfo, magma_tally3_int_t *ipiv, 
                      magma_tally3_int_t m, magma_tally3_int_t nb, 
                      magma_tally3_queue_t queue);


void
magma_tally3blas_cgeadd_batched_q(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex alpha,
    magma_tally3FloatComplex_const_ptr  const dAarray[], magma_tally3_int_t ldda,
    magma_tally3FloatComplex_ptr              dBarray[], magma_tally3_int_t lddb,
    magma_tally3_int_t batchCount, magma_tally3_queue_t queue );

void
magma_tally3blas_clacpy_batched(
    magma_tally3_uplo_t uplo, magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex_const_ptr  const dAarray[], magma_tally3_int_t ldda,
    magma_tally3FloatComplex_ptr              dBarray[], magma_tally3_int_t lddb,
    magma_tally3_int_t batchCount, magma_tally3_queue_t queue );

void
magma_tally3blas_cgeadd_batched(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex alpha,
    magma_tally3FloatComplex_const_ptr  const dAarray[], magma_tally3_int_t ldda,
    magma_tally3FloatComplex_ptr              dBarray[], magma_tally3_int_t lddb,
    magma_tally3_int_t batchCount );


void
magma_tally3blas_cgemv_batched(
    magma_tally3_trans_t trans, magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex alpha,
    magma_tally3FloatComplex_ptr dA_array[], magma_tally3_int_t ldda,
    magma_tally3FloatComplex_ptr dx_array[], magma_tally3_int_t incx,
    magma_tally3FloatComplex beta,
    magma_tally3FloatComplex_ptr dy_array[], magma_tally3_int_t incy,
    magma_tally3_int_t batchCount, magma_tally3_queue_t queue);


magma_tally3_int_t 
magma_tally3_cgeqrf_batched(
    magma_tally3_int_t m, magma_tally3_int_t n, 
    magma_tally3FloatComplex **dA_array,
    magma_tally3_int_t lda, 
    magma_tally3FloatComplex **tau_array,
    magma_tally3_int_t *info_array, magma_tally3_int_t batchCount, magma_tally3_queue_t queue);

magma_tally3_int_t 
magma_tally3_cgeqrf_batched_v4(
    magma_tally3_int_t m, magma_tally3_int_t n, 
    magma_tally3FloatComplex **dA_array,
    magma_tally3_int_t lda, 
    magma_tally3FloatComplex **tau_array,
    magma_tally3_int_t *info_array, magma_tally3_int_t batchCount);

magma_tally3_int_t
magma_tally3_cgeqrf_panel_batched(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t nb,    
    magma_tally3FloatComplex** dA_array,    magma_tally3_int_t ldda,
    magma_tally3FloatComplex** tau_array, 
    magma_tally3FloatComplex** dT_array, magma_tally3_int_t ldt, 
    magma_tally3FloatComplex** dR_array, magma_tally3_int_t ldr,
    magma_tally3FloatComplex** dW0_displ, 
    magma_tally3FloatComplex** dW1_displ,
    magma_tally3FloatComplex *dwork,  
    magma_tally3FloatComplex** W_array, 
    magma_tally3FloatComplex** W2_array,
    magma_tally3_int_t *info_array,
    magma_tally3_int_t batchCount, cublasHandle_t myhandle, magma_tally3_queue_t queue);

magma_tally3_int_t
magma_tally3_cgeqrf_panel_batched_v4(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t nb,    
    magma_tally3FloatComplex** dA_array,    magma_tally3_int_t ldda,
    magma_tally3FloatComplex** tau_array, 
    magma_tally3FloatComplex** dT_array, magma_tally3_int_t ldt, 
    magma_tally3FloatComplex** dR_array, magma_tally3_int_t ldr,
    float** dnorm_array,  
    magma_tally3FloatComplex** dW0_displ, 
    magma_tally3FloatComplex** dW1_displ,
    magma_tally3FloatComplex *dwork,  
    magma_tally3FloatComplex** W_array, 
    magma_tally3FloatComplex** W2_array,
    magma_tally3_int_t *info_array,
    magma_tally3_int_t batchCount, cublasHandle_t myhandle);

magma_tally3_int_t
magma_tally3_cgeqr2x_batched_v4(magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3FloatComplex **dA_array,
                  magma_tally3_int_t lda, 
                  magma_tally3FloatComplex **tau_array,
                  magma_tally3FloatComplex **dR_array, magma_tally3_int_t ldr,
                  float **dwork_array,  
                  magma_tally3_int_t *info, magma_tally3_int_t batchCount);

magma_tally3_int_t
magma_tally3_cgeqr2_batched(magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3FloatComplex **dA_array,
                  magma_tally3_int_t lda, 
                  magma_tally3FloatComplex **tau_array,
                  magma_tally3_int_t *info, magma_tally3_int_t batchCount, magma_tally3_queue_t queue);

magma_tally3_int_t
magma_tally3_clarfb_cgemm_batched(
                  cublasHandle_t myhandle,
                  magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
                  magma_tally3FloatComplex **dV_array,    magma_tally3_int_t ldv,
                  magma_tally3FloatComplex **dT_array,    magma_tally3_int_t ldt,
                  magma_tally3FloatComplex **dA_array,    magma_tally3_int_t lda,
                  magma_tally3FloatComplex **W_array,     magma_tally3_int_t ldw,
                  magma_tally3FloatComplex **W2_array,    magma_tally3_int_t ldw2,
                  magma_tally3_int_t batchCount, magma_tally3_queue_t queue);

magma_tally3_int_t
magma_tally3_clarfb_gemm_batched(
    magma_tally3_side_t side, magma_tally3_trans_t trans, magma_tally3_direct_t direct, magma_tally3_storev_t storev,
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    magma_tally3FloatComplex_const_ptr dV_array[],    magma_tally3_int_t lddv,
    magma_tally3FloatComplex_const_ptr dT_array[],    magma_tally3_int_t lddt,
    magma_tally3FloatComplex_ptr dC_array[],          magma_tally3_int_t lddc,
    magma_tally3FloatComplex_ptr dwork_array[],       magma_tally3_int_t ldwork,
    magma_tally3FloatComplex_ptr dworkvt_array[],     magma_tally3_int_t ldworkvt,
    magma_tally3_int_t batchCount, cublasHandle_t myhandle, magma_tally3_queue_t queue);




void
magma_tally3_clarft_batched_vold(magma_tally3_int_t n, magma_tally3_int_t k, magma_tally3FloatComplex **v_array, magma_tally3_int_t ldv,
                    magma_tally3FloatComplex **tau_array,
                    magma_tally3FloatComplex **T_array, magma_tally3_int_t ldt, 
                    magma_tally3_int_t batchCount);





magma_tally3_int_t
magma_tally3_clarft_batched(magma_tally3_int_t n, magma_tally3_int_t k, magma_tally3_int_t stair_T, 
                magma_tally3FloatComplex **v_array, magma_tally3_int_t ldv,
                magma_tally3FloatComplex **tau_array, magma_tally3FloatComplex **T_array, magma_tally3_int_t ldt, 
                magma_tally3FloatComplex **work_array, magma_tally3_int_t lwork, magma_tally3_int_t batchCount, cublasHandle_t myhandle, magma_tally3_queue_t queue);

void
magma_tally3_clarft_sm32x32_batched(magma_tally3_int_t n, magma_tally3_int_t k, magma_tally3FloatComplex **v_array, magma_tally3_int_t ldv,
                    magma_tally3FloatComplex **tau_array, magma_tally3FloatComplex **T_array, magma_tally3_int_t ldt, 
                    magma_tally3_int_t batchCount, cublasHandle_t myhandle, magma_tally3_queue_t queue);



void magma_tally3blas_clarft_recctrmv_sm32x32(
    magma_tally3_int_t m, magma_tally3_int_t n, 
    magma_tally3FloatComplex *tau, 
    magma_tally3FloatComplex *Trec, magma_tally3_int_t ldtrec, 
    magma_tally3FloatComplex *Ttri, magma_tally3_int_t ldttri);


void magma_tally3blas_clarft_recctrmv_sm32x32_batched(
    magma_tally3_int_t m, magma_tally3_int_t n, 
    magma_tally3FloatComplex **tau_array, 
    magma_tally3FloatComplex **Trec_array, magma_tally3_int_t ldtrec, 
    magma_tally3FloatComplex **Ttri_array, magma_tally3_int_t ldttri,
    magma_tally3_int_t batchCount, magma_tally3_queue_t queue);

void magma_tally3blas_clarft_ctrmv_sm32x32(
    magma_tally3_int_t m, magma_tally3_int_t n, 
    magma_tally3FloatComplex *tau, 
    magma_tally3FloatComplex *Tin, magma_tally3_int_t ldtin, 
    magma_tally3FloatComplex *Tout, magma_tally3_int_t ldtout);

void magma_tally3blas_clarft_ctrmv_sm32x32_batched(
    magma_tally3_int_t m, magma_tally3_int_t n, 
    magma_tally3FloatComplex **tau_array, 
    magma_tally3FloatComplex **Tin_array, magma_tally3_int_t ldtin, 
    magma_tally3FloatComplex **Tout_array, magma_tally3_int_t ldtout,
    magma_tally3_int_t batchCount, magma_tally3_queue_t queue);

void magma_tally3blas_clarft_gemv_loop_inside(
    int n, int k, 
    magma_tally3FloatComplex *tau, 
    magma_tally3FloatComplex *v, int ldv, 
    magma_tally3FloatComplex *T, int ldt);

void magma_tally3blas_clarft_gemv_loop_inside_batched(
    int n, int k, 
    magma_tally3FloatComplex **tau_array, 
    magma_tally3FloatComplex **v_array, int ldv, 
    magma_tally3FloatComplex **T_array, int ldt, magma_tally3_int_t batchCount, magma_tally3_queue_t queue);

void magma_tally3blas_clarft_gemvrowwise(
    magma_tally3_int_t m, magma_tally3_int_t i, 
    magma_tally3FloatComplex *tau, 
    magma_tally3FloatComplex *v, magma_tally3_int_t ldv, 
    magma_tally3FloatComplex *T, magma_tally3_int_t ldt,
    magma_tally3FloatComplex *W);

void magma_tally3blas_clarft_gemvrowwise_batched(
    magma_tally3_int_t m, magma_tally3_int_t i, 
    magma_tally3FloatComplex **tau_array, 
    magma_tally3FloatComplex **v_array, magma_tally3_int_t ldv, 
    magma_tally3FloatComplex **T_array, magma_tally3_int_t ldt,
    magma_tally3_int_t batchCount, magma_tally3_queue_t queue);

void magma_tally3blas_clarft_gemvcolwise(
    magma_tally3_int_t m,  magma_tally3_int_t step,
    magma_tally3FloatComplex *v, magma_tally3_int_t ldv, 
    magma_tally3FloatComplex *T,  magma_tally3_int_t ldt,
    magma_tally3FloatComplex *tau);

void magma_tally3blas_clarft_gemvcolwise_batched(
    magma_tally3_int_t m,  magma_tally3_int_t step,
    magma_tally3FloatComplex **v_array, magma_tally3_int_t ldv, 
    magma_tally3FloatComplex **T_array,  magma_tally3_int_t ldt,
    magma_tally3FloatComplex **tau_array, magma_tally3_int_t batchCount, magma_tally3_queue_t queue);




void cgeqrf_copy_upper_batched(                
                  magma_tally3_int_t n, magma_tally3_int_t nb,
                  magma_tally3FloatComplex **dV_array,    magma_tally3_int_t ldv,
                  magma_tally3FloatComplex **dR_array,    magma_tally3_int_t ldr,
          magma_tally3_int_t batchCount, magma_tally3_queue_t queue);



void 
magma_tally3blas_scnrm2_cols_batched(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex **dA_array, magma_tally3_int_t lda, 
    float **dxnorm_array, magma_tally3_int_t batchCount);
 
void 
magma_tally3_clarfgx_batched(magma_tally3_int_t n, magma_tally3FloatComplex **dx0_array, magma_tally3FloatComplex **dx_array, 
                  magma_tally3FloatComplex **dtau_array, float **dxnorm_array, 
                  magma_tally3FloatComplex **dR_array, magma_tally3_int_t it, magma_tally3_int_t batchCount);


void 
magma_tally3_clarfx_batched_v4(magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3FloatComplex **v_array, magma_tally3FloatComplex **tau_array,
                magma_tally3FloatComplex **C_array, magma_tally3_int_t ldc, float **xnorm_array, 
                magma_tally3_int_t step, 
                magma_tally3_int_t batchCount);


void 
magma_tally3blas_clarfg_batched(
    magma_tally3_int_t n,
    magma_tally3FloatComplex** dalpha_array, magma_tally3FloatComplex** dx_array, magma_tally3_int_t incx,
    magma_tally3FloatComplex** dtau_array, magma_tally3_int_t batchCount );





// for debugging purpose
void 
cset_stepinit_ipiv(magma_tally3_int_t **ipiv_array,
                 magma_tally3_int_t pm,
                 magma_tally3_int_t batchCount);



#ifdef __cplusplus
}
#endif

#undef COMPLEX

#endif  /* MAGMA_tally3_CBATCHED_H */
