/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Azzam Haidar
       @author Tingxing Dong

       @generated from magma_tally2_zbatched.h normal z -> c, Fri Jan 30 19:00:06 2015
*/

#ifndef MAGMA_tally2_CBATCHED_H
#define MAGMA_tally2_CBATCHED_H

#include "magma_tally2_types.h"

#define COMPLEX

#ifdef __cplusplus
extern "C" {
#endif
  /*
   *  local auxiliary routines
   */
void 
cset_pointer_int(magma_tally2_int_t **output_array,
        magma_tally2_int_t *input,
        magma_tally2_int_t lda,
        magma_tally2_int_t row, magma_tally2_int_t column, 
        magma_tally2_int_t batchSize,
        magma_tally2_int_t batchCount, magma_tally2_queue_t queue);

void 
cset_pointer(magma_tally2FloatComplex **output_array,
                 magma_tally2FloatComplex *input,
                 magma_tally2_int_t lda,
                 magma_tally2_int_t row, magma_tally2_int_t column,
                 magma_tally2_int_t batchSize,
                 magma_tally2_int_t batchCount, magma_tally2_queue_t queue);


void 
cset_array(magma_tally2FloatComplex **output_array,
               magma_tally2FloatComplex **input_array, magma_tally2_int_t lda,
               magma_tally2_int_t row, magma_tally2_int_t column,
               magma_tally2_int_t batchCount, magma_tally2_queue_t queue);

void 
magma_tally2_cdisplace_pointers(magma_tally2FloatComplex **output_array,
               magma_tally2FloatComplex **input_array, magma_tally2_int_t lda,
               magma_tally2_int_t row, magma_tally2_int_t column, 
               magma_tally2_int_t batchCount, magma_tally2_queue_t queue);

  /*
   *  LAPACK batched routines
   */

  /*
   *  BLAS batched routines
   */

void
magma_tally2blas_cgemm_batched(
    magma_tally2_trans_t transA, magma_tally2_trans_t transB, magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    magma_tally2FloatComplex alpha,
    magma_tally2FloatComplex const * const * dA_array, magma_tally2_int_t ldda,
    magma_tally2FloatComplex const * const * dB_array, magma_tally2_int_t lddb,
    magma_tally2FloatComplex beta,
    magma_tally2FloatComplex **dC_array, magma_tally2_int_t lddc, magma_tally2_int_t batchCount, magma_tally2_queue_t queue );

void
magma_tally2blas_cgemm_batched_lg(
    magma_tally2_trans_t transA, magma_tally2_trans_t transB, magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    magma_tally2FloatComplex alpha,
    magma_tally2FloatComplex const * const * dA_array, magma_tally2_int_t ldda,
    magma_tally2FloatComplex const * const * dB_array, magma_tally2_int_t lddb,
    magma_tally2FloatComplex beta,
    magma_tally2FloatComplex **dC_array, magma_tally2_int_t lddc, magma_tally2_int_t batchCount, magma_tally2_queue_t queue );
void
magma_tally2blas_cgemm_batched_k32(
    magma_tally2_trans_t transA, magma_tally2_trans_t transB, magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    magma_tally2FloatComplex alpha,
    magma_tally2FloatComplex const * const * dA_array, magma_tally2_int_t ldda,
    magma_tally2FloatComplex const * const * dB_array, magma_tally2_int_t lddb,
    magma_tally2FloatComplex beta,
    magma_tally2FloatComplex **dC_array, magma_tally2_int_t lddc, magma_tally2_int_t batchCount, magma_tally2_queue_t queue );



void 
magma_tally2blas_cherk_NC_batched( magma_tally2_trans_t TRANSA, magma_tally2_trans_t TRANSB, int m , int n , int k , 
                       magma_tally2FloatComplex alpha, magma_tally2FloatComplex **dA_array, int lda, 
                       magma_tally2FloatComplex **B_array, int ldb, 
                       magma_tally2FloatComplex beta,        magma_tally2FloatComplex **C_array, int ldc, 
                       magma_tally2_int_t batchCount);

void
magma_tally2blas_cherk_batched(
    magma_tally2_uplo_t uplo, magma_tally2_trans_t trans, magma_tally2_int_t n, magma_tally2_int_t k,
    float alpha,
    magma_tally2FloatComplex const * const * dA_array, magma_tally2_int_t ldda,
    float beta,
    magma_tally2FloatComplex **dC_array, magma_tally2_int_t lddc, magma_tally2_int_t batchCount, magma_tally2_queue_t queue );

void
magma_tally2blas_cherk_batched_lg(
    magma_tally2_uplo_t uplo, magma_tally2_trans_t trans, magma_tally2_int_t n, magma_tally2_int_t k,
    float alpha,
    magma_tally2FloatComplex const * const * dA_array, magma_tally2_int_t ldda,
    float beta,
    magma_tally2FloatComplex **dC_array, magma_tally2_int_t lddc, magma_tally2_int_t batchCount, magma_tally2_queue_t queue );

void
magma_tally2blas_cherk_batched_k32(
    magma_tally2_uplo_t uplo, magma_tally2_trans_t trans, magma_tally2_int_t n, magma_tally2_int_t k,
    float alpha,
    magma_tally2FloatComplex const * const * dA_array, magma_tally2_int_t ldda,
    float beta,
    magma_tally2FloatComplex **dC_array, magma_tally2_int_t lddc, magma_tally2_int_t batchCount, magma_tally2_queue_t queue );






magma_tally2_int_t 
magma_tally2_cpotf2_tile_batched(
    magma_tally2_uplo_t uplo, magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2FloatComplex **dA_array, magma_tally2_int_t lda,
    magma_tally2_int_t *info_array, magma_tally2_int_t gbstep, magma_tally2_int_t batchCount, magma_tally2_queue_t queue);

magma_tally2_int_t 
magma_tally2_cpotrf_panel(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n, magma_tally2_int_t nb,     
    magma_tally2FloatComplex *A, magma_tally2_int_t lda,
    magma_tally2FloatComplex* dinvA, magma_tally2FloatComplex** dinvA_array, magma_tally2_int_t invA_msize,
    magma_tally2FloatComplex* x, magma_tally2FloatComplex** x_array,  magma_tally2_int_t x_msize,
    magma_tally2_int_t *info_array,  magma_tally2_int_t batchCount, magma_tally2_int_t matrixSize, magma_tally2_queue_t queue);


void 
magma_tally2blas_ctrtri_diag_batched(
    magma_tally2_uplo_t uplo, magma_tally2_diag_t diag, magma_tally2_int_t n,
    magma_tally2FloatComplex const * const *dA_array, magma_tally2_int_t ldda,
    magma_tally2FloatComplex **dinvA_array, 
    magma_tally2_int_t resetozero, magma_tally2_int_t batchCount, magma_tally2_queue_t queue);



void 
magma_tally2blas_ctrsm_batched(
    magma_tally2_side_t side, magma_tally2_uplo_t uplo, magma_tally2_trans_t transA, magma_tally2_diag_t diag,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2FloatComplex alpha,
    magma_tally2FloatComplex** dA_array,    magma_tally2_int_t ldda,
    magma_tally2FloatComplex** dB_array,    magma_tally2_int_t lddb,
    magma_tally2_int_t batchCount, magma_tally2_queue_t queue);


void 
magma_tally2blas_ctrsm_work_batched(
    magma_tally2_side_t side, magma_tally2_uplo_t uplo, magma_tally2_trans_t transA, magma_tally2_diag_t diag,
    magma_tally2_int_t flag, magma_tally2_int_t m, magma_tally2_int_t n, 
    magma_tally2FloatComplex alpha, 
    magma_tally2FloatComplex** dA_array,    magma_tally2_int_t ldda,
    magma_tally2FloatComplex** dB_array,    magma_tally2_int_t lddb,
    magma_tally2FloatComplex** dX_array,    magma_tally2_int_t lddx, 
    magma_tally2FloatComplex** dinvA_array, magma_tally2_int_t dinvA_length,
    magma_tally2FloatComplex** dA_displ, magma_tally2FloatComplex** dB_displ, 
    magma_tally2FloatComplex** dX_displ, magma_tally2FloatComplex** dinvA_displ,
    magma_tally2_int_t resetozero, magma_tally2_int_t batchCount, magma_tally2_queue_t queue);

void 
magma_tally2blas_ctrsm_outofplace_batched(
    magma_tally2_side_t side, magma_tally2_uplo_t uplo, magma_tally2_trans_t transA, magma_tally2_diag_t diag,
    magma_tally2_int_t flag, magma_tally2_int_t m, magma_tally2_int_t n, 
    magma_tally2FloatComplex alpha, 
    magma_tally2FloatComplex** dA_array,    magma_tally2_int_t ldda,
    magma_tally2FloatComplex** dB_array,    magma_tally2_int_t lddb,
    magma_tally2FloatComplex** dX_array,    magma_tally2_int_t lddx, 
    magma_tally2FloatComplex** dinvA_array, magma_tally2_int_t dinvA_length,
    magma_tally2FloatComplex** dA_displ, magma_tally2FloatComplex** dB_displ, 
    magma_tally2FloatComplex** dX_displ, magma_tally2FloatComplex** dinvA_displ,
    magma_tally2_int_t resetozero, magma_tally2_int_t batchCount, magma_tally2_queue_t queue);

magma_tally2_int_t 
magma_tally2_cpotrf_batched(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n,
    magma_tally2FloatComplex **dA_array, magma_tally2_int_t lda,
    magma_tally2_int_t *info_array,  magma_tally2_int_t batchCount, magma_tally2_queue_t queue);


magma_tally2_int_t 
magma_tally2_cpotf2_batched(
    magma_tally2_uplo_t uplo, magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2FloatComplex **dA_array, magma_tally2_int_t lda,
    magma_tally2FloatComplex **dA_displ, 
    magma_tally2FloatComplex **dW_displ,
    magma_tally2FloatComplex **dB_displ, 
    magma_tally2FloatComplex **dC_displ, 
    magma_tally2_int_t *info_array, magma_tally2_int_t gbstep, 
    magma_tally2_int_t batchCount, cublasHandle_t myhandle, magma_tally2_queue_t queue);

magma_tally2_int_t 
magma_tally2_cpotrf_panel_batched(
    magma_tally2_uplo_t uplo, magma_tally2_int_t n, magma_tally2_int_t nb,     
    magma_tally2FloatComplex** dA_array,    magma_tally2_int_t ldda,
    magma_tally2FloatComplex** dX_array,    magma_tally2_int_t dX_length,
    magma_tally2FloatComplex** dinvA_array, magma_tally2_int_t dinvA_length,
    magma_tally2FloatComplex** dW0_displ, magma_tally2FloatComplex** dW1_displ, 
    magma_tally2FloatComplex** dW2_displ, magma_tally2FloatComplex** dW3_displ,
    magma_tally2FloatComplex** dW4_displ, 
    magma_tally2_int_t *info_array, magma_tally2_int_t gbstep,
    magma_tally2_int_t batchCount, cublasHandle_t myhandle, magma_tally2_queue_t queue);

magma_tally2_int_t 
magma_tally2_cpotrf_recpanel_batched(
    magma_tally2_uplo_t uplo, magma_tally2_int_t m, magma_tally2_int_t n, 
    magma_tally2_int_t min_recpnb,    
    magma_tally2FloatComplex** dA_array,    magma_tally2_int_t ldda,
    magma_tally2FloatComplex** dX_array,    magma_tally2_int_t dX_length,
    magma_tally2FloatComplex** dinvA_array, magma_tally2_int_t dinvA_length,
    magma_tally2FloatComplex** dW0_displ, magma_tally2FloatComplex** dW1_displ,  
    magma_tally2FloatComplex** dW2_displ, magma_tally2FloatComplex** dW3_displ,
    magma_tally2FloatComplex** dW4_displ,
    magma_tally2_int_t *info_array, magma_tally2_int_t gbstep, 
    magma_tally2_int_t batchCount, cublasHandle_t myhandle, magma_tally2_queue_t queue);

magma_tally2_int_t 
magma_tally2_cpotrf_rectile_batched(
    magma_tally2_uplo_t uplo, magma_tally2_int_t m, magma_tally2_int_t n, 
    magma_tally2_int_t min_recpnb,    
    magma_tally2FloatComplex** dA_array,    magma_tally2_int_t ldda,
    magma_tally2FloatComplex** dX_array,    magma_tally2_int_t dX_length,
    magma_tally2FloatComplex** dinvA_array, magma_tally2_int_t dinvA_length,
    magma_tally2FloatComplex** dW0_displ, magma_tally2FloatComplex** dW1_displ,  
    magma_tally2FloatComplex** dW2_displ, magma_tally2FloatComplex** dW3_displ,
    magma_tally2FloatComplex** dW4_displ,
    magma_tally2_int_t *info_array, magma_tally2_int_t gbstep,
    magma_tally2_int_t batchCount, cublasHandle_t myhandle, magma_tally2_queue_t queue);

magma_tally2_int_t 
magma_tally2_cpotrs_batched(
                  magma_tally2_uplo_t uplo, magma_tally2_int_t n, magma_tally2_int_t nrhs,
                  magma_tally2FloatComplex **dA_array, magma_tally2_int_t ldda,
                  magma_tally2FloatComplex **dB_array, magma_tally2_int_t lddb,
                  magma_tally2_int_t batchCount, magma_tally2_queue_t queue);

magma_tally2_int_t 
magma_tally2_cposv_batched(
                  magma_tally2_uplo_t uplo, magma_tally2_int_t n, magma_tally2_int_t nrhs,
                  magma_tally2FloatComplex **dA_array, magma_tally2_int_t ldda,
                  magma_tally2FloatComplex **dB_array, magma_tally2_int_t lddb,
                  magma_tally2_int_t *dinfo_array,
                  magma_tally2_int_t batchCount, magma_tally2_queue_t queue);

magma_tally2_int_t 
magma_tally2_cgetrs_batched(
                  magma_tally2_trans_t trans, magma_tally2_int_t n, magma_tally2_int_t nrhs,
                  magma_tally2FloatComplex **dA_array, magma_tally2_int_t ldda,
                  magma_tally2_int_t **dipiv_array, 
                  magma_tally2FloatComplex **dB_array, magma_tally2_int_t lddb,
                  magma_tally2_int_t batchCount, magma_tally2_queue_t queue);



void 
magma_tally2_claswp_rowparallel_batched( magma_tally2_int_t n, magma_tally2FloatComplex** input_array, magma_tally2_int_t ldi,
                   magma_tally2FloatComplex** output_array, magma_tally2_int_t ldo,
                   magma_tally2_int_t k1, magma_tally2_int_t k2,
                   magma_tally2_int_t **pivinfo_array, 
                   magma_tally2_int_t batchCount, magma_tally2_queue_t queue );

void 
magma_tally2_claswp_rowserial_batched(magma_tally2_int_t n, magma_tally2FloatComplex** dA_array, magma_tally2_int_t lda,
                   magma_tally2_int_t k1, magma_tally2_int_t k2,
                   magma_tally2_int_t **ipiv_array, 
                   magma_tally2_int_t batchCount, magma_tally2_queue_t queue);

void 
magma_tally2_claswp_columnserial_batched(magma_tally2_int_t n, magma_tally2FloatComplex** dA_array, magma_tally2_int_t lda,
                   magma_tally2_int_t k1, magma_tally2_int_t k2,
                   magma_tally2_int_t **ipiv_array, 
                   magma_tally2_int_t batchCount, magma_tally2_queue_t queue);

void 
magma_tally2blas_ctranspose_batched(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2FloatComplex **dA_array,  magma_tally2_int_t ldda,
    magma_tally2FloatComplex **dAT_array, magma_tally2_int_t lddat, magma_tally2_int_t batchCount );


void 
magma_tally2blas_claset_batched(
    magma_tally2_uplo_t uplo, magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2FloatComplex offdiag, magma_tally2FloatComplex diag,
    magma_tally2FloatComplex_ptr dAarray[], magma_tally2_int_t ldda,
    magma_tally2_int_t batchCount, magma_tally2_queue_t queue);

void 
magma_tally2blas_cmemset_batched(magma_tally2_int_t length, 
        magma_tally2FloatComplex_ptr dAarray[], magma_tally2FloatComplex val, 
        magma_tally2_int_t batchCount, magma_tally2_queue_t queue);

magma_tally2_int_t 
magma_tally2_cgetf2_batched(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2FloatComplex **dA_array, magma_tally2_int_t lda,
    magma_tally2FloatComplex **GERA_array,
    magma_tally2FloatComplex **GERB_array,
    magma_tally2FloatComplex **GERC_array,
    magma_tally2_int_t **ipiv_array,
    magma_tally2_int_t *info_array, 
    magma_tally2_int_t gbstep,            
    magma_tally2_int_t batchCount,
    cublasHandle_t myhandle, magma_tally2_queue_t queue);

magma_tally2_int_t 
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
    magma_tally2_int_t batchCount, cublasHandle_t myhandle, magma_tally2_queue_t queue);

magma_tally2_int_t 
magma_tally2_cgetrf_batched(
        magma_tally2_int_t m, magma_tally2_int_t n,
        magma_tally2FloatComplex **dA_array, 
        magma_tally2_int_t lda,
        magma_tally2_int_t **ipiv_array, 
        magma_tally2_int_t *info_array, 
        magma_tally2_int_t batchCount, magma_tally2_queue_t queue);

magma_tally2_int_t 
magma_tally2_cgetri_outofplace_batched( magma_tally2_int_t n, 
                  magma_tally2FloatComplex **dA_array, magma_tally2_int_t ldda,
                  magma_tally2_int_t **dipiv_array, 
                  magma_tally2FloatComplex **dinvA_array, magma_tally2_int_t lddia,
                  magma_tally2_int_t *info_array,
                  magma_tally2_int_t batchCount, magma_tally2_queue_t queue);

void 
magma_tally2_cdisplace_intpointers(magma_tally2_int_t **output_array,
               magma_tally2_int_t **input_array, magma_tally2_int_t lda,
               magma_tally2_int_t row, magma_tally2_int_t column, 
               magma_tally2_int_t batchCount, magma_tally2_queue_t queue);





void 
magma_tally2blas_icamax_atomic_batched(magma_tally2_int_t n, magma_tally2FloatComplex** x_array, magma_tally2_int_t incx, magma_tally2_int_t **max_id_array, magma_tally2_int_t batchCount);

void 
magma_tally2blas_icamax_tree_batched(magma_tally2_int_t n, magma_tally2FloatComplex** x_array, magma_tally2_int_t incx, magma_tally2_int_t **max_id_array, magma_tally2_int_t batchCount);



void 
magma_tally2blas_icamax_batched(magma_tally2_int_t n, magma_tally2FloatComplex** x_array, magma_tally2_int_t incx, magma_tally2_int_t **max_id_array, magma_tally2_int_t batchCount);

void 
magma_tally2blas_icamax(magma_tally2_int_t n, magma_tally2FloatComplex* x, magma_tally2_int_t incx, magma_tally2_int_t *max_id);


magma_tally2_int_t 
magma_tally2_icamax_batched(magma_tally2_int_t length, 
        magma_tally2FloatComplex **x_array, magma_tally2_int_t incx, magma_tally2_int_t step,  magma_tally2_int_t lda,
        magma_tally2_int_t** ipiv_array, magma_tally2_int_t *info_array, magma_tally2_int_t gbstep, magma_tally2_int_t batchCount, magma_tally2_queue_t queue);

magma_tally2_int_t 
magma_tally2_cswap_batched(magma_tally2_int_t n, magma_tally2FloatComplex **x_array, magma_tally2_int_t incx, magma_tally2_int_t j, 
                 magma_tally2_int_t** ipiv_array, magma_tally2_int_t batchCount, magma_tally2_queue_t queue);

magma_tally2_int_t 
magma_tally2_cscal_cgeru_batched(magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t step,
                                      magma_tally2FloatComplex **dA_array, magma_tally2_int_t lda,
                                      magma_tally2_int_t *info_array, magma_tally2_int_t gbstep, 
                                      magma_tally2_int_t batchCount, magma_tally2_queue_t queue);

magma_tally2_int_t 
magma_tally2_ccomputecolumn_batched(magma_tally2_int_t m, magma_tally2_int_t paneloffset, magma_tally2_int_t step, 
                                        magma_tally2FloatComplex **dA_array,  magma_tally2_int_t lda,
                                        magma_tally2_int_t **ipiv_array, 
                                        magma_tally2_int_t *info_array, magma_tally2_int_t gbstep, 
                                        magma_tally2_int_t batchCount, magma_tally2_queue_t queue);

void 
magma_tally2_cgetf2trsm_batched(magma_tally2_int_t ib, magma_tally2_int_t n, magma_tally2FloatComplex **dA_array,  magma_tally2_int_t j, magma_tally2_int_t lda,
                       magma_tally2_int_t batchCount, magma_tally2_queue_t queue);


magma_tally2_int_t 
magma_tally2_cgetf2_nopiv_batched(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2FloatComplex **dA_array, magma_tally2_int_t lda,
    magma_tally2FloatComplex **dW0_displ,
    magma_tally2FloatComplex **dW1_displ,
    magma_tally2FloatComplex **dW2_displ,
    magma_tally2_int_t *info_array,            
    magma_tally2_int_t gbstep, 
    magma_tally2_int_t batchCount,
    cublasHandle_t myhandle, magma_tally2_queue_t queue);

magma_tally2_int_t 
magma_tally2_cgetrf_panel_nopiv_batched(
    magma_tally2_int_t m, magma_tally2_int_t nb,    
    magma_tally2FloatComplex** dA_array,    magma_tally2_int_t ldda,
    magma_tally2FloatComplex** dX_array,    magma_tally2_int_t dX_length,
    magma_tally2FloatComplex** dinvA_array, magma_tally2_int_t dinvA_length,
    magma_tally2FloatComplex** dW0_displ, magma_tally2FloatComplex** dW1_displ,  
    magma_tally2FloatComplex** dW2_displ, magma_tally2FloatComplex** dW3_displ,
    magma_tally2FloatComplex** dW4_displ,     
    magma_tally2_int_t *info_array, magma_tally2_int_t gbstep, 
    magma_tally2_int_t batchCount, cublasHandle_t myhandle, magma_tally2_queue_t queue);

magma_tally2_int_t 
magma_tally2_cgetrf_recpanel_nopiv_batched(
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t min_recpnb,    
    magma_tally2FloatComplex** dA_array,    magma_tally2_int_t ldda,
    magma_tally2FloatComplex** dX_array,    magma_tally2_int_t dX_length,
    magma_tally2FloatComplex** dinvA_array, magma_tally2_int_t dinvA_length,
    magma_tally2FloatComplex** dW1_displ, magma_tally2FloatComplex** dW2_displ,  
    magma_tally2FloatComplex** dW3_displ, magma_tally2FloatComplex** dW4_displ,
    magma_tally2FloatComplex** dW5_displ, 
    magma_tally2_int_t *info_array, magma_tally2_int_t gbstep,
    magma_tally2_int_t batchCount, cublasHandle_t myhandle, magma_tally2_queue_t queue);


magma_tally2_int_t 
magma_tally2_cgetrf_nopiv_batched(
        magma_tally2_int_t m, magma_tally2_int_t n,
        magma_tally2FloatComplex **dA_array, 
        magma_tally2_int_t lda,
        magma_tally2_int_t *info_array, 
        magma_tally2_int_t batchCount, magma_tally2_queue_t queue);

magma_tally2_int_t 
magma_tally2_cgetrs_nopiv_batched(
                  magma_tally2_trans_t trans, magma_tally2_int_t n, magma_tally2_int_t nrhs,
                  magma_tally2FloatComplex **dA_array, magma_tally2_int_t ldda,
                  magma_tally2FloatComplex **dB_array, magma_tally2_int_t lddb,
                  magma_tally2_int_t *info_array,
                  magma_tally2_int_t batchCount, magma_tally2_queue_t queue);

magma_tally2_int_t 
magma_tally2_cgesv_nopiv_batched(
                  magma_tally2_int_t n, magma_tally2_int_t nrhs,
                  magma_tally2FloatComplex **dA_array, magma_tally2_int_t ldda,
                  magma_tally2FloatComplex **dB_array, magma_tally2_int_t lddb,
                  magma_tally2_int_t *info_array,
                  magma_tally2_int_t batchCount, magma_tally2_queue_t queue);

magma_tally2_int_t 
magma_tally2_cgesv_rbt_batched(
                  magma_tally2_int_t n, magma_tally2_int_t nrhs,
                  magma_tally2FloatComplex **dA_array, magma_tally2_int_t ldda,
                  magma_tally2FloatComplex **dB_array, magma_tally2_int_t lddb,
                  magma_tally2_int_t *info_array,
                  magma_tally2_int_t batchCount, magma_tally2_queue_t queue);

magma_tally2_int_t 
magma_tally2_cgesv_batched(
                  magma_tally2_int_t n, magma_tally2_int_t nrhs,
                  magma_tally2FloatComplex **dA_array, magma_tally2_int_t ldda,
                  magma_tally2_int_t **dipiv_array, 
                  magma_tally2FloatComplex **dB_array, magma_tally2_int_t lddb,
                  magma_tally2_int_t *dinfo_array,
                  magma_tally2_int_t batchCount, magma_tally2_queue_t queue);

magma_tally2_int_t
magma_tally2_cgerbt_batched(
    magma_tally2_bool_t gen, magma_tally2_int_t n, magma_tally2_int_t nrhs,
    magma_tally2FloatComplex **dA_array, magma_tally2_int_t ldda,
    magma_tally2FloatComplex **dB_array, magma_tally2_int_t lddb,
    magma_tally2FloatComplex *U, magma_tally2FloatComplex *V,
    magma_tally2_int_t *info, magma_tally2_int_t batchCount, magma_tally2_queue_t queue);

void 
magma_tally2blas_cprbt_batched(
    magma_tally2_int_t n, 
    magma_tally2FloatComplex **dA_array, magma_tally2_int_t ldda, 
    magma_tally2FloatComplex *du, magma_tally2FloatComplex *dv,
    magma_tally2_int_t batchCount, magma_tally2_queue_t queue);

void
magma_tally2blas_cprbt_mv_batched(
    magma_tally2_int_t n, 
    magma_tally2FloatComplex *dv, magma_tally2FloatComplex **db_array, magma_tally2_int_t batchCount, magma_tally2_queue_t queue);


void
magma_tally2blas_cprbt_mtv_batched(
    magma_tally2_int_t n, 
    magma_tally2FloatComplex *du, magma_tally2FloatComplex **db_array, magma_tally2_int_t batchCount, magma_tally2_queue_t queue);





void 
magma_tally2_clacgv_batched(magma_tally2_int_t n, magma_tally2FloatComplex **x_array, magma_tally2_int_t incx, int offset, int batchCount, magma_tally2_queue_t queue);

void 
magma_tally2_cpotf2_csscal_batched(magma_tally2_int_t n, magma_tally2FloatComplex **x_array, magma_tally2_int_t incx, magma_tally2_int_t offset, magma_tally2_int_t *info_array, magma_tally2_int_t batchCount, magma_tally2_queue_t queue);

void 
magma_tally2_cpotf2_cdotc_batched(magma_tally2_int_t n, magma_tally2FloatComplex **x_array, magma_tally2_int_t incx, magma_tally2_int_t offset, magma_tally2_int_t *info_array, magma_tally2_int_t gbstep, magma_tally2_int_t batchCount, magma_tally2_queue_t queue);


void 
setup_pivinfo( magma_tally2_int_t *pivinfo, magma_tally2_int_t *ipiv, 
                      magma_tally2_int_t m, magma_tally2_int_t nb, 
                      magma_tally2_queue_t queue);


void
magma_tally2blas_cgeadd_batched_q(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2FloatComplex alpha,
    magma_tally2FloatComplex_const_ptr  const dAarray[], magma_tally2_int_t ldda,
    magma_tally2FloatComplex_ptr              dBarray[], magma_tally2_int_t lddb,
    magma_tally2_int_t batchCount, magma_tally2_queue_t queue );

void
magma_tally2blas_clacpy_batched(
    magma_tally2_uplo_t uplo, magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2FloatComplex_const_ptr  const dAarray[], magma_tally2_int_t ldda,
    magma_tally2FloatComplex_ptr              dBarray[], magma_tally2_int_t lddb,
    magma_tally2_int_t batchCount, magma_tally2_queue_t queue );

void
magma_tally2blas_cgeadd_batched(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2FloatComplex alpha,
    magma_tally2FloatComplex_const_ptr  const dAarray[], magma_tally2_int_t ldda,
    magma_tally2FloatComplex_ptr              dBarray[], magma_tally2_int_t lddb,
    magma_tally2_int_t batchCount );


void
magma_tally2blas_cgemv_batched(
    magma_tally2_trans_t trans, magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2FloatComplex alpha,
    magma_tally2FloatComplex_ptr dA_array[], magma_tally2_int_t ldda,
    magma_tally2FloatComplex_ptr dx_array[], magma_tally2_int_t incx,
    magma_tally2FloatComplex beta,
    magma_tally2FloatComplex_ptr dy_array[], magma_tally2_int_t incy,
    magma_tally2_int_t batchCount, magma_tally2_queue_t queue);


magma_tally2_int_t 
magma_tally2_cgeqrf_batched(
    magma_tally2_int_t m, magma_tally2_int_t n, 
    magma_tally2FloatComplex **dA_array,
    magma_tally2_int_t lda, 
    magma_tally2FloatComplex **tau_array,
    magma_tally2_int_t *info_array, magma_tally2_int_t batchCount, magma_tally2_queue_t queue);

magma_tally2_int_t 
magma_tally2_cgeqrf_batched_v4(
    magma_tally2_int_t m, magma_tally2_int_t n, 
    magma_tally2FloatComplex **dA_array,
    magma_tally2_int_t lda, 
    magma_tally2FloatComplex **tau_array,
    magma_tally2_int_t *info_array, magma_tally2_int_t batchCount);

magma_tally2_int_t
magma_tally2_cgeqrf_panel_batched(
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t nb,    
    magma_tally2FloatComplex** dA_array,    magma_tally2_int_t ldda,
    magma_tally2FloatComplex** tau_array, 
    magma_tally2FloatComplex** dT_array, magma_tally2_int_t ldt, 
    magma_tally2FloatComplex** dR_array, magma_tally2_int_t ldr,
    magma_tally2FloatComplex** dW0_displ, 
    magma_tally2FloatComplex** dW1_displ,
    magma_tally2FloatComplex *dwork,  
    magma_tally2FloatComplex** W_array, 
    magma_tally2FloatComplex** W2_array,
    magma_tally2_int_t *info_array,
    magma_tally2_int_t batchCount, cublasHandle_t myhandle, magma_tally2_queue_t queue);

magma_tally2_int_t
magma_tally2_cgeqrf_panel_batched_v4(
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t nb,    
    magma_tally2FloatComplex** dA_array,    magma_tally2_int_t ldda,
    magma_tally2FloatComplex** tau_array, 
    magma_tally2FloatComplex** dT_array, magma_tally2_int_t ldt, 
    magma_tally2FloatComplex** dR_array, magma_tally2_int_t ldr,
    float** dnorm_array,  
    magma_tally2FloatComplex** dW0_displ, 
    magma_tally2FloatComplex** dW1_displ,
    magma_tally2FloatComplex *dwork,  
    magma_tally2FloatComplex** W_array, 
    magma_tally2FloatComplex** W2_array,
    magma_tally2_int_t *info_array,
    magma_tally2_int_t batchCount, cublasHandle_t myhandle);

magma_tally2_int_t
magma_tally2_cgeqr2x_batched_v4(magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2FloatComplex **dA_array,
                  magma_tally2_int_t lda, 
                  magma_tally2FloatComplex **tau_array,
                  magma_tally2FloatComplex **dR_array, magma_tally2_int_t ldr,
                  float **dwork_array,  
                  magma_tally2_int_t *info, magma_tally2_int_t batchCount);

magma_tally2_int_t
magma_tally2_cgeqr2_batched(magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2FloatComplex **dA_array,
                  magma_tally2_int_t lda, 
                  magma_tally2FloatComplex **tau_array,
                  magma_tally2_int_t *info, magma_tally2_int_t batchCount, magma_tally2_queue_t queue);

magma_tally2_int_t
magma_tally2_clarfb_cgemm_batched(
                  cublasHandle_t myhandle,
                  magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
                  magma_tally2FloatComplex **dV_array,    magma_tally2_int_t ldv,
                  magma_tally2FloatComplex **dT_array,    magma_tally2_int_t ldt,
                  magma_tally2FloatComplex **dA_array,    magma_tally2_int_t lda,
                  magma_tally2FloatComplex **W_array,     magma_tally2_int_t ldw,
                  magma_tally2FloatComplex **W2_array,    magma_tally2_int_t ldw2,
                  magma_tally2_int_t batchCount, magma_tally2_queue_t queue);

magma_tally2_int_t
magma_tally2_clarfb_gemm_batched(
    magma_tally2_side_t side, magma_tally2_trans_t trans, magma_tally2_direct_t direct, magma_tally2_storev_t storev,
    magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k,
    magma_tally2FloatComplex_const_ptr dV_array[],    magma_tally2_int_t lddv,
    magma_tally2FloatComplex_const_ptr dT_array[],    magma_tally2_int_t lddt,
    magma_tally2FloatComplex_ptr dC_array[],          magma_tally2_int_t lddc,
    magma_tally2FloatComplex_ptr dwork_array[],       magma_tally2_int_t ldwork,
    magma_tally2FloatComplex_ptr dworkvt_array[],     magma_tally2_int_t ldworkvt,
    magma_tally2_int_t batchCount, cublasHandle_t myhandle, magma_tally2_queue_t queue);




void
magma_tally2_clarft_batched_vold(magma_tally2_int_t n, magma_tally2_int_t k, magma_tally2FloatComplex **v_array, magma_tally2_int_t ldv,
                    magma_tally2FloatComplex **tau_array,
                    magma_tally2FloatComplex **T_array, magma_tally2_int_t ldt, 
                    magma_tally2_int_t batchCount);





magma_tally2_int_t
magma_tally2_clarft_batched(magma_tally2_int_t n, magma_tally2_int_t k, magma_tally2_int_t stair_T, 
                magma_tally2FloatComplex **v_array, magma_tally2_int_t ldv,
                magma_tally2FloatComplex **tau_array, magma_tally2FloatComplex **T_array, magma_tally2_int_t ldt, 
                magma_tally2FloatComplex **work_array, magma_tally2_int_t lwork, magma_tally2_int_t batchCount, cublasHandle_t myhandle, magma_tally2_queue_t queue);

void
magma_tally2_clarft_sm32x32_batched(magma_tally2_int_t n, magma_tally2_int_t k, magma_tally2FloatComplex **v_array, magma_tally2_int_t ldv,
                    magma_tally2FloatComplex **tau_array, magma_tally2FloatComplex **T_array, magma_tally2_int_t ldt, 
                    magma_tally2_int_t batchCount, cublasHandle_t myhandle, magma_tally2_queue_t queue);



void magma_tally2blas_clarft_recctrmv_sm32x32(
    magma_tally2_int_t m, magma_tally2_int_t n, 
    magma_tally2FloatComplex *tau, 
    magma_tally2FloatComplex *Trec, magma_tally2_int_t ldtrec, 
    magma_tally2FloatComplex *Ttri, magma_tally2_int_t ldttri);


void magma_tally2blas_clarft_recctrmv_sm32x32_batched(
    magma_tally2_int_t m, magma_tally2_int_t n, 
    magma_tally2FloatComplex **tau_array, 
    magma_tally2FloatComplex **Trec_array, magma_tally2_int_t ldtrec, 
    magma_tally2FloatComplex **Ttri_array, magma_tally2_int_t ldttri,
    magma_tally2_int_t batchCount, magma_tally2_queue_t queue);

void magma_tally2blas_clarft_ctrmv_sm32x32(
    magma_tally2_int_t m, magma_tally2_int_t n, 
    magma_tally2FloatComplex *tau, 
    magma_tally2FloatComplex *Tin, magma_tally2_int_t ldtin, 
    magma_tally2FloatComplex *Tout, magma_tally2_int_t ldtout);

void magma_tally2blas_clarft_ctrmv_sm32x32_batched(
    magma_tally2_int_t m, magma_tally2_int_t n, 
    magma_tally2FloatComplex **tau_array, 
    magma_tally2FloatComplex **Tin_array, magma_tally2_int_t ldtin, 
    magma_tally2FloatComplex **Tout_array, magma_tally2_int_t ldtout,
    magma_tally2_int_t batchCount, magma_tally2_queue_t queue);

void magma_tally2blas_clarft_gemv_loop_inside(
    int n, int k, 
    magma_tally2FloatComplex *tau, 
    magma_tally2FloatComplex *v, int ldv, 
    magma_tally2FloatComplex *T, int ldt);

void magma_tally2blas_clarft_gemv_loop_inside_batched(
    int n, int k, 
    magma_tally2FloatComplex **tau_array, 
    magma_tally2FloatComplex **v_array, int ldv, 
    magma_tally2FloatComplex **T_array, int ldt, magma_tally2_int_t batchCount, magma_tally2_queue_t queue);

void magma_tally2blas_clarft_gemvrowwise(
    magma_tally2_int_t m, magma_tally2_int_t i, 
    magma_tally2FloatComplex *tau, 
    magma_tally2FloatComplex *v, magma_tally2_int_t ldv, 
    magma_tally2FloatComplex *T, magma_tally2_int_t ldt,
    magma_tally2FloatComplex *W);

void magma_tally2blas_clarft_gemvrowwise_batched(
    magma_tally2_int_t m, magma_tally2_int_t i, 
    magma_tally2FloatComplex **tau_array, 
    magma_tally2FloatComplex **v_array, magma_tally2_int_t ldv, 
    magma_tally2FloatComplex **T_array, magma_tally2_int_t ldt,
    magma_tally2_int_t batchCount, magma_tally2_queue_t queue);

void magma_tally2blas_clarft_gemvcolwise(
    magma_tally2_int_t m,  magma_tally2_int_t step,
    magma_tally2FloatComplex *v, magma_tally2_int_t ldv, 
    magma_tally2FloatComplex *T,  magma_tally2_int_t ldt,
    magma_tally2FloatComplex *tau);

void magma_tally2blas_clarft_gemvcolwise_batched(
    magma_tally2_int_t m,  magma_tally2_int_t step,
    magma_tally2FloatComplex **v_array, magma_tally2_int_t ldv, 
    magma_tally2FloatComplex **T_array,  magma_tally2_int_t ldt,
    magma_tally2FloatComplex **tau_array, magma_tally2_int_t batchCount, magma_tally2_queue_t queue);




void cgeqrf_copy_upper_batched(                
                  magma_tally2_int_t n, magma_tally2_int_t nb,
                  magma_tally2FloatComplex **dV_array,    magma_tally2_int_t ldv,
                  magma_tally2FloatComplex **dR_array,    magma_tally2_int_t ldr,
          magma_tally2_int_t batchCount, magma_tally2_queue_t queue);



void 
magma_tally2blas_scnrm2_cols_batched(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2FloatComplex **dA_array, magma_tally2_int_t lda, 
    float **dxnorm_array, magma_tally2_int_t batchCount);
 
void 
magma_tally2_clarfgx_batched(magma_tally2_int_t n, magma_tally2FloatComplex **dx0_array, magma_tally2FloatComplex **dx_array, 
                  magma_tally2FloatComplex **dtau_array, float **dxnorm_array, 
                  magma_tally2FloatComplex **dR_array, magma_tally2_int_t it, magma_tally2_int_t batchCount);


void 
magma_tally2_clarfx_batched_v4(magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2FloatComplex **v_array, magma_tally2FloatComplex **tau_array,
                magma_tally2FloatComplex **C_array, magma_tally2_int_t ldc, float **xnorm_array, 
                magma_tally2_int_t step, 
                magma_tally2_int_t batchCount);


void 
magma_tally2blas_clarfg_batched(
    magma_tally2_int_t n,
    magma_tally2FloatComplex** dalpha_array, magma_tally2FloatComplex** dx_array, magma_tally2_int_t incx,
    magma_tally2FloatComplex** dtau_array, magma_tally2_int_t batchCount );





// for debugging purpose
void 
cset_stepinit_ipiv(magma_tally2_int_t **ipiv_array,
                 magma_tally2_int_t pm,
                 magma_tally2_int_t batchCount);



#ifdef __cplusplus
}
#endif

#undef COMPLEX

#endif  /* MAGMA_tally2_CBATCHED_H */
