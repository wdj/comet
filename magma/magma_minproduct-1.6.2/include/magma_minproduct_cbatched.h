/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Azzam Haidar
       @author Tingxing Dong

       @generated from magma_minproduct_zbatched.h normal z -> c, Fri Jan 30 19:00:06 2015
*/

#ifndef MAGMA_minproduct_CBATCHED_H
#define MAGMA_minproduct_CBATCHED_H

#include "magma_minproduct_types.h"

#define COMPLEX

#ifdef __cplusplus
extern "C" {
#endif
  /*
   *  local auxiliary routines
   */
void 
cset_pointer_int(magma_minproduct_int_t **output_array,
        magma_minproduct_int_t *input,
        magma_minproduct_int_t lda,
        magma_minproduct_int_t row, magma_minproduct_int_t column, 
        magma_minproduct_int_t batchSize,
        magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

void 
cset_pointer(magma_minproductFloatComplex **output_array,
                 magma_minproductFloatComplex *input,
                 magma_minproduct_int_t lda,
                 magma_minproduct_int_t row, magma_minproduct_int_t column,
                 magma_minproduct_int_t batchSize,
                 magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);


void 
cset_array(magma_minproductFloatComplex **output_array,
               magma_minproductFloatComplex **input_array, magma_minproduct_int_t lda,
               magma_minproduct_int_t row, magma_minproduct_int_t column,
               magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

void 
magma_minproduct_cdisplace_pointers(magma_minproductFloatComplex **output_array,
               magma_minproductFloatComplex **input_array, magma_minproduct_int_t lda,
               magma_minproduct_int_t row, magma_minproduct_int_t column, 
               magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

  /*
   *  LAPACK batched routines
   */

  /*
   *  BLAS batched routines
   */

void
magma_minproductblas_cgemm_batched(
    magma_minproduct_trans_t transA, magma_minproduct_trans_t transB, magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductFloatComplex alpha,
    magma_minproductFloatComplex const * const * dA_array, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex const * const * dB_array, magma_minproduct_int_t lddb,
    magma_minproductFloatComplex beta,
    magma_minproductFloatComplex **dC_array, magma_minproduct_int_t lddc, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue );

void
magma_minproductblas_cgemm_batched_lg(
    magma_minproduct_trans_t transA, magma_minproduct_trans_t transB, magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductFloatComplex alpha,
    magma_minproductFloatComplex const * const * dA_array, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex const * const * dB_array, magma_minproduct_int_t lddb,
    magma_minproductFloatComplex beta,
    magma_minproductFloatComplex **dC_array, magma_minproduct_int_t lddc, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue );
void
magma_minproductblas_cgemm_batched_k32(
    magma_minproduct_trans_t transA, magma_minproduct_trans_t transB, magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductFloatComplex alpha,
    magma_minproductFloatComplex const * const * dA_array, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex const * const * dB_array, magma_minproduct_int_t lddb,
    magma_minproductFloatComplex beta,
    magma_minproductFloatComplex **dC_array, magma_minproduct_int_t lddc, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue );



void 
magma_minproductblas_cherk_NC_batched( magma_minproduct_trans_t TRANSA, magma_minproduct_trans_t TRANSB, int m , int n , int k , 
                       magma_minproductFloatComplex alpha, magma_minproductFloatComplex **dA_array, int lda, 
                       magma_minproductFloatComplex **B_array, int ldb, 
                       magma_minproductFloatComplex beta,        magma_minproductFloatComplex **C_array, int ldc, 
                       magma_minproduct_int_t batchCount);

void
magma_minproductblas_cherk_batched(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans, magma_minproduct_int_t n, magma_minproduct_int_t k,
    float alpha,
    magma_minproductFloatComplex const * const * dA_array, magma_minproduct_int_t ldda,
    float beta,
    magma_minproductFloatComplex **dC_array, magma_minproduct_int_t lddc, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue );

void
magma_minproductblas_cherk_batched_lg(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans, magma_minproduct_int_t n, magma_minproduct_int_t k,
    float alpha,
    magma_minproductFloatComplex const * const * dA_array, magma_minproduct_int_t ldda,
    float beta,
    magma_minproductFloatComplex **dC_array, magma_minproduct_int_t lddc, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue );

void
magma_minproductblas_cherk_batched_k32(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans, magma_minproduct_int_t n, magma_minproduct_int_t k,
    float alpha,
    magma_minproductFloatComplex const * const * dA_array, magma_minproduct_int_t ldda,
    float beta,
    magma_minproductFloatComplex **dC_array, magma_minproduct_int_t lddc, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue );






magma_minproduct_int_t 
magma_minproduct_cpotf2_tile_batched(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex **dA_array, magma_minproduct_int_t lda,
    magma_minproduct_int_t *info_array, magma_minproduct_int_t gbstep, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_cpotrf_panel(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t nb,     
    magma_minproductFloatComplex *A, magma_minproduct_int_t lda,
    magma_minproductFloatComplex* dinvA, magma_minproductFloatComplex** dinvA_array, magma_minproduct_int_t invA_msize,
    magma_minproductFloatComplex* x, magma_minproductFloatComplex** x_array,  magma_minproduct_int_t x_msize,
    magma_minproduct_int_t *info_array,  magma_minproduct_int_t batchCount, magma_minproduct_int_t matrixSize, magma_minproduct_queue_t queue);


void 
magma_minproductblas_ctrtri_diag_batched(
    magma_minproduct_uplo_t uplo, magma_minproduct_diag_t diag, magma_minproduct_int_t n,
    magma_minproductFloatComplex const * const *dA_array, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex **dinvA_array, 
    magma_minproduct_int_t resetozero, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);



void 
magma_minproductblas_ctrsm_batched(
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo, magma_minproduct_trans_t transA, magma_minproduct_diag_t diag,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex alpha,
    magma_minproductFloatComplex** dA_array,    magma_minproduct_int_t ldda,
    magma_minproductFloatComplex** dB_array,    magma_minproduct_int_t lddb,
    magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);


void 
magma_minproductblas_ctrsm_work_batched(
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo, magma_minproduct_trans_t transA, magma_minproduct_diag_t diag,
    magma_minproduct_int_t flag, magma_minproduct_int_t m, magma_minproduct_int_t n, 
    magma_minproductFloatComplex alpha, 
    magma_minproductFloatComplex** dA_array,    magma_minproduct_int_t ldda,
    magma_minproductFloatComplex** dB_array,    magma_minproduct_int_t lddb,
    magma_minproductFloatComplex** dX_array,    magma_minproduct_int_t lddx, 
    magma_minproductFloatComplex** dinvA_array, magma_minproduct_int_t dinvA_length,
    magma_minproductFloatComplex** dA_displ, magma_minproductFloatComplex** dB_displ, 
    magma_minproductFloatComplex** dX_displ, magma_minproductFloatComplex** dinvA_displ,
    magma_minproduct_int_t resetozero, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

void 
magma_minproductblas_ctrsm_outofplace_batched(
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo, magma_minproduct_trans_t transA, magma_minproduct_diag_t diag,
    magma_minproduct_int_t flag, magma_minproduct_int_t m, magma_minproduct_int_t n, 
    magma_minproductFloatComplex alpha, 
    magma_minproductFloatComplex** dA_array,    magma_minproduct_int_t ldda,
    magma_minproductFloatComplex** dB_array,    magma_minproduct_int_t lddb,
    magma_minproductFloatComplex** dX_array,    magma_minproduct_int_t lddx, 
    magma_minproductFloatComplex** dinvA_array, magma_minproduct_int_t dinvA_length,
    magma_minproductFloatComplex** dA_displ, magma_minproductFloatComplex** dB_displ, 
    magma_minproductFloatComplex** dX_displ, magma_minproductFloatComplex** dinvA_displ,
    magma_minproduct_int_t resetozero, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_cpotrf_batched(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductFloatComplex **dA_array, magma_minproduct_int_t lda,
    magma_minproduct_int_t *info_array,  magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);


magma_minproduct_int_t 
magma_minproduct_cpotf2_batched(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex **dA_array, magma_minproduct_int_t lda,
    magma_minproductFloatComplex **dA_displ, 
    magma_minproductFloatComplex **dW_displ,
    magma_minproductFloatComplex **dB_displ, 
    magma_minproductFloatComplex **dC_displ, 
    magma_minproduct_int_t *info_array, magma_minproduct_int_t gbstep, 
    magma_minproduct_int_t batchCount, cublasHandle_t myhandle, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_cpotrf_panel_batched(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t nb,     
    magma_minproductFloatComplex** dA_array,    magma_minproduct_int_t ldda,
    magma_minproductFloatComplex** dX_array,    magma_minproduct_int_t dX_length,
    magma_minproductFloatComplex** dinvA_array, magma_minproduct_int_t dinvA_length,
    magma_minproductFloatComplex** dW0_displ, magma_minproductFloatComplex** dW1_displ, 
    magma_minproductFloatComplex** dW2_displ, magma_minproductFloatComplex** dW3_displ,
    magma_minproductFloatComplex** dW4_displ, 
    magma_minproduct_int_t *info_array, magma_minproduct_int_t gbstep,
    magma_minproduct_int_t batchCount, cublasHandle_t myhandle, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_cpotrf_recpanel_batched(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n, 
    magma_minproduct_int_t min_recpnb,    
    magma_minproductFloatComplex** dA_array,    magma_minproduct_int_t ldda,
    magma_minproductFloatComplex** dX_array,    magma_minproduct_int_t dX_length,
    magma_minproductFloatComplex** dinvA_array, magma_minproduct_int_t dinvA_length,
    magma_minproductFloatComplex** dW0_displ, magma_minproductFloatComplex** dW1_displ,  
    magma_minproductFloatComplex** dW2_displ, magma_minproductFloatComplex** dW3_displ,
    magma_minproductFloatComplex** dW4_displ,
    magma_minproduct_int_t *info_array, magma_minproduct_int_t gbstep, 
    magma_minproduct_int_t batchCount, cublasHandle_t myhandle, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_cpotrf_rectile_batched(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n, 
    magma_minproduct_int_t min_recpnb,    
    magma_minproductFloatComplex** dA_array,    magma_minproduct_int_t ldda,
    magma_minproductFloatComplex** dX_array,    magma_minproduct_int_t dX_length,
    magma_minproductFloatComplex** dinvA_array, magma_minproduct_int_t dinvA_length,
    magma_minproductFloatComplex** dW0_displ, magma_minproductFloatComplex** dW1_displ,  
    magma_minproductFloatComplex** dW2_displ, magma_minproductFloatComplex** dW3_displ,
    magma_minproductFloatComplex** dW4_displ,
    magma_minproduct_int_t *info_array, magma_minproduct_int_t gbstep,
    magma_minproduct_int_t batchCount, cublasHandle_t myhandle, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_cpotrs_batched(
                  magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
                  magma_minproductFloatComplex **dA_array, magma_minproduct_int_t ldda,
                  magma_minproductFloatComplex **dB_array, magma_minproduct_int_t lddb,
                  magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_cposv_batched(
                  magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
                  magma_minproductFloatComplex **dA_array, magma_minproduct_int_t ldda,
                  magma_minproductFloatComplex **dB_array, magma_minproduct_int_t lddb,
                  magma_minproduct_int_t *dinfo_array,
                  magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_cgetrs_batched(
                  magma_minproduct_trans_t trans, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
                  magma_minproductFloatComplex **dA_array, magma_minproduct_int_t ldda,
                  magma_minproduct_int_t **dipiv_array, 
                  magma_minproductFloatComplex **dB_array, magma_minproduct_int_t lddb,
                  magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);



void 
magma_minproduct_claswp_rowparallel_batched( magma_minproduct_int_t n, magma_minproductFloatComplex** input_array, magma_minproduct_int_t ldi,
                   magma_minproductFloatComplex** output_array, magma_minproduct_int_t ldo,
                   magma_minproduct_int_t k1, magma_minproduct_int_t k2,
                   magma_minproduct_int_t **pivinfo_array, 
                   magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue );

void 
magma_minproduct_claswp_rowserial_batched(magma_minproduct_int_t n, magma_minproductFloatComplex** dA_array, magma_minproduct_int_t lda,
                   magma_minproduct_int_t k1, magma_minproduct_int_t k2,
                   magma_minproduct_int_t **ipiv_array, 
                   magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

void 
magma_minproduct_claswp_columnserial_batched(magma_minproduct_int_t n, magma_minproductFloatComplex** dA_array, magma_minproduct_int_t lda,
                   magma_minproduct_int_t k1, magma_minproduct_int_t k2,
                   magma_minproduct_int_t **ipiv_array, 
                   magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

void 
magma_minproductblas_ctranspose_batched(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex **dA_array,  magma_minproduct_int_t ldda,
    magma_minproductFloatComplex **dAT_array, magma_minproduct_int_t lddat, magma_minproduct_int_t batchCount );


void 
magma_minproductblas_claset_batched(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex offdiag, magma_minproductFloatComplex diag,
    magma_minproductFloatComplex_ptr dAarray[], magma_minproduct_int_t ldda,
    magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

void 
magma_minproductblas_cmemset_batched(magma_minproduct_int_t length, 
        magma_minproductFloatComplex_ptr dAarray[], magma_minproductFloatComplex val, 
        magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_cgetf2_batched(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex **dA_array, magma_minproduct_int_t lda,
    magma_minproductFloatComplex **GERA_array,
    magma_minproductFloatComplex **GERB_array,
    magma_minproductFloatComplex **GERC_array,
    magma_minproduct_int_t **ipiv_array,
    magma_minproduct_int_t *info_array, 
    magma_minproduct_int_t gbstep,            
    magma_minproduct_int_t batchCount,
    cublasHandle_t myhandle, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_cgetrf_recpanel_batched(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t min_recpnb,    
    magma_minproductFloatComplex** dA_array,    magma_minproduct_int_t ldda,
    magma_minproduct_int_t** dipiv_array, magma_minproduct_int_t** dpivinfo_array,
    magma_minproductFloatComplex** dX_array,    magma_minproduct_int_t dX_length,
    magma_minproductFloatComplex** dinvA_array, magma_minproduct_int_t dinvA_length,
    magma_minproductFloatComplex** dW1_displ, magma_minproductFloatComplex** dW2_displ,  
    magma_minproductFloatComplex** dW3_displ, magma_minproductFloatComplex** dW4_displ,
    magma_minproductFloatComplex** dW5_displ,
    magma_minproduct_int_t *info_array, magma_minproduct_int_t gbstep,
    magma_minproduct_int_t batchCount, cublasHandle_t myhandle, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_cgetrf_batched(
        magma_minproduct_int_t m, magma_minproduct_int_t n,
        magma_minproductFloatComplex **dA_array, 
        magma_minproduct_int_t lda,
        magma_minproduct_int_t **ipiv_array, 
        magma_minproduct_int_t *info_array, 
        magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_cgetri_outofplace_batched( magma_minproduct_int_t n, 
                  magma_minproductFloatComplex **dA_array, magma_minproduct_int_t ldda,
                  magma_minproduct_int_t **dipiv_array, 
                  magma_minproductFloatComplex **dinvA_array, magma_minproduct_int_t lddia,
                  magma_minproduct_int_t *info_array,
                  magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

void 
magma_minproduct_cdisplace_intpointers(magma_minproduct_int_t **output_array,
               magma_minproduct_int_t **input_array, magma_minproduct_int_t lda,
               magma_minproduct_int_t row, magma_minproduct_int_t column, 
               magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);





void 
magma_minproductblas_icamax_atomic_batched(magma_minproduct_int_t n, magma_minproductFloatComplex** x_array, magma_minproduct_int_t incx, magma_minproduct_int_t **max_id_array, magma_minproduct_int_t batchCount);

void 
magma_minproductblas_icamax_tree_batched(magma_minproduct_int_t n, magma_minproductFloatComplex** x_array, magma_minproduct_int_t incx, magma_minproduct_int_t **max_id_array, magma_minproduct_int_t batchCount);



void 
magma_minproductblas_icamax_batched(magma_minproduct_int_t n, magma_minproductFloatComplex** x_array, magma_minproduct_int_t incx, magma_minproduct_int_t **max_id_array, magma_minproduct_int_t batchCount);

void 
magma_minproductblas_icamax(magma_minproduct_int_t n, magma_minproductFloatComplex* x, magma_minproduct_int_t incx, magma_minproduct_int_t *max_id);


magma_minproduct_int_t 
magma_minproduct_icamax_batched(magma_minproduct_int_t length, 
        magma_minproductFloatComplex **x_array, magma_minproduct_int_t incx, magma_minproduct_int_t step,  magma_minproduct_int_t lda,
        magma_minproduct_int_t** ipiv_array, magma_minproduct_int_t *info_array, magma_minproduct_int_t gbstep, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_cswap_batched(magma_minproduct_int_t n, magma_minproductFloatComplex **x_array, magma_minproduct_int_t incx, magma_minproduct_int_t j, 
                 magma_minproduct_int_t** ipiv_array, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_cscal_cgeru_batched(magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t step,
                                      magma_minproductFloatComplex **dA_array, magma_minproduct_int_t lda,
                                      magma_minproduct_int_t *info_array, magma_minproduct_int_t gbstep, 
                                      magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_ccomputecolumn_batched(magma_minproduct_int_t m, magma_minproduct_int_t paneloffset, magma_minproduct_int_t step, 
                                        magma_minproductFloatComplex **dA_array,  magma_minproduct_int_t lda,
                                        magma_minproduct_int_t **ipiv_array, 
                                        magma_minproduct_int_t *info_array, magma_minproduct_int_t gbstep, 
                                        magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

void 
magma_minproduct_cgetf2trsm_batched(magma_minproduct_int_t ib, magma_minproduct_int_t n, magma_minproductFloatComplex **dA_array,  magma_minproduct_int_t j, magma_minproduct_int_t lda,
                       magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);


magma_minproduct_int_t 
magma_minproduct_cgetf2_nopiv_batched(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex **dA_array, magma_minproduct_int_t lda,
    magma_minproductFloatComplex **dW0_displ,
    magma_minproductFloatComplex **dW1_displ,
    magma_minproductFloatComplex **dW2_displ,
    magma_minproduct_int_t *info_array,            
    magma_minproduct_int_t gbstep, 
    magma_minproduct_int_t batchCount,
    cublasHandle_t myhandle, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_cgetrf_panel_nopiv_batched(
    magma_minproduct_int_t m, magma_minproduct_int_t nb,    
    magma_minproductFloatComplex** dA_array,    magma_minproduct_int_t ldda,
    magma_minproductFloatComplex** dX_array,    magma_minproduct_int_t dX_length,
    magma_minproductFloatComplex** dinvA_array, magma_minproduct_int_t dinvA_length,
    magma_minproductFloatComplex** dW0_displ, magma_minproductFloatComplex** dW1_displ,  
    magma_minproductFloatComplex** dW2_displ, magma_minproductFloatComplex** dW3_displ,
    magma_minproductFloatComplex** dW4_displ,     
    magma_minproduct_int_t *info_array, magma_minproduct_int_t gbstep, 
    magma_minproduct_int_t batchCount, cublasHandle_t myhandle, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_cgetrf_recpanel_nopiv_batched(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t min_recpnb,    
    magma_minproductFloatComplex** dA_array,    magma_minproduct_int_t ldda,
    magma_minproductFloatComplex** dX_array,    magma_minproduct_int_t dX_length,
    magma_minproductFloatComplex** dinvA_array, magma_minproduct_int_t dinvA_length,
    magma_minproductFloatComplex** dW1_displ, magma_minproductFloatComplex** dW2_displ,  
    magma_minproductFloatComplex** dW3_displ, magma_minproductFloatComplex** dW4_displ,
    magma_minproductFloatComplex** dW5_displ, 
    magma_minproduct_int_t *info_array, magma_minproduct_int_t gbstep,
    magma_minproduct_int_t batchCount, cublasHandle_t myhandle, magma_minproduct_queue_t queue);


magma_minproduct_int_t 
magma_minproduct_cgetrf_nopiv_batched(
        magma_minproduct_int_t m, magma_minproduct_int_t n,
        magma_minproductFloatComplex **dA_array, 
        magma_minproduct_int_t lda,
        magma_minproduct_int_t *info_array, 
        magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_cgetrs_nopiv_batched(
                  magma_minproduct_trans_t trans, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
                  magma_minproductFloatComplex **dA_array, magma_minproduct_int_t ldda,
                  magma_minproductFloatComplex **dB_array, magma_minproduct_int_t lddb,
                  magma_minproduct_int_t *info_array,
                  magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_cgesv_nopiv_batched(
                  magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
                  magma_minproductFloatComplex **dA_array, magma_minproduct_int_t ldda,
                  magma_minproductFloatComplex **dB_array, magma_minproduct_int_t lddb,
                  magma_minproduct_int_t *info_array,
                  magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_cgesv_rbt_batched(
                  magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
                  magma_minproductFloatComplex **dA_array, magma_minproduct_int_t ldda,
                  magma_minproductFloatComplex **dB_array, magma_minproduct_int_t lddb,
                  magma_minproduct_int_t *info_array,
                  magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_cgesv_batched(
                  magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
                  magma_minproductFloatComplex **dA_array, magma_minproduct_int_t ldda,
                  magma_minproduct_int_t **dipiv_array, 
                  magma_minproductFloatComplex **dB_array, magma_minproduct_int_t lddb,
                  magma_minproduct_int_t *dinfo_array,
                  magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

magma_minproduct_int_t
magma_minproduct_cgerbt_batched(
    magma_minproduct_bool_t gen, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
    magma_minproductFloatComplex **dA_array, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex **dB_array, magma_minproduct_int_t lddb,
    magma_minproductFloatComplex *U, magma_minproductFloatComplex *V,
    magma_minproduct_int_t *info, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

void 
magma_minproductblas_cprbt_batched(
    magma_minproduct_int_t n, 
    magma_minproductFloatComplex **dA_array, magma_minproduct_int_t ldda, 
    magma_minproductFloatComplex *du, magma_minproductFloatComplex *dv,
    magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

void
magma_minproductblas_cprbt_mv_batched(
    magma_minproduct_int_t n, 
    magma_minproductFloatComplex *dv, magma_minproductFloatComplex **db_array, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);


void
magma_minproductblas_cprbt_mtv_batched(
    magma_minproduct_int_t n, 
    magma_minproductFloatComplex *du, magma_minproductFloatComplex **db_array, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);





void 
magma_minproduct_clacgv_batched(magma_minproduct_int_t n, magma_minproductFloatComplex **x_array, magma_minproduct_int_t incx, int offset, int batchCount, magma_minproduct_queue_t queue);

void 
magma_minproduct_cpotf2_csscal_batched(magma_minproduct_int_t n, magma_minproductFloatComplex **x_array, magma_minproduct_int_t incx, magma_minproduct_int_t offset, magma_minproduct_int_t *info_array, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

void 
magma_minproduct_cpotf2_cdotc_batched(magma_minproduct_int_t n, magma_minproductFloatComplex **x_array, magma_minproduct_int_t incx, magma_minproduct_int_t offset, magma_minproduct_int_t *info_array, magma_minproduct_int_t gbstep, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);


void 
setup_pivinfo( magma_minproduct_int_t *pivinfo, magma_minproduct_int_t *ipiv, 
                      magma_minproduct_int_t m, magma_minproduct_int_t nb, 
                      magma_minproduct_queue_t queue);


void
magma_minproductblas_cgeadd_batched_q(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex alpha,
    magma_minproductFloatComplex_const_ptr  const dAarray[], magma_minproduct_int_t ldda,
    magma_minproductFloatComplex_ptr              dBarray[], magma_minproduct_int_t lddb,
    magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue );

void
magma_minproductblas_clacpy_batched(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex_const_ptr  const dAarray[], magma_minproduct_int_t ldda,
    magma_minproductFloatComplex_ptr              dBarray[], magma_minproduct_int_t lddb,
    magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue );

void
magma_minproductblas_cgeadd_batched(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex alpha,
    magma_minproductFloatComplex_const_ptr  const dAarray[], magma_minproduct_int_t ldda,
    magma_minproductFloatComplex_ptr              dBarray[], magma_minproduct_int_t lddb,
    magma_minproduct_int_t batchCount );


void
magma_minproductblas_cgemv_batched(
    magma_minproduct_trans_t trans, magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex alpha,
    magma_minproductFloatComplex_ptr dA_array[], magma_minproduct_int_t ldda,
    magma_minproductFloatComplex_ptr dx_array[], magma_minproduct_int_t incx,
    magma_minproductFloatComplex beta,
    magma_minproductFloatComplex_ptr dy_array[], magma_minproduct_int_t incy,
    magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);


magma_minproduct_int_t 
magma_minproduct_cgeqrf_batched(
    magma_minproduct_int_t m, magma_minproduct_int_t n, 
    magma_minproductFloatComplex **dA_array,
    magma_minproduct_int_t lda, 
    magma_minproductFloatComplex **tau_array,
    magma_minproduct_int_t *info_array, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_cgeqrf_batched_v4(
    magma_minproduct_int_t m, magma_minproduct_int_t n, 
    magma_minproductFloatComplex **dA_array,
    magma_minproduct_int_t lda, 
    magma_minproductFloatComplex **tau_array,
    magma_minproduct_int_t *info_array, magma_minproduct_int_t batchCount);

magma_minproduct_int_t
magma_minproduct_cgeqrf_panel_batched(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t nb,    
    magma_minproductFloatComplex** dA_array,    magma_minproduct_int_t ldda,
    magma_minproductFloatComplex** tau_array, 
    magma_minproductFloatComplex** dT_array, magma_minproduct_int_t ldt, 
    magma_minproductFloatComplex** dR_array, magma_minproduct_int_t ldr,
    magma_minproductFloatComplex** dW0_displ, 
    magma_minproductFloatComplex** dW1_displ,
    magma_minproductFloatComplex *dwork,  
    magma_minproductFloatComplex** W_array, 
    magma_minproductFloatComplex** W2_array,
    magma_minproduct_int_t *info_array,
    magma_minproduct_int_t batchCount, cublasHandle_t myhandle, magma_minproduct_queue_t queue);

magma_minproduct_int_t
magma_minproduct_cgeqrf_panel_batched_v4(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t nb,    
    magma_minproductFloatComplex** dA_array,    magma_minproduct_int_t ldda,
    magma_minproductFloatComplex** tau_array, 
    magma_minproductFloatComplex** dT_array, magma_minproduct_int_t ldt, 
    magma_minproductFloatComplex** dR_array, magma_minproduct_int_t ldr,
    float** dnorm_array,  
    magma_minproductFloatComplex** dW0_displ, 
    magma_minproductFloatComplex** dW1_displ,
    magma_minproductFloatComplex *dwork,  
    magma_minproductFloatComplex** W_array, 
    magma_minproductFloatComplex** W2_array,
    magma_minproduct_int_t *info_array,
    magma_minproduct_int_t batchCount, cublasHandle_t myhandle);

magma_minproduct_int_t
magma_minproduct_cgeqr2x_batched_v4(magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproductFloatComplex **dA_array,
                  magma_minproduct_int_t lda, 
                  magma_minproductFloatComplex **tau_array,
                  magma_minproductFloatComplex **dR_array, magma_minproduct_int_t ldr,
                  float **dwork_array,  
                  magma_minproduct_int_t *info, magma_minproduct_int_t batchCount);

magma_minproduct_int_t
magma_minproduct_cgeqr2_batched(magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproductFloatComplex **dA_array,
                  magma_minproduct_int_t lda, 
                  magma_minproductFloatComplex **tau_array,
                  magma_minproduct_int_t *info, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

magma_minproduct_int_t
magma_minproduct_clarfb_cgemm_batched(
                  cublasHandle_t myhandle,
                  magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
                  magma_minproductFloatComplex **dV_array,    magma_minproduct_int_t ldv,
                  magma_minproductFloatComplex **dT_array,    magma_minproduct_int_t ldt,
                  magma_minproductFloatComplex **dA_array,    magma_minproduct_int_t lda,
                  magma_minproductFloatComplex **W_array,     magma_minproduct_int_t ldw,
                  magma_minproductFloatComplex **W2_array,    magma_minproduct_int_t ldw2,
                  magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

magma_minproduct_int_t
magma_minproduct_clarfb_gemm_batched(
    magma_minproduct_side_t side, magma_minproduct_trans_t trans, magma_minproduct_direct_t direct, magma_minproduct_storev_t storev,
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductFloatComplex_const_ptr dV_array[],    magma_minproduct_int_t lddv,
    magma_minproductFloatComplex_const_ptr dT_array[],    magma_minproduct_int_t lddt,
    magma_minproductFloatComplex_ptr dC_array[],          magma_minproduct_int_t lddc,
    magma_minproductFloatComplex_ptr dwork_array[],       magma_minproduct_int_t ldwork,
    magma_minproductFloatComplex_ptr dworkvt_array[],     magma_minproduct_int_t ldworkvt,
    magma_minproduct_int_t batchCount, cublasHandle_t myhandle, magma_minproduct_queue_t queue);




void
magma_minproduct_clarft_batched_vold(magma_minproduct_int_t n, magma_minproduct_int_t k, magma_minproductFloatComplex **v_array, magma_minproduct_int_t ldv,
                    magma_minproductFloatComplex **tau_array,
                    magma_minproductFloatComplex **T_array, magma_minproduct_int_t ldt, 
                    magma_minproduct_int_t batchCount);





magma_minproduct_int_t
magma_minproduct_clarft_batched(magma_minproduct_int_t n, magma_minproduct_int_t k, magma_minproduct_int_t stair_T, 
                magma_minproductFloatComplex **v_array, magma_minproduct_int_t ldv,
                magma_minproductFloatComplex **tau_array, magma_minproductFloatComplex **T_array, magma_minproduct_int_t ldt, 
                magma_minproductFloatComplex **work_array, magma_minproduct_int_t lwork, magma_minproduct_int_t batchCount, cublasHandle_t myhandle, magma_minproduct_queue_t queue);

void
magma_minproduct_clarft_sm32x32_batched(magma_minproduct_int_t n, magma_minproduct_int_t k, magma_minproductFloatComplex **v_array, magma_minproduct_int_t ldv,
                    magma_minproductFloatComplex **tau_array, magma_minproductFloatComplex **T_array, magma_minproduct_int_t ldt, 
                    magma_minproduct_int_t batchCount, cublasHandle_t myhandle, magma_minproduct_queue_t queue);



void magma_minproductblas_clarft_recctrmv_sm32x32(
    magma_minproduct_int_t m, magma_minproduct_int_t n, 
    magma_minproductFloatComplex *tau, 
    magma_minproductFloatComplex *Trec, magma_minproduct_int_t ldtrec, 
    magma_minproductFloatComplex *Ttri, magma_minproduct_int_t ldttri);


void magma_minproductblas_clarft_recctrmv_sm32x32_batched(
    magma_minproduct_int_t m, magma_minproduct_int_t n, 
    magma_minproductFloatComplex **tau_array, 
    magma_minproductFloatComplex **Trec_array, magma_minproduct_int_t ldtrec, 
    magma_minproductFloatComplex **Ttri_array, magma_minproduct_int_t ldttri,
    magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

void magma_minproductblas_clarft_ctrmv_sm32x32(
    magma_minproduct_int_t m, magma_minproduct_int_t n, 
    magma_minproductFloatComplex *tau, 
    magma_minproductFloatComplex *Tin, magma_minproduct_int_t ldtin, 
    magma_minproductFloatComplex *Tout, magma_minproduct_int_t ldtout);

void magma_minproductblas_clarft_ctrmv_sm32x32_batched(
    magma_minproduct_int_t m, magma_minproduct_int_t n, 
    magma_minproductFloatComplex **tau_array, 
    magma_minproductFloatComplex **Tin_array, magma_minproduct_int_t ldtin, 
    magma_minproductFloatComplex **Tout_array, magma_minproduct_int_t ldtout,
    magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

void magma_minproductblas_clarft_gemv_loop_inside(
    int n, int k, 
    magma_minproductFloatComplex *tau, 
    magma_minproductFloatComplex *v, int ldv, 
    magma_minproductFloatComplex *T, int ldt);

void magma_minproductblas_clarft_gemv_loop_inside_batched(
    int n, int k, 
    magma_minproductFloatComplex **tau_array, 
    magma_minproductFloatComplex **v_array, int ldv, 
    magma_minproductFloatComplex **T_array, int ldt, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

void magma_minproductblas_clarft_gemvrowwise(
    magma_minproduct_int_t m, magma_minproduct_int_t i, 
    magma_minproductFloatComplex *tau, 
    magma_minproductFloatComplex *v, magma_minproduct_int_t ldv, 
    magma_minproductFloatComplex *T, magma_minproduct_int_t ldt,
    magma_minproductFloatComplex *W);

void magma_minproductblas_clarft_gemvrowwise_batched(
    magma_minproduct_int_t m, magma_minproduct_int_t i, 
    magma_minproductFloatComplex **tau_array, 
    magma_minproductFloatComplex **v_array, magma_minproduct_int_t ldv, 
    magma_minproductFloatComplex **T_array, magma_minproduct_int_t ldt,
    magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

void magma_minproductblas_clarft_gemvcolwise(
    magma_minproduct_int_t m,  magma_minproduct_int_t step,
    magma_minproductFloatComplex *v, magma_minproduct_int_t ldv, 
    magma_minproductFloatComplex *T,  magma_minproduct_int_t ldt,
    magma_minproductFloatComplex *tau);

void magma_minproductblas_clarft_gemvcolwise_batched(
    magma_minproduct_int_t m,  magma_minproduct_int_t step,
    magma_minproductFloatComplex **v_array, magma_minproduct_int_t ldv, 
    magma_minproductFloatComplex **T_array,  magma_minproduct_int_t ldt,
    magma_minproductFloatComplex **tau_array, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);




void cgeqrf_copy_upper_batched(                
                  magma_minproduct_int_t n, magma_minproduct_int_t nb,
                  magma_minproductFloatComplex **dV_array,    magma_minproduct_int_t ldv,
                  magma_minproductFloatComplex **dR_array,    magma_minproduct_int_t ldr,
          magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);



void 
magma_minproductblas_scnrm2_cols_batched(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex **dA_array, magma_minproduct_int_t lda, 
    float **dxnorm_array, magma_minproduct_int_t batchCount);
 
void 
magma_minproduct_clarfgx_batched(magma_minproduct_int_t n, magma_minproductFloatComplex **dx0_array, magma_minproductFloatComplex **dx_array, 
                  magma_minproductFloatComplex **dtau_array, float **dxnorm_array, 
                  magma_minproductFloatComplex **dR_array, magma_minproduct_int_t it, magma_minproduct_int_t batchCount);


void 
magma_minproduct_clarfx_batched_v4(magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproductFloatComplex **v_array, magma_minproductFloatComplex **tau_array,
                magma_minproductFloatComplex **C_array, magma_minproduct_int_t ldc, float **xnorm_array, 
                magma_minproduct_int_t step, 
                magma_minproduct_int_t batchCount);


void 
magma_minproductblas_clarfg_batched(
    magma_minproduct_int_t n,
    magma_minproductFloatComplex** dalpha_array, magma_minproductFloatComplex** dx_array, magma_minproduct_int_t incx,
    magma_minproductFloatComplex** dtau_array, magma_minproduct_int_t batchCount );





// for debugging purpose
void 
cset_stepinit_ipiv(magma_minproduct_int_t **ipiv_array,
                 magma_minproduct_int_t pm,
                 magma_minproduct_int_t batchCount);



#ifdef __cplusplus
}
#endif

#undef COMPLEX

#endif  /* MAGMA_minproduct_CBATCHED_H */
