/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Azzam Haidar
       @author Tingxing Dong

       @generated from magma_minproduct_zbatched.h normal z -> s, Fri Jan 30 19:00:06 2015
*/

#ifndef MAGMA_minproduct_SBATCHED_H
#define MAGMA_minproduct_SBATCHED_H

#include "magma_minproduct_types.h"

#define REAL

#ifdef __cplusplus
extern "C" {
#endif
  /*
   *  local auxiliary routines
   */
void 
sset_pointer_int(magma_minproduct_int_t **output_array,
        magma_minproduct_int_t *input,
        magma_minproduct_int_t lda,
        magma_minproduct_int_t row, magma_minproduct_int_t column, 
        magma_minproduct_int_t batchSize,
        magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

void 
sset_pointer(float **output_array,
                 float *input,
                 magma_minproduct_int_t lda,
                 magma_minproduct_int_t row, magma_minproduct_int_t column,
                 magma_minproduct_int_t batchSize,
                 magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);


void 
sset_array(float **output_array,
               float **input_array, magma_minproduct_int_t lda,
               magma_minproduct_int_t row, magma_minproduct_int_t column,
               magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

void 
magma_minproduct_sdisplace_pointers(float **output_array,
               float **input_array, magma_minproduct_int_t lda,
               magma_minproduct_int_t row, magma_minproduct_int_t column, 
               magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

  /*
   *  LAPACK batched routines
   */

  /*
   *  BLAS batched routines
   */

void
magma_minproductblas_sgemm_batched(
    magma_minproduct_trans_t transA, magma_minproduct_trans_t transB, magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    float alpha,
    float const * const * dA_array, magma_minproduct_int_t ldda,
    float const * const * dB_array, magma_minproduct_int_t lddb,
    float beta,
    float **dC_array, magma_minproduct_int_t lddc, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue );

void
magma_minproductblas_sgemm_batched_lg(
    magma_minproduct_trans_t transA, magma_minproduct_trans_t transB, magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    float alpha,
    float const * const * dA_array, magma_minproduct_int_t ldda,
    float const * const * dB_array, magma_minproduct_int_t lddb,
    float beta,
    float **dC_array, magma_minproduct_int_t lddc, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue );
void
magma_minproductblas_sgemm_batched_k32(
    magma_minproduct_trans_t transA, magma_minproduct_trans_t transB, magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    float alpha,
    float const * const * dA_array, magma_minproduct_int_t ldda,
    float const * const * dB_array, magma_minproduct_int_t lddb,
    float beta,
    float **dC_array, magma_minproduct_int_t lddc, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue );



void 
magma_minproductblas_ssyrk_NC_batched( magma_minproduct_trans_t TRANSA, magma_minproduct_trans_t TRANSB, int m , int n , int k , 
                       float alpha, float **dA_array, int lda, 
                       float **B_array, int ldb, 
                       float beta,        float **C_array, int ldc, 
                       magma_minproduct_int_t batchCount);

void
magma_minproductblas_ssyrk_batched(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans, magma_minproduct_int_t n, magma_minproduct_int_t k,
    float alpha,
    float const * const * dA_array, magma_minproduct_int_t ldda,
    float beta,
    float **dC_array, magma_minproduct_int_t lddc, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue );

void
magma_minproductblas_ssyrk_batched_lg(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans, magma_minproduct_int_t n, magma_minproduct_int_t k,
    float alpha,
    float const * const * dA_array, magma_minproduct_int_t ldda,
    float beta,
    float **dC_array, magma_minproduct_int_t lddc, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue );

void
magma_minproductblas_ssyrk_batched_k32(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans, magma_minproduct_int_t n, magma_minproduct_int_t k,
    float alpha,
    float const * const * dA_array, magma_minproduct_int_t ldda,
    float beta,
    float **dC_array, magma_minproduct_int_t lddc, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue );






magma_minproduct_int_t 
magma_minproduct_spotf2_tile_batched(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n,
    float **dA_array, magma_minproduct_int_t lda,
    magma_minproduct_int_t *info_array, magma_minproduct_int_t gbstep, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_spotrf_panel(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t nb,     
    float *A, magma_minproduct_int_t lda,
    float* dinvA, float** dinvA_array, magma_minproduct_int_t invA_msize,
    float* x, float** x_array,  magma_minproduct_int_t x_msize,
    magma_minproduct_int_t *info_array,  magma_minproduct_int_t batchCount, magma_minproduct_int_t matrixSize, magma_minproduct_queue_t queue);


void 
magma_minproductblas_strtri_diag_batched(
    magma_minproduct_uplo_t uplo, magma_minproduct_diag_t diag, magma_minproduct_int_t n,
    float const * const *dA_array, magma_minproduct_int_t ldda,
    float **dinvA_array, 
    magma_minproduct_int_t resetozero, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);



void 
magma_minproductblas_strsm_batched(
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo, magma_minproduct_trans_t transA, magma_minproduct_diag_t diag,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    float alpha,
    float** dA_array,    magma_minproduct_int_t ldda,
    float** dB_array,    magma_minproduct_int_t lddb,
    magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);


void 
magma_minproductblas_strsm_work_batched(
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo, magma_minproduct_trans_t transA, magma_minproduct_diag_t diag,
    magma_minproduct_int_t flag, magma_minproduct_int_t m, magma_minproduct_int_t n, 
    float alpha, 
    float** dA_array,    magma_minproduct_int_t ldda,
    float** dB_array,    magma_minproduct_int_t lddb,
    float** dX_array,    magma_minproduct_int_t lddx, 
    float** dinvA_array, magma_minproduct_int_t dinvA_length,
    float** dA_displ, float** dB_displ, 
    float** dX_displ, float** dinvA_displ,
    magma_minproduct_int_t resetozero, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

void 
magma_minproductblas_strsm_outofplace_batched(
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo, magma_minproduct_trans_t transA, magma_minproduct_diag_t diag,
    magma_minproduct_int_t flag, magma_minproduct_int_t m, magma_minproduct_int_t n, 
    float alpha, 
    float** dA_array,    magma_minproduct_int_t ldda,
    float** dB_array,    magma_minproduct_int_t lddb,
    float** dX_array,    magma_minproduct_int_t lddx, 
    float** dinvA_array, magma_minproduct_int_t dinvA_length,
    float** dA_displ, float** dB_displ, 
    float** dX_displ, float** dinvA_displ,
    magma_minproduct_int_t resetozero, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_spotrf_batched(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    float **dA_array, magma_minproduct_int_t lda,
    magma_minproduct_int_t *info_array,  magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);


magma_minproduct_int_t 
magma_minproduct_spotf2_batched(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n,
    float **dA_array, magma_minproduct_int_t lda,
    float **dA_displ, 
    float **dW_displ,
    float **dB_displ, 
    float **dC_displ, 
    magma_minproduct_int_t *info_array, magma_minproduct_int_t gbstep, 
    magma_minproduct_int_t batchCount, cublasHandle_t myhandle, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_spotrf_panel_batched(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t nb,     
    float** dA_array,    magma_minproduct_int_t ldda,
    float** dX_array,    magma_minproduct_int_t dX_length,
    float** dinvA_array, magma_minproduct_int_t dinvA_length,
    float** dW0_displ, float** dW1_displ, 
    float** dW2_displ, float** dW3_displ,
    float** dW4_displ, 
    magma_minproduct_int_t *info_array, magma_minproduct_int_t gbstep,
    magma_minproduct_int_t batchCount, cublasHandle_t myhandle, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_spotrf_recpanel_batched(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n, 
    magma_minproduct_int_t min_recpnb,    
    float** dA_array,    magma_minproduct_int_t ldda,
    float** dX_array,    magma_minproduct_int_t dX_length,
    float** dinvA_array, magma_minproduct_int_t dinvA_length,
    float** dW0_displ, float** dW1_displ,  
    float** dW2_displ, float** dW3_displ,
    float** dW4_displ,
    magma_minproduct_int_t *info_array, magma_minproduct_int_t gbstep, 
    magma_minproduct_int_t batchCount, cublasHandle_t myhandle, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_spotrf_rectile_batched(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n, 
    magma_minproduct_int_t min_recpnb,    
    float** dA_array,    magma_minproduct_int_t ldda,
    float** dX_array,    magma_minproduct_int_t dX_length,
    float** dinvA_array, magma_minproduct_int_t dinvA_length,
    float** dW0_displ, float** dW1_displ,  
    float** dW2_displ, float** dW3_displ,
    float** dW4_displ,
    magma_minproduct_int_t *info_array, magma_minproduct_int_t gbstep,
    magma_minproduct_int_t batchCount, cublasHandle_t myhandle, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_spotrs_batched(
                  magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
                  float **dA_array, magma_minproduct_int_t ldda,
                  float **dB_array, magma_minproduct_int_t lddb,
                  magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_sposv_batched(
                  magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
                  float **dA_array, magma_minproduct_int_t ldda,
                  float **dB_array, magma_minproduct_int_t lddb,
                  magma_minproduct_int_t *dinfo_array,
                  magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_sgetrs_batched(
                  magma_minproduct_trans_t trans, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
                  float **dA_array, magma_minproduct_int_t ldda,
                  magma_minproduct_int_t **dipiv_array, 
                  float **dB_array, magma_minproduct_int_t lddb,
                  magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);



void 
magma_minproduct_slaswp_rowparallel_batched( magma_minproduct_int_t n, float** input_array, magma_minproduct_int_t ldi,
                   float** output_array, magma_minproduct_int_t ldo,
                   magma_minproduct_int_t k1, magma_minproduct_int_t k2,
                   magma_minproduct_int_t **pivinfo_array, 
                   magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue );

void 
magma_minproduct_slaswp_rowserial_batched(magma_minproduct_int_t n, float** dA_array, magma_minproduct_int_t lda,
                   magma_minproduct_int_t k1, magma_minproduct_int_t k2,
                   magma_minproduct_int_t **ipiv_array, 
                   magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

void 
magma_minproduct_slaswp_columnserial_batched(magma_minproduct_int_t n, float** dA_array, magma_minproduct_int_t lda,
                   magma_minproduct_int_t k1, magma_minproduct_int_t k2,
                   magma_minproduct_int_t **ipiv_array, 
                   magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

void 
magma_minproductblas_stranspose_batched(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    float **dA_array,  magma_minproduct_int_t ldda,
    float **dAT_array, magma_minproduct_int_t lddat, magma_minproduct_int_t batchCount );


void 
magma_minproductblas_slaset_batched(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n,
    float offdiag, float diag,
    magma_minproductFloat_ptr dAarray[], magma_minproduct_int_t ldda,
    magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

void 
magma_minproductblas_smemset_batched(magma_minproduct_int_t length, 
        magma_minproductFloat_ptr dAarray[], float val, 
        magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_sgetf2_batched(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    float **dA_array, magma_minproduct_int_t lda,
    float **GERA_array,
    float **GERB_array,
    float **GERC_array,
    magma_minproduct_int_t **ipiv_array,
    magma_minproduct_int_t *info_array, 
    magma_minproduct_int_t gbstep,            
    magma_minproduct_int_t batchCount,
    cublasHandle_t myhandle, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_sgetrf_recpanel_batched(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t min_recpnb,    
    float** dA_array,    magma_minproduct_int_t ldda,
    magma_minproduct_int_t** dipiv_array, magma_minproduct_int_t** dpivinfo_array,
    float** dX_array,    magma_minproduct_int_t dX_length,
    float** dinvA_array, magma_minproduct_int_t dinvA_length,
    float** dW1_displ, float** dW2_displ,  
    float** dW3_displ, float** dW4_displ,
    float** dW5_displ,
    magma_minproduct_int_t *info_array, magma_minproduct_int_t gbstep,
    magma_minproduct_int_t batchCount, cublasHandle_t myhandle, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_sgetrf_batched(
        magma_minproduct_int_t m, magma_minproduct_int_t n,
        float **dA_array, 
        magma_minproduct_int_t lda,
        magma_minproduct_int_t **ipiv_array, 
        magma_minproduct_int_t *info_array, 
        magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_sgetri_outofplace_batched( magma_minproduct_int_t n, 
                  float **dA_array, magma_minproduct_int_t ldda,
                  magma_minproduct_int_t **dipiv_array, 
                  float **dinvA_array, magma_minproduct_int_t lddia,
                  magma_minproduct_int_t *info_array,
                  magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

void 
magma_minproduct_sdisplace_intpointers(magma_minproduct_int_t **output_array,
               magma_minproduct_int_t **input_array, magma_minproduct_int_t lda,
               magma_minproduct_int_t row, magma_minproduct_int_t column, 
               magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);





void 
magma_minproductblas_isamax_atomic_batched(magma_minproduct_int_t n, float** x_array, magma_minproduct_int_t incx, magma_minproduct_int_t **max_id_array, magma_minproduct_int_t batchCount);

void 
magma_minproductblas_isamax_tree_batched(magma_minproduct_int_t n, float** x_array, magma_minproduct_int_t incx, magma_minproduct_int_t **max_id_array, magma_minproduct_int_t batchCount);



void 
magma_minproductblas_isamax_batched(magma_minproduct_int_t n, float** x_array, magma_minproduct_int_t incx, magma_minproduct_int_t **max_id_array, magma_minproduct_int_t batchCount);

void 
magma_minproductblas_isamax(magma_minproduct_int_t n, float* x, magma_minproduct_int_t incx, magma_minproduct_int_t *max_id);


magma_minproduct_int_t 
magma_minproduct_isamax_batched(magma_minproduct_int_t length, 
        float **x_array, magma_minproduct_int_t incx, magma_minproduct_int_t step,  magma_minproduct_int_t lda,
        magma_minproduct_int_t** ipiv_array, magma_minproduct_int_t *info_array, magma_minproduct_int_t gbstep, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_sswap_batched(magma_minproduct_int_t n, float **x_array, magma_minproduct_int_t incx, magma_minproduct_int_t j, 
                 magma_minproduct_int_t** ipiv_array, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_sscal_sger_batched(magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t step,
                                      float **dA_array, magma_minproduct_int_t lda,
                                      magma_minproduct_int_t *info_array, magma_minproduct_int_t gbstep, 
                                      magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_scomputecolumn_batched(magma_minproduct_int_t m, magma_minproduct_int_t paneloffset, magma_minproduct_int_t step, 
                                        float **dA_array,  magma_minproduct_int_t lda,
                                        magma_minproduct_int_t **ipiv_array, 
                                        magma_minproduct_int_t *info_array, magma_minproduct_int_t gbstep, 
                                        magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

void 
magma_minproduct_sgetf2trsm_batched(magma_minproduct_int_t ib, magma_minproduct_int_t n, float **dA_array,  magma_minproduct_int_t j, magma_minproduct_int_t lda,
                       magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);


magma_minproduct_int_t 
magma_minproduct_sgetf2_nopiv_batched(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    float **dA_array, magma_minproduct_int_t lda,
    float **dW0_displ,
    float **dW1_displ,
    float **dW2_displ,
    magma_minproduct_int_t *info_array,            
    magma_minproduct_int_t gbstep, 
    magma_minproduct_int_t batchCount,
    cublasHandle_t myhandle, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_sgetrf_panel_nopiv_batched(
    magma_minproduct_int_t m, magma_minproduct_int_t nb,    
    float** dA_array,    magma_minproduct_int_t ldda,
    float** dX_array,    magma_minproduct_int_t dX_length,
    float** dinvA_array, magma_minproduct_int_t dinvA_length,
    float** dW0_displ, float** dW1_displ,  
    float** dW2_displ, float** dW3_displ,
    float** dW4_displ,     
    magma_minproduct_int_t *info_array, magma_minproduct_int_t gbstep, 
    magma_minproduct_int_t batchCount, cublasHandle_t myhandle, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_sgetrf_recpanel_nopiv_batched(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t min_recpnb,    
    float** dA_array,    magma_minproduct_int_t ldda,
    float** dX_array,    magma_minproduct_int_t dX_length,
    float** dinvA_array, magma_minproduct_int_t dinvA_length,
    float** dW1_displ, float** dW2_displ,  
    float** dW3_displ, float** dW4_displ,
    float** dW5_displ, 
    magma_minproduct_int_t *info_array, magma_minproduct_int_t gbstep,
    magma_minproduct_int_t batchCount, cublasHandle_t myhandle, magma_minproduct_queue_t queue);


magma_minproduct_int_t 
magma_minproduct_sgetrf_nopiv_batched(
        magma_minproduct_int_t m, magma_minproduct_int_t n,
        float **dA_array, 
        magma_minproduct_int_t lda,
        magma_minproduct_int_t *info_array, 
        magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_sgetrs_nopiv_batched(
                  magma_minproduct_trans_t trans, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
                  float **dA_array, magma_minproduct_int_t ldda,
                  float **dB_array, magma_minproduct_int_t lddb,
                  magma_minproduct_int_t *info_array,
                  magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_sgesv_nopiv_batched(
                  magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
                  float **dA_array, magma_minproduct_int_t ldda,
                  float **dB_array, magma_minproduct_int_t lddb,
                  magma_minproduct_int_t *info_array,
                  magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_sgesv_rbt_batched(
                  magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
                  float **dA_array, magma_minproduct_int_t ldda,
                  float **dB_array, magma_minproduct_int_t lddb,
                  magma_minproduct_int_t *info_array,
                  magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_sgesv_batched(
                  magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
                  float **dA_array, magma_minproduct_int_t ldda,
                  magma_minproduct_int_t **dipiv_array, 
                  float **dB_array, magma_minproduct_int_t lddb,
                  magma_minproduct_int_t *dinfo_array,
                  magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

magma_minproduct_int_t
magma_minproduct_sgerbt_batched(
    magma_minproduct_bool_t gen, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
    float **dA_array, magma_minproduct_int_t ldda,
    float **dB_array, magma_minproduct_int_t lddb,
    float *U, float *V,
    magma_minproduct_int_t *info, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

void 
magma_minproductblas_sprbt_batched(
    magma_minproduct_int_t n, 
    float **dA_array, magma_minproduct_int_t ldda, 
    float *du, float *dv,
    magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

void
magma_minproductblas_sprbt_mv_batched(
    magma_minproduct_int_t n, 
    float *dv, float **db_array, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);


void
magma_minproductblas_sprbt_mtv_batched(
    magma_minproduct_int_t n, 
    float *du, float **db_array, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);





void 
magma_minproduct_slacgv_batched(magma_minproduct_int_t n, float **x_array, magma_minproduct_int_t incx, int offset, int batchCount, magma_minproduct_queue_t queue);

void 
magma_minproduct_spotf2_sscal_batched(magma_minproduct_int_t n, float **x_array, magma_minproduct_int_t incx, magma_minproduct_int_t offset, magma_minproduct_int_t *info_array, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

void 
magma_minproduct_spotf2_sdot_batched(magma_minproduct_int_t n, float **x_array, magma_minproduct_int_t incx, magma_minproduct_int_t offset, magma_minproduct_int_t *info_array, magma_minproduct_int_t gbstep, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);


void 
setup_pivinfo( magma_minproduct_int_t *pivinfo, magma_minproduct_int_t *ipiv, 
                      magma_minproduct_int_t m, magma_minproduct_int_t nb, 
                      magma_minproduct_queue_t queue);


void
magma_minproductblas_sgeadd_batched_q(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    float alpha,
    magma_minproductFloat_const_ptr  const dAarray[], magma_minproduct_int_t ldda,
    magma_minproductFloat_ptr              dBarray[], magma_minproduct_int_t lddb,
    magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue );

void
magma_minproductblas_slacpy_batched(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloat_const_ptr  const dAarray[], magma_minproduct_int_t ldda,
    magma_minproductFloat_ptr              dBarray[], magma_minproduct_int_t lddb,
    magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue );

void
magma_minproductblas_sgeadd_batched(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    float alpha,
    magma_minproductFloat_const_ptr  const dAarray[], magma_minproduct_int_t ldda,
    magma_minproductFloat_ptr              dBarray[], magma_minproduct_int_t lddb,
    magma_minproduct_int_t batchCount );


void
magma_minproductblas_sgemv_batched(
    magma_minproduct_trans_t trans, magma_minproduct_int_t m, magma_minproduct_int_t n,
    float alpha,
    magma_minproductFloat_ptr dA_array[], magma_minproduct_int_t ldda,
    magma_minproductFloat_ptr dx_array[], magma_minproduct_int_t incx,
    float beta,
    magma_minproductFloat_ptr dy_array[], magma_minproduct_int_t incy,
    magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);


magma_minproduct_int_t 
magma_minproduct_sgeqrf_batched(
    magma_minproduct_int_t m, magma_minproduct_int_t n, 
    float **dA_array,
    magma_minproduct_int_t lda, 
    float **tau_array,
    magma_minproduct_int_t *info_array, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_sgeqrf_batched_v4(
    magma_minproduct_int_t m, magma_minproduct_int_t n, 
    float **dA_array,
    magma_minproduct_int_t lda, 
    float **tau_array,
    magma_minproduct_int_t *info_array, magma_minproduct_int_t batchCount);

magma_minproduct_int_t
magma_minproduct_sgeqrf_panel_batched(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t nb,    
    float** dA_array,    magma_minproduct_int_t ldda,
    float** tau_array, 
    float** dT_array, magma_minproduct_int_t ldt, 
    float** dR_array, magma_minproduct_int_t ldr,
    float** dW0_displ, 
    float** dW1_displ,
    float *dwork,  
    float** W_array, 
    float** W2_array,
    magma_minproduct_int_t *info_array,
    magma_minproduct_int_t batchCount, cublasHandle_t myhandle, magma_minproduct_queue_t queue);

magma_minproduct_int_t
magma_minproduct_sgeqrf_panel_batched_v4(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t nb,    
    float** dA_array,    magma_minproduct_int_t ldda,
    float** tau_array, 
    float** dT_array, magma_minproduct_int_t ldt, 
    float** dR_array, magma_minproduct_int_t ldr,
    float** dnorm_array,  
    float** dW0_displ, 
    float** dW1_displ,
    float *dwork,  
    float** W_array, 
    float** W2_array,
    magma_minproduct_int_t *info_array,
    magma_minproduct_int_t batchCount, cublasHandle_t myhandle);

magma_minproduct_int_t
magma_minproduct_sgeqr2x_batched_v4(magma_minproduct_int_t m, magma_minproduct_int_t n, float **dA_array,
                  magma_minproduct_int_t lda, 
                  float **tau_array,
                  float **dR_array, magma_minproduct_int_t ldr,
                  float **dwork_array,  
                  magma_minproduct_int_t *info, magma_minproduct_int_t batchCount);

magma_minproduct_int_t
magma_minproduct_sgeqr2_batched(magma_minproduct_int_t m, magma_minproduct_int_t n, float **dA_array,
                  magma_minproduct_int_t lda, 
                  float **tau_array,
                  magma_minproduct_int_t *info, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

magma_minproduct_int_t
magma_minproduct_slarfb_sgemm_batched(
                  cublasHandle_t myhandle,
                  magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
                  float **dV_array,    magma_minproduct_int_t ldv,
                  float **dT_array,    magma_minproduct_int_t ldt,
                  float **dA_array,    magma_minproduct_int_t lda,
                  float **W_array,     magma_minproduct_int_t ldw,
                  float **W2_array,    magma_minproduct_int_t ldw2,
                  magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

magma_minproduct_int_t
magma_minproduct_slarfb_gemm_batched(
    magma_minproduct_side_t side, magma_minproduct_trans_t trans, magma_minproduct_direct_t direct, magma_minproduct_storev_t storev,
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductFloat_const_ptr dV_array[],    magma_minproduct_int_t lddv,
    magma_minproductFloat_const_ptr dT_array[],    magma_minproduct_int_t lddt,
    magma_minproductFloat_ptr dC_array[],          magma_minproduct_int_t lddc,
    magma_minproductFloat_ptr dwork_array[],       magma_minproduct_int_t ldwork,
    magma_minproductFloat_ptr dworkvt_array[],     magma_minproduct_int_t ldworkvt,
    magma_minproduct_int_t batchCount, cublasHandle_t myhandle, magma_minproduct_queue_t queue);




void
magma_minproduct_slarft_batched_vold(magma_minproduct_int_t n, magma_minproduct_int_t k, float **v_array, magma_minproduct_int_t ldv,
                    float **tau_array,
                    float **T_array, magma_minproduct_int_t ldt, 
                    magma_minproduct_int_t batchCount);





magma_minproduct_int_t
magma_minproduct_slarft_batched(magma_minproduct_int_t n, magma_minproduct_int_t k, magma_minproduct_int_t stair_T, 
                float **v_array, magma_minproduct_int_t ldv,
                float **tau_array, float **T_array, magma_minproduct_int_t ldt, 
                float **work_array, magma_minproduct_int_t lwork, magma_minproduct_int_t batchCount, cublasHandle_t myhandle, magma_minproduct_queue_t queue);

void
magma_minproduct_slarft_sm32x32_batched(magma_minproduct_int_t n, magma_minproduct_int_t k, float **v_array, magma_minproduct_int_t ldv,
                    float **tau_array, float **T_array, magma_minproduct_int_t ldt, 
                    magma_minproduct_int_t batchCount, cublasHandle_t myhandle, magma_minproduct_queue_t queue);



void magma_minproductblas_slarft_recstrmv_sm32x32(
    magma_minproduct_int_t m, magma_minproduct_int_t n, 
    float *tau, 
    float *Trec, magma_minproduct_int_t ldtrec, 
    float *Ttri, magma_minproduct_int_t ldttri);


void magma_minproductblas_slarft_recstrmv_sm32x32_batched(
    magma_minproduct_int_t m, magma_minproduct_int_t n, 
    float **tau_array, 
    float **Trec_array, magma_minproduct_int_t ldtrec, 
    float **Ttri_array, magma_minproduct_int_t ldttri,
    magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

void magma_minproductblas_slarft_strmv_sm32x32(
    magma_minproduct_int_t m, magma_minproduct_int_t n, 
    float *tau, 
    float *Tin, magma_minproduct_int_t ldtin, 
    float *Tout, magma_minproduct_int_t ldtout);

void magma_minproductblas_slarft_strmv_sm32x32_batched(
    magma_minproduct_int_t m, magma_minproduct_int_t n, 
    float **tau_array, 
    float **Tin_array, magma_minproduct_int_t ldtin, 
    float **Tout_array, magma_minproduct_int_t ldtout,
    magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

void magma_minproductblas_slarft_gemv_loop_inside(
    int n, int k, 
    float *tau, 
    float *v, int ldv, 
    float *T, int ldt);

void magma_minproductblas_slarft_gemv_loop_inside_batched(
    int n, int k, 
    float **tau_array, 
    float **v_array, int ldv, 
    float **T_array, int ldt, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

void magma_minproductblas_slarft_gemvrowwise(
    magma_minproduct_int_t m, magma_minproduct_int_t i, 
    float *tau, 
    float *v, magma_minproduct_int_t ldv, 
    float *T, magma_minproduct_int_t ldt,
    float *W);

void magma_minproductblas_slarft_gemvrowwise_batched(
    magma_minproduct_int_t m, magma_minproduct_int_t i, 
    float **tau_array, 
    float **v_array, magma_minproduct_int_t ldv, 
    float **T_array, magma_minproduct_int_t ldt,
    magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

void magma_minproductblas_slarft_gemvcolwise(
    magma_minproduct_int_t m,  magma_minproduct_int_t step,
    float *v, magma_minproduct_int_t ldv, 
    float *T,  magma_minproduct_int_t ldt,
    float *tau);

void magma_minproductblas_slarft_gemvcolwise_batched(
    magma_minproduct_int_t m,  magma_minproduct_int_t step,
    float **v_array, magma_minproduct_int_t ldv, 
    float **T_array,  magma_minproduct_int_t ldt,
    float **tau_array, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);




void sgeqrf_copy_upper_batched(                
                  magma_minproduct_int_t n, magma_minproduct_int_t nb,
                  float **dV_array,    magma_minproduct_int_t ldv,
                  float **dR_array,    magma_minproduct_int_t ldr,
          magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);



void 
magma_minproductblas_snrm2_cols_batched(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    float **dA_array, magma_minproduct_int_t lda, 
    float **dxnorm_array, magma_minproduct_int_t batchCount);
 
void 
magma_minproduct_slarfgx_batched(magma_minproduct_int_t n, float **dx0_array, float **dx_array, 
                  float **dtau_array, float **dxnorm_array, 
                  float **dR_array, magma_minproduct_int_t it, magma_minproduct_int_t batchCount);


void 
magma_minproduct_slarfx_batched_v4(magma_minproduct_int_t m, magma_minproduct_int_t n, float **v_array, float **tau_array,
                float **C_array, magma_minproduct_int_t ldc, float **xnorm_array, 
                magma_minproduct_int_t step, 
                magma_minproduct_int_t batchCount);


void 
magma_minproductblas_slarfg_batched(
    magma_minproduct_int_t n,
    float** dalpha_array, float** dx_array, magma_minproduct_int_t incx,
    float** dtau_array, magma_minproduct_int_t batchCount );





// for debugging purpose
void 
sset_stepinit_ipiv(magma_minproduct_int_t **ipiv_array,
                 magma_minproduct_int_t pm,
                 magma_minproduct_int_t batchCount);



#ifdef __cplusplus
}
#endif

#undef REAL

#endif  /* MAGMA_minproduct_SBATCHED_H */
