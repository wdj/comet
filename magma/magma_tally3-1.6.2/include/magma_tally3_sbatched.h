/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Azzam Haidar
       @author Tingxing Dong

       @generated from magma_tally3_zbatched.h normal z -> s, Fri Jan 30 19:00:06 2015
*/

#ifndef MAGMA_tally3_SBATCHED_H
#define MAGMA_tally3_SBATCHED_H

#include "magma_tally3_types.h"

#define REAL

#ifdef __cplusplus
extern "C" {
#endif
  /*
   *  local auxiliary routines
   */
void 
sset_pointer_int(magma_tally3_int_t **output_array,
        magma_tally3_int_t *input,
        magma_tally3_int_t lda,
        magma_tally3_int_t row, magma_tally3_int_t column, 
        magma_tally3_int_t batchSize,
        magma_tally3_int_t batchCount, magma_tally3_queue_t queue);

void 
sset_pointer(float **output_array,
                 float *input,
                 magma_tally3_int_t lda,
                 magma_tally3_int_t row, magma_tally3_int_t column,
                 magma_tally3_int_t batchSize,
                 magma_tally3_int_t batchCount, magma_tally3_queue_t queue);


void 
sset_array(float **output_array,
               float **input_array, magma_tally3_int_t lda,
               magma_tally3_int_t row, magma_tally3_int_t column,
               magma_tally3_int_t batchCount, magma_tally3_queue_t queue);

void 
magma_tally3_sdisplace_pointers(float **output_array,
               float **input_array, magma_tally3_int_t lda,
               magma_tally3_int_t row, magma_tally3_int_t column, 
               magma_tally3_int_t batchCount, magma_tally3_queue_t queue);

  /*
   *  LAPACK batched routines
   */

  /*
   *  BLAS batched routines
   */

void
magma_tally3blas_sgemm_batched(
    magma_tally3_trans_t transA, magma_tally3_trans_t transB, magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    float alpha,
    float const * const * dA_array, magma_tally3_int_t ldda,
    float const * const * dB_array, magma_tally3_int_t lddb,
    float beta,
    float **dC_array, magma_tally3_int_t lddc, magma_tally3_int_t batchCount, magma_tally3_queue_t queue );

void
magma_tally3blas_sgemm_batched_lg(
    magma_tally3_trans_t transA, magma_tally3_trans_t transB, magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    float alpha,
    float const * const * dA_array, magma_tally3_int_t ldda,
    float const * const * dB_array, magma_tally3_int_t lddb,
    float beta,
    float **dC_array, magma_tally3_int_t lddc, magma_tally3_int_t batchCount, magma_tally3_queue_t queue );
void
magma_tally3blas_sgemm_batched_k32(
    magma_tally3_trans_t transA, magma_tally3_trans_t transB, magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    float alpha,
    float const * const * dA_array, magma_tally3_int_t ldda,
    float const * const * dB_array, magma_tally3_int_t lddb,
    float beta,
    float **dC_array, magma_tally3_int_t lddc, magma_tally3_int_t batchCount, magma_tally3_queue_t queue );



void 
magma_tally3blas_ssyrk_NC_batched( magma_tally3_trans_t TRANSA, magma_tally3_trans_t TRANSB, int m , int n , int k , 
                       float alpha, float **dA_array, int lda, 
                       float **B_array, int ldb, 
                       float beta,        float **C_array, int ldc, 
                       magma_tally3_int_t batchCount);

void
magma_tally3blas_ssyrk_batched(
    magma_tally3_uplo_t uplo, magma_tally3_trans_t trans, magma_tally3_int_t n, magma_tally3_int_t k,
    float alpha,
    float const * const * dA_array, magma_tally3_int_t ldda,
    float beta,
    float **dC_array, magma_tally3_int_t lddc, magma_tally3_int_t batchCount, magma_tally3_queue_t queue );

void
magma_tally3blas_ssyrk_batched_lg(
    magma_tally3_uplo_t uplo, magma_tally3_trans_t trans, magma_tally3_int_t n, magma_tally3_int_t k,
    float alpha,
    float const * const * dA_array, magma_tally3_int_t ldda,
    float beta,
    float **dC_array, magma_tally3_int_t lddc, magma_tally3_int_t batchCount, magma_tally3_queue_t queue );

void
magma_tally3blas_ssyrk_batched_k32(
    magma_tally3_uplo_t uplo, magma_tally3_trans_t trans, magma_tally3_int_t n, magma_tally3_int_t k,
    float alpha,
    float const * const * dA_array, magma_tally3_int_t ldda,
    float beta,
    float **dC_array, magma_tally3_int_t lddc, magma_tally3_int_t batchCount, magma_tally3_queue_t queue );






magma_tally3_int_t 
magma_tally3_spotf2_tile_batched(
    magma_tally3_uplo_t uplo, magma_tally3_int_t m, magma_tally3_int_t n,
    float **dA_array, magma_tally3_int_t lda,
    magma_tally3_int_t *info_array, magma_tally3_int_t gbstep, magma_tally3_int_t batchCount, magma_tally3_queue_t queue);

magma_tally3_int_t 
magma_tally3_spotrf_panel(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n, magma_tally3_int_t nb,     
    float *A, magma_tally3_int_t lda,
    float* dinvA, float** dinvA_array, magma_tally3_int_t invA_msize,
    float* x, float** x_array,  magma_tally3_int_t x_msize,
    magma_tally3_int_t *info_array,  magma_tally3_int_t batchCount, magma_tally3_int_t matrixSize, magma_tally3_queue_t queue);


void 
magma_tally3blas_strtri_diag_batched(
    magma_tally3_uplo_t uplo, magma_tally3_diag_t diag, magma_tally3_int_t n,
    float const * const *dA_array, magma_tally3_int_t ldda,
    float **dinvA_array, 
    magma_tally3_int_t resetozero, magma_tally3_int_t batchCount, magma_tally3_queue_t queue);



void 
magma_tally3blas_strsm_batched(
    magma_tally3_side_t side, magma_tally3_uplo_t uplo, magma_tally3_trans_t transA, magma_tally3_diag_t diag,
    magma_tally3_int_t m, magma_tally3_int_t n,
    float alpha,
    float** dA_array,    magma_tally3_int_t ldda,
    float** dB_array,    magma_tally3_int_t lddb,
    magma_tally3_int_t batchCount, magma_tally3_queue_t queue);


void 
magma_tally3blas_strsm_work_batched(
    magma_tally3_side_t side, magma_tally3_uplo_t uplo, magma_tally3_trans_t transA, magma_tally3_diag_t diag,
    magma_tally3_int_t flag, magma_tally3_int_t m, magma_tally3_int_t n, 
    float alpha, 
    float** dA_array,    magma_tally3_int_t ldda,
    float** dB_array,    magma_tally3_int_t lddb,
    float** dX_array,    magma_tally3_int_t lddx, 
    float** dinvA_array, magma_tally3_int_t dinvA_length,
    float** dA_displ, float** dB_displ, 
    float** dX_displ, float** dinvA_displ,
    magma_tally3_int_t resetozero, magma_tally3_int_t batchCount, magma_tally3_queue_t queue);

void 
magma_tally3blas_strsm_outofplace_batched(
    magma_tally3_side_t side, magma_tally3_uplo_t uplo, magma_tally3_trans_t transA, magma_tally3_diag_t diag,
    magma_tally3_int_t flag, magma_tally3_int_t m, magma_tally3_int_t n, 
    float alpha, 
    float** dA_array,    magma_tally3_int_t ldda,
    float** dB_array,    magma_tally3_int_t lddb,
    float** dX_array,    magma_tally3_int_t lddx, 
    float** dinvA_array, magma_tally3_int_t dinvA_length,
    float** dA_displ, float** dB_displ, 
    float** dX_displ, float** dinvA_displ,
    magma_tally3_int_t resetozero, magma_tally3_int_t batchCount, magma_tally3_queue_t queue);

magma_tally3_int_t 
magma_tally3_spotrf_batched(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    float **dA_array, magma_tally3_int_t lda,
    magma_tally3_int_t *info_array,  magma_tally3_int_t batchCount, magma_tally3_queue_t queue);


magma_tally3_int_t 
magma_tally3_spotf2_batched(
    magma_tally3_uplo_t uplo, magma_tally3_int_t m, magma_tally3_int_t n,
    float **dA_array, magma_tally3_int_t lda,
    float **dA_displ, 
    float **dW_displ,
    float **dB_displ, 
    float **dC_displ, 
    magma_tally3_int_t *info_array, magma_tally3_int_t gbstep, 
    magma_tally3_int_t batchCount, cublasHandle_t myhandle, magma_tally3_queue_t queue);

magma_tally3_int_t 
magma_tally3_spotrf_panel_batched(
    magma_tally3_uplo_t uplo, magma_tally3_int_t n, magma_tally3_int_t nb,     
    float** dA_array,    magma_tally3_int_t ldda,
    float** dX_array,    magma_tally3_int_t dX_length,
    float** dinvA_array, magma_tally3_int_t dinvA_length,
    float** dW0_displ, float** dW1_displ, 
    float** dW2_displ, float** dW3_displ,
    float** dW4_displ, 
    magma_tally3_int_t *info_array, magma_tally3_int_t gbstep,
    magma_tally3_int_t batchCount, cublasHandle_t myhandle, magma_tally3_queue_t queue);

magma_tally3_int_t 
magma_tally3_spotrf_recpanel_batched(
    magma_tally3_uplo_t uplo, magma_tally3_int_t m, magma_tally3_int_t n, 
    magma_tally3_int_t min_recpnb,    
    float** dA_array,    magma_tally3_int_t ldda,
    float** dX_array,    magma_tally3_int_t dX_length,
    float** dinvA_array, magma_tally3_int_t dinvA_length,
    float** dW0_displ, float** dW1_displ,  
    float** dW2_displ, float** dW3_displ,
    float** dW4_displ,
    magma_tally3_int_t *info_array, magma_tally3_int_t gbstep, 
    magma_tally3_int_t batchCount, cublasHandle_t myhandle, magma_tally3_queue_t queue);

magma_tally3_int_t 
magma_tally3_spotrf_rectile_batched(
    magma_tally3_uplo_t uplo, magma_tally3_int_t m, magma_tally3_int_t n, 
    magma_tally3_int_t min_recpnb,    
    float** dA_array,    magma_tally3_int_t ldda,
    float** dX_array,    magma_tally3_int_t dX_length,
    float** dinvA_array, magma_tally3_int_t dinvA_length,
    float** dW0_displ, float** dW1_displ,  
    float** dW2_displ, float** dW3_displ,
    float** dW4_displ,
    magma_tally3_int_t *info_array, magma_tally3_int_t gbstep,
    magma_tally3_int_t batchCount, cublasHandle_t myhandle, magma_tally3_queue_t queue);

magma_tally3_int_t 
magma_tally3_spotrs_batched(
                  magma_tally3_uplo_t uplo, magma_tally3_int_t n, magma_tally3_int_t nrhs,
                  float **dA_array, magma_tally3_int_t ldda,
                  float **dB_array, magma_tally3_int_t lddb,
                  magma_tally3_int_t batchCount, magma_tally3_queue_t queue);

magma_tally3_int_t 
magma_tally3_sposv_batched(
                  magma_tally3_uplo_t uplo, magma_tally3_int_t n, magma_tally3_int_t nrhs,
                  float **dA_array, magma_tally3_int_t ldda,
                  float **dB_array, magma_tally3_int_t lddb,
                  magma_tally3_int_t *dinfo_array,
                  magma_tally3_int_t batchCount, magma_tally3_queue_t queue);

magma_tally3_int_t 
magma_tally3_sgetrs_batched(
                  magma_tally3_trans_t trans, magma_tally3_int_t n, magma_tally3_int_t nrhs,
                  float **dA_array, magma_tally3_int_t ldda,
                  magma_tally3_int_t **dipiv_array, 
                  float **dB_array, magma_tally3_int_t lddb,
                  magma_tally3_int_t batchCount, magma_tally3_queue_t queue);



void 
magma_tally3_slaswp_rowparallel_batched( magma_tally3_int_t n, float** input_array, magma_tally3_int_t ldi,
                   float** output_array, magma_tally3_int_t ldo,
                   magma_tally3_int_t k1, magma_tally3_int_t k2,
                   magma_tally3_int_t **pivinfo_array, 
                   magma_tally3_int_t batchCount, magma_tally3_queue_t queue );

void 
magma_tally3_slaswp_rowserial_batched(magma_tally3_int_t n, float** dA_array, magma_tally3_int_t lda,
                   magma_tally3_int_t k1, magma_tally3_int_t k2,
                   magma_tally3_int_t **ipiv_array, 
                   magma_tally3_int_t batchCount, magma_tally3_queue_t queue);

void 
magma_tally3_slaswp_columnserial_batched(magma_tally3_int_t n, float** dA_array, magma_tally3_int_t lda,
                   magma_tally3_int_t k1, magma_tally3_int_t k2,
                   magma_tally3_int_t **ipiv_array, 
                   magma_tally3_int_t batchCount, magma_tally3_queue_t queue);

void 
magma_tally3blas_stranspose_batched(
    magma_tally3_int_t m, magma_tally3_int_t n,
    float **dA_array,  magma_tally3_int_t ldda,
    float **dAT_array, magma_tally3_int_t lddat, magma_tally3_int_t batchCount );


void 
magma_tally3blas_slaset_batched(
    magma_tally3_uplo_t uplo, magma_tally3_int_t m, magma_tally3_int_t n,
    float offdiag, float diag,
    magma_tally3Float_ptr dAarray[], magma_tally3_int_t ldda,
    magma_tally3_int_t batchCount, magma_tally3_queue_t queue);

void 
magma_tally3blas_smemset_batched(magma_tally3_int_t length, 
        magma_tally3Float_ptr dAarray[], float val, 
        magma_tally3_int_t batchCount, magma_tally3_queue_t queue);

magma_tally3_int_t 
magma_tally3_sgetf2_batched(
    magma_tally3_int_t m, magma_tally3_int_t n,
    float **dA_array, magma_tally3_int_t lda,
    float **GERA_array,
    float **GERB_array,
    float **GERC_array,
    magma_tally3_int_t **ipiv_array,
    magma_tally3_int_t *info_array, 
    magma_tally3_int_t gbstep,            
    magma_tally3_int_t batchCount,
    cublasHandle_t myhandle, magma_tally3_queue_t queue);

magma_tally3_int_t 
magma_tally3_sgetrf_recpanel_batched(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t min_recpnb,    
    float** dA_array,    magma_tally3_int_t ldda,
    magma_tally3_int_t** dipiv_array, magma_tally3_int_t** dpivinfo_array,
    float** dX_array,    magma_tally3_int_t dX_length,
    float** dinvA_array, magma_tally3_int_t dinvA_length,
    float** dW1_displ, float** dW2_displ,  
    float** dW3_displ, float** dW4_displ,
    float** dW5_displ,
    magma_tally3_int_t *info_array, magma_tally3_int_t gbstep,
    magma_tally3_int_t batchCount, cublasHandle_t myhandle, magma_tally3_queue_t queue);

magma_tally3_int_t 
magma_tally3_sgetrf_batched(
        magma_tally3_int_t m, magma_tally3_int_t n,
        float **dA_array, 
        magma_tally3_int_t lda,
        magma_tally3_int_t **ipiv_array, 
        magma_tally3_int_t *info_array, 
        magma_tally3_int_t batchCount, magma_tally3_queue_t queue);

magma_tally3_int_t 
magma_tally3_sgetri_outofplace_batched( magma_tally3_int_t n, 
                  float **dA_array, magma_tally3_int_t ldda,
                  magma_tally3_int_t **dipiv_array, 
                  float **dinvA_array, magma_tally3_int_t lddia,
                  magma_tally3_int_t *info_array,
                  magma_tally3_int_t batchCount, magma_tally3_queue_t queue);

void 
magma_tally3_sdisplace_intpointers(magma_tally3_int_t **output_array,
               magma_tally3_int_t **input_array, magma_tally3_int_t lda,
               magma_tally3_int_t row, magma_tally3_int_t column, 
               magma_tally3_int_t batchCount, magma_tally3_queue_t queue);





void 
magma_tally3blas_isamax_atomic_batched(magma_tally3_int_t n, float** x_array, magma_tally3_int_t incx, magma_tally3_int_t **max_id_array, magma_tally3_int_t batchCount);

void 
magma_tally3blas_isamax_tree_batched(magma_tally3_int_t n, float** x_array, magma_tally3_int_t incx, magma_tally3_int_t **max_id_array, magma_tally3_int_t batchCount);



void 
magma_tally3blas_isamax_batched(magma_tally3_int_t n, float** x_array, magma_tally3_int_t incx, magma_tally3_int_t **max_id_array, magma_tally3_int_t batchCount);

void 
magma_tally3blas_isamax(magma_tally3_int_t n, float* x, magma_tally3_int_t incx, magma_tally3_int_t *max_id);


magma_tally3_int_t 
magma_tally3_isamax_batched(magma_tally3_int_t length, 
        float **x_array, magma_tally3_int_t incx, magma_tally3_int_t step,  magma_tally3_int_t lda,
        magma_tally3_int_t** ipiv_array, magma_tally3_int_t *info_array, magma_tally3_int_t gbstep, magma_tally3_int_t batchCount, magma_tally3_queue_t queue);

magma_tally3_int_t 
magma_tally3_sswap_batched(magma_tally3_int_t n, float **x_array, magma_tally3_int_t incx, magma_tally3_int_t j, 
                 magma_tally3_int_t** ipiv_array, magma_tally3_int_t batchCount, magma_tally3_queue_t queue);

magma_tally3_int_t 
magma_tally3_sscal_sger_batched(magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t step,
                                      float **dA_array, magma_tally3_int_t lda,
                                      magma_tally3_int_t *info_array, magma_tally3_int_t gbstep, 
                                      magma_tally3_int_t batchCount, magma_tally3_queue_t queue);

magma_tally3_int_t 
magma_tally3_scomputecolumn_batched(magma_tally3_int_t m, magma_tally3_int_t paneloffset, magma_tally3_int_t step, 
                                        float **dA_array,  magma_tally3_int_t lda,
                                        magma_tally3_int_t **ipiv_array, 
                                        magma_tally3_int_t *info_array, magma_tally3_int_t gbstep, 
                                        magma_tally3_int_t batchCount, magma_tally3_queue_t queue);

void 
magma_tally3_sgetf2trsm_batched(magma_tally3_int_t ib, magma_tally3_int_t n, float **dA_array,  magma_tally3_int_t j, magma_tally3_int_t lda,
                       magma_tally3_int_t batchCount, magma_tally3_queue_t queue);


magma_tally3_int_t 
magma_tally3_sgetf2_nopiv_batched(
    magma_tally3_int_t m, magma_tally3_int_t n,
    float **dA_array, magma_tally3_int_t lda,
    float **dW0_displ,
    float **dW1_displ,
    float **dW2_displ,
    magma_tally3_int_t *info_array,            
    magma_tally3_int_t gbstep, 
    magma_tally3_int_t batchCount,
    cublasHandle_t myhandle, magma_tally3_queue_t queue);

magma_tally3_int_t 
magma_tally3_sgetrf_panel_nopiv_batched(
    magma_tally3_int_t m, magma_tally3_int_t nb,    
    float** dA_array,    magma_tally3_int_t ldda,
    float** dX_array,    magma_tally3_int_t dX_length,
    float** dinvA_array, magma_tally3_int_t dinvA_length,
    float** dW0_displ, float** dW1_displ,  
    float** dW2_displ, float** dW3_displ,
    float** dW4_displ,     
    magma_tally3_int_t *info_array, magma_tally3_int_t gbstep, 
    magma_tally3_int_t batchCount, cublasHandle_t myhandle, magma_tally3_queue_t queue);

magma_tally3_int_t 
magma_tally3_sgetrf_recpanel_nopiv_batched(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t min_recpnb,    
    float** dA_array,    magma_tally3_int_t ldda,
    float** dX_array,    magma_tally3_int_t dX_length,
    float** dinvA_array, magma_tally3_int_t dinvA_length,
    float** dW1_displ, float** dW2_displ,  
    float** dW3_displ, float** dW4_displ,
    float** dW5_displ, 
    magma_tally3_int_t *info_array, magma_tally3_int_t gbstep,
    magma_tally3_int_t batchCount, cublasHandle_t myhandle, magma_tally3_queue_t queue);


magma_tally3_int_t 
magma_tally3_sgetrf_nopiv_batched(
        magma_tally3_int_t m, magma_tally3_int_t n,
        float **dA_array, 
        magma_tally3_int_t lda,
        magma_tally3_int_t *info_array, 
        magma_tally3_int_t batchCount, magma_tally3_queue_t queue);

magma_tally3_int_t 
magma_tally3_sgetrs_nopiv_batched(
                  magma_tally3_trans_t trans, magma_tally3_int_t n, magma_tally3_int_t nrhs,
                  float **dA_array, magma_tally3_int_t ldda,
                  float **dB_array, magma_tally3_int_t lddb,
                  magma_tally3_int_t *info_array,
                  magma_tally3_int_t batchCount, magma_tally3_queue_t queue);

magma_tally3_int_t 
magma_tally3_sgesv_nopiv_batched(
                  magma_tally3_int_t n, magma_tally3_int_t nrhs,
                  float **dA_array, magma_tally3_int_t ldda,
                  float **dB_array, magma_tally3_int_t lddb,
                  magma_tally3_int_t *info_array,
                  magma_tally3_int_t batchCount, magma_tally3_queue_t queue);

magma_tally3_int_t 
magma_tally3_sgesv_rbt_batched(
                  magma_tally3_int_t n, magma_tally3_int_t nrhs,
                  float **dA_array, magma_tally3_int_t ldda,
                  float **dB_array, magma_tally3_int_t lddb,
                  magma_tally3_int_t *info_array,
                  magma_tally3_int_t batchCount, magma_tally3_queue_t queue);

magma_tally3_int_t 
magma_tally3_sgesv_batched(
                  magma_tally3_int_t n, magma_tally3_int_t nrhs,
                  float **dA_array, magma_tally3_int_t ldda,
                  magma_tally3_int_t **dipiv_array, 
                  float **dB_array, magma_tally3_int_t lddb,
                  magma_tally3_int_t *dinfo_array,
                  magma_tally3_int_t batchCount, magma_tally3_queue_t queue);

magma_tally3_int_t
magma_tally3_sgerbt_batched(
    magma_tally3_bool_t gen, magma_tally3_int_t n, magma_tally3_int_t nrhs,
    float **dA_array, magma_tally3_int_t ldda,
    float **dB_array, magma_tally3_int_t lddb,
    float *U, float *V,
    magma_tally3_int_t *info, magma_tally3_int_t batchCount, magma_tally3_queue_t queue);

void 
magma_tally3blas_sprbt_batched(
    magma_tally3_int_t n, 
    float **dA_array, magma_tally3_int_t ldda, 
    float *du, float *dv,
    magma_tally3_int_t batchCount, magma_tally3_queue_t queue);

void
magma_tally3blas_sprbt_mv_batched(
    magma_tally3_int_t n, 
    float *dv, float **db_array, magma_tally3_int_t batchCount, magma_tally3_queue_t queue);


void
magma_tally3blas_sprbt_mtv_batched(
    magma_tally3_int_t n, 
    float *du, float **db_array, magma_tally3_int_t batchCount, magma_tally3_queue_t queue);





void 
magma_tally3_slacgv_batched(magma_tally3_int_t n, float **x_array, magma_tally3_int_t incx, int offset, int batchCount, magma_tally3_queue_t queue);

void 
magma_tally3_spotf2_sscal_batched(magma_tally3_int_t n, float **x_array, magma_tally3_int_t incx, magma_tally3_int_t offset, magma_tally3_int_t *info_array, magma_tally3_int_t batchCount, magma_tally3_queue_t queue);

void 
magma_tally3_spotf2_sdot_batched(magma_tally3_int_t n, float **x_array, magma_tally3_int_t incx, magma_tally3_int_t offset, magma_tally3_int_t *info_array, magma_tally3_int_t gbstep, magma_tally3_int_t batchCount, magma_tally3_queue_t queue);


void 
setup_pivinfo( magma_tally3_int_t *pivinfo, magma_tally3_int_t *ipiv, 
                      magma_tally3_int_t m, magma_tally3_int_t nb, 
                      magma_tally3_queue_t queue);


void
magma_tally3blas_sgeadd_batched_q(
    magma_tally3_int_t m, magma_tally3_int_t n,
    float alpha,
    magma_tally3Float_const_ptr  const dAarray[], magma_tally3_int_t ldda,
    magma_tally3Float_ptr              dBarray[], magma_tally3_int_t lddb,
    magma_tally3_int_t batchCount, magma_tally3_queue_t queue );

void
magma_tally3blas_slacpy_batched(
    magma_tally3_uplo_t uplo, magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3Float_const_ptr  const dAarray[], magma_tally3_int_t ldda,
    magma_tally3Float_ptr              dBarray[], magma_tally3_int_t lddb,
    magma_tally3_int_t batchCount, magma_tally3_queue_t queue );

void
magma_tally3blas_sgeadd_batched(
    magma_tally3_int_t m, magma_tally3_int_t n,
    float alpha,
    magma_tally3Float_const_ptr  const dAarray[], magma_tally3_int_t ldda,
    magma_tally3Float_ptr              dBarray[], magma_tally3_int_t lddb,
    magma_tally3_int_t batchCount );


void
magma_tally3blas_sgemv_batched(
    magma_tally3_trans_t trans, magma_tally3_int_t m, magma_tally3_int_t n,
    float alpha,
    magma_tally3Float_ptr dA_array[], magma_tally3_int_t ldda,
    magma_tally3Float_ptr dx_array[], magma_tally3_int_t incx,
    float beta,
    magma_tally3Float_ptr dy_array[], magma_tally3_int_t incy,
    magma_tally3_int_t batchCount, magma_tally3_queue_t queue);


magma_tally3_int_t 
magma_tally3_sgeqrf_batched(
    magma_tally3_int_t m, magma_tally3_int_t n, 
    float **dA_array,
    magma_tally3_int_t lda, 
    float **tau_array,
    magma_tally3_int_t *info_array, magma_tally3_int_t batchCount, magma_tally3_queue_t queue);

magma_tally3_int_t 
magma_tally3_sgeqrf_batched_v4(
    magma_tally3_int_t m, magma_tally3_int_t n, 
    float **dA_array,
    magma_tally3_int_t lda, 
    float **tau_array,
    magma_tally3_int_t *info_array, magma_tally3_int_t batchCount);

magma_tally3_int_t
magma_tally3_sgeqrf_panel_batched(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t nb,    
    float** dA_array,    magma_tally3_int_t ldda,
    float** tau_array, 
    float** dT_array, magma_tally3_int_t ldt, 
    float** dR_array, magma_tally3_int_t ldr,
    float** dW0_displ, 
    float** dW1_displ,
    float *dwork,  
    float** W_array, 
    float** W2_array,
    magma_tally3_int_t *info_array,
    magma_tally3_int_t batchCount, cublasHandle_t myhandle, magma_tally3_queue_t queue);

magma_tally3_int_t
magma_tally3_sgeqrf_panel_batched_v4(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t nb,    
    float** dA_array,    magma_tally3_int_t ldda,
    float** tau_array, 
    float** dT_array, magma_tally3_int_t ldt, 
    float** dR_array, magma_tally3_int_t ldr,
    float** dnorm_array,  
    float** dW0_displ, 
    float** dW1_displ,
    float *dwork,  
    float** W_array, 
    float** W2_array,
    magma_tally3_int_t *info_array,
    magma_tally3_int_t batchCount, cublasHandle_t myhandle);

magma_tally3_int_t
magma_tally3_sgeqr2x_batched_v4(magma_tally3_int_t m, magma_tally3_int_t n, float **dA_array,
                  magma_tally3_int_t lda, 
                  float **tau_array,
                  float **dR_array, magma_tally3_int_t ldr,
                  float **dwork_array,  
                  magma_tally3_int_t *info, magma_tally3_int_t batchCount);

magma_tally3_int_t
magma_tally3_sgeqr2_batched(magma_tally3_int_t m, magma_tally3_int_t n, float **dA_array,
                  magma_tally3_int_t lda, 
                  float **tau_array,
                  magma_tally3_int_t *info, magma_tally3_int_t batchCount, magma_tally3_queue_t queue);

magma_tally3_int_t
magma_tally3_slarfb_sgemm_batched(
                  cublasHandle_t myhandle,
                  magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
                  float **dV_array,    magma_tally3_int_t ldv,
                  float **dT_array,    magma_tally3_int_t ldt,
                  float **dA_array,    magma_tally3_int_t lda,
                  float **W_array,     magma_tally3_int_t ldw,
                  float **W2_array,    magma_tally3_int_t ldw2,
                  magma_tally3_int_t batchCount, magma_tally3_queue_t queue);

magma_tally3_int_t
magma_tally3_slarfb_gemm_batched(
    magma_tally3_side_t side, magma_tally3_trans_t trans, magma_tally3_direct_t direct, magma_tally3_storev_t storev,
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k,
    magma_tally3Float_const_ptr dV_array[],    magma_tally3_int_t lddv,
    magma_tally3Float_const_ptr dT_array[],    magma_tally3_int_t lddt,
    magma_tally3Float_ptr dC_array[],          magma_tally3_int_t lddc,
    magma_tally3Float_ptr dwork_array[],       magma_tally3_int_t ldwork,
    magma_tally3Float_ptr dworkvt_array[],     magma_tally3_int_t ldworkvt,
    magma_tally3_int_t batchCount, cublasHandle_t myhandle, magma_tally3_queue_t queue);




void
magma_tally3_slarft_batched_vold(magma_tally3_int_t n, magma_tally3_int_t k, float **v_array, magma_tally3_int_t ldv,
                    float **tau_array,
                    float **T_array, magma_tally3_int_t ldt, 
                    magma_tally3_int_t batchCount);





magma_tally3_int_t
magma_tally3_slarft_batched(magma_tally3_int_t n, magma_tally3_int_t k, magma_tally3_int_t stair_T, 
                float **v_array, magma_tally3_int_t ldv,
                float **tau_array, float **T_array, magma_tally3_int_t ldt, 
                float **work_array, magma_tally3_int_t lwork, magma_tally3_int_t batchCount, cublasHandle_t myhandle, magma_tally3_queue_t queue);

void
magma_tally3_slarft_sm32x32_batched(magma_tally3_int_t n, magma_tally3_int_t k, float **v_array, magma_tally3_int_t ldv,
                    float **tau_array, float **T_array, magma_tally3_int_t ldt, 
                    magma_tally3_int_t batchCount, cublasHandle_t myhandle, magma_tally3_queue_t queue);



void magma_tally3blas_slarft_recstrmv_sm32x32(
    magma_tally3_int_t m, magma_tally3_int_t n, 
    float *tau, 
    float *Trec, magma_tally3_int_t ldtrec, 
    float *Ttri, magma_tally3_int_t ldttri);


void magma_tally3blas_slarft_recstrmv_sm32x32_batched(
    magma_tally3_int_t m, magma_tally3_int_t n, 
    float **tau_array, 
    float **Trec_array, magma_tally3_int_t ldtrec, 
    float **Ttri_array, magma_tally3_int_t ldttri,
    magma_tally3_int_t batchCount, magma_tally3_queue_t queue);

void magma_tally3blas_slarft_strmv_sm32x32(
    magma_tally3_int_t m, magma_tally3_int_t n, 
    float *tau, 
    float *Tin, magma_tally3_int_t ldtin, 
    float *Tout, magma_tally3_int_t ldtout);

void magma_tally3blas_slarft_strmv_sm32x32_batched(
    magma_tally3_int_t m, magma_tally3_int_t n, 
    float **tau_array, 
    float **Tin_array, magma_tally3_int_t ldtin, 
    float **Tout_array, magma_tally3_int_t ldtout,
    magma_tally3_int_t batchCount, magma_tally3_queue_t queue);

void magma_tally3blas_slarft_gemv_loop_inside(
    int n, int k, 
    float *tau, 
    float *v, int ldv, 
    float *T, int ldt);

void magma_tally3blas_slarft_gemv_loop_inside_batched(
    int n, int k, 
    float **tau_array, 
    float **v_array, int ldv, 
    float **T_array, int ldt, magma_tally3_int_t batchCount, magma_tally3_queue_t queue);

void magma_tally3blas_slarft_gemvrowwise(
    magma_tally3_int_t m, magma_tally3_int_t i, 
    float *tau, 
    float *v, magma_tally3_int_t ldv, 
    float *T, magma_tally3_int_t ldt,
    float *W);

void magma_tally3blas_slarft_gemvrowwise_batched(
    magma_tally3_int_t m, magma_tally3_int_t i, 
    float **tau_array, 
    float **v_array, magma_tally3_int_t ldv, 
    float **T_array, magma_tally3_int_t ldt,
    magma_tally3_int_t batchCount, magma_tally3_queue_t queue);

void magma_tally3blas_slarft_gemvcolwise(
    magma_tally3_int_t m,  magma_tally3_int_t step,
    float *v, magma_tally3_int_t ldv, 
    float *T,  magma_tally3_int_t ldt,
    float *tau);

void magma_tally3blas_slarft_gemvcolwise_batched(
    magma_tally3_int_t m,  magma_tally3_int_t step,
    float **v_array, magma_tally3_int_t ldv, 
    float **T_array,  magma_tally3_int_t ldt,
    float **tau_array, magma_tally3_int_t batchCount, magma_tally3_queue_t queue);




void sgeqrf_copy_upper_batched(                
                  magma_tally3_int_t n, magma_tally3_int_t nb,
                  float **dV_array,    magma_tally3_int_t ldv,
                  float **dR_array,    magma_tally3_int_t ldr,
          magma_tally3_int_t batchCount, magma_tally3_queue_t queue);



void 
magma_tally3blas_snrm2_cols_batched(
    magma_tally3_int_t m, magma_tally3_int_t n,
    float **dA_array, magma_tally3_int_t lda, 
    float **dxnorm_array, magma_tally3_int_t batchCount);
 
void 
magma_tally3_slarfgx_batched(magma_tally3_int_t n, float **dx0_array, float **dx_array, 
                  float **dtau_array, float **dxnorm_array, 
                  float **dR_array, magma_tally3_int_t it, magma_tally3_int_t batchCount);


void 
magma_tally3_slarfx_batched_v4(magma_tally3_int_t m, magma_tally3_int_t n, float **v_array, float **tau_array,
                float **C_array, magma_tally3_int_t ldc, float **xnorm_array, 
                magma_tally3_int_t step, 
                magma_tally3_int_t batchCount);


void 
magma_tally3blas_slarfg_batched(
    magma_tally3_int_t n,
    float** dalpha_array, float** dx_array, magma_tally3_int_t incx,
    float** dtau_array, magma_tally3_int_t batchCount );





// for debugging purpose
void 
sset_stepinit_ipiv(magma_tally3_int_t **ipiv_array,
                 magma_tally3_int_t pm,
                 magma_tally3_int_t batchCount);



#ifdef __cplusplus
}
#endif

#undef REAL

#endif  /* MAGMA_tally3_SBATCHED_H */
