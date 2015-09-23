/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Azzam Haidar
       @author Tingxing Dong

       @precisions normal z -> s d c
*/

#ifndef MAGMA_minproduct_ZBATCHED_H
#define MAGMA_minproduct_ZBATCHED_H

#include "magma_minproduct_types.h"

#define COMPLEX

#ifdef __cplusplus
extern "C" {
#endif
  /*
   *  local auxiliary routines
   */
void 
zset_pointer_int(magma_minproduct_int_t **output_array,
        magma_minproduct_int_t *input,
        magma_minproduct_int_t lda,
        magma_minproduct_int_t row, magma_minproduct_int_t column, 
        magma_minproduct_int_t batchSize,
        magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

void 
zset_pointer(magma_minproductDoubleComplex **output_array,
                 magma_minproductDoubleComplex *input,
                 magma_minproduct_int_t lda,
                 magma_minproduct_int_t row, magma_minproduct_int_t column,
                 magma_minproduct_int_t batchSize,
                 magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);


void 
zset_array(magma_minproductDoubleComplex **output_array,
               magma_minproductDoubleComplex **input_array, magma_minproduct_int_t lda,
               magma_minproduct_int_t row, magma_minproduct_int_t column,
               magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

void 
magma_minproduct_zdisplace_pointers(magma_minproductDoubleComplex **output_array,
               magma_minproductDoubleComplex **input_array, magma_minproduct_int_t lda,
               magma_minproduct_int_t row, magma_minproduct_int_t column, 
               magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

  /*
   *  LAPACK batched routines
   */

  /*
   *  BLAS batched routines
   */

void
magma_minproductblas_zgemm_batched(
    magma_minproduct_trans_t transA, magma_minproduct_trans_t transB, magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductDoubleComplex alpha,
    magma_minproductDoubleComplex const * const * dA_array, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex const * const * dB_array, magma_minproduct_int_t lddb,
    magma_minproductDoubleComplex beta,
    magma_minproductDoubleComplex **dC_array, magma_minproduct_int_t lddc, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue );

void
magma_minproductblas_zgemm_batched_lg(
    magma_minproduct_trans_t transA, magma_minproduct_trans_t transB, magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductDoubleComplex alpha,
    magma_minproductDoubleComplex const * const * dA_array, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex const * const * dB_array, magma_minproduct_int_t lddb,
    magma_minproductDoubleComplex beta,
    magma_minproductDoubleComplex **dC_array, magma_minproduct_int_t lddc, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue );
void
magma_minproductblas_zgemm_batched_k32(
    magma_minproduct_trans_t transA, magma_minproduct_trans_t transB, magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductDoubleComplex alpha,
    magma_minproductDoubleComplex const * const * dA_array, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex const * const * dB_array, magma_minproduct_int_t lddb,
    magma_minproductDoubleComplex beta,
    magma_minproductDoubleComplex **dC_array, magma_minproduct_int_t lddc, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue );



void 
magma_minproductblas_zherk_NC_batched( magma_minproduct_trans_t TRANSA, magma_minproduct_trans_t TRANSB, int m , int n , int k , 
                       magma_minproductDoubleComplex alpha, magma_minproductDoubleComplex **dA_array, int lda, 
                       magma_minproductDoubleComplex **B_array, int ldb, 
                       magma_minproductDoubleComplex beta,        magma_minproductDoubleComplex **C_array, int ldc, 
                       magma_minproduct_int_t batchCount);

void
magma_minproductblas_zherk_batched(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans, magma_minproduct_int_t n, magma_minproduct_int_t k,
    double alpha,
    magma_minproductDoubleComplex const * const * dA_array, magma_minproduct_int_t ldda,
    double beta,
    magma_minproductDoubleComplex **dC_array, magma_minproduct_int_t lddc, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue );

void
magma_minproductblas_zherk_batched_lg(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans, magma_minproduct_int_t n, magma_minproduct_int_t k,
    double alpha,
    magma_minproductDoubleComplex const * const * dA_array, magma_minproduct_int_t ldda,
    double beta,
    magma_minproductDoubleComplex **dC_array, magma_minproduct_int_t lddc, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue );

void
magma_minproductblas_zherk_batched_k32(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans, magma_minproduct_int_t n, magma_minproduct_int_t k,
    double alpha,
    magma_minproductDoubleComplex const * const * dA_array, magma_minproduct_int_t ldda,
    double beta,
    magma_minproductDoubleComplex **dC_array, magma_minproduct_int_t lddc, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue );






magma_minproduct_int_t 
magma_minproduct_zpotf2_tile_batched(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex **dA_array, magma_minproduct_int_t lda,
    magma_minproduct_int_t *info_array, magma_minproduct_int_t gbstep, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_zpotrf_panel(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t nb,     
    magma_minproductDoubleComplex *A, magma_minproduct_int_t lda,
    magma_minproductDoubleComplex* dinvA, magma_minproductDoubleComplex** dinvA_array, magma_minproduct_int_t invA_msize,
    magma_minproductDoubleComplex* x, magma_minproductDoubleComplex** x_array,  magma_minproduct_int_t x_msize,
    magma_minproduct_int_t *info_array,  magma_minproduct_int_t batchCount, magma_minproduct_int_t matrixSize, magma_minproduct_queue_t queue);


void 
magma_minproductblas_ztrtri_diag_batched(
    magma_minproduct_uplo_t uplo, magma_minproduct_diag_t diag, magma_minproduct_int_t n,
    magma_minproductDoubleComplex const * const *dA_array, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex **dinvA_array, 
    magma_minproduct_int_t resetozero, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);



void 
magma_minproductblas_ztrsm_batched(
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo, magma_minproduct_trans_t transA, magma_minproduct_diag_t diag,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex alpha,
    magma_minproductDoubleComplex** dA_array,    magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex** dB_array,    magma_minproduct_int_t lddb,
    magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);


void 
magma_minproductblas_ztrsm_work_batched(
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo, magma_minproduct_trans_t transA, magma_minproduct_diag_t diag,
    magma_minproduct_int_t flag, magma_minproduct_int_t m, magma_minproduct_int_t n, 
    magma_minproductDoubleComplex alpha, 
    magma_minproductDoubleComplex** dA_array,    magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex** dB_array,    magma_minproduct_int_t lddb,
    magma_minproductDoubleComplex** dX_array,    magma_minproduct_int_t lddx, 
    magma_minproductDoubleComplex** dinvA_array, magma_minproduct_int_t dinvA_length,
    magma_minproductDoubleComplex** dA_displ, magma_minproductDoubleComplex** dB_displ, 
    magma_minproductDoubleComplex** dX_displ, magma_minproductDoubleComplex** dinvA_displ,
    magma_minproduct_int_t resetozero, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

void 
magma_minproductblas_ztrsm_outofplace_batched(
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo, magma_minproduct_trans_t transA, magma_minproduct_diag_t diag,
    magma_minproduct_int_t flag, magma_minproduct_int_t m, magma_minproduct_int_t n, 
    magma_minproductDoubleComplex alpha, 
    magma_minproductDoubleComplex** dA_array,    magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex** dB_array,    magma_minproduct_int_t lddb,
    magma_minproductDoubleComplex** dX_array,    magma_minproduct_int_t lddx, 
    magma_minproductDoubleComplex** dinvA_array, magma_minproduct_int_t dinvA_length,
    magma_minproductDoubleComplex** dA_displ, magma_minproductDoubleComplex** dB_displ, 
    magma_minproductDoubleComplex** dX_displ, magma_minproductDoubleComplex** dinvA_displ,
    magma_minproduct_int_t resetozero, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_zpotrf_batched(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductDoubleComplex **dA_array, magma_minproduct_int_t lda,
    magma_minproduct_int_t *info_array,  magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);


magma_minproduct_int_t 
magma_minproduct_zpotf2_batched(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex **dA_array, magma_minproduct_int_t lda,
    magma_minproductDoubleComplex **dA_displ, 
    magma_minproductDoubleComplex **dW_displ,
    magma_minproductDoubleComplex **dB_displ, 
    magma_minproductDoubleComplex **dC_displ, 
    magma_minproduct_int_t *info_array, magma_minproduct_int_t gbstep, 
    magma_minproduct_int_t batchCount, cublasHandle_t myhandle, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_zpotrf_panel_batched(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t nb,     
    magma_minproductDoubleComplex** dA_array,    magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex** dX_array,    magma_minproduct_int_t dX_length,
    magma_minproductDoubleComplex** dinvA_array, magma_minproduct_int_t dinvA_length,
    magma_minproductDoubleComplex** dW0_displ, magma_minproductDoubleComplex** dW1_displ, 
    magma_minproductDoubleComplex** dW2_displ, magma_minproductDoubleComplex** dW3_displ,
    magma_minproductDoubleComplex** dW4_displ, 
    magma_minproduct_int_t *info_array, magma_minproduct_int_t gbstep,
    magma_minproduct_int_t batchCount, cublasHandle_t myhandle, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_zpotrf_recpanel_batched(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n, 
    magma_minproduct_int_t min_recpnb,    
    magma_minproductDoubleComplex** dA_array,    magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex** dX_array,    magma_minproduct_int_t dX_length,
    magma_minproductDoubleComplex** dinvA_array, magma_minproduct_int_t dinvA_length,
    magma_minproductDoubleComplex** dW0_displ, magma_minproductDoubleComplex** dW1_displ,  
    magma_minproductDoubleComplex** dW2_displ, magma_minproductDoubleComplex** dW3_displ,
    magma_minproductDoubleComplex** dW4_displ,
    magma_minproduct_int_t *info_array, magma_minproduct_int_t gbstep, 
    magma_minproduct_int_t batchCount, cublasHandle_t myhandle, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_zpotrf_rectile_batched(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n, 
    magma_minproduct_int_t min_recpnb,    
    magma_minproductDoubleComplex** dA_array,    magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex** dX_array,    magma_minproduct_int_t dX_length,
    magma_minproductDoubleComplex** dinvA_array, magma_minproduct_int_t dinvA_length,
    magma_minproductDoubleComplex** dW0_displ, magma_minproductDoubleComplex** dW1_displ,  
    magma_minproductDoubleComplex** dW2_displ, magma_minproductDoubleComplex** dW3_displ,
    magma_minproductDoubleComplex** dW4_displ,
    magma_minproduct_int_t *info_array, magma_minproduct_int_t gbstep,
    magma_minproduct_int_t batchCount, cublasHandle_t myhandle, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_zpotrs_batched(
                  magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
                  magma_minproductDoubleComplex **dA_array, magma_minproduct_int_t ldda,
                  magma_minproductDoubleComplex **dB_array, magma_minproduct_int_t lddb,
                  magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_zposv_batched(
                  magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
                  magma_minproductDoubleComplex **dA_array, magma_minproduct_int_t ldda,
                  magma_minproductDoubleComplex **dB_array, magma_minproduct_int_t lddb,
                  magma_minproduct_int_t *dinfo_array,
                  magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_zgetrs_batched(
                  magma_minproduct_trans_t trans, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
                  magma_minproductDoubleComplex **dA_array, magma_minproduct_int_t ldda,
                  magma_minproduct_int_t **dipiv_array, 
                  magma_minproductDoubleComplex **dB_array, magma_minproduct_int_t lddb,
                  magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);



void 
magma_minproduct_zlaswp_rowparallel_batched( magma_minproduct_int_t n, magma_minproductDoubleComplex** input_array, magma_minproduct_int_t ldi,
                   magma_minproductDoubleComplex** output_array, magma_minproduct_int_t ldo,
                   magma_minproduct_int_t k1, magma_minproduct_int_t k2,
                   magma_minproduct_int_t **pivinfo_array, 
                   magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue );

void 
magma_minproduct_zlaswp_rowserial_batched(magma_minproduct_int_t n, magma_minproductDoubleComplex** dA_array, magma_minproduct_int_t lda,
                   magma_minproduct_int_t k1, magma_minproduct_int_t k2,
                   magma_minproduct_int_t **ipiv_array, 
                   magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

void 
magma_minproduct_zlaswp_columnserial_batched(magma_minproduct_int_t n, magma_minproductDoubleComplex** dA_array, magma_minproduct_int_t lda,
                   magma_minproduct_int_t k1, magma_minproduct_int_t k2,
                   magma_minproduct_int_t **ipiv_array, 
                   magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

void 
magma_minproductblas_ztranspose_batched(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex **dA_array,  magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex **dAT_array, magma_minproduct_int_t lddat, magma_minproduct_int_t batchCount );


void 
magma_minproductblas_zlaset_batched(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex offdiag, magma_minproductDoubleComplex diag,
    magma_minproductDoubleComplex_ptr dAarray[], magma_minproduct_int_t ldda,
    magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

void 
magma_minproductblas_zmemset_batched(magma_minproduct_int_t length, 
        magma_minproductDoubleComplex_ptr dAarray[], magma_minproductDoubleComplex val, 
        magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_zgetf2_batched(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex **dA_array, magma_minproduct_int_t lda,
    magma_minproductDoubleComplex **GERA_array,
    magma_minproductDoubleComplex **GERB_array,
    magma_minproductDoubleComplex **GERC_array,
    magma_minproduct_int_t **ipiv_array,
    magma_minproduct_int_t *info_array, 
    magma_minproduct_int_t gbstep,            
    magma_minproduct_int_t batchCount,
    cublasHandle_t myhandle, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
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
    magma_minproduct_int_t batchCount, cublasHandle_t myhandle, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_zgetrf_batched(
        magma_minproduct_int_t m, magma_minproduct_int_t n,
        magma_minproductDoubleComplex **dA_array, 
        magma_minproduct_int_t lda,
        magma_minproduct_int_t **ipiv_array, 
        magma_minproduct_int_t *info_array, 
        magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_zgetri_outofplace_batched( magma_minproduct_int_t n, 
                  magma_minproductDoubleComplex **dA_array, magma_minproduct_int_t ldda,
                  magma_minproduct_int_t **dipiv_array, 
                  magma_minproductDoubleComplex **dinvA_array, magma_minproduct_int_t lddia,
                  magma_minproduct_int_t *info_array,
                  magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

void 
magma_minproduct_zdisplace_intpointers(magma_minproduct_int_t **output_array,
               magma_minproduct_int_t **input_array, magma_minproduct_int_t lda,
               magma_minproduct_int_t row, magma_minproduct_int_t column, 
               magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);





void 
magma_minproductblas_izamax_atomic_batched(magma_minproduct_int_t n, magma_minproductDoubleComplex** x_array, magma_minproduct_int_t incx, magma_minproduct_int_t **max_id_array, magma_minproduct_int_t batchCount);

void 
magma_minproductblas_izamax_tree_batched(magma_minproduct_int_t n, magma_minproductDoubleComplex** x_array, magma_minproduct_int_t incx, magma_minproduct_int_t **max_id_array, magma_minproduct_int_t batchCount);



void 
magma_minproductblas_izamax_batched(magma_minproduct_int_t n, magma_minproductDoubleComplex** x_array, magma_minproduct_int_t incx, magma_minproduct_int_t **max_id_array, magma_minproduct_int_t batchCount);

void 
magma_minproductblas_izamax(magma_minproduct_int_t n, magma_minproductDoubleComplex* x, magma_minproduct_int_t incx, magma_minproduct_int_t *max_id);


magma_minproduct_int_t 
magma_minproduct_izamax_batched(magma_minproduct_int_t length, 
        magma_minproductDoubleComplex **x_array, magma_minproduct_int_t incx, magma_minproduct_int_t step,  magma_minproduct_int_t lda,
        magma_minproduct_int_t** ipiv_array, magma_minproduct_int_t *info_array, magma_minproduct_int_t gbstep, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_zswap_batched(magma_minproduct_int_t n, magma_minproductDoubleComplex **x_array, magma_minproduct_int_t incx, magma_minproduct_int_t j, 
                 magma_minproduct_int_t** ipiv_array, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_zscal_zgeru_batched(magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t step,
                                      magma_minproductDoubleComplex **dA_array, magma_minproduct_int_t lda,
                                      magma_minproduct_int_t *info_array, magma_minproduct_int_t gbstep, 
                                      magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_zcomputecolumn_batched(magma_minproduct_int_t m, magma_minproduct_int_t paneloffset, magma_minproduct_int_t step, 
                                        magma_minproductDoubleComplex **dA_array,  magma_minproduct_int_t lda,
                                        magma_minproduct_int_t **ipiv_array, 
                                        magma_minproduct_int_t *info_array, magma_minproduct_int_t gbstep, 
                                        magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

void 
magma_minproduct_zgetf2trsm_batched(magma_minproduct_int_t ib, magma_minproduct_int_t n, magma_minproductDoubleComplex **dA_array,  magma_minproduct_int_t j, magma_minproduct_int_t lda,
                       magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);


magma_minproduct_int_t 
magma_minproduct_zgetf2_nopiv_batched(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex **dA_array, magma_minproduct_int_t lda,
    magma_minproductDoubleComplex **dW0_displ,
    magma_minproductDoubleComplex **dW1_displ,
    magma_minproductDoubleComplex **dW2_displ,
    magma_minproduct_int_t *info_array,            
    magma_minproduct_int_t gbstep, 
    magma_minproduct_int_t batchCount,
    cublasHandle_t myhandle, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_zgetrf_panel_nopiv_batched(
    magma_minproduct_int_t m, magma_minproduct_int_t nb,    
    magma_minproductDoubleComplex** dA_array,    magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex** dX_array,    magma_minproduct_int_t dX_length,
    magma_minproductDoubleComplex** dinvA_array, magma_minproduct_int_t dinvA_length,
    magma_minproductDoubleComplex** dW0_displ, magma_minproductDoubleComplex** dW1_displ,  
    magma_minproductDoubleComplex** dW2_displ, magma_minproductDoubleComplex** dW3_displ,
    magma_minproductDoubleComplex** dW4_displ,     
    magma_minproduct_int_t *info_array, magma_minproduct_int_t gbstep, 
    magma_minproduct_int_t batchCount, cublasHandle_t myhandle, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_zgetrf_recpanel_nopiv_batched(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t min_recpnb,    
    magma_minproductDoubleComplex** dA_array,    magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex** dX_array,    magma_minproduct_int_t dX_length,
    magma_minproductDoubleComplex** dinvA_array, magma_minproduct_int_t dinvA_length,
    magma_minproductDoubleComplex** dW1_displ, magma_minproductDoubleComplex** dW2_displ,  
    magma_minproductDoubleComplex** dW3_displ, magma_minproductDoubleComplex** dW4_displ,
    magma_minproductDoubleComplex** dW5_displ, 
    magma_minproduct_int_t *info_array, magma_minproduct_int_t gbstep,
    magma_minproduct_int_t batchCount, cublasHandle_t myhandle, magma_minproduct_queue_t queue);


magma_minproduct_int_t 
magma_minproduct_zgetrf_nopiv_batched(
        magma_minproduct_int_t m, magma_minproduct_int_t n,
        magma_minproductDoubleComplex **dA_array, 
        magma_minproduct_int_t lda,
        magma_minproduct_int_t *info_array, 
        magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_zgetrs_nopiv_batched(
                  magma_minproduct_trans_t trans, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
                  magma_minproductDoubleComplex **dA_array, magma_minproduct_int_t ldda,
                  magma_minproductDoubleComplex **dB_array, magma_minproduct_int_t lddb,
                  magma_minproduct_int_t *info_array,
                  magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_zgesv_nopiv_batched(
                  magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
                  magma_minproductDoubleComplex **dA_array, magma_minproduct_int_t ldda,
                  magma_minproductDoubleComplex **dB_array, magma_minproduct_int_t lddb,
                  magma_minproduct_int_t *info_array,
                  magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_zgesv_rbt_batched(
                  magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
                  magma_minproductDoubleComplex **dA_array, magma_minproduct_int_t ldda,
                  magma_minproductDoubleComplex **dB_array, magma_minproduct_int_t lddb,
                  magma_minproduct_int_t *info_array,
                  magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_zgesv_batched(
                  magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
                  magma_minproductDoubleComplex **dA_array, magma_minproduct_int_t ldda,
                  magma_minproduct_int_t **dipiv_array, 
                  magma_minproductDoubleComplex **dB_array, magma_minproduct_int_t lddb,
                  magma_minproduct_int_t *dinfo_array,
                  magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

magma_minproduct_int_t
magma_minproduct_zgerbt_batched(
    magma_minproduct_bool_t gen, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
    magma_minproductDoubleComplex **dA_array, magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex **dB_array, magma_minproduct_int_t lddb,
    magma_minproductDoubleComplex *U, magma_minproductDoubleComplex *V,
    magma_minproduct_int_t *info, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

void 
magma_minproductblas_zprbt_batched(
    magma_minproduct_int_t n, 
    magma_minproductDoubleComplex **dA_array, magma_minproduct_int_t ldda, 
    magma_minproductDoubleComplex *du, magma_minproductDoubleComplex *dv,
    magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

void
magma_minproductblas_zprbt_mv_batched(
    magma_minproduct_int_t n, 
    magma_minproductDoubleComplex *dv, magma_minproductDoubleComplex **db_array, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);


void
magma_minproductblas_zprbt_mtv_batched(
    magma_minproduct_int_t n, 
    magma_minproductDoubleComplex *du, magma_minproductDoubleComplex **db_array, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);





void 
magma_minproduct_zlacgv_batched(magma_minproduct_int_t n, magma_minproductDoubleComplex **x_array, magma_minproduct_int_t incx, int offset, int batchCount, magma_minproduct_queue_t queue);

void 
magma_minproduct_zpotf2_zdscal_batched(magma_minproduct_int_t n, magma_minproductDoubleComplex **x_array, magma_minproduct_int_t incx, magma_minproduct_int_t offset, magma_minproduct_int_t *info_array, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

void 
magma_minproduct_zpotf2_zdotc_batched(magma_minproduct_int_t n, magma_minproductDoubleComplex **x_array, magma_minproduct_int_t incx, magma_minproduct_int_t offset, magma_minproduct_int_t *info_array, magma_minproduct_int_t gbstep, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);


void 
setup_pivinfo( magma_minproduct_int_t *pivinfo, magma_minproduct_int_t *ipiv, 
                      magma_minproduct_int_t m, magma_minproduct_int_t nb, 
                      magma_minproduct_queue_t queue);


void
magma_minproductblas_zgeadd_batched_q(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex alpha,
    magma_minproductDoubleComplex_const_ptr  const dAarray[], magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex_ptr              dBarray[], magma_minproduct_int_t lddb,
    magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue );

void
magma_minproductblas_zlacpy_batched(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex_const_ptr  const dAarray[], magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex_ptr              dBarray[], magma_minproduct_int_t lddb,
    magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue );

void
magma_minproductblas_zgeadd_batched(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex alpha,
    magma_minproductDoubleComplex_const_ptr  const dAarray[], magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex_ptr              dBarray[], magma_minproduct_int_t lddb,
    magma_minproduct_int_t batchCount );


void
magma_minproductblas_zgemv_batched(
    magma_minproduct_trans_t trans, magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex alpha,
    magma_minproductDoubleComplex_ptr dA_array[], magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex_ptr dx_array[], magma_minproduct_int_t incx,
    magma_minproductDoubleComplex beta,
    magma_minproductDoubleComplex_ptr dy_array[], magma_minproduct_int_t incy,
    magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);


magma_minproduct_int_t 
magma_minproduct_zgeqrf_batched(
    magma_minproduct_int_t m, magma_minproduct_int_t n, 
    magma_minproductDoubleComplex **dA_array,
    magma_minproduct_int_t lda, 
    magma_minproductDoubleComplex **tau_array,
    magma_minproduct_int_t *info_array, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_zgeqrf_batched_v4(
    magma_minproduct_int_t m, magma_minproduct_int_t n, 
    magma_minproductDoubleComplex **dA_array,
    magma_minproduct_int_t lda, 
    magma_minproductDoubleComplex **tau_array,
    magma_minproduct_int_t *info_array, magma_minproduct_int_t batchCount);

magma_minproduct_int_t
magma_minproduct_zgeqrf_panel_batched(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t nb,    
    magma_minproductDoubleComplex** dA_array,    magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex** tau_array, 
    magma_minproductDoubleComplex** dT_array, magma_minproduct_int_t ldt, 
    magma_minproductDoubleComplex** dR_array, magma_minproduct_int_t ldr,
    magma_minproductDoubleComplex** dW0_displ, 
    magma_minproductDoubleComplex** dW1_displ,
    magma_minproductDoubleComplex *dwork,  
    magma_minproductDoubleComplex** W_array, 
    magma_minproductDoubleComplex** W2_array,
    magma_minproduct_int_t *info_array,
    magma_minproduct_int_t batchCount, cublasHandle_t myhandle, magma_minproduct_queue_t queue);

magma_minproduct_int_t
magma_minproduct_zgeqrf_panel_batched_v4(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t nb,    
    magma_minproductDoubleComplex** dA_array,    magma_minproduct_int_t ldda,
    magma_minproductDoubleComplex** tau_array, 
    magma_minproductDoubleComplex** dT_array, magma_minproduct_int_t ldt, 
    magma_minproductDoubleComplex** dR_array, magma_minproduct_int_t ldr,
    double** dnorm_array,  
    magma_minproductDoubleComplex** dW0_displ, 
    magma_minproductDoubleComplex** dW1_displ,
    magma_minproductDoubleComplex *dwork,  
    magma_minproductDoubleComplex** W_array, 
    magma_minproductDoubleComplex** W2_array,
    magma_minproduct_int_t *info_array,
    magma_minproduct_int_t batchCount, cublasHandle_t myhandle);

magma_minproduct_int_t
magma_minproduct_zgeqr2x_batched_v4(magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproductDoubleComplex **dA_array,
                  magma_minproduct_int_t lda, 
                  magma_minproductDoubleComplex **tau_array,
                  magma_minproductDoubleComplex **dR_array, magma_minproduct_int_t ldr,
                  double **dwork_array,  
                  magma_minproduct_int_t *info, magma_minproduct_int_t batchCount);

magma_minproduct_int_t
magma_minproduct_zgeqr2_batched(magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproductDoubleComplex **dA_array,
                  magma_minproduct_int_t lda, 
                  magma_minproductDoubleComplex **tau_array,
                  magma_minproduct_int_t *info, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

magma_minproduct_int_t
magma_minproduct_zlarfb_zgemm_batched(
                  cublasHandle_t myhandle,
                  magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
                  magma_minproductDoubleComplex **dV_array,    magma_minproduct_int_t ldv,
                  magma_minproductDoubleComplex **dT_array,    magma_minproduct_int_t ldt,
                  magma_minproductDoubleComplex **dA_array,    magma_minproduct_int_t lda,
                  magma_minproductDoubleComplex **W_array,     magma_minproduct_int_t ldw,
                  magma_minproductDoubleComplex **W2_array,    magma_minproduct_int_t ldw2,
                  magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

magma_minproduct_int_t
magma_minproduct_zlarfb_gemm_batched(
    magma_minproduct_side_t side, magma_minproduct_trans_t trans, magma_minproduct_direct_t direct, magma_minproduct_storev_t storev,
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductDoubleComplex_const_ptr dV_array[],    magma_minproduct_int_t lddv,
    magma_minproductDoubleComplex_const_ptr dT_array[],    magma_minproduct_int_t lddt,
    magma_minproductDoubleComplex_ptr dC_array[],          magma_minproduct_int_t lddc,
    magma_minproductDoubleComplex_ptr dwork_array[],       magma_minproduct_int_t ldwork,
    magma_minproductDoubleComplex_ptr dworkvt_array[],     magma_minproduct_int_t ldworkvt,
    magma_minproduct_int_t batchCount, cublasHandle_t myhandle, magma_minproduct_queue_t queue);




void
magma_minproduct_zlarft_batched_vold(magma_minproduct_int_t n, magma_minproduct_int_t k, magma_minproductDoubleComplex **v_array, magma_minproduct_int_t ldv,
                    magma_minproductDoubleComplex **tau_array,
                    magma_minproductDoubleComplex **T_array, magma_minproduct_int_t ldt, 
                    magma_minproduct_int_t batchCount);





magma_minproduct_int_t
magma_minproduct_zlarft_batched(magma_minproduct_int_t n, magma_minproduct_int_t k, magma_minproduct_int_t stair_T, 
                magma_minproductDoubleComplex **v_array, magma_minproduct_int_t ldv,
                magma_minproductDoubleComplex **tau_array, magma_minproductDoubleComplex **T_array, magma_minproduct_int_t ldt, 
                magma_minproductDoubleComplex **work_array, magma_minproduct_int_t lwork, magma_minproduct_int_t batchCount, cublasHandle_t myhandle, magma_minproduct_queue_t queue);

void
magma_minproduct_zlarft_sm32x32_batched(magma_minproduct_int_t n, magma_minproduct_int_t k, magma_minproductDoubleComplex **v_array, magma_minproduct_int_t ldv,
                    magma_minproductDoubleComplex **tau_array, magma_minproductDoubleComplex **T_array, magma_minproduct_int_t ldt, 
                    magma_minproduct_int_t batchCount, cublasHandle_t myhandle, magma_minproduct_queue_t queue);



void magma_minproductblas_zlarft_recztrmv_sm32x32(
    magma_minproduct_int_t m, magma_minproduct_int_t n, 
    magma_minproductDoubleComplex *tau, 
    magma_minproductDoubleComplex *Trec, magma_minproduct_int_t ldtrec, 
    magma_minproductDoubleComplex *Ttri, magma_minproduct_int_t ldttri);


void magma_minproductblas_zlarft_recztrmv_sm32x32_batched(
    magma_minproduct_int_t m, magma_minproduct_int_t n, 
    magma_minproductDoubleComplex **tau_array, 
    magma_minproductDoubleComplex **Trec_array, magma_minproduct_int_t ldtrec, 
    magma_minproductDoubleComplex **Ttri_array, magma_minproduct_int_t ldttri,
    magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

void magma_minproductblas_zlarft_ztrmv_sm32x32(
    magma_minproduct_int_t m, magma_minproduct_int_t n, 
    magma_minproductDoubleComplex *tau, 
    magma_minproductDoubleComplex *Tin, magma_minproduct_int_t ldtin, 
    magma_minproductDoubleComplex *Tout, magma_minproduct_int_t ldtout);

void magma_minproductblas_zlarft_ztrmv_sm32x32_batched(
    magma_minproduct_int_t m, magma_minproduct_int_t n, 
    magma_minproductDoubleComplex **tau_array, 
    magma_minproductDoubleComplex **Tin_array, magma_minproduct_int_t ldtin, 
    magma_minproductDoubleComplex **Tout_array, magma_minproduct_int_t ldtout,
    magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

void magma_minproductblas_zlarft_gemv_loop_inside(
    int n, int k, 
    magma_minproductDoubleComplex *tau, 
    magma_minproductDoubleComplex *v, int ldv, 
    magma_minproductDoubleComplex *T, int ldt);

void magma_minproductblas_zlarft_gemv_loop_inside_batched(
    int n, int k, 
    magma_minproductDoubleComplex **tau_array, 
    magma_minproductDoubleComplex **v_array, int ldv, 
    magma_minproductDoubleComplex **T_array, int ldt, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

void magma_minproductblas_zlarft_gemvrowwise(
    magma_minproduct_int_t m, magma_minproduct_int_t i, 
    magma_minproductDoubleComplex *tau, 
    magma_minproductDoubleComplex *v, magma_minproduct_int_t ldv, 
    magma_minproductDoubleComplex *T, magma_minproduct_int_t ldt,
    magma_minproductDoubleComplex *W);

void magma_minproductblas_zlarft_gemvrowwise_batched(
    magma_minproduct_int_t m, magma_minproduct_int_t i, 
    magma_minproductDoubleComplex **tau_array, 
    magma_minproductDoubleComplex **v_array, magma_minproduct_int_t ldv, 
    magma_minproductDoubleComplex **T_array, magma_minproduct_int_t ldt,
    magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

void magma_minproductblas_zlarft_gemvcolwise(
    magma_minproduct_int_t m,  magma_minproduct_int_t step,
    magma_minproductDoubleComplex *v, magma_minproduct_int_t ldv, 
    magma_minproductDoubleComplex *T,  magma_minproduct_int_t ldt,
    magma_minproductDoubleComplex *tau);

void magma_minproductblas_zlarft_gemvcolwise_batched(
    magma_minproduct_int_t m,  magma_minproduct_int_t step,
    magma_minproductDoubleComplex **v_array, magma_minproduct_int_t ldv, 
    magma_minproductDoubleComplex **T_array,  magma_minproduct_int_t ldt,
    magma_minproductDoubleComplex **tau_array, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);




void zgeqrf_copy_upper_batched(                
                  magma_minproduct_int_t n, magma_minproduct_int_t nb,
                  magma_minproductDoubleComplex **dV_array,    magma_minproduct_int_t ldv,
                  magma_minproductDoubleComplex **dR_array,    magma_minproduct_int_t ldr,
          magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);



void 
magma_minproductblas_dznrm2_cols_batched(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex **dA_array, magma_minproduct_int_t lda, 
    double **dxnorm_array, magma_minproduct_int_t batchCount);
 
void 
magma_minproduct_zlarfgx_batched(magma_minproduct_int_t n, magma_minproductDoubleComplex **dx0_array, magma_minproductDoubleComplex **dx_array, 
                  magma_minproductDoubleComplex **dtau_array, double **dxnorm_array, 
                  magma_minproductDoubleComplex **dR_array, magma_minproduct_int_t it, magma_minproduct_int_t batchCount);


void 
magma_minproduct_zlarfx_batched_v4(magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproductDoubleComplex **v_array, magma_minproductDoubleComplex **tau_array,
                magma_minproductDoubleComplex **C_array, magma_minproduct_int_t ldc, double **xnorm_array, 
                magma_minproduct_int_t step, 
                magma_minproduct_int_t batchCount);


void 
magma_minproductblas_zlarfg_batched(
    magma_minproduct_int_t n,
    magma_minproductDoubleComplex** dalpha_array, magma_minproductDoubleComplex** dx_array, magma_minproduct_int_t incx,
    magma_minproductDoubleComplex** dtau_array, magma_minproduct_int_t batchCount );





// for debugging purpose
void 
zset_stepinit_ipiv(magma_minproduct_int_t **ipiv_array,
                 magma_minproduct_int_t pm,
                 magma_minproduct_int_t batchCount);



#ifdef __cplusplus
}
#endif

#undef COMPLEX

#endif  /* MAGMA_minproduct_ZBATCHED_H */
