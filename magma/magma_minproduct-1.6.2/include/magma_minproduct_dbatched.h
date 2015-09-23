/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Azzam Haidar
       @author Tingxing Dong

       @generated from magma_minproduct_zbatched.h normal z -> d, Fri Jan 30 19:00:06 2015
*/

#ifndef MAGMA_minproduct_DBATCHED_H
#define MAGMA_minproduct_DBATCHED_H

#include "magma_minproduct_types.h"

#define REAL

#ifdef __cplusplus
extern "C" {
#endif
  /*
   *  local auxiliary routines
   */
void 
dset_pointer_int(magma_minproduct_int_t **output_array,
        magma_minproduct_int_t *input,
        magma_minproduct_int_t lda,
        magma_minproduct_int_t row, magma_minproduct_int_t column, 
        magma_minproduct_int_t batchSize,
        magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

void 
dset_pointer(double **output_array,
                 double *input,
                 magma_minproduct_int_t lda,
                 magma_minproduct_int_t row, magma_minproduct_int_t column,
                 magma_minproduct_int_t batchSize,
                 magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);


void 
dset_array(double **output_array,
               double **input_array, magma_minproduct_int_t lda,
               magma_minproduct_int_t row, magma_minproduct_int_t column,
               magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

void 
magma_minproduct_ddisplace_pointers(double **output_array,
               double **input_array, magma_minproduct_int_t lda,
               magma_minproduct_int_t row, magma_minproduct_int_t column, 
               magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

  /*
   *  LAPACK batched routines
   */

  /*
   *  BLAS batched routines
   */

void
magma_minproductblas_dgemm_batched(
    magma_minproduct_trans_t transA, magma_minproduct_trans_t transB, magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    double alpha,
    double const * const * dA_array, magma_minproduct_int_t ldda,
    double const * const * dB_array, magma_minproduct_int_t lddb,
    double beta,
    double **dC_array, magma_minproduct_int_t lddc, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue );

void
magma_minproductblas_dgemm_batched_lg(
    magma_minproduct_trans_t transA, magma_minproduct_trans_t transB, magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    double alpha,
    double const * const * dA_array, magma_minproduct_int_t ldda,
    double const * const * dB_array, magma_minproduct_int_t lddb,
    double beta,
    double **dC_array, magma_minproduct_int_t lddc, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue );
void
magma_minproductblas_dgemm_batched_k32(
    magma_minproduct_trans_t transA, magma_minproduct_trans_t transB, magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    double alpha,
    double const * const * dA_array, magma_minproduct_int_t ldda,
    double const * const * dB_array, magma_minproduct_int_t lddb,
    double beta,
    double **dC_array, magma_minproduct_int_t lddc, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue );



void 
magma_minproductblas_dsyrk_NC_batched( magma_minproduct_trans_t TRANSA, magma_minproduct_trans_t TRANSB, int m , int n , int k , 
                       double alpha, double **dA_array, int lda, 
                       double **B_array, int ldb, 
                       double beta,        double **C_array, int ldc, 
                       magma_minproduct_int_t batchCount);

void
magma_minproductblas_dsyrk_batched(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans, magma_minproduct_int_t n, magma_minproduct_int_t k,
    double alpha,
    double const * const * dA_array, magma_minproduct_int_t ldda,
    double beta,
    double **dC_array, magma_minproduct_int_t lddc, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue );

void
magma_minproductblas_dsyrk_batched_lg(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans, magma_minproduct_int_t n, magma_minproduct_int_t k,
    double alpha,
    double const * const * dA_array, magma_minproduct_int_t ldda,
    double beta,
    double **dC_array, magma_minproduct_int_t lddc, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue );

void
magma_minproductblas_dsyrk_batched_k32(
    magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans, magma_minproduct_int_t n, magma_minproduct_int_t k,
    double alpha,
    double const * const * dA_array, magma_minproduct_int_t ldda,
    double beta,
    double **dC_array, magma_minproduct_int_t lddc, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue );






magma_minproduct_int_t 
magma_minproduct_dpotf2_tile_batched(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n,
    double **dA_array, magma_minproduct_int_t lda,
    magma_minproduct_int_t *info_array, magma_minproduct_int_t gbstep, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_dpotrf_panel(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t nb,     
    double *A, magma_minproduct_int_t lda,
    double* dinvA, double** dinvA_array, magma_minproduct_int_t invA_msize,
    double* x, double** x_array,  magma_minproduct_int_t x_msize,
    magma_minproduct_int_t *info_array,  magma_minproduct_int_t batchCount, magma_minproduct_int_t matrixSize, magma_minproduct_queue_t queue);


void 
magma_minproductblas_dtrtri_diag_batched(
    magma_minproduct_uplo_t uplo, magma_minproduct_diag_t diag, magma_minproduct_int_t n,
    double const * const *dA_array, magma_minproduct_int_t ldda,
    double **dinvA_array, 
    magma_minproduct_int_t resetozero, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);



void 
magma_minproductblas_dtrsm_batched(
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo, magma_minproduct_trans_t transA, magma_minproduct_diag_t diag,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    double alpha,
    double** dA_array,    magma_minproduct_int_t ldda,
    double** dB_array,    magma_minproduct_int_t lddb,
    magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);


void 
magma_minproductblas_dtrsm_work_batched(
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo, magma_minproduct_trans_t transA, magma_minproduct_diag_t diag,
    magma_minproduct_int_t flag, magma_minproduct_int_t m, magma_minproduct_int_t n, 
    double alpha, 
    double** dA_array,    magma_minproduct_int_t ldda,
    double** dB_array,    magma_minproduct_int_t lddb,
    double** dX_array,    magma_minproduct_int_t lddx, 
    double** dinvA_array, magma_minproduct_int_t dinvA_length,
    double** dA_displ, double** dB_displ, 
    double** dX_displ, double** dinvA_displ,
    magma_minproduct_int_t resetozero, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

void 
magma_minproductblas_dtrsm_outofplace_batched(
    magma_minproduct_side_t side, magma_minproduct_uplo_t uplo, magma_minproduct_trans_t transA, magma_minproduct_diag_t diag,
    magma_minproduct_int_t flag, magma_minproduct_int_t m, magma_minproduct_int_t n, 
    double alpha, 
    double** dA_array,    magma_minproduct_int_t ldda,
    double** dB_array,    magma_minproduct_int_t lddb,
    double** dX_array,    magma_minproduct_int_t lddx, 
    double** dinvA_array, magma_minproduct_int_t dinvA_length,
    double** dA_displ, double** dB_displ, 
    double** dX_displ, double** dinvA_displ,
    magma_minproduct_int_t resetozero, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_dpotrf_batched(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    double **dA_array, magma_minproduct_int_t lda,
    magma_minproduct_int_t *info_array,  magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);


magma_minproduct_int_t 
magma_minproduct_dpotf2_batched(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n,
    double **dA_array, magma_minproduct_int_t lda,
    double **dA_displ, 
    double **dW_displ,
    double **dB_displ, 
    double **dC_displ, 
    magma_minproduct_int_t *info_array, magma_minproduct_int_t gbstep, 
    magma_minproduct_int_t batchCount, cublasHandle_t myhandle, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_dpotrf_panel_batched(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t nb,     
    double** dA_array,    magma_minproduct_int_t ldda,
    double** dX_array,    magma_minproduct_int_t dX_length,
    double** dinvA_array, magma_minproduct_int_t dinvA_length,
    double** dW0_displ, double** dW1_displ, 
    double** dW2_displ, double** dW3_displ,
    double** dW4_displ, 
    magma_minproduct_int_t *info_array, magma_minproduct_int_t gbstep,
    magma_minproduct_int_t batchCount, cublasHandle_t myhandle, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_dpotrf_recpanel_batched(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n, 
    magma_minproduct_int_t min_recpnb,    
    double** dA_array,    magma_minproduct_int_t ldda,
    double** dX_array,    magma_minproduct_int_t dX_length,
    double** dinvA_array, magma_minproduct_int_t dinvA_length,
    double** dW0_displ, double** dW1_displ,  
    double** dW2_displ, double** dW3_displ,
    double** dW4_displ,
    magma_minproduct_int_t *info_array, magma_minproduct_int_t gbstep, 
    magma_minproduct_int_t batchCount, cublasHandle_t myhandle, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_dpotrf_rectile_batched(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n, 
    magma_minproduct_int_t min_recpnb,    
    double** dA_array,    magma_minproduct_int_t ldda,
    double** dX_array,    magma_minproduct_int_t dX_length,
    double** dinvA_array, magma_minproduct_int_t dinvA_length,
    double** dW0_displ, double** dW1_displ,  
    double** dW2_displ, double** dW3_displ,
    double** dW4_displ,
    magma_minproduct_int_t *info_array, magma_minproduct_int_t gbstep,
    magma_minproduct_int_t batchCount, cublasHandle_t myhandle, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_dpotrs_batched(
                  magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
                  double **dA_array, magma_minproduct_int_t ldda,
                  double **dB_array, magma_minproduct_int_t lddb,
                  magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_dposv_batched(
                  magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
                  double **dA_array, magma_minproduct_int_t ldda,
                  double **dB_array, magma_minproduct_int_t lddb,
                  magma_minproduct_int_t *dinfo_array,
                  magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_dgetrs_batched(
                  magma_minproduct_trans_t trans, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
                  double **dA_array, magma_minproduct_int_t ldda,
                  magma_minproduct_int_t **dipiv_array, 
                  double **dB_array, magma_minproduct_int_t lddb,
                  magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);



void 
magma_minproduct_dlaswp_rowparallel_batched( magma_minproduct_int_t n, double** input_array, magma_minproduct_int_t ldi,
                   double** output_array, magma_minproduct_int_t ldo,
                   magma_minproduct_int_t k1, magma_minproduct_int_t k2,
                   magma_minproduct_int_t **pivinfo_array, 
                   magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue );

void 
magma_minproduct_dlaswp_rowserial_batched(magma_minproduct_int_t n, double** dA_array, magma_minproduct_int_t lda,
                   magma_minproduct_int_t k1, magma_minproduct_int_t k2,
                   magma_minproduct_int_t **ipiv_array, 
                   magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

void 
magma_minproduct_dlaswp_columnserial_batched(magma_minproduct_int_t n, double** dA_array, magma_minproduct_int_t lda,
                   magma_minproduct_int_t k1, magma_minproduct_int_t k2,
                   magma_minproduct_int_t **ipiv_array, 
                   magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

void 
magma_minproductblas_dtranspose_batched(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    double **dA_array,  magma_minproduct_int_t ldda,
    double **dAT_array, magma_minproduct_int_t lddat, magma_minproduct_int_t batchCount );


void 
magma_minproductblas_dlaset_batched(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n,
    double offdiag, double diag,
    magma_minproductDouble_ptr dAarray[], magma_minproduct_int_t ldda,
    magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

void 
magma_minproductblas_dmemset_batched(magma_minproduct_int_t length, 
        magma_minproductDouble_ptr dAarray[], double val, 
        magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_dgetf2_batched(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    double **dA_array, magma_minproduct_int_t lda,
    double **GERA_array,
    double **GERB_array,
    double **GERC_array,
    magma_minproduct_int_t **ipiv_array,
    magma_minproduct_int_t *info_array, 
    magma_minproduct_int_t gbstep,            
    magma_minproduct_int_t batchCount,
    cublasHandle_t myhandle, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_dgetrf_recpanel_batched(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t min_recpnb,    
    double** dA_array,    magma_minproduct_int_t ldda,
    magma_minproduct_int_t** dipiv_array, magma_minproduct_int_t** dpivinfo_array,
    double** dX_array,    magma_minproduct_int_t dX_length,
    double** dinvA_array, magma_minproduct_int_t dinvA_length,
    double** dW1_displ, double** dW2_displ,  
    double** dW3_displ, double** dW4_displ,
    double** dW5_displ,
    magma_minproduct_int_t *info_array, magma_minproduct_int_t gbstep,
    magma_minproduct_int_t batchCount, cublasHandle_t myhandle, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_dgetrf_batched(
        magma_minproduct_int_t m, magma_minproduct_int_t n,
        double **dA_array, 
        magma_minproduct_int_t lda,
        magma_minproduct_int_t **ipiv_array, 
        magma_minproduct_int_t *info_array, 
        magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_dgetri_outofplace_batched( magma_minproduct_int_t n, 
                  double **dA_array, magma_minproduct_int_t ldda,
                  magma_minproduct_int_t **dipiv_array, 
                  double **dinvA_array, magma_minproduct_int_t lddia,
                  magma_minproduct_int_t *info_array,
                  magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

void 
magma_minproduct_ddisplace_intpointers(magma_minproduct_int_t **output_array,
               magma_minproduct_int_t **input_array, magma_minproduct_int_t lda,
               magma_minproduct_int_t row, magma_minproduct_int_t column, 
               magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);





void 
magma_minproductblas_idamax_atomic_batched(magma_minproduct_int_t n, double** x_array, magma_minproduct_int_t incx, magma_minproduct_int_t **max_id_array, magma_minproduct_int_t batchCount);

void 
magma_minproductblas_idamax_tree_batched(magma_minproduct_int_t n, double** x_array, magma_minproduct_int_t incx, magma_minproduct_int_t **max_id_array, magma_minproduct_int_t batchCount);



void 
magma_minproductblas_idamax_batched(magma_minproduct_int_t n, double** x_array, magma_minproduct_int_t incx, magma_minproduct_int_t **max_id_array, magma_minproduct_int_t batchCount);

void 
magma_minproductblas_idamax(magma_minproduct_int_t n, double* x, magma_minproduct_int_t incx, magma_minproduct_int_t *max_id);


magma_minproduct_int_t 
magma_minproduct_idamax_batched(magma_minproduct_int_t length, 
        double **x_array, magma_minproduct_int_t incx, magma_minproduct_int_t step,  magma_minproduct_int_t lda,
        magma_minproduct_int_t** ipiv_array, magma_minproduct_int_t *info_array, magma_minproduct_int_t gbstep, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_dswap_batched(magma_minproduct_int_t n, double **x_array, magma_minproduct_int_t incx, magma_minproduct_int_t j, 
                 magma_minproduct_int_t** ipiv_array, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_dscal_dger_batched(magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t step,
                                      double **dA_array, magma_minproduct_int_t lda,
                                      magma_minproduct_int_t *info_array, magma_minproduct_int_t gbstep, 
                                      magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_dcomputecolumn_batched(magma_minproduct_int_t m, magma_minproduct_int_t paneloffset, magma_minproduct_int_t step, 
                                        double **dA_array,  magma_minproduct_int_t lda,
                                        magma_minproduct_int_t **ipiv_array, 
                                        magma_minproduct_int_t *info_array, magma_minproduct_int_t gbstep, 
                                        magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

void 
magma_minproduct_dgetf2trsm_batched(magma_minproduct_int_t ib, magma_minproduct_int_t n, double **dA_array,  magma_minproduct_int_t j, magma_minproduct_int_t lda,
                       magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);


magma_minproduct_int_t 
magma_minproduct_dgetf2_nopiv_batched(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    double **dA_array, magma_minproduct_int_t lda,
    double **dW0_displ,
    double **dW1_displ,
    double **dW2_displ,
    magma_minproduct_int_t *info_array,            
    magma_minproduct_int_t gbstep, 
    magma_minproduct_int_t batchCount,
    cublasHandle_t myhandle, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_dgetrf_panel_nopiv_batched(
    magma_minproduct_int_t m, magma_minproduct_int_t nb,    
    double** dA_array,    magma_minproduct_int_t ldda,
    double** dX_array,    magma_minproduct_int_t dX_length,
    double** dinvA_array, magma_minproduct_int_t dinvA_length,
    double** dW0_displ, double** dW1_displ,  
    double** dW2_displ, double** dW3_displ,
    double** dW4_displ,     
    magma_minproduct_int_t *info_array, magma_minproduct_int_t gbstep, 
    magma_minproduct_int_t batchCount, cublasHandle_t myhandle, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_dgetrf_recpanel_nopiv_batched(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t min_recpnb,    
    double** dA_array,    magma_minproduct_int_t ldda,
    double** dX_array,    magma_minproduct_int_t dX_length,
    double** dinvA_array, magma_minproduct_int_t dinvA_length,
    double** dW1_displ, double** dW2_displ,  
    double** dW3_displ, double** dW4_displ,
    double** dW5_displ, 
    magma_minproduct_int_t *info_array, magma_minproduct_int_t gbstep,
    magma_minproduct_int_t batchCount, cublasHandle_t myhandle, magma_minproduct_queue_t queue);


magma_minproduct_int_t 
magma_minproduct_dgetrf_nopiv_batched(
        magma_minproduct_int_t m, magma_minproduct_int_t n,
        double **dA_array, 
        magma_minproduct_int_t lda,
        magma_minproduct_int_t *info_array, 
        magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_dgetrs_nopiv_batched(
                  magma_minproduct_trans_t trans, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
                  double **dA_array, magma_minproduct_int_t ldda,
                  double **dB_array, magma_minproduct_int_t lddb,
                  magma_minproduct_int_t *info_array,
                  magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_dgesv_nopiv_batched(
                  magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
                  double **dA_array, magma_minproduct_int_t ldda,
                  double **dB_array, magma_minproduct_int_t lddb,
                  magma_minproduct_int_t *info_array,
                  magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_dgesv_rbt_batched(
                  magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
                  double **dA_array, magma_minproduct_int_t ldda,
                  double **dB_array, magma_minproduct_int_t lddb,
                  magma_minproduct_int_t *info_array,
                  magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_dgesv_batched(
                  magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
                  double **dA_array, magma_minproduct_int_t ldda,
                  magma_minproduct_int_t **dipiv_array, 
                  double **dB_array, magma_minproduct_int_t lddb,
                  magma_minproduct_int_t *dinfo_array,
                  magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

magma_minproduct_int_t
magma_minproduct_dgerbt_batched(
    magma_minproduct_bool_t gen, magma_minproduct_int_t n, magma_minproduct_int_t nrhs,
    double **dA_array, magma_minproduct_int_t ldda,
    double **dB_array, magma_minproduct_int_t lddb,
    double *U, double *V,
    magma_minproduct_int_t *info, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

void 
magma_minproductblas_dprbt_batched(
    magma_minproduct_int_t n, 
    double **dA_array, magma_minproduct_int_t ldda, 
    double *du, double *dv,
    magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

void
magma_minproductblas_dprbt_mv_batched(
    magma_minproduct_int_t n, 
    double *dv, double **db_array, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);


void
magma_minproductblas_dprbt_mtv_batched(
    magma_minproduct_int_t n, 
    double *du, double **db_array, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);





void 
magma_minproduct_dlacgv_batched(magma_minproduct_int_t n, double **x_array, magma_minproduct_int_t incx, int offset, int batchCount, magma_minproduct_queue_t queue);

void 
magma_minproduct_dpotf2_dscal_batched(magma_minproduct_int_t n, double **x_array, magma_minproduct_int_t incx, magma_minproduct_int_t offset, magma_minproduct_int_t *info_array, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

void 
magma_minproduct_dpotf2_ddot_batched(magma_minproduct_int_t n, double **x_array, magma_minproduct_int_t incx, magma_minproduct_int_t offset, magma_minproduct_int_t *info_array, magma_minproduct_int_t gbstep, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);


void 
setup_pivinfo( magma_minproduct_int_t *pivinfo, magma_minproduct_int_t *ipiv, 
                      magma_minproduct_int_t m, magma_minproduct_int_t nb, 
                      magma_minproduct_queue_t queue);


void
magma_minproductblas_dgeadd_batched_q(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    double alpha,
    magma_minproductDouble_const_ptr  const dAarray[], magma_minproduct_int_t ldda,
    magma_minproductDouble_ptr              dBarray[], magma_minproduct_int_t lddb,
    magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue );

void
magma_minproductblas_dlacpy_batched(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDouble_const_ptr  const dAarray[], magma_minproduct_int_t ldda,
    magma_minproductDouble_ptr              dBarray[], magma_minproduct_int_t lddb,
    magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue );

void
magma_minproductblas_dgeadd_batched(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    double alpha,
    magma_minproductDouble_const_ptr  const dAarray[], magma_minproduct_int_t ldda,
    magma_minproductDouble_ptr              dBarray[], magma_minproduct_int_t lddb,
    magma_minproduct_int_t batchCount );


void
magma_minproductblas_dgemv_batched(
    magma_minproduct_trans_t trans, magma_minproduct_int_t m, magma_minproduct_int_t n,
    double alpha,
    magma_minproductDouble_ptr dA_array[], magma_minproduct_int_t ldda,
    magma_minproductDouble_ptr dx_array[], magma_minproduct_int_t incx,
    double beta,
    magma_minproductDouble_ptr dy_array[], magma_minproduct_int_t incy,
    magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);


magma_minproduct_int_t 
magma_minproduct_dgeqrf_batched(
    magma_minproduct_int_t m, magma_minproduct_int_t n, 
    double **dA_array,
    magma_minproduct_int_t lda, 
    double **tau_array,
    magma_minproduct_int_t *info_array, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

magma_minproduct_int_t 
magma_minproduct_dgeqrf_batched_v4(
    magma_minproduct_int_t m, magma_minproduct_int_t n, 
    double **dA_array,
    magma_minproduct_int_t lda, 
    double **tau_array,
    magma_minproduct_int_t *info_array, magma_minproduct_int_t batchCount);

magma_minproduct_int_t
magma_minproduct_dgeqrf_panel_batched(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t nb,    
    double** dA_array,    magma_minproduct_int_t ldda,
    double** tau_array, 
    double** dT_array, magma_minproduct_int_t ldt, 
    double** dR_array, magma_minproduct_int_t ldr,
    double** dW0_displ, 
    double** dW1_displ,
    double *dwork,  
    double** W_array, 
    double** W2_array,
    magma_minproduct_int_t *info_array,
    magma_minproduct_int_t batchCount, cublasHandle_t myhandle, magma_minproduct_queue_t queue);

magma_minproduct_int_t
magma_minproduct_dgeqrf_panel_batched_v4(
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t nb,    
    double** dA_array,    magma_minproduct_int_t ldda,
    double** tau_array, 
    double** dT_array, magma_minproduct_int_t ldt, 
    double** dR_array, magma_minproduct_int_t ldr,
    double** dnorm_array,  
    double** dW0_displ, 
    double** dW1_displ,
    double *dwork,  
    double** W_array, 
    double** W2_array,
    magma_minproduct_int_t *info_array,
    magma_minproduct_int_t batchCount, cublasHandle_t myhandle);

magma_minproduct_int_t
magma_minproduct_dgeqr2x_batched_v4(magma_minproduct_int_t m, magma_minproduct_int_t n, double **dA_array,
                  magma_minproduct_int_t lda, 
                  double **tau_array,
                  double **dR_array, magma_minproduct_int_t ldr,
                  double **dwork_array,  
                  magma_minproduct_int_t *info, magma_minproduct_int_t batchCount);

magma_minproduct_int_t
magma_minproduct_dgeqr2_batched(magma_minproduct_int_t m, magma_minproduct_int_t n, double **dA_array,
                  magma_minproduct_int_t lda, 
                  double **tau_array,
                  magma_minproduct_int_t *info, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

magma_minproduct_int_t
magma_minproduct_dlarfb_dgemm_batched(
                  cublasHandle_t myhandle,
                  magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
                  double **dV_array,    magma_minproduct_int_t ldv,
                  double **dT_array,    magma_minproduct_int_t ldt,
                  double **dA_array,    magma_minproduct_int_t lda,
                  double **W_array,     magma_minproduct_int_t ldw,
                  double **W2_array,    magma_minproduct_int_t ldw2,
                  magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

magma_minproduct_int_t
magma_minproduct_dlarfb_gemm_batched(
    magma_minproduct_side_t side, magma_minproduct_trans_t trans, magma_minproduct_direct_t direct, magma_minproduct_storev_t storev,
    magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k,
    magma_minproductDouble_const_ptr dV_array[],    magma_minproduct_int_t lddv,
    magma_minproductDouble_const_ptr dT_array[],    magma_minproduct_int_t lddt,
    magma_minproductDouble_ptr dC_array[],          magma_minproduct_int_t lddc,
    magma_minproductDouble_ptr dwork_array[],       magma_minproduct_int_t ldwork,
    magma_minproductDouble_ptr dworkvt_array[],     magma_minproduct_int_t ldworkvt,
    magma_minproduct_int_t batchCount, cublasHandle_t myhandle, magma_minproduct_queue_t queue);




void
magma_minproduct_dlarft_batched_vold(magma_minproduct_int_t n, magma_minproduct_int_t k, double **v_array, magma_minproduct_int_t ldv,
                    double **tau_array,
                    double **T_array, magma_minproduct_int_t ldt, 
                    magma_minproduct_int_t batchCount);





magma_minproduct_int_t
magma_minproduct_dlarft_batched(magma_minproduct_int_t n, magma_minproduct_int_t k, magma_minproduct_int_t stair_T, 
                double **v_array, magma_minproduct_int_t ldv,
                double **tau_array, double **T_array, magma_minproduct_int_t ldt, 
                double **work_array, magma_minproduct_int_t lwork, magma_minproduct_int_t batchCount, cublasHandle_t myhandle, magma_minproduct_queue_t queue);

void
magma_minproduct_dlarft_sm32x32_batched(magma_minproduct_int_t n, magma_minproduct_int_t k, double **v_array, magma_minproduct_int_t ldv,
                    double **tau_array, double **T_array, magma_minproduct_int_t ldt, 
                    magma_minproduct_int_t batchCount, cublasHandle_t myhandle, magma_minproduct_queue_t queue);



void magma_minproductblas_dlarft_recdtrmv_sm32x32(
    magma_minproduct_int_t m, magma_minproduct_int_t n, 
    double *tau, 
    double *Trec, magma_minproduct_int_t ldtrec, 
    double *Ttri, magma_minproduct_int_t ldttri);


void magma_minproductblas_dlarft_recdtrmv_sm32x32_batched(
    magma_minproduct_int_t m, magma_minproduct_int_t n, 
    double **tau_array, 
    double **Trec_array, magma_minproduct_int_t ldtrec, 
    double **Ttri_array, magma_minproduct_int_t ldttri,
    magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

void magma_minproductblas_dlarft_dtrmv_sm32x32(
    magma_minproduct_int_t m, magma_minproduct_int_t n, 
    double *tau, 
    double *Tin, magma_minproduct_int_t ldtin, 
    double *Tout, magma_minproduct_int_t ldtout);

void magma_minproductblas_dlarft_dtrmv_sm32x32_batched(
    magma_minproduct_int_t m, magma_minproduct_int_t n, 
    double **tau_array, 
    double **Tin_array, magma_minproduct_int_t ldtin, 
    double **Tout_array, magma_minproduct_int_t ldtout,
    magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

void magma_minproductblas_dlarft_gemv_loop_inside(
    int n, int k, 
    double *tau, 
    double *v, int ldv, 
    double *T, int ldt);

void magma_minproductblas_dlarft_gemv_loop_inside_batched(
    int n, int k, 
    double **tau_array, 
    double **v_array, int ldv, 
    double **T_array, int ldt, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

void magma_minproductblas_dlarft_gemvrowwise(
    magma_minproduct_int_t m, magma_minproduct_int_t i, 
    double *tau, 
    double *v, magma_minproduct_int_t ldv, 
    double *T, magma_minproduct_int_t ldt,
    double *W);

void magma_minproductblas_dlarft_gemvrowwise_batched(
    magma_minproduct_int_t m, magma_minproduct_int_t i, 
    double **tau_array, 
    double **v_array, magma_minproduct_int_t ldv, 
    double **T_array, magma_minproduct_int_t ldt,
    magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);

void magma_minproductblas_dlarft_gemvcolwise(
    magma_minproduct_int_t m,  magma_minproduct_int_t step,
    double *v, magma_minproduct_int_t ldv, 
    double *T,  magma_minproduct_int_t ldt,
    double *tau);

void magma_minproductblas_dlarft_gemvcolwise_batched(
    magma_minproduct_int_t m,  magma_minproduct_int_t step,
    double **v_array, magma_minproduct_int_t ldv, 
    double **T_array,  magma_minproduct_int_t ldt,
    double **tau_array, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);




void dgeqrf_copy_upper_batched(                
                  magma_minproduct_int_t n, magma_minproduct_int_t nb,
                  double **dV_array,    magma_minproduct_int_t ldv,
                  double **dR_array,    magma_minproduct_int_t ldr,
          magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue);



void 
magma_minproductblas_dnrm2_cols_batched(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    double **dA_array, magma_minproduct_int_t lda, 
    double **dxnorm_array, magma_minproduct_int_t batchCount);
 
void 
magma_minproduct_dlarfgx_batched(magma_minproduct_int_t n, double **dx0_array, double **dx_array, 
                  double **dtau_array, double **dxnorm_array, 
                  double **dR_array, magma_minproduct_int_t it, magma_minproduct_int_t batchCount);


void 
magma_minproduct_dlarfx_batched_v4(magma_minproduct_int_t m, magma_minproduct_int_t n, double **v_array, double **tau_array,
                double **C_array, magma_minproduct_int_t ldc, double **xnorm_array, 
                magma_minproduct_int_t step, 
                magma_minproduct_int_t batchCount);


void 
magma_minproductblas_dlarfg_batched(
    magma_minproduct_int_t n,
    double** dalpha_array, double** dx_array, magma_minproduct_int_t incx,
    double** dtau_array, magma_minproduct_int_t batchCount );





// for debugging purpose
void 
dset_stepinit_ipiv(magma_minproduct_int_t **ipiv_array,
                 magma_minproduct_int_t pm,
                 magma_minproduct_int_t batchCount);



#ifdef __cplusplus
}
#endif

#undef REAL

#endif  /* MAGMA_minproduct_DBATCHED_H */
