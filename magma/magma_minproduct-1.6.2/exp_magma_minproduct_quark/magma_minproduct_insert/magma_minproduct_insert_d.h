/* 
    -- MAGMA_minproduct (version 1.6.1) -- 
       Univ. of Tennessee, Knoxville 
       Univ. of California, Berkeley 
       Univ. of Colorado, Denver 
       May 2013 
 
       @author: Simplice Donfack 
*/
#ifndef MAGMA_minproduct_INSERT_D
#define MAGMA_minproduct_INSERT_D

#include "common_magma_minproduct.h"


 
 /*GPU task wrapper*/
void magma_minproduct_insert_dmalloc_pinned(magma_minproduct_int_t size, double **A, void *A_dep_ptr);
void magma_minproduct_insert_dfree_pinned(double *A, void *A_dep_ptr);
void magma_minproduct_insert_dfree_pinned_index(double **A, int index, void *A_dep_ptr);

void magma_minproduct_insert_dsetmatrix(magma_minproduct_int_t m, magma_minproduct_int_t nb, double *A_src, magma_minproduct_int_t LDA, double *dA_dst, magma_minproduct_int_t dA_LD);

void magma_minproduct_insert_dgetmatrix(magma_minproduct_int_t m, magma_minproduct_int_t nb, double *dA_src, magma_minproduct_int_t dA_LD, double *A_dst, magma_minproduct_int_t LDA);

void magma_minproduct_insert_dsetmatrix_transpose(magma_minproduct_int_t m, magma_minproduct_int_t nb, double *A_src, magma_minproduct_int_t LDA, double *dA_dst, magma_minproduct_int_t dA_LD, double *dwork, magma_minproduct_int_t dwork_LD, void *A_src_dep_ptr, void *dA_dst_dep_ptr);

void magma_minproduct_insert_dgetmatrix_transpose(magma_minproduct_int_t m, magma_minproduct_int_t nb, double *dA_src, magma_minproduct_int_t dA_LD, double *A_dst, magma_minproduct_int_t LDA, double *dwork, magma_minproduct_int_t dwork_LD, void *A_dst_dep_ptr);

void magma_minproduct_insert_dlaswp( magma_minproduct_int_t n, double *dA, magma_minproduct_int_t lda, magma_minproduct_int_t i1, magma_minproduct_int_t i2, magma_minproduct_int_t *ipiv, magma_minproduct_int_t inci, void *dA_dep_ptr);
void magma_minproduct_insert_dtrsm(magma_minproduct_side_t side, magma_minproduct_uplo_t uplo, magma_minproduct_trans_t trans, magma_minproduct_diag_t diag, 
                                 magma_minproduct_int_t m, magma_minproduct_int_t n, double alpha, 
                                 double *dA, magma_minproduct_int_t lda, double *dB, magma_minproduct_int_t ldb );
void magma_minproduct_insert_dgemm(magma_minproduct_trans_t transA, magma_minproduct_trans_t transB, 
                                 magma_minproduct_int_t m, magma_minproduct_int_t n, magma_minproduct_int_t k, double alpha, 
                                 double *dA, magma_minproduct_int_t lda, double *dB, magma_minproduct_int_t ldb, double beta, double *dC, magma_minproduct_int_t ldc );
#endif

