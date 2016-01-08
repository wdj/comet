/* 
    -- MAGMA_tally4 (version 1.6.1) -- 
       Univ. of Tennessee, Knoxville 
       Univ. of California, Berkeley 
       Univ. of Colorado, Denver 
       May 2013 
 
       @author: Simplice Donfack 
*/
#ifndef MAGMA_tally4_INSERT_D
#define MAGMA_tally4_INSERT_D

#include "common_magma_tally4.h"


 
 /*GPU task wrapper*/
void magma_tally4_insert_dmalloc_pinned(magma_tally4_int_t size, double **A, void *A_dep_ptr);
void magma_tally4_insert_dfree_pinned(double *A, void *A_dep_ptr);
void magma_tally4_insert_dfree_pinned_index(double **A, int index, void *A_dep_ptr);

void magma_tally4_insert_dsetmatrix(magma_tally4_int_t m, magma_tally4_int_t nb, double *A_src, magma_tally4_int_t LDA, double *dA_dst, magma_tally4_int_t dA_LD);

void magma_tally4_insert_dgetmatrix(magma_tally4_int_t m, magma_tally4_int_t nb, double *dA_src, magma_tally4_int_t dA_LD, double *A_dst, magma_tally4_int_t LDA);

void magma_tally4_insert_dsetmatrix_transpose(magma_tally4_int_t m, magma_tally4_int_t nb, double *A_src, magma_tally4_int_t LDA, double *dA_dst, magma_tally4_int_t dA_LD, double *dwork, magma_tally4_int_t dwork_LD, void *A_src_dep_ptr, void *dA_dst_dep_ptr);

void magma_tally4_insert_dgetmatrix_transpose(magma_tally4_int_t m, magma_tally4_int_t nb, double *dA_src, magma_tally4_int_t dA_LD, double *A_dst, magma_tally4_int_t LDA, double *dwork, magma_tally4_int_t dwork_LD, void *A_dst_dep_ptr);

void magma_tally4_insert_dlaswp( magma_tally4_int_t n, double *dA, magma_tally4_int_t lda, magma_tally4_int_t i1, magma_tally4_int_t i2, magma_tally4_int_t *ipiv, magma_tally4_int_t inci, void *dA_dep_ptr);
void magma_tally4_insert_dtrsm(magma_tally4_side_t side, magma_tally4_uplo_t uplo, magma_tally4_trans_t trans, magma_tally4_diag_t diag, 
                                 magma_tally4_int_t m, magma_tally4_int_t n, double alpha, 
                                 double *dA, magma_tally4_int_t lda, double *dB, magma_tally4_int_t ldb );
void magma_tally4_insert_dgemm(magma_tally4_trans_t transA, magma_tally4_trans_t transB, 
                                 magma_tally4_int_t m, magma_tally4_int_t n, magma_tally4_int_t k, double alpha, 
                                 double *dA, magma_tally4_int_t lda, double *dB, magma_tally4_int_t ldb, double beta, double *dC, magma_tally4_int_t ldc );
#endif

