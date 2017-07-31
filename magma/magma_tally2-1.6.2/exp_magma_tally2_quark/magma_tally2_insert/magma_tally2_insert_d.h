/* 
    -- MAGMA_tally2 (version 1.6.1) -- 
       Univ. of Tennessee, Knoxville 
       Univ. of California, Berkeley 
       Univ. of Colorado, Denver 
       May 2013 
 
       @author: Simplice Donfack 
*/
#ifndef MAGMA_tally2_INSERT_D
#define MAGMA_tally2_INSERT_D

#include "common_magma_tally2.h"


 
 /*GPU task wrapper*/
void magma_tally2_insert_dmalloc_pinned(magma_tally2_int_t size, double **A, void *A_dep_ptr);
void magma_tally2_insert_dfree_pinned(double *A, void *A_dep_ptr);
void magma_tally2_insert_dfree_pinned_index(double **A, int index, void *A_dep_ptr);

void magma_tally2_insert_dsetmatrix(magma_tally2_int_t m, magma_tally2_int_t nb, double *A_src, magma_tally2_int_t LDA, double *dA_dst, magma_tally2_int_t dA_LD);

void magma_tally2_insert_dgetmatrix(magma_tally2_int_t m, magma_tally2_int_t nb, double *dA_src, magma_tally2_int_t dA_LD, double *A_dst, magma_tally2_int_t LDA);

void magma_tally2_insert_dsetmatrix_transpose(magma_tally2_int_t m, magma_tally2_int_t nb, double *A_src, magma_tally2_int_t LDA, double *dA_dst, magma_tally2_int_t dA_LD, double *dwork, magma_tally2_int_t dwork_LD, void *A_src_dep_ptr, void *dA_dst_dep_ptr);

void magma_tally2_insert_dgetmatrix_transpose(magma_tally2_int_t m, magma_tally2_int_t nb, double *dA_src, magma_tally2_int_t dA_LD, double *A_dst, magma_tally2_int_t LDA, double *dwork, magma_tally2_int_t dwork_LD, void *A_dst_dep_ptr);

void magma_tally2_insert_dlaswp( magma_tally2_int_t n, double *dA, magma_tally2_int_t lda, magma_tally2_int_t i1, magma_tally2_int_t i2, magma_tally2_int_t *ipiv, magma_tally2_int_t inci, void *dA_dep_ptr);
void magma_tally2_insert_dtrsm(magma_tally2_side_t side, magma_tally2_uplo_t uplo, magma_tally2_trans_t trans, magma_tally2_diag_t diag, 
                                 magma_tally2_int_t m, magma_tally2_int_t n, double alpha, 
                                 double *dA, magma_tally2_int_t lda, double *dB, magma_tally2_int_t ldb );
void magma_tally2_insert_dgemm(magma_tally2_trans_t transA, magma_tally2_trans_t transB, 
                                 magma_tally2_int_t m, magma_tally2_int_t n, magma_tally2_int_t k, double alpha, 
                                 double *dA, magma_tally2_int_t lda, double *dB, magma_tally2_int_t ldb, double beta, double *dC, magma_tally2_int_t ldc );
#endif

