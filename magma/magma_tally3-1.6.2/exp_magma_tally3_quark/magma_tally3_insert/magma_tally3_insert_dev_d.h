/* 
    -- MAGMA_tally3 (version 1.6.1) -- 
       Univ. of Tennessee, Knoxville 
       Univ. of California, Berkeley 
       Univ. of Colorado, Denver 
       May 2013 
 
       @author: Simplice Donfack 
*/
#ifndef MAGMA_tally3_INSERT_DEV_D_H
#define MAGMA_tally3_INSERT_DEV_D_H



void magma_tally3_insert_dev_dmalloc_pinned(magma_tally3_int_t deviceID, magma_tally3_int_t size, double **A, void *A_dep_ptr);
void magma_tally3_insert_dev_dfree_pinned(magma_tally3_int_t deviceID, double *A, void *A_dep_ptr);
void magma_tally3_insert_dev_dfree_pinned_index(magma_tally3_int_t deviceID, double **A, int index, void *A_dep_ptr);


void magma_tally3_insert_dev_queue_sync(magma_tally3_int_t deviceID, magma_tally3_queue_t stream1, void *dep_ptr);

void magma_tally3_insert_dev_dsetmatrix(magma_tally3_int_t deviceID, magma_tally3_int_t m, magma_tally3_int_t nb, double *A_src, magma_tally3_int_t LDA, double *dA_dst, magma_tally3_int_t dA_LD);

void magma_tally3_insert_dev_dgetmatrix(magma_tally3_int_t deviceID, magma_tally3_int_t m, magma_tally3_int_t nb, double *dA_src, magma_tally3_int_t dA_LD, double *A_dst, magma_tally3_int_t LDA);

void magma_tally3_insert_dev_dsetmatrix_transpose(magma_tally3_int_t deviceID, magma_tally3_int_t m, magma_tally3_int_t nb, double *A_src, magma_tally3_int_t LDA, double *dA_dst, magma_tally3_int_t dA_LD, double *dwork, magma_tally3_int_t dwork_LD, void *A_src_dep_ptr, void *dA_dst_dep_ptr);
void magma_tally3_insert_dev_dsetmatrix_async_transpose(magma_tally3_int_t deviceID, magma_tally3_int_t m, magma_tally3_int_t nb, double *A_src, magma_tally3_int_t LDA, double *dA_dst, magma_tally3_int_t dA_LD, magma_tally3_queue_t stream1, double *dwork, magma_tally3_int_t dwork_LD, void *A_src_dep_ptr, void *dA_dst_dep_ptr);

void magma_tally3_insert_dev_dgetmatrix_transpose(magma_tally3_int_t deviceID, magma_tally3_int_t m, magma_tally3_int_t nb, double *dA_src, magma_tally3_int_t dA_LD, double *A_dst, magma_tally3_int_t LDA, double *dwork, magma_tally3_int_t dwork_LD, void *A_dst_dep_ptr);
void magma_tally3_insert_dev_dgetmatrix_async_transpose(magma_tally3_int_t deviceID, magma_tally3_int_t m, magma_tally3_int_t nb, double *dA_src, magma_tally3_int_t dA_LD, double *A_dst, magma_tally3_int_t LDA, magma_tally3_queue_t stream1, double *dwork, magma_tally3_int_t dwork_LD, void *A_dst_dep_ptr);

void magma_tally3_insert_dev_dlaswp(magma_tally3_int_t deviceID,  magma_tally3_int_t n, double *dA, magma_tally3_int_t lda, magma_tally3_int_t i1, magma_tally3_int_t i2, magma_tally3_int_t *ipiv, magma_tally3_int_t inci, void *dA_dep_ptr);
void magma_tally3_insert_dev_dtrsm(magma_tally3_int_t deviceID, magma_tally3_side_t side, magma_tally3_uplo_t uplo, magma_tally3_trans_t trans, magma_tally3_diag_t diag, 
                                 magma_tally3_int_t m, magma_tally3_int_t n, double alpha, 
                                 double *dA, magma_tally3_int_t lda, double *dB, magma_tally3_int_t ldb );
void magma_tally3_insert_dev_dgemm(magma_tally3_int_t deviceID, magma_tally3_trans_t transA, magma_tally3_trans_t transB, 
                                 magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t k, double alpha, 
                                 double *dA, magma_tally3_int_t lda, double *dB, magma_tally3_int_t ldb, double beta, double *dC, magma_tally3_int_t ldc );
#endif
