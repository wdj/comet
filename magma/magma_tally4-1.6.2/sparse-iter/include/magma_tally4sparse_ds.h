/*
-- MAGMA_tally4 (version 1.1) --
 Univ. of Tennessee, Knoxville
 Univ. of California, Berkeley
 Univ. of Colorado, Denver
 @date

 @generated from magma_tally4sparse_zc.h mixed zc -> ds, Fri May  1 21:08:05 2015
 @author Hartwig Anzt
*/

#ifndef MAGMA_tally4SPARSE_DS_H
#define MAGMA_tally4SPARSE_DS_H

#include "magma_tally4_types.h"
#include "magma_tally4sparse_types.h"

#define PRECISION_d


/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally4_SPARSE Matrix Descriptors
*/


#ifdef __cplusplus
extern "C" {
#endif

/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally4_SPARSE Auxiliary functions
*/
magma_tally4_int_t
magma_tally4_vector_dlag2s(
    magma_tally4_d_matrix x,
    magma_tally4_s_matrix *y,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_sparse_matrix_dlag2s(
    magma_tally4_d_matrix A,
    magma_tally4_s_matrix *B,
    magma_tally4_queue_t queue );


magma_tally4_int_t
magma_tally4_vector_slag2d(
    magma_tally4_s_matrix x,
    magma_tally4_d_matrix *y,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_sparse_matrix_slag2d(
    magma_tally4_s_matrix A,
    magma_tally4_d_matrix *B,
    magma_tally4_queue_t queue );

void
magma_tally4blas_dlag2s_sparse(
    magma_tally4_int_t M, 
    magma_tally4_int_t N , 
    magma_tally4Double_const_ptr dA, 
    magma_tally4_int_t lda, 
    magma_tally4Float_ptr dSA, 
    magma_tally4_int_t ldsa,
    magma_tally4_int_t *info,
    magma_tally4_queue_t queue );

void 
magma_tally4blas_slag2d_sparse(
    magma_tally4_int_t M, 
    magma_tally4_int_t N , 
    magma_tally4Float_const_ptr dSA, 
    magma_tally4_int_t ldsa, 
    magma_tally4Double_ptr dA, 
    magma_tally4_int_t lda,
    magma_tally4_int_t *info,
    magma_tally4_queue_t queue );

void 
magma_tally4_dlag2s_CSR_DENSE(
    magma_tally4_d_matrix A,
    magma_tally4_s_matrix *B,
    magma_tally4_queue_t queue );

void 
magma_tally4_dlag2s_CSR_DENSE_alloc(
    magma_tally4_d_matrix A,
    magma_tally4_s_matrix *B,
    magma_tally4_queue_t queue );

void 
magma_tally4_dlag2s_CSR_DENSE_convert(
    magma_tally4_d_matrix A,
    magma_tally4_s_matrix *B,
    magma_tally4_queue_t queue );

/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally4_SPARSE function definitions / Data on CPU
*/


/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally4_SPARSE function definitions / Data on CPU / Multi-GPU
*/

/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally4_SPARSE function definitions / Data on GPU
*/
magma_tally4_int_t
magma_tally4_dspgmres(
    magma_tally4_d_matrix A, 
    magma_tally4_d_matrix b, 
    magma_tally4_d_matrix *x,
    magma_tally4_d_solver_par *solver_par,
    magma_tally4_d_preconditioner *precond_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dspbicgstab(
    magma_tally4_d_matrix A, 
    magma_tally4_d_matrix b, 
    magma_tally4_d_matrix *x,
    magma_tally4_d_solver_par *solver_par,
    magma_tally4_d_preconditioner *precond_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dsir(
    magma_tally4_d_matrix A, 
    magma_tally4_d_matrix b, 
    magma_tally4_d_matrix *x,
    magma_tally4_d_solver_par *solver_par, 
    magma_tally4_d_preconditioner *precond_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dspir(
    magma_tally4_d_matrix A, 
    magma_tally4_d_matrix b, 
    magma_tally4_d_matrix *x,
    magma_tally4_d_solver_par *solver_par, 
    magma_tally4_d_preconditioner *precond_par,
    magma_tally4_queue_t queue );


/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally4_SPARSE utility function definitions
*/



/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally4_SPARSE BLAS function definitions
*/



#ifdef __cplusplus
}
#endif

#undef PRECISION_d
#endif /* MAGMA_tally4SPARSE_DS_H */
