/*
-- MAGMA_tally3 (version 1.1) --
 Univ. of Tennessee, Knoxville
 Univ. of California, Berkeley
 Univ. of Colorado, Denver
 @date

 @generated from magma_tally3sparse_zc.h mixed zc -> ds, Fri May  1 21:08:05 2015
 @author Hartwig Anzt
*/

#ifndef MAGMA_tally3SPARSE_DS_H
#define MAGMA_tally3SPARSE_DS_H

#include "magma_tally3_types.h"
#include "magma_tally3sparse_types.h"

#define PRECISION_d


/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally3_SPARSE Matrix Descriptors
*/


#ifdef __cplusplus
extern "C" {
#endif

/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally3_SPARSE Auxiliary functions
*/
magma_tally3_int_t
magma_tally3_vector_dlag2s(
    magma_tally3_d_matrix x,
    magma_tally3_s_matrix *y,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_sparse_matrix_dlag2s(
    magma_tally3_d_matrix A,
    magma_tally3_s_matrix *B,
    magma_tally3_queue_t queue );


magma_tally3_int_t
magma_tally3_vector_slag2d(
    magma_tally3_s_matrix x,
    magma_tally3_d_matrix *y,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_sparse_matrix_slag2d(
    magma_tally3_s_matrix A,
    magma_tally3_d_matrix *B,
    magma_tally3_queue_t queue );

void
magma_tally3blas_dlag2s_sparse(
    magma_tally3_int_t M, 
    magma_tally3_int_t N , 
    magma_tally3Double_const_ptr dA, 
    magma_tally3_int_t lda, 
    magma_tally3Float_ptr dSA, 
    magma_tally3_int_t ldsa,
    magma_tally3_int_t *info,
    magma_tally3_queue_t queue );

void 
magma_tally3blas_slag2d_sparse(
    magma_tally3_int_t M, 
    magma_tally3_int_t N , 
    magma_tally3Float_const_ptr dSA, 
    magma_tally3_int_t ldsa, 
    magma_tally3Double_ptr dA, 
    magma_tally3_int_t lda,
    magma_tally3_int_t *info,
    magma_tally3_queue_t queue );

void 
magma_tally3_dlag2s_CSR_DENSE(
    magma_tally3_d_matrix A,
    magma_tally3_s_matrix *B,
    magma_tally3_queue_t queue );

void 
magma_tally3_dlag2s_CSR_DENSE_alloc(
    magma_tally3_d_matrix A,
    magma_tally3_s_matrix *B,
    magma_tally3_queue_t queue );

void 
magma_tally3_dlag2s_CSR_DENSE_convert(
    magma_tally3_d_matrix A,
    magma_tally3_s_matrix *B,
    magma_tally3_queue_t queue );

/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally3_SPARSE function definitions / Data on CPU
*/


/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally3_SPARSE function definitions / Data on CPU / Multi-GPU
*/

/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally3_SPARSE function definitions / Data on GPU
*/
magma_tally3_int_t
magma_tally3_dspgmres(
    magma_tally3_d_matrix A, 
    magma_tally3_d_matrix b, 
    magma_tally3_d_matrix *x,
    magma_tally3_d_solver_par *solver_par,
    magma_tally3_d_preconditioner *precond_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dspbicgstab(
    magma_tally3_d_matrix A, 
    magma_tally3_d_matrix b, 
    magma_tally3_d_matrix *x,
    magma_tally3_d_solver_par *solver_par,
    magma_tally3_d_preconditioner *precond_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dsir(
    magma_tally3_d_matrix A, 
    magma_tally3_d_matrix b, 
    magma_tally3_d_matrix *x,
    magma_tally3_d_solver_par *solver_par, 
    magma_tally3_d_preconditioner *precond_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dspir(
    magma_tally3_d_matrix A, 
    magma_tally3_d_matrix b, 
    magma_tally3_d_matrix *x,
    magma_tally3_d_solver_par *solver_par, 
    magma_tally3_d_preconditioner *precond_par,
    magma_tally3_queue_t queue );


/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally3_SPARSE utility function definitions
*/



/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally3_SPARSE BLAS function definitions
*/



#ifdef __cplusplus
}
#endif

#undef PRECISION_d
#endif /* MAGMA_tally3SPARSE_DS_H */
