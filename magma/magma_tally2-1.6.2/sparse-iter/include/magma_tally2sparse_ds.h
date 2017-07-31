/*
-- MAGMA_tally2 (version 1.1) --
 Univ. of Tennessee, Knoxville
 Univ. of California, Berkeley
 Univ. of Colorado, Denver
 @date

 @generated from magma_tally2sparse_zc.h mixed zc -> ds, Fri May  1 21:08:05 2015
 @author Hartwig Anzt
*/

#ifndef MAGMA_tally2SPARSE_DS_H
#define MAGMA_tally2SPARSE_DS_H

#include "magma_tally2_types.h"
#include "magma_tally2sparse_types.h"

#define PRECISION_d


/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally2_SPARSE Matrix Descriptors
*/


#ifdef __cplusplus
extern "C" {
#endif

/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally2_SPARSE Auxiliary functions
*/
magma_tally2_int_t
magma_tally2_vector_dlag2s(
    magma_tally2_d_matrix x,
    magma_tally2_s_matrix *y,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_sparse_matrix_dlag2s(
    magma_tally2_d_matrix A,
    magma_tally2_s_matrix *B,
    magma_tally2_queue_t queue );


magma_tally2_int_t
magma_tally2_vector_slag2d(
    magma_tally2_s_matrix x,
    magma_tally2_d_matrix *y,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_sparse_matrix_slag2d(
    magma_tally2_s_matrix A,
    magma_tally2_d_matrix *B,
    magma_tally2_queue_t queue );

void
magma_tally2blas_dlag2s_sparse(
    magma_tally2_int_t M, 
    magma_tally2_int_t N , 
    magma_tally2Double_const_ptr dA, 
    magma_tally2_int_t lda, 
    magma_tally2Float_ptr dSA, 
    magma_tally2_int_t ldsa,
    magma_tally2_int_t *info,
    magma_tally2_queue_t queue );

void 
magma_tally2blas_slag2d_sparse(
    magma_tally2_int_t M, 
    magma_tally2_int_t N , 
    magma_tally2Float_const_ptr dSA, 
    magma_tally2_int_t ldsa, 
    magma_tally2Double_ptr dA, 
    magma_tally2_int_t lda,
    magma_tally2_int_t *info,
    magma_tally2_queue_t queue );

void 
magma_tally2_dlag2s_CSR_DENSE(
    magma_tally2_d_matrix A,
    magma_tally2_s_matrix *B,
    magma_tally2_queue_t queue );

void 
magma_tally2_dlag2s_CSR_DENSE_alloc(
    magma_tally2_d_matrix A,
    magma_tally2_s_matrix *B,
    magma_tally2_queue_t queue );

void 
magma_tally2_dlag2s_CSR_DENSE_convert(
    magma_tally2_d_matrix A,
    magma_tally2_s_matrix *B,
    magma_tally2_queue_t queue );

/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally2_SPARSE function definitions / Data on CPU
*/


/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally2_SPARSE function definitions / Data on CPU / Multi-GPU
*/

/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally2_SPARSE function definitions / Data on GPU
*/
magma_tally2_int_t
magma_tally2_dspgmres(
    magma_tally2_d_matrix A, 
    magma_tally2_d_matrix b, 
    magma_tally2_d_matrix *x,
    magma_tally2_d_solver_par *solver_par,
    magma_tally2_d_preconditioner *precond_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dspbicgstab(
    magma_tally2_d_matrix A, 
    magma_tally2_d_matrix b, 
    magma_tally2_d_matrix *x,
    magma_tally2_d_solver_par *solver_par,
    magma_tally2_d_preconditioner *precond_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dsir(
    magma_tally2_d_matrix A, 
    magma_tally2_d_matrix b, 
    magma_tally2_d_matrix *x,
    magma_tally2_d_solver_par *solver_par, 
    magma_tally2_d_preconditioner *precond_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dspir(
    magma_tally2_d_matrix A, 
    magma_tally2_d_matrix b, 
    magma_tally2_d_matrix *x,
    magma_tally2_d_solver_par *solver_par, 
    magma_tally2_d_preconditioner *precond_par,
    magma_tally2_queue_t queue );


/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally2_SPARSE utility function definitions
*/



/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally2_SPARSE BLAS function definitions
*/



#ifdef __cplusplus
}
#endif

#undef PRECISION_d
#endif /* MAGMA_tally2SPARSE_DS_H */
