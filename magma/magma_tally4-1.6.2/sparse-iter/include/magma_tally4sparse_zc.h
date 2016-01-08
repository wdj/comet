/*
-- MAGMA_tally4 (version 1.1) --
 Univ. of Tennessee, Knoxville
 Univ. of California, Berkeley
 Univ. of Colorado, Denver
 @date

 @precisions mixed zc -> ds
 @author Hartwig Anzt
*/

#ifndef MAGMA_tally4SPARSE_ZC_H
#define MAGMA_tally4SPARSE_ZC_H

#include "magma_tally4_types.h"
#include "magma_tally4sparse_types.h"

#define PRECISION_z


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
magma_tally4_vector_zlag2c(
    magma_tally4_z_matrix x,
    magma_tally4_c_matrix *y,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_sparse_matrix_zlag2c(
    magma_tally4_z_matrix A,
    magma_tally4_c_matrix *B,
    magma_tally4_queue_t queue );


magma_tally4_int_t
magma_tally4_vector_clag2z(
    magma_tally4_c_matrix x,
    magma_tally4_z_matrix *y,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_sparse_matrix_clag2z(
    magma_tally4_c_matrix A,
    magma_tally4_z_matrix *B,
    magma_tally4_queue_t queue );

void
magma_tally4blas_zlag2c_sparse(
    magma_tally4_int_t M, 
    magma_tally4_int_t N , 
    magma_tally4DoubleComplex_const_ptr dA, 
    magma_tally4_int_t lda, 
    magma_tally4FloatComplex_ptr dSA, 
    magma_tally4_int_t ldsa,
    magma_tally4_int_t *info,
    magma_tally4_queue_t queue );

void 
magma_tally4blas_clag2z_sparse(
    magma_tally4_int_t M, 
    magma_tally4_int_t N , 
    magma_tally4FloatComplex_const_ptr dSA, 
    magma_tally4_int_t ldsa, 
    magma_tally4DoubleComplex_ptr dA, 
    magma_tally4_int_t lda,
    magma_tally4_int_t *info,
    magma_tally4_queue_t queue );

void 
magma_tally4_zlag2c_CSR_DENSE(
    magma_tally4_z_matrix A,
    magma_tally4_c_matrix *B,
    magma_tally4_queue_t queue );

void 
magma_tally4_zlag2c_CSR_DENSE_alloc(
    magma_tally4_z_matrix A,
    magma_tally4_c_matrix *B,
    magma_tally4_queue_t queue );

void 
magma_tally4_zlag2c_CSR_DENSE_convert(
    magma_tally4_z_matrix A,
    magma_tally4_c_matrix *B,
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
magma_tally4_zcpgmres(
    magma_tally4_z_matrix A, 
    magma_tally4_z_matrix b, 
    magma_tally4_z_matrix *x,
    magma_tally4_z_solver_par *solver_par,
    magma_tally4_z_preconditioner *precond_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zcpbicgstab(
    magma_tally4_z_matrix A, 
    magma_tally4_z_matrix b, 
    magma_tally4_z_matrix *x,
    magma_tally4_z_solver_par *solver_par,
    magma_tally4_z_preconditioner *precond_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zcir(
    magma_tally4_z_matrix A, 
    magma_tally4_z_matrix b, 
    magma_tally4_z_matrix *x,
    magma_tally4_z_solver_par *solver_par, 
    magma_tally4_z_preconditioner *precond_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zcpir(
    magma_tally4_z_matrix A, 
    magma_tally4_z_matrix b, 
    magma_tally4_z_matrix *x,
    magma_tally4_z_solver_par *solver_par, 
    magma_tally4_z_preconditioner *precond_par,
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

#undef PRECISION_z
#endif /* MAGMA_tally4SPARSE_ZC_H */
