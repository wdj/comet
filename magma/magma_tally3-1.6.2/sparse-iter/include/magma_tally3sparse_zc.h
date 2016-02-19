/*
-- MAGMA_tally3 (version 1.1) --
 Univ. of Tennessee, Knoxville
 Univ. of California, Berkeley
 Univ. of Colorado, Denver
 @date

 @precisions mixed zc -> ds
 @author Hartwig Anzt
*/

#ifndef MAGMA_tally3SPARSE_ZC_H
#define MAGMA_tally3SPARSE_ZC_H

#include "magma_tally3_types.h"
#include "magma_tally3sparse_types.h"

#define PRECISION_z


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
magma_tally3_vector_zlag2c(
    magma_tally3_z_matrix x,
    magma_tally3_c_matrix *y,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_sparse_matrix_zlag2c(
    magma_tally3_z_matrix A,
    magma_tally3_c_matrix *B,
    magma_tally3_queue_t queue );


magma_tally3_int_t
magma_tally3_vector_clag2z(
    magma_tally3_c_matrix x,
    magma_tally3_z_matrix *y,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_sparse_matrix_clag2z(
    magma_tally3_c_matrix A,
    magma_tally3_z_matrix *B,
    magma_tally3_queue_t queue );

void
magma_tally3blas_zlag2c_sparse(
    magma_tally3_int_t M, 
    magma_tally3_int_t N , 
    magma_tally3DoubleComplex_const_ptr dA, 
    magma_tally3_int_t lda, 
    magma_tally3FloatComplex_ptr dSA, 
    magma_tally3_int_t ldsa,
    magma_tally3_int_t *info,
    magma_tally3_queue_t queue );

void 
magma_tally3blas_clag2z_sparse(
    magma_tally3_int_t M, 
    magma_tally3_int_t N , 
    magma_tally3FloatComplex_const_ptr dSA, 
    magma_tally3_int_t ldsa, 
    magma_tally3DoubleComplex_ptr dA, 
    magma_tally3_int_t lda,
    magma_tally3_int_t *info,
    magma_tally3_queue_t queue );

void 
magma_tally3_zlag2c_CSR_DENSE(
    magma_tally3_z_matrix A,
    magma_tally3_c_matrix *B,
    magma_tally3_queue_t queue );

void 
magma_tally3_zlag2c_CSR_DENSE_alloc(
    magma_tally3_z_matrix A,
    magma_tally3_c_matrix *B,
    magma_tally3_queue_t queue );

void 
magma_tally3_zlag2c_CSR_DENSE_convert(
    magma_tally3_z_matrix A,
    magma_tally3_c_matrix *B,
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
magma_tally3_zcpgmres(
    magma_tally3_z_matrix A, 
    magma_tally3_z_matrix b, 
    magma_tally3_z_matrix *x,
    magma_tally3_z_solver_par *solver_par,
    magma_tally3_z_preconditioner *precond_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zcpbicgstab(
    magma_tally3_z_matrix A, 
    magma_tally3_z_matrix b, 
    magma_tally3_z_matrix *x,
    magma_tally3_z_solver_par *solver_par,
    magma_tally3_z_preconditioner *precond_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zcir(
    magma_tally3_z_matrix A, 
    magma_tally3_z_matrix b, 
    magma_tally3_z_matrix *x,
    magma_tally3_z_solver_par *solver_par, 
    magma_tally3_z_preconditioner *precond_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zcpir(
    magma_tally3_z_matrix A, 
    magma_tally3_z_matrix b, 
    magma_tally3_z_matrix *x,
    magma_tally3_z_solver_par *solver_par, 
    magma_tally3_z_preconditioner *precond_par,
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

#undef PRECISION_z
#endif /* MAGMA_tally3SPARSE_ZC_H */
