/*
-- MAGMA_tally2 (version 1.1) --
 Univ. of Tennessee, Knoxville
 Univ. of California, Berkeley
 Univ. of Colorado, Denver
 @date

 @precisions mixed zc -> ds
 @author Hartwig Anzt
*/

#ifndef MAGMA_tally2SPARSE_ZC_H
#define MAGMA_tally2SPARSE_ZC_H

#include "magma_tally2_types.h"
#include "magma_tally2sparse_types.h"

#define PRECISION_z


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
magma_tally2_vector_zlag2c(
    magma_tally2_z_matrix x,
    magma_tally2_c_matrix *y,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_sparse_matrix_zlag2c(
    magma_tally2_z_matrix A,
    magma_tally2_c_matrix *B,
    magma_tally2_queue_t queue );


magma_tally2_int_t
magma_tally2_vector_clag2z(
    magma_tally2_c_matrix x,
    magma_tally2_z_matrix *y,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_sparse_matrix_clag2z(
    magma_tally2_c_matrix A,
    magma_tally2_z_matrix *B,
    magma_tally2_queue_t queue );

void
magma_tally2blas_zlag2c_sparse(
    magma_tally2_int_t M, 
    magma_tally2_int_t N , 
    magma_tally2DoubleComplex_const_ptr dA, 
    magma_tally2_int_t lda, 
    magma_tally2FloatComplex_ptr dSA, 
    magma_tally2_int_t ldsa,
    magma_tally2_int_t *info,
    magma_tally2_queue_t queue );

void 
magma_tally2blas_clag2z_sparse(
    magma_tally2_int_t M, 
    magma_tally2_int_t N , 
    magma_tally2FloatComplex_const_ptr dSA, 
    magma_tally2_int_t ldsa, 
    magma_tally2DoubleComplex_ptr dA, 
    magma_tally2_int_t lda,
    magma_tally2_int_t *info,
    magma_tally2_queue_t queue );

void 
magma_tally2_zlag2c_CSR_DENSE(
    magma_tally2_z_matrix A,
    magma_tally2_c_matrix *B,
    magma_tally2_queue_t queue );

void 
magma_tally2_zlag2c_CSR_DENSE_alloc(
    magma_tally2_z_matrix A,
    magma_tally2_c_matrix *B,
    magma_tally2_queue_t queue );

void 
magma_tally2_zlag2c_CSR_DENSE_convert(
    magma_tally2_z_matrix A,
    magma_tally2_c_matrix *B,
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
magma_tally2_zcpgmres(
    magma_tally2_z_matrix A, 
    magma_tally2_z_matrix b, 
    magma_tally2_z_matrix *x,
    magma_tally2_z_solver_par *solver_par,
    magma_tally2_z_preconditioner *precond_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zcpbicgstab(
    magma_tally2_z_matrix A, 
    magma_tally2_z_matrix b, 
    magma_tally2_z_matrix *x,
    magma_tally2_z_solver_par *solver_par,
    magma_tally2_z_preconditioner *precond_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zcir(
    magma_tally2_z_matrix A, 
    magma_tally2_z_matrix b, 
    magma_tally2_z_matrix *x,
    magma_tally2_z_solver_par *solver_par, 
    magma_tally2_z_preconditioner *precond_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zcpir(
    magma_tally2_z_matrix A, 
    magma_tally2_z_matrix b, 
    magma_tally2_z_matrix *x,
    magma_tally2_z_solver_par *solver_par, 
    magma_tally2_z_preconditioner *precond_par,
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

#undef PRECISION_z
#endif /* MAGMA_tally2SPARSE_ZC_H */
