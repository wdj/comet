/*
-- MAGMA_minproduct (version 1.1) --
 Univ. of Tennessee, Knoxville
 Univ. of California, Berkeley
 Univ. of Colorado, Denver
 @date

 @precisions mixed zc -> ds
 @author Hartwig Anzt
*/

#ifndef MAGMA_minproductSPARSE_ZC_H
#define MAGMA_minproductSPARSE_ZC_H

#include "magma_minproduct_types.h"
#include "magma_minproductsparse_types.h"

#define PRECISION_z


/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_minproduct_SPARSE Matrix Descriptors
*/


#ifdef __cplusplus
extern "C" {
#endif

/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_minproduct_SPARSE Auxiliary functions
*/
magma_minproduct_int_t
magma_minproduct_vector_zlag2c(
    magma_minproduct_z_matrix x,
    magma_minproduct_c_matrix *y,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_sparse_matrix_zlag2c(
    magma_minproduct_z_matrix A,
    magma_minproduct_c_matrix *B,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t
magma_minproduct_vector_clag2z(
    magma_minproduct_c_matrix x,
    magma_minproduct_z_matrix *y,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_sparse_matrix_clag2z(
    magma_minproduct_c_matrix A,
    magma_minproduct_z_matrix *B,
    magma_minproduct_queue_t queue );

void
magma_minproductblas_zlag2c_sparse(
    magma_minproduct_int_t M, 
    magma_minproduct_int_t N , 
    magma_minproductDoubleComplex_const_ptr dA, 
    magma_minproduct_int_t lda, 
    magma_minproductFloatComplex_ptr dSA, 
    magma_minproduct_int_t ldsa,
    magma_minproduct_int_t *info,
    magma_minproduct_queue_t queue );

void 
magma_minproductblas_clag2z_sparse(
    magma_minproduct_int_t M, 
    magma_minproduct_int_t N , 
    magma_minproductFloatComplex_const_ptr dSA, 
    magma_minproduct_int_t ldsa, 
    magma_minproductDoubleComplex_ptr dA, 
    magma_minproduct_int_t lda,
    magma_minproduct_int_t *info,
    magma_minproduct_queue_t queue );

void 
magma_minproduct_zlag2c_CSR_DENSE(
    magma_minproduct_z_matrix A,
    magma_minproduct_c_matrix *B,
    magma_minproduct_queue_t queue );

void 
magma_minproduct_zlag2c_CSR_DENSE_alloc(
    magma_minproduct_z_matrix A,
    magma_minproduct_c_matrix *B,
    magma_minproduct_queue_t queue );

void 
magma_minproduct_zlag2c_CSR_DENSE_convert(
    magma_minproduct_z_matrix A,
    magma_minproduct_c_matrix *B,
    magma_minproduct_queue_t queue );

/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_minproduct_SPARSE function definitions / Data on CPU
*/


/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_minproduct_SPARSE function definitions / Data on CPU / Multi-GPU
*/

/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_minproduct_SPARSE function definitions / Data on GPU
*/
magma_minproduct_int_t
magma_minproduct_zcpgmres(
    magma_minproduct_z_matrix A, 
    magma_minproduct_z_matrix b, 
    magma_minproduct_z_matrix *x,
    magma_minproduct_z_solver_par *solver_par,
    magma_minproduct_z_preconditioner *precond_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zcpbicgstab(
    magma_minproduct_z_matrix A, 
    magma_minproduct_z_matrix b, 
    magma_minproduct_z_matrix *x,
    magma_minproduct_z_solver_par *solver_par,
    magma_minproduct_z_preconditioner *precond_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zcir(
    magma_minproduct_z_matrix A, 
    magma_minproduct_z_matrix b, 
    magma_minproduct_z_matrix *x,
    magma_minproduct_z_solver_par *solver_par, 
    magma_minproduct_z_preconditioner *precond_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zcpir(
    magma_minproduct_z_matrix A, 
    magma_minproduct_z_matrix b, 
    magma_minproduct_z_matrix *x,
    magma_minproduct_z_solver_par *solver_par, 
    magma_minproduct_z_preconditioner *precond_par,
    magma_minproduct_queue_t queue );


/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_minproduct_SPARSE utility function definitions
*/



/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_minproduct_SPARSE BLAS function definitions
*/



#ifdef __cplusplus
}
#endif

#undef PRECISION_z
#endif /* MAGMA_minproductSPARSE_ZC_H */
