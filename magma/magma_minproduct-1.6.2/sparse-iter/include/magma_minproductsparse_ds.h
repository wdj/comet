/*
-- MAGMA_minproduct (version 1.1) --
 Univ. of Tennessee, Knoxville
 Univ. of California, Berkeley
 Univ. of Colorado, Denver
 @date

 @generated from magma_minproductsparse_zc.h mixed zc -> ds, Fri May  1 21:08:05 2015
 @author Hartwig Anzt
*/

#ifndef MAGMA_minproductSPARSE_DS_H
#define MAGMA_minproductSPARSE_DS_H

#include "magma_minproduct_types.h"
#include "magma_minproductsparse_types.h"

#define PRECISION_d


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
magma_minproduct_vector_dlag2s(
    magma_minproduct_d_matrix x,
    magma_minproduct_s_matrix *y,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_sparse_matrix_dlag2s(
    magma_minproduct_d_matrix A,
    magma_minproduct_s_matrix *B,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t
magma_minproduct_vector_slag2d(
    magma_minproduct_s_matrix x,
    magma_minproduct_d_matrix *y,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_sparse_matrix_slag2d(
    magma_minproduct_s_matrix A,
    magma_minproduct_d_matrix *B,
    magma_minproduct_queue_t queue );

void
magma_minproductblas_dlag2s_sparse(
    magma_minproduct_int_t M, 
    magma_minproduct_int_t N , 
    magma_minproductDouble_const_ptr dA, 
    magma_minproduct_int_t lda, 
    magma_minproductFloat_ptr dSA, 
    magma_minproduct_int_t ldsa,
    magma_minproduct_int_t *info,
    magma_minproduct_queue_t queue );

void 
magma_minproductblas_slag2d_sparse(
    magma_minproduct_int_t M, 
    magma_minproduct_int_t N , 
    magma_minproductFloat_const_ptr dSA, 
    magma_minproduct_int_t ldsa, 
    magma_minproductDouble_ptr dA, 
    magma_minproduct_int_t lda,
    magma_minproduct_int_t *info,
    magma_minproduct_queue_t queue );

void 
magma_minproduct_dlag2s_CSR_DENSE(
    magma_minproduct_d_matrix A,
    magma_minproduct_s_matrix *B,
    magma_minproduct_queue_t queue );

void 
magma_minproduct_dlag2s_CSR_DENSE_alloc(
    magma_minproduct_d_matrix A,
    magma_minproduct_s_matrix *B,
    magma_minproduct_queue_t queue );

void 
magma_minproduct_dlag2s_CSR_DENSE_convert(
    magma_minproduct_d_matrix A,
    magma_minproduct_s_matrix *B,
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
magma_minproduct_dspgmres(
    magma_minproduct_d_matrix A, 
    magma_minproduct_d_matrix b, 
    magma_minproduct_d_matrix *x,
    magma_minproduct_d_solver_par *solver_par,
    magma_minproduct_d_preconditioner *precond_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dspbicgstab(
    magma_minproduct_d_matrix A, 
    magma_minproduct_d_matrix b, 
    magma_minproduct_d_matrix *x,
    magma_minproduct_d_solver_par *solver_par,
    magma_minproduct_d_preconditioner *precond_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dsir(
    magma_minproduct_d_matrix A, 
    magma_minproduct_d_matrix b, 
    magma_minproduct_d_matrix *x,
    magma_minproduct_d_solver_par *solver_par, 
    magma_minproduct_d_preconditioner *precond_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dspir(
    magma_minproduct_d_matrix A, 
    magma_minproduct_d_matrix b, 
    magma_minproduct_d_matrix *x,
    magma_minproduct_d_solver_par *solver_par, 
    magma_minproduct_d_preconditioner *precond_par,
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

#undef PRECISION_d
#endif /* MAGMA_minproductSPARSE_DS_H */
