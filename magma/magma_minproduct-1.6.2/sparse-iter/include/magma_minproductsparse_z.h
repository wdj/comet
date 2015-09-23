/*
 -- MAGMA_minproduct (version 1.1) --
 Univ. of Tennessee, Knoxville
 Univ. of California, Berkeley
 Univ. of Colorado, Denver
 @date

 @precisions normal z -> s d c
 @author Hartwig Anzt
*/

#ifndef MAGMA_minproductSPARSE_Z_H
#define MAGMA_minproductSPARSE_Z_H

#include "magma_minproduct_types.h"
#include "magma_minproductsparse_types.h"

#define PRECISION_z


#ifdef __cplusplus
extern "C" {
#endif


/* ////////////////////////////////////////////////////////////////////////////
 -- For backwards compatability, map old (1.6.1) to new (1.6.2) function names
*/

#define magma_minproduct_z_mtranspose  magma_minproduct_zmtranspose
#define magma_minproduct_z_mtransfer   magma_minproduct_zmtransfer
#define magma_minproduct_z_vtransfer   magma_minproduct_zmtransfer
#define magma_minproduct_z_mconvert    magma_minproduct_zmconvert
#define magma_minproduct_z_vinit       magma_minproduct_zvinit
#define magma_minproduct_z_vvisu       magma_minproduct_zprint_vector
#define magma_minproduct_z_vread       magma_minproduct_zvread
#define magma_minproduct_z_vspread     magma_minproduct_zvspread
#define magma_minproduct_z_mvisu       magma_minproduct_zprint_matrix
#define magma_minproduct_z_mfree       magma_minproduct_zmfree
#define magma_minproduct_z_vfree       magma_minproduct_zmfree
#define write_z_csr_mtx     magma_minproduct_zwrite_csr_mtx
#define write_z_csrtomtx    magma_minproduct_zwrite_csrtomtx
#define print_z_csr         magma_minproduct_zprint_csr_mtx


/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_minproduct_SPARSE Auxiliary functions
*/


magma_minproduct_int_t
magma_minproduct_zparse_opts( 
    int argc, 
    char** argv, 
    magma_minproduct_zopts *opts, 
    int *matrices, 
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
read_z_csr_from_binary( 
    magma_minproduct_int_t* n_row, 
    magma_minproduct_int_t* n_col, 
    magma_minproduct_int_t* nnz, 
    magma_minproductDoubleComplex **val, 
    magma_minproduct_index_t **row, 
    magma_minproduct_index_t **col,
    const char * filename,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
read_z_csr_from_mtx( 
    magma_minproduct_storage_t *type, 
    magma_minproduct_location_t *location,
    magma_minproduct_int_t* n_row, 
    magma_minproduct_int_t* n_col, 
    magma_minproduct_int_t* nnz, 
    magma_minproductDoubleComplex **val, 
    magma_minproduct_index_t **row, 
    magma_minproduct_index_t **col, 
    const char *filename,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_z_csr_mtx( 
    magma_minproduct_z_matrix *A, 
    const char *filename,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_zcsrset( 
    magma_minproduct_int_t m, 
    magma_minproduct_int_t n, 
    magma_minproduct_index_t *row, 
    magma_minproduct_index_t *col, 
    magma_minproductDoubleComplex *val,
    magma_minproduct_z_matrix *A,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_zcsrget( 
    magma_minproduct_z_matrix A,
    magma_minproduct_int_t *m, 
    magma_minproduct_int_t *n, 
    magma_minproduct_index_t **row, 
    magma_minproduct_index_t **col, 
    magma_minproductDoubleComplex **val,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t 
magma_minproduct_zvset( 
    magma_minproduct_int_t m, 
    magma_minproduct_int_t n, 
    magma_minproductDoubleComplex *val,
    magma_minproduct_z_matrix *v,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_zvget( 
    magma_minproduct_z_matrix v,
    magma_minproduct_int_t *m, 
    magma_minproduct_int_t *n, 
    magma_minproductDoubleComplex **val,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_zvset_dev( 
    magma_minproduct_int_t m, 
    magma_minproduct_int_t n, 
    magma_minproductDoubleComplex_ptr val,
    magma_minproduct_z_matrix *v,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_zvget_dev( 
    magma_minproduct_z_matrix v,
    magma_minproduct_int_t *m, 
    magma_minproduct_int_t *n, 
    magma_minproductDoubleComplex_ptr *val,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t 
magma_minproduct_z_csr_mtxsymm( 
    magma_minproduct_z_matrix *A, 
    const char *filename,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_z_csr_compressor( 
    magma_minproductDoubleComplex ** val, 
    magma_minproduct_index_t ** row, 
    magma_minproduct_index_t ** col, 
    magma_minproductDoubleComplex ** valn, 
    magma_minproduct_index_t ** rown, 
    magma_minproduct_index_t ** coln, 
    magma_minproduct_int_t *n,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zmcsrcompressor( 
    magma_minproduct_z_matrix *A,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zmcsrcompressor_gpu( 
    magma_minproduct_z_matrix *A,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zvtranspose( 
    magma_minproduct_z_matrix x,
    magma_minproduct_z_matrix *y,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_z_cucsrtranspose( 
    magma_minproduct_z_matrix A, 
    magma_minproduct_z_matrix *B,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
z_transpose_csr( 
    magma_minproduct_int_t n_rows, 
    magma_minproduct_int_t n_cols, 
    magma_minproduct_int_t nnz,
    magma_minproductDoubleComplex *val, 
    magma_minproduct_index_t *row, 
    magma_minproduct_index_t *col, 
    magma_minproduct_int_t *new_n_rows, 
    magma_minproduct_int_t *new_n_cols, 
    magma_minproduct_int_t *new_nnz, 
    magma_minproductDoubleComplex **new_val, 
    magma_minproduct_index_t **new_row, 
    magma_minproduct_index_t **new_col,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zcsrsplit( 
    magma_minproduct_int_t bsize,
    magma_minproduct_z_matrix A,
    magma_minproduct_z_matrix *D,
    magma_minproduct_z_matrix *R,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zmscale( 
    magma_minproduct_z_matrix *A, 
    magma_minproduct_scale_t scaling,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_zmdiff( 
    magma_minproduct_z_matrix A, 
    magma_minproduct_z_matrix B, 
 real_Double_t *res,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zmdiagadd( 
    magma_minproduct_z_matrix *A, 
    magma_minproductDoubleComplex add,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_zmsort( 
    magma_minproduct_z_matrix *A,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zindexsort(
    magma_minproduct_index_t *x, 
    magma_minproduct_int_t first,
    magma_minproduct_int_t last,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zdomainoverlap(
    magma_minproduct_index_t num_rows,
    magma_minproduct_index_t *num_indices,
    magma_minproduct_index_t *rowptr,
    magma_minproduct_index_t *colidx,
    magma_minproduct_index_t *x,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zsymbilu( 
    magma_minproduct_z_matrix *A, 
    magma_minproduct_int_t levels,
    magma_minproduct_z_matrix *L,
    magma_minproduct_z_matrix *U,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t 
magma_minproduct_zwrite_csr_mtx( 
    magma_minproduct_z_matrix A,
    magma_minproduct_order_t MajorType,
 const char *filename,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_zwrite_csrtomtx( 
    magma_minproduct_z_matrix A,
    const char *filename,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_zprint_csr( 
    magma_minproduct_int_t n_row, 
    magma_minproduct_int_t n_col, 
    magma_minproduct_int_t nnz, 
    magma_minproductDoubleComplex **val, 
    magma_minproduct_index_t **row, 
    magma_minproduct_index_t **col,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_zprint_csr_mtx( 
    magma_minproduct_int_t n_row, 
    magma_minproduct_int_t n_col, 
    magma_minproduct_int_t nnz, 
    magma_minproductDoubleComplex **val, 
    magma_minproduct_index_t **row, 
    magma_minproduct_index_t **col, 
    magma_minproduct_order_t MajorType,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t 
magma_minproduct_zmtranspose(
    magma_minproduct_z_matrix A, 
    magma_minproduct_z_matrix *B,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t 
magma_minproduct_zmtransfer(
    magma_minproduct_z_matrix A, 
    magma_minproduct_z_matrix *B, 
    magma_minproduct_location_t src, 
    magma_minproduct_location_t dst,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_zmconvert(
    magma_minproduct_z_matrix A, 
    magma_minproduct_z_matrix *B, 
    magma_minproduct_storage_t old_format, 
    magma_minproduct_storage_t new_format,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t
magma_minproduct_zvinit(
    magma_minproduct_z_matrix *x, 
    magma_minproduct_location_t memory_location,
    magma_minproduct_int_t num_rows, 
    magma_minproduct_int_t num_cols,
    magma_minproductDoubleComplex values,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zprint_vector(
    magma_minproduct_z_matrix x, 
    magma_minproduct_int_t offset, 
    magma_minproduct_int_t displaylength,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zvread(
    magma_minproduct_z_matrix *x, 
    magma_minproduct_int_t length,
    char * filename,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zvspread(
    magma_minproduct_z_matrix *x, 
    const char * filename,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zprint_matrix(
    magma_minproduct_z_matrix A,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_zdiameter(
    magma_minproduct_z_matrix *A,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_zrowentries(
    magma_minproduct_z_matrix *A,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zmfree(
    magma_minproduct_z_matrix *A,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zresidual(
    magma_minproduct_z_matrix A, 
    magma_minproduct_z_matrix b, 
    magma_minproduct_z_matrix x, 
    double *res,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zresidualvec(
    magma_minproduct_z_matrix A,
    magma_minproduct_z_matrix b,
    magma_minproduct_z_matrix x,
    magma_minproduct_z_matrix *r,
    double *res,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zmgenerator(
    magma_minproduct_int_t n,
    magma_minproduct_int_t offdiags,
    magma_minproduct_index_t *diag_offset,
    magma_minproductDoubleComplex *diag_vals,
    magma_minproduct_z_matrix *A,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zm_27stencil(
    magma_minproduct_int_t n,
    magma_minproduct_z_matrix *A,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zm_5stencil(
    magma_minproduct_int_t n,
    magma_minproduct_z_matrix *A,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zsolverinfo(
    magma_minproduct_z_solver_par *solver_par, 
    magma_minproduct_z_preconditioner *precond_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zsolverinfo_init(
    magma_minproduct_z_solver_par *solver_par, 
    magma_minproduct_z_preconditioner *precond,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zeigensolverinfo_init(
    magma_minproduct_z_solver_par *solver_par,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t
magma_minproduct_zsolverinfo_free(
    magma_minproduct_z_solver_par *solver_par, 
    magma_minproduct_z_preconditioner *precond,
    magma_minproduct_queue_t queue );



/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_minproduct_SPARSE iterative incomplete factorizations
*/


magma_minproduct_int_t
magma_minproduct_ziterilusetup( 
    magma_minproduct_z_matrix A, 
    magma_minproduct_z_matrix b,                                 
    magma_minproduct_z_preconditioner *precond,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zitericsetup( 
    magma_minproduct_z_matrix A, 
    magma_minproduct_z_matrix b, 
    magma_minproduct_z_preconditioner *precond,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zitericupdate( 
    magma_minproduct_z_matrix A, 
    magma_minproduct_z_preconditioner *precond, 
    magma_minproduct_int_t updates,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zapplyiteric_l( 
    magma_minproduct_z_matrix b, 
    magma_minproduct_z_matrix *x, 
    magma_minproduct_z_preconditioner *precond,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zapplyiteric_r( 
    magma_minproduct_z_matrix b, 
    magma_minproduct_z_matrix *x, 
    magma_minproduct_z_preconditioner *precond,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t
magma_minproduct_ziterilu_csr( 
    magma_minproduct_z_matrix A,
    magma_minproduct_z_matrix L,
    magma_minproduct_z_matrix U,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_ziteric_csr( 
    magma_minproduct_z_matrix A,
    magma_minproduct_z_matrix A_CSR,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_zfrobenius( 
    magma_minproduct_z_matrix A, 
    magma_minproduct_z_matrix B, 
    real_Double_t *res,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_znonlinres(   
    magma_minproduct_z_matrix A, 
    magma_minproduct_z_matrix L,
    magma_minproduct_z_matrix U, 
    magma_minproduct_z_matrix *LU, 
    real_Double_t *res,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_zilures(   
    magma_minproduct_z_matrix A, 
    magma_minproduct_z_matrix L,
    magma_minproduct_z_matrix U, 
    magma_minproduct_z_matrix *LU, 
    real_Double_t *res,
    real_Double_t *nonlinres,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_zicres(       
    magma_minproduct_z_matrix A, 
    magma_minproduct_z_matrix C,
    magma_minproduct_z_matrix CT, 
    magma_minproduct_z_matrix *LU, 
    real_Double_t *res,
    real_Double_t *nonlinres,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_zinitguess( 
    magma_minproduct_z_matrix A, 
    magma_minproduct_z_matrix *L, 
    magma_minproduct_z_matrix *U,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_zinitrecursiveLU( 
    magma_minproduct_z_matrix A, 
    magma_minproduct_z_matrix *B,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_zmLdiagadd( 
    magma_minproduct_z_matrix *L,
    magma_minproduct_queue_t queue );




/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_minproduct_SPARSE function definitions / Data on CPU
*/




/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_minproduct_SPARSE function definitions / Data on CPU / Multi-GPU
*/

/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_minproduct_SPARSE solvers (Data on GPU)
*/

magma_minproduct_int_t 
magma_minproduct_zcg(
    magma_minproduct_z_matrix A, magma_minproduct_z_matrix b, 
    magma_minproduct_z_matrix *x, magma_minproduct_z_solver_par *solver_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_zcg_res(
    magma_minproduct_z_matrix A, magma_minproduct_z_matrix b, 
    magma_minproduct_z_matrix *x, magma_minproduct_z_solver_par *solver_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_zcg_merge(
    magma_minproduct_z_matrix A, magma_minproduct_z_matrix b, 
    magma_minproduct_z_matrix *x, magma_minproduct_z_solver_par *solver_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_zgmres(
    magma_minproduct_z_matrix A, magma_minproduct_z_matrix b, 
    magma_minproduct_z_matrix *x, magma_minproduct_z_solver_par *solver_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zbicgstab(
    magma_minproduct_z_matrix A, magma_minproduct_z_matrix b, magma_minproduct_z_matrix *x, 
    magma_minproduct_z_solver_par *solver_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zbicgstab_merge(
    magma_minproduct_z_matrix A, magma_minproduct_z_matrix b, 
    magma_minproduct_z_matrix *x, magma_minproduct_z_solver_par *solver_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zbicgstab_merge2(
    magma_minproduct_z_matrix A, magma_minproduct_z_matrix b, 
    magma_minproduct_z_matrix *x, magma_minproduct_z_solver_par *solver_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zpcg(
    magma_minproduct_z_matrix A, magma_minproduct_z_matrix b, 
    magma_minproduct_z_matrix *x, magma_minproduct_z_solver_par *solver_par, 
    magma_minproduct_z_preconditioner *precond_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zbpcg(
    magma_minproduct_z_matrix A, magma_minproduct_z_matrix b, 
    magma_minproduct_z_matrix *x, magma_minproduct_z_solver_par *solver_par, 
    magma_minproduct_z_preconditioner *precond_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zpbicgstab(
    magma_minproduct_z_matrix A, magma_minproduct_z_matrix b, 
    magma_minproduct_z_matrix *x, magma_minproduct_z_solver_par *solver_par, 
    magma_minproduct_z_preconditioner *precond_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zpgmres(
    magma_minproduct_z_matrix A, magma_minproduct_z_matrix b, 
    magma_minproduct_z_matrix *x, magma_minproduct_z_solver_par *solver_par, 
    magma_minproduct_z_preconditioner *precond_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zfgmres(
    magma_minproduct_z_matrix A, magma_minproduct_z_matrix b, 
    magma_minproduct_z_matrix *x, magma_minproduct_z_solver_par *solver_par, 
    magma_minproduct_z_preconditioner *precond_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zjacobi(
    magma_minproduct_z_matrix A, magma_minproduct_z_matrix b, 
    magma_minproduct_z_matrix *x, magma_minproduct_z_solver_par *solver_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zjacobidomainoverlap(
    magma_minproduct_z_matrix A, 
    magma_minproduct_z_matrix b, 
    magma_minproduct_z_matrix *x,  
    magma_minproduct_z_solver_par *solver_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zbaiter(
    magma_minproduct_z_matrix A, magma_minproduct_z_matrix b, 
    magma_minproduct_z_matrix *x, magma_minproduct_z_solver_par *solver_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_ziterref(
    magma_minproduct_z_matrix A, magma_minproduct_z_matrix b, 
    magma_minproduct_z_matrix *x, magma_minproduct_z_solver_par *solver_par, 
    magma_minproduct_z_preconditioner *precond_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zilu(
    magma_minproduct_z_matrix A, magma_minproduct_z_matrix b, 
    magma_minproduct_z_matrix *x, magma_minproduct_z_solver_par *solver_par, 
    magma_minproduct_int_t *ipiv,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zbcsrlu(
    magma_minproduct_z_matrix A, magma_minproduct_z_matrix b, 
    magma_minproduct_z_matrix *x, magma_minproduct_z_solver_par *solver_par,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t
magma_minproduct_zbcsrlutrf(
    magma_minproduct_z_matrix A, 
    magma_minproduct_z_matrix *M,
    magma_minproduct_int_t *ipiv, 
    magma_minproduct_int_t version,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zbcsrlusv(
    magma_minproduct_z_matrix A, magma_minproduct_z_matrix b, 
    magma_minproduct_z_matrix *x, magma_minproduct_z_solver_par *solver_par, 
    magma_minproduct_int_t *ipiv,
    magma_minproduct_queue_t queue );



magma_minproduct_int_t
magma_minproduct_zilucg(
    magma_minproduct_z_matrix A, magma_minproduct_z_matrix b, 
    magma_minproduct_z_matrix *x, magma_minproduct_z_solver_par *solver_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zilugmres(
    magma_minproduct_z_matrix A, magma_minproduct_z_matrix b, 
    magma_minproduct_z_matrix *x, magma_minproduct_z_solver_par *solver_par,
    magma_minproduct_queue_t queue ); 


magma_minproduct_int_t
magma_minproduct_zlobpcg_shift(
    magma_minproduct_int_t num_rows,
    magma_minproduct_int_t num_vecs, 
    magma_minproduct_int_t shift,
    magma_minproductDoubleComplex_ptr x,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zlobpcg_res(
    magma_minproduct_int_t num_rows,
    magma_minproduct_int_t num_vecs, 
    double *evalues, 
    magma_minproductDoubleComplex_ptr X,
    magma_minproductDoubleComplex_ptr R, 
    double *res,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zlobpcg_maxpy(
    magma_minproduct_int_t num_rows,
    magma_minproduct_int_t num_vecs, 
    magma_minproductDoubleComplex_ptr X,
    magma_minproductDoubleComplex_ptr Y,
    magma_minproduct_queue_t queue );


/*/////////////////////////////////////////////////////////////////////////////
 -- MAGMA_minproduct_SPARSE eigensolvers (Data on GPU)
*/
magma_minproduct_int_t
magma_minproduct_zlobpcg(
    magma_minproduct_z_matrix A, 
    magma_minproduct_z_solver_par *solver_par,
    magma_minproduct_z_preconditioner *precond_par, 
    magma_minproduct_queue_t queue );




/*/////////////////////////////////////////////////////////////////////////////
 -- MAGMA_minproduct_SPARSE preconditioners (Data on GPU)
*/

magma_minproduct_int_t
magma_minproduct_zjacobisetup(
    magma_minproduct_z_matrix A, 
    magma_minproduct_z_matrix b, 
    magma_minproduct_z_matrix *M, 
    magma_minproduct_z_matrix *c,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zjacobisetup_matrix(
    magma_minproduct_z_matrix A, 
    magma_minproduct_z_matrix *M, 
    magma_minproduct_z_matrix *d,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zjacobisetup_vector(
    magma_minproduct_z_matrix b, 
    magma_minproduct_z_matrix d, 
    magma_minproduct_z_matrix *c,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zjacobiiter(
    magma_minproduct_z_matrix M, 
    magma_minproduct_z_matrix c, 
    magma_minproduct_z_matrix *x, 
    magma_minproduct_z_solver_par *solver_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zjacobiiter_precond( 
    magma_minproduct_z_matrix M, 
    magma_minproduct_z_matrix *x, 
    magma_minproduct_z_solver_par *solver_par, 
    magma_minproduct_z_preconditioner *precond,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zjacobiiter_sys(
    magma_minproduct_z_matrix A, 
    magma_minproduct_z_matrix b, 
    magma_minproduct_z_matrix d, 
    magma_minproduct_z_matrix t, 
    magma_minproduct_z_matrix *x,  
    magma_minproduct_z_solver_par *solver_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zpastixsetup(
    magma_minproduct_z_matrix A, magma_minproduct_z_matrix b,
    magma_minproduct_z_preconditioner *precond,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t
magma_minproduct_zapplypastix(
    magma_minproduct_z_matrix b, magma_minproduct_z_matrix *x, 
    magma_minproduct_z_preconditioner *precond,
    magma_minproduct_queue_t queue );


// custom preconditioner
magma_minproduct_int_t
magma_minproduct_zapplycustomprecond_l(
    magma_minproduct_z_matrix b, 
    magma_minproduct_z_matrix *x, 
    magma_minproduct_z_preconditioner *precond,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zapplycustomprecond_r(
    magma_minproduct_z_matrix b, 
    magma_minproduct_z_matrix *x, 
    magma_minproduct_z_preconditioner *precond,
    magma_minproduct_queue_t queue );


// CUSPARSE preconditioner

magma_minproduct_int_t
magma_minproduct_zcuilusetup(
    magma_minproduct_z_matrix A, magma_minproduct_z_preconditioner *precond,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zapplycuilu_l(
    magma_minproduct_z_matrix b, magma_minproduct_z_matrix *x, 
    magma_minproduct_z_preconditioner *precond,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zapplycuilu_r(
    magma_minproduct_z_matrix b, magma_minproduct_z_matrix *x, 
    magma_minproduct_z_preconditioner *precond,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t
magma_minproduct_zcuiccsetup(
    magma_minproduct_z_matrix A, magma_minproduct_z_preconditioner *precond,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zapplycuicc_l(
    magma_minproduct_z_matrix b, magma_minproduct_z_matrix *x, 
    magma_minproduct_z_preconditioner *precond,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zapplycuicc_r(
    magma_minproduct_z_matrix b, magma_minproduct_z_matrix *x, 
    magma_minproduct_z_preconditioner *precond,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t
magma_minproduct_zcumilusetup(
    magma_minproduct_z_matrix A, 
    magma_minproduct_z_preconditioner *precond,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zcumilugeneratesolverinfo(
    magma_minproduct_z_preconditioner *precond,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zapplycumilu_l(
    magma_minproduct_z_matrix b, 
    magma_minproduct_z_matrix *x, 
    magma_minproduct_z_preconditioner *precond,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zapplycumilu_r(
    magma_minproduct_z_matrix b, 
    magma_minproduct_z_matrix *x, 
    magma_minproduct_z_preconditioner *precond,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t
magma_minproduct_zcumiccsetup(
    magma_minproduct_z_matrix A, 
    magma_minproduct_z_preconditioner *precond,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zcumicgeneratesolverinfo(
    magma_minproduct_z_preconditioner *precond,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zapplycumicc_l(
    magma_minproduct_z_matrix b, 
    magma_minproduct_z_matrix *x, 
    magma_minproduct_z_preconditioner *precond,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zapplycumicc_r(
    magma_minproduct_z_matrix b, 
    magma_minproduct_z_matrix *x, 
    magma_minproduct_z_preconditioner *precond,
    magma_minproduct_queue_t queue );


// block-asynchronous iteration

magma_minproduct_int_t
magma_minproduct_zbajac_csr(
    magma_minproduct_int_t localiters,
    magma_minproduct_z_matrix D,
    magma_minproduct_z_matrix R,
    magma_minproduct_z_matrix b,
    magma_minproduct_z_matrix *x,
    magma_minproduct_queue_t queue );

/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_minproduct_SPARSE utility function definitions
*/

magma_minproduct_int_t
magma_minproduct_z_spmv(
    magma_minproductDoubleComplex alpha, 
    magma_minproduct_z_matrix A, 
    magma_minproduct_z_matrix x, 
    magma_minproductDoubleComplex beta, 
    magma_minproduct_z_matrix y,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zcustomspmv(
    magma_minproductDoubleComplex alpha, 
    magma_minproduct_z_matrix x, 
    magma_minproductDoubleComplex beta, 
    magma_minproduct_z_matrix y,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_z_spmv_shift(
    magma_minproductDoubleComplex alpha, 
    magma_minproduct_z_matrix A, 
    magma_minproductDoubleComplex lambda,
    magma_minproduct_z_matrix x, 
    magma_minproductDoubleComplex beta, 
    magma_minproduct_int_t offset, 
    magma_minproduct_int_t blocksize,
    magma_minproductIndex_ptr dadd_vecs, 
    magma_minproduct_z_matrix y,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zcuspmm(
    magma_minproduct_z_matrix A, 
    magma_minproduct_z_matrix B, 
    magma_minproduct_z_matrix *AB,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_z_spmm(
    magma_minproductDoubleComplex alpha, 
    magma_minproduct_z_matrix A,
    magma_minproduct_z_matrix B,
    magma_minproduct_z_matrix *C,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zsymbilu( 
    magma_minproduct_z_matrix *A, 
    magma_minproduct_int_t levels, 
    magma_minproduct_z_matrix *L, 
    magma_minproduct_z_matrix *U,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zcuspaxpy(
    magma_minproductDoubleComplex_ptr alpha, magma_minproduct_z_matrix A, 
    magma_minproductDoubleComplex_ptr beta, magma_minproduct_z_matrix B, 
    magma_minproduct_z_matrix *AB,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_z_precond(
    magma_minproduct_z_matrix A, 
    magma_minproduct_z_matrix b, magma_minproduct_z_matrix *x,
    magma_minproduct_z_preconditioner *precond,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_z_solver(
    magma_minproduct_z_matrix A, magma_minproduct_z_matrix b, 
    magma_minproduct_z_matrix *x, magma_minproduct_zopts *zopts,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_z_precondsetup(
    magma_minproduct_z_matrix A, magma_minproduct_z_matrix b, 
    magma_minproduct_z_preconditioner *precond,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_z_applyprecond(
    magma_minproduct_z_matrix A, magma_minproduct_z_matrix b, 
    magma_minproduct_z_matrix *x, magma_minproduct_z_preconditioner *precond,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t
magma_minproduct_z_applyprecond_left(
    magma_minproduct_z_matrix A, magma_minproduct_z_matrix b, 
    magma_minproduct_z_matrix *x, magma_minproduct_z_preconditioner *precond,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t
magma_minproduct_z_applyprecond_right(
    magma_minproduct_z_matrix A, magma_minproduct_z_matrix b, 
    magma_minproduct_z_matrix *x, magma_minproduct_z_preconditioner *precond,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_z_initP2P(
    magma_minproduct_int_t *bandwidth_benchmark,
    magma_minproduct_int_t *num_gpus,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zcompact(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda,
    double *dnorms, double tol, 
    magma_minproduct_int_t *activeMask, magma_minproduct_int_t *cBlockSize,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zcompactActive(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex_ptr dA, magma_minproduct_int_t ldda, 
    magma_minproduct_int_t *active,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zmlumerge(    
    magma_minproduct_z_matrix L, 
    magma_minproduct_z_matrix U,
    magma_minproduct_z_matrix *A, 
    magma_minproduct_queue_t queue );


/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_minproduct_SPARSE BLAS function definitions
*/
magma_minproduct_int_t 
magma_minproduct_zgecsrmv(
    magma_minproduct_trans_t transA,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex alpha,
    magma_minproductDoubleComplex_ptr dval,
    magma_minproductIndex_ptr drowptr,
    magma_minproductIndex_ptr dcolind,
    magma_minproductDoubleComplex_ptr dx,
    magma_minproductDoubleComplex beta,
    magma_minproductDoubleComplex_ptr dy,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_zgecsrmv_shift(
    magma_minproduct_trans_t transA,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex alpha,
    magma_minproductDoubleComplex lambda,
    magma_minproductDoubleComplex_ptr dval,
    magma_minproductIndex_ptr drowptr,
    magma_minproductIndex_ptr dcolind,
    magma_minproductDoubleComplex_ptr dx,
    magma_minproductDoubleComplex beta,
    int offset,
    int blocksize,
    magma_minproductIndex_ptr dadd_rows,
    magma_minproductDoubleComplex_ptr dy,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_zmgecsrmv(
    magma_minproduct_trans_t transA,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproduct_int_t num_vecs,
    magma_minproductDoubleComplex alpha,
    magma_minproductDoubleComplex_ptr dval,
    magma_minproductIndex_ptr drowptr,
    magma_minproductIndex_ptr dcolind,
    magma_minproductDoubleComplex_ptr dx,
    magma_minproductDoubleComplex beta,
    magma_minproductDoubleComplex_ptr dy,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_zgeellmv(
    magma_minproduct_trans_t transA,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproduct_int_t nnz_per_row,
    magma_minproductDoubleComplex alpha,
    magma_minproductDoubleComplex_ptr dval,
    magma_minproductIndex_ptr dcolind,
    magma_minproductDoubleComplex_ptr dx,
    magma_minproductDoubleComplex beta,
    magma_minproductDoubleComplex_ptr dy,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_zgeellmv_shift(
    magma_minproduct_trans_t transA,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproduct_int_t nnz_per_row,
    magma_minproductDoubleComplex alpha,
    magma_minproductDoubleComplex lambda,
    magma_minproductDoubleComplex_ptr dval,
    magma_minproductIndex_ptr dcolind,
    magma_minproductDoubleComplex_ptr dx,
    magma_minproductDoubleComplex beta,
    int offset,
    int blocksize,
    magma_minproductIndex_ptr dadd_rows,
    magma_minproductDoubleComplex_ptr dy,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t 
magma_minproduct_zmgeellmv(
    magma_minproduct_trans_t transA,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproduct_int_t num_vecs,
    magma_minproduct_int_t nnz_per_row,
    magma_minproductDoubleComplex alpha,
    magma_minproductDoubleComplex_ptr dval,
    magma_minproductIndex_ptr dcolind,
    magma_minproductDoubleComplex_ptr dx,
    magma_minproductDoubleComplex beta,
    magma_minproductDoubleComplex_ptr dy,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t 
magma_minproduct_zgeelltmv(
    magma_minproduct_trans_t transA,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproduct_int_t nnz_per_row,
    magma_minproductDoubleComplex alpha,
    magma_minproductDoubleComplex_ptr dval,
    magma_minproductIndex_ptr dcolind,
    magma_minproductDoubleComplex_ptr dx,
    magma_minproductDoubleComplex beta,
    magma_minproductDoubleComplex_ptr dy,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_zgeelltmv_shift(
    magma_minproduct_trans_t transA,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproduct_int_t nnz_per_row,
    magma_minproductDoubleComplex alpha,
    magma_minproductDoubleComplex lambda,
    magma_minproductDoubleComplex_ptr dval,
    magma_minproductIndex_ptr dcolind,
    magma_minproductDoubleComplex_ptr dx,
    magma_minproductDoubleComplex beta,
    int offset,
    int blocksize,
    magma_minproductIndex_ptr dadd_rows,
    magma_minproductDoubleComplex_ptr dy,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t 
magma_minproduct_zmgeelltmv(
    magma_minproduct_trans_t transA,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproduct_int_t num_vecs,
    magma_minproduct_int_t nnz_per_row,
    magma_minproductDoubleComplex alpha,
    magma_minproductDoubleComplex_ptr dval,
    magma_minproductIndex_ptr dcolind,
    magma_minproductDoubleComplex_ptr dx,
    magma_minproductDoubleComplex beta,
    magma_minproductDoubleComplex_ptr dy,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_zgeellrtmv(
    magma_minproduct_trans_t transA,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproduct_int_t nnz_per_row,
    magma_minproductDoubleComplex alpha,
    magma_minproductDoubleComplex_ptr dval,
    magma_minproductIndex_ptr dcolind,
    magma_minproductIndex_ptr drowlength,
    magma_minproductDoubleComplex_ptr dx,
    magma_minproductDoubleComplex beta,
    magma_minproductDoubleComplex_ptr dy,
    magma_minproduct_int_t num_threads,
    magma_minproduct_int_t threads_per_row,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_zgesellcmv(
    magma_minproduct_trans_t transA,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproduct_int_t blocksize,
    magma_minproduct_int_t slices,
    magma_minproduct_int_t alignment,
    magma_minproductDoubleComplex alpha,
    magma_minproductDoubleComplex_ptr dval,
    magma_minproductIndex_ptr dcolind,
    magma_minproductIndex_ptr drowptr,
    magma_minproductDoubleComplex_ptr dx,
    magma_minproductDoubleComplex beta,
    magma_minproductDoubleComplex_ptr dy,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zgesellpmv(
    magma_minproduct_trans_t transA,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproduct_int_t blocksize,
    magma_minproduct_int_t slices,
    magma_minproduct_int_t alignment,
    magma_minproductDoubleComplex alpha,
    magma_minproductDoubleComplex_ptr dval,
    magma_minproductIndex_ptr dcolind,
    magma_minproductIndex_ptr drowptr,
    magma_minproductDoubleComplex_ptr dx,
    magma_minproductDoubleComplex beta,
    magma_minproductDoubleComplex_ptr dy,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zmgesellpmv(
    magma_minproduct_trans_t transA,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproduct_int_t num_vecs,
    magma_minproduct_int_t blocksize,
    magma_minproduct_int_t slices,
    magma_minproduct_int_t alignment,
    magma_minproductDoubleComplex alpha,
    magma_minproductDoubleComplex_ptr dval,
    magma_minproductIndex_ptr dcolind,
    magma_minproductIndex_ptr drowptr,
    magma_minproductDoubleComplex_ptr dx,
    magma_minproductDoubleComplex beta,
    magma_minproductDoubleComplex_ptr dy,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zmgesellpmv_blocked(
    magma_minproduct_trans_t transA,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproduct_int_t num_vecs,
    magma_minproduct_int_t blocksize,
    magma_minproduct_int_t slices,
    magma_minproduct_int_t alignment,
    magma_minproductDoubleComplex alpha,
    magma_minproductDoubleComplex_ptr dval,
    magma_minproductIndex_ptr dcolind,
    magma_minproductIndex_ptr drowptr,
    magma_minproductDoubleComplex_ptr dx,
    magma_minproductDoubleComplex beta,
    magma_minproductDoubleComplex_ptr dy,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t
magma_minproduct_zmergedgs(
    magma_minproduct_int_t n, 
    magma_minproduct_int_t ldh,
    magma_minproduct_int_t k, 
    magma_minproductDoubleComplex_ptr dv, 
    magma_minproductDoubleComplex_ptr dr,
    magma_minproductDoubleComplex_ptr dskp,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zcopyscale(    
    int n, 
    int k,
    magma_minproductDoubleComplex_ptr dr, 
    magma_minproductDoubleComplex_ptr dv,
    magma_minproductDoubleComplex_ptr dskp,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dznrm2scale(    
    int m, 
    magma_minproductDoubleComplex_ptr dr,    
    int lddr, 
    magma_minproductDoubleComplex *drnorm,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t
magma_minproduct_zjacobisetup_vector_gpu(
    int num_rows, 
    magma_minproduct_z_matrix b, 
    magma_minproduct_z_matrix d, 
    magma_minproduct_z_matrix c,
    magma_minproduct_z_matrix *x,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t
magma_minproduct_zjacobi_diagscal(    
    int num_rows, 
    magma_minproduct_z_matrix d, 
    magma_minproduct_z_matrix b, 
    magma_minproduct_z_matrix *c,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zjacobiupdate(
    magma_minproduct_z_matrix t, 
    magma_minproduct_z_matrix b, 
    magma_minproduct_z_matrix d, 
    magma_minproduct_z_matrix *x,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zjacobispmvupdate(
    magma_minproduct_int_t maxiter,
    magma_minproduct_z_matrix A, 
    magma_minproduct_z_matrix t, 
    magma_minproduct_z_matrix b, 
    magma_minproduct_z_matrix d, 
    magma_minproduct_z_matrix *x,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zjacobispmvupdate_bw(
    magma_minproduct_int_t maxiter,
    magma_minproduct_z_matrix A, 
    magma_minproduct_z_matrix t, 
    magma_minproduct_z_matrix b, 
    magma_minproduct_z_matrix d, 
    magma_minproduct_z_matrix *x,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zjacobispmvupdateselect(
    magma_minproduct_int_t maxiter,
    magma_minproduct_int_t num_updates,
    magma_minproduct_index_t *indices,
    magma_minproduct_z_matrix A,
    magma_minproduct_z_matrix t, 
    magma_minproduct_z_matrix b, 
    magma_minproduct_z_matrix d, 
    magma_minproduct_z_matrix tmp, 
    magma_minproduct_z_matrix *x,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zjacobisetup_diagscal(
    magma_minproduct_z_matrix A, magma_minproduct_z_matrix *d,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t
magma_minproduct_zbicgmerge1(    
    int n, 
    magma_minproductDoubleComplex_ptr dskp,
    magma_minproductDoubleComplex_ptr dv, 
    magma_minproductDoubleComplex_ptr dr, 
    magma_minproductDoubleComplex_ptr dp,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t
magma_minproduct_zbicgmerge2(
    int n, 
    magma_minproductDoubleComplex_ptr dskp, 
    magma_minproductDoubleComplex_ptr dr,
    magma_minproductDoubleComplex_ptr dv, 
    magma_minproductDoubleComplex_ptr ds,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zbicgmerge3(
    int n, 
    magma_minproductDoubleComplex_ptr dskp, 
    magma_minproductDoubleComplex_ptr dp,
    magma_minproductDoubleComplex_ptr ds,
    magma_minproductDoubleComplex_ptr dt,
    magma_minproductDoubleComplex_ptr dx, 
    magma_minproductDoubleComplex_ptr dr,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zbicgmerge4(
    int type, 
    magma_minproductDoubleComplex_ptr dskp,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zcgmerge_spmv1( 
    magma_minproduct_z_matrix A,
    magma_minproductDoubleComplex_ptr d1,
    magma_minproductDoubleComplex_ptr d2,
    magma_minproductDoubleComplex_ptr dd,
    magma_minproductDoubleComplex_ptr dz,
    magma_minproductDoubleComplex_ptr dskp,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zcgmerge_xrbeta( 
    int n,
    magma_minproductDoubleComplex_ptr d1,
    magma_minproductDoubleComplex_ptr d2,
    magma_minproductDoubleComplex_ptr dx,
    magma_minproductDoubleComplex_ptr dr,
    magma_minproductDoubleComplex_ptr dd,
    magma_minproductDoubleComplex_ptr dz, 
    magma_minproductDoubleComplex_ptr dskp,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zmdotc(
    magma_minproduct_int_t n, 
    magma_minproduct_int_t k, 
    magma_minproductDoubleComplex_ptr dv, 
    magma_minproductDoubleComplex_ptr dr,
    magma_minproductDoubleComplex_ptr dd1,
    magma_minproductDoubleComplex_ptr dd2,
    magma_minproductDoubleComplex_ptr dskp,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zgemvmdot(
    int n, 
    int k, 
    magma_minproductDoubleComplex_ptr dv, 
    magma_minproductDoubleComplex_ptr dr,
    magma_minproductDoubleComplex_ptr dd1,
    magma_minproductDoubleComplex_ptr dd2,
    magma_minproductDoubleComplex_ptr dskp,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zbicgmerge_spmv1( 
    magma_minproduct_z_matrix A,
    magma_minproductDoubleComplex_ptr dd1,
    magma_minproductDoubleComplex_ptr dd2,
    magma_minproductDoubleComplex_ptr dp,
    magma_minproductDoubleComplex_ptr dr,
    magma_minproductDoubleComplex_ptr dv,
    magma_minproductDoubleComplex_ptr dskp,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zbicgmerge_spmv2( 
    magma_minproduct_z_matrix A,
    magma_minproductDoubleComplex_ptr dd1,
    magma_minproductDoubleComplex_ptr dd2,
    magma_minproductDoubleComplex_ptr ds,
    magma_minproductDoubleComplex_ptr dt,
    magma_minproductDoubleComplex_ptr dskp,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zbicgmerge_xrbeta( 
    int n,
    magma_minproductDoubleComplex_ptr dd1,
    magma_minproductDoubleComplex_ptr dd2,
    magma_minproductDoubleComplex_ptr drr,
    magma_minproductDoubleComplex_ptr dr,
    magma_minproductDoubleComplex_ptr dp,
    magma_minproductDoubleComplex_ptr ds,
    magma_minproductDoubleComplex_ptr dt,
    magma_minproductDoubleComplex_ptr dx, 
    magma_minproductDoubleComplex_ptr dskp,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zbcsrswp(
    magma_minproduct_int_t n,
    magma_minproduct_int_t size_b, 
    magma_minproduct_int_t *ipiv,
    magma_minproductDoubleComplex_ptr dx,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zbcsrtrsv(
    magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t r_blocks,
    magma_minproduct_int_t c_blocks,
    magma_minproduct_int_t size_b, 
    magma_minproductDoubleComplex_ptr dA,
    magma_minproduct_index_t *blockinfo, 
    magma_minproductDoubleComplex_ptr dx,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zbcsrvalcpy(
    magma_minproduct_int_t size_b, 
    magma_minproduct_int_t num_blocks, 
    magma_minproduct_int_t num_zero_blocks, 
    magma_minproductDoubleComplex_ptr *dAval, 
    magma_minproductDoubleComplex_ptr *dBval,
    magma_minproductDoubleComplex_ptr *dBval2,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zbcsrluegemm(
    magma_minproduct_int_t size_b, 
    magma_minproduct_int_t num_block_rows,
    magma_minproduct_int_t kblocks,
    magma_minproductDoubleComplex_ptr *dA, 
    magma_minproductDoubleComplex_ptr *dB, 
    magma_minproductDoubleComplex_ptr *dC,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zbcsrlupivloc(
    magma_minproduct_int_t size_b, 
    magma_minproduct_int_t kblocks,
    magma_minproductDoubleComplex_ptr *dA, 
    magma_minproduct_int_t *ipiv,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_zbcsrblockinfo5(
    magma_minproduct_int_t lustep,
    magma_minproduct_int_t num_blocks, 
    magma_minproduct_int_t c_blocks, 
    magma_minproduct_int_t size_b,
    magma_minproduct_index_t *blockinfo,
    magma_minproductDoubleComplex_ptr dval,
    magma_minproductDoubleComplex_ptr *AII,
    magma_minproduct_queue_t queue );

 
#ifdef __cplusplus
}
#endif

#undef PRECISION_z
#endif /* MAGMA_minproductSPARSE_Z_H */
