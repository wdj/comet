/*
 -- MAGMA_minproduct (version 1.1) --
 Univ. of Tennessee, Knoxville
 Univ. of California, Berkeley
 Univ. of Colorado, Denver
 @date

 @generated from magma_minproductsparse_z.h normal z -> c, Mon May  4 17:11:25 2015
 @author Hartwig Anzt
*/

#ifndef MAGMA_minproductSPARSE_C_H
#define MAGMA_minproductSPARSE_C_H

#include "magma_minproduct_types.h"
#include "magma_minproductsparse_types.h"

#define PRECISION_c


#ifdef __cplusplus
extern "C" {
#endif


/* ////////////////////////////////////////////////////////////////////////////
 -- For backwards compatability, map old (1.6.1) to new (1.6.2) function names
*/

#define magma_minproduct_c_mtranspose  magma_minproduct_cmtranspose
#define magma_minproduct_c_mtransfer   magma_minproduct_cmtransfer
#define magma_minproduct_c_vtransfer   magma_minproduct_cmtransfer
#define magma_minproduct_c_mconvert    magma_minproduct_cmconvert
#define magma_minproduct_c_vinit       magma_minproduct_cvinit
#define magma_minproduct_c_vvisu       magma_minproduct_cprint_vector
#define magma_minproduct_c_vread       magma_minproduct_cvread
#define magma_minproduct_c_vspread     magma_minproduct_cvspread
#define magma_minproduct_c_mvisu       magma_minproduct_cprint_matrix
#define magma_minproduct_c_mfree       magma_minproduct_cmfree
#define magma_minproduct_c_vfree       magma_minproduct_cmfree
#define write_c_csr_mtx     magma_minproduct_cwrite_csr_mtx
#define write_c_csrtomtx    magma_minproduct_cwrite_csrtomtx
#define print_c_csr         magma_minproduct_cprint_csr_mtx


/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_minproduct_SPARSE Auxiliary functions
*/


magma_minproduct_int_t
magma_minproduct_cparse_opts( 
    int argc, 
    char** argv, 
    magma_minproduct_copts *opts, 
    int *matrices, 
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
read_c_csr_from_binary( 
    magma_minproduct_int_t* n_row, 
    magma_minproduct_int_t* n_col, 
    magma_minproduct_int_t* nnz, 
    magma_minproductFloatComplex **val, 
    magma_minproduct_index_t **row, 
    magma_minproduct_index_t **col,
    const char * filename,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
read_c_csr_from_mtx( 
    magma_minproduct_storage_t *type, 
    magma_minproduct_location_t *location,
    magma_minproduct_int_t* n_row, 
    magma_minproduct_int_t* n_col, 
    magma_minproduct_int_t* nnz, 
    magma_minproductFloatComplex **val, 
    magma_minproduct_index_t **row, 
    magma_minproduct_index_t **col, 
    const char *filename,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_c_csr_mtx( 
    magma_minproduct_c_matrix *A, 
    const char *filename,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_ccsrset( 
    magma_minproduct_int_t m, 
    magma_minproduct_int_t n, 
    magma_minproduct_index_t *row, 
    magma_minproduct_index_t *col, 
    magma_minproductFloatComplex *val,
    magma_minproduct_c_matrix *A,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_ccsrget( 
    magma_minproduct_c_matrix A,
    magma_minproduct_int_t *m, 
    magma_minproduct_int_t *n, 
    magma_minproduct_index_t **row, 
    magma_minproduct_index_t **col, 
    magma_minproductFloatComplex **val,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t 
magma_minproduct_cvset( 
    magma_minproduct_int_t m, 
    magma_minproduct_int_t n, 
    magma_minproductFloatComplex *val,
    magma_minproduct_c_matrix *v,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_cvget( 
    magma_minproduct_c_matrix v,
    magma_minproduct_int_t *m, 
    magma_minproduct_int_t *n, 
    magma_minproductFloatComplex **val,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_cvset_dev( 
    magma_minproduct_int_t m, 
    magma_minproduct_int_t n, 
    magma_minproductFloatComplex_ptr val,
    magma_minproduct_c_matrix *v,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_cvget_dev( 
    magma_minproduct_c_matrix v,
    magma_minproduct_int_t *m, 
    magma_minproduct_int_t *n, 
    magma_minproductFloatComplex_ptr *val,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t 
magma_minproduct_c_csr_mtxsymm( 
    magma_minproduct_c_matrix *A, 
    const char *filename,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_c_csr_compressor( 
    magma_minproductFloatComplex ** val, 
    magma_minproduct_index_t ** row, 
    magma_minproduct_index_t ** col, 
    magma_minproductFloatComplex ** valn, 
    magma_minproduct_index_t ** rown, 
    magma_minproduct_index_t ** coln, 
    magma_minproduct_int_t *n,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_cmcsrcompressor( 
    magma_minproduct_c_matrix *A,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_cmcsrcompressor_gpu( 
    magma_minproduct_c_matrix *A,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_cvtranspose( 
    magma_minproduct_c_matrix x,
    magma_minproduct_c_matrix *y,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_c_cucsrtranspose( 
    magma_minproduct_c_matrix A, 
    magma_minproduct_c_matrix *B,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
c_transpose_csr( 
    magma_minproduct_int_t n_rows, 
    magma_minproduct_int_t n_cols, 
    magma_minproduct_int_t nnz,
    magma_minproductFloatComplex *val, 
    magma_minproduct_index_t *row, 
    magma_minproduct_index_t *col, 
    magma_minproduct_int_t *new_n_rows, 
    magma_minproduct_int_t *new_n_cols, 
    magma_minproduct_int_t *new_nnz, 
    magma_minproductFloatComplex **new_val, 
    magma_minproduct_index_t **new_row, 
    magma_minproduct_index_t **new_col,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_ccsrsplit( 
    magma_minproduct_int_t bsize,
    magma_minproduct_c_matrix A,
    magma_minproduct_c_matrix *D,
    magma_minproduct_c_matrix *R,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_cmscale( 
    magma_minproduct_c_matrix *A, 
    magma_minproduct_scale_t scaling,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_cmdiff( 
    magma_minproduct_c_matrix A, 
    magma_minproduct_c_matrix B, 
 real_Double_t *res,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_cmdiagadd( 
    magma_minproduct_c_matrix *A, 
    magma_minproductFloatComplex add,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_cmsort( 
    magma_minproduct_c_matrix *A,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_cindexsort(
    magma_minproduct_index_t *x, 
    magma_minproduct_int_t first,
    magma_minproduct_int_t last,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_cdomainoverlap(
    magma_minproduct_index_t num_rows,
    magma_minproduct_index_t *num_indices,
    magma_minproduct_index_t *rowptr,
    magma_minproduct_index_t *colidx,
    magma_minproduct_index_t *x,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_csymbilu( 
    magma_minproduct_c_matrix *A, 
    magma_minproduct_int_t levels,
    magma_minproduct_c_matrix *L,
    magma_minproduct_c_matrix *U,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t 
magma_minproduct_cwrite_csr_mtx( 
    magma_minproduct_c_matrix A,
    magma_minproduct_order_t MajorType,
 const char *filename,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_cwrite_csrtomtx( 
    magma_minproduct_c_matrix A,
    const char *filename,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_cprint_csr( 
    magma_minproduct_int_t n_row, 
    magma_minproduct_int_t n_col, 
    magma_minproduct_int_t nnz, 
    magma_minproductFloatComplex **val, 
    magma_minproduct_index_t **row, 
    magma_minproduct_index_t **col,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_cprint_csr_mtx( 
    magma_minproduct_int_t n_row, 
    magma_minproduct_int_t n_col, 
    magma_minproduct_int_t nnz, 
    magma_minproductFloatComplex **val, 
    magma_minproduct_index_t **row, 
    magma_minproduct_index_t **col, 
    magma_minproduct_order_t MajorType,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t 
magma_minproduct_cmtranspose(
    magma_minproduct_c_matrix A, 
    magma_minproduct_c_matrix *B,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t 
magma_minproduct_cmtransfer(
    magma_minproduct_c_matrix A, 
    magma_minproduct_c_matrix *B, 
    magma_minproduct_location_t src, 
    magma_minproduct_location_t dst,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_cmconvert(
    magma_minproduct_c_matrix A, 
    magma_minproduct_c_matrix *B, 
    magma_minproduct_storage_t old_format, 
    magma_minproduct_storage_t new_format,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t
magma_minproduct_cvinit(
    magma_minproduct_c_matrix *x, 
    magma_minproduct_location_t memory_location,
    magma_minproduct_int_t num_rows, 
    magma_minproduct_int_t num_cols,
    magma_minproductFloatComplex values,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_cprint_vector(
    magma_minproduct_c_matrix x, 
    magma_minproduct_int_t offset, 
    magma_minproduct_int_t displaylength,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_cvread(
    magma_minproduct_c_matrix *x, 
    magma_minproduct_int_t length,
    char * filename,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_cvspread(
    magma_minproduct_c_matrix *x, 
    const char * filename,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_cprint_matrix(
    magma_minproduct_c_matrix A,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_cdiameter(
    magma_minproduct_c_matrix *A,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_crowentries(
    magma_minproduct_c_matrix *A,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_cmfree(
    magma_minproduct_c_matrix *A,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_cresidual(
    magma_minproduct_c_matrix A, 
    magma_minproduct_c_matrix b, 
    magma_minproduct_c_matrix x, 
    float *res,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_cresidualvec(
    magma_minproduct_c_matrix A,
    magma_minproduct_c_matrix b,
    magma_minproduct_c_matrix x,
    magma_minproduct_c_matrix *r,
    float *res,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_cmgenerator(
    magma_minproduct_int_t n,
    magma_minproduct_int_t offdiags,
    magma_minproduct_index_t *diag_offset,
    magma_minproductFloatComplex *diag_vals,
    magma_minproduct_c_matrix *A,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_cm_27stencil(
    magma_minproduct_int_t n,
    magma_minproduct_c_matrix *A,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_cm_5stencil(
    magma_minproduct_int_t n,
    magma_minproduct_c_matrix *A,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_csolverinfo(
    magma_minproduct_c_solver_par *solver_par, 
    magma_minproduct_c_preconditioner *precond_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_csolverinfo_init(
    magma_minproduct_c_solver_par *solver_par, 
    magma_minproduct_c_preconditioner *precond,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_ceigensolverinfo_init(
    magma_minproduct_c_solver_par *solver_par,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t
magma_minproduct_csolverinfo_free(
    magma_minproduct_c_solver_par *solver_par, 
    magma_minproduct_c_preconditioner *precond,
    magma_minproduct_queue_t queue );



/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_minproduct_SPARSE iterative incomplete factorizations
*/


magma_minproduct_int_t
magma_minproduct_citerilusetup( 
    magma_minproduct_c_matrix A, 
    magma_minproduct_c_matrix b,                                 
    magma_minproduct_c_preconditioner *precond,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_citericsetup( 
    magma_minproduct_c_matrix A, 
    magma_minproduct_c_matrix b, 
    magma_minproduct_c_preconditioner *precond,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_citericupdate( 
    magma_minproduct_c_matrix A, 
    magma_minproduct_c_preconditioner *precond, 
    magma_minproduct_int_t updates,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_capplyiteric_l( 
    magma_minproduct_c_matrix b, 
    magma_minproduct_c_matrix *x, 
    magma_minproduct_c_preconditioner *precond,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_capplyiteric_r( 
    magma_minproduct_c_matrix b, 
    magma_minproduct_c_matrix *x, 
    magma_minproduct_c_preconditioner *precond,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t
magma_minproduct_citerilu_csr( 
    magma_minproduct_c_matrix A,
    magma_minproduct_c_matrix L,
    magma_minproduct_c_matrix U,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_citeric_csr( 
    magma_minproduct_c_matrix A,
    magma_minproduct_c_matrix A_CSR,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_cfrobenius( 
    magma_minproduct_c_matrix A, 
    magma_minproduct_c_matrix B, 
    real_Double_t *res,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_cnonlinres(   
    magma_minproduct_c_matrix A, 
    magma_minproduct_c_matrix L,
    magma_minproduct_c_matrix U, 
    magma_minproduct_c_matrix *LU, 
    real_Double_t *res,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_cilures(   
    magma_minproduct_c_matrix A, 
    magma_minproduct_c_matrix L,
    magma_minproduct_c_matrix U, 
    magma_minproduct_c_matrix *LU, 
    real_Double_t *res,
    real_Double_t *nonlinres,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_cicres(       
    magma_minproduct_c_matrix A, 
    magma_minproduct_c_matrix C,
    magma_minproduct_c_matrix CT, 
    magma_minproduct_c_matrix *LU, 
    real_Double_t *res,
    real_Double_t *nonlinres,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_cinitguess( 
    magma_minproduct_c_matrix A, 
    magma_minproduct_c_matrix *L, 
    magma_minproduct_c_matrix *U,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_cinitrecursiveLU( 
    magma_minproduct_c_matrix A, 
    magma_minproduct_c_matrix *B,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_cmLdiagadd( 
    magma_minproduct_c_matrix *L,
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
magma_minproduct_ccg(
    magma_minproduct_c_matrix A, magma_minproduct_c_matrix b, 
    magma_minproduct_c_matrix *x, magma_minproduct_c_solver_par *solver_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_ccg_res(
    magma_minproduct_c_matrix A, magma_minproduct_c_matrix b, 
    magma_minproduct_c_matrix *x, magma_minproduct_c_solver_par *solver_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_ccg_merge(
    magma_minproduct_c_matrix A, magma_minproduct_c_matrix b, 
    magma_minproduct_c_matrix *x, magma_minproduct_c_solver_par *solver_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_cgmres(
    magma_minproduct_c_matrix A, magma_minproduct_c_matrix b, 
    magma_minproduct_c_matrix *x, magma_minproduct_c_solver_par *solver_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_cbicgstab(
    magma_minproduct_c_matrix A, magma_minproduct_c_matrix b, magma_minproduct_c_matrix *x, 
    magma_minproduct_c_solver_par *solver_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_cbicgstab_merge(
    magma_minproduct_c_matrix A, magma_minproduct_c_matrix b, 
    magma_minproduct_c_matrix *x, magma_minproduct_c_solver_par *solver_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_cbicgstab_merge2(
    magma_minproduct_c_matrix A, magma_minproduct_c_matrix b, 
    magma_minproduct_c_matrix *x, magma_minproduct_c_solver_par *solver_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_cpcg(
    magma_minproduct_c_matrix A, magma_minproduct_c_matrix b, 
    magma_minproduct_c_matrix *x, magma_minproduct_c_solver_par *solver_par, 
    magma_minproduct_c_preconditioner *precond_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_cbpcg(
    magma_minproduct_c_matrix A, magma_minproduct_c_matrix b, 
    magma_minproduct_c_matrix *x, magma_minproduct_c_solver_par *solver_par, 
    magma_minproduct_c_preconditioner *precond_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_cpbicgstab(
    magma_minproduct_c_matrix A, magma_minproduct_c_matrix b, 
    magma_minproduct_c_matrix *x, magma_minproduct_c_solver_par *solver_par, 
    magma_minproduct_c_preconditioner *precond_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_cpgmres(
    magma_minproduct_c_matrix A, magma_minproduct_c_matrix b, 
    magma_minproduct_c_matrix *x, magma_minproduct_c_solver_par *solver_par, 
    magma_minproduct_c_preconditioner *precond_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_cfgmres(
    magma_minproduct_c_matrix A, magma_minproduct_c_matrix b, 
    magma_minproduct_c_matrix *x, magma_minproduct_c_solver_par *solver_par, 
    magma_minproduct_c_preconditioner *precond_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_cjacobi(
    magma_minproduct_c_matrix A, magma_minproduct_c_matrix b, 
    magma_minproduct_c_matrix *x, magma_minproduct_c_solver_par *solver_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_cjacobidomainoverlap(
    magma_minproduct_c_matrix A, 
    magma_minproduct_c_matrix b, 
    magma_minproduct_c_matrix *x,  
    magma_minproduct_c_solver_par *solver_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_cbaiter(
    magma_minproduct_c_matrix A, magma_minproduct_c_matrix b, 
    magma_minproduct_c_matrix *x, magma_minproduct_c_solver_par *solver_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_citerref(
    magma_minproduct_c_matrix A, magma_minproduct_c_matrix b, 
    magma_minproduct_c_matrix *x, magma_minproduct_c_solver_par *solver_par, 
    magma_minproduct_c_preconditioner *precond_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_cilu(
    magma_minproduct_c_matrix A, magma_minproduct_c_matrix b, 
    magma_minproduct_c_matrix *x, magma_minproduct_c_solver_par *solver_par, 
    magma_minproduct_int_t *ipiv,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_cbcsrlu(
    magma_minproduct_c_matrix A, magma_minproduct_c_matrix b, 
    magma_minproduct_c_matrix *x, magma_minproduct_c_solver_par *solver_par,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t
magma_minproduct_cbcsrlutrf(
    magma_minproduct_c_matrix A, 
    magma_minproduct_c_matrix *M,
    magma_minproduct_int_t *ipiv, 
    magma_minproduct_int_t version,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_cbcsrlusv(
    magma_minproduct_c_matrix A, magma_minproduct_c_matrix b, 
    magma_minproduct_c_matrix *x, magma_minproduct_c_solver_par *solver_par, 
    magma_minproduct_int_t *ipiv,
    magma_minproduct_queue_t queue );



magma_minproduct_int_t
magma_minproduct_cilucg(
    magma_minproduct_c_matrix A, magma_minproduct_c_matrix b, 
    magma_minproduct_c_matrix *x, magma_minproduct_c_solver_par *solver_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_cilugmres(
    magma_minproduct_c_matrix A, magma_minproduct_c_matrix b, 
    magma_minproduct_c_matrix *x, magma_minproduct_c_solver_par *solver_par,
    magma_minproduct_queue_t queue ); 


magma_minproduct_int_t
magma_minproduct_clobpcg_shift(
    magma_minproduct_int_t num_rows,
    magma_minproduct_int_t num_vecs, 
    magma_minproduct_int_t shift,
    magma_minproductFloatComplex_ptr x,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_clobpcg_res(
    magma_minproduct_int_t num_rows,
    magma_minproduct_int_t num_vecs, 
    float *evalues, 
    magma_minproductFloatComplex_ptr X,
    magma_minproductFloatComplex_ptr R, 
    float *res,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_clobpcg_maxpy(
    magma_minproduct_int_t num_rows,
    magma_minproduct_int_t num_vecs, 
    magma_minproductFloatComplex_ptr X,
    magma_minproductFloatComplex_ptr Y,
    magma_minproduct_queue_t queue );


/*/////////////////////////////////////////////////////////////////////////////
 -- MAGMA_minproduct_SPARSE eigensolvers (Data on GPU)
*/
magma_minproduct_int_t
magma_minproduct_clobpcg(
    magma_minproduct_c_matrix A, 
    magma_minproduct_c_solver_par *solver_par,
    magma_minproduct_c_preconditioner *precond_par, 
    magma_minproduct_queue_t queue );




/*/////////////////////////////////////////////////////////////////////////////
 -- MAGMA_minproduct_SPARSE preconditioners (Data on GPU)
*/

magma_minproduct_int_t
magma_minproduct_cjacobisetup(
    magma_minproduct_c_matrix A, 
    magma_minproduct_c_matrix b, 
    magma_minproduct_c_matrix *M, 
    magma_minproduct_c_matrix *c,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_cjacobisetup_matrix(
    magma_minproduct_c_matrix A, 
    magma_minproduct_c_matrix *M, 
    magma_minproduct_c_matrix *d,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_cjacobisetup_vector(
    magma_minproduct_c_matrix b, 
    magma_minproduct_c_matrix d, 
    magma_minproduct_c_matrix *c,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_cjacobiiter(
    magma_minproduct_c_matrix M, 
    magma_minproduct_c_matrix c, 
    magma_minproduct_c_matrix *x, 
    magma_minproduct_c_solver_par *solver_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_cjacobiiter_precond( 
    magma_minproduct_c_matrix M, 
    magma_minproduct_c_matrix *x, 
    magma_minproduct_c_solver_par *solver_par, 
    magma_minproduct_c_preconditioner *precond,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_cjacobiiter_sys(
    magma_minproduct_c_matrix A, 
    magma_minproduct_c_matrix b, 
    magma_minproduct_c_matrix d, 
    magma_minproduct_c_matrix t, 
    magma_minproduct_c_matrix *x,  
    magma_minproduct_c_solver_par *solver_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_cpastixsetup(
    magma_minproduct_c_matrix A, magma_minproduct_c_matrix b,
    magma_minproduct_c_preconditioner *precond,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t
magma_minproduct_capplypastix(
    magma_minproduct_c_matrix b, magma_minproduct_c_matrix *x, 
    magma_minproduct_c_preconditioner *precond,
    magma_minproduct_queue_t queue );


// custom preconditioner
magma_minproduct_int_t
magma_minproduct_capplycustomprecond_l(
    magma_minproduct_c_matrix b, 
    magma_minproduct_c_matrix *x, 
    magma_minproduct_c_preconditioner *precond,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_capplycustomprecond_r(
    magma_minproduct_c_matrix b, 
    magma_minproduct_c_matrix *x, 
    magma_minproduct_c_preconditioner *precond,
    magma_minproduct_queue_t queue );


// CUSPARSE preconditioner

magma_minproduct_int_t
magma_minproduct_ccuilusetup(
    magma_minproduct_c_matrix A, magma_minproduct_c_preconditioner *precond,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_capplycuilu_l(
    magma_minproduct_c_matrix b, magma_minproduct_c_matrix *x, 
    magma_minproduct_c_preconditioner *precond,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_capplycuilu_r(
    magma_minproduct_c_matrix b, magma_minproduct_c_matrix *x, 
    magma_minproduct_c_preconditioner *precond,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t
magma_minproduct_ccuiccsetup(
    magma_minproduct_c_matrix A, magma_minproduct_c_preconditioner *precond,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_capplycuicc_l(
    magma_minproduct_c_matrix b, magma_minproduct_c_matrix *x, 
    magma_minproduct_c_preconditioner *precond,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_capplycuicc_r(
    magma_minproduct_c_matrix b, magma_minproduct_c_matrix *x, 
    magma_minproduct_c_preconditioner *precond,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t
magma_minproduct_ccumilusetup(
    magma_minproduct_c_matrix A, 
    magma_minproduct_c_preconditioner *precond,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_ccumilugeneratesolverinfo(
    magma_minproduct_c_preconditioner *precond,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_capplycumilu_l(
    magma_minproduct_c_matrix b, 
    magma_minproduct_c_matrix *x, 
    magma_minproduct_c_preconditioner *precond,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_capplycumilu_r(
    magma_minproduct_c_matrix b, 
    magma_minproduct_c_matrix *x, 
    magma_minproduct_c_preconditioner *precond,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t
magma_minproduct_ccumiccsetup(
    magma_minproduct_c_matrix A, 
    magma_minproduct_c_preconditioner *precond,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_ccumicgeneratesolverinfo(
    magma_minproduct_c_preconditioner *precond,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_capplycumicc_l(
    magma_minproduct_c_matrix b, 
    magma_minproduct_c_matrix *x, 
    magma_minproduct_c_preconditioner *precond,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_capplycumicc_r(
    magma_minproduct_c_matrix b, 
    magma_minproduct_c_matrix *x, 
    magma_minproduct_c_preconditioner *precond,
    magma_minproduct_queue_t queue );


// block-asynchronous iteration

magma_minproduct_int_t
magma_minproduct_cbajac_csr(
    magma_minproduct_int_t localiters,
    magma_minproduct_c_matrix D,
    magma_minproduct_c_matrix R,
    magma_minproduct_c_matrix b,
    magma_minproduct_c_matrix *x,
    magma_minproduct_queue_t queue );

/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_minproduct_SPARSE utility function definitions
*/

magma_minproduct_int_t
magma_minproduct_c_spmv(
    magma_minproductFloatComplex alpha, 
    magma_minproduct_c_matrix A, 
    magma_minproduct_c_matrix x, 
    magma_minproductFloatComplex beta, 
    magma_minproduct_c_matrix y,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_ccustomspmv(
    magma_minproductFloatComplex alpha, 
    magma_minproduct_c_matrix x, 
    magma_minproductFloatComplex beta, 
    magma_minproduct_c_matrix y,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_c_spmv_shift(
    magma_minproductFloatComplex alpha, 
    magma_minproduct_c_matrix A, 
    magma_minproductFloatComplex lambda,
    magma_minproduct_c_matrix x, 
    magma_minproductFloatComplex beta, 
    magma_minproduct_int_t offset, 
    magma_minproduct_int_t blocksize,
    magma_minproductIndex_ptr dadd_vecs, 
    magma_minproduct_c_matrix y,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_ccuspmm(
    magma_minproduct_c_matrix A, 
    magma_minproduct_c_matrix B, 
    magma_minproduct_c_matrix *AB,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_c_spmm(
    magma_minproductFloatComplex alpha, 
    magma_minproduct_c_matrix A,
    magma_minproduct_c_matrix B,
    magma_minproduct_c_matrix *C,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_csymbilu( 
    magma_minproduct_c_matrix *A, 
    magma_minproduct_int_t levels, 
    magma_minproduct_c_matrix *L, 
    magma_minproduct_c_matrix *U,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_ccuspaxpy(
    magma_minproductFloatComplex_ptr alpha, magma_minproduct_c_matrix A, 
    magma_minproductFloatComplex_ptr beta, magma_minproduct_c_matrix B, 
    magma_minproduct_c_matrix *AB,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_c_precond(
    magma_minproduct_c_matrix A, 
    magma_minproduct_c_matrix b, magma_minproduct_c_matrix *x,
    magma_minproduct_c_preconditioner *precond,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_c_solver(
    magma_minproduct_c_matrix A, magma_minproduct_c_matrix b, 
    magma_minproduct_c_matrix *x, magma_minproduct_copts *zopts,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_c_precondsetup(
    magma_minproduct_c_matrix A, magma_minproduct_c_matrix b, 
    magma_minproduct_c_preconditioner *precond,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_c_applyprecond(
    magma_minproduct_c_matrix A, magma_minproduct_c_matrix b, 
    magma_minproduct_c_matrix *x, magma_minproduct_c_preconditioner *precond,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t
magma_minproduct_c_applyprecond_left(
    magma_minproduct_c_matrix A, magma_minproduct_c_matrix b, 
    magma_minproduct_c_matrix *x, magma_minproduct_c_preconditioner *precond,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t
magma_minproduct_c_applyprecond_right(
    magma_minproduct_c_matrix A, magma_minproduct_c_matrix b, 
    magma_minproduct_c_matrix *x, magma_minproduct_c_preconditioner *precond,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_c_initP2P(
    magma_minproduct_int_t *bandwidth_benchmark,
    magma_minproduct_int_t *num_gpus,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_ccompact(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda,
    float *dnorms, float tol, 
    magma_minproduct_int_t *activeMask, magma_minproduct_int_t *cBlockSize,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_ccompactActive(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda, 
    magma_minproduct_int_t *active,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_cmlumerge(    
    magma_minproduct_c_matrix L, 
    magma_minproduct_c_matrix U,
    magma_minproduct_c_matrix *A, 
    magma_minproduct_queue_t queue );


/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_minproduct_SPARSE BLAS function definitions
*/
magma_minproduct_int_t 
magma_minproduct_cgecsrmv(
    magma_minproduct_trans_t transA,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex alpha,
    magma_minproductFloatComplex_ptr dval,
    magma_minproductIndex_ptr drowptr,
    magma_minproductIndex_ptr dcolind,
    magma_minproductFloatComplex_ptr dx,
    magma_minproductFloatComplex beta,
    magma_minproductFloatComplex_ptr dy,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_cgecsrmv_shift(
    magma_minproduct_trans_t transA,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex alpha,
    magma_minproductFloatComplex lambda,
    magma_minproductFloatComplex_ptr dval,
    magma_minproductIndex_ptr drowptr,
    magma_minproductIndex_ptr dcolind,
    magma_minproductFloatComplex_ptr dx,
    magma_minproductFloatComplex beta,
    int offset,
    int blocksize,
    magma_minproductIndex_ptr dadd_rows,
    magma_minproductFloatComplex_ptr dy,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_cmgecsrmv(
    magma_minproduct_trans_t transA,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproduct_int_t num_vecs,
    magma_minproductFloatComplex alpha,
    magma_minproductFloatComplex_ptr dval,
    magma_minproductIndex_ptr drowptr,
    magma_minproductIndex_ptr dcolind,
    magma_minproductFloatComplex_ptr dx,
    magma_minproductFloatComplex beta,
    magma_minproductFloatComplex_ptr dy,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_cgeellmv(
    magma_minproduct_trans_t transA,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproduct_int_t nnz_per_row,
    magma_minproductFloatComplex alpha,
    magma_minproductFloatComplex_ptr dval,
    magma_minproductIndex_ptr dcolind,
    magma_minproductFloatComplex_ptr dx,
    magma_minproductFloatComplex beta,
    magma_minproductFloatComplex_ptr dy,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_cgeellmv_shift(
    magma_minproduct_trans_t transA,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproduct_int_t nnz_per_row,
    magma_minproductFloatComplex alpha,
    magma_minproductFloatComplex lambda,
    magma_minproductFloatComplex_ptr dval,
    magma_minproductIndex_ptr dcolind,
    magma_minproductFloatComplex_ptr dx,
    magma_minproductFloatComplex beta,
    int offset,
    int blocksize,
    magma_minproductIndex_ptr dadd_rows,
    magma_minproductFloatComplex_ptr dy,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t 
magma_minproduct_cmgeellmv(
    magma_minproduct_trans_t transA,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproduct_int_t num_vecs,
    magma_minproduct_int_t nnz_per_row,
    magma_minproductFloatComplex alpha,
    magma_minproductFloatComplex_ptr dval,
    magma_minproductIndex_ptr dcolind,
    magma_minproductFloatComplex_ptr dx,
    magma_minproductFloatComplex beta,
    magma_minproductFloatComplex_ptr dy,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t 
magma_minproduct_cgeelltmv(
    magma_minproduct_trans_t transA,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproduct_int_t nnz_per_row,
    magma_minproductFloatComplex alpha,
    magma_minproductFloatComplex_ptr dval,
    magma_minproductIndex_ptr dcolind,
    magma_minproductFloatComplex_ptr dx,
    magma_minproductFloatComplex beta,
    magma_minproductFloatComplex_ptr dy,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_cgeelltmv_shift(
    magma_minproduct_trans_t transA,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproduct_int_t nnz_per_row,
    magma_minproductFloatComplex alpha,
    magma_minproductFloatComplex lambda,
    magma_minproductFloatComplex_ptr dval,
    magma_minproductIndex_ptr dcolind,
    magma_minproductFloatComplex_ptr dx,
    magma_minproductFloatComplex beta,
    int offset,
    int blocksize,
    magma_minproductIndex_ptr dadd_rows,
    magma_minproductFloatComplex_ptr dy,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t 
magma_minproduct_cmgeelltmv(
    magma_minproduct_trans_t transA,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproduct_int_t num_vecs,
    magma_minproduct_int_t nnz_per_row,
    magma_minproductFloatComplex alpha,
    magma_minproductFloatComplex_ptr dval,
    magma_minproductIndex_ptr dcolind,
    magma_minproductFloatComplex_ptr dx,
    magma_minproductFloatComplex beta,
    magma_minproductFloatComplex_ptr dy,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_cgeellrtmv(
    magma_minproduct_trans_t transA,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproduct_int_t nnz_per_row,
    magma_minproductFloatComplex alpha,
    magma_minproductFloatComplex_ptr dval,
    magma_minproductIndex_ptr dcolind,
    magma_minproductIndex_ptr drowlength,
    magma_minproductFloatComplex_ptr dx,
    magma_minproductFloatComplex beta,
    magma_minproductFloatComplex_ptr dy,
    magma_minproduct_int_t num_threads,
    magma_minproduct_int_t threads_per_row,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_cgesellcmv(
    magma_minproduct_trans_t transA,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproduct_int_t blocksize,
    magma_minproduct_int_t slices,
    magma_minproduct_int_t alignment,
    magma_minproductFloatComplex alpha,
    magma_minproductFloatComplex_ptr dval,
    magma_minproductIndex_ptr dcolind,
    magma_minproductIndex_ptr drowptr,
    magma_minproductFloatComplex_ptr dx,
    magma_minproductFloatComplex beta,
    magma_minproductFloatComplex_ptr dy,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_cgesellpmv(
    magma_minproduct_trans_t transA,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproduct_int_t blocksize,
    magma_minproduct_int_t slices,
    magma_minproduct_int_t alignment,
    magma_minproductFloatComplex alpha,
    magma_minproductFloatComplex_ptr dval,
    magma_minproductIndex_ptr dcolind,
    magma_minproductIndex_ptr drowptr,
    magma_minproductFloatComplex_ptr dx,
    magma_minproductFloatComplex beta,
    magma_minproductFloatComplex_ptr dy,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_cmgesellpmv(
    magma_minproduct_trans_t transA,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproduct_int_t num_vecs,
    magma_minproduct_int_t blocksize,
    magma_minproduct_int_t slices,
    magma_minproduct_int_t alignment,
    magma_minproductFloatComplex alpha,
    magma_minproductFloatComplex_ptr dval,
    magma_minproductIndex_ptr dcolind,
    magma_minproductIndex_ptr drowptr,
    magma_minproductFloatComplex_ptr dx,
    magma_minproductFloatComplex beta,
    magma_minproductFloatComplex_ptr dy,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_cmgesellpmv_blocked(
    magma_minproduct_trans_t transA,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproduct_int_t num_vecs,
    magma_minproduct_int_t blocksize,
    magma_minproduct_int_t slices,
    magma_minproduct_int_t alignment,
    magma_minproductFloatComplex alpha,
    magma_minproductFloatComplex_ptr dval,
    magma_minproductIndex_ptr dcolind,
    magma_minproductIndex_ptr drowptr,
    magma_minproductFloatComplex_ptr dx,
    magma_minproductFloatComplex beta,
    magma_minproductFloatComplex_ptr dy,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t
magma_minproduct_cmergedgs(
    magma_minproduct_int_t n, 
    magma_minproduct_int_t ldh,
    magma_minproduct_int_t k, 
    magma_minproductFloatComplex_ptr dv, 
    magma_minproductFloatComplex_ptr dr,
    magma_minproductFloatComplex_ptr dskp,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_ccopyscale(    
    int n, 
    int k,
    magma_minproductFloatComplex_ptr dr, 
    magma_minproductFloatComplex_ptr dv,
    magma_minproductFloatComplex_ptr dskp,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_scnrm2scale(    
    int m, 
    magma_minproductFloatComplex_ptr dr,    
    int lddr, 
    magma_minproductFloatComplex *drnorm,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t
magma_minproduct_cjacobisetup_vector_gpu(
    int num_rows, 
    magma_minproduct_c_matrix b, 
    magma_minproduct_c_matrix d, 
    magma_minproduct_c_matrix c,
    magma_minproduct_c_matrix *x,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t
magma_minproduct_cjacobi_diagscal(    
    int num_rows, 
    magma_minproduct_c_matrix d, 
    magma_minproduct_c_matrix b, 
    magma_minproduct_c_matrix *c,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_cjacobiupdate(
    magma_minproduct_c_matrix t, 
    magma_minproduct_c_matrix b, 
    magma_minproduct_c_matrix d, 
    magma_minproduct_c_matrix *x,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_cjacobispmvupdate(
    magma_minproduct_int_t maxiter,
    magma_minproduct_c_matrix A, 
    magma_minproduct_c_matrix t, 
    magma_minproduct_c_matrix b, 
    magma_minproduct_c_matrix d, 
    magma_minproduct_c_matrix *x,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_cjacobispmvupdate_bw(
    magma_minproduct_int_t maxiter,
    magma_minproduct_c_matrix A, 
    magma_minproduct_c_matrix t, 
    magma_minproduct_c_matrix b, 
    magma_minproduct_c_matrix d, 
    magma_minproduct_c_matrix *x,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_cjacobispmvupdateselect(
    magma_minproduct_int_t maxiter,
    magma_minproduct_int_t num_updates,
    magma_minproduct_index_t *indices,
    magma_minproduct_c_matrix A,
    magma_minproduct_c_matrix t, 
    magma_minproduct_c_matrix b, 
    magma_minproduct_c_matrix d, 
    magma_minproduct_c_matrix tmp, 
    magma_minproduct_c_matrix *x,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_cjacobisetup_diagscal(
    magma_minproduct_c_matrix A, magma_minproduct_c_matrix *d,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t
magma_minproduct_cbicgmerge1(    
    int n, 
    magma_minproductFloatComplex_ptr dskp,
    magma_minproductFloatComplex_ptr dv, 
    magma_minproductFloatComplex_ptr dr, 
    magma_minproductFloatComplex_ptr dp,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t
magma_minproduct_cbicgmerge2(
    int n, 
    magma_minproductFloatComplex_ptr dskp, 
    magma_minproductFloatComplex_ptr dr,
    magma_minproductFloatComplex_ptr dv, 
    magma_minproductFloatComplex_ptr ds,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_cbicgmerge3(
    int n, 
    magma_minproductFloatComplex_ptr dskp, 
    magma_minproductFloatComplex_ptr dp,
    magma_minproductFloatComplex_ptr ds,
    magma_minproductFloatComplex_ptr dt,
    magma_minproductFloatComplex_ptr dx, 
    magma_minproductFloatComplex_ptr dr,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_cbicgmerge4(
    int type, 
    magma_minproductFloatComplex_ptr dskp,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_ccgmerge_spmv1( 
    magma_minproduct_c_matrix A,
    magma_minproductFloatComplex_ptr d1,
    magma_minproductFloatComplex_ptr d2,
    magma_minproductFloatComplex_ptr dd,
    magma_minproductFloatComplex_ptr dz,
    magma_minproductFloatComplex_ptr dskp,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_ccgmerge_xrbeta( 
    int n,
    magma_minproductFloatComplex_ptr d1,
    magma_minproductFloatComplex_ptr d2,
    magma_minproductFloatComplex_ptr dx,
    magma_minproductFloatComplex_ptr dr,
    magma_minproductFloatComplex_ptr dd,
    magma_minproductFloatComplex_ptr dz, 
    magma_minproductFloatComplex_ptr dskp,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_cmdotc(
    magma_minproduct_int_t n, 
    magma_minproduct_int_t k, 
    magma_minproductFloatComplex_ptr dv, 
    magma_minproductFloatComplex_ptr dr,
    magma_minproductFloatComplex_ptr dd1,
    magma_minproductFloatComplex_ptr dd2,
    magma_minproductFloatComplex_ptr dskp,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_cgemvmdot(
    int n, 
    int k, 
    magma_minproductFloatComplex_ptr dv, 
    magma_minproductFloatComplex_ptr dr,
    magma_minproductFloatComplex_ptr dd1,
    magma_minproductFloatComplex_ptr dd2,
    magma_minproductFloatComplex_ptr dskp,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_cbicgmerge_spmv1( 
    magma_minproduct_c_matrix A,
    magma_minproductFloatComplex_ptr dd1,
    magma_minproductFloatComplex_ptr dd2,
    magma_minproductFloatComplex_ptr dp,
    magma_minproductFloatComplex_ptr dr,
    magma_minproductFloatComplex_ptr dv,
    magma_minproductFloatComplex_ptr dskp,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_cbicgmerge_spmv2( 
    magma_minproduct_c_matrix A,
    magma_minproductFloatComplex_ptr dd1,
    magma_minproductFloatComplex_ptr dd2,
    magma_minproductFloatComplex_ptr ds,
    magma_minproductFloatComplex_ptr dt,
    magma_minproductFloatComplex_ptr dskp,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_cbicgmerge_xrbeta( 
    int n,
    magma_minproductFloatComplex_ptr dd1,
    magma_minproductFloatComplex_ptr dd2,
    magma_minproductFloatComplex_ptr drr,
    magma_minproductFloatComplex_ptr dr,
    magma_minproductFloatComplex_ptr dp,
    magma_minproductFloatComplex_ptr ds,
    magma_minproductFloatComplex_ptr dt,
    magma_minproductFloatComplex_ptr dx, 
    magma_minproductFloatComplex_ptr dskp,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_cbcsrswp(
    magma_minproduct_int_t n,
    magma_minproduct_int_t size_b, 
    magma_minproduct_int_t *ipiv,
    magma_minproductFloatComplex_ptr dx,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_cbcsrtrsv(
    magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t r_blocks,
    magma_minproduct_int_t c_blocks,
    magma_minproduct_int_t size_b, 
    magma_minproductFloatComplex_ptr dA,
    magma_minproduct_index_t *blockinfo, 
    magma_minproductFloatComplex_ptr dx,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_cbcsrvalcpy(
    magma_minproduct_int_t size_b, 
    magma_minproduct_int_t num_blocks, 
    magma_minproduct_int_t num_zero_blocks, 
    magma_minproductFloatComplex_ptr *dAval, 
    magma_minproductFloatComplex_ptr *dBval,
    magma_minproductFloatComplex_ptr *dBval2,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_cbcsrluegemm(
    magma_minproduct_int_t size_b, 
    magma_minproduct_int_t num_block_rows,
    magma_minproduct_int_t kblocks,
    magma_minproductFloatComplex_ptr *dA, 
    magma_minproductFloatComplex_ptr *dB, 
    magma_minproductFloatComplex_ptr *dC,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_cbcsrlupivloc(
    magma_minproduct_int_t size_b, 
    magma_minproduct_int_t kblocks,
    magma_minproductFloatComplex_ptr *dA, 
    magma_minproduct_int_t *ipiv,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_cbcsrblockinfo5(
    magma_minproduct_int_t lustep,
    magma_minproduct_int_t num_blocks, 
    magma_minproduct_int_t c_blocks, 
    magma_minproduct_int_t size_b,
    magma_minproduct_index_t *blockinfo,
    magma_minproductFloatComplex_ptr dval,
    magma_minproductFloatComplex_ptr *AII,
    magma_minproduct_queue_t queue );

 
#ifdef __cplusplus
}
#endif

#undef PRECISION_c
#endif /* MAGMA_minproductSPARSE_C_H */
