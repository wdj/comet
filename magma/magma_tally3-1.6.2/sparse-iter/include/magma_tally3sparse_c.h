/*
 -- MAGMA_tally3 (version 1.1) --
 Univ. of Tennessee, Knoxville
 Univ. of California, Berkeley
 Univ. of Colorado, Denver
 @date

 @generated from magma_tally3sparse_z.h normal z -> c, Mon May  4 17:11:25 2015
 @author Hartwig Anzt
*/

#ifndef MAGMA_tally3SPARSE_C_H
#define MAGMA_tally3SPARSE_C_H

#include "magma_tally3_types.h"
#include "magma_tally3sparse_types.h"

#define PRECISION_c


#ifdef __cplusplus
extern "C" {
#endif


/* ////////////////////////////////////////////////////////////////////////////
 -- For backwards compatability, map old (1.6.1) to new (1.6.2) function names
*/

#define magma_tally3_c_mtranspose  magma_tally3_cmtranspose
#define magma_tally3_c_mtransfer   magma_tally3_cmtransfer
#define magma_tally3_c_vtransfer   magma_tally3_cmtransfer
#define magma_tally3_c_mconvert    magma_tally3_cmconvert
#define magma_tally3_c_vinit       magma_tally3_cvinit
#define magma_tally3_c_vvisu       magma_tally3_cprint_vector
#define magma_tally3_c_vread       magma_tally3_cvread
#define magma_tally3_c_vspread     magma_tally3_cvspread
#define magma_tally3_c_mvisu       magma_tally3_cprint_matrix
#define magma_tally3_c_mfree       magma_tally3_cmfree
#define magma_tally3_c_vfree       magma_tally3_cmfree
#define write_c_csr_mtx     magma_tally3_cwrite_csr_mtx
#define write_c_csrtomtx    magma_tally3_cwrite_csrtomtx
#define print_c_csr         magma_tally3_cprint_csr_mtx


/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally3_SPARSE Auxiliary functions
*/


magma_tally3_int_t
magma_tally3_cparse_opts( 
    int argc, 
    char** argv, 
    magma_tally3_copts *opts, 
    int *matrices, 
    magma_tally3_queue_t queue );

magma_tally3_int_t 
read_c_csr_from_binary( 
    magma_tally3_int_t* n_row, 
    magma_tally3_int_t* n_col, 
    magma_tally3_int_t* nnz, 
    magma_tally3FloatComplex **val, 
    magma_tally3_index_t **row, 
    magma_tally3_index_t **col,
    const char * filename,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
read_c_csr_from_mtx( 
    magma_tally3_storage_t *type, 
    magma_tally3_location_t *location,
    magma_tally3_int_t* n_row, 
    magma_tally3_int_t* n_col, 
    magma_tally3_int_t* nnz, 
    magma_tally3FloatComplex **val, 
    magma_tally3_index_t **row, 
    magma_tally3_index_t **col, 
    const char *filename,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_c_csr_mtx( 
    magma_tally3_c_matrix *A, 
    const char *filename,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_ccsrset( 
    magma_tally3_int_t m, 
    magma_tally3_int_t n, 
    magma_tally3_index_t *row, 
    magma_tally3_index_t *col, 
    magma_tally3FloatComplex *val,
    magma_tally3_c_matrix *A,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_ccsrget( 
    magma_tally3_c_matrix A,
    magma_tally3_int_t *m, 
    magma_tally3_int_t *n, 
    magma_tally3_index_t **row, 
    magma_tally3_index_t **col, 
    magma_tally3FloatComplex **val,
    magma_tally3_queue_t queue );


magma_tally3_int_t 
magma_tally3_cvset( 
    magma_tally3_int_t m, 
    magma_tally3_int_t n, 
    magma_tally3FloatComplex *val,
    magma_tally3_c_matrix *v,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_cvget( 
    magma_tally3_c_matrix v,
    magma_tally3_int_t *m, 
    magma_tally3_int_t *n, 
    magma_tally3FloatComplex **val,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_cvset_dev( 
    magma_tally3_int_t m, 
    magma_tally3_int_t n, 
    magma_tally3FloatComplex_ptr val,
    magma_tally3_c_matrix *v,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_cvget_dev( 
    magma_tally3_c_matrix v,
    magma_tally3_int_t *m, 
    magma_tally3_int_t *n, 
    magma_tally3FloatComplex_ptr *val,
    magma_tally3_queue_t queue );


magma_tally3_int_t 
magma_tally3_c_csr_mtxsymm( 
    magma_tally3_c_matrix *A, 
    const char *filename,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_c_csr_compressor( 
    magma_tally3FloatComplex ** val, 
    magma_tally3_index_t ** row, 
    magma_tally3_index_t ** col, 
    magma_tally3FloatComplex ** valn, 
    magma_tally3_index_t ** rown, 
    magma_tally3_index_t ** coln, 
    magma_tally3_int_t *n,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_cmcsrcompressor( 
    magma_tally3_c_matrix *A,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_cmcsrcompressor_gpu( 
    magma_tally3_c_matrix *A,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_cvtranspose( 
    magma_tally3_c_matrix x,
    magma_tally3_c_matrix *y,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_c_cucsrtranspose( 
    magma_tally3_c_matrix A, 
    magma_tally3_c_matrix *B,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
c_transpose_csr( 
    magma_tally3_int_t n_rows, 
    magma_tally3_int_t n_cols, 
    magma_tally3_int_t nnz,
    magma_tally3FloatComplex *val, 
    magma_tally3_index_t *row, 
    magma_tally3_index_t *col, 
    magma_tally3_int_t *new_n_rows, 
    magma_tally3_int_t *new_n_cols, 
    magma_tally3_int_t *new_nnz, 
    magma_tally3FloatComplex **new_val, 
    magma_tally3_index_t **new_row, 
    magma_tally3_index_t **new_col,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_ccsrsplit( 
    magma_tally3_int_t bsize,
    magma_tally3_c_matrix A,
    magma_tally3_c_matrix *D,
    magma_tally3_c_matrix *R,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_cmscale( 
    magma_tally3_c_matrix *A, 
    magma_tally3_scale_t scaling,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_cmdiff( 
    magma_tally3_c_matrix A, 
    magma_tally3_c_matrix B, 
 real_Double_t *res,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_cmdiagadd( 
    magma_tally3_c_matrix *A, 
    magma_tally3FloatComplex add,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_cmsort( 
    magma_tally3_c_matrix *A,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_cindexsort(
    magma_tally3_index_t *x, 
    magma_tally3_int_t first,
    magma_tally3_int_t last,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_cdomainoverlap(
    magma_tally3_index_t num_rows,
    magma_tally3_index_t *num_indices,
    magma_tally3_index_t *rowptr,
    magma_tally3_index_t *colidx,
    magma_tally3_index_t *x,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_csymbilu( 
    magma_tally3_c_matrix *A, 
    magma_tally3_int_t levels,
    magma_tally3_c_matrix *L,
    magma_tally3_c_matrix *U,
    magma_tally3_queue_t queue );


magma_tally3_int_t 
magma_tally3_cwrite_csr_mtx( 
    magma_tally3_c_matrix A,
    magma_tally3_order_t MajorType,
 const char *filename,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_cwrite_csrtomtx( 
    magma_tally3_c_matrix A,
    const char *filename,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_cprint_csr( 
    magma_tally3_int_t n_row, 
    magma_tally3_int_t n_col, 
    magma_tally3_int_t nnz, 
    magma_tally3FloatComplex **val, 
    magma_tally3_index_t **row, 
    magma_tally3_index_t **col,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_cprint_csr_mtx( 
    magma_tally3_int_t n_row, 
    magma_tally3_int_t n_col, 
    magma_tally3_int_t nnz, 
    magma_tally3FloatComplex **val, 
    magma_tally3_index_t **row, 
    magma_tally3_index_t **col, 
    magma_tally3_order_t MajorType,
    magma_tally3_queue_t queue );


magma_tally3_int_t 
magma_tally3_cmtranspose(
    magma_tally3_c_matrix A, 
    magma_tally3_c_matrix *B,
    magma_tally3_queue_t queue );


magma_tally3_int_t 
magma_tally3_cmtransfer(
    magma_tally3_c_matrix A, 
    magma_tally3_c_matrix *B, 
    magma_tally3_location_t src, 
    magma_tally3_location_t dst,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_cmconvert(
    magma_tally3_c_matrix A, 
    magma_tally3_c_matrix *B, 
    magma_tally3_storage_t old_format, 
    magma_tally3_storage_t new_format,
    magma_tally3_queue_t queue );


magma_tally3_int_t
magma_tally3_cvinit(
    magma_tally3_c_matrix *x, 
    magma_tally3_location_t memory_location,
    magma_tally3_int_t num_rows, 
    magma_tally3_int_t num_cols,
    magma_tally3FloatComplex values,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_cprint_vector(
    magma_tally3_c_matrix x, 
    magma_tally3_int_t offset, 
    magma_tally3_int_t displaylength,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_cvread(
    magma_tally3_c_matrix *x, 
    magma_tally3_int_t length,
    char * filename,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_cvspread(
    magma_tally3_c_matrix *x, 
    const char * filename,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_cprint_matrix(
    magma_tally3_c_matrix A,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_cdiameter(
    magma_tally3_c_matrix *A,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_crowentries(
    magma_tally3_c_matrix *A,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_cmfree(
    magma_tally3_c_matrix *A,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_cresidual(
    magma_tally3_c_matrix A, 
    magma_tally3_c_matrix b, 
    magma_tally3_c_matrix x, 
    float *res,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_cresidualvec(
    magma_tally3_c_matrix A,
    magma_tally3_c_matrix b,
    magma_tally3_c_matrix x,
    magma_tally3_c_matrix *r,
    float *res,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_cmgenerator(
    magma_tally3_int_t n,
    magma_tally3_int_t offdiags,
    magma_tally3_index_t *diag_offset,
    magma_tally3FloatComplex *diag_vals,
    magma_tally3_c_matrix *A,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_cm_27stencil(
    magma_tally3_int_t n,
    magma_tally3_c_matrix *A,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_cm_5stencil(
    magma_tally3_int_t n,
    magma_tally3_c_matrix *A,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_csolverinfo(
    magma_tally3_c_solver_par *solver_par, 
    magma_tally3_c_preconditioner *precond_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_csolverinfo_init(
    magma_tally3_c_solver_par *solver_par, 
    magma_tally3_c_preconditioner *precond,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_ceigensolverinfo_init(
    magma_tally3_c_solver_par *solver_par,
    magma_tally3_queue_t queue );


magma_tally3_int_t
magma_tally3_csolverinfo_free(
    magma_tally3_c_solver_par *solver_par, 
    magma_tally3_c_preconditioner *precond,
    magma_tally3_queue_t queue );



/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally3_SPARSE iterative incomplete factorizations
*/


magma_tally3_int_t
magma_tally3_citerilusetup( 
    magma_tally3_c_matrix A, 
    magma_tally3_c_matrix b,                                 
    magma_tally3_c_preconditioner *precond,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_citericsetup( 
    magma_tally3_c_matrix A, 
    magma_tally3_c_matrix b, 
    magma_tally3_c_preconditioner *precond,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_citericupdate( 
    magma_tally3_c_matrix A, 
    magma_tally3_c_preconditioner *precond, 
    magma_tally3_int_t updates,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_capplyiteric_l( 
    magma_tally3_c_matrix b, 
    magma_tally3_c_matrix *x, 
    magma_tally3_c_preconditioner *precond,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_capplyiteric_r( 
    magma_tally3_c_matrix b, 
    magma_tally3_c_matrix *x, 
    magma_tally3_c_preconditioner *precond,
    magma_tally3_queue_t queue );


magma_tally3_int_t
magma_tally3_citerilu_csr( 
    magma_tally3_c_matrix A,
    magma_tally3_c_matrix L,
    magma_tally3_c_matrix U,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_citeric_csr( 
    magma_tally3_c_matrix A,
    magma_tally3_c_matrix A_CSR,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_cfrobenius( 
    magma_tally3_c_matrix A, 
    magma_tally3_c_matrix B, 
    real_Double_t *res,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_cnonlinres(   
    magma_tally3_c_matrix A, 
    magma_tally3_c_matrix L,
    magma_tally3_c_matrix U, 
    magma_tally3_c_matrix *LU, 
    real_Double_t *res,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_cilures(   
    magma_tally3_c_matrix A, 
    magma_tally3_c_matrix L,
    magma_tally3_c_matrix U, 
    magma_tally3_c_matrix *LU, 
    real_Double_t *res,
    real_Double_t *nonlinres,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_cicres(       
    magma_tally3_c_matrix A, 
    magma_tally3_c_matrix C,
    magma_tally3_c_matrix CT, 
    magma_tally3_c_matrix *LU, 
    real_Double_t *res,
    real_Double_t *nonlinres,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_cinitguess( 
    magma_tally3_c_matrix A, 
    magma_tally3_c_matrix *L, 
    magma_tally3_c_matrix *U,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_cinitrecursiveLU( 
    magma_tally3_c_matrix A, 
    magma_tally3_c_matrix *B,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_cmLdiagadd( 
    magma_tally3_c_matrix *L,
    magma_tally3_queue_t queue );




/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally3_SPARSE function definitions / Data on CPU
*/




/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally3_SPARSE function definitions / Data on CPU / Multi-GPU
*/

/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally3_SPARSE solvers (Data on GPU)
*/

magma_tally3_int_t 
magma_tally3_ccg(
    magma_tally3_c_matrix A, magma_tally3_c_matrix b, 
    magma_tally3_c_matrix *x, magma_tally3_c_solver_par *solver_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_ccg_res(
    magma_tally3_c_matrix A, magma_tally3_c_matrix b, 
    magma_tally3_c_matrix *x, magma_tally3_c_solver_par *solver_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_ccg_merge(
    magma_tally3_c_matrix A, magma_tally3_c_matrix b, 
    magma_tally3_c_matrix *x, magma_tally3_c_solver_par *solver_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_cgmres(
    magma_tally3_c_matrix A, magma_tally3_c_matrix b, 
    magma_tally3_c_matrix *x, magma_tally3_c_solver_par *solver_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_cbicgstab(
    magma_tally3_c_matrix A, magma_tally3_c_matrix b, magma_tally3_c_matrix *x, 
    magma_tally3_c_solver_par *solver_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_cbicgstab_merge(
    magma_tally3_c_matrix A, magma_tally3_c_matrix b, 
    magma_tally3_c_matrix *x, magma_tally3_c_solver_par *solver_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_cbicgstab_merge2(
    magma_tally3_c_matrix A, magma_tally3_c_matrix b, 
    magma_tally3_c_matrix *x, magma_tally3_c_solver_par *solver_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_cpcg(
    magma_tally3_c_matrix A, magma_tally3_c_matrix b, 
    magma_tally3_c_matrix *x, magma_tally3_c_solver_par *solver_par, 
    magma_tally3_c_preconditioner *precond_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_cbpcg(
    magma_tally3_c_matrix A, magma_tally3_c_matrix b, 
    magma_tally3_c_matrix *x, magma_tally3_c_solver_par *solver_par, 
    magma_tally3_c_preconditioner *precond_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_cpbicgstab(
    magma_tally3_c_matrix A, magma_tally3_c_matrix b, 
    magma_tally3_c_matrix *x, magma_tally3_c_solver_par *solver_par, 
    magma_tally3_c_preconditioner *precond_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_cpgmres(
    magma_tally3_c_matrix A, magma_tally3_c_matrix b, 
    magma_tally3_c_matrix *x, magma_tally3_c_solver_par *solver_par, 
    magma_tally3_c_preconditioner *precond_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_cfgmres(
    magma_tally3_c_matrix A, magma_tally3_c_matrix b, 
    magma_tally3_c_matrix *x, magma_tally3_c_solver_par *solver_par, 
    magma_tally3_c_preconditioner *precond_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_cjacobi(
    magma_tally3_c_matrix A, magma_tally3_c_matrix b, 
    magma_tally3_c_matrix *x, magma_tally3_c_solver_par *solver_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_cjacobidomainoverlap(
    magma_tally3_c_matrix A, 
    magma_tally3_c_matrix b, 
    magma_tally3_c_matrix *x,  
    magma_tally3_c_solver_par *solver_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_cbaiter(
    magma_tally3_c_matrix A, magma_tally3_c_matrix b, 
    magma_tally3_c_matrix *x, magma_tally3_c_solver_par *solver_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_citerref(
    magma_tally3_c_matrix A, magma_tally3_c_matrix b, 
    magma_tally3_c_matrix *x, magma_tally3_c_solver_par *solver_par, 
    magma_tally3_c_preconditioner *precond_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_cilu(
    magma_tally3_c_matrix A, magma_tally3_c_matrix b, 
    magma_tally3_c_matrix *x, magma_tally3_c_solver_par *solver_par, 
    magma_tally3_int_t *ipiv,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_cbcsrlu(
    magma_tally3_c_matrix A, magma_tally3_c_matrix b, 
    magma_tally3_c_matrix *x, magma_tally3_c_solver_par *solver_par,
    magma_tally3_queue_t queue );


magma_tally3_int_t
magma_tally3_cbcsrlutrf(
    magma_tally3_c_matrix A, 
    magma_tally3_c_matrix *M,
    magma_tally3_int_t *ipiv, 
    magma_tally3_int_t version,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_cbcsrlusv(
    magma_tally3_c_matrix A, magma_tally3_c_matrix b, 
    magma_tally3_c_matrix *x, magma_tally3_c_solver_par *solver_par, 
    magma_tally3_int_t *ipiv,
    magma_tally3_queue_t queue );



magma_tally3_int_t
magma_tally3_cilucg(
    magma_tally3_c_matrix A, magma_tally3_c_matrix b, 
    magma_tally3_c_matrix *x, magma_tally3_c_solver_par *solver_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_cilugmres(
    magma_tally3_c_matrix A, magma_tally3_c_matrix b, 
    magma_tally3_c_matrix *x, magma_tally3_c_solver_par *solver_par,
    magma_tally3_queue_t queue ); 


magma_tally3_int_t
magma_tally3_clobpcg_shift(
    magma_tally3_int_t num_rows,
    magma_tally3_int_t num_vecs, 
    magma_tally3_int_t shift,
    magma_tally3FloatComplex_ptr x,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_clobpcg_res(
    magma_tally3_int_t num_rows,
    magma_tally3_int_t num_vecs, 
    float *evalues, 
    magma_tally3FloatComplex_ptr X,
    magma_tally3FloatComplex_ptr R, 
    float *res,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_clobpcg_maxpy(
    magma_tally3_int_t num_rows,
    magma_tally3_int_t num_vecs, 
    magma_tally3FloatComplex_ptr X,
    magma_tally3FloatComplex_ptr Y,
    magma_tally3_queue_t queue );


/*/////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally3_SPARSE eigensolvers (Data on GPU)
*/
magma_tally3_int_t
magma_tally3_clobpcg(
    magma_tally3_c_matrix A, 
    magma_tally3_c_solver_par *solver_par,
    magma_tally3_c_preconditioner *precond_par, 
    magma_tally3_queue_t queue );




/*/////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally3_SPARSE preconditioners (Data on GPU)
*/

magma_tally3_int_t
magma_tally3_cjacobisetup(
    magma_tally3_c_matrix A, 
    magma_tally3_c_matrix b, 
    magma_tally3_c_matrix *M, 
    magma_tally3_c_matrix *c,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_cjacobisetup_matrix(
    magma_tally3_c_matrix A, 
    magma_tally3_c_matrix *M, 
    magma_tally3_c_matrix *d,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_cjacobisetup_vector(
    magma_tally3_c_matrix b, 
    magma_tally3_c_matrix d, 
    magma_tally3_c_matrix *c,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_cjacobiiter(
    magma_tally3_c_matrix M, 
    magma_tally3_c_matrix c, 
    magma_tally3_c_matrix *x, 
    magma_tally3_c_solver_par *solver_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_cjacobiiter_precond( 
    magma_tally3_c_matrix M, 
    magma_tally3_c_matrix *x, 
    magma_tally3_c_solver_par *solver_par, 
    magma_tally3_c_preconditioner *precond,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_cjacobiiter_sys(
    magma_tally3_c_matrix A, 
    magma_tally3_c_matrix b, 
    magma_tally3_c_matrix d, 
    magma_tally3_c_matrix t, 
    magma_tally3_c_matrix *x,  
    magma_tally3_c_solver_par *solver_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_cpastixsetup(
    magma_tally3_c_matrix A, magma_tally3_c_matrix b,
    magma_tally3_c_preconditioner *precond,
    magma_tally3_queue_t queue );


magma_tally3_int_t
magma_tally3_capplypastix(
    magma_tally3_c_matrix b, magma_tally3_c_matrix *x, 
    magma_tally3_c_preconditioner *precond,
    magma_tally3_queue_t queue );


// custom preconditioner
magma_tally3_int_t
magma_tally3_capplycustomprecond_l(
    magma_tally3_c_matrix b, 
    magma_tally3_c_matrix *x, 
    magma_tally3_c_preconditioner *precond,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_capplycustomprecond_r(
    magma_tally3_c_matrix b, 
    magma_tally3_c_matrix *x, 
    magma_tally3_c_preconditioner *precond,
    magma_tally3_queue_t queue );


// CUSPARSE preconditioner

magma_tally3_int_t
magma_tally3_ccuilusetup(
    magma_tally3_c_matrix A, magma_tally3_c_preconditioner *precond,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_capplycuilu_l(
    magma_tally3_c_matrix b, magma_tally3_c_matrix *x, 
    magma_tally3_c_preconditioner *precond,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_capplycuilu_r(
    magma_tally3_c_matrix b, magma_tally3_c_matrix *x, 
    magma_tally3_c_preconditioner *precond,
    magma_tally3_queue_t queue );


magma_tally3_int_t
magma_tally3_ccuiccsetup(
    magma_tally3_c_matrix A, magma_tally3_c_preconditioner *precond,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_capplycuicc_l(
    magma_tally3_c_matrix b, magma_tally3_c_matrix *x, 
    magma_tally3_c_preconditioner *precond,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_capplycuicc_r(
    magma_tally3_c_matrix b, magma_tally3_c_matrix *x, 
    magma_tally3_c_preconditioner *precond,
    magma_tally3_queue_t queue );


magma_tally3_int_t
magma_tally3_ccumilusetup(
    magma_tally3_c_matrix A, 
    magma_tally3_c_preconditioner *precond,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_ccumilugeneratesolverinfo(
    magma_tally3_c_preconditioner *precond,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_capplycumilu_l(
    magma_tally3_c_matrix b, 
    magma_tally3_c_matrix *x, 
    magma_tally3_c_preconditioner *precond,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_capplycumilu_r(
    magma_tally3_c_matrix b, 
    magma_tally3_c_matrix *x, 
    magma_tally3_c_preconditioner *precond,
    magma_tally3_queue_t queue );


magma_tally3_int_t
magma_tally3_ccumiccsetup(
    magma_tally3_c_matrix A, 
    magma_tally3_c_preconditioner *precond,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_ccumicgeneratesolverinfo(
    magma_tally3_c_preconditioner *precond,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_capplycumicc_l(
    magma_tally3_c_matrix b, 
    magma_tally3_c_matrix *x, 
    magma_tally3_c_preconditioner *precond,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_capplycumicc_r(
    magma_tally3_c_matrix b, 
    magma_tally3_c_matrix *x, 
    magma_tally3_c_preconditioner *precond,
    magma_tally3_queue_t queue );


// block-asynchronous iteration

magma_tally3_int_t
magma_tally3_cbajac_csr(
    magma_tally3_int_t localiters,
    magma_tally3_c_matrix D,
    magma_tally3_c_matrix R,
    magma_tally3_c_matrix b,
    magma_tally3_c_matrix *x,
    magma_tally3_queue_t queue );

/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally3_SPARSE utility function definitions
*/

magma_tally3_int_t
magma_tally3_c_spmv(
    magma_tally3FloatComplex alpha, 
    magma_tally3_c_matrix A, 
    magma_tally3_c_matrix x, 
    magma_tally3FloatComplex beta, 
    magma_tally3_c_matrix y,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_ccustomspmv(
    magma_tally3FloatComplex alpha, 
    magma_tally3_c_matrix x, 
    magma_tally3FloatComplex beta, 
    magma_tally3_c_matrix y,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_c_spmv_shift(
    magma_tally3FloatComplex alpha, 
    magma_tally3_c_matrix A, 
    magma_tally3FloatComplex lambda,
    magma_tally3_c_matrix x, 
    magma_tally3FloatComplex beta, 
    magma_tally3_int_t offset, 
    magma_tally3_int_t blocksize,
    magma_tally3Index_ptr dadd_vecs, 
    magma_tally3_c_matrix y,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_ccuspmm(
    magma_tally3_c_matrix A, 
    magma_tally3_c_matrix B, 
    magma_tally3_c_matrix *AB,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_c_spmm(
    magma_tally3FloatComplex alpha, 
    magma_tally3_c_matrix A,
    magma_tally3_c_matrix B,
    magma_tally3_c_matrix *C,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_csymbilu( 
    magma_tally3_c_matrix *A, 
    magma_tally3_int_t levels, 
    magma_tally3_c_matrix *L, 
    magma_tally3_c_matrix *U,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_ccuspaxpy(
    magma_tally3FloatComplex_ptr alpha, magma_tally3_c_matrix A, 
    magma_tally3FloatComplex_ptr beta, magma_tally3_c_matrix B, 
    magma_tally3_c_matrix *AB,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_c_precond(
    magma_tally3_c_matrix A, 
    magma_tally3_c_matrix b, magma_tally3_c_matrix *x,
    magma_tally3_c_preconditioner *precond,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_c_solver(
    magma_tally3_c_matrix A, magma_tally3_c_matrix b, 
    magma_tally3_c_matrix *x, magma_tally3_copts *zopts,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_c_precondsetup(
    magma_tally3_c_matrix A, magma_tally3_c_matrix b, 
    magma_tally3_c_preconditioner *precond,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_c_applyprecond(
    magma_tally3_c_matrix A, magma_tally3_c_matrix b, 
    magma_tally3_c_matrix *x, magma_tally3_c_preconditioner *precond,
    magma_tally3_queue_t queue );


magma_tally3_int_t
magma_tally3_c_applyprecond_left(
    magma_tally3_c_matrix A, magma_tally3_c_matrix b, 
    magma_tally3_c_matrix *x, magma_tally3_c_preconditioner *precond,
    magma_tally3_queue_t queue );


magma_tally3_int_t
magma_tally3_c_applyprecond_right(
    magma_tally3_c_matrix A, magma_tally3_c_matrix b, 
    magma_tally3_c_matrix *x, magma_tally3_c_preconditioner *precond,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_c_initP2P(
    magma_tally3_int_t *bandwidth_benchmark,
    magma_tally3_int_t *num_gpus,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_ccompact(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex_ptr dA, magma_tally3_int_t ldda,
    float *dnorms, float tol, 
    magma_tally3_int_t *activeMask, magma_tally3_int_t *cBlockSize,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_ccompactActive(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex_ptr dA, magma_tally3_int_t ldda, 
    magma_tally3_int_t *active,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_cmlumerge(    
    magma_tally3_c_matrix L, 
    magma_tally3_c_matrix U,
    magma_tally3_c_matrix *A, 
    magma_tally3_queue_t queue );


/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally3_SPARSE BLAS function definitions
*/
magma_tally3_int_t 
magma_tally3_cgecsrmv(
    magma_tally3_trans_t transA,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex alpha,
    magma_tally3FloatComplex_ptr dval,
    magma_tally3Index_ptr drowptr,
    magma_tally3Index_ptr dcolind,
    magma_tally3FloatComplex_ptr dx,
    magma_tally3FloatComplex beta,
    magma_tally3FloatComplex_ptr dy,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_cgecsrmv_shift(
    magma_tally3_trans_t transA,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3FloatComplex alpha,
    magma_tally3FloatComplex lambda,
    magma_tally3FloatComplex_ptr dval,
    magma_tally3Index_ptr drowptr,
    magma_tally3Index_ptr dcolind,
    magma_tally3FloatComplex_ptr dx,
    magma_tally3FloatComplex beta,
    int offset,
    int blocksize,
    magma_tally3Index_ptr dadd_rows,
    magma_tally3FloatComplex_ptr dy,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_cmgecsrmv(
    magma_tally3_trans_t transA,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3_int_t num_vecs,
    magma_tally3FloatComplex alpha,
    magma_tally3FloatComplex_ptr dval,
    magma_tally3Index_ptr drowptr,
    magma_tally3Index_ptr dcolind,
    magma_tally3FloatComplex_ptr dx,
    magma_tally3FloatComplex beta,
    magma_tally3FloatComplex_ptr dy,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_cgeellmv(
    magma_tally3_trans_t transA,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3_int_t nnz_per_row,
    magma_tally3FloatComplex alpha,
    magma_tally3FloatComplex_ptr dval,
    magma_tally3Index_ptr dcolind,
    magma_tally3FloatComplex_ptr dx,
    magma_tally3FloatComplex beta,
    magma_tally3FloatComplex_ptr dy,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_cgeellmv_shift(
    magma_tally3_trans_t transA,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3_int_t nnz_per_row,
    magma_tally3FloatComplex alpha,
    magma_tally3FloatComplex lambda,
    magma_tally3FloatComplex_ptr dval,
    magma_tally3Index_ptr dcolind,
    magma_tally3FloatComplex_ptr dx,
    magma_tally3FloatComplex beta,
    int offset,
    int blocksize,
    magma_tally3Index_ptr dadd_rows,
    magma_tally3FloatComplex_ptr dy,
    magma_tally3_queue_t queue );


magma_tally3_int_t 
magma_tally3_cmgeellmv(
    magma_tally3_trans_t transA,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3_int_t num_vecs,
    magma_tally3_int_t nnz_per_row,
    magma_tally3FloatComplex alpha,
    magma_tally3FloatComplex_ptr dval,
    magma_tally3Index_ptr dcolind,
    magma_tally3FloatComplex_ptr dx,
    magma_tally3FloatComplex beta,
    magma_tally3FloatComplex_ptr dy,
    magma_tally3_queue_t queue );


magma_tally3_int_t 
magma_tally3_cgeelltmv(
    magma_tally3_trans_t transA,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3_int_t nnz_per_row,
    magma_tally3FloatComplex alpha,
    magma_tally3FloatComplex_ptr dval,
    magma_tally3Index_ptr dcolind,
    magma_tally3FloatComplex_ptr dx,
    magma_tally3FloatComplex beta,
    magma_tally3FloatComplex_ptr dy,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_cgeelltmv_shift(
    magma_tally3_trans_t transA,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3_int_t nnz_per_row,
    magma_tally3FloatComplex alpha,
    magma_tally3FloatComplex lambda,
    magma_tally3FloatComplex_ptr dval,
    magma_tally3Index_ptr dcolind,
    magma_tally3FloatComplex_ptr dx,
    magma_tally3FloatComplex beta,
    int offset,
    int blocksize,
    magma_tally3Index_ptr dadd_rows,
    magma_tally3FloatComplex_ptr dy,
    magma_tally3_queue_t queue );


magma_tally3_int_t 
magma_tally3_cmgeelltmv(
    magma_tally3_trans_t transA,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3_int_t num_vecs,
    magma_tally3_int_t nnz_per_row,
    magma_tally3FloatComplex alpha,
    magma_tally3FloatComplex_ptr dval,
    magma_tally3Index_ptr dcolind,
    magma_tally3FloatComplex_ptr dx,
    magma_tally3FloatComplex beta,
    magma_tally3FloatComplex_ptr dy,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_cgeellrtmv(
    magma_tally3_trans_t transA,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3_int_t nnz_per_row,
    magma_tally3FloatComplex alpha,
    magma_tally3FloatComplex_ptr dval,
    magma_tally3Index_ptr dcolind,
    magma_tally3Index_ptr drowlength,
    magma_tally3FloatComplex_ptr dx,
    magma_tally3FloatComplex beta,
    magma_tally3FloatComplex_ptr dy,
    magma_tally3_int_t num_threads,
    magma_tally3_int_t threads_per_row,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_cgesellcmv(
    magma_tally3_trans_t transA,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3_int_t blocksize,
    magma_tally3_int_t slices,
    magma_tally3_int_t alignment,
    magma_tally3FloatComplex alpha,
    magma_tally3FloatComplex_ptr dval,
    magma_tally3Index_ptr dcolind,
    magma_tally3Index_ptr drowptr,
    magma_tally3FloatComplex_ptr dx,
    magma_tally3FloatComplex beta,
    magma_tally3FloatComplex_ptr dy,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_cgesellpmv(
    magma_tally3_trans_t transA,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3_int_t blocksize,
    magma_tally3_int_t slices,
    magma_tally3_int_t alignment,
    magma_tally3FloatComplex alpha,
    magma_tally3FloatComplex_ptr dval,
    magma_tally3Index_ptr dcolind,
    magma_tally3Index_ptr drowptr,
    magma_tally3FloatComplex_ptr dx,
    magma_tally3FloatComplex beta,
    magma_tally3FloatComplex_ptr dy,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_cmgesellpmv(
    magma_tally3_trans_t transA,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3_int_t num_vecs,
    magma_tally3_int_t blocksize,
    magma_tally3_int_t slices,
    magma_tally3_int_t alignment,
    magma_tally3FloatComplex alpha,
    magma_tally3FloatComplex_ptr dval,
    magma_tally3Index_ptr dcolind,
    magma_tally3Index_ptr drowptr,
    magma_tally3FloatComplex_ptr dx,
    magma_tally3FloatComplex beta,
    magma_tally3FloatComplex_ptr dy,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_cmgesellpmv_blocked(
    magma_tally3_trans_t transA,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3_int_t num_vecs,
    magma_tally3_int_t blocksize,
    magma_tally3_int_t slices,
    magma_tally3_int_t alignment,
    magma_tally3FloatComplex alpha,
    magma_tally3FloatComplex_ptr dval,
    magma_tally3Index_ptr dcolind,
    magma_tally3Index_ptr drowptr,
    magma_tally3FloatComplex_ptr dx,
    magma_tally3FloatComplex beta,
    magma_tally3FloatComplex_ptr dy,
    magma_tally3_queue_t queue );


magma_tally3_int_t
magma_tally3_cmergedgs(
    magma_tally3_int_t n, 
    magma_tally3_int_t ldh,
    magma_tally3_int_t k, 
    magma_tally3FloatComplex_ptr dv, 
    magma_tally3FloatComplex_ptr dr,
    magma_tally3FloatComplex_ptr dskp,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_ccopyscale(    
    int n, 
    int k,
    magma_tally3FloatComplex_ptr dr, 
    magma_tally3FloatComplex_ptr dv,
    magma_tally3FloatComplex_ptr dskp,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_scnrm2scale(    
    int m, 
    magma_tally3FloatComplex_ptr dr,    
    int lddr, 
    magma_tally3FloatComplex *drnorm,
    magma_tally3_queue_t queue );


magma_tally3_int_t
magma_tally3_cjacobisetup_vector_gpu(
    int num_rows, 
    magma_tally3_c_matrix b, 
    magma_tally3_c_matrix d, 
    magma_tally3_c_matrix c,
    magma_tally3_c_matrix *x,
    magma_tally3_queue_t queue );


magma_tally3_int_t
magma_tally3_cjacobi_diagscal(    
    int num_rows, 
    magma_tally3_c_matrix d, 
    magma_tally3_c_matrix b, 
    magma_tally3_c_matrix *c,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_cjacobiupdate(
    magma_tally3_c_matrix t, 
    magma_tally3_c_matrix b, 
    magma_tally3_c_matrix d, 
    magma_tally3_c_matrix *x,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_cjacobispmvupdate(
    magma_tally3_int_t maxiter,
    magma_tally3_c_matrix A, 
    magma_tally3_c_matrix t, 
    magma_tally3_c_matrix b, 
    magma_tally3_c_matrix d, 
    magma_tally3_c_matrix *x,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_cjacobispmvupdate_bw(
    magma_tally3_int_t maxiter,
    magma_tally3_c_matrix A, 
    magma_tally3_c_matrix t, 
    magma_tally3_c_matrix b, 
    magma_tally3_c_matrix d, 
    magma_tally3_c_matrix *x,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_cjacobispmvupdateselect(
    magma_tally3_int_t maxiter,
    magma_tally3_int_t num_updates,
    magma_tally3_index_t *indices,
    magma_tally3_c_matrix A,
    magma_tally3_c_matrix t, 
    magma_tally3_c_matrix b, 
    magma_tally3_c_matrix d, 
    magma_tally3_c_matrix tmp, 
    magma_tally3_c_matrix *x,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_cjacobisetup_diagscal(
    magma_tally3_c_matrix A, magma_tally3_c_matrix *d,
    magma_tally3_queue_t queue );


magma_tally3_int_t
magma_tally3_cbicgmerge1(    
    int n, 
    magma_tally3FloatComplex_ptr dskp,
    magma_tally3FloatComplex_ptr dv, 
    magma_tally3FloatComplex_ptr dr, 
    magma_tally3FloatComplex_ptr dp,
    magma_tally3_queue_t queue );


magma_tally3_int_t
magma_tally3_cbicgmerge2(
    int n, 
    magma_tally3FloatComplex_ptr dskp, 
    magma_tally3FloatComplex_ptr dr,
    magma_tally3FloatComplex_ptr dv, 
    magma_tally3FloatComplex_ptr ds,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_cbicgmerge3(
    int n, 
    magma_tally3FloatComplex_ptr dskp, 
    magma_tally3FloatComplex_ptr dp,
    magma_tally3FloatComplex_ptr ds,
    magma_tally3FloatComplex_ptr dt,
    magma_tally3FloatComplex_ptr dx, 
    magma_tally3FloatComplex_ptr dr,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_cbicgmerge4(
    int type, 
    magma_tally3FloatComplex_ptr dskp,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_ccgmerge_spmv1( 
    magma_tally3_c_matrix A,
    magma_tally3FloatComplex_ptr d1,
    magma_tally3FloatComplex_ptr d2,
    magma_tally3FloatComplex_ptr dd,
    magma_tally3FloatComplex_ptr dz,
    magma_tally3FloatComplex_ptr dskp,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_ccgmerge_xrbeta( 
    int n,
    magma_tally3FloatComplex_ptr d1,
    magma_tally3FloatComplex_ptr d2,
    magma_tally3FloatComplex_ptr dx,
    magma_tally3FloatComplex_ptr dr,
    magma_tally3FloatComplex_ptr dd,
    magma_tally3FloatComplex_ptr dz, 
    magma_tally3FloatComplex_ptr dskp,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_cmdotc(
    magma_tally3_int_t n, 
    magma_tally3_int_t k, 
    magma_tally3FloatComplex_ptr dv, 
    magma_tally3FloatComplex_ptr dr,
    magma_tally3FloatComplex_ptr dd1,
    magma_tally3FloatComplex_ptr dd2,
    magma_tally3FloatComplex_ptr dskp,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_cgemvmdot(
    int n, 
    int k, 
    magma_tally3FloatComplex_ptr dv, 
    magma_tally3FloatComplex_ptr dr,
    magma_tally3FloatComplex_ptr dd1,
    magma_tally3FloatComplex_ptr dd2,
    magma_tally3FloatComplex_ptr dskp,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_cbicgmerge_spmv1( 
    magma_tally3_c_matrix A,
    magma_tally3FloatComplex_ptr dd1,
    magma_tally3FloatComplex_ptr dd2,
    magma_tally3FloatComplex_ptr dp,
    magma_tally3FloatComplex_ptr dr,
    magma_tally3FloatComplex_ptr dv,
    magma_tally3FloatComplex_ptr dskp,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_cbicgmerge_spmv2( 
    magma_tally3_c_matrix A,
    magma_tally3FloatComplex_ptr dd1,
    magma_tally3FloatComplex_ptr dd2,
    magma_tally3FloatComplex_ptr ds,
    magma_tally3FloatComplex_ptr dt,
    magma_tally3FloatComplex_ptr dskp,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_cbicgmerge_xrbeta( 
    int n,
    magma_tally3FloatComplex_ptr dd1,
    magma_tally3FloatComplex_ptr dd2,
    magma_tally3FloatComplex_ptr drr,
    magma_tally3FloatComplex_ptr dr,
    magma_tally3FloatComplex_ptr dp,
    magma_tally3FloatComplex_ptr ds,
    magma_tally3FloatComplex_ptr dt,
    magma_tally3FloatComplex_ptr dx, 
    magma_tally3FloatComplex_ptr dskp,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_cbcsrswp(
    magma_tally3_int_t n,
    magma_tally3_int_t size_b, 
    magma_tally3_int_t *ipiv,
    magma_tally3FloatComplex_ptr dx,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_cbcsrtrsv(
    magma_tally3_uplo_t uplo,
    magma_tally3_int_t r_blocks,
    magma_tally3_int_t c_blocks,
    magma_tally3_int_t size_b, 
    magma_tally3FloatComplex_ptr dA,
    magma_tally3_index_t *blockinfo, 
    magma_tally3FloatComplex_ptr dx,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_cbcsrvalcpy(
    magma_tally3_int_t size_b, 
    magma_tally3_int_t num_blocks, 
    magma_tally3_int_t num_zero_blocks, 
    magma_tally3FloatComplex_ptr *dAval, 
    magma_tally3FloatComplex_ptr *dBval,
    magma_tally3FloatComplex_ptr *dBval2,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_cbcsrluegemm(
    magma_tally3_int_t size_b, 
    magma_tally3_int_t num_block_rows,
    magma_tally3_int_t kblocks,
    magma_tally3FloatComplex_ptr *dA, 
    magma_tally3FloatComplex_ptr *dB, 
    magma_tally3FloatComplex_ptr *dC,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_cbcsrlupivloc(
    magma_tally3_int_t size_b, 
    magma_tally3_int_t kblocks,
    magma_tally3FloatComplex_ptr *dA, 
    magma_tally3_int_t *ipiv,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_cbcsrblockinfo5(
    magma_tally3_int_t lustep,
    magma_tally3_int_t num_blocks, 
    magma_tally3_int_t c_blocks, 
    magma_tally3_int_t size_b,
    magma_tally3_index_t *blockinfo,
    magma_tally3FloatComplex_ptr dval,
    magma_tally3FloatComplex_ptr *AII,
    magma_tally3_queue_t queue );

 
#ifdef __cplusplus
}
#endif

#undef PRECISION_c
#endif /* MAGMA_tally3SPARSE_C_H */
