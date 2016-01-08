/*
 -- MAGMA_tally4 (version 1.1) --
 Univ. of Tennessee, Knoxville
 Univ. of California, Berkeley
 Univ. of Colorado, Denver
 @date

 @generated from magma_tally4sparse_z.h normal z -> c, Mon May  4 17:11:25 2015
 @author Hartwig Anzt
*/

#ifndef MAGMA_tally4SPARSE_C_H
#define MAGMA_tally4SPARSE_C_H

#include "magma_tally4_types.h"
#include "magma_tally4sparse_types.h"

#define PRECISION_c


#ifdef __cplusplus
extern "C" {
#endif


/* ////////////////////////////////////////////////////////////////////////////
 -- For backwards compatability, map old (1.6.1) to new (1.6.2) function names
*/

#define magma_tally4_c_mtranspose  magma_tally4_cmtranspose
#define magma_tally4_c_mtransfer   magma_tally4_cmtransfer
#define magma_tally4_c_vtransfer   magma_tally4_cmtransfer
#define magma_tally4_c_mconvert    magma_tally4_cmconvert
#define magma_tally4_c_vinit       magma_tally4_cvinit
#define magma_tally4_c_vvisu       magma_tally4_cprint_vector
#define magma_tally4_c_vread       magma_tally4_cvread
#define magma_tally4_c_vspread     magma_tally4_cvspread
#define magma_tally4_c_mvisu       magma_tally4_cprint_matrix
#define magma_tally4_c_mfree       magma_tally4_cmfree
#define magma_tally4_c_vfree       magma_tally4_cmfree
#define write_c_csr_mtx     magma_tally4_cwrite_csr_mtx
#define write_c_csrtomtx    magma_tally4_cwrite_csrtomtx
#define print_c_csr         magma_tally4_cprint_csr_mtx


/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally4_SPARSE Auxiliary functions
*/


magma_tally4_int_t
magma_tally4_cparse_opts( 
    int argc, 
    char** argv, 
    magma_tally4_copts *opts, 
    int *matrices, 
    magma_tally4_queue_t queue );

magma_tally4_int_t 
read_c_csr_from_binary( 
    magma_tally4_int_t* n_row, 
    magma_tally4_int_t* n_col, 
    magma_tally4_int_t* nnz, 
    magma_tally4FloatComplex **val, 
    magma_tally4_index_t **row, 
    magma_tally4_index_t **col,
    const char * filename,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
read_c_csr_from_mtx( 
    magma_tally4_storage_t *type, 
    magma_tally4_location_t *location,
    magma_tally4_int_t* n_row, 
    magma_tally4_int_t* n_col, 
    magma_tally4_int_t* nnz, 
    magma_tally4FloatComplex **val, 
    magma_tally4_index_t **row, 
    magma_tally4_index_t **col, 
    const char *filename,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_c_csr_mtx( 
    magma_tally4_c_matrix *A, 
    const char *filename,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_ccsrset( 
    magma_tally4_int_t m, 
    magma_tally4_int_t n, 
    magma_tally4_index_t *row, 
    magma_tally4_index_t *col, 
    magma_tally4FloatComplex *val,
    magma_tally4_c_matrix *A,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_ccsrget( 
    magma_tally4_c_matrix A,
    magma_tally4_int_t *m, 
    magma_tally4_int_t *n, 
    magma_tally4_index_t **row, 
    magma_tally4_index_t **col, 
    magma_tally4FloatComplex **val,
    magma_tally4_queue_t queue );


magma_tally4_int_t 
magma_tally4_cvset( 
    magma_tally4_int_t m, 
    magma_tally4_int_t n, 
    magma_tally4FloatComplex *val,
    magma_tally4_c_matrix *v,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_cvget( 
    magma_tally4_c_matrix v,
    magma_tally4_int_t *m, 
    magma_tally4_int_t *n, 
    magma_tally4FloatComplex **val,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_cvset_dev( 
    magma_tally4_int_t m, 
    magma_tally4_int_t n, 
    magma_tally4FloatComplex_ptr val,
    magma_tally4_c_matrix *v,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_cvget_dev( 
    magma_tally4_c_matrix v,
    magma_tally4_int_t *m, 
    magma_tally4_int_t *n, 
    magma_tally4FloatComplex_ptr *val,
    magma_tally4_queue_t queue );


magma_tally4_int_t 
magma_tally4_c_csr_mtxsymm( 
    magma_tally4_c_matrix *A, 
    const char *filename,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_c_csr_compressor( 
    magma_tally4FloatComplex ** val, 
    magma_tally4_index_t ** row, 
    magma_tally4_index_t ** col, 
    magma_tally4FloatComplex ** valn, 
    magma_tally4_index_t ** rown, 
    magma_tally4_index_t ** coln, 
    magma_tally4_int_t *n,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_cmcsrcompressor( 
    magma_tally4_c_matrix *A,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_cmcsrcompressor_gpu( 
    magma_tally4_c_matrix *A,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_cvtranspose( 
    magma_tally4_c_matrix x,
    magma_tally4_c_matrix *y,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_c_cucsrtranspose( 
    magma_tally4_c_matrix A, 
    magma_tally4_c_matrix *B,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
c_transpose_csr( 
    magma_tally4_int_t n_rows, 
    magma_tally4_int_t n_cols, 
    magma_tally4_int_t nnz,
    magma_tally4FloatComplex *val, 
    magma_tally4_index_t *row, 
    magma_tally4_index_t *col, 
    magma_tally4_int_t *new_n_rows, 
    magma_tally4_int_t *new_n_cols, 
    magma_tally4_int_t *new_nnz, 
    magma_tally4FloatComplex **new_val, 
    magma_tally4_index_t **new_row, 
    magma_tally4_index_t **new_col,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_ccsrsplit( 
    magma_tally4_int_t bsize,
    magma_tally4_c_matrix A,
    magma_tally4_c_matrix *D,
    magma_tally4_c_matrix *R,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_cmscale( 
    magma_tally4_c_matrix *A, 
    magma_tally4_scale_t scaling,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_cmdiff( 
    magma_tally4_c_matrix A, 
    magma_tally4_c_matrix B, 
 real_Double_t *res,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_cmdiagadd( 
    magma_tally4_c_matrix *A, 
    magma_tally4FloatComplex add,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_cmsort( 
    magma_tally4_c_matrix *A,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_cindexsort(
    magma_tally4_index_t *x, 
    magma_tally4_int_t first,
    magma_tally4_int_t last,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_cdomainoverlap(
    magma_tally4_index_t num_rows,
    magma_tally4_index_t *num_indices,
    magma_tally4_index_t *rowptr,
    magma_tally4_index_t *colidx,
    magma_tally4_index_t *x,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_csymbilu( 
    magma_tally4_c_matrix *A, 
    magma_tally4_int_t levels,
    magma_tally4_c_matrix *L,
    magma_tally4_c_matrix *U,
    magma_tally4_queue_t queue );


magma_tally4_int_t 
magma_tally4_cwrite_csr_mtx( 
    magma_tally4_c_matrix A,
    magma_tally4_order_t MajorType,
 const char *filename,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_cwrite_csrtomtx( 
    magma_tally4_c_matrix A,
    const char *filename,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_cprint_csr( 
    magma_tally4_int_t n_row, 
    magma_tally4_int_t n_col, 
    magma_tally4_int_t nnz, 
    magma_tally4FloatComplex **val, 
    magma_tally4_index_t **row, 
    magma_tally4_index_t **col,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_cprint_csr_mtx( 
    magma_tally4_int_t n_row, 
    magma_tally4_int_t n_col, 
    magma_tally4_int_t nnz, 
    magma_tally4FloatComplex **val, 
    magma_tally4_index_t **row, 
    magma_tally4_index_t **col, 
    magma_tally4_order_t MajorType,
    magma_tally4_queue_t queue );


magma_tally4_int_t 
magma_tally4_cmtranspose(
    magma_tally4_c_matrix A, 
    magma_tally4_c_matrix *B,
    magma_tally4_queue_t queue );


magma_tally4_int_t 
magma_tally4_cmtransfer(
    magma_tally4_c_matrix A, 
    magma_tally4_c_matrix *B, 
    magma_tally4_location_t src, 
    magma_tally4_location_t dst,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_cmconvert(
    magma_tally4_c_matrix A, 
    magma_tally4_c_matrix *B, 
    magma_tally4_storage_t old_format, 
    magma_tally4_storage_t new_format,
    magma_tally4_queue_t queue );


magma_tally4_int_t
magma_tally4_cvinit(
    magma_tally4_c_matrix *x, 
    magma_tally4_location_t memory_location,
    magma_tally4_int_t num_rows, 
    magma_tally4_int_t num_cols,
    magma_tally4FloatComplex values,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_cprint_vector(
    magma_tally4_c_matrix x, 
    magma_tally4_int_t offset, 
    magma_tally4_int_t displaylength,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_cvread(
    magma_tally4_c_matrix *x, 
    magma_tally4_int_t length,
    char * filename,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_cvspread(
    magma_tally4_c_matrix *x, 
    const char * filename,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_cprint_matrix(
    magma_tally4_c_matrix A,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_cdiameter(
    magma_tally4_c_matrix *A,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_crowentries(
    magma_tally4_c_matrix *A,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_cmfree(
    magma_tally4_c_matrix *A,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_cresidual(
    magma_tally4_c_matrix A, 
    magma_tally4_c_matrix b, 
    magma_tally4_c_matrix x, 
    float *res,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_cresidualvec(
    magma_tally4_c_matrix A,
    magma_tally4_c_matrix b,
    magma_tally4_c_matrix x,
    magma_tally4_c_matrix *r,
    float *res,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_cmgenerator(
    magma_tally4_int_t n,
    magma_tally4_int_t offdiags,
    magma_tally4_index_t *diag_offset,
    magma_tally4FloatComplex *diag_vals,
    magma_tally4_c_matrix *A,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_cm_27stencil(
    magma_tally4_int_t n,
    magma_tally4_c_matrix *A,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_cm_5stencil(
    magma_tally4_int_t n,
    magma_tally4_c_matrix *A,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_csolverinfo(
    magma_tally4_c_solver_par *solver_par, 
    magma_tally4_c_preconditioner *precond_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_csolverinfo_init(
    magma_tally4_c_solver_par *solver_par, 
    magma_tally4_c_preconditioner *precond,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_ceigensolverinfo_init(
    magma_tally4_c_solver_par *solver_par,
    magma_tally4_queue_t queue );


magma_tally4_int_t
magma_tally4_csolverinfo_free(
    magma_tally4_c_solver_par *solver_par, 
    magma_tally4_c_preconditioner *precond,
    magma_tally4_queue_t queue );



/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally4_SPARSE iterative incomplete factorizations
*/


magma_tally4_int_t
magma_tally4_citerilusetup( 
    magma_tally4_c_matrix A, 
    magma_tally4_c_matrix b,                                 
    magma_tally4_c_preconditioner *precond,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_citericsetup( 
    magma_tally4_c_matrix A, 
    magma_tally4_c_matrix b, 
    magma_tally4_c_preconditioner *precond,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_citericupdate( 
    magma_tally4_c_matrix A, 
    magma_tally4_c_preconditioner *precond, 
    magma_tally4_int_t updates,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_capplyiteric_l( 
    magma_tally4_c_matrix b, 
    magma_tally4_c_matrix *x, 
    magma_tally4_c_preconditioner *precond,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_capplyiteric_r( 
    magma_tally4_c_matrix b, 
    magma_tally4_c_matrix *x, 
    magma_tally4_c_preconditioner *precond,
    magma_tally4_queue_t queue );


magma_tally4_int_t
magma_tally4_citerilu_csr( 
    magma_tally4_c_matrix A,
    magma_tally4_c_matrix L,
    magma_tally4_c_matrix U,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_citeric_csr( 
    magma_tally4_c_matrix A,
    magma_tally4_c_matrix A_CSR,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_cfrobenius( 
    magma_tally4_c_matrix A, 
    magma_tally4_c_matrix B, 
    real_Double_t *res,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_cnonlinres(   
    magma_tally4_c_matrix A, 
    magma_tally4_c_matrix L,
    magma_tally4_c_matrix U, 
    magma_tally4_c_matrix *LU, 
    real_Double_t *res,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_cilures(   
    magma_tally4_c_matrix A, 
    magma_tally4_c_matrix L,
    magma_tally4_c_matrix U, 
    magma_tally4_c_matrix *LU, 
    real_Double_t *res,
    real_Double_t *nonlinres,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_cicres(       
    magma_tally4_c_matrix A, 
    magma_tally4_c_matrix C,
    magma_tally4_c_matrix CT, 
    magma_tally4_c_matrix *LU, 
    real_Double_t *res,
    real_Double_t *nonlinres,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_cinitguess( 
    magma_tally4_c_matrix A, 
    magma_tally4_c_matrix *L, 
    magma_tally4_c_matrix *U,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_cinitrecursiveLU( 
    magma_tally4_c_matrix A, 
    magma_tally4_c_matrix *B,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_cmLdiagadd( 
    magma_tally4_c_matrix *L,
    magma_tally4_queue_t queue );




/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally4_SPARSE function definitions / Data on CPU
*/




/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally4_SPARSE function definitions / Data on CPU / Multi-GPU
*/

/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally4_SPARSE solvers (Data on GPU)
*/

magma_tally4_int_t 
magma_tally4_ccg(
    magma_tally4_c_matrix A, magma_tally4_c_matrix b, 
    magma_tally4_c_matrix *x, magma_tally4_c_solver_par *solver_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_ccg_res(
    magma_tally4_c_matrix A, magma_tally4_c_matrix b, 
    magma_tally4_c_matrix *x, magma_tally4_c_solver_par *solver_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_ccg_merge(
    magma_tally4_c_matrix A, magma_tally4_c_matrix b, 
    magma_tally4_c_matrix *x, magma_tally4_c_solver_par *solver_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_cgmres(
    magma_tally4_c_matrix A, magma_tally4_c_matrix b, 
    magma_tally4_c_matrix *x, magma_tally4_c_solver_par *solver_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_cbicgstab(
    magma_tally4_c_matrix A, magma_tally4_c_matrix b, magma_tally4_c_matrix *x, 
    magma_tally4_c_solver_par *solver_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_cbicgstab_merge(
    magma_tally4_c_matrix A, magma_tally4_c_matrix b, 
    magma_tally4_c_matrix *x, magma_tally4_c_solver_par *solver_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_cbicgstab_merge2(
    magma_tally4_c_matrix A, magma_tally4_c_matrix b, 
    magma_tally4_c_matrix *x, magma_tally4_c_solver_par *solver_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_cpcg(
    magma_tally4_c_matrix A, magma_tally4_c_matrix b, 
    magma_tally4_c_matrix *x, magma_tally4_c_solver_par *solver_par, 
    magma_tally4_c_preconditioner *precond_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_cbpcg(
    magma_tally4_c_matrix A, magma_tally4_c_matrix b, 
    magma_tally4_c_matrix *x, magma_tally4_c_solver_par *solver_par, 
    magma_tally4_c_preconditioner *precond_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_cpbicgstab(
    magma_tally4_c_matrix A, magma_tally4_c_matrix b, 
    magma_tally4_c_matrix *x, magma_tally4_c_solver_par *solver_par, 
    magma_tally4_c_preconditioner *precond_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_cpgmres(
    magma_tally4_c_matrix A, magma_tally4_c_matrix b, 
    magma_tally4_c_matrix *x, magma_tally4_c_solver_par *solver_par, 
    magma_tally4_c_preconditioner *precond_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_cfgmres(
    magma_tally4_c_matrix A, magma_tally4_c_matrix b, 
    magma_tally4_c_matrix *x, magma_tally4_c_solver_par *solver_par, 
    magma_tally4_c_preconditioner *precond_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_cjacobi(
    magma_tally4_c_matrix A, magma_tally4_c_matrix b, 
    magma_tally4_c_matrix *x, magma_tally4_c_solver_par *solver_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_cjacobidomainoverlap(
    magma_tally4_c_matrix A, 
    magma_tally4_c_matrix b, 
    magma_tally4_c_matrix *x,  
    magma_tally4_c_solver_par *solver_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_cbaiter(
    magma_tally4_c_matrix A, magma_tally4_c_matrix b, 
    magma_tally4_c_matrix *x, magma_tally4_c_solver_par *solver_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_citerref(
    magma_tally4_c_matrix A, magma_tally4_c_matrix b, 
    magma_tally4_c_matrix *x, magma_tally4_c_solver_par *solver_par, 
    magma_tally4_c_preconditioner *precond_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_cilu(
    magma_tally4_c_matrix A, magma_tally4_c_matrix b, 
    magma_tally4_c_matrix *x, magma_tally4_c_solver_par *solver_par, 
    magma_tally4_int_t *ipiv,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_cbcsrlu(
    magma_tally4_c_matrix A, magma_tally4_c_matrix b, 
    magma_tally4_c_matrix *x, magma_tally4_c_solver_par *solver_par,
    magma_tally4_queue_t queue );


magma_tally4_int_t
magma_tally4_cbcsrlutrf(
    magma_tally4_c_matrix A, 
    magma_tally4_c_matrix *M,
    magma_tally4_int_t *ipiv, 
    magma_tally4_int_t version,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_cbcsrlusv(
    magma_tally4_c_matrix A, magma_tally4_c_matrix b, 
    magma_tally4_c_matrix *x, magma_tally4_c_solver_par *solver_par, 
    magma_tally4_int_t *ipiv,
    magma_tally4_queue_t queue );



magma_tally4_int_t
magma_tally4_cilucg(
    magma_tally4_c_matrix A, magma_tally4_c_matrix b, 
    magma_tally4_c_matrix *x, magma_tally4_c_solver_par *solver_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_cilugmres(
    magma_tally4_c_matrix A, magma_tally4_c_matrix b, 
    magma_tally4_c_matrix *x, magma_tally4_c_solver_par *solver_par,
    magma_tally4_queue_t queue ); 


magma_tally4_int_t
magma_tally4_clobpcg_shift(
    magma_tally4_int_t num_rows,
    magma_tally4_int_t num_vecs, 
    magma_tally4_int_t shift,
    magma_tally4FloatComplex_ptr x,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_clobpcg_res(
    magma_tally4_int_t num_rows,
    magma_tally4_int_t num_vecs, 
    float *evalues, 
    magma_tally4FloatComplex_ptr X,
    magma_tally4FloatComplex_ptr R, 
    float *res,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_clobpcg_maxpy(
    magma_tally4_int_t num_rows,
    magma_tally4_int_t num_vecs, 
    magma_tally4FloatComplex_ptr X,
    magma_tally4FloatComplex_ptr Y,
    magma_tally4_queue_t queue );


/*/////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally4_SPARSE eigensolvers (Data on GPU)
*/
magma_tally4_int_t
magma_tally4_clobpcg(
    magma_tally4_c_matrix A, 
    magma_tally4_c_solver_par *solver_par,
    magma_tally4_c_preconditioner *precond_par, 
    magma_tally4_queue_t queue );




/*/////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally4_SPARSE preconditioners (Data on GPU)
*/

magma_tally4_int_t
magma_tally4_cjacobisetup(
    magma_tally4_c_matrix A, 
    magma_tally4_c_matrix b, 
    magma_tally4_c_matrix *M, 
    magma_tally4_c_matrix *c,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_cjacobisetup_matrix(
    magma_tally4_c_matrix A, 
    magma_tally4_c_matrix *M, 
    magma_tally4_c_matrix *d,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_cjacobisetup_vector(
    magma_tally4_c_matrix b, 
    magma_tally4_c_matrix d, 
    magma_tally4_c_matrix *c,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_cjacobiiter(
    magma_tally4_c_matrix M, 
    magma_tally4_c_matrix c, 
    magma_tally4_c_matrix *x, 
    magma_tally4_c_solver_par *solver_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_cjacobiiter_precond( 
    magma_tally4_c_matrix M, 
    magma_tally4_c_matrix *x, 
    magma_tally4_c_solver_par *solver_par, 
    magma_tally4_c_preconditioner *precond,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_cjacobiiter_sys(
    magma_tally4_c_matrix A, 
    magma_tally4_c_matrix b, 
    magma_tally4_c_matrix d, 
    magma_tally4_c_matrix t, 
    magma_tally4_c_matrix *x,  
    magma_tally4_c_solver_par *solver_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_cpastixsetup(
    magma_tally4_c_matrix A, magma_tally4_c_matrix b,
    magma_tally4_c_preconditioner *precond,
    magma_tally4_queue_t queue );


magma_tally4_int_t
magma_tally4_capplypastix(
    magma_tally4_c_matrix b, magma_tally4_c_matrix *x, 
    magma_tally4_c_preconditioner *precond,
    magma_tally4_queue_t queue );


// custom preconditioner
magma_tally4_int_t
magma_tally4_capplycustomprecond_l(
    magma_tally4_c_matrix b, 
    magma_tally4_c_matrix *x, 
    magma_tally4_c_preconditioner *precond,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_capplycustomprecond_r(
    magma_tally4_c_matrix b, 
    magma_tally4_c_matrix *x, 
    magma_tally4_c_preconditioner *precond,
    magma_tally4_queue_t queue );


// CUSPARSE preconditioner

magma_tally4_int_t
magma_tally4_ccuilusetup(
    magma_tally4_c_matrix A, magma_tally4_c_preconditioner *precond,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_capplycuilu_l(
    magma_tally4_c_matrix b, magma_tally4_c_matrix *x, 
    magma_tally4_c_preconditioner *precond,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_capplycuilu_r(
    magma_tally4_c_matrix b, magma_tally4_c_matrix *x, 
    magma_tally4_c_preconditioner *precond,
    magma_tally4_queue_t queue );


magma_tally4_int_t
magma_tally4_ccuiccsetup(
    magma_tally4_c_matrix A, magma_tally4_c_preconditioner *precond,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_capplycuicc_l(
    magma_tally4_c_matrix b, magma_tally4_c_matrix *x, 
    magma_tally4_c_preconditioner *precond,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_capplycuicc_r(
    magma_tally4_c_matrix b, magma_tally4_c_matrix *x, 
    magma_tally4_c_preconditioner *precond,
    magma_tally4_queue_t queue );


magma_tally4_int_t
magma_tally4_ccumilusetup(
    magma_tally4_c_matrix A, 
    magma_tally4_c_preconditioner *precond,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_ccumilugeneratesolverinfo(
    magma_tally4_c_preconditioner *precond,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_capplycumilu_l(
    magma_tally4_c_matrix b, 
    magma_tally4_c_matrix *x, 
    magma_tally4_c_preconditioner *precond,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_capplycumilu_r(
    magma_tally4_c_matrix b, 
    magma_tally4_c_matrix *x, 
    magma_tally4_c_preconditioner *precond,
    magma_tally4_queue_t queue );


magma_tally4_int_t
magma_tally4_ccumiccsetup(
    magma_tally4_c_matrix A, 
    magma_tally4_c_preconditioner *precond,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_ccumicgeneratesolverinfo(
    magma_tally4_c_preconditioner *precond,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_capplycumicc_l(
    magma_tally4_c_matrix b, 
    magma_tally4_c_matrix *x, 
    magma_tally4_c_preconditioner *precond,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_capplycumicc_r(
    magma_tally4_c_matrix b, 
    magma_tally4_c_matrix *x, 
    magma_tally4_c_preconditioner *precond,
    magma_tally4_queue_t queue );


// block-asynchronous iteration

magma_tally4_int_t
magma_tally4_cbajac_csr(
    magma_tally4_int_t localiters,
    magma_tally4_c_matrix D,
    magma_tally4_c_matrix R,
    magma_tally4_c_matrix b,
    magma_tally4_c_matrix *x,
    magma_tally4_queue_t queue );

/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally4_SPARSE utility function definitions
*/

magma_tally4_int_t
magma_tally4_c_spmv(
    magma_tally4FloatComplex alpha, 
    magma_tally4_c_matrix A, 
    magma_tally4_c_matrix x, 
    magma_tally4FloatComplex beta, 
    magma_tally4_c_matrix y,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_ccustomspmv(
    magma_tally4FloatComplex alpha, 
    magma_tally4_c_matrix x, 
    magma_tally4FloatComplex beta, 
    magma_tally4_c_matrix y,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_c_spmv_shift(
    magma_tally4FloatComplex alpha, 
    magma_tally4_c_matrix A, 
    magma_tally4FloatComplex lambda,
    magma_tally4_c_matrix x, 
    magma_tally4FloatComplex beta, 
    magma_tally4_int_t offset, 
    magma_tally4_int_t blocksize,
    magma_tally4Index_ptr dadd_vecs, 
    magma_tally4_c_matrix y,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_ccuspmm(
    magma_tally4_c_matrix A, 
    magma_tally4_c_matrix B, 
    magma_tally4_c_matrix *AB,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_c_spmm(
    magma_tally4FloatComplex alpha, 
    magma_tally4_c_matrix A,
    magma_tally4_c_matrix B,
    magma_tally4_c_matrix *C,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_csymbilu( 
    magma_tally4_c_matrix *A, 
    magma_tally4_int_t levels, 
    magma_tally4_c_matrix *L, 
    magma_tally4_c_matrix *U,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_ccuspaxpy(
    magma_tally4FloatComplex_ptr alpha, magma_tally4_c_matrix A, 
    magma_tally4FloatComplex_ptr beta, magma_tally4_c_matrix B, 
    magma_tally4_c_matrix *AB,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_c_precond(
    magma_tally4_c_matrix A, 
    magma_tally4_c_matrix b, magma_tally4_c_matrix *x,
    magma_tally4_c_preconditioner *precond,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_c_solver(
    magma_tally4_c_matrix A, magma_tally4_c_matrix b, 
    magma_tally4_c_matrix *x, magma_tally4_copts *zopts,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_c_precondsetup(
    magma_tally4_c_matrix A, magma_tally4_c_matrix b, 
    magma_tally4_c_preconditioner *precond,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_c_applyprecond(
    magma_tally4_c_matrix A, magma_tally4_c_matrix b, 
    magma_tally4_c_matrix *x, magma_tally4_c_preconditioner *precond,
    magma_tally4_queue_t queue );


magma_tally4_int_t
magma_tally4_c_applyprecond_left(
    magma_tally4_c_matrix A, magma_tally4_c_matrix b, 
    magma_tally4_c_matrix *x, magma_tally4_c_preconditioner *precond,
    magma_tally4_queue_t queue );


magma_tally4_int_t
magma_tally4_c_applyprecond_right(
    magma_tally4_c_matrix A, magma_tally4_c_matrix b, 
    magma_tally4_c_matrix *x, magma_tally4_c_preconditioner *precond,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_c_initP2P(
    magma_tally4_int_t *bandwidth_benchmark,
    magma_tally4_int_t *num_gpus,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_ccompact(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4FloatComplex_ptr dA, magma_tally4_int_t ldda,
    float *dnorms, float tol, 
    magma_tally4_int_t *activeMask, magma_tally4_int_t *cBlockSize,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_ccompactActive(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4FloatComplex_ptr dA, magma_tally4_int_t ldda, 
    magma_tally4_int_t *active,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_cmlumerge(    
    magma_tally4_c_matrix L, 
    magma_tally4_c_matrix U,
    magma_tally4_c_matrix *A, 
    magma_tally4_queue_t queue );


/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally4_SPARSE BLAS function definitions
*/
magma_tally4_int_t 
magma_tally4_cgecsrmv(
    magma_tally4_trans_t transA,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4FloatComplex alpha,
    magma_tally4FloatComplex_ptr dval,
    magma_tally4Index_ptr drowptr,
    magma_tally4Index_ptr dcolind,
    magma_tally4FloatComplex_ptr dx,
    magma_tally4FloatComplex beta,
    magma_tally4FloatComplex_ptr dy,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_cgecsrmv_shift(
    magma_tally4_trans_t transA,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4FloatComplex alpha,
    magma_tally4FloatComplex lambda,
    magma_tally4FloatComplex_ptr dval,
    magma_tally4Index_ptr drowptr,
    magma_tally4Index_ptr dcolind,
    magma_tally4FloatComplex_ptr dx,
    magma_tally4FloatComplex beta,
    int offset,
    int blocksize,
    magma_tally4Index_ptr dadd_rows,
    magma_tally4FloatComplex_ptr dy,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_cmgecsrmv(
    magma_tally4_trans_t transA,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4_int_t num_vecs,
    magma_tally4FloatComplex alpha,
    magma_tally4FloatComplex_ptr dval,
    magma_tally4Index_ptr drowptr,
    magma_tally4Index_ptr dcolind,
    magma_tally4FloatComplex_ptr dx,
    magma_tally4FloatComplex beta,
    magma_tally4FloatComplex_ptr dy,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_cgeellmv(
    magma_tally4_trans_t transA,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4_int_t nnz_per_row,
    magma_tally4FloatComplex alpha,
    magma_tally4FloatComplex_ptr dval,
    magma_tally4Index_ptr dcolind,
    magma_tally4FloatComplex_ptr dx,
    magma_tally4FloatComplex beta,
    magma_tally4FloatComplex_ptr dy,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_cgeellmv_shift(
    magma_tally4_trans_t transA,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4_int_t nnz_per_row,
    magma_tally4FloatComplex alpha,
    magma_tally4FloatComplex lambda,
    magma_tally4FloatComplex_ptr dval,
    magma_tally4Index_ptr dcolind,
    magma_tally4FloatComplex_ptr dx,
    magma_tally4FloatComplex beta,
    int offset,
    int blocksize,
    magma_tally4Index_ptr dadd_rows,
    magma_tally4FloatComplex_ptr dy,
    magma_tally4_queue_t queue );


magma_tally4_int_t 
magma_tally4_cmgeellmv(
    magma_tally4_trans_t transA,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4_int_t num_vecs,
    magma_tally4_int_t nnz_per_row,
    magma_tally4FloatComplex alpha,
    magma_tally4FloatComplex_ptr dval,
    magma_tally4Index_ptr dcolind,
    magma_tally4FloatComplex_ptr dx,
    magma_tally4FloatComplex beta,
    magma_tally4FloatComplex_ptr dy,
    magma_tally4_queue_t queue );


magma_tally4_int_t 
magma_tally4_cgeelltmv(
    magma_tally4_trans_t transA,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4_int_t nnz_per_row,
    magma_tally4FloatComplex alpha,
    magma_tally4FloatComplex_ptr dval,
    magma_tally4Index_ptr dcolind,
    magma_tally4FloatComplex_ptr dx,
    magma_tally4FloatComplex beta,
    magma_tally4FloatComplex_ptr dy,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_cgeelltmv_shift(
    magma_tally4_trans_t transA,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4_int_t nnz_per_row,
    magma_tally4FloatComplex alpha,
    magma_tally4FloatComplex lambda,
    magma_tally4FloatComplex_ptr dval,
    magma_tally4Index_ptr dcolind,
    magma_tally4FloatComplex_ptr dx,
    magma_tally4FloatComplex beta,
    int offset,
    int blocksize,
    magma_tally4Index_ptr dadd_rows,
    magma_tally4FloatComplex_ptr dy,
    magma_tally4_queue_t queue );


magma_tally4_int_t 
magma_tally4_cmgeelltmv(
    magma_tally4_trans_t transA,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4_int_t num_vecs,
    magma_tally4_int_t nnz_per_row,
    magma_tally4FloatComplex alpha,
    magma_tally4FloatComplex_ptr dval,
    magma_tally4Index_ptr dcolind,
    magma_tally4FloatComplex_ptr dx,
    magma_tally4FloatComplex beta,
    magma_tally4FloatComplex_ptr dy,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_cgeellrtmv(
    magma_tally4_trans_t transA,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4_int_t nnz_per_row,
    magma_tally4FloatComplex alpha,
    magma_tally4FloatComplex_ptr dval,
    magma_tally4Index_ptr dcolind,
    magma_tally4Index_ptr drowlength,
    magma_tally4FloatComplex_ptr dx,
    magma_tally4FloatComplex beta,
    magma_tally4FloatComplex_ptr dy,
    magma_tally4_int_t num_threads,
    magma_tally4_int_t threads_per_row,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_cgesellcmv(
    magma_tally4_trans_t transA,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4_int_t blocksize,
    magma_tally4_int_t slices,
    magma_tally4_int_t alignment,
    magma_tally4FloatComplex alpha,
    magma_tally4FloatComplex_ptr dval,
    magma_tally4Index_ptr dcolind,
    magma_tally4Index_ptr drowptr,
    magma_tally4FloatComplex_ptr dx,
    magma_tally4FloatComplex beta,
    magma_tally4FloatComplex_ptr dy,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_cgesellpmv(
    magma_tally4_trans_t transA,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4_int_t blocksize,
    magma_tally4_int_t slices,
    magma_tally4_int_t alignment,
    magma_tally4FloatComplex alpha,
    magma_tally4FloatComplex_ptr dval,
    magma_tally4Index_ptr dcolind,
    magma_tally4Index_ptr drowptr,
    magma_tally4FloatComplex_ptr dx,
    magma_tally4FloatComplex beta,
    magma_tally4FloatComplex_ptr dy,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_cmgesellpmv(
    magma_tally4_trans_t transA,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4_int_t num_vecs,
    magma_tally4_int_t blocksize,
    magma_tally4_int_t slices,
    magma_tally4_int_t alignment,
    magma_tally4FloatComplex alpha,
    magma_tally4FloatComplex_ptr dval,
    magma_tally4Index_ptr dcolind,
    magma_tally4Index_ptr drowptr,
    magma_tally4FloatComplex_ptr dx,
    magma_tally4FloatComplex beta,
    magma_tally4FloatComplex_ptr dy,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_cmgesellpmv_blocked(
    magma_tally4_trans_t transA,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4_int_t num_vecs,
    magma_tally4_int_t blocksize,
    magma_tally4_int_t slices,
    magma_tally4_int_t alignment,
    magma_tally4FloatComplex alpha,
    magma_tally4FloatComplex_ptr dval,
    magma_tally4Index_ptr dcolind,
    magma_tally4Index_ptr drowptr,
    magma_tally4FloatComplex_ptr dx,
    magma_tally4FloatComplex beta,
    magma_tally4FloatComplex_ptr dy,
    magma_tally4_queue_t queue );


magma_tally4_int_t
magma_tally4_cmergedgs(
    magma_tally4_int_t n, 
    magma_tally4_int_t ldh,
    magma_tally4_int_t k, 
    magma_tally4FloatComplex_ptr dv, 
    magma_tally4FloatComplex_ptr dr,
    magma_tally4FloatComplex_ptr dskp,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_ccopyscale(    
    int n, 
    int k,
    magma_tally4FloatComplex_ptr dr, 
    magma_tally4FloatComplex_ptr dv,
    magma_tally4FloatComplex_ptr dskp,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_scnrm2scale(    
    int m, 
    magma_tally4FloatComplex_ptr dr,    
    int lddr, 
    magma_tally4FloatComplex *drnorm,
    magma_tally4_queue_t queue );


magma_tally4_int_t
magma_tally4_cjacobisetup_vector_gpu(
    int num_rows, 
    magma_tally4_c_matrix b, 
    magma_tally4_c_matrix d, 
    magma_tally4_c_matrix c,
    magma_tally4_c_matrix *x,
    magma_tally4_queue_t queue );


magma_tally4_int_t
magma_tally4_cjacobi_diagscal(    
    int num_rows, 
    magma_tally4_c_matrix d, 
    magma_tally4_c_matrix b, 
    magma_tally4_c_matrix *c,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_cjacobiupdate(
    magma_tally4_c_matrix t, 
    magma_tally4_c_matrix b, 
    magma_tally4_c_matrix d, 
    magma_tally4_c_matrix *x,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_cjacobispmvupdate(
    magma_tally4_int_t maxiter,
    magma_tally4_c_matrix A, 
    magma_tally4_c_matrix t, 
    magma_tally4_c_matrix b, 
    magma_tally4_c_matrix d, 
    magma_tally4_c_matrix *x,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_cjacobispmvupdate_bw(
    magma_tally4_int_t maxiter,
    magma_tally4_c_matrix A, 
    magma_tally4_c_matrix t, 
    magma_tally4_c_matrix b, 
    magma_tally4_c_matrix d, 
    magma_tally4_c_matrix *x,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_cjacobispmvupdateselect(
    magma_tally4_int_t maxiter,
    magma_tally4_int_t num_updates,
    magma_tally4_index_t *indices,
    magma_tally4_c_matrix A,
    magma_tally4_c_matrix t, 
    magma_tally4_c_matrix b, 
    magma_tally4_c_matrix d, 
    magma_tally4_c_matrix tmp, 
    magma_tally4_c_matrix *x,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_cjacobisetup_diagscal(
    magma_tally4_c_matrix A, magma_tally4_c_matrix *d,
    magma_tally4_queue_t queue );


magma_tally4_int_t
magma_tally4_cbicgmerge1(    
    int n, 
    magma_tally4FloatComplex_ptr dskp,
    magma_tally4FloatComplex_ptr dv, 
    magma_tally4FloatComplex_ptr dr, 
    magma_tally4FloatComplex_ptr dp,
    magma_tally4_queue_t queue );


magma_tally4_int_t
magma_tally4_cbicgmerge2(
    int n, 
    magma_tally4FloatComplex_ptr dskp, 
    magma_tally4FloatComplex_ptr dr,
    magma_tally4FloatComplex_ptr dv, 
    magma_tally4FloatComplex_ptr ds,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_cbicgmerge3(
    int n, 
    magma_tally4FloatComplex_ptr dskp, 
    magma_tally4FloatComplex_ptr dp,
    magma_tally4FloatComplex_ptr ds,
    magma_tally4FloatComplex_ptr dt,
    magma_tally4FloatComplex_ptr dx, 
    magma_tally4FloatComplex_ptr dr,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_cbicgmerge4(
    int type, 
    magma_tally4FloatComplex_ptr dskp,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_ccgmerge_spmv1( 
    magma_tally4_c_matrix A,
    magma_tally4FloatComplex_ptr d1,
    magma_tally4FloatComplex_ptr d2,
    magma_tally4FloatComplex_ptr dd,
    magma_tally4FloatComplex_ptr dz,
    magma_tally4FloatComplex_ptr dskp,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_ccgmerge_xrbeta( 
    int n,
    magma_tally4FloatComplex_ptr d1,
    magma_tally4FloatComplex_ptr d2,
    magma_tally4FloatComplex_ptr dx,
    magma_tally4FloatComplex_ptr dr,
    magma_tally4FloatComplex_ptr dd,
    magma_tally4FloatComplex_ptr dz, 
    magma_tally4FloatComplex_ptr dskp,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_cmdotc(
    magma_tally4_int_t n, 
    magma_tally4_int_t k, 
    magma_tally4FloatComplex_ptr dv, 
    magma_tally4FloatComplex_ptr dr,
    magma_tally4FloatComplex_ptr dd1,
    magma_tally4FloatComplex_ptr dd2,
    magma_tally4FloatComplex_ptr dskp,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_cgemvmdot(
    int n, 
    int k, 
    magma_tally4FloatComplex_ptr dv, 
    magma_tally4FloatComplex_ptr dr,
    magma_tally4FloatComplex_ptr dd1,
    magma_tally4FloatComplex_ptr dd2,
    magma_tally4FloatComplex_ptr dskp,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_cbicgmerge_spmv1( 
    magma_tally4_c_matrix A,
    magma_tally4FloatComplex_ptr dd1,
    magma_tally4FloatComplex_ptr dd2,
    magma_tally4FloatComplex_ptr dp,
    magma_tally4FloatComplex_ptr dr,
    magma_tally4FloatComplex_ptr dv,
    magma_tally4FloatComplex_ptr dskp,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_cbicgmerge_spmv2( 
    magma_tally4_c_matrix A,
    magma_tally4FloatComplex_ptr dd1,
    magma_tally4FloatComplex_ptr dd2,
    magma_tally4FloatComplex_ptr ds,
    magma_tally4FloatComplex_ptr dt,
    magma_tally4FloatComplex_ptr dskp,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_cbicgmerge_xrbeta( 
    int n,
    magma_tally4FloatComplex_ptr dd1,
    magma_tally4FloatComplex_ptr dd2,
    magma_tally4FloatComplex_ptr drr,
    magma_tally4FloatComplex_ptr dr,
    magma_tally4FloatComplex_ptr dp,
    magma_tally4FloatComplex_ptr ds,
    magma_tally4FloatComplex_ptr dt,
    magma_tally4FloatComplex_ptr dx, 
    magma_tally4FloatComplex_ptr dskp,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_cbcsrswp(
    magma_tally4_int_t n,
    magma_tally4_int_t size_b, 
    magma_tally4_int_t *ipiv,
    magma_tally4FloatComplex_ptr dx,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_cbcsrtrsv(
    magma_tally4_uplo_t uplo,
    magma_tally4_int_t r_blocks,
    magma_tally4_int_t c_blocks,
    magma_tally4_int_t size_b, 
    magma_tally4FloatComplex_ptr dA,
    magma_tally4_index_t *blockinfo, 
    magma_tally4FloatComplex_ptr dx,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_cbcsrvalcpy(
    magma_tally4_int_t size_b, 
    magma_tally4_int_t num_blocks, 
    magma_tally4_int_t num_zero_blocks, 
    magma_tally4FloatComplex_ptr *dAval, 
    magma_tally4FloatComplex_ptr *dBval,
    magma_tally4FloatComplex_ptr *dBval2,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_cbcsrluegemm(
    magma_tally4_int_t size_b, 
    magma_tally4_int_t num_block_rows,
    magma_tally4_int_t kblocks,
    magma_tally4FloatComplex_ptr *dA, 
    magma_tally4FloatComplex_ptr *dB, 
    magma_tally4FloatComplex_ptr *dC,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_cbcsrlupivloc(
    magma_tally4_int_t size_b, 
    magma_tally4_int_t kblocks,
    magma_tally4FloatComplex_ptr *dA, 
    magma_tally4_int_t *ipiv,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_cbcsrblockinfo5(
    magma_tally4_int_t lustep,
    magma_tally4_int_t num_blocks, 
    magma_tally4_int_t c_blocks, 
    magma_tally4_int_t size_b,
    magma_tally4_index_t *blockinfo,
    magma_tally4FloatComplex_ptr dval,
    magma_tally4FloatComplex_ptr *AII,
    magma_tally4_queue_t queue );

 
#ifdef __cplusplus
}
#endif

#undef PRECISION_c
#endif /* MAGMA_tally4SPARSE_C_H */
