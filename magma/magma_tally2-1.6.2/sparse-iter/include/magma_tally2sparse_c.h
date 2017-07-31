/*
 -- MAGMA_tally2 (version 1.1) --
 Univ. of Tennessee, Knoxville
 Univ. of California, Berkeley
 Univ. of Colorado, Denver
 @date

 @generated from magma_tally2sparse_z.h normal z -> c, Mon May  4 17:11:25 2015
 @author Hartwig Anzt
*/

#ifndef MAGMA_tally2SPARSE_C_H
#define MAGMA_tally2SPARSE_C_H

#include "magma_tally2_types.h"
#include "magma_tally2sparse_types.h"

#define PRECISION_c


#ifdef __cplusplus
extern "C" {
#endif


/* ////////////////////////////////////////////////////////////////////////////
 -- For backwards compatability, map old (1.6.1) to new (1.6.2) function names
*/

#define magma_tally2_c_mtranspose  magma_tally2_cmtranspose
#define magma_tally2_c_mtransfer   magma_tally2_cmtransfer
#define magma_tally2_c_vtransfer   magma_tally2_cmtransfer
#define magma_tally2_c_mconvert    magma_tally2_cmconvert
#define magma_tally2_c_vinit       magma_tally2_cvinit
#define magma_tally2_c_vvisu       magma_tally2_cprint_vector
#define magma_tally2_c_vread       magma_tally2_cvread
#define magma_tally2_c_vspread     magma_tally2_cvspread
#define magma_tally2_c_mvisu       magma_tally2_cprint_matrix
#define magma_tally2_c_mfree       magma_tally2_cmfree
#define magma_tally2_c_vfree       magma_tally2_cmfree
#define write_c_csr_mtx     magma_tally2_cwrite_csr_mtx
#define write_c_csrtomtx    magma_tally2_cwrite_csrtomtx
#define print_c_csr         magma_tally2_cprint_csr_mtx


/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally2_SPARSE Auxiliary functions
*/


magma_tally2_int_t
magma_tally2_cparse_opts( 
    int argc, 
    char** argv, 
    magma_tally2_copts *opts, 
    int *matrices, 
    magma_tally2_queue_t queue );

magma_tally2_int_t 
read_c_csr_from_binary( 
    magma_tally2_int_t* n_row, 
    magma_tally2_int_t* n_col, 
    magma_tally2_int_t* nnz, 
    magma_tally2FloatComplex **val, 
    magma_tally2_index_t **row, 
    magma_tally2_index_t **col,
    const char * filename,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
read_c_csr_from_mtx( 
    magma_tally2_storage_t *type, 
    magma_tally2_location_t *location,
    magma_tally2_int_t* n_row, 
    magma_tally2_int_t* n_col, 
    magma_tally2_int_t* nnz, 
    magma_tally2FloatComplex **val, 
    magma_tally2_index_t **row, 
    magma_tally2_index_t **col, 
    const char *filename,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_c_csr_mtx( 
    magma_tally2_c_matrix *A, 
    const char *filename,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_ccsrset( 
    magma_tally2_int_t m, 
    magma_tally2_int_t n, 
    magma_tally2_index_t *row, 
    magma_tally2_index_t *col, 
    magma_tally2FloatComplex *val,
    magma_tally2_c_matrix *A,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_ccsrget( 
    magma_tally2_c_matrix A,
    magma_tally2_int_t *m, 
    magma_tally2_int_t *n, 
    magma_tally2_index_t **row, 
    magma_tally2_index_t **col, 
    magma_tally2FloatComplex **val,
    magma_tally2_queue_t queue );


magma_tally2_int_t 
magma_tally2_cvset( 
    magma_tally2_int_t m, 
    magma_tally2_int_t n, 
    magma_tally2FloatComplex *val,
    magma_tally2_c_matrix *v,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_cvget( 
    magma_tally2_c_matrix v,
    magma_tally2_int_t *m, 
    magma_tally2_int_t *n, 
    magma_tally2FloatComplex **val,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_cvset_dev( 
    magma_tally2_int_t m, 
    magma_tally2_int_t n, 
    magma_tally2FloatComplex_ptr val,
    magma_tally2_c_matrix *v,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_cvget_dev( 
    magma_tally2_c_matrix v,
    magma_tally2_int_t *m, 
    magma_tally2_int_t *n, 
    magma_tally2FloatComplex_ptr *val,
    magma_tally2_queue_t queue );


magma_tally2_int_t 
magma_tally2_c_csr_mtxsymm( 
    magma_tally2_c_matrix *A, 
    const char *filename,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_c_csr_compressor( 
    magma_tally2FloatComplex ** val, 
    magma_tally2_index_t ** row, 
    magma_tally2_index_t ** col, 
    magma_tally2FloatComplex ** valn, 
    magma_tally2_index_t ** rown, 
    magma_tally2_index_t ** coln, 
    magma_tally2_int_t *n,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_cmcsrcompressor( 
    magma_tally2_c_matrix *A,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_cmcsrcompressor_gpu( 
    magma_tally2_c_matrix *A,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_cvtranspose( 
    magma_tally2_c_matrix x,
    magma_tally2_c_matrix *y,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_c_cucsrtranspose( 
    magma_tally2_c_matrix A, 
    magma_tally2_c_matrix *B,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
c_transpose_csr( 
    magma_tally2_int_t n_rows, 
    magma_tally2_int_t n_cols, 
    magma_tally2_int_t nnz,
    magma_tally2FloatComplex *val, 
    magma_tally2_index_t *row, 
    magma_tally2_index_t *col, 
    magma_tally2_int_t *new_n_rows, 
    magma_tally2_int_t *new_n_cols, 
    magma_tally2_int_t *new_nnz, 
    magma_tally2FloatComplex **new_val, 
    magma_tally2_index_t **new_row, 
    magma_tally2_index_t **new_col,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_ccsrsplit( 
    magma_tally2_int_t bsize,
    magma_tally2_c_matrix A,
    magma_tally2_c_matrix *D,
    magma_tally2_c_matrix *R,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_cmscale( 
    magma_tally2_c_matrix *A, 
    magma_tally2_scale_t scaling,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_cmdiff( 
    magma_tally2_c_matrix A, 
    magma_tally2_c_matrix B, 
 real_Double_t *res,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_cmdiagadd( 
    magma_tally2_c_matrix *A, 
    magma_tally2FloatComplex add,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_cmsort( 
    magma_tally2_c_matrix *A,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_cindexsort(
    magma_tally2_index_t *x, 
    magma_tally2_int_t first,
    magma_tally2_int_t last,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_cdomainoverlap(
    magma_tally2_index_t num_rows,
    magma_tally2_index_t *num_indices,
    magma_tally2_index_t *rowptr,
    magma_tally2_index_t *colidx,
    magma_tally2_index_t *x,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_csymbilu( 
    magma_tally2_c_matrix *A, 
    magma_tally2_int_t levels,
    magma_tally2_c_matrix *L,
    magma_tally2_c_matrix *U,
    magma_tally2_queue_t queue );


magma_tally2_int_t 
magma_tally2_cwrite_csr_mtx( 
    magma_tally2_c_matrix A,
    magma_tally2_order_t MajorType,
 const char *filename,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_cwrite_csrtomtx( 
    magma_tally2_c_matrix A,
    const char *filename,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_cprint_csr( 
    magma_tally2_int_t n_row, 
    magma_tally2_int_t n_col, 
    magma_tally2_int_t nnz, 
    magma_tally2FloatComplex **val, 
    magma_tally2_index_t **row, 
    magma_tally2_index_t **col,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_cprint_csr_mtx( 
    magma_tally2_int_t n_row, 
    magma_tally2_int_t n_col, 
    magma_tally2_int_t nnz, 
    magma_tally2FloatComplex **val, 
    magma_tally2_index_t **row, 
    magma_tally2_index_t **col, 
    magma_tally2_order_t MajorType,
    magma_tally2_queue_t queue );


magma_tally2_int_t 
magma_tally2_cmtranspose(
    magma_tally2_c_matrix A, 
    magma_tally2_c_matrix *B,
    magma_tally2_queue_t queue );


magma_tally2_int_t 
magma_tally2_cmtransfer(
    magma_tally2_c_matrix A, 
    magma_tally2_c_matrix *B, 
    magma_tally2_location_t src, 
    magma_tally2_location_t dst,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_cmconvert(
    magma_tally2_c_matrix A, 
    magma_tally2_c_matrix *B, 
    magma_tally2_storage_t old_format, 
    magma_tally2_storage_t new_format,
    magma_tally2_queue_t queue );


magma_tally2_int_t
magma_tally2_cvinit(
    magma_tally2_c_matrix *x, 
    magma_tally2_location_t memory_location,
    magma_tally2_int_t num_rows, 
    magma_tally2_int_t num_cols,
    magma_tally2FloatComplex values,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_cprint_vector(
    magma_tally2_c_matrix x, 
    magma_tally2_int_t offset, 
    magma_tally2_int_t displaylength,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_cvread(
    magma_tally2_c_matrix *x, 
    magma_tally2_int_t length,
    char * filename,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_cvspread(
    magma_tally2_c_matrix *x, 
    const char * filename,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_cprint_matrix(
    magma_tally2_c_matrix A,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_cdiameter(
    magma_tally2_c_matrix *A,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_crowentries(
    magma_tally2_c_matrix *A,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_cmfree(
    magma_tally2_c_matrix *A,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_cresidual(
    magma_tally2_c_matrix A, 
    magma_tally2_c_matrix b, 
    magma_tally2_c_matrix x, 
    float *res,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_cresidualvec(
    magma_tally2_c_matrix A,
    magma_tally2_c_matrix b,
    magma_tally2_c_matrix x,
    magma_tally2_c_matrix *r,
    float *res,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_cmgenerator(
    magma_tally2_int_t n,
    magma_tally2_int_t offdiags,
    magma_tally2_index_t *diag_offset,
    magma_tally2FloatComplex *diag_vals,
    magma_tally2_c_matrix *A,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_cm_27stencil(
    magma_tally2_int_t n,
    magma_tally2_c_matrix *A,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_cm_5stencil(
    magma_tally2_int_t n,
    magma_tally2_c_matrix *A,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_csolverinfo(
    magma_tally2_c_solver_par *solver_par, 
    magma_tally2_c_preconditioner *precond_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_csolverinfo_init(
    magma_tally2_c_solver_par *solver_par, 
    magma_tally2_c_preconditioner *precond,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_ceigensolverinfo_init(
    magma_tally2_c_solver_par *solver_par,
    magma_tally2_queue_t queue );


magma_tally2_int_t
magma_tally2_csolverinfo_free(
    magma_tally2_c_solver_par *solver_par, 
    magma_tally2_c_preconditioner *precond,
    magma_tally2_queue_t queue );



/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally2_SPARSE iterative incomplete factorizations
*/


magma_tally2_int_t
magma_tally2_citerilusetup( 
    magma_tally2_c_matrix A, 
    magma_tally2_c_matrix b,                                 
    magma_tally2_c_preconditioner *precond,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_citericsetup( 
    magma_tally2_c_matrix A, 
    magma_tally2_c_matrix b, 
    magma_tally2_c_preconditioner *precond,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_citericupdate( 
    magma_tally2_c_matrix A, 
    magma_tally2_c_preconditioner *precond, 
    magma_tally2_int_t updates,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_capplyiteric_l( 
    magma_tally2_c_matrix b, 
    magma_tally2_c_matrix *x, 
    magma_tally2_c_preconditioner *precond,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_capplyiteric_r( 
    magma_tally2_c_matrix b, 
    magma_tally2_c_matrix *x, 
    magma_tally2_c_preconditioner *precond,
    magma_tally2_queue_t queue );


magma_tally2_int_t
magma_tally2_citerilu_csr( 
    magma_tally2_c_matrix A,
    magma_tally2_c_matrix L,
    magma_tally2_c_matrix U,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_citeric_csr( 
    magma_tally2_c_matrix A,
    magma_tally2_c_matrix A_CSR,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_cfrobenius( 
    magma_tally2_c_matrix A, 
    magma_tally2_c_matrix B, 
    real_Double_t *res,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_cnonlinres(   
    magma_tally2_c_matrix A, 
    magma_tally2_c_matrix L,
    magma_tally2_c_matrix U, 
    magma_tally2_c_matrix *LU, 
    real_Double_t *res,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_cilures(   
    magma_tally2_c_matrix A, 
    magma_tally2_c_matrix L,
    magma_tally2_c_matrix U, 
    magma_tally2_c_matrix *LU, 
    real_Double_t *res,
    real_Double_t *nonlinres,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_cicres(       
    magma_tally2_c_matrix A, 
    magma_tally2_c_matrix C,
    magma_tally2_c_matrix CT, 
    magma_tally2_c_matrix *LU, 
    real_Double_t *res,
    real_Double_t *nonlinres,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_cinitguess( 
    magma_tally2_c_matrix A, 
    magma_tally2_c_matrix *L, 
    magma_tally2_c_matrix *U,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_cinitrecursiveLU( 
    magma_tally2_c_matrix A, 
    magma_tally2_c_matrix *B,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_cmLdiagadd( 
    magma_tally2_c_matrix *L,
    magma_tally2_queue_t queue );




/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally2_SPARSE function definitions / Data on CPU
*/




/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally2_SPARSE function definitions / Data on CPU / Multi-GPU
*/

/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally2_SPARSE solvers (Data on GPU)
*/

magma_tally2_int_t 
magma_tally2_ccg(
    magma_tally2_c_matrix A, magma_tally2_c_matrix b, 
    magma_tally2_c_matrix *x, magma_tally2_c_solver_par *solver_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_ccg_res(
    magma_tally2_c_matrix A, magma_tally2_c_matrix b, 
    magma_tally2_c_matrix *x, magma_tally2_c_solver_par *solver_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_ccg_merge(
    magma_tally2_c_matrix A, magma_tally2_c_matrix b, 
    magma_tally2_c_matrix *x, magma_tally2_c_solver_par *solver_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_cgmres(
    magma_tally2_c_matrix A, magma_tally2_c_matrix b, 
    magma_tally2_c_matrix *x, magma_tally2_c_solver_par *solver_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_cbicgstab(
    magma_tally2_c_matrix A, magma_tally2_c_matrix b, magma_tally2_c_matrix *x, 
    magma_tally2_c_solver_par *solver_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_cbicgstab_merge(
    magma_tally2_c_matrix A, magma_tally2_c_matrix b, 
    magma_tally2_c_matrix *x, magma_tally2_c_solver_par *solver_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_cbicgstab_merge2(
    magma_tally2_c_matrix A, magma_tally2_c_matrix b, 
    magma_tally2_c_matrix *x, magma_tally2_c_solver_par *solver_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_cpcg(
    magma_tally2_c_matrix A, magma_tally2_c_matrix b, 
    magma_tally2_c_matrix *x, magma_tally2_c_solver_par *solver_par, 
    magma_tally2_c_preconditioner *precond_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_cbpcg(
    magma_tally2_c_matrix A, magma_tally2_c_matrix b, 
    magma_tally2_c_matrix *x, magma_tally2_c_solver_par *solver_par, 
    magma_tally2_c_preconditioner *precond_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_cpbicgstab(
    magma_tally2_c_matrix A, magma_tally2_c_matrix b, 
    magma_tally2_c_matrix *x, magma_tally2_c_solver_par *solver_par, 
    magma_tally2_c_preconditioner *precond_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_cpgmres(
    magma_tally2_c_matrix A, magma_tally2_c_matrix b, 
    magma_tally2_c_matrix *x, magma_tally2_c_solver_par *solver_par, 
    magma_tally2_c_preconditioner *precond_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_cfgmres(
    magma_tally2_c_matrix A, magma_tally2_c_matrix b, 
    magma_tally2_c_matrix *x, magma_tally2_c_solver_par *solver_par, 
    magma_tally2_c_preconditioner *precond_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_cjacobi(
    magma_tally2_c_matrix A, magma_tally2_c_matrix b, 
    magma_tally2_c_matrix *x, magma_tally2_c_solver_par *solver_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_cjacobidomainoverlap(
    magma_tally2_c_matrix A, 
    magma_tally2_c_matrix b, 
    magma_tally2_c_matrix *x,  
    magma_tally2_c_solver_par *solver_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_cbaiter(
    magma_tally2_c_matrix A, magma_tally2_c_matrix b, 
    magma_tally2_c_matrix *x, magma_tally2_c_solver_par *solver_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_citerref(
    magma_tally2_c_matrix A, magma_tally2_c_matrix b, 
    magma_tally2_c_matrix *x, magma_tally2_c_solver_par *solver_par, 
    magma_tally2_c_preconditioner *precond_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_cilu(
    magma_tally2_c_matrix A, magma_tally2_c_matrix b, 
    magma_tally2_c_matrix *x, magma_tally2_c_solver_par *solver_par, 
    magma_tally2_int_t *ipiv,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_cbcsrlu(
    magma_tally2_c_matrix A, magma_tally2_c_matrix b, 
    magma_tally2_c_matrix *x, magma_tally2_c_solver_par *solver_par,
    magma_tally2_queue_t queue );


magma_tally2_int_t
magma_tally2_cbcsrlutrf(
    magma_tally2_c_matrix A, 
    magma_tally2_c_matrix *M,
    magma_tally2_int_t *ipiv, 
    magma_tally2_int_t version,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_cbcsrlusv(
    magma_tally2_c_matrix A, magma_tally2_c_matrix b, 
    magma_tally2_c_matrix *x, magma_tally2_c_solver_par *solver_par, 
    magma_tally2_int_t *ipiv,
    magma_tally2_queue_t queue );



magma_tally2_int_t
magma_tally2_cilucg(
    magma_tally2_c_matrix A, magma_tally2_c_matrix b, 
    magma_tally2_c_matrix *x, magma_tally2_c_solver_par *solver_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_cilugmres(
    magma_tally2_c_matrix A, magma_tally2_c_matrix b, 
    magma_tally2_c_matrix *x, magma_tally2_c_solver_par *solver_par,
    magma_tally2_queue_t queue ); 


magma_tally2_int_t
magma_tally2_clobpcg_shift(
    magma_tally2_int_t num_rows,
    magma_tally2_int_t num_vecs, 
    magma_tally2_int_t shift,
    magma_tally2FloatComplex_ptr x,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_clobpcg_res(
    magma_tally2_int_t num_rows,
    magma_tally2_int_t num_vecs, 
    float *evalues, 
    magma_tally2FloatComplex_ptr X,
    magma_tally2FloatComplex_ptr R, 
    float *res,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_clobpcg_maxpy(
    magma_tally2_int_t num_rows,
    magma_tally2_int_t num_vecs, 
    magma_tally2FloatComplex_ptr X,
    magma_tally2FloatComplex_ptr Y,
    magma_tally2_queue_t queue );


/*/////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally2_SPARSE eigensolvers (Data on GPU)
*/
magma_tally2_int_t
magma_tally2_clobpcg(
    magma_tally2_c_matrix A, 
    magma_tally2_c_solver_par *solver_par,
    magma_tally2_c_preconditioner *precond_par, 
    magma_tally2_queue_t queue );




/*/////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally2_SPARSE preconditioners (Data on GPU)
*/

magma_tally2_int_t
magma_tally2_cjacobisetup(
    magma_tally2_c_matrix A, 
    magma_tally2_c_matrix b, 
    magma_tally2_c_matrix *M, 
    magma_tally2_c_matrix *c,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_cjacobisetup_matrix(
    magma_tally2_c_matrix A, 
    magma_tally2_c_matrix *M, 
    magma_tally2_c_matrix *d,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_cjacobisetup_vector(
    magma_tally2_c_matrix b, 
    magma_tally2_c_matrix d, 
    magma_tally2_c_matrix *c,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_cjacobiiter(
    magma_tally2_c_matrix M, 
    magma_tally2_c_matrix c, 
    magma_tally2_c_matrix *x, 
    magma_tally2_c_solver_par *solver_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_cjacobiiter_precond( 
    magma_tally2_c_matrix M, 
    magma_tally2_c_matrix *x, 
    magma_tally2_c_solver_par *solver_par, 
    magma_tally2_c_preconditioner *precond,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_cjacobiiter_sys(
    magma_tally2_c_matrix A, 
    magma_tally2_c_matrix b, 
    magma_tally2_c_matrix d, 
    magma_tally2_c_matrix t, 
    magma_tally2_c_matrix *x,  
    magma_tally2_c_solver_par *solver_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_cpastixsetup(
    magma_tally2_c_matrix A, magma_tally2_c_matrix b,
    magma_tally2_c_preconditioner *precond,
    magma_tally2_queue_t queue );


magma_tally2_int_t
magma_tally2_capplypastix(
    magma_tally2_c_matrix b, magma_tally2_c_matrix *x, 
    magma_tally2_c_preconditioner *precond,
    magma_tally2_queue_t queue );


// custom preconditioner
magma_tally2_int_t
magma_tally2_capplycustomprecond_l(
    magma_tally2_c_matrix b, 
    magma_tally2_c_matrix *x, 
    magma_tally2_c_preconditioner *precond,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_capplycustomprecond_r(
    magma_tally2_c_matrix b, 
    magma_tally2_c_matrix *x, 
    magma_tally2_c_preconditioner *precond,
    magma_tally2_queue_t queue );


// CUSPARSE preconditioner

magma_tally2_int_t
magma_tally2_ccuilusetup(
    magma_tally2_c_matrix A, magma_tally2_c_preconditioner *precond,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_capplycuilu_l(
    magma_tally2_c_matrix b, magma_tally2_c_matrix *x, 
    magma_tally2_c_preconditioner *precond,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_capplycuilu_r(
    magma_tally2_c_matrix b, magma_tally2_c_matrix *x, 
    magma_tally2_c_preconditioner *precond,
    magma_tally2_queue_t queue );


magma_tally2_int_t
magma_tally2_ccuiccsetup(
    magma_tally2_c_matrix A, magma_tally2_c_preconditioner *precond,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_capplycuicc_l(
    magma_tally2_c_matrix b, magma_tally2_c_matrix *x, 
    magma_tally2_c_preconditioner *precond,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_capplycuicc_r(
    magma_tally2_c_matrix b, magma_tally2_c_matrix *x, 
    magma_tally2_c_preconditioner *precond,
    magma_tally2_queue_t queue );


magma_tally2_int_t
magma_tally2_ccumilusetup(
    magma_tally2_c_matrix A, 
    magma_tally2_c_preconditioner *precond,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_ccumilugeneratesolverinfo(
    magma_tally2_c_preconditioner *precond,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_capplycumilu_l(
    magma_tally2_c_matrix b, 
    magma_tally2_c_matrix *x, 
    magma_tally2_c_preconditioner *precond,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_capplycumilu_r(
    magma_tally2_c_matrix b, 
    magma_tally2_c_matrix *x, 
    magma_tally2_c_preconditioner *precond,
    magma_tally2_queue_t queue );


magma_tally2_int_t
magma_tally2_ccumiccsetup(
    magma_tally2_c_matrix A, 
    magma_tally2_c_preconditioner *precond,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_ccumicgeneratesolverinfo(
    magma_tally2_c_preconditioner *precond,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_capplycumicc_l(
    magma_tally2_c_matrix b, 
    magma_tally2_c_matrix *x, 
    magma_tally2_c_preconditioner *precond,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_capplycumicc_r(
    magma_tally2_c_matrix b, 
    magma_tally2_c_matrix *x, 
    magma_tally2_c_preconditioner *precond,
    magma_tally2_queue_t queue );


// block-asynchronous iteration

magma_tally2_int_t
magma_tally2_cbajac_csr(
    magma_tally2_int_t localiters,
    magma_tally2_c_matrix D,
    magma_tally2_c_matrix R,
    magma_tally2_c_matrix b,
    magma_tally2_c_matrix *x,
    magma_tally2_queue_t queue );

/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally2_SPARSE utility function definitions
*/

magma_tally2_int_t
magma_tally2_c_spmv(
    magma_tally2FloatComplex alpha, 
    magma_tally2_c_matrix A, 
    magma_tally2_c_matrix x, 
    magma_tally2FloatComplex beta, 
    magma_tally2_c_matrix y,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_ccustomspmv(
    magma_tally2FloatComplex alpha, 
    magma_tally2_c_matrix x, 
    magma_tally2FloatComplex beta, 
    magma_tally2_c_matrix y,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_c_spmv_shift(
    magma_tally2FloatComplex alpha, 
    magma_tally2_c_matrix A, 
    magma_tally2FloatComplex lambda,
    magma_tally2_c_matrix x, 
    magma_tally2FloatComplex beta, 
    magma_tally2_int_t offset, 
    magma_tally2_int_t blocksize,
    magma_tally2Index_ptr dadd_vecs, 
    magma_tally2_c_matrix y,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_ccuspmm(
    magma_tally2_c_matrix A, 
    magma_tally2_c_matrix B, 
    magma_tally2_c_matrix *AB,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_c_spmm(
    magma_tally2FloatComplex alpha, 
    magma_tally2_c_matrix A,
    magma_tally2_c_matrix B,
    magma_tally2_c_matrix *C,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_csymbilu( 
    magma_tally2_c_matrix *A, 
    magma_tally2_int_t levels, 
    magma_tally2_c_matrix *L, 
    magma_tally2_c_matrix *U,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_ccuspaxpy(
    magma_tally2FloatComplex_ptr alpha, magma_tally2_c_matrix A, 
    magma_tally2FloatComplex_ptr beta, magma_tally2_c_matrix B, 
    magma_tally2_c_matrix *AB,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_c_precond(
    magma_tally2_c_matrix A, 
    magma_tally2_c_matrix b, magma_tally2_c_matrix *x,
    magma_tally2_c_preconditioner *precond,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_c_solver(
    magma_tally2_c_matrix A, magma_tally2_c_matrix b, 
    magma_tally2_c_matrix *x, magma_tally2_copts *zopts,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_c_precondsetup(
    magma_tally2_c_matrix A, magma_tally2_c_matrix b, 
    magma_tally2_c_preconditioner *precond,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_c_applyprecond(
    magma_tally2_c_matrix A, magma_tally2_c_matrix b, 
    magma_tally2_c_matrix *x, magma_tally2_c_preconditioner *precond,
    magma_tally2_queue_t queue );


magma_tally2_int_t
magma_tally2_c_applyprecond_left(
    magma_tally2_c_matrix A, magma_tally2_c_matrix b, 
    magma_tally2_c_matrix *x, magma_tally2_c_preconditioner *precond,
    magma_tally2_queue_t queue );


magma_tally2_int_t
magma_tally2_c_applyprecond_right(
    magma_tally2_c_matrix A, magma_tally2_c_matrix b, 
    magma_tally2_c_matrix *x, magma_tally2_c_preconditioner *precond,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_c_initP2P(
    magma_tally2_int_t *bandwidth_benchmark,
    magma_tally2_int_t *num_gpus,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_ccompact(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2FloatComplex_ptr dA, magma_tally2_int_t ldda,
    float *dnorms, float tol, 
    magma_tally2_int_t *activeMask, magma_tally2_int_t *cBlockSize,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_ccompactActive(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2FloatComplex_ptr dA, magma_tally2_int_t ldda, 
    magma_tally2_int_t *active,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_cmlumerge(    
    magma_tally2_c_matrix L, 
    magma_tally2_c_matrix U,
    magma_tally2_c_matrix *A, 
    magma_tally2_queue_t queue );


/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally2_SPARSE BLAS function definitions
*/
magma_tally2_int_t 
magma_tally2_cgecsrmv(
    magma_tally2_trans_t transA,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2FloatComplex alpha,
    magma_tally2FloatComplex_ptr dval,
    magma_tally2Index_ptr drowptr,
    magma_tally2Index_ptr dcolind,
    magma_tally2FloatComplex_ptr dx,
    magma_tally2FloatComplex beta,
    magma_tally2FloatComplex_ptr dy,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_cgecsrmv_shift(
    magma_tally2_trans_t transA,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2FloatComplex alpha,
    magma_tally2FloatComplex lambda,
    magma_tally2FloatComplex_ptr dval,
    magma_tally2Index_ptr drowptr,
    magma_tally2Index_ptr dcolind,
    magma_tally2FloatComplex_ptr dx,
    magma_tally2FloatComplex beta,
    int offset,
    int blocksize,
    magma_tally2Index_ptr dadd_rows,
    magma_tally2FloatComplex_ptr dy,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_cmgecsrmv(
    magma_tally2_trans_t transA,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2_int_t num_vecs,
    magma_tally2FloatComplex alpha,
    magma_tally2FloatComplex_ptr dval,
    magma_tally2Index_ptr drowptr,
    magma_tally2Index_ptr dcolind,
    magma_tally2FloatComplex_ptr dx,
    magma_tally2FloatComplex beta,
    magma_tally2FloatComplex_ptr dy,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_cgeellmv(
    magma_tally2_trans_t transA,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2_int_t nnz_per_row,
    magma_tally2FloatComplex alpha,
    magma_tally2FloatComplex_ptr dval,
    magma_tally2Index_ptr dcolind,
    magma_tally2FloatComplex_ptr dx,
    magma_tally2FloatComplex beta,
    magma_tally2FloatComplex_ptr dy,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_cgeellmv_shift(
    magma_tally2_trans_t transA,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2_int_t nnz_per_row,
    magma_tally2FloatComplex alpha,
    magma_tally2FloatComplex lambda,
    magma_tally2FloatComplex_ptr dval,
    magma_tally2Index_ptr dcolind,
    magma_tally2FloatComplex_ptr dx,
    magma_tally2FloatComplex beta,
    int offset,
    int blocksize,
    magma_tally2Index_ptr dadd_rows,
    magma_tally2FloatComplex_ptr dy,
    magma_tally2_queue_t queue );


magma_tally2_int_t 
magma_tally2_cmgeellmv(
    magma_tally2_trans_t transA,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2_int_t num_vecs,
    magma_tally2_int_t nnz_per_row,
    magma_tally2FloatComplex alpha,
    magma_tally2FloatComplex_ptr dval,
    magma_tally2Index_ptr dcolind,
    magma_tally2FloatComplex_ptr dx,
    magma_tally2FloatComplex beta,
    magma_tally2FloatComplex_ptr dy,
    magma_tally2_queue_t queue );


magma_tally2_int_t 
magma_tally2_cgeelltmv(
    magma_tally2_trans_t transA,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2_int_t nnz_per_row,
    magma_tally2FloatComplex alpha,
    magma_tally2FloatComplex_ptr dval,
    magma_tally2Index_ptr dcolind,
    magma_tally2FloatComplex_ptr dx,
    magma_tally2FloatComplex beta,
    magma_tally2FloatComplex_ptr dy,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_cgeelltmv_shift(
    magma_tally2_trans_t transA,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2_int_t nnz_per_row,
    magma_tally2FloatComplex alpha,
    magma_tally2FloatComplex lambda,
    magma_tally2FloatComplex_ptr dval,
    magma_tally2Index_ptr dcolind,
    magma_tally2FloatComplex_ptr dx,
    magma_tally2FloatComplex beta,
    int offset,
    int blocksize,
    magma_tally2Index_ptr dadd_rows,
    magma_tally2FloatComplex_ptr dy,
    magma_tally2_queue_t queue );


magma_tally2_int_t 
magma_tally2_cmgeelltmv(
    magma_tally2_trans_t transA,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2_int_t num_vecs,
    magma_tally2_int_t nnz_per_row,
    magma_tally2FloatComplex alpha,
    magma_tally2FloatComplex_ptr dval,
    magma_tally2Index_ptr dcolind,
    magma_tally2FloatComplex_ptr dx,
    magma_tally2FloatComplex beta,
    magma_tally2FloatComplex_ptr dy,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_cgeellrtmv(
    magma_tally2_trans_t transA,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2_int_t nnz_per_row,
    magma_tally2FloatComplex alpha,
    magma_tally2FloatComplex_ptr dval,
    magma_tally2Index_ptr dcolind,
    magma_tally2Index_ptr drowlength,
    magma_tally2FloatComplex_ptr dx,
    magma_tally2FloatComplex beta,
    magma_tally2FloatComplex_ptr dy,
    magma_tally2_int_t num_threads,
    magma_tally2_int_t threads_per_row,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_cgesellcmv(
    magma_tally2_trans_t transA,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2_int_t blocksize,
    magma_tally2_int_t slices,
    magma_tally2_int_t alignment,
    magma_tally2FloatComplex alpha,
    magma_tally2FloatComplex_ptr dval,
    magma_tally2Index_ptr dcolind,
    magma_tally2Index_ptr drowptr,
    magma_tally2FloatComplex_ptr dx,
    magma_tally2FloatComplex beta,
    magma_tally2FloatComplex_ptr dy,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_cgesellpmv(
    magma_tally2_trans_t transA,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2_int_t blocksize,
    magma_tally2_int_t slices,
    magma_tally2_int_t alignment,
    magma_tally2FloatComplex alpha,
    magma_tally2FloatComplex_ptr dval,
    magma_tally2Index_ptr dcolind,
    magma_tally2Index_ptr drowptr,
    magma_tally2FloatComplex_ptr dx,
    magma_tally2FloatComplex beta,
    magma_tally2FloatComplex_ptr dy,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_cmgesellpmv(
    magma_tally2_trans_t transA,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2_int_t num_vecs,
    magma_tally2_int_t blocksize,
    magma_tally2_int_t slices,
    magma_tally2_int_t alignment,
    magma_tally2FloatComplex alpha,
    magma_tally2FloatComplex_ptr dval,
    magma_tally2Index_ptr dcolind,
    magma_tally2Index_ptr drowptr,
    magma_tally2FloatComplex_ptr dx,
    magma_tally2FloatComplex beta,
    magma_tally2FloatComplex_ptr dy,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_cmgesellpmv_blocked(
    magma_tally2_trans_t transA,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2_int_t num_vecs,
    magma_tally2_int_t blocksize,
    magma_tally2_int_t slices,
    magma_tally2_int_t alignment,
    magma_tally2FloatComplex alpha,
    magma_tally2FloatComplex_ptr dval,
    magma_tally2Index_ptr dcolind,
    magma_tally2Index_ptr drowptr,
    magma_tally2FloatComplex_ptr dx,
    magma_tally2FloatComplex beta,
    magma_tally2FloatComplex_ptr dy,
    magma_tally2_queue_t queue );


magma_tally2_int_t
magma_tally2_cmergedgs(
    magma_tally2_int_t n, 
    magma_tally2_int_t ldh,
    magma_tally2_int_t k, 
    magma_tally2FloatComplex_ptr dv, 
    magma_tally2FloatComplex_ptr dr,
    magma_tally2FloatComplex_ptr dskp,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_ccopyscale(    
    int n, 
    int k,
    magma_tally2FloatComplex_ptr dr, 
    magma_tally2FloatComplex_ptr dv,
    magma_tally2FloatComplex_ptr dskp,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_scnrm2scale(    
    int m, 
    magma_tally2FloatComplex_ptr dr,    
    int lddr, 
    magma_tally2FloatComplex *drnorm,
    magma_tally2_queue_t queue );


magma_tally2_int_t
magma_tally2_cjacobisetup_vector_gpu(
    int num_rows, 
    magma_tally2_c_matrix b, 
    magma_tally2_c_matrix d, 
    magma_tally2_c_matrix c,
    magma_tally2_c_matrix *x,
    magma_tally2_queue_t queue );


magma_tally2_int_t
magma_tally2_cjacobi_diagscal(    
    int num_rows, 
    magma_tally2_c_matrix d, 
    magma_tally2_c_matrix b, 
    magma_tally2_c_matrix *c,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_cjacobiupdate(
    magma_tally2_c_matrix t, 
    magma_tally2_c_matrix b, 
    magma_tally2_c_matrix d, 
    magma_tally2_c_matrix *x,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_cjacobispmvupdate(
    magma_tally2_int_t maxiter,
    magma_tally2_c_matrix A, 
    magma_tally2_c_matrix t, 
    magma_tally2_c_matrix b, 
    magma_tally2_c_matrix d, 
    magma_tally2_c_matrix *x,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_cjacobispmvupdate_bw(
    magma_tally2_int_t maxiter,
    magma_tally2_c_matrix A, 
    magma_tally2_c_matrix t, 
    magma_tally2_c_matrix b, 
    magma_tally2_c_matrix d, 
    magma_tally2_c_matrix *x,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_cjacobispmvupdateselect(
    magma_tally2_int_t maxiter,
    magma_tally2_int_t num_updates,
    magma_tally2_index_t *indices,
    magma_tally2_c_matrix A,
    magma_tally2_c_matrix t, 
    magma_tally2_c_matrix b, 
    magma_tally2_c_matrix d, 
    magma_tally2_c_matrix tmp, 
    magma_tally2_c_matrix *x,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_cjacobisetup_diagscal(
    magma_tally2_c_matrix A, magma_tally2_c_matrix *d,
    magma_tally2_queue_t queue );


magma_tally2_int_t
magma_tally2_cbicgmerge1(    
    int n, 
    magma_tally2FloatComplex_ptr dskp,
    magma_tally2FloatComplex_ptr dv, 
    magma_tally2FloatComplex_ptr dr, 
    magma_tally2FloatComplex_ptr dp,
    magma_tally2_queue_t queue );


magma_tally2_int_t
magma_tally2_cbicgmerge2(
    int n, 
    magma_tally2FloatComplex_ptr dskp, 
    magma_tally2FloatComplex_ptr dr,
    magma_tally2FloatComplex_ptr dv, 
    magma_tally2FloatComplex_ptr ds,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_cbicgmerge3(
    int n, 
    magma_tally2FloatComplex_ptr dskp, 
    magma_tally2FloatComplex_ptr dp,
    magma_tally2FloatComplex_ptr ds,
    magma_tally2FloatComplex_ptr dt,
    magma_tally2FloatComplex_ptr dx, 
    magma_tally2FloatComplex_ptr dr,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_cbicgmerge4(
    int type, 
    magma_tally2FloatComplex_ptr dskp,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_ccgmerge_spmv1( 
    magma_tally2_c_matrix A,
    magma_tally2FloatComplex_ptr d1,
    magma_tally2FloatComplex_ptr d2,
    magma_tally2FloatComplex_ptr dd,
    magma_tally2FloatComplex_ptr dz,
    magma_tally2FloatComplex_ptr dskp,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_ccgmerge_xrbeta( 
    int n,
    magma_tally2FloatComplex_ptr d1,
    magma_tally2FloatComplex_ptr d2,
    magma_tally2FloatComplex_ptr dx,
    magma_tally2FloatComplex_ptr dr,
    magma_tally2FloatComplex_ptr dd,
    magma_tally2FloatComplex_ptr dz, 
    magma_tally2FloatComplex_ptr dskp,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_cmdotc(
    magma_tally2_int_t n, 
    magma_tally2_int_t k, 
    magma_tally2FloatComplex_ptr dv, 
    magma_tally2FloatComplex_ptr dr,
    magma_tally2FloatComplex_ptr dd1,
    magma_tally2FloatComplex_ptr dd2,
    magma_tally2FloatComplex_ptr dskp,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_cgemvmdot(
    int n, 
    int k, 
    magma_tally2FloatComplex_ptr dv, 
    magma_tally2FloatComplex_ptr dr,
    magma_tally2FloatComplex_ptr dd1,
    magma_tally2FloatComplex_ptr dd2,
    magma_tally2FloatComplex_ptr dskp,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_cbicgmerge_spmv1( 
    magma_tally2_c_matrix A,
    magma_tally2FloatComplex_ptr dd1,
    magma_tally2FloatComplex_ptr dd2,
    magma_tally2FloatComplex_ptr dp,
    magma_tally2FloatComplex_ptr dr,
    magma_tally2FloatComplex_ptr dv,
    magma_tally2FloatComplex_ptr dskp,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_cbicgmerge_spmv2( 
    magma_tally2_c_matrix A,
    magma_tally2FloatComplex_ptr dd1,
    magma_tally2FloatComplex_ptr dd2,
    magma_tally2FloatComplex_ptr ds,
    magma_tally2FloatComplex_ptr dt,
    magma_tally2FloatComplex_ptr dskp,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_cbicgmerge_xrbeta( 
    int n,
    magma_tally2FloatComplex_ptr dd1,
    magma_tally2FloatComplex_ptr dd2,
    magma_tally2FloatComplex_ptr drr,
    magma_tally2FloatComplex_ptr dr,
    magma_tally2FloatComplex_ptr dp,
    magma_tally2FloatComplex_ptr ds,
    magma_tally2FloatComplex_ptr dt,
    magma_tally2FloatComplex_ptr dx, 
    magma_tally2FloatComplex_ptr dskp,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_cbcsrswp(
    magma_tally2_int_t n,
    magma_tally2_int_t size_b, 
    magma_tally2_int_t *ipiv,
    magma_tally2FloatComplex_ptr dx,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_cbcsrtrsv(
    magma_tally2_uplo_t uplo,
    magma_tally2_int_t r_blocks,
    magma_tally2_int_t c_blocks,
    magma_tally2_int_t size_b, 
    magma_tally2FloatComplex_ptr dA,
    magma_tally2_index_t *blockinfo, 
    magma_tally2FloatComplex_ptr dx,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_cbcsrvalcpy(
    magma_tally2_int_t size_b, 
    magma_tally2_int_t num_blocks, 
    magma_tally2_int_t num_zero_blocks, 
    magma_tally2FloatComplex_ptr *dAval, 
    magma_tally2FloatComplex_ptr *dBval,
    magma_tally2FloatComplex_ptr *dBval2,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_cbcsrluegemm(
    magma_tally2_int_t size_b, 
    magma_tally2_int_t num_block_rows,
    magma_tally2_int_t kblocks,
    magma_tally2FloatComplex_ptr *dA, 
    magma_tally2FloatComplex_ptr *dB, 
    magma_tally2FloatComplex_ptr *dC,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_cbcsrlupivloc(
    magma_tally2_int_t size_b, 
    magma_tally2_int_t kblocks,
    magma_tally2FloatComplex_ptr *dA, 
    magma_tally2_int_t *ipiv,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_cbcsrblockinfo5(
    magma_tally2_int_t lustep,
    magma_tally2_int_t num_blocks, 
    magma_tally2_int_t c_blocks, 
    magma_tally2_int_t size_b,
    magma_tally2_index_t *blockinfo,
    magma_tally2FloatComplex_ptr dval,
    magma_tally2FloatComplex_ptr *AII,
    magma_tally2_queue_t queue );

 
#ifdef __cplusplus
}
#endif

#undef PRECISION_c
#endif /* MAGMA_tally2SPARSE_C_H */
