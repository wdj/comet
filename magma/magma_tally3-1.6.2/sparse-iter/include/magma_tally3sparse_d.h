/*
 -- MAGMA_tally3 (version 1.1) --
 Univ. of Tennessee, Knoxville
 Univ. of California, Berkeley
 Univ. of Colorado, Denver
 @date

 @generated from magma_tally3sparse_z.h normal z -> d, Mon May  4 17:11:25 2015
 @author Hartwig Anzt
*/

#ifndef MAGMA_tally3SPARSE_D_H
#define MAGMA_tally3SPARSE_D_H

#include "magma_tally3_types.h"
#include "magma_tally3sparse_types.h"

#define PRECISION_d


#ifdef __cplusplus
extern "C" {
#endif


/* ////////////////////////////////////////////////////////////////////////////
 -- For backwards compatability, map old (1.6.1) to new (1.6.2) function names
*/

#define magma_tally3_d_mtranspose  magma_tally3_dmtranspose
#define magma_tally3_d_mtransfer   magma_tally3_dmtransfer
#define magma_tally3_d_vtransfer   magma_tally3_dmtransfer
#define magma_tally3_d_mconvert    magma_tally3_dmconvert
#define magma_tally3_d_vinit       magma_tally3_dvinit
#define magma_tally3_d_vvisu       magma_tally3_dprint_vector
#define magma_tally3_d_vread       magma_tally3_dvread
#define magma_tally3_d_vspread     magma_tally3_dvspread
#define magma_tally3_d_mvisu       magma_tally3_dprint_matrix
#define magma_tally3_d_mfree       magma_tally3_dmfree
#define magma_tally3_d_vfree       magma_tally3_dmfree
#define write_d_csr_mtx     magma_tally3_dwrite_csr_mtx
#define write_d_csrtomtx    magma_tally3_dwrite_csrtomtx
#define print_d_csr         magma_tally3_dprint_csr_mtx


/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally3_SPARSE Auxiliary functions
*/


magma_tally3_int_t
magma_tally3_dparse_opts( 
    int argc, 
    char** argv, 
    magma_tally3_dopts *opts, 
    int *matrices, 
    magma_tally3_queue_t queue );

magma_tally3_int_t 
read_d_csr_from_binary( 
    magma_tally3_int_t* n_row, 
    magma_tally3_int_t* n_col, 
    magma_tally3_int_t* nnz, 
    double **val, 
    magma_tally3_index_t **row, 
    magma_tally3_index_t **col,
    const char * filename,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
read_d_csr_from_mtx( 
    magma_tally3_storage_t *type, 
    magma_tally3_location_t *location,
    magma_tally3_int_t* n_row, 
    magma_tally3_int_t* n_col, 
    magma_tally3_int_t* nnz, 
    double **val, 
    magma_tally3_index_t **row, 
    magma_tally3_index_t **col, 
    const char *filename,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_d_csr_mtx( 
    magma_tally3_d_matrix *A, 
    const char *filename,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_dcsrset( 
    magma_tally3_int_t m, 
    magma_tally3_int_t n, 
    magma_tally3_index_t *row, 
    magma_tally3_index_t *col, 
    double *val,
    magma_tally3_d_matrix *A,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_dcsrget( 
    magma_tally3_d_matrix A,
    magma_tally3_int_t *m, 
    magma_tally3_int_t *n, 
    magma_tally3_index_t **row, 
    magma_tally3_index_t **col, 
    double **val,
    magma_tally3_queue_t queue );


magma_tally3_int_t 
magma_tally3_dvset( 
    magma_tally3_int_t m, 
    magma_tally3_int_t n, 
    double *val,
    magma_tally3_d_matrix *v,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_dvget( 
    magma_tally3_d_matrix v,
    magma_tally3_int_t *m, 
    magma_tally3_int_t *n, 
    double **val,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_dvset_dev( 
    magma_tally3_int_t m, 
    magma_tally3_int_t n, 
    magma_tally3Double_ptr val,
    magma_tally3_d_matrix *v,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_dvget_dev( 
    magma_tally3_d_matrix v,
    magma_tally3_int_t *m, 
    magma_tally3_int_t *n, 
    magma_tally3Double_ptr *val,
    magma_tally3_queue_t queue );


magma_tally3_int_t 
magma_tally3_d_csr_mtxsymm( 
    magma_tally3_d_matrix *A, 
    const char *filename,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_d_csr_compressor( 
    double ** val, 
    magma_tally3_index_t ** row, 
    magma_tally3_index_t ** col, 
    double ** valn, 
    magma_tally3_index_t ** rown, 
    magma_tally3_index_t ** coln, 
    magma_tally3_int_t *n,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dmcsrcompressor( 
    magma_tally3_d_matrix *A,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dmcsrcompressor_gpu( 
    magma_tally3_d_matrix *A,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dvtranspose( 
    magma_tally3_d_matrix x,
    magma_tally3_d_matrix *y,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_d_cucsrtranspose( 
    magma_tally3_d_matrix A, 
    magma_tally3_d_matrix *B,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
d_transpose_csr( 
    magma_tally3_int_t n_rows, 
    magma_tally3_int_t n_cols, 
    magma_tally3_int_t nnz,
    double *val, 
    magma_tally3_index_t *row, 
    magma_tally3_index_t *col, 
    magma_tally3_int_t *new_n_rows, 
    magma_tally3_int_t *new_n_cols, 
    magma_tally3_int_t *new_nnz, 
    double **new_val, 
    magma_tally3_index_t **new_row, 
    magma_tally3_index_t **new_col,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dcsrsplit( 
    magma_tally3_int_t bsize,
    magma_tally3_d_matrix A,
    magma_tally3_d_matrix *D,
    magma_tally3_d_matrix *R,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dmscale( 
    magma_tally3_d_matrix *A, 
    magma_tally3_scale_t scaling,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_dmdiff( 
    magma_tally3_d_matrix A, 
    magma_tally3_d_matrix B, 
 real_Double_t *res,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dmdiagadd( 
    magma_tally3_d_matrix *A, 
    double add,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_dmsort( 
    magma_tally3_d_matrix *A,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dindexsort(
    magma_tally3_index_t *x, 
    magma_tally3_int_t first,
    magma_tally3_int_t last,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_ddomainoverlap(
    magma_tally3_index_t num_rows,
    magma_tally3_index_t *num_indices,
    magma_tally3_index_t *rowptr,
    magma_tally3_index_t *colidx,
    magma_tally3_index_t *x,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dsymbilu( 
    magma_tally3_d_matrix *A, 
    magma_tally3_int_t levels,
    magma_tally3_d_matrix *L,
    magma_tally3_d_matrix *U,
    magma_tally3_queue_t queue );


magma_tally3_int_t 
magma_tally3_dwrite_csr_mtx( 
    magma_tally3_d_matrix A,
    magma_tally3_order_t MajorType,
 const char *filename,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_dwrite_csrtomtx( 
    magma_tally3_d_matrix A,
    const char *filename,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_dprint_csr( 
    magma_tally3_int_t n_row, 
    magma_tally3_int_t n_col, 
    magma_tally3_int_t nnz, 
    double **val, 
    magma_tally3_index_t **row, 
    magma_tally3_index_t **col,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_dprint_csr_mtx( 
    magma_tally3_int_t n_row, 
    magma_tally3_int_t n_col, 
    magma_tally3_int_t nnz, 
    double **val, 
    magma_tally3_index_t **row, 
    magma_tally3_index_t **col, 
    magma_tally3_order_t MajorType,
    magma_tally3_queue_t queue );


magma_tally3_int_t 
magma_tally3_dmtranspose(
    magma_tally3_d_matrix A, 
    magma_tally3_d_matrix *B,
    magma_tally3_queue_t queue );


magma_tally3_int_t 
magma_tally3_dmtransfer(
    magma_tally3_d_matrix A, 
    magma_tally3_d_matrix *B, 
    magma_tally3_location_t src, 
    magma_tally3_location_t dst,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_dmconvert(
    magma_tally3_d_matrix A, 
    magma_tally3_d_matrix *B, 
    magma_tally3_storage_t old_format, 
    magma_tally3_storage_t new_format,
    magma_tally3_queue_t queue );


magma_tally3_int_t
magma_tally3_dvinit(
    magma_tally3_d_matrix *x, 
    magma_tally3_location_t memory_location,
    magma_tally3_int_t num_rows, 
    magma_tally3_int_t num_cols,
    double values,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dprint_vector(
    magma_tally3_d_matrix x, 
    magma_tally3_int_t offset, 
    magma_tally3_int_t displaylength,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dvread(
    magma_tally3_d_matrix *x, 
    magma_tally3_int_t length,
    char * filename,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dvspread(
    magma_tally3_d_matrix *x, 
    const char * filename,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dprint_matrix(
    magma_tally3_d_matrix A,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_ddiameter(
    magma_tally3_d_matrix *A,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_drowentries(
    magma_tally3_d_matrix *A,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dmfree(
    magma_tally3_d_matrix *A,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dresidual(
    magma_tally3_d_matrix A, 
    magma_tally3_d_matrix b, 
    magma_tally3_d_matrix x, 
    double *res,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dresidualvec(
    magma_tally3_d_matrix A,
    magma_tally3_d_matrix b,
    magma_tally3_d_matrix x,
    magma_tally3_d_matrix *r,
    double *res,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dmgenerator(
    magma_tally3_int_t n,
    magma_tally3_int_t offdiags,
    magma_tally3_index_t *diag_offset,
    double *diag_vals,
    magma_tally3_d_matrix *A,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dm_27stencil(
    magma_tally3_int_t n,
    magma_tally3_d_matrix *A,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dm_5stencil(
    magma_tally3_int_t n,
    magma_tally3_d_matrix *A,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dsolverinfo(
    magma_tally3_d_solver_par *solver_par, 
    magma_tally3_d_preconditioner *precond_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dsolverinfo_init(
    magma_tally3_d_solver_par *solver_par, 
    magma_tally3_d_preconditioner *precond,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_deigensolverinfo_init(
    magma_tally3_d_solver_par *solver_par,
    magma_tally3_queue_t queue );


magma_tally3_int_t
magma_tally3_dsolverinfo_free(
    magma_tally3_d_solver_par *solver_par, 
    magma_tally3_d_preconditioner *precond,
    magma_tally3_queue_t queue );



/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally3_SPARSE iterative incomplete factorizations
*/


magma_tally3_int_t
magma_tally3_diterilusetup( 
    magma_tally3_d_matrix A, 
    magma_tally3_d_matrix b,                                 
    magma_tally3_d_preconditioner *precond,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_ditericsetup( 
    magma_tally3_d_matrix A, 
    magma_tally3_d_matrix b, 
    magma_tally3_d_preconditioner *precond,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_ditericupdate( 
    magma_tally3_d_matrix A, 
    magma_tally3_d_preconditioner *precond, 
    magma_tally3_int_t updates,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dapplyiteric_l( 
    magma_tally3_d_matrix b, 
    magma_tally3_d_matrix *x, 
    magma_tally3_d_preconditioner *precond,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dapplyiteric_r( 
    magma_tally3_d_matrix b, 
    magma_tally3_d_matrix *x, 
    magma_tally3_d_preconditioner *precond,
    magma_tally3_queue_t queue );


magma_tally3_int_t
magma_tally3_diterilu_csr( 
    magma_tally3_d_matrix A,
    magma_tally3_d_matrix L,
    magma_tally3_d_matrix U,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_diteric_csr( 
    magma_tally3_d_matrix A,
    magma_tally3_d_matrix A_CSR,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_dfrobenius( 
    magma_tally3_d_matrix A, 
    magma_tally3_d_matrix B, 
    real_Double_t *res,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_dnonlinres(   
    magma_tally3_d_matrix A, 
    magma_tally3_d_matrix L,
    magma_tally3_d_matrix U, 
    magma_tally3_d_matrix *LU, 
    real_Double_t *res,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_dilures(   
    magma_tally3_d_matrix A, 
    magma_tally3_d_matrix L,
    magma_tally3_d_matrix U, 
    magma_tally3_d_matrix *LU, 
    real_Double_t *res,
    real_Double_t *nonlinres,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_dicres(       
    magma_tally3_d_matrix A, 
    magma_tally3_d_matrix C,
    magma_tally3_d_matrix CT, 
    magma_tally3_d_matrix *LU, 
    real_Double_t *res,
    real_Double_t *nonlinres,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_dinitguess( 
    magma_tally3_d_matrix A, 
    magma_tally3_d_matrix *L, 
    magma_tally3_d_matrix *U,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_dinitrecursiveLU( 
    magma_tally3_d_matrix A, 
    magma_tally3_d_matrix *B,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_dmLdiagadd( 
    magma_tally3_d_matrix *L,
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
magma_tally3_dcg(
    magma_tally3_d_matrix A, magma_tally3_d_matrix b, 
    magma_tally3_d_matrix *x, magma_tally3_d_solver_par *solver_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_dcg_res(
    magma_tally3_d_matrix A, magma_tally3_d_matrix b, 
    magma_tally3_d_matrix *x, magma_tally3_d_solver_par *solver_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_dcg_merge(
    magma_tally3_d_matrix A, magma_tally3_d_matrix b, 
    magma_tally3_d_matrix *x, magma_tally3_d_solver_par *solver_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_dgmres(
    magma_tally3_d_matrix A, magma_tally3_d_matrix b, 
    magma_tally3_d_matrix *x, magma_tally3_d_solver_par *solver_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dbicgstab(
    magma_tally3_d_matrix A, magma_tally3_d_matrix b, magma_tally3_d_matrix *x, 
    magma_tally3_d_solver_par *solver_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dbicgstab_merge(
    magma_tally3_d_matrix A, magma_tally3_d_matrix b, 
    magma_tally3_d_matrix *x, magma_tally3_d_solver_par *solver_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dbicgstab_merge2(
    magma_tally3_d_matrix A, magma_tally3_d_matrix b, 
    magma_tally3_d_matrix *x, magma_tally3_d_solver_par *solver_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dpcg(
    magma_tally3_d_matrix A, magma_tally3_d_matrix b, 
    magma_tally3_d_matrix *x, magma_tally3_d_solver_par *solver_par, 
    magma_tally3_d_preconditioner *precond_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dbpcg(
    magma_tally3_d_matrix A, magma_tally3_d_matrix b, 
    magma_tally3_d_matrix *x, magma_tally3_d_solver_par *solver_par, 
    magma_tally3_d_preconditioner *precond_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dpbicgstab(
    magma_tally3_d_matrix A, magma_tally3_d_matrix b, 
    magma_tally3_d_matrix *x, magma_tally3_d_solver_par *solver_par, 
    magma_tally3_d_preconditioner *precond_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dpgmres(
    magma_tally3_d_matrix A, magma_tally3_d_matrix b, 
    magma_tally3_d_matrix *x, magma_tally3_d_solver_par *solver_par, 
    magma_tally3_d_preconditioner *precond_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dfgmres(
    magma_tally3_d_matrix A, magma_tally3_d_matrix b, 
    magma_tally3_d_matrix *x, magma_tally3_d_solver_par *solver_par, 
    magma_tally3_d_preconditioner *precond_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_djacobi(
    magma_tally3_d_matrix A, magma_tally3_d_matrix b, 
    magma_tally3_d_matrix *x, magma_tally3_d_solver_par *solver_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_djacobidomainoverlap(
    magma_tally3_d_matrix A, 
    magma_tally3_d_matrix b, 
    magma_tally3_d_matrix *x,  
    magma_tally3_d_solver_par *solver_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dbaiter(
    magma_tally3_d_matrix A, magma_tally3_d_matrix b, 
    magma_tally3_d_matrix *x, magma_tally3_d_solver_par *solver_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_diterref(
    magma_tally3_d_matrix A, magma_tally3_d_matrix b, 
    magma_tally3_d_matrix *x, magma_tally3_d_solver_par *solver_par, 
    magma_tally3_d_preconditioner *precond_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dilu(
    magma_tally3_d_matrix A, magma_tally3_d_matrix b, 
    magma_tally3_d_matrix *x, magma_tally3_d_solver_par *solver_par, 
    magma_tally3_int_t *ipiv,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dbcsrlu(
    magma_tally3_d_matrix A, magma_tally3_d_matrix b, 
    magma_tally3_d_matrix *x, magma_tally3_d_solver_par *solver_par,
    magma_tally3_queue_t queue );


magma_tally3_int_t
magma_tally3_dbcsrlutrf(
    magma_tally3_d_matrix A, 
    magma_tally3_d_matrix *M,
    magma_tally3_int_t *ipiv, 
    magma_tally3_int_t version,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dbcsrlusv(
    magma_tally3_d_matrix A, magma_tally3_d_matrix b, 
    magma_tally3_d_matrix *x, magma_tally3_d_solver_par *solver_par, 
    magma_tally3_int_t *ipiv,
    magma_tally3_queue_t queue );



magma_tally3_int_t
magma_tally3_dilucg(
    magma_tally3_d_matrix A, magma_tally3_d_matrix b, 
    magma_tally3_d_matrix *x, magma_tally3_d_solver_par *solver_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dilugmres(
    magma_tally3_d_matrix A, magma_tally3_d_matrix b, 
    magma_tally3_d_matrix *x, magma_tally3_d_solver_par *solver_par,
    magma_tally3_queue_t queue ); 


magma_tally3_int_t
magma_tally3_dlobpcg_shift(
    magma_tally3_int_t num_rows,
    magma_tally3_int_t num_vecs, 
    magma_tally3_int_t shift,
    magma_tally3Double_ptr x,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dlobpcg_res(
    magma_tally3_int_t num_rows,
    magma_tally3_int_t num_vecs, 
    double *evalues, 
    magma_tally3Double_ptr X,
    magma_tally3Double_ptr R, 
    double *res,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dlobpcg_maxpy(
    magma_tally3_int_t num_rows,
    magma_tally3_int_t num_vecs, 
    magma_tally3Double_ptr X,
    magma_tally3Double_ptr Y,
    magma_tally3_queue_t queue );


/*/////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally3_SPARSE eigensolvers (Data on GPU)
*/
magma_tally3_int_t
magma_tally3_dlobpcg(
    magma_tally3_d_matrix A, 
    magma_tally3_d_solver_par *solver_par,
    magma_tally3_d_preconditioner *precond_par, 
    magma_tally3_queue_t queue );




/*/////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally3_SPARSE preconditioners (Data on GPU)
*/

magma_tally3_int_t
magma_tally3_djacobisetup(
    magma_tally3_d_matrix A, 
    magma_tally3_d_matrix b, 
    magma_tally3_d_matrix *M, 
    magma_tally3_d_matrix *c,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_djacobisetup_matrix(
    magma_tally3_d_matrix A, 
    magma_tally3_d_matrix *M, 
    magma_tally3_d_matrix *d,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_djacobisetup_vector(
    magma_tally3_d_matrix b, 
    magma_tally3_d_matrix d, 
    magma_tally3_d_matrix *c,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_djacobiiter(
    magma_tally3_d_matrix M, 
    magma_tally3_d_matrix c, 
    magma_tally3_d_matrix *x, 
    magma_tally3_d_solver_par *solver_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_djacobiiter_precond( 
    magma_tally3_d_matrix M, 
    magma_tally3_d_matrix *x, 
    magma_tally3_d_solver_par *solver_par, 
    magma_tally3_d_preconditioner *precond,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_djacobiiter_sys(
    magma_tally3_d_matrix A, 
    magma_tally3_d_matrix b, 
    magma_tally3_d_matrix d, 
    magma_tally3_d_matrix t, 
    magma_tally3_d_matrix *x,  
    magma_tally3_d_solver_par *solver_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dpastixsetup(
    magma_tally3_d_matrix A, magma_tally3_d_matrix b,
    magma_tally3_d_preconditioner *precond,
    magma_tally3_queue_t queue );


magma_tally3_int_t
magma_tally3_dapplypastix(
    magma_tally3_d_matrix b, magma_tally3_d_matrix *x, 
    magma_tally3_d_preconditioner *precond,
    magma_tally3_queue_t queue );


// custom preconditioner
magma_tally3_int_t
magma_tally3_dapplycustomprecond_l(
    magma_tally3_d_matrix b, 
    magma_tally3_d_matrix *x, 
    magma_tally3_d_preconditioner *precond,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dapplycustomprecond_r(
    magma_tally3_d_matrix b, 
    magma_tally3_d_matrix *x, 
    magma_tally3_d_preconditioner *precond,
    magma_tally3_queue_t queue );


// CUSPARSE preconditioner

magma_tally3_int_t
magma_tally3_dcuilusetup(
    magma_tally3_d_matrix A, magma_tally3_d_preconditioner *precond,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dapplycuilu_l(
    magma_tally3_d_matrix b, magma_tally3_d_matrix *x, 
    magma_tally3_d_preconditioner *precond,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dapplycuilu_r(
    magma_tally3_d_matrix b, magma_tally3_d_matrix *x, 
    magma_tally3_d_preconditioner *precond,
    magma_tally3_queue_t queue );


magma_tally3_int_t
magma_tally3_dcuiccsetup(
    magma_tally3_d_matrix A, magma_tally3_d_preconditioner *precond,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dapplycuicc_l(
    magma_tally3_d_matrix b, magma_tally3_d_matrix *x, 
    magma_tally3_d_preconditioner *precond,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dapplycuicc_r(
    magma_tally3_d_matrix b, magma_tally3_d_matrix *x, 
    magma_tally3_d_preconditioner *precond,
    magma_tally3_queue_t queue );


magma_tally3_int_t
magma_tally3_dcumilusetup(
    magma_tally3_d_matrix A, 
    magma_tally3_d_preconditioner *precond,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dcumilugeneratesolverinfo(
    magma_tally3_d_preconditioner *precond,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dapplycumilu_l(
    magma_tally3_d_matrix b, 
    magma_tally3_d_matrix *x, 
    magma_tally3_d_preconditioner *precond,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dapplycumilu_r(
    magma_tally3_d_matrix b, 
    magma_tally3_d_matrix *x, 
    magma_tally3_d_preconditioner *precond,
    magma_tally3_queue_t queue );


magma_tally3_int_t
magma_tally3_dcumiccsetup(
    magma_tally3_d_matrix A, 
    magma_tally3_d_preconditioner *precond,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dcumicgeneratesolverinfo(
    magma_tally3_d_preconditioner *precond,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dapplycumicc_l(
    magma_tally3_d_matrix b, 
    magma_tally3_d_matrix *x, 
    magma_tally3_d_preconditioner *precond,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dapplycumicc_r(
    magma_tally3_d_matrix b, 
    magma_tally3_d_matrix *x, 
    magma_tally3_d_preconditioner *precond,
    magma_tally3_queue_t queue );


// block-asynchronous iteration

magma_tally3_int_t
magma_tally3_dbajac_csr(
    magma_tally3_int_t localiters,
    magma_tally3_d_matrix D,
    magma_tally3_d_matrix R,
    magma_tally3_d_matrix b,
    magma_tally3_d_matrix *x,
    magma_tally3_queue_t queue );

/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally3_SPARSE utility function definitions
*/

magma_tally3_int_t
magma_tally3_d_spmv(
    double alpha, 
    magma_tally3_d_matrix A, 
    magma_tally3_d_matrix x, 
    double beta, 
    magma_tally3_d_matrix y,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dcustomspmv(
    double alpha, 
    magma_tally3_d_matrix x, 
    double beta, 
    magma_tally3_d_matrix y,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_d_spmv_shift(
    double alpha, 
    magma_tally3_d_matrix A, 
    double lambda,
    magma_tally3_d_matrix x, 
    double beta, 
    magma_tally3_int_t offset, 
    magma_tally3_int_t blocksize,
    magma_tally3Index_ptr dadd_vecs, 
    magma_tally3_d_matrix y,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dcuspmm(
    magma_tally3_d_matrix A, 
    magma_tally3_d_matrix B, 
    magma_tally3_d_matrix *AB,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_d_spmm(
    double alpha, 
    magma_tally3_d_matrix A,
    magma_tally3_d_matrix B,
    magma_tally3_d_matrix *C,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dsymbilu( 
    magma_tally3_d_matrix *A, 
    magma_tally3_int_t levels, 
    magma_tally3_d_matrix *L, 
    magma_tally3_d_matrix *U,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dcuspaxpy(
    magma_tally3Double_ptr alpha, magma_tally3_d_matrix A, 
    magma_tally3Double_ptr beta, magma_tally3_d_matrix B, 
    magma_tally3_d_matrix *AB,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_d_precond(
    magma_tally3_d_matrix A, 
    magma_tally3_d_matrix b, magma_tally3_d_matrix *x,
    magma_tally3_d_preconditioner *precond,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_d_solver(
    magma_tally3_d_matrix A, magma_tally3_d_matrix b, 
    magma_tally3_d_matrix *x, magma_tally3_dopts *zopts,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_d_precondsetup(
    magma_tally3_d_matrix A, magma_tally3_d_matrix b, 
    magma_tally3_d_preconditioner *precond,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_d_applyprecond(
    magma_tally3_d_matrix A, magma_tally3_d_matrix b, 
    magma_tally3_d_matrix *x, magma_tally3_d_preconditioner *precond,
    magma_tally3_queue_t queue );


magma_tally3_int_t
magma_tally3_d_applyprecond_left(
    magma_tally3_d_matrix A, magma_tally3_d_matrix b, 
    magma_tally3_d_matrix *x, magma_tally3_d_preconditioner *precond,
    magma_tally3_queue_t queue );


magma_tally3_int_t
magma_tally3_d_applyprecond_right(
    magma_tally3_d_matrix A, magma_tally3_d_matrix b, 
    magma_tally3_d_matrix *x, magma_tally3_d_preconditioner *precond,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_d_initP2P(
    magma_tally3_int_t *bandwidth_benchmark,
    magma_tally3_int_t *num_gpus,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dcompact(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3Double_ptr dA, magma_tally3_int_t ldda,
    double *dnorms, double tol, 
    magma_tally3_int_t *activeMask, magma_tally3_int_t *cBlockSize,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dcompactActive(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3Double_ptr dA, magma_tally3_int_t ldda, 
    magma_tally3_int_t *active,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dmlumerge(    
    magma_tally3_d_matrix L, 
    magma_tally3_d_matrix U,
    magma_tally3_d_matrix *A, 
    magma_tally3_queue_t queue );


/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally3_SPARSE BLAS function definitions
*/
magma_tally3_int_t 
magma_tally3_dgecsrmv(
    magma_tally3_trans_t transA,
    magma_tally3_int_t m, magma_tally3_int_t n,
    double alpha,
    magma_tally3Double_ptr dval,
    magma_tally3Index_ptr drowptr,
    magma_tally3Index_ptr dcolind,
    magma_tally3Double_ptr dx,
    double beta,
    magma_tally3Double_ptr dy,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_dgecsrmv_shift(
    magma_tally3_trans_t transA,
    magma_tally3_int_t m, magma_tally3_int_t n,
    double alpha,
    double lambda,
    magma_tally3Double_ptr dval,
    magma_tally3Index_ptr drowptr,
    magma_tally3Index_ptr dcolind,
    magma_tally3Double_ptr dx,
    double beta,
    int offset,
    int blocksize,
    magma_tally3Index_ptr dadd_rows,
    magma_tally3Double_ptr dy,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_dmgecsrmv(
    magma_tally3_trans_t transA,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3_int_t num_vecs,
    double alpha,
    magma_tally3Double_ptr dval,
    magma_tally3Index_ptr drowptr,
    magma_tally3Index_ptr dcolind,
    magma_tally3Double_ptr dx,
    double beta,
    magma_tally3Double_ptr dy,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_dgeellmv(
    magma_tally3_trans_t transA,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3_int_t nnz_per_row,
    double alpha,
    magma_tally3Double_ptr dval,
    magma_tally3Index_ptr dcolind,
    magma_tally3Double_ptr dx,
    double beta,
    magma_tally3Double_ptr dy,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_dgeellmv_shift(
    magma_tally3_trans_t transA,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3_int_t nnz_per_row,
    double alpha,
    double lambda,
    magma_tally3Double_ptr dval,
    magma_tally3Index_ptr dcolind,
    magma_tally3Double_ptr dx,
    double beta,
    int offset,
    int blocksize,
    magma_tally3Index_ptr dadd_rows,
    magma_tally3Double_ptr dy,
    magma_tally3_queue_t queue );


magma_tally3_int_t 
magma_tally3_dmgeellmv(
    magma_tally3_trans_t transA,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3_int_t num_vecs,
    magma_tally3_int_t nnz_per_row,
    double alpha,
    magma_tally3Double_ptr dval,
    magma_tally3Index_ptr dcolind,
    magma_tally3Double_ptr dx,
    double beta,
    magma_tally3Double_ptr dy,
    magma_tally3_queue_t queue );


magma_tally3_int_t 
magma_tally3_dgeelltmv(
    magma_tally3_trans_t transA,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3_int_t nnz_per_row,
    double alpha,
    magma_tally3Double_ptr dval,
    magma_tally3Index_ptr dcolind,
    magma_tally3Double_ptr dx,
    double beta,
    magma_tally3Double_ptr dy,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_dgeelltmv_shift(
    magma_tally3_trans_t transA,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3_int_t nnz_per_row,
    double alpha,
    double lambda,
    magma_tally3Double_ptr dval,
    magma_tally3Index_ptr dcolind,
    magma_tally3Double_ptr dx,
    double beta,
    int offset,
    int blocksize,
    magma_tally3Index_ptr dadd_rows,
    magma_tally3Double_ptr dy,
    magma_tally3_queue_t queue );


magma_tally3_int_t 
magma_tally3_dmgeelltmv(
    magma_tally3_trans_t transA,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3_int_t num_vecs,
    magma_tally3_int_t nnz_per_row,
    double alpha,
    magma_tally3Double_ptr dval,
    magma_tally3Index_ptr dcolind,
    magma_tally3Double_ptr dx,
    double beta,
    magma_tally3Double_ptr dy,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_dgeellrtmv(
    magma_tally3_trans_t transA,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3_int_t nnz_per_row,
    double alpha,
    magma_tally3Double_ptr dval,
    magma_tally3Index_ptr dcolind,
    magma_tally3Index_ptr drowlength,
    magma_tally3Double_ptr dx,
    double beta,
    magma_tally3Double_ptr dy,
    magma_tally3_int_t num_threads,
    magma_tally3_int_t threads_per_row,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_dgesellcmv(
    magma_tally3_trans_t transA,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3_int_t blocksize,
    magma_tally3_int_t slices,
    magma_tally3_int_t alignment,
    double alpha,
    magma_tally3Double_ptr dval,
    magma_tally3Index_ptr dcolind,
    magma_tally3Index_ptr drowptr,
    magma_tally3Double_ptr dx,
    double beta,
    magma_tally3Double_ptr dy,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dgesellpmv(
    magma_tally3_trans_t transA,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3_int_t blocksize,
    magma_tally3_int_t slices,
    magma_tally3_int_t alignment,
    double alpha,
    magma_tally3Double_ptr dval,
    magma_tally3Index_ptr dcolind,
    magma_tally3Index_ptr drowptr,
    magma_tally3Double_ptr dx,
    double beta,
    magma_tally3Double_ptr dy,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dmgesellpmv(
    magma_tally3_trans_t transA,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3_int_t num_vecs,
    magma_tally3_int_t blocksize,
    magma_tally3_int_t slices,
    magma_tally3_int_t alignment,
    double alpha,
    magma_tally3Double_ptr dval,
    magma_tally3Index_ptr dcolind,
    magma_tally3Index_ptr drowptr,
    magma_tally3Double_ptr dx,
    double beta,
    magma_tally3Double_ptr dy,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dmgesellpmv_blocked(
    magma_tally3_trans_t transA,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3_int_t num_vecs,
    magma_tally3_int_t blocksize,
    magma_tally3_int_t slices,
    magma_tally3_int_t alignment,
    double alpha,
    magma_tally3Double_ptr dval,
    magma_tally3Index_ptr dcolind,
    magma_tally3Index_ptr drowptr,
    magma_tally3Double_ptr dx,
    double beta,
    magma_tally3Double_ptr dy,
    magma_tally3_queue_t queue );


magma_tally3_int_t
magma_tally3_dmergedgs(
    magma_tally3_int_t n, 
    magma_tally3_int_t ldh,
    magma_tally3_int_t k, 
    magma_tally3Double_ptr dv, 
    magma_tally3Double_ptr dr,
    magma_tally3Double_ptr dskp,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dcopyscale(    
    int n, 
    int k,
    magma_tally3Double_ptr dr, 
    magma_tally3Double_ptr dv,
    magma_tally3Double_ptr dskp,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dnrm2scale(    
    int m, 
    magma_tally3Double_ptr dr,    
    int lddr, 
    double *drnorm,
    magma_tally3_queue_t queue );


magma_tally3_int_t
magma_tally3_djacobisetup_vector_gpu(
    int num_rows, 
    magma_tally3_d_matrix b, 
    magma_tally3_d_matrix d, 
    magma_tally3_d_matrix c,
    magma_tally3_d_matrix *x,
    magma_tally3_queue_t queue );


magma_tally3_int_t
magma_tally3_djacobi_diagscal(    
    int num_rows, 
    magma_tally3_d_matrix d, 
    magma_tally3_d_matrix b, 
    magma_tally3_d_matrix *c,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_djacobiupdate(
    magma_tally3_d_matrix t, 
    magma_tally3_d_matrix b, 
    magma_tally3_d_matrix d, 
    magma_tally3_d_matrix *x,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_djacobispmvupdate(
    magma_tally3_int_t maxiter,
    magma_tally3_d_matrix A, 
    magma_tally3_d_matrix t, 
    magma_tally3_d_matrix b, 
    magma_tally3_d_matrix d, 
    magma_tally3_d_matrix *x,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_djacobispmvupdate_bw(
    magma_tally3_int_t maxiter,
    magma_tally3_d_matrix A, 
    magma_tally3_d_matrix t, 
    magma_tally3_d_matrix b, 
    magma_tally3_d_matrix d, 
    magma_tally3_d_matrix *x,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_djacobispmvupdateselect(
    magma_tally3_int_t maxiter,
    magma_tally3_int_t num_updates,
    magma_tally3_index_t *indices,
    magma_tally3_d_matrix A,
    magma_tally3_d_matrix t, 
    magma_tally3_d_matrix b, 
    magma_tally3_d_matrix d, 
    magma_tally3_d_matrix tmp, 
    magma_tally3_d_matrix *x,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_djacobisetup_diagscal(
    magma_tally3_d_matrix A, magma_tally3_d_matrix *d,
    magma_tally3_queue_t queue );


magma_tally3_int_t
magma_tally3_dbicgmerge1(    
    int n, 
    magma_tally3Double_ptr dskp,
    magma_tally3Double_ptr dv, 
    magma_tally3Double_ptr dr, 
    magma_tally3Double_ptr dp,
    magma_tally3_queue_t queue );


magma_tally3_int_t
magma_tally3_dbicgmerge2(
    int n, 
    magma_tally3Double_ptr dskp, 
    magma_tally3Double_ptr dr,
    magma_tally3Double_ptr dv, 
    magma_tally3Double_ptr ds,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dbicgmerge3(
    int n, 
    magma_tally3Double_ptr dskp, 
    magma_tally3Double_ptr dp,
    magma_tally3Double_ptr ds,
    magma_tally3Double_ptr dt,
    magma_tally3Double_ptr dx, 
    magma_tally3Double_ptr dr,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dbicgmerge4(
    int type, 
    magma_tally3Double_ptr dskp,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dcgmerge_spmv1( 
    magma_tally3_d_matrix A,
    magma_tally3Double_ptr d1,
    magma_tally3Double_ptr d2,
    magma_tally3Double_ptr dd,
    magma_tally3Double_ptr dz,
    magma_tally3Double_ptr dskp,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dcgmerge_xrbeta( 
    int n,
    magma_tally3Double_ptr d1,
    magma_tally3Double_ptr d2,
    magma_tally3Double_ptr dx,
    magma_tally3Double_ptr dr,
    magma_tally3Double_ptr dd,
    magma_tally3Double_ptr dz, 
    magma_tally3Double_ptr dskp,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dmdotc(
    magma_tally3_int_t n, 
    magma_tally3_int_t k, 
    magma_tally3Double_ptr dv, 
    magma_tally3Double_ptr dr,
    magma_tally3Double_ptr dd1,
    magma_tally3Double_ptr dd2,
    magma_tally3Double_ptr dskp,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dgemvmdot(
    int n, 
    int k, 
    magma_tally3Double_ptr dv, 
    magma_tally3Double_ptr dr,
    magma_tally3Double_ptr dd1,
    magma_tally3Double_ptr dd2,
    magma_tally3Double_ptr dskp,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dbicgmerge_spmv1( 
    magma_tally3_d_matrix A,
    magma_tally3Double_ptr dd1,
    magma_tally3Double_ptr dd2,
    magma_tally3Double_ptr dp,
    magma_tally3Double_ptr dr,
    magma_tally3Double_ptr dv,
    magma_tally3Double_ptr dskp,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dbicgmerge_spmv2( 
    magma_tally3_d_matrix A,
    magma_tally3Double_ptr dd1,
    magma_tally3Double_ptr dd2,
    magma_tally3Double_ptr ds,
    magma_tally3Double_ptr dt,
    magma_tally3Double_ptr dskp,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dbicgmerge_xrbeta( 
    int n,
    magma_tally3Double_ptr dd1,
    magma_tally3Double_ptr dd2,
    magma_tally3Double_ptr drr,
    magma_tally3Double_ptr dr,
    magma_tally3Double_ptr dp,
    magma_tally3Double_ptr ds,
    magma_tally3Double_ptr dt,
    magma_tally3Double_ptr dx, 
    magma_tally3Double_ptr dskp,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dbcsrswp(
    magma_tally3_int_t n,
    magma_tally3_int_t size_b, 
    magma_tally3_int_t *ipiv,
    magma_tally3Double_ptr dx,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dbcsrtrsv(
    magma_tally3_uplo_t uplo,
    magma_tally3_int_t r_blocks,
    magma_tally3_int_t c_blocks,
    magma_tally3_int_t size_b, 
    magma_tally3Double_ptr dA,
    magma_tally3_index_t *blockinfo, 
    magma_tally3Double_ptr dx,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dbcsrvalcpy(
    magma_tally3_int_t size_b, 
    magma_tally3_int_t num_blocks, 
    magma_tally3_int_t num_zero_blocks, 
    magma_tally3Double_ptr *dAval, 
    magma_tally3Double_ptr *dBval,
    magma_tally3Double_ptr *dBval2,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dbcsrluegemm(
    magma_tally3_int_t size_b, 
    magma_tally3_int_t num_block_rows,
    magma_tally3_int_t kblocks,
    magma_tally3Double_ptr *dA, 
    magma_tally3Double_ptr *dB, 
    magma_tally3Double_ptr *dC,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dbcsrlupivloc(
    magma_tally3_int_t size_b, 
    magma_tally3_int_t kblocks,
    magma_tally3Double_ptr *dA, 
    magma_tally3_int_t *ipiv,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dbcsrblockinfo5(
    magma_tally3_int_t lustep,
    magma_tally3_int_t num_blocks, 
    magma_tally3_int_t c_blocks, 
    magma_tally3_int_t size_b,
    magma_tally3_index_t *blockinfo,
    magma_tally3Double_ptr dval,
    magma_tally3Double_ptr *AII,
    magma_tally3_queue_t queue );

 
#ifdef __cplusplus
}
#endif

#undef PRECISION_d
#endif /* MAGMA_tally3SPARSE_D_H */
