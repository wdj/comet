/*
 -- MAGMA_tally4 (version 1.1) --
 Univ. of Tennessee, Knoxville
 Univ. of California, Berkeley
 Univ. of Colorado, Denver
 @date

 @generated from magma_tally4sparse_z.h normal z -> d, Mon May  4 17:11:25 2015
 @author Hartwig Anzt
*/

#ifndef MAGMA_tally4SPARSE_D_H
#define MAGMA_tally4SPARSE_D_H

#include "magma_tally4_types.h"
#include "magma_tally4sparse_types.h"

#define PRECISION_d


#ifdef __cplusplus
extern "C" {
#endif


/* ////////////////////////////////////////////////////////////////////////////
 -- For backwards compatability, map old (1.6.1) to new (1.6.2) function names
*/

#define magma_tally4_d_mtranspose  magma_tally4_dmtranspose
#define magma_tally4_d_mtransfer   magma_tally4_dmtransfer
#define magma_tally4_d_vtransfer   magma_tally4_dmtransfer
#define magma_tally4_d_mconvert    magma_tally4_dmconvert
#define magma_tally4_d_vinit       magma_tally4_dvinit
#define magma_tally4_d_vvisu       magma_tally4_dprint_vector
#define magma_tally4_d_vread       magma_tally4_dvread
#define magma_tally4_d_vspread     magma_tally4_dvspread
#define magma_tally4_d_mvisu       magma_tally4_dprint_matrix
#define magma_tally4_d_mfree       magma_tally4_dmfree
#define magma_tally4_d_vfree       magma_tally4_dmfree
#define write_d_csr_mtx     magma_tally4_dwrite_csr_mtx
#define write_d_csrtomtx    magma_tally4_dwrite_csrtomtx
#define print_d_csr         magma_tally4_dprint_csr_mtx


/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally4_SPARSE Auxiliary functions
*/


magma_tally4_int_t
magma_tally4_dparse_opts( 
    int argc, 
    char** argv, 
    magma_tally4_dopts *opts, 
    int *matrices, 
    magma_tally4_queue_t queue );

magma_tally4_int_t 
read_d_csr_from_binary( 
    magma_tally4_int_t* n_row, 
    magma_tally4_int_t* n_col, 
    magma_tally4_int_t* nnz, 
    double **val, 
    magma_tally4_index_t **row, 
    magma_tally4_index_t **col,
    const char * filename,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
read_d_csr_from_mtx( 
    magma_tally4_storage_t *type, 
    magma_tally4_location_t *location,
    magma_tally4_int_t* n_row, 
    magma_tally4_int_t* n_col, 
    magma_tally4_int_t* nnz, 
    double **val, 
    magma_tally4_index_t **row, 
    magma_tally4_index_t **col, 
    const char *filename,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_d_csr_mtx( 
    magma_tally4_d_matrix *A, 
    const char *filename,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_dcsrset( 
    magma_tally4_int_t m, 
    magma_tally4_int_t n, 
    magma_tally4_index_t *row, 
    magma_tally4_index_t *col, 
    double *val,
    magma_tally4_d_matrix *A,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_dcsrget( 
    magma_tally4_d_matrix A,
    magma_tally4_int_t *m, 
    magma_tally4_int_t *n, 
    magma_tally4_index_t **row, 
    magma_tally4_index_t **col, 
    double **val,
    magma_tally4_queue_t queue );


magma_tally4_int_t 
magma_tally4_dvset( 
    magma_tally4_int_t m, 
    magma_tally4_int_t n, 
    double *val,
    magma_tally4_d_matrix *v,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_dvget( 
    magma_tally4_d_matrix v,
    magma_tally4_int_t *m, 
    magma_tally4_int_t *n, 
    double **val,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_dvset_dev( 
    magma_tally4_int_t m, 
    magma_tally4_int_t n, 
    magma_tally4Double_ptr val,
    magma_tally4_d_matrix *v,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_dvget_dev( 
    magma_tally4_d_matrix v,
    magma_tally4_int_t *m, 
    magma_tally4_int_t *n, 
    magma_tally4Double_ptr *val,
    magma_tally4_queue_t queue );


magma_tally4_int_t 
magma_tally4_d_csr_mtxsymm( 
    magma_tally4_d_matrix *A, 
    const char *filename,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_d_csr_compressor( 
    double ** val, 
    magma_tally4_index_t ** row, 
    magma_tally4_index_t ** col, 
    double ** valn, 
    magma_tally4_index_t ** rown, 
    magma_tally4_index_t ** coln, 
    magma_tally4_int_t *n,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dmcsrcompressor( 
    magma_tally4_d_matrix *A,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dmcsrcompressor_gpu( 
    magma_tally4_d_matrix *A,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dvtranspose( 
    magma_tally4_d_matrix x,
    magma_tally4_d_matrix *y,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_d_cucsrtranspose( 
    magma_tally4_d_matrix A, 
    magma_tally4_d_matrix *B,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
d_transpose_csr( 
    magma_tally4_int_t n_rows, 
    magma_tally4_int_t n_cols, 
    magma_tally4_int_t nnz,
    double *val, 
    magma_tally4_index_t *row, 
    magma_tally4_index_t *col, 
    magma_tally4_int_t *new_n_rows, 
    magma_tally4_int_t *new_n_cols, 
    magma_tally4_int_t *new_nnz, 
    double **new_val, 
    magma_tally4_index_t **new_row, 
    magma_tally4_index_t **new_col,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dcsrsplit( 
    magma_tally4_int_t bsize,
    magma_tally4_d_matrix A,
    magma_tally4_d_matrix *D,
    magma_tally4_d_matrix *R,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dmscale( 
    magma_tally4_d_matrix *A, 
    magma_tally4_scale_t scaling,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_dmdiff( 
    magma_tally4_d_matrix A, 
    magma_tally4_d_matrix B, 
 real_Double_t *res,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dmdiagadd( 
    magma_tally4_d_matrix *A, 
    double add,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_dmsort( 
    magma_tally4_d_matrix *A,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dindexsort(
    magma_tally4_index_t *x, 
    magma_tally4_int_t first,
    magma_tally4_int_t last,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_ddomainoverlap(
    magma_tally4_index_t num_rows,
    magma_tally4_index_t *num_indices,
    magma_tally4_index_t *rowptr,
    magma_tally4_index_t *colidx,
    magma_tally4_index_t *x,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dsymbilu( 
    magma_tally4_d_matrix *A, 
    magma_tally4_int_t levels,
    magma_tally4_d_matrix *L,
    magma_tally4_d_matrix *U,
    magma_tally4_queue_t queue );


magma_tally4_int_t 
magma_tally4_dwrite_csr_mtx( 
    magma_tally4_d_matrix A,
    magma_tally4_order_t MajorType,
 const char *filename,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_dwrite_csrtomtx( 
    magma_tally4_d_matrix A,
    const char *filename,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_dprint_csr( 
    magma_tally4_int_t n_row, 
    magma_tally4_int_t n_col, 
    magma_tally4_int_t nnz, 
    double **val, 
    magma_tally4_index_t **row, 
    magma_tally4_index_t **col,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_dprint_csr_mtx( 
    magma_tally4_int_t n_row, 
    magma_tally4_int_t n_col, 
    magma_tally4_int_t nnz, 
    double **val, 
    magma_tally4_index_t **row, 
    magma_tally4_index_t **col, 
    magma_tally4_order_t MajorType,
    magma_tally4_queue_t queue );


magma_tally4_int_t 
magma_tally4_dmtranspose(
    magma_tally4_d_matrix A, 
    magma_tally4_d_matrix *B,
    magma_tally4_queue_t queue );


magma_tally4_int_t 
magma_tally4_dmtransfer(
    magma_tally4_d_matrix A, 
    magma_tally4_d_matrix *B, 
    magma_tally4_location_t src, 
    magma_tally4_location_t dst,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_dmconvert(
    magma_tally4_d_matrix A, 
    magma_tally4_d_matrix *B, 
    magma_tally4_storage_t old_format, 
    magma_tally4_storage_t new_format,
    magma_tally4_queue_t queue );


magma_tally4_int_t
magma_tally4_dvinit(
    magma_tally4_d_matrix *x, 
    magma_tally4_location_t memory_location,
    magma_tally4_int_t num_rows, 
    magma_tally4_int_t num_cols,
    double values,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dprint_vector(
    magma_tally4_d_matrix x, 
    magma_tally4_int_t offset, 
    magma_tally4_int_t displaylength,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dvread(
    magma_tally4_d_matrix *x, 
    magma_tally4_int_t length,
    char * filename,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dvspread(
    magma_tally4_d_matrix *x, 
    const char * filename,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dprint_matrix(
    magma_tally4_d_matrix A,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_ddiameter(
    magma_tally4_d_matrix *A,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_drowentries(
    magma_tally4_d_matrix *A,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dmfree(
    magma_tally4_d_matrix *A,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dresidual(
    magma_tally4_d_matrix A, 
    magma_tally4_d_matrix b, 
    magma_tally4_d_matrix x, 
    double *res,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dresidualvec(
    magma_tally4_d_matrix A,
    magma_tally4_d_matrix b,
    magma_tally4_d_matrix x,
    magma_tally4_d_matrix *r,
    double *res,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dmgenerator(
    magma_tally4_int_t n,
    magma_tally4_int_t offdiags,
    magma_tally4_index_t *diag_offset,
    double *diag_vals,
    magma_tally4_d_matrix *A,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dm_27stencil(
    magma_tally4_int_t n,
    magma_tally4_d_matrix *A,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dm_5stencil(
    magma_tally4_int_t n,
    magma_tally4_d_matrix *A,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dsolverinfo(
    magma_tally4_d_solver_par *solver_par, 
    magma_tally4_d_preconditioner *precond_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dsolverinfo_init(
    magma_tally4_d_solver_par *solver_par, 
    magma_tally4_d_preconditioner *precond,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_deigensolverinfo_init(
    magma_tally4_d_solver_par *solver_par,
    magma_tally4_queue_t queue );


magma_tally4_int_t
magma_tally4_dsolverinfo_free(
    magma_tally4_d_solver_par *solver_par, 
    magma_tally4_d_preconditioner *precond,
    magma_tally4_queue_t queue );



/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally4_SPARSE iterative incomplete factorizations
*/


magma_tally4_int_t
magma_tally4_diterilusetup( 
    magma_tally4_d_matrix A, 
    magma_tally4_d_matrix b,                                 
    magma_tally4_d_preconditioner *precond,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_ditericsetup( 
    magma_tally4_d_matrix A, 
    magma_tally4_d_matrix b, 
    magma_tally4_d_preconditioner *precond,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_ditericupdate( 
    magma_tally4_d_matrix A, 
    magma_tally4_d_preconditioner *precond, 
    magma_tally4_int_t updates,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dapplyiteric_l( 
    magma_tally4_d_matrix b, 
    magma_tally4_d_matrix *x, 
    magma_tally4_d_preconditioner *precond,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dapplyiteric_r( 
    magma_tally4_d_matrix b, 
    magma_tally4_d_matrix *x, 
    magma_tally4_d_preconditioner *precond,
    magma_tally4_queue_t queue );


magma_tally4_int_t
magma_tally4_diterilu_csr( 
    magma_tally4_d_matrix A,
    magma_tally4_d_matrix L,
    magma_tally4_d_matrix U,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_diteric_csr( 
    magma_tally4_d_matrix A,
    magma_tally4_d_matrix A_CSR,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_dfrobenius( 
    magma_tally4_d_matrix A, 
    magma_tally4_d_matrix B, 
    real_Double_t *res,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_dnonlinres(   
    magma_tally4_d_matrix A, 
    magma_tally4_d_matrix L,
    magma_tally4_d_matrix U, 
    magma_tally4_d_matrix *LU, 
    real_Double_t *res,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_dilures(   
    magma_tally4_d_matrix A, 
    magma_tally4_d_matrix L,
    magma_tally4_d_matrix U, 
    magma_tally4_d_matrix *LU, 
    real_Double_t *res,
    real_Double_t *nonlinres,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_dicres(       
    magma_tally4_d_matrix A, 
    magma_tally4_d_matrix C,
    magma_tally4_d_matrix CT, 
    magma_tally4_d_matrix *LU, 
    real_Double_t *res,
    real_Double_t *nonlinres,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_dinitguess( 
    magma_tally4_d_matrix A, 
    magma_tally4_d_matrix *L, 
    magma_tally4_d_matrix *U,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_dinitrecursiveLU( 
    magma_tally4_d_matrix A, 
    magma_tally4_d_matrix *B,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_dmLdiagadd( 
    magma_tally4_d_matrix *L,
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
magma_tally4_dcg(
    magma_tally4_d_matrix A, magma_tally4_d_matrix b, 
    magma_tally4_d_matrix *x, magma_tally4_d_solver_par *solver_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_dcg_res(
    magma_tally4_d_matrix A, magma_tally4_d_matrix b, 
    magma_tally4_d_matrix *x, magma_tally4_d_solver_par *solver_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_dcg_merge(
    magma_tally4_d_matrix A, magma_tally4_d_matrix b, 
    magma_tally4_d_matrix *x, magma_tally4_d_solver_par *solver_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_dgmres(
    magma_tally4_d_matrix A, magma_tally4_d_matrix b, 
    magma_tally4_d_matrix *x, magma_tally4_d_solver_par *solver_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dbicgstab(
    magma_tally4_d_matrix A, magma_tally4_d_matrix b, magma_tally4_d_matrix *x, 
    magma_tally4_d_solver_par *solver_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dbicgstab_merge(
    magma_tally4_d_matrix A, magma_tally4_d_matrix b, 
    magma_tally4_d_matrix *x, magma_tally4_d_solver_par *solver_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dbicgstab_merge2(
    magma_tally4_d_matrix A, magma_tally4_d_matrix b, 
    magma_tally4_d_matrix *x, magma_tally4_d_solver_par *solver_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dpcg(
    magma_tally4_d_matrix A, magma_tally4_d_matrix b, 
    magma_tally4_d_matrix *x, magma_tally4_d_solver_par *solver_par, 
    magma_tally4_d_preconditioner *precond_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dbpcg(
    magma_tally4_d_matrix A, magma_tally4_d_matrix b, 
    magma_tally4_d_matrix *x, magma_tally4_d_solver_par *solver_par, 
    magma_tally4_d_preconditioner *precond_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dpbicgstab(
    magma_tally4_d_matrix A, magma_tally4_d_matrix b, 
    magma_tally4_d_matrix *x, magma_tally4_d_solver_par *solver_par, 
    magma_tally4_d_preconditioner *precond_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dpgmres(
    magma_tally4_d_matrix A, magma_tally4_d_matrix b, 
    magma_tally4_d_matrix *x, magma_tally4_d_solver_par *solver_par, 
    magma_tally4_d_preconditioner *precond_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dfgmres(
    magma_tally4_d_matrix A, magma_tally4_d_matrix b, 
    magma_tally4_d_matrix *x, magma_tally4_d_solver_par *solver_par, 
    magma_tally4_d_preconditioner *precond_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_djacobi(
    magma_tally4_d_matrix A, magma_tally4_d_matrix b, 
    magma_tally4_d_matrix *x, magma_tally4_d_solver_par *solver_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_djacobidomainoverlap(
    magma_tally4_d_matrix A, 
    magma_tally4_d_matrix b, 
    magma_tally4_d_matrix *x,  
    magma_tally4_d_solver_par *solver_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dbaiter(
    magma_tally4_d_matrix A, magma_tally4_d_matrix b, 
    magma_tally4_d_matrix *x, magma_tally4_d_solver_par *solver_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_diterref(
    magma_tally4_d_matrix A, magma_tally4_d_matrix b, 
    magma_tally4_d_matrix *x, magma_tally4_d_solver_par *solver_par, 
    magma_tally4_d_preconditioner *precond_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dilu(
    magma_tally4_d_matrix A, magma_tally4_d_matrix b, 
    magma_tally4_d_matrix *x, magma_tally4_d_solver_par *solver_par, 
    magma_tally4_int_t *ipiv,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dbcsrlu(
    magma_tally4_d_matrix A, magma_tally4_d_matrix b, 
    magma_tally4_d_matrix *x, magma_tally4_d_solver_par *solver_par,
    magma_tally4_queue_t queue );


magma_tally4_int_t
magma_tally4_dbcsrlutrf(
    magma_tally4_d_matrix A, 
    magma_tally4_d_matrix *M,
    magma_tally4_int_t *ipiv, 
    magma_tally4_int_t version,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dbcsrlusv(
    magma_tally4_d_matrix A, magma_tally4_d_matrix b, 
    magma_tally4_d_matrix *x, magma_tally4_d_solver_par *solver_par, 
    magma_tally4_int_t *ipiv,
    magma_tally4_queue_t queue );



magma_tally4_int_t
magma_tally4_dilucg(
    magma_tally4_d_matrix A, magma_tally4_d_matrix b, 
    magma_tally4_d_matrix *x, magma_tally4_d_solver_par *solver_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dilugmres(
    magma_tally4_d_matrix A, magma_tally4_d_matrix b, 
    magma_tally4_d_matrix *x, magma_tally4_d_solver_par *solver_par,
    magma_tally4_queue_t queue ); 


magma_tally4_int_t
magma_tally4_dlobpcg_shift(
    magma_tally4_int_t num_rows,
    magma_tally4_int_t num_vecs, 
    magma_tally4_int_t shift,
    magma_tally4Double_ptr x,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dlobpcg_res(
    magma_tally4_int_t num_rows,
    magma_tally4_int_t num_vecs, 
    double *evalues, 
    magma_tally4Double_ptr X,
    magma_tally4Double_ptr R, 
    double *res,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dlobpcg_maxpy(
    magma_tally4_int_t num_rows,
    magma_tally4_int_t num_vecs, 
    magma_tally4Double_ptr X,
    magma_tally4Double_ptr Y,
    magma_tally4_queue_t queue );


/*/////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally4_SPARSE eigensolvers (Data on GPU)
*/
magma_tally4_int_t
magma_tally4_dlobpcg(
    magma_tally4_d_matrix A, 
    magma_tally4_d_solver_par *solver_par,
    magma_tally4_d_preconditioner *precond_par, 
    magma_tally4_queue_t queue );




/*/////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally4_SPARSE preconditioners (Data on GPU)
*/

magma_tally4_int_t
magma_tally4_djacobisetup(
    magma_tally4_d_matrix A, 
    magma_tally4_d_matrix b, 
    magma_tally4_d_matrix *M, 
    magma_tally4_d_matrix *c,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_djacobisetup_matrix(
    magma_tally4_d_matrix A, 
    magma_tally4_d_matrix *M, 
    magma_tally4_d_matrix *d,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_djacobisetup_vector(
    magma_tally4_d_matrix b, 
    magma_tally4_d_matrix d, 
    magma_tally4_d_matrix *c,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_djacobiiter(
    magma_tally4_d_matrix M, 
    magma_tally4_d_matrix c, 
    magma_tally4_d_matrix *x, 
    magma_tally4_d_solver_par *solver_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_djacobiiter_precond( 
    magma_tally4_d_matrix M, 
    magma_tally4_d_matrix *x, 
    magma_tally4_d_solver_par *solver_par, 
    magma_tally4_d_preconditioner *precond,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_djacobiiter_sys(
    magma_tally4_d_matrix A, 
    magma_tally4_d_matrix b, 
    magma_tally4_d_matrix d, 
    magma_tally4_d_matrix t, 
    magma_tally4_d_matrix *x,  
    magma_tally4_d_solver_par *solver_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dpastixsetup(
    magma_tally4_d_matrix A, magma_tally4_d_matrix b,
    magma_tally4_d_preconditioner *precond,
    magma_tally4_queue_t queue );


magma_tally4_int_t
magma_tally4_dapplypastix(
    magma_tally4_d_matrix b, magma_tally4_d_matrix *x, 
    magma_tally4_d_preconditioner *precond,
    magma_tally4_queue_t queue );


// custom preconditioner
magma_tally4_int_t
magma_tally4_dapplycustomprecond_l(
    magma_tally4_d_matrix b, 
    magma_tally4_d_matrix *x, 
    magma_tally4_d_preconditioner *precond,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dapplycustomprecond_r(
    magma_tally4_d_matrix b, 
    magma_tally4_d_matrix *x, 
    magma_tally4_d_preconditioner *precond,
    magma_tally4_queue_t queue );


// CUSPARSE preconditioner

magma_tally4_int_t
magma_tally4_dcuilusetup(
    magma_tally4_d_matrix A, magma_tally4_d_preconditioner *precond,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dapplycuilu_l(
    magma_tally4_d_matrix b, magma_tally4_d_matrix *x, 
    magma_tally4_d_preconditioner *precond,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dapplycuilu_r(
    magma_tally4_d_matrix b, magma_tally4_d_matrix *x, 
    magma_tally4_d_preconditioner *precond,
    magma_tally4_queue_t queue );


magma_tally4_int_t
magma_tally4_dcuiccsetup(
    magma_tally4_d_matrix A, magma_tally4_d_preconditioner *precond,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dapplycuicc_l(
    magma_tally4_d_matrix b, magma_tally4_d_matrix *x, 
    magma_tally4_d_preconditioner *precond,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dapplycuicc_r(
    magma_tally4_d_matrix b, magma_tally4_d_matrix *x, 
    magma_tally4_d_preconditioner *precond,
    magma_tally4_queue_t queue );


magma_tally4_int_t
magma_tally4_dcumilusetup(
    magma_tally4_d_matrix A, 
    magma_tally4_d_preconditioner *precond,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dcumilugeneratesolverinfo(
    magma_tally4_d_preconditioner *precond,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dapplycumilu_l(
    magma_tally4_d_matrix b, 
    magma_tally4_d_matrix *x, 
    magma_tally4_d_preconditioner *precond,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dapplycumilu_r(
    magma_tally4_d_matrix b, 
    magma_tally4_d_matrix *x, 
    magma_tally4_d_preconditioner *precond,
    magma_tally4_queue_t queue );


magma_tally4_int_t
magma_tally4_dcumiccsetup(
    magma_tally4_d_matrix A, 
    magma_tally4_d_preconditioner *precond,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dcumicgeneratesolverinfo(
    magma_tally4_d_preconditioner *precond,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dapplycumicc_l(
    magma_tally4_d_matrix b, 
    magma_tally4_d_matrix *x, 
    magma_tally4_d_preconditioner *precond,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dapplycumicc_r(
    magma_tally4_d_matrix b, 
    magma_tally4_d_matrix *x, 
    magma_tally4_d_preconditioner *precond,
    magma_tally4_queue_t queue );


// block-asynchronous iteration

magma_tally4_int_t
magma_tally4_dbajac_csr(
    magma_tally4_int_t localiters,
    magma_tally4_d_matrix D,
    magma_tally4_d_matrix R,
    magma_tally4_d_matrix b,
    magma_tally4_d_matrix *x,
    magma_tally4_queue_t queue );

/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally4_SPARSE utility function definitions
*/

magma_tally4_int_t
magma_tally4_d_spmv(
    double alpha, 
    magma_tally4_d_matrix A, 
    magma_tally4_d_matrix x, 
    double beta, 
    magma_tally4_d_matrix y,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dcustomspmv(
    double alpha, 
    magma_tally4_d_matrix x, 
    double beta, 
    magma_tally4_d_matrix y,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_d_spmv_shift(
    double alpha, 
    magma_tally4_d_matrix A, 
    double lambda,
    magma_tally4_d_matrix x, 
    double beta, 
    magma_tally4_int_t offset, 
    magma_tally4_int_t blocksize,
    magma_tally4Index_ptr dadd_vecs, 
    magma_tally4_d_matrix y,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dcuspmm(
    magma_tally4_d_matrix A, 
    magma_tally4_d_matrix B, 
    magma_tally4_d_matrix *AB,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_d_spmm(
    double alpha, 
    magma_tally4_d_matrix A,
    magma_tally4_d_matrix B,
    magma_tally4_d_matrix *C,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dsymbilu( 
    magma_tally4_d_matrix *A, 
    magma_tally4_int_t levels, 
    magma_tally4_d_matrix *L, 
    magma_tally4_d_matrix *U,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dcuspaxpy(
    magma_tally4Double_ptr alpha, magma_tally4_d_matrix A, 
    magma_tally4Double_ptr beta, magma_tally4_d_matrix B, 
    magma_tally4_d_matrix *AB,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_d_precond(
    magma_tally4_d_matrix A, 
    magma_tally4_d_matrix b, magma_tally4_d_matrix *x,
    magma_tally4_d_preconditioner *precond,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_d_solver(
    magma_tally4_d_matrix A, magma_tally4_d_matrix b, 
    magma_tally4_d_matrix *x, magma_tally4_dopts *zopts,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_d_precondsetup(
    magma_tally4_d_matrix A, magma_tally4_d_matrix b, 
    magma_tally4_d_preconditioner *precond,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_d_applyprecond(
    magma_tally4_d_matrix A, magma_tally4_d_matrix b, 
    magma_tally4_d_matrix *x, magma_tally4_d_preconditioner *precond,
    magma_tally4_queue_t queue );


magma_tally4_int_t
magma_tally4_d_applyprecond_left(
    magma_tally4_d_matrix A, magma_tally4_d_matrix b, 
    magma_tally4_d_matrix *x, magma_tally4_d_preconditioner *precond,
    magma_tally4_queue_t queue );


magma_tally4_int_t
magma_tally4_d_applyprecond_right(
    magma_tally4_d_matrix A, magma_tally4_d_matrix b, 
    magma_tally4_d_matrix *x, magma_tally4_d_preconditioner *precond,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_d_initP2P(
    magma_tally4_int_t *bandwidth_benchmark,
    magma_tally4_int_t *num_gpus,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dcompact(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4Double_ptr dA, magma_tally4_int_t ldda,
    double *dnorms, double tol, 
    magma_tally4_int_t *activeMask, magma_tally4_int_t *cBlockSize,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dcompactActive(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4Double_ptr dA, magma_tally4_int_t ldda, 
    magma_tally4_int_t *active,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dmlumerge(    
    magma_tally4_d_matrix L, 
    magma_tally4_d_matrix U,
    magma_tally4_d_matrix *A, 
    magma_tally4_queue_t queue );


/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally4_SPARSE BLAS function definitions
*/
magma_tally4_int_t 
magma_tally4_dgecsrmv(
    magma_tally4_trans_t transA,
    magma_tally4_int_t m, magma_tally4_int_t n,
    double alpha,
    magma_tally4Double_ptr dval,
    magma_tally4Index_ptr drowptr,
    magma_tally4Index_ptr dcolind,
    magma_tally4Double_ptr dx,
    double beta,
    magma_tally4Double_ptr dy,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_dgecsrmv_shift(
    magma_tally4_trans_t transA,
    magma_tally4_int_t m, magma_tally4_int_t n,
    double alpha,
    double lambda,
    magma_tally4Double_ptr dval,
    magma_tally4Index_ptr drowptr,
    magma_tally4Index_ptr dcolind,
    magma_tally4Double_ptr dx,
    double beta,
    int offset,
    int blocksize,
    magma_tally4Index_ptr dadd_rows,
    magma_tally4Double_ptr dy,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_dmgecsrmv(
    magma_tally4_trans_t transA,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4_int_t num_vecs,
    double alpha,
    magma_tally4Double_ptr dval,
    magma_tally4Index_ptr drowptr,
    magma_tally4Index_ptr dcolind,
    magma_tally4Double_ptr dx,
    double beta,
    magma_tally4Double_ptr dy,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_dgeellmv(
    magma_tally4_trans_t transA,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4_int_t nnz_per_row,
    double alpha,
    magma_tally4Double_ptr dval,
    magma_tally4Index_ptr dcolind,
    magma_tally4Double_ptr dx,
    double beta,
    magma_tally4Double_ptr dy,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_dgeellmv_shift(
    magma_tally4_trans_t transA,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4_int_t nnz_per_row,
    double alpha,
    double lambda,
    magma_tally4Double_ptr dval,
    magma_tally4Index_ptr dcolind,
    magma_tally4Double_ptr dx,
    double beta,
    int offset,
    int blocksize,
    magma_tally4Index_ptr dadd_rows,
    magma_tally4Double_ptr dy,
    magma_tally4_queue_t queue );


magma_tally4_int_t 
magma_tally4_dmgeellmv(
    magma_tally4_trans_t transA,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4_int_t num_vecs,
    magma_tally4_int_t nnz_per_row,
    double alpha,
    magma_tally4Double_ptr dval,
    magma_tally4Index_ptr dcolind,
    magma_tally4Double_ptr dx,
    double beta,
    magma_tally4Double_ptr dy,
    magma_tally4_queue_t queue );


magma_tally4_int_t 
magma_tally4_dgeelltmv(
    magma_tally4_trans_t transA,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4_int_t nnz_per_row,
    double alpha,
    magma_tally4Double_ptr dval,
    magma_tally4Index_ptr dcolind,
    magma_tally4Double_ptr dx,
    double beta,
    magma_tally4Double_ptr dy,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_dgeelltmv_shift(
    magma_tally4_trans_t transA,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4_int_t nnz_per_row,
    double alpha,
    double lambda,
    magma_tally4Double_ptr dval,
    magma_tally4Index_ptr dcolind,
    magma_tally4Double_ptr dx,
    double beta,
    int offset,
    int blocksize,
    magma_tally4Index_ptr dadd_rows,
    magma_tally4Double_ptr dy,
    magma_tally4_queue_t queue );


magma_tally4_int_t 
magma_tally4_dmgeelltmv(
    magma_tally4_trans_t transA,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4_int_t num_vecs,
    magma_tally4_int_t nnz_per_row,
    double alpha,
    magma_tally4Double_ptr dval,
    magma_tally4Index_ptr dcolind,
    magma_tally4Double_ptr dx,
    double beta,
    magma_tally4Double_ptr dy,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_dgeellrtmv(
    magma_tally4_trans_t transA,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4_int_t nnz_per_row,
    double alpha,
    magma_tally4Double_ptr dval,
    magma_tally4Index_ptr dcolind,
    magma_tally4Index_ptr drowlength,
    magma_tally4Double_ptr dx,
    double beta,
    magma_tally4Double_ptr dy,
    magma_tally4_int_t num_threads,
    magma_tally4_int_t threads_per_row,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_dgesellcmv(
    magma_tally4_trans_t transA,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4_int_t blocksize,
    magma_tally4_int_t slices,
    magma_tally4_int_t alignment,
    double alpha,
    magma_tally4Double_ptr dval,
    magma_tally4Index_ptr dcolind,
    magma_tally4Index_ptr drowptr,
    magma_tally4Double_ptr dx,
    double beta,
    magma_tally4Double_ptr dy,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dgesellpmv(
    magma_tally4_trans_t transA,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4_int_t blocksize,
    magma_tally4_int_t slices,
    magma_tally4_int_t alignment,
    double alpha,
    magma_tally4Double_ptr dval,
    magma_tally4Index_ptr dcolind,
    magma_tally4Index_ptr drowptr,
    magma_tally4Double_ptr dx,
    double beta,
    magma_tally4Double_ptr dy,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dmgesellpmv(
    magma_tally4_trans_t transA,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4_int_t num_vecs,
    magma_tally4_int_t blocksize,
    magma_tally4_int_t slices,
    magma_tally4_int_t alignment,
    double alpha,
    magma_tally4Double_ptr dval,
    magma_tally4Index_ptr dcolind,
    magma_tally4Index_ptr drowptr,
    magma_tally4Double_ptr dx,
    double beta,
    magma_tally4Double_ptr dy,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dmgesellpmv_blocked(
    magma_tally4_trans_t transA,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4_int_t num_vecs,
    magma_tally4_int_t blocksize,
    magma_tally4_int_t slices,
    magma_tally4_int_t alignment,
    double alpha,
    magma_tally4Double_ptr dval,
    magma_tally4Index_ptr dcolind,
    magma_tally4Index_ptr drowptr,
    magma_tally4Double_ptr dx,
    double beta,
    magma_tally4Double_ptr dy,
    magma_tally4_queue_t queue );


magma_tally4_int_t
magma_tally4_dmergedgs(
    magma_tally4_int_t n, 
    magma_tally4_int_t ldh,
    magma_tally4_int_t k, 
    magma_tally4Double_ptr dv, 
    magma_tally4Double_ptr dr,
    magma_tally4Double_ptr dskp,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dcopyscale(    
    int n, 
    int k,
    magma_tally4Double_ptr dr, 
    magma_tally4Double_ptr dv,
    magma_tally4Double_ptr dskp,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dnrm2scale(    
    int m, 
    magma_tally4Double_ptr dr,    
    int lddr, 
    double *drnorm,
    magma_tally4_queue_t queue );


magma_tally4_int_t
magma_tally4_djacobisetup_vector_gpu(
    int num_rows, 
    magma_tally4_d_matrix b, 
    magma_tally4_d_matrix d, 
    magma_tally4_d_matrix c,
    magma_tally4_d_matrix *x,
    magma_tally4_queue_t queue );


magma_tally4_int_t
magma_tally4_djacobi_diagscal(    
    int num_rows, 
    magma_tally4_d_matrix d, 
    magma_tally4_d_matrix b, 
    magma_tally4_d_matrix *c,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_djacobiupdate(
    magma_tally4_d_matrix t, 
    magma_tally4_d_matrix b, 
    magma_tally4_d_matrix d, 
    magma_tally4_d_matrix *x,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_djacobispmvupdate(
    magma_tally4_int_t maxiter,
    magma_tally4_d_matrix A, 
    magma_tally4_d_matrix t, 
    magma_tally4_d_matrix b, 
    magma_tally4_d_matrix d, 
    magma_tally4_d_matrix *x,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_djacobispmvupdate_bw(
    magma_tally4_int_t maxiter,
    magma_tally4_d_matrix A, 
    magma_tally4_d_matrix t, 
    magma_tally4_d_matrix b, 
    magma_tally4_d_matrix d, 
    magma_tally4_d_matrix *x,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_djacobispmvupdateselect(
    magma_tally4_int_t maxiter,
    magma_tally4_int_t num_updates,
    magma_tally4_index_t *indices,
    magma_tally4_d_matrix A,
    magma_tally4_d_matrix t, 
    magma_tally4_d_matrix b, 
    magma_tally4_d_matrix d, 
    magma_tally4_d_matrix tmp, 
    magma_tally4_d_matrix *x,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_djacobisetup_diagscal(
    magma_tally4_d_matrix A, magma_tally4_d_matrix *d,
    magma_tally4_queue_t queue );


magma_tally4_int_t
magma_tally4_dbicgmerge1(    
    int n, 
    magma_tally4Double_ptr dskp,
    magma_tally4Double_ptr dv, 
    magma_tally4Double_ptr dr, 
    magma_tally4Double_ptr dp,
    magma_tally4_queue_t queue );


magma_tally4_int_t
magma_tally4_dbicgmerge2(
    int n, 
    magma_tally4Double_ptr dskp, 
    magma_tally4Double_ptr dr,
    magma_tally4Double_ptr dv, 
    magma_tally4Double_ptr ds,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dbicgmerge3(
    int n, 
    magma_tally4Double_ptr dskp, 
    magma_tally4Double_ptr dp,
    magma_tally4Double_ptr ds,
    magma_tally4Double_ptr dt,
    magma_tally4Double_ptr dx, 
    magma_tally4Double_ptr dr,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dbicgmerge4(
    int type, 
    magma_tally4Double_ptr dskp,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dcgmerge_spmv1( 
    magma_tally4_d_matrix A,
    magma_tally4Double_ptr d1,
    magma_tally4Double_ptr d2,
    magma_tally4Double_ptr dd,
    magma_tally4Double_ptr dz,
    magma_tally4Double_ptr dskp,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dcgmerge_xrbeta( 
    int n,
    magma_tally4Double_ptr d1,
    magma_tally4Double_ptr d2,
    magma_tally4Double_ptr dx,
    magma_tally4Double_ptr dr,
    magma_tally4Double_ptr dd,
    magma_tally4Double_ptr dz, 
    magma_tally4Double_ptr dskp,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dmdotc(
    magma_tally4_int_t n, 
    magma_tally4_int_t k, 
    magma_tally4Double_ptr dv, 
    magma_tally4Double_ptr dr,
    magma_tally4Double_ptr dd1,
    magma_tally4Double_ptr dd2,
    magma_tally4Double_ptr dskp,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dgemvmdot(
    int n, 
    int k, 
    magma_tally4Double_ptr dv, 
    magma_tally4Double_ptr dr,
    magma_tally4Double_ptr dd1,
    magma_tally4Double_ptr dd2,
    magma_tally4Double_ptr dskp,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dbicgmerge_spmv1( 
    magma_tally4_d_matrix A,
    magma_tally4Double_ptr dd1,
    magma_tally4Double_ptr dd2,
    magma_tally4Double_ptr dp,
    magma_tally4Double_ptr dr,
    magma_tally4Double_ptr dv,
    magma_tally4Double_ptr dskp,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dbicgmerge_spmv2( 
    magma_tally4_d_matrix A,
    magma_tally4Double_ptr dd1,
    magma_tally4Double_ptr dd2,
    magma_tally4Double_ptr ds,
    magma_tally4Double_ptr dt,
    magma_tally4Double_ptr dskp,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dbicgmerge_xrbeta( 
    int n,
    magma_tally4Double_ptr dd1,
    magma_tally4Double_ptr dd2,
    magma_tally4Double_ptr drr,
    magma_tally4Double_ptr dr,
    magma_tally4Double_ptr dp,
    magma_tally4Double_ptr ds,
    magma_tally4Double_ptr dt,
    magma_tally4Double_ptr dx, 
    magma_tally4Double_ptr dskp,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dbcsrswp(
    magma_tally4_int_t n,
    magma_tally4_int_t size_b, 
    magma_tally4_int_t *ipiv,
    magma_tally4Double_ptr dx,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dbcsrtrsv(
    magma_tally4_uplo_t uplo,
    magma_tally4_int_t r_blocks,
    magma_tally4_int_t c_blocks,
    magma_tally4_int_t size_b, 
    magma_tally4Double_ptr dA,
    magma_tally4_index_t *blockinfo, 
    magma_tally4Double_ptr dx,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dbcsrvalcpy(
    magma_tally4_int_t size_b, 
    magma_tally4_int_t num_blocks, 
    magma_tally4_int_t num_zero_blocks, 
    magma_tally4Double_ptr *dAval, 
    magma_tally4Double_ptr *dBval,
    magma_tally4Double_ptr *dBval2,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dbcsrluegemm(
    magma_tally4_int_t size_b, 
    magma_tally4_int_t num_block_rows,
    magma_tally4_int_t kblocks,
    magma_tally4Double_ptr *dA, 
    magma_tally4Double_ptr *dB, 
    magma_tally4Double_ptr *dC,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dbcsrlupivloc(
    magma_tally4_int_t size_b, 
    magma_tally4_int_t kblocks,
    magma_tally4Double_ptr *dA, 
    magma_tally4_int_t *ipiv,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dbcsrblockinfo5(
    magma_tally4_int_t lustep,
    magma_tally4_int_t num_blocks, 
    magma_tally4_int_t c_blocks, 
    magma_tally4_int_t size_b,
    magma_tally4_index_t *blockinfo,
    magma_tally4Double_ptr dval,
    magma_tally4Double_ptr *AII,
    magma_tally4_queue_t queue );

 
#ifdef __cplusplus
}
#endif

#undef PRECISION_d
#endif /* MAGMA_tally4SPARSE_D_H */
