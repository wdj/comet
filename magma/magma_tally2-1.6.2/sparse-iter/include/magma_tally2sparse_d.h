/*
 -- MAGMA_tally2 (version 1.1) --
 Univ. of Tennessee, Knoxville
 Univ. of California, Berkeley
 Univ. of Colorado, Denver
 @date

 @generated from magma_tally2sparse_z.h normal z -> d, Mon May  4 17:11:25 2015
 @author Hartwig Anzt
*/

#ifndef MAGMA_tally2SPARSE_D_H
#define MAGMA_tally2SPARSE_D_H

#include "magma_tally2_types.h"
#include "magma_tally2sparse_types.h"

#define PRECISION_d


#ifdef __cplusplus
extern "C" {
#endif


/* ////////////////////////////////////////////////////////////////////////////
 -- For backwards compatability, map old (1.6.1) to new (1.6.2) function names
*/

#define magma_tally2_d_mtranspose  magma_tally2_dmtranspose
#define magma_tally2_d_mtransfer   magma_tally2_dmtransfer
#define magma_tally2_d_vtransfer   magma_tally2_dmtransfer
#define magma_tally2_d_mconvert    magma_tally2_dmconvert
#define magma_tally2_d_vinit       magma_tally2_dvinit
#define magma_tally2_d_vvisu       magma_tally2_dprint_vector
#define magma_tally2_d_vread       magma_tally2_dvread
#define magma_tally2_d_vspread     magma_tally2_dvspread
#define magma_tally2_d_mvisu       magma_tally2_dprint_matrix
#define magma_tally2_d_mfree       magma_tally2_dmfree
#define magma_tally2_d_vfree       magma_tally2_dmfree
#define write_d_csr_mtx     magma_tally2_dwrite_csr_mtx
#define write_d_csrtomtx    magma_tally2_dwrite_csrtomtx
#define print_d_csr         magma_tally2_dprint_csr_mtx


/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally2_SPARSE Auxiliary functions
*/


magma_tally2_int_t
magma_tally2_dparse_opts( 
    int argc, 
    char** argv, 
    magma_tally2_dopts *opts, 
    int *matrices, 
    magma_tally2_queue_t queue );

magma_tally2_int_t 
read_d_csr_from_binary( 
    magma_tally2_int_t* n_row, 
    magma_tally2_int_t* n_col, 
    magma_tally2_int_t* nnz, 
    double **val, 
    magma_tally2_index_t **row, 
    magma_tally2_index_t **col,
    const char * filename,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
read_d_csr_from_mtx( 
    magma_tally2_storage_t *type, 
    magma_tally2_location_t *location,
    magma_tally2_int_t* n_row, 
    magma_tally2_int_t* n_col, 
    magma_tally2_int_t* nnz, 
    double **val, 
    magma_tally2_index_t **row, 
    magma_tally2_index_t **col, 
    const char *filename,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_d_csr_mtx( 
    magma_tally2_d_matrix *A, 
    const char *filename,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_dcsrset( 
    magma_tally2_int_t m, 
    magma_tally2_int_t n, 
    magma_tally2_index_t *row, 
    magma_tally2_index_t *col, 
    double *val,
    magma_tally2_d_matrix *A,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_dcsrget( 
    magma_tally2_d_matrix A,
    magma_tally2_int_t *m, 
    magma_tally2_int_t *n, 
    magma_tally2_index_t **row, 
    magma_tally2_index_t **col, 
    double **val,
    magma_tally2_queue_t queue );


magma_tally2_int_t 
magma_tally2_dvset( 
    magma_tally2_int_t m, 
    magma_tally2_int_t n, 
    double *val,
    magma_tally2_d_matrix *v,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_dvget( 
    magma_tally2_d_matrix v,
    magma_tally2_int_t *m, 
    magma_tally2_int_t *n, 
    double **val,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_dvset_dev( 
    magma_tally2_int_t m, 
    magma_tally2_int_t n, 
    magma_tally2Double_ptr val,
    magma_tally2_d_matrix *v,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_dvget_dev( 
    magma_tally2_d_matrix v,
    magma_tally2_int_t *m, 
    magma_tally2_int_t *n, 
    magma_tally2Double_ptr *val,
    magma_tally2_queue_t queue );


magma_tally2_int_t 
magma_tally2_d_csr_mtxsymm( 
    magma_tally2_d_matrix *A, 
    const char *filename,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_d_csr_compressor( 
    double ** val, 
    magma_tally2_index_t ** row, 
    magma_tally2_index_t ** col, 
    double ** valn, 
    magma_tally2_index_t ** rown, 
    magma_tally2_index_t ** coln, 
    magma_tally2_int_t *n,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dmcsrcompressor( 
    magma_tally2_d_matrix *A,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dmcsrcompressor_gpu( 
    magma_tally2_d_matrix *A,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dvtranspose( 
    magma_tally2_d_matrix x,
    magma_tally2_d_matrix *y,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_d_cucsrtranspose( 
    magma_tally2_d_matrix A, 
    magma_tally2_d_matrix *B,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
d_transpose_csr( 
    magma_tally2_int_t n_rows, 
    magma_tally2_int_t n_cols, 
    magma_tally2_int_t nnz,
    double *val, 
    magma_tally2_index_t *row, 
    magma_tally2_index_t *col, 
    magma_tally2_int_t *new_n_rows, 
    magma_tally2_int_t *new_n_cols, 
    magma_tally2_int_t *new_nnz, 
    double **new_val, 
    magma_tally2_index_t **new_row, 
    magma_tally2_index_t **new_col,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dcsrsplit( 
    magma_tally2_int_t bsize,
    magma_tally2_d_matrix A,
    magma_tally2_d_matrix *D,
    magma_tally2_d_matrix *R,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dmscale( 
    magma_tally2_d_matrix *A, 
    magma_tally2_scale_t scaling,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_dmdiff( 
    magma_tally2_d_matrix A, 
    magma_tally2_d_matrix B, 
 real_Double_t *res,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dmdiagadd( 
    magma_tally2_d_matrix *A, 
    double add,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_dmsort( 
    magma_tally2_d_matrix *A,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dindexsort(
    magma_tally2_index_t *x, 
    magma_tally2_int_t first,
    magma_tally2_int_t last,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_ddomainoverlap(
    magma_tally2_index_t num_rows,
    magma_tally2_index_t *num_indices,
    magma_tally2_index_t *rowptr,
    magma_tally2_index_t *colidx,
    magma_tally2_index_t *x,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dsymbilu( 
    magma_tally2_d_matrix *A, 
    magma_tally2_int_t levels,
    magma_tally2_d_matrix *L,
    magma_tally2_d_matrix *U,
    magma_tally2_queue_t queue );


magma_tally2_int_t 
magma_tally2_dwrite_csr_mtx( 
    magma_tally2_d_matrix A,
    magma_tally2_order_t MajorType,
 const char *filename,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_dwrite_csrtomtx( 
    magma_tally2_d_matrix A,
    const char *filename,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_dprint_csr( 
    magma_tally2_int_t n_row, 
    magma_tally2_int_t n_col, 
    magma_tally2_int_t nnz, 
    double **val, 
    magma_tally2_index_t **row, 
    magma_tally2_index_t **col,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_dprint_csr_mtx( 
    magma_tally2_int_t n_row, 
    magma_tally2_int_t n_col, 
    magma_tally2_int_t nnz, 
    double **val, 
    magma_tally2_index_t **row, 
    magma_tally2_index_t **col, 
    magma_tally2_order_t MajorType,
    magma_tally2_queue_t queue );


magma_tally2_int_t 
magma_tally2_dmtranspose(
    magma_tally2_d_matrix A, 
    magma_tally2_d_matrix *B,
    magma_tally2_queue_t queue );


magma_tally2_int_t 
magma_tally2_dmtransfer(
    magma_tally2_d_matrix A, 
    magma_tally2_d_matrix *B, 
    magma_tally2_location_t src, 
    magma_tally2_location_t dst,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_dmconvert(
    magma_tally2_d_matrix A, 
    magma_tally2_d_matrix *B, 
    magma_tally2_storage_t old_format, 
    magma_tally2_storage_t new_format,
    magma_tally2_queue_t queue );


magma_tally2_int_t
magma_tally2_dvinit(
    magma_tally2_d_matrix *x, 
    magma_tally2_location_t memory_location,
    magma_tally2_int_t num_rows, 
    magma_tally2_int_t num_cols,
    double values,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dprint_vector(
    magma_tally2_d_matrix x, 
    magma_tally2_int_t offset, 
    magma_tally2_int_t displaylength,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dvread(
    magma_tally2_d_matrix *x, 
    magma_tally2_int_t length,
    char * filename,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dvspread(
    magma_tally2_d_matrix *x, 
    const char * filename,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dprint_matrix(
    magma_tally2_d_matrix A,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_ddiameter(
    magma_tally2_d_matrix *A,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_drowentries(
    magma_tally2_d_matrix *A,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dmfree(
    magma_tally2_d_matrix *A,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dresidual(
    magma_tally2_d_matrix A, 
    magma_tally2_d_matrix b, 
    magma_tally2_d_matrix x, 
    double *res,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dresidualvec(
    magma_tally2_d_matrix A,
    magma_tally2_d_matrix b,
    magma_tally2_d_matrix x,
    magma_tally2_d_matrix *r,
    double *res,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dmgenerator(
    magma_tally2_int_t n,
    magma_tally2_int_t offdiags,
    magma_tally2_index_t *diag_offset,
    double *diag_vals,
    magma_tally2_d_matrix *A,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dm_27stencil(
    magma_tally2_int_t n,
    magma_tally2_d_matrix *A,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dm_5stencil(
    magma_tally2_int_t n,
    magma_tally2_d_matrix *A,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dsolverinfo(
    magma_tally2_d_solver_par *solver_par, 
    magma_tally2_d_preconditioner *precond_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dsolverinfo_init(
    magma_tally2_d_solver_par *solver_par, 
    magma_tally2_d_preconditioner *precond,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_deigensolverinfo_init(
    magma_tally2_d_solver_par *solver_par,
    magma_tally2_queue_t queue );


magma_tally2_int_t
magma_tally2_dsolverinfo_free(
    magma_tally2_d_solver_par *solver_par, 
    magma_tally2_d_preconditioner *precond,
    magma_tally2_queue_t queue );



/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally2_SPARSE iterative incomplete factorizations
*/


magma_tally2_int_t
magma_tally2_diterilusetup( 
    magma_tally2_d_matrix A, 
    magma_tally2_d_matrix b,                                 
    magma_tally2_d_preconditioner *precond,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_ditericsetup( 
    magma_tally2_d_matrix A, 
    magma_tally2_d_matrix b, 
    magma_tally2_d_preconditioner *precond,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_ditericupdate( 
    magma_tally2_d_matrix A, 
    magma_tally2_d_preconditioner *precond, 
    magma_tally2_int_t updates,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dapplyiteric_l( 
    magma_tally2_d_matrix b, 
    magma_tally2_d_matrix *x, 
    magma_tally2_d_preconditioner *precond,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dapplyiteric_r( 
    magma_tally2_d_matrix b, 
    magma_tally2_d_matrix *x, 
    magma_tally2_d_preconditioner *precond,
    magma_tally2_queue_t queue );


magma_tally2_int_t
magma_tally2_diterilu_csr( 
    magma_tally2_d_matrix A,
    magma_tally2_d_matrix L,
    magma_tally2_d_matrix U,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_diteric_csr( 
    magma_tally2_d_matrix A,
    magma_tally2_d_matrix A_CSR,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_dfrobenius( 
    magma_tally2_d_matrix A, 
    magma_tally2_d_matrix B, 
    real_Double_t *res,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_dnonlinres(   
    magma_tally2_d_matrix A, 
    magma_tally2_d_matrix L,
    magma_tally2_d_matrix U, 
    magma_tally2_d_matrix *LU, 
    real_Double_t *res,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_dilures(   
    magma_tally2_d_matrix A, 
    magma_tally2_d_matrix L,
    magma_tally2_d_matrix U, 
    magma_tally2_d_matrix *LU, 
    real_Double_t *res,
    real_Double_t *nonlinres,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_dicres(       
    magma_tally2_d_matrix A, 
    magma_tally2_d_matrix C,
    magma_tally2_d_matrix CT, 
    magma_tally2_d_matrix *LU, 
    real_Double_t *res,
    real_Double_t *nonlinres,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_dinitguess( 
    magma_tally2_d_matrix A, 
    magma_tally2_d_matrix *L, 
    magma_tally2_d_matrix *U,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_dinitrecursiveLU( 
    magma_tally2_d_matrix A, 
    magma_tally2_d_matrix *B,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_dmLdiagadd( 
    magma_tally2_d_matrix *L,
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
magma_tally2_dcg(
    magma_tally2_d_matrix A, magma_tally2_d_matrix b, 
    magma_tally2_d_matrix *x, magma_tally2_d_solver_par *solver_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_dcg_res(
    magma_tally2_d_matrix A, magma_tally2_d_matrix b, 
    magma_tally2_d_matrix *x, magma_tally2_d_solver_par *solver_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_dcg_merge(
    magma_tally2_d_matrix A, magma_tally2_d_matrix b, 
    magma_tally2_d_matrix *x, magma_tally2_d_solver_par *solver_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_dgmres(
    magma_tally2_d_matrix A, magma_tally2_d_matrix b, 
    magma_tally2_d_matrix *x, magma_tally2_d_solver_par *solver_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dbicgstab(
    magma_tally2_d_matrix A, magma_tally2_d_matrix b, magma_tally2_d_matrix *x, 
    magma_tally2_d_solver_par *solver_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dbicgstab_merge(
    magma_tally2_d_matrix A, magma_tally2_d_matrix b, 
    magma_tally2_d_matrix *x, magma_tally2_d_solver_par *solver_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dbicgstab_merge2(
    magma_tally2_d_matrix A, magma_tally2_d_matrix b, 
    magma_tally2_d_matrix *x, magma_tally2_d_solver_par *solver_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dpcg(
    magma_tally2_d_matrix A, magma_tally2_d_matrix b, 
    magma_tally2_d_matrix *x, magma_tally2_d_solver_par *solver_par, 
    magma_tally2_d_preconditioner *precond_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dbpcg(
    magma_tally2_d_matrix A, magma_tally2_d_matrix b, 
    magma_tally2_d_matrix *x, magma_tally2_d_solver_par *solver_par, 
    magma_tally2_d_preconditioner *precond_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dpbicgstab(
    magma_tally2_d_matrix A, magma_tally2_d_matrix b, 
    magma_tally2_d_matrix *x, magma_tally2_d_solver_par *solver_par, 
    magma_tally2_d_preconditioner *precond_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dpgmres(
    magma_tally2_d_matrix A, magma_tally2_d_matrix b, 
    magma_tally2_d_matrix *x, magma_tally2_d_solver_par *solver_par, 
    magma_tally2_d_preconditioner *precond_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dfgmres(
    magma_tally2_d_matrix A, magma_tally2_d_matrix b, 
    magma_tally2_d_matrix *x, magma_tally2_d_solver_par *solver_par, 
    magma_tally2_d_preconditioner *precond_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_djacobi(
    magma_tally2_d_matrix A, magma_tally2_d_matrix b, 
    magma_tally2_d_matrix *x, magma_tally2_d_solver_par *solver_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_djacobidomainoverlap(
    magma_tally2_d_matrix A, 
    magma_tally2_d_matrix b, 
    magma_tally2_d_matrix *x,  
    magma_tally2_d_solver_par *solver_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dbaiter(
    magma_tally2_d_matrix A, magma_tally2_d_matrix b, 
    magma_tally2_d_matrix *x, magma_tally2_d_solver_par *solver_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_diterref(
    magma_tally2_d_matrix A, magma_tally2_d_matrix b, 
    magma_tally2_d_matrix *x, magma_tally2_d_solver_par *solver_par, 
    magma_tally2_d_preconditioner *precond_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dilu(
    magma_tally2_d_matrix A, magma_tally2_d_matrix b, 
    magma_tally2_d_matrix *x, magma_tally2_d_solver_par *solver_par, 
    magma_tally2_int_t *ipiv,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dbcsrlu(
    magma_tally2_d_matrix A, magma_tally2_d_matrix b, 
    magma_tally2_d_matrix *x, magma_tally2_d_solver_par *solver_par,
    magma_tally2_queue_t queue );


magma_tally2_int_t
magma_tally2_dbcsrlutrf(
    magma_tally2_d_matrix A, 
    magma_tally2_d_matrix *M,
    magma_tally2_int_t *ipiv, 
    magma_tally2_int_t version,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dbcsrlusv(
    magma_tally2_d_matrix A, magma_tally2_d_matrix b, 
    magma_tally2_d_matrix *x, magma_tally2_d_solver_par *solver_par, 
    magma_tally2_int_t *ipiv,
    magma_tally2_queue_t queue );



magma_tally2_int_t
magma_tally2_dilucg(
    magma_tally2_d_matrix A, magma_tally2_d_matrix b, 
    magma_tally2_d_matrix *x, magma_tally2_d_solver_par *solver_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dilugmres(
    magma_tally2_d_matrix A, magma_tally2_d_matrix b, 
    magma_tally2_d_matrix *x, magma_tally2_d_solver_par *solver_par,
    magma_tally2_queue_t queue ); 


magma_tally2_int_t
magma_tally2_dlobpcg_shift(
    magma_tally2_int_t num_rows,
    magma_tally2_int_t num_vecs, 
    magma_tally2_int_t shift,
    magma_tally2Double_ptr x,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dlobpcg_res(
    magma_tally2_int_t num_rows,
    magma_tally2_int_t num_vecs, 
    double *evalues, 
    magma_tally2Double_ptr X,
    magma_tally2Double_ptr R, 
    double *res,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dlobpcg_maxpy(
    magma_tally2_int_t num_rows,
    magma_tally2_int_t num_vecs, 
    magma_tally2Double_ptr X,
    magma_tally2Double_ptr Y,
    magma_tally2_queue_t queue );


/*/////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally2_SPARSE eigensolvers (Data on GPU)
*/
magma_tally2_int_t
magma_tally2_dlobpcg(
    magma_tally2_d_matrix A, 
    magma_tally2_d_solver_par *solver_par,
    magma_tally2_d_preconditioner *precond_par, 
    magma_tally2_queue_t queue );




/*/////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally2_SPARSE preconditioners (Data on GPU)
*/

magma_tally2_int_t
magma_tally2_djacobisetup(
    magma_tally2_d_matrix A, 
    magma_tally2_d_matrix b, 
    magma_tally2_d_matrix *M, 
    magma_tally2_d_matrix *c,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_djacobisetup_matrix(
    magma_tally2_d_matrix A, 
    magma_tally2_d_matrix *M, 
    magma_tally2_d_matrix *d,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_djacobisetup_vector(
    magma_tally2_d_matrix b, 
    magma_tally2_d_matrix d, 
    magma_tally2_d_matrix *c,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_djacobiiter(
    magma_tally2_d_matrix M, 
    magma_tally2_d_matrix c, 
    magma_tally2_d_matrix *x, 
    magma_tally2_d_solver_par *solver_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_djacobiiter_precond( 
    magma_tally2_d_matrix M, 
    magma_tally2_d_matrix *x, 
    magma_tally2_d_solver_par *solver_par, 
    magma_tally2_d_preconditioner *precond,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_djacobiiter_sys(
    magma_tally2_d_matrix A, 
    magma_tally2_d_matrix b, 
    magma_tally2_d_matrix d, 
    magma_tally2_d_matrix t, 
    magma_tally2_d_matrix *x,  
    magma_tally2_d_solver_par *solver_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dpastixsetup(
    magma_tally2_d_matrix A, magma_tally2_d_matrix b,
    magma_tally2_d_preconditioner *precond,
    magma_tally2_queue_t queue );


magma_tally2_int_t
magma_tally2_dapplypastix(
    magma_tally2_d_matrix b, magma_tally2_d_matrix *x, 
    magma_tally2_d_preconditioner *precond,
    magma_tally2_queue_t queue );


// custom preconditioner
magma_tally2_int_t
magma_tally2_dapplycustomprecond_l(
    magma_tally2_d_matrix b, 
    magma_tally2_d_matrix *x, 
    magma_tally2_d_preconditioner *precond,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dapplycustomprecond_r(
    magma_tally2_d_matrix b, 
    magma_tally2_d_matrix *x, 
    magma_tally2_d_preconditioner *precond,
    magma_tally2_queue_t queue );


// CUSPARSE preconditioner

magma_tally2_int_t
magma_tally2_dcuilusetup(
    magma_tally2_d_matrix A, magma_tally2_d_preconditioner *precond,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dapplycuilu_l(
    magma_tally2_d_matrix b, magma_tally2_d_matrix *x, 
    magma_tally2_d_preconditioner *precond,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dapplycuilu_r(
    magma_tally2_d_matrix b, magma_tally2_d_matrix *x, 
    magma_tally2_d_preconditioner *precond,
    magma_tally2_queue_t queue );


magma_tally2_int_t
magma_tally2_dcuiccsetup(
    magma_tally2_d_matrix A, magma_tally2_d_preconditioner *precond,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dapplycuicc_l(
    magma_tally2_d_matrix b, magma_tally2_d_matrix *x, 
    magma_tally2_d_preconditioner *precond,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dapplycuicc_r(
    magma_tally2_d_matrix b, magma_tally2_d_matrix *x, 
    magma_tally2_d_preconditioner *precond,
    magma_tally2_queue_t queue );


magma_tally2_int_t
magma_tally2_dcumilusetup(
    magma_tally2_d_matrix A, 
    magma_tally2_d_preconditioner *precond,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dcumilugeneratesolverinfo(
    magma_tally2_d_preconditioner *precond,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dapplycumilu_l(
    magma_tally2_d_matrix b, 
    magma_tally2_d_matrix *x, 
    magma_tally2_d_preconditioner *precond,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dapplycumilu_r(
    magma_tally2_d_matrix b, 
    magma_tally2_d_matrix *x, 
    magma_tally2_d_preconditioner *precond,
    magma_tally2_queue_t queue );


magma_tally2_int_t
magma_tally2_dcumiccsetup(
    magma_tally2_d_matrix A, 
    magma_tally2_d_preconditioner *precond,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dcumicgeneratesolverinfo(
    magma_tally2_d_preconditioner *precond,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dapplycumicc_l(
    magma_tally2_d_matrix b, 
    magma_tally2_d_matrix *x, 
    magma_tally2_d_preconditioner *precond,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dapplycumicc_r(
    magma_tally2_d_matrix b, 
    magma_tally2_d_matrix *x, 
    magma_tally2_d_preconditioner *precond,
    magma_tally2_queue_t queue );


// block-asynchronous iteration

magma_tally2_int_t
magma_tally2_dbajac_csr(
    magma_tally2_int_t localiters,
    magma_tally2_d_matrix D,
    magma_tally2_d_matrix R,
    magma_tally2_d_matrix b,
    magma_tally2_d_matrix *x,
    magma_tally2_queue_t queue );

/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally2_SPARSE utility function definitions
*/

magma_tally2_int_t
magma_tally2_d_spmv(
    double alpha, 
    magma_tally2_d_matrix A, 
    magma_tally2_d_matrix x, 
    double beta, 
    magma_tally2_d_matrix y,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dcustomspmv(
    double alpha, 
    magma_tally2_d_matrix x, 
    double beta, 
    magma_tally2_d_matrix y,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_d_spmv_shift(
    double alpha, 
    magma_tally2_d_matrix A, 
    double lambda,
    magma_tally2_d_matrix x, 
    double beta, 
    magma_tally2_int_t offset, 
    magma_tally2_int_t blocksize,
    magma_tally2Index_ptr dadd_vecs, 
    magma_tally2_d_matrix y,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dcuspmm(
    magma_tally2_d_matrix A, 
    magma_tally2_d_matrix B, 
    magma_tally2_d_matrix *AB,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_d_spmm(
    double alpha, 
    magma_tally2_d_matrix A,
    magma_tally2_d_matrix B,
    magma_tally2_d_matrix *C,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dsymbilu( 
    magma_tally2_d_matrix *A, 
    magma_tally2_int_t levels, 
    magma_tally2_d_matrix *L, 
    magma_tally2_d_matrix *U,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dcuspaxpy(
    magma_tally2Double_ptr alpha, magma_tally2_d_matrix A, 
    magma_tally2Double_ptr beta, magma_tally2_d_matrix B, 
    magma_tally2_d_matrix *AB,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_d_precond(
    magma_tally2_d_matrix A, 
    magma_tally2_d_matrix b, magma_tally2_d_matrix *x,
    magma_tally2_d_preconditioner *precond,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_d_solver(
    magma_tally2_d_matrix A, magma_tally2_d_matrix b, 
    magma_tally2_d_matrix *x, magma_tally2_dopts *zopts,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_d_precondsetup(
    magma_tally2_d_matrix A, magma_tally2_d_matrix b, 
    magma_tally2_d_preconditioner *precond,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_d_applyprecond(
    magma_tally2_d_matrix A, magma_tally2_d_matrix b, 
    magma_tally2_d_matrix *x, magma_tally2_d_preconditioner *precond,
    magma_tally2_queue_t queue );


magma_tally2_int_t
magma_tally2_d_applyprecond_left(
    magma_tally2_d_matrix A, magma_tally2_d_matrix b, 
    magma_tally2_d_matrix *x, magma_tally2_d_preconditioner *precond,
    magma_tally2_queue_t queue );


magma_tally2_int_t
magma_tally2_d_applyprecond_right(
    magma_tally2_d_matrix A, magma_tally2_d_matrix b, 
    magma_tally2_d_matrix *x, magma_tally2_d_preconditioner *precond,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_d_initP2P(
    magma_tally2_int_t *bandwidth_benchmark,
    magma_tally2_int_t *num_gpus,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dcompact(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda,
    double *dnorms, double tol, 
    magma_tally2_int_t *activeMask, magma_tally2_int_t *cBlockSize,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dcompactActive(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2Double_ptr dA, magma_tally2_int_t ldda, 
    magma_tally2_int_t *active,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dmlumerge(    
    magma_tally2_d_matrix L, 
    magma_tally2_d_matrix U,
    magma_tally2_d_matrix *A, 
    magma_tally2_queue_t queue );


/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally2_SPARSE BLAS function definitions
*/
magma_tally2_int_t 
magma_tally2_dgecsrmv(
    magma_tally2_trans_t transA,
    magma_tally2_int_t m, magma_tally2_int_t n,
    double alpha,
    magma_tally2Double_ptr dval,
    magma_tally2Index_ptr drowptr,
    magma_tally2Index_ptr dcolind,
    magma_tally2Double_ptr dx,
    double beta,
    magma_tally2Double_ptr dy,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_dgecsrmv_shift(
    magma_tally2_trans_t transA,
    magma_tally2_int_t m, magma_tally2_int_t n,
    double alpha,
    double lambda,
    magma_tally2Double_ptr dval,
    magma_tally2Index_ptr drowptr,
    magma_tally2Index_ptr dcolind,
    magma_tally2Double_ptr dx,
    double beta,
    int offset,
    int blocksize,
    magma_tally2Index_ptr dadd_rows,
    magma_tally2Double_ptr dy,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_dmgecsrmv(
    magma_tally2_trans_t transA,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2_int_t num_vecs,
    double alpha,
    magma_tally2Double_ptr dval,
    magma_tally2Index_ptr drowptr,
    magma_tally2Index_ptr dcolind,
    magma_tally2Double_ptr dx,
    double beta,
    magma_tally2Double_ptr dy,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_dgeellmv(
    magma_tally2_trans_t transA,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2_int_t nnz_per_row,
    double alpha,
    magma_tally2Double_ptr dval,
    magma_tally2Index_ptr dcolind,
    magma_tally2Double_ptr dx,
    double beta,
    magma_tally2Double_ptr dy,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_dgeellmv_shift(
    magma_tally2_trans_t transA,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2_int_t nnz_per_row,
    double alpha,
    double lambda,
    magma_tally2Double_ptr dval,
    magma_tally2Index_ptr dcolind,
    magma_tally2Double_ptr dx,
    double beta,
    int offset,
    int blocksize,
    magma_tally2Index_ptr dadd_rows,
    magma_tally2Double_ptr dy,
    magma_tally2_queue_t queue );


magma_tally2_int_t 
magma_tally2_dmgeellmv(
    magma_tally2_trans_t transA,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2_int_t num_vecs,
    magma_tally2_int_t nnz_per_row,
    double alpha,
    magma_tally2Double_ptr dval,
    magma_tally2Index_ptr dcolind,
    magma_tally2Double_ptr dx,
    double beta,
    magma_tally2Double_ptr dy,
    magma_tally2_queue_t queue );


magma_tally2_int_t 
magma_tally2_dgeelltmv(
    magma_tally2_trans_t transA,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2_int_t nnz_per_row,
    double alpha,
    magma_tally2Double_ptr dval,
    magma_tally2Index_ptr dcolind,
    magma_tally2Double_ptr dx,
    double beta,
    magma_tally2Double_ptr dy,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_dgeelltmv_shift(
    magma_tally2_trans_t transA,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2_int_t nnz_per_row,
    double alpha,
    double lambda,
    magma_tally2Double_ptr dval,
    magma_tally2Index_ptr dcolind,
    magma_tally2Double_ptr dx,
    double beta,
    int offset,
    int blocksize,
    magma_tally2Index_ptr dadd_rows,
    magma_tally2Double_ptr dy,
    magma_tally2_queue_t queue );


magma_tally2_int_t 
magma_tally2_dmgeelltmv(
    magma_tally2_trans_t transA,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2_int_t num_vecs,
    magma_tally2_int_t nnz_per_row,
    double alpha,
    magma_tally2Double_ptr dval,
    magma_tally2Index_ptr dcolind,
    magma_tally2Double_ptr dx,
    double beta,
    magma_tally2Double_ptr dy,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_dgeellrtmv(
    magma_tally2_trans_t transA,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2_int_t nnz_per_row,
    double alpha,
    magma_tally2Double_ptr dval,
    magma_tally2Index_ptr dcolind,
    magma_tally2Index_ptr drowlength,
    magma_tally2Double_ptr dx,
    double beta,
    magma_tally2Double_ptr dy,
    magma_tally2_int_t num_threads,
    magma_tally2_int_t threads_per_row,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_dgesellcmv(
    magma_tally2_trans_t transA,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2_int_t blocksize,
    magma_tally2_int_t slices,
    magma_tally2_int_t alignment,
    double alpha,
    magma_tally2Double_ptr dval,
    magma_tally2Index_ptr dcolind,
    magma_tally2Index_ptr drowptr,
    magma_tally2Double_ptr dx,
    double beta,
    magma_tally2Double_ptr dy,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dgesellpmv(
    magma_tally2_trans_t transA,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2_int_t blocksize,
    magma_tally2_int_t slices,
    magma_tally2_int_t alignment,
    double alpha,
    magma_tally2Double_ptr dval,
    magma_tally2Index_ptr dcolind,
    magma_tally2Index_ptr drowptr,
    magma_tally2Double_ptr dx,
    double beta,
    magma_tally2Double_ptr dy,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dmgesellpmv(
    magma_tally2_trans_t transA,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2_int_t num_vecs,
    magma_tally2_int_t blocksize,
    magma_tally2_int_t slices,
    magma_tally2_int_t alignment,
    double alpha,
    magma_tally2Double_ptr dval,
    magma_tally2Index_ptr dcolind,
    magma_tally2Index_ptr drowptr,
    magma_tally2Double_ptr dx,
    double beta,
    magma_tally2Double_ptr dy,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dmgesellpmv_blocked(
    magma_tally2_trans_t transA,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2_int_t num_vecs,
    magma_tally2_int_t blocksize,
    magma_tally2_int_t slices,
    magma_tally2_int_t alignment,
    double alpha,
    magma_tally2Double_ptr dval,
    magma_tally2Index_ptr dcolind,
    magma_tally2Index_ptr drowptr,
    magma_tally2Double_ptr dx,
    double beta,
    magma_tally2Double_ptr dy,
    magma_tally2_queue_t queue );


magma_tally2_int_t
magma_tally2_dmergedgs(
    magma_tally2_int_t n, 
    magma_tally2_int_t ldh,
    magma_tally2_int_t k, 
    magma_tally2Double_ptr dv, 
    magma_tally2Double_ptr dr,
    magma_tally2Double_ptr dskp,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dcopyscale(    
    int n, 
    int k,
    magma_tally2Double_ptr dr, 
    magma_tally2Double_ptr dv,
    magma_tally2Double_ptr dskp,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dnrm2scale(    
    int m, 
    magma_tally2Double_ptr dr,    
    int lddr, 
    double *drnorm,
    magma_tally2_queue_t queue );


magma_tally2_int_t
magma_tally2_djacobisetup_vector_gpu(
    int num_rows, 
    magma_tally2_d_matrix b, 
    magma_tally2_d_matrix d, 
    magma_tally2_d_matrix c,
    magma_tally2_d_matrix *x,
    magma_tally2_queue_t queue );


magma_tally2_int_t
magma_tally2_djacobi_diagscal(    
    int num_rows, 
    magma_tally2_d_matrix d, 
    magma_tally2_d_matrix b, 
    magma_tally2_d_matrix *c,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_djacobiupdate(
    magma_tally2_d_matrix t, 
    magma_tally2_d_matrix b, 
    magma_tally2_d_matrix d, 
    magma_tally2_d_matrix *x,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_djacobispmvupdate(
    magma_tally2_int_t maxiter,
    magma_tally2_d_matrix A, 
    magma_tally2_d_matrix t, 
    magma_tally2_d_matrix b, 
    magma_tally2_d_matrix d, 
    magma_tally2_d_matrix *x,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_djacobispmvupdate_bw(
    magma_tally2_int_t maxiter,
    magma_tally2_d_matrix A, 
    magma_tally2_d_matrix t, 
    magma_tally2_d_matrix b, 
    magma_tally2_d_matrix d, 
    magma_tally2_d_matrix *x,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_djacobispmvupdateselect(
    magma_tally2_int_t maxiter,
    magma_tally2_int_t num_updates,
    magma_tally2_index_t *indices,
    magma_tally2_d_matrix A,
    magma_tally2_d_matrix t, 
    magma_tally2_d_matrix b, 
    magma_tally2_d_matrix d, 
    magma_tally2_d_matrix tmp, 
    magma_tally2_d_matrix *x,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_djacobisetup_diagscal(
    magma_tally2_d_matrix A, magma_tally2_d_matrix *d,
    magma_tally2_queue_t queue );


magma_tally2_int_t
magma_tally2_dbicgmerge1(    
    int n, 
    magma_tally2Double_ptr dskp,
    magma_tally2Double_ptr dv, 
    magma_tally2Double_ptr dr, 
    magma_tally2Double_ptr dp,
    magma_tally2_queue_t queue );


magma_tally2_int_t
magma_tally2_dbicgmerge2(
    int n, 
    magma_tally2Double_ptr dskp, 
    magma_tally2Double_ptr dr,
    magma_tally2Double_ptr dv, 
    magma_tally2Double_ptr ds,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dbicgmerge3(
    int n, 
    magma_tally2Double_ptr dskp, 
    magma_tally2Double_ptr dp,
    magma_tally2Double_ptr ds,
    magma_tally2Double_ptr dt,
    magma_tally2Double_ptr dx, 
    magma_tally2Double_ptr dr,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dbicgmerge4(
    int type, 
    magma_tally2Double_ptr dskp,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dcgmerge_spmv1( 
    magma_tally2_d_matrix A,
    magma_tally2Double_ptr d1,
    magma_tally2Double_ptr d2,
    magma_tally2Double_ptr dd,
    magma_tally2Double_ptr dz,
    magma_tally2Double_ptr dskp,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dcgmerge_xrbeta( 
    int n,
    magma_tally2Double_ptr d1,
    magma_tally2Double_ptr d2,
    magma_tally2Double_ptr dx,
    magma_tally2Double_ptr dr,
    magma_tally2Double_ptr dd,
    magma_tally2Double_ptr dz, 
    magma_tally2Double_ptr dskp,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dmdotc(
    magma_tally2_int_t n, 
    magma_tally2_int_t k, 
    magma_tally2Double_ptr dv, 
    magma_tally2Double_ptr dr,
    magma_tally2Double_ptr dd1,
    magma_tally2Double_ptr dd2,
    magma_tally2Double_ptr dskp,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dgemvmdot(
    int n, 
    int k, 
    magma_tally2Double_ptr dv, 
    magma_tally2Double_ptr dr,
    magma_tally2Double_ptr dd1,
    magma_tally2Double_ptr dd2,
    magma_tally2Double_ptr dskp,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dbicgmerge_spmv1( 
    magma_tally2_d_matrix A,
    magma_tally2Double_ptr dd1,
    magma_tally2Double_ptr dd2,
    magma_tally2Double_ptr dp,
    magma_tally2Double_ptr dr,
    magma_tally2Double_ptr dv,
    magma_tally2Double_ptr dskp,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dbicgmerge_spmv2( 
    magma_tally2_d_matrix A,
    magma_tally2Double_ptr dd1,
    magma_tally2Double_ptr dd2,
    magma_tally2Double_ptr ds,
    magma_tally2Double_ptr dt,
    magma_tally2Double_ptr dskp,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dbicgmerge_xrbeta( 
    int n,
    magma_tally2Double_ptr dd1,
    magma_tally2Double_ptr dd2,
    magma_tally2Double_ptr drr,
    magma_tally2Double_ptr dr,
    magma_tally2Double_ptr dp,
    magma_tally2Double_ptr ds,
    magma_tally2Double_ptr dt,
    magma_tally2Double_ptr dx, 
    magma_tally2Double_ptr dskp,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dbcsrswp(
    magma_tally2_int_t n,
    magma_tally2_int_t size_b, 
    magma_tally2_int_t *ipiv,
    magma_tally2Double_ptr dx,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dbcsrtrsv(
    magma_tally2_uplo_t uplo,
    magma_tally2_int_t r_blocks,
    magma_tally2_int_t c_blocks,
    magma_tally2_int_t size_b, 
    magma_tally2Double_ptr dA,
    magma_tally2_index_t *blockinfo, 
    magma_tally2Double_ptr dx,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dbcsrvalcpy(
    magma_tally2_int_t size_b, 
    magma_tally2_int_t num_blocks, 
    magma_tally2_int_t num_zero_blocks, 
    magma_tally2Double_ptr *dAval, 
    magma_tally2Double_ptr *dBval,
    magma_tally2Double_ptr *dBval2,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dbcsrluegemm(
    magma_tally2_int_t size_b, 
    magma_tally2_int_t num_block_rows,
    magma_tally2_int_t kblocks,
    magma_tally2Double_ptr *dA, 
    magma_tally2Double_ptr *dB, 
    magma_tally2Double_ptr *dC,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dbcsrlupivloc(
    magma_tally2_int_t size_b, 
    magma_tally2_int_t kblocks,
    magma_tally2Double_ptr *dA, 
    magma_tally2_int_t *ipiv,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dbcsrblockinfo5(
    magma_tally2_int_t lustep,
    magma_tally2_int_t num_blocks, 
    magma_tally2_int_t c_blocks, 
    magma_tally2_int_t size_b,
    magma_tally2_index_t *blockinfo,
    magma_tally2Double_ptr dval,
    magma_tally2Double_ptr *AII,
    magma_tally2_queue_t queue );

 
#ifdef __cplusplus
}
#endif

#undef PRECISION_d
#endif /* MAGMA_tally2SPARSE_D_H */
