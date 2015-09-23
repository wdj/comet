/*
 -- MAGMA_minproduct (version 1.1) --
 Univ. of Tennessee, Knoxville
 Univ. of California, Berkeley
 Univ. of Colorado, Denver
 @date

 @generated from magma_minproductsparse_z.h normal z -> d, Mon May  4 17:11:25 2015
 @author Hartwig Anzt
*/

#ifndef MAGMA_minproductSPARSE_D_H
#define MAGMA_minproductSPARSE_D_H

#include "magma_minproduct_types.h"
#include "magma_minproductsparse_types.h"

#define PRECISION_d


#ifdef __cplusplus
extern "C" {
#endif


/* ////////////////////////////////////////////////////////////////////////////
 -- For backwards compatability, map old (1.6.1) to new (1.6.2) function names
*/

#define magma_minproduct_d_mtranspose  magma_minproduct_dmtranspose
#define magma_minproduct_d_mtransfer   magma_minproduct_dmtransfer
#define magma_minproduct_d_vtransfer   magma_minproduct_dmtransfer
#define magma_minproduct_d_mconvert    magma_minproduct_dmconvert
#define magma_minproduct_d_vinit       magma_minproduct_dvinit
#define magma_minproduct_d_vvisu       magma_minproduct_dprint_vector
#define magma_minproduct_d_vread       magma_minproduct_dvread
#define magma_minproduct_d_vspread     magma_minproduct_dvspread
#define magma_minproduct_d_mvisu       magma_minproduct_dprint_matrix
#define magma_minproduct_d_mfree       magma_minproduct_dmfree
#define magma_minproduct_d_vfree       magma_minproduct_dmfree
#define write_d_csr_mtx     magma_minproduct_dwrite_csr_mtx
#define write_d_csrtomtx    magma_minproduct_dwrite_csrtomtx
#define print_d_csr         magma_minproduct_dprint_csr_mtx


/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_minproduct_SPARSE Auxiliary functions
*/


magma_minproduct_int_t
magma_minproduct_dparse_opts( 
    int argc, 
    char** argv, 
    magma_minproduct_dopts *opts, 
    int *matrices, 
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
read_d_csr_from_binary( 
    magma_minproduct_int_t* n_row, 
    magma_minproduct_int_t* n_col, 
    magma_minproduct_int_t* nnz, 
    double **val, 
    magma_minproduct_index_t **row, 
    magma_minproduct_index_t **col,
    const char * filename,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
read_d_csr_from_mtx( 
    magma_minproduct_storage_t *type, 
    magma_minproduct_location_t *location,
    magma_minproduct_int_t* n_row, 
    magma_minproduct_int_t* n_col, 
    magma_minproduct_int_t* nnz, 
    double **val, 
    magma_minproduct_index_t **row, 
    magma_minproduct_index_t **col, 
    const char *filename,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_d_csr_mtx( 
    magma_minproduct_d_matrix *A, 
    const char *filename,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_dcsrset( 
    magma_minproduct_int_t m, 
    magma_minproduct_int_t n, 
    magma_minproduct_index_t *row, 
    magma_minproduct_index_t *col, 
    double *val,
    magma_minproduct_d_matrix *A,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_dcsrget( 
    magma_minproduct_d_matrix A,
    magma_minproduct_int_t *m, 
    magma_minproduct_int_t *n, 
    magma_minproduct_index_t **row, 
    magma_minproduct_index_t **col, 
    double **val,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t 
magma_minproduct_dvset( 
    magma_minproduct_int_t m, 
    magma_minproduct_int_t n, 
    double *val,
    magma_minproduct_d_matrix *v,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_dvget( 
    magma_minproduct_d_matrix v,
    magma_minproduct_int_t *m, 
    magma_minproduct_int_t *n, 
    double **val,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_dvset_dev( 
    magma_minproduct_int_t m, 
    magma_minproduct_int_t n, 
    magma_minproductDouble_ptr val,
    magma_minproduct_d_matrix *v,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_dvget_dev( 
    magma_minproduct_d_matrix v,
    magma_minproduct_int_t *m, 
    magma_minproduct_int_t *n, 
    magma_minproductDouble_ptr *val,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t 
magma_minproduct_d_csr_mtxsymm( 
    magma_minproduct_d_matrix *A, 
    const char *filename,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_d_csr_compressor( 
    double ** val, 
    magma_minproduct_index_t ** row, 
    magma_minproduct_index_t ** col, 
    double ** valn, 
    magma_minproduct_index_t ** rown, 
    magma_minproduct_index_t ** coln, 
    magma_minproduct_int_t *n,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dmcsrcompressor( 
    magma_minproduct_d_matrix *A,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dmcsrcompressor_gpu( 
    magma_minproduct_d_matrix *A,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dvtranspose( 
    magma_minproduct_d_matrix x,
    magma_minproduct_d_matrix *y,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_d_cucsrtranspose( 
    magma_minproduct_d_matrix A, 
    magma_minproduct_d_matrix *B,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
d_transpose_csr( 
    magma_minproduct_int_t n_rows, 
    magma_minproduct_int_t n_cols, 
    magma_minproduct_int_t nnz,
    double *val, 
    magma_minproduct_index_t *row, 
    magma_minproduct_index_t *col, 
    magma_minproduct_int_t *new_n_rows, 
    magma_minproduct_int_t *new_n_cols, 
    magma_minproduct_int_t *new_nnz, 
    double **new_val, 
    magma_minproduct_index_t **new_row, 
    magma_minproduct_index_t **new_col,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dcsrsplit( 
    magma_minproduct_int_t bsize,
    magma_minproduct_d_matrix A,
    magma_minproduct_d_matrix *D,
    magma_minproduct_d_matrix *R,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dmscale( 
    magma_minproduct_d_matrix *A, 
    magma_minproduct_scale_t scaling,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_dmdiff( 
    magma_minproduct_d_matrix A, 
    magma_minproduct_d_matrix B, 
 real_Double_t *res,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dmdiagadd( 
    magma_minproduct_d_matrix *A, 
    double add,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_dmsort( 
    magma_minproduct_d_matrix *A,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dindexsort(
    magma_minproduct_index_t *x, 
    magma_minproduct_int_t first,
    magma_minproduct_int_t last,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_ddomainoverlap(
    magma_minproduct_index_t num_rows,
    magma_minproduct_index_t *num_indices,
    magma_minproduct_index_t *rowptr,
    magma_minproduct_index_t *colidx,
    magma_minproduct_index_t *x,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dsymbilu( 
    magma_minproduct_d_matrix *A, 
    magma_minproduct_int_t levels,
    magma_minproduct_d_matrix *L,
    magma_minproduct_d_matrix *U,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t 
magma_minproduct_dwrite_csr_mtx( 
    magma_minproduct_d_matrix A,
    magma_minproduct_order_t MajorType,
 const char *filename,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_dwrite_csrtomtx( 
    magma_minproduct_d_matrix A,
    const char *filename,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_dprint_csr( 
    magma_minproduct_int_t n_row, 
    magma_minproduct_int_t n_col, 
    magma_minproduct_int_t nnz, 
    double **val, 
    magma_minproduct_index_t **row, 
    magma_minproduct_index_t **col,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_dprint_csr_mtx( 
    magma_minproduct_int_t n_row, 
    magma_minproduct_int_t n_col, 
    magma_minproduct_int_t nnz, 
    double **val, 
    magma_minproduct_index_t **row, 
    magma_minproduct_index_t **col, 
    magma_minproduct_order_t MajorType,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t 
magma_minproduct_dmtranspose(
    magma_minproduct_d_matrix A, 
    magma_minproduct_d_matrix *B,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t 
magma_minproduct_dmtransfer(
    magma_minproduct_d_matrix A, 
    magma_minproduct_d_matrix *B, 
    magma_minproduct_location_t src, 
    magma_minproduct_location_t dst,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_dmconvert(
    magma_minproduct_d_matrix A, 
    magma_minproduct_d_matrix *B, 
    magma_minproduct_storage_t old_format, 
    magma_minproduct_storage_t new_format,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t
magma_minproduct_dvinit(
    magma_minproduct_d_matrix *x, 
    magma_minproduct_location_t memory_location,
    magma_minproduct_int_t num_rows, 
    magma_minproduct_int_t num_cols,
    double values,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dprint_vector(
    magma_minproduct_d_matrix x, 
    magma_minproduct_int_t offset, 
    magma_minproduct_int_t displaylength,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dvread(
    magma_minproduct_d_matrix *x, 
    magma_minproduct_int_t length,
    char * filename,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dvspread(
    magma_minproduct_d_matrix *x, 
    const char * filename,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dprint_matrix(
    magma_minproduct_d_matrix A,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_ddiameter(
    magma_minproduct_d_matrix *A,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_drowentries(
    magma_minproduct_d_matrix *A,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dmfree(
    magma_minproduct_d_matrix *A,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dresidual(
    magma_minproduct_d_matrix A, 
    magma_minproduct_d_matrix b, 
    magma_minproduct_d_matrix x, 
    double *res,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dresidualvec(
    magma_minproduct_d_matrix A,
    magma_minproduct_d_matrix b,
    magma_minproduct_d_matrix x,
    magma_minproduct_d_matrix *r,
    double *res,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dmgenerator(
    magma_minproduct_int_t n,
    magma_minproduct_int_t offdiags,
    magma_minproduct_index_t *diag_offset,
    double *diag_vals,
    magma_minproduct_d_matrix *A,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dm_27stencil(
    magma_minproduct_int_t n,
    magma_minproduct_d_matrix *A,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dm_5stencil(
    magma_minproduct_int_t n,
    magma_minproduct_d_matrix *A,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dsolverinfo(
    magma_minproduct_d_solver_par *solver_par, 
    magma_minproduct_d_preconditioner *precond_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dsolverinfo_init(
    magma_minproduct_d_solver_par *solver_par, 
    magma_minproduct_d_preconditioner *precond,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_deigensolverinfo_init(
    magma_minproduct_d_solver_par *solver_par,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t
magma_minproduct_dsolverinfo_free(
    magma_minproduct_d_solver_par *solver_par, 
    magma_minproduct_d_preconditioner *precond,
    magma_minproduct_queue_t queue );



/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_minproduct_SPARSE iterative incomplete factorizations
*/


magma_minproduct_int_t
magma_minproduct_diterilusetup( 
    magma_minproduct_d_matrix A, 
    magma_minproduct_d_matrix b,                                 
    magma_minproduct_d_preconditioner *precond,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_ditericsetup( 
    magma_minproduct_d_matrix A, 
    magma_minproduct_d_matrix b, 
    magma_minproduct_d_preconditioner *precond,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_ditericupdate( 
    magma_minproduct_d_matrix A, 
    magma_minproduct_d_preconditioner *precond, 
    magma_minproduct_int_t updates,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dapplyiteric_l( 
    magma_minproduct_d_matrix b, 
    magma_minproduct_d_matrix *x, 
    magma_minproduct_d_preconditioner *precond,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dapplyiteric_r( 
    magma_minproduct_d_matrix b, 
    magma_minproduct_d_matrix *x, 
    magma_minproduct_d_preconditioner *precond,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t
magma_minproduct_diterilu_csr( 
    magma_minproduct_d_matrix A,
    magma_minproduct_d_matrix L,
    magma_minproduct_d_matrix U,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_diteric_csr( 
    magma_minproduct_d_matrix A,
    magma_minproduct_d_matrix A_CSR,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_dfrobenius( 
    magma_minproduct_d_matrix A, 
    magma_minproduct_d_matrix B, 
    real_Double_t *res,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_dnonlinres(   
    magma_minproduct_d_matrix A, 
    magma_minproduct_d_matrix L,
    magma_minproduct_d_matrix U, 
    magma_minproduct_d_matrix *LU, 
    real_Double_t *res,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_dilures(   
    magma_minproduct_d_matrix A, 
    magma_minproduct_d_matrix L,
    magma_minproduct_d_matrix U, 
    magma_minproduct_d_matrix *LU, 
    real_Double_t *res,
    real_Double_t *nonlinres,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_dicres(       
    magma_minproduct_d_matrix A, 
    magma_minproduct_d_matrix C,
    magma_minproduct_d_matrix CT, 
    magma_minproduct_d_matrix *LU, 
    real_Double_t *res,
    real_Double_t *nonlinres,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_dinitguess( 
    magma_minproduct_d_matrix A, 
    magma_minproduct_d_matrix *L, 
    magma_minproduct_d_matrix *U,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_dinitrecursiveLU( 
    magma_minproduct_d_matrix A, 
    magma_minproduct_d_matrix *B,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_dmLdiagadd( 
    magma_minproduct_d_matrix *L,
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
magma_minproduct_dcg(
    magma_minproduct_d_matrix A, magma_minproduct_d_matrix b, 
    magma_minproduct_d_matrix *x, magma_minproduct_d_solver_par *solver_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_dcg_res(
    magma_minproduct_d_matrix A, magma_minproduct_d_matrix b, 
    magma_minproduct_d_matrix *x, magma_minproduct_d_solver_par *solver_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_dcg_merge(
    magma_minproduct_d_matrix A, magma_minproduct_d_matrix b, 
    magma_minproduct_d_matrix *x, magma_minproduct_d_solver_par *solver_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_dgmres(
    magma_minproduct_d_matrix A, magma_minproduct_d_matrix b, 
    magma_minproduct_d_matrix *x, magma_minproduct_d_solver_par *solver_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dbicgstab(
    magma_minproduct_d_matrix A, magma_minproduct_d_matrix b, magma_minproduct_d_matrix *x, 
    magma_minproduct_d_solver_par *solver_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dbicgstab_merge(
    magma_minproduct_d_matrix A, magma_minproduct_d_matrix b, 
    magma_minproduct_d_matrix *x, magma_minproduct_d_solver_par *solver_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dbicgstab_merge2(
    magma_minproduct_d_matrix A, magma_minproduct_d_matrix b, 
    magma_minproduct_d_matrix *x, magma_minproduct_d_solver_par *solver_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dpcg(
    magma_minproduct_d_matrix A, magma_minproduct_d_matrix b, 
    magma_minproduct_d_matrix *x, magma_minproduct_d_solver_par *solver_par, 
    magma_minproduct_d_preconditioner *precond_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dbpcg(
    magma_minproduct_d_matrix A, magma_minproduct_d_matrix b, 
    magma_minproduct_d_matrix *x, magma_minproduct_d_solver_par *solver_par, 
    magma_minproduct_d_preconditioner *precond_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dpbicgstab(
    magma_minproduct_d_matrix A, magma_minproduct_d_matrix b, 
    magma_minproduct_d_matrix *x, magma_minproduct_d_solver_par *solver_par, 
    magma_minproduct_d_preconditioner *precond_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dpgmres(
    magma_minproduct_d_matrix A, magma_minproduct_d_matrix b, 
    magma_minproduct_d_matrix *x, magma_minproduct_d_solver_par *solver_par, 
    magma_minproduct_d_preconditioner *precond_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dfgmres(
    magma_minproduct_d_matrix A, magma_minproduct_d_matrix b, 
    magma_minproduct_d_matrix *x, magma_minproduct_d_solver_par *solver_par, 
    magma_minproduct_d_preconditioner *precond_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_djacobi(
    magma_minproduct_d_matrix A, magma_minproduct_d_matrix b, 
    magma_minproduct_d_matrix *x, magma_minproduct_d_solver_par *solver_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_djacobidomainoverlap(
    magma_minproduct_d_matrix A, 
    magma_minproduct_d_matrix b, 
    magma_minproduct_d_matrix *x,  
    magma_minproduct_d_solver_par *solver_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dbaiter(
    magma_minproduct_d_matrix A, magma_minproduct_d_matrix b, 
    magma_minproduct_d_matrix *x, magma_minproduct_d_solver_par *solver_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_diterref(
    magma_minproduct_d_matrix A, magma_minproduct_d_matrix b, 
    magma_minproduct_d_matrix *x, magma_minproduct_d_solver_par *solver_par, 
    magma_minproduct_d_preconditioner *precond_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dilu(
    magma_minproduct_d_matrix A, magma_minproduct_d_matrix b, 
    magma_minproduct_d_matrix *x, magma_minproduct_d_solver_par *solver_par, 
    magma_minproduct_int_t *ipiv,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dbcsrlu(
    magma_minproduct_d_matrix A, magma_minproduct_d_matrix b, 
    magma_minproduct_d_matrix *x, magma_minproduct_d_solver_par *solver_par,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t
magma_minproduct_dbcsrlutrf(
    magma_minproduct_d_matrix A, 
    magma_minproduct_d_matrix *M,
    magma_minproduct_int_t *ipiv, 
    magma_minproduct_int_t version,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dbcsrlusv(
    magma_minproduct_d_matrix A, magma_minproduct_d_matrix b, 
    magma_minproduct_d_matrix *x, magma_minproduct_d_solver_par *solver_par, 
    magma_minproduct_int_t *ipiv,
    magma_minproduct_queue_t queue );



magma_minproduct_int_t
magma_minproduct_dilucg(
    magma_minproduct_d_matrix A, magma_minproduct_d_matrix b, 
    magma_minproduct_d_matrix *x, magma_minproduct_d_solver_par *solver_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dilugmres(
    magma_minproduct_d_matrix A, magma_minproduct_d_matrix b, 
    magma_minproduct_d_matrix *x, magma_minproduct_d_solver_par *solver_par,
    magma_minproduct_queue_t queue ); 


magma_minproduct_int_t
magma_minproduct_dlobpcg_shift(
    magma_minproduct_int_t num_rows,
    magma_minproduct_int_t num_vecs, 
    magma_minproduct_int_t shift,
    magma_minproductDouble_ptr x,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dlobpcg_res(
    magma_minproduct_int_t num_rows,
    magma_minproduct_int_t num_vecs, 
    double *evalues, 
    magma_minproductDouble_ptr X,
    magma_minproductDouble_ptr R, 
    double *res,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dlobpcg_maxpy(
    magma_minproduct_int_t num_rows,
    magma_minproduct_int_t num_vecs, 
    magma_minproductDouble_ptr X,
    magma_minproductDouble_ptr Y,
    magma_minproduct_queue_t queue );


/*/////////////////////////////////////////////////////////////////////////////
 -- MAGMA_minproduct_SPARSE eigensolvers (Data on GPU)
*/
magma_minproduct_int_t
magma_minproduct_dlobpcg(
    magma_minproduct_d_matrix A, 
    magma_minproduct_d_solver_par *solver_par,
    magma_minproduct_d_preconditioner *precond_par, 
    magma_minproduct_queue_t queue );




/*/////////////////////////////////////////////////////////////////////////////
 -- MAGMA_minproduct_SPARSE preconditioners (Data on GPU)
*/

magma_minproduct_int_t
magma_minproduct_djacobisetup(
    magma_minproduct_d_matrix A, 
    magma_minproduct_d_matrix b, 
    magma_minproduct_d_matrix *M, 
    magma_minproduct_d_matrix *c,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_djacobisetup_matrix(
    magma_minproduct_d_matrix A, 
    magma_minproduct_d_matrix *M, 
    magma_minproduct_d_matrix *d,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_djacobisetup_vector(
    magma_minproduct_d_matrix b, 
    magma_minproduct_d_matrix d, 
    magma_minproduct_d_matrix *c,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_djacobiiter(
    magma_minproduct_d_matrix M, 
    magma_minproduct_d_matrix c, 
    magma_minproduct_d_matrix *x, 
    magma_minproduct_d_solver_par *solver_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_djacobiiter_precond( 
    magma_minproduct_d_matrix M, 
    magma_minproduct_d_matrix *x, 
    magma_minproduct_d_solver_par *solver_par, 
    magma_minproduct_d_preconditioner *precond,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_djacobiiter_sys(
    magma_minproduct_d_matrix A, 
    magma_minproduct_d_matrix b, 
    magma_minproduct_d_matrix d, 
    magma_minproduct_d_matrix t, 
    magma_minproduct_d_matrix *x,  
    magma_minproduct_d_solver_par *solver_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dpastixsetup(
    magma_minproduct_d_matrix A, magma_minproduct_d_matrix b,
    magma_minproduct_d_preconditioner *precond,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t
magma_minproduct_dapplypastix(
    magma_minproduct_d_matrix b, magma_minproduct_d_matrix *x, 
    magma_minproduct_d_preconditioner *precond,
    magma_minproduct_queue_t queue );


// custom preconditioner
magma_minproduct_int_t
magma_minproduct_dapplycustomprecond_l(
    magma_minproduct_d_matrix b, 
    magma_minproduct_d_matrix *x, 
    magma_minproduct_d_preconditioner *precond,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dapplycustomprecond_r(
    magma_minproduct_d_matrix b, 
    magma_minproduct_d_matrix *x, 
    magma_minproduct_d_preconditioner *precond,
    magma_minproduct_queue_t queue );


// CUSPARSE preconditioner

magma_minproduct_int_t
magma_minproduct_dcuilusetup(
    magma_minproduct_d_matrix A, magma_minproduct_d_preconditioner *precond,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dapplycuilu_l(
    magma_minproduct_d_matrix b, magma_minproduct_d_matrix *x, 
    magma_minproduct_d_preconditioner *precond,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dapplycuilu_r(
    magma_minproduct_d_matrix b, magma_minproduct_d_matrix *x, 
    magma_minproduct_d_preconditioner *precond,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t
magma_minproduct_dcuiccsetup(
    magma_minproduct_d_matrix A, magma_minproduct_d_preconditioner *precond,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dapplycuicc_l(
    magma_minproduct_d_matrix b, magma_minproduct_d_matrix *x, 
    magma_minproduct_d_preconditioner *precond,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dapplycuicc_r(
    magma_minproduct_d_matrix b, magma_minproduct_d_matrix *x, 
    magma_minproduct_d_preconditioner *precond,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t
magma_minproduct_dcumilusetup(
    magma_minproduct_d_matrix A, 
    magma_minproduct_d_preconditioner *precond,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dcumilugeneratesolverinfo(
    magma_minproduct_d_preconditioner *precond,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dapplycumilu_l(
    magma_minproduct_d_matrix b, 
    magma_minproduct_d_matrix *x, 
    magma_minproduct_d_preconditioner *precond,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dapplycumilu_r(
    magma_minproduct_d_matrix b, 
    magma_minproduct_d_matrix *x, 
    magma_minproduct_d_preconditioner *precond,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t
magma_minproduct_dcumiccsetup(
    magma_minproduct_d_matrix A, 
    magma_minproduct_d_preconditioner *precond,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dcumicgeneratesolverinfo(
    magma_minproduct_d_preconditioner *precond,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dapplycumicc_l(
    magma_minproduct_d_matrix b, 
    magma_minproduct_d_matrix *x, 
    magma_minproduct_d_preconditioner *precond,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dapplycumicc_r(
    magma_minproduct_d_matrix b, 
    magma_minproduct_d_matrix *x, 
    magma_minproduct_d_preconditioner *precond,
    magma_minproduct_queue_t queue );


// block-asynchronous iteration

magma_minproduct_int_t
magma_minproduct_dbajac_csr(
    magma_minproduct_int_t localiters,
    magma_minproduct_d_matrix D,
    magma_minproduct_d_matrix R,
    magma_minproduct_d_matrix b,
    magma_minproduct_d_matrix *x,
    magma_minproduct_queue_t queue );

/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_minproduct_SPARSE utility function definitions
*/

magma_minproduct_int_t
magma_minproduct_d_spmv(
    double alpha, 
    magma_minproduct_d_matrix A, 
    magma_minproduct_d_matrix x, 
    double beta, 
    magma_minproduct_d_matrix y,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dcustomspmv(
    double alpha, 
    magma_minproduct_d_matrix x, 
    double beta, 
    magma_minproduct_d_matrix y,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_d_spmv_shift(
    double alpha, 
    magma_minproduct_d_matrix A, 
    double lambda,
    magma_minproduct_d_matrix x, 
    double beta, 
    magma_minproduct_int_t offset, 
    magma_minproduct_int_t blocksize,
    magma_minproductIndex_ptr dadd_vecs, 
    magma_minproduct_d_matrix y,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dcuspmm(
    magma_minproduct_d_matrix A, 
    magma_minproduct_d_matrix B, 
    magma_minproduct_d_matrix *AB,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_d_spmm(
    double alpha, 
    magma_minproduct_d_matrix A,
    magma_minproduct_d_matrix B,
    magma_minproduct_d_matrix *C,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dsymbilu( 
    magma_minproduct_d_matrix *A, 
    magma_minproduct_int_t levels, 
    magma_minproduct_d_matrix *L, 
    magma_minproduct_d_matrix *U,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dcuspaxpy(
    magma_minproductDouble_ptr alpha, magma_minproduct_d_matrix A, 
    magma_minproductDouble_ptr beta, magma_minproduct_d_matrix B, 
    magma_minproduct_d_matrix *AB,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_d_precond(
    magma_minproduct_d_matrix A, 
    magma_minproduct_d_matrix b, magma_minproduct_d_matrix *x,
    magma_minproduct_d_preconditioner *precond,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_d_solver(
    magma_minproduct_d_matrix A, magma_minproduct_d_matrix b, 
    magma_minproduct_d_matrix *x, magma_minproduct_dopts *zopts,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_d_precondsetup(
    magma_minproduct_d_matrix A, magma_minproduct_d_matrix b, 
    magma_minproduct_d_preconditioner *precond,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_d_applyprecond(
    magma_minproduct_d_matrix A, magma_minproduct_d_matrix b, 
    magma_minproduct_d_matrix *x, magma_minproduct_d_preconditioner *precond,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t
magma_minproduct_d_applyprecond_left(
    magma_minproduct_d_matrix A, magma_minproduct_d_matrix b, 
    magma_minproduct_d_matrix *x, magma_minproduct_d_preconditioner *precond,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t
magma_minproduct_d_applyprecond_right(
    magma_minproduct_d_matrix A, magma_minproduct_d_matrix b, 
    magma_minproduct_d_matrix *x, magma_minproduct_d_preconditioner *precond,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_d_initP2P(
    magma_minproduct_int_t *bandwidth_benchmark,
    magma_minproduct_int_t *num_gpus,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dcompact(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda,
    double *dnorms, double tol, 
    magma_minproduct_int_t *activeMask, magma_minproduct_int_t *cBlockSize,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dcompactActive(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDouble_ptr dA, magma_minproduct_int_t ldda, 
    magma_minproduct_int_t *active,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dmlumerge(    
    magma_minproduct_d_matrix L, 
    magma_minproduct_d_matrix U,
    magma_minproduct_d_matrix *A, 
    magma_minproduct_queue_t queue );


/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_minproduct_SPARSE BLAS function definitions
*/
magma_minproduct_int_t 
magma_minproduct_dgecsrmv(
    magma_minproduct_trans_t transA,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    double alpha,
    magma_minproductDouble_ptr dval,
    magma_minproductIndex_ptr drowptr,
    magma_minproductIndex_ptr dcolind,
    magma_minproductDouble_ptr dx,
    double beta,
    magma_minproductDouble_ptr dy,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_dgecsrmv_shift(
    magma_minproduct_trans_t transA,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    double alpha,
    double lambda,
    magma_minproductDouble_ptr dval,
    magma_minproductIndex_ptr drowptr,
    magma_minproductIndex_ptr dcolind,
    magma_minproductDouble_ptr dx,
    double beta,
    int offset,
    int blocksize,
    magma_minproductIndex_ptr dadd_rows,
    magma_minproductDouble_ptr dy,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_dmgecsrmv(
    magma_minproduct_trans_t transA,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproduct_int_t num_vecs,
    double alpha,
    magma_minproductDouble_ptr dval,
    magma_minproductIndex_ptr drowptr,
    magma_minproductIndex_ptr dcolind,
    magma_minproductDouble_ptr dx,
    double beta,
    magma_minproductDouble_ptr dy,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_dgeellmv(
    magma_minproduct_trans_t transA,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproduct_int_t nnz_per_row,
    double alpha,
    magma_minproductDouble_ptr dval,
    magma_minproductIndex_ptr dcolind,
    magma_minproductDouble_ptr dx,
    double beta,
    magma_minproductDouble_ptr dy,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_dgeellmv_shift(
    magma_minproduct_trans_t transA,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproduct_int_t nnz_per_row,
    double alpha,
    double lambda,
    magma_minproductDouble_ptr dval,
    magma_minproductIndex_ptr dcolind,
    magma_minproductDouble_ptr dx,
    double beta,
    int offset,
    int blocksize,
    magma_minproductIndex_ptr dadd_rows,
    magma_minproductDouble_ptr dy,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t 
magma_minproduct_dmgeellmv(
    magma_minproduct_trans_t transA,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproduct_int_t num_vecs,
    magma_minproduct_int_t nnz_per_row,
    double alpha,
    magma_minproductDouble_ptr dval,
    magma_minproductIndex_ptr dcolind,
    magma_minproductDouble_ptr dx,
    double beta,
    magma_minproductDouble_ptr dy,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t 
magma_minproduct_dgeelltmv(
    magma_minproduct_trans_t transA,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproduct_int_t nnz_per_row,
    double alpha,
    magma_minproductDouble_ptr dval,
    magma_minproductIndex_ptr dcolind,
    magma_minproductDouble_ptr dx,
    double beta,
    magma_minproductDouble_ptr dy,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_dgeelltmv_shift(
    magma_minproduct_trans_t transA,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproduct_int_t nnz_per_row,
    double alpha,
    double lambda,
    magma_minproductDouble_ptr dval,
    magma_minproductIndex_ptr dcolind,
    magma_minproductDouble_ptr dx,
    double beta,
    int offset,
    int blocksize,
    magma_minproductIndex_ptr dadd_rows,
    magma_minproductDouble_ptr dy,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t 
magma_minproduct_dmgeelltmv(
    magma_minproduct_trans_t transA,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproduct_int_t num_vecs,
    magma_minproduct_int_t nnz_per_row,
    double alpha,
    magma_minproductDouble_ptr dval,
    magma_minproductIndex_ptr dcolind,
    magma_minproductDouble_ptr dx,
    double beta,
    magma_minproductDouble_ptr dy,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_dgeellrtmv(
    magma_minproduct_trans_t transA,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproduct_int_t nnz_per_row,
    double alpha,
    magma_minproductDouble_ptr dval,
    magma_minproductIndex_ptr dcolind,
    magma_minproductIndex_ptr drowlength,
    magma_minproductDouble_ptr dx,
    double beta,
    magma_minproductDouble_ptr dy,
    magma_minproduct_int_t num_threads,
    magma_minproduct_int_t threads_per_row,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_dgesellcmv(
    magma_minproduct_trans_t transA,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproduct_int_t blocksize,
    magma_minproduct_int_t slices,
    magma_minproduct_int_t alignment,
    double alpha,
    magma_minproductDouble_ptr dval,
    magma_minproductIndex_ptr dcolind,
    magma_minproductIndex_ptr drowptr,
    magma_minproductDouble_ptr dx,
    double beta,
    magma_minproductDouble_ptr dy,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dgesellpmv(
    magma_minproduct_trans_t transA,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproduct_int_t blocksize,
    magma_minproduct_int_t slices,
    magma_minproduct_int_t alignment,
    double alpha,
    magma_minproductDouble_ptr dval,
    magma_minproductIndex_ptr dcolind,
    magma_minproductIndex_ptr drowptr,
    magma_minproductDouble_ptr dx,
    double beta,
    magma_minproductDouble_ptr dy,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dmgesellpmv(
    magma_minproduct_trans_t transA,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproduct_int_t num_vecs,
    magma_minproduct_int_t blocksize,
    magma_minproduct_int_t slices,
    magma_minproduct_int_t alignment,
    double alpha,
    magma_minproductDouble_ptr dval,
    magma_minproductIndex_ptr dcolind,
    magma_minproductIndex_ptr drowptr,
    magma_minproductDouble_ptr dx,
    double beta,
    magma_minproductDouble_ptr dy,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dmgesellpmv_blocked(
    magma_minproduct_trans_t transA,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproduct_int_t num_vecs,
    magma_minproduct_int_t blocksize,
    magma_minproduct_int_t slices,
    magma_minproduct_int_t alignment,
    double alpha,
    magma_minproductDouble_ptr dval,
    magma_minproductIndex_ptr dcolind,
    magma_minproductIndex_ptr drowptr,
    magma_minproductDouble_ptr dx,
    double beta,
    magma_minproductDouble_ptr dy,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t
magma_minproduct_dmergedgs(
    magma_minproduct_int_t n, 
    magma_minproduct_int_t ldh,
    magma_minproduct_int_t k, 
    magma_minproductDouble_ptr dv, 
    magma_minproductDouble_ptr dr,
    magma_minproductDouble_ptr dskp,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dcopyscale(    
    int n, 
    int k,
    magma_minproductDouble_ptr dr, 
    magma_minproductDouble_ptr dv,
    magma_minproductDouble_ptr dskp,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dnrm2scale(    
    int m, 
    magma_minproductDouble_ptr dr,    
    int lddr, 
    double *drnorm,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t
magma_minproduct_djacobisetup_vector_gpu(
    int num_rows, 
    magma_minproduct_d_matrix b, 
    magma_minproduct_d_matrix d, 
    magma_minproduct_d_matrix c,
    magma_minproduct_d_matrix *x,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t
magma_minproduct_djacobi_diagscal(    
    int num_rows, 
    magma_minproduct_d_matrix d, 
    magma_minproduct_d_matrix b, 
    magma_minproduct_d_matrix *c,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_djacobiupdate(
    magma_minproduct_d_matrix t, 
    magma_minproduct_d_matrix b, 
    magma_minproduct_d_matrix d, 
    magma_minproduct_d_matrix *x,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_djacobispmvupdate(
    magma_minproduct_int_t maxiter,
    magma_minproduct_d_matrix A, 
    magma_minproduct_d_matrix t, 
    magma_minproduct_d_matrix b, 
    magma_minproduct_d_matrix d, 
    magma_minproduct_d_matrix *x,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_djacobispmvupdate_bw(
    magma_minproduct_int_t maxiter,
    magma_minproduct_d_matrix A, 
    magma_minproduct_d_matrix t, 
    magma_minproduct_d_matrix b, 
    magma_minproduct_d_matrix d, 
    magma_minproduct_d_matrix *x,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_djacobispmvupdateselect(
    magma_minproduct_int_t maxiter,
    magma_minproduct_int_t num_updates,
    magma_minproduct_index_t *indices,
    magma_minproduct_d_matrix A,
    magma_minproduct_d_matrix t, 
    magma_minproduct_d_matrix b, 
    magma_minproduct_d_matrix d, 
    magma_minproduct_d_matrix tmp, 
    magma_minproduct_d_matrix *x,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_djacobisetup_diagscal(
    magma_minproduct_d_matrix A, magma_minproduct_d_matrix *d,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t
magma_minproduct_dbicgmerge1(    
    int n, 
    magma_minproductDouble_ptr dskp,
    magma_minproductDouble_ptr dv, 
    magma_minproductDouble_ptr dr, 
    magma_minproductDouble_ptr dp,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t
magma_minproduct_dbicgmerge2(
    int n, 
    magma_minproductDouble_ptr dskp, 
    magma_minproductDouble_ptr dr,
    magma_minproductDouble_ptr dv, 
    magma_minproductDouble_ptr ds,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dbicgmerge3(
    int n, 
    magma_minproductDouble_ptr dskp, 
    magma_minproductDouble_ptr dp,
    magma_minproductDouble_ptr ds,
    magma_minproductDouble_ptr dt,
    magma_minproductDouble_ptr dx, 
    magma_minproductDouble_ptr dr,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dbicgmerge4(
    int type, 
    magma_minproductDouble_ptr dskp,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dcgmerge_spmv1( 
    magma_minproduct_d_matrix A,
    magma_minproductDouble_ptr d1,
    magma_minproductDouble_ptr d2,
    magma_minproductDouble_ptr dd,
    magma_minproductDouble_ptr dz,
    magma_minproductDouble_ptr dskp,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dcgmerge_xrbeta( 
    int n,
    magma_minproductDouble_ptr d1,
    magma_minproductDouble_ptr d2,
    magma_minproductDouble_ptr dx,
    magma_minproductDouble_ptr dr,
    magma_minproductDouble_ptr dd,
    magma_minproductDouble_ptr dz, 
    magma_minproductDouble_ptr dskp,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dmdotc(
    magma_minproduct_int_t n, 
    magma_minproduct_int_t k, 
    magma_minproductDouble_ptr dv, 
    magma_minproductDouble_ptr dr,
    magma_minproductDouble_ptr dd1,
    magma_minproductDouble_ptr dd2,
    magma_minproductDouble_ptr dskp,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dgemvmdot(
    int n, 
    int k, 
    magma_minproductDouble_ptr dv, 
    magma_minproductDouble_ptr dr,
    magma_minproductDouble_ptr dd1,
    magma_minproductDouble_ptr dd2,
    magma_minproductDouble_ptr dskp,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dbicgmerge_spmv1( 
    magma_minproduct_d_matrix A,
    magma_minproductDouble_ptr dd1,
    magma_minproductDouble_ptr dd2,
    magma_minproductDouble_ptr dp,
    magma_minproductDouble_ptr dr,
    magma_minproductDouble_ptr dv,
    magma_minproductDouble_ptr dskp,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dbicgmerge_spmv2( 
    magma_minproduct_d_matrix A,
    magma_minproductDouble_ptr dd1,
    magma_minproductDouble_ptr dd2,
    magma_minproductDouble_ptr ds,
    magma_minproductDouble_ptr dt,
    magma_minproductDouble_ptr dskp,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dbicgmerge_xrbeta( 
    int n,
    magma_minproductDouble_ptr dd1,
    magma_minproductDouble_ptr dd2,
    magma_minproductDouble_ptr drr,
    magma_minproductDouble_ptr dr,
    magma_minproductDouble_ptr dp,
    magma_minproductDouble_ptr ds,
    magma_minproductDouble_ptr dt,
    magma_minproductDouble_ptr dx, 
    magma_minproductDouble_ptr dskp,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dbcsrswp(
    magma_minproduct_int_t n,
    magma_minproduct_int_t size_b, 
    magma_minproduct_int_t *ipiv,
    magma_minproductDouble_ptr dx,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dbcsrtrsv(
    magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t r_blocks,
    magma_minproduct_int_t c_blocks,
    magma_minproduct_int_t size_b, 
    magma_minproductDouble_ptr dA,
    magma_minproduct_index_t *blockinfo, 
    magma_minproductDouble_ptr dx,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dbcsrvalcpy(
    magma_minproduct_int_t size_b, 
    magma_minproduct_int_t num_blocks, 
    magma_minproduct_int_t num_zero_blocks, 
    magma_minproductDouble_ptr *dAval, 
    magma_minproductDouble_ptr *dBval,
    magma_minproductDouble_ptr *dBval2,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dbcsrluegemm(
    magma_minproduct_int_t size_b, 
    magma_minproduct_int_t num_block_rows,
    magma_minproduct_int_t kblocks,
    magma_minproductDouble_ptr *dA, 
    magma_minproductDouble_ptr *dB, 
    magma_minproductDouble_ptr *dC,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dbcsrlupivloc(
    magma_minproduct_int_t size_b, 
    magma_minproduct_int_t kblocks,
    magma_minproductDouble_ptr *dA, 
    magma_minproduct_int_t *ipiv,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_dbcsrblockinfo5(
    magma_minproduct_int_t lustep,
    magma_minproduct_int_t num_blocks, 
    magma_minproduct_int_t c_blocks, 
    magma_minproduct_int_t size_b,
    magma_minproduct_index_t *blockinfo,
    magma_minproductDouble_ptr dval,
    magma_minproductDouble_ptr *AII,
    magma_minproduct_queue_t queue );

 
#ifdef __cplusplus
}
#endif

#undef PRECISION_d
#endif /* MAGMA_minproductSPARSE_D_H */
