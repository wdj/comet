/*
 -- MAGMA_minproduct (version 1.1) --
 Univ. of Tennessee, Knoxville
 Univ. of California, Berkeley
 Univ. of Colorado, Denver
 @date

 @generated from magma_minproductsparse_z.h normal z -> s, Mon May  4 17:11:25 2015
 @author Hartwig Anzt
*/

#ifndef MAGMA_minproductSPARSE_S_H
#define MAGMA_minproductSPARSE_S_H

#include "magma_minproduct_types.h"
#include "magma_minproductsparse_types.h"

#define PRECISION_s


#ifdef __cplusplus
extern "C" {
#endif


/* ////////////////////////////////////////////////////////////////////////////
 -- For backwards compatability, map old (1.6.1) to new (1.6.2) function names
*/

#define magma_minproduct_s_mtranspose  magma_minproduct_smtranspose
#define magma_minproduct_s_mtransfer   magma_minproduct_smtransfer
#define magma_minproduct_s_vtransfer   magma_minproduct_smtransfer
#define magma_minproduct_s_mconvert    magma_minproduct_smconvert
#define magma_minproduct_s_vinit       magma_minproduct_svinit
#define magma_minproduct_s_vvisu       magma_minproduct_sprint_vector
#define magma_minproduct_s_vread       magma_minproduct_svread
#define magma_minproduct_s_vspread     magma_minproduct_svspread
#define magma_minproduct_s_mvisu       magma_minproduct_sprint_matrix
#define magma_minproduct_s_mfree       magma_minproduct_smfree
#define magma_minproduct_s_vfree       magma_minproduct_smfree
#define write_s_csr_mtx     magma_minproduct_swrite_csr_mtx
#define write_s_csrtomtx    magma_minproduct_swrite_csrtomtx
#define print_s_csr         magma_minproduct_sprint_csr_mtx


/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_minproduct_SPARSE Auxiliary functions
*/


magma_minproduct_int_t
magma_minproduct_sparse_opts( 
    int argc, 
    char** argv, 
    magma_minproduct_sopts *opts, 
    int *matrices, 
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
read_s_csr_from_binary( 
    magma_minproduct_int_t* n_row, 
    magma_minproduct_int_t* n_col, 
    magma_minproduct_int_t* nnz, 
    float **val, 
    magma_minproduct_index_t **row, 
    magma_minproduct_index_t **col,
    const char * filename,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
read_s_csr_from_mtx( 
    magma_minproduct_storage_t *type, 
    magma_minproduct_location_t *location,
    magma_minproduct_int_t* n_row, 
    magma_minproduct_int_t* n_col, 
    magma_minproduct_int_t* nnz, 
    float **val, 
    magma_minproduct_index_t **row, 
    magma_minproduct_index_t **col, 
    const char *filename,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_s_csr_mtx( 
    magma_minproduct_s_matrix *A, 
    const char *filename,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_scsrset( 
    magma_minproduct_int_t m, 
    magma_minproduct_int_t n, 
    magma_minproduct_index_t *row, 
    magma_minproduct_index_t *col, 
    float *val,
    magma_minproduct_s_matrix *A,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_scsrget( 
    magma_minproduct_s_matrix A,
    magma_minproduct_int_t *m, 
    magma_minproduct_int_t *n, 
    magma_minproduct_index_t **row, 
    magma_minproduct_index_t **col, 
    float **val,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t 
magma_minproduct_svset( 
    magma_minproduct_int_t m, 
    magma_minproduct_int_t n, 
    float *val,
    magma_minproduct_s_matrix *v,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_svget( 
    magma_minproduct_s_matrix v,
    magma_minproduct_int_t *m, 
    magma_minproduct_int_t *n, 
    float **val,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_svset_dev( 
    magma_minproduct_int_t m, 
    magma_minproduct_int_t n, 
    magma_minproductFloat_ptr val,
    magma_minproduct_s_matrix *v,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_svget_dev( 
    magma_minproduct_s_matrix v,
    magma_minproduct_int_t *m, 
    magma_minproduct_int_t *n, 
    magma_minproductFloat_ptr *val,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t 
magma_minproduct_s_csr_mtxsymm( 
    magma_minproduct_s_matrix *A, 
    const char *filename,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_s_csr_compressor( 
    float ** val, 
    magma_minproduct_index_t ** row, 
    magma_minproduct_index_t ** col, 
    float ** valn, 
    magma_minproduct_index_t ** rown, 
    magma_minproduct_index_t ** coln, 
    magma_minproduct_int_t *n,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_smcsrcompressor( 
    magma_minproduct_s_matrix *A,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_smcsrcompressor_gpu( 
    magma_minproduct_s_matrix *A,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_svtranspose( 
    magma_minproduct_s_matrix x,
    magma_minproduct_s_matrix *y,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_s_cucsrtranspose( 
    magma_minproduct_s_matrix A, 
    magma_minproduct_s_matrix *B,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
s_transpose_csr( 
    magma_minproduct_int_t n_rows, 
    magma_minproduct_int_t n_cols, 
    magma_minproduct_int_t nnz,
    float *val, 
    magma_minproduct_index_t *row, 
    magma_minproduct_index_t *col, 
    magma_minproduct_int_t *new_n_rows, 
    magma_minproduct_int_t *new_n_cols, 
    magma_minproduct_int_t *new_nnz, 
    float **new_val, 
    magma_minproduct_index_t **new_row, 
    magma_minproduct_index_t **new_col,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_scsrsplit( 
    magma_minproduct_int_t bsize,
    magma_minproduct_s_matrix A,
    magma_minproduct_s_matrix *D,
    magma_minproduct_s_matrix *R,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_smscale( 
    magma_minproduct_s_matrix *A, 
    magma_minproduct_scale_t scaling,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_smdiff( 
    magma_minproduct_s_matrix A, 
    magma_minproduct_s_matrix B, 
 real_Double_t *res,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_smdiagadd( 
    magma_minproduct_s_matrix *A, 
    float add,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_smsort( 
    magma_minproduct_s_matrix *A,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_sindexsort(
    magma_minproduct_index_t *x, 
    magma_minproduct_int_t first,
    magma_minproduct_int_t last,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_sdomainoverlap(
    magma_minproduct_index_t num_rows,
    magma_minproduct_index_t *num_indices,
    magma_minproduct_index_t *rowptr,
    magma_minproduct_index_t *colidx,
    magma_minproduct_index_t *x,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_ssymbilu( 
    magma_minproduct_s_matrix *A, 
    magma_minproduct_int_t levels,
    magma_minproduct_s_matrix *L,
    magma_minproduct_s_matrix *U,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t 
magma_minproduct_swrite_csr_mtx( 
    magma_minproduct_s_matrix A,
    magma_minproduct_order_t MajorType,
 const char *filename,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_swrite_csrtomtx( 
    magma_minproduct_s_matrix A,
    const char *filename,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_sprint_csr( 
    magma_minproduct_int_t n_row, 
    magma_minproduct_int_t n_col, 
    magma_minproduct_int_t nnz, 
    float **val, 
    magma_minproduct_index_t **row, 
    magma_minproduct_index_t **col,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_sprint_csr_mtx( 
    magma_minproduct_int_t n_row, 
    magma_minproduct_int_t n_col, 
    magma_minproduct_int_t nnz, 
    float **val, 
    magma_minproduct_index_t **row, 
    magma_minproduct_index_t **col, 
    magma_minproduct_order_t MajorType,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t 
magma_minproduct_smtranspose(
    magma_minproduct_s_matrix A, 
    magma_minproduct_s_matrix *B,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t 
magma_minproduct_smtransfer(
    magma_minproduct_s_matrix A, 
    magma_minproduct_s_matrix *B, 
    magma_minproduct_location_t src, 
    magma_minproduct_location_t dst,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_smconvert(
    magma_minproduct_s_matrix A, 
    magma_minproduct_s_matrix *B, 
    magma_minproduct_storage_t old_format, 
    magma_minproduct_storage_t new_format,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t
magma_minproduct_svinit(
    magma_minproduct_s_matrix *x, 
    magma_minproduct_location_t memory_location,
    magma_minproduct_int_t num_rows, 
    magma_minproduct_int_t num_cols,
    float values,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_sprint_vector(
    magma_minproduct_s_matrix x, 
    magma_minproduct_int_t offset, 
    magma_minproduct_int_t displaylength,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_svread(
    magma_minproduct_s_matrix *x, 
    magma_minproduct_int_t length,
    char * filename,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_svspread(
    magma_minproduct_s_matrix *x, 
    const char * filename,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_sprint_matrix(
    magma_minproduct_s_matrix A,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_sdiameter(
    magma_minproduct_s_matrix *A,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_srowentries(
    magma_minproduct_s_matrix *A,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_smfree(
    magma_minproduct_s_matrix *A,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_sresidual(
    magma_minproduct_s_matrix A, 
    magma_minproduct_s_matrix b, 
    magma_minproduct_s_matrix x, 
    float *res,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_sresidualvec(
    magma_minproduct_s_matrix A,
    magma_minproduct_s_matrix b,
    magma_minproduct_s_matrix x,
    magma_minproduct_s_matrix *r,
    float *res,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_smgenerator(
    magma_minproduct_int_t n,
    magma_minproduct_int_t offdiags,
    magma_minproduct_index_t *diag_offset,
    float *diag_vals,
    magma_minproduct_s_matrix *A,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_sm_27stencil(
    magma_minproduct_int_t n,
    magma_minproduct_s_matrix *A,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_sm_5stencil(
    magma_minproduct_int_t n,
    magma_minproduct_s_matrix *A,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_ssolverinfo(
    magma_minproduct_s_solver_par *solver_par, 
    magma_minproduct_s_preconditioner *precond_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_ssolverinfo_init(
    magma_minproduct_s_solver_par *solver_par, 
    magma_minproduct_s_preconditioner *precond,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_seigensolverinfo_init(
    magma_minproduct_s_solver_par *solver_par,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t
magma_minproduct_ssolverinfo_free(
    magma_minproduct_s_solver_par *solver_par, 
    magma_minproduct_s_preconditioner *precond,
    magma_minproduct_queue_t queue );



/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_minproduct_SPARSE iterative incomplete factorizations
*/


magma_minproduct_int_t
magma_minproduct_siterilusetup( 
    magma_minproduct_s_matrix A, 
    magma_minproduct_s_matrix b,                                 
    magma_minproduct_s_preconditioner *precond,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_sitericsetup( 
    magma_minproduct_s_matrix A, 
    magma_minproduct_s_matrix b, 
    magma_minproduct_s_preconditioner *precond,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_sitericupdate( 
    magma_minproduct_s_matrix A, 
    magma_minproduct_s_preconditioner *precond, 
    magma_minproduct_int_t updates,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_sapplyiteric_l( 
    magma_minproduct_s_matrix b, 
    magma_minproduct_s_matrix *x, 
    magma_minproduct_s_preconditioner *precond,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_sapplyiteric_r( 
    magma_minproduct_s_matrix b, 
    magma_minproduct_s_matrix *x, 
    magma_minproduct_s_preconditioner *precond,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t
magma_minproduct_siterilu_csr( 
    magma_minproduct_s_matrix A,
    magma_minproduct_s_matrix L,
    magma_minproduct_s_matrix U,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_siteric_csr( 
    magma_minproduct_s_matrix A,
    magma_minproduct_s_matrix A_CSR,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_sfrobenius( 
    magma_minproduct_s_matrix A, 
    magma_minproduct_s_matrix B, 
    real_Double_t *res,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_snonlinres(   
    magma_minproduct_s_matrix A, 
    magma_minproduct_s_matrix L,
    magma_minproduct_s_matrix U, 
    magma_minproduct_s_matrix *LU, 
    real_Double_t *res,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_silures(   
    magma_minproduct_s_matrix A, 
    magma_minproduct_s_matrix L,
    magma_minproduct_s_matrix U, 
    magma_minproduct_s_matrix *LU, 
    real_Double_t *res,
    real_Double_t *nonlinres,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_sicres(       
    magma_minproduct_s_matrix A, 
    magma_minproduct_s_matrix C,
    magma_minproduct_s_matrix CT, 
    magma_minproduct_s_matrix *LU, 
    real_Double_t *res,
    real_Double_t *nonlinres,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_sinitguess( 
    magma_minproduct_s_matrix A, 
    magma_minproduct_s_matrix *L, 
    magma_minproduct_s_matrix *U,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_sinitrecursiveLU( 
    magma_minproduct_s_matrix A, 
    magma_minproduct_s_matrix *B,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_smLdiagadd( 
    magma_minproduct_s_matrix *L,
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
magma_minproduct_scg(
    magma_minproduct_s_matrix A, magma_minproduct_s_matrix b, 
    magma_minproduct_s_matrix *x, magma_minproduct_s_solver_par *solver_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_scg_res(
    magma_minproduct_s_matrix A, magma_minproduct_s_matrix b, 
    magma_minproduct_s_matrix *x, magma_minproduct_s_solver_par *solver_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_scg_merge(
    magma_minproduct_s_matrix A, magma_minproduct_s_matrix b, 
    magma_minproduct_s_matrix *x, magma_minproduct_s_solver_par *solver_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_sgmres(
    magma_minproduct_s_matrix A, magma_minproduct_s_matrix b, 
    magma_minproduct_s_matrix *x, magma_minproduct_s_solver_par *solver_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_sbicgstab(
    magma_minproduct_s_matrix A, magma_minproduct_s_matrix b, magma_minproduct_s_matrix *x, 
    magma_minproduct_s_solver_par *solver_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_sbicgstab_merge(
    magma_minproduct_s_matrix A, magma_minproduct_s_matrix b, 
    magma_minproduct_s_matrix *x, magma_minproduct_s_solver_par *solver_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_sbicgstab_merge2(
    magma_minproduct_s_matrix A, magma_minproduct_s_matrix b, 
    magma_minproduct_s_matrix *x, magma_minproduct_s_solver_par *solver_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_spcg(
    magma_minproduct_s_matrix A, magma_minproduct_s_matrix b, 
    magma_minproduct_s_matrix *x, magma_minproduct_s_solver_par *solver_par, 
    magma_minproduct_s_preconditioner *precond_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_sbpcg(
    magma_minproduct_s_matrix A, magma_minproduct_s_matrix b, 
    magma_minproduct_s_matrix *x, magma_minproduct_s_solver_par *solver_par, 
    magma_minproduct_s_preconditioner *precond_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_spbicgstab(
    magma_minproduct_s_matrix A, magma_minproduct_s_matrix b, 
    magma_minproduct_s_matrix *x, magma_minproduct_s_solver_par *solver_par, 
    magma_minproduct_s_preconditioner *precond_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_spgmres(
    magma_minproduct_s_matrix A, magma_minproduct_s_matrix b, 
    magma_minproduct_s_matrix *x, magma_minproduct_s_solver_par *solver_par, 
    magma_minproduct_s_preconditioner *precond_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_sfgmres(
    magma_minproduct_s_matrix A, magma_minproduct_s_matrix b, 
    magma_minproduct_s_matrix *x, magma_minproduct_s_solver_par *solver_par, 
    magma_minproduct_s_preconditioner *precond_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_sjacobi(
    magma_minproduct_s_matrix A, magma_minproduct_s_matrix b, 
    magma_minproduct_s_matrix *x, magma_minproduct_s_solver_par *solver_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_sjacobidomainoverlap(
    magma_minproduct_s_matrix A, 
    magma_minproduct_s_matrix b, 
    magma_minproduct_s_matrix *x,  
    magma_minproduct_s_solver_par *solver_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_sbaiter(
    magma_minproduct_s_matrix A, magma_minproduct_s_matrix b, 
    magma_minproduct_s_matrix *x, magma_minproduct_s_solver_par *solver_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_siterref(
    magma_minproduct_s_matrix A, magma_minproduct_s_matrix b, 
    magma_minproduct_s_matrix *x, magma_minproduct_s_solver_par *solver_par, 
    magma_minproduct_s_preconditioner *precond_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_silu(
    magma_minproduct_s_matrix A, magma_minproduct_s_matrix b, 
    magma_minproduct_s_matrix *x, magma_minproduct_s_solver_par *solver_par, 
    magma_minproduct_int_t *ipiv,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_sbcsrlu(
    magma_minproduct_s_matrix A, magma_minproduct_s_matrix b, 
    magma_minproduct_s_matrix *x, magma_minproduct_s_solver_par *solver_par,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t
magma_minproduct_sbcsrlutrf(
    magma_minproduct_s_matrix A, 
    magma_minproduct_s_matrix *M,
    magma_minproduct_int_t *ipiv, 
    magma_minproduct_int_t version,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_sbcsrlusv(
    magma_minproduct_s_matrix A, magma_minproduct_s_matrix b, 
    magma_minproduct_s_matrix *x, magma_minproduct_s_solver_par *solver_par, 
    magma_minproduct_int_t *ipiv,
    magma_minproduct_queue_t queue );



magma_minproduct_int_t
magma_minproduct_silucg(
    magma_minproduct_s_matrix A, magma_minproduct_s_matrix b, 
    magma_minproduct_s_matrix *x, magma_minproduct_s_solver_par *solver_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_silugmres(
    magma_minproduct_s_matrix A, magma_minproduct_s_matrix b, 
    magma_minproduct_s_matrix *x, magma_minproduct_s_solver_par *solver_par,
    magma_minproduct_queue_t queue ); 


magma_minproduct_int_t
magma_minproduct_slobpcg_shift(
    magma_minproduct_int_t num_rows,
    magma_minproduct_int_t num_vecs, 
    magma_minproduct_int_t shift,
    magma_minproductFloat_ptr x,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_slobpcg_res(
    magma_minproduct_int_t num_rows,
    magma_minproduct_int_t num_vecs, 
    float *evalues, 
    magma_minproductFloat_ptr X,
    magma_minproductFloat_ptr R, 
    float *res,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_slobpcg_maxpy(
    magma_minproduct_int_t num_rows,
    magma_minproduct_int_t num_vecs, 
    magma_minproductFloat_ptr X,
    magma_minproductFloat_ptr Y,
    magma_minproduct_queue_t queue );


/*/////////////////////////////////////////////////////////////////////////////
 -- MAGMA_minproduct_SPARSE eigensolvers (Data on GPU)
*/
magma_minproduct_int_t
magma_minproduct_slobpcg(
    magma_minproduct_s_matrix A, 
    magma_minproduct_s_solver_par *solver_par,
    magma_minproduct_s_preconditioner *precond_par, 
    magma_minproduct_queue_t queue );




/*/////////////////////////////////////////////////////////////////////////////
 -- MAGMA_minproduct_SPARSE preconditioners (Data on GPU)
*/

magma_minproduct_int_t
magma_minproduct_sjacobisetup(
    magma_minproduct_s_matrix A, 
    magma_minproduct_s_matrix b, 
    magma_minproduct_s_matrix *M, 
    magma_minproduct_s_matrix *c,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_sjacobisetup_matrix(
    magma_minproduct_s_matrix A, 
    magma_minproduct_s_matrix *M, 
    magma_minproduct_s_matrix *d,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_sjacobisetup_vector(
    magma_minproduct_s_matrix b, 
    magma_minproduct_s_matrix d, 
    magma_minproduct_s_matrix *c,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_sjacobiiter(
    magma_minproduct_s_matrix M, 
    magma_minproduct_s_matrix c, 
    magma_minproduct_s_matrix *x, 
    magma_minproduct_s_solver_par *solver_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_sjacobiiter_precond( 
    magma_minproduct_s_matrix M, 
    magma_minproduct_s_matrix *x, 
    magma_minproduct_s_solver_par *solver_par, 
    magma_minproduct_s_preconditioner *precond,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_sjacobiiter_sys(
    magma_minproduct_s_matrix A, 
    magma_minproduct_s_matrix b, 
    magma_minproduct_s_matrix d, 
    magma_minproduct_s_matrix t, 
    magma_minproduct_s_matrix *x,  
    magma_minproduct_s_solver_par *solver_par,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_spastixsetup(
    magma_minproduct_s_matrix A, magma_minproduct_s_matrix b,
    magma_minproduct_s_preconditioner *precond,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t
magma_minproduct_sapplypastix(
    magma_minproduct_s_matrix b, magma_minproduct_s_matrix *x, 
    magma_minproduct_s_preconditioner *precond,
    magma_minproduct_queue_t queue );


// custom preconditioner
magma_minproduct_int_t
magma_minproduct_sapplycustomprecond_l(
    magma_minproduct_s_matrix b, 
    magma_minproduct_s_matrix *x, 
    magma_minproduct_s_preconditioner *precond,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_sapplycustomprecond_r(
    magma_minproduct_s_matrix b, 
    magma_minproduct_s_matrix *x, 
    magma_minproduct_s_preconditioner *precond,
    magma_minproduct_queue_t queue );


// CUSPARSE preconditioner

magma_minproduct_int_t
magma_minproduct_scuilusetup(
    magma_minproduct_s_matrix A, magma_minproduct_s_preconditioner *precond,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_sapplycuilu_l(
    magma_minproduct_s_matrix b, magma_minproduct_s_matrix *x, 
    magma_minproduct_s_preconditioner *precond,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_sapplycuilu_r(
    magma_minproduct_s_matrix b, magma_minproduct_s_matrix *x, 
    magma_minproduct_s_preconditioner *precond,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t
magma_minproduct_scuiccsetup(
    magma_minproduct_s_matrix A, magma_minproduct_s_preconditioner *precond,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_sapplycuicc_l(
    magma_minproduct_s_matrix b, magma_minproduct_s_matrix *x, 
    magma_minproduct_s_preconditioner *precond,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_sapplycuicc_r(
    magma_minproduct_s_matrix b, magma_minproduct_s_matrix *x, 
    magma_minproduct_s_preconditioner *precond,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t
magma_minproduct_scumilusetup(
    magma_minproduct_s_matrix A, 
    magma_minproduct_s_preconditioner *precond,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_scumilugeneratesolverinfo(
    magma_minproduct_s_preconditioner *precond,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_sapplycumilu_l(
    magma_minproduct_s_matrix b, 
    magma_minproduct_s_matrix *x, 
    magma_minproduct_s_preconditioner *precond,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_sapplycumilu_r(
    magma_minproduct_s_matrix b, 
    magma_minproduct_s_matrix *x, 
    magma_minproduct_s_preconditioner *precond,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t
magma_minproduct_scumiccsetup(
    magma_minproduct_s_matrix A, 
    magma_minproduct_s_preconditioner *precond,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_scumicgeneratesolverinfo(
    magma_minproduct_s_preconditioner *precond,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_sapplycumicc_l(
    magma_minproduct_s_matrix b, 
    magma_minproduct_s_matrix *x, 
    magma_minproduct_s_preconditioner *precond,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_sapplycumicc_r(
    magma_minproduct_s_matrix b, 
    magma_minproduct_s_matrix *x, 
    magma_minproduct_s_preconditioner *precond,
    magma_minproduct_queue_t queue );


// block-asynchronous iteration

magma_minproduct_int_t
magma_minproduct_sbajac_csr(
    magma_minproduct_int_t localiters,
    magma_minproduct_s_matrix D,
    magma_minproduct_s_matrix R,
    magma_minproduct_s_matrix b,
    magma_minproduct_s_matrix *x,
    magma_minproduct_queue_t queue );

/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_minproduct_SPARSE utility function definitions
*/

magma_minproduct_int_t
magma_minproduct_s_spmv(
    float alpha, 
    magma_minproduct_s_matrix A, 
    magma_minproduct_s_matrix x, 
    float beta, 
    magma_minproduct_s_matrix y,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_scustomspmv(
    float alpha, 
    magma_minproduct_s_matrix x, 
    float beta, 
    magma_minproduct_s_matrix y,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_s_spmv_shift(
    float alpha, 
    magma_minproduct_s_matrix A, 
    float lambda,
    magma_minproduct_s_matrix x, 
    float beta, 
    magma_minproduct_int_t offset, 
    magma_minproduct_int_t blocksize,
    magma_minproductIndex_ptr dadd_vecs, 
    magma_minproduct_s_matrix y,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_scuspmm(
    magma_minproduct_s_matrix A, 
    magma_minproduct_s_matrix B, 
    magma_minproduct_s_matrix *AB,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_s_spmm(
    float alpha, 
    magma_minproduct_s_matrix A,
    magma_minproduct_s_matrix B,
    magma_minproduct_s_matrix *C,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_ssymbilu( 
    magma_minproduct_s_matrix *A, 
    magma_minproduct_int_t levels, 
    magma_minproduct_s_matrix *L, 
    magma_minproduct_s_matrix *U,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_scuspaxpy(
    magma_minproductFloat_ptr alpha, magma_minproduct_s_matrix A, 
    magma_minproductFloat_ptr beta, magma_minproduct_s_matrix B, 
    magma_minproduct_s_matrix *AB,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_s_precond(
    magma_minproduct_s_matrix A, 
    magma_minproduct_s_matrix b, magma_minproduct_s_matrix *x,
    magma_minproduct_s_preconditioner *precond,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_s_solver(
    magma_minproduct_s_matrix A, magma_minproduct_s_matrix b, 
    magma_minproduct_s_matrix *x, magma_minproduct_sopts *zopts,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_s_precondsetup(
    magma_minproduct_s_matrix A, magma_minproduct_s_matrix b, 
    magma_minproduct_s_preconditioner *precond,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_s_applyprecond(
    magma_minproduct_s_matrix A, magma_minproduct_s_matrix b, 
    magma_minproduct_s_matrix *x, magma_minproduct_s_preconditioner *precond,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t
magma_minproduct_s_applyprecond_left(
    magma_minproduct_s_matrix A, magma_minproduct_s_matrix b, 
    magma_minproduct_s_matrix *x, magma_minproduct_s_preconditioner *precond,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t
magma_minproduct_s_applyprecond_right(
    magma_minproduct_s_matrix A, magma_minproduct_s_matrix b, 
    magma_minproduct_s_matrix *x, magma_minproduct_s_preconditioner *precond,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_s_initP2P(
    magma_minproduct_int_t *bandwidth_benchmark,
    magma_minproduct_int_t *num_gpus,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_scompact(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda,
    float *dnorms, float tol, 
    magma_minproduct_int_t *activeMask, magma_minproduct_int_t *cBlockSize,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_scompactActive(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloat_ptr dA, magma_minproduct_int_t ldda, 
    magma_minproduct_int_t *active,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_smlumerge(    
    magma_minproduct_s_matrix L, 
    magma_minproduct_s_matrix U,
    magma_minproduct_s_matrix *A, 
    magma_minproduct_queue_t queue );


/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_minproduct_SPARSE BLAS function definitions
*/
magma_minproduct_int_t 
magma_minproduct_sgecsrmv(
    magma_minproduct_trans_t transA,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    float alpha,
    magma_minproductFloat_ptr dval,
    magma_minproductIndex_ptr drowptr,
    magma_minproductIndex_ptr dcolind,
    magma_minproductFloat_ptr dx,
    float beta,
    magma_minproductFloat_ptr dy,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_sgecsrmv_shift(
    magma_minproduct_trans_t transA,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    float alpha,
    float lambda,
    magma_minproductFloat_ptr dval,
    magma_minproductIndex_ptr drowptr,
    magma_minproductIndex_ptr dcolind,
    magma_minproductFloat_ptr dx,
    float beta,
    int offset,
    int blocksize,
    magma_minproductIndex_ptr dadd_rows,
    magma_minproductFloat_ptr dy,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_smgecsrmv(
    magma_minproduct_trans_t transA,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproduct_int_t num_vecs,
    float alpha,
    magma_minproductFloat_ptr dval,
    magma_minproductIndex_ptr drowptr,
    magma_minproductIndex_ptr dcolind,
    magma_minproductFloat_ptr dx,
    float beta,
    magma_minproductFloat_ptr dy,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_sgeellmv(
    magma_minproduct_trans_t transA,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproduct_int_t nnz_per_row,
    float alpha,
    magma_minproductFloat_ptr dval,
    magma_minproductIndex_ptr dcolind,
    magma_minproductFloat_ptr dx,
    float beta,
    magma_minproductFloat_ptr dy,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_sgeellmv_shift(
    magma_minproduct_trans_t transA,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproduct_int_t nnz_per_row,
    float alpha,
    float lambda,
    magma_minproductFloat_ptr dval,
    magma_minproductIndex_ptr dcolind,
    magma_minproductFloat_ptr dx,
    float beta,
    int offset,
    int blocksize,
    magma_minproductIndex_ptr dadd_rows,
    magma_minproductFloat_ptr dy,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t 
magma_minproduct_smgeellmv(
    magma_minproduct_trans_t transA,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproduct_int_t num_vecs,
    magma_minproduct_int_t nnz_per_row,
    float alpha,
    magma_minproductFloat_ptr dval,
    magma_minproductIndex_ptr dcolind,
    magma_minproductFloat_ptr dx,
    float beta,
    magma_minproductFloat_ptr dy,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t 
magma_minproduct_sgeelltmv(
    magma_minproduct_trans_t transA,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproduct_int_t nnz_per_row,
    float alpha,
    magma_minproductFloat_ptr dval,
    magma_minproductIndex_ptr dcolind,
    magma_minproductFloat_ptr dx,
    float beta,
    magma_minproductFloat_ptr dy,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_sgeelltmv_shift(
    magma_minproduct_trans_t transA,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproduct_int_t nnz_per_row,
    float alpha,
    float lambda,
    magma_minproductFloat_ptr dval,
    magma_minproductIndex_ptr dcolind,
    magma_minproductFloat_ptr dx,
    float beta,
    int offset,
    int blocksize,
    magma_minproductIndex_ptr dadd_rows,
    magma_minproductFloat_ptr dy,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t 
magma_minproduct_smgeelltmv(
    magma_minproduct_trans_t transA,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproduct_int_t num_vecs,
    magma_minproduct_int_t nnz_per_row,
    float alpha,
    magma_minproductFloat_ptr dval,
    magma_minproductIndex_ptr dcolind,
    magma_minproductFloat_ptr dx,
    float beta,
    magma_minproductFloat_ptr dy,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_sgeellrtmv(
    magma_minproduct_trans_t transA,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproduct_int_t nnz_per_row,
    float alpha,
    magma_minproductFloat_ptr dval,
    magma_minproductIndex_ptr dcolind,
    magma_minproductIndex_ptr drowlength,
    magma_minproductFloat_ptr dx,
    float beta,
    magma_minproductFloat_ptr dy,
    magma_minproduct_int_t num_threads,
    magma_minproduct_int_t threads_per_row,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t 
magma_minproduct_sgesellcmv(
    magma_minproduct_trans_t transA,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproduct_int_t blocksize,
    magma_minproduct_int_t slices,
    magma_minproduct_int_t alignment,
    float alpha,
    magma_minproductFloat_ptr dval,
    magma_minproductIndex_ptr dcolind,
    magma_minproductIndex_ptr drowptr,
    magma_minproductFloat_ptr dx,
    float beta,
    magma_minproductFloat_ptr dy,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_sgesellpmv(
    magma_minproduct_trans_t transA,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproduct_int_t blocksize,
    magma_minproduct_int_t slices,
    magma_minproduct_int_t alignment,
    float alpha,
    magma_minproductFloat_ptr dval,
    magma_minproductIndex_ptr dcolind,
    magma_minproductIndex_ptr drowptr,
    magma_minproductFloat_ptr dx,
    float beta,
    magma_minproductFloat_ptr dy,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_smgesellpmv(
    magma_minproduct_trans_t transA,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproduct_int_t num_vecs,
    magma_minproduct_int_t blocksize,
    magma_minproduct_int_t slices,
    magma_minproduct_int_t alignment,
    float alpha,
    magma_minproductFloat_ptr dval,
    magma_minproductIndex_ptr dcolind,
    magma_minproductIndex_ptr drowptr,
    magma_minproductFloat_ptr dx,
    float beta,
    magma_minproductFloat_ptr dy,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_smgesellpmv_blocked(
    magma_minproduct_trans_t transA,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproduct_int_t num_vecs,
    magma_minproduct_int_t blocksize,
    magma_minproduct_int_t slices,
    magma_minproduct_int_t alignment,
    float alpha,
    magma_minproductFloat_ptr dval,
    magma_minproductIndex_ptr dcolind,
    magma_minproductIndex_ptr drowptr,
    magma_minproductFloat_ptr dx,
    float beta,
    magma_minproductFloat_ptr dy,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t
magma_minproduct_smergedgs(
    magma_minproduct_int_t n, 
    magma_minproduct_int_t ldh,
    magma_minproduct_int_t k, 
    magma_minproductFloat_ptr dv, 
    magma_minproductFloat_ptr dr,
    magma_minproductFloat_ptr dskp,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_scopyscale(    
    int n, 
    int k,
    magma_minproductFloat_ptr dr, 
    magma_minproductFloat_ptr dv,
    magma_minproductFloat_ptr dskp,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_snrm2scale(    
    int m, 
    magma_minproductFloat_ptr dr,    
    int lddr, 
    float *drnorm,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t
magma_minproduct_sjacobisetup_vector_gpu(
    int num_rows, 
    magma_minproduct_s_matrix b, 
    magma_minproduct_s_matrix d, 
    magma_minproduct_s_matrix c,
    magma_minproduct_s_matrix *x,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t
magma_minproduct_sjacobi_diagscal(    
    int num_rows, 
    magma_minproduct_s_matrix d, 
    magma_minproduct_s_matrix b, 
    magma_minproduct_s_matrix *c,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_sjacobiupdate(
    magma_minproduct_s_matrix t, 
    magma_minproduct_s_matrix b, 
    magma_minproduct_s_matrix d, 
    magma_minproduct_s_matrix *x,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_sjacobispmvupdate(
    magma_minproduct_int_t maxiter,
    magma_minproduct_s_matrix A, 
    magma_minproduct_s_matrix t, 
    magma_minproduct_s_matrix b, 
    magma_minproduct_s_matrix d, 
    magma_minproduct_s_matrix *x,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_sjacobispmvupdate_bw(
    magma_minproduct_int_t maxiter,
    magma_minproduct_s_matrix A, 
    magma_minproduct_s_matrix t, 
    magma_minproduct_s_matrix b, 
    magma_minproduct_s_matrix d, 
    magma_minproduct_s_matrix *x,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_sjacobispmvupdateselect(
    magma_minproduct_int_t maxiter,
    magma_minproduct_int_t num_updates,
    magma_minproduct_index_t *indices,
    magma_minproduct_s_matrix A,
    magma_minproduct_s_matrix t, 
    magma_minproduct_s_matrix b, 
    magma_minproduct_s_matrix d, 
    magma_minproduct_s_matrix tmp, 
    magma_minproduct_s_matrix *x,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_sjacobisetup_diagscal(
    magma_minproduct_s_matrix A, magma_minproduct_s_matrix *d,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t
magma_minproduct_sbicgmerge1(    
    int n, 
    magma_minproductFloat_ptr dskp,
    magma_minproductFloat_ptr dv, 
    magma_minproductFloat_ptr dr, 
    magma_minproductFloat_ptr dp,
    magma_minproduct_queue_t queue );


magma_minproduct_int_t
magma_minproduct_sbicgmerge2(
    int n, 
    magma_minproductFloat_ptr dskp, 
    magma_minproductFloat_ptr dr,
    magma_minproductFloat_ptr dv, 
    magma_minproductFloat_ptr ds,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_sbicgmerge3(
    int n, 
    magma_minproductFloat_ptr dskp, 
    magma_minproductFloat_ptr dp,
    magma_minproductFloat_ptr ds,
    magma_minproductFloat_ptr dt,
    magma_minproductFloat_ptr dx, 
    magma_minproductFloat_ptr dr,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_sbicgmerge4(
    int type, 
    magma_minproductFloat_ptr dskp,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_scgmerge_spmv1( 
    magma_minproduct_s_matrix A,
    magma_minproductFloat_ptr d1,
    magma_minproductFloat_ptr d2,
    magma_minproductFloat_ptr dd,
    magma_minproductFloat_ptr dz,
    magma_minproductFloat_ptr dskp,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_scgmerge_xrbeta( 
    int n,
    magma_minproductFloat_ptr d1,
    magma_minproductFloat_ptr d2,
    magma_minproductFloat_ptr dx,
    magma_minproductFloat_ptr dr,
    magma_minproductFloat_ptr dd,
    magma_minproductFloat_ptr dz, 
    magma_minproductFloat_ptr dskp,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_smdotc(
    magma_minproduct_int_t n, 
    magma_minproduct_int_t k, 
    magma_minproductFloat_ptr dv, 
    magma_minproductFloat_ptr dr,
    magma_minproductFloat_ptr dd1,
    magma_minproductFloat_ptr dd2,
    magma_minproductFloat_ptr dskp,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_sgemvmdot(
    int n, 
    int k, 
    magma_minproductFloat_ptr dv, 
    magma_minproductFloat_ptr dr,
    magma_minproductFloat_ptr dd1,
    magma_minproductFloat_ptr dd2,
    magma_minproductFloat_ptr dskp,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_sbicgmerge_spmv1( 
    magma_minproduct_s_matrix A,
    magma_minproductFloat_ptr dd1,
    magma_minproductFloat_ptr dd2,
    magma_minproductFloat_ptr dp,
    magma_minproductFloat_ptr dr,
    magma_minproductFloat_ptr dv,
    magma_minproductFloat_ptr dskp,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_sbicgmerge_spmv2( 
    magma_minproduct_s_matrix A,
    magma_minproductFloat_ptr dd1,
    magma_minproductFloat_ptr dd2,
    magma_minproductFloat_ptr ds,
    magma_minproductFloat_ptr dt,
    magma_minproductFloat_ptr dskp,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_sbicgmerge_xrbeta( 
    int n,
    magma_minproductFloat_ptr dd1,
    magma_minproductFloat_ptr dd2,
    magma_minproductFloat_ptr drr,
    magma_minproductFloat_ptr dr,
    magma_minproductFloat_ptr dp,
    magma_minproductFloat_ptr ds,
    magma_minproductFloat_ptr dt,
    magma_minproductFloat_ptr dx, 
    magma_minproductFloat_ptr dskp,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_sbcsrswp(
    magma_minproduct_int_t n,
    magma_minproduct_int_t size_b, 
    magma_minproduct_int_t *ipiv,
    magma_minproductFloat_ptr dx,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_sbcsrtrsv(
    magma_minproduct_uplo_t uplo,
    magma_minproduct_int_t r_blocks,
    magma_minproduct_int_t c_blocks,
    magma_minproduct_int_t size_b, 
    magma_minproductFloat_ptr dA,
    magma_minproduct_index_t *blockinfo, 
    magma_minproductFloat_ptr dx,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_sbcsrvalcpy(
    magma_minproduct_int_t size_b, 
    magma_minproduct_int_t num_blocks, 
    magma_minproduct_int_t num_zero_blocks, 
    magma_minproductFloat_ptr *dAval, 
    magma_minproductFloat_ptr *dBval,
    magma_minproductFloat_ptr *dBval2,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_sbcsrluegemm(
    magma_minproduct_int_t size_b, 
    magma_minproduct_int_t num_block_rows,
    magma_minproduct_int_t kblocks,
    magma_minproductFloat_ptr *dA, 
    magma_minproductFloat_ptr *dB, 
    magma_minproductFloat_ptr *dC,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_sbcsrlupivloc(
    magma_minproduct_int_t size_b, 
    magma_minproduct_int_t kblocks,
    magma_minproductFloat_ptr *dA, 
    magma_minproduct_int_t *ipiv,
    magma_minproduct_queue_t queue );

magma_minproduct_int_t
magma_minproduct_sbcsrblockinfo5(
    magma_minproduct_int_t lustep,
    magma_minproduct_int_t num_blocks, 
    magma_minproduct_int_t c_blocks, 
    magma_minproduct_int_t size_b,
    magma_minproduct_index_t *blockinfo,
    magma_minproductFloat_ptr dval,
    magma_minproductFloat_ptr *AII,
    magma_minproduct_queue_t queue );

 
#ifdef __cplusplus
}
#endif

#undef PRECISION_s
#endif /* MAGMA_minproductSPARSE_S_H */
