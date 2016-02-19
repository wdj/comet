/*
 -- MAGMA_tally3 (version 1.1) --
 Univ. of Tennessee, Knoxville
 Univ. of California, Berkeley
 Univ. of Colorado, Denver
 @date

 @generated from magma_tally3sparse_z.h normal z -> s, Mon May  4 17:11:25 2015
 @author Hartwig Anzt
*/

#ifndef MAGMA_tally3SPARSE_S_H
#define MAGMA_tally3SPARSE_S_H

#include "magma_tally3_types.h"
#include "magma_tally3sparse_types.h"

#define PRECISION_s


#ifdef __cplusplus
extern "C" {
#endif


/* ////////////////////////////////////////////////////////////////////////////
 -- For backwards compatability, map old (1.6.1) to new (1.6.2) function names
*/

#define magma_tally3_s_mtranspose  magma_tally3_smtranspose
#define magma_tally3_s_mtransfer   magma_tally3_smtransfer
#define magma_tally3_s_vtransfer   magma_tally3_smtransfer
#define magma_tally3_s_mconvert    magma_tally3_smconvert
#define magma_tally3_s_vinit       magma_tally3_svinit
#define magma_tally3_s_vvisu       magma_tally3_sprint_vector
#define magma_tally3_s_vread       magma_tally3_svread
#define magma_tally3_s_vspread     magma_tally3_svspread
#define magma_tally3_s_mvisu       magma_tally3_sprint_matrix
#define magma_tally3_s_mfree       magma_tally3_smfree
#define magma_tally3_s_vfree       magma_tally3_smfree
#define write_s_csr_mtx     magma_tally3_swrite_csr_mtx
#define write_s_csrtomtx    magma_tally3_swrite_csrtomtx
#define print_s_csr         magma_tally3_sprint_csr_mtx


/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally3_SPARSE Auxiliary functions
*/


magma_tally3_int_t
magma_tally3_sparse_opts( 
    int argc, 
    char** argv, 
    magma_tally3_sopts *opts, 
    int *matrices, 
    magma_tally3_queue_t queue );

magma_tally3_int_t 
read_s_csr_from_binary( 
    magma_tally3_int_t* n_row, 
    magma_tally3_int_t* n_col, 
    magma_tally3_int_t* nnz, 
    float **val, 
    magma_tally3_index_t **row, 
    magma_tally3_index_t **col,
    const char * filename,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
read_s_csr_from_mtx( 
    magma_tally3_storage_t *type, 
    magma_tally3_location_t *location,
    magma_tally3_int_t* n_row, 
    magma_tally3_int_t* n_col, 
    magma_tally3_int_t* nnz, 
    float **val, 
    magma_tally3_index_t **row, 
    magma_tally3_index_t **col, 
    const char *filename,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_s_csr_mtx( 
    magma_tally3_s_matrix *A, 
    const char *filename,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_scsrset( 
    magma_tally3_int_t m, 
    magma_tally3_int_t n, 
    magma_tally3_index_t *row, 
    magma_tally3_index_t *col, 
    float *val,
    magma_tally3_s_matrix *A,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_scsrget( 
    magma_tally3_s_matrix A,
    magma_tally3_int_t *m, 
    magma_tally3_int_t *n, 
    magma_tally3_index_t **row, 
    magma_tally3_index_t **col, 
    float **val,
    magma_tally3_queue_t queue );


magma_tally3_int_t 
magma_tally3_svset( 
    magma_tally3_int_t m, 
    magma_tally3_int_t n, 
    float *val,
    magma_tally3_s_matrix *v,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_svget( 
    magma_tally3_s_matrix v,
    magma_tally3_int_t *m, 
    magma_tally3_int_t *n, 
    float **val,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_svset_dev( 
    magma_tally3_int_t m, 
    magma_tally3_int_t n, 
    magma_tally3Float_ptr val,
    magma_tally3_s_matrix *v,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_svget_dev( 
    magma_tally3_s_matrix v,
    magma_tally3_int_t *m, 
    magma_tally3_int_t *n, 
    magma_tally3Float_ptr *val,
    magma_tally3_queue_t queue );


magma_tally3_int_t 
magma_tally3_s_csr_mtxsymm( 
    magma_tally3_s_matrix *A, 
    const char *filename,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_s_csr_compressor( 
    float ** val, 
    magma_tally3_index_t ** row, 
    magma_tally3_index_t ** col, 
    float ** valn, 
    magma_tally3_index_t ** rown, 
    magma_tally3_index_t ** coln, 
    magma_tally3_int_t *n,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_smcsrcompressor( 
    magma_tally3_s_matrix *A,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_smcsrcompressor_gpu( 
    magma_tally3_s_matrix *A,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_svtranspose( 
    magma_tally3_s_matrix x,
    magma_tally3_s_matrix *y,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_s_cucsrtranspose( 
    magma_tally3_s_matrix A, 
    magma_tally3_s_matrix *B,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
s_transpose_csr( 
    magma_tally3_int_t n_rows, 
    magma_tally3_int_t n_cols, 
    magma_tally3_int_t nnz,
    float *val, 
    magma_tally3_index_t *row, 
    magma_tally3_index_t *col, 
    magma_tally3_int_t *new_n_rows, 
    magma_tally3_int_t *new_n_cols, 
    magma_tally3_int_t *new_nnz, 
    float **new_val, 
    magma_tally3_index_t **new_row, 
    magma_tally3_index_t **new_col,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_scsrsplit( 
    magma_tally3_int_t bsize,
    magma_tally3_s_matrix A,
    magma_tally3_s_matrix *D,
    magma_tally3_s_matrix *R,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_smscale( 
    magma_tally3_s_matrix *A, 
    magma_tally3_scale_t scaling,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_smdiff( 
    magma_tally3_s_matrix A, 
    magma_tally3_s_matrix B, 
 real_Double_t *res,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_smdiagadd( 
    magma_tally3_s_matrix *A, 
    float add,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_smsort( 
    magma_tally3_s_matrix *A,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_sindexsort(
    magma_tally3_index_t *x, 
    magma_tally3_int_t first,
    magma_tally3_int_t last,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_sdomainoverlap(
    magma_tally3_index_t num_rows,
    magma_tally3_index_t *num_indices,
    magma_tally3_index_t *rowptr,
    magma_tally3_index_t *colidx,
    magma_tally3_index_t *x,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_ssymbilu( 
    magma_tally3_s_matrix *A, 
    magma_tally3_int_t levels,
    magma_tally3_s_matrix *L,
    magma_tally3_s_matrix *U,
    magma_tally3_queue_t queue );


magma_tally3_int_t 
magma_tally3_swrite_csr_mtx( 
    magma_tally3_s_matrix A,
    magma_tally3_order_t MajorType,
 const char *filename,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_swrite_csrtomtx( 
    magma_tally3_s_matrix A,
    const char *filename,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_sprint_csr( 
    magma_tally3_int_t n_row, 
    magma_tally3_int_t n_col, 
    magma_tally3_int_t nnz, 
    float **val, 
    magma_tally3_index_t **row, 
    magma_tally3_index_t **col,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_sprint_csr_mtx( 
    magma_tally3_int_t n_row, 
    magma_tally3_int_t n_col, 
    magma_tally3_int_t nnz, 
    float **val, 
    magma_tally3_index_t **row, 
    magma_tally3_index_t **col, 
    magma_tally3_order_t MajorType,
    magma_tally3_queue_t queue );


magma_tally3_int_t 
magma_tally3_smtranspose(
    magma_tally3_s_matrix A, 
    magma_tally3_s_matrix *B,
    magma_tally3_queue_t queue );


magma_tally3_int_t 
magma_tally3_smtransfer(
    magma_tally3_s_matrix A, 
    magma_tally3_s_matrix *B, 
    magma_tally3_location_t src, 
    magma_tally3_location_t dst,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_smconvert(
    magma_tally3_s_matrix A, 
    magma_tally3_s_matrix *B, 
    magma_tally3_storage_t old_format, 
    magma_tally3_storage_t new_format,
    magma_tally3_queue_t queue );


magma_tally3_int_t
magma_tally3_svinit(
    magma_tally3_s_matrix *x, 
    magma_tally3_location_t memory_location,
    magma_tally3_int_t num_rows, 
    magma_tally3_int_t num_cols,
    float values,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_sprint_vector(
    magma_tally3_s_matrix x, 
    magma_tally3_int_t offset, 
    magma_tally3_int_t displaylength,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_svread(
    magma_tally3_s_matrix *x, 
    magma_tally3_int_t length,
    char * filename,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_svspread(
    magma_tally3_s_matrix *x, 
    const char * filename,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_sprint_matrix(
    magma_tally3_s_matrix A,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_sdiameter(
    magma_tally3_s_matrix *A,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_srowentries(
    magma_tally3_s_matrix *A,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_smfree(
    magma_tally3_s_matrix *A,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_sresidual(
    magma_tally3_s_matrix A, 
    magma_tally3_s_matrix b, 
    magma_tally3_s_matrix x, 
    float *res,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_sresidualvec(
    magma_tally3_s_matrix A,
    magma_tally3_s_matrix b,
    magma_tally3_s_matrix x,
    magma_tally3_s_matrix *r,
    float *res,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_smgenerator(
    magma_tally3_int_t n,
    magma_tally3_int_t offdiags,
    magma_tally3_index_t *diag_offset,
    float *diag_vals,
    magma_tally3_s_matrix *A,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_sm_27stencil(
    magma_tally3_int_t n,
    magma_tally3_s_matrix *A,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_sm_5stencil(
    magma_tally3_int_t n,
    magma_tally3_s_matrix *A,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_ssolverinfo(
    magma_tally3_s_solver_par *solver_par, 
    magma_tally3_s_preconditioner *precond_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_ssolverinfo_init(
    magma_tally3_s_solver_par *solver_par, 
    magma_tally3_s_preconditioner *precond,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_seigensolverinfo_init(
    magma_tally3_s_solver_par *solver_par,
    magma_tally3_queue_t queue );


magma_tally3_int_t
magma_tally3_ssolverinfo_free(
    magma_tally3_s_solver_par *solver_par, 
    magma_tally3_s_preconditioner *precond,
    magma_tally3_queue_t queue );



/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally3_SPARSE iterative incomplete factorizations
*/


magma_tally3_int_t
magma_tally3_siterilusetup( 
    magma_tally3_s_matrix A, 
    magma_tally3_s_matrix b,                                 
    magma_tally3_s_preconditioner *precond,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_sitericsetup( 
    magma_tally3_s_matrix A, 
    magma_tally3_s_matrix b, 
    magma_tally3_s_preconditioner *precond,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_sitericupdate( 
    magma_tally3_s_matrix A, 
    magma_tally3_s_preconditioner *precond, 
    magma_tally3_int_t updates,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_sapplyiteric_l( 
    magma_tally3_s_matrix b, 
    magma_tally3_s_matrix *x, 
    magma_tally3_s_preconditioner *precond,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_sapplyiteric_r( 
    magma_tally3_s_matrix b, 
    magma_tally3_s_matrix *x, 
    magma_tally3_s_preconditioner *precond,
    magma_tally3_queue_t queue );


magma_tally3_int_t
magma_tally3_siterilu_csr( 
    magma_tally3_s_matrix A,
    magma_tally3_s_matrix L,
    magma_tally3_s_matrix U,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_siteric_csr( 
    magma_tally3_s_matrix A,
    magma_tally3_s_matrix A_CSR,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_sfrobenius( 
    magma_tally3_s_matrix A, 
    magma_tally3_s_matrix B, 
    real_Double_t *res,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_snonlinres(   
    magma_tally3_s_matrix A, 
    magma_tally3_s_matrix L,
    magma_tally3_s_matrix U, 
    magma_tally3_s_matrix *LU, 
    real_Double_t *res,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_silures(   
    magma_tally3_s_matrix A, 
    magma_tally3_s_matrix L,
    magma_tally3_s_matrix U, 
    magma_tally3_s_matrix *LU, 
    real_Double_t *res,
    real_Double_t *nonlinres,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_sicres(       
    magma_tally3_s_matrix A, 
    magma_tally3_s_matrix C,
    magma_tally3_s_matrix CT, 
    magma_tally3_s_matrix *LU, 
    real_Double_t *res,
    real_Double_t *nonlinres,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_sinitguess( 
    magma_tally3_s_matrix A, 
    magma_tally3_s_matrix *L, 
    magma_tally3_s_matrix *U,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_sinitrecursiveLU( 
    magma_tally3_s_matrix A, 
    magma_tally3_s_matrix *B,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_smLdiagadd( 
    magma_tally3_s_matrix *L,
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
magma_tally3_scg(
    magma_tally3_s_matrix A, magma_tally3_s_matrix b, 
    magma_tally3_s_matrix *x, magma_tally3_s_solver_par *solver_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_scg_res(
    magma_tally3_s_matrix A, magma_tally3_s_matrix b, 
    magma_tally3_s_matrix *x, magma_tally3_s_solver_par *solver_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_scg_merge(
    magma_tally3_s_matrix A, magma_tally3_s_matrix b, 
    magma_tally3_s_matrix *x, magma_tally3_s_solver_par *solver_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_sgmres(
    magma_tally3_s_matrix A, magma_tally3_s_matrix b, 
    magma_tally3_s_matrix *x, magma_tally3_s_solver_par *solver_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_sbicgstab(
    magma_tally3_s_matrix A, magma_tally3_s_matrix b, magma_tally3_s_matrix *x, 
    magma_tally3_s_solver_par *solver_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_sbicgstab_merge(
    magma_tally3_s_matrix A, magma_tally3_s_matrix b, 
    magma_tally3_s_matrix *x, magma_tally3_s_solver_par *solver_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_sbicgstab_merge2(
    magma_tally3_s_matrix A, magma_tally3_s_matrix b, 
    magma_tally3_s_matrix *x, magma_tally3_s_solver_par *solver_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_spcg(
    magma_tally3_s_matrix A, magma_tally3_s_matrix b, 
    magma_tally3_s_matrix *x, magma_tally3_s_solver_par *solver_par, 
    magma_tally3_s_preconditioner *precond_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_sbpcg(
    magma_tally3_s_matrix A, magma_tally3_s_matrix b, 
    magma_tally3_s_matrix *x, magma_tally3_s_solver_par *solver_par, 
    magma_tally3_s_preconditioner *precond_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_spbicgstab(
    magma_tally3_s_matrix A, magma_tally3_s_matrix b, 
    magma_tally3_s_matrix *x, magma_tally3_s_solver_par *solver_par, 
    magma_tally3_s_preconditioner *precond_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_spgmres(
    magma_tally3_s_matrix A, magma_tally3_s_matrix b, 
    magma_tally3_s_matrix *x, magma_tally3_s_solver_par *solver_par, 
    magma_tally3_s_preconditioner *precond_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_sfgmres(
    magma_tally3_s_matrix A, magma_tally3_s_matrix b, 
    magma_tally3_s_matrix *x, magma_tally3_s_solver_par *solver_par, 
    magma_tally3_s_preconditioner *precond_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_sjacobi(
    magma_tally3_s_matrix A, magma_tally3_s_matrix b, 
    magma_tally3_s_matrix *x, magma_tally3_s_solver_par *solver_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_sjacobidomainoverlap(
    magma_tally3_s_matrix A, 
    magma_tally3_s_matrix b, 
    magma_tally3_s_matrix *x,  
    magma_tally3_s_solver_par *solver_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_sbaiter(
    magma_tally3_s_matrix A, magma_tally3_s_matrix b, 
    magma_tally3_s_matrix *x, magma_tally3_s_solver_par *solver_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_siterref(
    magma_tally3_s_matrix A, magma_tally3_s_matrix b, 
    magma_tally3_s_matrix *x, magma_tally3_s_solver_par *solver_par, 
    magma_tally3_s_preconditioner *precond_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_silu(
    magma_tally3_s_matrix A, magma_tally3_s_matrix b, 
    magma_tally3_s_matrix *x, magma_tally3_s_solver_par *solver_par, 
    magma_tally3_int_t *ipiv,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_sbcsrlu(
    magma_tally3_s_matrix A, magma_tally3_s_matrix b, 
    magma_tally3_s_matrix *x, magma_tally3_s_solver_par *solver_par,
    magma_tally3_queue_t queue );


magma_tally3_int_t
magma_tally3_sbcsrlutrf(
    magma_tally3_s_matrix A, 
    magma_tally3_s_matrix *M,
    magma_tally3_int_t *ipiv, 
    magma_tally3_int_t version,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_sbcsrlusv(
    magma_tally3_s_matrix A, magma_tally3_s_matrix b, 
    magma_tally3_s_matrix *x, magma_tally3_s_solver_par *solver_par, 
    magma_tally3_int_t *ipiv,
    magma_tally3_queue_t queue );



magma_tally3_int_t
magma_tally3_silucg(
    magma_tally3_s_matrix A, magma_tally3_s_matrix b, 
    magma_tally3_s_matrix *x, magma_tally3_s_solver_par *solver_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_silugmres(
    magma_tally3_s_matrix A, magma_tally3_s_matrix b, 
    magma_tally3_s_matrix *x, magma_tally3_s_solver_par *solver_par,
    magma_tally3_queue_t queue ); 


magma_tally3_int_t
magma_tally3_slobpcg_shift(
    magma_tally3_int_t num_rows,
    magma_tally3_int_t num_vecs, 
    magma_tally3_int_t shift,
    magma_tally3Float_ptr x,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_slobpcg_res(
    magma_tally3_int_t num_rows,
    magma_tally3_int_t num_vecs, 
    float *evalues, 
    magma_tally3Float_ptr X,
    magma_tally3Float_ptr R, 
    float *res,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_slobpcg_maxpy(
    magma_tally3_int_t num_rows,
    magma_tally3_int_t num_vecs, 
    magma_tally3Float_ptr X,
    magma_tally3Float_ptr Y,
    magma_tally3_queue_t queue );


/*/////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally3_SPARSE eigensolvers (Data on GPU)
*/
magma_tally3_int_t
magma_tally3_slobpcg(
    magma_tally3_s_matrix A, 
    magma_tally3_s_solver_par *solver_par,
    magma_tally3_s_preconditioner *precond_par, 
    magma_tally3_queue_t queue );




/*/////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally3_SPARSE preconditioners (Data on GPU)
*/

magma_tally3_int_t
magma_tally3_sjacobisetup(
    magma_tally3_s_matrix A, 
    magma_tally3_s_matrix b, 
    magma_tally3_s_matrix *M, 
    magma_tally3_s_matrix *c,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_sjacobisetup_matrix(
    magma_tally3_s_matrix A, 
    magma_tally3_s_matrix *M, 
    magma_tally3_s_matrix *d,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_sjacobisetup_vector(
    magma_tally3_s_matrix b, 
    magma_tally3_s_matrix d, 
    magma_tally3_s_matrix *c,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_sjacobiiter(
    magma_tally3_s_matrix M, 
    magma_tally3_s_matrix c, 
    magma_tally3_s_matrix *x, 
    magma_tally3_s_solver_par *solver_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_sjacobiiter_precond( 
    magma_tally3_s_matrix M, 
    magma_tally3_s_matrix *x, 
    magma_tally3_s_solver_par *solver_par, 
    magma_tally3_s_preconditioner *precond,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_sjacobiiter_sys(
    magma_tally3_s_matrix A, 
    magma_tally3_s_matrix b, 
    magma_tally3_s_matrix d, 
    magma_tally3_s_matrix t, 
    magma_tally3_s_matrix *x,  
    magma_tally3_s_solver_par *solver_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_spastixsetup(
    magma_tally3_s_matrix A, magma_tally3_s_matrix b,
    magma_tally3_s_preconditioner *precond,
    magma_tally3_queue_t queue );


magma_tally3_int_t
magma_tally3_sapplypastix(
    magma_tally3_s_matrix b, magma_tally3_s_matrix *x, 
    magma_tally3_s_preconditioner *precond,
    magma_tally3_queue_t queue );


// custom preconditioner
magma_tally3_int_t
magma_tally3_sapplycustomprecond_l(
    magma_tally3_s_matrix b, 
    magma_tally3_s_matrix *x, 
    magma_tally3_s_preconditioner *precond,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_sapplycustomprecond_r(
    magma_tally3_s_matrix b, 
    magma_tally3_s_matrix *x, 
    magma_tally3_s_preconditioner *precond,
    magma_tally3_queue_t queue );


// CUSPARSE preconditioner

magma_tally3_int_t
magma_tally3_scuilusetup(
    magma_tally3_s_matrix A, magma_tally3_s_preconditioner *precond,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_sapplycuilu_l(
    magma_tally3_s_matrix b, magma_tally3_s_matrix *x, 
    magma_tally3_s_preconditioner *precond,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_sapplycuilu_r(
    magma_tally3_s_matrix b, magma_tally3_s_matrix *x, 
    magma_tally3_s_preconditioner *precond,
    magma_tally3_queue_t queue );


magma_tally3_int_t
magma_tally3_scuiccsetup(
    magma_tally3_s_matrix A, magma_tally3_s_preconditioner *precond,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_sapplycuicc_l(
    magma_tally3_s_matrix b, magma_tally3_s_matrix *x, 
    magma_tally3_s_preconditioner *precond,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_sapplycuicc_r(
    magma_tally3_s_matrix b, magma_tally3_s_matrix *x, 
    magma_tally3_s_preconditioner *precond,
    magma_tally3_queue_t queue );


magma_tally3_int_t
magma_tally3_scumilusetup(
    magma_tally3_s_matrix A, 
    magma_tally3_s_preconditioner *precond,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_scumilugeneratesolverinfo(
    magma_tally3_s_preconditioner *precond,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_sapplycumilu_l(
    magma_tally3_s_matrix b, 
    magma_tally3_s_matrix *x, 
    magma_tally3_s_preconditioner *precond,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_sapplycumilu_r(
    magma_tally3_s_matrix b, 
    magma_tally3_s_matrix *x, 
    magma_tally3_s_preconditioner *precond,
    magma_tally3_queue_t queue );


magma_tally3_int_t
magma_tally3_scumiccsetup(
    magma_tally3_s_matrix A, 
    magma_tally3_s_preconditioner *precond,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_scumicgeneratesolverinfo(
    magma_tally3_s_preconditioner *precond,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_sapplycumicc_l(
    magma_tally3_s_matrix b, 
    magma_tally3_s_matrix *x, 
    magma_tally3_s_preconditioner *precond,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_sapplycumicc_r(
    magma_tally3_s_matrix b, 
    magma_tally3_s_matrix *x, 
    magma_tally3_s_preconditioner *precond,
    magma_tally3_queue_t queue );


// block-asynchronous iteration

magma_tally3_int_t
magma_tally3_sbajac_csr(
    magma_tally3_int_t localiters,
    magma_tally3_s_matrix D,
    magma_tally3_s_matrix R,
    magma_tally3_s_matrix b,
    magma_tally3_s_matrix *x,
    magma_tally3_queue_t queue );

/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally3_SPARSE utility function definitions
*/

magma_tally3_int_t
magma_tally3_s_spmv(
    float alpha, 
    magma_tally3_s_matrix A, 
    magma_tally3_s_matrix x, 
    float beta, 
    magma_tally3_s_matrix y,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_scustomspmv(
    float alpha, 
    magma_tally3_s_matrix x, 
    float beta, 
    magma_tally3_s_matrix y,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_s_spmv_shift(
    float alpha, 
    magma_tally3_s_matrix A, 
    float lambda,
    magma_tally3_s_matrix x, 
    float beta, 
    magma_tally3_int_t offset, 
    magma_tally3_int_t blocksize,
    magma_tally3Index_ptr dadd_vecs, 
    magma_tally3_s_matrix y,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_scuspmm(
    magma_tally3_s_matrix A, 
    magma_tally3_s_matrix B, 
    magma_tally3_s_matrix *AB,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_s_spmm(
    float alpha, 
    magma_tally3_s_matrix A,
    magma_tally3_s_matrix B,
    magma_tally3_s_matrix *C,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_ssymbilu( 
    magma_tally3_s_matrix *A, 
    magma_tally3_int_t levels, 
    magma_tally3_s_matrix *L, 
    magma_tally3_s_matrix *U,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_scuspaxpy(
    magma_tally3Float_ptr alpha, magma_tally3_s_matrix A, 
    magma_tally3Float_ptr beta, magma_tally3_s_matrix B, 
    magma_tally3_s_matrix *AB,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_s_precond(
    magma_tally3_s_matrix A, 
    magma_tally3_s_matrix b, magma_tally3_s_matrix *x,
    magma_tally3_s_preconditioner *precond,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_s_solver(
    magma_tally3_s_matrix A, magma_tally3_s_matrix b, 
    magma_tally3_s_matrix *x, magma_tally3_sopts *zopts,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_s_precondsetup(
    magma_tally3_s_matrix A, magma_tally3_s_matrix b, 
    magma_tally3_s_preconditioner *precond,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_s_applyprecond(
    magma_tally3_s_matrix A, magma_tally3_s_matrix b, 
    magma_tally3_s_matrix *x, magma_tally3_s_preconditioner *precond,
    magma_tally3_queue_t queue );


magma_tally3_int_t
magma_tally3_s_applyprecond_left(
    magma_tally3_s_matrix A, magma_tally3_s_matrix b, 
    magma_tally3_s_matrix *x, magma_tally3_s_preconditioner *precond,
    magma_tally3_queue_t queue );


magma_tally3_int_t
magma_tally3_s_applyprecond_right(
    magma_tally3_s_matrix A, magma_tally3_s_matrix b, 
    magma_tally3_s_matrix *x, magma_tally3_s_preconditioner *precond,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_s_initP2P(
    magma_tally3_int_t *bandwidth_benchmark,
    magma_tally3_int_t *num_gpus,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_scompact(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda,
    float *dnorms, float tol, 
    magma_tally3_int_t *activeMask, magma_tally3_int_t *cBlockSize,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_scompactActive(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda, 
    magma_tally3_int_t *active,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_smlumerge(    
    magma_tally3_s_matrix L, 
    magma_tally3_s_matrix U,
    magma_tally3_s_matrix *A, 
    magma_tally3_queue_t queue );


/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally3_SPARSE BLAS function definitions
*/
magma_tally3_int_t 
magma_tally3_sgecsrmv(
    magma_tally3_trans_t transA,
    magma_tally3_int_t m, magma_tally3_int_t n,
    float alpha,
    magma_tally3Float_ptr dval,
    magma_tally3Index_ptr drowptr,
    magma_tally3Index_ptr dcolind,
    magma_tally3Float_ptr dx,
    float beta,
    magma_tally3Float_ptr dy,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_sgecsrmv_shift(
    magma_tally3_trans_t transA,
    magma_tally3_int_t m, magma_tally3_int_t n,
    float alpha,
    float lambda,
    magma_tally3Float_ptr dval,
    magma_tally3Index_ptr drowptr,
    magma_tally3Index_ptr dcolind,
    magma_tally3Float_ptr dx,
    float beta,
    int offset,
    int blocksize,
    magma_tally3Index_ptr dadd_rows,
    magma_tally3Float_ptr dy,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_smgecsrmv(
    magma_tally3_trans_t transA,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3_int_t num_vecs,
    float alpha,
    magma_tally3Float_ptr dval,
    magma_tally3Index_ptr drowptr,
    magma_tally3Index_ptr dcolind,
    magma_tally3Float_ptr dx,
    float beta,
    magma_tally3Float_ptr dy,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_sgeellmv(
    magma_tally3_trans_t transA,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3_int_t nnz_per_row,
    float alpha,
    magma_tally3Float_ptr dval,
    magma_tally3Index_ptr dcolind,
    magma_tally3Float_ptr dx,
    float beta,
    magma_tally3Float_ptr dy,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_sgeellmv_shift(
    magma_tally3_trans_t transA,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3_int_t nnz_per_row,
    float alpha,
    float lambda,
    magma_tally3Float_ptr dval,
    magma_tally3Index_ptr dcolind,
    magma_tally3Float_ptr dx,
    float beta,
    int offset,
    int blocksize,
    magma_tally3Index_ptr dadd_rows,
    magma_tally3Float_ptr dy,
    magma_tally3_queue_t queue );


magma_tally3_int_t 
magma_tally3_smgeellmv(
    magma_tally3_trans_t transA,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3_int_t num_vecs,
    magma_tally3_int_t nnz_per_row,
    float alpha,
    magma_tally3Float_ptr dval,
    magma_tally3Index_ptr dcolind,
    magma_tally3Float_ptr dx,
    float beta,
    magma_tally3Float_ptr dy,
    magma_tally3_queue_t queue );


magma_tally3_int_t 
magma_tally3_sgeelltmv(
    magma_tally3_trans_t transA,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3_int_t nnz_per_row,
    float alpha,
    magma_tally3Float_ptr dval,
    magma_tally3Index_ptr dcolind,
    magma_tally3Float_ptr dx,
    float beta,
    magma_tally3Float_ptr dy,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_sgeelltmv_shift(
    magma_tally3_trans_t transA,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3_int_t nnz_per_row,
    float alpha,
    float lambda,
    magma_tally3Float_ptr dval,
    magma_tally3Index_ptr dcolind,
    magma_tally3Float_ptr dx,
    float beta,
    int offset,
    int blocksize,
    magma_tally3Index_ptr dadd_rows,
    magma_tally3Float_ptr dy,
    magma_tally3_queue_t queue );


magma_tally3_int_t 
magma_tally3_smgeelltmv(
    magma_tally3_trans_t transA,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3_int_t num_vecs,
    magma_tally3_int_t nnz_per_row,
    float alpha,
    magma_tally3Float_ptr dval,
    magma_tally3Index_ptr dcolind,
    magma_tally3Float_ptr dx,
    float beta,
    magma_tally3Float_ptr dy,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_sgeellrtmv(
    magma_tally3_trans_t transA,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3_int_t nnz_per_row,
    float alpha,
    magma_tally3Float_ptr dval,
    magma_tally3Index_ptr dcolind,
    magma_tally3Index_ptr drowlength,
    magma_tally3Float_ptr dx,
    float beta,
    magma_tally3Float_ptr dy,
    magma_tally3_int_t num_threads,
    magma_tally3_int_t threads_per_row,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_sgesellcmv(
    magma_tally3_trans_t transA,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3_int_t blocksize,
    magma_tally3_int_t slices,
    magma_tally3_int_t alignment,
    float alpha,
    magma_tally3Float_ptr dval,
    magma_tally3Index_ptr dcolind,
    magma_tally3Index_ptr drowptr,
    magma_tally3Float_ptr dx,
    float beta,
    magma_tally3Float_ptr dy,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_sgesellpmv(
    magma_tally3_trans_t transA,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3_int_t blocksize,
    magma_tally3_int_t slices,
    magma_tally3_int_t alignment,
    float alpha,
    magma_tally3Float_ptr dval,
    magma_tally3Index_ptr dcolind,
    magma_tally3Index_ptr drowptr,
    magma_tally3Float_ptr dx,
    float beta,
    magma_tally3Float_ptr dy,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_smgesellpmv(
    magma_tally3_trans_t transA,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3_int_t num_vecs,
    magma_tally3_int_t blocksize,
    magma_tally3_int_t slices,
    magma_tally3_int_t alignment,
    float alpha,
    magma_tally3Float_ptr dval,
    magma_tally3Index_ptr dcolind,
    magma_tally3Index_ptr drowptr,
    magma_tally3Float_ptr dx,
    float beta,
    magma_tally3Float_ptr dy,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_smgesellpmv_blocked(
    magma_tally3_trans_t transA,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3_int_t num_vecs,
    magma_tally3_int_t blocksize,
    magma_tally3_int_t slices,
    magma_tally3_int_t alignment,
    float alpha,
    magma_tally3Float_ptr dval,
    magma_tally3Index_ptr dcolind,
    magma_tally3Index_ptr drowptr,
    magma_tally3Float_ptr dx,
    float beta,
    magma_tally3Float_ptr dy,
    magma_tally3_queue_t queue );


magma_tally3_int_t
magma_tally3_smergedgs(
    magma_tally3_int_t n, 
    magma_tally3_int_t ldh,
    magma_tally3_int_t k, 
    magma_tally3Float_ptr dv, 
    magma_tally3Float_ptr dr,
    magma_tally3Float_ptr dskp,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_scopyscale(    
    int n, 
    int k,
    magma_tally3Float_ptr dr, 
    magma_tally3Float_ptr dv,
    magma_tally3Float_ptr dskp,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_snrm2scale(    
    int m, 
    magma_tally3Float_ptr dr,    
    int lddr, 
    float *drnorm,
    magma_tally3_queue_t queue );


magma_tally3_int_t
magma_tally3_sjacobisetup_vector_gpu(
    int num_rows, 
    magma_tally3_s_matrix b, 
    magma_tally3_s_matrix d, 
    magma_tally3_s_matrix c,
    magma_tally3_s_matrix *x,
    magma_tally3_queue_t queue );


magma_tally3_int_t
magma_tally3_sjacobi_diagscal(    
    int num_rows, 
    magma_tally3_s_matrix d, 
    magma_tally3_s_matrix b, 
    magma_tally3_s_matrix *c,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_sjacobiupdate(
    magma_tally3_s_matrix t, 
    magma_tally3_s_matrix b, 
    magma_tally3_s_matrix d, 
    magma_tally3_s_matrix *x,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_sjacobispmvupdate(
    magma_tally3_int_t maxiter,
    magma_tally3_s_matrix A, 
    magma_tally3_s_matrix t, 
    magma_tally3_s_matrix b, 
    magma_tally3_s_matrix d, 
    magma_tally3_s_matrix *x,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_sjacobispmvupdate_bw(
    magma_tally3_int_t maxiter,
    magma_tally3_s_matrix A, 
    magma_tally3_s_matrix t, 
    magma_tally3_s_matrix b, 
    magma_tally3_s_matrix d, 
    magma_tally3_s_matrix *x,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_sjacobispmvupdateselect(
    magma_tally3_int_t maxiter,
    magma_tally3_int_t num_updates,
    magma_tally3_index_t *indices,
    magma_tally3_s_matrix A,
    magma_tally3_s_matrix t, 
    magma_tally3_s_matrix b, 
    magma_tally3_s_matrix d, 
    magma_tally3_s_matrix tmp, 
    magma_tally3_s_matrix *x,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_sjacobisetup_diagscal(
    magma_tally3_s_matrix A, magma_tally3_s_matrix *d,
    magma_tally3_queue_t queue );


magma_tally3_int_t
magma_tally3_sbicgmerge1(    
    int n, 
    magma_tally3Float_ptr dskp,
    magma_tally3Float_ptr dv, 
    magma_tally3Float_ptr dr, 
    magma_tally3Float_ptr dp,
    magma_tally3_queue_t queue );


magma_tally3_int_t
magma_tally3_sbicgmerge2(
    int n, 
    magma_tally3Float_ptr dskp, 
    magma_tally3Float_ptr dr,
    magma_tally3Float_ptr dv, 
    magma_tally3Float_ptr ds,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_sbicgmerge3(
    int n, 
    magma_tally3Float_ptr dskp, 
    magma_tally3Float_ptr dp,
    magma_tally3Float_ptr ds,
    magma_tally3Float_ptr dt,
    magma_tally3Float_ptr dx, 
    magma_tally3Float_ptr dr,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_sbicgmerge4(
    int type, 
    magma_tally3Float_ptr dskp,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_scgmerge_spmv1( 
    magma_tally3_s_matrix A,
    magma_tally3Float_ptr d1,
    magma_tally3Float_ptr d2,
    magma_tally3Float_ptr dd,
    magma_tally3Float_ptr dz,
    magma_tally3Float_ptr dskp,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_scgmerge_xrbeta( 
    int n,
    magma_tally3Float_ptr d1,
    magma_tally3Float_ptr d2,
    magma_tally3Float_ptr dx,
    magma_tally3Float_ptr dr,
    magma_tally3Float_ptr dd,
    magma_tally3Float_ptr dz, 
    magma_tally3Float_ptr dskp,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_smdotc(
    magma_tally3_int_t n, 
    magma_tally3_int_t k, 
    magma_tally3Float_ptr dv, 
    magma_tally3Float_ptr dr,
    magma_tally3Float_ptr dd1,
    magma_tally3Float_ptr dd2,
    magma_tally3Float_ptr dskp,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_sgemvmdot(
    int n, 
    int k, 
    magma_tally3Float_ptr dv, 
    magma_tally3Float_ptr dr,
    magma_tally3Float_ptr dd1,
    magma_tally3Float_ptr dd2,
    magma_tally3Float_ptr dskp,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_sbicgmerge_spmv1( 
    magma_tally3_s_matrix A,
    magma_tally3Float_ptr dd1,
    magma_tally3Float_ptr dd2,
    magma_tally3Float_ptr dp,
    magma_tally3Float_ptr dr,
    magma_tally3Float_ptr dv,
    magma_tally3Float_ptr dskp,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_sbicgmerge_spmv2( 
    magma_tally3_s_matrix A,
    magma_tally3Float_ptr dd1,
    magma_tally3Float_ptr dd2,
    magma_tally3Float_ptr ds,
    magma_tally3Float_ptr dt,
    magma_tally3Float_ptr dskp,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_sbicgmerge_xrbeta( 
    int n,
    magma_tally3Float_ptr dd1,
    magma_tally3Float_ptr dd2,
    magma_tally3Float_ptr drr,
    magma_tally3Float_ptr dr,
    magma_tally3Float_ptr dp,
    magma_tally3Float_ptr ds,
    magma_tally3Float_ptr dt,
    magma_tally3Float_ptr dx, 
    magma_tally3Float_ptr dskp,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_sbcsrswp(
    magma_tally3_int_t n,
    magma_tally3_int_t size_b, 
    magma_tally3_int_t *ipiv,
    magma_tally3Float_ptr dx,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_sbcsrtrsv(
    magma_tally3_uplo_t uplo,
    magma_tally3_int_t r_blocks,
    magma_tally3_int_t c_blocks,
    magma_tally3_int_t size_b, 
    magma_tally3Float_ptr dA,
    magma_tally3_index_t *blockinfo, 
    magma_tally3Float_ptr dx,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_sbcsrvalcpy(
    magma_tally3_int_t size_b, 
    magma_tally3_int_t num_blocks, 
    magma_tally3_int_t num_zero_blocks, 
    magma_tally3Float_ptr *dAval, 
    magma_tally3Float_ptr *dBval,
    magma_tally3Float_ptr *dBval2,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_sbcsrluegemm(
    magma_tally3_int_t size_b, 
    magma_tally3_int_t num_block_rows,
    magma_tally3_int_t kblocks,
    magma_tally3Float_ptr *dA, 
    magma_tally3Float_ptr *dB, 
    magma_tally3Float_ptr *dC,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_sbcsrlupivloc(
    magma_tally3_int_t size_b, 
    magma_tally3_int_t kblocks,
    magma_tally3Float_ptr *dA, 
    magma_tally3_int_t *ipiv,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_sbcsrblockinfo5(
    magma_tally3_int_t lustep,
    magma_tally3_int_t num_blocks, 
    magma_tally3_int_t c_blocks, 
    magma_tally3_int_t size_b,
    magma_tally3_index_t *blockinfo,
    magma_tally3Float_ptr dval,
    magma_tally3Float_ptr *AII,
    magma_tally3_queue_t queue );

 
#ifdef __cplusplus
}
#endif

#undef PRECISION_s
#endif /* MAGMA_tally3SPARSE_S_H */
