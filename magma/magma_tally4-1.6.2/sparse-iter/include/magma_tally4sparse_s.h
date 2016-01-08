/*
 -- MAGMA_tally4 (version 1.1) --
 Univ. of Tennessee, Knoxville
 Univ. of California, Berkeley
 Univ. of Colorado, Denver
 @date

 @generated from magma_tally4sparse_z.h normal z -> s, Mon May  4 17:11:25 2015
 @author Hartwig Anzt
*/

#ifndef MAGMA_tally4SPARSE_S_H
#define MAGMA_tally4SPARSE_S_H

#include "magma_tally4_types.h"
#include "magma_tally4sparse_types.h"

#define PRECISION_s


#ifdef __cplusplus
extern "C" {
#endif


/* ////////////////////////////////////////////////////////////////////////////
 -- For backwards compatability, map old (1.6.1) to new (1.6.2) function names
*/

#define magma_tally4_s_mtranspose  magma_tally4_smtranspose
#define magma_tally4_s_mtransfer   magma_tally4_smtransfer
#define magma_tally4_s_vtransfer   magma_tally4_smtransfer
#define magma_tally4_s_mconvert    magma_tally4_smconvert
#define magma_tally4_s_vinit       magma_tally4_svinit
#define magma_tally4_s_vvisu       magma_tally4_sprint_vector
#define magma_tally4_s_vread       magma_tally4_svread
#define magma_tally4_s_vspread     magma_tally4_svspread
#define magma_tally4_s_mvisu       magma_tally4_sprint_matrix
#define magma_tally4_s_mfree       magma_tally4_smfree
#define magma_tally4_s_vfree       magma_tally4_smfree
#define write_s_csr_mtx     magma_tally4_swrite_csr_mtx
#define write_s_csrtomtx    magma_tally4_swrite_csrtomtx
#define print_s_csr         magma_tally4_sprint_csr_mtx


/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally4_SPARSE Auxiliary functions
*/


magma_tally4_int_t
magma_tally4_sparse_opts( 
    int argc, 
    char** argv, 
    magma_tally4_sopts *opts, 
    int *matrices, 
    magma_tally4_queue_t queue );

magma_tally4_int_t 
read_s_csr_from_binary( 
    magma_tally4_int_t* n_row, 
    magma_tally4_int_t* n_col, 
    magma_tally4_int_t* nnz, 
    float **val, 
    magma_tally4_index_t **row, 
    magma_tally4_index_t **col,
    const char * filename,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
read_s_csr_from_mtx( 
    magma_tally4_storage_t *type, 
    magma_tally4_location_t *location,
    magma_tally4_int_t* n_row, 
    magma_tally4_int_t* n_col, 
    magma_tally4_int_t* nnz, 
    float **val, 
    magma_tally4_index_t **row, 
    magma_tally4_index_t **col, 
    const char *filename,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_s_csr_mtx( 
    magma_tally4_s_matrix *A, 
    const char *filename,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_scsrset( 
    magma_tally4_int_t m, 
    magma_tally4_int_t n, 
    magma_tally4_index_t *row, 
    magma_tally4_index_t *col, 
    float *val,
    magma_tally4_s_matrix *A,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_scsrget( 
    magma_tally4_s_matrix A,
    magma_tally4_int_t *m, 
    magma_tally4_int_t *n, 
    magma_tally4_index_t **row, 
    magma_tally4_index_t **col, 
    float **val,
    magma_tally4_queue_t queue );


magma_tally4_int_t 
magma_tally4_svset( 
    magma_tally4_int_t m, 
    magma_tally4_int_t n, 
    float *val,
    magma_tally4_s_matrix *v,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_svget( 
    magma_tally4_s_matrix v,
    magma_tally4_int_t *m, 
    magma_tally4_int_t *n, 
    float **val,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_svset_dev( 
    magma_tally4_int_t m, 
    magma_tally4_int_t n, 
    magma_tally4Float_ptr val,
    magma_tally4_s_matrix *v,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_svget_dev( 
    magma_tally4_s_matrix v,
    magma_tally4_int_t *m, 
    magma_tally4_int_t *n, 
    magma_tally4Float_ptr *val,
    magma_tally4_queue_t queue );


magma_tally4_int_t 
magma_tally4_s_csr_mtxsymm( 
    magma_tally4_s_matrix *A, 
    const char *filename,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_s_csr_compressor( 
    float ** val, 
    magma_tally4_index_t ** row, 
    magma_tally4_index_t ** col, 
    float ** valn, 
    magma_tally4_index_t ** rown, 
    magma_tally4_index_t ** coln, 
    magma_tally4_int_t *n,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_smcsrcompressor( 
    magma_tally4_s_matrix *A,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_smcsrcompressor_gpu( 
    magma_tally4_s_matrix *A,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_svtranspose( 
    magma_tally4_s_matrix x,
    magma_tally4_s_matrix *y,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_s_cucsrtranspose( 
    magma_tally4_s_matrix A, 
    magma_tally4_s_matrix *B,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
s_transpose_csr( 
    magma_tally4_int_t n_rows, 
    magma_tally4_int_t n_cols, 
    magma_tally4_int_t nnz,
    float *val, 
    magma_tally4_index_t *row, 
    magma_tally4_index_t *col, 
    magma_tally4_int_t *new_n_rows, 
    magma_tally4_int_t *new_n_cols, 
    magma_tally4_int_t *new_nnz, 
    float **new_val, 
    magma_tally4_index_t **new_row, 
    magma_tally4_index_t **new_col,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_scsrsplit( 
    magma_tally4_int_t bsize,
    magma_tally4_s_matrix A,
    magma_tally4_s_matrix *D,
    magma_tally4_s_matrix *R,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_smscale( 
    magma_tally4_s_matrix *A, 
    magma_tally4_scale_t scaling,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_smdiff( 
    magma_tally4_s_matrix A, 
    magma_tally4_s_matrix B, 
 real_Double_t *res,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_smdiagadd( 
    magma_tally4_s_matrix *A, 
    float add,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_smsort( 
    magma_tally4_s_matrix *A,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_sindexsort(
    magma_tally4_index_t *x, 
    magma_tally4_int_t first,
    magma_tally4_int_t last,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_sdomainoverlap(
    magma_tally4_index_t num_rows,
    magma_tally4_index_t *num_indices,
    magma_tally4_index_t *rowptr,
    magma_tally4_index_t *colidx,
    magma_tally4_index_t *x,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_ssymbilu( 
    magma_tally4_s_matrix *A, 
    magma_tally4_int_t levels,
    magma_tally4_s_matrix *L,
    magma_tally4_s_matrix *U,
    magma_tally4_queue_t queue );


magma_tally4_int_t 
magma_tally4_swrite_csr_mtx( 
    magma_tally4_s_matrix A,
    magma_tally4_order_t MajorType,
 const char *filename,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_swrite_csrtomtx( 
    magma_tally4_s_matrix A,
    const char *filename,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_sprint_csr( 
    magma_tally4_int_t n_row, 
    magma_tally4_int_t n_col, 
    magma_tally4_int_t nnz, 
    float **val, 
    magma_tally4_index_t **row, 
    magma_tally4_index_t **col,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_sprint_csr_mtx( 
    magma_tally4_int_t n_row, 
    magma_tally4_int_t n_col, 
    magma_tally4_int_t nnz, 
    float **val, 
    magma_tally4_index_t **row, 
    magma_tally4_index_t **col, 
    magma_tally4_order_t MajorType,
    magma_tally4_queue_t queue );


magma_tally4_int_t 
magma_tally4_smtranspose(
    magma_tally4_s_matrix A, 
    magma_tally4_s_matrix *B,
    magma_tally4_queue_t queue );


magma_tally4_int_t 
magma_tally4_smtransfer(
    magma_tally4_s_matrix A, 
    magma_tally4_s_matrix *B, 
    magma_tally4_location_t src, 
    magma_tally4_location_t dst,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_smconvert(
    magma_tally4_s_matrix A, 
    magma_tally4_s_matrix *B, 
    magma_tally4_storage_t old_format, 
    magma_tally4_storage_t new_format,
    magma_tally4_queue_t queue );


magma_tally4_int_t
magma_tally4_svinit(
    magma_tally4_s_matrix *x, 
    magma_tally4_location_t memory_location,
    magma_tally4_int_t num_rows, 
    magma_tally4_int_t num_cols,
    float values,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_sprint_vector(
    magma_tally4_s_matrix x, 
    magma_tally4_int_t offset, 
    magma_tally4_int_t displaylength,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_svread(
    magma_tally4_s_matrix *x, 
    magma_tally4_int_t length,
    char * filename,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_svspread(
    magma_tally4_s_matrix *x, 
    const char * filename,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_sprint_matrix(
    magma_tally4_s_matrix A,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_sdiameter(
    magma_tally4_s_matrix *A,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_srowentries(
    magma_tally4_s_matrix *A,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_smfree(
    magma_tally4_s_matrix *A,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_sresidual(
    magma_tally4_s_matrix A, 
    magma_tally4_s_matrix b, 
    magma_tally4_s_matrix x, 
    float *res,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_sresidualvec(
    magma_tally4_s_matrix A,
    magma_tally4_s_matrix b,
    magma_tally4_s_matrix x,
    magma_tally4_s_matrix *r,
    float *res,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_smgenerator(
    magma_tally4_int_t n,
    magma_tally4_int_t offdiags,
    magma_tally4_index_t *diag_offset,
    float *diag_vals,
    magma_tally4_s_matrix *A,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_sm_27stencil(
    magma_tally4_int_t n,
    magma_tally4_s_matrix *A,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_sm_5stencil(
    magma_tally4_int_t n,
    magma_tally4_s_matrix *A,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_ssolverinfo(
    magma_tally4_s_solver_par *solver_par, 
    magma_tally4_s_preconditioner *precond_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_ssolverinfo_init(
    magma_tally4_s_solver_par *solver_par, 
    magma_tally4_s_preconditioner *precond,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_seigensolverinfo_init(
    magma_tally4_s_solver_par *solver_par,
    magma_tally4_queue_t queue );


magma_tally4_int_t
magma_tally4_ssolverinfo_free(
    magma_tally4_s_solver_par *solver_par, 
    magma_tally4_s_preconditioner *precond,
    magma_tally4_queue_t queue );



/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally4_SPARSE iterative incomplete factorizations
*/


magma_tally4_int_t
magma_tally4_siterilusetup( 
    magma_tally4_s_matrix A, 
    magma_tally4_s_matrix b,                                 
    magma_tally4_s_preconditioner *precond,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_sitericsetup( 
    magma_tally4_s_matrix A, 
    magma_tally4_s_matrix b, 
    magma_tally4_s_preconditioner *precond,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_sitericupdate( 
    magma_tally4_s_matrix A, 
    magma_tally4_s_preconditioner *precond, 
    magma_tally4_int_t updates,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_sapplyiteric_l( 
    magma_tally4_s_matrix b, 
    magma_tally4_s_matrix *x, 
    magma_tally4_s_preconditioner *precond,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_sapplyiteric_r( 
    magma_tally4_s_matrix b, 
    magma_tally4_s_matrix *x, 
    magma_tally4_s_preconditioner *precond,
    magma_tally4_queue_t queue );


magma_tally4_int_t
magma_tally4_siterilu_csr( 
    magma_tally4_s_matrix A,
    magma_tally4_s_matrix L,
    magma_tally4_s_matrix U,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_siteric_csr( 
    magma_tally4_s_matrix A,
    magma_tally4_s_matrix A_CSR,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_sfrobenius( 
    magma_tally4_s_matrix A, 
    magma_tally4_s_matrix B, 
    real_Double_t *res,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_snonlinres(   
    magma_tally4_s_matrix A, 
    magma_tally4_s_matrix L,
    magma_tally4_s_matrix U, 
    magma_tally4_s_matrix *LU, 
    real_Double_t *res,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_silures(   
    magma_tally4_s_matrix A, 
    magma_tally4_s_matrix L,
    magma_tally4_s_matrix U, 
    magma_tally4_s_matrix *LU, 
    real_Double_t *res,
    real_Double_t *nonlinres,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_sicres(       
    magma_tally4_s_matrix A, 
    magma_tally4_s_matrix C,
    magma_tally4_s_matrix CT, 
    magma_tally4_s_matrix *LU, 
    real_Double_t *res,
    real_Double_t *nonlinres,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_sinitguess( 
    magma_tally4_s_matrix A, 
    magma_tally4_s_matrix *L, 
    magma_tally4_s_matrix *U,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_sinitrecursiveLU( 
    magma_tally4_s_matrix A, 
    magma_tally4_s_matrix *B,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_smLdiagadd( 
    magma_tally4_s_matrix *L,
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
magma_tally4_scg(
    magma_tally4_s_matrix A, magma_tally4_s_matrix b, 
    magma_tally4_s_matrix *x, magma_tally4_s_solver_par *solver_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_scg_res(
    magma_tally4_s_matrix A, magma_tally4_s_matrix b, 
    magma_tally4_s_matrix *x, magma_tally4_s_solver_par *solver_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_scg_merge(
    magma_tally4_s_matrix A, magma_tally4_s_matrix b, 
    magma_tally4_s_matrix *x, magma_tally4_s_solver_par *solver_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_sgmres(
    magma_tally4_s_matrix A, magma_tally4_s_matrix b, 
    magma_tally4_s_matrix *x, magma_tally4_s_solver_par *solver_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_sbicgstab(
    magma_tally4_s_matrix A, magma_tally4_s_matrix b, magma_tally4_s_matrix *x, 
    magma_tally4_s_solver_par *solver_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_sbicgstab_merge(
    magma_tally4_s_matrix A, magma_tally4_s_matrix b, 
    magma_tally4_s_matrix *x, magma_tally4_s_solver_par *solver_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_sbicgstab_merge2(
    magma_tally4_s_matrix A, magma_tally4_s_matrix b, 
    magma_tally4_s_matrix *x, magma_tally4_s_solver_par *solver_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_spcg(
    magma_tally4_s_matrix A, magma_tally4_s_matrix b, 
    magma_tally4_s_matrix *x, magma_tally4_s_solver_par *solver_par, 
    magma_tally4_s_preconditioner *precond_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_sbpcg(
    magma_tally4_s_matrix A, magma_tally4_s_matrix b, 
    magma_tally4_s_matrix *x, magma_tally4_s_solver_par *solver_par, 
    magma_tally4_s_preconditioner *precond_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_spbicgstab(
    magma_tally4_s_matrix A, magma_tally4_s_matrix b, 
    magma_tally4_s_matrix *x, magma_tally4_s_solver_par *solver_par, 
    magma_tally4_s_preconditioner *precond_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_spgmres(
    magma_tally4_s_matrix A, magma_tally4_s_matrix b, 
    magma_tally4_s_matrix *x, magma_tally4_s_solver_par *solver_par, 
    magma_tally4_s_preconditioner *precond_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_sfgmres(
    magma_tally4_s_matrix A, magma_tally4_s_matrix b, 
    magma_tally4_s_matrix *x, magma_tally4_s_solver_par *solver_par, 
    magma_tally4_s_preconditioner *precond_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_sjacobi(
    magma_tally4_s_matrix A, magma_tally4_s_matrix b, 
    magma_tally4_s_matrix *x, magma_tally4_s_solver_par *solver_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_sjacobidomainoverlap(
    magma_tally4_s_matrix A, 
    magma_tally4_s_matrix b, 
    magma_tally4_s_matrix *x,  
    magma_tally4_s_solver_par *solver_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_sbaiter(
    magma_tally4_s_matrix A, magma_tally4_s_matrix b, 
    magma_tally4_s_matrix *x, magma_tally4_s_solver_par *solver_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_siterref(
    magma_tally4_s_matrix A, magma_tally4_s_matrix b, 
    magma_tally4_s_matrix *x, magma_tally4_s_solver_par *solver_par, 
    magma_tally4_s_preconditioner *precond_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_silu(
    magma_tally4_s_matrix A, magma_tally4_s_matrix b, 
    magma_tally4_s_matrix *x, magma_tally4_s_solver_par *solver_par, 
    magma_tally4_int_t *ipiv,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_sbcsrlu(
    magma_tally4_s_matrix A, magma_tally4_s_matrix b, 
    magma_tally4_s_matrix *x, magma_tally4_s_solver_par *solver_par,
    magma_tally4_queue_t queue );


magma_tally4_int_t
magma_tally4_sbcsrlutrf(
    magma_tally4_s_matrix A, 
    magma_tally4_s_matrix *M,
    magma_tally4_int_t *ipiv, 
    magma_tally4_int_t version,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_sbcsrlusv(
    magma_tally4_s_matrix A, magma_tally4_s_matrix b, 
    magma_tally4_s_matrix *x, magma_tally4_s_solver_par *solver_par, 
    magma_tally4_int_t *ipiv,
    magma_tally4_queue_t queue );



magma_tally4_int_t
magma_tally4_silucg(
    magma_tally4_s_matrix A, magma_tally4_s_matrix b, 
    magma_tally4_s_matrix *x, magma_tally4_s_solver_par *solver_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_silugmres(
    magma_tally4_s_matrix A, magma_tally4_s_matrix b, 
    magma_tally4_s_matrix *x, magma_tally4_s_solver_par *solver_par,
    magma_tally4_queue_t queue ); 


magma_tally4_int_t
magma_tally4_slobpcg_shift(
    magma_tally4_int_t num_rows,
    magma_tally4_int_t num_vecs, 
    magma_tally4_int_t shift,
    magma_tally4Float_ptr x,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_slobpcg_res(
    magma_tally4_int_t num_rows,
    magma_tally4_int_t num_vecs, 
    float *evalues, 
    magma_tally4Float_ptr X,
    magma_tally4Float_ptr R, 
    float *res,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_slobpcg_maxpy(
    magma_tally4_int_t num_rows,
    magma_tally4_int_t num_vecs, 
    magma_tally4Float_ptr X,
    magma_tally4Float_ptr Y,
    magma_tally4_queue_t queue );


/*/////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally4_SPARSE eigensolvers (Data on GPU)
*/
magma_tally4_int_t
magma_tally4_slobpcg(
    magma_tally4_s_matrix A, 
    magma_tally4_s_solver_par *solver_par,
    magma_tally4_s_preconditioner *precond_par, 
    magma_tally4_queue_t queue );




/*/////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally4_SPARSE preconditioners (Data on GPU)
*/

magma_tally4_int_t
magma_tally4_sjacobisetup(
    magma_tally4_s_matrix A, 
    magma_tally4_s_matrix b, 
    magma_tally4_s_matrix *M, 
    magma_tally4_s_matrix *c,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_sjacobisetup_matrix(
    magma_tally4_s_matrix A, 
    magma_tally4_s_matrix *M, 
    magma_tally4_s_matrix *d,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_sjacobisetup_vector(
    magma_tally4_s_matrix b, 
    magma_tally4_s_matrix d, 
    magma_tally4_s_matrix *c,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_sjacobiiter(
    magma_tally4_s_matrix M, 
    magma_tally4_s_matrix c, 
    magma_tally4_s_matrix *x, 
    magma_tally4_s_solver_par *solver_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_sjacobiiter_precond( 
    magma_tally4_s_matrix M, 
    magma_tally4_s_matrix *x, 
    magma_tally4_s_solver_par *solver_par, 
    magma_tally4_s_preconditioner *precond,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_sjacobiiter_sys(
    magma_tally4_s_matrix A, 
    magma_tally4_s_matrix b, 
    magma_tally4_s_matrix d, 
    magma_tally4_s_matrix t, 
    magma_tally4_s_matrix *x,  
    magma_tally4_s_solver_par *solver_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_spastixsetup(
    magma_tally4_s_matrix A, magma_tally4_s_matrix b,
    magma_tally4_s_preconditioner *precond,
    magma_tally4_queue_t queue );


magma_tally4_int_t
magma_tally4_sapplypastix(
    magma_tally4_s_matrix b, magma_tally4_s_matrix *x, 
    magma_tally4_s_preconditioner *precond,
    magma_tally4_queue_t queue );


// custom preconditioner
magma_tally4_int_t
magma_tally4_sapplycustomprecond_l(
    magma_tally4_s_matrix b, 
    magma_tally4_s_matrix *x, 
    magma_tally4_s_preconditioner *precond,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_sapplycustomprecond_r(
    magma_tally4_s_matrix b, 
    magma_tally4_s_matrix *x, 
    magma_tally4_s_preconditioner *precond,
    magma_tally4_queue_t queue );


// CUSPARSE preconditioner

magma_tally4_int_t
magma_tally4_scuilusetup(
    magma_tally4_s_matrix A, magma_tally4_s_preconditioner *precond,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_sapplycuilu_l(
    magma_tally4_s_matrix b, magma_tally4_s_matrix *x, 
    magma_tally4_s_preconditioner *precond,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_sapplycuilu_r(
    magma_tally4_s_matrix b, magma_tally4_s_matrix *x, 
    magma_tally4_s_preconditioner *precond,
    magma_tally4_queue_t queue );


magma_tally4_int_t
magma_tally4_scuiccsetup(
    magma_tally4_s_matrix A, magma_tally4_s_preconditioner *precond,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_sapplycuicc_l(
    magma_tally4_s_matrix b, magma_tally4_s_matrix *x, 
    magma_tally4_s_preconditioner *precond,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_sapplycuicc_r(
    magma_tally4_s_matrix b, magma_tally4_s_matrix *x, 
    magma_tally4_s_preconditioner *precond,
    magma_tally4_queue_t queue );


magma_tally4_int_t
magma_tally4_scumilusetup(
    magma_tally4_s_matrix A, 
    magma_tally4_s_preconditioner *precond,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_scumilugeneratesolverinfo(
    magma_tally4_s_preconditioner *precond,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_sapplycumilu_l(
    magma_tally4_s_matrix b, 
    magma_tally4_s_matrix *x, 
    magma_tally4_s_preconditioner *precond,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_sapplycumilu_r(
    magma_tally4_s_matrix b, 
    magma_tally4_s_matrix *x, 
    magma_tally4_s_preconditioner *precond,
    magma_tally4_queue_t queue );


magma_tally4_int_t
magma_tally4_scumiccsetup(
    magma_tally4_s_matrix A, 
    magma_tally4_s_preconditioner *precond,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_scumicgeneratesolverinfo(
    magma_tally4_s_preconditioner *precond,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_sapplycumicc_l(
    magma_tally4_s_matrix b, 
    magma_tally4_s_matrix *x, 
    magma_tally4_s_preconditioner *precond,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_sapplycumicc_r(
    magma_tally4_s_matrix b, 
    magma_tally4_s_matrix *x, 
    magma_tally4_s_preconditioner *precond,
    magma_tally4_queue_t queue );


// block-asynchronous iteration

magma_tally4_int_t
magma_tally4_sbajac_csr(
    magma_tally4_int_t localiters,
    magma_tally4_s_matrix D,
    magma_tally4_s_matrix R,
    magma_tally4_s_matrix b,
    magma_tally4_s_matrix *x,
    magma_tally4_queue_t queue );

/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally4_SPARSE utility function definitions
*/

magma_tally4_int_t
magma_tally4_s_spmv(
    float alpha, 
    magma_tally4_s_matrix A, 
    magma_tally4_s_matrix x, 
    float beta, 
    magma_tally4_s_matrix y,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_scustomspmv(
    float alpha, 
    magma_tally4_s_matrix x, 
    float beta, 
    magma_tally4_s_matrix y,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_s_spmv_shift(
    float alpha, 
    magma_tally4_s_matrix A, 
    float lambda,
    magma_tally4_s_matrix x, 
    float beta, 
    magma_tally4_int_t offset, 
    magma_tally4_int_t blocksize,
    magma_tally4Index_ptr dadd_vecs, 
    magma_tally4_s_matrix y,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_scuspmm(
    magma_tally4_s_matrix A, 
    magma_tally4_s_matrix B, 
    magma_tally4_s_matrix *AB,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_s_spmm(
    float alpha, 
    magma_tally4_s_matrix A,
    magma_tally4_s_matrix B,
    magma_tally4_s_matrix *C,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_ssymbilu( 
    magma_tally4_s_matrix *A, 
    magma_tally4_int_t levels, 
    magma_tally4_s_matrix *L, 
    magma_tally4_s_matrix *U,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_scuspaxpy(
    magma_tally4Float_ptr alpha, magma_tally4_s_matrix A, 
    magma_tally4Float_ptr beta, magma_tally4_s_matrix B, 
    magma_tally4_s_matrix *AB,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_s_precond(
    magma_tally4_s_matrix A, 
    magma_tally4_s_matrix b, magma_tally4_s_matrix *x,
    magma_tally4_s_preconditioner *precond,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_s_solver(
    magma_tally4_s_matrix A, magma_tally4_s_matrix b, 
    magma_tally4_s_matrix *x, magma_tally4_sopts *zopts,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_s_precondsetup(
    magma_tally4_s_matrix A, magma_tally4_s_matrix b, 
    magma_tally4_s_preconditioner *precond,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_s_applyprecond(
    magma_tally4_s_matrix A, magma_tally4_s_matrix b, 
    magma_tally4_s_matrix *x, magma_tally4_s_preconditioner *precond,
    magma_tally4_queue_t queue );


magma_tally4_int_t
magma_tally4_s_applyprecond_left(
    magma_tally4_s_matrix A, magma_tally4_s_matrix b, 
    magma_tally4_s_matrix *x, magma_tally4_s_preconditioner *precond,
    magma_tally4_queue_t queue );


magma_tally4_int_t
magma_tally4_s_applyprecond_right(
    magma_tally4_s_matrix A, magma_tally4_s_matrix b, 
    magma_tally4_s_matrix *x, magma_tally4_s_preconditioner *precond,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_s_initP2P(
    magma_tally4_int_t *bandwidth_benchmark,
    magma_tally4_int_t *num_gpus,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_scompact(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4Float_ptr dA, magma_tally4_int_t ldda,
    float *dnorms, float tol, 
    magma_tally4_int_t *activeMask, magma_tally4_int_t *cBlockSize,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_scompactActive(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4Float_ptr dA, magma_tally4_int_t ldda, 
    magma_tally4_int_t *active,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_smlumerge(    
    magma_tally4_s_matrix L, 
    magma_tally4_s_matrix U,
    magma_tally4_s_matrix *A, 
    magma_tally4_queue_t queue );


/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally4_SPARSE BLAS function definitions
*/
magma_tally4_int_t 
magma_tally4_sgecsrmv(
    magma_tally4_trans_t transA,
    magma_tally4_int_t m, magma_tally4_int_t n,
    float alpha,
    magma_tally4Float_ptr dval,
    magma_tally4Index_ptr drowptr,
    magma_tally4Index_ptr dcolind,
    magma_tally4Float_ptr dx,
    float beta,
    magma_tally4Float_ptr dy,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_sgecsrmv_shift(
    magma_tally4_trans_t transA,
    magma_tally4_int_t m, magma_tally4_int_t n,
    float alpha,
    float lambda,
    magma_tally4Float_ptr dval,
    magma_tally4Index_ptr drowptr,
    magma_tally4Index_ptr dcolind,
    magma_tally4Float_ptr dx,
    float beta,
    int offset,
    int blocksize,
    magma_tally4Index_ptr dadd_rows,
    magma_tally4Float_ptr dy,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_smgecsrmv(
    magma_tally4_trans_t transA,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4_int_t num_vecs,
    float alpha,
    magma_tally4Float_ptr dval,
    magma_tally4Index_ptr drowptr,
    magma_tally4Index_ptr dcolind,
    magma_tally4Float_ptr dx,
    float beta,
    magma_tally4Float_ptr dy,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_sgeellmv(
    magma_tally4_trans_t transA,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4_int_t nnz_per_row,
    float alpha,
    magma_tally4Float_ptr dval,
    magma_tally4Index_ptr dcolind,
    magma_tally4Float_ptr dx,
    float beta,
    magma_tally4Float_ptr dy,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_sgeellmv_shift(
    magma_tally4_trans_t transA,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4_int_t nnz_per_row,
    float alpha,
    float lambda,
    magma_tally4Float_ptr dval,
    magma_tally4Index_ptr dcolind,
    magma_tally4Float_ptr dx,
    float beta,
    int offset,
    int blocksize,
    magma_tally4Index_ptr dadd_rows,
    magma_tally4Float_ptr dy,
    magma_tally4_queue_t queue );


magma_tally4_int_t 
magma_tally4_smgeellmv(
    magma_tally4_trans_t transA,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4_int_t num_vecs,
    magma_tally4_int_t nnz_per_row,
    float alpha,
    magma_tally4Float_ptr dval,
    magma_tally4Index_ptr dcolind,
    magma_tally4Float_ptr dx,
    float beta,
    magma_tally4Float_ptr dy,
    magma_tally4_queue_t queue );


magma_tally4_int_t 
magma_tally4_sgeelltmv(
    magma_tally4_trans_t transA,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4_int_t nnz_per_row,
    float alpha,
    magma_tally4Float_ptr dval,
    magma_tally4Index_ptr dcolind,
    magma_tally4Float_ptr dx,
    float beta,
    magma_tally4Float_ptr dy,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_sgeelltmv_shift(
    magma_tally4_trans_t transA,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4_int_t nnz_per_row,
    float alpha,
    float lambda,
    magma_tally4Float_ptr dval,
    magma_tally4Index_ptr dcolind,
    magma_tally4Float_ptr dx,
    float beta,
    int offset,
    int blocksize,
    magma_tally4Index_ptr dadd_rows,
    magma_tally4Float_ptr dy,
    magma_tally4_queue_t queue );


magma_tally4_int_t 
magma_tally4_smgeelltmv(
    magma_tally4_trans_t transA,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4_int_t num_vecs,
    magma_tally4_int_t nnz_per_row,
    float alpha,
    magma_tally4Float_ptr dval,
    magma_tally4Index_ptr dcolind,
    magma_tally4Float_ptr dx,
    float beta,
    magma_tally4Float_ptr dy,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_sgeellrtmv(
    magma_tally4_trans_t transA,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4_int_t nnz_per_row,
    float alpha,
    magma_tally4Float_ptr dval,
    magma_tally4Index_ptr dcolind,
    magma_tally4Index_ptr drowlength,
    magma_tally4Float_ptr dx,
    float beta,
    magma_tally4Float_ptr dy,
    magma_tally4_int_t num_threads,
    magma_tally4_int_t threads_per_row,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_sgesellcmv(
    magma_tally4_trans_t transA,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4_int_t blocksize,
    magma_tally4_int_t slices,
    magma_tally4_int_t alignment,
    float alpha,
    magma_tally4Float_ptr dval,
    magma_tally4Index_ptr dcolind,
    magma_tally4Index_ptr drowptr,
    magma_tally4Float_ptr dx,
    float beta,
    magma_tally4Float_ptr dy,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_sgesellpmv(
    magma_tally4_trans_t transA,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4_int_t blocksize,
    magma_tally4_int_t slices,
    magma_tally4_int_t alignment,
    float alpha,
    magma_tally4Float_ptr dval,
    magma_tally4Index_ptr dcolind,
    magma_tally4Index_ptr drowptr,
    magma_tally4Float_ptr dx,
    float beta,
    magma_tally4Float_ptr dy,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_smgesellpmv(
    magma_tally4_trans_t transA,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4_int_t num_vecs,
    magma_tally4_int_t blocksize,
    magma_tally4_int_t slices,
    magma_tally4_int_t alignment,
    float alpha,
    magma_tally4Float_ptr dval,
    magma_tally4Index_ptr dcolind,
    magma_tally4Index_ptr drowptr,
    magma_tally4Float_ptr dx,
    float beta,
    magma_tally4Float_ptr dy,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_smgesellpmv_blocked(
    magma_tally4_trans_t transA,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4_int_t num_vecs,
    magma_tally4_int_t blocksize,
    magma_tally4_int_t slices,
    magma_tally4_int_t alignment,
    float alpha,
    magma_tally4Float_ptr dval,
    magma_tally4Index_ptr dcolind,
    magma_tally4Index_ptr drowptr,
    magma_tally4Float_ptr dx,
    float beta,
    magma_tally4Float_ptr dy,
    magma_tally4_queue_t queue );


magma_tally4_int_t
magma_tally4_smergedgs(
    magma_tally4_int_t n, 
    magma_tally4_int_t ldh,
    magma_tally4_int_t k, 
    magma_tally4Float_ptr dv, 
    magma_tally4Float_ptr dr,
    magma_tally4Float_ptr dskp,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_scopyscale(    
    int n, 
    int k,
    magma_tally4Float_ptr dr, 
    magma_tally4Float_ptr dv,
    magma_tally4Float_ptr dskp,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_snrm2scale(    
    int m, 
    magma_tally4Float_ptr dr,    
    int lddr, 
    float *drnorm,
    magma_tally4_queue_t queue );


magma_tally4_int_t
magma_tally4_sjacobisetup_vector_gpu(
    int num_rows, 
    magma_tally4_s_matrix b, 
    magma_tally4_s_matrix d, 
    magma_tally4_s_matrix c,
    magma_tally4_s_matrix *x,
    magma_tally4_queue_t queue );


magma_tally4_int_t
magma_tally4_sjacobi_diagscal(    
    int num_rows, 
    magma_tally4_s_matrix d, 
    magma_tally4_s_matrix b, 
    magma_tally4_s_matrix *c,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_sjacobiupdate(
    magma_tally4_s_matrix t, 
    magma_tally4_s_matrix b, 
    magma_tally4_s_matrix d, 
    magma_tally4_s_matrix *x,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_sjacobispmvupdate(
    magma_tally4_int_t maxiter,
    magma_tally4_s_matrix A, 
    magma_tally4_s_matrix t, 
    magma_tally4_s_matrix b, 
    magma_tally4_s_matrix d, 
    magma_tally4_s_matrix *x,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_sjacobispmvupdate_bw(
    magma_tally4_int_t maxiter,
    magma_tally4_s_matrix A, 
    magma_tally4_s_matrix t, 
    magma_tally4_s_matrix b, 
    magma_tally4_s_matrix d, 
    magma_tally4_s_matrix *x,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_sjacobispmvupdateselect(
    magma_tally4_int_t maxiter,
    magma_tally4_int_t num_updates,
    magma_tally4_index_t *indices,
    magma_tally4_s_matrix A,
    magma_tally4_s_matrix t, 
    magma_tally4_s_matrix b, 
    magma_tally4_s_matrix d, 
    magma_tally4_s_matrix tmp, 
    magma_tally4_s_matrix *x,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_sjacobisetup_diagscal(
    magma_tally4_s_matrix A, magma_tally4_s_matrix *d,
    magma_tally4_queue_t queue );


magma_tally4_int_t
magma_tally4_sbicgmerge1(    
    int n, 
    magma_tally4Float_ptr dskp,
    magma_tally4Float_ptr dv, 
    magma_tally4Float_ptr dr, 
    magma_tally4Float_ptr dp,
    magma_tally4_queue_t queue );


magma_tally4_int_t
magma_tally4_sbicgmerge2(
    int n, 
    magma_tally4Float_ptr dskp, 
    magma_tally4Float_ptr dr,
    magma_tally4Float_ptr dv, 
    magma_tally4Float_ptr ds,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_sbicgmerge3(
    int n, 
    magma_tally4Float_ptr dskp, 
    magma_tally4Float_ptr dp,
    magma_tally4Float_ptr ds,
    magma_tally4Float_ptr dt,
    magma_tally4Float_ptr dx, 
    magma_tally4Float_ptr dr,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_sbicgmerge4(
    int type, 
    magma_tally4Float_ptr dskp,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_scgmerge_spmv1( 
    magma_tally4_s_matrix A,
    magma_tally4Float_ptr d1,
    magma_tally4Float_ptr d2,
    magma_tally4Float_ptr dd,
    magma_tally4Float_ptr dz,
    magma_tally4Float_ptr dskp,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_scgmerge_xrbeta( 
    int n,
    magma_tally4Float_ptr d1,
    magma_tally4Float_ptr d2,
    magma_tally4Float_ptr dx,
    magma_tally4Float_ptr dr,
    magma_tally4Float_ptr dd,
    magma_tally4Float_ptr dz, 
    magma_tally4Float_ptr dskp,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_smdotc(
    magma_tally4_int_t n, 
    magma_tally4_int_t k, 
    magma_tally4Float_ptr dv, 
    magma_tally4Float_ptr dr,
    magma_tally4Float_ptr dd1,
    magma_tally4Float_ptr dd2,
    magma_tally4Float_ptr dskp,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_sgemvmdot(
    int n, 
    int k, 
    magma_tally4Float_ptr dv, 
    magma_tally4Float_ptr dr,
    magma_tally4Float_ptr dd1,
    magma_tally4Float_ptr dd2,
    magma_tally4Float_ptr dskp,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_sbicgmerge_spmv1( 
    magma_tally4_s_matrix A,
    magma_tally4Float_ptr dd1,
    magma_tally4Float_ptr dd2,
    magma_tally4Float_ptr dp,
    magma_tally4Float_ptr dr,
    magma_tally4Float_ptr dv,
    magma_tally4Float_ptr dskp,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_sbicgmerge_spmv2( 
    magma_tally4_s_matrix A,
    magma_tally4Float_ptr dd1,
    magma_tally4Float_ptr dd2,
    magma_tally4Float_ptr ds,
    magma_tally4Float_ptr dt,
    magma_tally4Float_ptr dskp,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_sbicgmerge_xrbeta( 
    int n,
    magma_tally4Float_ptr dd1,
    magma_tally4Float_ptr dd2,
    magma_tally4Float_ptr drr,
    magma_tally4Float_ptr dr,
    magma_tally4Float_ptr dp,
    magma_tally4Float_ptr ds,
    magma_tally4Float_ptr dt,
    magma_tally4Float_ptr dx, 
    magma_tally4Float_ptr dskp,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_sbcsrswp(
    magma_tally4_int_t n,
    magma_tally4_int_t size_b, 
    magma_tally4_int_t *ipiv,
    magma_tally4Float_ptr dx,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_sbcsrtrsv(
    magma_tally4_uplo_t uplo,
    magma_tally4_int_t r_blocks,
    magma_tally4_int_t c_blocks,
    magma_tally4_int_t size_b, 
    magma_tally4Float_ptr dA,
    magma_tally4_index_t *blockinfo, 
    magma_tally4Float_ptr dx,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_sbcsrvalcpy(
    magma_tally4_int_t size_b, 
    magma_tally4_int_t num_blocks, 
    magma_tally4_int_t num_zero_blocks, 
    magma_tally4Float_ptr *dAval, 
    magma_tally4Float_ptr *dBval,
    magma_tally4Float_ptr *dBval2,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_sbcsrluegemm(
    magma_tally4_int_t size_b, 
    magma_tally4_int_t num_block_rows,
    magma_tally4_int_t kblocks,
    magma_tally4Float_ptr *dA, 
    magma_tally4Float_ptr *dB, 
    magma_tally4Float_ptr *dC,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_sbcsrlupivloc(
    magma_tally4_int_t size_b, 
    magma_tally4_int_t kblocks,
    magma_tally4Float_ptr *dA, 
    magma_tally4_int_t *ipiv,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_sbcsrblockinfo5(
    magma_tally4_int_t lustep,
    magma_tally4_int_t num_blocks, 
    magma_tally4_int_t c_blocks, 
    magma_tally4_int_t size_b,
    magma_tally4_index_t *blockinfo,
    magma_tally4Float_ptr dval,
    magma_tally4Float_ptr *AII,
    magma_tally4_queue_t queue );

 
#ifdef __cplusplus
}
#endif

#undef PRECISION_s
#endif /* MAGMA_tally4SPARSE_S_H */
