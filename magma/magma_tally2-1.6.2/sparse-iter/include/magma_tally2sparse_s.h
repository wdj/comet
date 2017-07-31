/*
 -- MAGMA_tally2 (version 1.1) --
 Univ. of Tennessee, Knoxville
 Univ. of California, Berkeley
 Univ. of Colorado, Denver
 @date

 @generated from magma_tally2sparse_z.h normal z -> s, Mon May  4 17:11:25 2015
 @author Hartwig Anzt
*/

#ifndef MAGMA_tally2SPARSE_S_H
#define MAGMA_tally2SPARSE_S_H

#include "magma_tally2_types.h"
#include "magma_tally2sparse_types.h"

#define PRECISION_s


#ifdef __cplusplus
extern "C" {
#endif


/* ////////////////////////////////////////////////////////////////////////////
 -- For backwards compatability, map old (1.6.1) to new (1.6.2) function names
*/

#define magma_tally2_s_mtranspose  magma_tally2_smtranspose
#define magma_tally2_s_mtransfer   magma_tally2_smtransfer
#define magma_tally2_s_vtransfer   magma_tally2_smtransfer
#define magma_tally2_s_mconvert    magma_tally2_smconvert
#define magma_tally2_s_vinit       magma_tally2_svinit
#define magma_tally2_s_vvisu       magma_tally2_sprint_vector
#define magma_tally2_s_vread       magma_tally2_svread
#define magma_tally2_s_vspread     magma_tally2_svspread
#define magma_tally2_s_mvisu       magma_tally2_sprint_matrix
#define magma_tally2_s_mfree       magma_tally2_smfree
#define magma_tally2_s_vfree       magma_tally2_smfree
#define write_s_csr_mtx     magma_tally2_swrite_csr_mtx
#define write_s_csrtomtx    magma_tally2_swrite_csrtomtx
#define print_s_csr         magma_tally2_sprint_csr_mtx


/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally2_SPARSE Auxiliary functions
*/


magma_tally2_int_t
magma_tally2_sparse_opts( 
    int argc, 
    char** argv, 
    magma_tally2_sopts *opts, 
    int *matrices, 
    magma_tally2_queue_t queue );

magma_tally2_int_t 
read_s_csr_from_binary( 
    magma_tally2_int_t* n_row, 
    magma_tally2_int_t* n_col, 
    magma_tally2_int_t* nnz, 
    float **val, 
    magma_tally2_index_t **row, 
    magma_tally2_index_t **col,
    const char * filename,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
read_s_csr_from_mtx( 
    magma_tally2_storage_t *type, 
    magma_tally2_location_t *location,
    magma_tally2_int_t* n_row, 
    magma_tally2_int_t* n_col, 
    magma_tally2_int_t* nnz, 
    float **val, 
    magma_tally2_index_t **row, 
    magma_tally2_index_t **col, 
    const char *filename,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_s_csr_mtx( 
    magma_tally2_s_matrix *A, 
    const char *filename,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_scsrset( 
    magma_tally2_int_t m, 
    magma_tally2_int_t n, 
    magma_tally2_index_t *row, 
    magma_tally2_index_t *col, 
    float *val,
    magma_tally2_s_matrix *A,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_scsrget( 
    magma_tally2_s_matrix A,
    magma_tally2_int_t *m, 
    magma_tally2_int_t *n, 
    magma_tally2_index_t **row, 
    magma_tally2_index_t **col, 
    float **val,
    magma_tally2_queue_t queue );


magma_tally2_int_t 
magma_tally2_svset( 
    magma_tally2_int_t m, 
    magma_tally2_int_t n, 
    float *val,
    magma_tally2_s_matrix *v,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_svget( 
    magma_tally2_s_matrix v,
    magma_tally2_int_t *m, 
    magma_tally2_int_t *n, 
    float **val,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_svset_dev( 
    magma_tally2_int_t m, 
    magma_tally2_int_t n, 
    magma_tally2Float_ptr val,
    magma_tally2_s_matrix *v,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_svget_dev( 
    magma_tally2_s_matrix v,
    magma_tally2_int_t *m, 
    magma_tally2_int_t *n, 
    magma_tally2Float_ptr *val,
    magma_tally2_queue_t queue );


magma_tally2_int_t 
magma_tally2_s_csr_mtxsymm( 
    magma_tally2_s_matrix *A, 
    const char *filename,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_s_csr_compressor( 
    float ** val, 
    magma_tally2_index_t ** row, 
    magma_tally2_index_t ** col, 
    float ** valn, 
    magma_tally2_index_t ** rown, 
    magma_tally2_index_t ** coln, 
    magma_tally2_int_t *n,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_smcsrcompressor( 
    magma_tally2_s_matrix *A,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_smcsrcompressor_gpu( 
    magma_tally2_s_matrix *A,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_svtranspose( 
    magma_tally2_s_matrix x,
    magma_tally2_s_matrix *y,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_s_cucsrtranspose( 
    magma_tally2_s_matrix A, 
    magma_tally2_s_matrix *B,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
s_transpose_csr( 
    magma_tally2_int_t n_rows, 
    magma_tally2_int_t n_cols, 
    magma_tally2_int_t nnz,
    float *val, 
    magma_tally2_index_t *row, 
    magma_tally2_index_t *col, 
    magma_tally2_int_t *new_n_rows, 
    magma_tally2_int_t *new_n_cols, 
    magma_tally2_int_t *new_nnz, 
    float **new_val, 
    magma_tally2_index_t **new_row, 
    magma_tally2_index_t **new_col,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_scsrsplit( 
    magma_tally2_int_t bsize,
    magma_tally2_s_matrix A,
    magma_tally2_s_matrix *D,
    magma_tally2_s_matrix *R,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_smscale( 
    magma_tally2_s_matrix *A, 
    magma_tally2_scale_t scaling,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_smdiff( 
    magma_tally2_s_matrix A, 
    magma_tally2_s_matrix B, 
 real_Double_t *res,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_smdiagadd( 
    magma_tally2_s_matrix *A, 
    float add,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_smsort( 
    magma_tally2_s_matrix *A,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_sindexsort(
    magma_tally2_index_t *x, 
    magma_tally2_int_t first,
    magma_tally2_int_t last,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_sdomainoverlap(
    magma_tally2_index_t num_rows,
    magma_tally2_index_t *num_indices,
    magma_tally2_index_t *rowptr,
    magma_tally2_index_t *colidx,
    magma_tally2_index_t *x,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_ssymbilu( 
    magma_tally2_s_matrix *A, 
    magma_tally2_int_t levels,
    magma_tally2_s_matrix *L,
    magma_tally2_s_matrix *U,
    magma_tally2_queue_t queue );


magma_tally2_int_t 
magma_tally2_swrite_csr_mtx( 
    magma_tally2_s_matrix A,
    magma_tally2_order_t MajorType,
 const char *filename,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_swrite_csrtomtx( 
    magma_tally2_s_matrix A,
    const char *filename,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_sprint_csr( 
    magma_tally2_int_t n_row, 
    magma_tally2_int_t n_col, 
    magma_tally2_int_t nnz, 
    float **val, 
    magma_tally2_index_t **row, 
    magma_tally2_index_t **col,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_sprint_csr_mtx( 
    magma_tally2_int_t n_row, 
    magma_tally2_int_t n_col, 
    magma_tally2_int_t nnz, 
    float **val, 
    magma_tally2_index_t **row, 
    magma_tally2_index_t **col, 
    magma_tally2_order_t MajorType,
    magma_tally2_queue_t queue );


magma_tally2_int_t 
magma_tally2_smtranspose(
    magma_tally2_s_matrix A, 
    magma_tally2_s_matrix *B,
    magma_tally2_queue_t queue );


magma_tally2_int_t 
magma_tally2_smtransfer(
    magma_tally2_s_matrix A, 
    magma_tally2_s_matrix *B, 
    magma_tally2_location_t src, 
    magma_tally2_location_t dst,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_smconvert(
    magma_tally2_s_matrix A, 
    magma_tally2_s_matrix *B, 
    magma_tally2_storage_t old_format, 
    magma_tally2_storage_t new_format,
    magma_tally2_queue_t queue );


magma_tally2_int_t
magma_tally2_svinit(
    magma_tally2_s_matrix *x, 
    magma_tally2_location_t memory_location,
    magma_tally2_int_t num_rows, 
    magma_tally2_int_t num_cols,
    float values,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_sprint_vector(
    magma_tally2_s_matrix x, 
    magma_tally2_int_t offset, 
    magma_tally2_int_t displaylength,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_svread(
    magma_tally2_s_matrix *x, 
    magma_tally2_int_t length,
    char * filename,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_svspread(
    magma_tally2_s_matrix *x, 
    const char * filename,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_sprint_matrix(
    magma_tally2_s_matrix A,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_sdiameter(
    magma_tally2_s_matrix *A,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_srowentries(
    magma_tally2_s_matrix *A,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_smfree(
    magma_tally2_s_matrix *A,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_sresidual(
    magma_tally2_s_matrix A, 
    magma_tally2_s_matrix b, 
    magma_tally2_s_matrix x, 
    float *res,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_sresidualvec(
    magma_tally2_s_matrix A,
    magma_tally2_s_matrix b,
    magma_tally2_s_matrix x,
    magma_tally2_s_matrix *r,
    float *res,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_smgenerator(
    magma_tally2_int_t n,
    magma_tally2_int_t offdiags,
    magma_tally2_index_t *diag_offset,
    float *diag_vals,
    magma_tally2_s_matrix *A,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_sm_27stencil(
    magma_tally2_int_t n,
    magma_tally2_s_matrix *A,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_sm_5stencil(
    magma_tally2_int_t n,
    magma_tally2_s_matrix *A,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_ssolverinfo(
    magma_tally2_s_solver_par *solver_par, 
    magma_tally2_s_preconditioner *precond_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_ssolverinfo_init(
    magma_tally2_s_solver_par *solver_par, 
    magma_tally2_s_preconditioner *precond,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_seigensolverinfo_init(
    magma_tally2_s_solver_par *solver_par,
    magma_tally2_queue_t queue );


magma_tally2_int_t
magma_tally2_ssolverinfo_free(
    magma_tally2_s_solver_par *solver_par, 
    magma_tally2_s_preconditioner *precond,
    magma_tally2_queue_t queue );



/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally2_SPARSE iterative incomplete factorizations
*/


magma_tally2_int_t
magma_tally2_siterilusetup( 
    magma_tally2_s_matrix A, 
    magma_tally2_s_matrix b,                                 
    magma_tally2_s_preconditioner *precond,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_sitericsetup( 
    magma_tally2_s_matrix A, 
    magma_tally2_s_matrix b, 
    magma_tally2_s_preconditioner *precond,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_sitericupdate( 
    magma_tally2_s_matrix A, 
    magma_tally2_s_preconditioner *precond, 
    magma_tally2_int_t updates,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_sapplyiteric_l( 
    magma_tally2_s_matrix b, 
    magma_tally2_s_matrix *x, 
    magma_tally2_s_preconditioner *precond,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_sapplyiteric_r( 
    magma_tally2_s_matrix b, 
    magma_tally2_s_matrix *x, 
    magma_tally2_s_preconditioner *precond,
    magma_tally2_queue_t queue );


magma_tally2_int_t
magma_tally2_siterilu_csr( 
    magma_tally2_s_matrix A,
    magma_tally2_s_matrix L,
    magma_tally2_s_matrix U,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_siteric_csr( 
    magma_tally2_s_matrix A,
    magma_tally2_s_matrix A_CSR,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_sfrobenius( 
    magma_tally2_s_matrix A, 
    magma_tally2_s_matrix B, 
    real_Double_t *res,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_snonlinres(   
    magma_tally2_s_matrix A, 
    magma_tally2_s_matrix L,
    magma_tally2_s_matrix U, 
    magma_tally2_s_matrix *LU, 
    real_Double_t *res,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_silures(   
    magma_tally2_s_matrix A, 
    magma_tally2_s_matrix L,
    magma_tally2_s_matrix U, 
    magma_tally2_s_matrix *LU, 
    real_Double_t *res,
    real_Double_t *nonlinres,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_sicres(       
    magma_tally2_s_matrix A, 
    magma_tally2_s_matrix C,
    magma_tally2_s_matrix CT, 
    magma_tally2_s_matrix *LU, 
    real_Double_t *res,
    real_Double_t *nonlinres,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_sinitguess( 
    magma_tally2_s_matrix A, 
    magma_tally2_s_matrix *L, 
    magma_tally2_s_matrix *U,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_sinitrecursiveLU( 
    magma_tally2_s_matrix A, 
    magma_tally2_s_matrix *B,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_smLdiagadd( 
    magma_tally2_s_matrix *L,
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
magma_tally2_scg(
    magma_tally2_s_matrix A, magma_tally2_s_matrix b, 
    magma_tally2_s_matrix *x, magma_tally2_s_solver_par *solver_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_scg_res(
    magma_tally2_s_matrix A, magma_tally2_s_matrix b, 
    magma_tally2_s_matrix *x, magma_tally2_s_solver_par *solver_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_scg_merge(
    magma_tally2_s_matrix A, magma_tally2_s_matrix b, 
    magma_tally2_s_matrix *x, magma_tally2_s_solver_par *solver_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_sgmres(
    magma_tally2_s_matrix A, magma_tally2_s_matrix b, 
    magma_tally2_s_matrix *x, magma_tally2_s_solver_par *solver_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_sbicgstab(
    magma_tally2_s_matrix A, magma_tally2_s_matrix b, magma_tally2_s_matrix *x, 
    magma_tally2_s_solver_par *solver_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_sbicgstab_merge(
    magma_tally2_s_matrix A, magma_tally2_s_matrix b, 
    magma_tally2_s_matrix *x, magma_tally2_s_solver_par *solver_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_sbicgstab_merge2(
    magma_tally2_s_matrix A, magma_tally2_s_matrix b, 
    magma_tally2_s_matrix *x, magma_tally2_s_solver_par *solver_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_spcg(
    magma_tally2_s_matrix A, magma_tally2_s_matrix b, 
    magma_tally2_s_matrix *x, magma_tally2_s_solver_par *solver_par, 
    magma_tally2_s_preconditioner *precond_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_sbpcg(
    magma_tally2_s_matrix A, magma_tally2_s_matrix b, 
    magma_tally2_s_matrix *x, magma_tally2_s_solver_par *solver_par, 
    magma_tally2_s_preconditioner *precond_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_spbicgstab(
    magma_tally2_s_matrix A, magma_tally2_s_matrix b, 
    magma_tally2_s_matrix *x, magma_tally2_s_solver_par *solver_par, 
    magma_tally2_s_preconditioner *precond_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_spgmres(
    magma_tally2_s_matrix A, magma_tally2_s_matrix b, 
    magma_tally2_s_matrix *x, magma_tally2_s_solver_par *solver_par, 
    magma_tally2_s_preconditioner *precond_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_sfgmres(
    magma_tally2_s_matrix A, magma_tally2_s_matrix b, 
    magma_tally2_s_matrix *x, magma_tally2_s_solver_par *solver_par, 
    magma_tally2_s_preconditioner *precond_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_sjacobi(
    magma_tally2_s_matrix A, magma_tally2_s_matrix b, 
    magma_tally2_s_matrix *x, magma_tally2_s_solver_par *solver_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_sjacobidomainoverlap(
    magma_tally2_s_matrix A, 
    magma_tally2_s_matrix b, 
    magma_tally2_s_matrix *x,  
    magma_tally2_s_solver_par *solver_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_sbaiter(
    magma_tally2_s_matrix A, magma_tally2_s_matrix b, 
    magma_tally2_s_matrix *x, magma_tally2_s_solver_par *solver_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_siterref(
    magma_tally2_s_matrix A, magma_tally2_s_matrix b, 
    magma_tally2_s_matrix *x, magma_tally2_s_solver_par *solver_par, 
    magma_tally2_s_preconditioner *precond_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_silu(
    magma_tally2_s_matrix A, magma_tally2_s_matrix b, 
    magma_tally2_s_matrix *x, magma_tally2_s_solver_par *solver_par, 
    magma_tally2_int_t *ipiv,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_sbcsrlu(
    magma_tally2_s_matrix A, magma_tally2_s_matrix b, 
    magma_tally2_s_matrix *x, magma_tally2_s_solver_par *solver_par,
    magma_tally2_queue_t queue );


magma_tally2_int_t
magma_tally2_sbcsrlutrf(
    magma_tally2_s_matrix A, 
    magma_tally2_s_matrix *M,
    magma_tally2_int_t *ipiv, 
    magma_tally2_int_t version,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_sbcsrlusv(
    magma_tally2_s_matrix A, magma_tally2_s_matrix b, 
    magma_tally2_s_matrix *x, magma_tally2_s_solver_par *solver_par, 
    magma_tally2_int_t *ipiv,
    magma_tally2_queue_t queue );



magma_tally2_int_t
magma_tally2_silucg(
    magma_tally2_s_matrix A, magma_tally2_s_matrix b, 
    magma_tally2_s_matrix *x, magma_tally2_s_solver_par *solver_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_silugmres(
    magma_tally2_s_matrix A, magma_tally2_s_matrix b, 
    magma_tally2_s_matrix *x, magma_tally2_s_solver_par *solver_par,
    magma_tally2_queue_t queue ); 


magma_tally2_int_t
magma_tally2_slobpcg_shift(
    magma_tally2_int_t num_rows,
    magma_tally2_int_t num_vecs, 
    magma_tally2_int_t shift,
    magma_tally2Float_ptr x,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_slobpcg_res(
    magma_tally2_int_t num_rows,
    magma_tally2_int_t num_vecs, 
    float *evalues, 
    magma_tally2Float_ptr X,
    magma_tally2Float_ptr R, 
    float *res,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_slobpcg_maxpy(
    magma_tally2_int_t num_rows,
    magma_tally2_int_t num_vecs, 
    magma_tally2Float_ptr X,
    magma_tally2Float_ptr Y,
    magma_tally2_queue_t queue );


/*/////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally2_SPARSE eigensolvers (Data on GPU)
*/
magma_tally2_int_t
magma_tally2_slobpcg(
    magma_tally2_s_matrix A, 
    magma_tally2_s_solver_par *solver_par,
    magma_tally2_s_preconditioner *precond_par, 
    magma_tally2_queue_t queue );




/*/////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally2_SPARSE preconditioners (Data on GPU)
*/

magma_tally2_int_t
magma_tally2_sjacobisetup(
    magma_tally2_s_matrix A, 
    magma_tally2_s_matrix b, 
    magma_tally2_s_matrix *M, 
    magma_tally2_s_matrix *c,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_sjacobisetup_matrix(
    magma_tally2_s_matrix A, 
    magma_tally2_s_matrix *M, 
    magma_tally2_s_matrix *d,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_sjacobisetup_vector(
    magma_tally2_s_matrix b, 
    magma_tally2_s_matrix d, 
    magma_tally2_s_matrix *c,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_sjacobiiter(
    magma_tally2_s_matrix M, 
    magma_tally2_s_matrix c, 
    magma_tally2_s_matrix *x, 
    magma_tally2_s_solver_par *solver_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_sjacobiiter_precond( 
    magma_tally2_s_matrix M, 
    magma_tally2_s_matrix *x, 
    magma_tally2_s_solver_par *solver_par, 
    magma_tally2_s_preconditioner *precond,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_sjacobiiter_sys(
    magma_tally2_s_matrix A, 
    magma_tally2_s_matrix b, 
    magma_tally2_s_matrix d, 
    magma_tally2_s_matrix t, 
    magma_tally2_s_matrix *x,  
    magma_tally2_s_solver_par *solver_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_spastixsetup(
    magma_tally2_s_matrix A, magma_tally2_s_matrix b,
    magma_tally2_s_preconditioner *precond,
    magma_tally2_queue_t queue );


magma_tally2_int_t
magma_tally2_sapplypastix(
    magma_tally2_s_matrix b, magma_tally2_s_matrix *x, 
    magma_tally2_s_preconditioner *precond,
    magma_tally2_queue_t queue );


// custom preconditioner
magma_tally2_int_t
magma_tally2_sapplycustomprecond_l(
    magma_tally2_s_matrix b, 
    magma_tally2_s_matrix *x, 
    magma_tally2_s_preconditioner *precond,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_sapplycustomprecond_r(
    magma_tally2_s_matrix b, 
    magma_tally2_s_matrix *x, 
    magma_tally2_s_preconditioner *precond,
    magma_tally2_queue_t queue );


// CUSPARSE preconditioner

magma_tally2_int_t
magma_tally2_scuilusetup(
    magma_tally2_s_matrix A, magma_tally2_s_preconditioner *precond,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_sapplycuilu_l(
    magma_tally2_s_matrix b, magma_tally2_s_matrix *x, 
    magma_tally2_s_preconditioner *precond,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_sapplycuilu_r(
    magma_tally2_s_matrix b, magma_tally2_s_matrix *x, 
    magma_tally2_s_preconditioner *precond,
    magma_tally2_queue_t queue );


magma_tally2_int_t
magma_tally2_scuiccsetup(
    magma_tally2_s_matrix A, magma_tally2_s_preconditioner *precond,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_sapplycuicc_l(
    magma_tally2_s_matrix b, magma_tally2_s_matrix *x, 
    magma_tally2_s_preconditioner *precond,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_sapplycuicc_r(
    magma_tally2_s_matrix b, magma_tally2_s_matrix *x, 
    magma_tally2_s_preconditioner *precond,
    magma_tally2_queue_t queue );


magma_tally2_int_t
magma_tally2_scumilusetup(
    magma_tally2_s_matrix A, 
    magma_tally2_s_preconditioner *precond,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_scumilugeneratesolverinfo(
    magma_tally2_s_preconditioner *precond,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_sapplycumilu_l(
    magma_tally2_s_matrix b, 
    magma_tally2_s_matrix *x, 
    magma_tally2_s_preconditioner *precond,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_sapplycumilu_r(
    magma_tally2_s_matrix b, 
    magma_tally2_s_matrix *x, 
    magma_tally2_s_preconditioner *precond,
    magma_tally2_queue_t queue );


magma_tally2_int_t
magma_tally2_scumiccsetup(
    magma_tally2_s_matrix A, 
    magma_tally2_s_preconditioner *precond,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_scumicgeneratesolverinfo(
    magma_tally2_s_preconditioner *precond,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_sapplycumicc_l(
    magma_tally2_s_matrix b, 
    magma_tally2_s_matrix *x, 
    magma_tally2_s_preconditioner *precond,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_sapplycumicc_r(
    magma_tally2_s_matrix b, 
    magma_tally2_s_matrix *x, 
    magma_tally2_s_preconditioner *precond,
    magma_tally2_queue_t queue );


// block-asynchronous iteration

magma_tally2_int_t
magma_tally2_sbajac_csr(
    magma_tally2_int_t localiters,
    magma_tally2_s_matrix D,
    magma_tally2_s_matrix R,
    magma_tally2_s_matrix b,
    magma_tally2_s_matrix *x,
    magma_tally2_queue_t queue );

/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally2_SPARSE utility function definitions
*/

magma_tally2_int_t
magma_tally2_s_spmv(
    float alpha, 
    magma_tally2_s_matrix A, 
    magma_tally2_s_matrix x, 
    float beta, 
    magma_tally2_s_matrix y,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_scustomspmv(
    float alpha, 
    magma_tally2_s_matrix x, 
    float beta, 
    magma_tally2_s_matrix y,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_s_spmv_shift(
    float alpha, 
    magma_tally2_s_matrix A, 
    float lambda,
    magma_tally2_s_matrix x, 
    float beta, 
    magma_tally2_int_t offset, 
    magma_tally2_int_t blocksize,
    magma_tally2Index_ptr dadd_vecs, 
    magma_tally2_s_matrix y,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_scuspmm(
    magma_tally2_s_matrix A, 
    magma_tally2_s_matrix B, 
    magma_tally2_s_matrix *AB,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_s_spmm(
    float alpha, 
    magma_tally2_s_matrix A,
    magma_tally2_s_matrix B,
    magma_tally2_s_matrix *C,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_ssymbilu( 
    magma_tally2_s_matrix *A, 
    magma_tally2_int_t levels, 
    magma_tally2_s_matrix *L, 
    magma_tally2_s_matrix *U,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_scuspaxpy(
    magma_tally2Float_ptr alpha, magma_tally2_s_matrix A, 
    magma_tally2Float_ptr beta, magma_tally2_s_matrix B, 
    magma_tally2_s_matrix *AB,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_s_precond(
    magma_tally2_s_matrix A, 
    magma_tally2_s_matrix b, magma_tally2_s_matrix *x,
    magma_tally2_s_preconditioner *precond,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_s_solver(
    magma_tally2_s_matrix A, magma_tally2_s_matrix b, 
    magma_tally2_s_matrix *x, magma_tally2_sopts *zopts,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_s_precondsetup(
    magma_tally2_s_matrix A, magma_tally2_s_matrix b, 
    magma_tally2_s_preconditioner *precond,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_s_applyprecond(
    magma_tally2_s_matrix A, magma_tally2_s_matrix b, 
    magma_tally2_s_matrix *x, magma_tally2_s_preconditioner *precond,
    magma_tally2_queue_t queue );


magma_tally2_int_t
magma_tally2_s_applyprecond_left(
    magma_tally2_s_matrix A, magma_tally2_s_matrix b, 
    magma_tally2_s_matrix *x, magma_tally2_s_preconditioner *precond,
    magma_tally2_queue_t queue );


magma_tally2_int_t
magma_tally2_s_applyprecond_right(
    magma_tally2_s_matrix A, magma_tally2_s_matrix b, 
    magma_tally2_s_matrix *x, magma_tally2_s_preconditioner *precond,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_s_initP2P(
    magma_tally2_int_t *bandwidth_benchmark,
    magma_tally2_int_t *num_gpus,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_scompact(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2Float_ptr dA, magma_tally2_int_t ldda,
    float *dnorms, float tol, 
    magma_tally2_int_t *activeMask, magma_tally2_int_t *cBlockSize,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_scompactActive(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2Float_ptr dA, magma_tally2_int_t ldda, 
    magma_tally2_int_t *active,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_smlumerge(    
    magma_tally2_s_matrix L, 
    magma_tally2_s_matrix U,
    magma_tally2_s_matrix *A, 
    magma_tally2_queue_t queue );


/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally2_SPARSE BLAS function definitions
*/
magma_tally2_int_t 
magma_tally2_sgecsrmv(
    magma_tally2_trans_t transA,
    magma_tally2_int_t m, magma_tally2_int_t n,
    float alpha,
    magma_tally2Float_ptr dval,
    magma_tally2Index_ptr drowptr,
    magma_tally2Index_ptr dcolind,
    magma_tally2Float_ptr dx,
    float beta,
    magma_tally2Float_ptr dy,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_sgecsrmv_shift(
    magma_tally2_trans_t transA,
    magma_tally2_int_t m, magma_tally2_int_t n,
    float alpha,
    float lambda,
    magma_tally2Float_ptr dval,
    magma_tally2Index_ptr drowptr,
    magma_tally2Index_ptr dcolind,
    magma_tally2Float_ptr dx,
    float beta,
    int offset,
    int blocksize,
    magma_tally2Index_ptr dadd_rows,
    magma_tally2Float_ptr dy,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_smgecsrmv(
    magma_tally2_trans_t transA,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2_int_t num_vecs,
    float alpha,
    magma_tally2Float_ptr dval,
    magma_tally2Index_ptr drowptr,
    magma_tally2Index_ptr dcolind,
    magma_tally2Float_ptr dx,
    float beta,
    magma_tally2Float_ptr dy,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_sgeellmv(
    magma_tally2_trans_t transA,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2_int_t nnz_per_row,
    float alpha,
    magma_tally2Float_ptr dval,
    magma_tally2Index_ptr dcolind,
    magma_tally2Float_ptr dx,
    float beta,
    magma_tally2Float_ptr dy,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_sgeellmv_shift(
    magma_tally2_trans_t transA,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2_int_t nnz_per_row,
    float alpha,
    float lambda,
    magma_tally2Float_ptr dval,
    magma_tally2Index_ptr dcolind,
    magma_tally2Float_ptr dx,
    float beta,
    int offset,
    int blocksize,
    magma_tally2Index_ptr dadd_rows,
    magma_tally2Float_ptr dy,
    magma_tally2_queue_t queue );


magma_tally2_int_t 
magma_tally2_smgeellmv(
    magma_tally2_trans_t transA,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2_int_t num_vecs,
    magma_tally2_int_t nnz_per_row,
    float alpha,
    magma_tally2Float_ptr dval,
    magma_tally2Index_ptr dcolind,
    magma_tally2Float_ptr dx,
    float beta,
    magma_tally2Float_ptr dy,
    magma_tally2_queue_t queue );


magma_tally2_int_t 
magma_tally2_sgeelltmv(
    magma_tally2_trans_t transA,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2_int_t nnz_per_row,
    float alpha,
    magma_tally2Float_ptr dval,
    magma_tally2Index_ptr dcolind,
    magma_tally2Float_ptr dx,
    float beta,
    magma_tally2Float_ptr dy,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_sgeelltmv_shift(
    magma_tally2_trans_t transA,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2_int_t nnz_per_row,
    float alpha,
    float lambda,
    magma_tally2Float_ptr dval,
    magma_tally2Index_ptr dcolind,
    magma_tally2Float_ptr dx,
    float beta,
    int offset,
    int blocksize,
    magma_tally2Index_ptr dadd_rows,
    magma_tally2Float_ptr dy,
    magma_tally2_queue_t queue );


magma_tally2_int_t 
magma_tally2_smgeelltmv(
    magma_tally2_trans_t transA,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2_int_t num_vecs,
    magma_tally2_int_t nnz_per_row,
    float alpha,
    magma_tally2Float_ptr dval,
    magma_tally2Index_ptr dcolind,
    magma_tally2Float_ptr dx,
    float beta,
    magma_tally2Float_ptr dy,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_sgeellrtmv(
    magma_tally2_trans_t transA,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2_int_t nnz_per_row,
    float alpha,
    magma_tally2Float_ptr dval,
    magma_tally2Index_ptr dcolind,
    magma_tally2Index_ptr drowlength,
    magma_tally2Float_ptr dx,
    float beta,
    magma_tally2Float_ptr dy,
    magma_tally2_int_t num_threads,
    magma_tally2_int_t threads_per_row,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_sgesellcmv(
    magma_tally2_trans_t transA,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2_int_t blocksize,
    magma_tally2_int_t slices,
    magma_tally2_int_t alignment,
    float alpha,
    magma_tally2Float_ptr dval,
    magma_tally2Index_ptr dcolind,
    magma_tally2Index_ptr drowptr,
    magma_tally2Float_ptr dx,
    float beta,
    magma_tally2Float_ptr dy,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_sgesellpmv(
    magma_tally2_trans_t transA,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2_int_t blocksize,
    magma_tally2_int_t slices,
    magma_tally2_int_t alignment,
    float alpha,
    magma_tally2Float_ptr dval,
    magma_tally2Index_ptr dcolind,
    magma_tally2Index_ptr drowptr,
    magma_tally2Float_ptr dx,
    float beta,
    magma_tally2Float_ptr dy,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_smgesellpmv(
    magma_tally2_trans_t transA,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2_int_t num_vecs,
    magma_tally2_int_t blocksize,
    magma_tally2_int_t slices,
    magma_tally2_int_t alignment,
    float alpha,
    magma_tally2Float_ptr dval,
    magma_tally2Index_ptr dcolind,
    magma_tally2Index_ptr drowptr,
    magma_tally2Float_ptr dx,
    float beta,
    magma_tally2Float_ptr dy,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_smgesellpmv_blocked(
    magma_tally2_trans_t transA,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2_int_t num_vecs,
    magma_tally2_int_t blocksize,
    magma_tally2_int_t slices,
    magma_tally2_int_t alignment,
    float alpha,
    magma_tally2Float_ptr dval,
    magma_tally2Index_ptr dcolind,
    magma_tally2Index_ptr drowptr,
    magma_tally2Float_ptr dx,
    float beta,
    magma_tally2Float_ptr dy,
    magma_tally2_queue_t queue );


magma_tally2_int_t
magma_tally2_smergedgs(
    magma_tally2_int_t n, 
    magma_tally2_int_t ldh,
    magma_tally2_int_t k, 
    magma_tally2Float_ptr dv, 
    magma_tally2Float_ptr dr,
    magma_tally2Float_ptr dskp,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_scopyscale(    
    int n, 
    int k,
    magma_tally2Float_ptr dr, 
    magma_tally2Float_ptr dv,
    magma_tally2Float_ptr dskp,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_snrm2scale(    
    int m, 
    magma_tally2Float_ptr dr,    
    int lddr, 
    float *drnorm,
    magma_tally2_queue_t queue );


magma_tally2_int_t
magma_tally2_sjacobisetup_vector_gpu(
    int num_rows, 
    magma_tally2_s_matrix b, 
    magma_tally2_s_matrix d, 
    magma_tally2_s_matrix c,
    magma_tally2_s_matrix *x,
    magma_tally2_queue_t queue );


magma_tally2_int_t
magma_tally2_sjacobi_diagscal(    
    int num_rows, 
    magma_tally2_s_matrix d, 
    magma_tally2_s_matrix b, 
    magma_tally2_s_matrix *c,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_sjacobiupdate(
    magma_tally2_s_matrix t, 
    magma_tally2_s_matrix b, 
    magma_tally2_s_matrix d, 
    magma_tally2_s_matrix *x,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_sjacobispmvupdate(
    magma_tally2_int_t maxiter,
    magma_tally2_s_matrix A, 
    magma_tally2_s_matrix t, 
    magma_tally2_s_matrix b, 
    magma_tally2_s_matrix d, 
    magma_tally2_s_matrix *x,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_sjacobispmvupdate_bw(
    magma_tally2_int_t maxiter,
    magma_tally2_s_matrix A, 
    magma_tally2_s_matrix t, 
    magma_tally2_s_matrix b, 
    magma_tally2_s_matrix d, 
    magma_tally2_s_matrix *x,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_sjacobispmvupdateselect(
    magma_tally2_int_t maxiter,
    magma_tally2_int_t num_updates,
    magma_tally2_index_t *indices,
    magma_tally2_s_matrix A,
    magma_tally2_s_matrix t, 
    magma_tally2_s_matrix b, 
    magma_tally2_s_matrix d, 
    magma_tally2_s_matrix tmp, 
    magma_tally2_s_matrix *x,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_sjacobisetup_diagscal(
    magma_tally2_s_matrix A, magma_tally2_s_matrix *d,
    magma_tally2_queue_t queue );


magma_tally2_int_t
magma_tally2_sbicgmerge1(    
    int n, 
    magma_tally2Float_ptr dskp,
    magma_tally2Float_ptr dv, 
    magma_tally2Float_ptr dr, 
    magma_tally2Float_ptr dp,
    magma_tally2_queue_t queue );


magma_tally2_int_t
magma_tally2_sbicgmerge2(
    int n, 
    magma_tally2Float_ptr dskp, 
    magma_tally2Float_ptr dr,
    magma_tally2Float_ptr dv, 
    magma_tally2Float_ptr ds,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_sbicgmerge3(
    int n, 
    magma_tally2Float_ptr dskp, 
    magma_tally2Float_ptr dp,
    magma_tally2Float_ptr ds,
    magma_tally2Float_ptr dt,
    magma_tally2Float_ptr dx, 
    magma_tally2Float_ptr dr,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_sbicgmerge4(
    int type, 
    magma_tally2Float_ptr dskp,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_scgmerge_spmv1( 
    magma_tally2_s_matrix A,
    magma_tally2Float_ptr d1,
    magma_tally2Float_ptr d2,
    magma_tally2Float_ptr dd,
    magma_tally2Float_ptr dz,
    magma_tally2Float_ptr dskp,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_scgmerge_xrbeta( 
    int n,
    magma_tally2Float_ptr d1,
    magma_tally2Float_ptr d2,
    magma_tally2Float_ptr dx,
    magma_tally2Float_ptr dr,
    magma_tally2Float_ptr dd,
    magma_tally2Float_ptr dz, 
    magma_tally2Float_ptr dskp,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_smdotc(
    magma_tally2_int_t n, 
    magma_tally2_int_t k, 
    magma_tally2Float_ptr dv, 
    magma_tally2Float_ptr dr,
    magma_tally2Float_ptr dd1,
    magma_tally2Float_ptr dd2,
    magma_tally2Float_ptr dskp,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_sgemvmdot(
    int n, 
    int k, 
    magma_tally2Float_ptr dv, 
    magma_tally2Float_ptr dr,
    magma_tally2Float_ptr dd1,
    magma_tally2Float_ptr dd2,
    magma_tally2Float_ptr dskp,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_sbicgmerge_spmv1( 
    magma_tally2_s_matrix A,
    magma_tally2Float_ptr dd1,
    magma_tally2Float_ptr dd2,
    magma_tally2Float_ptr dp,
    magma_tally2Float_ptr dr,
    magma_tally2Float_ptr dv,
    magma_tally2Float_ptr dskp,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_sbicgmerge_spmv2( 
    magma_tally2_s_matrix A,
    magma_tally2Float_ptr dd1,
    magma_tally2Float_ptr dd2,
    magma_tally2Float_ptr ds,
    magma_tally2Float_ptr dt,
    magma_tally2Float_ptr dskp,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_sbicgmerge_xrbeta( 
    int n,
    magma_tally2Float_ptr dd1,
    magma_tally2Float_ptr dd2,
    magma_tally2Float_ptr drr,
    magma_tally2Float_ptr dr,
    magma_tally2Float_ptr dp,
    magma_tally2Float_ptr ds,
    magma_tally2Float_ptr dt,
    magma_tally2Float_ptr dx, 
    magma_tally2Float_ptr dskp,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_sbcsrswp(
    magma_tally2_int_t n,
    magma_tally2_int_t size_b, 
    magma_tally2_int_t *ipiv,
    magma_tally2Float_ptr dx,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_sbcsrtrsv(
    magma_tally2_uplo_t uplo,
    magma_tally2_int_t r_blocks,
    magma_tally2_int_t c_blocks,
    magma_tally2_int_t size_b, 
    magma_tally2Float_ptr dA,
    magma_tally2_index_t *blockinfo, 
    magma_tally2Float_ptr dx,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_sbcsrvalcpy(
    magma_tally2_int_t size_b, 
    magma_tally2_int_t num_blocks, 
    magma_tally2_int_t num_zero_blocks, 
    magma_tally2Float_ptr *dAval, 
    magma_tally2Float_ptr *dBval,
    magma_tally2Float_ptr *dBval2,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_sbcsrluegemm(
    magma_tally2_int_t size_b, 
    magma_tally2_int_t num_block_rows,
    magma_tally2_int_t kblocks,
    magma_tally2Float_ptr *dA, 
    magma_tally2Float_ptr *dB, 
    magma_tally2Float_ptr *dC,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_sbcsrlupivloc(
    magma_tally2_int_t size_b, 
    magma_tally2_int_t kblocks,
    magma_tally2Float_ptr *dA, 
    magma_tally2_int_t *ipiv,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_sbcsrblockinfo5(
    magma_tally2_int_t lustep,
    magma_tally2_int_t num_blocks, 
    magma_tally2_int_t c_blocks, 
    magma_tally2_int_t size_b,
    magma_tally2_index_t *blockinfo,
    magma_tally2Float_ptr dval,
    magma_tally2Float_ptr *AII,
    magma_tally2_queue_t queue );

 
#ifdef __cplusplus
}
#endif

#undef PRECISION_s
#endif /* MAGMA_tally2SPARSE_S_H */
