/*
 -- MAGMA_tally3 (version 1.1) --
 Univ. of Tennessee, Knoxville
 Univ. of California, Berkeley
 Univ. of Colorado, Denver
 @date

 @precisions normal z -> s d c
 @author Hartwig Anzt
*/

#ifndef MAGMA_tally3SPARSE_Z_H
#define MAGMA_tally3SPARSE_Z_H

#include "magma_tally3_types.h"
#include "magma_tally3sparse_types.h"

#define PRECISION_z


#ifdef __cplusplus
extern "C" {
#endif


/* ////////////////////////////////////////////////////////////////////////////
 -- For backwards compatability, map old (1.6.1) to new (1.6.2) function names
*/

#define magma_tally3_z_mtranspose  magma_tally3_zmtranspose
#define magma_tally3_z_mtransfer   magma_tally3_zmtransfer
#define magma_tally3_z_vtransfer   magma_tally3_zmtransfer
#define magma_tally3_z_mconvert    magma_tally3_zmconvert
#define magma_tally3_z_vinit       magma_tally3_zvinit
#define magma_tally3_z_vvisu       magma_tally3_zprint_vector
#define magma_tally3_z_vread       magma_tally3_zvread
#define magma_tally3_z_vspread     magma_tally3_zvspread
#define magma_tally3_z_mvisu       magma_tally3_zprint_matrix
#define magma_tally3_z_mfree       magma_tally3_zmfree
#define magma_tally3_z_vfree       magma_tally3_zmfree
#define write_z_csr_mtx     magma_tally3_zwrite_csr_mtx
#define write_z_csrtomtx    magma_tally3_zwrite_csrtomtx
#define print_z_csr         magma_tally3_zprint_csr_mtx


/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally3_SPARSE Auxiliary functions
*/


magma_tally3_int_t
magma_tally3_zparse_opts( 
    int argc, 
    char** argv, 
    magma_tally3_zopts *opts, 
    int *matrices, 
    magma_tally3_queue_t queue );

magma_tally3_int_t 
read_z_csr_from_binary( 
    magma_tally3_int_t* n_row, 
    magma_tally3_int_t* n_col, 
    magma_tally3_int_t* nnz, 
    magma_tally3DoubleComplex **val, 
    magma_tally3_index_t **row, 
    magma_tally3_index_t **col,
    const char * filename,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
read_z_csr_from_mtx( 
    magma_tally3_storage_t *type, 
    magma_tally3_location_t *location,
    magma_tally3_int_t* n_row, 
    magma_tally3_int_t* n_col, 
    magma_tally3_int_t* nnz, 
    magma_tally3DoubleComplex **val, 
    magma_tally3_index_t **row, 
    magma_tally3_index_t **col, 
    const char *filename,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_z_csr_mtx( 
    magma_tally3_z_matrix *A, 
    const char *filename,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_zcsrset( 
    magma_tally3_int_t m, 
    magma_tally3_int_t n, 
    magma_tally3_index_t *row, 
    magma_tally3_index_t *col, 
    magma_tally3DoubleComplex *val,
    magma_tally3_z_matrix *A,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_zcsrget( 
    magma_tally3_z_matrix A,
    magma_tally3_int_t *m, 
    magma_tally3_int_t *n, 
    magma_tally3_index_t **row, 
    magma_tally3_index_t **col, 
    magma_tally3DoubleComplex **val,
    magma_tally3_queue_t queue );


magma_tally3_int_t 
magma_tally3_zvset( 
    magma_tally3_int_t m, 
    magma_tally3_int_t n, 
    magma_tally3DoubleComplex *val,
    magma_tally3_z_matrix *v,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_zvget( 
    magma_tally3_z_matrix v,
    magma_tally3_int_t *m, 
    magma_tally3_int_t *n, 
    magma_tally3DoubleComplex **val,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_zvset_dev( 
    magma_tally3_int_t m, 
    magma_tally3_int_t n, 
    magma_tally3DoubleComplex_ptr val,
    magma_tally3_z_matrix *v,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_zvget_dev( 
    magma_tally3_z_matrix v,
    magma_tally3_int_t *m, 
    magma_tally3_int_t *n, 
    magma_tally3DoubleComplex_ptr *val,
    magma_tally3_queue_t queue );


magma_tally3_int_t 
magma_tally3_z_csr_mtxsymm( 
    magma_tally3_z_matrix *A, 
    const char *filename,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_z_csr_compressor( 
    magma_tally3DoubleComplex ** val, 
    magma_tally3_index_t ** row, 
    magma_tally3_index_t ** col, 
    magma_tally3DoubleComplex ** valn, 
    magma_tally3_index_t ** rown, 
    magma_tally3_index_t ** coln, 
    magma_tally3_int_t *n,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zmcsrcompressor( 
    magma_tally3_z_matrix *A,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zmcsrcompressor_gpu( 
    magma_tally3_z_matrix *A,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zvtranspose( 
    magma_tally3_z_matrix x,
    magma_tally3_z_matrix *y,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_z_cucsrtranspose( 
    magma_tally3_z_matrix A, 
    magma_tally3_z_matrix *B,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
z_transpose_csr( 
    magma_tally3_int_t n_rows, 
    magma_tally3_int_t n_cols, 
    magma_tally3_int_t nnz,
    magma_tally3DoubleComplex *val, 
    magma_tally3_index_t *row, 
    magma_tally3_index_t *col, 
    magma_tally3_int_t *new_n_rows, 
    magma_tally3_int_t *new_n_cols, 
    magma_tally3_int_t *new_nnz, 
    magma_tally3DoubleComplex **new_val, 
    magma_tally3_index_t **new_row, 
    magma_tally3_index_t **new_col,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zcsrsplit( 
    magma_tally3_int_t bsize,
    magma_tally3_z_matrix A,
    magma_tally3_z_matrix *D,
    magma_tally3_z_matrix *R,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zmscale( 
    magma_tally3_z_matrix *A, 
    magma_tally3_scale_t scaling,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_zmdiff( 
    magma_tally3_z_matrix A, 
    magma_tally3_z_matrix B, 
 real_Double_t *res,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zmdiagadd( 
    magma_tally3_z_matrix *A, 
    magma_tally3DoubleComplex add,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_zmsort( 
    magma_tally3_z_matrix *A,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zindexsort(
    magma_tally3_index_t *x, 
    magma_tally3_int_t first,
    magma_tally3_int_t last,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zdomainoverlap(
    magma_tally3_index_t num_rows,
    magma_tally3_index_t *num_indices,
    magma_tally3_index_t *rowptr,
    magma_tally3_index_t *colidx,
    magma_tally3_index_t *x,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zsymbilu( 
    magma_tally3_z_matrix *A, 
    magma_tally3_int_t levels,
    magma_tally3_z_matrix *L,
    magma_tally3_z_matrix *U,
    magma_tally3_queue_t queue );


magma_tally3_int_t 
magma_tally3_zwrite_csr_mtx( 
    magma_tally3_z_matrix A,
    magma_tally3_order_t MajorType,
 const char *filename,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_zwrite_csrtomtx( 
    magma_tally3_z_matrix A,
    const char *filename,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_zprint_csr( 
    magma_tally3_int_t n_row, 
    magma_tally3_int_t n_col, 
    magma_tally3_int_t nnz, 
    magma_tally3DoubleComplex **val, 
    magma_tally3_index_t **row, 
    magma_tally3_index_t **col,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_zprint_csr_mtx( 
    magma_tally3_int_t n_row, 
    magma_tally3_int_t n_col, 
    magma_tally3_int_t nnz, 
    magma_tally3DoubleComplex **val, 
    magma_tally3_index_t **row, 
    magma_tally3_index_t **col, 
    magma_tally3_order_t MajorType,
    magma_tally3_queue_t queue );


magma_tally3_int_t 
magma_tally3_zmtranspose(
    magma_tally3_z_matrix A, 
    magma_tally3_z_matrix *B,
    magma_tally3_queue_t queue );


magma_tally3_int_t 
magma_tally3_zmtransfer(
    magma_tally3_z_matrix A, 
    magma_tally3_z_matrix *B, 
    magma_tally3_location_t src, 
    magma_tally3_location_t dst,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_zmconvert(
    magma_tally3_z_matrix A, 
    magma_tally3_z_matrix *B, 
    magma_tally3_storage_t old_format, 
    magma_tally3_storage_t new_format,
    magma_tally3_queue_t queue );


magma_tally3_int_t
magma_tally3_zvinit(
    magma_tally3_z_matrix *x, 
    magma_tally3_location_t memory_location,
    magma_tally3_int_t num_rows, 
    magma_tally3_int_t num_cols,
    magma_tally3DoubleComplex values,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zprint_vector(
    magma_tally3_z_matrix x, 
    magma_tally3_int_t offset, 
    magma_tally3_int_t displaylength,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zvread(
    magma_tally3_z_matrix *x, 
    magma_tally3_int_t length,
    char * filename,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zvspread(
    magma_tally3_z_matrix *x, 
    const char * filename,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zprint_matrix(
    magma_tally3_z_matrix A,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_zdiameter(
    magma_tally3_z_matrix *A,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_zrowentries(
    magma_tally3_z_matrix *A,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zmfree(
    magma_tally3_z_matrix *A,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zresidual(
    magma_tally3_z_matrix A, 
    magma_tally3_z_matrix b, 
    magma_tally3_z_matrix x, 
    double *res,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zresidualvec(
    magma_tally3_z_matrix A,
    magma_tally3_z_matrix b,
    magma_tally3_z_matrix x,
    magma_tally3_z_matrix *r,
    double *res,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zmgenerator(
    magma_tally3_int_t n,
    magma_tally3_int_t offdiags,
    magma_tally3_index_t *diag_offset,
    magma_tally3DoubleComplex *diag_vals,
    magma_tally3_z_matrix *A,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zm_27stencil(
    magma_tally3_int_t n,
    magma_tally3_z_matrix *A,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zm_5stencil(
    magma_tally3_int_t n,
    magma_tally3_z_matrix *A,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zsolverinfo(
    magma_tally3_z_solver_par *solver_par, 
    magma_tally3_z_preconditioner *precond_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zsolverinfo_init(
    magma_tally3_z_solver_par *solver_par, 
    magma_tally3_z_preconditioner *precond,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zeigensolverinfo_init(
    magma_tally3_z_solver_par *solver_par,
    magma_tally3_queue_t queue );


magma_tally3_int_t
magma_tally3_zsolverinfo_free(
    magma_tally3_z_solver_par *solver_par, 
    magma_tally3_z_preconditioner *precond,
    magma_tally3_queue_t queue );



/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally3_SPARSE iterative incomplete factorizations
*/


magma_tally3_int_t
magma_tally3_ziterilusetup( 
    magma_tally3_z_matrix A, 
    magma_tally3_z_matrix b,                                 
    magma_tally3_z_preconditioner *precond,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zitericsetup( 
    magma_tally3_z_matrix A, 
    magma_tally3_z_matrix b, 
    magma_tally3_z_preconditioner *precond,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zitericupdate( 
    magma_tally3_z_matrix A, 
    magma_tally3_z_preconditioner *precond, 
    magma_tally3_int_t updates,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zapplyiteric_l( 
    magma_tally3_z_matrix b, 
    magma_tally3_z_matrix *x, 
    magma_tally3_z_preconditioner *precond,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zapplyiteric_r( 
    magma_tally3_z_matrix b, 
    magma_tally3_z_matrix *x, 
    magma_tally3_z_preconditioner *precond,
    magma_tally3_queue_t queue );


magma_tally3_int_t
magma_tally3_ziterilu_csr( 
    magma_tally3_z_matrix A,
    magma_tally3_z_matrix L,
    magma_tally3_z_matrix U,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_ziteric_csr( 
    magma_tally3_z_matrix A,
    magma_tally3_z_matrix A_CSR,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_zfrobenius( 
    magma_tally3_z_matrix A, 
    magma_tally3_z_matrix B, 
    real_Double_t *res,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_znonlinres(   
    magma_tally3_z_matrix A, 
    magma_tally3_z_matrix L,
    magma_tally3_z_matrix U, 
    magma_tally3_z_matrix *LU, 
    real_Double_t *res,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_zilures(   
    magma_tally3_z_matrix A, 
    magma_tally3_z_matrix L,
    magma_tally3_z_matrix U, 
    magma_tally3_z_matrix *LU, 
    real_Double_t *res,
    real_Double_t *nonlinres,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_zicres(       
    magma_tally3_z_matrix A, 
    magma_tally3_z_matrix C,
    magma_tally3_z_matrix CT, 
    magma_tally3_z_matrix *LU, 
    real_Double_t *res,
    real_Double_t *nonlinres,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_zinitguess( 
    magma_tally3_z_matrix A, 
    magma_tally3_z_matrix *L, 
    magma_tally3_z_matrix *U,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_zinitrecursiveLU( 
    magma_tally3_z_matrix A, 
    magma_tally3_z_matrix *B,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_zmLdiagadd( 
    magma_tally3_z_matrix *L,
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
magma_tally3_zcg(
    magma_tally3_z_matrix A, magma_tally3_z_matrix b, 
    magma_tally3_z_matrix *x, magma_tally3_z_solver_par *solver_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_zcg_res(
    magma_tally3_z_matrix A, magma_tally3_z_matrix b, 
    magma_tally3_z_matrix *x, magma_tally3_z_solver_par *solver_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_zcg_merge(
    magma_tally3_z_matrix A, magma_tally3_z_matrix b, 
    magma_tally3_z_matrix *x, magma_tally3_z_solver_par *solver_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_zgmres(
    magma_tally3_z_matrix A, magma_tally3_z_matrix b, 
    magma_tally3_z_matrix *x, magma_tally3_z_solver_par *solver_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zbicgstab(
    magma_tally3_z_matrix A, magma_tally3_z_matrix b, magma_tally3_z_matrix *x, 
    magma_tally3_z_solver_par *solver_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zbicgstab_merge(
    magma_tally3_z_matrix A, magma_tally3_z_matrix b, 
    magma_tally3_z_matrix *x, magma_tally3_z_solver_par *solver_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zbicgstab_merge2(
    magma_tally3_z_matrix A, magma_tally3_z_matrix b, 
    magma_tally3_z_matrix *x, magma_tally3_z_solver_par *solver_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zpcg(
    magma_tally3_z_matrix A, magma_tally3_z_matrix b, 
    magma_tally3_z_matrix *x, magma_tally3_z_solver_par *solver_par, 
    magma_tally3_z_preconditioner *precond_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zbpcg(
    magma_tally3_z_matrix A, magma_tally3_z_matrix b, 
    magma_tally3_z_matrix *x, magma_tally3_z_solver_par *solver_par, 
    magma_tally3_z_preconditioner *precond_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zpbicgstab(
    magma_tally3_z_matrix A, magma_tally3_z_matrix b, 
    magma_tally3_z_matrix *x, magma_tally3_z_solver_par *solver_par, 
    magma_tally3_z_preconditioner *precond_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zpgmres(
    magma_tally3_z_matrix A, magma_tally3_z_matrix b, 
    magma_tally3_z_matrix *x, magma_tally3_z_solver_par *solver_par, 
    magma_tally3_z_preconditioner *precond_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zfgmres(
    magma_tally3_z_matrix A, magma_tally3_z_matrix b, 
    magma_tally3_z_matrix *x, magma_tally3_z_solver_par *solver_par, 
    magma_tally3_z_preconditioner *precond_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zjacobi(
    magma_tally3_z_matrix A, magma_tally3_z_matrix b, 
    magma_tally3_z_matrix *x, magma_tally3_z_solver_par *solver_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zjacobidomainoverlap(
    magma_tally3_z_matrix A, 
    magma_tally3_z_matrix b, 
    magma_tally3_z_matrix *x,  
    magma_tally3_z_solver_par *solver_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zbaiter(
    magma_tally3_z_matrix A, magma_tally3_z_matrix b, 
    magma_tally3_z_matrix *x, magma_tally3_z_solver_par *solver_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_ziterref(
    magma_tally3_z_matrix A, magma_tally3_z_matrix b, 
    magma_tally3_z_matrix *x, magma_tally3_z_solver_par *solver_par, 
    magma_tally3_z_preconditioner *precond_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zilu(
    magma_tally3_z_matrix A, magma_tally3_z_matrix b, 
    magma_tally3_z_matrix *x, magma_tally3_z_solver_par *solver_par, 
    magma_tally3_int_t *ipiv,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zbcsrlu(
    magma_tally3_z_matrix A, magma_tally3_z_matrix b, 
    magma_tally3_z_matrix *x, magma_tally3_z_solver_par *solver_par,
    magma_tally3_queue_t queue );


magma_tally3_int_t
magma_tally3_zbcsrlutrf(
    magma_tally3_z_matrix A, 
    magma_tally3_z_matrix *M,
    magma_tally3_int_t *ipiv, 
    magma_tally3_int_t version,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zbcsrlusv(
    magma_tally3_z_matrix A, magma_tally3_z_matrix b, 
    magma_tally3_z_matrix *x, magma_tally3_z_solver_par *solver_par, 
    magma_tally3_int_t *ipiv,
    magma_tally3_queue_t queue );



magma_tally3_int_t
magma_tally3_zilucg(
    magma_tally3_z_matrix A, magma_tally3_z_matrix b, 
    magma_tally3_z_matrix *x, magma_tally3_z_solver_par *solver_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zilugmres(
    magma_tally3_z_matrix A, magma_tally3_z_matrix b, 
    magma_tally3_z_matrix *x, magma_tally3_z_solver_par *solver_par,
    magma_tally3_queue_t queue ); 


magma_tally3_int_t
magma_tally3_zlobpcg_shift(
    magma_tally3_int_t num_rows,
    magma_tally3_int_t num_vecs, 
    magma_tally3_int_t shift,
    magma_tally3DoubleComplex_ptr x,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zlobpcg_res(
    magma_tally3_int_t num_rows,
    magma_tally3_int_t num_vecs, 
    double *evalues, 
    magma_tally3DoubleComplex_ptr X,
    magma_tally3DoubleComplex_ptr R, 
    double *res,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zlobpcg_maxpy(
    magma_tally3_int_t num_rows,
    magma_tally3_int_t num_vecs, 
    magma_tally3DoubleComplex_ptr X,
    magma_tally3DoubleComplex_ptr Y,
    magma_tally3_queue_t queue );


/*/////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally3_SPARSE eigensolvers (Data on GPU)
*/
magma_tally3_int_t
magma_tally3_zlobpcg(
    magma_tally3_z_matrix A, 
    magma_tally3_z_solver_par *solver_par,
    magma_tally3_z_preconditioner *precond_par, 
    magma_tally3_queue_t queue );




/*/////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally3_SPARSE preconditioners (Data on GPU)
*/

magma_tally3_int_t
magma_tally3_zjacobisetup(
    magma_tally3_z_matrix A, 
    magma_tally3_z_matrix b, 
    magma_tally3_z_matrix *M, 
    magma_tally3_z_matrix *c,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zjacobisetup_matrix(
    magma_tally3_z_matrix A, 
    magma_tally3_z_matrix *M, 
    magma_tally3_z_matrix *d,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zjacobisetup_vector(
    magma_tally3_z_matrix b, 
    magma_tally3_z_matrix d, 
    magma_tally3_z_matrix *c,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zjacobiiter(
    magma_tally3_z_matrix M, 
    magma_tally3_z_matrix c, 
    magma_tally3_z_matrix *x, 
    magma_tally3_z_solver_par *solver_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zjacobiiter_precond( 
    magma_tally3_z_matrix M, 
    magma_tally3_z_matrix *x, 
    magma_tally3_z_solver_par *solver_par, 
    magma_tally3_z_preconditioner *precond,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zjacobiiter_sys(
    magma_tally3_z_matrix A, 
    magma_tally3_z_matrix b, 
    magma_tally3_z_matrix d, 
    magma_tally3_z_matrix t, 
    magma_tally3_z_matrix *x,  
    magma_tally3_z_solver_par *solver_par,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zpastixsetup(
    magma_tally3_z_matrix A, magma_tally3_z_matrix b,
    magma_tally3_z_preconditioner *precond,
    magma_tally3_queue_t queue );


magma_tally3_int_t
magma_tally3_zapplypastix(
    magma_tally3_z_matrix b, magma_tally3_z_matrix *x, 
    magma_tally3_z_preconditioner *precond,
    magma_tally3_queue_t queue );


// custom preconditioner
magma_tally3_int_t
magma_tally3_zapplycustomprecond_l(
    magma_tally3_z_matrix b, 
    magma_tally3_z_matrix *x, 
    magma_tally3_z_preconditioner *precond,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zapplycustomprecond_r(
    magma_tally3_z_matrix b, 
    magma_tally3_z_matrix *x, 
    magma_tally3_z_preconditioner *precond,
    magma_tally3_queue_t queue );


// CUSPARSE preconditioner

magma_tally3_int_t
magma_tally3_zcuilusetup(
    magma_tally3_z_matrix A, magma_tally3_z_preconditioner *precond,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zapplycuilu_l(
    magma_tally3_z_matrix b, magma_tally3_z_matrix *x, 
    magma_tally3_z_preconditioner *precond,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zapplycuilu_r(
    magma_tally3_z_matrix b, magma_tally3_z_matrix *x, 
    magma_tally3_z_preconditioner *precond,
    magma_tally3_queue_t queue );


magma_tally3_int_t
magma_tally3_zcuiccsetup(
    magma_tally3_z_matrix A, magma_tally3_z_preconditioner *precond,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zapplycuicc_l(
    magma_tally3_z_matrix b, magma_tally3_z_matrix *x, 
    magma_tally3_z_preconditioner *precond,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zapplycuicc_r(
    magma_tally3_z_matrix b, magma_tally3_z_matrix *x, 
    magma_tally3_z_preconditioner *precond,
    magma_tally3_queue_t queue );


magma_tally3_int_t
magma_tally3_zcumilusetup(
    magma_tally3_z_matrix A, 
    magma_tally3_z_preconditioner *precond,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zcumilugeneratesolverinfo(
    magma_tally3_z_preconditioner *precond,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zapplycumilu_l(
    magma_tally3_z_matrix b, 
    magma_tally3_z_matrix *x, 
    magma_tally3_z_preconditioner *precond,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zapplycumilu_r(
    magma_tally3_z_matrix b, 
    magma_tally3_z_matrix *x, 
    magma_tally3_z_preconditioner *precond,
    magma_tally3_queue_t queue );


magma_tally3_int_t
magma_tally3_zcumiccsetup(
    magma_tally3_z_matrix A, 
    magma_tally3_z_preconditioner *precond,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zcumicgeneratesolverinfo(
    magma_tally3_z_preconditioner *precond,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zapplycumicc_l(
    magma_tally3_z_matrix b, 
    magma_tally3_z_matrix *x, 
    magma_tally3_z_preconditioner *precond,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zapplycumicc_r(
    magma_tally3_z_matrix b, 
    magma_tally3_z_matrix *x, 
    magma_tally3_z_preconditioner *precond,
    magma_tally3_queue_t queue );


// block-asynchronous iteration

magma_tally3_int_t
magma_tally3_zbajac_csr(
    magma_tally3_int_t localiters,
    magma_tally3_z_matrix D,
    magma_tally3_z_matrix R,
    magma_tally3_z_matrix b,
    magma_tally3_z_matrix *x,
    magma_tally3_queue_t queue );

/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally3_SPARSE utility function definitions
*/

magma_tally3_int_t
magma_tally3_z_spmv(
    magma_tally3DoubleComplex alpha, 
    magma_tally3_z_matrix A, 
    magma_tally3_z_matrix x, 
    magma_tally3DoubleComplex beta, 
    magma_tally3_z_matrix y,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zcustomspmv(
    magma_tally3DoubleComplex alpha, 
    magma_tally3_z_matrix x, 
    magma_tally3DoubleComplex beta, 
    magma_tally3_z_matrix y,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_z_spmv_shift(
    magma_tally3DoubleComplex alpha, 
    magma_tally3_z_matrix A, 
    magma_tally3DoubleComplex lambda,
    magma_tally3_z_matrix x, 
    magma_tally3DoubleComplex beta, 
    magma_tally3_int_t offset, 
    magma_tally3_int_t blocksize,
    magma_tally3Index_ptr dadd_vecs, 
    magma_tally3_z_matrix y,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zcuspmm(
    magma_tally3_z_matrix A, 
    magma_tally3_z_matrix B, 
    magma_tally3_z_matrix *AB,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_z_spmm(
    magma_tally3DoubleComplex alpha, 
    magma_tally3_z_matrix A,
    magma_tally3_z_matrix B,
    magma_tally3_z_matrix *C,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zsymbilu( 
    magma_tally3_z_matrix *A, 
    magma_tally3_int_t levels, 
    magma_tally3_z_matrix *L, 
    magma_tally3_z_matrix *U,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zcuspaxpy(
    magma_tally3DoubleComplex_ptr alpha, magma_tally3_z_matrix A, 
    magma_tally3DoubleComplex_ptr beta, magma_tally3_z_matrix B, 
    magma_tally3_z_matrix *AB,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_z_precond(
    magma_tally3_z_matrix A, 
    magma_tally3_z_matrix b, magma_tally3_z_matrix *x,
    magma_tally3_z_preconditioner *precond,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_z_solver(
    magma_tally3_z_matrix A, magma_tally3_z_matrix b, 
    magma_tally3_z_matrix *x, magma_tally3_zopts *zopts,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_z_precondsetup(
    magma_tally3_z_matrix A, magma_tally3_z_matrix b, 
    magma_tally3_z_preconditioner *precond,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_z_applyprecond(
    magma_tally3_z_matrix A, magma_tally3_z_matrix b, 
    magma_tally3_z_matrix *x, magma_tally3_z_preconditioner *precond,
    magma_tally3_queue_t queue );


magma_tally3_int_t
magma_tally3_z_applyprecond_left(
    magma_tally3_z_matrix A, magma_tally3_z_matrix b, 
    magma_tally3_z_matrix *x, magma_tally3_z_preconditioner *precond,
    magma_tally3_queue_t queue );


magma_tally3_int_t
magma_tally3_z_applyprecond_right(
    magma_tally3_z_matrix A, magma_tally3_z_matrix b, 
    magma_tally3_z_matrix *x, magma_tally3_z_preconditioner *precond,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_z_initP2P(
    magma_tally3_int_t *bandwidth_benchmark,
    magma_tally3_int_t *num_gpus,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zcompact(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3DoubleComplex_ptr dA, magma_tally3_int_t ldda,
    double *dnorms, double tol, 
    magma_tally3_int_t *activeMask, magma_tally3_int_t *cBlockSize,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zcompactActive(
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3DoubleComplex_ptr dA, magma_tally3_int_t ldda, 
    magma_tally3_int_t *active,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zmlumerge(    
    magma_tally3_z_matrix L, 
    magma_tally3_z_matrix U,
    magma_tally3_z_matrix *A, 
    magma_tally3_queue_t queue );


/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally3_SPARSE BLAS function definitions
*/
magma_tally3_int_t 
magma_tally3_zgecsrmv(
    magma_tally3_trans_t transA,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3DoubleComplex alpha,
    magma_tally3DoubleComplex_ptr dval,
    magma_tally3Index_ptr drowptr,
    magma_tally3Index_ptr dcolind,
    magma_tally3DoubleComplex_ptr dx,
    magma_tally3DoubleComplex beta,
    magma_tally3DoubleComplex_ptr dy,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_zgecsrmv_shift(
    magma_tally3_trans_t transA,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3DoubleComplex alpha,
    magma_tally3DoubleComplex lambda,
    magma_tally3DoubleComplex_ptr dval,
    magma_tally3Index_ptr drowptr,
    magma_tally3Index_ptr dcolind,
    magma_tally3DoubleComplex_ptr dx,
    magma_tally3DoubleComplex beta,
    int offset,
    int blocksize,
    magma_tally3Index_ptr dadd_rows,
    magma_tally3DoubleComplex_ptr dy,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_zmgecsrmv(
    magma_tally3_trans_t transA,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3_int_t num_vecs,
    magma_tally3DoubleComplex alpha,
    magma_tally3DoubleComplex_ptr dval,
    magma_tally3Index_ptr drowptr,
    magma_tally3Index_ptr dcolind,
    magma_tally3DoubleComplex_ptr dx,
    magma_tally3DoubleComplex beta,
    magma_tally3DoubleComplex_ptr dy,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_zgeellmv(
    magma_tally3_trans_t transA,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3_int_t nnz_per_row,
    magma_tally3DoubleComplex alpha,
    magma_tally3DoubleComplex_ptr dval,
    magma_tally3Index_ptr dcolind,
    magma_tally3DoubleComplex_ptr dx,
    magma_tally3DoubleComplex beta,
    magma_tally3DoubleComplex_ptr dy,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_zgeellmv_shift(
    magma_tally3_trans_t transA,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3_int_t nnz_per_row,
    magma_tally3DoubleComplex alpha,
    magma_tally3DoubleComplex lambda,
    magma_tally3DoubleComplex_ptr dval,
    magma_tally3Index_ptr dcolind,
    magma_tally3DoubleComplex_ptr dx,
    magma_tally3DoubleComplex beta,
    int offset,
    int blocksize,
    magma_tally3Index_ptr dadd_rows,
    magma_tally3DoubleComplex_ptr dy,
    magma_tally3_queue_t queue );


magma_tally3_int_t 
magma_tally3_zmgeellmv(
    magma_tally3_trans_t transA,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3_int_t num_vecs,
    magma_tally3_int_t nnz_per_row,
    magma_tally3DoubleComplex alpha,
    magma_tally3DoubleComplex_ptr dval,
    magma_tally3Index_ptr dcolind,
    magma_tally3DoubleComplex_ptr dx,
    magma_tally3DoubleComplex beta,
    magma_tally3DoubleComplex_ptr dy,
    magma_tally3_queue_t queue );


magma_tally3_int_t 
magma_tally3_zgeelltmv(
    magma_tally3_trans_t transA,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3_int_t nnz_per_row,
    magma_tally3DoubleComplex alpha,
    magma_tally3DoubleComplex_ptr dval,
    magma_tally3Index_ptr dcolind,
    magma_tally3DoubleComplex_ptr dx,
    magma_tally3DoubleComplex beta,
    magma_tally3DoubleComplex_ptr dy,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_zgeelltmv_shift(
    magma_tally3_trans_t transA,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3_int_t nnz_per_row,
    magma_tally3DoubleComplex alpha,
    magma_tally3DoubleComplex lambda,
    magma_tally3DoubleComplex_ptr dval,
    magma_tally3Index_ptr dcolind,
    magma_tally3DoubleComplex_ptr dx,
    magma_tally3DoubleComplex beta,
    int offset,
    int blocksize,
    magma_tally3Index_ptr dadd_rows,
    magma_tally3DoubleComplex_ptr dy,
    magma_tally3_queue_t queue );


magma_tally3_int_t 
magma_tally3_zmgeelltmv(
    magma_tally3_trans_t transA,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3_int_t num_vecs,
    magma_tally3_int_t nnz_per_row,
    magma_tally3DoubleComplex alpha,
    magma_tally3DoubleComplex_ptr dval,
    magma_tally3Index_ptr dcolind,
    magma_tally3DoubleComplex_ptr dx,
    magma_tally3DoubleComplex beta,
    magma_tally3DoubleComplex_ptr dy,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_zgeellrtmv(
    magma_tally3_trans_t transA,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3_int_t nnz_per_row,
    magma_tally3DoubleComplex alpha,
    magma_tally3DoubleComplex_ptr dval,
    magma_tally3Index_ptr dcolind,
    magma_tally3Index_ptr drowlength,
    magma_tally3DoubleComplex_ptr dx,
    magma_tally3DoubleComplex beta,
    magma_tally3DoubleComplex_ptr dy,
    magma_tally3_int_t num_threads,
    magma_tally3_int_t threads_per_row,
    magma_tally3_queue_t queue );

magma_tally3_int_t 
magma_tally3_zgesellcmv(
    magma_tally3_trans_t transA,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3_int_t blocksize,
    magma_tally3_int_t slices,
    magma_tally3_int_t alignment,
    magma_tally3DoubleComplex alpha,
    magma_tally3DoubleComplex_ptr dval,
    magma_tally3Index_ptr dcolind,
    magma_tally3Index_ptr drowptr,
    magma_tally3DoubleComplex_ptr dx,
    magma_tally3DoubleComplex beta,
    magma_tally3DoubleComplex_ptr dy,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zgesellpmv(
    magma_tally3_trans_t transA,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3_int_t blocksize,
    magma_tally3_int_t slices,
    magma_tally3_int_t alignment,
    magma_tally3DoubleComplex alpha,
    magma_tally3DoubleComplex_ptr dval,
    magma_tally3Index_ptr dcolind,
    magma_tally3Index_ptr drowptr,
    magma_tally3DoubleComplex_ptr dx,
    magma_tally3DoubleComplex beta,
    magma_tally3DoubleComplex_ptr dy,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zmgesellpmv(
    magma_tally3_trans_t transA,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3_int_t num_vecs,
    magma_tally3_int_t blocksize,
    magma_tally3_int_t slices,
    magma_tally3_int_t alignment,
    magma_tally3DoubleComplex alpha,
    magma_tally3DoubleComplex_ptr dval,
    magma_tally3Index_ptr dcolind,
    magma_tally3Index_ptr drowptr,
    magma_tally3DoubleComplex_ptr dx,
    magma_tally3DoubleComplex beta,
    magma_tally3DoubleComplex_ptr dy,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zmgesellpmv_blocked(
    magma_tally3_trans_t transA,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3_int_t num_vecs,
    magma_tally3_int_t blocksize,
    magma_tally3_int_t slices,
    magma_tally3_int_t alignment,
    magma_tally3DoubleComplex alpha,
    magma_tally3DoubleComplex_ptr dval,
    magma_tally3Index_ptr dcolind,
    magma_tally3Index_ptr drowptr,
    magma_tally3DoubleComplex_ptr dx,
    magma_tally3DoubleComplex beta,
    magma_tally3DoubleComplex_ptr dy,
    magma_tally3_queue_t queue );


magma_tally3_int_t
magma_tally3_zmergedgs(
    magma_tally3_int_t n, 
    magma_tally3_int_t ldh,
    magma_tally3_int_t k, 
    magma_tally3DoubleComplex_ptr dv, 
    magma_tally3DoubleComplex_ptr dr,
    magma_tally3DoubleComplex_ptr dskp,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zcopyscale(    
    int n, 
    int k,
    magma_tally3DoubleComplex_ptr dr, 
    magma_tally3DoubleComplex_ptr dv,
    magma_tally3DoubleComplex_ptr dskp,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_dznrm2scale(    
    int m, 
    magma_tally3DoubleComplex_ptr dr,    
    int lddr, 
    magma_tally3DoubleComplex *drnorm,
    magma_tally3_queue_t queue );


magma_tally3_int_t
magma_tally3_zjacobisetup_vector_gpu(
    int num_rows, 
    magma_tally3_z_matrix b, 
    magma_tally3_z_matrix d, 
    magma_tally3_z_matrix c,
    magma_tally3_z_matrix *x,
    magma_tally3_queue_t queue );


magma_tally3_int_t
magma_tally3_zjacobi_diagscal(    
    int num_rows, 
    magma_tally3_z_matrix d, 
    magma_tally3_z_matrix b, 
    magma_tally3_z_matrix *c,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zjacobiupdate(
    magma_tally3_z_matrix t, 
    magma_tally3_z_matrix b, 
    magma_tally3_z_matrix d, 
    magma_tally3_z_matrix *x,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zjacobispmvupdate(
    magma_tally3_int_t maxiter,
    magma_tally3_z_matrix A, 
    magma_tally3_z_matrix t, 
    magma_tally3_z_matrix b, 
    magma_tally3_z_matrix d, 
    magma_tally3_z_matrix *x,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zjacobispmvupdate_bw(
    magma_tally3_int_t maxiter,
    magma_tally3_z_matrix A, 
    magma_tally3_z_matrix t, 
    magma_tally3_z_matrix b, 
    magma_tally3_z_matrix d, 
    magma_tally3_z_matrix *x,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zjacobispmvupdateselect(
    magma_tally3_int_t maxiter,
    magma_tally3_int_t num_updates,
    magma_tally3_index_t *indices,
    magma_tally3_z_matrix A,
    magma_tally3_z_matrix t, 
    magma_tally3_z_matrix b, 
    magma_tally3_z_matrix d, 
    magma_tally3_z_matrix tmp, 
    magma_tally3_z_matrix *x,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zjacobisetup_diagscal(
    magma_tally3_z_matrix A, magma_tally3_z_matrix *d,
    magma_tally3_queue_t queue );


magma_tally3_int_t
magma_tally3_zbicgmerge1(    
    int n, 
    magma_tally3DoubleComplex_ptr dskp,
    magma_tally3DoubleComplex_ptr dv, 
    magma_tally3DoubleComplex_ptr dr, 
    magma_tally3DoubleComplex_ptr dp,
    magma_tally3_queue_t queue );


magma_tally3_int_t
magma_tally3_zbicgmerge2(
    int n, 
    magma_tally3DoubleComplex_ptr dskp, 
    magma_tally3DoubleComplex_ptr dr,
    magma_tally3DoubleComplex_ptr dv, 
    magma_tally3DoubleComplex_ptr ds,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zbicgmerge3(
    int n, 
    magma_tally3DoubleComplex_ptr dskp, 
    magma_tally3DoubleComplex_ptr dp,
    magma_tally3DoubleComplex_ptr ds,
    magma_tally3DoubleComplex_ptr dt,
    magma_tally3DoubleComplex_ptr dx, 
    magma_tally3DoubleComplex_ptr dr,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zbicgmerge4(
    int type, 
    magma_tally3DoubleComplex_ptr dskp,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zcgmerge_spmv1( 
    magma_tally3_z_matrix A,
    magma_tally3DoubleComplex_ptr d1,
    magma_tally3DoubleComplex_ptr d2,
    magma_tally3DoubleComplex_ptr dd,
    magma_tally3DoubleComplex_ptr dz,
    magma_tally3DoubleComplex_ptr dskp,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zcgmerge_xrbeta( 
    int n,
    magma_tally3DoubleComplex_ptr d1,
    magma_tally3DoubleComplex_ptr d2,
    magma_tally3DoubleComplex_ptr dx,
    magma_tally3DoubleComplex_ptr dr,
    magma_tally3DoubleComplex_ptr dd,
    magma_tally3DoubleComplex_ptr dz, 
    magma_tally3DoubleComplex_ptr dskp,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zmdotc(
    magma_tally3_int_t n, 
    magma_tally3_int_t k, 
    magma_tally3DoubleComplex_ptr dv, 
    magma_tally3DoubleComplex_ptr dr,
    magma_tally3DoubleComplex_ptr dd1,
    magma_tally3DoubleComplex_ptr dd2,
    magma_tally3DoubleComplex_ptr dskp,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zgemvmdot(
    int n, 
    int k, 
    magma_tally3DoubleComplex_ptr dv, 
    magma_tally3DoubleComplex_ptr dr,
    magma_tally3DoubleComplex_ptr dd1,
    magma_tally3DoubleComplex_ptr dd2,
    magma_tally3DoubleComplex_ptr dskp,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zbicgmerge_spmv1( 
    magma_tally3_z_matrix A,
    magma_tally3DoubleComplex_ptr dd1,
    magma_tally3DoubleComplex_ptr dd2,
    magma_tally3DoubleComplex_ptr dp,
    magma_tally3DoubleComplex_ptr dr,
    magma_tally3DoubleComplex_ptr dv,
    magma_tally3DoubleComplex_ptr dskp,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zbicgmerge_spmv2( 
    magma_tally3_z_matrix A,
    magma_tally3DoubleComplex_ptr dd1,
    magma_tally3DoubleComplex_ptr dd2,
    magma_tally3DoubleComplex_ptr ds,
    magma_tally3DoubleComplex_ptr dt,
    magma_tally3DoubleComplex_ptr dskp,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zbicgmerge_xrbeta( 
    int n,
    magma_tally3DoubleComplex_ptr dd1,
    magma_tally3DoubleComplex_ptr dd2,
    magma_tally3DoubleComplex_ptr drr,
    magma_tally3DoubleComplex_ptr dr,
    magma_tally3DoubleComplex_ptr dp,
    magma_tally3DoubleComplex_ptr ds,
    magma_tally3DoubleComplex_ptr dt,
    magma_tally3DoubleComplex_ptr dx, 
    magma_tally3DoubleComplex_ptr dskp,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zbcsrswp(
    magma_tally3_int_t n,
    magma_tally3_int_t size_b, 
    magma_tally3_int_t *ipiv,
    magma_tally3DoubleComplex_ptr dx,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zbcsrtrsv(
    magma_tally3_uplo_t uplo,
    magma_tally3_int_t r_blocks,
    magma_tally3_int_t c_blocks,
    magma_tally3_int_t size_b, 
    magma_tally3DoubleComplex_ptr dA,
    magma_tally3_index_t *blockinfo, 
    magma_tally3DoubleComplex_ptr dx,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zbcsrvalcpy(
    magma_tally3_int_t size_b, 
    magma_tally3_int_t num_blocks, 
    magma_tally3_int_t num_zero_blocks, 
    magma_tally3DoubleComplex_ptr *dAval, 
    magma_tally3DoubleComplex_ptr *dBval,
    magma_tally3DoubleComplex_ptr *dBval2,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zbcsrluegemm(
    magma_tally3_int_t size_b, 
    magma_tally3_int_t num_block_rows,
    magma_tally3_int_t kblocks,
    magma_tally3DoubleComplex_ptr *dA, 
    magma_tally3DoubleComplex_ptr *dB, 
    magma_tally3DoubleComplex_ptr *dC,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zbcsrlupivloc(
    magma_tally3_int_t size_b, 
    magma_tally3_int_t kblocks,
    magma_tally3DoubleComplex_ptr *dA, 
    magma_tally3_int_t *ipiv,
    magma_tally3_queue_t queue );

magma_tally3_int_t
magma_tally3_zbcsrblockinfo5(
    magma_tally3_int_t lustep,
    magma_tally3_int_t num_blocks, 
    magma_tally3_int_t c_blocks, 
    magma_tally3_int_t size_b,
    magma_tally3_index_t *blockinfo,
    magma_tally3DoubleComplex_ptr dval,
    magma_tally3DoubleComplex_ptr *AII,
    magma_tally3_queue_t queue );

 
#ifdef __cplusplus
}
#endif

#undef PRECISION_z
#endif /* MAGMA_tally3SPARSE_Z_H */
