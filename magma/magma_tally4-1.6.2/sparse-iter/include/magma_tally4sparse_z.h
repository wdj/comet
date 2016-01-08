/*
 -- MAGMA_tally4 (version 1.1) --
 Univ. of Tennessee, Knoxville
 Univ. of California, Berkeley
 Univ. of Colorado, Denver
 @date

 @precisions normal z -> s d c
 @author Hartwig Anzt
*/

#ifndef MAGMA_tally4SPARSE_Z_H
#define MAGMA_tally4SPARSE_Z_H

#include "magma_tally4_types.h"
#include "magma_tally4sparse_types.h"

#define PRECISION_z


#ifdef __cplusplus
extern "C" {
#endif


/* ////////////////////////////////////////////////////////////////////////////
 -- For backwards compatability, map old (1.6.1) to new (1.6.2) function names
*/

#define magma_tally4_z_mtranspose  magma_tally4_zmtranspose
#define magma_tally4_z_mtransfer   magma_tally4_zmtransfer
#define magma_tally4_z_vtransfer   magma_tally4_zmtransfer
#define magma_tally4_z_mconvert    magma_tally4_zmconvert
#define magma_tally4_z_vinit       magma_tally4_zvinit
#define magma_tally4_z_vvisu       magma_tally4_zprint_vector
#define magma_tally4_z_vread       magma_tally4_zvread
#define magma_tally4_z_vspread     magma_tally4_zvspread
#define magma_tally4_z_mvisu       magma_tally4_zprint_matrix
#define magma_tally4_z_mfree       magma_tally4_zmfree
#define magma_tally4_z_vfree       magma_tally4_zmfree
#define write_z_csr_mtx     magma_tally4_zwrite_csr_mtx
#define write_z_csrtomtx    magma_tally4_zwrite_csrtomtx
#define print_z_csr         magma_tally4_zprint_csr_mtx


/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally4_SPARSE Auxiliary functions
*/


magma_tally4_int_t
magma_tally4_zparse_opts( 
    int argc, 
    char** argv, 
    magma_tally4_zopts *opts, 
    int *matrices, 
    magma_tally4_queue_t queue );

magma_tally4_int_t 
read_z_csr_from_binary( 
    magma_tally4_int_t* n_row, 
    magma_tally4_int_t* n_col, 
    magma_tally4_int_t* nnz, 
    magma_tally4DoubleComplex **val, 
    magma_tally4_index_t **row, 
    magma_tally4_index_t **col,
    const char * filename,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
read_z_csr_from_mtx( 
    magma_tally4_storage_t *type, 
    magma_tally4_location_t *location,
    magma_tally4_int_t* n_row, 
    magma_tally4_int_t* n_col, 
    magma_tally4_int_t* nnz, 
    magma_tally4DoubleComplex **val, 
    magma_tally4_index_t **row, 
    magma_tally4_index_t **col, 
    const char *filename,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_z_csr_mtx( 
    magma_tally4_z_matrix *A, 
    const char *filename,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_zcsrset( 
    magma_tally4_int_t m, 
    magma_tally4_int_t n, 
    magma_tally4_index_t *row, 
    magma_tally4_index_t *col, 
    magma_tally4DoubleComplex *val,
    magma_tally4_z_matrix *A,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_zcsrget( 
    magma_tally4_z_matrix A,
    magma_tally4_int_t *m, 
    magma_tally4_int_t *n, 
    magma_tally4_index_t **row, 
    magma_tally4_index_t **col, 
    magma_tally4DoubleComplex **val,
    magma_tally4_queue_t queue );


magma_tally4_int_t 
magma_tally4_zvset( 
    magma_tally4_int_t m, 
    magma_tally4_int_t n, 
    magma_tally4DoubleComplex *val,
    magma_tally4_z_matrix *v,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_zvget( 
    magma_tally4_z_matrix v,
    magma_tally4_int_t *m, 
    magma_tally4_int_t *n, 
    magma_tally4DoubleComplex **val,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_zvset_dev( 
    magma_tally4_int_t m, 
    magma_tally4_int_t n, 
    magma_tally4DoubleComplex_ptr val,
    magma_tally4_z_matrix *v,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_zvget_dev( 
    magma_tally4_z_matrix v,
    magma_tally4_int_t *m, 
    magma_tally4_int_t *n, 
    magma_tally4DoubleComplex_ptr *val,
    magma_tally4_queue_t queue );


magma_tally4_int_t 
magma_tally4_z_csr_mtxsymm( 
    magma_tally4_z_matrix *A, 
    const char *filename,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_z_csr_compressor( 
    magma_tally4DoubleComplex ** val, 
    magma_tally4_index_t ** row, 
    magma_tally4_index_t ** col, 
    magma_tally4DoubleComplex ** valn, 
    magma_tally4_index_t ** rown, 
    magma_tally4_index_t ** coln, 
    magma_tally4_int_t *n,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zmcsrcompressor( 
    magma_tally4_z_matrix *A,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zmcsrcompressor_gpu( 
    magma_tally4_z_matrix *A,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zvtranspose( 
    magma_tally4_z_matrix x,
    magma_tally4_z_matrix *y,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_z_cucsrtranspose( 
    magma_tally4_z_matrix A, 
    magma_tally4_z_matrix *B,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
z_transpose_csr( 
    magma_tally4_int_t n_rows, 
    magma_tally4_int_t n_cols, 
    magma_tally4_int_t nnz,
    magma_tally4DoubleComplex *val, 
    magma_tally4_index_t *row, 
    magma_tally4_index_t *col, 
    magma_tally4_int_t *new_n_rows, 
    magma_tally4_int_t *new_n_cols, 
    magma_tally4_int_t *new_nnz, 
    magma_tally4DoubleComplex **new_val, 
    magma_tally4_index_t **new_row, 
    magma_tally4_index_t **new_col,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zcsrsplit( 
    magma_tally4_int_t bsize,
    magma_tally4_z_matrix A,
    magma_tally4_z_matrix *D,
    magma_tally4_z_matrix *R,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zmscale( 
    magma_tally4_z_matrix *A, 
    magma_tally4_scale_t scaling,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_zmdiff( 
    magma_tally4_z_matrix A, 
    magma_tally4_z_matrix B, 
 real_Double_t *res,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zmdiagadd( 
    magma_tally4_z_matrix *A, 
    magma_tally4DoubleComplex add,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_zmsort( 
    magma_tally4_z_matrix *A,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zindexsort(
    magma_tally4_index_t *x, 
    magma_tally4_int_t first,
    magma_tally4_int_t last,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zdomainoverlap(
    magma_tally4_index_t num_rows,
    magma_tally4_index_t *num_indices,
    magma_tally4_index_t *rowptr,
    magma_tally4_index_t *colidx,
    magma_tally4_index_t *x,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zsymbilu( 
    magma_tally4_z_matrix *A, 
    magma_tally4_int_t levels,
    magma_tally4_z_matrix *L,
    magma_tally4_z_matrix *U,
    magma_tally4_queue_t queue );


magma_tally4_int_t 
magma_tally4_zwrite_csr_mtx( 
    magma_tally4_z_matrix A,
    magma_tally4_order_t MajorType,
 const char *filename,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_zwrite_csrtomtx( 
    magma_tally4_z_matrix A,
    const char *filename,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_zprint_csr( 
    magma_tally4_int_t n_row, 
    magma_tally4_int_t n_col, 
    magma_tally4_int_t nnz, 
    magma_tally4DoubleComplex **val, 
    magma_tally4_index_t **row, 
    magma_tally4_index_t **col,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_zprint_csr_mtx( 
    magma_tally4_int_t n_row, 
    magma_tally4_int_t n_col, 
    magma_tally4_int_t nnz, 
    magma_tally4DoubleComplex **val, 
    magma_tally4_index_t **row, 
    magma_tally4_index_t **col, 
    magma_tally4_order_t MajorType,
    magma_tally4_queue_t queue );


magma_tally4_int_t 
magma_tally4_zmtranspose(
    magma_tally4_z_matrix A, 
    magma_tally4_z_matrix *B,
    magma_tally4_queue_t queue );


magma_tally4_int_t 
magma_tally4_zmtransfer(
    magma_tally4_z_matrix A, 
    magma_tally4_z_matrix *B, 
    magma_tally4_location_t src, 
    magma_tally4_location_t dst,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_zmconvert(
    magma_tally4_z_matrix A, 
    magma_tally4_z_matrix *B, 
    magma_tally4_storage_t old_format, 
    magma_tally4_storage_t new_format,
    magma_tally4_queue_t queue );


magma_tally4_int_t
magma_tally4_zvinit(
    magma_tally4_z_matrix *x, 
    magma_tally4_location_t memory_location,
    magma_tally4_int_t num_rows, 
    magma_tally4_int_t num_cols,
    magma_tally4DoubleComplex values,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zprint_vector(
    magma_tally4_z_matrix x, 
    magma_tally4_int_t offset, 
    magma_tally4_int_t displaylength,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zvread(
    magma_tally4_z_matrix *x, 
    magma_tally4_int_t length,
    char * filename,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zvspread(
    magma_tally4_z_matrix *x, 
    const char * filename,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zprint_matrix(
    magma_tally4_z_matrix A,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_zdiameter(
    magma_tally4_z_matrix *A,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_zrowentries(
    magma_tally4_z_matrix *A,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zmfree(
    magma_tally4_z_matrix *A,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zresidual(
    magma_tally4_z_matrix A, 
    magma_tally4_z_matrix b, 
    magma_tally4_z_matrix x, 
    double *res,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zresidualvec(
    magma_tally4_z_matrix A,
    magma_tally4_z_matrix b,
    magma_tally4_z_matrix x,
    magma_tally4_z_matrix *r,
    double *res,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zmgenerator(
    magma_tally4_int_t n,
    magma_tally4_int_t offdiags,
    magma_tally4_index_t *diag_offset,
    magma_tally4DoubleComplex *diag_vals,
    magma_tally4_z_matrix *A,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zm_27stencil(
    magma_tally4_int_t n,
    magma_tally4_z_matrix *A,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zm_5stencil(
    magma_tally4_int_t n,
    magma_tally4_z_matrix *A,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zsolverinfo(
    magma_tally4_z_solver_par *solver_par, 
    magma_tally4_z_preconditioner *precond_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zsolverinfo_init(
    magma_tally4_z_solver_par *solver_par, 
    magma_tally4_z_preconditioner *precond,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zeigensolverinfo_init(
    magma_tally4_z_solver_par *solver_par,
    magma_tally4_queue_t queue );


magma_tally4_int_t
magma_tally4_zsolverinfo_free(
    magma_tally4_z_solver_par *solver_par, 
    magma_tally4_z_preconditioner *precond,
    magma_tally4_queue_t queue );



/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally4_SPARSE iterative incomplete factorizations
*/


magma_tally4_int_t
magma_tally4_ziterilusetup( 
    magma_tally4_z_matrix A, 
    magma_tally4_z_matrix b,                                 
    magma_tally4_z_preconditioner *precond,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zitericsetup( 
    magma_tally4_z_matrix A, 
    magma_tally4_z_matrix b, 
    magma_tally4_z_preconditioner *precond,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zitericupdate( 
    magma_tally4_z_matrix A, 
    magma_tally4_z_preconditioner *precond, 
    magma_tally4_int_t updates,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zapplyiteric_l( 
    magma_tally4_z_matrix b, 
    magma_tally4_z_matrix *x, 
    magma_tally4_z_preconditioner *precond,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zapplyiteric_r( 
    magma_tally4_z_matrix b, 
    magma_tally4_z_matrix *x, 
    magma_tally4_z_preconditioner *precond,
    magma_tally4_queue_t queue );


magma_tally4_int_t
magma_tally4_ziterilu_csr( 
    magma_tally4_z_matrix A,
    magma_tally4_z_matrix L,
    magma_tally4_z_matrix U,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_ziteric_csr( 
    magma_tally4_z_matrix A,
    magma_tally4_z_matrix A_CSR,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_zfrobenius( 
    magma_tally4_z_matrix A, 
    magma_tally4_z_matrix B, 
    real_Double_t *res,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_znonlinres(   
    magma_tally4_z_matrix A, 
    magma_tally4_z_matrix L,
    magma_tally4_z_matrix U, 
    magma_tally4_z_matrix *LU, 
    real_Double_t *res,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_zilures(   
    magma_tally4_z_matrix A, 
    magma_tally4_z_matrix L,
    magma_tally4_z_matrix U, 
    magma_tally4_z_matrix *LU, 
    real_Double_t *res,
    real_Double_t *nonlinres,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_zicres(       
    magma_tally4_z_matrix A, 
    magma_tally4_z_matrix C,
    magma_tally4_z_matrix CT, 
    magma_tally4_z_matrix *LU, 
    real_Double_t *res,
    real_Double_t *nonlinres,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_zinitguess( 
    magma_tally4_z_matrix A, 
    magma_tally4_z_matrix *L, 
    magma_tally4_z_matrix *U,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_zinitrecursiveLU( 
    magma_tally4_z_matrix A, 
    magma_tally4_z_matrix *B,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_zmLdiagadd( 
    magma_tally4_z_matrix *L,
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
magma_tally4_zcg(
    magma_tally4_z_matrix A, magma_tally4_z_matrix b, 
    magma_tally4_z_matrix *x, magma_tally4_z_solver_par *solver_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_zcg_res(
    magma_tally4_z_matrix A, magma_tally4_z_matrix b, 
    magma_tally4_z_matrix *x, magma_tally4_z_solver_par *solver_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_zcg_merge(
    magma_tally4_z_matrix A, magma_tally4_z_matrix b, 
    magma_tally4_z_matrix *x, magma_tally4_z_solver_par *solver_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_zgmres(
    magma_tally4_z_matrix A, magma_tally4_z_matrix b, 
    magma_tally4_z_matrix *x, magma_tally4_z_solver_par *solver_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zbicgstab(
    magma_tally4_z_matrix A, magma_tally4_z_matrix b, magma_tally4_z_matrix *x, 
    magma_tally4_z_solver_par *solver_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zbicgstab_merge(
    magma_tally4_z_matrix A, magma_tally4_z_matrix b, 
    magma_tally4_z_matrix *x, magma_tally4_z_solver_par *solver_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zbicgstab_merge2(
    magma_tally4_z_matrix A, magma_tally4_z_matrix b, 
    magma_tally4_z_matrix *x, magma_tally4_z_solver_par *solver_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zpcg(
    magma_tally4_z_matrix A, magma_tally4_z_matrix b, 
    magma_tally4_z_matrix *x, magma_tally4_z_solver_par *solver_par, 
    magma_tally4_z_preconditioner *precond_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zbpcg(
    magma_tally4_z_matrix A, magma_tally4_z_matrix b, 
    magma_tally4_z_matrix *x, magma_tally4_z_solver_par *solver_par, 
    magma_tally4_z_preconditioner *precond_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zpbicgstab(
    magma_tally4_z_matrix A, magma_tally4_z_matrix b, 
    magma_tally4_z_matrix *x, magma_tally4_z_solver_par *solver_par, 
    magma_tally4_z_preconditioner *precond_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zpgmres(
    magma_tally4_z_matrix A, magma_tally4_z_matrix b, 
    magma_tally4_z_matrix *x, magma_tally4_z_solver_par *solver_par, 
    magma_tally4_z_preconditioner *precond_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zfgmres(
    magma_tally4_z_matrix A, magma_tally4_z_matrix b, 
    magma_tally4_z_matrix *x, magma_tally4_z_solver_par *solver_par, 
    magma_tally4_z_preconditioner *precond_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zjacobi(
    magma_tally4_z_matrix A, magma_tally4_z_matrix b, 
    magma_tally4_z_matrix *x, magma_tally4_z_solver_par *solver_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zjacobidomainoverlap(
    magma_tally4_z_matrix A, 
    magma_tally4_z_matrix b, 
    magma_tally4_z_matrix *x,  
    magma_tally4_z_solver_par *solver_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zbaiter(
    magma_tally4_z_matrix A, magma_tally4_z_matrix b, 
    magma_tally4_z_matrix *x, magma_tally4_z_solver_par *solver_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_ziterref(
    magma_tally4_z_matrix A, magma_tally4_z_matrix b, 
    magma_tally4_z_matrix *x, magma_tally4_z_solver_par *solver_par, 
    magma_tally4_z_preconditioner *precond_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zilu(
    magma_tally4_z_matrix A, magma_tally4_z_matrix b, 
    magma_tally4_z_matrix *x, magma_tally4_z_solver_par *solver_par, 
    magma_tally4_int_t *ipiv,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zbcsrlu(
    magma_tally4_z_matrix A, magma_tally4_z_matrix b, 
    magma_tally4_z_matrix *x, magma_tally4_z_solver_par *solver_par,
    magma_tally4_queue_t queue );


magma_tally4_int_t
magma_tally4_zbcsrlutrf(
    magma_tally4_z_matrix A, 
    magma_tally4_z_matrix *M,
    magma_tally4_int_t *ipiv, 
    magma_tally4_int_t version,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zbcsrlusv(
    magma_tally4_z_matrix A, magma_tally4_z_matrix b, 
    magma_tally4_z_matrix *x, magma_tally4_z_solver_par *solver_par, 
    magma_tally4_int_t *ipiv,
    magma_tally4_queue_t queue );



magma_tally4_int_t
magma_tally4_zilucg(
    magma_tally4_z_matrix A, magma_tally4_z_matrix b, 
    magma_tally4_z_matrix *x, magma_tally4_z_solver_par *solver_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zilugmres(
    magma_tally4_z_matrix A, magma_tally4_z_matrix b, 
    magma_tally4_z_matrix *x, magma_tally4_z_solver_par *solver_par,
    magma_tally4_queue_t queue ); 


magma_tally4_int_t
magma_tally4_zlobpcg_shift(
    magma_tally4_int_t num_rows,
    magma_tally4_int_t num_vecs, 
    magma_tally4_int_t shift,
    magma_tally4DoubleComplex_ptr x,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zlobpcg_res(
    magma_tally4_int_t num_rows,
    magma_tally4_int_t num_vecs, 
    double *evalues, 
    magma_tally4DoubleComplex_ptr X,
    magma_tally4DoubleComplex_ptr R, 
    double *res,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zlobpcg_maxpy(
    magma_tally4_int_t num_rows,
    magma_tally4_int_t num_vecs, 
    magma_tally4DoubleComplex_ptr X,
    magma_tally4DoubleComplex_ptr Y,
    magma_tally4_queue_t queue );


/*/////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally4_SPARSE eigensolvers (Data on GPU)
*/
magma_tally4_int_t
magma_tally4_zlobpcg(
    magma_tally4_z_matrix A, 
    magma_tally4_z_solver_par *solver_par,
    magma_tally4_z_preconditioner *precond_par, 
    magma_tally4_queue_t queue );




/*/////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally4_SPARSE preconditioners (Data on GPU)
*/

magma_tally4_int_t
magma_tally4_zjacobisetup(
    magma_tally4_z_matrix A, 
    magma_tally4_z_matrix b, 
    magma_tally4_z_matrix *M, 
    magma_tally4_z_matrix *c,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zjacobisetup_matrix(
    magma_tally4_z_matrix A, 
    magma_tally4_z_matrix *M, 
    magma_tally4_z_matrix *d,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zjacobisetup_vector(
    magma_tally4_z_matrix b, 
    magma_tally4_z_matrix d, 
    magma_tally4_z_matrix *c,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zjacobiiter(
    magma_tally4_z_matrix M, 
    magma_tally4_z_matrix c, 
    magma_tally4_z_matrix *x, 
    magma_tally4_z_solver_par *solver_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zjacobiiter_precond( 
    magma_tally4_z_matrix M, 
    magma_tally4_z_matrix *x, 
    magma_tally4_z_solver_par *solver_par, 
    magma_tally4_z_preconditioner *precond,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zjacobiiter_sys(
    magma_tally4_z_matrix A, 
    magma_tally4_z_matrix b, 
    magma_tally4_z_matrix d, 
    magma_tally4_z_matrix t, 
    magma_tally4_z_matrix *x,  
    magma_tally4_z_solver_par *solver_par,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zpastixsetup(
    magma_tally4_z_matrix A, magma_tally4_z_matrix b,
    magma_tally4_z_preconditioner *precond,
    magma_tally4_queue_t queue );


magma_tally4_int_t
magma_tally4_zapplypastix(
    magma_tally4_z_matrix b, magma_tally4_z_matrix *x, 
    magma_tally4_z_preconditioner *precond,
    magma_tally4_queue_t queue );


// custom preconditioner
magma_tally4_int_t
magma_tally4_zapplycustomprecond_l(
    magma_tally4_z_matrix b, 
    magma_tally4_z_matrix *x, 
    magma_tally4_z_preconditioner *precond,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zapplycustomprecond_r(
    magma_tally4_z_matrix b, 
    magma_tally4_z_matrix *x, 
    magma_tally4_z_preconditioner *precond,
    magma_tally4_queue_t queue );


// CUSPARSE preconditioner

magma_tally4_int_t
magma_tally4_zcuilusetup(
    magma_tally4_z_matrix A, magma_tally4_z_preconditioner *precond,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zapplycuilu_l(
    magma_tally4_z_matrix b, magma_tally4_z_matrix *x, 
    magma_tally4_z_preconditioner *precond,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zapplycuilu_r(
    magma_tally4_z_matrix b, magma_tally4_z_matrix *x, 
    magma_tally4_z_preconditioner *precond,
    magma_tally4_queue_t queue );


magma_tally4_int_t
magma_tally4_zcuiccsetup(
    magma_tally4_z_matrix A, magma_tally4_z_preconditioner *precond,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zapplycuicc_l(
    magma_tally4_z_matrix b, magma_tally4_z_matrix *x, 
    magma_tally4_z_preconditioner *precond,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zapplycuicc_r(
    magma_tally4_z_matrix b, magma_tally4_z_matrix *x, 
    magma_tally4_z_preconditioner *precond,
    magma_tally4_queue_t queue );


magma_tally4_int_t
magma_tally4_zcumilusetup(
    magma_tally4_z_matrix A, 
    magma_tally4_z_preconditioner *precond,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zcumilugeneratesolverinfo(
    magma_tally4_z_preconditioner *precond,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zapplycumilu_l(
    magma_tally4_z_matrix b, 
    magma_tally4_z_matrix *x, 
    magma_tally4_z_preconditioner *precond,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zapplycumilu_r(
    magma_tally4_z_matrix b, 
    magma_tally4_z_matrix *x, 
    magma_tally4_z_preconditioner *precond,
    magma_tally4_queue_t queue );


magma_tally4_int_t
magma_tally4_zcumiccsetup(
    magma_tally4_z_matrix A, 
    magma_tally4_z_preconditioner *precond,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zcumicgeneratesolverinfo(
    magma_tally4_z_preconditioner *precond,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zapplycumicc_l(
    magma_tally4_z_matrix b, 
    magma_tally4_z_matrix *x, 
    magma_tally4_z_preconditioner *precond,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zapplycumicc_r(
    magma_tally4_z_matrix b, 
    magma_tally4_z_matrix *x, 
    magma_tally4_z_preconditioner *precond,
    magma_tally4_queue_t queue );


// block-asynchronous iteration

magma_tally4_int_t
magma_tally4_zbajac_csr(
    magma_tally4_int_t localiters,
    magma_tally4_z_matrix D,
    magma_tally4_z_matrix R,
    magma_tally4_z_matrix b,
    magma_tally4_z_matrix *x,
    magma_tally4_queue_t queue );

/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally4_SPARSE utility function definitions
*/

magma_tally4_int_t
magma_tally4_z_spmv(
    magma_tally4DoubleComplex alpha, 
    magma_tally4_z_matrix A, 
    magma_tally4_z_matrix x, 
    magma_tally4DoubleComplex beta, 
    magma_tally4_z_matrix y,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zcustomspmv(
    magma_tally4DoubleComplex alpha, 
    magma_tally4_z_matrix x, 
    magma_tally4DoubleComplex beta, 
    magma_tally4_z_matrix y,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_z_spmv_shift(
    magma_tally4DoubleComplex alpha, 
    magma_tally4_z_matrix A, 
    magma_tally4DoubleComplex lambda,
    magma_tally4_z_matrix x, 
    magma_tally4DoubleComplex beta, 
    magma_tally4_int_t offset, 
    magma_tally4_int_t blocksize,
    magma_tally4Index_ptr dadd_vecs, 
    magma_tally4_z_matrix y,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zcuspmm(
    magma_tally4_z_matrix A, 
    magma_tally4_z_matrix B, 
    magma_tally4_z_matrix *AB,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_z_spmm(
    magma_tally4DoubleComplex alpha, 
    magma_tally4_z_matrix A,
    magma_tally4_z_matrix B,
    magma_tally4_z_matrix *C,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zsymbilu( 
    magma_tally4_z_matrix *A, 
    magma_tally4_int_t levels, 
    magma_tally4_z_matrix *L, 
    magma_tally4_z_matrix *U,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zcuspaxpy(
    magma_tally4DoubleComplex_ptr alpha, magma_tally4_z_matrix A, 
    magma_tally4DoubleComplex_ptr beta, magma_tally4_z_matrix B, 
    magma_tally4_z_matrix *AB,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_z_precond(
    magma_tally4_z_matrix A, 
    magma_tally4_z_matrix b, magma_tally4_z_matrix *x,
    magma_tally4_z_preconditioner *precond,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_z_solver(
    magma_tally4_z_matrix A, magma_tally4_z_matrix b, 
    magma_tally4_z_matrix *x, magma_tally4_zopts *zopts,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_z_precondsetup(
    magma_tally4_z_matrix A, magma_tally4_z_matrix b, 
    magma_tally4_z_preconditioner *precond,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_z_applyprecond(
    magma_tally4_z_matrix A, magma_tally4_z_matrix b, 
    magma_tally4_z_matrix *x, magma_tally4_z_preconditioner *precond,
    magma_tally4_queue_t queue );


magma_tally4_int_t
magma_tally4_z_applyprecond_left(
    magma_tally4_z_matrix A, magma_tally4_z_matrix b, 
    magma_tally4_z_matrix *x, magma_tally4_z_preconditioner *precond,
    magma_tally4_queue_t queue );


magma_tally4_int_t
magma_tally4_z_applyprecond_right(
    magma_tally4_z_matrix A, magma_tally4_z_matrix b, 
    magma_tally4_z_matrix *x, magma_tally4_z_preconditioner *precond,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_z_initP2P(
    magma_tally4_int_t *bandwidth_benchmark,
    magma_tally4_int_t *num_gpus,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zcompact(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda,
    double *dnorms, double tol, 
    magma_tally4_int_t *activeMask, magma_tally4_int_t *cBlockSize,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zcompactActive(
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex_ptr dA, magma_tally4_int_t ldda, 
    magma_tally4_int_t *active,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zmlumerge(    
    magma_tally4_z_matrix L, 
    magma_tally4_z_matrix U,
    magma_tally4_z_matrix *A, 
    magma_tally4_queue_t queue );


/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally4_SPARSE BLAS function definitions
*/
magma_tally4_int_t 
magma_tally4_zgecsrmv(
    magma_tally4_trans_t transA,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex alpha,
    magma_tally4DoubleComplex_ptr dval,
    magma_tally4Index_ptr drowptr,
    magma_tally4Index_ptr dcolind,
    magma_tally4DoubleComplex_ptr dx,
    magma_tally4DoubleComplex beta,
    magma_tally4DoubleComplex_ptr dy,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_zgecsrmv_shift(
    magma_tally4_trans_t transA,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4DoubleComplex alpha,
    magma_tally4DoubleComplex lambda,
    magma_tally4DoubleComplex_ptr dval,
    magma_tally4Index_ptr drowptr,
    magma_tally4Index_ptr dcolind,
    magma_tally4DoubleComplex_ptr dx,
    magma_tally4DoubleComplex beta,
    int offset,
    int blocksize,
    magma_tally4Index_ptr dadd_rows,
    magma_tally4DoubleComplex_ptr dy,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_zmgecsrmv(
    magma_tally4_trans_t transA,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4_int_t num_vecs,
    magma_tally4DoubleComplex alpha,
    magma_tally4DoubleComplex_ptr dval,
    magma_tally4Index_ptr drowptr,
    magma_tally4Index_ptr dcolind,
    magma_tally4DoubleComplex_ptr dx,
    magma_tally4DoubleComplex beta,
    magma_tally4DoubleComplex_ptr dy,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_zgeellmv(
    magma_tally4_trans_t transA,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4_int_t nnz_per_row,
    magma_tally4DoubleComplex alpha,
    magma_tally4DoubleComplex_ptr dval,
    magma_tally4Index_ptr dcolind,
    magma_tally4DoubleComplex_ptr dx,
    magma_tally4DoubleComplex beta,
    magma_tally4DoubleComplex_ptr dy,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_zgeellmv_shift(
    magma_tally4_trans_t transA,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4_int_t nnz_per_row,
    magma_tally4DoubleComplex alpha,
    magma_tally4DoubleComplex lambda,
    magma_tally4DoubleComplex_ptr dval,
    magma_tally4Index_ptr dcolind,
    magma_tally4DoubleComplex_ptr dx,
    magma_tally4DoubleComplex beta,
    int offset,
    int blocksize,
    magma_tally4Index_ptr dadd_rows,
    magma_tally4DoubleComplex_ptr dy,
    magma_tally4_queue_t queue );


magma_tally4_int_t 
magma_tally4_zmgeellmv(
    magma_tally4_trans_t transA,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4_int_t num_vecs,
    magma_tally4_int_t nnz_per_row,
    magma_tally4DoubleComplex alpha,
    magma_tally4DoubleComplex_ptr dval,
    magma_tally4Index_ptr dcolind,
    magma_tally4DoubleComplex_ptr dx,
    magma_tally4DoubleComplex beta,
    magma_tally4DoubleComplex_ptr dy,
    magma_tally4_queue_t queue );


magma_tally4_int_t 
magma_tally4_zgeelltmv(
    magma_tally4_trans_t transA,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4_int_t nnz_per_row,
    magma_tally4DoubleComplex alpha,
    magma_tally4DoubleComplex_ptr dval,
    magma_tally4Index_ptr dcolind,
    magma_tally4DoubleComplex_ptr dx,
    magma_tally4DoubleComplex beta,
    magma_tally4DoubleComplex_ptr dy,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_zgeelltmv_shift(
    magma_tally4_trans_t transA,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4_int_t nnz_per_row,
    magma_tally4DoubleComplex alpha,
    magma_tally4DoubleComplex lambda,
    magma_tally4DoubleComplex_ptr dval,
    magma_tally4Index_ptr dcolind,
    magma_tally4DoubleComplex_ptr dx,
    magma_tally4DoubleComplex beta,
    int offset,
    int blocksize,
    magma_tally4Index_ptr dadd_rows,
    magma_tally4DoubleComplex_ptr dy,
    magma_tally4_queue_t queue );


magma_tally4_int_t 
magma_tally4_zmgeelltmv(
    magma_tally4_trans_t transA,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4_int_t num_vecs,
    magma_tally4_int_t nnz_per_row,
    magma_tally4DoubleComplex alpha,
    magma_tally4DoubleComplex_ptr dval,
    magma_tally4Index_ptr dcolind,
    magma_tally4DoubleComplex_ptr dx,
    magma_tally4DoubleComplex beta,
    magma_tally4DoubleComplex_ptr dy,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_zgeellrtmv(
    magma_tally4_trans_t transA,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4_int_t nnz_per_row,
    magma_tally4DoubleComplex alpha,
    magma_tally4DoubleComplex_ptr dval,
    magma_tally4Index_ptr dcolind,
    magma_tally4Index_ptr drowlength,
    magma_tally4DoubleComplex_ptr dx,
    magma_tally4DoubleComplex beta,
    magma_tally4DoubleComplex_ptr dy,
    magma_tally4_int_t num_threads,
    magma_tally4_int_t threads_per_row,
    magma_tally4_queue_t queue );

magma_tally4_int_t 
magma_tally4_zgesellcmv(
    magma_tally4_trans_t transA,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4_int_t blocksize,
    magma_tally4_int_t slices,
    magma_tally4_int_t alignment,
    magma_tally4DoubleComplex alpha,
    magma_tally4DoubleComplex_ptr dval,
    magma_tally4Index_ptr dcolind,
    magma_tally4Index_ptr drowptr,
    magma_tally4DoubleComplex_ptr dx,
    magma_tally4DoubleComplex beta,
    magma_tally4DoubleComplex_ptr dy,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zgesellpmv(
    magma_tally4_trans_t transA,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4_int_t blocksize,
    magma_tally4_int_t slices,
    magma_tally4_int_t alignment,
    magma_tally4DoubleComplex alpha,
    magma_tally4DoubleComplex_ptr dval,
    magma_tally4Index_ptr dcolind,
    magma_tally4Index_ptr drowptr,
    magma_tally4DoubleComplex_ptr dx,
    magma_tally4DoubleComplex beta,
    magma_tally4DoubleComplex_ptr dy,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zmgesellpmv(
    magma_tally4_trans_t transA,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4_int_t num_vecs,
    magma_tally4_int_t blocksize,
    magma_tally4_int_t slices,
    magma_tally4_int_t alignment,
    magma_tally4DoubleComplex alpha,
    magma_tally4DoubleComplex_ptr dval,
    magma_tally4Index_ptr dcolind,
    magma_tally4Index_ptr drowptr,
    magma_tally4DoubleComplex_ptr dx,
    magma_tally4DoubleComplex beta,
    magma_tally4DoubleComplex_ptr dy,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zmgesellpmv_blocked(
    magma_tally4_trans_t transA,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4_int_t num_vecs,
    magma_tally4_int_t blocksize,
    magma_tally4_int_t slices,
    magma_tally4_int_t alignment,
    magma_tally4DoubleComplex alpha,
    magma_tally4DoubleComplex_ptr dval,
    magma_tally4Index_ptr dcolind,
    magma_tally4Index_ptr drowptr,
    magma_tally4DoubleComplex_ptr dx,
    magma_tally4DoubleComplex beta,
    magma_tally4DoubleComplex_ptr dy,
    magma_tally4_queue_t queue );


magma_tally4_int_t
magma_tally4_zmergedgs(
    magma_tally4_int_t n, 
    magma_tally4_int_t ldh,
    magma_tally4_int_t k, 
    magma_tally4DoubleComplex_ptr dv, 
    magma_tally4DoubleComplex_ptr dr,
    magma_tally4DoubleComplex_ptr dskp,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zcopyscale(    
    int n, 
    int k,
    magma_tally4DoubleComplex_ptr dr, 
    magma_tally4DoubleComplex_ptr dv,
    magma_tally4DoubleComplex_ptr dskp,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_dznrm2scale(    
    int m, 
    magma_tally4DoubleComplex_ptr dr,    
    int lddr, 
    magma_tally4DoubleComplex *drnorm,
    magma_tally4_queue_t queue );


magma_tally4_int_t
magma_tally4_zjacobisetup_vector_gpu(
    int num_rows, 
    magma_tally4_z_matrix b, 
    magma_tally4_z_matrix d, 
    magma_tally4_z_matrix c,
    magma_tally4_z_matrix *x,
    magma_tally4_queue_t queue );


magma_tally4_int_t
magma_tally4_zjacobi_diagscal(    
    int num_rows, 
    magma_tally4_z_matrix d, 
    magma_tally4_z_matrix b, 
    magma_tally4_z_matrix *c,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zjacobiupdate(
    magma_tally4_z_matrix t, 
    magma_tally4_z_matrix b, 
    magma_tally4_z_matrix d, 
    magma_tally4_z_matrix *x,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zjacobispmvupdate(
    magma_tally4_int_t maxiter,
    magma_tally4_z_matrix A, 
    magma_tally4_z_matrix t, 
    magma_tally4_z_matrix b, 
    magma_tally4_z_matrix d, 
    magma_tally4_z_matrix *x,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zjacobispmvupdate_bw(
    magma_tally4_int_t maxiter,
    magma_tally4_z_matrix A, 
    magma_tally4_z_matrix t, 
    magma_tally4_z_matrix b, 
    magma_tally4_z_matrix d, 
    magma_tally4_z_matrix *x,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zjacobispmvupdateselect(
    magma_tally4_int_t maxiter,
    magma_tally4_int_t num_updates,
    magma_tally4_index_t *indices,
    magma_tally4_z_matrix A,
    magma_tally4_z_matrix t, 
    magma_tally4_z_matrix b, 
    magma_tally4_z_matrix d, 
    magma_tally4_z_matrix tmp, 
    magma_tally4_z_matrix *x,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zjacobisetup_diagscal(
    magma_tally4_z_matrix A, magma_tally4_z_matrix *d,
    magma_tally4_queue_t queue );


magma_tally4_int_t
magma_tally4_zbicgmerge1(    
    int n, 
    magma_tally4DoubleComplex_ptr dskp,
    magma_tally4DoubleComplex_ptr dv, 
    magma_tally4DoubleComplex_ptr dr, 
    magma_tally4DoubleComplex_ptr dp,
    magma_tally4_queue_t queue );


magma_tally4_int_t
magma_tally4_zbicgmerge2(
    int n, 
    magma_tally4DoubleComplex_ptr dskp, 
    magma_tally4DoubleComplex_ptr dr,
    magma_tally4DoubleComplex_ptr dv, 
    magma_tally4DoubleComplex_ptr ds,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zbicgmerge3(
    int n, 
    magma_tally4DoubleComplex_ptr dskp, 
    magma_tally4DoubleComplex_ptr dp,
    magma_tally4DoubleComplex_ptr ds,
    magma_tally4DoubleComplex_ptr dt,
    magma_tally4DoubleComplex_ptr dx, 
    magma_tally4DoubleComplex_ptr dr,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zbicgmerge4(
    int type, 
    magma_tally4DoubleComplex_ptr dskp,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zcgmerge_spmv1( 
    magma_tally4_z_matrix A,
    magma_tally4DoubleComplex_ptr d1,
    magma_tally4DoubleComplex_ptr d2,
    magma_tally4DoubleComplex_ptr dd,
    magma_tally4DoubleComplex_ptr dz,
    magma_tally4DoubleComplex_ptr dskp,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zcgmerge_xrbeta( 
    int n,
    magma_tally4DoubleComplex_ptr d1,
    magma_tally4DoubleComplex_ptr d2,
    magma_tally4DoubleComplex_ptr dx,
    magma_tally4DoubleComplex_ptr dr,
    magma_tally4DoubleComplex_ptr dd,
    magma_tally4DoubleComplex_ptr dz, 
    magma_tally4DoubleComplex_ptr dskp,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zmdotc(
    magma_tally4_int_t n, 
    magma_tally4_int_t k, 
    magma_tally4DoubleComplex_ptr dv, 
    magma_tally4DoubleComplex_ptr dr,
    magma_tally4DoubleComplex_ptr dd1,
    magma_tally4DoubleComplex_ptr dd2,
    magma_tally4DoubleComplex_ptr dskp,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zgemvmdot(
    int n, 
    int k, 
    magma_tally4DoubleComplex_ptr dv, 
    magma_tally4DoubleComplex_ptr dr,
    magma_tally4DoubleComplex_ptr dd1,
    magma_tally4DoubleComplex_ptr dd2,
    magma_tally4DoubleComplex_ptr dskp,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zbicgmerge_spmv1( 
    magma_tally4_z_matrix A,
    magma_tally4DoubleComplex_ptr dd1,
    magma_tally4DoubleComplex_ptr dd2,
    magma_tally4DoubleComplex_ptr dp,
    magma_tally4DoubleComplex_ptr dr,
    magma_tally4DoubleComplex_ptr dv,
    magma_tally4DoubleComplex_ptr dskp,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zbicgmerge_spmv2( 
    magma_tally4_z_matrix A,
    magma_tally4DoubleComplex_ptr dd1,
    magma_tally4DoubleComplex_ptr dd2,
    magma_tally4DoubleComplex_ptr ds,
    magma_tally4DoubleComplex_ptr dt,
    magma_tally4DoubleComplex_ptr dskp,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zbicgmerge_xrbeta( 
    int n,
    magma_tally4DoubleComplex_ptr dd1,
    magma_tally4DoubleComplex_ptr dd2,
    magma_tally4DoubleComplex_ptr drr,
    magma_tally4DoubleComplex_ptr dr,
    magma_tally4DoubleComplex_ptr dp,
    magma_tally4DoubleComplex_ptr ds,
    magma_tally4DoubleComplex_ptr dt,
    magma_tally4DoubleComplex_ptr dx, 
    magma_tally4DoubleComplex_ptr dskp,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zbcsrswp(
    magma_tally4_int_t n,
    magma_tally4_int_t size_b, 
    magma_tally4_int_t *ipiv,
    magma_tally4DoubleComplex_ptr dx,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zbcsrtrsv(
    magma_tally4_uplo_t uplo,
    magma_tally4_int_t r_blocks,
    magma_tally4_int_t c_blocks,
    magma_tally4_int_t size_b, 
    magma_tally4DoubleComplex_ptr dA,
    magma_tally4_index_t *blockinfo, 
    magma_tally4DoubleComplex_ptr dx,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zbcsrvalcpy(
    magma_tally4_int_t size_b, 
    magma_tally4_int_t num_blocks, 
    magma_tally4_int_t num_zero_blocks, 
    magma_tally4DoubleComplex_ptr *dAval, 
    magma_tally4DoubleComplex_ptr *dBval,
    magma_tally4DoubleComplex_ptr *dBval2,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zbcsrluegemm(
    magma_tally4_int_t size_b, 
    magma_tally4_int_t num_block_rows,
    magma_tally4_int_t kblocks,
    magma_tally4DoubleComplex_ptr *dA, 
    magma_tally4DoubleComplex_ptr *dB, 
    magma_tally4DoubleComplex_ptr *dC,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zbcsrlupivloc(
    magma_tally4_int_t size_b, 
    magma_tally4_int_t kblocks,
    magma_tally4DoubleComplex_ptr *dA, 
    magma_tally4_int_t *ipiv,
    magma_tally4_queue_t queue );

magma_tally4_int_t
magma_tally4_zbcsrblockinfo5(
    magma_tally4_int_t lustep,
    magma_tally4_int_t num_blocks, 
    magma_tally4_int_t c_blocks, 
    magma_tally4_int_t size_b,
    magma_tally4_index_t *blockinfo,
    magma_tally4DoubleComplex_ptr dval,
    magma_tally4DoubleComplex_ptr *AII,
    magma_tally4_queue_t queue );

 
#ifdef __cplusplus
}
#endif

#undef PRECISION_z
#endif /* MAGMA_tally4SPARSE_Z_H */
