/*
 -- MAGMA_tally2 (version 1.1) --
 Univ. of Tennessee, Knoxville
 Univ. of California, Berkeley
 Univ. of Colorado, Denver
 @date

 @precisions normal z -> s d c
 @author Hartwig Anzt
*/

#ifndef MAGMA_tally2SPARSE_Z_H
#define MAGMA_tally2SPARSE_Z_H

#include "magma_tally2_types.h"
#include "magma_tally2sparse_types.h"

#define PRECISION_z


#ifdef __cplusplus
extern "C" {
#endif


/* ////////////////////////////////////////////////////////////////////////////
 -- For backwards compatability, map old (1.6.1) to new (1.6.2) function names
*/

#define magma_tally2_z_mtranspose  magma_tally2_zmtranspose
#define magma_tally2_z_mtransfer   magma_tally2_zmtransfer
#define magma_tally2_z_vtransfer   magma_tally2_zmtransfer
#define magma_tally2_z_mconvert    magma_tally2_zmconvert
#define magma_tally2_z_vinit       magma_tally2_zvinit
#define magma_tally2_z_vvisu       magma_tally2_zprint_vector
#define magma_tally2_z_vread       magma_tally2_zvread
#define magma_tally2_z_vspread     magma_tally2_zvspread
#define magma_tally2_z_mvisu       magma_tally2_zprint_matrix
#define magma_tally2_z_mfree       magma_tally2_zmfree
#define magma_tally2_z_vfree       magma_tally2_zmfree
#define write_z_csr_mtx     magma_tally2_zwrite_csr_mtx
#define write_z_csrtomtx    magma_tally2_zwrite_csrtomtx
#define print_z_csr         magma_tally2_zprint_csr_mtx


/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally2_SPARSE Auxiliary functions
*/


magma_tally2_int_t
magma_tally2_zparse_opts( 
    int argc, 
    char** argv, 
    magma_tally2_zopts *opts, 
    int *matrices, 
    magma_tally2_queue_t queue );

magma_tally2_int_t 
read_z_csr_from_binary( 
    magma_tally2_int_t* n_row, 
    magma_tally2_int_t* n_col, 
    magma_tally2_int_t* nnz, 
    magma_tally2DoubleComplex **val, 
    magma_tally2_index_t **row, 
    magma_tally2_index_t **col,
    const char * filename,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
read_z_csr_from_mtx( 
    magma_tally2_storage_t *type, 
    magma_tally2_location_t *location,
    magma_tally2_int_t* n_row, 
    magma_tally2_int_t* n_col, 
    magma_tally2_int_t* nnz, 
    magma_tally2DoubleComplex **val, 
    magma_tally2_index_t **row, 
    magma_tally2_index_t **col, 
    const char *filename,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_z_csr_mtx( 
    magma_tally2_z_matrix *A, 
    const char *filename,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_zcsrset( 
    magma_tally2_int_t m, 
    magma_tally2_int_t n, 
    magma_tally2_index_t *row, 
    magma_tally2_index_t *col, 
    magma_tally2DoubleComplex *val,
    magma_tally2_z_matrix *A,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_zcsrget( 
    magma_tally2_z_matrix A,
    magma_tally2_int_t *m, 
    magma_tally2_int_t *n, 
    magma_tally2_index_t **row, 
    magma_tally2_index_t **col, 
    magma_tally2DoubleComplex **val,
    magma_tally2_queue_t queue );


magma_tally2_int_t 
magma_tally2_zvset( 
    magma_tally2_int_t m, 
    magma_tally2_int_t n, 
    magma_tally2DoubleComplex *val,
    magma_tally2_z_matrix *v,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_zvget( 
    magma_tally2_z_matrix v,
    magma_tally2_int_t *m, 
    magma_tally2_int_t *n, 
    magma_tally2DoubleComplex **val,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_zvset_dev( 
    magma_tally2_int_t m, 
    magma_tally2_int_t n, 
    magma_tally2DoubleComplex_ptr val,
    magma_tally2_z_matrix *v,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_zvget_dev( 
    magma_tally2_z_matrix v,
    magma_tally2_int_t *m, 
    magma_tally2_int_t *n, 
    magma_tally2DoubleComplex_ptr *val,
    magma_tally2_queue_t queue );


magma_tally2_int_t 
magma_tally2_z_csr_mtxsymm( 
    magma_tally2_z_matrix *A, 
    const char *filename,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_z_csr_compressor( 
    magma_tally2DoubleComplex ** val, 
    magma_tally2_index_t ** row, 
    magma_tally2_index_t ** col, 
    magma_tally2DoubleComplex ** valn, 
    magma_tally2_index_t ** rown, 
    magma_tally2_index_t ** coln, 
    magma_tally2_int_t *n,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zmcsrcompressor( 
    magma_tally2_z_matrix *A,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zmcsrcompressor_gpu( 
    magma_tally2_z_matrix *A,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zvtranspose( 
    magma_tally2_z_matrix x,
    magma_tally2_z_matrix *y,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_z_cucsrtranspose( 
    magma_tally2_z_matrix A, 
    magma_tally2_z_matrix *B,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
z_transpose_csr( 
    magma_tally2_int_t n_rows, 
    magma_tally2_int_t n_cols, 
    magma_tally2_int_t nnz,
    magma_tally2DoubleComplex *val, 
    magma_tally2_index_t *row, 
    magma_tally2_index_t *col, 
    magma_tally2_int_t *new_n_rows, 
    magma_tally2_int_t *new_n_cols, 
    magma_tally2_int_t *new_nnz, 
    magma_tally2DoubleComplex **new_val, 
    magma_tally2_index_t **new_row, 
    magma_tally2_index_t **new_col,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zcsrsplit( 
    magma_tally2_int_t bsize,
    magma_tally2_z_matrix A,
    magma_tally2_z_matrix *D,
    magma_tally2_z_matrix *R,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zmscale( 
    magma_tally2_z_matrix *A, 
    magma_tally2_scale_t scaling,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_zmdiff( 
    magma_tally2_z_matrix A, 
    magma_tally2_z_matrix B, 
 real_Double_t *res,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zmdiagadd( 
    magma_tally2_z_matrix *A, 
    magma_tally2DoubleComplex add,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_zmsort( 
    magma_tally2_z_matrix *A,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zindexsort(
    magma_tally2_index_t *x, 
    magma_tally2_int_t first,
    magma_tally2_int_t last,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zdomainoverlap(
    magma_tally2_index_t num_rows,
    magma_tally2_index_t *num_indices,
    magma_tally2_index_t *rowptr,
    magma_tally2_index_t *colidx,
    magma_tally2_index_t *x,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zsymbilu( 
    magma_tally2_z_matrix *A, 
    magma_tally2_int_t levels,
    magma_tally2_z_matrix *L,
    magma_tally2_z_matrix *U,
    magma_tally2_queue_t queue );


magma_tally2_int_t 
magma_tally2_zwrite_csr_mtx( 
    magma_tally2_z_matrix A,
    magma_tally2_order_t MajorType,
 const char *filename,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_zwrite_csrtomtx( 
    magma_tally2_z_matrix A,
    const char *filename,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_zprint_csr( 
    magma_tally2_int_t n_row, 
    magma_tally2_int_t n_col, 
    magma_tally2_int_t nnz, 
    magma_tally2DoubleComplex **val, 
    magma_tally2_index_t **row, 
    magma_tally2_index_t **col,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_zprint_csr_mtx( 
    magma_tally2_int_t n_row, 
    magma_tally2_int_t n_col, 
    magma_tally2_int_t nnz, 
    magma_tally2DoubleComplex **val, 
    magma_tally2_index_t **row, 
    magma_tally2_index_t **col, 
    magma_tally2_order_t MajorType,
    magma_tally2_queue_t queue );


magma_tally2_int_t 
magma_tally2_zmtranspose(
    magma_tally2_z_matrix A, 
    magma_tally2_z_matrix *B,
    magma_tally2_queue_t queue );


magma_tally2_int_t 
magma_tally2_zmtransfer(
    magma_tally2_z_matrix A, 
    magma_tally2_z_matrix *B, 
    magma_tally2_location_t src, 
    magma_tally2_location_t dst,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_zmconvert(
    magma_tally2_z_matrix A, 
    magma_tally2_z_matrix *B, 
    magma_tally2_storage_t old_format, 
    magma_tally2_storage_t new_format,
    magma_tally2_queue_t queue );


magma_tally2_int_t
magma_tally2_zvinit(
    magma_tally2_z_matrix *x, 
    magma_tally2_location_t memory_location,
    magma_tally2_int_t num_rows, 
    magma_tally2_int_t num_cols,
    magma_tally2DoubleComplex values,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zprint_vector(
    magma_tally2_z_matrix x, 
    magma_tally2_int_t offset, 
    magma_tally2_int_t displaylength,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zvread(
    magma_tally2_z_matrix *x, 
    magma_tally2_int_t length,
    char * filename,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zvspread(
    magma_tally2_z_matrix *x, 
    const char * filename,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zprint_matrix(
    magma_tally2_z_matrix A,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_zdiameter(
    magma_tally2_z_matrix *A,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_zrowentries(
    magma_tally2_z_matrix *A,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zmfree(
    magma_tally2_z_matrix *A,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zresidual(
    magma_tally2_z_matrix A, 
    magma_tally2_z_matrix b, 
    magma_tally2_z_matrix x, 
    double *res,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zresidualvec(
    magma_tally2_z_matrix A,
    magma_tally2_z_matrix b,
    magma_tally2_z_matrix x,
    magma_tally2_z_matrix *r,
    double *res,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zmgenerator(
    magma_tally2_int_t n,
    magma_tally2_int_t offdiags,
    magma_tally2_index_t *diag_offset,
    magma_tally2DoubleComplex *diag_vals,
    magma_tally2_z_matrix *A,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zm_27stencil(
    magma_tally2_int_t n,
    magma_tally2_z_matrix *A,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zm_5stencil(
    magma_tally2_int_t n,
    magma_tally2_z_matrix *A,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zsolverinfo(
    magma_tally2_z_solver_par *solver_par, 
    magma_tally2_z_preconditioner *precond_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zsolverinfo_init(
    magma_tally2_z_solver_par *solver_par, 
    magma_tally2_z_preconditioner *precond,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zeigensolverinfo_init(
    magma_tally2_z_solver_par *solver_par,
    magma_tally2_queue_t queue );


magma_tally2_int_t
magma_tally2_zsolverinfo_free(
    magma_tally2_z_solver_par *solver_par, 
    magma_tally2_z_preconditioner *precond,
    magma_tally2_queue_t queue );



/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally2_SPARSE iterative incomplete factorizations
*/


magma_tally2_int_t
magma_tally2_ziterilusetup( 
    magma_tally2_z_matrix A, 
    magma_tally2_z_matrix b,                                 
    magma_tally2_z_preconditioner *precond,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zitericsetup( 
    magma_tally2_z_matrix A, 
    magma_tally2_z_matrix b, 
    magma_tally2_z_preconditioner *precond,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zitericupdate( 
    magma_tally2_z_matrix A, 
    magma_tally2_z_preconditioner *precond, 
    magma_tally2_int_t updates,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zapplyiteric_l( 
    magma_tally2_z_matrix b, 
    magma_tally2_z_matrix *x, 
    magma_tally2_z_preconditioner *precond,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zapplyiteric_r( 
    magma_tally2_z_matrix b, 
    magma_tally2_z_matrix *x, 
    magma_tally2_z_preconditioner *precond,
    magma_tally2_queue_t queue );


magma_tally2_int_t
magma_tally2_ziterilu_csr( 
    magma_tally2_z_matrix A,
    magma_tally2_z_matrix L,
    magma_tally2_z_matrix U,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_ziteric_csr( 
    magma_tally2_z_matrix A,
    magma_tally2_z_matrix A_CSR,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_zfrobenius( 
    magma_tally2_z_matrix A, 
    magma_tally2_z_matrix B, 
    real_Double_t *res,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_znonlinres(   
    magma_tally2_z_matrix A, 
    magma_tally2_z_matrix L,
    magma_tally2_z_matrix U, 
    magma_tally2_z_matrix *LU, 
    real_Double_t *res,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_zilures(   
    magma_tally2_z_matrix A, 
    magma_tally2_z_matrix L,
    magma_tally2_z_matrix U, 
    magma_tally2_z_matrix *LU, 
    real_Double_t *res,
    real_Double_t *nonlinres,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_zicres(       
    magma_tally2_z_matrix A, 
    magma_tally2_z_matrix C,
    magma_tally2_z_matrix CT, 
    magma_tally2_z_matrix *LU, 
    real_Double_t *res,
    real_Double_t *nonlinres,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_zinitguess( 
    magma_tally2_z_matrix A, 
    magma_tally2_z_matrix *L, 
    magma_tally2_z_matrix *U,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_zinitrecursiveLU( 
    magma_tally2_z_matrix A, 
    magma_tally2_z_matrix *B,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_zmLdiagadd( 
    magma_tally2_z_matrix *L,
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
magma_tally2_zcg(
    magma_tally2_z_matrix A, magma_tally2_z_matrix b, 
    magma_tally2_z_matrix *x, magma_tally2_z_solver_par *solver_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_zcg_res(
    magma_tally2_z_matrix A, magma_tally2_z_matrix b, 
    magma_tally2_z_matrix *x, magma_tally2_z_solver_par *solver_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_zcg_merge(
    magma_tally2_z_matrix A, magma_tally2_z_matrix b, 
    magma_tally2_z_matrix *x, magma_tally2_z_solver_par *solver_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_zgmres(
    magma_tally2_z_matrix A, magma_tally2_z_matrix b, 
    magma_tally2_z_matrix *x, magma_tally2_z_solver_par *solver_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zbicgstab(
    magma_tally2_z_matrix A, magma_tally2_z_matrix b, magma_tally2_z_matrix *x, 
    magma_tally2_z_solver_par *solver_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zbicgstab_merge(
    magma_tally2_z_matrix A, magma_tally2_z_matrix b, 
    magma_tally2_z_matrix *x, magma_tally2_z_solver_par *solver_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zbicgstab_merge2(
    magma_tally2_z_matrix A, magma_tally2_z_matrix b, 
    magma_tally2_z_matrix *x, magma_tally2_z_solver_par *solver_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zpcg(
    magma_tally2_z_matrix A, magma_tally2_z_matrix b, 
    magma_tally2_z_matrix *x, magma_tally2_z_solver_par *solver_par, 
    magma_tally2_z_preconditioner *precond_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zbpcg(
    magma_tally2_z_matrix A, magma_tally2_z_matrix b, 
    magma_tally2_z_matrix *x, magma_tally2_z_solver_par *solver_par, 
    magma_tally2_z_preconditioner *precond_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zpbicgstab(
    magma_tally2_z_matrix A, magma_tally2_z_matrix b, 
    magma_tally2_z_matrix *x, magma_tally2_z_solver_par *solver_par, 
    magma_tally2_z_preconditioner *precond_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zpgmres(
    magma_tally2_z_matrix A, magma_tally2_z_matrix b, 
    magma_tally2_z_matrix *x, magma_tally2_z_solver_par *solver_par, 
    magma_tally2_z_preconditioner *precond_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zfgmres(
    magma_tally2_z_matrix A, magma_tally2_z_matrix b, 
    magma_tally2_z_matrix *x, magma_tally2_z_solver_par *solver_par, 
    magma_tally2_z_preconditioner *precond_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zjacobi(
    magma_tally2_z_matrix A, magma_tally2_z_matrix b, 
    magma_tally2_z_matrix *x, magma_tally2_z_solver_par *solver_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zjacobidomainoverlap(
    magma_tally2_z_matrix A, 
    magma_tally2_z_matrix b, 
    magma_tally2_z_matrix *x,  
    magma_tally2_z_solver_par *solver_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zbaiter(
    magma_tally2_z_matrix A, magma_tally2_z_matrix b, 
    magma_tally2_z_matrix *x, magma_tally2_z_solver_par *solver_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_ziterref(
    magma_tally2_z_matrix A, magma_tally2_z_matrix b, 
    magma_tally2_z_matrix *x, magma_tally2_z_solver_par *solver_par, 
    magma_tally2_z_preconditioner *precond_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zilu(
    magma_tally2_z_matrix A, magma_tally2_z_matrix b, 
    magma_tally2_z_matrix *x, magma_tally2_z_solver_par *solver_par, 
    magma_tally2_int_t *ipiv,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zbcsrlu(
    magma_tally2_z_matrix A, magma_tally2_z_matrix b, 
    magma_tally2_z_matrix *x, magma_tally2_z_solver_par *solver_par,
    magma_tally2_queue_t queue );


magma_tally2_int_t
magma_tally2_zbcsrlutrf(
    magma_tally2_z_matrix A, 
    magma_tally2_z_matrix *M,
    magma_tally2_int_t *ipiv, 
    magma_tally2_int_t version,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zbcsrlusv(
    magma_tally2_z_matrix A, magma_tally2_z_matrix b, 
    magma_tally2_z_matrix *x, magma_tally2_z_solver_par *solver_par, 
    magma_tally2_int_t *ipiv,
    magma_tally2_queue_t queue );



magma_tally2_int_t
magma_tally2_zilucg(
    magma_tally2_z_matrix A, magma_tally2_z_matrix b, 
    magma_tally2_z_matrix *x, magma_tally2_z_solver_par *solver_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zilugmres(
    magma_tally2_z_matrix A, magma_tally2_z_matrix b, 
    magma_tally2_z_matrix *x, magma_tally2_z_solver_par *solver_par,
    magma_tally2_queue_t queue ); 


magma_tally2_int_t
magma_tally2_zlobpcg_shift(
    magma_tally2_int_t num_rows,
    magma_tally2_int_t num_vecs, 
    magma_tally2_int_t shift,
    magma_tally2DoubleComplex_ptr x,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zlobpcg_res(
    magma_tally2_int_t num_rows,
    magma_tally2_int_t num_vecs, 
    double *evalues, 
    magma_tally2DoubleComplex_ptr X,
    magma_tally2DoubleComplex_ptr R, 
    double *res,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zlobpcg_maxpy(
    magma_tally2_int_t num_rows,
    magma_tally2_int_t num_vecs, 
    magma_tally2DoubleComplex_ptr X,
    magma_tally2DoubleComplex_ptr Y,
    magma_tally2_queue_t queue );


/*/////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally2_SPARSE eigensolvers (Data on GPU)
*/
magma_tally2_int_t
magma_tally2_zlobpcg(
    magma_tally2_z_matrix A, 
    magma_tally2_z_solver_par *solver_par,
    magma_tally2_z_preconditioner *precond_par, 
    magma_tally2_queue_t queue );




/*/////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally2_SPARSE preconditioners (Data on GPU)
*/

magma_tally2_int_t
magma_tally2_zjacobisetup(
    magma_tally2_z_matrix A, 
    magma_tally2_z_matrix b, 
    magma_tally2_z_matrix *M, 
    magma_tally2_z_matrix *c,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zjacobisetup_matrix(
    magma_tally2_z_matrix A, 
    magma_tally2_z_matrix *M, 
    magma_tally2_z_matrix *d,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zjacobisetup_vector(
    magma_tally2_z_matrix b, 
    magma_tally2_z_matrix d, 
    magma_tally2_z_matrix *c,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zjacobiiter(
    magma_tally2_z_matrix M, 
    magma_tally2_z_matrix c, 
    magma_tally2_z_matrix *x, 
    magma_tally2_z_solver_par *solver_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zjacobiiter_precond( 
    magma_tally2_z_matrix M, 
    magma_tally2_z_matrix *x, 
    magma_tally2_z_solver_par *solver_par, 
    magma_tally2_z_preconditioner *precond,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zjacobiiter_sys(
    magma_tally2_z_matrix A, 
    magma_tally2_z_matrix b, 
    magma_tally2_z_matrix d, 
    magma_tally2_z_matrix t, 
    magma_tally2_z_matrix *x,  
    magma_tally2_z_solver_par *solver_par,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zpastixsetup(
    magma_tally2_z_matrix A, magma_tally2_z_matrix b,
    magma_tally2_z_preconditioner *precond,
    magma_tally2_queue_t queue );


magma_tally2_int_t
magma_tally2_zapplypastix(
    magma_tally2_z_matrix b, magma_tally2_z_matrix *x, 
    magma_tally2_z_preconditioner *precond,
    magma_tally2_queue_t queue );


// custom preconditioner
magma_tally2_int_t
magma_tally2_zapplycustomprecond_l(
    magma_tally2_z_matrix b, 
    magma_tally2_z_matrix *x, 
    magma_tally2_z_preconditioner *precond,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zapplycustomprecond_r(
    magma_tally2_z_matrix b, 
    magma_tally2_z_matrix *x, 
    magma_tally2_z_preconditioner *precond,
    magma_tally2_queue_t queue );


// CUSPARSE preconditioner

magma_tally2_int_t
magma_tally2_zcuilusetup(
    magma_tally2_z_matrix A, magma_tally2_z_preconditioner *precond,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zapplycuilu_l(
    magma_tally2_z_matrix b, magma_tally2_z_matrix *x, 
    magma_tally2_z_preconditioner *precond,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zapplycuilu_r(
    magma_tally2_z_matrix b, magma_tally2_z_matrix *x, 
    magma_tally2_z_preconditioner *precond,
    magma_tally2_queue_t queue );


magma_tally2_int_t
magma_tally2_zcuiccsetup(
    magma_tally2_z_matrix A, magma_tally2_z_preconditioner *precond,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zapplycuicc_l(
    magma_tally2_z_matrix b, magma_tally2_z_matrix *x, 
    magma_tally2_z_preconditioner *precond,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zapplycuicc_r(
    magma_tally2_z_matrix b, magma_tally2_z_matrix *x, 
    magma_tally2_z_preconditioner *precond,
    magma_tally2_queue_t queue );


magma_tally2_int_t
magma_tally2_zcumilusetup(
    magma_tally2_z_matrix A, 
    magma_tally2_z_preconditioner *precond,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zcumilugeneratesolverinfo(
    magma_tally2_z_preconditioner *precond,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zapplycumilu_l(
    magma_tally2_z_matrix b, 
    magma_tally2_z_matrix *x, 
    magma_tally2_z_preconditioner *precond,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zapplycumilu_r(
    magma_tally2_z_matrix b, 
    magma_tally2_z_matrix *x, 
    magma_tally2_z_preconditioner *precond,
    magma_tally2_queue_t queue );


magma_tally2_int_t
magma_tally2_zcumiccsetup(
    magma_tally2_z_matrix A, 
    magma_tally2_z_preconditioner *precond,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zcumicgeneratesolverinfo(
    magma_tally2_z_preconditioner *precond,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zapplycumicc_l(
    magma_tally2_z_matrix b, 
    magma_tally2_z_matrix *x, 
    magma_tally2_z_preconditioner *precond,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zapplycumicc_r(
    magma_tally2_z_matrix b, 
    magma_tally2_z_matrix *x, 
    magma_tally2_z_preconditioner *precond,
    magma_tally2_queue_t queue );


// block-asynchronous iteration

magma_tally2_int_t
magma_tally2_zbajac_csr(
    magma_tally2_int_t localiters,
    magma_tally2_z_matrix D,
    magma_tally2_z_matrix R,
    magma_tally2_z_matrix b,
    magma_tally2_z_matrix *x,
    magma_tally2_queue_t queue );

/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally2_SPARSE utility function definitions
*/

magma_tally2_int_t
magma_tally2_z_spmv(
    magma_tally2DoubleComplex alpha, 
    magma_tally2_z_matrix A, 
    magma_tally2_z_matrix x, 
    magma_tally2DoubleComplex beta, 
    magma_tally2_z_matrix y,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zcustomspmv(
    magma_tally2DoubleComplex alpha, 
    magma_tally2_z_matrix x, 
    magma_tally2DoubleComplex beta, 
    magma_tally2_z_matrix y,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_z_spmv_shift(
    magma_tally2DoubleComplex alpha, 
    magma_tally2_z_matrix A, 
    magma_tally2DoubleComplex lambda,
    magma_tally2_z_matrix x, 
    magma_tally2DoubleComplex beta, 
    magma_tally2_int_t offset, 
    magma_tally2_int_t blocksize,
    magma_tally2Index_ptr dadd_vecs, 
    magma_tally2_z_matrix y,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zcuspmm(
    magma_tally2_z_matrix A, 
    magma_tally2_z_matrix B, 
    magma_tally2_z_matrix *AB,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_z_spmm(
    magma_tally2DoubleComplex alpha, 
    magma_tally2_z_matrix A,
    magma_tally2_z_matrix B,
    magma_tally2_z_matrix *C,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zsymbilu( 
    magma_tally2_z_matrix *A, 
    magma_tally2_int_t levels, 
    magma_tally2_z_matrix *L, 
    magma_tally2_z_matrix *U,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zcuspaxpy(
    magma_tally2DoubleComplex_ptr alpha, magma_tally2_z_matrix A, 
    magma_tally2DoubleComplex_ptr beta, magma_tally2_z_matrix B, 
    magma_tally2_z_matrix *AB,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_z_precond(
    magma_tally2_z_matrix A, 
    magma_tally2_z_matrix b, magma_tally2_z_matrix *x,
    magma_tally2_z_preconditioner *precond,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_z_solver(
    magma_tally2_z_matrix A, magma_tally2_z_matrix b, 
    magma_tally2_z_matrix *x, magma_tally2_zopts *zopts,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_z_precondsetup(
    magma_tally2_z_matrix A, magma_tally2_z_matrix b, 
    magma_tally2_z_preconditioner *precond,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_z_applyprecond(
    magma_tally2_z_matrix A, magma_tally2_z_matrix b, 
    magma_tally2_z_matrix *x, magma_tally2_z_preconditioner *precond,
    magma_tally2_queue_t queue );


magma_tally2_int_t
magma_tally2_z_applyprecond_left(
    magma_tally2_z_matrix A, magma_tally2_z_matrix b, 
    magma_tally2_z_matrix *x, magma_tally2_z_preconditioner *precond,
    magma_tally2_queue_t queue );


magma_tally2_int_t
magma_tally2_z_applyprecond_right(
    magma_tally2_z_matrix A, magma_tally2_z_matrix b, 
    magma_tally2_z_matrix *x, magma_tally2_z_preconditioner *precond,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_z_initP2P(
    magma_tally2_int_t *bandwidth_benchmark,
    magma_tally2_int_t *num_gpus,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zcompact(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda,
    double *dnorms, double tol, 
    magma_tally2_int_t *activeMask, magma_tally2_int_t *cBlockSize,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zcompactActive(
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex_ptr dA, magma_tally2_int_t ldda, 
    magma_tally2_int_t *active,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zmlumerge(    
    magma_tally2_z_matrix L, 
    magma_tally2_z_matrix U,
    magma_tally2_z_matrix *A, 
    magma_tally2_queue_t queue );


/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA_tally2_SPARSE BLAS function definitions
*/
magma_tally2_int_t 
magma_tally2_zgecsrmv(
    magma_tally2_trans_t transA,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex alpha,
    magma_tally2DoubleComplex_ptr dval,
    magma_tally2Index_ptr drowptr,
    magma_tally2Index_ptr dcolind,
    magma_tally2DoubleComplex_ptr dx,
    magma_tally2DoubleComplex beta,
    magma_tally2DoubleComplex_ptr dy,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_zgecsrmv_shift(
    magma_tally2_trans_t transA,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2DoubleComplex alpha,
    magma_tally2DoubleComplex lambda,
    magma_tally2DoubleComplex_ptr dval,
    magma_tally2Index_ptr drowptr,
    magma_tally2Index_ptr dcolind,
    magma_tally2DoubleComplex_ptr dx,
    magma_tally2DoubleComplex beta,
    int offset,
    int blocksize,
    magma_tally2Index_ptr dadd_rows,
    magma_tally2DoubleComplex_ptr dy,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_zmgecsrmv(
    magma_tally2_trans_t transA,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2_int_t num_vecs,
    magma_tally2DoubleComplex alpha,
    magma_tally2DoubleComplex_ptr dval,
    magma_tally2Index_ptr drowptr,
    magma_tally2Index_ptr dcolind,
    magma_tally2DoubleComplex_ptr dx,
    magma_tally2DoubleComplex beta,
    magma_tally2DoubleComplex_ptr dy,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_zgeellmv(
    magma_tally2_trans_t transA,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2_int_t nnz_per_row,
    magma_tally2DoubleComplex alpha,
    magma_tally2DoubleComplex_ptr dval,
    magma_tally2Index_ptr dcolind,
    magma_tally2DoubleComplex_ptr dx,
    magma_tally2DoubleComplex beta,
    magma_tally2DoubleComplex_ptr dy,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_zgeellmv_shift(
    magma_tally2_trans_t transA,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2_int_t nnz_per_row,
    magma_tally2DoubleComplex alpha,
    magma_tally2DoubleComplex lambda,
    magma_tally2DoubleComplex_ptr dval,
    magma_tally2Index_ptr dcolind,
    magma_tally2DoubleComplex_ptr dx,
    magma_tally2DoubleComplex beta,
    int offset,
    int blocksize,
    magma_tally2Index_ptr dadd_rows,
    magma_tally2DoubleComplex_ptr dy,
    magma_tally2_queue_t queue );


magma_tally2_int_t 
magma_tally2_zmgeellmv(
    magma_tally2_trans_t transA,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2_int_t num_vecs,
    magma_tally2_int_t nnz_per_row,
    magma_tally2DoubleComplex alpha,
    magma_tally2DoubleComplex_ptr dval,
    magma_tally2Index_ptr dcolind,
    magma_tally2DoubleComplex_ptr dx,
    magma_tally2DoubleComplex beta,
    magma_tally2DoubleComplex_ptr dy,
    magma_tally2_queue_t queue );


magma_tally2_int_t 
magma_tally2_zgeelltmv(
    magma_tally2_trans_t transA,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2_int_t nnz_per_row,
    magma_tally2DoubleComplex alpha,
    magma_tally2DoubleComplex_ptr dval,
    magma_tally2Index_ptr dcolind,
    magma_tally2DoubleComplex_ptr dx,
    magma_tally2DoubleComplex beta,
    magma_tally2DoubleComplex_ptr dy,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_zgeelltmv_shift(
    magma_tally2_trans_t transA,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2_int_t nnz_per_row,
    magma_tally2DoubleComplex alpha,
    magma_tally2DoubleComplex lambda,
    magma_tally2DoubleComplex_ptr dval,
    magma_tally2Index_ptr dcolind,
    magma_tally2DoubleComplex_ptr dx,
    magma_tally2DoubleComplex beta,
    int offset,
    int blocksize,
    magma_tally2Index_ptr dadd_rows,
    magma_tally2DoubleComplex_ptr dy,
    magma_tally2_queue_t queue );


magma_tally2_int_t 
magma_tally2_zmgeelltmv(
    magma_tally2_trans_t transA,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2_int_t num_vecs,
    magma_tally2_int_t nnz_per_row,
    magma_tally2DoubleComplex alpha,
    magma_tally2DoubleComplex_ptr dval,
    magma_tally2Index_ptr dcolind,
    magma_tally2DoubleComplex_ptr dx,
    magma_tally2DoubleComplex beta,
    magma_tally2DoubleComplex_ptr dy,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_zgeellrtmv(
    magma_tally2_trans_t transA,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2_int_t nnz_per_row,
    magma_tally2DoubleComplex alpha,
    magma_tally2DoubleComplex_ptr dval,
    magma_tally2Index_ptr dcolind,
    magma_tally2Index_ptr drowlength,
    magma_tally2DoubleComplex_ptr dx,
    magma_tally2DoubleComplex beta,
    magma_tally2DoubleComplex_ptr dy,
    magma_tally2_int_t num_threads,
    magma_tally2_int_t threads_per_row,
    magma_tally2_queue_t queue );

magma_tally2_int_t 
magma_tally2_zgesellcmv(
    magma_tally2_trans_t transA,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2_int_t blocksize,
    magma_tally2_int_t slices,
    magma_tally2_int_t alignment,
    magma_tally2DoubleComplex alpha,
    magma_tally2DoubleComplex_ptr dval,
    magma_tally2Index_ptr dcolind,
    magma_tally2Index_ptr drowptr,
    magma_tally2DoubleComplex_ptr dx,
    magma_tally2DoubleComplex beta,
    magma_tally2DoubleComplex_ptr dy,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zgesellpmv(
    magma_tally2_trans_t transA,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2_int_t blocksize,
    magma_tally2_int_t slices,
    magma_tally2_int_t alignment,
    magma_tally2DoubleComplex alpha,
    magma_tally2DoubleComplex_ptr dval,
    magma_tally2Index_ptr dcolind,
    magma_tally2Index_ptr drowptr,
    magma_tally2DoubleComplex_ptr dx,
    magma_tally2DoubleComplex beta,
    magma_tally2DoubleComplex_ptr dy,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zmgesellpmv(
    magma_tally2_trans_t transA,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2_int_t num_vecs,
    magma_tally2_int_t blocksize,
    magma_tally2_int_t slices,
    magma_tally2_int_t alignment,
    magma_tally2DoubleComplex alpha,
    magma_tally2DoubleComplex_ptr dval,
    magma_tally2Index_ptr dcolind,
    magma_tally2Index_ptr drowptr,
    magma_tally2DoubleComplex_ptr dx,
    magma_tally2DoubleComplex beta,
    magma_tally2DoubleComplex_ptr dy,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zmgesellpmv_blocked(
    magma_tally2_trans_t transA,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2_int_t num_vecs,
    magma_tally2_int_t blocksize,
    magma_tally2_int_t slices,
    magma_tally2_int_t alignment,
    magma_tally2DoubleComplex alpha,
    magma_tally2DoubleComplex_ptr dval,
    magma_tally2Index_ptr dcolind,
    magma_tally2Index_ptr drowptr,
    magma_tally2DoubleComplex_ptr dx,
    magma_tally2DoubleComplex beta,
    magma_tally2DoubleComplex_ptr dy,
    magma_tally2_queue_t queue );


magma_tally2_int_t
magma_tally2_zmergedgs(
    magma_tally2_int_t n, 
    magma_tally2_int_t ldh,
    magma_tally2_int_t k, 
    magma_tally2DoubleComplex_ptr dv, 
    magma_tally2DoubleComplex_ptr dr,
    magma_tally2DoubleComplex_ptr dskp,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zcopyscale(    
    int n, 
    int k,
    magma_tally2DoubleComplex_ptr dr, 
    magma_tally2DoubleComplex_ptr dv,
    magma_tally2DoubleComplex_ptr dskp,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_dznrm2scale(    
    int m, 
    magma_tally2DoubleComplex_ptr dr,    
    int lddr, 
    magma_tally2DoubleComplex *drnorm,
    magma_tally2_queue_t queue );


magma_tally2_int_t
magma_tally2_zjacobisetup_vector_gpu(
    int num_rows, 
    magma_tally2_z_matrix b, 
    magma_tally2_z_matrix d, 
    magma_tally2_z_matrix c,
    magma_tally2_z_matrix *x,
    magma_tally2_queue_t queue );


magma_tally2_int_t
magma_tally2_zjacobi_diagscal(    
    int num_rows, 
    magma_tally2_z_matrix d, 
    magma_tally2_z_matrix b, 
    magma_tally2_z_matrix *c,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zjacobiupdate(
    magma_tally2_z_matrix t, 
    magma_tally2_z_matrix b, 
    magma_tally2_z_matrix d, 
    magma_tally2_z_matrix *x,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zjacobispmvupdate(
    magma_tally2_int_t maxiter,
    magma_tally2_z_matrix A, 
    magma_tally2_z_matrix t, 
    magma_tally2_z_matrix b, 
    magma_tally2_z_matrix d, 
    magma_tally2_z_matrix *x,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zjacobispmvupdate_bw(
    magma_tally2_int_t maxiter,
    magma_tally2_z_matrix A, 
    magma_tally2_z_matrix t, 
    magma_tally2_z_matrix b, 
    magma_tally2_z_matrix d, 
    magma_tally2_z_matrix *x,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zjacobispmvupdateselect(
    magma_tally2_int_t maxiter,
    magma_tally2_int_t num_updates,
    magma_tally2_index_t *indices,
    magma_tally2_z_matrix A,
    magma_tally2_z_matrix t, 
    magma_tally2_z_matrix b, 
    magma_tally2_z_matrix d, 
    magma_tally2_z_matrix tmp, 
    magma_tally2_z_matrix *x,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zjacobisetup_diagscal(
    magma_tally2_z_matrix A, magma_tally2_z_matrix *d,
    magma_tally2_queue_t queue );


magma_tally2_int_t
magma_tally2_zbicgmerge1(    
    int n, 
    magma_tally2DoubleComplex_ptr dskp,
    magma_tally2DoubleComplex_ptr dv, 
    magma_tally2DoubleComplex_ptr dr, 
    magma_tally2DoubleComplex_ptr dp,
    magma_tally2_queue_t queue );


magma_tally2_int_t
magma_tally2_zbicgmerge2(
    int n, 
    magma_tally2DoubleComplex_ptr dskp, 
    magma_tally2DoubleComplex_ptr dr,
    magma_tally2DoubleComplex_ptr dv, 
    magma_tally2DoubleComplex_ptr ds,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zbicgmerge3(
    int n, 
    magma_tally2DoubleComplex_ptr dskp, 
    magma_tally2DoubleComplex_ptr dp,
    magma_tally2DoubleComplex_ptr ds,
    magma_tally2DoubleComplex_ptr dt,
    magma_tally2DoubleComplex_ptr dx, 
    magma_tally2DoubleComplex_ptr dr,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zbicgmerge4(
    int type, 
    magma_tally2DoubleComplex_ptr dskp,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zcgmerge_spmv1( 
    magma_tally2_z_matrix A,
    magma_tally2DoubleComplex_ptr d1,
    magma_tally2DoubleComplex_ptr d2,
    magma_tally2DoubleComplex_ptr dd,
    magma_tally2DoubleComplex_ptr dz,
    magma_tally2DoubleComplex_ptr dskp,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zcgmerge_xrbeta( 
    int n,
    magma_tally2DoubleComplex_ptr d1,
    magma_tally2DoubleComplex_ptr d2,
    magma_tally2DoubleComplex_ptr dx,
    magma_tally2DoubleComplex_ptr dr,
    magma_tally2DoubleComplex_ptr dd,
    magma_tally2DoubleComplex_ptr dz, 
    magma_tally2DoubleComplex_ptr dskp,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zmdotc(
    magma_tally2_int_t n, 
    magma_tally2_int_t k, 
    magma_tally2DoubleComplex_ptr dv, 
    magma_tally2DoubleComplex_ptr dr,
    magma_tally2DoubleComplex_ptr dd1,
    magma_tally2DoubleComplex_ptr dd2,
    magma_tally2DoubleComplex_ptr dskp,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zgemvmdot(
    int n, 
    int k, 
    magma_tally2DoubleComplex_ptr dv, 
    magma_tally2DoubleComplex_ptr dr,
    magma_tally2DoubleComplex_ptr dd1,
    magma_tally2DoubleComplex_ptr dd2,
    magma_tally2DoubleComplex_ptr dskp,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zbicgmerge_spmv1( 
    magma_tally2_z_matrix A,
    magma_tally2DoubleComplex_ptr dd1,
    magma_tally2DoubleComplex_ptr dd2,
    magma_tally2DoubleComplex_ptr dp,
    magma_tally2DoubleComplex_ptr dr,
    magma_tally2DoubleComplex_ptr dv,
    magma_tally2DoubleComplex_ptr dskp,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zbicgmerge_spmv2( 
    magma_tally2_z_matrix A,
    magma_tally2DoubleComplex_ptr dd1,
    magma_tally2DoubleComplex_ptr dd2,
    magma_tally2DoubleComplex_ptr ds,
    magma_tally2DoubleComplex_ptr dt,
    magma_tally2DoubleComplex_ptr dskp,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zbicgmerge_xrbeta( 
    int n,
    magma_tally2DoubleComplex_ptr dd1,
    magma_tally2DoubleComplex_ptr dd2,
    magma_tally2DoubleComplex_ptr drr,
    magma_tally2DoubleComplex_ptr dr,
    magma_tally2DoubleComplex_ptr dp,
    magma_tally2DoubleComplex_ptr ds,
    magma_tally2DoubleComplex_ptr dt,
    magma_tally2DoubleComplex_ptr dx, 
    magma_tally2DoubleComplex_ptr dskp,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zbcsrswp(
    magma_tally2_int_t n,
    magma_tally2_int_t size_b, 
    magma_tally2_int_t *ipiv,
    magma_tally2DoubleComplex_ptr dx,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zbcsrtrsv(
    magma_tally2_uplo_t uplo,
    magma_tally2_int_t r_blocks,
    magma_tally2_int_t c_blocks,
    magma_tally2_int_t size_b, 
    magma_tally2DoubleComplex_ptr dA,
    magma_tally2_index_t *blockinfo, 
    magma_tally2DoubleComplex_ptr dx,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zbcsrvalcpy(
    magma_tally2_int_t size_b, 
    magma_tally2_int_t num_blocks, 
    magma_tally2_int_t num_zero_blocks, 
    magma_tally2DoubleComplex_ptr *dAval, 
    magma_tally2DoubleComplex_ptr *dBval,
    magma_tally2DoubleComplex_ptr *dBval2,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zbcsrluegemm(
    magma_tally2_int_t size_b, 
    magma_tally2_int_t num_block_rows,
    magma_tally2_int_t kblocks,
    magma_tally2DoubleComplex_ptr *dA, 
    magma_tally2DoubleComplex_ptr *dB, 
    magma_tally2DoubleComplex_ptr *dC,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zbcsrlupivloc(
    magma_tally2_int_t size_b, 
    magma_tally2_int_t kblocks,
    magma_tally2DoubleComplex_ptr *dA, 
    magma_tally2_int_t *ipiv,
    magma_tally2_queue_t queue );

magma_tally2_int_t
magma_tally2_zbcsrblockinfo5(
    magma_tally2_int_t lustep,
    magma_tally2_int_t num_blocks, 
    magma_tally2_int_t c_blocks, 
    magma_tally2_int_t size_b,
    magma_tally2_index_t *blockinfo,
    magma_tally2DoubleComplex_ptr dval,
    magma_tally2DoubleComplex_ptr *AII,
    magma_tally2_queue_t queue );

 
#ifdef __cplusplus
}
#endif

#undef PRECISION_z
#endif /* MAGMA_tally2SPARSE_Z_H */
