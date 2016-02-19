/*
    -- MAGMA_tally3 (version 1.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date
*/

#ifndef MAGMA_tally3SPARSE_TYPES_H
#define MAGMA_tally3SPARSE_TYPES_H


#if defined(HAVE_PASTIX)
//PaStiX include
#include <stdint.h>
/* to access functions from the libpastix, respect this order */
#include <pastix.h>
#include <read_matrix.h>
#include <get_options.h>
#include <assert.h>
#endif

// includes CUDA
#include <cusparse_v2.h>



#ifdef __cplusplus
extern "C" {
#endif




typedef struct magma_tally3_z_matrix{

    magma_tally3_storage_t    storage_type;            // matrix format - CSR, ELL, SELL-P
    magma_tally3_location_t   memory_location;         // CPU or DEV
    magma_tally3_symmetry_t   sym;                     // opt: indicate symmetry
    magma_tally3_diagorder_t  diagorder_type;          // opt: only needed for factorization matrices
    magma_tally3_fillmode_t   fill_mode;               // fill mode full/lower/upper
    magma_tally3_int_t        num_rows;                // number of rows
    magma_tally3_int_t        num_cols;                // number of columns
    magma_tally3_int_t        nnz;                     // opt: number of nonzeros
    magma_tally3_int_t        max_nnz_row;             // opt: max number of nonzeros in one row
    magma_tally3_int_t        diameter;                // opt: max distance of entry from main diagonal
    union {
        magma_tally3DoubleComplex      *val;           // array containing values in CPU case
        magma_tally3DoubleComplex_ptr  dval;           // array containing values in DEV case
    };
    union {
        magma_tally3DoubleComplex      *diag;          // opt: diagonal entries in CPU case
        magma_tally3DoubleComplex_ptr  ddiag;          // opt: diagonal entries in DEV case
    };
    union {
        magma_tally3_index_t           *row;           // row pointer CPU case
        magma_tally3Index_ptr          drow;           // row pointer DEV case
    };
    union {
        magma_tally3_index_t           *rowidx;        // opt: array containing row indices CPU case
        magma_tally3Index_ptr          drowidx;        // opt: array containing row indices DEV case
    };
    union {
        magma_tally3_index_t           *col;           // array containing col indices CPU case
        magma_tally3Index_ptr          dcol;           // array containing col indices DEV case
    };
    magma_tally3_index_t      *blockinfo;              // opt: for BCSR format CPU case
    magma_tally3_int_t        blocksize;               // opt: info for SELL-P/BCSR
    magma_tally3_int_t        numblocks;               // opt: info for SELL-P/BCSR
    magma_tally3_int_t        alignment;               // opt: info for SELL-P/BCSR
    magma_tally3_order_t      major;                   // opt: row/col major for dense matrices
    magma_tally3_int_t        ld;                      // opt: leading dimension for dense

}magma_tally3_z_matrix;

typedef struct magma_tally3_c_matrix{

    magma_tally3_storage_t    storage_type;            // matrix format - CSR, ELL, SELL-P
    magma_tally3_location_t   memory_location;         // CPU or DEV
    magma_tally3_symmetry_t   sym;                     // opt: indicate symmetry
    magma_tally3_diagorder_t  diagorder_type;          // opt: only needed for factorization matrices
    magma_tally3_fillmode_t   fill_mode;               // fill mode full/lower/upper
    magma_tally3_int_t        num_rows;                // number of rows
    magma_tally3_int_t        num_cols;                // number of columns
    magma_tally3_int_t        nnz;                     // opt: number of nonzeros
    magma_tally3_int_t        max_nnz_row;             // opt: max number of nonzeros in one row
    magma_tally3_int_t        diameter;                // opt: max distance of entry from main diagonal
    union {
        magma_tally3FloatComplex       *val;           // array containing values in CPU case
        magma_tally3FloatComplex_ptr   dval;           // array containing values in DEV case
    };
    union {
        magma_tally3FloatComplex       *diag;          // opt: diagonal entries in CPU case
        magma_tally3FloatComplex_ptr   ddiag;          // opt: diagonal entries in DEV case
    };
    union {
        magma_tally3_index_t           *row;           // row pointer CPU case
        magma_tally3Index_ptr          drow;           // row pointer DEV case
    };
    union {
        magma_tally3_index_t           *rowidx;        // opt: array containing row indices CPU case
        magma_tally3Index_ptr          drowidx;        // opt: array containing row indices DEV case
    };
    union {
        magma_tally3_index_t           *col;           // array containing col indices CPU case
        magma_tally3Index_ptr          dcol;           // array containing col indices DEV case
    };
    magma_tally3_index_t      *blockinfo;              // opt: for BCSR format CPU case
    magma_tally3_int_t        blocksize;               // opt: info for SELL-P/BCSR
    magma_tally3_int_t        numblocks;               // opt: info for SELL-P/BCSR
    magma_tally3_int_t        alignment;               // opt: info for SELL-P/BCSR
    magma_tally3_order_t      major;                   // opt: row/col major for dense matrices
    magma_tally3_int_t        ld;                      // opt: leading dimension for dense

}magma_tally3_c_matrix;


typedef struct magma_tally3_d_matrix{

    magma_tally3_storage_t    storage_type;            // matrix format - CSR, ELL, SELL-P
    magma_tally3_location_t   memory_location;         // CPU or DEV
    magma_tally3_symmetry_t   sym;                     // opt: indicate symmetry
    magma_tally3_diagorder_t  diagorder_type;          // opt: only needed for factorization matrices
    magma_tally3_fillmode_t   fill_mode;               // fill mode full/lower/upper
    magma_tally3_int_t        num_rows;                // number of rows
    magma_tally3_int_t        num_cols;                // number of columns
    magma_tally3_int_t        nnz;                     // opt: number of nonzeros
    magma_tally3_int_t        max_nnz_row;             // opt: max number of nonzeros in one row
    magma_tally3_int_t        diameter;                // opt: max distance of entry from main diagonal
    union {
        double                  *val;           // array containing values in CPU case
        magma_tally3Double_ptr         dval;           // array containing values in DEV case
    };
    union {
        double                  *diag;          // opt: diagonal entries in CPU case
        magma_tally3Double_ptr         ddiag;          // opt: diagonal entries in DEV case
    };
    union {
        magma_tally3_index_t           *row;           // row pointer CPU case
        magma_tally3Index_ptr          drow;           // row pointer DEV case
    };
    union {
        magma_tally3_index_t           *rowidx;        // opt: array containing row indices CPU case
        magma_tally3Index_ptr          drowidx;        // opt: array containing row indices DEV case
    };
    union {
        magma_tally3_index_t           *col;           // array containing col indices CPU case
        magma_tally3Index_ptr          dcol;           // array containing col indices DEV case
    };
    magma_tally3_index_t      *blockinfo;              // opt: for BCSR format CPU case
    magma_tally3_int_t        blocksize;               // opt: info for SELL-P/BCSR
    magma_tally3_int_t        numblocks;               // opt: info for SELL-P/BCSR
    magma_tally3_int_t        alignment;               // opt: info for SELL-P/BCSR
    magma_tally3_order_t      major;                   // opt: row/col major for dense matrices
    magma_tally3_int_t        ld;                      // opt: leading dimension for dense

}magma_tally3_d_matrix;


typedef struct magma_tally3_s_matrix{

    magma_tally3_storage_t    storage_type;            // matrix format - CSR, ELL, SELL-P
    magma_tally3_location_t   memory_location;         // CPU or DEV
    magma_tally3_symmetry_t   sym;                     // opt: indicate symmetry
    magma_tally3_diagorder_t  diagorder_type;          // opt: only needed for factorization matrices
    magma_tally3_fillmode_t   fill_mode;               // fill mode full/lower/upper
    magma_tally3_int_t        num_rows;                // number of rows
    magma_tally3_int_t        num_cols;                // number of columns
    magma_tally3_int_t        nnz;                     // opt: number of nonzeros
    magma_tally3_int_t        max_nnz_row;             // opt: max number of nonzeros in one row
    magma_tally3_int_t        diameter;                // opt: max distance of entry from main diagonal
    union {
        float                   *val;           // array containing values in CPU case
        magma_tally3Float_ptr          dval;           // array containing values in DEV case
    };
    union {
        float                   *diag;          // opt: diagonal entries in CPU case
        magma_tally3Float_ptr          ddiag;          // opt: diagonal entries in DEV case
    };
    union {
        magma_tally3_index_t           *row;           // row pointer CPU case
        magma_tally3Index_ptr          drow;           // row pointer DEV case
    };
    union {
        magma_tally3_index_t           *rowidx;        // opt: array containing row indices CPU case
        magma_tally3Index_ptr          drowidx;        // opt: array containing row indices DEV case
    };
    union {
        magma_tally3_index_t           *col;           // array containing col indices CPU case
        magma_tally3Index_ptr          dcol;           // array containing col indices DEV case
    };
    magma_tally3_index_t      *blockinfo;              // opt: for BCSR format CPU case
    magma_tally3_int_t        blocksize;               // opt: info for SELL-P/BCSR
    magma_tally3_int_t        numblocks;               // opt: info for SELL-P/BCSR
    magma_tally3_int_t        alignment;               // opt: info for SELL-P/BCSR
    magma_tally3_order_t      major;                   // opt: row/col major for dense matrices
    magma_tally3_int_t        ld;                      // opt: leading dimension for dense
    
}magma_tally3_s_matrix;


// for backwards compatability, make these aliases.
typedef magma_tally3_s_matrix magma_tally3_s_sparse_matrix;
typedef magma_tally3_d_matrix magma_tally3_d_sparse_matrix;
typedef magma_tally3_c_matrix magma_tally3_c_sparse_matrix;
typedef magma_tally3_z_matrix magma_tally3_z_sparse_matrix;

typedef magma_tally3_s_matrix magma_tally3_s_vector;
typedef magma_tally3_d_matrix magma_tally3_d_vector;
typedef magma_tally3_c_matrix magma_tally3_c_vector;
typedef magma_tally3_z_matrix magma_tally3_z_vector;

/*
typedef struct magma_tally3_z_vector{

    magma_tally3_location_t   memory_location;         // CPU or DEV
    magma_tally3_int_t        num_rows;                // number of rows
    magma_tally3_int_t        num_cols;                // number of columns (in case of a block of vectors)
    magma_tally3_int_t        nnz;                     // opt: number of nonzeros
    union {
        magma_tally3DoubleComplex      *val;           // array containing values in CPU case
        magma_tally3DoubleComplex_ptr  dval;           // array containing values in DEV case
    };
    magma_tally3_order_t      major;                   // storage type:Row/Column-Major

}magma_tally3_z_vector;

typedef struct magma_tally3_c_vector{

    magma_tally3_location_t   memory_location;         // CPU or DEV
    magma_tally3_int_t        num_rows;                // number of rows
    magma_tally3_int_t        num_cols;                // number of columns (in case of a block of vectors)
    magma_tally3_int_t        nnz;                     // opt: number of nonzeros
    union {
        magma_tally3FloatComplex       *val;           // array containing values in CPU case
        magma_tally3FloatComplex_ptr   dval;           // array containing values in DEV case
    };
    magma_tally3_order_t      major;                   // storage type:Row/Column-Major

}magma_tally3_c_vector;


typedef struct magma_tally3_d_vector{

    magma_tally3_location_t   memory_location;         // CPU or DEV
    magma_tally3_int_t        num_rows;                // number of rows
    magma_tally3_int_t        num_cols;                // number of columns (in case of a block of vectors)
    magma_tally3_int_t        nnz;                     // opt: number of nonzeros
    union {
        double                  *val;           // array containing values in CPU case
        magma_tally3Double_ptr         dval;           // array containing values in DEV case
    };
    magma_tally3_order_t      major;                   // storage type:Row/Column-Major

}magma_tally3_d_vector;


typedef struct magma_tally3_s_vector{

    magma_tally3_location_t   memory_location;         // CPU or DEV
    magma_tally3_int_t        num_rows;                // number of rows
    magma_tally3_int_t        num_cols;                // number of columns (in case of a block of vectors)
    magma_tally3_int_t        nnz;                     // opt: number of nonzeros
    union {
        float                   *val;           // array containing values in CPU case
        magma_tally3Float_ptr          dval;           // array containing values in DEV case
    };
    magma_tally3_order_t      major;                   // storage type:Row/Column-Major

}magma_tally3_s_vector;
*/


//*****************     solver parameters     ********************************//

typedef struct magma_tally3_z_solver_par{

    magma_tally3_solver_type  solver;                  // solver type
    magma_tally3_int_t        version;                 // sometimes there are different versions
    double             epsilon;                 // relative residual stopping criterion
    magma_tally3_int_t        maxiter;                 // upper iteration limit
    magma_tally3_int_t        restart;                 // for GMRES
    magma_tally3_ortho_t      ortho;                   // for GMRES
    magma_tally3_int_t        numiter;                 // feedback: number of needed iterations
    double             init_res;                // feedback: initial residual
    double             final_res;               // feedback: final residual
    double             iter_res;                // feedback: iteratively computed residual
    real_Double_t      runtime;                 // feedback: runtime needed
    real_Double_t      *res_vec;                // feedback: array containing residuals
    real_Double_t      *timing;                 // feedback: detailed timing
    magma_tally3_int_t        verbose;                 // print residual ever 'verbose' iterations
    magma_tally3_int_t        num_eigenvalues;         // number of EV for eigensolvers
    magma_tally3_int_t        ev_length;               // needed for framework
    double             *eigenvalues;            // feedback: array containing eigenvalues
    magma_tally3DoubleComplex_ptr      eigenvectors;   // feedback: array containing eigenvectors on DEV
    magma_tally3_int_t        info;                    // feedback: did the solver converge etc.

//---------------------------------
// the input for verbose is:
// 0 = production mode
// k>0 = convergence and timing is monitored in *res_vec and *timeing every  
// k-th iteration 
//
// the output of info is:
//  0 = convergence (stopping criterion met)
// -1 = no convergence
// -2 = convergence but stopping criterion not met within maxiter
//--------------------------------

}magma_tally3_z_solver_par;



typedef struct magma_tally3_c_solver_par{

    magma_tally3_solver_type  solver;                  // solver type
    magma_tally3_int_t        version;                 // sometimes there are different versions
    float              epsilon;                 // relative residual stopping criterion
    magma_tally3_int_t        maxiter;                 // upper iteration limit
    magma_tally3_int_t        restart;                 // for GMRES
    magma_tally3_ortho_t      ortho;                   // for GMRES
    magma_tally3_int_t        numiter;                 // feedback: number of needed iterations
    float              init_res;                // feedback: initial residual
    float              final_res;               // feedback: final residual
    float              iter_res;                // feedback: iteratively computed residual
    real_Double_t      runtime;                 // feedback: runtime needed
    real_Double_t      *res_vec;                // feedback: array containing residuals
    real_Double_t      *timing;                 // feedback: detailed timing
    magma_tally3_int_t        verbose;                 // print residual ever 'verbose' iterations
    magma_tally3_int_t        num_eigenvalues;         // number of EV for eigensolvers
    magma_tally3_int_t        ev_length;               // needed for framework
    float              *eigenvalues;            // feedback: array containing eigenvalues
    magma_tally3FloatComplex_ptr       eigenvectors;   // feedback: array containing eigenvectors on DEV
    magma_tally3_int_t        info;                    // feedback: did the solver converge etc.

//---------------------------------
// the input for verbose is:
// 0 = production mode
// k>0 = convergence and timing is monitored in *res_vec and *timeing every  
// k-th iteration 
//
// the output of info is:
//  0 = convergence (stopping criterion met)
// -1 = no convergence
// -2 = convergence but stopping criterion not met within maxiter
//--------------------------------

}magma_tally3_c_solver_par;



typedef struct magma_tally3_d_solver_par{

    magma_tally3_solver_type  solver;                  // solver type
    magma_tally3_int_t        version;                 // sometimes there are different versions
    double             epsilon;                 // relative residual stopping criterion
    magma_tally3_int_t        maxiter;                 // upper iteration limit
    magma_tally3_int_t        restart;                 // for GMRES
    magma_tally3_ortho_t      ortho;                   // for GMRES
    magma_tally3_int_t        numiter;                 // feedback: number of needed iterations
    double             init_res;                // feedback: initial residual
    double             final_res;               // feedback: final residual
    double             iter_res;                // feedback: iteratively computed residual
    real_Double_t      runtime;                 // feedback: runtime needed
    real_Double_t      *res_vec;                // feedback: array containing residuals
    real_Double_t      *timing;                 // feedback: detailed timing
    magma_tally3_int_t        verbose;                 // print residual ever 'verbose' iterations
    magma_tally3_int_t        num_eigenvalues;         // number of EV for eigensolvers
    magma_tally3_int_t        ev_length;               // needed for framework
    double             *eigenvalues;            // feedback: array containing eigenvalues
    magma_tally3Double_ptr             eigenvectors;   // feedback: array containing eigenvectors on DEV
    magma_tally3_int_t        info;                    // feedback: did the solver converge etc.

//---------------------------------
// the input for verbose is:
// 0 = production mode
// k>0 = convergence and timing is monitored in *res_vec and *timeing every  
// k-th iteration 
//
// the output of info is:
//  0 = convergence (stopping criterion met)
// -1 = no convergence
// -2 = convergence but stopping criterion not met within maxiter
//--------------------------------

}magma_tally3_d_solver_par;



typedef struct magma_tally3_s_solver_par{

    magma_tally3_solver_type  solver;                  // solver type
    magma_tally3_int_t        version;                 // sometimes there are different versions
    float              epsilon;                 // relative residual stopping criterion
    magma_tally3_int_t        maxiter;                 // upper iteration limit
    magma_tally3_int_t        restart;                 // for GMRES
    magma_tally3_ortho_t      ortho;                   // for GMRES
    magma_tally3_int_t        numiter;                 // feedback: number of needed iterations
    float              init_res;                // feedback: initial residual
    float              final_res;               // feedback: final residual
    float              iter_res;                // feedback: iteratively computed residual
    real_Double_t      runtime;                 // feedback: runtime needed
    real_Double_t      *res_vec;                // feedback: array containing residuals
    real_Double_t      *timing;                 // feedback: detailed timing
    magma_tally3_int_t        verbose;                 // print residual ever 'verbose' iterations
    magma_tally3_int_t        num_eigenvalues;         // number of EV for eigensolvers
    magma_tally3_int_t        ev_length;               // needed for framework
    float              *eigenvalues;            // feedback: array containing eigenvalues
    magma_tally3Float_ptr              eigenvectors;   // feedback: array containing eigenvectors on DEV
    magma_tally3_int_t        info;                    // feedback: did the solver converge etc.

//---------------------------------
// the input for verbose is:
// 0 = production mode
// k>0 = convergence and timing is monitored in *res_vec and *timeing every  
// k-th iteration 
//
// the output of info is:
//       0          Success.
//      -117        Not supported.
//      -201        No convergence within iteration limit. 
//      -202        No convergence.
//      -203        Operator A is not positive definite.
//--------------------------------

}magma_tally3_s_solver_par;



//************            preconditioner parameters       ********************//

typedef struct magma_tally3_z_preconditioner{

    magma_tally3_solver_type       solver;
    magma_tally3_int_t             levels;
    magma_tally3_int_t             sweeps;
    magma_tally3_precision         format;
    double                  epsilon;  
    magma_tally3_int_t             maxiter;
    magma_tally3_int_t             restart; 
    magma_tally3_int_t             numiter;
    double                  init_res;
    double                  final_res;
    magma_tally3_z_matrix   M;
    magma_tally3_z_matrix   L;
    magma_tally3_z_matrix   U;
    magma_tally3_z_matrix   LD;
    magma_tally3_z_matrix   UD;
    magma_tally3_z_matrix          d;
    magma_tally3_z_matrix          d2;
    magma_tally3_z_matrix          work1;
    magma_tally3_z_matrix          work2;
    magma_tally3_int_t*            int_array_1;
    magma_tally3_int_t*            int_array_2;
    cusparseSolveAnalysisInfo_t cuinfo;
    cusparseSolveAnalysisInfo_t cuinfoL;
    cusparseSolveAnalysisInfo_t cuinfoU;
#if defined(HAVE_PASTIX)
    pastix_data_t*          pastix_data;
    magma_tally3_int_t*            iparm;
    double*                 dparm;
#endif

}magma_tally3_z_preconditioner;

typedef struct magma_tally3_c_preconditioner{

    magma_tally3_solver_type       solver;
    magma_tally3_int_t             levels;
    magma_tally3_int_t             sweeps;
    magma_tally3_precision         format;
    float                   epsilon;  
    magma_tally3_int_t             maxiter;
    magma_tally3_int_t             restart; 
    magma_tally3_int_t             numiter;
    float                   init_res;
    float                   final_res;
    magma_tally3_c_matrix   M;
    magma_tally3_c_matrix   L;
    magma_tally3_c_matrix   U;
    magma_tally3_c_matrix   LD;
    magma_tally3_c_matrix   UD;
    magma_tally3_c_matrix          d;
    magma_tally3_c_matrix          d2;
    magma_tally3_c_matrix          work1;
    magma_tally3_c_matrix          work2;
    magma_tally3_int_t*            int_array_1;
    magma_tally3_int_t*            int_array_2;
    cusparseSolveAnalysisInfo_t cuinfo;
    cusparseSolveAnalysisInfo_t cuinfoL;
    cusparseSolveAnalysisInfo_t cuinfoU;
#if defined(HAVE_PASTIX)
    pastix_data_t*          pastix_data;
    magma_tally3_int_t*            iparm;
    float*                  dparm;
#endif

}magma_tally3_c_preconditioner;


typedef struct magma_tally3_d_preconditioner{

    magma_tally3_solver_type       solver;
    magma_tally3_int_t             levels;
    magma_tally3_int_t             sweeps;
    magma_tally3_precision         format;
    double                  epsilon;  
    magma_tally3_int_t             maxiter;
    magma_tally3_int_t             restart; 
    magma_tally3_int_t             numiter;
    double                  init_res;
    double                  final_res;
    magma_tally3_d_matrix   M;
    magma_tally3_d_matrix   L;
    magma_tally3_d_matrix   U;
    magma_tally3_d_matrix   LD;
    magma_tally3_d_matrix   UD;
    magma_tally3_d_matrix          d;
    magma_tally3_d_matrix          d2;
    magma_tally3_d_matrix          work1;
    magma_tally3_d_matrix          work2;
    magma_tally3_int_t*            int_array_1;
    magma_tally3_int_t*            int_array_2;
    cusparseSolveAnalysisInfo_t cuinfo;
    cusparseSolveAnalysisInfo_t cuinfoL;
    cusparseSolveAnalysisInfo_t cuinfoU;
#if defined(HAVE_PASTIX)
    pastix_data_t*          pastix_data;
    magma_tally3_int_t*            iparm;
    double*                 dparm;
#endif

}magma_tally3_d_preconditioner;


typedef struct magma_tally3_s_preconditioner{

    magma_tally3_solver_type       solver;
    magma_tally3_int_t             levels;
    magma_tally3_int_t             sweeps;
    magma_tally3_precision         format;
    float                   epsilon;  
    magma_tally3_int_t             maxiter;
    magma_tally3_int_t             restart; 
    magma_tally3_int_t             numiter;
    float                   init_res;
    float                   final_res;
    magma_tally3_s_matrix   M;
    magma_tally3_s_matrix   L;
    magma_tally3_s_matrix   U;
    magma_tally3_s_matrix   LD;
    magma_tally3_s_matrix   UD;
    magma_tally3_s_matrix          d;
    magma_tally3_s_matrix          d2;
    magma_tally3_s_matrix          work1;
    magma_tally3_s_matrix          work2;
    magma_tally3_int_t*            int_array_1;
    magma_tally3_int_t*            int_array_2;
    cusparseSolveAnalysisInfo_t cuinfo;
    cusparseSolveAnalysisInfo_t cuinfoL;
    cusparseSolveAnalysisInfo_t cuinfoU;
#if defined(HAVE_PASTIX)
    pastix_data_t*          pastix_data;
    magma_tally3_int_t*            iparm;
    float*                  dparm;
#endif

}magma_tally3_s_preconditioner;


//##############################################################################
//
//              opts for the testers
//
//##############################################################################

typedef struct magma_tally3_zopts{

    magma_tally3_z_solver_par      solver_par;
    magma_tally3_z_preconditioner  precond_par;
    magma_tally3_storage_t         input_format;
    int                     blocksize;
    int                     alignment;
    magma_tally3_storage_t         output_format;
    magma_tally3_location_t        input_location;
    magma_tally3_location_t        output_location;
    magma_tally3_scale_t           scaling;

}magma_tally3_zopts;

typedef struct magma_tally3_copts{

    magma_tally3_c_solver_par      solver_par;
    magma_tally3_c_preconditioner  precond_par;
    magma_tally3_storage_t         input_format;
    int                     blocksize;
    int                     alignment;
    magma_tally3_storage_t         output_format;
    magma_tally3_location_t        input_location;
    magma_tally3_location_t        output_location;
    magma_tally3_scale_t           scaling;

}magma_tally3_copts;

typedef struct magma_tally3_dopts{

    magma_tally3_d_solver_par      solver_par;
    magma_tally3_d_preconditioner  precond_par;
    magma_tally3_storage_t         input_format;
    int                     blocksize;
    int                     alignment;
    magma_tally3_storage_t         output_format;
    magma_tally3_location_t        input_location;
    magma_tally3_location_t        output_location;
    magma_tally3_scale_t           scaling;

}magma_tally3_dopts;

typedef struct magma_tally3_sopts{

    magma_tally3_s_solver_par      solver_par;
    magma_tally3_s_preconditioner  precond_par;
    magma_tally3_storage_t         input_format;
    int                     blocksize;
    int                     alignment;
    magma_tally3_storage_t         output_format;
    magma_tally3_location_t        input_location;
    magma_tally3_location_t        output_location;
    magma_tally3_scale_t           scaling;

}magma_tally3_sopts;

#ifdef __cplusplus
}
#endif

#endif        //  #ifndef MAGMA_tally3SPARSE_TYPES_H
