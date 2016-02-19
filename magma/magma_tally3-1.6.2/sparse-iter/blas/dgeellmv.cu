/*
    -- MAGMA_tally3 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from zgeellmv.cu normal z -> d, Sun May  3 11:22:58 2015

*/

#include "common_magma_tally3.h"

#define BLOCK_SIZE 512


// ELLPACK SpMV kernel
//Michael Garland
__global__ void 
dgeellmv_kernel( 
    int num_rows, 
    int num_cols,
    int num_cols_per_row,
    double alpha, 
    double * dval, 
    magma_tally3_index_t * dcolind,
    double * dx,
    double beta, 
    double * dy)
{
int row = blockDim.x * blockIdx.x + threadIdx.x ;
    if(row < num_rows ){
        double dot = MAGMA_tally3_D_MAKE(0.0, 0.0);
        for ( int n = 0; n < num_cols_per_row ; n ++){
            int col = dcolind [ num_cols_per_row * row + n ];
            double val = dval [ num_cols_per_row * row + n ];
            if( val != 0)
                dot += val * dx[col ];
        }
        dy[ row ] = dot * alpha + beta * dy [ row ];
    }
}

// shifted ELLPACK SpMV kernel
//Michael Garland
__global__ void 
dgeellmv_kernel_shift( 
    int num_rows, 
    int num_cols,
    int num_cols_per_row,
    double alpha, 
    double lambda, 
    double * dval, 
    magma_tally3_index_t * dcolind,
    double * dx,
    double beta, 
    int offset,
    int blocksize,
    magma_tally3_index_t * addrows,
    double * dy)
{
int row = blockDim.x * blockIdx.x + threadIdx.x ;
    if(row < num_rows ){
        double dot = MAGMA_tally3_D_MAKE(0.0, 0.0);
        for ( int n = 0; n < num_cols_per_row ; n ++){
            int col = dcolind [ num_cols_per_row * row + n ];
            double val = dval [ num_cols_per_row * row + n ];
            if( val != 0)
                dot += val * dx[col ];
        }
        if( row<blocksize )
            dy[ row ] = dot * alpha - lambda * dx[ offset+row ] + beta * dy [ row ];
        else
            dy[ row ] = dot * alpha - lambda * dx[ addrows[row-blocksize] ] + beta * dy [ row ];   
    }
}





/**
    Purpose
    -------
    
    This routine computes y = alpha *  A *  x + beta * y on the GPU.
    Input format is ELLPACK.
    
    Arguments
    ---------
    
    @param[in]
    transA      magma_tally3_trans_t
                transposition parameter for A
                
    @param[in]
    m           magma_tally3_int_t
                number of rows in A

    @param[in]
    n           magma_tally3_int_t
                number of columns in A 
                
    @param[in]
    nnz_per_row magma_tally3_int_t
                number of elements in the longest row 

    @param[in]
    alpha       double
                scalar multiplier

    @param[in]
    dval        magma_tally3Double_ptr
                array containing values of A in ELLPACK

    @param[in]
    dcolind     magma_tally3Index_ptr
                columnindices of A in ELLPACK

    @param[in]
    dx          magma_tally3Double_ptr
                input vector x

    @param[in]
    beta        double
                scalar multiplier

    @param[out]
    dy          magma_tally3Double_ptr
                input/output vector y

    @param[in]
    queue       magma_tally3_queue_t
                Queue to execute in.

    @ingroup magma_tally3sparse_dblas
    ********************************************************************/

extern "C" magma_tally3_int_t
magma_tally3_dgeellmv(
    magma_tally3_trans_t transA,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3_int_t nnz_per_row,
    double alpha,
    magma_tally3Double_ptr dval,
    magma_tally3Index_ptr dcolind,
    magma_tally3Double_ptr dx,
    double beta,
    magma_tally3Double_ptr dy,
    magma_tally3_queue_t queue )
{
    dim3 grid( magma_tally3_ceildiv( m, BLOCK_SIZE ) );
    magma_tally3_int_t threads = BLOCK_SIZE;
   dgeellmv_kernel<<< grid, threads, 0, queue >>>
                  ( m, n, nnz_per_row, alpha, dval, dcolind, dx, beta, dy );


   return MAGMA_tally3_SUCCESS;
}



/**
    Purpose
    -------
    
    This routine computes y = alpha *( A - lambda I ) * x + beta * y on the GPU.
    Input format is ELLPACK.
    It is the shifted version of the ELLPACK SpMV.
    
    Arguments
    ---------
    
    @param[in]
    transA      magma_tally3_trans_t
                transposition parameter for A

    @param[in]
    m           magma_tally3_int_t
                number of rows in A

    @param[in]
    n           magma_tally3_int_t
                number of columns in A 
    @param[in]
    nnz_per_row magma_tally3_int_t
                number of elements in the longest row 
                
    @param[in]
    alpha       double
                scalar multiplier
                
    @param[in]
    lambda      double
                scalar multiplier

    @param[in]
    dval        magma_tally3Double_ptr
                array containing values of A in ELLPACK

    @param[in]
    dcolind     magma_tally3Index_ptr
                columnindices of A in ELLPACK

    @param[in]
    dx          magma_tally3Double_ptr
                input vector x

    @param[in]
    beta        double
                scalar multiplier
                
    @param[in]
    offset      magma_tally3_int_t 
                in case not the main diagonal is scaled
                
    @param[in]
    blocksize   magma_tally3_int_t 
                in case of processing multiple vectors  
                
    @param[in]
    addrows     magma_tally3Index_ptr
                in case the matrixpowerskernel is used

    @param[out]
    dy          magma_tally3Double_ptr
                input/output vector y

    @param[in]
    queue       magma_tally3_queue_t
                Queue to execute in.

    @ingroup magma_tally3sparse_dblas
    ********************************************************************/

extern "C" magma_tally3_int_t
magma_tally3_dgeellmv_shift(
    magma_tally3_trans_t transA,
    magma_tally3_int_t m, magma_tally3_int_t n,
    magma_tally3_int_t nnz_per_row,
    double alpha,
    double lambda,
    magma_tally3Double_ptr dval,
    magma_tally3Index_ptr dcolind,
    magma_tally3Double_ptr dx,
    double beta,
    int offset,
    int blocksize,
    magma_tally3Index_ptr addrows,
    magma_tally3Double_ptr dy,
    magma_tally3_queue_t queue )
{
    dim3 grid( magma_tally3_ceildiv( m, BLOCK_SIZE ) );
    magma_tally3_int_t threads = BLOCK_SIZE;
   dgeellmv_kernel_shift<<< grid, threads, 0, queue >>>
                  ( m, n, nnz_per_row, alpha, lambda, dval, dcolind, dx, 
                                    beta, offset, blocksize, addrows, dy );


   return MAGMA_tally3_SUCCESS;
}



