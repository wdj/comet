/*
    -- MAGMA_tally3 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @precisions normal z -> c d s

*/

#include "common_magma_tally3.h"

#define BLOCK_SIZE 512


// ELLPACK SpMV kernel
//Michael Garland
__global__ void 
zgeellmv_kernel( 
    int num_rows, 
    int num_cols,
    int num_cols_per_row,
    magma_tally3DoubleComplex alpha, 
    magma_tally3DoubleComplex * dval, 
    magma_tally3_index_t * dcolind,
    magma_tally3DoubleComplex * dx,
    magma_tally3DoubleComplex beta, 
    magma_tally3DoubleComplex * dy)
{
int row = blockDim.x * blockIdx.x + threadIdx.x ;
    if(row < num_rows ){
        magma_tally3DoubleComplex dot = MAGMA_tally3_Z_MAKE(0.0, 0.0);
        for ( int n = 0; n < num_cols_per_row ; n ++){
            int col = dcolind [ num_cols_per_row * row + n ];
            magma_tally3DoubleComplex val = dval [ num_cols_per_row * row + n ];
            if( val != 0)
                dot += val * dx[col ];
        }
        dy[ row ] = dot * alpha + beta * dy [ row ];
    }
}

// shifted ELLPACK SpMV kernel
//Michael Garland
__global__ void 
zgeellmv_kernel_shift( 
    int num_rows, 
    int num_cols,
    int num_cols_per_row,
    magma_tally3DoubleComplex alpha, 
    magma_tally3DoubleComplex lambda, 
    magma_tally3DoubleComplex * dval, 
    magma_tally3_index_t * dcolind,
    magma_tally3DoubleComplex * dx,
    magma_tally3DoubleComplex beta, 
    int offset,
    int blocksize,
    magma_tally3_index_t * addrows,
    magma_tally3DoubleComplex * dy)
{
int row = blockDim.x * blockIdx.x + threadIdx.x ;
    if(row < num_rows ){
        magma_tally3DoubleComplex dot = MAGMA_tally3_Z_MAKE(0.0, 0.0);
        for ( int n = 0; n < num_cols_per_row ; n ++){
            int col = dcolind [ num_cols_per_row * row + n ];
            magma_tally3DoubleComplex val = dval [ num_cols_per_row * row + n ];
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
    alpha       magma_tally3DoubleComplex
                scalar multiplier

    @param[in]
    dval        magma_tally3DoubleComplex_ptr
                array containing values of A in ELLPACK

    @param[in]
    dcolind     magma_tally3Index_ptr
                columnindices of A in ELLPACK

    @param[in]
    dx          magma_tally3DoubleComplex_ptr
                input vector x

    @param[in]
    beta        magma_tally3DoubleComplex
                scalar multiplier

    @param[out]
    dy          magma_tally3DoubleComplex_ptr
                input/output vector y

    @param[in]
    queue       magma_tally3_queue_t
                Queue to execute in.

    @ingroup magma_tally3sparse_zblas
    ********************************************************************/

extern "C" magma_tally3_int_t
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
    magma_tally3_queue_t queue )
{
    dim3 grid( magma_tally3_ceildiv( m, BLOCK_SIZE ) );
    magma_tally3_int_t threads = BLOCK_SIZE;
   zgeellmv_kernel<<< grid, threads, 0, queue >>>
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
    alpha       magma_tally3DoubleComplex
                scalar multiplier
                
    @param[in]
    lambda      magma_tally3DoubleComplex
                scalar multiplier

    @param[in]
    dval        magma_tally3DoubleComplex_ptr
                array containing values of A in ELLPACK

    @param[in]
    dcolind     magma_tally3Index_ptr
                columnindices of A in ELLPACK

    @param[in]
    dx          magma_tally3DoubleComplex_ptr
                input vector x

    @param[in]
    beta        magma_tally3DoubleComplex
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
    dy          magma_tally3DoubleComplex_ptr
                input/output vector y

    @param[in]
    queue       magma_tally3_queue_t
                Queue to execute in.

    @ingroup magma_tally3sparse_zblas
    ********************************************************************/

extern "C" magma_tally3_int_t
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
    magma_tally3Index_ptr addrows,
    magma_tally3DoubleComplex_ptr dy,
    magma_tally3_queue_t queue )
{
    dim3 grid( magma_tally3_ceildiv( m, BLOCK_SIZE ) );
    magma_tally3_int_t threads = BLOCK_SIZE;
   zgeellmv_kernel_shift<<< grid, threads, 0, queue >>>
                  ( m, n, nnz_per_row, alpha, lambda, dval, dcolind, dx, 
                                    beta, offset, blocksize, addrows, dy );


   return MAGMA_tally3_SUCCESS;
}



