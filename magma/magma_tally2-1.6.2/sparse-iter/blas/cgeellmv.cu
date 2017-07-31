/*
    -- MAGMA_tally2 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from zgeellmv.cu normal z -> c, Sun May  3 11:22:58 2015

*/

#include "common_magma_tally2.h"

#define BLOCK_SIZE 512


// ELLPACK SpMV kernel
//Michael Garland
__global__ void 
cgeellmv_kernel( 
    int num_rows, 
    int num_cols,
    int num_cols_per_row,
    magma_tally2FloatComplex alpha, 
    magma_tally2FloatComplex * dval, 
    magma_tally2_index_t * dcolind,
    magma_tally2FloatComplex * dx,
    magma_tally2FloatComplex beta, 
    magma_tally2FloatComplex * dy)
{
int row = blockDim.x * blockIdx.x + threadIdx.x ;
    if(row < num_rows ){
        magma_tally2FloatComplex dot = MAGMA_tally2_C_MAKE(0.0, 0.0);
        for ( int n = 0; n < num_cols_per_row ; n ++){
            int col = dcolind [ num_cols_per_row * row + n ];
            magma_tally2FloatComplex val = dval [ num_cols_per_row * row + n ];
            if( val != 0)
                dot += val * dx[col ];
        }
        dy[ row ] = dot * alpha + beta * dy [ row ];
    }
}

// shifted ELLPACK SpMV kernel
//Michael Garland
__global__ void 
cgeellmv_kernel_shift( 
    int num_rows, 
    int num_cols,
    int num_cols_per_row,
    magma_tally2FloatComplex alpha, 
    magma_tally2FloatComplex lambda, 
    magma_tally2FloatComplex * dval, 
    magma_tally2_index_t * dcolind,
    magma_tally2FloatComplex * dx,
    magma_tally2FloatComplex beta, 
    int offset,
    int blocksize,
    magma_tally2_index_t * addrows,
    magma_tally2FloatComplex * dy)
{
int row = blockDim.x * blockIdx.x + threadIdx.x ;
    if(row < num_rows ){
        magma_tally2FloatComplex dot = MAGMA_tally2_C_MAKE(0.0, 0.0);
        for ( int n = 0; n < num_cols_per_row ; n ++){
            int col = dcolind [ num_cols_per_row * row + n ];
            magma_tally2FloatComplex val = dval [ num_cols_per_row * row + n ];
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
    transA      magma_tally2_trans_t
                transposition parameter for A
                
    @param[in]
    m           magma_tally2_int_t
                number of rows in A

    @param[in]
    n           magma_tally2_int_t
                number of columns in A 
                
    @param[in]
    nnz_per_row magma_tally2_int_t
                number of elements in the longest row 

    @param[in]
    alpha       magma_tally2FloatComplex
                scalar multiplier

    @param[in]
    dval        magma_tally2FloatComplex_ptr
                array containing values of A in ELLPACK

    @param[in]
    dcolind     magma_tally2Index_ptr
                columnindices of A in ELLPACK

    @param[in]
    dx          magma_tally2FloatComplex_ptr
                input vector x

    @param[in]
    beta        magma_tally2FloatComplex
                scalar multiplier

    @param[out]
    dy          magma_tally2FloatComplex_ptr
                input/output vector y

    @param[in]
    queue       magma_tally2_queue_t
                Queue to execute in.

    @ingroup magma_tally2sparse_cblas
    ********************************************************************/

extern "C" magma_tally2_int_t
magma_tally2_cgeellmv(
    magma_tally2_trans_t transA,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2_int_t nnz_per_row,
    magma_tally2FloatComplex alpha,
    magma_tally2FloatComplex_ptr dval,
    magma_tally2Index_ptr dcolind,
    magma_tally2FloatComplex_ptr dx,
    magma_tally2FloatComplex beta,
    magma_tally2FloatComplex_ptr dy,
    magma_tally2_queue_t queue )
{
    dim3 grid( magma_tally2_ceildiv( m, BLOCK_SIZE ) );
    magma_tally2_int_t threads = BLOCK_SIZE;
   cgeellmv_kernel<<< grid, threads, 0, queue >>>
                  ( m, n, nnz_per_row, alpha, dval, dcolind, dx, beta, dy );


   return MAGMA_tally2_SUCCESS;
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
    transA      magma_tally2_trans_t
                transposition parameter for A

    @param[in]
    m           magma_tally2_int_t
                number of rows in A

    @param[in]
    n           magma_tally2_int_t
                number of columns in A 
    @param[in]
    nnz_per_row magma_tally2_int_t
                number of elements in the longest row 
                
    @param[in]
    alpha       magma_tally2FloatComplex
                scalar multiplier
                
    @param[in]
    lambda      magma_tally2FloatComplex
                scalar multiplier

    @param[in]
    dval        magma_tally2FloatComplex_ptr
                array containing values of A in ELLPACK

    @param[in]
    dcolind     magma_tally2Index_ptr
                columnindices of A in ELLPACK

    @param[in]
    dx          magma_tally2FloatComplex_ptr
                input vector x

    @param[in]
    beta        magma_tally2FloatComplex
                scalar multiplier
                
    @param[in]
    offset      magma_tally2_int_t 
                in case not the main diagonal is scaled
                
    @param[in]
    blocksize   magma_tally2_int_t 
                in case of processing multiple vectors  
                
    @param[in]
    addrows     magma_tally2Index_ptr
                in case the matrixpowerskernel is used

    @param[out]
    dy          magma_tally2FloatComplex_ptr
                input/output vector y

    @param[in]
    queue       magma_tally2_queue_t
                Queue to execute in.

    @ingroup magma_tally2sparse_cblas
    ********************************************************************/

extern "C" magma_tally2_int_t
magma_tally2_cgeellmv_shift(
    magma_tally2_trans_t transA,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2_int_t nnz_per_row,
    magma_tally2FloatComplex alpha,
    magma_tally2FloatComplex lambda,
    magma_tally2FloatComplex_ptr dval,
    magma_tally2Index_ptr dcolind,
    magma_tally2FloatComplex_ptr dx,
    magma_tally2FloatComplex beta,
    int offset,
    int blocksize,
    magma_tally2Index_ptr addrows,
    magma_tally2FloatComplex_ptr dy,
    magma_tally2_queue_t queue )
{
    dim3 grid( magma_tally2_ceildiv( m, BLOCK_SIZE ) );
    magma_tally2_int_t threads = BLOCK_SIZE;
   cgeellmv_kernel_shift<<< grid, threads, 0, queue >>>
                  ( m, n, nnz_per_row, alpha, lambda, dval, dcolind, dx, 
                                    beta, offset, blocksize, addrows, dy );


   return MAGMA_tally2_SUCCESS;
}



