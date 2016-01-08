/*
    -- MAGMA_tally4 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @precisions normal z -> c d s

*/

#include "common_magma_tally4.h"

#define BLOCK_SIZE 512


// ELL SpMV kernel
//Michael Garland
__global__ void 
zgeelltmv_kernel( 
    int num_rows, 
    int num_cols,
    int num_cols_per_row,
    magma_tally4DoubleComplex alpha, 
    magma_tally4DoubleComplex * dval, 
    magma_tally4_index_t * dcolind,
    magma_tally4DoubleComplex * dx,
    magma_tally4DoubleComplex beta, 
    magma_tally4DoubleComplex * dy)
{
    int row = blockDim.x * blockIdx.x + threadIdx.x ;
    if(row < num_rows ){
        magma_tally4DoubleComplex dot = MAGMA_tally4_Z_MAKE(0.0, 0.0);
        for ( int n = 0; n < num_cols_per_row ; n ++){
            int col = dcolind [ num_rows * n + row ];
            magma_tally4DoubleComplex val = dval [ num_rows * n + row ];
            if( val != 0)
                dot += val * dx[col ];
        }
        dy[ row ] = dot * alpha + beta * dy [ row ];
    }
}

// shifted ELL SpMV kernel
//Michael Garland
__global__ void 
zgeelltmv_kernel_shift( 
    int num_rows, 
    int num_cols,
    int num_cols_per_row,
    magma_tally4DoubleComplex alpha, 
    magma_tally4DoubleComplex lambda, 
    magma_tally4DoubleComplex * dval, 
    magma_tally4_index_t * dcolind,
    magma_tally4DoubleComplex * dx,
    magma_tally4DoubleComplex beta, 
    int offset,
    int blocksize,
    magma_tally4_index_t * addrows,
    magma_tally4DoubleComplex * dy)
{

    int row = blockDim.x * blockIdx.x + threadIdx.x ;
    if(row < num_rows ){
        magma_tally4DoubleComplex dot = MAGMA_tally4_Z_MAKE(0.0, 0.0);
        for ( int n = 0; n < num_cols_per_row ; n ++){
            int col = dcolind [ num_rows * n + row ];
            magma_tally4DoubleComplex val = dval [ num_rows * n + row ];
            if( val != 0)
                dot += val * dx[col ];
        }
        if( row<blocksize )
            dy[ row ] = dot * alpha - lambda 
                    * dx[ offset+row ] + beta * dy [ row ];
        else
            dy[ row ] = dot * alpha - lambda 
                    * dx[ addrows[row-blocksize] ] + beta * dy [ row ];            
    }
}




/**
    Purpose
    -------
    
    This routine computes y = alpha *  A^t *  x + beta * y on the GPU.
    Input format is ELL.
    
    Arguments
    ---------
    
    @param[in]
    transA      magma_tally4_trans_t
                transposition parameter for A
                
    @param[in]
    m           magma_tally4_int_t
                number of rows in A

    @param[in]
    n           magma_tally4_int_t
                number of columns in A 
                
    @param[in]
    nnz_per_row magma_tally4_int_t
                number of elements in the longest row 

    @param[in]
    alpha       magma_tally4DoubleComplex
                scalar multiplier

    @param[in]
    dval        magma_tally4DoubleComplex_ptr
                array containing values of A in ELL

    @param[in]
    dcolind     magma_tally4Index_ptr
                columnindices of A in ELL

    @param[in]
    dx          magma_tally4DoubleComplex_ptr
                input vector x

    @param[in]
    beta        magma_tally4DoubleComplex
                scalar multiplier

    @param[out]
    dy          magma_tally4DoubleComplex_ptr
                input/output vector y

    @param[in]
    queue       magma_tally4_queue_t
                Queue to execute in.

    @ingroup magma_tally4sparse_z
    ********************************************************************/

extern "C" magma_tally4_int_t
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
    magma_tally4_queue_t queue )
{
    dim3 grid( magma_tally4_ceildiv( m, BLOCK_SIZE ) );
    magma_tally4_int_t threads = BLOCK_SIZE;
    zgeelltmv_kernel<<< grid, threads, 0, queue >>>
                  ( m, n, nnz_per_row, alpha, dval, dcolind, dx, beta, dy );


   return MAGMA_tally4_SUCCESS;
}


/**
    Purpose
    -------
    
    This routine computes y = alpha *( A - lambda I ) * x + beta * y on the GPU.
    Input format is ELL.
    
    Arguments
    ---------

    @param[in]
    transA      magma_tally4_trans_t
                transposition parameter for A    

    @param[in]
    m           magma_tally4_int_t
                number of rows in A

    @param[in]
    n           magma_tally4_int_t
                number of columns in A 
                
    @param[in]
    nnz_per_row magma_tally4_int_t
                number of elements in the longest row 

    @param[in]
    alpha       magma_tally4DoubleComplex
                scalar multiplier

    @param[in]
    lambda      magma_tally4DoubleComplex
                scalar multiplier

    @param[in]
    dval        magma_tally4DoubleComplex_ptr
                array containing values of A in ELL

    @param[in]
    dcolind     magma_tally4Index_ptr
                columnindices of A in ELL

    @param[in]
    dx          magma_tally4DoubleComplex_ptr
                input vector x

    @param[in]
    beta        magma_tally4DoubleComplex
                scalar multiplier
                
    @param[in]
    offset      magma_tally4_int_t 
                in case not the main diagonal is scaled
                
    @param[in]
    blocksize   magma_tally4_int_t 
                in case of processing multiple vectors  
                
    @param[in]
    addrows     magma_tally4Index_ptr
                in case the matrixpowerskernel is used

    @param[out]
    dy          magma_tally4DoubleComplex_ptr
                input/output vector y

    @param[in]
    queue       magma_tally4_queue_t
                Queue to execute in.

    @ingroup magma_tally4sparse_zblas
    ********************************************************************/

extern "C" magma_tally4_int_t
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
    magma_tally4Index_ptr addrows,
    magma_tally4DoubleComplex_ptr dy,
    magma_tally4_queue_t queue )
{
    dim3 grid( magma_tally4_ceildiv( m, BLOCK_SIZE ) );
    magma_tally4_int_t threads = BLOCK_SIZE;
    magma_tally4DoubleComplex tmp_shift;
    //magma_tally4_zsetvector(1,&lambda,1,&tmp_shift,1); 
    tmp_shift = lambda;
    zgeelltmv_kernel_shift<<< grid, threads, 0, queue >>>
                  ( m, n, nnz_per_row, alpha, tmp_shift, dval, dcolind, dx, 
                            beta, offset, blocksize, addrows, dy );


   return MAGMA_tally4_SUCCESS;
}



