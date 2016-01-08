/*
    -- MAGMA_tally4 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @precisions normal z -> c d s

*/

#include "common_magma_tally4.h"

// 512 is maximum number of threads for CUDA capability 1.x
#define BLOCK_SIZE  512



__global__ void
magma_tally4_zlobpcg_maxpy_kernel( 
    magma_tally4_int_t num_rows, 
    magma_tally4_int_t num_vecs, 
    magma_tally4DoubleComplex * X, 
    magma_tally4DoubleComplex * Y)
{

    int row = blockIdx.x * blockDim.x + threadIdx.x; // global row index

    if( row<num_rows ){
        for( int i=0; i<num_vecs; i++ ){ 

            Y[ row + i*num_rows ] += X[ row + i*num_rows ];
        }
    }
}




/**
    Purpose
    -------
    
    This routine computes a axpy for a mxn matrix:
        
        Y = X + Y
        
    It replaces:
            magma_tally4_zaxpy(m*n, c_one, Y, 1, X, 1);


        / x1[0] x2[0] x3[0] \
        | x1[1] x2[1] x3[1] |
    X = | x1[2] x2[2] x3[2] | = x1[0] x1[1] x1[2] x1[3] x1[4] x2[0] x2[1] .
        | x1[3] x2[3] x3[3] |
        \ x1[4] x2[4] x3[4] /
    
    Arguments
    ---------

    @param[in]
    num_rows    magma_tally4_int_t
                number of rows

    @param[in]
    num_vecs    magma_tally4_int_t
                number of vectors

    @param[in]
    X           magma_tally4DoubleComplex_ptr 
                input vector X

    @param[in/out]
    Y           magma_tally4DoubleComplex_ptr 
                input/output vector Y

    @param[in]
    queue       magma_tally4_queue_t
                Queue to execute in.

    @ingroup magma_tally4sparse_zgegpuk
    ********************************************************************/

extern "C" magma_tally4_int_t
magma_tally4_zlobpcg_maxpy(
    magma_tally4_int_t num_rows,
    magma_tally4_int_t num_vecs, 
    magma_tally4DoubleComplex_ptr X,
    magma_tally4DoubleComplex_ptr Y,
    magma_tally4_queue_t queue )
{
    // every thread handles one row

    magma_tally4_int_t block_size = BLOCK_SIZE;
     magma_tally4_int_t threads = BLOCK_SIZE;
    dim3 block( block_size );
    dim3 grid( magma_tally4_ceildiv( num_rows, block_size ) );

    magma_tally4_zlobpcg_maxpy_kernel<<< grid, threads, 0, queue >>>
                                ( num_rows, num_vecs, X, Y );


    return MAGMA_tally4_SUCCESS;
}



