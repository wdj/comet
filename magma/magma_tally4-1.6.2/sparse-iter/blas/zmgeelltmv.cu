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


__global__ void 
zmgeelltmv_kernel( 
        int num_rows, 
        int num_cols,
        int num_vecs,
        int num_cols_per_row,
        magma_tally4DoubleComplex alpha, 
        magma_tally4DoubleComplex * dval, 
        magma_tally4_index_t * dcolind,
        magma_tally4DoubleComplex * dx,
        magma_tally4DoubleComplex beta, 
        magma_tally4DoubleComplex * dy)
{
    extern __shared__ magma_tally4DoubleComplex dot[];
    int row = blockDim.x * blockIdx.x + threadIdx.x ;
    if(row < num_rows ){
        for( int i=0; i<num_vecs; i++ )
                dot[ threadIdx.x+ i*blockDim.x ] = MAGMA_tally4_Z_MAKE(0.0, 0.0);
        for ( int n = 0; n < num_cols_per_row ; n ++){
            int col = dcolind [ num_rows * n + row ];
            magma_tally4DoubleComplex val = dval [ num_rows * n + row ];
            if( val != 0){
                for( int i=0; i<num_vecs; i++ )
                    dot[ threadIdx.x + i*blockDim.x ] += 
                                        val * dx[col + i * num_cols ];
            }
        }
        for( int i=0; i<num_vecs; i++ )
                dy[ row + i*num_cols ] = dot[ threadIdx.x + i*blockDim.x ] 
                                * alpha + beta * dy [ row + i*num_cols ];
    }
}





/**
    Purpose
    -------
    
    This routine computes Y = alpha *  A *  X + beta * Y for X and Y sets of 
    num_vec vectors on the GPU. Input format is ELL. 
    
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
    num_vecs    mama_int_t
                number of vectors
                
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

    @ingroup magma_tally4sparse_zblas
    ********************************************************************/

extern "C" magma_tally4_int_t
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
    magma_tally4_queue_t queue )
{
    dim3 grid( magma_tally4_ceildiv( m, BLOCK_SIZE ) );
    magma_tally4_int_t threads = BLOCK_SIZE;
    unsigned int MEM_SIZE =  num_vecs* BLOCK_SIZE 
                * sizeof( magma_tally4DoubleComplex ); // num_vecs vectors 
    zmgeelltmv_kernel<<< grid, threads, MEM_SIZE, queue >>>
        ( m, n, num_vecs, nnz_per_row, alpha, dval, dcolind, dx, beta, dy );


    return MAGMA_tally4_SUCCESS;
}



