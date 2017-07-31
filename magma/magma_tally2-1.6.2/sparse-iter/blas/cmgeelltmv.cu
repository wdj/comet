/*
    -- MAGMA_tally2 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from zmgeelltmv.cu normal z -> c, Sun May  3 11:22:58 2015

*/

#include "common_magma_tally2.h"

#define BLOCK_SIZE 512


__global__ void 
cmgeelltmv_kernel( 
        int num_rows, 
        int num_cols,
        int num_vecs,
        int num_cols_per_row,
        magma_tally2FloatComplex alpha, 
        magma_tally2FloatComplex * dval, 
        magma_tally2_index_t * dcolind,
        magma_tally2FloatComplex * dx,
        magma_tally2FloatComplex beta, 
        magma_tally2FloatComplex * dy)
{
    extern __shared__ magma_tally2FloatComplex dot[];
    int row = blockDim.x * blockIdx.x + threadIdx.x ;
    if(row < num_rows ){
        for( int i=0; i<num_vecs; i++ )
                dot[ threadIdx.x+ i*blockDim.x ] = MAGMA_tally2_C_MAKE(0.0, 0.0);
        for ( int n = 0; n < num_cols_per_row ; n ++){
            int col = dcolind [ num_rows * n + row ];
            magma_tally2FloatComplex val = dval [ num_rows * n + row ];
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
    transA      magma_tally2_trans_t
                transposition parameter for A

    @param[in]
    m           magma_tally2_int_t
                number of rows in A

    @param[in]
    n           magma_tally2_int_t
                number of columns in A 
                
    @param[in]
    num_vecs    mama_int_t
                number of vectors
                
    @param[in]
    nnz_per_row magma_tally2_int_t
                number of elements in the longest row 
                
    @param[in]
    alpha       magma_tally2FloatComplex
                scalar multiplier

    @param[in]
    dval        magma_tally2FloatComplex_ptr
                array containing values of A in ELL

    @param[in]
    dcolind     magma_tally2Index_ptr
                columnindices of A in ELL

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
magma_tally2_cmgeelltmv(
    magma_tally2_trans_t transA,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2_int_t num_vecs,
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
    unsigned int MEM_SIZE =  num_vecs* BLOCK_SIZE 
                * sizeof( magma_tally2FloatComplex ); // num_vecs vectors 
    cmgeelltmv_kernel<<< grid, threads, MEM_SIZE, queue >>>
        ( m, n, num_vecs, nnz_per_row, alpha, dval, dcolind, dx, beta, dy );


    return MAGMA_tally2_SUCCESS;
}



