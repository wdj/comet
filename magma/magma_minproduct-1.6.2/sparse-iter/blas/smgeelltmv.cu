/*
    -- MAGMA_minproduct (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from zmgeelltmv.cu normal z -> s, Sun May  3 11:22:58 2015

*/

#include "common_magma_minproduct.h"

#define BLOCK_SIZE 512


__global__ void 
smgeelltmv_kernel( 
        int num_rows, 
        int num_cols,
        int num_vecs,
        int num_cols_per_row,
        float alpha, 
        float * dval, 
        magma_minproduct_index_t * dcolind,
        float * dx,
        float beta, 
        float * dy)
{
    extern __shared__ float dot[];
    int row = blockDim.x * blockIdx.x + threadIdx.x ;
    if(row < num_rows ){
        for( int i=0; i<num_vecs; i++ )
                dot[ threadIdx.x+ i*blockDim.x ] = MAGMA_minproduct_S_MAKE(0.0, 0.0);
        for ( int n = 0; n < num_cols_per_row ; n ++){
            int col = dcolind [ num_rows * n + row ];
            float val = dval [ num_rows * n + row ];
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
    transA      magma_minproduct_trans_t
                transposition parameter for A

    @param[in]
    m           magma_minproduct_int_t
                number of rows in A

    @param[in]
    n           magma_minproduct_int_t
                number of columns in A 
                
    @param[in]
    num_vecs    mama_int_t
                number of vectors
                
    @param[in]
    nnz_per_row magma_minproduct_int_t
                number of elements in the longest row 
                
    @param[in]
    alpha       float
                scalar multiplier

    @param[in]
    dval        magma_minproductFloat_ptr
                array containing values of A in ELL

    @param[in]
    dcolind     magma_minproductIndex_ptr
                columnindices of A in ELL

    @param[in]
    dx          magma_minproductFloat_ptr
                input vector x

    @param[in]
    beta        float
                scalar multiplier

    @param[out]
    dy          magma_minproductFloat_ptr
                input/output vector y

    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_sblas
    ********************************************************************/

extern "C" magma_minproduct_int_t
magma_minproduct_smgeelltmv(
    magma_minproduct_trans_t transA,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproduct_int_t num_vecs,
    magma_minproduct_int_t nnz_per_row,
    float alpha,
    magma_minproductFloat_ptr dval,
    magma_minproductIndex_ptr dcolind,
    magma_minproductFloat_ptr dx,
    float beta,
    magma_minproductFloat_ptr dy,
    magma_minproduct_queue_t queue )
{
    dim3 grid( magma_minproduct_ceildiv( m, BLOCK_SIZE ) );
    magma_minproduct_int_t threads = BLOCK_SIZE;
    unsigned int MEM_SIZE =  num_vecs* BLOCK_SIZE 
                * sizeof( float ); // num_vecs vectors 
    smgeelltmv_kernel<<< grid, threads, MEM_SIZE, queue >>>
        ( m, n, num_vecs, nnz_per_row, alpha, dval, dcolind, dx, beta, dy );


    return MAGMA_minproduct_SUCCESS;
}



