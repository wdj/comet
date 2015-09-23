/*
    -- MAGMA_minproduct (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @precisions normal z -> c d s

*/

#include "common_magma_minproduct.h"

#define BLOCK_SIZE 512


__global__ void 
zmgeelltmv_kernel( 
        int num_rows, 
        int num_cols,
        int num_vecs,
        int num_cols_per_row,
        magma_minproductDoubleComplex alpha, 
        magma_minproductDoubleComplex * dval, 
        magma_minproduct_index_t * dcolind,
        magma_minproductDoubleComplex * dx,
        magma_minproductDoubleComplex beta, 
        magma_minproductDoubleComplex * dy)
{
    extern __shared__ magma_minproductDoubleComplex dot[];
    int row = blockDim.x * blockIdx.x + threadIdx.x ;
    if(row < num_rows ){
        for( int i=0; i<num_vecs; i++ )
                dot[ threadIdx.x+ i*blockDim.x ] = MAGMA_minproduct_Z_MAKE(0.0, 0.0);
        for ( int n = 0; n < num_cols_per_row ; n ++){
            int col = dcolind [ num_rows * n + row ];
            magma_minproductDoubleComplex val = dval [ num_rows * n + row ];
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
    alpha       magma_minproductDoubleComplex
                scalar multiplier

    @param[in]
    dval        magma_minproductDoubleComplex_ptr
                array containing values of A in ELL

    @param[in]
    dcolind     magma_minproductIndex_ptr
                columnindices of A in ELL

    @param[in]
    dx          magma_minproductDoubleComplex_ptr
                input vector x

    @param[in]
    beta        magma_minproductDoubleComplex
                scalar multiplier

    @param[out]
    dy          magma_minproductDoubleComplex_ptr
                input/output vector y

    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_zblas
    ********************************************************************/

extern "C" magma_minproduct_int_t
magma_minproduct_zmgeelltmv(
    magma_minproduct_trans_t transA,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproduct_int_t num_vecs,
    magma_minproduct_int_t nnz_per_row,
    magma_minproductDoubleComplex alpha,
    magma_minproductDoubleComplex_ptr dval,
    magma_minproductIndex_ptr dcolind,
    magma_minproductDoubleComplex_ptr dx,
    magma_minproductDoubleComplex beta,
    magma_minproductDoubleComplex_ptr dy,
    magma_minproduct_queue_t queue )
{
    dim3 grid( magma_minproduct_ceildiv( m, BLOCK_SIZE ) );
    magma_minproduct_int_t threads = BLOCK_SIZE;
    unsigned int MEM_SIZE =  num_vecs* BLOCK_SIZE 
                * sizeof( magma_minproductDoubleComplex ); // num_vecs vectors 
    zmgeelltmv_kernel<<< grid, threads, MEM_SIZE, queue >>>
        ( m, n, num_vecs, nnz_per_row, alpha, dval, dcolind, dx, beta, dy );


    return MAGMA_minproduct_SUCCESS;
}



