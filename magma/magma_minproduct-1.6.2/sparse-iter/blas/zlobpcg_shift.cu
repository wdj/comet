/*
    -- MAGMA_minproduct (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @precisions normal z -> c d s

*/

#include "common_magma_minproduct.h"

__global__ void
magma_minproduct_zlobpcg_shift_kernel( 
    magma_minproduct_int_t num_rows, 
    magma_minproduct_int_t num_vecs, 
    magma_minproduct_int_t shift, 
    magma_minproductDoubleComplex * x )
{

    int idx = threadIdx.x ;     // thread in row
    int row = blockIdx.y * gridDim.x + blockIdx.x; // global block index

    if( row<num_rows){
        magma_minproductDoubleComplex tmp = x[idx];
        __syncthreads();

        if( idx > shift-1 ){
            idx-=shift;
            x[idx] = tmp;
            __syncthreads();
        }

    }
}




/**
    Purpose
    -------
    
    For a Block-LOBPCG, the set of residuals (entries consecutive in memory)  
    shrinks and the vectors are shifted in case shift residuals drop below 
    threshold. The memory layout of x is:

        / x1[0] x2[0] x3[0] \
        | x1[1] x2[1] x3[1] |
    x = | x1[2] x2[2] x3[2] | = x1[0] x2[0] x3[0] x1[1] x2[1] x3[1] x1[2] .
        | x1[3] x2[3] x3[3] |
        \ x1[4] x2[4] x3[4] /
    
    Arguments
    ---------

    @param[in]
    num_rows    magma_minproduct_int_t
                number of rows

    @param[in]
    num_vecs    magma_minproduct_int_t
                number of vectors

    @param[in]
    shift       magma_minproduct_int_t
                shift number

    @param[in/out]
    x           magma_minproductDoubleComplex_ptr 
                input/output vector x

    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_zaux
    ********************************************************************/

extern "C" magma_minproduct_int_t
magma_minproduct_zlobpcg_shift(
    magma_minproduct_int_t num_rows,
    magma_minproduct_int_t num_vecs, 
    magma_minproduct_int_t shift,
    magma_minproductDoubleComplex_ptr x,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t num_threads = num_vecs;
    // every thread handles one row containing the 
    if (  num_threads > 1024 )
        printf("error: too many threads requested.\n");

    int Ms = num_threads * sizeof( magma_minproductDoubleComplex );
    if (  Ms > 1024*8 )
        printf("error: too much shared memory requested.\n");

    dim3 block( num_threads, 1, 1 );

    int dimgrid1 = (int) sqrt( (double) num_rows);
    int dimgrid2 = magma_minproduct_ceildiv( num_rows, dimgrid1 );

    dim3 grid( dimgrid1, dimgrid2, 1);

    magma_minproduct_zlobpcg_shift_kernel<<< grid, block, Ms, queue >>>
            ( num_rows, num_vecs, shift, x );


    return MAGMA_minproduct_SUCCESS;
}



