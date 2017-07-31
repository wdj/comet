/*
    -- MAGMA_tally2 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @precisions normal z -> c d s

*/
#include "common_magma_tally2.h"

#define BLOCK_SIZE 512


__global__ void 
zmgecsrmv_kernel( 
    int num_rows, 
    int num_cols, 
    int num_vecs,
    magma_tally2DoubleComplex alpha, 
    magma_tally2DoubleComplex * dval, 
    magma_tally2_index_t * drowptr, 
    magma_tally2_index_t * dcolind,
    magma_tally2DoubleComplex * dx,
    magma_tally2DoubleComplex beta, 
    magma_tally2DoubleComplex * dy)
{

    int row = blockIdx.x*blockDim.x+threadIdx.x;
    int j;
    extern __shared__ magma_tally2DoubleComplex dot[];

    if( row<num_rows ){
        for( int i=0; i<num_vecs; i++ )
                dot[ threadIdx.x+ i*blockDim.x ] = MAGMA_tally2_Z_MAKE(0.0, 0.0);
        int start = drowptr[ row ] ;
        int end = drowptr[ row+1 ];
        for( j=start; j<end; j++ ){
            int col = dcolind [ j ];
            magma_tally2DoubleComplex val = dval[ j ];
            for( int i=0; i<num_vecs; i++ )
                dot[ threadIdx.x + i*blockDim.x ] += 
                                    val * dx[ col + i*num_cols ];
        }
        for( int i=0; i<num_vecs; i++ )
            dy[ row +i*num_cols ] = alpha * dot[ threadIdx.x + i*blockDim.x ] 
                                             + beta * dy[ row + i*num_cols ];
    }
}



/**
    Purpose
    -------
    
    This routine computes Y = alpha *  A *  X + beta * Y for X and Y sets of 
    num_vec vectors on the GPU. Input format is CSR. 
    
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
    alpha       magma_tally2DoubleComplex
                scalar multiplier

    @param[in]
    dval        magma_tally2DoubleComplex_ptr
                array containing values of A in CSR

    @param[in]
    drowptr     magma_tally2Index_ptr
                rowpointer of A in CSR

    @param[in]
    dcolind     magma_tally2Index_ptr
                columnindices of A in CSR

    @param[in]
    dx          magma_tally2DoubleComplex_ptr
                input vector x

    @param[in]
    beta        magma_tally2DoubleComplex
                scalar multiplier

    @param[out]
    dy          magma_tally2DoubleComplex_ptr
                input/output vector y

    @param[in]
    queue       magma_tally2_queue_t
                Queue to execute in.

    @ingroup magma_tally2sparse_zblas
    ********************************************************************/

extern "C" magma_tally2_int_t
magma_tally2_zmgecsrmv(
    magma_tally2_trans_t transA,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2_int_t num_vecs, 
    magma_tally2DoubleComplex alpha,
    magma_tally2DoubleComplex_ptr dval,
    magma_tally2Index_ptr drowptr,
    magma_tally2Index_ptr dcolind,
    magma_tally2DoubleComplex_ptr dx,
    magma_tally2DoubleComplex beta,
    magma_tally2DoubleComplex_ptr dy,
    magma_tally2_queue_t queue )
{
    dim3 grid( magma_tally2_ceildiv( m, BLOCK_SIZE ), 1, 1);
    magma_tally2_int_t threads = BLOCK_SIZE;
    unsigned int MEM_SIZE =  num_vecs* BLOCK_SIZE 
                    * sizeof( magma_tally2DoubleComplex ); // num_vecs vectors 
    zmgecsrmv_kernel<<< grid, threads, MEM_SIZE >>>
            (m, n, num_vecs, alpha, dval, drowptr, dcolind, dx, beta, dy);

   return MAGMA_tally2_SUCCESS;
}



