/*
    -- MAGMA_tally4 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from zmgecsrmv.cu normal z -> c, Sun May  3 11:22:58 2015

*/
#include "common_magma_tally4.h"

#define BLOCK_SIZE 512


__global__ void 
cmgecsrmv_kernel( 
    int num_rows, 
    int num_cols, 
    int num_vecs,
    magma_tally4FloatComplex alpha, 
    magma_tally4FloatComplex * dval, 
    magma_tally4_index_t * drowptr, 
    magma_tally4_index_t * dcolind,
    magma_tally4FloatComplex * dx,
    magma_tally4FloatComplex beta, 
    magma_tally4FloatComplex * dy)
{

    int row = blockIdx.x*blockDim.x+threadIdx.x;
    int j;
    extern __shared__ magma_tally4FloatComplex dot[];

    if( row<num_rows ){
        for( int i=0; i<num_vecs; i++ )
                dot[ threadIdx.x+ i*blockDim.x ] = MAGMA_tally4_C_MAKE(0.0, 0.0);
        int start = drowptr[ row ] ;
        int end = drowptr[ row+1 ];
        for( j=start; j<end; j++ ){
            int col = dcolind [ j ];
            magma_tally4FloatComplex val = dval[ j ];
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
    alpha       magma_tally4FloatComplex
                scalar multiplier

    @param[in]
    dval        magma_tally4FloatComplex_ptr
                array containing values of A in CSR

    @param[in]
    drowptr     magma_tally4Index_ptr
                rowpointer of A in CSR

    @param[in]
    dcolind     magma_tally4Index_ptr
                columnindices of A in CSR

    @param[in]
    dx          magma_tally4FloatComplex_ptr
                input vector x

    @param[in]
    beta        magma_tally4FloatComplex
                scalar multiplier

    @param[out]
    dy          magma_tally4FloatComplex_ptr
                input/output vector y

    @param[in]
    queue       magma_tally4_queue_t
                Queue to execute in.

    @ingroup magma_tally4sparse_cblas
    ********************************************************************/

extern "C" magma_tally4_int_t
magma_tally4_cmgecsrmv(
    magma_tally4_trans_t transA,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4_int_t num_vecs, 
    magma_tally4FloatComplex alpha,
    magma_tally4FloatComplex_ptr dval,
    magma_tally4Index_ptr drowptr,
    magma_tally4Index_ptr dcolind,
    magma_tally4FloatComplex_ptr dx,
    magma_tally4FloatComplex beta,
    magma_tally4FloatComplex_ptr dy,
    magma_tally4_queue_t queue )
{
    dim3 grid( magma_tally4_ceildiv( m, BLOCK_SIZE ), 1, 1);
    magma_tally4_int_t threads = BLOCK_SIZE;
    unsigned int MEM_SIZE =  num_vecs* BLOCK_SIZE 
                    * sizeof( magma_tally4FloatComplex ); // num_vecs vectors 
    cmgecsrmv_kernel<<< grid, threads, MEM_SIZE >>>
            (m, n, num_vecs, alpha, dval, drowptr, dcolind, dx, beta, dy);

   return MAGMA_tally4_SUCCESS;
}



