/*
    -- MAGMA_tally4 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from zgecsrmv.cu normal z -> c, Sun May  3 11:22:58 2015

*/
#include "common_magma_tally4.h"

#define BLOCK_SIZE 256


// CSR-SpMV kernel
__global__ void 
cgecsrmv_kernel( 
    int num_rows, 
    int num_cols, 
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

    if(row<num_rows){
        magma_tally4FloatComplex dot = MAGMA_tally4_C_ZERO;
        int start = drowptr[ row ];
        int end = drowptr[ row+1 ];
        for( j=start; j<end; j++)
            dot += dval[ j ] * dx[ dcolind[j] ];
        dy[ row ] =  dot *alpha + beta * dy[ row ];
    }
}

// shifted CSR-SpMV kernel
__global__ void 
cgecsrmv_kernel_shift( 
    int num_rows, 
    int num_cols, 
    magma_tally4FloatComplex alpha, 
    magma_tally4FloatComplex lambda, 
    magma_tally4FloatComplex * dval, 
    magma_tally4_index_t * drowptr, 
    magma_tally4_index_t * dcolind,
    magma_tally4FloatComplex * dx,
    magma_tally4FloatComplex beta, 
    int offset,
    int blocksize,
    magma_tally4_index_t * addrows,
    magma_tally4FloatComplex * dy)
{

    int row = blockIdx.x*blockDim.x+threadIdx.x;
    int j;

    if(row<num_rows){
        magma_tally4FloatComplex dot = MAGMA_tally4_C_ZERO;
        int start = drowptr[ row ];
        int end = drowptr[ row+1 ];
        for( j=start; j<end; j++)
            dot += dval[ j ] * dx[ dcolind[j] ];
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
    
    This routine computes y = alpha *  A *  x + beta * y on the GPU.
    The input format is CSR (val, row, col).
    
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
magma_tally4_cgecsrmv(
    magma_tally4_trans_t transA,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4FloatComplex alpha,
    magma_tally4FloatComplex_ptr dval,
    magma_tally4Index_ptr drowptr,
    magma_tally4Index_ptr dcolind,
    magma_tally4FloatComplex_ptr dx,
    magma_tally4FloatComplex beta,
    magma_tally4FloatComplex_ptr dy,
    magma_tally4_queue_t queue )
{
    dim3 grid( magma_tally4_ceildiv( m, BLOCK_SIZE ) );
    magma_tally4_int_t threads = BLOCK_SIZE;
    cgecsrmv_kernel<<< grid, threads, 0, queue >>>
                    (m, n, alpha, dval, drowptr, dcolind, dx, beta, dy);

    return MAGMA_tally4_SUCCESS;
}



/**
    Purpose
    -------
    
    This routine computes y = alpha * ( A -lambda I ) * x + beta * y on the GPU.
    It is a shifted version of the CSR-SpMV.
    
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
    alpha       magma_tally4FloatComplex
                scalar multiplier

    @param[in]
    lambda      magma_tally4FloatComplex
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
    dy          magma_tally4FloatComplex_ptr
                output vector y  
    @param[in]
    queue       magma_tally4_queue_t
                Queue to execute in.

    @ingroup magma_tally4sparse_cblas
    ********************************************************************/

extern "C" magma_tally4_int_t
magma_tally4_cgecsrmv_shift(
    magma_tally4_trans_t transA,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4FloatComplex alpha,
    magma_tally4FloatComplex lambda,
    magma_tally4FloatComplex_ptr dval,
    magma_tally4Index_ptr drowptr,
    magma_tally4Index_ptr dcolind,
    magma_tally4FloatComplex_ptr dx,
    magma_tally4FloatComplex beta,
    int offset,
    int blocksize,
    magma_tally4_index_t * addrows,
    magma_tally4FloatComplex_ptr dy,
    magma_tally4_queue_t queue )
{
    dim3 grid( magma_tally4_ceildiv( m, BLOCK_SIZE ) );
    magma_tally4_int_t threads = BLOCK_SIZE;
    cgecsrmv_kernel_shift<<< grid, threads, 0, queue >>>
                         (m, n, alpha, lambda, dval, drowptr, dcolind, dx, 
                                    beta, offset, blocksize, addrows, dy);

    return MAGMA_tally4_SUCCESS;
}



