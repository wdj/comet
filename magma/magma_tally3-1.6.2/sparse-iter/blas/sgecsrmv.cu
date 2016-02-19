/*
    -- MAGMA_tally3 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from zgecsrmv.cu normal z -> s, Sun May  3 11:22:58 2015

*/
#include "common_magma_tally3.h"

#define BLOCK_SIZE 256


// CSR-SpMV kernel
__global__ void 
sgecsrmv_kernel( 
    int num_rows, 
    int num_cols, 
    float alpha, 
    float * dval, 
    magma_tally3_index_t * drowptr, 
    magma_tally3_index_t * dcolind,
    float * dx,
    float beta, 
    float * dy)
{

    int row = blockIdx.x*blockDim.x+threadIdx.x;
    int j;

    if(row<num_rows){
        float dot = MAGMA_tally3_S_ZERO;
        int start = drowptr[ row ];
        int end = drowptr[ row+1 ];
        for( j=start; j<end; j++)
            dot += dval[ j ] * dx[ dcolind[j] ];
        dy[ row ] =  dot *alpha + beta * dy[ row ];
    }
}

// shifted CSR-SpMV kernel
__global__ void 
sgecsrmv_kernel_shift( 
    int num_rows, 
    int num_cols, 
    float alpha, 
    float lambda, 
    float * dval, 
    magma_tally3_index_t * drowptr, 
    magma_tally3_index_t * dcolind,
    float * dx,
    float beta, 
    int offset,
    int blocksize,
    magma_tally3_index_t * addrows,
    float * dy)
{

    int row = blockIdx.x*blockDim.x+threadIdx.x;
    int j;

    if(row<num_rows){
        float dot = MAGMA_tally3_S_ZERO;
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
    transA      magma_tally3_trans_t
                transposition parameter for A
                
    @param[in]
    m           magma_tally3_int_t
                number of rows in A

    @param[in]
    n           magma_tally3_int_t
                number of columns in A 

    @param[in]
    alpha       float
                scalar multiplier

    @param[in]
    dval        magma_tally3Float_ptr
                array containing values of A in CSR

    @param[in]
    drowptr     magma_tally3Index_ptr
                rowpointer of A in CSR

    @param[in]
    dcolind     magma_tally3Index_ptr
                columnindices of A in CSR

    @param[in]
    dx          magma_tally3Float_ptr
                input vector x

    @param[in]
    beta        float
                scalar multiplier

    @param[out]
    dy          magma_tally3Float_ptr
                input/output vector y

    @param[in]
    queue       magma_tally3_queue_t
                Queue to execute in.

    @ingroup magma_tally3sparse_sblas
    ********************************************************************/

extern "C" magma_tally3_int_t
magma_tally3_sgecsrmv(
    magma_tally3_trans_t transA,
    magma_tally3_int_t m, magma_tally3_int_t n,
    float alpha,
    magma_tally3Float_ptr dval,
    magma_tally3Index_ptr drowptr,
    magma_tally3Index_ptr dcolind,
    magma_tally3Float_ptr dx,
    float beta,
    magma_tally3Float_ptr dy,
    magma_tally3_queue_t queue )
{
    dim3 grid( magma_tally3_ceildiv( m, BLOCK_SIZE ) );
    magma_tally3_int_t threads = BLOCK_SIZE;
    sgecsrmv_kernel<<< grid, threads, 0, queue >>>
                    (m, n, alpha, dval, drowptr, dcolind, dx, beta, dy);

    return MAGMA_tally3_SUCCESS;
}



/**
    Purpose
    -------
    
    This routine computes y = alpha * ( A -lambda I ) * x + beta * y on the GPU.
    It is a shifted version of the CSR-SpMV.
    
    Arguments
    ---------
    
    @param[in]
    transA      magma_tally3_trans_t
                transposition parameter for A

    @param[in]
    m           magma_tally3_int_t
                number of rows in A

    @param[in]
    n           magma_tally3_int_t
                number of columns in A 

    @param[in]
    alpha       float
                scalar multiplier

    @param[in]
    lambda      float
                scalar multiplier

    @param[in]
    dval        magma_tally3Float_ptr
                array containing values of A in CSR

    @param[in]
    drowptr     magma_tally3Index_ptr
                rowpointer of A in CSR

    @param[in]
    dcolind     magma_tally3Index_ptr
                columnindices of A in CSR

    @param[in]
    dx          magma_tally3Float_ptr
                input vector x

    @param[in]
    beta        float
                scalar multiplier

    @param[in]
    offset      magma_tally3_int_t 
                in case not the main diagonal is scaled
                
    @param[in]
    blocksize   magma_tally3_int_t 
                in case of processing multiple vectors  
                
    @param[in]
    addrows     magma_tally3Index_ptr
                in case the matrixpowerskernel is used
                
    @param[out]
    dy          magma_tally3Float_ptr
                output vector y  
    @param[in]
    queue       magma_tally3_queue_t
                Queue to execute in.

    @ingroup magma_tally3sparse_sblas
    ********************************************************************/

extern "C" magma_tally3_int_t
magma_tally3_sgecsrmv_shift(
    magma_tally3_trans_t transA,
    magma_tally3_int_t m, magma_tally3_int_t n,
    float alpha,
    float lambda,
    magma_tally3Float_ptr dval,
    magma_tally3Index_ptr drowptr,
    magma_tally3Index_ptr dcolind,
    magma_tally3Float_ptr dx,
    float beta,
    int offset,
    int blocksize,
    magma_tally3_index_t * addrows,
    magma_tally3Float_ptr dy,
    magma_tally3_queue_t queue )
{
    dim3 grid( magma_tally3_ceildiv( m, BLOCK_SIZE ) );
    magma_tally3_int_t threads = BLOCK_SIZE;
    sgecsrmv_kernel_shift<<< grid, threads, 0, queue >>>
                         (m, n, alpha, lambda, dval, drowptr, dcolind, dx, 
                                    beta, offset, blocksize, addrows, dy);

    return MAGMA_tally3_SUCCESS;
}



