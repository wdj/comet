/*
    -- MAGMA_minproduct (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from zgecsrmv.cu normal z -> d, Sun May  3 11:22:58 2015

*/
#include "common_magma_minproduct.h"

#define BLOCK_SIZE 256


// CSR-SpMV kernel
__global__ void 
dgecsrmv_kernel( 
    int num_rows, 
    int num_cols, 
    double alpha, 
    double * dval, 
    magma_minproduct_index_t * drowptr, 
    magma_minproduct_index_t * dcolind,
    double * dx,
    double beta, 
    double * dy)
{

    int row = blockIdx.x*blockDim.x+threadIdx.x;
    int j;

    if(row<num_rows){
        double dot = MAGMA_minproduct_D_ZERO;
        int start = drowptr[ row ];
        int end = drowptr[ row+1 ];
        for( j=start; j<end; j++)
            dot += dval[ j ] * dx[ dcolind[j] ];
        dy[ row ] =  dot *alpha + beta * dy[ row ];
    }
}

// shifted CSR-SpMV kernel
__global__ void 
dgecsrmv_kernel_shift( 
    int num_rows, 
    int num_cols, 
    double alpha, 
    double lambda, 
    double * dval, 
    magma_minproduct_index_t * drowptr, 
    magma_minproduct_index_t * dcolind,
    double * dx,
    double beta, 
    int offset,
    int blocksize,
    magma_minproduct_index_t * addrows,
    double * dy)
{

    int row = blockIdx.x*blockDim.x+threadIdx.x;
    int j;

    if(row<num_rows){
        double dot = MAGMA_minproduct_D_ZERO;
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
    transA      magma_minproduct_trans_t
                transposition parameter for A
                
    @param[in]
    m           magma_minproduct_int_t
                number of rows in A

    @param[in]
    n           magma_minproduct_int_t
                number of columns in A 

    @param[in]
    alpha       double
                scalar multiplier

    @param[in]
    dval        magma_minproductDouble_ptr
                array containing values of A in CSR

    @param[in]
    drowptr     magma_minproductIndex_ptr
                rowpointer of A in CSR

    @param[in]
    dcolind     magma_minproductIndex_ptr
                columnindices of A in CSR

    @param[in]
    dx          magma_minproductDouble_ptr
                input vector x

    @param[in]
    beta        double
                scalar multiplier

    @param[out]
    dy          magma_minproductDouble_ptr
                input/output vector y

    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_dblas
    ********************************************************************/

extern "C" magma_minproduct_int_t
magma_minproduct_dgecsrmv(
    magma_minproduct_trans_t transA,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    double alpha,
    magma_minproductDouble_ptr dval,
    magma_minproductIndex_ptr drowptr,
    magma_minproductIndex_ptr dcolind,
    magma_minproductDouble_ptr dx,
    double beta,
    magma_minproductDouble_ptr dy,
    magma_minproduct_queue_t queue )
{
    dim3 grid( magma_minproduct_ceildiv( m, BLOCK_SIZE ) );
    magma_minproduct_int_t threads = BLOCK_SIZE;
    dgecsrmv_kernel<<< grid, threads, 0, queue >>>
                    (m, n, alpha, dval, drowptr, dcolind, dx, beta, dy);

    return MAGMA_minproduct_SUCCESS;
}



/**
    Purpose
    -------
    
    This routine computes y = alpha * ( A -lambda I ) * x + beta * y on the GPU.
    It is a shifted version of the CSR-SpMV.
    
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
    alpha       double
                scalar multiplier

    @param[in]
    lambda      double
                scalar multiplier

    @param[in]
    dval        magma_minproductDouble_ptr
                array containing values of A in CSR

    @param[in]
    drowptr     magma_minproductIndex_ptr
                rowpointer of A in CSR

    @param[in]
    dcolind     magma_minproductIndex_ptr
                columnindices of A in CSR

    @param[in]
    dx          magma_minproductDouble_ptr
                input vector x

    @param[in]
    beta        double
                scalar multiplier

    @param[in]
    offset      magma_minproduct_int_t 
                in case not the main diagonal is scaled
                
    @param[in]
    blocksize   magma_minproduct_int_t 
                in case of processing multiple vectors  
                
    @param[in]
    addrows     magma_minproductIndex_ptr
                in case the matrixpowerskernel is used
                
    @param[out]
    dy          magma_minproductDouble_ptr
                output vector y  
    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_dblas
    ********************************************************************/

extern "C" magma_minproduct_int_t
magma_minproduct_dgecsrmv_shift(
    magma_minproduct_trans_t transA,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    double alpha,
    double lambda,
    magma_minproductDouble_ptr dval,
    magma_minproductIndex_ptr drowptr,
    magma_minproductIndex_ptr dcolind,
    magma_minproductDouble_ptr dx,
    double beta,
    int offset,
    int blocksize,
    magma_minproduct_index_t * addrows,
    magma_minproductDouble_ptr dy,
    magma_minproduct_queue_t queue )
{
    dim3 grid( magma_minproduct_ceildiv( m, BLOCK_SIZE ) );
    magma_minproduct_int_t threads = BLOCK_SIZE;
    dgecsrmv_kernel_shift<<< grid, threads, 0, queue >>>
                         (m, n, alpha, lambda, dval, drowptr, dcolind, dx, 
                                    beta, offset, blocksize, addrows, dy);

    return MAGMA_minproduct_SUCCESS;
}



