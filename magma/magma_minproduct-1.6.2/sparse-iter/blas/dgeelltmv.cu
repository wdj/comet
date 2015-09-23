/*
    -- MAGMA_minproduct (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from zgeelltmv.cu normal z -> d, Sun May  3 11:22:58 2015

*/

#include "common_magma_minproduct.h"

#define BLOCK_SIZE 512


// ELL SpMV kernel
//Michael Garland
__global__ void 
dgeelltmv_kernel( 
    int num_rows, 
    int num_cols,
    int num_cols_per_row,
    double alpha, 
    double * dval, 
    magma_minproduct_index_t * dcolind,
    double * dx,
    double beta, 
    double * dy)
{
    int row = blockDim.x * blockIdx.x + threadIdx.x ;
    if(row < num_rows ){
        double dot = MAGMA_minproduct_D_MAKE(0.0, 0.0);
        for ( int n = 0; n < num_cols_per_row ; n ++){
            int col = dcolind [ num_rows * n + row ];
            double val = dval [ num_rows * n + row ];
            if( val != 0)
                dot += val * dx[col ];
        }
        dy[ row ] = dot * alpha + beta * dy [ row ];
    }
}

// shifted ELL SpMV kernel
//Michael Garland
__global__ void 
dgeelltmv_kernel_shift( 
    int num_rows, 
    int num_cols,
    int num_cols_per_row,
    double alpha, 
    double lambda, 
    double * dval, 
    magma_minproduct_index_t * dcolind,
    double * dx,
    double beta, 
    int offset,
    int blocksize,
    magma_minproduct_index_t * addrows,
    double * dy)
{

    int row = blockDim.x * blockIdx.x + threadIdx.x ;
    if(row < num_rows ){
        double dot = MAGMA_minproduct_D_MAKE(0.0, 0.0);
        for ( int n = 0; n < num_cols_per_row ; n ++){
            int col = dcolind [ num_rows * n + row ];
            double val = dval [ num_rows * n + row ];
            if( val != 0)
                dot += val * dx[col ];
        }
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
    
    This routine computes y = alpha *  A^t *  x + beta * y on the GPU.
    Input format is ELL.
    
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
    nnz_per_row magma_minproduct_int_t
                number of elements in the longest row 

    @param[in]
    alpha       double
                scalar multiplier

    @param[in]
    dval        magma_minproductDouble_ptr
                array containing values of A in ELL

    @param[in]
    dcolind     magma_minproductIndex_ptr
                columnindices of A in ELL

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

    @ingroup magma_minproductsparse_d
    ********************************************************************/

extern "C" magma_minproduct_int_t
magma_minproduct_dgeelltmv(
    magma_minproduct_trans_t transA,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproduct_int_t nnz_per_row,
    double alpha,
    magma_minproductDouble_ptr dval,
    magma_minproductIndex_ptr dcolind,
    magma_minproductDouble_ptr dx,
    double beta,
    magma_minproductDouble_ptr dy,
    magma_minproduct_queue_t queue )
{
    dim3 grid( magma_minproduct_ceildiv( m, BLOCK_SIZE ) );
    magma_minproduct_int_t threads = BLOCK_SIZE;
    dgeelltmv_kernel<<< grid, threads, 0, queue >>>
                  ( m, n, nnz_per_row, alpha, dval, dcolind, dx, beta, dy );


   return MAGMA_minproduct_SUCCESS;
}


/**
    Purpose
    -------
    
    This routine computes y = alpha *( A - lambda I ) * x + beta * y on the GPU.
    Input format is ELL.
    
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
    nnz_per_row magma_minproduct_int_t
                number of elements in the longest row 

    @param[in]
    alpha       double
                scalar multiplier

    @param[in]
    lambda      double
                scalar multiplier

    @param[in]
    dval        magma_minproductDouble_ptr
                array containing values of A in ELL

    @param[in]
    dcolind     magma_minproductIndex_ptr
                columnindices of A in ELL

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
                input/output vector y

    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_dblas
    ********************************************************************/

extern "C" magma_minproduct_int_t
magma_minproduct_dgeelltmv_shift(
    magma_minproduct_trans_t transA,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproduct_int_t nnz_per_row,
    double alpha,
    double lambda,
    magma_minproductDouble_ptr dval,
    magma_minproductIndex_ptr dcolind,
    magma_minproductDouble_ptr dx,
    double beta,
    int offset,
    int blocksize,
    magma_minproductIndex_ptr addrows,
    magma_minproductDouble_ptr dy,
    magma_minproduct_queue_t queue )
{
    dim3 grid( magma_minproduct_ceildiv( m, BLOCK_SIZE ) );
    magma_minproduct_int_t threads = BLOCK_SIZE;
    double tmp_shift;
    //magma_minproduct_dsetvector(1,&lambda,1,&tmp_shift,1); 
    tmp_shift = lambda;
    dgeelltmv_kernel_shift<<< grid, threads, 0, queue >>>
                  ( m, n, nnz_per_row, alpha, tmp_shift, dval, dcolind, dx, 
                            beta, offset, blocksize, addrows, dy );


   return MAGMA_minproduct_SUCCESS;
}



