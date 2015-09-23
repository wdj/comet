/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions normal z -> s d c

       @author Adrien Remy
       @author Azzam Haidar
*/
#include "common_magma_minproduct.h"
#include "zgerbt.h"


#define block_height  32
#define block_width  4
#define block_length 256
#define NB 64
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
    Purpose
    -------
    ZPRBT_MVT compute B = UTB to randomize B
    
    Arguments
    ---------
    @param[in]
    n       INTEGER
            The number of values of db.  n >= 0.

    @param[in]
    du     COMPLEX_16 array, dimension (n,2)
            The 2*n vector representing the random butterfly matrix V
    
    @param[in,out]
    db     COMPLEX_16 array, dimension (n)
            The n vector db computed by ZGESV_NOPIV_GPU
            On exit db = du*db
    
    @param[in]
    queue   magma_minproduct_queue_t
            Queue to execute in.
    ********************************************************************/
extern "C" void
magma_minproductblas_zprbt_mtv_batched(
    magma_minproduct_int_t n, 
    magma_minproductDoubleComplex *du, magma_minproductDoubleComplex **db_array,
    magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue)
{
    /*

     */
    magma_minproduct_int_t threads = block_length;
    dim3 grid ( n/(4*block_length) + ((n%(4*block_length))!=0), batchCount);

    magma_minproductblas_zapply_transpose_vector_kernel_batched<<< grid, threads, 0, queue >>>(n/2, du, n, db_array, 0);
    magma_minproductblas_zapply_transpose_vector_kernel_batched<<< grid, threads, 0, queue >>>(n/2, du, n+n/2, db_array, n/2);

    threads = block_length;
    grid = n/(2*block_length) + ((n%(2*block_length))!=0);
    magma_minproductblas_zapply_transpose_vector_kernel_batched<<< grid, threads, 0, queue >>>(n, du, 0, db_array, 0);
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
    Purpose
    -------
    ZPRBT_MV compute B = VB to obtain the non randomized solution
    
    Arguments
    ---------
    @param[in]
    n       INTEGER
            The number of values of db.  n >= 0.
    
    @param[in,out]
    db      COMPLEX_16 array, dimension (n)
            The n vector db computed by ZGESV_NOPIV_GPU
            On exit db = dv*db
    
    @param[in]
    dv      COMPLEX_16 array, dimension (n,2)
            The 2*n vector representing the random butterfly matrix V
    
    @param[in]
    queue   magma_minproduct_queue_t
            Queue to execute in.
    ********************************************************************/
extern "C" void
magma_minproductblas_zprbt_mv_batched(
    magma_minproduct_int_t n, 
    magma_minproductDoubleComplex *dv, magma_minproductDoubleComplex **db_array,
    magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue)
{

    magma_minproduct_int_t threads = block_length;
    dim3 grid ( n/(2*block_length) + ((n%(2*block_length))!=0), batchCount);

    magma_minproductblas_zapply_vector_kernel_batched<<< grid, threads, 0, queue >>>(n, dv, 0, db_array, 0);


    threads = block_length;
    grid = n/(4*block_length) + ((n%(4*block_length))!=0);

    magma_minproductblas_zapply_vector_kernel_batched<<< grid, threads, 0, queue >>>(n/2, dv, n, db_array, 0);
    magma_minproductblas_zapply_vector_kernel_batched<<< grid, threads, 0, queue >>>(n/2, dv, n+n/2, db_array, n/2);


}




/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
    Purpose
    -------
    ZPRBT randomize a square general matrix using partial randomized transformation
    
    Arguments
    ---------
    @param[in]
    n       INTEGER
            The number of columns and rows of the matrix dA.  n >= 0.
    
    @param[in,out]
    dA      COMPLEX_16 array, dimension (n,ldda)
            The n-by-n matrix dA
            On exit dA = duT*dA*d_V
    
    @param[in]
    ldda    INTEGER
            The leading dimension of the array dA.  LDA >= max(1,n).
    
    @param[in]
    du      COMPLEX_16 array, dimension (n,2)
            The 2*n vector representing the random butterfly matrix U
    
    @param[in]
    dv      COMPLEX_16 array, dimension (n,2)
            The 2*n vector representing the random butterfly matrix V
    
    @param[in]
    queue   magma_minproduct_queue_t
            Queue to execute in.

    ********************************************************************/
extern "C" void 
magma_minproductblas_zprbt_batched(
    magma_minproduct_int_t n, 
    magma_minproductDoubleComplex **dA_array, magma_minproduct_int_t ldda, 
    magma_minproductDoubleComplex *du, magma_minproductDoubleComplex *dv,
    magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue)
{
    du += ldda;
    dv += ldda;

    dim3 threads(block_height, block_width);
    dim3 grid(n/(4*block_height) + ((n%(4*block_height))!=0), 
            n/(4*block_width)  + ((n%(4*block_width))!=0),
            batchCount);

    magma_minproductblas_zelementary_multiplication_kernel_batched<<< grid, threads, 0, queue >>>(n/2, dA_array,            0, ldda, du,   0, dv,   0);
    magma_minproductblas_zelementary_multiplication_kernel_batched<<< grid, threads, 0, queue >>>(n/2, dA_array,     ldda*n/2, ldda, du,   0, dv, n/2);
    magma_minproductblas_zelementary_multiplication_kernel_batched<<< grid, threads, 0, queue >>>(n/2, dA_array,          n/2, ldda, du, n/2, dv,   0);
    magma_minproductblas_zelementary_multiplication_kernel_batched<<< grid, threads, 0, queue >>>(n/2, dA_array, ldda*n/2+n/2, ldda, du, n/2, dv, n/2);

    dim3 threads2(block_height, block_width);
    dim3 grid2(n/(2*block_height) + ((n%(2*block_height))!=0), 
            n/(2*block_width)  + ((n%(2*block_width))!=0),
            batchCount);
    magma_minproductblas_zelementary_multiplication_kernel_batched<<< grid2, threads2, 0, queue >>>(n, dA_array, 0, ldda, du, -ldda, dv, -ldda);
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////


