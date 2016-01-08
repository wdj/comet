/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions normal z -> c d s

       @author Adrien Remy
       @author Azzam Haidar
       
       Definitions used in zgerbt.cu zgerbt_batched.cu
*/

#ifndef ZGERBT_H
#define ZGERBT_H
/////////////////////////////////////
// classical prototypes
/////////////////////////////////////

__global__ void 
magma_tally4blas_zelementary_multiplication_kernel(
    magma_tally4_int_t n,
    magma_tally4DoubleComplex *dA, magma_tally4_int_t offsetA, magma_tally4_int_t ldda, 
    magma_tally4DoubleComplex *du, magma_tally4_int_t offsetu, 
    magma_tally4DoubleComplex *dv, magma_tally4_int_t offsetv);

__global__ void 
magma_tally4blas_zapply_vector_kernel(
    magma_tally4_int_t n,
    magma_tally4DoubleComplex *du, magma_tally4_int_t offsetu,  magma_tally4DoubleComplex *db, magma_tally4_int_t offsetb );

__global__ void 
magma_tally4blas_zapply_transpose_vector_kernel(
    magma_tally4_int_t n,
    magma_tally4DoubleComplex *du, magma_tally4_int_t offsetu, magma_tally4DoubleComplex *db, magma_tally4_int_t offsetb );
/////////////////////////////////////
// batched prototypes
/////////////////////////////////////
__global__ void 
magma_tally4blas_zelementary_multiplication_kernel_batched(
    magma_tally4_int_t n,
    magma_tally4DoubleComplex **dA_array, magma_tally4_int_t offsetA, magma_tally4_int_t ldda, 
    magma_tally4DoubleComplex *du, magma_tally4_int_t offsetu, 
    magma_tally4DoubleComplex *dv, magma_tally4_int_t offsetv);

__global__ void 
magma_tally4blas_zapply_vector_kernel_batched(
    magma_tally4_int_t n,
    magma_tally4DoubleComplex *du, magma_tally4_int_t offsetu, magma_tally4DoubleComplex **db_array, magma_tally4_int_t offsetb );

__global__ void 
magma_tally4blas_zapply_transpose_vector_kernel_batched(
    magma_tally4_int_t n,
    magma_tally4DoubleComplex *du, magma_tally4_int_t offsetu, magma_tally4DoubleComplex **db_array, magma_tally4_int_t offsetb );

#endif        //  #ifndef ZGERBT_H
