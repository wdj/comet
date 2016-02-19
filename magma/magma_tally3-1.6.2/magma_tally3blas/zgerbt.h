/*
    -- MAGMA_tally3 (version 1.6.1) --
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
magma_tally3blas_zelementary_multiplication_kernel(
    magma_tally3_int_t n,
    magma_tally3DoubleComplex *dA, magma_tally3_int_t offsetA, magma_tally3_int_t ldda, 
    magma_tally3DoubleComplex *du, magma_tally3_int_t offsetu, 
    magma_tally3DoubleComplex *dv, magma_tally3_int_t offsetv);

__global__ void 
magma_tally3blas_zapply_vector_kernel(
    magma_tally3_int_t n,
    magma_tally3DoubleComplex *du, magma_tally3_int_t offsetu,  magma_tally3DoubleComplex *db, magma_tally3_int_t offsetb );

__global__ void 
magma_tally3blas_zapply_transpose_vector_kernel(
    magma_tally3_int_t n,
    magma_tally3DoubleComplex *du, magma_tally3_int_t offsetu, magma_tally3DoubleComplex *db, magma_tally3_int_t offsetb );
/////////////////////////////////////
// batched prototypes
/////////////////////////////////////
__global__ void 
magma_tally3blas_zelementary_multiplication_kernel_batched(
    magma_tally3_int_t n,
    magma_tally3DoubleComplex **dA_array, magma_tally3_int_t offsetA, magma_tally3_int_t ldda, 
    magma_tally3DoubleComplex *du, magma_tally3_int_t offsetu, 
    magma_tally3DoubleComplex *dv, magma_tally3_int_t offsetv);

__global__ void 
magma_tally3blas_zapply_vector_kernel_batched(
    magma_tally3_int_t n,
    magma_tally3DoubleComplex *du, magma_tally3_int_t offsetu, magma_tally3DoubleComplex **db_array, magma_tally3_int_t offsetb );

__global__ void 
magma_tally3blas_zapply_transpose_vector_kernel_batched(
    magma_tally3_int_t n,
    magma_tally3DoubleComplex *du, magma_tally3_int_t offsetu, magma_tally3DoubleComplex **db_array, magma_tally3_int_t offsetb );

#endif        //  #ifndef ZGERBT_H
