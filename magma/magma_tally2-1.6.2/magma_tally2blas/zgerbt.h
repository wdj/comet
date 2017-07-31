/*
    -- MAGMA_tally2 (version 1.6.1) --
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
magma_tally2blas_zelementary_multiplication_kernel(
    magma_tally2_int_t n,
    magma_tally2DoubleComplex *dA, magma_tally2_int_t offsetA, magma_tally2_int_t ldda, 
    magma_tally2DoubleComplex *du, magma_tally2_int_t offsetu, 
    magma_tally2DoubleComplex *dv, magma_tally2_int_t offsetv);

__global__ void 
magma_tally2blas_zapply_vector_kernel(
    magma_tally2_int_t n,
    magma_tally2DoubleComplex *du, magma_tally2_int_t offsetu,  magma_tally2DoubleComplex *db, magma_tally2_int_t offsetb );

__global__ void 
magma_tally2blas_zapply_transpose_vector_kernel(
    magma_tally2_int_t n,
    magma_tally2DoubleComplex *du, magma_tally2_int_t offsetu, magma_tally2DoubleComplex *db, magma_tally2_int_t offsetb );
/////////////////////////////////////
// batched prototypes
/////////////////////////////////////
__global__ void 
magma_tally2blas_zelementary_multiplication_kernel_batched(
    magma_tally2_int_t n,
    magma_tally2DoubleComplex **dA_array, magma_tally2_int_t offsetA, magma_tally2_int_t ldda, 
    magma_tally2DoubleComplex *du, magma_tally2_int_t offsetu, 
    magma_tally2DoubleComplex *dv, magma_tally2_int_t offsetv);

__global__ void 
magma_tally2blas_zapply_vector_kernel_batched(
    magma_tally2_int_t n,
    magma_tally2DoubleComplex *du, magma_tally2_int_t offsetu, magma_tally2DoubleComplex **db_array, magma_tally2_int_t offsetb );

__global__ void 
magma_tally2blas_zapply_transpose_vector_kernel_batched(
    magma_tally2_int_t n,
    magma_tally2DoubleComplex *du, magma_tally2_int_t offsetu, magma_tally2DoubleComplex **db_array, magma_tally2_int_t offsetb );

#endif        //  #ifndef ZGERBT_H
