/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from zgerbt.h normal z -> c, Fri Jan 30 19:00:10 2015

       @author Adrien Remy
       @author Azzam Haidar
       
       Definitions used in cgerbt.cu cgerbt_batched.cu
*/

#ifndef CGERBT_H
#define CGERBT_H
/////////////////////////////////////
// classical prototypes
/////////////////////////////////////

__global__ void 
magma_tally4blas_celementary_multiplication_kernel(
    magma_tally4_int_t n,
    magma_tally4FloatComplex *dA, magma_tally4_int_t offsetA, magma_tally4_int_t ldda, 
    magma_tally4FloatComplex *du, magma_tally4_int_t offsetu, 
    magma_tally4FloatComplex *dv, magma_tally4_int_t offsetv);

__global__ void 
magma_tally4blas_capply_vector_kernel(
    magma_tally4_int_t n,
    magma_tally4FloatComplex *du, magma_tally4_int_t offsetu,  magma_tally4FloatComplex *db, magma_tally4_int_t offsetb );

__global__ void 
magma_tally4blas_capply_transpose_vector_kernel(
    magma_tally4_int_t n,
    magma_tally4FloatComplex *du, magma_tally4_int_t offsetu, magma_tally4FloatComplex *db, magma_tally4_int_t offsetb );
/////////////////////////////////////
// batched prototypes
/////////////////////////////////////
__global__ void 
magma_tally4blas_celementary_multiplication_kernel_batched(
    magma_tally4_int_t n,
    magma_tally4FloatComplex **dA_array, magma_tally4_int_t offsetA, magma_tally4_int_t ldda, 
    magma_tally4FloatComplex *du, magma_tally4_int_t offsetu, 
    magma_tally4FloatComplex *dv, magma_tally4_int_t offsetv);

__global__ void 
magma_tally4blas_capply_vector_kernel_batched(
    magma_tally4_int_t n,
    magma_tally4FloatComplex *du, magma_tally4_int_t offsetu, magma_tally4FloatComplex **db_array, magma_tally4_int_t offsetb );

__global__ void 
magma_tally4blas_capply_transpose_vector_kernel_batched(
    magma_tally4_int_t n,
    magma_tally4FloatComplex *du, magma_tally4_int_t offsetu, magma_tally4FloatComplex **db_array, magma_tally4_int_t offsetb );

#endif        //  #ifndef CGERBT_H
