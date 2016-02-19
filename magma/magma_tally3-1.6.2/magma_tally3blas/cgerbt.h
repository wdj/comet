/*
    -- MAGMA_tally3 (version 1.6.1) --
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
magma_tally3blas_celementary_multiplication_kernel(
    magma_tally3_int_t n,
    magma_tally3FloatComplex *dA, magma_tally3_int_t offsetA, magma_tally3_int_t ldda, 
    magma_tally3FloatComplex *du, magma_tally3_int_t offsetu, 
    magma_tally3FloatComplex *dv, magma_tally3_int_t offsetv);

__global__ void 
magma_tally3blas_capply_vector_kernel(
    magma_tally3_int_t n,
    magma_tally3FloatComplex *du, magma_tally3_int_t offsetu,  magma_tally3FloatComplex *db, magma_tally3_int_t offsetb );

__global__ void 
magma_tally3blas_capply_transpose_vector_kernel(
    magma_tally3_int_t n,
    magma_tally3FloatComplex *du, magma_tally3_int_t offsetu, magma_tally3FloatComplex *db, magma_tally3_int_t offsetb );
/////////////////////////////////////
// batched prototypes
/////////////////////////////////////
__global__ void 
magma_tally3blas_celementary_multiplication_kernel_batched(
    magma_tally3_int_t n,
    magma_tally3FloatComplex **dA_array, magma_tally3_int_t offsetA, magma_tally3_int_t ldda, 
    magma_tally3FloatComplex *du, magma_tally3_int_t offsetu, 
    magma_tally3FloatComplex *dv, magma_tally3_int_t offsetv);

__global__ void 
magma_tally3blas_capply_vector_kernel_batched(
    magma_tally3_int_t n,
    magma_tally3FloatComplex *du, magma_tally3_int_t offsetu, magma_tally3FloatComplex **db_array, magma_tally3_int_t offsetb );

__global__ void 
magma_tally3blas_capply_transpose_vector_kernel_batched(
    magma_tally3_int_t n,
    magma_tally3FloatComplex *du, magma_tally3_int_t offsetu, magma_tally3FloatComplex **db_array, magma_tally3_int_t offsetb );

#endif        //  #ifndef CGERBT_H
