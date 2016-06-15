/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from zgerbt.h normal z -> s, Fri Jan 30 19:00:10 2015

       @author Adrien Remy
       @author Azzam Haidar
       
       Definitions used in sgerbt.cu sgerbt_batched.cu
*/

#ifndef SGERBT_H
#define SGERBT_H
/////////////////////////////////////
// classical prototypes
/////////////////////////////////////

__global__ void 
magma_tally3blas_selementary_multiplication_kernel(
    magma_tally3_int_t n,
    float *dA, magma_tally3_int_t offsetA, magma_tally3_int_t ldda, 
    float *du, magma_tally3_int_t offsetu, 
    float *dv, magma_tally3_int_t offsetv);

__global__ void 
magma_tally3blas_sapply_vector_kernel(
    magma_tally3_int_t n,
    float *du, magma_tally3_int_t offsetu,  float *db, magma_tally3_int_t offsetb );

__global__ void 
magma_tally3blas_sapply_transpose_vector_kernel(
    magma_tally3_int_t n,
    float *du, magma_tally3_int_t offsetu, float *db, magma_tally3_int_t offsetb );
/////////////////////////////////////
// batched prototypes
/////////////////////////////////////
__global__ void 
magma_tally3blas_selementary_multiplication_kernel_batched(
    magma_tally3_int_t n,
    float **dA_array, magma_tally3_int_t offsetA, magma_tally3_int_t ldda, 
    float *du, magma_tally3_int_t offsetu, 
    float *dv, magma_tally3_int_t offsetv);

__global__ void 
magma_tally3blas_sapply_vector_kernel_batched(
    magma_tally3_int_t n,
    float *du, magma_tally3_int_t offsetu, float **db_array, magma_tally3_int_t offsetb );

__global__ void 
magma_tally3blas_sapply_transpose_vector_kernel_batched(
    magma_tally3_int_t n,
    float *du, magma_tally3_int_t offsetu, float **db_array, magma_tally3_int_t offsetb );

#endif        //  #ifndef SGERBT_H