/*
    -- MAGMA_minproduct (version 1.6.1) --
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
magma_minproductblas_celementary_multiplication_kernel(
    magma_minproduct_int_t n,
    magma_minproductFloatComplex *dA, magma_minproduct_int_t offsetA, magma_minproduct_int_t ldda, 
    magma_minproductFloatComplex *du, magma_minproduct_int_t offsetu, 
    magma_minproductFloatComplex *dv, magma_minproduct_int_t offsetv);

__global__ void 
magma_minproductblas_capply_vector_kernel(
    magma_minproduct_int_t n,
    magma_minproductFloatComplex *du, magma_minproduct_int_t offsetu,  magma_minproductFloatComplex *db, magma_minproduct_int_t offsetb );

__global__ void 
magma_minproductblas_capply_transpose_vector_kernel(
    magma_minproduct_int_t n,
    magma_minproductFloatComplex *du, magma_minproduct_int_t offsetu, magma_minproductFloatComplex *db, magma_minproduct_int_t offsetb );
/////////////////////////////////////
// batched prototypes
/////////////////////////////////////
__global__ void 
magma_minproductblas_celementary_multiplication_kernel_batched(
    magma_minproduct_int_t n,
    magma_minproductFloatComplex **dA_array, magma_minproduct_int_t offsetA, magma_minproduct_int_t ldda, 
    magma_minproductFloatComplex *du, magma_minproduct_int_t offsetu, 
    magma_minproductFloatComplex *dv, magma_minproduct_int_t offsetv);

__global__ void 
magma_minproductblas_capply_vector_kernel_batched(
    magma_minproduct_int_t n,
    magma_minproductFloatComplex *du, magma_minproduct_int_t offsetu, magma_minproductFloatComplex **db_array, magma_minproduct_int_t offsetb );

__global__ void 
magma_minproductblas_capply_transpose_vector_kernel_batched(
    magma_minproduct_int_t n,
    magma_minproductFloatComplex *du, magma_minproduct_int_t offsetu, magma_minproductFloatComplex **db_array, magma_minproduct_int_t offsetb );

#endif        //  #ifndef CGERBT_H
