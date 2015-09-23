/*
    -- MAGMA_minproduct (version 1.6.1) --
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
magma_minproductblas_zelementary_multiplication_kernel(
    magma_minproduct_int_t n,
    magma_minproductDoubleComplex *dA, magma_minproduct_int_t offsetA, magma_minproduct_int_t ldda, 
    magma_minproductDoubleComplex *du, magma_minproduct_int_t offsetu, 
    magma_minproductDoubleComplex *dv, magma_minproduct_int_t offsetv);

__global__ void 
magma_minproductblas_zapply_vector_kernel(
    magma_minproduct_int_t n,
    magma_minproductDoubleComplex *du, magma_minproduct_int_t offsetu,  magma_minproductDoubleComplex *db, magma_minproduct_int_t offsetb );

__global__ void 
magma_minproductblas_zapply_transpose_vector_kernel(
    magma_minproduct_int_t n,
    magma_minproductDoubleComplex *du, magma_minproduct_int_t offsetu, magma_minproductDoubleComplex *db, magma_minproduct_int_t offsetb );
/////////////////////////////////////
// batched prototypes
/////////////////////////////////////
__global__ void 
magma_minproductblas_zelementary_multiplication_kernel_batched(
    magma_minproduct_int_t n,
    magma_minproductDoubleComplex **dA_array, magma_minproduct_int_t offsetA, magma_minproduct_int_t ldda, 
    magma_minproductDoubleComplex *du, magma_minproduct_int_t offsetu, 
    magma_minproductDoubleComplex *dv, magma_minproduct_int_t offsetv);

__global__ void 
magma_minproductblas_zapply_vector_kernel_batched(
    magma_minproduct_int_t n,
    magma_minproductDoubleComplex *du, magma_minproduct_int_t offsetu, magma_minproductDoubleComplex **db_array, magma_minproduct_int_t offsetb );

__global__ void 
magma_minproductblas_zapply_transpose_vector_kernel_batched(
    magma_minproduct_int_t n,
    magma_minproductDoubleComplex *du, magma_minproduct_int_t offsetu, magma_minproductDoubleComplex **db_array, magma_minproduct_int_t offsetb );

#endif        //  #ifndef ZGERBT_H
