/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from zgerbt_kernels.cu normal z -> c, Fri Jan 30 19:00:08 2015


       @author Adrien REMY
*/
#include "common_magma_tally2.h"


#define block_height  32
#define block_width  4
#define block_length 256
#define NB 64
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
static __device__ void 
magma_tally2blas_celementary_multiplication_devfunc(
    magma_tally2_int_t n,
    magma_tally2FloatComplex *dA, magma_tally2_int_t ldda, 
    magma_tally2FloatComplex *du, 
    magma_tally2FloatComplex *dv)
{    
    magma_tally2_int_t idx, idy;

    idx = blockIdx.x * blockDim.x + threadIdx.x;
    idy = blockIdx.y * blockDim.y + threadIdx.y;

    if ((idx < n/2)&&(idy < n/2)){

        dA += idx + idy * ldda;

        magma_tally2FloatComplex a00, a10, a01, a11, b1, b2, b3, b4;
        __shared__ magma_tally2FloatComplex u1[block_height], u2[block_height], v1[block_width], v2[block_width];

        du += idx;
        dv += idy;

        u1[threadIdx.x]=du[0];
        u2[threadIdx.x]=du[n/2];
        v1[threadIdx.y]=dv[0];
        v2[threadIdx.y]=dv[n/2];

        __syncthreads();

        a00 = dA[0];
        a01 = dA[ldda*n/2];
        a10 = dA[n/2];
        a11 = dA[ldda*n/2+n/2];

        b1 = a00 + a01;
        b2 = a10 + a11;
        b3 = a00 - a01;
        b4 = a10 - a11;

        dA[0] = u1[threadIdx.x] * v1[threadIdx.y] * (b1 + b2);
        dA[ldda*n/2] = u1[threadIdx.x] * v2[threadIdx.y] * (b3 + b4);
        dA[n/2] = u2[threadIdx.x] * v1[threadIdx.y] * (b1 - b2);
        dA[ldda*n/2+n/2] = u2[threadIdx.x] * v2[threadIdx.y] *(b3 - b4);
    }
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
__global__ void 
magma_tally2blas_celementary_multiplication_kernel(
    magma_tally2_int_t n,
    magma_tally2FloatComplex *dA, magma_tally2_int_t offsetA, magma_tally2_int_t ldda, 
    magma_tally2FloatComplex *du, magma_tally2_int_t offsetu, 
    magma_tally2FloatComplex *dv, magma_tally2_int_t offsetv)
{    
    magma_tally2blas_celementary_multiplication_devfunc( n, dA+offsetA, ldda, du+offsetu, dv+offsetv);
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
__global__ void 
magma_tally2blas_celementary_multiplication_kernel_batched(
    magma_tally2_int_t n,
    magma_tally2FloatComplex **dA_array, magma_tally2_int_t offsetA, magma_tally2_int_t ldda, 
    magma_tally2FloatComplex *du, magma_tally2_int_t offsetu, 
    magma_tally2FloatComplex *dv, magma_tally2_int_t offsetv)
{    
    int batchid = blockIdx.z;
    magma_tally2blas_celementary_multiplication_devfunc( n, dA_array[batchid]+offsetA, ldda, du+offsetu, dv+offsetv);
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
static __device__ void 
magma_tally2blas_capply_vector_devfunc(
    magma_tally2_int_t n,
    magma_tally2FloatComplex *du, magma_tally2FloatComplex *db)
{
    magma_tally2_int_t idx;

    idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (idx < n/2){

        du += idx;
        db += idx;

        magma_tally2FloatComplex a1,a2;

        a1 = du[0]*db[0];
        a2 = du[n/2]*db[n/2];

        db[0] = a1 + a2;
        db[n/2] = a1 -a2;
    }
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
__global__ void 
magma_tally2blas_capply_vector_kernel(
    magma_tally2_int_t n,
    magma_tally2FloatComplex *du, magma_tally2_int_t offsetu,  magma_tally2FloatComplex *db, magma_tally2_int_t offsetb )
{
    magma_tally2blas_capply_vector_devfunc(n, du+offsetu, db+offsetb);
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
__global__ void 
magma_tally2blas_capply_vector_kernel_batched(
    magma_tally2_int_t n,
    magma_tally2FloatComplex *du, magma_tally2_int_t offsetu, magma_tally2FloatComplex **db_array, magma_tally2_int_t offsetb )
{

    int batchid = blockIdx.y;
    magma_tally2blas_capply_vector_devfunc(n, du+offsetu, db_array[batchid]+offsetb);
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

static __device__ void 
magma_tally2blas_capply_transpose_vector_devfunc(
    magma_tally2_int_t n,
    magma_tally2FloatComplex *du,magma_tally2FloatComplex *db )
{
    magma_tally2_int_t idx;

    idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (idx < n/2){

        du += idx;
        db += idx;

        magma_tally2FloatComplex a1,a2;

        a1 = db[0] + db[n/2];
        a2 = db[0] - db[n/2];

        db[0] = du[0]*a1;
        db[n/2] = du[n/2]*a2;
    }
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
__global__ void 
magma_tally2blas_capply_transpose_vector_kernel(
    magma_tally2_int_t n,
    magma_tally2FloatComplex *du, magma_tally2_int_t offsetu, magma_tally2FloatComplex *db, magma_tally2_int_t offsetb )
{
    magma_tally2blas_capply_transpose_vector_devfunc(n, du+offsetu, db+offsetb);
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
__global__ void 
magma_tally2blas_capply_transpose_vector_kernel_batched(
    magma_tally2_int_t n,
    magma_tally2FloatComplex *du, magma_tally2_int_t offsetu, magma_tally2FloatComplex **db_array, magma_tally2_int_t offsetb )
{

    int batchid = blockIdx.y;
    magma_tally2blas_capply_transpose_vector_devfunc(n, du+offsetu, db_array[batchid]+offsetb);
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////


