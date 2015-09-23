/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2013
       
       @author Azzam Haidar
       @author Tingxing Dong

       @generated from zpotf2_kernels.cu normal z -> c, Fri Jan 30 19:00:10 2015
*/

#include "common_magma_minproduct.h"
#include "batched_kernel_param.h"
#include "magma_minproduct_templates.h"

#define PRECISION_c


#define A(i, j)  (A + (i) + (j)*lda)   // A(i, j) means at i row, j column


// dynamically allocated shared memory, set to size number of threads when the kernel is launched.
// See CUDA Guide B.2.3
extern __shared__ magma_minproductFloatComplex shared_data[];


// dynamically allocated shared memory, set to size number of threads when the kernel is launched.
// See CUDA Guide B.2.3
extern __shared__ float dble_shared_data[];

/////////////////////////////////////////////////////////////////////////////////////////////////
__global__ void cdotc_kernel_batched(int n, magma_minproductFloatComplex **x_array, int incx, int offset, magma_minproduct_int_t *info_array, int gbstep)
{
    int tx = threadIdx.x;

    magma_minproductFloatComplex *x = x_array[blockIdx.z]+offset;

    float *sdata = dble_shared_data;

    magma_minproductFloatComplex res = MAGMA_minproduct_C_ZERO;

    if (tx < n) {
       res = x[tx*incx];
    }

    sdata[tx] = MAGMA_minproduct_C_REAL(res * MAGMA_minproduct_C_CNJG(res));

    __syncthreads();

    for(int s = blockDim.x/2; s > 32; s >>= 1 ) {
        if (tx < s) {
            sdata[tx] += sdata[tx+s];
        }
        __syncthreads();
    }

    if (tx < 32) {
        volatile float* smem = sdata;
        smem[tx] += smem[tx+32];
        smem[tx] += smem[tx+16];
        smem[tx] += smem[tx+8];
        smem[tx] += smem[tx+4];
        smem[tx] += smem[tx+2];
        smem[tx] += smem[tx+1];
    }

    if (tx == 0) {
        float xreal = MAGMA_minproduct_C_REAL(x[n*incx]);        
        //MAGMA_minproduct_C_SET2REAL(x[n*incx], sqrt(xreal - sdata[0]));
        x[n*incx] = MAGMA_minproduct_C_MAKE(sqrt(xreal - sdata[0]), 0);
        if(x[n*incx] == MAGMA_minproduct_C_ZERO){
            info_array[blockIdx.z] = offset + gbstep + 1;
        }
    }
}


void magma_minproduct_cpotf2_cdotc_batched(magma_minproduct_int_t n, magma_minproductFloatComplex **x_array, magma_minproduct_int_t incx, magma_minproduct_int_t offset, magma_minproduct_int_t *info_array, magma_minproduct_int_t gbstep, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue)
{
/*
    Specialized Cdotc
    1) performs cdotc sum = x[0:n-1]*conj(x[0:n-1])
    2) updates x[n] = sqrt(x[n]-sum);

*/
    if (n > MAX_NTHREADS) {
        printf("n = %d > %d is not supported in cpotf2_cdotc\n", (int) n, (int) MAX_NTHREADS);
        
    }
    int threadSize;

    if (n <= 1024 && n > 512) {
        threadSize = 1024;
    }
    else if (n <= 512 && n > 256 ) {
        threadSize = 512;
    }
    else if (n <= 256 && n > 128) {
        threadSize = 256;
    }
    else if (n <= 128 && n > 64) {
        threadSize = 128;
    }
    else {
        threadSize = 64;
    }

    
    dim3 grid(1, 1, batchCount);
    cdotc_kernel_batched<<< grid, threadSize, 
                  threadSize * sizeof(float), queue>>> (n, x_array, incx, offset, info_array, gbstep);
}

/////////////////////////////////////////////////////////////////////////////////////////////////
__global__ void csscal_kernel_batched(int n, magma_minproductFloatComplex **x_array, int incx, int offset, magma_minproduct_int_t *info_array)
{
    // checkinfo to avoid computation of the singular matrix
    if(info_array[blockIdx.z] != 0 ) return;

    int id = threadIdx.x;
    magma_minproductFloatComplex *x = x_array[blockIdx.z]+offset;

    __shared__ magma_minproductFloatComplex factor;

    if (threadIdx.x == 0) {
        factor = MAGMA_minproduct_C_MAKE(1.0/MAGMA_minproduct_C_REAL(x[0]), 0.0);
    }

    __syncthreads();

    if ( id < n && id >0) {
        x[id*incx] = x[id*incx] * factor;
        //printf("x=%f", x[id*incx]);
    }
}


void magma_minproduct_cpotf2_csscal_batched(magma_minproduct_int_t n, magma_minproductFloatComplex **x_array, magma_minproduct_int_t incx, magma_minproduct_int_t offset, magma_minproduct_int_t *info_array, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue)
{
/*
    Specialized Csscal perform x[1:n-1]/x[0]

*/
    dim3 grid(1, 1, batchCount);
    dim3 threads(n, 1, 1); 

    csscal_kernel_batched<<< grid, threads, 0, queue >>> (n, x_array, incx, offset, info_array);
}

/////////////////////////////////////////////////////////////////////////////////////////////////


#if defined(PRECISION_z) || defined(PRECISION_c)

__global__ void clacgv_kernel_batched(int n, magma_minproductFloatComplex **x_array, int incx, int offset)
{
    int id = threadIdx.x;

    magma_minproductFloatComplex *x = x_array[blockIdx.z]+offset;

    if ( id < n ) {
        x[id*incx] = MAGMA_minproduct_C_CNJG(x[id*incx]);
    }
}

void magma_minproduct_clacgv_batched(magma_minproduct_int_t n, magma_minproductFloatComplex **x_array, magma_minproduct_int_t incx, int offset, int batchCount, magma_minproduct_queue_t queue)
{
/*
    Purpose
    =======

    CLACGV conjugates a complex vector of length N.

    Arguments
    =========

    N       (input) INTEGER
            The length of the vector X.  N >= 0.

    X       (input/output) COMPLEX array, dimension
                           (1+(N-1)*abs(INCX))
            On entry, the vector of length N to be conjugated.
            On exit, X is overwritten with conjg(X).

    INCX    (input) INTEGER
            The spacing between successive elements of X.

    ===================================================================== */

    dim3 grid(1, 1, batchCount);
    dim3 threads(n, 1, 1);
   
    clacgv_kernel_batched<<< grid, threads, 0, queue >>> (n, x_array, incx, offset);
}

#endif // defined(PRECISION_z) || defined(PRECISION_c)



/////////////////////////////////////////////////////////////////////////////////////////////////
static __device__ void cpotf2_device(int m, int n, 
                              magma_minproductFloatComplex *A, int lda, 
                              magma_minproductFloatComplex alpha, 
                              magma_minproductFloatComplex beta, magma_minproduct_int_t *info, int gbstep)
{
/*
    Each thread block load entire A into shared memory
    factorize it and copy back. n must be small enough to fit shared memory.
    n is checked by a macro POTF2_TILE_SIZE before the kernel. 
*/
    // checkinfo to avoid computation of the singular matrix
    if(*info != 0 ) return;

    int tx = threadIdx.x;
    magma_minproductFloatComplex *sdata_A = shared_data;
    __shared__ magma_minproductFloatComplex factor;
    __shared__ float sum[POTF2_TILE_SIZE];

    // load A into sdata_A
    if(tx < m)
    {
        for(int i=0; i<n; i++)
        {  
             sdata_A[tx + i * m] =  A[tx + i * lda];
        }
    }
    __syncthreads();

    for(int iter=0; iter<n; iter++)
    {
        float res = MAGMA_minproduct_D_ZERO;
        magma_minproductFloatComplex res1 = MAGMA_minproduct_C_ZERO;

        //1) performs cdotc sum = A[iter, 0:iter-1]*conj(A[iter, 0:iter-1])
        //2) updates A[iter,iter] = sqrt(A[iter,iter]-sum);
        if(tx<iter)
        {
            res = MAGMA_minproduct_C_REAL (sdata_A[iter + tx * m] * MAGMA_minproduct_C_CNJG(sdata_A[iter + tx * m]));         
            sum[tx] = res;
        }
        else
        {
            sum[tx] = 0.0;
        }
        __syncthreads();
        magma_minproduct_sum_reduce<POTF2_TILE_SIZE>(tx, sum);//tried on K40: if m=32 n=32 the overall cpotf2_device routine time is 60ms n=16 time=25 n=8 time=20ms 
        //magma_minproduct_sum_reduce_n(iter, tx, sum); //tried on K40: if m=32 n=32 the time went from 61ms to 70ms when switching to reduce_n. n=16 time=28.
        //magma_minproduct_sum_reduce_inlined(iter, tx, sum); //tried on K40: similar to magma_minproduct_sum_reduce<POTF2_TILE_SIZE>(tx, sum);
        
        if (tx == 0) {
              float xreal = MAGMA_minproduct_C_REAL(sdata_A[iter + iter * m]);        
              sdata_A[iter + iter * m] = MAGMA_minproduct_C_MAKE(sqrt(xreal - sum[0]), 0);
              if(sdata_A[iter + iter * m] == MAGMA_minproduct_C_ZERO){
                  *info = iter + gbstep + 1;
              }
        }
        __syncthreads();
        if(sdata_A[iter + iter * m] == MAGMA_minproduct_C_ZERO) return;
        __syncthreads();

        //clacgv conjugates a complex vector of length iter. //TODO
        #if defined(PRECISION_z) || defined(PRECISION_c)
        if(tx < iter)
        {
             sdata_A[iter + tx * m] = MAGMA_minproduct_C_CNJG(sdata_A[iter + tx * m]);
        }
        __syncthreads();  
        #endif
  
        // cgemv  
        // Compute elements iter:n-1 of column iter = A(iter:n,0:iter-1) * A(iter-1,0:iter-1) (row).
        if(tx < m && tx > iter)
        {
            for(int j=0; j < iter; j++)
            {
                res1 += sdata_A[tx + j * m]  *  sdata_A[iter + j * m]; // TODO move the clacgv conj to be done automatically here implicitly.
            }   
            sdata_A [tx + iter * m] = alpha * res1 + sdata_A [tx + iter * m] * beta;   
        }
        __syncthreads();  

        //clacgv conjugates a complex vector of length iter.
        #if defined(PRECISION_z) || defined(PRECISION_c)
        if(tx < iter)
        {
             sdata_A[iter + tx * m] = MAGMA_minproduct_C_CNJG(sdata_A[iter + tx * m]);
        }
        __syncthreads();  
        #endif

        // csscal perform A[iter:n-1, iter]/A[iter,iter];
        if (tx == 0) {
            factor = MAGMA_minproduct_C_MAKE(1.0/MAGMA_minproduct_C_REAL(sdata_A[iter + iter * m]), 0.0);
        }
        __syncthreads();

        if ( tx < m && tx > iter) {
            sdata_A[ tx + iter * m ]  *= factor;
        }
        __syncthreads();
    }// end of iter

    //copy sdata_A to A
    if(tx < m)
    {
        for(int i=0; i<n; i++)
        {  
             A[tx + i * lda] = sdata_A[tx + i * m];
        }
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////
__global__ void cpotf2_kernel_batched(int m, int n, 
                              magma_minproductFloatComplex **dA_array, int lda, 
                              magma_minproductFloatComplex alpha, 
                              magma_minproductFloatComplex beta, 
                              magma_minproduct_int_t *info_array, int gbstep)
{
/*
    Each thread block load entire dA_array[blockIdx.z] into shared memory
    factorize it and copy back. n must be small enough to fit shared memory.
    n is checked by a macro POTF2_TILE_SIZE before the kernel. 
*/
    int batchid = blockIdx.z;
    cpotf2_device(m, n, dA_array[batchid], lda, alpha, beta, &(info_array[batchid]), gbstep);
}
/////////////////////////////////////////////////////////////////////////////////////////////////
__global__ void cpotf2_kernel(int m, int n, 
                              magma_minproductFloatComplex *dA, int lda, 
                              magma_minproductFloatComplex alpha, 
                              magma_minproductFloatComplex beta,
                              magma_minproduct_int_t *info)
{
    cpotf2_device(m, n, dA, lda, alpha, beta, info, 0);
}
/////////////////////////////////////////////////////////////////////////////////////////////////
/**
    Purpose
    -------

    cpotf2 computes the Cholesky factorization of a real symmetric
    positive definite matrix A.

    The factorization has the form
        A = U**H * U,  if UPLO = Magma_minproductUpper, or
        A = L  * L**H, if UPLO = Magma_minproductLower,
    where U is an upper triangular matrix and L is lower triangular.

    This is the unblocked version of the algorithm, calling Level 2 BLAS.

    Arguments
    ---------

    @param[in]
    uplo    magma_minproduct_uplo_t
            Specifies whether the upper or lower triangular part of the
            symmetric matrix A is stored.
      -     = Magma_minproductUpper:  Upper triangular
      -     = Magma_minproductLower:  Lower triangular

    @param[in]
    n       INTEGER
            The order of the matrix A.  N >= 0 and N <= 512.

    @param[in,out]
    dA      COMPLEX array, dimension (LDDA,N)
            On entry, the symmetric matrix A.  If UPLO = Magma_minproductUpper, the leading
            n by n upper triangular part of A contains the upper
            triangular part of the matrix A, and the strictly lower
            triangular part of A is not referenced.  If UPLO = Magma_minproductLower, the
            leading n by n lower triangular part of A contains the lower
            triangular part of the matrix A, and the strictly upper
            triangular part of A is not referenced.
    \n
            On exit, if INFO = 0, the factor U or L from the Cholesky
            factorization A = U**H * U  or A = L * L**H.

    @param[in]
    ldda    INTEGER
            The leading dimension of the array A.  LDDA >= max(1,N).

    @param[out]
    info    INTEGER
      -     = 0: successful exit
      -     < 0: if INFO = -k, the k-th argument had an illegal value
      -     > 0: if INFO = k, the leading minor of order k is not
                 positive definite, and the factorization could not be
                 completed.

    @ingroup magma_minproduct_cposv_aux
    ********************************************************************/
extern "C" magma_minproduct_int_t
magma_minproduct_cpotf2_tile_batched(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex **dA_array, magma_minproduct_int_t lda,
    magma_minproduct_int_t *info_array, magma_minproduct_int_t gbstep, magma_minproduct_int_t batchCount, magma_minproduct_queue_t queue)
{

    magma_minproduct_int_t arginfo = 0;
    
    if ( uplo != Magma_minproductUpper && uplo != Magma_minproductLower) {
        arginfo = -1;
    } else if (m < 0 || n < 0 || m > POTF2_TILE_SIZE || n > POTF2_TILE_SIZE) {
        arginfo = -2;
    } else if (lda < max(1,m)) {
        arginfo = -4;
    } else if (m < n) {
        arginfo = -10;
    }
    if (uplo == Magma_minproductUpper) {
        printf("Upper side is unavailable \n");
        arginfo = -1;
    }

    if (arginfo != 0) {
        magma_minproduct_xerbla( __func__, -(arginfo) );
        return arginfo;
    }
    
    // Quick return if possible
    if (m == 0 || n == 0) {
        return arginfo;
    }

    magma_minproductFloatComplex alpha = MAGMA_minproduct_C_NEG_ONE;
    magma_minproductFloatComplex beta  = MAGMA_minproduct_C_ONE;

    dim3 dimGrid(1, 1, batchCount);
    dim3 threads(POTF2_TILE_SIZE, 1);
    int shared_mem_size = sizeof(magma_minproductFloatComplex)*m*n; // + sizeof(float)*(POTF2_TILE_SIZE+1);

    cpotf2_kernel_batched<<<dimGrid, threads, shared_mem_size, queue >>>(m, n, dA_array, lda, alpha, beta, info_array, gbstep);

    return arginfo;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////
extern "C" magma_minproduct_int_t
magma_minproduct_cpotf2_tile(
    magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductFloatComplex *dA, magma_minproduct_int_t lda,
    magma_minproduct_int_t *info)
{

    *info = 0;
    if ( uplo != Magma_minproductUpper && uplo != Magma_minproductLower) {
        *info = -1;
    } else if (m < 0 || n < 0 || m > POTF2_TILE_SIZE) {
        *info = -2;
    } else if (lda < max(1,m)) {
        *info = -4;
    } else if (m < n) {
        *info = -10;
    }
    if (uplo == Magma_minproductUpper) {
        printf("Upper side is unavailable \n");
        *info = -1;
    }


    if (*info != 0) {
        magma_minproduct_xerbla( __func__, -(*info) );
        return *info;
    }

    // Quick return if possible
    if (m == 0 || n == 0) {
        return *info;
    }

    magma_minproductFloatComplex alpha = MAGMA_minproduct_C_NEG_ONE;
    magma_minproductFloatComplex beta  = MAGMA_minproduct_C_ONE;

    dim3 dimGrid(1);
    dim3 threads(POTF2_TILE_SIZE, 1);
    int shared_mem_size = sizeof(magma_minproductFloatComplex)*m*n; // + sizeof(float)*(POTF2_TILE_SIZE+1);

    cpotf2_kernel<<<dimGrid, threads, shared_mem_size, magma_minproduct_stream >>>(m, n, dA, lda, alpha, beta, info);

    return *info;
}

