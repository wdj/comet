/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @generated from zlaswp_batched.cu normal z -> d, Fri Jan 30 19:00:10 2015
       
       @author Azzam Haidar
       @author Tingxing Dong
*/
#include "common_magma_tally2.h"
#include "batched_kernel_param.h"

#define BLK_SIZE 256
// SWP_WIDTH is number of threads in a block
// 64 and 256 are better on Kepler; 
extern __shared__ double shared_data[];


//=================================================================================================
static __device__ 
void dlaswp_rowparallel_devfunc(  
                              int n, int width, int height,
                              double *dA, int lda, 
                              double *dout, int ldo,
                              magma_tally2_int_t* pivinfo)
{

    //int height = k2- k1;
    //int height = blockDim.x;
    unsigned int tid = threadIdx.x;
    dA   += SWP_WIDTH * blockIdx.x * lda;
    dout += SWP_WIDTH * blockIdx.x * ldo;
    double *sdata = shared_data;

    if(blockIdx.x == gridDim.x -1)
    {
       width = n - blockIdx.x * SWP_WIDTH;
    }

    if(tid < height)
    {
        int mynewroworig = pivinfo[tid]-1; //-1 to get the index in C
        int itsreplacement = pivinfo[mynewroworig] -1 ; //-1 to get the index in C
        #pragma unroll
        for(int i=0; i<width; i++)
        {
          sdata[ tid + i * height ]    = dA[ mynewroworig + i * lda ];
          dA[ mynewroworig + i * lda ] = dA[ itsreplacement + i * lda ];
        }
    }
    __syncthreads();

    if(tid < height)
    {
        // copy back the upper swapped portion of A to dout 
        #pragma unroll
        for(int i=0; i<width; i++)
        {
           dout[tid + i * ldo] = sdata[tid + i * height];
        }
    }
}

//=================================================================================================
// parallel swap the swaped dA(1:nb,i:n) is stored in dout 
//=================================================================================================
__global__ 
void dlaswp_rowparallel_kernel( 
                                int n, int width, int height,
                                double *dinput, int ldi, 
                                double *doutput, int ldo,
                                magma_tally2_int_t*  pivinfo)
{

    dlaswp_rowparallel_devfunc(n, width, height, dinput, ldi, doutput, ldo, pivinfo);

}
//=================================================================================================

__global__ 
void dlaswp_rowparallel_kernel_batched(
                                int n, int width, int height,
                                double **input_array, int ldi, 
                                double **output_array, int ldo,
                                magma_tally2_int_t** pivinfo_array)
{
    int batchid = blockIdx.z;
    dlaswp_rowparallel_devfunc(n, width, height, input_array[batchid], ldi, output_array[batchid], ldo, pivinfo_array[batchid]);
}


//=================================================================================================

//=================================================================================================
extern "C" void
magma_tally2_dlaswp_rowparallel_batched( magma_tally2_int_t n, 
                       double** input_array, magma_tally2_int_t ldi,
                       double** output_array, magma_tally2_int_t ldo,
                       magma_tally2_int_t k1, magma_tally2_int_t k2,
                       magma_tally2_int_t **pivinfo_array, 
                       magma_tally2_int_t batchCount, magma_tally2_queue_t queue)
{

    if(n == 0 ) return ;
    int height = k2-k1;
    if( height  > 1024) 
    {
       printf(" n=%d > 1024, not supported \n", n);

    }

    int blocks =  (n-1)/ SWP_WIDTH + 1;
    dim3  grid(blocks, 1, batchCount);

    if( n < SWP_WIDTH)
    {
        dlaswp_rowparallel_kernel_batched<<<grid, height, sizeof(double) * height * n, queue >>>
                                           ( n, n, height, input_array, ldi, output_array, ldo, pivinfo_array ); 
    }
    else
    {
        dlaswp_rowparallel_kernel_batched<<< grid, height, sizeof(double) * height * SWP_WIDTH , queue >>>
                                            (n, SWP_WIDTH, height, input_array, ldi, output_array, ldo, pivinfo_array ); 
 
    }
}

//=================================================================================================




//=================================================================================================
extern "C" void
magma_tally2_dlaswp_rowparallel_q( magma_tally2_int_t n, 
                       double* input, magma_tally2_int_t ldi,
                       double* output, magma_tally2_int_t ldo,
                       magma_tally2_int_t k1, magma_tally2_int_t k2,
                       magma_tally2_int_t *pivinfo, 
                       magma_tally2_queue_t queue)
{
    if(n == 0 ) return ;
    int height = k2-k1;
    if( height  > MAX_NTHREADS) 
    {
       printf(" height=%d > %d, magma_tally2_dlaswp_rowparallel_q not supported \n", n,MAX_NTHREADS);

    }

    int blocks =  (n-1)/ SWP_WIDTH + 1;
    dim3  grid(blocks, 1, 1);

    if( n < SWP_WIDTH)
    {
        dlaswp_rowparallel_kernel<<<grid, height, sizeof(double) * height * n, queue >>>
                                   ( n, n, height, input, ldi, output, ldo, pivinfo ); 
    }
    else
    {
        dlaswp_rowparallel_kernel<<< grid, height, sizeof(double) * height * SWP_WIDTH , queue >>>
                                    (n, SWP_WIDTH, height, input, ldi, output, ldo, pivinfo ); 
    }
}


//=================================================================================================

extern "C" void
magma_tally2_dlaswp_rowparallel( magma_tally2_int_t n, double* input, magma_tally2_int_t ldi,
                   double* output, magma_tally2_int_t ldo,
                   magma_tally2_int_t k1, magma_tally2_int_t k2,
                   magma_tally2_int_t *pivinfo)
{
    magma_tally2_dlaswp_rowparallel_q(n, input, ldi, output, ldo, k1, k2, pivinfo, magma_tally2_stream);
}

//=================================================================================================





//=================================================================================================
//  serial swap that does swapping one row by one row
//=================================================================================================
__global__ void dlaswp_rowserial_kernel_batched( int n, double **dA_array, int lda, int k1, int k2, magma_tally2_int_t** ipiv_array )
{
    double* dA = dA_array[blockIdx.z];
    magma_tally2_int_t *d_ipiv = ipiv_array[blockIdx.z];
    
    unsigned int tid = threadIdx.x + blockDim.x*blockIdx.x;
    
    k1--;
    k2--;

    if( tid < n) {

        double A1;

        for( int i1 = k1; i1 < k2; i1++ ) 
        {
            int i2 = d_ipiv[i1] - 1;  // Fortran index, switch i1 and i2
            if( i2 != i1)
            {
                A1 = dA[i1 + tid * lda];
                dA[i1 + tid * lda] = dA[i2 + tid * lda];
                dA[i2 + tid * lda] = A1;
            }
        }
    }
}

//=================================================================================================
//  serial swap that does swapping one row by one row, similar to LAPACK
//  K1, K2 are in Fortran indexing  
//=================================================================================================
extern "C" void
magma_tally2_dlaswp_rowserial_batched(magma_tally2_int_t n, double** dA_array, magma_tally2_int_t lda,
                   magma_tally2_int_t k1, magma_tally2_int_t k2,
                   magma_tally2_int_t **ipiv_array, 
                   magma_tally2_int_t batchCount, magma_tally2_queue_t queue)
{

    if(n == 0 ) return ;

    int blocks =  (n-1)/ BLK_SIZE + 1;
    dim3  grid(blocks, 1, batchCount);

    dlaswp_rowserial_kernel_batched<<< grid, max(BLK_SIZE, n), 0, queue >>>(
        n, dA_array, lda, k1, k2, ipiv_array); 

}



//=================================================================================================
//  serial swap that does swapping one column by one column
//=================================================================================================
__global__ void dlaswp_columnserial_kernel_batched( int n, double **dA_array, int lda, int k1, int k2, magma_tally2_int_t** ipiv_array )
{
    double* dA = dA_array[blockIdx.z];
    magma_tally2_int_t *d_ipiv = ipiv_array[blockIdx.z];
    
    unsigned int tid = threadIdx.x + blockDim.x*blockIdx.x;
    k1--;
    k2--;
    if( k1 < 0 || k2 < 0 ) return;


    if( tid < n) {
        double A1;
        if(k1 <= k2)
        {
            for( int i1 = k1; i1 <= k2; i1++ ) 
            {
                int i2 = d_ipiv[i1] - 1;  // Fortran index, switch i1 and i2
                if( i2 != i1)
                {
                    A1 = dA[i1 * lda + tid];
                    dA[i1 * lda + tid] = dA[i2 * lda + tid];
                    dA[i2 * lda + tid] = A1;
                }
            }
        }else
        {
            for( int i1 = k1; i1 >= k2; i1-- ) 
            {
                int i2 = d_ipiv[i1] - 1;  // Fortran index, switch i1 and i2
                if( i2 != i1)
                {
                    A1 = dA[i1 * lda + tid];
                    dA[i1 * lda + tid] = dA[i2 * lda + tid];
                    dA[i2 * lda + tid] = A1;
                }
            }
        }
    }
}

//=================================================================================================
//  serial swap that does swapping one column by one column
//  K1, K2 are in Fortran indexing  
//=================================================================================================
extern "C" void
magma_tally2_dlaswp_columnserial_batched(magma_tally2_int_t n, double** dA_array, magma_tally2_int_t lda,
                   magma_tally2_int_t k1, magma_tally2_int_t k2,
                   magma_tally2_int_t **ipiv_array, 
                   magma_tally2_int_t batchCount, magma_tally2_queue_t queue)
{

    if(n == 0 ) return ;

    int blocks =  (n-1)/ BLK_SIZE + 1;
    dim3  grid(blocks, 1, batchCount);

    dlaswp_columnserial_kernel_batched<<< grid, min(BLK_SIZE, n), 0, queue >>>(
        n, dA_array, lda, k1, k2, ipiv_array); 

}

