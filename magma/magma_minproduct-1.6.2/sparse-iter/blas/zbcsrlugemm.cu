/*
    -- MAGMA_minproduct (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @precisions normal z -> c d s

*/
#include "common_magma_minproductsparse.h"

#define BLOCK_SIZE 512

#define PRECISION_z

#define  Ablockinfo(i,j)  Ablockinfo[(i)*c_blocks   + (j)]
#define  Bblockinfo(i,j)  Bblockinfo[(i)*c_blocks   + (j)]
#define A(i,j) ((Ablockinfo(i,j)-1)*size_b*size_b)
#define B(i,j) ((Bblockinfo(i,j)-1)*size_b*size_b)

//============================================================

#define ldb m
#define lda m
#define ldc m


#define fetch_x_A(i) (((i)<m*m)?Aval[i]:0)
#define fetch_x_B(i) (((i)<m*m)?B[i]:0)


// every multiprocessor handles one BCSR-block
__global__ void 
zbcsr_gemm_kernel32( 
    int m,
    int n,
    int kblocks,   
    mdouble **Avals, 
    double **Bval,
    double **Cval)
{
#if (__CUDA_ARCH__ >= 200)

#if defined(PRECISION_d)
    const  int tx = threadIdx.x;
    const  int ty = threadIdx.y;
  
    const int idt = ty * 64 + tx;

    const int tx2 = idt%16;
    const int ty2 = idt/16;

    double xxB[4];
    magma_minproductDouble_ptr B;

    int trackA = __mul24( ty2, lda) + tx2 ;
    magma_minproductDouble_ptr Aval = Avals[blockIdx.z];

    __shared__ double Abs[64][65];
    __shared__ double  Bb[16][65];


    for(int j=ty2; j<64; j+=16){
        for(int y=tx2; y<64; y+=16){
           Abs[y][j] = fetch_x_A(trackA + y-tx2) ;
            }
        trackA += __mul24( 16, m);
    }

    for(int k=0; k<kblocks; k++){
        B = Bval[k];
        int trackB = tx2+ __mul24( ty2 * 16, ldb );

        // Prefetch part of B
          #pragma unroll
          for(int y=0; y<4; y++){
                 Bb[tx2][ty2*4+y] = fetch_x_B( trackB + y * ldb) ;
          }
        __syncthreads();    // this is necessary!!!

        double Axs[4];
        double Bxp[4];
        double Cb[16] = {0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0};

        int k1;
        for(k1=0; k1<m-16; k1+=16)
        {
                trackB += 16;

                #pragma unroll
                for( int y=0; y<4; y++)
                        xxB[y] = fetch_x_B( trackB + y*ldb);
                #pragma unroll
                for( int j1=0;j1<16;j1++)
                {
                        #pragma unroll
                        for( int y=0; y<4; y++){
                                Axs[y] =  Abs[tx2+y*16][j1+k1] ;
                        }

                        #pragma unroll
                        for( int y=0; y<4; y++){
                                Bxp[y]= Bb[j1][ty2+y*16];
                        }

                        #pragma unroll
                        for( int x=0; x<4; x++)
                        {
                                #pragma unroll
                                for( int y=0; y<4; y++)
                                {
                                        Cb[x*4+y]  += Axs[x]*Bxp[y];
                                }
                        }

                }
                #pragma unroll
                for(int y=0; y<4; y++)
                        Bb[tx2][ty2*4 + y] = xxB[y];

                __syncthreads();     // this is necessary!!!
        }
        // Prepare where to write the result
        magma_minproductDouble_ptr C = Cval[blockIdx.z * kblocks + k];
        C += tx2 + __mul24 (ty2 ,ldc);

        #pragma unroll
        for(int j1=0;j1<16;j1++)
        {

                #pragma unroll
                for( int y=0; y<4; y++)
                        Axs[y] =  Abs[tx2 + y*16][j1+k1] ;

                #pragma unroll
                for( int y=0; y<4; y++)
                        Bxp[y]= Bb[j1][ty2 + y*16];

                #pragma unroll
                for( int x=0; x<4; x++)
                {
                        #pragma unroll
                        for( int y=0;y<4; y++)
                        {
                                Cb[x*4 + y]  += Axs[x]*Bxp[y];
                        }
                }
        }   
        int gy = ty2;
        #pragma unroll
        for( int y=0;y<4;y++, gy+=16)
        {
                int gx = tx2;
        #pragma unroll
                for(int x=0;x<4;x++, gx+=16)
                {
                        if (gx < m && gy < n){
                              C[x*16] -= Cb[y+x*4];
                       }
                }
                C += ldc*16;
        }
      }
#endif

#endif
}

// every multiprocessor handles one BCSR-block
__global__ void 
zbcsr_gemm_kernel64( 
    int m,
    int n,
    int kblocks,   
    double **Avals, 
    double **Bval,
    double **Cval)
{
#if (__CUDA_ARCH__ >= 200)

#if defined(PRECISION_d)
    const  int tx = threadIdx.x;
    const  int ty = threadIdx.y;
  
    const int idt = ty * 64 + tx;

    const int tx2 = idt%16;
    const int ty2 = idt/16;

    double xxB[4];

    magma_minproductDouble_ptr B;

    int trackA = __mul24( ty2, lda) + tx2 ;
    magma_minproductDouble_ptr Aval = Avals[blockIdx.z];

    __shared__ double Abs[64][65];
    __shared__ double  Bb[16][65];


    for(int j=ty2; j<64; j+=16){
        for(int y=tx2; y<64; y+=16){
           Abs[y][j] = fetch_x_A(trackA + y-tx2) ;
            }
        trackA += __mul24( 16, m);
    }


    for(int k=0; k<kblocks; k++){

        B = Bval[k];
        int trackB = tx2+ __mul24( ty2 * 4, ldb );

        // Prefetch part of B
          #pragma unroll
          for(int y=0; y<4; y++){
                 Bb[tx2][ty2*4+y] = fetch_x_B( trackB + y * ldb) ;
          }

        __syncthreads();    // this is necessary!!!

        double Axs[4];
        double Bxp[4];

        double Cb[16] = {0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0};

        int k1;
        for(k1=0; k1<m-16; k1+=16)
        {
                trackB += 16;

                #pragma unroll
                for( int y=0; y<4; y++)
                        xxB[y] = fetch_x_B( trackB + y*ldb);

                #pragma unroll
                for( int j1=0;j1<16;j1++)
                {
                        #pragma unroll
                        for( int y=0; y<4; y++){
                                Axs[y] =  Abs[tx2+y*16][j1+k1] ;
                        }

                        #pragma unroll
                        for( int y=0; y<4; y++){
                                Bxp[y]= Bb[j1][ty2+y*16];
                        }

                        #pragma unroll
                        for( int x=0; x<4; x++)
                        {
                                #pragma unroll
                                for( int y=0; y<4; y++)
                                {
                                        Cb[x*4+y]  += Axs[x]*Bxp[y];
                                }
                        }

                }

                __syncthreads();
                #pragma unroll
                for(int y=0; y<4; y++)
                        Bb[tx2][ty2*4 + y] = xxB[y];

                __syncthreads();     // this is necessary!!!

        }
        // Prepare where to write the result
        magma_minproductDouble_ptr C = Cval[blockIdx.z * kblocks + k];
        C += tx2 + __mul24 (ty2 ,ldc);

        #pragma unroll
        for(int j1=0;j1<16;j1++)
        {

                #pragma unroll
                for( int y=0; y<4; y++)
                        Axs[y] =  Abs[tx2 + y*16][j1+k1] ;

                #pragma unroll
                for( int y=0; y<4; y++)
                        Bxp[y]= Bb[j1][ty2 + y*16];

                #pragma unroll
                for( int x=0; x<4; x++)
                {
                        #pragma unroll
                        for( int y=0;y<4; y++)
                        {
                                Cb[x*4 + y]  += Axs[x]*Bxp[y];
                        }
                }
        }   

        int gy = ty2;
        #pragma unroll
        for( int y=0;y<4;y++, gy+=16)
        {
                int gx = tx2;
        #pragma unroll
                for(int x=0;x<4;x++, gx+=16)
                {
                        if (gx < m && gy < n){
                              C[x*16] -= Cb[y+x*4];
                       }
                }

                C += ldc*16;
        }

      }
#endif

#endif
}





/**
    Purpose
    -------
    
    For a Block-CSR ILU factorization, this routine updates all blocks in
    the trailing matrix.
    
    Arguments
    ---------

    @param[in]
    size_b      magma_minproduct_int_t
                blocksize in BCSR

    @param[in]
    num_brows   magma_minproduct_int_t
                number of block rows

    @param[in]
    kblocks     magma_minproduct_int_t
                number of blocks in row

    @param[in]
    dA          magma_minproductDoubleComplex**
                input blocks of matrix A
                
    @param[in]
    dB          magma_minproductDoubleComplex**
                input blocks of matrix B
                
    @param[in]
    dC          magma_minproductDoubleComplex**
                output blocks of matrix C
    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_zgegpuk
    ********************************************************************/

extern "C" magma_minproduct_int_t
magma_minproduct_zbcsrluegemm(
    magma_minproduct_int_t size_b, 
    magma_minproduct_int_t num_brows,
    magma_minproduct_int_t kblocks,
    magma_minproductDoubleComplex_ptr *dA,  
    magma_minproductDoubleComplex_ptr *dB,  
    magma_minproductDoubleComplex_ptr *dC,
    magma_minproduct_queue_t queue )
{
    #if defined(PRECISION_d)

    magma_minproduct_int_t arch = magma_minproduct_getdevice_arch();

    if ( arch < 200  ) {
        printf("error: magma_minproduct_zbcsrluegemm needs a CUDA architecture"
               " with at least 48K shared memory (Fermi +).\n"
               "Please run zbcsrlu.cpp using CUBLAS batched.\n");
    
    }
    else {

    dim3 threads( 64, 4 );

    dim3 grid(1, 1, num_brows);
    zbcsr_gemm_kernel64<<< grid, threads, 0, queue >>>( 
                  size_b, size_b, kblocks, dA, dB, dC );

    }

#else
    printf("error: currently only supported for double precision.\n"
           "Please run zbcsrlu.cpp using CUBLAS batched.\n");
#endif

    return MAGMA_minproduct_SUCCESS;
}



