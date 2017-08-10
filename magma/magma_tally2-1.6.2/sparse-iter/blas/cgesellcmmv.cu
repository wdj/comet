/*
    -- MAGMA_tally2 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from zgesellcmmv.cu normal z -> c, Sun May  3 11:22:58 2015

*/
#include "common_magma_tally2sparse.h"

#define PRECISION_c

//#define TEXTURE

/*
// SELLP SpMV kernel
// see paper by M. KREUTZER, G. HAGER, G WELLEIN, H. FEHSKE A. BISHOP
// A UNIFIED SPARSE MATRIX DATA FORMAT 
// FOR MODERN PROCESSORS WITH WIDE SIMD UNITS
// SELLC SpMV kernel modified assigning multiple threads to each row - 2D kernel
__global__ void 
zgesellptmv2d_kernel_4_ldg( 
    int num_rows, 
    int num_cols,
    int blocksize,
    int T,
    magma_tally2FloatComplex alpha, 
    magma_tally2FloatComplex * dval, 
    magma_tally2_index_t * dcolind,
    magma_tally2_index_t * drowptr,
    const magma_tally2FloatComplex *  __restrict__ dx,
    magma_tally2FloatComplex beta, 
    magma_tally2FloatComplex * dy)
{

#if defined(TEXTURE) && (__CUDA_ARCH__ >= 300)
   // T threads assigned to each row
    int idx = threadIdx.y ;     // thread in row
    int idy = threadIdx.x;      // local row
    int ldx = idx * blocksize + idy;
    int bdx = blockIdx.y * gridDim.x + blockIdx.x; // global block index
    int row = bdx * blocksize + idy;  // global row index

    extern __shared__ magma_tally2FloatComplex shared[];

    if(row < num_rows ){
        magma_tally2FloatComplex dot = MAGMA_tally2_C_MAKE(0.0, 0.0);
        int offset = drowptr[ bdx ];
        int block = blocksize * T; // total number of threads

        int max_ = (drowptr[ bdx+1 ]-offset)/block;  
            // number of elements each thread handles

        int kk, i1, i2;
        magma_tally2FloatComplex x1, x2, v1, v2;
        dcolind += offset + ldx ;
        dval += offset + ldx;
        for ( kk = 0; kk < max_-1 ; kk+=2 ){
            i1 = dcolind[ block*kk];
            i2 = dcolind[ block*kk + block];

            x1 = __ldg( dx+ i1  );   
            x2 = __ldg( dx+ i2  ); 

            v1 = dval[ block*kk ];
            v2 = dval[ block*kk + block];

            dot += v1 * x1;
            dot += v2 * x2;
        }
  
        if (kk<max_){
           x1 = __ldg( dx + dcolind[ block*kk]  );            
           v1 = dval[ block*kk ];

            dot += v1 * x1;
        }

        shared[ldx]  = dot;

        __syncthreads();
        if( idx < 2 ){
            shared[ldx]+=shared[ldx+blocksize*2];              
            __syncthreads();
            if( idx == 0 ) {
                dy[row] = 
                (shared[ldx]+shared[ldx+blocksize*1])*alpha + beta*dy [row];
            }

        }

    }
#endif
}
*/


// SELLP SpMV kernel
// see paper by M. KREUTZER, G. HAGER, G WELLEIN, H. FEHSKE A. BISHOP
// A UNIFIED SPARSE MATRIX DATA FORMAT 
// FOR MODERN PROCESSORS WITH WIDE SIMD UNITS
// SELLC SpMV kernel modified assigning one thread to each row - 1D kernel
__global__ void 
zgesellptmv2d_kernel_1( 
    int num_rows, 
    int num_cols,
    int blocksize,
    int T,
    magma_tally2FloatComplex alpha, 
    magma_tally2FloatComplex * dval, 
    magma_tally2_index_t * dcolind,
    magma_tally2_index_t * drowptr,
    magma_tally2FloatComplex *  dx,
    magma_tally2FloatComplex beta, 
    magma_tally2FloatComplex * dy)
{

    // threads assigned to rows
    int Idx = blockDim.x * blockIdx.x + threadIdx.x ;
    int offset = drowptr[ blockIdx.x ];
    int border = (drowptr[ blockIdx.x+1 ]-offset)/blocksize;
    if(Idx < num_rows ){
        magma_tally2FloatComplex dot = MAGMA_tally2_C_MAKE(0.0, 0.0);
        for ( int n = 0; n < border; n++){ 
            int col = dcolind [offset+ blocksize * n + threadIdx.x ];
            magma_tally2FloatComplex val = dval[offset+ blocksize * n + threadIdx.x];
            if( val != 0){
                  dot=dot+val*dx[col];
            }
        }

        dy[ Idx ] = dot * alpha + beta * dy [ Idx ];
    }
}


// SELLP SpMV kernel
// see paper by M. KREUTZER, G. HAGER, G WELLEIN, H. FEHSKE A. BISHOP
// A UNIFIED SPARSE MATRIX DATA FORMAT 
// FOR MODERN PROCESSORS WITH WIDE SIMD UNITS
// SELLC SpMV kernel modified assigning multiple threads to each row - 2D kernel
__global__ void 
zgesellptmv2d_kernel_4( 
    int num_rows, 
    int num_cols,
    int blocksize,
    int T,
    magma_tally2FloatComplex alpha, 
    magma_tally2FloatComplex * dval, 
    magma_tally2_index_t * dcolind,
    magma_tally2_index_t * drowptr,
    magma_tally2FloatComplex *  dx,
    magma_tally2FloatComplex beta, 
    magma_tally2FloatComplex * dy)
{
   // T threads assigned to each row
    int idx = threadIdx.y ;     // thread in row
    int idy = threadIdx.x;      // local row
    int ldx = idx * blocksize + idy;
    int bdx = blockIdx.y * gridDim.x + blockIdx.x; // global block index
    int row = bdx * blocksize + idy;  // global row index

    extern __shared__ magma_tally2FloatComplex shared[];

    if(row < num_rows ){
        magma_tally2FloatComplex dot = MAGMA_tally2_C_MAKE(0.0, 0.0);
        int offset = drowptr[ bdx ];
        int block = blocksize * T; // total number of threads

        int max_ = (drowptr[ bdx+1 ]-offset)/block;  
            // number of elements each thread handles

        int kk, i1, i2;
        magma_tally2FloatComplex x1, x2, v1, v2;
        dcolind += offset + ldx ;
        dval += offset + ldx;
        for ( kk = 0; kk < max_-1 ; kk+=2 ){
            i1 = dcolind[ block*kk];
            i2 = dcolind[ block*kk + block];

            x1 = dx[ i1 ];   
            x2 = dx[ i2 ]; 

            v1 = dval[ block*kk ];
            v2 = dval[ block*kk + block];

            dot += v1 * x1;
            dot += v2 * x2;
        }
  
        if (kk<max_){
           x1 = dx[ dcolind[ block*kk] ];            
           v1 = dval[ block*kk ];

            dot += v1 * x1;
        }

        shared[ldx]  = dot;

        __syncthreads();
        if( idx < 2 ){
            shared[ldx]+=shared[ldx+blocksize*2];              
            __syncthreads();
            if( idx == 0 ) {
                dy[row] = 
                (shared[ldx]+shared[ldx+blocksize*1])*alpha + beta*dy [row];
            }

        }

    }
}


// SELLP SpMV kernel
// see paper by M. KREUTZER, G. HAGER, G WELLEIN, H. FEHSKE A. BISHOP
// A UNIFIED SPARSE MATRIX DATA FORMAT 
// FOR MODERN PROCESSORS WITH WIDE SIMD UNITS
// SELLC SpMV kernel modified assigning multiple threads to each row - 2D kernel
__global__ void 
zgesellptmv2d_kernel_8( 
    int num_rows, 
    int num_cols,
    int blocksize,
    int T,
    magma_tally2FloatComplex alpha, 
    magma_tally2FloatComplex * dval, 
    magma_tally2_index_t * dcolind,
    magma_tally2_index_t * drowptr,
    magma_tally2FloatComplex *  dx,
    magma_tally2FloatComplex beta, 
    magma_tally2FloatComplex * dy)
{
   // T threads assigned to each row
    int idx = threadIdx.y ;     // thread in row
    int idy = threadIdx.x;      // local row
    int ldx = idx * blocksize + idy;
    int bdx = blockIdx.y * gridDim.x + blockIdx.x; // global block index
    int row = bdx * blocksize + idy;  // global row index

    extern __shared__ magma_tally2FloatComplex shared[];

    if(row < num_rows ){
        magma_tally2FloatComplex dot = MAGMA_tally2_C_MAKE(0.0, 0.0);
        int offset = drowptr[ bdx ];
        int block = blocksize * T; // total number of threads

        int max_ = (drowptr[ bdx+1 ]-offset)/block;  
            // number of elements each thread handles

        int kk, i1, i2;
        magma_tally2FloatComplex x1, x2, v1, v2;
        dcolind += offset + ldx ;
        dval += offset + ldx;
        for ( kk = 0; kk < max_-1 ; kk+=2 ){
            i1 = dcolind[ block*kk];
            i2 = dcolind[ block*kk + block];

            x1 = dx[ i1 ];   
            x2 = dx[ i2 ]; 

            v1 = dval[ block*kk ];
            v2 = dval[ block*kk + block];

            dot += v1 * x1;
            dot += v2 * x2;
        }
  
        if (kk<max_){
           x1 = dx[ dcolind[ block*kk] ];            
           v1 = dval[ block*kk ];

            dot += v1 * x1;
        }

        shared[ldx]  = dot;

        __syncthreads();
        if( idx < 4 ){
            shared[ldx]+=shared[ldx+blocksize*4];              
            __syncthreads();
            if( idx < 2 ) shared[ldx]+=shared[ldx+blocksize*2];   
            __syncthreads();
            if( idx == 0 ) {
                dy[row] = 
                (shared[ldx]+shared[ldx+blocksize*1])*alpha + beta*dy [row];
            }

        }

    }
}


// SELLP SpMV kernel
// see paper by M. KREUTZER, G. HAGER, G WELLEIN, H. FEHSKE A. BISHOP
// A UNIFIED SPARSE MATRIX DATA FORMAT 
// FOR MODERN PROCESSORS WITH WIDE SIMD UNITS
// SELLC SpMV kernel modified assigning multiple threads to each row - 2D kernel
__global__ void 
zgesellptmv2d_kernel_16( 
    int num_rows, 
    int num_cols,
    int blocksize,
    int T,
    magma_tally2FloatComplex alpha, 
    magma_tally2FloatComplex * dval, 
    magma_tally2_index_t * dcolind,
    magma_tally2_index_t * drowptr,
    magma_tally2FloatComplex *  dx,
    magma_tally2FloatComplex beta, 
    magma_tally2FloatComplex * dy)
{
   // T threads assigned to each row
    int idx = threadIdx.y ;     // thread in row
    int idy = threadIdx.x;      // local row
    int ldx = idx * blocksize + idy;
    int bdx = blockIdx.y * gridDim.x + blockIdx.x; // global block index
    int row = bdx * blocksize + idy;  // global row index

    extern __shared__ magma_tally2FloatComplex shared[];

    if(row < num_rows ){
        magma_tally2FloatComplex dot = MAGMA_tally2_C_MAKE(0.0, 0.0);
        int offset = drowptr[ bdx ];
        int block = blocksize * T; // total number of threads

        int max_ = (drowptr[ bdx+1 ]-offset)/block;  
            // number of elements each thread handles

        for ( int k = 0; k < max_ ; k++ ){
            magma_tally2FloatComplex val = 
                        dval[ offset + ldx + block*k ];
            int col = 
                    dcolind[ offset + ldx + block*k ];

            dot += val * dx[ col ];
        }
        shared[ldx]  = dot;

        __syncthreads();
        if( idx < 8 ){
            shared[ldx]+=shared[ldx+blocksize*8];              
            __syncthreads();
            if( idx < 4 ) shared[ldx]+=shared[ldx+blocksize*4];   
            __syncthreads();
            if( idx < 2 ) shared[ldx]+=shared[ldx+blocksize*2];   
            __syncthreads();
            if( idx == 0 ) {
                dy[row] = 
                (shared[ldx]+shared[ldx+blocksize*1])*alpha + beta*dy [row];
            }

        }

    }
}


// SELLP SpMV kernel
// see paper by M. KREUTZER, G. HAGER, G WELLEIN, H. FEHSKE A. BISHOP
// A UNIFIED SPARSE MATRIX DATA FORMAT 
// FOR MODERN PROCESSORS WITH WIDE SIMD UNITS
// SELLC SpMV kernel modified assigning multiple threads to each row - 2D kernel
__global__ void 
zgesellptmv2d_kernel_32( 
    int num_rows, 
    int num_cols,
    int blocksize,
    int T,
    magma_tally2FloatComplex alpha, 
    magma_tally2FloatComplex * dval, 
    magma_tally2_index_t * dcolind,
    magma_tally2_index_t * drowptr,
    magma_tally2FloatComplex *  dx,
    magma_tally2FloatComplex beta, 
    magma_tally2FloatComplex * dy)
{
   // T threads assigned to each row
    int idx = threadIdx.y ;     // thread in row
    int idy = threadIdx.x;      // local row
    int ldx = idx * blocksize + idy;
    int bdx = blockIdx.y * gridDim.x + blockIdx.x; // global block index
    int row = bdx * blocksize + idy;  // global row index

    extern __shared__ magma_tally2FloatComplex shared[];

    if(row < num_rows ){
        magma_tally2FloatComplex dot = MAGMA_tally2_C_MAKE(0.0, 0.0);
        int offset = drowptr[ bdx ];
        int block = blocksize * T; // total number of threads

        int max_ = (drowptr[ bdx+1 ]-offset)/block;  
            // number of elements each thread handles
        for ( int k = 0; k < max_ ; k++ ){
            magma_tally2FloatComplex val = 
                        dval[ offset + ldx + block*k ];
            int col = 
                    dcolind[ offset + ldx + block*k ];

            dot += val * dx[ col ];
        }
        shared[ldx]  = dot;

        __syncthreads();
        if( idx < 16 ){
            shared[ldx]+=shared[ldx+blocksize*16];              
            __syncthreads();
            if( idx < 8 ) shared[ldx]+=shared[ldx+blocksize*8];  
            __syncthreads();
            if( idx < 4 ) shared[ldx]+=shared[ldx+blocksize*4];   
            __syncthreads();
            if( idx < 2 ) shared[ldx]+=shared[ldx+blocksize*2];   
            __syncthreads();
            if( idx == 0 ) {
                dy[row] = 
                (shared[ldx]+shared[ldx+blocksize*1])*alpha + beta*dy [row];
            }

        }

    }
}



/************************* same but using texture mem *************************/

#if defined(PRECISION_d) && defined(TEXTURE)

__inline__ __device__ float 
read_from_tex( cudaTextureObject_t texdx, const int& i){
  int2 temp = tex1Dfetch<int2>( texdx, i ); 
  return __hiloint2float(temp.y,temp.x);
}

// SELLP SpMV kernel
// see paper by M. KREUTZER, G. HAGER, G WELLEIN, H. FEHSKE A. BISHOP
// A UNIFIED SPARSE MATRIX DATA FORMAT 
// FOR MODERN PROCESSORS WITH WIDE SIMD UNITS
// SELLC SpMV kernel modified assigning multiple threads to each row - 2D kernel
__global__ void 
zgesellptmv2d_kernel_4_tex( 
    int num_rows, 
    int num_cols,
    int blocksize,
    int T,
    magma_tally2FloatComplex alpha, 
    magma_tally2FloatComplex * dval, 
    magma_tally2_index_t * dcolind,
    magma_tally2_index_t * drowptr,
    cudaTextureObject_t texdx,
    magma_tally2FloatComplex beta, 
    magma_tally2FloatComplex * dy)
{
   // T threads assigned to each row
    int idx = threadIdx.y ;     // thread in row
    int idy = threadIdx.x;      // local row
    int ldx = idx * blocksize + idy;
    int bdx = blockIdx.y * gridDim.x + blockIdx.x; // global block index
    int row = bdx * blocksize + idy;  // global row index

    extern __shared__ magma_tally2FloatComplex shared[];

    if(row < num_rows ){
        magma_tally2FloatComplex dot = MAGMA_tally2_C_MAKE(0.0, 0.0);
        int offset = drowptr[ bdx ];
        int block = blocksize * T; // total number of threads

        int max_ = (drowptr[ bdx+1 ]-offset)/block;  
            // number of elements each thread handles

        int kk, i1, i2;
        magma_tally2FloatComplex x1, x2, v1, v2;
        dcolind += offset + ldx ;
        dval += offset + ldx;
        for ( kk = 0; kk < max_-1 ; kk+=2 ){
            i1 = dcolind[ block*kk];
            i2 = dcolind[ block*kk + block];

            x1 = read_from_tex( texdx, i1 );
            x2 = read_from_tex( texdx, i2 );

            v1 = dval[ block*kk ];
            v2 = dval[ block*kk + block];

            dot += v1 * x1;
            dot += v2 * x2;
        }
  
        if (kk<max_){
           x1 = read_from_tex( texdx, dcolind[ block*kk] );
           v1 = dval[ block*kk ];

            dot += v1 * x1;
        }

        shared[ldx]  = dot;

        __syncthreads();
        if( idx < 2 ){
            shared[ldx]+=shared[ldx+blocksize*2];              
            __syncthreads();
            if( idx == 0 ) {
                dy[row] = 
                (shared[ldx]+shared[ldx+blocksize*1])*alpha + beta*dy [row];
            }

        }

    }
}


// SELLP SpMV kernel
// see paper by M. KREUTZER, G. HAGER, G WELLEIN, H. FEHSKE A. BISHOP
// A UNIFIED SPARSE MATRIX DATA FORMAT 
// FOR MODERN PROCESSORS WITH WIDE SIMD UNITS
// SELLC SpMV kernel modified assigning multiple threads to each row - 2D kernel
__global__ void 
zgesellptmv2d_kernel_8_tex( 
    int num_rows, 
    int num_cols,
    int blocksize,
    int T,
    magma_tally2FloatComplex alpha, 
    magma_tally2FloatComplex * dval, 
    magma_tally2_index_t * dcolind,
    magma_tally2_index_t * drowptr,
    cudaTextureObject_t texdx,
    magma_tally2FloatComplex beta, 
    magma_tally2FloatComplex * dy)
{
   // T threads assigned to each row
    int idx = threadIdx.y ;     // thread in row
    int idy = threadIdx.x;      // local row
    int ldx = idx * blocksize + idy;
    int bdx = blockIdx.y * gridDim.x + blockIdx.x; // global block index
    int row = bdx * blocksize + idy;  // global row index

    extern __shared__ magma_tally2FloatComplex shared[];

    if(row < num_rows ){
        magma_tally2FloatComplex dot = MAGMA_tally2_C_MAKE(0.0, 0.0);
        int offset = drowptr[ bdx ];
        int block = blocksize * T; // total number of threads

        int max_ = (drowptr[ bdx+1 ]-offset)/block;  
            // number of elements each thread handles

        int kk, i1, i2;
        magma_tally2FloatComplex x1, x2, v1, v2;
        dcolind += offset + ldx ;
        dval += offset + ldx;
        for ( kk = 0; kk < max_-1 ; kk+=2 ){
            i1 = dcolind[ block*kk];
            i2 = dcolind[ block*kk + block];

            x1 = read_from_tex( texdx, i1 );
            x2 = read_from_tex( texdx, i2 );

            v1 = dval[ block*kk ];
            v2 = dval[ block*kk + block];

            dot += v1 * x1;
            dot += v2 * x2;
        }
  
        if (kk<max_){
           x1 = read_from_tex( texdx, dcolind[ block*kk] );
           v1 = dval[ block*kk ];

            dot += v1 * x1;
        }

        shared[ldx]  = dot;

        __syncthreads();
        if( idx < 4 ){
            shared[ldx]+=shared[ldx+blocksize*4];              
            __syncthreads();
            if( idx < 2 ) shared[ldx]+=shared[ldx+blocksize*2];   
            __syncthreads();
            if( idx == 0 ) {
                dy[row] = 
                (shared[ldx]+shared[ldx+blocksize*1])*alpha + beta*dy [row];
            }

        }

    }
}


// SELLP SpMV kernel
// see paper by M. KREUTZER, G. HAGER, G WELLEIN, H. FEHSKE A. BISHOP
// A UNIFIED SPARSE MATRIX DATA FORMAT 
// FOR MODERN PROCESSORS WITH WIDE SIMD UNITS
// SELLC SpMV kernel modified assigning multiple threads to each row - 2D kernel
__global__ void 
zgesellptmv2d_kernel_16_tex( 
    int num_rows, 
    int num_cols,
    int blocksize,
    int T,
    magma_tally2FloatComplex alpha, 
    magma_tally2FloatComplex * dval, 
    magma_tally2_index_t * dcolind,
    magma_tally2_index_t * drowptr,
    cudaTextureObject_t texdx,
    magma_tally2FloatComplex beta, 
    magma_tally2FloatComplex * dy)
{
   // T threads assigned to each row
    int idx = threadIdx.y ;     // thread in row
    int idy = threadIdx.x;      // local row
    int ldx = idx * blocksize + idy;
    int bdx = blockIdx.y * gridDim.x + blockIdx.x; // global block index
    int row = bdx * blocksize + idy;  // global row index

    extern __shared__ magma_tally2FloatComplex shared[];

    if(row < num_rows ){
        magma_tally2FloatComplex dot = MAGMA_tally2_C_MAKE(0.0, 0.0);
        int offset = drowptr[ bdx ];
        int block = blocksize * T; // total number of threads

        int max_ = (drowptr[ bdx+1 ]-offset)/block;  
            // number of elements each thread handles

        for ( int k = 0; k < max_ ; k++ ){
            magma_tally2FloatComplex val = 
                        dval[ offset + ldx + block*k ];
            int col = 
                    dcolind[ offset + ldx + block*k ];

            dot += val * read_from_tex( texdx, col );
        }
        shared[ldx]  = dot;

        __syncthreads();
        if( idx < 8 ){
            shared[ldx]+=shared[ldx+blocksize*8];              
            __syncthreads();
            if( idx < 4 ) shared[ldx]+=shared[ldx+blocksize*4];   
            __syncthreads();
            if( idx < 2 ) shared[ldx]+=shared[ldx+blocksize*2];   
            __syncthreads();
            if( idx == 0 ) {
                dy[row] = 
                (shared[ldx]+shared[ldx+blocksize*1])*alpha + beta*dy [row];
            }

        }

    }
}


// SELLP SpMV kernel
// see paper by M. KREUTZER, G. HAGER, G WELLEIN, H. FEHSKE A. BISHOP
// A UNIFIED SPARSE MATRIX DATA FORMAT 
// FOR MODERN PROCESSORS WITH WIDE SIMD UNITS
// SELLC SpMV kernel modified assigning multiple threads to each row - 2D kernel
__global__ void 
zgesellptmv2d_kernel_32_tex( 
    int num_rows, 
    int num_cols,
    int blocksize,
    int T,
    magma_tally2FloatComplex alpha, 
    magma_tally2FloatComplex * dval, 
    magma_tally2_index_t * dcolind,
    magma_tally2_index_t * drowptr,
    cudaTextureObject_t texdx,
    magma_tally2FloatComplex beta, 
    magma_tally2FloatComplex * dy)
{
   // T threads assigned to each row
    int idx = threadIdx.y ;     // thread in row
    int idy = threadIdx.x;      // local row
    int ldx = idx * blocksize + idy;
    int bdx = blockIdx.y * gridDim.x + blockIdx.x; // global block index
    int row = bdx * blocksize + idy;  // global row index

    extern __shared__ magma_tally2FloatComplex shared[];

    if(row < num_rows ){
        magma_tally2FloatComplex dot = MAGMA_tally2_C_MAKE(0.0, 0.0);
        int offset = drowptr[ bdx ];
        int block = blocksize * T; // total number of threads

        int max_ = (drowptr[ bdx+1 ]-offset)/block;  
            // number of elements each thread handles
        for ( int k = 0; k < max_ ; k++ ){
            magma_tally2FloatComplex val = 
                        dval[ offset + ldx + block*k ];
            int col = 
                    dcolind[ offset + ldx + block*k ];

            dot += val * read_from_tex( texdx, col );
        }
        shared[ldx]  = dot;

        __syncthreads();
        if( idx < 16 ){
            shared[ldx]+=shared[ldx+blocksize*16];              
            __syncthreads();
            if( idx < 8 ) shared[ldx]+=shared[ldx+blocksize*8];  
            __syncthreads();
            if( idx < 4 ) shared[ldx]+=shared[ldx+blocksize*4];   
            __syncthreads();
            if( idx < 2 ) shared[ldx]+=shared[ldx+blocksize*2];   
            __syncthreads();
            if( idx == 0 ) {
                dy[row] = 
                (shared[ldx]+shared[ldx+blocksize*1])*alpha + beta*dy [row];
            }

        }

    }
}

#endif

/*********************     end of texture versions   **************************/

/**
    Purpose
    -------
    
    This routine computes y = alpha *  A^t *  x + beta * y on the GPU.
    Input format is SELLP.
    
    Arguments
    ---------

    @param[in]
    transA      magma_tally2_trans_t
                transposition parameter for A

    @param[in]
    m           magma_tally2_int_t
                number of rows in A

    @param[in]
    n           magma_tally2_int_t
                number of columns in A 

    @param[in]
    blocksize   magma_tally2_int_t
                number of rows in one ELL-slice

    @param[in]
    slices      magma_tally2_int_t
                number of slices in matrix

    @param[in]
    alignment   magma_tally2_int_t
                number of threads assigned to one row

    @param[in]
    alpha       magma_tally2FloatComplex
                scalar multiplier

    @param[in]
    dval        magma_tally2FloatComplex_ptr
                array containing values of A in SELLP

    @param[in]
    dcolind     magma_tally2Index_ptr
                columnindices of A in SELLP

    @param[in]
    drowptr     magma_tally2Index_ptr
                rowpointer of SELLP

    @param[in]
    dx          magma_tally2FloatComplex_ptr
                input vector x

    @param[in]
    beta        magma_tally2FloatComplex
                scalar multiplier

    @param[out]
    dy          magma_tally2FloatComplex_ptr
                input/output vector y

    @param[in]
    queue       magma_tally2_queue_t
                Queue to execute in.

    @ingroup magma_tally2sparse_cblas
    ********************************************************************/

extern "C" magma_tally2_int_t
magma_tally2_cgesellpmv(
    magma_tally2_trans_t transA,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2_int_t blocksize,
    magma_tally2_int_t slices,
    magma_tally2_int_t alignment,
    magma_tally2FloatComplex alpha,
    magma_tally2FloatComplex_ptr dval,
    magma_tally2Index_ptr dcolind,
    magma_tally2Index_ptr drowptr,
    magma_tally2FloatComplex_ptr dx,
    magma_tally2FloatComplex beta,
    magma_tally2FloatComplex_ptr dy,
    magma_tally2_queue_t queue )
{
    // using a 2D thread grid

    int num_threads = blocksize*alignment;
    magma_tally2_int_t arch = magma_tally2_getdevice_arch();
    if ( arch < 200 && num_threads > 256 )
        printf("error: too much shared memory requested.\n");

    dim3 block( blocksize, alignment, 1);

    int dimgrid1 = (int) sqrt( (float)slices );
    int dimgrid2 = magma_tally2_ceildiv( slices, dimgrid1 );

    dim3 grid( dimgrid1, dimgrid2, 1);
    int Ms = num_threads * sizeof( magma_tally2FloatComplex );

    #if defined(PRECISION_d) && defined(TEXTURE)

        // Create channel.
        cudaChannelFormatDesc channel_desc;
        channel_desc = 
            cudaCreateChannelDesc(32, 32, 0, 0, cudaChannelFormatKindSigned);

        // Create resource descriptor.
        struct cudaResourceDesc resDescdx;
        memset(&resDescdx, 0, sizeof(resDescdx));
        resDescdx.resType = cudaResourceTypeLinear;
        resDescdx.res.linear.devPtr = (void*)dx;
        resDescdx.res.linear.desc = channel_desc;
        resDescdx.res.linear.sizeInBytes = m*sizeof(float);

        // Specify texture object parameters.
        struct cudaTextureDesc texDesc;
        memset(&texDesc, 0, sizeof(texDesc));
        texDesc.addressMode[0] = cudaAddressModeClamp;
        texDesc.filterMode     = cudaFilterModePoint;
        texDesc.readMode       = cudaReadModeElementType;

        // Create texture object.
        cudaTextureObject_t texdx = 0;
        cudaCreateTextureObject(&texdx, &resDescdx, &texDesc, NULL);

        cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeEightByte);

        if ( alignment == 4)
            zgesellptmv2d_kernel_4_tex<<< grid, block, Ms, queue >>>
            ( m, n, blocksize, alignment, alpha,
                dval, dcolind, drowptr, texdx, beta, dy );

        else if ( alignment == 8)
            zgesellptmv2d_kernel_8_tex<<< grid, block, Ms, queue >>>
            ( m, n, blocksize, alignment, alpha,
                dval, dcolind, drowptr, texdx, beta, dy );

        else if ( alignment == 16)
            zgesellptmv2d_kernel_16_tex<<< grid, block, Ms, queue >>>
            ( m, n, blocksize, alignment, alpha,
                dval, dcolind, drowptr, texdx, beta, dy );

        else if ( alignment == 32)
            zgesellptmv2d_kernel_32_tex<<< grid, block, Ms, queue >>>
            ( m, n, blocksize, alignment, alpha,
                dval, dcolind, drowptr, texdx, beta, dy );

        else {
            printf("error: alignment %d not supported.\n", alignment);
            return MAGMA_tally2_ERR_NOT_SUPPORTED;
        }

        cudaDestroyTextureObject(texdx);

    #else 
        if ( alignment == 1)
            zgesellptmv2d_kernel_1<<< grid, block, Ms, queue >>>
            ( m, n, blocksize, alignment, alpha,
                dval, dcolind, drowptr, dx, beta, dy );

        else if ( alignment == 4)
            zgesellptmv2d_kernel_4<<< grid, block, Ms, queue >>>
            ( m, n, blocksize, alignment, alpha,
                dval, dcolind, drowptr, dx, beta, dy );

        else if ( alignment == 8)
            zgesellptmv2d_kernel_8<<< grid, block, Ms, queue >>>
            ( m, n, blocksize, alignment, alpha,
                dval, dcolind, drowptr, dx, beta, dy );

        else if ( alignment == 16)
            zgesellptmv2d_kernel_16<<< grid, block, Ms, queue >>>
            ( m, n, blocksize, alignment, alpha,
                dval, dcolind, drowptr, dx, beta, dy );

        else if ( alignment == 32)
            zgesellptmv2d_kernel_32<<< grid, block, Ms, queue >>>
            ( m, n, blocksize, alignment, alpha,
                dval, dcolind, drowptr, dx, beta, dy );

        else {
            printf("error: alignment %d not supported.\n", alignment);
            return MAGMA_tally2_ERR_NOT_SUPPORTED;
        }
    #endif

   return MAGMA_tally2_SUCCESS;
}
