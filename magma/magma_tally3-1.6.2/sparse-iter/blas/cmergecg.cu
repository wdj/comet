/*
    -- MAGMA_tally3 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from zmergecg.cu normal z -> c, Sun May  3 11:22:58 2015
       @author Hartwig Anzt

*/
#include "common_magma_tally3sparse.h"

#define BLOCK_SIZE 512

#define PRECISION_c


// These routines merge multiple kernels from cmergecg into one
// for a description see 
// "Reformulated Conjugate Gradient for the Energy-Aware 
// Solution of Linear Systems on GPUs (ICPP '13)

// accelerated reduction for one vector
__global__ void
magma_tally3_ccgreduce_kernel_spmv1( 
    int Gs,
    int n, 
    magma_tally3FloatComplex * vtmp,
    magma_tally3FloatComplex * vtmp2 )
{

    extern __shared__ magma_tally3FloatComplex temp[];    
    int Idx = threadIdx.x;
    int blockSize = 128;
    int gridSize = blockSize  * 2 * gridDim.x; 
    temp[Idx] = MAGMA_tally3_C_MAKE( 0.0, 0.0);
    int i = blockIdx.x * ( blockSize * 2 ) + Idx;   
    while (i < Gs ) {
        temp[ Idx  ] += vtmp[ i ]; 
        temp[ Idx  ] += ( i + blockSize < Gs ) ? vtmp[ i + blockSize ] 
                                                : MAGMA_tally3_C_MAKE( 0.0, 0.0); 
        i += gridSize;
    }
    __syncthreads();
    if ( Idx < 64 ){
        temp[ Idx ] += temp[ Idx + 64 ];
    }
    __syncthreads();
    #if defined(PRECISION_z) || defined(PRECISION_c)
        if( Idx < 32 ){
            temp[ Idx ] += temp[ Idx + 32 ];__syncthreads();
            temp[ Idx ] += temp[ Idx + 16 ];__syncthreads();
            temp[ Idx ] += temp[ Idx + 8 ];__syncthreads();
            temp[ Idx ] += temp[ Idx + 4 ];__syncthreads();
            temp[ Idx ] += temp[ Idx + 2 ];__syncthreads();
            temp[ Idx ] += temp[ Idx + 1 ];__syncthreads();
        }
    #endif
    #if defined(PRECISION_d)
        if( Idx < 32 ){
            volatile float *temp2 = temp;
            temp2[ Idx ] += temp2[ Idx + 32 ];
            temp2[ Idx ] += temp2[ Idx + 16 ];
            temp2[ Idx ] += temp2[ Idx + 8 ];
            temp2[ Idx ] += temp2[ Idx + 4 ];
            temp2[ Idx ] += temp2[ Idx + 2 ];
            temp2[ Idx ] += temp2[ Idx + 1 ];
        }
    #endif
    #if defined(PRECISION_s)
        if( Idx < 32 ){
            volatile float *temp2 = temp;
            temp2[ Idx ] += temp2[ Idx + 32 ];
            temp2[ Idx ] += temp2[ Idx + 16 ];
            temp2[ Idx ] += temp2[ Idx + 8 ];
            temp2[ Idx ] += temp2[ Idx + 4 ];
            temp2[ Idx ] += temp2[ Idx + 2 ];
            temp2[ Idx ] += temp2[ Idx + 1 ];
        }
    #endif
    if ( Idx == 0 ){
        vtmp2[ blockIdx.x ] = temp[ 0 ];
    }
}

// computes the SpMV using CSR and the first step of the reduction
__global__ void
magma_tally3_ccgmerge_spmvcsr_kernel(  
    int n,
    magma_tally3FloatComplex * dval, 
    magma_tally3_index_t * drowptr, 
    magma_tally3_index_t * dcolind,
    magma_tally3FloatComplex * d,
    magma_tally3FloatComplex * z,
    magma_tally3FloatComplex * vtmp )
{

    extern __shared__ magma_tally3FloatComplex temp[]; 
    int Idx = threadIdx.x;   
    int i   = blockIdx.x * blockDim.x + Idx;
    int j;

    temp[ Idx ] = MAGMA_tally3_C_MAKE( 0.0, 0.0);

    if( i<n ){
        magma_tally3FloatComplex dot = MAGMA_tally3_C_ZERO;
        int start = drowptr[ i ];
        int end = drowptr[ i+1 ];
        for( j=start; j<end; j++)
            dot += dval[ j ] * d[ dcolind[j] ];
        z[ i ] =  dot;
        temp[ Idx ] =  d[ i ] * dot;
    }

    __syncthreads();
    if ( Idx < 128 ){
        temp[ Idx ] += temp[ Idx + 128 ];
    }
    __syncthreads();
    if ( Idx < 64 ){
        temp[ Idx ] += temp[ Idx + 64 ];
    }
    __syncthreads();
    #if defined(PRECISION_z) || defined(PRECISION_c)
        if( Idx < 32 ){
            temp[ Idx ] += temp[ Idx + 32 ];__syncthreads();
            temp[ Idx ] += temp[ Idx + 16 ];__syncthreads();
            temp[ Idx ] += temp[ Idx + 8 ];__syncthreads();
            temp[ Idx ] += temp[ Idx + 4 ];__syncthreads();
            temp[ Idx ] += temp[ Idx + 2 ];__syncthreads();
            temp[ Idx ] += temp[ Idx + 1 ];__syncthreads();
        }
    #endif
    #if defined(PRECISION_d)
        if( Idx < 32 ){
            volatile float *temp2 = temp;
            temp2[ Idx ] += temp2[ Idx + 32 ];
            temp2[ Idx ] += temp2[ Idx + 16 ];
            temp2[ Idx ] += temp2[ Idx + 8 ];
            temp2[ Idx ] += temp2[ Idx + 4 ];
            temp2[ Idx ] += temp2[ Idx + 2 ];
            temp2[ Idx ] += temp2[ Idx + 1 ];
        }
    #endif
    #if defined(PRECISION_s)
        if( Idx < 32 ){
            volatile float *temp2 = temp;
            temp2[ Idx ] += temp2[ Idx + 32 ];
            temp2[ Idx ] += temp2[ Idx + 16 ];
            temp2[ Idx ] += temp2[ Idx + 8 ];
            temp2[ Idx ] += temp2[ Idx + 4 ];
            temp2[ Idx ] += temp2[ Idx + 2 ];
            temp2[ Idx ] += temp2[ Idx + 1 ];
        }
    #endif

    if ( Idx == 0 ){
            vtmp[ blockIdx.x ] = temp[ 0 ];
    }
}

// computes the SpMV using ELL and the first step of the reduction
__global__ void
magma_tally3_ccgmerge_spmvell_kernel(  
    int n,
    int num_cols_per_row,
    magma_tally3FloatComplex * dval, 
    magma_tally3_index_t * dcolind,
    magma_tally3FloatComplex * d,
    magma_tally3FloatComplex * z,
    magma_tally3FloatComplex * vtmp )
{

    extern __shared__ magma_tally3FloatComplex temp[]; 
    int Idx = threadIdx.x;   
    int i   = blockIdx.x * blockDim.x + Idx;

    temp[ Idx ] = MAGMA_tally3_C_MAKE( 0.0, 0.0);

    if(i < n ){
        magma_tally3FloatComplex dot = MAGMA_tally3_C_MAKE(0.0, 0.0);
        for ( int k = 0; k < num_cols_per_row ; k ++){
            int col = dcolind [ n * k + i ];
            magma_tally3FloatComplex val = dval [ n * k + i ];
            if( val != 0)
                dot += val * d[ col ];
        }
        z[ i ] =  dot;
        temp[ Idx ] = d[ i ] * dot;
    }

    __syncthreads();
    if ( Idx < 128 ){
        temp[ Idx ] += temp[ Idx + 128 ];
    }
    __syncthreads();
    if ( Idx < 64 ){
        temp[ Idx ] += temp[ Idx + 64 ];
    }
    __syncthreads();
    #if defined(PRECISION_z) || defined(PRECISION_c)
        if( Idx < 32 ){
            temp[ Idx ] += temp[ Idx + 32 ];__syncthreads();
            temp[ Idx ] += temp[ Idx + 16 ];__syncthreads();
            temp[ Idx ] += temp[ Idx + 8 ];__syncthreads();
            temp[ Idx ] += temp[ Idx + 4 ];__syncthreads();
            temp[ Idx ] += temp[ Idx + 2 ];__syncthreads();
            temp[ Idx ] += temp[ Idx + 1 ];__syncthreads();
        }
    #endif
    #if defined(PRECISION_d)
        if( Idx < 32 ){
            volatile float *temp2 = temp;
            temp2[ Idx ] += temp2[ Idx + 32 ];
            temp2[ Idx ] += temp2[ Idx + 16 ];
            temp2[ Idx ] += temp2[ Idx + 8 ];
            temp2[ Idx ] += temp2[ Idx + 4 ];
            temp2[ Idx ] += temp2[ Idx + 2 ];
            temp2[ Idx ] += temp2[ Idx + 1 ];
        }
    #endif
    #if defined(PRECISION_s)
        if( Idx < 32 ){
            volatile float *temp2 = temp;
            temp2[ Idx ] += temp2[ Idx + 32 ];
            temp2[ Idx ] += temp2[ Idx + 16 ];
            temp2[ Idx ] += temp2[ Idx + 8 ];
            temp2[ Idx ] += temp2[ Idx + 4 ];
            temp2[ Idx ] += temp2[ Idx + 2 ];
            temp2[ Idx ] += temp2[ Idx + 1 ];
        }
    #endif

    if ( Idx == 0 ){
            vtmp[ blockIdx.x ] = temp[ 0 ];
    }
}

// computes the SpMV using ELLPACK and the first step of the reduction
__global__ void
magma_tally3_ccgmerge_spmvellpack_kernel(  
    int n,
    int num_cols_per_row,
    magma_tally3FloatComplex * dval, 
    magma_tally3_index_t * dcolind,
    magma_tally3FloatComplex * d,
    magma_tally3FloatComplex * z,
    magma_tally3FloatComplex * vtmp )
{

    extern __shared__ magma_tally3FloatComplex temp[]; 
    int Idx = threadIdx.x;   
    int i   = blockIdx.x * blockDim.x + Idx;

    temp[ Idx ] = MAGMA_tally3_C_MAKE( 0.0, 0.0);

    if(i < n ){
        magma_tally3FloatComplex dot = MAGMA_tally3_C_MAKE(0.0, 0.0);
        for ( int k = 0; k < num_cols_per_row ; k ++){
            int col = dcolind [ num_cols_per_row * i + k ];
            magma_tally3FloatComplex val = dval [ num_cols_per_row * i + k ];
            if( val != 0)
                dot += val * d[ col ];
        }
        z[ i ] =  dot;
        temp[ Idx ] = d[ i ] * dot;
    }

    __syncthreads();
    if ( Idx < 128 ){
        temp[ Idx ] += temp[ Idx + 128 ];
    }
    __syncthreads();
    if ( Idx < 64 ){
        temp[ Idx ] += temp[ Idx + 64 ];
    }
    __syncthreads();
    #if defined(PRECISION_z) || defined(PRECISION_c)
        if( Idx < 32 ){
            temp[ Idx ] += temp[ Idx + 32 ];__syncthreads();
            temp[ Idx ] += temp[ Idx + 16 ];__syncthreads();
            temp[ Idx ] += temp[ Idx + 8 ];__syncthreads();
            temp[ Idx ] += temp[ Idx + 4 ];__syncthreads();
            temp[ Idx ] += temp[ Idx + 2 ];__syncthreads();
            temp[ Idx ] += temp[ Idx + 1 ];__syncthreads();
        }
    #endif
    #if defined(PRECISION_d)
        if( Idx < 32 ){
            volatile float *temp2 = temp;
            temp2[ Idx ] += temp2[ Idx + 32 ];
            temp2[ Idx ] += temp2[ Idx + 16 ];
            temp2[ Idx ] += temp2[ Idx + 8 ];
            temp2[ Idx ] += temp2[ Idx + 4 ];
            temp2[ Idx ] += temp2[ Idx + 2 ];
            temp2[ Idx ] += temp2[ Idx + 1 ];
        }
    #endif
    #if defined(PRECISION_s)
        if( Idx < 32 ){
            volatile float *temp2 = temp;
            temp2[ Idx ] += temp2[ Idx + 32 ];
            temp2[ Idx ] += temp2[ Idx + 16 ];
            temp2[ Idx ] += temp2[ Idx + 8 ];
            temp2[ Idx ] += temp2[ Idx + 4 ];
            temp2[ Idx ] += temp2[ Idx + 2 ];
            temp2[ Idx ] += temp2[ Idx + 1 ];
        }
    #endif

    if ( Idx == 0 ){
            vtmp[ blockIdx.x ] = temp[ 0 ];
    }
}

// computes the SpMV using ELLRT 8 threads per row
__global__ void
magma_tally3_ccgmerge_spmvellpackrt_kernel_8(  
    int n,
    magma_tally3FloatComplex * dval, 
    magma_tally3_index_t * dcolind,
    magma_tally3_index_t * drowlength,
    magma_tally3FloatComplex * d,
    magma_tally3FloatComplex * z,
    magma_tally3FloatComplex * vtmp,
    magma_tally3_int_t T, 
    magma_tally3_int_t alignment  )
{

int idx = blockIdx.y * gridDim.x * blockDim.x + 
          blockDim.x * blockIdx.x + threadIdx.x ; // global thread index
int idb = threadIdx.x ;  // local thread index
int idp = idb%T;  // number of threads assigned to one row
int i = idx/T;  // row index

extern __shared__ magma_tally3FloatComplex shared[];

    if(i < n ){
        magma_tally3FloatComplex dot = MAGMA_tally3_C_MAKE(0.0, 0.0);
        int max_ = magma_tally3_ceildiv( drowlength[i], T );  
            // number of elements each thread handles

        for ( int k = 0; k < max_ ; k++ ){

            // original code in paper (not working for me)
            //magma_tally3FloatComplex val = dval[ k*(T*alignment)+(i*T)+idp ];  
            //int col = dcolind [ k*(T*alignment)+(i*T)+idp ];    

            // new code (working for me)        
            magma_tally3FloatComplex val = dval[ k*(T)+(i*alignment)+idp ];
            int col = dcolind [ k*(T)+(i*alignment)+idp ];

            dot += val * d[ col ];
        }
        shared[idb]  = dot;
        if( idp < 4 ){
            shared[idb]+=shared[idb+4];
            if( idp < 2 ) shared[idb]+=shared[idb+2];
            if( idp == 0 ) {
                z[i] = (shared[idb]+shared[idb+1]);
            }

        }
    }

}

// computes the SpMV using ELLRT 8 threads per row
__global__ void
magma_tally3_ccgmerge_spmvellpackrt_kernel_16(  
    int n,
    magma_tally3FloatComplex * dval, 
    magma_tally3_index_t * dcolind,
    magma_tally3_index_t * drowlength,
    magma_tally3FloatComplex * d,
    magma_tally3FloatComplex * z,
    magma_tally3FloatComplex * vtmp,
    magma_tally3_int_t T, 
    magma_tally3_int_t alignment  )
{

int idx = blockIdx.y * gridDim.x * blockDim.x + 
          blockDim.x * blockIdx.x + threadIdx.x ; // global thread index
int idb = threadIdx.x ;  // local thread index
int idp = idb%T;  // number of threads assigned to one row
int i = idx/T;  // row index

extern __shared__ magma_tally3FloatComplex shared[];

    if(i < n ){
        magma_tally3FloatComplex dot = MAGMA_tally3_C_MAKE(0.0, 0.0);
        int max_ = magma_tally3_ceildiv( drowlength[i], T );  
            // number of elements each thread handles

        for ( int k = 0; k < max_ ; k++ ){

            // original code in paper (not working for me)
            //magma_tally3FloatComplex val = dval[ k*(T*alignment)+(i*T)+idp ];  
            //int col = dcolind [ k*(T*alignment)+(i*T)+idp ];    

            // new code (working for me)        
            magma_tally3FloatComplex val = dval[ k*(T)+(i*alignment)+idp ];
            int col = dcolind [ k*(T)+(i*alignment)+idp ];

            dot += val * d[ col ];
        }
        shared[idb]  = dot;
        if( idp < 8 ){
            shared[idb]+=shared[idb+8];
            if( idp < 4 ) shared[idb]+=shared[idb+4];
            if( idp < 2 ) shared[idb]+=shared[idb+2];
            if( idp == 0 ) {
                z[i] = (shared[idb]+shared[idb+1]);
            }

        }
    }

}

// computes the SpMV using ELLRT 8 threads per row
__global__ void
magma_tally3_ccgmerge_spmvellpackrt_kernel_32(  
    int n,
    magma_tally3FloatComplex * dval, 
    magma_tally3_index_t * dcolind,
    magma_tally3_index_t * drowlength,
    magma_tally3FloatComplex * d,
    magma_tally3FloatComplex * z,
    magma_tally3FloatComplex * vtmp,
    magma_tally3_int_t T, 
    magma_tally3_int_t alignment  )
{

int idx = blockIdx.y * gridDim.x * blockDim.x + 
          blockDim.x * blockIdx.x + threadIdx.x ; // global thread index
int idb = threadIdx.x ;  // local thread index
int idp = idb%T;  // number of threads assigned to one row
int i = idx/T;  // row index

extern __shared__ magma_tally3FloatComplex shared[];

    if(i < n ){
        magma_tally3FloatComplex dot = MAGMA_tally3_C_MAKE(0.0, 0.0);
        int max_ = magma_tally3_ceildiv( drowlength[i], T );  
            // number of elements each thread handles

        for ( int k = 0; k < max_ ; k++ ){

            // original code in paper (not working for me)
            //magma_tally3FloatComplex val = dval[ k*(T*alignment)+(i*T)+idp ];  
            //int col = dcolind [ k*(T*alignment)+(i*T)+idp ];    

            // new code (working for me)        
            magma_tally3FloatComplex val = dval[ k*(T)+(i*alignment)+idp ];
            int col = dcolind [ k*(T)+(i*alignment)+idp ];

            dot += val * d[ col ];
        }
        shared[idb]  = dot;
        if( idp < 16 ){
            shared[idb]+=shared[idb+16];
            if( idp < 8 ) shared[idb]+=shared[idb+8];
            if( idp < 4 ) shared[idb]+=shared[idb+4];
            if( idp < 2 ) shared[idb]+=shared[idb+2];
            if( idp == 0 ) {
                z[i] = (shared[idb]+shared[idb+1]);
            }

        }
    }

}





// additional kernel necessary to compute first reduction step
__global__ void
magma_tally3_ccgmerge_spmvellpackrt_kernel2(  
    int n,
    magma_tally3FloatComplex * z,
    magma_tally3FloatComplex * d,
    magma_tally3FloatComplex * vtmp2 )
{

    extern __shared__ magma_tally3FloatComplex temp[]; 
    int Idx = threadIdx.x;   
    int i   = blockIdx.x * blockDim.x + Idx;
    

    temp[ Idx ] = ( i < n ) ? z[i]*d[i] : MAGMA_tally3_C_MAKE(0.0, 0.0);
    __syncthreads();
    if ( Idx < 128 ){
        temp[ Idx ] += temp[ Idx + 128 ];
    }
    __syncthreads();
    if ( Idx < 64 ){
        temp[ Idx ] += temp[ Idx + 64 ];
    }
    __syncthreads();
    #if defined(PRECISION_z) || defined(PRECISION_c)
        if( Idx < 32 ){
            temp[ Idx ] += temp[ Idx + 32 ];__syncthreads();
            temp[ Idx ] += temp[ Idx + 16 ];__syncthreads();
            temp[ Idx ] += temp[ Idx + 8 ];__syncthreads();
            temp[ Idx ] += temp[ Idx + 4 ];__syncthreads();
            temp[ Idx ] += temp[ Idx + 2 ];__syncthreads();
            temp[ Idx ] += temp[ Idx + 1 ];__syncthreads();
        }
    #endif
    #if defined(PRECISION_d)
        if( Idx < 32 ){
            volatile float *temp2 = temp;
            temp2[ Idx ] += temp2[ Idx + 32 ];
            temp2[ Idx ] += temp2[ Idx + 16 ];
            temp2[ Idx ] += temp2[ Idx + 8 ];
            temp2[ Idx ] += temp2[ Idx + 4 ];
            temp2[ Idx ] += temp2[ Idx + 2 ];
            temp2[ Idx ] += temp2[ Idx + 1 ];
        }
    #endif
    #if defined(PRECISION_s)
        if( Idx < 32 ){
            volatile float *temp2 = temp;
            temp2[ Idx ] += temp2[ Idx + 32 ];
            temp2[ Idx ] += temp2[ Idx + 16 ];
            temp2[ Idx ] += temp2[ Idx + 8 ];
            temp2[ Idx ] += temp2[ Idx + 4 ];
            temp2[ Idx ] += temp2[ Idx + 2 ];
            temp2[ Idx ] += temp2[ Idx + 1 ];
        }
    #endif

    if ( Idx == 0 ){
            vtmp2[ blockIdx.x ] = temp[ 0 ];
    }
}



// computes the SpMV using SELLC
__global__ void
magma_tally3_ccgmerge_spmvsellc_kernel(   
    int num_rows, 
    int blocksize,
    magma_tally3FloatComplex * dval, 
    magma_tally3_index_t * dcolind,
    magma_tally3_index_t * drowptr,
    magma_tally3FloatComplex * d,
    magma_tally3FloatComplex * z,
    magma_tally3FloatComplex * vtmp)
{

    extern __shared__ magma_tally3FloatComplex temp[]; 
    int Idx = threadIdx.x;   
    int i   = blockIdx.x * blockDim.x + Idx;
    int offset = drowptr[ blockIdx.x ];
    int border = (drowptr[ blockIdx.x+1 ]-offset)/blocksize;

 temp[ Idx ] = MAGMA_tally3_C_MAKE( 0.0, 0.0);


    if(i < num_rows ){
        magma_tally3FloatComplex dot = MAGMA_tally3_C_MAKE(0.0, 0.0);
        for ( int n = 0; n < border; n ++){
            int col = dcolind [offset+ blocksize * n + Idx ];
            magma_tally3FloatComplex val = dval[offset+ blocksize * n + Idx];
            if( val != 0){
                  dot=dot+val*d[col];
            }
        }
        z[ i ] =  dot;
        temp[ Idx ] = d[ i ] * dot;
    }
    __syncthreads();
    if ( Idx < 128 ){
        temp[ Idx ] += temp[ Idx + 128 ];
    }
    __syncthreads();
    if ( Idx < 64 ){
        temp[ Idx ] += temp[ Idx + 64 ];
    }
    __syncthreads();
    #if defined(PRECISION_z) || defined(PRECISION_c)
        if( Idx < 32 ){
            temp[ Idx ] += temp[ Idx + 32 ];__syncthreads();
            temp[ Idx ] += temp[ Idx + 16 ];__syncthreads();
            temp[ Idx ] += temp[ Idx + 8 ];__syncthreads();
            temp[ Idx ] += temp[ Idx + 4 ];__syncthreads();
            temp[ Idx ] += temp[ Idx + 2 ];__syncthreads();
            temp[ Idx ] += temp[ Idx + 1 ];__syncthreads();
        }
    #endif
    #if defined(PRECISION_d)
        if( Idx < 32 ){
            volatile float *temp2 = temp;
            temp2[ Idx ] += temp2[ Idx + 32 ];
            temp2[ Idx ] += temp2[ Idx + 16 ];
            temp2[ Idx ] += temp2[ Idx + 8 ];
            temp2[ Idx ] += temp2[ Idx + 4 ];
            temp2[ Idx ] += temp2[ Idx + 2 ];
            temp2[ Idx ] += temp2[ Idx + 1 ];
        }
    #endif
    #if defined(PRECISION_s)
        if( Idx < 32 ){
            volatile float *temp2 = temp;
            temp2[ Idx ] += temp2[ Idx + 32 ];
            temp2[ Idx ] += temp2[ Idx + 16 ];
            temp2[ Idx ] += temp2[ Idx + 8 ];
            temp2[ Idx ] += temp2[ Idx + 4 ];
            temp2[ Idx ] += temp2[ Idx + 2 ];
            temp2[ Idx ] += temp2[ Idx + 1 ];
        }
    #endif

    if ( Idx == 0 ){
            vtmp[ blockIdx.x ] = temp[ 0 ];
    }
}


// SELLP SpMV kernel
// see paper by M. KREUTZER, G. HAGER, G WELLEIN, H. FEHSKE A. BISHOP
// A UNIFIED SPARSE MATRIX DATA FORMAT 
// FOR MODERN PROCESSORS WITH WIDE SIMD UNITS
// SELLC SpMV kernel modified assigning multiple threads to each row - 2D kernel
__global__ void
magma_tally3_ccgmerge_spmvsellpt_kernel_8( 
    int num_rows, 
    int blocksize,
    int T,
    magma_tally3FloatComplex * dval, 
    magma_tally3_index_t * dcolind,
    magma_tally3_index_t * drowptr,
    magma_tally3FloatComplex * d,
    magma_tally3FloatComplex * z)
{
   // T threads assigned to each row
    int idx = threadIdx.y ;     // thread in row
    int idy = threadIdx.x;      // local row
    int ldx = idx * blocksize + idy;
    int bdx = blockIdx.y * gridDim.x + blockIdx.x; // global block index
    int row = bdx * blocksize + idy;  // global row index

    extern __shared__ magma_tally3FloatComplex shared[];

    if(row < num_rows ){
        magma_tally3FloatComplex dot = MAGMA_tally3_C_MAKE(0.0, 0.0);
        int offset = drowptr[ bdx ];
        int block = blocksize * T; // total number of threads

        int max_ = (drowptr[ bdx+1 ]-offset)/block;  
            // number of elements each thread handles
        for ( int k = 0; k < max_ ; k++ ){
            magma_tally3FloatComplex val = 
                        dval[ offset + ldx + block*k ];
            int col = 
                    dcolind[ offset + ldx + block*k ];
            dot += val * d[ col ];
        }
        shared[ldx]  = dot;

        __syncthreads();
        if( idx < 4 ){
            shared[ldx]+=shared[ldx+blocksize*4];              
            __syncthreads();
            if( idx < 2 ) shared[ldx]+=shared[ldx+blocksize*2];   
            __syncthreads();
            if( idx == 0 ) {
                z[row] = 
                (shared[ldx]+shared[ldx+blocksize*1]);
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
magma_tally3_ccgmerge_spmvsellpt_kernel_16( 
    int num_rows, 
    int blocksize,
    int T,
    magma_tally3FloatComplex * dval, 
    magma_tally3_index_t * dcolind,
    magma_tally3_index_t * drowptr,
    magma_tally3FloatComplex * d,
    magma_tally3FloatComplex * z)
{
   // T threads assigned to each row
    int idx = threadIdx.y ;     // thread in row
    int idy = threadIdx.x;      // local row
    int ldx = idx * blocksize + idy;
    int bdx = blockIdx.y * gridDim.x + blockIdx.x; // global block index
    int row = bdx * blocksize + idy;  // global row index

    extern __shared__ magma_tally3FloatComplex shared[];

    if(row < num_rows ){
        magma_tally3FloatComplex dot = MAGMA_tally3_C_MAKE(0.0, 0.0);
        int offset = drowptr[ bdx ];
        int block = blocksize * T; // total number of threads

        int max_ = (drowptr[ bdx+1 ]-offset)/block;  
            // number of elements each thread handles
        for ( int k = 0; k < max_ ; k++ ){
            magma_tally3FloatComplex val = 
                        dval[ offset + ldx + block*k ];
            int col = 
                    dcolind[ offset + ldx + block*k ];
            dot += val * d[ col ];
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
                z[row] = 
                (shared[ldx]+shared[ldx+blocksize*1]);
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
magma_tally3_ccgmerge_spmvsellpt_kernel_32( 
    int num_rows, 
    int blocksize,
    int T,
    magma_tally3FloatComplex * dval, 
    magma_tally3_index_t * dcolind,
    magma_tally3_index_t * drowptr,
    magma_tally3FloatComplex * d,
    magma_tally3FloatComplex * z)
{
   // T threads assigned to each row
    int idx = threadIdx.y ;     // thread in row
    int idy = threadIdx.x;      // local row
    int ldx = idx * blocksize + idy;
    int bdx = blockIdx.y * gridDim.x + blockIdx.x; // global block index
    int row = bdx * blocksize + idy;  // global row index

    extern __shared__ magma_tally3FloatComplex shared[];

    if(row < num_rows ){
        magma_tally3FloatComplex dot = MAGMA_tally3_C_MAKE(0.0, 0.0);
        int offset = drowptr[ bdx ];
        int block = blocksize * T; // total number of threads

        int max_ = (drowptr[ bdx+1 ]-offset)/block;  
            // number of elements each thread handles
        for ( int k = 0; k < max_ ; k++ ){
            magma_tally3FloatComplex val = 
                        dval[ offset + ldx + block*k ];
            int col = 
                    dcolind[ offset + ldx + block*k ];
            dot += val * d[ col ];
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
                z[row] = 
                (shared[ldx]+shared[ldx+blocksize*1]);
            }

        }

    }
}


// kernel to handle scalars
__global__ void // rho = beta/tmp; gamma = beta;
magma_tally3_ccg_rhokernel(  
    magma_tally3FloatComplex * skp ){
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if( i==0 ){
        magma_tally3FloatComplex tmp = skp[1];
        skp[3] = tmp/skp[4];
        skp[2] = tmp;
    }
}

/**
    Purpose
    -------

    Merges the first SpmV using different formats with the dot product 
    and the computation of rho

    Arguments
    ---------

    @param[in]
    A           magma_tally3_c_matrix
                input matrix 

    @param[in]
    d1          magma_tally3FloatComplex_ptr 
                temporary vector

    @param[in]
    d2          magma_tally3FloatComplex_ptr 
                temporary vector

    @param[in]
    dd          magma_tally3FloatComplex_ptr 
                input vector d

    @param[out]
    dz          magma_tally3FloatComplex_ptr 
                input vector z

    @param[out]
    skp         magma_tally3FloatComplex_ptr 
                array for parameters ( skp[3]=rho )

    @param[in]
    queue       magma_tally3_queue_t
                Queue to execute in.

    @ingroup magma_tally3sparse_cgegpuk
    ********************************************************************/

extern "C" magma_tally3_int_t
magma_tally3_ccgmerge_spmv1(
    magma_tally3_c_matrix A,
    magma_tally3FloatComplex_ptr d1,
    magma_tally3FloatComplex_ptr d2,
    magma_tally3FloatComplex_ptr dd,
    magma_tally3FloatComplex_ptr dz,
    magma_tally3FloatComplex_ptr skp,
    magma_tally3_queue_t queue )
{
    // set queue for old dense routines
    magma_tally3_queue_t orig_queue;
    magma_tally3blasGetKernelStream( &orig_queue );

    int local_block_size=256;
    dim3 Bs( local_block_size );
    dim3 Gs( magma_tally3_ceildiv( A.num_rows, local_block_size ) );
    dim3 Gs_next;
    int Ms =  local_block_size * sizeof( magma_tally3FloatComplex ); 
    magma_tally3FloatComplex_ptr aux1 = d1, aux2 = d2;
    int b = 1;        

    if ( A.storage_type == Magma_tally3_CSR )
        magma_tally3_ccgmerge_spmvcsr_kernel<<<Gs, Bs, Ms, queue >>>
        ( A.num_rows, A.dval, A.drow, A.dcol, dd, dz, d1 );
    else if ( A.storage_type == Magma_tally3_ELLPACKT )
        magma_tally3_ccgmerge_spmvellpack_kernel<<<Gs, Bs, Ms, queue >>>
        ( A.num_rows, A.max_nnz_row, A.dval, A.dcol, dd, dz, d1 );
    else if ( A.storage_type == Magma_tally3_ELL )
        magma_tally3_ccgmerge_spmvell_kernel<<<Gs, Bs, Ms, queue >>>
        ( A.num_rows, A.max_nnz_row, A.dval, A.dcol, dd, dz, d1 );
    else if ( A.storage_type == Magma_tally3_SELLP ) {
            int num_threadssellp = A.blocksize*A.alignment;
            magma_tally3_int_t arch = magma_tally3_getdevice_arch();
            if ( arch < 200 && num_threadssellp > 256 )
                printf("error: too much shared memory requested.\n");

            dim3 block( A.blocksize, A.alignment, 1);
            int dimgrid1 = sqrt(A.numblocks);
            int dimgrid2 = magma_tally3_ceildiv( A.numblocks, dimgrid1 );

            dim3 gridsellp( dimgrid1, dimgrid2, 1);
            int Mssellp = num_threadssellp * sizeof( magma_tally3FloatComplex );

            if ( A.alignment == 8)
                magma_tally3_ccgmerge_spmvsellpt_kernel_8
                <<< gridsellp, block, Mssellp, queue >>>
                ( A.num_rows, A.blocksize, A.alignment, 
                    A.dval, A.dcol, A.drow, dd, dz);

            else if ( A.alignment == 16)
                magma_tally3_ccgmerge_spmvsellpt_kernel_16
                <<< gridsellp, block, Mssellp, queue >>>
                ( A.num_rows, A.blocksize, A.alignment, 
                    A.dval, A.dcol, A.drow, dd, dz);

            else if ( A.alignment == 32)
                magma_tally3_ccgmerge_spmvsellpt_kernel_32
                <<< gridsellp, block, Mssellp, queue >>>
                ( A.num_rows, A.blocksize, A.alignment, 
                    A.dval, A.dcol, A.drow, dd, dz);

            else
                printf("error: alignment not supported.\n");

        // in case of using SELLP, we can't efficiently merge the 
        // dot product and the first reduction loop into the SpMV kernel
        // as the SpMV grid would result in low occupancy.
        magma_tally3_ccgmerge_spmvellpackrt_kernel2<<<Gs, Bs, Ms, queue >>>
                              ( A.num_rows, dz, dd, d1 );

    }
    else if ( A.storage_type == Magma_tally3_ELLRT ) {
        // in case of using ELLRT, we need a different grid, assigning
        // threads_per_row processors to each row
        // the block size is num_threads
        // fixed values


    int num_blocks = ( (A.num_rows+A.blocksize-1)/A.blocksize);

    int num_threads = A.alignment*A.blocksize;

    int real_row_length = ((int)(A.max_nnz_row+A.alignment-1)/A.alignment)
                            *A.alignment;

    magma_tally3_int_t arch = magma_tally3_getdevice_arch();
    if ( arch < 200 && num_threads > 256 )
        printf("error: too much shared memory requested.\n");

    int dimgrid1 = sqrt(num_blocks);
    int dimgrid2 = magma_tally3_ceildiv( num_blocks, dimgrid1 );
    dim3 gridellrt( dimgrid1, dimgrid2, 1);

    int Mellrt = A.alignment * A.blocksize * sizeof( magma_tally3FloatComplex );
    // printf("launch kernel: %dx%d %d %d\n", grid.x, grid.y, num_threads , Ms);

    if ( A.alignment == 32 ) {
        magma_tally3_ccgmerge_spmvellpackrt_kernel_32
                <<< gridellrt, num_threads , Mellrt, queue >>>
                 ( A.num_rows, A.dval, A.dcol, A.drow, dd, dz, d1, 
                                                 A.alignment, real_row_length );
    }
    else if ( A.alignment == 16 ) {
        magma_tally3_ccgmerge_spmvellpackrt_kernel_16
                <<< gridellrt, num_threads , Mellrt, queue >>>
                 ( A.num_rows, A.dval, A.dcol, A.drow, dd, dz, d1, 
                                                 A.alignment, real_row_length );
    }
    else if ( A.alignment == 8 ) {
        magma_tally3_ccgmerge_spmvellpackrt_kernel_8
                <<< gridellrt, num_threads , Mellrt, queue >>>
                 ( A.num_rows, A.dval, A.dcol, A.drow, dd, dz, d1, 
                                                 A.alignment, real_row_length );
    }
    else {
        printf("error: alignment %d not supported.\n", A.alignment);
        return MAGMA_tally3_ERR_NOT_SUPPORTED;
    }
        // in case of using ELLRT, we can't efficiently merge the 
        // dot product and the first reduction loop into the SpMV kernel
        // as the SpMV grid would result in low occupancy.

        magma_tally3_ccgmerge_spmvellpackrt_kernel2<<<Gs, Bs, Ms, queue >>>
                              ( A.num_rows, dz, dd, d1 );
    }

    while( Gs.x > 1 ) {
        Gs_next.x = ( Gs.x+Bs.x-1 )/ Bs.x ;
        if ( Gs_next.x == 1 ) Gs_next.x = 2;
        magma_tally3_ccgreduce_kernel_spmv1<<< Gs_next.x/2, Bs.x/2, Ms/2 >>> 
                                        ( Gs.x,  A.num_rows, aux1, aux2 );
        Gs_next.x = Gs_next.x /2;
        Gs.x = Gs_next.x;
        b = 1 - b;
        if ( b ) { aux1 = d1; aux2 = d2; }
        else   { aux2 = d1; aux1 = d2; }
    }


    magma_tally3_ccopyvector( 1, aux1, 1, skp+4, 1 );
    dim3 Bs2( 2 );
    dim3 Gs2( 1 );
    magma_tally3_ccg_rhokernel<<<Gs2, Bs2, 0>>>( skp );

   magma_tally3blasSetKernelStream( orig_queue );
   return MAGMA_tally3_SUCCESS;
}


/* -------------------------------------------------------------------------- */

// updates x and r and computes the first part of the dot product r*r
__global__ void
magma_tally3_ccgmerge_xrbeta_kernel(  
    int n, 
    magma_tally3FloatComplex * x, 
    magma_tally3FloatComplex * r,
    magma_tally3FloatComplex * d,
    magma_tally3FloatComplex * z,
    magma_tally3FloatComplex * skp,
    magma_tally3FloatComplex * vtmp )
{

    extern __shared__ magma_tally3FloatComplex temp[]; 
    int Idx = threadIdx.x;   
    int i   = blockIdx.x * blockDim.x + Idx;

    magma_tally3FloatComplex rho = skp[3];
    magma_tally3FloatComplex mrho = MAGMA_tally3_C_MAKE( -1.0, 0.0)*rho;

    temp[ Idx ] = MAGMA_tally3_C_MAKE( 0.0, 0.0);

    if( i<n ){
        x[i] += rho * d[i] ;
        r[i] += mrho  * z[i];
        temp[ Idx ] = r[i] * r[i];
    }
    __syncthreads();
    if ( Idx < 128 ){
        temp[ Idx ] += temp[ Idx + 128 ];
    }
    __syncthreads();
    if ( Idx < 64 ){
        temp[ Idx ] += temp[ Idx + 64 ];
    }
    __syncthreads();
    #if defined(PRECISION_z) || defined(PRECISION_c)
        if( Idx < 32 ){
            temp[ Idx ] += temp[ Idx + 32 ];__syncthreads();
            temp[ Idx ] += temp[ Idx + 16 ];__syncthreads();
            temp[ Idx ] += temp[ Idx + 8 ];__syncthreads();
            temp[ Idx ] += temp[ Idx + 4 ];__syncthreads();
            temp[ Idx ] += temp[ Idx + 2 ];__syncthreads();
            temp[ Idx ] += temp[ Idx + 1 ];__syncthreads();
        }
    #endif
    #if defined(PRECISION_d)
        if( Idx < 32 ){
            volatile float *temp2 = temp;
            temp2[ Idx ] += temp2[ Idx + 32 ];
            temp2[ Idx ] += temp2[ Idx + 16 ];
            temp2[ Idx ] += temp2[ Idx + 8 ];
            temp2[ Idx ] += temp2[ Idx + 4 ];
            temp2[ Idx ] += temp2[ Idx + 2 ];
            temp2[ Idx ] += temp2[ Idx + 1 ];
        }
    #endif
    #if defined(PRECISION_s)
        if( Idx < 32 ){
            volatile float *temp2 = temp;
            temp2[ Idx ] += temp2[ Idx + 32 ];
            temp2[ Idx ] += temp2[ Idx + 16 ];
            temp2[ Idx ] += temp2[ Idx + 8 ];
            temp2[ Idx ] += temp2[ Idx + 4 ];
            temp2[ Idx ] += temp2[ Idx + 2 ];
            temp2[ Idx ] += temp2[ Idx + 1 ];
        }
    #endif

    if ( Idx == 0 ){
            vtmp[ blockIdx.x ] = temp[ 0 ];
    }

}

// kernel to handle scalars
__global__ void //alpha = beta / gamma
magma_tally3_ccg_alphabetakernel(  
    magma_tally3FloatComplex * skp )
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if( i==0 ){
        magma_tally3FloatComplex tmp1 = skp[1];
        skp[0] =  tmp1/skp[2];
        //printf("beta=%e\n", MAGMA_tally3_C_REAL(tmp1));
    }
}

// update search Krylov vector d
__global__ void
magma_tally3_ccg_d_kernel(  
    int n, 
    magma_tally3FloatComplex * skp,
    magma_tally3FloatComplex * r,
    magma_tally3FloatComplex * d )
{
  
    int i   = blockIdx.x * blockDim.x + threadIdx.x;

    magma_tally3FloatComplex alpha = skp[0];

    if( i<n ){
        d[i] = r[i] + alpha * d[i];
    }

}



/**
    Purpose
    -------

    Merges the update of r and x with the dot product and performs then 
    the update for the Krylov vector d

    Arguments
    ---------

    @param[in]
    n           int
                dimension n

    @param[in]
    d1          magma_tally3FloatComplex_ptr 
                temporary vector

    @param[in]
    d2          magma_tally3FloatComplex_ptr 
                temporary vector

    @param[in/out]
    dx          magma_tally3FloatComplex_ptr
                input vector x

    @param[in/out]
    dr          magma_tally3FloatComplex_ptr 
                input/output vector r

    @param[in]
    dd          magma_tally3FloatComplex_ptr 
                input vector d

    @param[in]
    dz          magma_tally3FloatComplex_ptr 
                input vector z
    @param[in]
    skp         magma_tally3FloatComplex_ptr 
                array for parameters

    @param[in]
    queue       magma_tally3_queue_t
                Queue to execute in.

    @ingroup magma_tally3sparse_csygpuk
    ********************************************************************/

extern "C" magma_tally3_int_t
magma_tally3_ccgmerge_xrbeta(
    int n,
    magma_tally3FloatComplex_ptr d1,
    magma_tally3FloatComplex_ptr d2,
    magma_tally3FloatComplex_ptr dx,
    magma_tally3FloatComplex_ptr dr,
    magma_tally3FloatComplex_ptr dd,
    magma_tally3FloatComplex_ptr dz, 
    magma_tally3FloatComplex_ptr skp,
    magma_tally3_queue_t queue )
{
    // set queue for old dense routines
    magma_tally3_queue_t orig_queue;
    magma_tally3blasGetKernelStream( &orig_queue );

    int local_block_size=256;
    dim3 Bs( local_block_size );
    dim3 Gs( magma_tally3_ceildiv( n, local_block_size ) );
    dim3 Gs_next;
    int Ms =  2*local_block_size * sizeof( magma_tally3FloatComplex ); 
    magma_tally3FloatComplex_ptr aux1 = d1, aux2 = d2;
    int b = 1;        
    magma_tally3_ccgmerge_xrbeta_kernel<<<Gs, Bs, Ms>>>
                                    ( n, dx, dr, dd, dz, skp, d1);  



    while( Gs.x > 1 ) {
        Gs_next.x = ( Gs.x+Bs.x-1 )/ Bs.x ;
        if ( Gs_next.x == 1 ) Gs_next.x = 2;
        magma_tally3_ccgreduce_kernel_spmv1<<< Gs_next.x/2, Bs.x/2, Ms/2 >>> 
                                    ( Gs.x, n, aux1, aux2 );
        Gs_next.x = Gs_next.x /2;
        Gs.x = Gs_next.x;
        b = 1 - b;
        if ( b ) { aux1 = d1; aux2 = d2; }
        else   { aux2 = d1; aux1 = d2; }
    }


    magma_tally3_ccopyvector( 1, aux1, 1, skp+1, 1 );
    dim3 Bs2( 2 );
    dim3 Gs2( 1 );
    magma_tally3_ccg_alphabetakernel<<<Gs2, Bs2, 0>>>( skp );

    dim3 Bs3( local_block_size );
    dim3 Gs3( magma_tally3_ceildiv( n, local_block_size ) );
    magma_tally3_ccg_d_kernel<<<Gs3, Bs3, 0>>>( n, skp, dr, dd );  

   magma_tally3blasSetKernelStream( orig_queue );
   return MAGMA_tally3_SUCCESS;
}

/* -------------------------------------------------------------------------- */

