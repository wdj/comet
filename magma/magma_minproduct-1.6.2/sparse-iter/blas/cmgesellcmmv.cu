/*
    -- MAGMA_minproduct (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from zmgesellcmmv.cu normal z -> c, Sun May  3 11:22:58 2015

*/
#include "common_magma_minproductsparse.h"

#define PRECISION_c

#define TEXTURE


// SELLP SpMV kernel 3D grid
// see paper by M. KREUTZER, G. HAGER, G WELLEIN, H. FEHSKE A. BISHOP
// A UNIFIED SPARSE MATRIX DATA FORMAT 
// FOR MODERN PROCESSORS WITH WIDE SIMD UNITS
// SELLC SpMV kernel modified assigning multiple threads to each row - 2D kernel
__global__ void 
zmgesellptmv_kernel_1_3D( 
    int num_rows, 
    int num_cols,
    int num_vecs,
    int blocksize,
    int T,
    magma_minproductFloatComplex alpha, 
    magma_minproductFloatComplex * dval, 
    magma_minproduct_index_t * dcolind,
    magma_minproduct_index_t * drowptr,
    magma_minproductFloatComplex * dx,
    magma_minproductFloatComplex beta, 
    magma_minproductFloatComplex * dy)
{
   // T threads assigned to each row
    int idx = threadIdx.x;      // local row
    int idy = threadIdx.y;      // vector
    int bdx = blockIdx.y * gridDim.x + blockIdx.x; // global block index
    int row = bdx * blocksize + idx;  // global row index


    if(row < num_rows ){
        magma_minproductFloatComplex dot = MAGMA_minproduct_C_MAKE(0.0, 0.0);
        int offset = drowptr[ bdx ];
        int max_ = (drowptr[ bdx+1 ]-offset)/blocksize;  
            // number of elements each thread handles

        for ( int k = 0; k < max_ ; k++ ){
            magma_minproductFloatComplex val = 
                        dval[ offset + idx + blocksize*k ];
            int col = 
                    dcolind[ offset + idx + blocksize*k ] ;

            dot += val * dx[ col*num_vecs+idy ];
        }
        dy[ row+idy*num_rows ] = dot*alpha + beta*dy [ row+idy*num_rows ];

    }

}


// SELLP SpMV kernel 3D grid
// see paper by M. KREUTZER, G. HAGER, G WELLEIN, H. FEHSKE A. BISHOP
// A UNIFIED SPARSE MATRIX DATA FORMAT 
// FOR MODERN PROCESSORS WITH WIDE SIMD UNITS
// SELLC SpMV kernel modified assigning multiple threads to each row - 2D kernel
__global__ void 
zmgesellptmv_kernel_4_3D( 
    int num_rows, 
    int num_cols,
    int num_vecs,
    int blocksize,
    int T,
    magma_minproductFloatComplex alpha, 
    magma_minproductFloatComplex * dval, 
    magma_minproduct_index_t * dcolind,
    magma_minproduct_index_t * drowptr,
    magma_minproductFloatComplex * dx,
    magma_minproductFloatComplex beta, 
    magma_minproductFloatComplex * dy)
{
   // T threads assigned to each row
    int idx = threadIdx.y ;     // thread in row
    int idy = threadIdx.x;      // local row
    int idz = threadIdx.z;      // vector
    int ldx = idx * blocksize + idy;
    int ldz = idz * blocksize * T + idx * blocksize + idy;
    int bdx = blockIdx.y * gridDim.x + blockIdx.x; // global block index
    int row = bdx * blocksize + idy;  // global row index
    int vec = idz*num_rows;

    extern __shared__ magma_minproductFloatComplex shared[];


    if(row < num_rows ){
        magma_minproductFloatComplex dot = MAGMA_minproduct_C_MAKE(0.0, 0.0);
        int offset = drowptr[ bdx ];
        int block = blocksize * T; // total number of threads

        int max_ = (drowptr[ bdx+1 ]-offset)/block;  
            // number of elements each thread handles



        for ( int k = 0; k < max_ ; k++ ){
            magma_minproductFloatComplex val = 
                        dval[ offset + ldx + block*k ];
            int col = 
                    dcolind[ offset + ldx + block*k ] ;

            dot += val * dx[ col+vec ];
        }
        shared[ldz]  = dot;

        __syncthreads();
        if( idx < 2 ){
            shared[ldz]+=shared[ldz+blocksize*2];               
            __syncthreads();
            if( idx == 0 ) {
                dy[row+vec] = 
                (shared[ldz]+shared[ldz+blocksize*1])*alpha 
                                            + beta*dy [row+vec];
            }

        }

    }

}


// SELLP SpMV kernel 3D grid
// see paper by M. KREUTZER, G. HAGER, G WELLEIN, H. FEHSKE A. BISHOP
// A UNIFIED SPARSE MATRIX DATA FORMAT 
// FOR MODERN PROCESSORS WITH WIDE SIMD UNITS
// SELLC SpMV kernel modified assigning multiple threads to each row - 2D kernel
__global__ void 
zmgesellptmv_kernel_8_3D( 
    int num_rows, 
    int num_cols,
    int num_vecs,
    int blocksize,
    int T,
    magma_minproductFloatComplex alpha, 
    magma_minproductFloatComplex * dval, 
    magma_minproduct_index_t * dcolind,
    magma_minproduct_index_t * drowptr,
    const magma_minproductFloatComplex * __restrict__ dx,
    magma_minproductFloatComplex beta, 
    magma_minproductFloatComplex * dy)
{
   // T threads assigned to each row
    int idx = threadIdx.y ;     // thread in row
    int idy = threadIdx.x;      // local row
    int idz = threadIdx.z;      // vector
    int ldx = idx * blocksize + idy;
    int ldz = idz * blocksize * T + idx * blocksize + idy;
    int bdx = blockIdx.y * gridDim.x + blockIdx.x; // global block index
    int row = bdx * blocksize + idy;  // global row index
    int vec = idz*num_rows;

    extern __shared__ magma_minproductFloatComplex shared[];


    if(row < num_rows ){
        magma_minproductFloatComplex dot = MAGMA_minproduct_C_MAKE(0.0, 0.0);
        int offset = drowptr[ bdx ];
        int block = blocksize * T; // total number of threads

        int max_ = (drowptr[ bdx+1 ]-offset)/block;  
            // number of elements each thread handles



        for ( int k = 0; k < max_ ; k++ ){
            magma_minproductFloatComplex val = 
                        dval[ offset + ldx + block*k ];
            int col = 
                    dcolind[ offset + ldx + block*k ] ;

            dot += val * dx[ col+vec ];
        }
        shared[ldz]  = dot;

        __syncthreads();
        if( idx < 4 ){
            shared[ldz]+=shared[ldz+blocksize*4];               
            __syncthreads();
            if( idx < 2 ) shared[ldz]+=shared[ldz+blocksize*2];   
            __syncthreads();
            if( idx == 0 ) {
                dy[row+vec] = 
                (shared[ldz]+shared[ldz+blocksize*1])*alpha 
                                            + beta*dy [row+vec];
            }

        }

    }

}


// SELLP SpMV kernel 3D grid
// see paper by M. KREUTZER, G. HAGER, G WELLEIN, H. FEHSKE A. BISHOP
// A UNIFIED SPARSE MATRIX DATA FORMAT 
// FOR MODERN PROCESSORS WITH WIDE SIMD UNITS
// SELLC SpMV kernel modified assigning multiple threads to each row - 2D kernel
__global__ void 
zmgesellptmv_kernel_16_3D( 
    int num_rows, 
    int num_cols,
    int num_vecs,
    int blocksize,
    int T,
    magma_minproductFloatComplex alpha, 
    magma_minproductFloatComplex * dval, 
    magma_minproduct_index_t * dcolind,
    magma_minproduct_index_t * drowptr,
    magma_minproductFloatComplex * dx,
    magma_minproductFloatComplex beta, 
    magma_minproductFloatComplex * dy)
{
   // T threads assigned to each row
    int idx = threadIdx.y ;     // thread in row
    int idy = threadIdx.x;      // local row
    int idz = threadIdx.z;      // vector
    int ldx = idx * blocksize + idy;
    int ldz = idz * blocksize * T + idx * blocksize + idy;
    int bdx = blockIdx.y * gridDim.x + blockIdx.x; // global block index
    int row = bdx * blocksize + idy;  // global row index
    int vec = idz*num_rows;

    extern __shared__ magma_minproductFloatComplex shared[];


    if(row < num_rows ){
        magma_minproductFloatComplex dot = MAGMA_minproduct_C_MAKE(0.0, 0.0);
        int offset = drowptr[ bdx ];
        int block = blocksize * T; // total number of threads

        int max_ = (drowptr[ bdx+1 ]-offset)/block;  
            // number of elements each thread handles
        for ( int k = 0; k < max_ ; k++ ){
            magma_minproductFloatComplex val = 
                        dval[ offset + ldx + block*k ];
            int col = 
                    dcolind[ offset + ldx + block*k ];

            dot += val * dx[ col+vec ];
        }
        shared[ldz]  = dot;

        __syncthreads();
        if( idx < 8 ){
            shared[ldz]+=shared[ldz+blocksize*8];              
            __syncthreads();
            if( idx < 4 ) shared[ldz]+=shared[ldz+blocksize*4];   
            __syncthreads();
            if( idx < 2 ) shared[ldz]+=shared[ldz+blocksize*2];   
            __syncthreads();
            if( idx == 0 ) {
                dy[row+vec] = 
                (shared[ldz]+shared[ldz+blocksize*1])*alpha 
                                            + beta*dy [row+vec];
            }

        }

    }

}


// SELLP SpMV kernel 3D grid
// see paper by M. KREUTZER, G. HAGER, G WELLEIN, H. FEHSKE A. BISHOP
// A UNIFIED SPARSE MATRIX DATA FORMAT 
// FOR MODERN PROCESSORS WITH WIDE SIMD UNITS
// SELLC SpMV kernel modified assigning multiple threads to each row - 2D kernel
__global__ void 
zmgesellptmv_kernel_32_3D( 
    int num_rows, 
    int num_cols,
    int num_vecs,
    int blocksize,
    int T,
    magma_minproductFloatComplex alpha, 
    magma_minproductFloatComplex * dval, 
    magma_minproduct_index_t * dcolind,
    magma_minproduct_index_t * drowptr,
    magma_minproductFloatComplex * dx,
    magma_minproductFloatComplex beta, 
    magma_minproductFloatComplex * dy)
{
   // T threads assigned to each row
    int idx = threadIdx.y ;     // thread in row
    int idy = threadIdx.x;      // local row
    int idz = threadIdx.z;      // vector
    int ldx = idx * blocksize + idy;
    int ldz = idz * blocksize * T + idx * blocksize + idy;
    int bdx = blockIdx.y * gridDim.x + blockIdx.x; // global block index
    int row = bdx * blocksize + idy;  // global row index
    int vec = idz*num_rows;

    extern __shared__ magma_minproductFloatComplex shared[];


    if(row < num_rows ){
        magma_minproductFloatComplex dot = MAGMA_minproduct_C_MAKE(0.0, 0.0);
        int offset = drowptr[ bdx ];
        int block = blocksize * T; // total number of threads

        int max_ = (drowptr[ bdx+1 ]-offset)/block;  
            // number of elements each thread handles
        for ( int k = 0; k < max_ ; k++ ){
            magma_minproductFloatComplex val = 
                        dval[ offset + ldx + block*k ];
            int col = 
                    dcolind[ offset + ldx + block*k ];

            dot += val * dx[ col+vec ];
        }
        shared[ldz]  = dot;

        __syncthreads();
        if( idx < 16 ){
            shared[ldz]+=shared[ldz+blocksize*16];              
            __syncthreads();
            if( idx < 8 ) shared[ldz]+=shared[ldz+blocksize*8];  
            __syncthreads();
            if( idx < 4 ) shared[ldz]+=shared[ldz+blocksize*4];   
            __syncthreads();
            if( idx < 2 ) shared[ldz]+=shared[ldz+blocksize*2];   
            __syncthreads();
            if( idx == 0 ) {
                dy[row+vec] = 
                (shared[ldz]+shared[ldz+blocksize*1])*alpha 
                                            + beta*dy [row+vec];
            }

        }

    }

}

/************************* same but using texture mem *************************/



// SELLP SpMV kernel 2D grid - for large number of vectors
// see paper by M. KREUTZER, G. HAGER, G WELLEIN, H. FEHSKE A. BISHOP
// A UNIFIED SPARSE MATRIX DATA FORMAT 
// FOR MODERN PROCESSORS WITH WIDE SIMD UNITS
// SELLC SpMV kernel modified assigning multiple threads to each row - 2D kernel
__global__ void 
zmgesellptmv_kernel_1_3D_tex( 
    int num_rows, 
    int num_cols,
    int num_vecs,
    int blocksize,
    int T,
    magma_minproductFloatComplex alpha, 
    magma_minproductFloatComplex * dval, 
    magma_minproduct_index_t * dcolind,
    magma_minproduct_index_t * drowptr,
    cudaTextureObject_t texdx,
    magma_minproductFloatComplex beta, 
    magma_minproductFloatComplex * dy)
{
#if defined(PRECISION_d) && defined(TEXTURE) && (__CUDA_ARCH__ >= 300)

    int idx = threadIdx.x;      // local row
    int idy = threadIdx.y;      // vector
    int bdx = blockIdx.y * gridDim.x + blockIdx.x; // global block index
    int row = bdx * blocksize + idx;  // global row index

    if(row < num_rows ){
        magma_minproductFloatComplex dot1 = MAGMA_minproduct_C_MAKE(0.0, 0.0);
        magma_minproductFloatComplex dot2 = MAGMA_minproduct_C_MAKE(0.0, 0.0);
        int offset = drowptr[ bdx ];
        int max_ = (drowptr[ bdx+1 ]-offset)/blocksize;  
            // number of elements each thread handles

        for ( int k = 0; k < max_ ; k++ ){
            magma_minproductFloatComplex val = 
                        dval[ offset + idx + blocksize*k ];
            int col = 
                    num_vecs * dcolind[ offset + idx + blocksize*k ] ;

            int4 v = tex1Dfetch<int4>(texdx, col/2 + idy );
            dot1 += val * __hiloint2float(v.y, v.x);
            dot2 += val * __hiloint2float(v.w, v.z);
        }
        dy[row+num_rows*idy*2] = 
                            dot1*alpha
                            + beta*dy [row*num_vecs+idy*2];
        dy[row+num_rows*idy*2+num_rows] = 
                            dot2*alpha
                            + beta*dy [row*num_vecs+idy*2+1];
    }
#endif
}


// SELLP SpMV kernel 3D grid
// see paper by M. KREUTZER, G. HAGER, G WELLEIN, H. FEHSKE A. BISHOP
// A UNIFIED SPARSE MATRIX DATA FORMAT 
// FOR MODERN PROCESSORS WITH WIDE SIMD UNITS
// SELLC SpMV kernel modified assigning multiple threads to each row - 2D kernel
__global__ void 
zmgesellptmv_kernel_4_3D_tex( 
    int num_rows, 
    int num_cols,
    int num_vecs,
    int blocksize,
    int T,
    magma_minproductFloatComplex alpha, 
    magma_minproductFloatComplex * dval, 
    magma_minproduct_index_t * dcolind,
    magma_minproduct_index_t * drowptr,
    cudaTextureObject_t texdx,
    magma_minproductFloatComplex beta, 
    magma_minproductFloatComplex * dy)
{
#if defined(PRECISION_d) && defined(TEXTURE) && (__CUDA_ARCH__ >= 300)
   // T threads assigned to each row
    int idx = threadIdx.y ;     // thread in row
    int idy = threadIdx.x;      // local row
    int idz = threadIdx.z;      // vector
    int ldx = idx * blocksize + idy;
    int ldz = idz * blocksize * T + idx * blocksize + idy;
    int bdx = blockIdx.y * gridDim.x + blockIdx.x; // global block index
    int row = bdx * blocksize + idy;  // global row index
    int sv = num_vecs/2 * blocksize * T;

    extern __shared__ magma_minproductFloatComplex shared[];


    if(row < num_rows ){
        magma_minproductFloatComplex dot1 = MAGMA_minproduct_C_MAKE(0.0, 0.0);
        magma_minproductFloatComplex dot2 = MAGMA_minproduct_C_MAKE(0.0, 0.0);
        int offset = drowptr[ bdx ];
        int block = blocksize * T; // total number of threads

        int max_ = (drowptr[ bdx+1 ]-offset)/block;  
            // number of elements each thread handles



        for ( int k = 0; k < max_ ; k++ ){
            magma_minproductFloatComplex val = 
                        dval[ offset + ldx + block*k ];
            int col = 
                    num_vecs * dcolind[ offset + ldx + block*k ] ;

            int4 v = tex1Dfetch<int4>(texdx, col/2 + idz );
            dot1 += val * __hiloint2float(v.y, v.x);
            dot2 += val * __hiloint2float(v.w, v.z);
        }
        shared[ldz]  = dot1;
        shared[ldz+sv]  = dot2;

        __syncthreads();
        if( idx < 2 ){
            shared[ldz]+=shared[ldz+blocksize*2];    
            shared[ldz+sv]+=shared[ldz+sv+blocksize*2];               
            __syncthreads();
            if( idx == 0 ) {
                dy[row+num_rows*idz*2] = 
                (shared[ldz]+shared[ldz+blocksize*1])*alpha
                                            + beta*dy [row*num_vecs+idz*2];
                dy[row+num_rows*idz*2+num_rows] = 
                (shared[ldz+sv]+shared[ldz+sv+blocksize*1])*alpha
                                            + beta*dy [row*num_vecs+idz*2+1];
            }

        }

    }
#endif
}
 

// SELLP SpMV kernel 3D grid
// see paper by M. KREUTZER, G. HAGER, G WELLEIN, H. FEHSKE A. BISHOP
// A UNIFIED SPARSE MATRIX DATA FORMAT 
// FOR MODERN PROCESSORS WITH WIDE SIMD UNITS
// SELLC SpMV kernel modified assigning multiple threads to each row - 2D kernel
__global__ void 
zmgesellptmv_kernel_8_3D_tex( 
    int num_rows, 
    int num_cols,
    int num_vecs,
    int blocksize,
    int T,
    magma_minproductFloatComplex alpha, 
    magma_minproductFloatComplex * dval, 
    magma_minproduct_index_t * dcolind,
    magma_minproduct_index_t * drowptr,
    cudaTextureObject_t texdx,
    magma_minproductFloatComplex beta, 
    magma_minproductFloatComplex * dy)
{
#if defined(PRECISION_d) && defined(TEXTURE) && (__CUDA_ARCH__ >= 300)
   // T threads assigned to each row
    int idx = threadIdx.y ;     // thread in row
    int idy = threadIdx.x;      // local row
    int idz = threadIdx.z;      // vector
    int ldx = idx * blocksize + idy;
    int ldz = idz * blocksize * T + idx * blocksize + idy;
    int bdx = blockIdx.y * gridDim.x + blockIdx.x; // global block index
    int row = bdx * blocksize + idy;  // global row index
    int sv = num_vecs/2 * blocksize * T;

    extern __shared__ magma_minproductFloatComplex shared[];


    if(row < num_rows ){
        magma_minproductFloatComplex dot1 = MAGMA_minproduct_C_MAKE(0.0, 0.0);
        magma_minproductFloatComplex dot2 = MAGMA_minproduct_C_MAKE(0.0, 0.0);
        int offset = drowptr[ bdx ];
        int block = blocksize * T; // total number of threads

        int max_ = (drowptr[ bdx+1 ]-offset)/block;  
            // number of elements each thread handles



        for ( int k = 0; k < max_ ; k++ ){
            magma_minproductFloatComplex val = 
                        dval[ offset + ldx + block*k ];
            int col = 
                    num_vecs * dcolind[ offset + ldx + block*k ] ;

            int4 v = tex1Dfetch<int4>(texdx, col/2 + idz );
            dot1 += val * __hiloint2float(v.y, v.x);
            dot2 += val * __hiloint2float(v.w, v.z);
        }
        shared[ldz]  = dot1;
        shared[ldz+sv]  = dot2;

        __syncthreads();
        if( idx < 4 ){
            shared[ldz]+=shared[ldz+blocksize*4];    
            shared[ldz+sv]+=shared[ldz+sv+blocksize*4];               
            __syncthreads();
            if( idx < 2 ){
                shared[ldz]+=shared[ldz+blocksize*2];   
                shared[ldz+sv]+=shared[ldz+sv+blocksize*2];   
            }
            __syncthreads();
            if( idx == 0 ) {
                dy[row+num_rows*idz*2] = 
                (shared[ldz]+shared[ldz+blocksize*1])*alpha
                                            + beta*dy [row*num_vecs+idz*2];
                dy[row+num_rows*idz*2+num_rows] = 
                (shared[ldz+sv]+shared[ldz+sv+blocksize*1])*alpha
                                            + beta*dy [row*num_vecs+idz*2+1];
            }

        }

    }
#endif
}


// SELLP SpMV kernel 3D grid
// see paper by M. KREUTZER, G. HAGER, G WELLEIN, H. FEHSKE A. BISHOP
// A UNIFIED SPARSE MATRIX DATA FORMAT 
// FOR MODERN PROCESSORS WITH WIDE SIMD UNITS
// SELLC SpMV kernel modified assigning multiple threads to each row - 2D kernel
__global__ void 
zmgesellptmv_kernel_16_3D_tex( 
    int num_rows, 
    int num_cols,
    int num_vecs,
    int blocksize,
    int T,
    magma_minproductFloatComplex alpha, 
    magma_minproductFloatComplex * dval, 
    magma_minproduct_index_t * dcolind,
    magma_minproduct_index_t * drowptr,
    cudaTextureObject_t texdx,
    magma_minproductFloatComplex beta, 
    magma_minproductFloatComplex * dy)
{
#if defined(PRECISION_d) && defined(TEXTURE) && (__CUDA_ARCH__ >= 300)
   // T threads assigned to each row
    int idx = threadIdx.y ;     // thread in row
    int idy = threadIdx.x;      // local row
    int idz = threadIdx.z;      // vector
    int ldx = idx * blocksize + idy;
    int ldz = idz * blocksize * T + idx * blocksize + idy;
    int bdx = blockIdx.y * gridDim.x + blockIdx.x; // global block index
    int row = bdx * blocksize + idy;  // global row index
    int sv = num_vecs/2 * blocksize * T;

    extern __shared__ magma_minproductFloatComplex shared[];


    if(row < num_rows ){
        magma_minproductFloatComplex dot1 = MAGMA_minproduct_C_MAKE(0.0, 0.0);
        magma_minproductFloatComplex dot2 = MAGMA_minproduct_C_MAKE(0.0, 0.0);
        int offset = drowptr[ bdx ];
        int block = blocksize * T; // total number of threads

        int max_ = (drowptr[ bdx+1 ]-offset)/block;  
            // number of elements each thread handles



        for ( int k = 0; k < max_ ; k++ ){
            magma_minproductFloatComplex val = 
                        dval[ offset + ldx + block*k ];
            int col = 
                    num_vecs * dcolind[ offset + ldx + block*k ] ;

            int4 v = tex1Dfetch<int4>(texdx, col/2 + idz );
            dot1 += val * __hiloint2float(v.y, v.x);
            dot2 += val * __hiloint2float(v.w, v.z);
        }
        shared[ldz]  = dot1;
        shared[ldz+sv]  = dot2;

        __syncthreads();
        if( idx < 8 ){
            shared[ldz]+=shared[ldz+blocksize*8];    
            shared[ldz+sv]+=shared[ldz+sv+blocksize*8];               
            __syncthreads();
            if( idx < 4 ){
                shared[ldz]+=shared[ldz+blocksize*4];   
                shared[ldz+sv]+=shared[ldz+sv+blocksize*4];   
            }
            if( idx < 2 ){
                shared[ldz]+=shared[ldz+blocksize*2];   
                shared[ldz+sv]+=shared[ldz+sv+blocksize*2];   
            }
            __syncthreads();
            if( idx == 0 ) {
                dy[row+num_rows*idz*2] = 
                (shared[ldz]+shared[ldz+blocksize*1])*alpha
                                            + beta*dy [row*num_vecs+idz*2];
                dy[row+num_rows*idz*2+num_rows] = 
                (shared[ldz+sv]+shared[ldz+sv+blocksize*1])*alpha
                                            + beta*dy [row*num_vecs+idz*2+1];
            }

        }

    }
#endif
}


// SELLP SpMV kernel 3D grid
// see paper by M. KREUTZER, G. HAGER, G WELLEIN, H. FEHSKE A. BISHOP
// A UNIFIED SPARSE MATRIX DATA FORMAT 
// FOR MODERN PROCESSORS WITH WIDE SIMD UNITS
// SELLC SpMV kernel modified assigning multiple threads to each row - 2D kernel
__global__ void 
zmgesellptmv_kernel_32_3D_tex( 
    int num_rows, 
    int num_cols,
    int num_vecs,
    int blocksize,
    int T,
    magma_minproductFloatComplex alpha, 
    magma_minproductFloatComplex * dval, 
    magma_minproduct_index_t * dcolind,
    magma_minproduct_index_t * drowptr,
    cudaTextureObject_t texdx,
    magma_minproductFloatComplex beta, 
    magma_minproductFloatComplex * dy)
{
#if defined(PRECISION_d) && defined(TEXTURE) && (__CUDA_ARCH__ >= 300)
   // T threads assigned to each row
    int idx = threadIdx.y ;     // thread in row
    int idy = threadIdx.x;      // local row
    int idz = threadIdx.z;      // vector
    int ldx = idx * blocksize + idy;
    int ldz = idz * blocksize * T + idx * blocksize + idy;
    int bdx = blockIdx.y * gridDim.x + blockIdx.x; // global block index
    int row = bdx * blocksize + idy;  // global row index
    int sv = num_vecs/2 * blocksize * T;

    extern __shared__ magma_minproductFloatComplex shared[];


    if(row < num_rows ){
        magma_minproductFloatComplex dot1 = MAGMA_minproduct_C_MAKE(0.0, 0.0);
        magma_minproductFloatComplex dot2 = MAGMA_minproduct_C_MAKE(0.0, 0.0);
        int offset = drowptr[ bdx ];
        int block = blocksize * T; // total number of threads

        int max_ = (drowptr[ bdx+1 ]-offset)/block;  
            // number of elements each thread handles



        for ( int k = 0; k < max_ ; k++ ){
            magma_minproductFloatComplex val = 
                        dval[ offset + ldx + block*k ];
            int col = 
                    num_vecs * dcolind[ offset + ldx + block*k ] ;

            int4 v = tex1Dfetch<int4>(texdx, col/2 + idz );
            dot1 += val * __hiloint2float(v.y, v.x);
            dot2 += val * __hiloint2float(v.w, v.z);
        }
        shared[ldz]  = dot1;
        shared[ldz+sv]  = dot2;

        __syncthreads();
        if( idx < 16 ){
            shared[ldz]+=shared[ldz+blocksize*16];    
            shared[ldz+sv]+=shared[ldz+sv+blocksize*16];               
            __syncthreads();
            if( idx < 8 ){
                shared[ldz]+=shared[ldz+blocksize*8];   
                shared[ldz+sv]+=shared[ldz+sv+blocksize*8];   
            }
            if( idx < 4 ){
                shared[ldz]+=shared[ldz+blocksize*4];   
                shared[ldz+sv]+=shared[ldz+sv+blocksize*4];   
            }
            if( idx < 2 ){
                shared[ldz]+=shared[ldz+blocksize*2];   
                shared[ldz+sv]+=shared[ldz+sv+blocksize*2];   
            }
            __syncthreads();
            if( idx == 0 ) {
                dy[row+num_rows*idz*2] = 
                (shared[ldz]+shared[ldz+blocksize*1])*alpha
                                            + beta*dy [row*num_vecs+idz*2];
                dy[row+num_rows*idz*2+num_rows] = 
                (shared[ldz+sv]+shared[ldz+sv+blocksize*1])*alpha
                                            + beta*dy [row*num_vecs+idz*2+1];
            }

        }

    }
#endif
}

//***************** routines for beta = 0 ************************************//


// SELLP SpMV kernel 2D grid - for large number of vectors
// see paper by M. KREUTZER, G. HAGER, G WELLEIN, H. FEHSKE A. BISHOP
// A UNIFIED SPARSE MATRIX DATA FORMAT 
// FOR MODERN PROCESSORS WITH WIDE SIMD UNITS
// SELLC SpMV kernel modified assigning multiple threads to each row - 2D kernel
__global__ void 
zmgesellptmv_kernel_1_3D_texb( 
    int num_rows, 
    int num_cols,
    int num_vecs,
    int blocksize,
    int T,
    magma_minproductFloatComplex alpha, 
    magma_minproductFloatComplex * dval, 
    magma_minproduct_index_t * dcolind,
    magma_minproduct_index_t * drowptr,
    cudaTextureObject_t texdx,
    magma_minproductFloatComplex * dy)
{
#if defined(PRECISION_d) && defined(TEXTURE) && (__CUDA_ARCH__ >= 300)

    int idx = threadIdx.x;      // local row
    int idy = threadIdx.y;      // vector
    int bdx = blockIdx.y * gridDim.x + blockIdx.x; // global block index
    int row = bdx * blocksize + idx;  // global row index

    if(row < num_rows ){
        magma_minproductFloatComplex dot1 = MAGMA_minproduct_C_MAKE(0.0, 0.0);
        magma_minproductFloatComplex dot2 = MAGMA_minproduct_C_MAKE(0.0, 0.0);
        int offset = drowptr[ bdx ];
        int block = blocksize * T; // total number of threads

        int max_ = (drowptr[ bdx+1 ]-offset)/block;  
            // number of elements each thread handles

        for ( int k = 0; k < max_ ; k++ ){
            magma_minproductFloatComplex val = 
                        dval[ offset + idx + blocksize*k ];
            int col = 
                    num_vecs * dcolind[ offset + idx + blocksize*k ] ;

            int4 v = tex1Dfetch<int4>(texdx, col/2 + idy );
            dot1 += val * __hiloint2float(v.y, v.x);
            dot2 += val * __hiloint2float(v.w, v.z);
        }
        dy[row+num_rows*idy*2] = 
                            dot1*alpha;
        dy[row+num_rows*idy*2+num_rows] = 
                            dot2*alpha;
    }
#endif
}


// SELLP SpMV kernel 3D grid
// see paper by M. KREUTZER, G. HAGER, G WELLEIN, H. FEHSKE A. BISHOP
// A UNIFIED SPARSE MATRIX DATA FORMAT 
// FOR MODERN PROCESSORS WITH WIDE SIMD UNITS
// SELLC SpMV kernel modified assigning multiple threads to each row - 2D kernel
__global__ void 
zmgesellptmv_kernel_4_3D_texb( 
    int num_rows, 
    int num_cols,
    int num_vecs,
    int blocksize,
    int T,
    magma_minproductFloatComplex alpha, 
    magma_minproductFloatComplex * dval, 
    magma_minproduct_index_t * dcolind,
    magma_minproduct_index_t * drowptr,
    cudaTextureObject_t texdx,
    magma_minproductFloatComplex * dy)
{
#if defined(PRECISION_d) && defined(TEXTURE) && (__CUDA_ARCH__ >= 300)
   // T threads assigned to each row
    int idx = threadIdx.y ;     // thread in row
    int idy = threadIdx.x;      // local row
    int idz = threadIdx.z;      // vector
    int ldx = idx * blocksize + idy;
    int ldz = idz * blocksize * T + idx * blocksize + idy;
    int bdx = blockIdx.y * gridDim.x + blockIdx.x; // global block index
    int row = bdx * blocksize + idy;  // global row index
    int sv = num_vecs/2 * blocksize * T;

    extern __shared__ magma_minproductFloatComplex shared[];


    if(row < num_rows ){
        magma_minproductFloatComplex dot1 = MAGMA_minproduct_C_MAKE(0.0, 0.0);
        magma_minproductFloatComplex dot2 = MAGMA_minproduct_C_MAKE(0.0, 0.0);
        int offset = drowptr[ bdx ];
        int block = blocksize * T; // total number of threads

        int max_ = (drowptr[ bdx+1 ]-offset)/block;  
            // number of elements each thread handles



        for ( int k = 0; k < max_ ; k++ ){
            magma_minproductFloatComplex val = 
                        dval[ offset + ldx + block*k ];
            int col = 
                    num_vecs * dcolind[ offset + ldx + block*k ] ;

            int4 v = tex1Dfetch<int4>(texdx, col/2 + idz );
            dot1 += val * __hiloint2float(v.y, v.x);
            dot2 += val * __hiloint2float(v.w, v.z);
        }
        shared[ldz]  = dot1;
        shared[ldz+sv]  = dot2;

        __syncthreads();
        if( idx < 2 ){
            shared[ldz]+=shared[ldz+blocksize*2];    
            shared[ldz+sv]+=shared[ldz+sv+blocksize*2];               
            __syncthreads();
            if( idx == 0 ) {
                dy[row+num_rows*idz*2] = 
                (shared[ldz]+shared[ldz+blocksize*1])*alpha;
                dy[row+num_rows*idz*2+num_rows] = 
                (shared[ldz+sv]+shared[ldz+sv+blocksize*1])*alpha;
            }

        }

    }
#endif
}
 

// SELLP SpMV kernel 3D grid
// see paper by M. KREUTZER, G. HAGER, G WELLEIN, H. FEHSKE A. BISHOP
// A UNIFIED SPARSE MATRIX DATA FORMAT 
// FOR MODERN PROCESSORS WITH WIDE SIMD UNITS
// SELLC SpMV kernel modified assigning multiple threads to each row - 2D kernel
__global__ void 
zmgesellptmv_kernel_8_3D_texb( 
    int num_rows, 
    int num_cols,
    int num_vecs,
    int blocksize,
    int T,
    magma_minproductFloatComplex alpha, 
    magma_minproductFloatComplex * dval, 
    magma_minproduct_index_t * dcolind,
    magma_minproduct_index_t * drowptr,
    cudaTextureObject_t texdx,
    magma_minproductFloatComplex * dy)
{
#if defined(PRECISION_d) && defined(TEXTURE) && (__CUDA_ARCH__ >= 300)
   // T threads assigned to each row
    int idx = threadIdx.y ;     // thread in row
    int idy = threadIdx.x;      // local row
    int idz = threadIdx.z;      // vector
    int ldx = idx * blocksize + idy;
    int ldz = idz * blocksize * T + idx * blocksize + idy;
    int bdx = blockIdx.y * gridDim.x + blockIdx.x; // global block index
    int row = bdx * blocksize + idy;  // global row index
    int sv = num_vecs/2 * blocksize * T;

    extern __shared__ magma_minproductFloatComplex shared[];


    if(row < num_rows ){
        magma_minproductFloatComplex dot1 = MAGMA_minproduct_C_MAKE(0.0, 0.0);
        magma_minproductFloatComplex dot2 = MAGMA_minproduct_C_MAKE(0.0, 0.0);
        int offset = drowptr[ bdx ];
        int block = blocksize * T; // total number of threads

        int max_ = (drowptr[ bdx+1 ]-offset)/block;  
            // number of elements each thread handles



        for ( int k = 0; k < max_ ; k++ ){
            magma_minproductFloatComplex val = 
                        dval[ offset + ldx + block*k ];
            int col = 
                    num_vecs * dcolind[ offset + ldx + block*k ] ;

            int4 v = tex1Dfetch<int4>(texdx, col/2 + idz );
            dot1 += val * __hiloint2float(v.y, v.x);
            dot2 += val * __hiloint2float(v.w, v.z);
        }
        shared[ldz]  = dot1;
        shared[ldz+sv]  = dot2;

        __syncthreads();
        if( idx < 4 ){
            shared[ldz]+=shared[ldz+blocksize*4];    
            shared[ldz+sv]+=shared[ldz+sv+blocksize*4];               
            __syncthreads();
            if( idx < 2 ){
                shared[ldz]+=shared[ldz+blocksize*2];   
                shared[ldz+sv]+=shared[ldz+sv+blocksize*2];   
            }
            __syncthreads();
            if( idx == 0 ) {
                dy[row+num_rows*idz*2] = 
                (shared[ldz]+shared[ldz+blocksize*1])*alpha;

                dy[row+num_rows*idz*2+num_rows] = 
                (shared[ldz+sv]+shared[ldz+sv+blocksize*1])*alpha;
            }

        }

    }
#endif
}


// SELLP SpMV kernel 3D grid
// see paper by M. KREUTZER, G. HAGER, G WELLEIN, H. FEHSKE A. BISHOP
// A UNIFIED SPARSE MATRIX DATA FORMAT 
// FOR MODERN PROCESSORS WITH WIDE SIMD UNITS
// SELLC SpMV kernel modified assigning multiple threads to each row - 2D kernel
__global__ void 
zmgesellptmv_kernel_16_3D_texb(
    int num_rows, 
    int num_cols,
    int num_vecs,
    int blocksize,
    int T,
    magma_minproductFloatComplex alpha, 
    magma_minproductFloatComplex * dval, 
    magma_minproduct_index_t * dcolind,
    magma_minproduct_index_t * drowptr,
    cudaTextureObject_t texdx,
    magma_minproductFloatComplex * dy)
{
#if defined(PRECISION_d) && defined(TEXTURE) && (__CUDA_ARCH__ >= 300)
   // T threads assigned to each row
    int idx = threadIdx.y ;     // thread in row
    int idy = threadIdx.x;      // local row
    int idz = threadIdx.z;      // vector
    int ldx = idx * blocksize + idy;
    int ldz = idz * blocksize * T + idx * blocksize + idy;
    int bdx = blockIdx.y * gridDim.x + blockIdx.x; // global block index
    int row = bdx * blocksize + idy;  // global row index
    int sv = num_vecs/2 * blocksize * T;

    extern __shared__ magma_minproductFloatComplex shared[];


    if(row < num_rows ){
        magma_minproductFloatComplex dot1 = MAGMA_minproduct_C_MAKE(0.0, 0.0);
        magma_minproductFloatComplex dot2 = MAGMA_minproduct_C_MAKE(0.0, 0.0);
        int offset = drowptr[ bdx ];
        int block = blocksize * T; // total number of threads

        int max_ = (drowptr[ bdx+1 ]-offset)/block;  
            // number of elements each thread handles



        for ( int k = 0; k < max_ ; k++ ){
            magma_minproductFloatComplex val = 
                        dval[ offset + ldx + block*k ];
            int col = 
                    num_vecs * dcolind[ offset + ldx + block*k ] ;

            int4 v = tex1Dfetch<int4>(texdx, col/2 + idz );
            dot1 += val * __hiloint2float(v.y, v.x);
            dot2 += val * __hiloint2float(v.w, v.z);
        }
        shared[ldz]  = dot1;
        shared[ldz+sv]  = dot2;

        __syncthreads();
        if( idx < 8 ){
            shared[ldz]+=shared[ldz+blocksize*8];    
            shared[ldz+sv]+=shared[ldz+sv+blocksize*8];               
            __syncthreads();
            if( idx < 4 ){
                shared[ldz]+=shared[ldz+blocksize*4];   
                shared[ldz+sv]+=shared[ldz+sv+blocksize*4];   
            }
            if( idx < 2 ){
                shared[ldz]+=shared[ldz+blocksize*2];   
                shared[ldz+sv]+=shared[ldz+sv+blocksize*2];   
            }
            __syncthreads();
            if( idx == 0 ) {
                dy[row+num_rows*idz*2] = 
                (shared[ldz]+shared[ldz+blocksize*1])*alpha;

                dy[row+num_rows*idz*2+num_rows] = 
                (shared[ldz+sv]+shared[ldz+sv+blocksize*1])*alpha;
            }

        }

    }
#endif
}


// SELLP SpMV kernel 3D grid
// see paper by M. KREUTZER, G. HAGER, G WELLEIN, H. FEHSKE A. BISHOP
// A UNIFIED SPARSE MATRIX DATA FORMAT 
// FOR MODERN PROCESSORS WITH WIDE SIMD UNITS
// SELLC SpMV kernel modified assigning multiple threads to each row - 2D kernel
__global__ void 
zmgesellptmv_kernel_32_3D_texb( 
    int num_rows, 
    int num_cols,
    int num_vecs,
    int blocksize,
    int T,
    magma_minproductFloatComplex alpha, 
    magma_minproductFloatComplex * dval, 
    magma_minproduct_index_t * dcolind,
    magma_minproduct_index_t * drowptr,
    cudaTextureObject_t texdx,
    magma_minproductFloatComplex * dy)
{
#if defined(PRECISION_d) && defined(TEXTURE) && (__CUDA_ARCH__ >= 300)
   // T threads assigned to each row
    int idx = threadIdx.y ;     // thread in row
    int idy = threadIdx.x;      // local row
    int idz = threadIdx.z;      // vector
    int ldx = idx * blocksize + idy;
    int ldz = idz * blocksize * T + idx * blocksize + idy;
    int bdx = blockIdx.y * gridDim.x + blockIdx.x; // global block index
    int row = bdx * blocksize + idy;  // global row index
    int sv = num_vecs/2 * blocksize * T;

    extern __shared__ magma_minproductFloatComplex shared[];


    if(row < num_rows ){
        magma_minproductFloatComplex dot1 = MAGMA_minproduct_C_MAKE(0.0, 0.0);
        magma_minproductFloatComplex dot2 = MAGMA_minproduct_C_MAKE(0.0, 0.0);
        int offset = drowptr[ bdx ];
        int block = blocksize * T; // total number of threads

        int max_ = (drowptr[ bdx+1 ]-offset)/block;  
            // number of elements each thread handles



        for ( int k = 0; k < max_ ; k++ ){
            magma_minproductFloatComplex val = 
                        dval[ offset + ldx + block*k ];
            int col = 
                    num_vecs * dcolind[ offset + ldx + block*k ] ;

            int4 v = tex1Dfetch<int4>(texdx, col/2 + idz );
            dot1 += val * __hiloint2float(v.y, v.x);
            dot2 += val * __hiloint2float(v.w, v.z);
        }
        shared[ldz]  = dot1;
        shared[ldz+sv]  = dot2;

        __syncthreads();
        if( idx < 16 ){
            shared[ldz]+=shared[ldz+blocksize*16];    
            shared[ldz+sv]+=shared[ldz+sv+blocksize*16];               
            __syncthreads();
            if( idx < 8 ){
                shared[ldz]+=shared[ldz+blocksize*8];   
                shared[ldz+sv]+=shared[ldz+sv+blocksize*8];   
            }
            if( idx < 4 ){
                shared[ldz]+=shared[ldz+blocksize*4];   
                shared[ldz+sv]+=shared[ldz+sv+blocksize*4];   
            }
            if( idx < 2 ){
                shared[ldz]+=shared[ldz+blocksize*2];   
                shared[ldz+sv]+=shared[ldz+sv+blocksize*2];   
            }
            __syncthreads();
            if( idx == 0 ) {
                dy[row+num_rows*idz*2] = 
                (shared[ldz]+shared[ldz+blocksize*1])*alpha;

                dy[row+num_rows*idz*2+num_rows] = 
                (shared[ldz+sv]+shared[ldz+sv+blocksize*1])*alpha;
            }

        }

    }
#endif
}


//*************************** end  kernels using texture  ********************//



/**
    Purpose
    -------
    
    This routine computes Y = alpha *  A^t *  X + beta * Y on the GPU.
    Input format is SELLP. Note, that the input format for X is row-major
    while the output format for Y is column major!
    
    Arguments
    ---------

    @param[in]
    transA      magma_minproduct_trans_t
                transpose A?

    @param[in]
    m           magma_minproduct_int_t
                number of rows in A

    @param[in]
    n           magma_minproduct_int_t
                number of columns in A 

    @param[in]
    num_vecs    magma_minproduct_int_t
                number of columns in X and Y

    @param[in]
    blocksize   magma_minproduct_int_t
                number of rows in one ELL-slice

    @param[in]
    slices      magma_minproduct_int_t
                number of slices in matrix

    @param[in]
    alignment   magma_minproduct_int_t
                number of threads assigned to one row

    @param[in]
    alpha       magma_minproductFloatComplex
                scalar multiplier

    @param[in]
    dval        magma_minproductFloatComplex_ptr
                array containing values of A in SELLP

    @param[in]
    dcolind     magma_minproductIndex_ptr
                columnindices of A in SELLP

    @param[in]
    drowptr     magma_minproductIndex_ptr
                rowpointer of SELLP

    @param[in]
    dx          magma_minproductFloatComplex_ptr
                input vector x

    @param[in]
    beta        magma_minproductFloatComplex
                scalar multiplier

    @param[out]
    dy          magma_minproductFloatComplex_ptr
                input/output vector y

    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_cblas
    ********************************************************************/

extern "C" magma_minproduct_int_t
magma_minproduct_cmgesellpmv(
    magma_minproduct_trans_t transA,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproduct_int_t num_vecs,
    magma_minproduct_int_t blocksize,
    magma_minproduct_int_t slices,
    magma_minproduct_int_t alignment,
    magma_minproductFloatComplex alpha,
    magma_minproductFloatComplex_ptr dval,
    magma_minproductIndex_ptr dcolind,
    magma_minproductIndex_ptr drowptr,
    magma_minproductFloatComplex_ptr dx,
    magma_minproductFloatComplex beta,
    magma_minproductFloatComplex_ptr dy,
    magma_minproduct_queue_t queue )
{
    // using a 3D thread grid for small num_vecs, a 2D grid otherwise

    int texture=0, kepler=0, precision=0;

    magma_minproduct_int_t arch = magma_minproduct_getdevice_arch();
    if ( arch > 300 )
        kepler = 1;
           
    #if defined(PRECISION_d)
        precision = 1;
    #endif

    #if defined(TEXTURE)
        texture = 1;
    #endif

    if ( (texture==1) && (precision==1) && (kepler==1) ) {

        // Create channel.
        cudaChannelFormatDesc channel_desc;
        channel_desc = cudaCreateChannelDesc(32, 32, 32, 32, 
                                        cudaChannelFormatKindSigned);

        // Create resource descriptor.
        struct cudaResourceDesc resDescdx;
        memset(&resDescdx, 0, sizeof(resDescdx));
        resDescdx.resType = cudaResourceTypeLinear;
        resDescdx.res.linear.devPtr = (void*)dx;
        resDescdx.res.linear.desc = channel_desc;
        resDescdx.res.linear.sizeInBytes = m * num_vecs * sizeof(float);

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

        if ( num_vecs%2 ==1 ) { // only multiple of 2 can be processed
            printf("error: number of vectors has to be multiple of 2.\n");
            return MAGMA_minproduct_ERR_NOT_SUPPORTED;
        }
        if ( num_vecs > 8 ) // avoid running into memory problems
            alignment = 1; 

        int num_threads = (num_vecs/2) * blocksize*alignment;   
        
        // every thread handles two vectors
        if (  num_threads > 1024 )
            printf("error: too many threads requested.\n");

        dim3 block( blocksize, alignment, num_vecs/2 );

        int dimgrid1 = sqrt(slices);
        int dimgrid2 = magma_minproduct_ceildiv( slices, dimgrid1 );

        dim3 grid( dimgrid1, dimgrid2, 1);
        int Ms = num_vecs * blocksize*alignment * sizeof( magma_minproductFloatComplex );


        if ( alignment == 1) {
            dim3 block( blocksize, num_vecs/2, 1 );
            if ( beta == MAGMA_minproduct_C_MAKE( 0.0, 0.0 ) )
            zmgesellptmv_kernel_1_3D_texb<<< grid, block, 0, queue >>>
            ( m, n, num_vecs, blocksize, alignment, alpha,
                dval, dcolind, drowptr, texdx, dy );
            else
            zmgesellptmv_kernel_1_3D_tex<<< grid, block, 0, queue >>>
            ( m, n, num_vecs, blocksize, alignment, alpha,
                dval, dcolind, drowptr, texdx, beta, dy );
        }
        else if ( alignment == 4) {
            dim3 block( blocksize, alignment, num_vecs/2 );
            if ( beta == MAGMA_minproduct_C_MAKE( 0.0, 0.0 ) )
            zmgesellptmv_kernel_4_3D_texb<<< grid, block, Ms, queue >>>
            ( m, n, num_vecs, blocksize, alignment, alpha,
                dval, dcolind, drowptr, texdx, dy );
            else
            zmgesellptmv_kernel_4_3D_tex<<< grid, block, Ms, queue >>>
            ( m, n, num_vecs, blocksize, alignment, alpha,
                dval, dcolind, drowptr, texdx, beta, dy );
        }
        else if ( alignment == 8) {
            dim3 block( blocksize, alignment, num_vecs/2 );
            if ( beta == MAGMA_minproduct_C_MAKE( 0.0, 0.0 ) )
            zmgesellptmv_kernel_8_3D_texb<<< grid, block, Ms, queue >>>
            ( m, n, num_vecs, blocksize, alignment, alpha,
                dval, dcolind, drowptr, texdx, dy );
            else
            zmgesellptmv_kernel_8_3D_tex<<< grid, block, Ms, queue >>>
            ( m, n, num_vecs, blocksize, alignment, alpha,
                dval, dcolind, drowptr, texdx, beta, dy );
        }
        else if ( alignment == 16) {
            dim3 block( blocksize, alignment, num_vecs/2 );
            if ( beta == MAGMA_minproduct_C_MAKE( 0.0, 0.0 ) )
            zmgesellptmv_kernel_16_3D_texb<<< grid, block, Ms, queue >>>
            ( m, n, num_vecs, blocksize, alignment, alpha,
                dval, dcolind, drowptr, texdx, dy );
            else
            zmgesellptmv_kernel_16_3D_tex<<< grid, block, Ms, queue >>>
            ( m, n, num_vecs, blocksize, alignment, alpha,
                dval, dcolind, drowptr, texdx, beta, dy );
        }
        else if ( alignment == 32) {
            dim3 block( blocksize, alignment, num_vecs/2 );
            if ( beta == MAGMA_minproduct_C_MAKE( 0.0, 0.0 ) )
            zmgesellptmv_kernel_32_3D_texb<<< grid, block, Ms, queue >>>
            ( m, n, num_vecs, blocksize, alignment, alpha,
                dval, dcolind, drowptr, texdx, dy );
            else
            zmgesellptmv_kernel_32_3D_tex<<< grid, block, Ms, queue >>>
            ( m, n, num_vecs, blocksize, alignment, alpha,
                dval, dcolind, drowptr, texdx, beta, dy );
        }
        else {
            printf("error: alignment %d not supported.\n", alignment);
            return MAGMA_minproduct_ERR_NOT_SUPPORTED;
        }

    } else {

        if ( num_vecs%2 ==1 ) { // only multiple of 2 can be processed
            printf("error: number of vectors has to be multiple of 2.\n");
            return MAGMA_minproduct_ERR_NOT_SUPPORTED;
        }

        if ( num_vecs > 8 ) // avoid running into memory problems
            alignment = 1;

        int num_threads = num_vecs * blocksize*alignment;

        // every thread handles two vectors
        if (  num_threads > 1024 )
            printf("error: too many threads requested.\n");

        int dimgrid1 = sqrt(slices);
        int dimgrid2 = magma_minproduct_ceildiv( slices, dimgrid1 );

        dim3 grid( dimgrid1, dimgrid2, 1);
        int Ms =  num_threads * sizeof( magma_minproductFloatComplex );

        if ( alignment == 1) {
            dim3 block( blocksize, num_vecs, 1 ); 
            zmgesellptmv_kernel_1_3D<<< grid, block, 0, queue >>>
            ( m, n, num_vecs, blocksize, alignment, alpha,
                dval, dcolind, drowptr, dx, beta, dy );
        }
        else if ( alignment == 4) {
            dim3 block( blocksize, alignment, num_vecs );
            zmgesellptmv_kernel_4_3D<<< grid, block, Ms, queue >>>
            ( m, n, num_vecs, blocksize, alignment, alpha,
                dval, dcolind, drowptr, dx, beta, dy );
        }
        else if ( alignment == 8) {
            dim3 block( blocksize, alignment, num_vecs );
            zmgesellptmv_kernel_8_3D<<< grid, block, Ms, queue >>>
            ( m, n, num_vecs, blocksize, alignment, alpha,
                dval, dcolind, drowptr, dx, beta, dy );
        }
        else if ( alignment == 16) {
            dim3 block( blocksize, alignment, num_vecs );
            zmgesellptmv_kernel_16_3D<<< grid, block, Ms, queue >>>
            ( m, n, num_vecs, blocksize, alignment, alpha,
                dval, dcolind, drowptr, dx, beta, dy );
        }
        else if ( alignment == 32) {
            dim3 block( blocksize, alignment, num_vecs );
            zmgesellptmv_kernel_32_3D<<< grid, block, Ms, queue >>>
            ( m, n, num_vecs, blocksize, alignment, alpha,
                dval, dcolind, drowptr, dx, beta, dy );
        }
        else {
            printf("error: alignment %d not supported.\n", alignment);
            return MAGMA_minproduct_ERR_NOT_SUPPORTED;
        }
    }

   return MAGMA_minproduct_SUCCESS;
}

