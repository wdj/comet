/*
    -- MAGMA_tally2 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from zgesellcmv.cu normal z -> s, Sun May  3 11:22:58 2015

*/
#include "common_magma_tally2sparse.h"

#define BLOCK_SIZE 512


#define PRECISION_s


// SELLC SpMV kernel
// see paper by M. KREUTZER, G. HAGER, G WELLEIN, H. FEHSKE A. BISHOP
// A UNIFIED SPARSE MATRIX DATA FORMAT 
// FOR MODERN PROCESSORS WITH WIDE SIMD UNITS
__global__ void 
sgesellcmv_kernel(   
    int num_rows, 
    int num_cols,
    int blocksize,
    float alpha, 
    float * dval, 
    magma_tally2_index_t * dcolind,
    magma_tally2_index_t * drowptr,
    float * dx,
    float beta, 
    float * dy)
{
    // threads assigned to rows
    int Idx = blockDim.x * blockIdx.x + threadIdx.x ;
    int offset = drowptr[ blockIdx.x ];
    int border = (drowptr[ blockIdx.x+1 ]-offset)/blocksize;
    if(Idx < num_rows ){
        float dot = MAGMA_tally2_S_MAKE(0.0, 0.0);
        for ( int n = 0; n < border; n++){ 
            int col = dcolind [offset+ blocksize * n + threadIdx.x ];
            float val = dval[offset+ blocksize * n + threadIdx.x];
            if( val != 0){
                  dot=dot+val*dx[col];
            }
        }

        dy[ Idx ] = dot * alpha + beta * dy [ Idx ];
    }
}


/**
    Purpose
    -------
    
    This routine computes y = alpha *  A^t *  x + beta * y on the GPU.
    Input format is SELLC/SELLP.
    
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
                number of threads assigned to one row (=1)

    @param[in]
    alpha       float
                scalar multiplier

    @param[in]
    dval        magma_tally2Float_ptr
                array containing values of A in SELLC/P

    @param[in]
    dcolind     magma_tally2Index_ptr
                columnindices of A in SELLC/P

    @param[in]
    drowptr     magma_tally2Index_ptr
                rowpointer of SELLP

    @param[in]
    dx          magma_tally2Float_ptr
                input vector x

    @param[in]
    beta        float
                scalar multiplier

    @param[out]
    dy          magma_tally2Float_ptr
                input/output vector y

    @param[in]
    queue       magma_tally2_queue_t
                Queue to execute in.

    @ingroup magma_tally2sparse_sblas
    ********************************************************************/

extern "C" magma_tally2_int_t
magma_tally2_sgesellcmv(
    magma_tally2_trans_t transA,
    magma_tally2_int_t m, magma_tally2_int_t n,
    magma_tally2_int_t blocksize,
    magma_tally2_int_t slices,
    magma_tally2_int_t alignment,
    float alpha,
    magma_tally2Float_ptr dval,
    magma_tally2Index_ptr dcolind,
    magma_tally2Index_ptr drowptr,
    magma_tally2Float_ptr dx,
    float beta,
    magma_tally2Float_ptr dy,
    magma_tally2_queue_t queue )
{
    // the kernel can only handle up to 65535 slices 
   // (~2M rows for blocksize 32)
   dim3 grid( slices, 1, 1);
   magma_tally2_int_t threads = blocksize;
   sgesellcmv_kernel<<< grid, threads, 0, queue >>>
   ( m, n, blocksize, alpha,
        dval, dcolind, drowptr, dx, beta, dy );

   return MAGMA_tally2_SUCCESS;
}

