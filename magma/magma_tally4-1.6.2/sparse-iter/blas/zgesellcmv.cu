/*
    -- MAGMA_tally4 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @precisions normal z -> c d s

*/
#include "common_magma_tally4sparse.h"

#define BLOCK_SIZE 512


#define PRECISION_z


// SELLC SpMV kernel
// see paper by M. KREUTZER, G. HAGER, G WELLEIN, H. FEHSKE A. BISHOP
// A UNIFIED SPARSE MATRIX DATA FORMAT 
// FOR MODERN PROCESSORS WITH WIDE SIMD UNITS
__global__ void 
zgesellcmv_kernel(   
    int num_rows, 
    int num_cols,
    int blocksize,
    magma_tally4DoubleComplex alpha, 
    magma_tally4DoubleComplex * dval, 
    magma_tally4_index_t * dcolind,
    magma_tally4_index_t * drowptr,
    magma_tally4DoubleComplex * dx,
    magma_tally4DoubleComplex beta, 
    magma_tally4DoubleComplex * dy)
{
    // threads assigned to rows
    int Idx = blockDim.x * blockIdx.x + threadIdx.x ;
    int offset = drowptr[ blockIdx.x ];
    int border = (drowptr[ blockIdx.x+1 ]-offset)/blocksize;
    if(Idx < num_rows ){
        magma_tally4DoubleComplex dot = MAGMA_tally4_Z_MAKE(0.0, 0.0);
        for ( int n = 0; n < border; n++){ 
            int col = dcolind [offset+ blocksize * n + threadIdx.x ];
            magma_tally4DoubleComplex val = dval[offset+ blocksize * n + threadIdx.x];
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
    transA      magma_tally4_trans_t
                transposition parameter for A

    @param[in]
    m           magma_tally4_int_t
                number of rows in A

    @param[in]
    n           magma_tally4_int_t
                number of columns in A 

    @param[in]
    blocksize   magma_tally4_int_t
                number of rows in one ELL-slice

    @param[in]
    slices      magma_tally4_int_t
                number of slices in matrix

    @param[in]
    alignment   magma_tally4_int_t
                number of threads assigned to one row (=1)

    @param[in]
    alpha       magma_tally4DoubleComplex
                scalar multiplier

    @param[in]
    dval        magma_tally4DoubleComplex_ptr
                array containing values of A in SELLC/P

    @param[in]
    dcolind     magma_tally4Index_ptr
                columnindices of A in SELLC/P

    @param[in]
    drowptr     magma_tally4Index_ptr
                rowpointer of SELLP

    @param[in]
    dx          magma_tally4DoubleComplex_ptr
                input vector x

    @param[in]
    beta        magma_tally4DoubleComplex
                scalar multiplier

    @param[out]
    dy          magma_tally4DoubleComplex_ptr
                input/output vector y

    @param[in]
    queue       magma_tally4_queue_t
                Queue to execute in.

    @ingroup magma_tally4sparse_zblas
    ********************************************************************/

extern "C" magma_tally4_int_t
magma_tally4_zgesellcmv(
    magma_tally4_trans_t transA,
    magma_tally4_int_t m, magma_tally4_int_t n,
    magma_tally4_int_t blocksize,
    magma_tally4_int_t slices,
    magma_tally4_int_t alignment,
    magma_tally4DoubleComplex alpha,
    magma_tally4DoubleComplex_ptr dval,
    magma_tally4Index_ptr dcolind,
    magma_tally4Index_ptr drowptr,
    magma_tally4DoubleComplex_ptr dx,
    magma_tally4DoubleComplex beta,
    magma_tally4DoubleComplex_ptr dy,
    magma_tally4_queue_t queue )
{
    // the kernel can only handle up to 65535 slices 
   // (~2M rows for blocksize 32)
   dim3 grid( slices, 1, 1);
   magma_tally4_int_t threads = blocksize;
   zgesellcmv_kernel<<< grid, threads, 0, queue >>>
   ( m, n, blocksize, alpha,
        dval, dcolind, drowptr, dx, beta, dy );

   return MAGMA_tally4_SUCCESS;
}

