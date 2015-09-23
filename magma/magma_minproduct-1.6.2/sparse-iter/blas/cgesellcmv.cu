/*
    -- MAGMA_minproduct (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from zgesellcmv.cu normal z -> c, Sun May  3 11:22:58 2015

*/
#include "common_magma_minproductsparse.h"

#define BLOCK_SIZE 512


#define PRECISION_c


// SELLC SpMV kernel
// see paper by M. KREUTZER, G. HAGER, G WELLEIN, H. FEHSKE A. BISHOP
// A UNIFIED SPARSE MATRIX DATA FORMAT 
// FOR MODERN PROCESSORS WITH WIDE SIMD UNITS
__global__ void 
cgesellcmv_kernel(   
    int num_rows, 
    int num_cols,
    int blocksize,
    magma_minproductFloatComplex alpha, 
    magma_minproductFloatComplex * dval, 
    magma_minproduct_index_t * dcolind,
    magma_minproduct_index_t * drowptr,
    magma_minproductFloatComplex * dx,
    magma_minproductFloatComplex beta, 
    magma_minproductFloatComplex * dy)
{
    // threads assigned to rows
    int Idx = blockDim.x * blockIdx.x + threadIdx.x ;
    int offset = drowptr[ blockIdx.x ];
    int border = (drowptr[ blockIdx.x+1 ]-offset)/blocksize;
    if(Idx < num_rows ){
        magma_minproductFloatComplex dot = MAGMA_minproduct_C_MAKE(0.0, 0.0);
        for ( int n = 0; n < border; n++){ 
            int col = dcolind [offset+ blocksize * n + threadIdx.x ];
            magma_minproductFloatComplex val = dval[offset+ blocksize * n + threadIdx.x];
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
    transA      magma_minproduct_trans_t
                transposition parameter for A

    @param[in]
    m           magma_minproduct_int_t
                number of rows in A

    @param[in]
    n           magma_minproduct_int_t
                number of columns in A 

    @param[in]
    blocksize   magma_minproduct_int_t
                number of rows in one ELL-slice

    @param[in]
    slices      magma_minproduct_int_t
                number of slices in matrix

    @param[in]
    alignment   magma_minproduct_int_t
                number of threads assigned to one row (=1)

    @param[in]
    alpha       magma_minproductFloatComplex
                scalar multiplier

    @param[in]
    dval        magma_minproductFloatComplex_ptr
                array containing values of A in SELLC/P

    @param[in]
    dcolind     magma_minproductIndex_ptr
                columnindices of A in SELLC/P

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
magma_minproduct_cgesellcmv(
    magma_minproduct_trans_t transA,
    magma_minproduct_int_t m, magma_minproduct_int_t n,
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
    // the kernel can only handle up to 65535 slices 
   // (~2M rows for blocksize 32)
   dim3 grid( slices, 1, 1);
   magma_minproduct_int_t threads = blocksize;
   cgesellcmv_kernel<<< grid, threads, 0, queue >>>
   ( m, n, blocksize, alpha,
        dval, dcolind, drowptr, dx, beta, dy );

   return MAGMA_minproduct_SUCCESS;
}

