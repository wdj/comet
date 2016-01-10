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

#define  Ablockinfo(i,j)  Ablockinfo[(i)*c_blocks   + (j)]
#define  Bblockinfo(i,j)  Bblockinfo[(i)*c_blocks   + (j)]
#define A(i,j) ((Ablockinfo(i,j)-1)*size_b*size_b)
#define B(i,j) ((Bblockinfo(i,j)-1)*size_b*size_b)

//============================================================

#define ldb m
#define lda m
#define ldc m


// every multiprocessor handles one BCSR-block
__global__ void 
zbcsrlupivloc_kernel( 
    int size_b,
    int kblocks,   
    double **A, 
    magma_tally4_int_t *ipiv)
{
    if( blockIdx.x < kblocks ) {
        if(threadIdx.x < size_b ){
            for( int i=0; i<size_b; i++){
                int dst = ipiv[i]-1;
                if( dst != i ){
                    double *A1 = A[blockIdx.x]+threadIdx.x*size_b+i;
                    double *A2 = A[blockIdx.x]+threadIdx.x*size_b+dst;
                    double tmp = *A2;
                    *A2 = *A1;
                    *A1 = tmp;
                }               
            }
            
        }
    }

}





/**
    Purpose
    -------
    
    For a Block-CSR ILU factorization, this routine updates all blocks in
    the trailing matrix.
    
    Arguments
    ---------

    @param[in]
    size_b      magma_tally4_int_t
                blocksize in BCSR
    
    @param[in]
    kblocks     magma_tally4_int_t
                number of blocks
                
    @param[in]
    dA          magma_tally4DoubleComplex_ptr *
                matrix in BCSR

    @param[in]
    ipiv        magma_tally4Int_ptr
                array containing pivots
    @param[in]
    queue       magma_tally4_queue_t
                Queue to execute in.

    @ingroup magma_tally4sparse_zgegpuk
    ********************************************************************/

extern "C" magma_tally4_int_t
magma_tally4_zbcsrlupivloc(
    magma_tally4_int_t size_b, 
    magma_tally4_int_t kblocks,
    magma_tally4DoubleComplex_ptr *dA,  
    magma_tally4Int_ptr ipiv,
    magma_tally4_queue_t queue )
{
    #if defined(PRECISION_d)
    dim3 threads( 64, 1 );

    dim3 grid(kblocks, 1, 1);
    zbcsrlupivloc_kernel<<< grid, threads, 0, queue >>>( 
                  size_b, kblocks, dA, ipiv );

#endif


    return MAGMA_tally4_SUCCESS;
}


