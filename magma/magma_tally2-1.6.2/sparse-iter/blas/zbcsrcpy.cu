/*
    -- MAGMA_tally2 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @precisions normal z -> c d s

*/

#include "common_magma_tally2.h"

#define BLOCK_SIZE 512


// every multiprocessor handles one BCSR-block to copy from A
__global__ void 
zbcsrvalcpy_kernel( 
    int size_b,
    magma_tally2_int_t num_blocks,
    magma_tally2DoubleComplex **Aval, 
    magma_tally2DoubleComplex **Bval )
{
    if(blockIdx.x*65535+blockIdx.y < num_blocks){
        magma_tally2DoubleComplex *dA = Aval[ blockIdx.x*65535+blockIdx.y ];
        magma_tally2DoubleComplex *dB = Bval[ blockIdx.x*65535+blockIdx.y ];
        int i = threadIdx.x;

        while( i<size_b*size_b ){
                dB[i] = dA[i];
                i+=BLOCK_SIZE;
        }
    }
}

// every multiprocessor handles one BCSR-block to initialize with 0
__global__ void 
zbcsrvalzro_kernel( 
    int size_b,
    magma_tally2_int_t num_blocks,
    magma_tally2DoubleComplex **Bval )
{
    if(blockIdx.x*65535+blockIdx.y < num_blocks){
        magma_tally2DoubleComplex *dB = Bval[ blockIdx.x*65535+blockIdx.y ];
        int i = threadIdx.x;
        //dB += i;

        while( i<size_b*size_b ){
                dB[i] = MAGMA_tally2_Z_MAKE(0.0, 0.0);
                i+=BLOCK_SIZE;
        }
    }

}



/**
    Purpose
    -------
    
    For a Block-CSR ILU factorization, this routine copies the filled blocks
    from the original matrix A and initializes the blocks that will later be 
    filled in the factorization process with zeros.
    
    Arguments
    ---------


    @param[in]
    size_b      magma_tally2_int_t
                blocksize in BCSR

    @param[in]
    num_blocks  magma_tally2_int_t
                number of nonzero blocks

    @param[in]
    num_zblocks magma_tally2_int_t
                number of zero-blocks (will later be filled)

    @param[in]
    Aval        magma_tally2DoubleComplex_ptr *
                pointers to the nonzero blocks in A

    @param[in]
    Bval        magma_tally2DoubleComplex_ptr *
                pointers to the nonzero blocks in B

    @param[in]
    Bval2       magma_tally2DoubleComplex_ptr *
                pointers to the zero blocks in B

    @param[in]
    queue       magma_tally2_queue_t
                Queue to execute in.

    @ingroup magma_tally2sparse_zgegpuk
    ********************************************************************/

extern "C" magma_tally2_int_t
magma_tally2_zbcsrvalcpy(
    magma_tally2_int_t size_b, 
    magma_tally2_int_t num_blocks, 
    magma_tally2_int_t num_zblocks, 
    magma_tally2DoubleComplex_ptr *Aval, 
    magma_tally2DoubleComplex_ptr *Bval,
    magma_tally2DoubleComplex_ptr *Bval2,
    magma_tally2_queue_t queue )
{
    dim3 dimBlock( BLOCK_SIZE, 1, 1 );

        // the grids are adapted to the number of nonzero/zero blocks 
        // the upper block-number the kernels can handle is 65535*65535
        int dimgrid1 = 65535;
        int dimgrid2 = magma_tally2_ceildiv( num_blocks, 65535 );
        int dimgrid3 = magma_tally2_ceildiv( num_zblocks, 65535 );
        dim3 dimGrid( dimgrid2, dimgrid1, 1 );

        zbcsrvalcpy_kernel<<<dimGrid,dimBlock, 0, queue >>>
                            ( size_b, num_blocks, Aval, Bval );

        dim3 dimGrid2( dimgrid3, dimgrid1, 1 );

        zbcsrvalzro_kernel<<<dimGrid2,dimBlock, 0, queue >>>
                            ( size_b, num_zblocks, Bval2 );

        return MAGMA_tally2_SUCCESS;
}



