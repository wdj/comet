/*
    -- MAGMA_tally3 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @precisions normal z -> c d s

*/

#include "common_magma_tally3.h"

#define BLOCK_SIZE 512

#define  blockinfo(i,j)  blockinfo[(i)*c_blocks   + (j)]
#define  val(i,j) val+((blockinfo(i,j)-1)*size_b*size_b)



// every thread initializes one entry
__global__ void 
zbcsrblockinfo5_kernel( 
    magma_tally3_int_t num_blocks,
    magma_tally3DoubleComplex * address,
    magma_tally3DoubleComplex **AII )
{
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if( i < num_blocks ){
        *AII[ i ] = *address;
        if(i==0)
        printf("address: %d\n", address);
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
    lustep      magma_tally3_int_t
                lustep

    @param[in]
    num_blocks  magma_tally3_int_t
                number of nonzero blocks

    @param[in]
    c_blocks    magma_tally3_int_t
                number of column-blocks
                
    @param[in]
    size_b      magma_tally3_int_t
                blocksize
                
    @param[in]
    blockinfo   magma_tally3_int_t*
                block filled? location?

    @param[in]
    val         magma_tally3DoubleComplex*
                pointers to the nonzero blocks in A

    @param[in]
    AII         magma_tally3DoubleComplex**
                pointers to the respective nonzero blocks in B

    @param[in]
    queue       magma_tally3_queue_t
                Queue to execute in.

    @ingroup magma_tally3sparse_zgegpuk
    ********************************************************************/

extern "C" magma_tally3_int_t
magma_tally3_zbcsrblockinfo5(
    magma_tally3_int_t lustep,
    magma_tally3_int_t num_blocks, 
    magma_tally3_int_t c_blocks, 
    magma_tally3_int_t size_b,
    magma_tally3_index_t *blockinfo,
    magma_tally3DoubleComplex_ptr val,
    magma_tally3DoubleComplex_ptr *AII,
    magma_tally3_queue_t queue )
{
    dim3 dimBlock( BLOCK_SIZE, 1, 1 );

        int dimgrid = magma_tally3_ceildiv( num_blocks, BLOCK_SIZE );
        dim3 dimGrid( dimgrid, 1, 1 );


        printf("dim grid: %d x %d", dimgrid, BLOCK_SIZE);
        magma_tally3DoubleComplex **hAII;
        magma_tally3_malloc((void **)&hAII, num_blocks*sizeof(magma_tally3DoubleComplex*));

        for(int i=0; i<num_blocks; i++) {
           hAII[i] = val(lustep,lustep);
        }
        magma_tally3_setvector( num_blocks, sizeof(magma_tally3DoubleComplex*), 
                                                            hAII, 1, AII, 1 );
/*
    magma_tally3_setvector( 1, sizeof(magma_tally3DoubleComplex*), address, 1, daddress, 1 );
    zbcsrblockinfo5_kernel<<<dimGrid,dimBlock, 0, queue >>>
                        ( num_blocks, daddress, AII );

*/
        return MAGMA_tally3_SUCCESS;
}



