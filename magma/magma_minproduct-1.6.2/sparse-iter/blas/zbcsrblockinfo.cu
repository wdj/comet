/*
    -- MAGMA_minproduct (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @precisions normal z -> c d s

*/

#include "common_magma_minproduct.h"

#define BLOCK_SIZE 512

#define  blockinfo(i,j)  blockinfo[(i)*c_blocks   + (j)]
#define  val(i,j) val+((blockinfo(i,j)-1)*size_b*size_b)



// every thread initializes one entry
__global__ void 
zbcsrblockinfo5_kernel( 
    magma_minproduct_int_t num_blocks,
    magma_minproductDoubleComplex * address,
    magma_minproductDoubleComplex **AII )
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
    lustep      magma_minproduct_int_t
                lustep

    @param[in]
    num_blocks  magma_minproduct_int_t
                number of nonzero blocks

    @param[in]
    c_blocks    magma_minproduct_int_t
                number of column-blocks
                
    @param[in]
    size_b      magma_minproduct_int_t
                blocksize
                
    @param[in]
    blockinfo   magma_minproduct_int_t*
                block filled? location?

    @param[in]
    val         magma_minproductDoubleComplex*
                pointers to the nonzero blocks in A

    @param[in]
    AII         magma_minproductDoubleComplex**
                pointers to the respective nonzero blocks in B

    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_zgegpuk
    ********************************************************************/

extern "C" magma_minproduct_int_t
magma_minproduct_zbcsrblockinfo5(
    magma_minproduct_int_t lustep,
    magma_minproduct_int_t num_blocks, 
    magma_minproduct_int_t c_blocks, 
    magma_minproduct_int_t size_b,
    magma_minproduct_index_t *blockinfo,
    magma_minproductDoubleComplex_ptr val,
    magma_minproductDoubleComplex_ptr *AII,
    magma_minproduct_queue_t queue )
{
    dim3 dimBlock( BLOCK_SIZE, 1, 1 );

        int dimgrid = magma_minproduct_ceildiv( num_blocks, BLOCK_SIZE );
        dim3 dimGrid( dimgrid, 1, 1 );


        printf("dim grid: %d x %d", dimgrid, BLOCK_SIZE);
        magma_minproductDoubleComplex **hAII;
        magma_minproduct_malloc((void **)&hAII, num_blocks*sizeof(magma_minproductDoubleComplex*));

        for(int i=0; i<num_blocks; i++) {
           hAII[i] = val(lustep,lustep);
        }
        magma_minproduct_setvector( num_blocks, sizeof(magma_minproductDoubleComplex*), 
                                                            hAII, 1, AII, 1 );
/*
    magma_minproduct_setvector( 1, sizeof(magma_minproductDoubleComplex*), address, 1, daddress, 1 );
    zbcsrblockinfo5_kernel<<<dimGrid,dimBlock, 0, queue >>>
                        ( num_blocks, daddress, AII );

*/
        return MAGMA_minproduct_SUCCESS;
}



