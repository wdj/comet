/*
    -- MAGMA_tally4 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @precisions normal z -> c d s

*/

#include "common_magma_tally4.h"


#define  blockinfo(i,j)  blockinfo[(i)*c_blocks   + (j)]
#define A(i,j) A+((blockinfo(i,j)-1)*size_b*size_b)
#define x(i) x+(i*size_b)


/**
    Purpose
    -------
    
    For a Block-CSR ILU factorization, this routine performs the triangular 
    solves.
    
    Arguments
    ---------

    @param[in]
    uplo        magma_tally4_uplo_t
                upper/lower fill structure

    @param[in]
    r_blocks    magma_tally4_int_t
                number of blocks in row
                
    @param[in]
    c_blocks    magma_tally4_int_t
                number of blocks in column    
                
    @param[in]
    size_b      magma_tally4_int_t
                blocksize in BCSR
 
    @param[in]
    A           magma_tally4DoubleComplex_ptr 
                upper/lower factor

    @param[in]
    blockinfo   magma_tally4_int_t*
                array containing matrix information

    @param[in]
    x           magma_tally4DoubleComplex_ptr 
                input/output vector x

    @param[in]
    queue       magma_tally4_queue_t
                Queue to execute in.

    @ingroup magma_tally4sparse_zgegpuk
    ********************************************************************/

extern "C" magma_tally4_int_t
magma_tally4_zbcsrtrsv(
    magma_tally4_uplo_t uplo,
    magma_tally4_int_t r_blocks,
    magma_tally4_int_t c_blocks,
    magma_tally4_int_t size_b, 
    magma_tally4DoubleComplex_ptr A,
    magma_tally4_index_t *blockinfo,   
    magma_tally4DoubleComplex_ptr x,
    magma_tally4_queue_t queue )
{
    // set queue for old dense routines
    magma_tally4_queue_t orig_queue;
    magma_tally4blasGetKernelStream( &orig_queue );

    // some useful variables
    magma_tally4DoubleComplex one = MAGMA_tally4_Z_MAKE(1.0, 0.0);
    magma_tally4DoubleComplex mone = MAGMA_tally4_Z_MAKE(-1.0, 0.0);
    magma_tally4_int_t j,k;

    if ( uplo==Magma_tally4Lower ) { 
        // forward solve
        for( k=0; k<r_blocks; k++) {
            // do the forward triangular solve for block M(k,k): L(k,k)y = b
            magma_tally4_ztrsv(Magma_tally4Lower, Magma_tally4NoTrans, Magma_tally4Unit, size_b, A(k,k), 
                                                             size_b, x(k), 1 );

             // update for all nonzero blocks below M(k,k) 
                    // the respective values of y
            for( j=k+1; j<c_blocks; j++ ) {
                if ( (blockinfo(j,k)!=0) ) {
                    magma_tally4blas_zgemv( Magma_tally4NoTrans, size_b, size_b, 
                                     mone, A(j,k), size_b,
                                     x(k), 1, one,  x(j), 1 );

                }
            }
        }
    }
    else if ( uplo==Magma_tally4Upper ) {
        // backward solve
        for( k=r_blocks-1; k>=0; k--) {
            // do the backward triangular solve for block M(k,k): U(k,k)x = y
            magma_tally4_ztrsv(Magma_tally4Upper, Magma_tally4NoTrans, Magma_tally4NonUnit, size_b, A(k,k), 
                                                             size_b, x(k), 1 );

            // update for all nonzero blocks above M(k,k) 
                    // the respective values of y
            for( j=k-1; j>=0; j-- ) {
                if ( (blockinfo(j,k)!=0) ) {
                    magma_tally4blas_zgemv( Magma_tally4NoTrans, size_b, size_b, 
                                     mone, A(j,k), size_b,
                                     x(k), 1, one,  x(j), 1 );

                }
            }
        }
    }

    magma_tally4blasSetKernelStream( orig_queue );
    return MAGMA_tally4_SUCCESS;
}



