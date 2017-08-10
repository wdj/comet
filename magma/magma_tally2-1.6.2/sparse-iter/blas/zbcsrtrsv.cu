/*
    -- MAGMA_tally2 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @precisions normal z -> c d s

*/

#include "common_magma_tally2.h"


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
    uplo        magma_tally2_uplo_t
                upper/lower fill structure

    @param[in]
    r_blocks    magma_tally2_int_t
                number of blocks in row
                
    @param[in]
    c_blocks    magma_tally2_int_t
                number of blocks in column    
                
    @param[in]
    size_b      magma_tally2_int_t
                blocksize in BCSR
 
    @param[in]
    A           magma_tally2DoubleComplex_ptr 
                upper/lower factor

    @param[in]
    blockinfo   magma_tally2_int_t*
                array containing matrix information

    @param[in]
    x           magma_tally2DoubleComplex_ptr 
                input/output vector x

    @param[in]
    queue       magma_tally2_queue_t
                Queue to execute in.

    @ingroup magma_tally2sparse_zgegpuk
    ********************************************************************/

extern "C" magma_tally2_int_t
magma_tally2_zbcsrtrsv(
    magma_tally2_uplo_t uplo,
    magma_tally2_int_t r_blocks,
    magma_tally2_int_t c_blocks,
    magma_tally2_int_t size_b, 
    magma_tally2DoubleComplex_ptr A,
    magma_tally2_index_t *blockinfo,   
    magma_tally2DoubleComplex_ptr x,
    magma_tally2_queue_t queue )
{
    // set queue for old dense routines
    magma_tally2_queue_t orig_queue;
    magma_tally2blasGetKernelStream( &orig_queue );

    // some useful variables
    magma_tally2DoubleComplex one = MAGMA_tally2_Z_MAKE(1.0, 0.0);
    magma_tally2DoubleComplex mone = MAGMA_tally2_Z_MAKE(-1.0, 0.0);
    magma_tally2_int_t j,k;

    if ( uplo==Magma_tally2Lower ) { 
        // forward solve
        for( k=0; k<r_blocks; k++) {
            // do the forward triangular solve for block M(k,k): L(k,k)y = b
            magma_tally2_ztrsv(Magma_tally2Lower, Magma_tally2NoTrans, Magma_tally2Unit, size_b, A(k,k), 
                                                             size_b, x(k), 1 );

             // update for all nonzero blocks below M(k,k) 
                    // the respective values of y
            for( j=k+1; j<c_blocks; j++ ) {
                if ( (blockinfo(j,k)!=0) ) {
                    magma_tally2blas_zgemv( Magma_tally2NoTrans, size_b, size_b, 
                                     mone, A(j,k), size_b,
                                     x(k), 1, one,  x(j), 1 );

                }
            }
        }
    }
    else if ( uplo==Magma_tally2Upper ) {
        // backward solve
        for( k=r_blocks-1; k>=0; k--) {
            // do the backward triangular solve for block M(k,k): U(k,k)x = y
            magma_tally2_ztrsv(Magma_tally2Upper, Magma_tally2NoTrans, Magma_tally2NonUnit, size_b, A(k,k), 
                                                             size_b, x(k), 1 );

            // update for all nonzero blocks above M(k,k) 
                    // the respective values of y
            for( j=k-1; j>=0; j-- ) {
                if ( (blockinfo(j,k)!=0) ) {
                    magma_tally2blas_zgemv( Magma_tally2NoTrans, size_b, size_b, 
                                     mone, A(j,k), size_b,
                                     x(k), 1, one,  x(j), 1 );

                }
            }
        }
    }

    magma_tally2blasSetKernelStream( orig_queue );
    return MAGMA_tally2_SUCCESS;
}


