/*
    -- MAGMA_tally3 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @precisions normal z -> c d s

*/

#include "common_magma_tally3.h"


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
    uplo        magma_tally3_uplo_t
                upper/lower fill structure

    @param[in]
    r_blocks    magma_tally3_int_t
                number of blocks in row
                
    @param[in]
    c_blocks    magma_tally3_int_t
                number of blocks in column    
                
    @param[in]
    size_b      magma_tally3_int_t
                blocksize in BCSR
 
    @param[in]
    A           magma_tally3DoubleComplex_ptr 
                upper/lower factor

    @param[in]
    blockinfo   magma_tally3_int_t*
                array containing matrix information

    @param[in]
    x           magma_tally3DoubleComplex_ptr 
                input/output vector x

    @param[in]
    queue       magma_tally3_queue_t
                Queue to execute in.

    @ingroup magma_tally3sparse_zgegpuk
    ********************************************************************/

extern "C" magma_tally3_int_t
magma_tally3_zbcsrtrsv(
    magma_tally3_uplo_t uplo,
    magma_tally3_int_t r_blocks,
    magma_tally3_int_t c_blocks,
    magma_tally3_int_t size_b, 
    magma_tally3DoubleComplex_ptr A,
    magma_tally3_index_t *blockinfo,   
    magma_tally3DoubleComplex_ptr x,
    magma_tally3_queue_t queue )
{
    // set queue for old dense routines
    magma_tally3_queue_t orig_queue;
    magma_tally3blasGetKernelStream( &orig_queue );

    // some useful variables
    magma_tally3DoubleComplex one = MAGMA_tally3_Z_MAKE(1.0, 0.0);
    magma_tally3DoubleComplex mone = MAGMA_tally3_Z_MAKE(-1.0, 0.0);
    magma_tally3_int_t j,k;

    if ( uplo==Magma_tally3Lower ) { 
        // forward solve
        for( k=0; k<r_blocks; k++) {
            // do the forward triangular solve for block M(k,k): L(k,k)y = b
            magma_tally3_ztrsv(Magma_tally3Lower, Magma_tally3NoTrans, Magma_tally3Unit, size_b, A(k,k), 
                                                             size_b, x(k), 1 );

             // update for all nonzero blocks below M(k,k) 
                    // the respective values of y
            for( j=k+1; j<c_blocks; j++ ) {
                if ( (blockinfo(j,k)!=0) ) {
                    magma_tally3blas_zgemv( Magma_tally3NoTrans, size_b, size_b, 
                                     mone, A(j,k), size_b,
                                     x(k), 1, one,  x(j), 1 );

                }
            }
        }
    }
    else if ( uplo==Magma_tally3Upper ) {
        // backward solve
        for( k=r_blocks-1; k>=0; k--) {
            // do the backward triangular solve for block M(k,k): U(k,k)x = y
            magma_tally3_ztrsv(Magma_tally3Upper, Magma_tally3NoTrans, Magma_tally3NonUnit, size_b, A(k,k), 
                                                             size_b, x(k), 1 );

            // update for all nonzero blocks above M(k,k) 
                    // the respective values of y
            for( j=k-1; j>=0; j-- ) {
                if ( (blockinfo(j,k)!=0) ) {
                    magma_tally3blas_zgemv( Magma_tally3NoTrans, size_b, size_b, 
                                     mone, A(j,k), size_b,
                                     x(k), 1, one,  x(j), 1 );

                }
            }
        }
    }

    magma_tally3blasSetKernelStream( orig_queue );
    return MAGMA_tally3_SUCCESS;
}



