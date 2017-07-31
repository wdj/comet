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


/**
    Purpose
    -------
    
    For a Block-CSR ILU factorization, this routine swaps rows in the vector *x
    according to the pivoting in *ipiv.
    
    Arguments
    ---------

    @param[in]
    r_blocks    magma_tally2_int_t
                number of blocks

    @param[in]
    size_b      magma_tally2_int_t
                blocksize in BCSR

    @param[in]
    ipiv        magma_tally2_int_t*
                array containing pivots

    @param[in]
    x           magma_tally2DoubleComplex_ptr 
                input/output vector x

    @param[in]
    queue       magma_tally2_queue_t
                Queue to execute in.

    @ingroup magma_tally2sparse_zgegpuk
    ********************************************************************/

extern "C" magma_tally2_int_t
magma_tally2_zbcsrswp(
    magma_tally2_int_t r_blocks,
    magma_tally2_int_t size_b, 
    magma_tally2Int_ptr ipiv,
    magma_tally2DoubleComplex_ptr x,
    magma_tally2_queue_t queue )
{
    const magma_tally2_int_t nrhs = 1, n = r_blocks*size_b, ione = 1, inc = 1;

   magma_tally2DoubleComplex_ptr work; 
    magma_tally2_zmalloc_cpu( &work, r_blocks*size_b );

    // first shift the pivot elements
    for( magma_tally2_int_t k=0; k<r_blocks; k++) {
            for( magma_tally2_int_t l=0; l<size_b; l++)
            ipiv[ k*size_b+l ] = ipiv[ k*size_b+l ] + k*size_b;
    }

    // now the usual pivoting
    magma_tally2_zgetmatrix(n, 1, x, n, work, n);
    lapackf77_zlaswp(&nrhs, work, &n, &ione, &n, ipiv, &inc);
    magma_tally2_zsetmatrix(n, 1, work, n, x, n);

    magma_tally2_free_cpu(work);

    return MAGMA_tally2_SUCCESS;
}



