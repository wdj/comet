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


/**
    Purpose
    -------
    
    For a Block-CSR ILU factorization, this routine swaps rows in the vector *x
    according to the pivoting in *ipiv.
    
    Arguments
    ---------

    @param[in]
    r_blocks    magma_minproduct_int_t
                number of blocks

    @param[in]
    size_b      magma_minproduct_int_t
                blocksize in BCSR

    @param[in]
    ipiv        magma_minproduct_int_t*
                array containing pivots

    @param[in]
    x           magma_minproductDoubleComplex_ptr 
                input/output vector x

    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_zgegpuk
    ********************************************************************/

extern "C" magma_minproduct_int_t
magma_minproduct_zbcsrswp(
    magma_minproduct_int_t r_blocks,
    magma_minproduct_int_t size_b, 
    magma_minproductInt_ptr ipiv,
    magma_minproductDoubleComplex_ptr x,
    magma_minproduct_queue_t queue )
{
    const magma_minproduct_int_t nrhs = 1, n = r_blocks*size_b, ione = 1, inc = 1;

   magma_minproductDoubleComplex_ptr work; 
    magma_minproduct_zmalloc_cpu( &work, r_blocks*size_b );

    // first shift the pivot elements
    for( magma_minproduct_int_t k=0; k<r_blocks; k++) {
            for( magma_minproduct_int_t l=0; l<size_b; l++)
            ipiv[ k*size_b+l ] = ipiv[ k*size_b+l ] + k*size_b;
    }

    // now the usual pivoting
    magma_minproduct_zgetmatrix(n, 1, x, n, work, n);
    lapackf77_zlaswp(&nrhs, work, &n, &ione, &n, ipiv, &inc);
    magma_minproduct_zsetmatrix(n, 1, work, n, x, n);

    magma_minproduct_free_cpu(work);

    return MAGMA_minproduct_SUCCESS;
}



