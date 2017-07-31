/*
    -- MAGMA_tally2 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @precisions normal z -> s d c
       @author Hartwig Anzt
*/
#include "common_magma_tally2sparse.h"


/**
    Purpose
    -------

    Allocates memory for magma_tally2_z_matrix and initializes it
    with the passed value.


    Arguments
    ---------

    @param[out]
    x           magma_tally2_z_matrix*
                vector to initialize

    @param[in]
    mem_loc     magma_tally2_location_t
                memory for vector

    @param[in]
    num_rows    magma_tally2_int_t
                desired length of vector

    @param[in]
    values      magma_tally2DoubleComplex
                entries in vector

    @param[in]
    queue       magma_tally2_queue_t
                Queue to execute in.

    @ingroup magma_tally2sparse_zaux
    ********************************************************************/

extern "C" magma_tally2_int_t
magma_tally2_zvinit(
    magma_tally2_z_matrix *x,
    magma_tally2_location_t mem_loc,
    magma_tally2_int_t num_rows,
    magma_tally2_int_t num_cols,
    magma_tally2DoubleComplex values,
    magma_tally2_queue_t queue )
{
    magma_tally2_int_t info = 0;
    
    // set queue for old dense routines
    magma_tally2_queue_t orig_queue=NULL;
    magma_tally2blasGetKernelStream( &orig_queue );

    x->memory_location = Magma_tally2_CPU;
    x->num_rows = num_rows;
    x->storage_type = Magma_tally2_DENSE;
    x->ld = num_rows;
    x->num_cols = num_cols;
    x->nnz = num_rows*num_cols;
    x->major = Magma_tally2ColMajor;
    if ( mem_loc == Magma_tally2_CPU ) {
        x->memory_location = Magma_tally2_CPU;
        CHECK( magma_tally2_zmalloc_cpu( &x->val, x->nnz ));
        for( magma_tally2_int_t i=0; i<x->nnz; i++)
             x->val[i] = values;
    }
    else if ( mem_loc == Magma_tally2_DEV ) {
        x->memory_location = Magma_tally2_DEV;
        CHECK( magma_tally2_zmalloc( &x->val, x->nnz ));
        magma_tally2blas_zlaset(Magma_tally2Full, x->num_rows, x->num_cols, values, values, x->val, x->num_rows);
    }
    
cleanup:
    magma_tally2blasSetKernelStream( orig_queue );
    return MAGMA_tally2_SUCCESS;
}



   


