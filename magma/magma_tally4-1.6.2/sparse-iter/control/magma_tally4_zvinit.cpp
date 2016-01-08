/*
    -- MAGMA_tally4 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @precisions normal z -> s d c
       @author Hartwig Anzt
*/
#include "common_magma_tally4sparse.h"


/**
    Purpose
    -------

    Allocates memory for magma_tally4_z_matrix and initializes it
    with the passed value.


    Arguments
    ---------

    @param[out]
    x           magma_tally4_z_matrix*
                vector to initialize

    @param[in]
    mem_loc     magma_tally4_location_t
                memory for vector

    @param[in]
    num_rows    magma_tally4_int_t
                desired length of vector

    @param[in]
    values      magma_tally4DoubleComplex
                entries in vector

    @param[in]
    queue       magma_tally4_queue_t
                Queue to execute in.

    @ingroup magma_tally4sparse_zaux
    ********************************************************************/

extern "C" magma_tally4_int_t
magma_tally4_zvinit(
    magma_tally4_z_matrix *x,
    magma_tally4_location_t mem_loc,
    magma_tally4_int_t num_rows,
    magma_tally4_int_t num_cols,
    magma_tally4DoubleComplex values,
    magma_tally4_queue_t queue )
{
    magma_tally4_int_t info = 0;
    
    // set queue for old dense routines
    magma_tally4_queue_t orig_queue=NULL;
    magma_tally4blasGetKernelStream( &orig_queue );

    x->memory_location = Magma_tally4_CPU;
    x->num_rows = num_rows;
    x->storage_type = Magma_tally4_DENSE;
    x->ld = num_rows;
    x->num_cols = num_cols;
    x->nnz = num_rows*num_cols;
    x->major = Magma_tally4ColMajor;
    if ( mem_loc == Magma_tally4_CPU ) {
        x->memory_location = Magma_tally4_CPU;
        CHECK( magma_tally4_zmalloc_cpu( &x->val, x->nnz ));
        for( magma_tally4_int_t i=0; i<x->nnz; i++)
             x->val[i] = values;
    }
    else if ( mem_loc == Magma_tally4_DEV ) {
        x->memory_location = Magma_tally4_DEV;
        CHECK( magma_tally4_zmalloc( &x->val, x->nnz ));
        magma_tally4blas_zlaset(Magma_tally4Full, x->num_rows, x->num_cols, values, values, x->val, x->num_rows);
    }
    
cleanup:
    magma_tally4blasSetKernelStream( orig_queue );
    return MAGMA_tally4_SUCCESS;
}



   


