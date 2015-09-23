/*
    -- MAGMA_minproduct (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @precisions normal z -> s d c
       @author Hartwig Anzt
*/
#include "common_magma_minproductsparse.h"


/**
    Purpose
    -------

    Allocates memory for magma_minproduct_z_matrix and initializes it
    with the passed value.


    Arguments
    ---------

    @param[out]
    x           magma_minproduct_z_matrix*
                vector to initialize

    @param[in]
    mem_loc     magma_minproduct_location_t
                memory for vector

    @param[in]
    num_rows    magma_minproduct_int_t
                desired length of vector

    @param[in]
    values      magma_minproductDoubleComplex
                entries in vector

    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_zaux
    ********************************************************************/

extern "C" magma_minproduct_int_t
magma_minproduct_zvinit(
    magma_minproduct_z_matrix *x,
    magma_minproduct_location_t mem_loc,
    magma_minproduct_int_t num_rows,
    magma_minproduct_int_t num_cols,
    magma_minproductDoubleComplex values,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t info = 0;
    
    // set queue for old dense routines
    magma_minproduct_queue_t orig_queue=NULL;
    magma_minproductblasGetKernelStream( &orig_queue );

    x->memory_location = Magma_minproduct_CPU;
    x->num_rows = num_rows;
    x->storage_type = Magma_minproduct_DENSE;
    x->ld = num_rows;
    x->num_cols = num_cols;
    x->nnz = num_rows*num_cols;
    x->major = Magma_minproductColMajor;
    if ( mem_loc == Magma_minproduct_CPU ) {
        x->memory_location = Magma_minproduct_CPU;
        CHECK( magma_minproduct_zmalloc_cpu( &x->val, x->nnz ));
        for( magma_minproduct_int_t i=0; i<x->nnz; i++)
             x->val[i] = values;
    }
    else if ( mem_loc == Magma_minproduct_DEV ) {
        x->memory_location = Magma_minproduct_DEV;
        CHECK( magma_minproduct_zmalloc( &x->val, x->nnz ));
        magma_minproductblas_zlaset(Magma_minproductFull, x->num_rows, x->num_cols, values, values, x->val, x->num_rows);
    }
    
cleanup:
    magma_minproductblasSetKernelStream( orig_queue );
    return MAGMA_minproduct_SUCCESS;
}



   


