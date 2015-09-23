/*
    -- MAGMA_minproduct (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @precisions normal z -> s d c
       @author Hartwig Anzt
*/

//  in this file, many routines are taken from
//  the IO functions provided by MatrixMarket

#include "common_magma_minproductsparse.h"


/**
    Purpose
    -------

    Passes a vector to MAGMA_minproduct.

    Arguments
    ---------

    @param[in]
    m           magma_minproduct_int_t
                number of rows

    @param[in]
    n           magma_minproduct_int_t
                number of columns

    @param[in]
    val         magma_minproductDoubleComplex*
                array containing vector entries

    @param[out]
    v           magma_minproduct_z_matrix*
                magma_minproduct vector
    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_zaux
    ********************************************************************/

extern "C"
magma_minproduct_int_t
magma_minproduct_zvset(
    magma_minproduct_int_t m, magma_minproduct_int_t n,
    magma_minproductDoubleComplex *val,
    magma_minproduct_z_matrix *v,
    magma_minproduct_queue_t queue )
{
    v->num_rows = m;
    v->num_cols = n;
    v->nnz = m*n;
    v->memory_location = Magma_minproduct_CPU;
    v->val = val;
    v->major = Magma_minproductColMajor;
    v->storage_type = Magma_minproduct_DENSE;

    return MAGMA_minproduct_SUCCESS;
}


/**
    Purpose
    -------

    Passes a MAGMA_minproduct vector back.

    Arguments
    ---------

    @param[in]
    v           magma_minproduct_z_matrix
                magma_minproduct vector

    @param[out]
    m           magma_minproduct_int_t
                number of rows

    @param[out]
    n           magma_minproduct_int_t
                number of columns

    @param[out]
    val         magma_minproductDoubleComplex*
                array containing vector entries

    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_zaux
    ********************************************************************/

extern "C"
magma_minproduct_int_t
magma_minproduct_zvget(
    magma_minproduct_z_matrix v,
    magma_minproduct_int_t *m, magma_minproduct_int_t *n,
    magma_minproductDoubleComplex **val,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_z_matrix v_CPU={Magma_minproduct_CSR};
    magma_minproduct_int_t info =0;
    
    if ( v.memory_location == Magma_minproduct_CPU ) {

        *m = v.num_rows;
        *n = v.num_cols;
        *val = v.val;
    } else {
        CHECK( magma_minproduct_zmtransfer( v, &v_CPU, v.memory_location, Magma_minproduct_CPU, queue ));
        CHECK( magma_minproduct_zvget( v_CPU, m, n, val, queue ));
    }
    
cleanup:
    magma_minproduct_zmfree( &v_CPU, queue );
    return info;
}


