/*
    -- MAGMA_tally3 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @precisions normal z -> s d c
       @author Hartwig Anzt

*/
#include "common_magma_tally3sparse.h"

#define  r(i_)  (r->dval)+i_*dofs
#define  b(i_)  (b.dval)+i_*dofs

/**
    Purpose
    -------

    Computes the residual r = b-Ax for a solution approximation x.
    It returns both, the actual residual and the residual vector

    Arguments
    ---------

    @param[in]
    A           magma_tally3_z_matrix
                input matrix A

    @param[in]
    b           magma_tally3_z_matrix
                RHS b

    @param[in]
    x           magma_tally3_z_matrix
                solution approximation

    @param[in,out]
    r           magma_tally3_z_matrix*
                residual vector 
                
    @param[out]
    res         magma_tally3DoubleComplex*
                return residual 

    @param[in]
    queue       magma_tally3_queue_t
                Queue to execute in.

    @ingroup magma_tally3sparse_zaux
    ********************************************************************/

extern "C" magma_tally3_int_t
magma_tally3_zresidualvec(
    magma_tally3_z_matrix A, magma_tally3_z_matrix b, magma_tally3_z_matrix x,
    magma_tally3_z_matrix *r, double *res,
    magma_tally3_queue_t queue )
{
    magma_tally3_int_t info =0;
    
    // set queue for old dense routines
    magma_tally3_queue_t orig_queue=NULL;
    magma_tally3blasGetKernelStream( &orig_queue );

    // some useful variables
    magma_tally3DoubleComplex zero = MAGMA_tally3_Z_ZERO, one = MAGMA_tally3_Z_ONE,
                                            mone = MAGMA_tally3_Z_NEG_ONE;
    magma_tally3_int_t dofs = A.num_rows;
    
    if ( A.num_rows == b.num_rows ) {

        CHECK( magma_tally3_z_spmv( mone, A, x, zero, *r, queue ));      // r = A x
        magma_tally3_zaxpy(dofs, one, b.dval, 1, r->dval, 1);          // r = r - b
        *res =  magma_tally3_dznrm2(dofs, r->dval, 1);            // res = ||r||
        //               /magma_tally3_dznrm2(dofs, b.dval, 1);               /||b||
        //printf( "relative residual: %e\n", *res );

    } else if ((b.num_rows*b.num_cols)%A.num_rows== 0 ) {
        magma_tally3_int_t num_vecs = b.num_rows*b.num_cols/A.num_rows;

        CHECK( magma_tally3_z_spmv( mone, A, x, zero, *r, queue ));           // r = A x

        for( magma_tally3_int_t i=0; i<num_vecs; i++) {
            magma_tally3_zaxpy(dofs, one, b(i), 1, r(i), 1);   // r = r - b
            res[i] =  magma_tally3_dznrm2(dofs, r(i), 1);        // res = ||r||
        }
        //               /magma_tally3_dznrm2(dofs, b.dval, 1);               /||b||
        //printf( "relative residual: %e\n", *res );
    } else {
        printf("error: dimensions do not match.\n");
        info = MAGMA_tally3_ERR_NOT_SUPPORTED;
    }

cleanup:
    magma_tally3blasSetKernelStream( orig_queue );
    return info;
}

