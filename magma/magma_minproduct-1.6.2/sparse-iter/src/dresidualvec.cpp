/*
    -- MAGMA_minproduct (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from zresidualvec.cpp normal z -> d, Sun May  3 11:22:59 2015
       @author Hartwig Anzt

*/
#include "common_magma_minproductsparse.h"

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
    A           magma_minproduct_d_matrix
                input matrix A

    @param[in]
    b           magma_minproduct_d_matrix
                RHS b

    @param[in]
    x           magma_minproduct_d_matrix
                solution approximation

    @param[in,out]
    r           magma_minproduct_d_matrix*
                residual vector 
                
    @param[out]
    res         double*
                return residual 

    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_daux
    ********************************************************************/

extern "C" magma_minproduct_int_t
magma_minproduct_dresidualvec(
    magma_minproduct_d_matrix A, magma_minproduct_d_matrix b, magma_minproduct_d_matrix x,
    magma_minproduct_d_matrix *r, double *res,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t info =0;
    
    // set queue for old dense routines
    magma_minproduct_queue_t orig_queue=NULL;
    magma_minproductblasGetKernelStream( &orig_queue );

    // some useful variables
    double zero = MAGMA_minproduct_D_ZERO, one = MAGMA_minproduct_D_ONE,
                                            mone = MAGMA_minproduct_D_NEG_ONE;
    magma_minproduct_int_t dofs = A.num_rows;
    
    if ( A.num_rows == b.num_rows ) {

        CHECK( magma_minproduct_d_spmv( mone, A, x, zero, *r, queue ));      // r = A x
        magma_minproduct_daxpy(dofs, one, b.dval, 1, r->dval, 1);          // r = r - b
        *res =  magma_minproduct_dnrm2(dofs, r->dval, 1);            // res = ||r||
        //               /magma_minproduct_dnrm2(dofs, b.dval, 1);               /||b||
        //printf( "relative residual: %e\n", *res );

    } else if ((b.num_rows*b.num_cols)%A.num_rows== 0 ) {
        magma_minproduct_int_t num_vecs = b.num_rows*b.num_cols/A.num_rows;

        CHECK( magma_minproduct_d_spmv( mone, A, x, zero, *r, queue ));           // r = A x

        for( magma_minproduct_int_t i=0; i<num_vecs; i++) {
            magma_minproduct_daxpy(dofs, one, b(i), 1, r(i), 1);   // r = r - b
            res[i] =  magma_minproduct_dnrm2(dofs, r(i), 1);        // res = ||r||
        }
        //               /magma_minproduct_dnrm2(dofs, b.dval, 1);               /||b||
        //printf( "relative residual: %e\n", *res );
    } else {
        printf("error: dimensions do not match.\n");
        info = MAGMA_minproduct_ERR_NOT_SUPPORTED;
    }

cleanup:
    magma_minproductblasSetKernelStream( orig_queue );
    return info;
}

