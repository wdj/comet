/*
    -- MAGMA_tally2 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from zresidual.cpp normal z -> c, Sun May  3 11:22:59 2015
       @author Hartwig Anzt

*/
#include "common_magma_tally2sparse.h"

#define  r(i)  r.dval+i*dofs
#define  b(i)  b.dval+i*dofs

/**
    Purpose
    -------

    Computes the residual ||b-Ax|| for a solution approximation x.

    Arguments
    ---------

    @param[in]
    A           magma_tally2_c_matrix
                input matrix A

    @param[in]
    b           magma_tally2_c_matrix
                RHS b

    @param[in]
    x           magma_tally2_c_matrix
                solution approximation

    @param[out]
    res         magma_tally2FloatComplex*
                return residual

    @param[in]
    queue       magma_tally2_queue_t
                Queue to execute in.

    @ingroup magma_tally2sparse_caux
    ********************************************************************/

extern "C" magma_tally2_int_t
magma_tally2_cresidual(
    magma_tally2_c_matrix A, magma_tally2_c_matrix b, magma_tally2_c_matrix x,
    float *res,
    magma_tally2_queue_t queue )
{
    magma_tally2_int_t info = 0;
    
    // set queue for old dense routines
    magma_tally2_queue_t orig_queue=NULL;
    magma_tally2blasGetKernelStream( &orig_queue );

    // some useful variables
    magma_tally2FloatComplex zero = MAGMA_tally2_C_ZERO, one = MAGMA_tally2_C_ONE,
                                            mone = MAGMA_tally2_C_NEG_ONE;
    magma_tally2_int_t dofs = A.num_rows;
    magma_tally2_int_t num_vecs = b.num_rows*b.num_cols/A.num_rows;
    
    magma_tally2_c_matrix r={Magma_tally2_CSR};
    
    if ( A.num_rows == b.num_rows ) {
        CHECK( magma_tally2_cvinit( &r, Magma_tally2_DEV, A.num_rows, b.num_cols, zero, queue ));

        CHECK( magma_tally2_c_spmv( one, A, x, zero, r, queue ));           // r = A x
        magma_tally2_caxpy(dofs, mone, b.dval, 1, r.dval, 1);          // r = r - b
        *res =  magma_tally2_scnrm2(dofs, r.dval, 1);            // res = ||r||
        //               /magma_tally2_scnrm2(dofs, b.dval, 1);               /||b||
        //printf( "relative residual: %e\n", *res );

    } else if ((b.num_rows*b.num_cols)%A.num_rows== 0 ) {
        
        CHECK( magma_tally2_cvinit( &r, Magma_tally2_DEV, b.num_rows,b.num_cols, zero, queue ));

        CHECK( magma_tally2_c_spmv( one, A, x, zero, r, queue ));           // r = A x

        for( magma_tally2_int_t i=0; i<num_vecs; i++) {
            magma_tally2_caxpy(dofs, mone, b(i), 1, r(i), 1);   // r = r - b
            res[i] =  magma_tally2_scnrm2(dofs, r(i), 1);        // res = ||r||
        }
        //               /magma_tally2_scnrm2(dofs, b.dval, 1);               /||b||
        //printf( "relative residual: %e\n", *res );

    } else {
        printf("error: dimensions do not match.\n");
        info = MAGMA_tally2_ERR_NOT_SUPPORTED;
    }
    
cleanup:
    magma_tally2_cmfree(&r, queue );
    magma_tally2blasSetKernelStream( orig_queue );
    return info;
}

