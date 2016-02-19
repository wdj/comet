/*
    -- MAGMA_tally3 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from zresidual.cpp normal z -> s, Sun May  3 11:22:59 2015
       @author Hartwig Anzt

*/
#include "common_magma_tally3sparse.h"

#define  r(i)  r.dval+i*dofs
#define  b(i)  b.dval+i*dofs

/**
    Purpose
    -------

    Computes the residual ||b-Ax|| for a solution approximation x.

    Arguments
    ---------

    @param[in]
    A           magma_tally3_s_matrix
                input matrix A

    @param[in]
    b           magma_tally3_s_matrix
                RHS b

    @param[in]
    x           magma_tally3_s_matrix
                solution approximation

    @param[out]
    res         float*
                return residual

    @param[in]
    queue       magma_tally3_queue_t
                Queue to execute in.

    @ingroup magma_tally3sparse_saux
    ********************************************************************/

extern "C" magma_tally3_int_t
magma_tally3_sresidual(
    magma_tally3_s_matrix A, magma_tally3_s_matrix b, magma_tally3_s_matrix x,
    float *res,
    magma_tally3_queue_t queue )
{
    magma_tally3_int_t info = 0;
    
    // set queue for old dense routines
    magma_tally3_queue_t orig_queue=NULL;
    magma_tally3blasGetKernelStream( &orig_queue );

    // some useful variables
    float zero = MAGMA_tally3_S_ZERO, one = MAGMA_tally3_S_ONE,
                                            mone = MAGMA_tally3_S_NEG_ONE;
    magma_tally3_int_t dofs = A.num_rows;
    magma_tally3_int_t num_vecs = b.num_rows*b.num_cols/A.num_rows;
    
    magma_tally3_s_matrix r={Magma_tally3_CSR};
    
    if ( A.num_rows == b.num_rows ) {
        CHECK( magma_tally3_svinit( &r, Magma_tally3_DEV, A.num_rows, b.num_cols, zero, queue ));

        CHECK( magma_tally3_s_spmv( one, A, x, zero, r, queue ));           // r = A x
        magma_tally3_saxpy(dofs, mone, b.dval, 1, r.dval, 1);          // r = r - b
        *res =  magma_tally3_snrm2(dofs, r.dval, 1);            // res = ||r||
        //               /magma_tally3_snrm2(dofs, b.dval, 1);               /||b||
        //printf( "relative residual: %e\n", *res );

    } else if ((b.num_rows*b.num_cols)%A.num_rows== 0 ) {
        
        CHECK( magma_tally3_svinit( &r, Magma_tally3_DEV, b.num_rows,b.num_cols, zero, queue ));

        CHECK( magma_tally3_s_spmv( one, A, x, zero, r, queue ));           // r = A x

        for( magma_tally3_int_t i=0; i<num_vecs; i++) {
            magma_tally3_saxpy(dofs, mone, b(i), 1, r(i), 1);   // r = r - b
            res[i] =  magma_tally3_snrm2(dofs, r(i), 1);        // res = ||r||
        }
        //               /magma_tally3_snrm2(dofs, b.dval, 1);               /||b||
        //printf( "relative residual: %e\n", *res );

    } else {
        printf("error: dimensions do not match.\n");
        info = MAGMA_tally3_ERR_NOT_SUPPORTED;
    }
    
cleanup:
    magma_tally3_smfree(&r, queue );
    magma_tally3blasSetKernelStream( orig_queue );
    return info;
}

