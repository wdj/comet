/*
    -- MAGMA_tally2 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @author Mark Gates
*/

#ifndef MAGMA_tally2SPARSE_COMMON_H
#define MAGMA_tally2SPARSE_COMMON_H

#include "common_magma_tally2.h"
#include "magma_tally2sparse.h"

#ifdef __cplusplus
extern "C" {
#endif


magma_tally2_int_t cusparse2magma_tally2_error( cusparseStatus_t status );


/**
    Macro checks the return code of a function;
    if non-zero, sets info to err, then does goto cleanup.
    err is evaluated only once.
    Assumes variable info and label cleanup exist.
    Usually, all paths (successful and error) exit through the cleanup code.
    Example:
    
        magma_tally2_int_t function()
        {
            magma_tally2_int_t info = 0;
            double *A=NULL, *B=NULL;
            CHECK( magma_tally2_malloc( &A, sizeA ));
            CHECK( magma_tally2_malloc( &B, sizeB ));
            ...
        cleanup:
            magma_tally2_free( A );
            magma_tally2_free( B );
            return info;
        }
    
    @ingroup internal
    ********************************************************************/
#define CHECK( err )             \
    do {                         \
        magma_tally2_int_t e_ = (err);  \
        if ( e_ != 0 ) {         \
            info = e_;           \
            goto cleanup;        \
        }                        \
    } while(0)


/**
    Macro checks the return code of a cusparse function;
    if non-zero, maps the cusparse error to a magma_tally2 error and sets info,
    then does goto cleanup.
    
    @see CHECK
    @ingroup internal
    ********************************************************************/
#define CHECK_CUSPARSE( err )                   \
    do {                                        \
        cusparseStatus_t e_ = (err);            \
        if ( e_ != 0 ) {                        \
            info = cusparse2magma_tally2_error( e_ );  \
            goto cleanup;                       \
        }                                       \
    } while(0)


#ifdef __cplusplus
} // extern C
#endif

#endif        //  #ifndef MAGMA_tally2SPARSE_COMMON_H
